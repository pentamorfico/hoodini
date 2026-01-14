import concurrent.futures
import itertools
import time
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path

import polars as pl
import requests
import taxoniq

from hoodini.utils.classes import IPGXMLFile
from hoodini.utils.logging_utils import error, info, warn


def nuc2ass(nucleotide_ids, apikey=None, temp_dir="temp", chunk_size=10, max_concurrent=9):
    """
    Fetch nucleotide summaries and link nucleotide IDs to assemblies.

    Parameters:
        nucleotide_ids (list): List of nucleotide IDs to query.
        api_key (str): NCBI API key.
        temp_dir (str, optional): Directory to store temporary XML files. Defaults to "temp".
        chunk_size (int, optional): Number of IDs to query per request. Defaults to 100.

    Returns:
        pl.DataFrame: DataFrame containing columns `AccessionVersion`, `AssemblyAccession`, `Taxid`, and `superkingdom`.
    """
    link_list = create_ncbi_links(
        chunk_list=nucleotide_ids,
        engine="efetch",
        chunk_size=chunk_size,
        db="nuccore",
        rettype="docsum",
        retmode="xml",
        apikey=apikey,
    )
    download_files(
        urls=link_list, folder=f"{temp_dir}/nucsum", max_concurrent_downloads=max_concurrent
    )

    df_nucsum = parseXML(f"{temp_dir}/nucsum", "nucsum")
    df_nucsum = df_nucsum[df_nucsum["doc_id"] != "0"]

    df_nucsum = df_nucsum[["doc_id", "AccessionVersion"]]

    link_list = create_ncbi_links(
        chunk_list=df_nucsum["AccessionVersion"],
        engine="elink",
        chunk_size=chunk_size,
        db="nuccore",
        dbto="assembly",
        retmode="xml",
        apikey=apikey,
    )
    download_files(
        urls=link_list, folder=f"{temp_dir}/elink", max_concurrent_downloads=max_concurrent
    )

    nuc2ass = parseXML(f"{temp_dir}/elink", "nuc2ass")
    df_nucsum = df_nucsum.join(
        nuc2ass[["id_list", "linked_id"]], left_on="doc_id", right_on="id_list", how="left"
    )

    assembly_ids = df_nucsum["linked_id"].drop_nulls().unique().to_list()
    link_list = create_ncbi_links(
        chunk_list=assembly_ids,
        engine="efetch",
        chunk_size=chunk_size,
        db="assembly",
        rettype="docsum",
        retmode="xml",
        apikey=apikey,
    )
    download_files(
        urls=link_list, folder=f"{temp_dir}/asssum", max_concurrent_downloads=max_concurrent
    )

    df_asssum = parseXML(f"{temp_dir}/asssum", "asssum")

    dicc_tax = {}
    for taxid in df_asssum["Taxid"].unique():
        t = taxoniq.Taxon(taxid)
        dicc_tax[taxid] = {t.rank.name: t.scientific_name for t in t.ranked_lineage}

    taxdf = pl.DataFrame(dicc_tax).T
    df_asssum = df_asssum.join(taxdf, left_on="Taxid", right_index=True, how="left")

    df_nucsum = df_nucsum.join(
        df_asssum[["uid", "AssemblyAccession", "Taxid", "superkingdom"]],
        left_on="linked_id",
        right_on="uid",
        how="left",
    )
    df_nucsum = df_nucsum.dropna(subset=["AssemblyAccession"])
    df_nucsum["AssemblyAccession"] = df_nucsum.apply(
        lambda x: (
            f"GCF_{x['AssemblyAccession'].split('_')[1]}"
            if "_" in x["AccessionVersion"]
            else f"GCA_{x['AssemblyAccession'].split('_')[1]}"
        ),
        how="horizontal",
    )

    return df_nucsum[["AccessionVersion", "AssemblyAccession", "Taxid", "superkingdom"]]


def nuc2len(nucleotide_ids, apikey, temp_dir="temp", chunk_size=100, max_concurrent=10):
    """
    Fetch nucleotide summaries and link nucleotide IDs to assemblies.

    Parameters:
        nucleotide_ids (list): List of nucleotide IDs to query.
        apikey (str): NCBI API key.
        temp_dir (str, optional): Directory to store temporary XML files. Defaults to "temp".
        chunk_size (int, optional): Number of IDs to query per request. Defaults to 100.

    Returns:
        pl.DataFrame: DataFrame containing columns `AccessionVersion`, `AssemblyAccession`, `Taxid`, and `superkingdom`.
    """
    link_list = create_ncbi_links(
        chunk_list=nucleotide_ids,
        engine="efetch",
        chunk_size=chunk_size,
        db="nuccore",
        rettype="docsum",
        retmode="xml",
        apikey=apikey,
    )
    download_files(
        urls=link_list, folder=f"{temp_dir}/nuclen", max_concurrent_downloads=max_concurrent
    )

    df_nucsum = parseXML(f"{temp_dir}/nuclen", "nucsum")
    df_nucsum = df_nucsum[df_nucsum["doc_id"] != "0"]

    df_nucsum = df_nucsum[["AccessionVersion", "Length"]]
    df_nucsum = df_nucsum.rename(
        {"AccessionVersion": "nucleotide_id", "Length": "nucleotide_length"}
    )
    return df_nucsum


def chunked_iterable(iterable, size):
    it = iter(iterable)
    while True:
        chunk = tuple(itertools.islice(it, size))
        if not chunk:
            break
        yield chunk


def download_file(url, index, folder):
    folder = Path(folder)
    filename = folder / f"{index}.txt"
    while True:
        try:
            response = requests.get(url, timeout=65)
            if response.content:
                if (
                    "error" not in response.content.decode("utf-8")
                    and "Error" not in response.content.decode("utf-8")
                    and "ERROR" not in response.content.decode("utf-8")
                ):
                    break
            else:
                warn("No response received, retrying...")
                time.sleep(5)
        except requests.exceptions.ChunkedEncodingError as e:
            warn(f"Download interrupted (ChunkedEncodingError): {e}, retrying...")
            time.sleep(5)
        except requests.RequestException as e:
            warn(f"General network error: {e}, retrying...")
            time.sleep(5)
        except Exception as e:
            error(f"Unexpected error: {e}, aborting...")
            break

    with open(filename, "wb") as file:
        file.write(response.content)
    time.sleep(2)


def download_files(urls, folder, max_concurrent_downloads):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_concurrent_downloads) as executor:
        for index, url in enumerate(urls):
            info(f"Downloading {url} to {folder / f'{index}.txt'}")
            executor.submit(download_file, url, index, folder)


def create_ncbi_links(
    chunk_list, chunk_size, db, retmode, apikey, engine=None, rettype=None, dbto=None
):
    assert engine in ["efetch", "elink"]
    link_list = []
    for c in chunked_iterable(chunk_list, size=chunk_size):
        base_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/{engine}.fcgi?"
        if engine == "efetch":
            chunk = ",".join([n.strip() for n in c if isinstance(n, str)])
            base_url += f"&db={db}&rettype={rettype}&id="
        elif engine == "elink":
            if dbto is None:
                raise ValueError("dbto parameter should be provided when using engine=elink")
            chunk = "&id=".join([n.strip() for n in c])
            base_url += f"&dbfrom={db}&db={dbto}&id="
        end_url = "&retmode=" + retmode
        if apikey:
            end_url += "&api_key=" + apikey
        url = base_url + chunk + end_url
        link_list.append(url)
    info(f"Created {len(link_list)} links for {engine} with chunk size {chunk_size}.")
    return link_list


def process_xml(file_path, mode):
    """
    Function to process a single file and return the resulting dataframe.

    Args:
    file_path (str): Path to the XML file to be processed.
    mode (str): The mode to use when processing the file ("ipg", "nucsum", etc.).

    Returns:
    DataFrame: Processed dataframe from the file.
    """
    nuc2ass_file = IPGXMLFile(file_path)
    return nuc2ass_file.to_dataframe(mode=mode)


def parseXML(folder_path, mode):
    """
    Parses XML files in a given folder based on the specified mode and returns a concatenated DataFrame.

    Args:
    folder_path (str): Path to the folder containing XML files.
    mode (str): The mode to use for parsing ("ipg", "nucsum", "asssum", etc.).

    Returns:
    DataFrame: A concatenated DataFrame with the contents of all processed XML files.
    """

    folder_path = Path(folder_path)
    all_file_paths = [p for p in folder_path.iterdir() if p.is_file()]

    final_df = pl.DataFrame()

    with ProcessPoolExecutor() as executor:
        results = executor.map(partial(process_xml, mode=mode), all_file_paths)

        for df in results:
            final_df = pl.concat([final_df, df], how="vertical") if final_df.height > 0 else df

    return final_df
