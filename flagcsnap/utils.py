import xml.etree.ElementTree as ET
import pandas as pd
import os
import re
import rich_click as click
import itertools
import requests
import concurrent.futures
import numpy as np
from pandas.core.base import PandasObject
import glob
from Bio import SeqIO
import subprocess
## Classes

class IPGXMLFile:
    def __init__(self, file_path):
        """
        Initialize the IPGXMLFile object with the path to an IPG XML file and parse it.

        Args:
        file_path (str): The path to the IPG XML file.
        """
        self.file_path = file_path
        self.xml_data = ET.parse(self.file_path)

    def to_dict(self):
        """
        Converts the IPG XML element tree into a dictionary.

        Returns:
        dict: A dictionary representation of the XML data.
        """
        root = self.xml_data.getroot()
        ipg_dict = {}

        for ipg_report in root.findall('IPGReport'):
            ipg_id = ipg_report.get('ipg')
            product_acc = ipg_report.get('product_acc')

            product = ipg_report.find('Product')
            product_details = {attr: product.get(attr) for attr in product.keys()} if product is not None else {}

            protein_list = ipg_report.find('ProteinList')
            proteins = []
            if protein_list is not None:
                for protein in protein_list.findall('Protein'):
                    protein_accver = protein.get('accver')
                    protein_details = {attr: protein.get(attr) for attr in protein.keys() if attr != 'accver'}
                    
                    cds_list = protein.find('CDSList')
                    cds_details = []
                    if cds_list is not None:
                        for cds in cds_list.findall('CDS'):
                            cds_info = {attr: cds.get(attr) for attr in cds.keys()}
                            cds_details.append(cds_info)

                    protein_details['protein_accver'] = protein_accver
                    protein_details['cds_list'] = cds_details
                    proteins.append(protein_details)

            ipg_dict[ipg_id] = {
                'product_acc': product_acc,
                'product_details': product_details,
                'proteins': proteins
            }

        return ipg_dict

    def to_dataframe(self):
        """
        Converts the IPG XML element tree into a pandas DataFrame.

        Returns:
        DataFrame: A pandas DataFrame representing the IPG data.
        """
        ipg_dict = self.to_dict()
        flattened_data = []

        for ipg_id, details in ipg_dict.items():
            basic_info = {'ipg_id': ipg_id, 'product_acc': details['product_acc']}
            product_info = details['product_details']

            if details['proteins']:
                for protein in details['proteins']:
                    protein_info = {k: v for k, v in protein.items() if k != 'cds_list'}
                    if protein['cds_list']:
                        for cds in protein['cds_list']:
                            # Create separate columns for protein and nucleotide accver
                            row = {
                                **basic_info, 
                                **product_info, 
                                **protein_info, 
                                'protein_accver': protein_info['protein_accver'], 
                                'nucleotide_accver': cds['accver'],
                                **{k: v for k, v in cds.items() if k != 'accver'}
                            }
                            flattened_data.append(row)
                    else:
                        # Include protein info even if there are no CDS elements
                        row = {**basic_info, **product_info, **protein_info}
                        flattened_data.append(row)
            else:
                # Include basic and product info if there are no proteins
                row = {**basic_info, **product_info}
                flattened_data.append(row)

        return pd.DataFrame(flattened_data)

## Inputs
def validate_input_file(ctx, param, value):
    """Validate the input file and return a list of its lines."""
    if not os.path.isfile(value):
        raise click.BadParameter(f"File '{value}' does not exist.")

    try:
        with open(value, 'r') as file:
            lines = [line.strip() for line in file.readlines() if line.strip()]
            if len(lines) <= 1:
                raise click.BadParameter("File must contain multiple lines.")

            # Check for single-column format
            for line in lines:
                if ',' in line or '\t' in line:
                    raise click.BadParameter("File must be a single-column text file without delimiters like commas or tabs.")

            # Check for special characters
            special_char_pattern = re.compile('[^A-Za-z0-9\s._]+')
            for i, line in enumerate(lines):
                match = special_char_pattern.search(line)
                if match:
                    raise click.BadParameter(f"Invalid character '{match.group()}' found in line {i+1}: \"{line}\"")

            return value
    except Exception as e:
        raise click.BadParameter(f"Error reading file: {e}")
    

# Parsing data

def download_assembly_summary(file_path):
    assembly_df = None
    for i in ["refseq","genbank","refseq_historical","genbank_historical"]:
        url =   "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_{}.txt".format(i)
        command = [
        "aria2c",  # The aria2 executable name
        "--dir", file_path,  # Output directory for downloaded files
        "-x", "8",  # Number of connections per download
        "-s", "8",  # Split the file into 8 segments for parallel download
        url  # The URL of the file you want to download
        ]
        result = subprocess.run(command,capture_output=True)
        custom_headers = ["assembly_accession","bioproject","biosample","wgs_master","refseq_category","taxid","species_taxid","organism_name","infraspecific_name","isolate","version_status","assembly_level","release_type","genome_rep","seq_rel_date","asm_name","asm_submitter","gbrs_paired_asm","paired_asm_comp","ftp_path","excluded_from_refseq","relation_to_type_material","asm_not_live_date","assembly_type","group","genome_size","genome_size_ungapped","gc_percent","replicon_count","scaffold_count","contig_count","annotation_provider","annotation_name","annotation_date","total_gene_count","protein_coding_gene_count","non_coding_gene_count","pubmed_id"]
        if isinstance(assembly_df, pd.DataFrame) :
            assembly_df = pd.concat([assembly_df,pd.read_csv(file_path+"/assembly_summary_{}.txt".format(i),sep="\t",engine="pyarrow",header=None,names=custom_headers,skiprows=2)], ignore_index=True)
        else:
            assembly_df = pd.read_csv(file_path+"/assembly_summary_{}.txt".format(i),sep="\t",engine="pyarrow",header=None,names=custom_headers,skiprows=2)
    assembly_df.to_csv(file_path+"/assembly_dataframe.tsv",sep="\t",index=False)

def chunked_iterable(iterable, size):
    it = iter(iterable)
    while True:
        chunk = tuple(itertools.islice(it, size))
        if not chunk:
            break
        yield chunk

def download_file(url, index, folder):
    response = requests.get(url,timeout=25)
    filename = f"{folder}/{index}.txt"  # Assign a value based on the index
    with open(filename, "wb") as file:
        file.write(response.content)

def download_files(urls,folder,max_concurrent_downloads):
    if not os.path.exists(folder):
        os.makedirs(folder)
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_concurrent_downloads) as executor:
        for index, url in enumerate(urls):  # Enumerate through the URLs to get the index
            executor.submit(download_file, url, index,folder)

def create_ncbi_links(chunk_list,chunk_size,db,rettype,retmode,apikey):
    link_list = []
    for c in chunked_iterable(chunk_list,size=chunk_size):
        base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db="+db+"&id="
        chunk =",".join([n.strip() for n in c if isinstance(n, str)]) # Get a comma-separated list from the chunk
        end_url = "&rettype="+rettype+"&retmode="+retmode+"&api_key="+apikey
        url = base_url + chunk + end_url
        link_list.append(url)
    return link_list

def assign_rank(row):
    category = row['source']
    completeness = row['assembly_level']
    type_strain = row["strain_number_header"]
    
    if isinstance(type_strain, str):
        return 1
    if category == 'RefSeq' and completeness == 'Chromosome':    
        return 2
    elif category == 'RefSeq' and completeness == 'Complete Genome':  
        return 3
    elif category != 'RefSeq' and (completeness == 'Chromosome' or completeness == 'Complete Genome'):  
        return 4
    else: 
        return 5
    

def choose_candidates(df_data,input_file,assembly,mode):
    if mode=="best_ipg":
        filt_key="ipg_id"
    elif mode=="best_id":
        filt_key="protein_accver"
        df_data = df_data[df_data['protein_accver'].isin(input_file)]
    elif mode=="one_id":
        filt_key="protein_accver"
        df_data = df_data[df_data['protein_accver'].isin(input_file)]
    elif mode=="best_cluster":
        #TODO add file representing the cluster>85%
        filt_key="Cluster"
    merged_df = pd.merge(df_data, assembly, left_on='assembly', right_on='assembly_accession',how="left",suffixes=('', '_df2'))
    merged_df['ranked'] = merged_df.apply(assign_rank, axis=1)
    grouped_df = merged_df.groupby(filt_key).agg(lambda x: x.tolist()).reset_index()
    grouped_df['condition'] = grouped_df['ranked'].apply(lambda x: 1 if any(val <= 4 for val in x) else 0)
    merged_df = pd.merge(merged_df, grouped_df[[filt_key,"condition"]], left_on=filt_key, right_on=filt_key,how="left")
    merged_df.sort_values('ranked', inplace=True)
    if mode=="one_id":
        filtered_df=merged_df.drop_duplicates(subset=[filt_key],keep="first")
    else:
        filtered_df = pd.concat([merged_df[merged_df['condition'] == 0],merged_df[merged_df['condition'] == 1].drop_duplicates(subset=[filt_key],keep="first")])
        filtered_df.drop_duplicates(subset=[filt_key,"start","stop"], inplace=True, keep="first")
        
    return filtered_df

def unwrap_attributes(df):
    attribute_names = set()
    for attributes in df['attributes']:
        attributes_list = attributes.split(';')
        for attr in attributes_list:
            attribute_name = attr.split('=')[0]
            attribute_names.add(attribute_name)

    # Create new columns with attribute names and default values
    for attribute_name in attribute_names:
        df = df.assign(**{attribute_name: None})

    # Function to extract attribute values and assign them to columns
    def extract_attribute_values(row):
        attributes_list = row['attributes'].split(';')
        for attr in attributes_list:
            attribute_name, attribute_value = attr.split('=')
            row[attribute_name] = attribute_value
        return row

    # Apply the function to each row
    df = df.apply(lambda row: extract_attribute_values(row), axis=1)

    # Drop the original "attributes" columns
    df = df.drop(columns=['attributes'])
    return df

def read_fasta(filename):
    with open(filename, "r") as file:
        records = file.read().split(">")[1:] # skip the first empty split
        records = [record.split("\n", 1) for record in records]
        records = [(t[0].split(" ")[0], "".join(t[1].split())) for t in records]
    return pd.DataFrame(records, columns=["id", "sequence"])

def to_fasta(df, id_col, seq_col, output_file):
    with open(output_file, 'w') as f:
        for _, row in df.iterrows():
            seq_id = row[id_col]
            sequence = row[seq_col]
            
            f.write(f'>{seq_id}\n')
            f.write(f'{sequence}\n')
PandasObject.to_fasta = to_fasta

def extract_neighborhood(protein,assembly,mod,wn): 
    custom_header = ['seqid', 'source', 'type','start','end','score','strand','phase','attributes']
    if isinstance(assembly, str):
        if os.path.exists(assembly+"/genomic.gff"):
            try:
                query = "protein_id="+protein
                gff = pd.read_csv(assembly+"/genomic.gff",sep="\t",comment="#",names=custom_header,engine="c")
                start = gff[gff['attributes'].str.contains(query)]["start"].tolist()[0]
                end = gff[gff['attributes'].str.contains(query)]["end"].tolist()[0]
                strand = gff[gff['attributes'].str.contains(query)]["strand"].tolist()[0]
                target_nuc = gff[gff['attributes'].str.contains(query)]["seqid"].tolist()[0]
                if mod=="win_nts":
                    start_win = start-wn
                    end_win = end+wn
                    subgff=gff.query("seqid == @target_nuc & type =='CDS' & start>=@start_win & end<=@end_win")
                elif mod=="win_ngen":
                    subgff=gff.query("seqid == @target_nuc & type =='CDS'").reset_index(drop=True)
                    prot_index = subgff[subgff['attributes'].str.contains(query)].index.tolist()[0]
                    subgff=subgff[prot_index-wn:prot_index+wn]
                subgff=unwrap_attributes(subgff)
                faa_df = read_fasta(assembly+"/protein.faa")
                subgff = pd.merge(subgff, faa_df[["id","sequence"]], left_on='protein_id', right_on='id',how="left")
                subgff["assembly_accession"]=assembly.rsplit("/",1)[1]   
                subgff["target_prot"]=protein   
                subgff["rel_start"]=subgff["start"]-start
                subgff["rel_end"]=subgff["end"]-start
                if strand =="-":
                    indices = subgff.index[subgff["id"] == protein]
                    relend_index = subgff.columns.get_loc("rel_end")
                    relstart_index = subgff.columns.get_loc("rel_start")
                    delta = subgff.iloc[indices, relend_index].tolist()
                    neg_strand = subgff['strand']=="-"
                    upstream = subgff['rel_start']>0
                    subgff['flip_strand']=np.where(neg_strand, "+", "-")
                    subgff['flipped']="flipped"
                    new_start=np.where(upstream, -(subgff["rel_end"]-delta), (-subgff["rel_end"])+delta)
                    new_end=np.where(upstream, -(subgff["rel_start"]-delta), (-subgff["rel_start"])+delta)
                    subgff['rel_start'],subgff['rel_end']=new_start,new_end
                    subgff.iloc[indices, relstart_index]=0
                    subgff.iloc[indices, relend_index]=delta
                else:
                    subgff['flipped']=""
                    subgff['flip_strand']=subgff['strand']
                return subgff
            except:
                return None

    else:
        return None

def merge_cluster_result(result_df,cluster_df):
    value_counts = cluster_df['clu_rep_seq'].value_counts()
    cluster_df['clu_size'] = cluster_df['clu_rep_seq'].map(value_counts)
    new_df = cluster_df.drop_duplicates(subset='clu_rep_seq')
    new_df = new_df.sort_values('clu_size', ascending=False).reset_index(drop=True)
    new_df['fam_cluster'] = new_df.index.astype(str)
    new_df = new_df.loc[new_df['clu_size'] >= 2]
    merged_df = pd.merge(cluster_df, new_df[["clu_rep_seq","fam_cluster"]], left_on="clu_rep_seq", right_on="clu_rep_seq",how="left")
    results = pd.merge(result_df, merged_df[["clu_rep_seq","fam_cluster","member"]], left_on="id", right_on="member",how="left").drop('member', axis=1)
    return results

def flat(lis):
    flatList = []
    # Iterate with outer list
    for element in lis:
        if type(element) is list:
            # Check if type is list than iterate through the sublist
            for item in element:
                flatList.append(item)
        else:
            flatList.append(element)
    return flatList

def read_contig_from_multifasta(file_path, contig_name):
    contig=None
    for record in SeqIO.parse(file_path, 'fasta'):
        if record.id == contig_name:
            contig = record.seq
            break
    if contig is None:
        return None  # Contig not found
    return contig

def extract_subsequence(file_path, contig_name, start, end, strand, prot_id, mode, win,flipped):
    if mode=="win_nts":
        start_win = start-win
        end_win = end+win
    elif mode=="win_ngen":
        start_win = start - win*1000
        start_win = end + win*1000
    if start_win<=0:
        start_win=1
    file_path = glob.glob(file_path+"/*.fna")[0]
    # Extract the subsequence based on the strand information
    sequence = read_contig_from_multifasta(file_path,contig_name)
    if sequence is not None:
        if end_win>len(sequence):
            end_win=int(len(sequence))
        subsequence = sequence[int(start_win)-1:int(end_win)]  # Adjust for 0-based indexing

        return pd.DataFrame({'contig_id':contig_name, 'sequence': str(subsequence),"start_win" : start_win,"end_win" : end_win,"strand":strand,"start_prot":start,"end_prot":end,"prot_id":prot_id,"flipped":flipped},index=[0])

# Styling

from rich.console import Console
from rich.logging import RichHandler
import logging

FORMAT = "%(message)s"
logging.basicConfig(
    level=logging.INFO, format=FORMAT, datefmt="[%X]", handlers=[RichHandler()]
)
log = logging.getLogger("rich")
console = Console()


def desaturate(color_list,f,alpha):
    r,g,b = color_list[0],color_list[1],color_list[2]
    if r>1 or g>1 or b>1:
        r,g,b = r/255,g/255,b/255
    L = 0.3*r + 0.6*g + 0.1*b
    new_r = (r + f * (L - r))*255
    new_g = (g + f * (L - g))*255
    new_b = (b + f * (L - b))*255
    new_alpha = alpha*255
    return [new_r,new_g,new_b,new_alpha] 

def darken(color_list, d, alpha):
    r, g, b = color_list[0], color_list[1], color_list[2]
    if r > 1 or g > 1 or b > 1:
        r, g, b = r/255, g/255, b/255
    new_r = r * (1 - d)
    new_g = g * (1 - d)
    new_b = b * (1 - d)
    new_alpha = alpha*255
    return [new_r*255, new_g*255, new_b*255, new_alpha]