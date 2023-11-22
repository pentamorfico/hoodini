import subprocess
import pandas as pd


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
    