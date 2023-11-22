from flagcsnap.utils import IPGXMLFile,create_ncbi_links, download_files, choose_candidates, console
import pandas as pd 
import os
import subprocess
from multiprocessing import Pool
import concurrent.futures

class IPGParser:
    def __init__(self, obj):
        self.master = obj
        for key, val in vars(obj).items():
            setattr(self, key, val)

    def run(self):
    
        with console.status("[bold green]Fetching data from NCBI ...") as status:
            if len(self.input_list)>0:
                link_list = create_ncbi_links(self.input_list,100,"ipg","ipg","xml",self.apikey)
                download_files(link_list,self.output+"/ipg",self.max_concurrent_downloads)

                #This can be done in parallel??
                final_df = False
                for filename in os.listdir(self.output+"/ipg/"):
                    file_path = os.path.join(self.output+"/ipg", filename)
                    if os.path.isfile(file_path):
                        # Process the file
                        ipg_file = IPGXMLFile(file_path)
                        ipg_dict = ipg_file.to_dict()
                        ipg_df = ipg_file.to_dataframe()
                        if isinstance(final_df, pd.DataFrame):
                            final_df = pd.concat([final_df, ipg_df], ignore_index=True)
                        else:
                            final_df = ipg_df
                #Filter out those that are not Bacteria or Viruses (kingdom column)
                df = final_df[final_df["kingdom"].isin(["Bacteria","Viruses"])]
                filt_ipg_df = pd.merge(df, self.assembly[["assembly_accession","species_taxid","organism_name","infraspecific_name"]], left_on='assembly',right_on="assembly_accession",how="left")
                filt_ipg_df["temp"] = filt_ipg_df["assembly_accession"].str.split(".").str[0]
                filt_ipg_df = pd.merge(filt_ipg_df, self.type_strains, left_on='temp',right_on="type_strain_assembly",how="left").drop(columns=["type_strain_assembly","temp"])
                filt_ipg_df = choose_candidates(filt_ipg_df,self.input_list,self.assembly,mode=self.cand_mode)
                filt_ipg_df['condition'] = filt_ipg_df.duplicated(subset='protein_accver', keep=False).astype(int)
                nucs=list(set(list(filt_ipg_df[filt_ipg_df['condition']==1]["nucleotide_accver"])))
                if len(nucs)>=1:
                    link_list = create_ncbi_links(nucs,200,"nuccore","","xml",self.apikey)
                    download_files(link_list,self.output+"/docsum",self.max_concurrent_downloads)
                    final_df = False
                    for filename in os.listdir(self.output+"/docsum"):
                        file_path = os.path.join(self.output+"/docsum", filename)
                        if os.path.isfile(file_path):
                            # Process the file
                            try:
                                docsum_df = pd.read_xml(open(file_path,"r"))
                                if isinstance(final_df, pd.DataFrame):
                                    final_df = pd.concat([final_df, docsum_df], ignore_index=True)
                                else:
                                    final_df = docsum_df
                            except:
                                continue
                    ipg_len_df = pd.merge(filt_ipg_df, final_df[["GBSeq_accession-version","GBSeq_length"]], left_on='nucleotide_accver', right_on='GBSeq_accession-version',how="left")
                    ipg_len_df.sort_values('GBSeq_length', inplace=True,ascending=False)
                    filt_ipg_df=pd.concat([ipg_len_df[ipg_len_df['condition'] == 0],ipg_len_df[ipg_len_df['condition'] == 1].drop_duplicates(subset=["protein_accver"],keep="first")])

                dropped=filt_ipg_df.dropna(subset=["protein_accver","assembly"])
                self.dropped = dropped
                filt_ipg_df["assembly"].dropna().to_csv(self.output+"/assembly_list.txt",index=False,header=False)
                self.ipg_df = filt_ipg_df
                self.valid_ipg = dropped["protein_accver"].unique()

                console.print(f"✔️ \t{len(self.valid_ipg)} out of {len(self.input_list)} IDs are valid and ready to be downloaded\n")

                
