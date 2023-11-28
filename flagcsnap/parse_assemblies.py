from multiprocessing import Pool
import subprocess
import pandas as pd
from flagcsnap.utils import console,extract_neighborhood

class AssemblyParser:
    def __init__(self, obj):
        self.master = obj
        for key, val in vars(obj).items():
            setattr(self, key, val)

    def run(self):

        with console.status("[bold green]Downloading data from NCBI ...") as status:

            if not self.assembly_folder:
                if self.ncrna or self.cctyper:
                    formatfile="protein,gff3,genome"
                else:
                    formatfile="protein,gff3"
                try:
                    command = ["datasets","download", "genome", "accession", "--inputfile", self.output+"/assembly_list.txt", 
                                "--dehydrated", "--filename", self.output+"/dehydrated.zip", "--api-key", self.apikey, 
                                "--include", formatfile, "--annotated"]
                    subprocess.run(command, check=True)
                    command = ["unzip",self.output+"/dehydrated.zip","-d",self.output+"/assembly_folder"]
                    subprocess.run(command, check=True)
                    command = ["datasets", "rehydrate", "--directory",  self.output+"/assembly_folder", "--api-key", self.apikey, "--max-workers", str(self.max_concurrent_downloads)]
                    subprocess.run(command, check=True)
                    assembly_folder = self.output+"/assembly_folder/ncbi_dataset/data/"
                    self.assembly_folder = assembly_folder
                except:
                    pass
            else:
                assembly_folder = self.assembly_folder

            console.print(f"✔️ \tDownloaded assemblies")

        
        with console.status("[bold green]Extracting neighborhoods from asssemblies ...") as status:

            pool = Pool(self.num_threads) # number of cores you want to use
            file_list = [[n,i,self.mod,self.wn]for n,i in zip(self.dropped['protein_accver'],assembly_folder+self.dropped['assembly'])] ##FisX THIS, weord way to pass arguments to a function
            df_list = pool.starmap(extract_neighborhood, file_list) #creates a list of the loaded df's
            results = pd.concat(df_list) # concatenates all the df's into a single df
            self.results = results

            results[["id","sequence"]].dropna().drop_duplicates(subset=["id"]).to_fasta("id","sequence", self.output+"/results.fasta")

            console.print(f"✔️ \t{len([n for n in df_list if n is not None])} out of {len(df_list)} assemblies were downloaded and parsed\n")














