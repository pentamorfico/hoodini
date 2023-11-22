from multiprocessing import Pool
import subprocess
import pandas as pd
from flagcsnap.utils import console,merge_cluster_result
import os 
from io import StringIO

class ProteinClusterer:
    def __init__(self, obj):
        self.master = obj
        for key, val in vars(obj).items():
            setattr(self, key, val)

    def run(self):

        with console.status("[bold green]Clustering neighbor proteins ...") as status:

            if self.clust_method =="diamond_deepclust":
                command = ["diamond", "deepclust", "-d",  self.output+"/results.fasta"]
                result = subprocess.run(command, check=True,capture_output=True)
                clusterdf = pd.read_csv(StringIO(str(result.stdout,"utf-8")),sep="\t",names=["clu_rep_seq","member"],header=None)
            elif self.clust_method =="deepmmseqs":
                #run mmseqs
                clusterdf = pd.read_csv(self.output+"/results.tsv",sep="\t",names=["clu_rep_seq","member"],header=None)
            elif self.clust_method =="jackhmmer":
                temp = None
                #run jackhmmer  
            elif self.clust_method=="custom":
                temp = None
                #run jackhmmer      

            results=merge_cluster_result(self.results,clusterdf)
            results  = results.dropna(subset=['sequence'])
            results = results.assign(ID = "ID=" + results["id"])
            results[["seqid","source","type","start","end","score","strand","phase","ID"]].drop_duplicates(subset=['ID']).drop_duplicates(subset=['start','end',"seqid"]).to_csv(self.output+"/results.gff",sep="\t",index=False,header=None)
            self.results=results
            console.print(f"✔️\tClustering done\n")

        


