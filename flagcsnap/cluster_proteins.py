import subprocess
import os
from io import StringIO
import pandas as pd
from flagcsnap.utils import console,merge_cluster_result
from rich.progress import Progress
from flagcsnap.deep_mmseqs import run_mmseqs_clustering
from flagcsnap.deep_jackhmmer import parallel_jackhmmer


class ProteinClusterer:
    def __init__(self, obj):
        self.master = obj
        for key, val in vars(obj).items():
            setattr(self, key, val)

    def run(self):

            if self.clust_method =="diamond_deepclust":

                with console.status("[bold green]Clustering neighbor proteins with Diamond deepclust...") as status:

                    command = ["diamond", "deepclust", "-d",  self.output+"/results.fasta", "--member-cover","0.7"]
                    result = subprocess.run(command, check=True,capture_output=True)
                    clusterdf = pd.read_csv(StringIO(str(result.stdout,"utf-8")),sep="\t",names=["clu_rep_seq","member"],header=None)

            elif self.clust_method =="deepmmseqs":

                with console.status("[bold green]Clustering neighbor proteins with MMseqs2...") as status:

                    if not os.path.exists(self.output+"/tmp_mmseqs"):
                        os.makedirs(self.output+"/tmp_mmseqs")
                    run_mmseqs_clustering(self.output+"/results.fasta",self.output+"/tmp_mmseqs",max_steps=5, sensitivity=15, cluster_mode=1, cluster_steps=9, cov_mode=0, coverage=0.7, output=self.output+"/deepmmseqs_results.tsv")
                    clusterdf = pd.read_csv(self.output+"/deepmmseqs_results.tsv",sep="\t",names=["clu_rep_seq","member"],header=None)


            elif self.clust_method =="jackhmmer":

                clusterdf = parallel_jackhmmer(self.output+"/results.fasta")  




            results=merge_cluster_result(self.results,clusterdf)
            results  = results.dropna(subset=['sequence'])
            results = results.assign(ID = "ID=" + results["id"])
            results[["seqid","source","type","start","end","score","strand","phase","ID"]].drop_duplicates(subset=['ID']).drop_duplicates(subset=['start','end',"seqid"]).to_csv(self.output+"/results.gff",sep="\t",index=False,header=None)
            self.results=results
            console.print(f"✔️\tClustering done\n")

        


