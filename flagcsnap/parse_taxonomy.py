from multiprocessing import Pool
import subprocess
import pandas as pd
from flagcsnap.utils import console,merge_cluster_result
import os 
from io import StringIO
from ete3 import NCBITaxa
from itertools import combinations_with_replacement
import numpy as np
from scipy.cluster import hierarchy
from scipy.spatial import distance
import treeswift

class TaxonomyParser:
    def __init__(self, obj):
        self.master = obj
        for key, val in vars(obj).items():
            setattr(self, key, val)

    def run(self):        

        with console.status("[bold green]Parsing taxonomy ...") as status:
            results=pd.merge(self.results, self.assembly[["assembly_accession","taxid","species_taxid","organism_name","infraspecific_name"]], on='assembly_accession',how="left")
            results["species_taxid"]=results["species_taxid"].fillna("32644").astype("int")
            results["species_taxid"].dropna().drop_duplicates().to_csv(self.output+"/taxa.tab",sep="\t",index=False,header=None)

            ps = subprocess.run(['taxonkit', 'lineage', self.output+"/taxa.tab","--data-dir",os.environ.get('CONDA_PREFIX')+"/db/taxonkit"], check=True,capture_output=True)
            ps2 = subprocess.run(['taxonkit', 'reformat',"--data-dir",os.environ.get('CONDA_PREFIX')+"/db/taxonkit"],input=ps.stdout, capture_output=True)
            ps3 = subprocess.run(['awk',"-v",'OFS="\t"','"{print $1,$3}"'],input=ps2.stdout, capture_output=True)

            taxdf = pd.read_csv(StringIO(str(ps3.stdout,"utf-8")),sep="\t",names=["taxid_taxonkit","ncbi_tax","clean_tax"],header=None) 
            taxdf[["kingdom","phylum","class","order","family","genus","species"]] = taxdf["clean_tax"].str.split(';',expand=True)
            results=pd.merge(results, taxdf, left_on='species_taxid', right_on='taxid_taxonkit',how="left").drop(columns='taxid_taxonkit')
            targets = results["target_prot"].drop_duplicates().tolist()
            results[results["id"].isin(targets)][["id","sequence"]].drop_duplicates().to_fasta("id","sequence",self.output+"/target_prots.fasta")
            console.print(f"✔️\tTaxonomy done\n")

        with console.status("[bold green]Making trees ...") as status:


            if self.tree_mode=="taxonomy":
                console.print(f"🌴\tCreating taxonomic tree")
                ncbi = NCBITaxa()
                tree = ncbi.get_topology(results["species_taxid"]).write()
                tree = treeswift.read_tree_newick(tree)
                tree.write_tree_newick(self.output+"/target_prots.nwk")
                pdm = tree.distance_matrix(leaf_labels=True)
                target_taxid = results[["target_prot","species_taxid"]].drop_duplicates(subset=['target_prot']).set_index('target_prot')['species_taxid'].to_dict()
                dispd = pd.DataFrame(index=range(len(targets)), columns=range(len(targets)))
                dispd.columns = dispd.index = targets
                for i in combinations_with_replacement(targets, 2):
                    if str(target_taxid[i[0]])==str(target_taxid[i[1]]):
                        dispd.loc[i[0], i[1]] = dispd.loc[i[1], i[0]] = 0
                    else:
                        dispd.loc[i[0], i[1]] = dispd.loc[i[1], i[0]]  = pdm[str(target_taxid[i[0]])][str(target_taxid[i[1]])]
                #save distance matrix 
                dispd.to_csv(self.output+"/distance_matrix.csv")
            elif self.tree_mode=="fast_phylo":
                console.print(f"🌴\tCreating quick protein tree")
                command = ["famsa", "-dist_export", "-square_matrix", self.output+"/target_prots.fasta", self.output+"/distance_matrix.csv"]
                subprocess.run(command, check=True,capture_output=True)
            elif self.tree_mode=="input_tree":
                tree = treeswift.read_tree(self.tree_file, schema="newick")
                pdm = tree.distance_matrix(leaf_labels=True)
                dispd = pd.DataFrame(index=range(len(targets)), columns=range(len(targets)))
                dispd.columns = dispd.index = targets
                for i in combinations_with_replacement(targets, 2):
                    if i[0]==i[1]:
                        dispd.loc[i[0], i[1]] = dispd.loc[i[1], i[0]] = 0
                    else:
                        dispd.loc[i[0], i[1]] = dispd.loc[i[1], i[0]]  = pdm[i[0]][i[1]]
                dispd.to_csv(self.output+"/distance_matrix.csv")
            elif self.tree_mode=="make_tree":
                console.print(f"🌴\tCreating fast phylogenetic tree")
                command = ["famsa", self.output+"/target_prots.fasta", self.output+"/target_prots.aln"]
                subprocess.run(command, check=True,stdout=subprocess.DEVNULL)
                command = ["FastTree", self.output+"/target_prots.aln"]
                with open(self.output+"/target_prots.nwk", 'w') as output_file:
                    subprocess.run(command, check=True, stdout=output_file)
                tree = treeswift.read_tree(self.output+"/target_prots.nwk", schema="newick")
                pdm = tree.distance_matrix(leaf_labels=True)
                dispd = pd.DataFrame(index=range(len(targets)), columns=range(len(targets)))
                dispd.columns = dispd.index = targets
                for i in combinations_with_replacement(targets, 2):
                    if i[0]==i[1]:
                        dispd.loc[i[0], i[1]] = dispd.loc[i[1], i[0]] = 0
                    else:
                        dispd.loc[i[0], i[1]] = dispd.loc[i[1], i[0]]  = pdm[i[0]][i[1]]
                dispd.to_csv(self.output+"/distance_matrix.csv")

            console.print(f"✔️\tTree done\n")

        dispd = pd.read_csv(self.output+"/distance_matrix.csv",index_col=0)
        array = dispd.values
        np.fill_diagonal(array, 0)
        dispd = pd.DataFrame(array, columns=dispd.columns, index=dispd.index)
        Z = hierarchy.linkage(distance.squareform(dispd), method = 'single')
        dendro = hierarchy.dendrogram(Z, no_plot=True, count_sort = 'descending',labels=targets)
        icoord, dcoord, leaf_labels = dendro['icoord'], dendro['dcoord'], dendro["ivl"]
        den_data = pd.DataFrame(leaf_labels,columns=["leaf_labels"])
        den_data['y'] = list(np.linspace(min([num for sublist in icoord for num in sublist]), max([num for sublist in icoord for num in sublist]), len(leaf_labels)))
        den_data['x'] = [1 for y in den_data['y']] 
        den_data = pd.merge(den_data,results[["kingdom","phylum","class","order","family","genus","species","id","assembly_accession"]],left_on="leaf_labels",right_on="id",how="left").drop(columns=["id"])
        results=pd.merge(results, den_data[["y","leaf_labels"]], left_on='target_prot', right_on='leaf_labels',how="left").drop(columns='leaf_labels')
        self.results=results
        self.den_data=den_data
        self.leaf_labels=leaf_labels
        self.dendro=dendro
