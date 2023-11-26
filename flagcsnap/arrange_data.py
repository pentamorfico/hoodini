from multiprocessing import Pool
import subprocess
import pandas as pd
from flagcsnap.utils import console,flat,desaturate
import os 
from io import StringIO
import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib
import random


class Arranger:
    def __init__(self, obj):
        self.master = obj
        for key, val in vars(obj).items():
            setattr(self, key, val)

    def run(self):       

        with console.status("[bold green]Structuring data for plotting ...") as status:
            height_factor=20
            den_data = self.den_data
            results = self.results
            y_range = [0, max(den_data['y'])+min(den_data['y'])+(den_data['y'][1]-den_data['y'][0])]
            if len(self.leaf_labels) < 5:
                height = int(len(self.leaf_labels)*height_factor)*3
            elif len(self.leaf_labels) < 10:
                height = int(len(self.leaf_labels)*height_factor)*2
            else:
                height = int(len(self.leaf_labels)*height_factor)

            # Find the step size for the y axis, that would be den_data['y'][n] - den_data['y'][n-1] for the firt n in which the difference is not 0
            ## TODO: see why there are duplicate values

            for i in range(1, len(den_data['y'])):
                if den_data['y'][i] - den_data['y'][i-1] != 0:
                    y_step = den_data['y'][i] - den_data['y'][i-1]
                    break

            y_half_height = y_step/6
            neg_strand = results['flip_strand']=="-"
            results["dx"] = results["rel_end"]-results["rel_start"]
            gene_x_tail = np.where(neg_strand, results['rel_end'], results['rel_start'])
            gene_dx = np.where(neg_strand, -1*results['dx'], results['dx'])
            gene_x_head = gene_x_tail + gene_dx
            gene_x_head_start = np.where(neg_strand, gene_x_head+100, gene_x_head-100)
            text_x = np.where(neg_strand,gene_x_tail - (gene_x_tail-gene_x_head_start)/2, gene_x_tail + (gene_x_head_start-gene_x_tail)/2)
            results['text_x']=text_x
            results['text_y']=results['y']+y_half_height
            results['xs']=[list(n) for n in zip(list(gene_x_tail.tolist())[:],list(gene_x_tail.tolist())[:],gene_x_head_start.tolist()[:],gene_x_head.tolist()[:],gene_x_head_start.tolist()[:])]
            def divide_by_10(lst):
                return [x / 100 for x in lst]

            # Apply the function to the column
            results['xs'] = results['xs'].apply(divide_by_10)
            results['ys']=[list(n) for n in zip(results['y']-y_half_height, results['y']+y_half_height, results['y']+y_half_height, results['y'], results['y']-y_half_height)]
            results["coordinates"]=[[n] for n in [list(list(x) for x in zip(n[0],n[1])) for n in list(zip(results["xs"],results["ys"]))]]

            min_x = min([n[0] for n in flat(flat(results["coordinates"]))])
            min_y = min([n[1] for n in flat(flat(results["coordinates"]))])
            max_x = max([n[0] for n in flat(flat(results["coordinates"]))])
            max_text_len = max([len(n) for n in den_data["species"]])*4
            path=[lst for lst in [[list(n) for n in zip([-n*math.sqrt(len(self.leaf_labels))+min_x-max_text_len for n in self.dendro["dcoord"][i]],[n for n in self.dendro["icoord"][i]])] for i,n in enumerate(self.dendro["icoord"])]]
            dendrogram = pd.DataFrame({'path': path,"color":[[0,0,0,100] for _ in path]})
            den_data['coordinates'] = [[x+min_x-max_text_len,y] for x,y in zip(den_data['x'],den_data['y'])]
            den_data["start_line"]=[[min_x,y] for y in den_data['y']]
            den_data["end_line"]=[[max_x,y] for y in den_data['y']]
            den_data = den_data.loc[den_data.astype(str).drop_duplicates().index]
            
            families = results["fam_cluster"].dropna().unique().tolist()
            colors_rgb = [plt.cm.rainbow(random.random()) for _ in families]
            colors_hex = [matplotlib.colors.to_hex(color) for color in colors_rgb]
            colors_dic = {num:desaturate([n for n in color],0.4,1) for num,color in zip(families,colors_rgb)}
            results["fillcolor"]=results["fam_cluster"].map(colors_dic).apply(lambda d: d if isinstance(d, list) else [230, 230, 230,255])
            results=results.drop_duplicates(subset=['start','end',"seqid","target_prot"])
            results[["id","target_prot","seqid","fam_cluster","start","end","rel_start","rel_end","strand","locus_tag","product","assembly_accession","species_taxid","kingdom","phylum","class","order","family","genus","species","sequence"]].drop_duplicates(subset=['start','end',"seqid"]).to_csv(self.output+"/results.txt",index=False,header=False)
            
            self.results=results
            self.den_data=den_data
            self.dendrogram=dendrogram