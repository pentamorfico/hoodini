from multiprocessing import Pool
import subprocess
import pandas as pd
from flagcsnap.utils import console, flat, desaturate, darken
import os
from io import StringIO
import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib
import random
import pandas as pd

class Arranger:
    def __init__(self, obj):
        self.master = obj
        for key, val in vars(obj).items():
            setattr(self, key, val)

    def run(self):
        with console.status("[bold green]Structuring data for plotting ...") as status:

            if self.domains:

                start_domains = np.where(self.domains_data["flip_strand"]=="+", self.domains_data["start"]*3 + self.domains_data["rel_start"], self.domains_data["rel_end"] - self.domains_data["end"]*3)
                end_domains = np.where(self.domains_data["flip_strand"]=="+", self.domains_data["end"]*3 + self.domains_data["rel_start"], self.domains_data["rel_end"] - self.domains_data["start"]*3)
                self.domains_data["start"] = start_domains
                self.domains_data["end"] = end_domains
                temp_evalue = self.domains_data["e_value"].copy()
                self.domains_data['e_value'] = self.domains_data['e_value'].apply(lambda x: 1e-100 if float(x) < 1e-100 else float(x))
                self.domains_data['log_values'] = np.log10(self.domains_data['e_value'])

                #change the 0 to 1e-200
                norm = matplotlib.colors.LogNorm(vmin=self.domains_data['e_value'].min(), vmax=self.domains_data['e_value'].max())
                cmap = plt.cm.plasma
                self.domains_data['colors'] = [list(cmap(norm(value))) for value in self.domains_data['e_value']]
                self.domains_data['colors'] = self.domains_data['colors'].apply(lambda x: list(int(255*i) for i in x))
                ## set the alpha (the last element of the list) to 100, no matter the value that was there before
                self.domains_data['colors'] = self.domains_data['colors'].apply(lambda x: x[:-1]+[100])
                self.domains_data['e_value'] = temp_evalue.astype(str)

            height_factor = 20
            y_range = [0, max(self.den_data['y'])+min(self.den_data['y']) +
                       (self.den_data['y'][1]-self.den_data['y'][0])]
            if len(self.leaf_labels) < 5:
                height = int(len(self.leaf_labels)*height_factor)*3
            elif len(self.leaf_labels) < 10:
                height = int(len(self.leaf_labels)*height_factor)*2
            else:
                height = int(len(self.leaf_labels)*height_factor)

            # Find the step size for the y axis, that would be den_data['y'][n] - den_data['y'][n-1] for the firt n in which the difference is not 0
            for i in range(1, len(self.den_data['y'])):
                if self.den_data['y'][i] - self.den_data['y'][i-1] != 0:
                    y_step = self.den_data['y'][i] - self.den_data['y'][i-1]
                    break
            y_half_height = y_step/6
            y_quarter_height = y_step/12


            neg_strand = self.results['flip_strand'] == "-"
            self.results["dx"] = self.results["rel_end"] - \
                self.results["rel_start"]
            gene_x_tail = np.where(
                neg_strand, self.results['rel_end'], self.results['rel_start'])
            gene_dx = np.where(
                neg_strand, -1*self.results['dx'], self.results['dx'])
            gene_x_head = gene_x_tail + gene_dx
            gene_x_head_start = np.where(
                neg_strand, gene_x_head+100, gene_x_head-150)
            
            text_x = np.where(neg_strand, gene_x_tail - (gene_x_tail-gene_x_head_start) /
                              2, gene_x_tail + (gene_x_head_start-gene_x_tail)/2)
            
            self.results['text_x'] = text_x
            self.results['text_y'] = self.results['y']+y_half_height
            self.results['xs'] = [list(n) for n in zip(list(gene_x_tail.tolist())[:], list(gene_x_tail.tolist())[
                :], gene_x_head_start.tolist()[:], gene_x_head.tolist()[:], gene_x_head_start.tolist()[:])]

            def divide_by_100(lst):
                return [x / 100 for x in lst]

            # Apply the function to the column
            self.results['xs'] = self.results['xs'].apply(divide_by_100)
            self.results["text_x"] = self.results['text_x']/100
            self.results['ys'] = [list(n) for n in zip(self.results['y']-y_half_height, self.results['y'] +
                                y_half_height, self.results['y']+y_half_height, self.results['y'], self.results['y']-y_half_height)]
            self.results["text_coordinates"] = [[x, y] for x, y in zip(
                self.results['text_x'], self.results['text_y'])]
            self.results["coordinates"] = [[n] for n in [list(list(x) for x in zip(
                n[0], n[1])) for n in list(zip(self.results["xs"], self.results["ys"]))]]

            if self.domains:
                self.domains_data['xs'] = [list(n) for n in zip(list(self.domains_data["start"].tolist())[:], list(self.domains_data["start"].tolist())[
                    :], self.domains_data["end"].tolist()[:], self.domains_data["end"].tolist()[:])]
                self.domains_data['xs'] = self.domains_data['xs'].apply(divide_by_100)
                self.domains_data['ys'] = [list(n) for n in zip(self.domains_data['y']-y_half_height-y_quarter_height*self.domains_data['y_pos'], 
                                                        (self.domains_data['y']-y_half_height-y_quarter_height*self.domains_data['y_pos'])-y_quarter_height,
                                                        (self.domains_data['y']-y_half_height-y_quarter_height*self.domains_data['y_pos'])-y_quarter_height,
                                                        self.domains_data['y']-y_half_height-y_quarter_height*self.domains_data['y_pos'])]
                self.domains_data["coordinates"] = [[n] for n in [list(list(x) for x in zip(
                    n[0], n[1])) for n in list(zip(self.domains_data["xs"], self.domains_data["ys"]))]]

            min_x = min([n[0]
                        for n in flat(flat(self.results["coordinates"]))])
            min_y = min([n[1]
                        for n in flat(flat(self.results["coordinates"]))])
            max_y = max([n[1]
                        for n in flat(flat(self.results["coordinates"]))])
            max_x = max([n[0]
                        for n in flat(flat(self.results["coordinates"]))])
            max_text_len = max([len(n) for n in self.den_data["species"]])*4
            path = [lst for lst in [[list(n) for n in zip([-n*math.sqrt(len(self.leaf_labels))+min_x-max_text_len for n in self.dendro["dcoord"][i]], [
                n for n in self.dendro["icoord"][i]])] for i, n in enumerate(self.dendro["icoord"])]]
            dendrogram = pd.DataFrame(
                {'path': path, "color": [[0, 0, 0, 100] for _ in path]})
            self.den_data['coordinates'] = [[x+min_x-max_text_len, y]
                                            for x, y in zip(self.den_data['x'], self.den_data['y'])]
            self.den_data["start_line"] = [[min_x, y]
                                           for y in self.den_data['y']]
            self.den_data["end_line"] = [[max_x, y]
                                         for y in self.den_data['y']]
            
           
            min_bp = self.results["rel_start"].min()
            max_bp = self.results["rel_end"].max()
            min_thousand = math.floor(abs(min_bp) / 1000) * 1000
            max_thousand = math.floor(max_bp / 1000) * 1000
            baseline_y = max_y + y_step
            lines_start = [[x/100,baseline_y-y_quarter_height] for x in range(-min_thousand, max_thousand + 1000, 1000)]
            lines_start.append([min_x, baseline_y])            
            lines_end = [[x/100,baseline_y+y_quarter_height] for x in range(-min_thousand, max_thousand + 1000, 1000)]
            lines_end.append([max_x, baseline_y])
            ticks = pd.DataFrame({'start': lines_start, 'end': lines_end})
            tick_text_coords = [[x/100, baseline_y+y_quarter_height*2] for x in range(-min_thousand, max_thousand + 1000, 1000)]
            tick_text = pd.DataFrame({'text': [str(x) for x in range(-min_thousand, max_thousand + 1000, 1000)], 'coordinates': tick_text_coords})
            self.ticks = ticks
            self.tick_text = tick_text



            


            self.den_data = self.den_data.loc[self.den_data.astype(
                str).drop_duplicates().index]
            
            families = self.results["fam_cluster"].dropna().unique().tolist()

            #compute the prevalence percentaje of each family. For that, count the unique assembly_accession for each family and divide by the total number of unique assembly_accession
            prevalence = self.results.groupby("fam_cluster")["assembly_accession"].nunique().reset_index()
            prevalence["prevalence"] = (prevalence["assembly_accession"]/self.results["assembly_accession"].nunique()).round(2)
            prevalence = prevalence.set_index("fam_cluster")["prevalence"].to_dict()
            self.results["prevalence"] = self.results["fam_cluster"].map(prevalence)

            colors_rgb = [plt.cm.gist_ncar(random.random()) for _ in families]


            colors_dic = {num: desaturate([n for n in color], 1-prevalence[num]*0.6, 1) if prevalence[num] >= self.min_prevalence else [230, 230, 230, 255] for num, color in zip(families, colors_rgb)}
            self.results["fillcolor"] = self.results["fam_cluster"].map(colors_dic).apply(
                lambda d: d if isinstance(d, list) else [230, 230, 230, 255])
            self.results['linecolor'] = self.results['fillcolor'].apply(
                lambda x: darken(x, 0.5, 1))

            self.results = self.results.drop_duplicates(
                subset=['start', 'end', "seqid", "target_prot"])
            self.results[["id", "target_prot", "seqid", "fam_cluster", "start", "end", "rel_start", "rel_end", "strand", "flip_strand", "locus_tag", "product", "assembly_accession", "species_taxid", "kingdom",
                          "phylum", "class", "order", "family", "genus", "species", "sequence"]].drop_duplicates(subset=['start', 'end', "seqid"]).to_csv(self.output+"/results.txt", index=False)

            self.dendrogram = dendrogram
