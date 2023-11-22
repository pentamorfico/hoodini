from multiprocessing import Pool
import subprocess
import pandas as pd
from utils import console,flat
import os 
from io import StringIO
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import random
from utils import desaturate,darken,to_fasta
import concurrent.futures
from utils import extract_subsequence
from Bio import SeqIO
from ast import literal_eval

class ExtraAnnotation:
    def __init__(self, obj):
        self.master = obj
        for key, val in vars(obj).items():
            setattr(self, key, val)

    def run(self):       

        results = self.results

        if self.padloc:
            with console.status("[bold green]Running PADLOC...") as status:
                console.print("🛡️\tRunning Padloc")

                try:
                    if not os.path.exists(self.output+'/padloc'):
                        os.makedirs(self.output+'/padloc')
                    command = ["padloc", "--faa", self.output+"/results.fasta", "--gff",self.output+"/results.gff","--cpu",str(self.num_threads),"--outdir", self.output+"/padloc"]
                    subprocess.run(command, check=True)

                    padloc_df = pd.read_csv(self.output+"/padloc/results.fasta_padloc.csv")
                    padloc_df= padloc_df.rename(columns={'system': 'padloc_system',"protein.name":"padloc_gene"})
                    results=pd.merge(results, padloc_df[["padloc_system","padloc_gene","target.description"]], left_on='id', right_on='target.description',how="left").drop(columns='target.description')
                except:
                    console.print("☹️\tPADLOC failed")
        if self.deffinder:
            with console.status("[bold green]Running DefenseFinder...") as status:
                console.print("🛡️\tRunning DefenseFinder")

                try:   
                    if not os.path.exists(self.output+"/defense_finder"):
                        os.makedirs(self.output+"/defense_finder")
                    deffinder_df = results[["locus_tag","sequence"]]
                    deffinder_df.dropna().drop_duplicates().to_fasta("locus_tag","sequence",self.output+"/defense_finder/proteome.fasta")
                    command = ["defense-finder", "run", self.output+"/defense_finder/proteome.fasta", "--db-type","gembase","-o",self.output+"/defense_finder"]
                    subprocess.run(command, check=True,stdout=subprocess.DEVNULL)
                    deffinder_df = pd.read_csv(self.output+"/defense_finder/defense_finder_genes.tsv",sep="\t")
                    deffinder_df.rename(columns={'gene_name': 'deffinder_system'})
                    deffinder_df[['deffinder_system', 'deffinder_gene']] = deffinder_df['gene_name'].str.split('__', expand=True)
                    results=pd.merge(results, deffinder_df[["deffinder_system","deffinder_gene","hit_id"]], left_on='locus_tag', right_on='hit_id',how="left").drop(columns='hit_id')
                except:
                    print("Defense finder failed")
                    results["deffinder_system"] = np.nan
                    results["deffinder_gene"] = np.nan


        with console.status("[bold green]Extracting nucleotide fasta...") as status:
            if self.ncrna or self.cctyper:
                if not os.path.exists(self.output+'/neighborhood'):
                    os.makedirs(self.output+'/neighborhood')
                neighs = None
                executor = concurrent.futures.ThreadPoolExecutor(max_workers=self.num_threads)
                with executor:
                    futures = []
                    dropped = pd.merge(self.dropped,
                    results[["assembly_accession","seqid","flipped"]],left_on="assembly",right_on="assembly_accession",suffixes=('', '_df2')).drop(columns="assembly_accession_df2").drop_duplicates(subset=["assembly_accession","seqid"],keep='last')
                    for index, row in dropped.iterrows():
                        #Check if prot is not IMGVR or Plasmid
                        future = executor.submit(extract_subsequence, self.assembly_folder + str(row['assembly']), row["seqid"], int(row['start']), int(row['stop']), row['strand'], row["protein_accver"], self.mod, self.wn,row["flipped"])
                        futures.append(future)
                    for i,future in enumerate(concurrent.futures.as_completed(futures)):
                        neigh = future.result()
                        if neighs is None:
                            neighs=neigh
                        else:
                            if not neigh is None:
                                neighs = pd.concat([neighs,neigh],ignore_index=True)

                neighs.to_fasta("contig_id","sequence",self.output+"/neighborhood/neighborhoods.fasta")

            

        if self.ncrna:
            with console.status("[bold green]Annotating ncRNAs...") as status:
                console.print("🧬\tRunning Infernal")
                if not os.path.exists(self.output+'/ncrna'):
                    os.makedirs(self.output+'/ncrna')
                command = ["cmsearch","--tblout", self.output+"/ncrna/results.txt","-A",self.output+"/ncrna/results.sto", "-E", "1", "--incE","10", "--cpu", 
                            str(self.num_threads), os.environ.get('CONDA_PREFIX')+"/db/cm_models/all.cm", self.output+"/neighborhood/neighborhoods.fasta"]
                subprocess.run(command, check=True,stdout=subprocess.DEVNULL)
                column_names = ['nucid', '-', 'nc_feature',"--","cm","mdlfrom","mdlto","seqfrom","seqto","strand_ncrna","trunc","pass","gc","bias","score","E-value","inc","desc"]
                if not os.path.getsize(self.output+"/ncrna/results.sto") == 0:
                    cmdf = pd.read_csv(self.output+'/ncrna/results.txt', sep=r'\s+', engine='python',comment="#",header=None,names=column_names)
                    for record in SeqIO.parse(self.output+"/ncrna/results.sto", "stockholm"):
                        seqfrom = record.id.split("/")[1].split("-")[0]
                        seqto = record.id.split("/")[1].split("-")[1]
                        seqid = record.id.split("/")[0]
                        sequence= str(record.seq).replace(".","").replace("-","")
                        criteria = {'nucid': str(seqid),
                                    'seqfrom': int(seqfrom),
                                    'seqto': int(seqto)}
                        cmdf.loc[(cmdf['nucid'] == criteria['nucid']) &
                                        (cmdf['seqfrom'] == criteria['seqfrom']) &
                                        (cmdf['seqto'] == criteria['seqto']),"sequence"]=sequence
                    results_ncRNA=pd.merge(neighs[["contig_id","start_win","end_win","strand","start_prot","end_prot","prot_id","flipped"]],cmdf,left_on="contig_id",right_on="nucid")

                    results_ncRNA["ncRNA_start"]=results_ncRNA["start_win"]+results_ncRNA["seqfrom"]
                    results_ncRNA["ncRNA_end"]=results_ncRNA["start_win"]+results_ncRNA["seqto"]
                    results_ncRNA["rel_start"]=results_ncRNA["ncRNA_start"]-results_ncRNA["start_prot"]
                    results_ncRNA["rel_end"]=results_ncRNA["ncRNA_end"]-results_ncRNA["start_prot"]
                    results_ncRNA["delta"]=results_ncRNA["end_prot"]-results_ncRNA["start_prot"]

                    upstream = results_ncRNA["rel_start"]>0
                    downstream = results_ncRNA["rel_start"]<0
                    flipped = results_ncRNA["flipped"]=="flipped"

                    new_start = np.where(upstream & flipped, -(results_ncRNA["rel_end"] - results_ncRNA["delta"]),
                                        np.where(downstream & flipped, -(results_ncRNA["rel_end"]) + results_ncRNA["delta"], results_ncRNA['rel_start']))

                    new_end = np.where(upstream & flipped, -(results_ncRNA["rel_start"] - results_ncRNA["delta"]),
                                    np.where(downstream & flipped, -(results_ncRNA["rel_start"]) + results_ncRNA["delta"], results_ncRNA['rel_end']))


                    results_ncRNA['rel_start'],results_ncRNA['rel_end']=new_start,new_end

                    for index, row in results_ncRNA.iterrows():
                        matching_rows = results[results['id'] == row["prot_id"]]
                        for i,j in matching_rows.iterrows():
                            if j["start"]-self.wn<row["ncRNA_end"]<j["end"]+self.wn:
                                results_ncRNA.at[index,"y"]=j["y"]
                                results_ncRNA.at[index,"product"]="ncRNA"
                                results_ncRNA.at[index,"target_prot"]=j["target_prot"]
                                results_ncRNA.at[index,"assembly_accession"]=j["assembly_accession"]
                                results_ncRNA.at[index,"species"]=j["species"]
                                results_ncRNA.at[index,"genus"]=j["genus"]
                                results_ncRNA.at[index,"family"]=j["family"]
                                results_ncRNA.at[index,"order"]=j["order"]
                                results_ncRNA.at[index,"class"]=j["class"]
                                results_ncRNA.at[index,"phylum"]=j["phylum"]
                                results_ncRNA.at[index,"kingdom"]=j["kingdom"]
                                results_ncRNA.at[index,"seqid"]=j["seqid"]

                    y_step = self.den_data['y'][1] - self.den_data['y'][0]
                    y_half_height = y_step/7
                    neg_strand = results_ncRNA['strand']=="-"
                    results_ncRNA["dx"] = results_ncRNA["rel_end"]-results_ncRNA["rel_start"]
                    gene_x_tail = results_ncRNA['rel_start']
                    gene_dx = results_ncRNA['dx']
                    gene_x_head = gene_x_tail + gene_dx
                    gene_x_head_start = gene_x_head
                    text_x = gene_x_tail + (gene_x_head_start-gene_x_tail)/2
                    results_ncRNA['text_x']=text_x
                    results_ncRNA['text_y']=results_ncRNA['y']+y_half_height
                    results_ncRNA['xs']=[list(n) for n in zip(list(gene_x_tail.tolist())[:],list(gene_x_tail.tolist())[:],gene_x_head_start.tolist()[:],gene_x_head.tolist()[:],gene_x_head_start.tolist()[:])]
                    def divide_by_10(lst):
                        return [x / 100 for x in lst]

                    # Apply the function to the column
                    results_ncRNA['xs'] = results_ncRNA['xs'].apply(divide_by_10)
                    results_ncRNA['ys']=[list(n) for n in zip(results_ncRNA['y']-y_half_height, results_ncRNA['y']+y_half_height, results_ncRNA['y']+y_half_height, results_ncRNA['y'], results_ncRNA['y']-y_half_height)]
                    results_ncRNA["coordinates"]=[[n] for n in [list(list(x) for x in zip(n[0],n[1])) for n in list(zip(results_ncRNA["xs"],results_ncRNA["ys"]))]]

                    families = results_ncRNA["nc_feature"].dropna().unique().tolist()
                    colors_rgb = [plt.cm.rainbow(random.random()) for _ in families]
                    colors_hex = [matplotlib.colors.to_hex(color) for color in colors_rgb]
                    colors_dic = {num:desaturate([n for n in color],0.6,1) for num,color in zip(families,colors_rgb)}
                    results_ncRNA["fillcolor"]=results_ncRNA["nc_feature"].map(colors_dic).apply(lambda d: d if isinstance(d, list) else [230, 230, 230,255])
                    #results_ncRNA = pd.concat([results,results_ncRNA],join="outer")   

        if self.cctyper:
            cctyper = self.cctyper
            with console.status("[bold green]Annotating CRISPRs...") as status:
                console.print("✂️\tRunning CRISPRCasTyper")
                cctyper_array=False
                if not os.path.exists(self.output+'/cctyper'):
                    os.makedirs(self.output+'/cctyper')            
                command = ["cctyper","--gff", self.output+"/results.gff","--prot",self.output+"/results.fasta", "-t",str(self.num_threads),self.output+"/neighborhood/neighborhoods.fasta",self.output+"/cctyper/output"]
                subprocess.run(command, check=True,stdout=subprocess.DEVNULL)


                if os.path.exists(self.output+"/cctyper/output/cas_operons.tab"):
                    cctyper_df = pd.read_csv(self.output+"/cctyper/output/cas_operons.tab",sep="\t")
                    cctyper_df['Genes']=cctyper_df['Genes'].apply(literal_eval)
                    cctyper_df['Prot_IDs']=cctyper_df['Prot_IDs'].apply(literal_eval)
                    split_data = {'Genes': [], 'Prot_IDs': [],'Best_type': []}

                    # Loop through the original DataFrame and split the lists
                    for _, row in cctyper_df.iterrows():
                        genes_list = row['Genes']
                        prot_ids_list = row['Prot_IDs']
                        best_type = row['Best_type']
                        for gene, prot_id in zip(genes_list, prot_ids_list):
                            split_data['Genes'].append(gene)
                            split_data['Prot_IDs'].append(prot_id)
                            split_data['Best_type'].append(best_type)
                    cctyper_df = pd.DataFrame(split_data)
                    cctyper_df= cctyper_df.rename(columns={'Best_type': 'cctyper_system',"Genes":"cctyper_gene"})
                    results=pd.merge(results, cctyper_df[["cctyper_system","cctyper_gene","Prot_IDs"]], left_on='id', right_on='Prot_IDs',how="left").drop(columns='Prot_IDs')
                    
                    if os.path.exists(self.output+"/cctyper/output/crisprs_all.tab"):
                        cctyper_array=True
                        crispr_df = pd.read_csv(self.output+'/cctyper/output/crisprs_all.tab', sep='\t', engine='python',header=0)
                        results_crispr=pd.merge(neighs[["contig_id","start_win","end_win","strand","start_prot","end_prot","prot_id","flipped"]],crispr_df,left_on="contig_id",right_on="Contig")


                        results_crispr["crispr_start"]=results_crispr["start_win"]+results_crispr["Start"]
                        results_crispr["crispr_end"]=results_crispr["start_win"]+results_crispr["End"]
                        results_crispr["rel_start"]=results_crispr["crispr_start"]-results_crispr["start_prot"]
                        results_crispr["rel_end"]=results_crispr["crispr_end"]-results_crispr["start_prot"]
                        results_crispr["delta"]=results_crispr["end_prot"]-results_crispr["start_prot"]

                        upstream = results_crispr["rel_start"]>0
                        downstream = results_crispr["rel_start"]<0
                        flipped = results_crispr["flipped"]=="flipped"

                        new_start = np.where(upstream & flipped, -(results_crispr["rel_end"] - results_crispr["delta"]),
                                            np.where(downstream & flipped, -(results_crispr["rel_end"]) + results_crispr["delta"], results_crispr['rel_start']))

                        new_end = np.where(upstream & flipped, -(results_crispr["rel_start"] - results_crispr["delta"]),
                                        np.where(downstream & flipped, -(results_crispr["rel_start"]) + results_crispr["delta"], results_crispr['rel_end']))


                        results_crispr['rel_start'],results_crispr['rel_end']=new_start,new_end



                        for index, row in results_crispr.iterrows():
                            matching_rows = results[results['id'] == row["prot_id"]]
                            for i,j in matching_rows.iterrows():
                                if j["start"]-self.wn<row["crispr_end"]<j["end"]+self.wn:
                                    results_crispr.at[index,"y"]=j["y"]
                                    results_crispr.at[index,"product"]="CRISPR array"
                                    results_crispr.at[index,"target_prot"]=j["target_prot"]
                                    results_crispr.at[index,"assembly_accession"]=j["assembly_accession"]
                                    results_crispr.at[index,"species"]=j["species"]
                                    results_crispr.at[index,"genus"]=j["genus"]
                                    results_crispr.at[index,"family"]=j["family"]
                                    results_crispr.at[index,"order"]=j["order"]
                                    results_crispr.at[index,"class"]=j["class"]
                                    results_crispr.at[index,"phylum"]=j["phylum"]
                                    results_crispr.at[index,"kingdom"]=j["kingdom"]
                                    results_crispr.at[index,"seqid"]=j["seqid"]

                        y_step = self.den_data['y'][1] - self.den_data['y'][0]
                        y_half_height = y_step/7
                        neg_strand = results_crispr['strand']=="-"
                        results_crispr["dx"] = results_crispr["rel_end"]-results_crispr["rel_start"]
                        gene_x_tail = results_crispr['rel_start']
                        gene_dx = results_crispr['dx']
                        gene_x_head = gene_x_tail + gene_dx
                        gene_x_head_start = gene_x_head
                        text_x = gene_x_tail + (gene_x_head_start-gene_x_tail)/2
                        results_crispr['text_x']=text_x
                        results_crispr['text_y']=results_crispr['y']+y_half_height
                        results_crispr['xs']=[list(n) for n in zip(list(gene_x_tail.tolist())[:],list(gene_x_tail.tolist())[:],gene_x_head_start.tolist()[:],gene_x_head.tolist()[:],gene_x_head_start.tolist()[:])]
                        def divide_by_10(lst):
                            return [x / 100 for x in lst]

                        # Apply the function to the column
                        results_crispr['xs'] = results_crispr['xs'].apply(divide_by_10)
                        results_crispr['ys']=[list(n) for n in zip(results_crispr['y']-y_half_height, results_crispr['y']+y_half_height, results_crispr['y']+y_half_height, results_crispr['y'], results_crispr['y']-y_half_height)]
                        results_crispr["coordinates"]=[[n] for n in [list(list(x) for x in zip(n[0],n[1])) for n in list(zip(results_crispr["xs"],results_crispr["ys"]))]]
                        results_crispr.rename(columns={'Subtype': 'nc_feature'},inplace=True)


        for n in ["padloc_gene", "padloc_system", "deffinder_gene", "deffinder_system", "cctyper_gene", "cctyper_system"]:
            if n not in results.columns:
                results[n] = np.nan

        results['defense'] = ""
        if self.padloc:
            results['defense'] = results['defense'] + results["padloc_gene"].fillna("")
        if self.deffinder:
            results['defense'] = results['defense'] + results["deffinder_gene"].fillna("") 
        if self.cctyper:
            results['defense'] = results['defense'] + results["cctyper_gene"].fillna("")

        results['defense'].replace('', np.nan, inplace=True)
        indexes_df = results['defense'].dropna().index
        indexes_nodf = results['defense'].isna().index
        results.loc[indexes_nodf, 'linecolor'] = results.loc[indexes_nodf, 'fillcolor'].apply(lambda x: darken(x, 0.5, 1))
        results.loc[indexes_df, 'linecolor'] = results.loc[indexes_df, 'fillcolor'].apply(lambda x: [0, 0, 0,255])
        #results.loc[indexes, 'fillcolor']=results.loc[indexes, 'fillcolor'].apply(lambda d: [0, 0, 0,255])
        results=results.drop_duplicates(subset=["id","target_prot"])
        results['linetype'] = np.where(pd.notna(results['defense']), 0.3, 0.2)

        results=results.drop_duplicates(subset=['start','end',"seqid","target_prot"])
        results[["id","target_prot","seqid","fam_cluster","padloc_system","padloc_gene","deffinder_system","deffinder_gene","start","end","rel_start","rel_end","strand","flip_strand","product","locus_tag","assembly_accession","species_taxid","kingdom","phylum","class","order","family","genus","species","linecolor","fillcolor"]].to_csv(self.output+"/results_defense.txt",index=False,header=False)

        results_ncfeatures =  pd.DataFrame()
        if self.ncrna or self.cctyper:
            if self.ncrna:
                results_ncRNA["nc_feature"] = "ncRNA " + results_ncRNA["nc_feature"]
                results_ncfeatures = pd.concat([results_ncfeatures, results_ncRNA], ignore_index=True)
            if self.cctyper:
                if cctyper_array:
                    results_crispr["nc_feature"] = "CRISPR array " + results_crispr["nc_feature"]
                    results_ncfeatures = pd.concat([results_ncfeatures, results_crispr], ignore_index=True)
                else:
                    cctyper=None
            
            if self.cctyper or self.ncrna:
                families = results_ncfeatures["nc_feature"].dropna().unique().tolist()
                colors_rgb = [plt.cm.rainbow(random.random()) for _ in families]
                colors_hex = [matplotlib.colors.to_hex(color) for color in colors_rgb]
                colors_dic = {num:desaturate([n for n in color],0.6,1) for num,color in zip(families,colors_rgb)}
                results_ncfeatures["fillcolor"]=results_ncfeatures["nc_feature"].map(colors_dic).apply(lambda d: d if isinstance(d, list) else [230, 230, 230,255])
        self.ncfeatures=results_ncfeatures
