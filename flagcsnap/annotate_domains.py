import subprocess
from io import StringIO
import pandas as pd
from flagcsnap.utils import console,merge_cluster_result
from rich.progress import Progress


class AnnotDomains:
    def __init__(self, obj):
        self.master = obj
        for key, val in vars(obj).items():
            setattr(self, key, val)

    def run(self):

        with console.status("[bold green]Annotating protein domains with MMseqs2...") as status:

            command = ["mmseqs", "easy-search", self.output+"/results.fasta", self.domains, self.output+"/domains.m8", self.output+"/tmp"]
            result = subprocess.run(command, check=True,capture_output=True)


            column_names = ['protein_id', 'domain_id', 'identity', 'alignment_length', 'mismatches', 
            'gap_openings', 'q_start', 'q_end', 's_start', 's_end', 'e_value', 'bit_score']

            domains_df = pd.read_csv(self.output+"/domains.m8",names=column_names,sep="\t")

            domains_df['bit_score'] = domains_df['bit_score'].astype(float)
            
            domains_df['alignment_length'] = domains_df['alignment_length'].astype(float)
            domains_df['bit_score*alignment_length'] = domains_df['bit_score']*domains_df['alignment_length']
            domains_df = domains_df.sort_values(by=['bit_score*alignment_length'],ascending=False)

            grouped = domains_df.groupby('protein_id')

            all_domains_df = pd.DataFrame(columns=['protein_id',"domain_id","e_value", 'start', 'end', 'y_pos', 'bit_score'])

            def find_next_y_pos(domain, df):
                y_pos = -1
                while True:
                    overlapping_domains = df[(df['y_pos'] == y_pos) &
                                            ~((df['q_end'] <= domain['q_start']) | (df['q_start'] >= domain['q_end']))]
                    if overlapping_domains.empty:
                        return y_pos
                    y_pos += 1

            for protein_id, group in grouped:
                group['y_pos'] = -1
                for index, domain in group.iterrows():
                    group.at[index, 'y_pos'] = find_next_y_pos(domain, group)
                new_df = group[['q_start', 'q_end', 'y_pos','bit_score',"domain_id","e_value"]].copy()
                new_df.rename(columns={'q_start': 'start', 'q_end': 'end'}, inplace=True)
                new_df['protein_id'] = protein_id
                all_domains_df = pd.concat([all_domains_df, new_df], ignore_index=True)

            all_domains_df = all_domains_df.merge(self.results[["id","rel_start","rel_end","flip_strand","y"]], right_on="id", left_on="protein_id", how="left")
            domain_metadata = pd.read_csv(self.domains_metadata,sep="\t",names=["domain_id","clan","family","short_name","description"])
            all_domains_df["domain_id"] = all_domains_df["domain_id"].str.split(".").str[0]
            all_domains_df = all_domains_df.merge(domain_metadata, on="domain_id", how="left") 

            self.domains_data = all_domains_df


            console.print(f"✔️\tDomain annotation done\n")

        


