import sys
import os
from time import sleep

import rich
import logging
from rich.logging import RichHandler
import rich_click as click
from rich.progress import Progress
from rich.console import Console
from rich.prompt import Prompt

from flagcsnap.utils import console, log, download_assembly_summary
import shutil
import pandas as pd


class Controller:
    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)

    def run(self):
        """
        Checks if the assembly database exists, if not, asks whether to download it or not
        """

        self.path_assembly_db = os.path.join(
            self.assembly_db, "assembly_dataframe.tsv")

        if not all(os.path.exists(path) for path in [self.assembly_db, self.path_assembly_db]):

            console.print(
                f"⚠️\t[bold][yellow]Assembly database not found at {self.assembly_db}[/bold]. Do you want to download it?")

            while True:
                user_input = Prompt.ask(
                    "⌨️\tPlease enter yes/no (y/N)").lower()

                if user_input in ["y", "yes"]:
                    console.print(
                        f"📁\tDownloading assembly database at {self.assembly_db}")
                    with console.status("[bold green]Downloading assembly summary data ...") as status:
                        download_assembly_summary(self.assembly_db)
                    break

                elif user_input in ["n", "no"]:
                    console.print(
                        f"🚫\t[bold][red]Aborting flagcsnap! [bold]{self.assembly_db}[/bold] folder won't be removed")
                    sys.exit()

                else:
                    # Prompt for a valid answer
                    console.print("❗\tPlease select a valid answer (y/N)")

            console.print(f"✔️\tAssembly database at {self.assembly_db}\n")

        else:

            console.print(f"✔️\tAssembly database at {self.assembly_db}\n")

        with console.status("[bold green]Reading assembly database ...") as status:

            self.assembly = pd.read_csv(
                self.path_assembly_db, sep="\t", engine="pyarrow", header=0)
            type_strains = pd.read_csv(
                self.assembly_db+"/type_strains_curated.txt", sep="\t", engine="pyarrow", header=0)
            self.type_strains = type_strains.rename(
                columns={'Genome seq. accession number': 'type_strain_assembly'})

        # Chech if outouf foder exists, if not, ask whether to remove it or not
        if self.output is None:
            self.output = self.input_path.rsplit("/", 1)[-1].rsplit(".", 1)[0]
            console.print(
                f"📥\tUsing [bold]{self.output}[/bold] as output folder {self.output}")

        if self.output is not None:
            if os.path.exists(self.output):
                console.print(
                    f"⚠️\t[bold][yellow]Output folder [bold]{self.output}[/bold] already exists. Do you want to remove it?")
                while True:
                    # Convert input to lowercase for consistency
                    user_input = Prompt.ask(
                        "⌨️\tPlease enter yes/no (y/N)").lower()
                    if user_input in ["y", "yes"]:
                        console.print(
                            f"🗑️\tRemoving content of output folder [bold]{self.output}[/bold]\n")
                        shutil.rmtree(self.output)
                        break  # Exit the loop after handling the input
                    elif user_input in ["n", "no"]:
                        console.print(
                            f"🚫\t[bold][red]Aborting flagcsnap! [bold]{self.output}[/bold] folder won't be removed")
                        sys.exit()
                    else:
                        # Prompt for a valid answer
                        console.print("❗\tPlease select a valid answer (y/N)")
            else:
                console.print(
                    f"✔️\tCreating output folder [bold]{self.output}[/bold]\n")
            os.makedirs(self.output)
        else:
            # If no output is declared, use the name of the input file

            self.output = os.getcwd()
            console.print(
                f"✔️\tUsing [bold]{self.output}[/bold] as output folder {self.output}\n")

        # Read input

        input_list = open(self.input_path, "r").read().splitlines()
        input_img = [n for n in input_list if n.startswith("IMG")]
        input_ncbi = [n for n in input_list if not n.startswith("IMG")]
        self.input_list = input_list
        self.input_img = input_img
        self.input_ncbi = input_ncbi

        console.print("✔️\tSuccessfuly read input file", style="bold green")
        console.print(
            f"➡️\t{len(input_list)} IDs have been read from the input file\n")
