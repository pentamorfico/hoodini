import sys, os
from time import sleep

import rich
import logging
from rich.logging import RichHandler
import rich_click as click
from rich.progress import Progress
from rich.console import Console
from rich.prompt import Prompt

from flagcsnap.stream import Controller
from flagcsnap.parse_ipg import IPGParser
from flagcsnap.parse_assemblies import AssemblyParser
from flagcsnap.cluster_proteins import ProteinClusterer
from flagcsnap.parse_taxonomy import TaxonomyParser
from flagcsnap.arrange_data import Arranger
from flagcsnap.basic_plot import Plotter
from flagcsnap.extra_tools import ExtraAnnotation
from flagcsnap.extra_plot import ExtraPlotter


from flagcsnap.utils import console,log,validate_input_file

# Reading parameters and inputs

click.rich_click.USE_RICH_MARKUP = True
click.rich_click.STYLE_ERRORS_SUGGESTIONS = "yellow"
click.rich_click.STYLE_HELP_OPTIONS = "bold cyan"
click.rich_click.STYLE_HELP_OPTIONS_DEFAULTS = "dim"

@click.command()
@click.option('--max-concurrent-downloads', default=10, type=int, help='Maximum concurrent downloads.')
@click.option('--api-key', "apikey", default='0dd59f676839567a5bb2eb826cff33c8fc09', help='NCBI API Key.')
@click.option('--input', 'input_path', type=str, help='Path to input file.',required=True, callback=validate_input_file)
@click.option('--num-threads',"num_threads", default=10, type=int, help='Number of threads.')
@click.option('--assembly-folder',"assembly_folder", default=None, help='Path to assembly folder.')
@click.option('--cand-mode',"cand_mode", default="best_id", help='Mode of chossing best Identical Protein Groups representative')
@click.option('--win-mode',"mod", default='win_nts', help='Genomic window in numer of genes or nucleotides')
@click.option('--win',"wn", default=20000, type=int, help='Window size')
@click.option('--height-factor', default=20, type=int, help='Height factor.')
@click.option('--padloc', is_flag=True, help='Run PADLOC for antiphage defense systems identification')
@click.option('--deffinder', is_flag=True, help='Run DefenseFinder for antiphage defense systems identification')
@click.option('--tree-mode', default='fast_phylo', help='Tree building method')
@click.option('--tree-file', default='target_prots.nwk', help='Path to tree file.')
@click.option('--output', default=None, help='Output folder name.')
@click.option('--ncRNA',"ncrna", is_flag=True, help='ncRNA prediction using Infernal')
@click.option('--cctyper', is_flag=True, help='CRISPR-Cas systems prediction using CCtyper.')
@click.option('--ngenes', default=10, type=int, help='Number of genes.')
@click.option('--clust-method', default='diamond_deepclust', help='Clustering method.')
@click.option('--assembly-db', default=os.environ.get('CONDA_PREFIX')+"/db/flagcsnap", help='Assembly DB name.')
@click.option('--img-db', default=None, help='Additional files from IMG databases search')
@click.option('--img-nuc', default=None, help='Additional files from IMG databases search')


# Main function


def main(**kwargs):
    """
    flagcsnap: A tool for gene-centric comparative genomic analysis using publicly available data

    Read the documentation at: 
    """    

    # Setting up logging
    FORMAT = "%(message)s"
    logging.basicConfig(
        level="NOTSET", format=FORMAT, datefmt="[%X]", handlers=[RichHandler()]
    )
    log = logging.getLogger("rich")
    console = Console()
    console.print("""_____________________________________________________________________________________________""")
    console.print("""\n[bold]flagcsnap[/bold]: A tool for gene-centric comparative genomic analysis using publicly available data\n""")
    console.print("""¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯""")

    ## Initialize
    

    console.print("🚀\tInitializing flagcsnap")

    controller = Controller(**kwargs)
    controller.run()


    console.print("📥\tFetching data from NCBI Identical Protein Groups")

    ipgparser = IPGParser(controller)
    ipgparser.run()

    console.print("📥\tDownloading assemblies from NCBI")

    assemblyparser = AssemblyParser(ipgparser)
    assemblyparser.run()

    console.print("✨\tClustering neighbor proteins")

    assemblyparser = ProteinClusterer(assemblyparser)
    assemblyparser.run()

    console.print("🦠\tExtracting taxonomic information")

    taxparser = TaxonomyParser(assemblyparser)
    taxparser.run()

    console.print("🧹\tArranging data")
    arranger = Arranger(taxparser)
    arranger.run()

    console.print("🖼️\tPlotting data")
    plotter = Plotter(arranger)
    plotter.run()

    console.print("🧬\tRunning extra annotation tools")
    plotter = ExtraAnnotation(arranger)
    plotter.run()

    console.print("🖼️\tPlotting extra data")
    plotter = ExtraPlotter(plotter)
    plotter.run()

if __name__ == '__main__':
    main()
