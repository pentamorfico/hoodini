import itertools
import sys

import ete3

from hoodini.utils.logging_utils import error, info, warn


def calculate_taxid_distances(taxids, update_db=False):
    """
    Calculate pairwise distances between given taxonomic IDs using ete3 by passing node objects.

    Parameters:
    - taxids (list): A list of taxonomic IDs (integers or strings).
    - update_db (bool): Whether to update the local NCBI taxonomy database. Default is False.

    Returns:
    - dict: A dictionary with keys as tuples of taxid pairs and values as their distances.

    Raises:
    - ValueError: If any of the taxids are missing from the taxonomy database.
    """

    taxids_str = [str(int(taxid)) for taxid in taxids]

    ncbi = ete3.NCBITaxa()

    if update_db:
        try:
            info("Updating NCBI taxonomy database. This may take a while...")
            ncbi.update_taxonomy_database()
            info("Taxonomy database updated successfully.")
        except Exception as e:
            error(f"Error updating taxonomy database: {e}")
            sys.exit(1)

    try:
        tree = ncbi.get_topology(taxids_str, intermediate_nodes=True)
    except Exception as e:
        error(f"Error retrieving topology: {e}")
        sys.exit(1)

    tree_taxids = set()
    taxid_to_node = {}
    for node in tree.traverse():
        try:
            taxid = int(node.name)
            tree_taxids.add(taxid)
            taxid_to_node[taxid] = node
        except ValueError:
            continue

    try:
        taxids_int = [int(taxid) for taxid in taxids_str]
    except ValueError as ve:
        error(f"Error converting taxids to integers: {ve}")
        sys.exit(1)

    missing_taxids = [taxid for taxid in taxids_int if taxid not in tree_taxids]
    if missing_taxids:
        warn("The following taxids are missing from the taxonomy tree:")
        for mtaxid in missing_taxids:
            try:
                name = ncbi.get_taxid_translator([mtaxid]).get(mtaxid, "Unknown")
            except Exception:
                name = "Unknown"
            warn(f" - {mtaxid} ({name})")
        raise ValueError(
            "Some taxids are missing from the taxonomy tree. Please verify their validity."
        )
    else:
        info("All taxids are present in the taxonomy tree.")

    distances = {}

    info("Calculating pairwise distances using node objects...")
    for taxid1, taxid2 in itertools.combinations(taxids_int, 2):
        try:
            node1 = taxid_to_node[taxid1]
            node2 = taxid_to_node[taxid2]
            distance = tree.get_distance(node1, node2)
            distances[(taxid1, taxid2)] = distance
        except Exception as e:
            warn(f"Error calculating distance between {taxid1} and {taxid2}: {e}")

    if not distances:
        raise ValueError("No distances were calculated. Please check the taxids and try again.")

    info("Pairwise distances calculated successfully.")

    return distances
