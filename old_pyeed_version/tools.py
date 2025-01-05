#from msilib import sequence
from typing import List
from numpy import diff
from pathlib import Path
import os

from pyeed.core import ProteinInfo


def find_family(
        blast_results: List[ProteinInfo], 
        word: str
        ) -> Tuple[List, List]:
    """
    Separates proteins into two lists based on whether their name contains the given word or not.
    
    Args:
        blast_results (List[ProteinInfo]): A list of ProteinInfo objects obtained from a blast search.
        word (str): The word to search for in the protein names.
        
    Returns:
        Tuple[List, List]: A tuple containing two lists: found_family and different.
            found_family: A list of source IDs of proteins whose name contains the given word.
            different: A list of source IDs of proteins whose name does not contain the given word.
    """
    found_family = []
    different = []
    nameless = []

    for protein in blast_results:
        if protein.name:
            if word in protein.name:
                found_family.append(protein.source_id)
            else:
                different.append(protein.source_id)
        elif isinstance(protein, dict) and 'source_id' in protein:
            nameless.append(protein['source_id'])
        else:
            print("Protein is not a dictionary.")

    return found_family, different



def categorize_organism(
        blast_results: List[ProteinInfo], 
        category: str
        ) -> Dict[str, List[str]]:
    """
    Categorizes proteins based on their organism information.

    Args:
        blast_results (List[ProteinInfo]): A list of blast results, where each result is an instance of the ProteinInfo class.
        category (str): The category to use for categorizing the proteins. Can be "kingdom", "domain", or "phylum".

    Returns:
        Dict[str, List[str]]: A dictionary where the keys are the categories (kingdom, domain, or phylum) and the values are lists of protein source IDs belonging to that category.
    """
    proteins = {}
    length = len(blast_results)

    if category == "kingdom":
        category_key = "organism.kingdom"
        print("Kingdoms:")
    elif category == "domain":
        category_key = "organism.domain"
        print("Domains:")
    elif category == "phylum":
        category_key = "organism.phylum"
        print("Phylum:")
    else:
        raise NameError("The given category is not valid. Try kingdom, domain, or phylum")

    for protein in blast_results:
        category_value = getattr(protein, category_key)
        if category_value not in proteins:
            proteins[category_value] = [protein.source_id]
        else:
            proteins[category_value].append(protein.source_id)

    for key in proteins:
        percentage = round(len(proteins[key]) / length * 100, 4)
        print(f"Percentage of {key}: {percentage}%")

    return proteins


def categorize_len(blast_results: List[ProteinInfo]) -> Tuple[List, List, List]:
    """
    Categorizes a list of ProteinInfo objects based on the length of their sequences into three categories: small, middle, and long.

    Args:
    - blast_results (List[ProteinInfo]): A list of ProteinInfo objects representing the results of a protein blast search.

    Returns:
    - Tuple[List, List, List]: A tuple containing three lists: small, middle, and long.
      - small (List): A list of source IDs of proteins with sequence length less than 200.
      - middle (List): A list of source IDs of proteins with sequence length between 200 and 299.
      - long (List): A list of source IDs of proteins with sequence length 300 or more.
    """

    small = []
    middle = []
    long = []

    for protein in blast_results:
        sequence_length = len(protein.sequence)
        if sequence_length < 200:
            small.append(protein.source_id)
        elif sequence_length < 300:
            middle.append(protein.source_id)
        else:
            long.append(protein.source_id)

    return small, middle, long



def remove_duplicates(): 
    """
    The `remove_duplicates` function removes duplicate files from a directory based on the `source_id` attribute of `ProteinInfo` objects.
    """

    proteins = {}

    for i in range(7):
        for path in Path("./blast_results/blast_results_" + str(i+1)).rglob("*.json"):
            with open(str(path.absolute())) as f:
                protein = ProteinInfo.from_json(f)
                if not protein.source_id in proteins: 
                    proteins[protein.source_id] = [path]
                else:
                    proteins[protein.source_id].extend([path])

    for file_path in proteins.values():
        if len(file_path) > 1: 
            print(f"Duplicate files found:\n{path}\n")
            for path in file_path[1:]:
                os.remove(path)
                print(f"{path} has been deleted.\n")