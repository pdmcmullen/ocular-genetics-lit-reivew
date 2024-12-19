"""Given a list of pubmed identifiers, determine whether the cooresponding paper
is available in PMC, and download if possible"""


import os
import requests
from Bio import Entrez

Entrez.email = 'pmcmullen@luriechildrens.org'

def parse_pmids_from_file(data_file_path: str) -> list[str]:
    """Gather a list of pmids from a file, one per line.

    Args:
        data_file_path (str): Path to a file containing a list of pubmed IDs.

    Returns:
        list[str]: List of pubmed IDs.
    """

    # TODO: open file, read data, check that pmids only contain digits, return as a list of strings

    return pmids


def query_pubmed_for_pmcid(pmid):
    pass

    # do something to throttle https requests so we don't overwhelm the system by building in a short pause (.2 seconds or something)

    # call Entrez.efetch to get the data corresponding to the pmid
    # https://biopython.org/docs/1.75/api/Bio.Entrez.html#Bio.Entrez.efetch

    # parse the record returned by efetch to identify the pmcid (if it exists)

    # verify that this is a string that starts with "PMC" and is followed by digits

    # return the pmcid if it exists; None if not


def download_paper(pmid, pmc_id, output_dir):
    pass

    # construct the url corresponding to the paper
    pdf_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmc_id}/pdf/"
    

    response = requests.get(pdf_url)
    if response.status_code == 200:
        pdf_path = os.path.join(output_dir, f"{pmid}-{pmc_id}.pdf") # we probably want to make these a little more human-readable. e.g., include the gene name, first author in the file name.
        with open(pdf_path, "wb") as pdf_file:
            pdf_file.write(response.content)




def main():
    pmids = parse_pmids_from_file('Data/gcep_publications_hgmd2024-1_unique_pmids.txt')
    for pmid in pmids:
        # call query_pubmed_for_pmcid, call download_paper if applicable
        pass


if __name__ == '__main__':
    main()