"""Given a list of pubmed identifiers, determine whether the cooresponding paper
is available in PMC, and download if possible"""


import os
import requests
import xml.etree.ElementTree as ET
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

    pmids = []

    with open(data_file_path, "r") as file:
        for pmid in file:
            pmid_str = str(pmid).strip()
            if pmid_str.isnumeric():
                pmids.append(pmid_str)

    return pmids

def query_pubmed_for_pmcid(pmid):
    # do something to throttle https requests so we don't overwhelm the system by building in a short pause (.2 seconds or something)

    # call Entrez.efetch to get the data corresponding to the pmid
    # https://biopython.org/docs/1.75/api/Bio.Entrez.html#Bio.Entrez.efetch

    # parse the record returned by efetch to identify the pmcid (if it exists)

    # verify that this is a string that starts with "PMC" and is followed by digits

    # return the pmcid if it exists; None if not
    
    pmcid = None
    
    # Get the result summary
    handle = Entrez.efetch(db = "pubmed", id = pmid)
    e_summary_result = handle.read().strip().decode('UTF-8')
    handle.close()

    # Convert to XML and traverse to find the pmc id
    tree = ET.fromstring(e_summary_result)
    for root_child in tree[0]:
        if root_child.tag == 'PubmedData':
            for pubmed_data_child in root_child:
                if pubmed_data_child.tag == 'ArticleIdList':
                    for aid_child in pubmed_data_child:
                        if aid_child.attrib['IdType'] == 'pmc':
                            pmcid = aid_child.text
                            break

    return pmcid

def download_paper(pmid, pmc_id, output_dir):
    # construct the url corresponding to the paper
    pdf_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmc_id}/pdf/"
    
    response = requests.get(pdf_url)
    if response.status_code == 200:
        pdf_path = os.path.join(output_dir, f"{pmid}-{pmc_id}.pdf") # we probably want to make these a little more human-readable. e.g., include the gene name, first author in the file name.
        with open(pdf_path, "wb") as pdf_file:
            pdf_file.write(response.content)

def main():
    id_list = []

    # pmids = [25596306, 20884774] #hard coding to test with only two
    pmids = parse_pmids_from_file('Data/gcep_publications_hgmd2024-1_unique_pmids.txt')    
    # print(f'Found {len(pmids)} pmids')

    for pmid in pmids:
        # call query_pubmed_for_pmcid, call download_paper if applicable
        pmcid = query_pubmed_for_pmcid(pmid)
        if pmcid is not None:
            download_paper(pmid, pmcid, os.path.join(os.path.dirname(os.path.realpath(__file__)), 'pdfs'))


if __name__ == '__main__':
    main()