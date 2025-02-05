"""Given a list of pubmed identifiers, determine whether the cooresponding paper
is available in PMC, and download if possible"""


import os
import requests
import xml.etree.ElementTree as ET
from Bio import Entrez
import gzip
import shutil
import tarfile
import csv

Entrez.email = 'pmcmullen@luriechildrens.org'

def parse_pmids_from_file(data_file_path: str) -> list[str]:
    """Gather a list of pmids from a file, one per line.

    Args:
        data_file_path (str): Path to a file containing a list of pubmed IDs.

    Returns:
        list[str]: List of pubmed IDs.
    """

    pmids = []

    with open(data_file_path, "r") as file:
        for pmid in file:
            pmid_str = str(pmid).strip()
            if pmid_str.isnumeric():
                pmids.append(pmid_str)

    return pmids

def query_pubmed_for_pmcid(pmid):

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

    #collecting metadata
    author_name = f'{tree.find(".//AuthorList/Author/LastName").text}' #, {tree.find(".//AuthorList/Author/ForeName").text}' (first name)
    title = tree.find(".//Article/ArticleTitle").text
    journal = tree.find(".//Journal/Title").text
    pub_year = tree.find(".//PubDate/Year").text

    return pmcid, author_name, title, journal, pub_year

def download_paper(pmid, pmc_id, temp_workspace, output_dir):
    # construct the url corresponding to the paper
    metadata_url = f"https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id={pmc_id}"
    folder_name = None
    
    response = requests.get(metadata_url)
    if response.status_code == 200:
        if "is not Open Access" in response.content.decode('UTF-8'):
            #add to a list of pmcids/pmids without open access
            return False
        else:
            tree = ET.fromstring(response.content.decode('UTF-8'))
            pdf_url = tree[2][0][0].attrib['href']
            gz_path = os.path.join(temp_workspace, f"{pmid}-{pmc_id}.tar.gz") 
            gz_url = pdf_url.replace('ftp:', 'https:')
            gz_response = requests.get(gz_url)
            with open(gz_path, "wb") as gz_file:
                gz_file.write(gz_response.content)

            with gzip.open(gz_path, 'rb') as gz_file:
                tar_path = os.path.join(temp_workspace, f"{pmid}-{pmc_id}.tar") 
                with open(tar_path, 'wb') as extracted_file:
                    shutil.copyfileobj(gz_file, extracted_file)

            with tarfile.open(tar_path, 'r') as tar:
                folder_name = os.path.join(output_dir, tar.getmembers()[0].path)
                tar.extractall(output_dir)
            return True

    return folder_name


def export_list(papers, csv_path):
    if len(papers) > 0:
        dict_keys = list(papers[0].keys())

        with open(csv_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=dict_keys)
            writer.writeheader()
            writer.writerows(papers)

def main():
    papers = []
    temp_workspace = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'temp')
    output_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'pdfs')

    #pmids = [25596306, 20884774] #hard coding to test with only two
    pmids = parse_pmids_from_file('Data/gcep_publications_hgmd2024-1_unique_pmids.txt')    
    # print(f'Found {len(pmids)} pmids')

    for pmid in pmids:
        try:
            paper = {
                'PMID': pmid,
                'First_Author': None,
                'Title': None,
                'Journal': None,
                'Publication_Year': None,
                'PMCID': None,
                'Open_Access': None,
                'File_Name': None
                }
            # call query_pubmed_for_pmcid, call download_paper if applicable
            pmcid, author_name, title, journal, pub_year = query_pubmed_for_pmcid(pmid)
            paper['First_Author'] = author_name
            paper['Title'] = title
            paper['Journal'] = journal
            paper['Publication_Year'] = pub_year

            if pmcid is not None:
                paper['PMCID'] = pmcid
                folder_name = download_paper(pmid, pmcid, temp_workspace, output_dir)
                #if folder_name is not None:
                    #{pmid}_{pmcid}_{first author surname}_{journal}_{year}.pdf
                    #folder_name = f'{pmid}_{pmcid}_{author_name}_{journal}_{pub_year}.pdf'
                    #pass
                paper['Open_Access'] = folder_name
                if folder_name == True:
                    paper["File_Name"] = f'{pmid}_{pmcid}_{author_name}_{journal}_{pub_year}.pdf'

        except:
            pass

        papers.append(paper)

    # Export the list of papers and metadata to csv
    export_list(papers, os.path.join(os.path.dirname(os.path.realpath(__file__)), 'papers.csv'))

    # Delete everything in temp workspace


if __name__ == '__main__':
    main()