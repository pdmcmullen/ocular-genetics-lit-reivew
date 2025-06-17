from azure.search.documents import SearchClient
from azure.identity import DefaultAzureCredential
import unicodedata


import openai
from openai import AzureOpenAI

import os

credential = DefaultAzureCredential()



gpt_prompt = """You are an assistant that summarizes evidence for genetic
variants in the scientific literature. This context comes from a single
publication describing a scientific study. Your goal is to extract and tabulate information about the participants in the present research study, 
ignoring previously existing evidence for variants that this publication cites or references. Some participants will be have a relevant phenotype and therefore should be denoted "affected", while otherse will not.
Please prepare a table with a row for each individual described in the
manuscript. Include the following
columns: gene, whether the individual is affected, variant (if any, in HGVS c-dot nomenclature), variant (if any, in HGVS p-dot
nomenclature), the rs-identifier of the variant, patient identifier used in the
paper, the relationship of the individual to the study proband, phenotype(s), molecular technology used to determine the variant, and
novel functional characterization of the variant, and the specific context you used to determine that this individual should be included in the table. Report only variants in the
gene CRYAA. Report the table as both human-readable and json formats."""



# Setup
search_endpoint = "https://rossenmanuscriptsearchcatalog.search.windows.net"
index_name = "manuscript-index"
vector_config = "default-vector-profile"

openai.api_type = "azure"
openai.api_base = "https://oric-rag-demo-resource.openai.azure.com/"
openai.api_version = "2024-02-15-preview"
openai.api_key = os.getenv("AZURE_OPENAI_API_KEY")  # Set this in your environment


# List of PMCIDs to loop over
pmcids = [
#"3159683", # focus on one paper for now Sun et al 2011
"3534140", # su et al 2012
#"10297118",
#"2886234",
#"4526073",
#"4787201",
#"4794711",
#"5633377",
#"5736644",
#"6194747",
#"7116826",
#"7427288",
#"7691105",
#"7909819"
]
search_client = SearchClient(
    endpoint=search_endpoint,
    index_name=index_name,
    credential=credential
)

for pmcid in pmcids:
    results = search_client.search(
        search_text="*",
        filter=f"pmcid eq '{pmcid}'",
        semantic_configuration_name="ai-search-1749232431805-semantic-configuration"
    )

    combined_chunks = "\n".join([doc["chunk"] for doc in results])
    safe_chunks = unicodedata.normalize("NFKC", combined_chunks)






    client = AzureOpenAI(
        azure_endpoint="https://oric-rag-demo-resource.openai.azure.com",
        api_key=os.getenv("AZURE_OPENAI_API_KEY"), 
        api_version="2024-02-15-preview"  
    )

    response = client.chat.completions.create(
        model="gpt-4.1",  # deployment ID
        messages=[
            {"role": "system", "content": gpt_prompt},
            {"role": "user", "content": safe_chunks}
        ]
    )
    print(f"\nResults for PMCID {pmcid}:\n")
    print(response.choices[0].message.content)