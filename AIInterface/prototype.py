from azure.search.documents import SearchClient
from azure.identity import DefaultAzureCredential

import openai
from openai import AzureOpenAI

import os

credential = DefaultAzureCredential()


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
"10297118",
"2886234",
"3159683",
"3534140",
"4526073",
"4787201",
"4794711",
"5633377",
"5736644",
"6194747",
"7116826",
"7427288",
"7691105",
"7909819"
]
search_client = SearchClient(
    endpoint=search_endpoint,
    index_name=index_name,
    credential=credential
)

for pmcid in pmcids:
    results = search_client.search(
        search_text="what variants are mentioned?",
        filter=f"pmcid eq '{pmcid}'",
        semantic_configuration_name="ai-search-1749232431805-semantic-configuration",
        top=5
    )

    combined_chunks = "\n".join([doc["chunk"] for doc in results])

    # Construct the prompt
    system_msg = {
        "role": "system",
        "content": "Extract all genetic variants mentioned and return a table with columns: gene, variant, phenotype, technology used, functional study, citation."
    }
    user_msg = {
        "role": "user",
        "content": combined_chunks
    }


    client = AzureOpenAI(
        azure_endpoint="https://oric-rag-demo-resource.openai.azure.com",
        api_key=os.getenv("AZURE_OPENAI_API_KEY"),  # or use DefaultAzureCredential if set up
        api_version="2024-02-15-preview"  # or whatever version you're using
    )

    response = client.chat.completions.create(
        model="gpt-4.1",  # deployment ID
        messages=[
            {"role": "system", "content": "You are an assistant that summarizes genetic variant information."},
            {"role": "user", "content": combined_chunks}
        ]
    )
    print(f"\nResults for PMCID {pmcid}:\n")
    print(response.choices[0].message.content)