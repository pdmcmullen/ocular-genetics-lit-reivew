from azure.search.documents import SearchClient
from azure.identity import DefaultAzureCredential
import unicodedata


import openai
from openai import AzureOpenAI

import os

credential = DefaultAzureCredential()



gpt_prompt = """Do something."""




# Setup
search_endpoint = "https://rossenmanuscriptsearchcatalog.search.windows.net"
index_name = "manuscript-index"
vector_config = "default-vector-profile"

openai.api_type = "azure"
openai.api_base = "https://oric-rag-demo-resource.openai.azure.com/"
openai.api_version = "2024-02-15-preview"
openai.api_key = os.getenv("AZURE_OPENAI_API_KEY")  # Set this in your environment






client = AzureOpenAI(
    azure_endpoint="https://oric-rag-demo-resource.openai.azure.com",
    api_key=os.getenv("AZURE_OPENAI_API_KEY"), 
    api_version="2024-02-15-preview"  
)

response = client.chat.completions.create(
    model="gpt-5",  # deployment ID
    messages=[
        {"role": "system", "content": gpt_prompt},
        {"role": "user", "content": 'blah'}
    ]
    )
print(response)