import pandas as pd
import requests
from bs4 import BeautifulSoup
import time

# function, fetch PRJNA IDs via BeautifulSoup
def fetch_sra_from_gse(gse):
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}"
    response = requests.get(url)
    soup = BeautifulSoup(response.content, 'html.parser')
    link_elements = soup.select('a[href^="https://www.ncbi.nlm.nih.gov/bioproject/PRJ"]')
    sra_ids = []
    for link in link_elements:
        sra_id = link.get_text()
        if sra_id.startswith("PRJ"):
            sra_ids.append(sra_id)
    return sra_ids

# built you dataFrame here 
results_df = pd.DataFrame(columns=["GSE_ID", "PRJNA_IDs"])

# input your GEO IDs here
gse_ids = pd.read_csv('YOU_gseid.csv')


for gse in gse_ids:
    # wait, to avoid exceeding limits
    time.sleep(0.2)
    prjn_ids = fetch_sra_from_gse(gse)
    results_df = results_df.append({"GSE_ID": gse, "PRJNA_IDs": ', '.join(prjn_ids)}, ignore_index=True)

# output you results
results_df.to_csv('YOU_PRJNA_IDs_output.csv', index=False)
