import csv
from Bio import Entrez

def setup_entrez(email, api_key):
    """Configure Entrez with user email and API key."""
    Entrez.email = email
    Entrez.api_key = api_key

def fetch_pmid_from_doi(doi, email, api_key):
    """Fetch PMID from a given DOI using the Entrez system."""
    setup_entrez(email, api_key)
    try:
        handle = Entrez.esearch(db="pubmed", term=doi, rettype="uilist")
        record = Entrez.read(handle)
        handle.close()
        if record['IdList']:
            return record['IdList'][0], True
        else:
            return "No PMID found", False
    except Exception as e:
        return str(e), False

def process_dois(input_csv, output_csv, email, api_key):
    """Process a CSV file containing DOIs and output a CSV with DOIs and corresponding PMIDs."""
    with open(input_csv, mode='r', newline='') as infile, open(output_csv, mode='w', newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = ['DOI', 'PMID']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()

        for row in reader:
            doi = row['DOI']
            pmid, success = fetch_pmid_from_doi(doi, email, api_key)
            writer.writerow({'DOI': doi, 'PMID': pmid if success else "Not Found"})

# Usage
email = "YOU_EMAIL"
api_key = "API-KEY"
input_csv = "YOU.csv"
output_csv = "YOU.csv"

process_dois(input_csv, output_csv, email, api_key)
