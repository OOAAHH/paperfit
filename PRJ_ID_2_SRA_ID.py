import csv
import argparse
from Bio import Entrez
from tqdm import tqdm

def fetch_sra_ids(project_id):
    """fetch SRA ID via project ID"""
    try:
        handle = Entrez.esearch(db="sra", term=project_id, retmax=10000)
        record = Entrez.read(handle)
        handle.close()
        return record["IdList"]
    except Exception as e:
        print(f"Error fetching SRA IDs for project {project_id}: {e}")
        return []

def main(email, api_key, input_csv_path, output_csv_path, view):
    Entrez.email = email
    Entrez.api_key = api_key

    # Read in IDs from csv files
    with open(input_csv_path, newline='') as csvfile:
        reader = csv.reader(csvfile)
        project_ids = [row[0] for row in reader]

    # Query SRA IDs and write down.
    with open(output_csv_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Project ID', 'SRA ID'])  # Built colunms
        for project_id in tqdm(project_ids, desc="Processing Project IDs"):
            sra_ids = fetch_sra_ids(project_id)
            for sra_id in sra_ids:
                full_sra_id = f"SRR{sra_id}"
                writer.writerow([project_id, full_sra_id])
                if view:  # Determines whether to print details according to the -view parameter
                    print(f"Written: {project_id}, {full_sra_id}")  # Print progress if need

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch SRA IDs for given project IDs.")
    parser.add_argument('-e', '--email', required=True, help="Your email address.")
    parser.add_argument('-k', '--api_key', default='', help="Your NCBI API key (optional).")
    parser.add_argument('-i', '--input_csv', required=True, help="Path to the input CSV file containing project IDs.")
    parser.add_argument('-o', '--output_csv', required=True, help="Path to the output CSV file.")
    parser.add_argument('-v', '--view', action='store_true', help="Display detailed progress information.")

    args = parser.parse_args()

    main(args.email, args.api_key, args.input_csv, args.output_csv, args.view)
