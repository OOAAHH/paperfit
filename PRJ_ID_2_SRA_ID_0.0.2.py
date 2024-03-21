import csv
import os
import concurrent.futures
import subprocess
import re
from Bio import Entrez
from tqdm import tqdm
import argparse
import xml.etree.ElementTree as ET
import time
from functools import wraps




#====================================================
def rate_limited(max_per_second):
    min_interval = 1.0 / float(max_per_second)
    def decorate(func):
        last_called = [0.0]
        @wraps(func)
        def rate_limited_function(*args, **kwargs):
            elapsed = time.clock() - last_called[0]
            wait_required = min_interval - elapsed
            if wait_required > 0:
                time.sleep(wait_required)
            ret = func(*args, **kwargs)
            last_called[0] = time.clock()
            return ret
        return rate_limited_function
    return decorate

def retry(attempts=3, delay=2):
    "
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            for attempt in range(attempts):
                try:
                    return func(*args, **kwargs)
                except Exception as e:
                    print(f"Attempt {attempt + 1}/{attempts} failed with error: {e}. Retrying in {delay} seconds...")
                    time.sleep(delay)
            raise Exception(f"All {attempts} attempts failed.")
        return wrapper
    return decorator

#====================================================

# function
# project ID 2 SRA ID
# input project_id
# return record["IdList"]
@rate_limited(10)
@retry(attempts=3, delay=2)
def fetch_sra_ids(project_id):
    try:
        print(f"Fetching data for {project_id} at {time.strftime('%X')}")
        handle = Entrez.esearch(db="sra", term=project_id, retmax=10000)
        record = Entrez.read(handle)
        handle.close()
        return record["IdList"]
    except Exception as e:
        print(f"Error fetching SRA IDs for project {project_id}: {e}")
        raise  

# function
# SRA ID to SRR ID(accession)
# return list
@rate_limited(10)
@retry(attempts=3, delay=2)
def get_srr_ids_from_sra(sra_id):
    try:
        print(f"Fetching data for {sra_id} at {time.strftime('%X')}")
        handle = Entrez.efetch(db="sra", id=sra_id)
        data = handle.read()
        handle.close()
        srr_ids = []
        root = ET.fromstring(data)
        for run in root.findall('.//RUN_SET/RUN'):
            srr_id = run.attrib['accession']
            srr_ids.append(srr_id)
        return srr_ids  # Added return statement
    except Exception as e:
        print(f"Error fetching SRA IDs for project {sra_id}: {e}")
        raise

# function
# Using sratoolkit::vdb-dump --info, to query files format via SRR IDs 
# Input srr_id
# Input SRR IDs
# Return line.split(':')[1].strip()
# Return files' format
@rate_limited(10)
@retry(attempts=3, delay=2)
def get_sra_info(srr_id):
    info = {"format": "Unknown", "size": 0} 
    try:
        result = subprocess.run(['vdb-dump', '--info', srr_id], capture_output=True, text=True)
        for line in result.stdout.split('\n'):
            
            fmt_match = re.search(r"FMT\s+:\s+(.+)", line)
            size_match = re.search(r"size\s+:\s+([\d,]+)", line)
            if fmt_match:
                info['format'] = fmt_match.group(1)
            if size_match:
                size_str = size_match.group(1).replace(',', '')
                try:
                    info['size'] = int(size_str)  
                except ValueError:
                    print(f"Error converting size for SRR ID {srr_id}: {size_str} is not a valid integer.")
                    info['size'] = 0  
    except Exception as e:
        print(f"Error fetching info for SRR ID {srr_id}: {e}")
    raise  


# function
# Writing SRA ID to CSV files
# Input (*)email, (*)api_key, input_csv_path, intermediate_csv_path
#   (*) The api_key may acclerate your query, but it not essential.
#   (*) The email tells NCBI who are you.
# Return intermediate_csv_path
# Return CSV files, contais Project IDs and SRA IDs. Some data may be lost depending on the network？It's weird.
def write_sra_ids_to_csv(email, api_key, input_csv_path, intermediate_csv_path):
    Entrez.email = email
    Entrez.api_key = api_key

    with open(input_csv_path, newline='') as csvfile, \
         open(intermediate_csv_path, 'w', newline='') as outcsv:
        reader = csv.DictReader(csvfile)  # 使用DictReader而不是reader
        writer = csv.writer(outcsv)
        writer.writerow(['Project ID', 'SRA ID'])

        for row in tqdm(reader, desc="Fetching SRA IDs"):
            project_id = row.get('Project ID')  
            if project_id: 
                project_ids = fetch_sra_ids(project_id)
                for sra_id in project_ids:
                    writer.writerow([project_id, sra_id])

# function
# Writing SRR IDs to CSV files.
# Input intermediate_csv_path, output_csv_path, view
# Return output_csv_path
# Return output CSV files, contais Project IDs and SRA IDs and SRR IDs.
def update_csv_with_srr_id(intermediate_csv_path, output_csv_path, view):
    with open(intermediate_csv_path, newline='') as incsv, open(output_csv_path, 'w', newline='') as outcsv:
        reader = csv.DictReader(incsv)
        fieldnames = ['Project ID', 'SRA ID', 'SRR ID'] 
        writer = csv.DictWriter(outcsv, fieldnames=fieldnames)
        writer.writeheader()

        sra_ids = [row for row in reader]

        progress = tqdm(total=len(sra_ids), desc="Extracting SRR IDs")

      
        with concurrent.futures.ThreadPoolExecutor(max_workers=1) as executor:
            future_to_row = {executor.submit(get_srr_ids_from_sra, row['SRA ID']): row for row in sra_ids}
            
            for future in concurrent.futures.as_completed(future_to_row):
                row = future_to_row[future]
                srr_ids = future.result()
                
                for srr_id in srr_ids:
                    new_row = {'Project ID': row['Project ID'], 'SRA ID': row['SRA ID'], 'SRR ID': srr_id}
                    writer.writerow(new_row)
                    
                progress.update(1)
                if view:
                    tqdm.write(f"Updated: {row['Project ID']}, {row['SRA ID']}, SRR IDs: {', '.join(srr_ids)}")
        
        progress.close()


# function
# Input output_csv_path, final_csv_path
# Input output CSV files, and new output files' path
# Return final_csv_path
# Return new output CSV files
def update_csv_with_file_format(output_csv_path, final_csv_path):
    with open(output_csv_path, newline='') as incsv, open(final_csv_path, 'w', newline='') as outcsv:
        reader = csv.DictReader(incsv)
        fieldnames = reader.fieldnames + ['File Format'] 
        writer = csv.DictWriter(outcsv, fieldnames=fieldnames)
        writer.writeheader()

        rows = list(reader)
        progress = tqdm(total=len(rows), desc="Fetching file formats")

        with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
            future_to_row = {executor.submit(get_sra_format, row['SRR ID']): row for row in rows}
            
            for future in concurrent.futures.as_completed(future_to_row):
                row = future_to_row[future]
                file_format = future.result()
                row['File Format'] = file_format
                writer.writerow(row)
                progress.update(1)
        
        progress.close()

# Main function
def main(email, api_key, input_csv_path, output_dir, view):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    intermediate_csv_path = os.path.join(output_dir, "intermediate_sra_ids.csv")
    final_csv_path = os.path.join(output_dir, "final_output.csv")

    write_sra_ids_to_csv(email, api_key, input_csv_path, intermediate_csv_path)
    output_csv_path = os.path.join(output_dir, "output_with_srr_ids.csv")
    update_csv_with_srr_id(intermediate_csv_path, output_csv_path, view)
    update_csv_with_file_format(output_csv_path, final_csv_path)

# Options
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch SRA and SRR IDs for given project IDs.")
    parser.add_argument('-e', '--email', required=True, help="Your email address.")
    parser.add_argument('-k', '--api_key', default='', help="Your NCBI API key (optional).")
    parser.add_argument('-i', '--input_csv', required=True, help="Path to the input CSV file containing project IDs.")
    parser.add_argument('-o', '--output_dir', required=True, help="Directory to store the output CSV files.")
    parser.add_argument('-v', '--view', action='store_true', help="Display detailed progress information.")

    args = parser.parse_args()

    main(args.email, args.api_key, args.input_csv, args.output_dir, args.view)

