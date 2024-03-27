import csv
from Bio import Entrez
from tqdm import tqdm
import os

# 配置你的Entrez访问
Entrez.email = "YOU_email"
Entrez.api_key = "YOU_apikey“

def check_open_access(pmid):
    try:
        #links = Entrez.elink(dbfrom="pubmed", db="pmc", id=pmid, linkname="pubmed_pmc_refs")
        links = Entrez.elink(dbfrom="pubmed", db="pmc", id=pmid)
        record = Entrez.read(links)
        if record[0]["LinkSetDb"]:
            pmc_id = record[0]["LinkSetDb"][0]["Link"][0]["Id"]
            return "Open Access", f"PMC{pmc_id}"  # 返回访问状态和PMC ID
        else:
            return "Paper not open", None  # 文献不是开放获取
    except Exception as e:
        return "Error checking access", None  # 检查访问状态时发生错误

def process_csv(input_csv, output_csv, output_dir):
    with open(input_csv, newline='') as infile, open(output_csv, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames + ['Access', 'PMC_ID', 'FullTextSaved', 'FullTextStatus']  # Add 'FullTextStatus'
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()

        for row in tqdm(reader, desc="Processing PMIDs"):
            pmid = row.get('PMID')
            if pmid:
                access_status, pmc_id = check_open_access(pmid)
                row['Access'] = access_status
                row['PMC_ID'] = pmc_id if pmc_id else "N/A"
                if pmc_id and pmc_id.startswith("PMC"):
                    success, status_message = fetch_full_text(pmc_id, output_dir)
                    row['FullTextSaved'] = "Yes" if success else "No"
                    row['FullTextStatus'] = status_message  # Update the status message
                else:
                    row['FullTextSaved'] = "No"
                    row['FullTextStatus'] = "PMC ID not found"  # Or any other appropriate message
                writer.writerow(row)




def fetch_full_text(pmc_id, output_dir, retries=3, delay=5):
    """
    Fetches the full text of a given PMC ID in XML format and saves it to a file,
    with specified retries upon failure.

    :param pmc_id: The PMC ID of the paper.
    :param output_dir: Directory where the full text XML will be saved.
    :param retries: Number of retries if fetching fails.
    :param delay: Delay between retries in seconds.
    :return: A tuple of (bool, str) indicating success/failure and status message.
    """
    attempt = 0
    while attempt < retries:
        try:
            handle = Entrez.efetch(db="pmc", id=pmc_id, retmode="xml")
            xml_data = handle.read()  # This is in bytes
            handle.close()

            if b"The publisher of this article does not allow downloading of the full text in XML form" in xml_data:
                return False, "Full text not allowed"

            # Build the output file path
            output_path = os.path.join(output_dir, f"{pmc_id}.xml")

            # Save the full text XML to a file in binary write mode
            with open(output_path, 'wb') as f:
                f.write(xml_data)
            return True, f"Saved full text for {pmc_id} to {output_path}"
        except Exception as e:
            print(f"Attempt {attempt + 1} failed: {e}")
            time.sleep(delay)  # Wait for a bit before retrying
            attempt += 1

    return False, f"Failed to fetch full text for {pmc_id} after {retries} attempts."



# 输入和输出CSV文件的路径
input_csv = "/content/PMID_list.csv"
output_csv = "/content/PMID_list_XML_paper_supple.csv"
output_dir = "/content/drive/MyDrive/request_list_XML_5"

process_csv(input_csv, output_csv, output_dir)
