import csv
from Bio import Entrez
from tqdm import tqdm
import argparse

def fetch_abstract_and_title(email, api_key, pmid):
    Entrez.email = email
    Entrez.api_key = api_key
    try:
        handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        title = records['PubmedArticle'][0]['MedlineCitation']['Article']['ArticleTitle']
        abstract_elements = records['PubmedArticle'][0]['MedlineCitation']['Article']['Abstract']['AbstractText']
        # 初始化摘要文本为空字符串
        abstract_text = ""
        # 遍历所有摘要元素
        for element in abstract_elements:
            # 检查摘要元素是否有Label属性，如果有，加上Label并换行
            if 'Label' in element.attributes:
                abstract_text += f"{element.attributes['Label']}: "
            abstract_text += element + " "  # 添加摘要文本和空格
        return title, abstract_text.strip(), True  # 返回标题、整理好的摘要文本，并标记成功
    except Exception as e:  # 处理其他错误，如请求失败等
        return "Request failed", str(e), False


def process_pmids(input_csv, output_csv, email, api_key, view):
    success_count = 0
    fail_count = 0
    with open(input_csv, newline='') as infile, open(output_csv, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = ['PMID', 'Title', 'Abstract']  # 定义输出CSV的列名
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()

        for row in tqdm(reader, desc="Processing PMIDs"):
            pmid = row['PMID']
            title, abstract, success = fetch_abstract_and_title(email, api_key, pmid)
            writer.writerow({'PMID': pmid, 'Title': title, 'Abstract': abstract})
            if success:
                success_count += 1
            else:
                fail_count += 1
            if view:  # 实时显示
                print(f"PMID: {pmid} | Status: {'Success' if success else 'Failed'}")

    print(f"Process completed. Success: {success_count}, Failed: {fail_count}")

# Main function
def main(input_csv, output_csv, email, api_key, view):
    process_pmids(input_csv, output_csv, email, api_key, view)

# Options
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch PubMed abstracts and titles for given PMIDs.")
    parser.add_argument('-e', '--email', required=True, help="Your email address.")
    parser.add_argument('-k', '--api_key', default='', help="Your NCBI API key (optional).")
    parser.add_argument('-i', '--input_csv', required=True, help="Path to the input CSV file containing PMIDs.")
    parser.add_argument('-o', '--output_csv', required=True, help="Path to store the output CSV file with titles and abstracts.")
    parser.add_argument('-v', '--view', action='store_true', help="Display detailed progress information.")

    args = parser.parse_args()

    main(args.input_csv, args.output_csv, args.email, args.api_key, args.view)
