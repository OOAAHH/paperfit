# rna_structure_workflow.py
"""Enhanced pipeline for triaging RNA‑structure–prediction papers
---------------------------------------------------------------
This version adds
1. **Metadata re‑use** – if the input CSV already contains *title*, *abstract*,
   or *doi*, we use them directly; PubMed fetch becomes a fallback.
2. **Cheap heuristic filter** – obvious non‑algorithm papers are discarded
   without an LLM call, saving tokens.
3. **Citation retrieval** – retrieves citations for each article and stores them separately.
4. **LLM-based analysis** – uses LLM to analyze title and abstract only, providing detailed reasoning.
5. **Web search enhancement** - uses local search implementation to gather additional context.

Run:
    python rna_structure_workflow.py articles.csv output.csv

Env vars (new):
    ENTREZ_EMAIL     – for NCBI Entrez
    OPENAI_API_KEY   – OpenAI key
    S2_API_KEY       – Semantic Scholar token (optional)
    OPENAI_BASE_URL  - 自定义API基础URL，如DeepSeek API
    USE_LOCAL_SEARCH - 使用本地搜索方案而非SerpAPI
"""
from __future__ import annotations
import os, re, sys, time, json, argparse, logging
from typing import Optional, Dict, Any, List, Tuple
import requests, pandas as pd
from tqdm import tqdm
from xml.etree import ElementTree as ET
from Bio import Entrez
from openai import OpenAI
from bs4 import BeautifulSoup
import urllib.parse

# ---------------- configuration ---------------- #
MODEL_NAME           = os.getenv("OPENAI_MODEL", "deepseek-reasoner")  # 使用deepseek-reasoner模型
OPENAI_TEMP          = float(os.getenv("OPENAI_TEMP", 0))  # 设置温度为0以保证结果一致性
RATE_LIMIT_PER_SEC   = 3          # 对NCBI请求的礼貌延迟，每秒最多3个请求
S2_API_KEY           = os.getenv("S2_API_KEY", "")  # Semantic Scholar API密钥
OPENAI_BASE_URL      = os.getenv("OPENAI_BASE_URL", "https://api.deepseek.com")  # DeepSeek API基础URL
SERPAPI_API_KEY      = os.getenv("SERPAPI_API_KEY", "")  # SerpAPI搜索引擎API密钥
USE_WEB_SEARCH       = os.getenv("USE_WEB_SEARCH", "false").lower() == "true"  # 是否启用网络搜索
USE_LOCAL_SEARCH     = os.getenv("USE_LOCAL_SEARCH", "false").lower() == "true"  # 是否使用本地搜索方案
LOCAL_EMBEDDINGS     = os.getenv("USE_LOCAL_EMBEDDINGS", "false").lower() == "true"  # 是否使用本地嵌入
DEBUG                = os.getenv("DEBUG", "false").lower() == "true"  # 是否启用详细调试模式

# 配置日志格式，使用INFO级别，如果DEBUG模式则使用DEBUG级别
logging.basicConfig(
    level=logging.DEBUG if DEBUG else logging.INFO, 
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)

# 设置NCBI Entrez邮箱（必需）
Entrez.email = os.getenv("ENTREZ_EMAIL", "oahsun@outlook.com")
if Entrez.email == "oahsun@outlook.com":
    logging.warning("Set ENTREZ_EMAIL to comply with NCBI policies.")

# ---------------- simple heuristics ---------------- #
# 定义正向关键词正则表达式，用于识别可能是结构预测方法的文章
KEY_POS = re.compile(r"\b(predict(?:ion|ive)?|model(?:ing)?|fold(?:ing)?)\b", re.I)
# 定义负向关键词正则表达式，用于排除明显不是算法的文章
KEY_NEG = re.compile(r"\b(review|database|protocol|benchmark|assessment|case study)\b", re.I)


def heuristic_is_structure_method(title: str, abstract: str) -> Tuple[bool, str]:
    """
    使用简单启发式规则进行快速筛选
    
    参数:
        title: 文章标题
        abstract: 文章摘要
        
    返回:
        Tuple[bool, str]: 
            - 如果包含正向关键词且不包含明显的负向关键词，返回(True, 原因)
            - 否则返回(False, 原因)
    """
    text = f"{title} {abstract}"
    has_pos = bool(KEY_POS.search(text))
    has_neg = bool(KEY_NEG.search(text))
    
    # 查找匹配的关键词，用于解释
    pos_matches = [m.group() for m in KEY_POS.finditer(text)]
    neg_matches = [m.group() for m in KEY_NEG.finditer(text)]
    
    # 构建解释
    if has_pos and not has_neg:
        return True, f"包含预测关键词: {', '.join(pos_matches)}"
    elif not has_pos:
        return False, "不包含预测相关关键词"
    else:  # has_pos and has_neg
        return False, f"虽然包含预测关键词({', '.join(pos_matches)})，但也包含排除关键词({', '.join(neg_matches)})"


# ---------------- Web Search Utilities ---------------- #

def local_web_search(query: str, num_results: int = 3) -> str:
    """
    使用本地化方案进行网络搜索，替代SerpAPI
    
    参数:
        query: 搜索查询词
        num_results: 返回结果数量
        
    返回:
        str: 搜索结果摘要
    """
    try:
        logging.debug(f"使用本地搜索方案，查询: {query}")
        
        # 对查询进行URL编码
        encoded_query = urllib.parse.quote(query)
        
        # 使用Google Scholar搜索
        url = f"https://scholar.google.com/scholar?q={encoded_query}&hl=en&as_sdt=0,5&num={num_results}"
        
        # 设置headers模拟浏览器
        headers = {
            "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36",
            "Accept-Language": "en-US,en;q=0.9",
            "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8",
            "Connection": "keep-alive",
            "Accept-Encoding": "gzip, deflate, br",
        }
        
        # 发送请求
        response = requests.get(url, headers=headers, timeout=20)
        
        if not response.ok:
            logging.warning(f"Google Scholar搜索返回错误: {response.status_code}")
            
            # 尝试使用普通Google搜索作为备选
            backup_url = f"https://www.google.com/search?q={encoded_query}&num={num_results}"
            response = requests.get(backup_url, headers=headers, timeout=20)
            
            if not response.ok:
                logging.warning(f"备用搜索也失败: {response.status_code}")
                return ""
        
        # 解析HTML
        soup = BeautifulSoup(response.text, 'html.parser')
        
        # 提取搜索结果
        search_results = []
        
        # 对于Google Scholar
        results = soup.select('.gs_ri')
        if results:
            for i, result in enumerate(results[:num_results]):
                title_elem = result.select_one('.gs_rt')
                snippet_elem = result.select_one('.gs_rs')
                
                title = title_elem.text if title_elem else ""
                snippet = snippet_elem.text if snippet_elem else ""
                
                # 获取链接
                link = ""
                if title_elem and title_elem.a:
                    link = title_elem.a.get('href', '')
                
                search_results.append(f"标题: {title}\n摘要: {snippet}\n链接: {link}\n")
        else:
            # 对于普通Google搜索
            results = soup.select('.g')
            for i, result in enumerate(results[:num_results]):
                title_elem = result.select_one('h3')
                snippet_elem = result.select_one('.VwiC3b')
                
                title = title_elem.text if title_elem else ""
                snippet = snippet_elem.text if snippet_elem else ""
                
                # 获取链接
                link = ""
                link_elem = result.select_one('a')
                if link_elem:
                    link = link_elem.get('href', '')
                    if link.startswith('/url?q='):
                        link = link[7:link.find('&sa=')]
                
                search_results.append(f"标题: {title}\n摘要: {snippet}\n链接: {link}\n")
        
        # 将所有结果合并为一个字符串
        combined_results = "\n".join(search_results)
        logging.debug(f"本地搜索返回{len(search_results)}条结果")
        
        # 等待一段时间，避免过快请求被封IP
        time.sleep(3)
        
        return combined_results
            
    except Exception as e:
        logging.warning(f"本地网络搜索失败: {e}")
        return ""


def search_semantic_scholar(query: str, num_results: int = 3) -> str:
    """
    使用Semantic Scholar API搜索学术文献
    
    参数:
        query: 搜索查询词
        num_results: 返回结果数量
        
    返回:
        str: 搜索结果摘要
    """
    try:
        logging.debug(f"使用Semantic Scholar搜索: {query}")
        
        # Semantic Scholar API的搜索端点
        url = "https://api.semanticscholar.org/graph/v1/paper/search"
        
        # 设置查询参数
        params = {
            "query": query,
            "limit": num_results,
            "fields": "title,abstract,url,year,authors,venue"
        }
        
        # 如果有API密钥，添加到headers中
        headers = {}
        if S2_API_KEY:
            headers["x-api-key"] = S2_API_KEY
        
        # 发送请求
        response = requests.get(url, params=params, headers=headers, timeout=20)
        
        if not response.ok:
            logging.warning(f"Semantic Scholar API返回错误: {response.status_code}")
            return ""
        
        # 解析结果
        data = response.json()
        search_results = []
        
        for paper in data.get('data', []):
            title = paper.get('title', '')
            abstract = paper.get('abstract', '')
            paper_url = paper.get('url', '')
            year = paper.get('year', '')
            
            # 获取作者信息
            authors = []
            for author in paper.get('authors', []):
                authors.append(author.get('name', ''))
            
            author_str = ', '.join(authors[:3])
            if len(authors) > 3:
                author_str += ' et al.'
                
            venue = paper.get('venue', '')
            
            result_str = f"标题: {title}\n"
            if abstract:
                result_str += f"摘要: {abstract}\n"
            result_str += f"作者: {author_str}\n"
            if year:
                result_str += f"年份: {year}\n"
            if venue:
                result_str += f"期刊/会议: {venue}\n"
            result_str += f"链接: {paper_url}\n"
            
            search_results.append(result_str)
        
        # 将所有结果合并为一个字符串
        combined_results = "\n".join(search_results)
        logging.debug(f"Semantic Scholar搜索返回{len(search_results)}条结果")
        
        return combined_results
        
    except Exception as e:
        logging.warning(f"Semantic Scholar搜索失败: {e}")
        return ""


def search_pubmed_articles(query: str, num_results: int = 3) -> str:
    """
    使用PubMed API搜索医学文献
    
    参数:
        query: 搜索查询词
        num_results: 返回结果数量
        
    返回:
        str: 搜索结果摘要
    """
    try:
        logging.debug(f"使用PubMed搜索: {query}")
        
        # 使用Entrez搜索PubMed
        handle = Entrez.esearch(db="pubmed", term=query, retmax=num_results)
        record = Entrez.read(handle)
        handle.close()
        
        if not record.get("IdList"):
            logging.warning("PubMed搜索未返回结果")
            return ""
        
        # 获取文章详情
        ids = record["IdList"]
        handle = Entrez.efetch(db="pubmed", id=",".join(ids), rettype="abstract", retmode="xml")
        articles = Entrez.read(handle)
        handle.close()
        
        search_results = []
        
        # 处理每篇文章
        for article in articles["PubmedArticle"]:
            try:
                # 获取文章信息
                medline_citation = article["MedlineCitation"]
                article_data = medline_citation["Article"]
                
                # 标题
                title = article_data.get("ArticleTitle", "")
                
                # 摘要
                abstract = ""
                if "Abstract" in article_data:
                    abstract_texts = article_data["Abstract"].get("AbstractText", [])
                    if isinstance(abstract_texts, list):
                        abstract = " ".join([str(text) for text in abstract_texts])
                    else:
                        abstract = str(abstract_texts)
                
                # 作者
                authors = []
                if "AuthorList" in article_data:
                    for author in article_data["AuthorList"]:
                        if "LastName" in author and "ForeName" in author:
                            authors.append(f"{author['LastName']} {author['ForeName']}")
                
                # 期刊信息
                journal = article_data.get("Journal", {})
                journal_title = journal.get("Title", "")
                
                # 发布年份
                pub_date = None
                if "PubDate" in journal.get("JournalIssue", {}):
                    pub_date = journal["JournalIssue"]["PubDate"]
                    if "Year" in pub_date:
                        pub_date = pub_date["Year"]
                    else:
                        pub_date = None
                
                # DOI
                doi = ""
                for id_elem in article_data.get("ELocationID", []):
                    if id_elem.attributes["EIdType"] == "doi":
                        doi = str(id_elem)
                        break
                
                # 构建结果字符串
                result_str = f"标题: {title}\n"
                if abstract:
                    result_str += f"摘要: {abstract}\n"
                if authors:
                    author_str = ", ".join(authors[:3])
                    if len(authors) > 3:
                        author_str += " et al."
                    result_str += f"作者: {author_str}\n"
                if pub_date:
                    result_str += f"年份: {pub_date}\n"
                if journal_title:
                    result_str += f"期刊: {journal_title}\n"
                if doi:
                    result_str += f"DOI: {doi}\n"
                
                result_str += f"PMID: {medline_citation.get('PMID', '')}\n"
                
                search_results.append(result_str)
                
            except Exception as e:
                logging.warning(f"处理PubMed文章时出错: {e}")
                continue
        
        # 合并结果
        combined_results = "\n".join(search_results)
        logging.debug(f"PubMed搜索返回{len(search_results)}条结果")
        
        return combined_results
        
    except Exception as e:
        logging.warning(f"PubMed搜索失败: {e}")
        return ""


def search_web(query: str, num_results: int = 3) -> str:
    """
    使用搜索引擎API搜索相关信息
    
    参数:
        query: 搜索查询词
        num_results: 返回结果数量
        
    返回:
        str: 搜索结果摘要
    """
    # 如果使用本地搜索方案，则不使用SerpAPI
    if USE_LOCAL_SEARCH:
        # 首先尝试PubMed搜索，因为我们主要处理生物医学文献
        pubmed_results = search_pubmed_articles(query, num_results)
        
        # 如果PubMed搜索成功，则返回结果
        if pubmed_results:
            return pubmed_results
            
        # 否则尝试Semantic Scholar搜索
        semantic_results = search_semantic_scholar(query, num_results)
        if semantic_results:
            return semantic_results
            
        # 如果前两种方法都失败，则使用本地网络搜索
        return local_web_search(query, num_results)
    
    # 如果不使用本地搜索方案，则使用SerpAPI
    if not SERPAPI_API_KEY:
        logging.warning("未设置SERPAPI_API_KEY，无法进行网络搜索")
        return ""
        
    try:
        logging.debug(f"搜索查询: {query}")
        
        # 使用SerpAPI搜索
        params = {
            "engine": "google",
            "q": query,
            "api_key": SERPAPI_API_KEY,
            "num": num_results
        }
        
        response = requests.get(
            "https://serpapi.com/search", 
            params=params,
            timeout=30
        )
        
        if not response.ok:
            logging.warning(f"搜索API返回错误: {response.status_code}")
            return ""
            
        results = response.json()
        
        # 提取搜索结果摘要
        search_results = []
        organic_results = results.get("organic_results", [])
        
        for result in organic_results[:num_results]:
            title = result.get("title", "")
            snippet = result.get("snippet", "")
            link = result.get("link", "")
            search_results.append(f"标题: {title}\n摘要: {snippet}\n链接: {link}\n")
        
        # 将所有结果合并为一个字符串
        combined_results = "\n".join(search_results)
        logging.debug(f"搜索返回{len(search_results)}条结果")
        
        return combined_results
            
    except Exception as e:
        logging.warning(f"网络搜索失败: {e}")
        return ""


def get_enhanced_context(title: str, abstract: str, doi: str = "", pmid: str = "") -> str:
    """
    获取增强的上下文信息，通过网络搜索相关内容
    
    参数:
        title: 文章标题
        abstract: 文章摘要
        doi: 文章DOI
        pmid: PubMed ID
        
    返回:
        str: 增强的上下文信息
    """
    if not USE_WEB_SEARCH:
        return ""
        
    # 构建搜索查询
    # 使用文章标题+RNA structure prediction关键词
    query = f"{title} RNA structure prediction algorithm method"
    
    # 如果有DOI，添加到查询中
    if doi:
        query += f" {doi}"
        
    # 获取搜索结果
    search_results = search_web(query)
    
    return search_results


# --------------- Entrez utilities --------------- #

def efetch_xml(db: str, id_: str, rettype: str = "abstract", retmode: str = "xml") -> ET.Element:
    """
    从NCBI的Entrez数据库获取XML格式数据，带有重试和延迟
    
    参数:
        db: 数据库名称（如'pubmed'）
        id_: 记录ID
        rettype: 返回类型，默认为'abstract'
        retmode: 返回模式，默认为'xml'
        
    返回:
        ET.Element: 解析后的XML元素
    """
    for attempt in range(3):
        try:
            logging.debug(f"尝试从{db}获取ID:{id_}的数据 (尝试 {attempt+1}/3)")
            with Entrez.efetch(db=db, id=id_, rettype=rettype, retmode=retmode) as h:
                xml_text = h.read()
            time.sleep(1 / RATE_LIMIT_PER_SEC)
            return ET.fromstring(xml_text)
        except Exception as e:
            logging.warning(f"获取{db}中ID:{id_}失败: {e}")
            if attempt == 2:
                raise
            time.sleep(2 ** attempt)
    return ET.Element("empty")


def fetch_pubmed_record(pmid: str) -> Dict[str, Any]:
    """
    从PubMed获取文章的基本元数据
    
    参数:
        pmid: PubMed ID
        
    返回:
        Dict: 包含标题、摘要和DOI的字典
    """
    logging.debug(f"获取PMID:{pmid}的PubMed记录")
    root = efetch_xml("pubmed", pmid)
    article = root.find(".//PubmedArticle/MedlineCitation/Article")
    if article is None:
        logging.warning(f"PMID:{pmid}没有找到文章信息")
        return {}
    title = " ".join(article.findtext("ArticleTitle", default="").split())
    abstract = " ".join(" ".join(p.text or "" for p in article.findall("Abstract/AbstractText")).split())
    doi = None
    for id_node in article.findall("ELocationID"):
        if id_node.get("EIdType") == "doi":
            doi = id_node.text
            break
    
    logging.debug(f"PMID:{pmid}获取成功 - 标题: {title[:30]}... DOI: {doi}")
    return {"title": title, "abstract": abstract, "doi": doi}


def is_valid_doi(doi: str) -> bool:
    """
    检查DOI是否有效
    
    参数:
        doi: 待检查的DOI字符串
        
    返回:
        bool: 如果DOI格式有效返回True
    """
    # 基本的DOI格式检查 - 大多数DOI以10.开头
    if not doi or not isinstance(doi, str):
        return False
        
    doi = doi.strip()
    # 检查基本格式
    if not doi.startswith('10.'):
        return False
        
    # 检查是否有合理的长度和格式
    if '/' not in doi or len(doi) < 8:
        return False
        
    return True


# ---------------- Citation retrieval functions ---------------- #

def get_pubmed_references(pmid: str) -> List[Dict[str, str]]:
    """
    获取PubMed文章引用的参考文献列表
    
    参数:
        pmid: PubMed ID
        
    返回:
        List[Dict[str, str]]: 参考文献列表，每个引文包含pmid, title, year等信息
    """
    references = []
    try:
        logging.debug(f"获取PMID:{pmid}的参考文献")
        # 使用elink获取该文章引用的文献（参考文献）
        with Entrez.elink(dbfrom="pubmed", id=pmid, linkname="pubmed_pubmed_refs") as h:
            data = h.read()
        root = ET.fromstring(data)
        
        # 提取参考文献的PMID列表
        ref_pmids = []
        for link in root.findall(".//LinkSetDb/Link/Id"):
            ref_pmids.append(link.text)
        
        if not ref_pmids:
            logging.debug(f"PMID:{pmid}没有找到参考文献")
            return []
            
        # 获取这些参考文献的基本信息
        logging.debug(f"找到{len(ref_pmids)}篇PMID:{pmid}引用的参考文献，获取详情")
        # 将PMID列表分批处理，避免请求过大
        batch_size = 50
        for i in range(0, len(ref_pmids), batch_size):
            batch = ref_pmids[i:i+batch_size]
            with Entrez.efetch(db="pubmed", id=",".join(batch), rettype="abstract", retmode="xml") as h:
                reference_data = h.read()
            reference_root = ET.fromstring(reference_data)
            
            # 处理每篇参考文献
            for article_elem in reference_root.findall(".//PubmedArticle"):
                try:
                    ref_pmid = article_elem.find(".//PMID").text
                    article = article_elem.find(".//Article")
                    if article is None:
                        continue
                        
                    title = article.findtext("ArticleTitle", "")
                    year = ""
                    date_elem = article_elem.find(".//PubDate/Year")
                    if date_elem is not None:
                        year = date_elem.text
                        
                    journal = article.findtext("Journal/Title", "")
                    
                    # 添加到参考文献列表
                    references.append({
                        "pmid": ref_pmid,
                        "title": title,
                        "year": year,
                        "journal": journal
                    })
                except Exception as e:
                    logging.warning(f"处理参考文献PMID时出错: {e}")
            
            # 添加延迟，避免过快请求
            time.sleep(1 / RATE_LIMIT_PER_SEC)
            
        logging.debug(f"成功获取PMID:{pmid}的{len(references)}篇参考文献详情")
    except Exception as e:
        logging.warning(f"获取参考文献失败: {e}")
    
    return references


def get_semantic_scholar_references(doi_or_pmid: str) -> List[Dict[str, str]]:
    """
    使用Semantic Scholar API获取文章的参考文献
    
    参数:
        doi_or_pmid: 文章的DOI或PubMed ID
        
    返回:
        List[Dict[str, str]]: 参考文献列表
    """
    if not S2_API_KEY:
        return []
        
    references = []
    try:
        logging.debug(f"通过Semantic Scholar获取{doi_or_pmid}的参考文献")
        
        # 确定是DOI还是PMID
        id_type = "DOI" if is_valid_doi(doi_or_pmid) else "PMID"
        paper_id = f"{id_type}:{doi_or_pmid}"
        
        # 构建API请求 - 注意使用references而不是citations
        headers = {"x-api-key": S2_API_KEY}
        url = f"https://api.semanticscholar.org/graph/v1/paper/{paper_id}/references?fields=title,year,journal,authors,externalIds"
        
        response = requests.get(url, headers=headers, timeout=20)
        if not response.ok:
            logging.warning(f"Semantic Scholar API返回错误: {response.status_code}")
            return []
            
        data = response.json()
        for item in data.get("data", []):
            try:
                # 在references接口中，引用的论文在citedPaper字段
                ref_data = item.get("citedPaper", {})
                if not ref_data:
                    continue
                    
                external_ids = ref_data.get("externalIds", {})
                reference = {
                    "title": ref_data.get("title", ""),
                    "year": str(ref_data.get("year", "")),
                    "journal": ref_data.get("journal", {}).get("name", ""),
                    "s2_id": ref_data.get("paperId", ""),
                    "pmid": external_ids.get("PubMed", ""),
                    "doi": external_ids.get("DOI", "")
                }
                references.append(reference)
            except Exception as e:
                logging.warning(f"处理Semantic Scholar参考文献项时出错: {e}")
                
        logging.debug(f"从Semantic Scholar获取到{len(references)}篇参考文献")
    except Exception as e:
        logging.warning(f"Semantic Scholar参考文献获取失败: {e}")
        
    return references


def get_article_references(pmid: str, doi: Optional[str] = None) -> List[Dict[str, str]]:
    """
    获取文章的参考文献，优先使用Semantic Scholar，失败则使用PubMed
    
    参数:
        pmid: PubMed ID
        doi: 文章DOI
        
    返回:
        List[Dict[str, str]]: 参考文献列表
    """
    # 先尝试使用Semantic Scholar（如果有API密钥）
    references = []
    if S2_API_KEY:
        # 优先使用DOI查询
        if doi and is_valid_doi(doi):
            references = get_semantic_scholar_references(doi)
        
        # 如果使用DOI失败，尝试使用PMID
        if not references and pmid:
            references = get_semantic_scholar_references(pmid)
    
    # 如果Semantic Scholar未获取到参考文献，使用PubMed
    if not references and pmid:
        references = get_pubmed_references(pmid)
    
    return references


# ---------------- OpenAI utilities ---------------- #
# 初始化OpenAI客户端，使用新版API
api_key = os.getenv("OPENAI_API_KEY")
if not api_key:
    sys.exit("Set OPENAI_API_KEY first.")

try:
    client = OpenAI(
        api_key=api_key,
        base_url=OPENAI_BASE_URL,  # 使用DeepSeek API基础URL
    )
    logging.info(f"成功初始化AI客户端，使用模型: {MODEL_NAME}，API基础URL: {OPENAI_BASE_URL}")
except Exception as e:
    logging.error(f"初始化AI客户端失败: {e}")
    sys.exit(f"初始化AI客户端失败: {e}")

# 定义系统提示，明确要求模型只回答"yes"或"no"并提供原因
SYSTEM_PROMPT = (
    "You are a biomedical literature triage assistant."
    " Given a TITLE and ABSTRACT, analyze if the paper describes an RNA structure prediction method."
    " Answer in JSON format: {\"decision\": \"yes/no\", \"reasoning\": \"your detailed explanation\"}"
    " Answer \"yes\" **only** if the paper presents an ALGORITHM, METHOD, or computational PIPELINE whose goals include PREDICTING RNA secondary or tertiary (3D) STRUCTURE."
    " Reviews, surveys, databases, experimental protocols, purely protein‑structure methods, or papers without a predictive component => \"no\"."
    " Be specific in your reasoning, mentioning key terms or concepts that influenced your decision."
    " Your reasoning should be comprehensive and explain why this is or is not an RNA structure prediction method."
)

# 增强版系统提示，包含网络搜索结果
ENHANCED_SYSTEM_PROMPT = (
    "You are a biomedical literature triage assistant with access to additional web search results."
    " Given a TITLE, ABSTRACT, and WEB SEARCH RESULTS, analyze if the paper describes an RNA structure prediction method."
    " Use the additional context from web search to make a more informed decision."
    " Answer in JSON format: {\"decision\": \"yes/no\", \"reasoning\": \"your detailed explanation\", \"search_insights\": \"how search results influenced your decision\"}"
    " Answer \"yes\" **only** if the paper presents an ALGORITHM, METHOD, or computational PIPELINE whose goals include PREDICTING RNA secondary or tertiary (3D) STRUCTURE."
    " Reviews, surveys, databases, experimental protocols, purely protein‑structure methods, or papers without a predictive component => \"no\"."
    " Be specific in your reasoning, mentioning key terms or concepts that influenced your decision."
    " Your reasoning should be comprehensive and explain why this is or is not an RNA structure prediction method."
    " If web search results provide additional relevant information that helped your decision, mention it in the search_insights field."
)


def llm_is_structure_method(title: str, abstract: str, search_results: str = "") -> Tuple[bool, str]:
    """
    使用大语言模型判断文章是否描述RNA结构预测方法
    
    参数:
        title: 文章标题
        abstract: 文章摘要
        search_results: 网络搜索结果（可选）
        
    返回:
        Tuple[bool, str]: (是否是结构预测方法, 原始响应内容)
    """
    try:
        logging.debug(f"发送标题'{title[:30]}...'至AI模型进行分析")
        
        # 根据是否有搜索结果选择不同的提示和内容
        if search_results:
            system_prompt = ENHANCED_SYSTEM_PROMPT
            user_content = f"TITLE: {title}\nABSTRACT: {abstract}\nWEB SEARCH RESULTS: {search_results}"
            logging.debug("使用增强版分析（包含网络搜索结果）")
        else:
            system_prompt = SYSTEM_PROMPT
            user_content = f"TITLE: {title}\nABSTRACT: {abstract}"
            logging.debug("使用标准版分析（仅标题和摘要）")
        
        # 准备消息
        messages = [
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": user_content},
        ]
        
        # 使用DeepSeek API
        response = client.chat.completions.create(
            model=MODEL_NAME,
            messages=messages,
            temperature=OPENAI_TEMP,
            max_tokens=800,  # 增加token上限以获取更详细的解释
            stream=False
        )
        
        content = response.choices[0].message.content.strip()
        logging.debug(f"AI模型响应: {content}")
        
        # 不再尝试解析JSON，而是直接从文本中寻找"yes"或"no"关键词
        is_yes = "yes" in content.lower()
        
        # 返回判断结果和原始响应内容
        return is_yes, content
            
    except Exception as e:
        error_msg = f"AI API错误: {e}"
        logging.error(error_msg)
        return False, error_msg


# ---------------- pipeline per row ---------------- #

def process_article(row: Dict[str, Any], force_llm: bool = False) -> Dict[str, Any]:
    """
    处理单篇文章，执行以下步骤:
    1. 获取/补充文章元数据
    2. 使用启发式规则进行初步筛选
    3. 根据需要调用LLM进行判断
    4. 获取文章参考文献
    
    参数:
        row: 包含文章信息的字典
        force_llm: 是否强制使用LLM判断（即使启发式规则判断为否）
        
    返回:
        Dict[str, Any]: 处理结果
    """
    # 获取基本信息
    pmid = str(row.get("pmid", "") or row.get("PubMed_ID", ""))
    title = row.get("title", "") or row.get("Title", "")
    abstract = row.get("abstract", "") or row.get("Abstract", "")
    doi = row.get("doi", "") or row.get("DOI", "")

    logging.info(f"处理PMID:{pmid}, 标题: {title[:50]}...")
    
    # 数据收集，用于诊断和调试
    diagnostic_info = {
        "pmid": pmid,
        "title_length": len(title),
        "abstract_length": len(abstract),
        "has_doi": bool(doi)
    }

    # 补充缺失字段
    if not title or not abstract or not doi:
        logging.debug(f"PMID:{pmid}的元数据不完整，尝试补充")
        fetched = fetch_pubmed_record(pmid)
        title = title or fetched.get("title", "")
        abstract = abstract or fetched.get("abstract", "")
        doi = doi or fetched.get("doi", "")
        diagnostic_info["fetched_from_pubmed"] = True
        diagnostic_info["pubmed_fetch_success"] = bool(title or abstract)

    # 如果标题和摘要都缺失，无法进行分析
    if not (title or abstract):
        error_msg = f"PMID:{pmid}无法获取标题和摘要，跳过"
        logging.warning(error_msg)
        return {
            "pmid": pmid,
            "error": "no_title_abstract",
            "decision_path": "error",
            "LLM_Response": error_msg
        }

    # 使用启发式规则初步筛选
    heuris_flag, heuris_reason = heuristic_is_structure_method(title, abstract)
    diagnostic_info["heuristic_result"] = heuris_flag
    diagnostic_info["heuristic_reason"] = heuris_reason
    
    logging.debug(f"启发式分析结果: {heuris_flag}, 原因: {heuris_reason}")

    # 根据启发式结果或强制标志决定是否调用LLM
    is_structure = False
    llm_response = heuris_reason
    decision_path = "heuristic_only"
    search_results = ""
    
    if heuris_flag or force_llm:
        # 对于可能的结构预测方法，尝试获取网络搜索结果作为额外上下文
        if USE_WEB_SEARCH:
            logging.debug(f"为PMID:{pmid}获取网络搜索结果")
            search_results = get_enhanced_context(title, abstract, doi, pmid)
            diagnostic_info["web_search_performed"] = True
            diagnostic_info["search_results_length"] = len(search_results)
        
        # 调用LLM进行判断，如果有搜索结果则一并使用
        is_structure, llm_response = llm_is_structure_method(title, abstract, search_results)
        
        # 设置决策路径
        if search_results:
            decision_path = "llm_with_search" if heuris_flag else "force_llm_with_search"
        else:
            decision_path = "llm_after_heuristic" if heuris_flag else "force_llm"
            
        diagnostic_info["llm_called"] = True
        diagnostic_info["llm_result"] = is_structure
        
        logging.debug(f"LLM分析结果: {is_structure}")
    else:
        logging.debug(f"跳过LLM调用，仅使用启发式结果: {heuris_flag}")
        diagnostic_info["llm_called"] = False
    
    # 获取文章参考文献
    references = []
    try:
        if is_structure:  # 只为结构预测方法的文章获取参考文献
            logging.debug(f"获取PMID:{pmid}的参考文献")
            references = get_article_references(pmid, doi)
            diagnostic_info["references_count"] = len(references)
            logging.debug(f"获取到{len(references)}篇参考文献")
    except Exception as e:
        logging.warning(f"获取参考文献时出错: {e}")
        diagnostic_info["reference_error"] = str(e)
    
    # 准备结果
    result = {
        "Structure_Prediction_Method": "Yes" if is_structure else "No",
        "LLM_Response": llm_response,  # 直接保存原始响应
        "Decision_Path": decision_path,
        "References_Count": len(references),
        "References": json.dumps(references) if references else "",
        "Web_Search_Used": "Yes" if search_results else "No",
        "Diagnostic_Info": json.dumps(diagnostic_info) if DEBUG else ""
    }
    
    # 记录处理结果
    status = "结构预测方法" if is_structure else "非结构预测方法"
    search_info = "使用网络搜索" if search_results else "无网络搜索"
    logging.info(f"PMID:{pmid}处理完成: {status}, {search_info}, 获取参考文献: {len(references)}篇")
    
    return result


# ---------------- main CLI ---------------- #

def main(in_csv: str, out_csv: str, force_llm: bool = False, inplace: bool = True):
    """
    主函数：读取CSV，处理每一行，输出结果
    
    参数:
        in_csv: 输入CSV文件路径
        out_csv: 输出CSV文件路径
        force_llm: 是否对所有文章强制使用LLM
        inplace: 是否直接在原始CSV中添加结果列
    """
    # 读取输入CSV
    try:
        df_in = pd.read_csv(in_csv)
        logging.info(f"成功读取输入文件 {in_csv}，共{len(df_in)}条记录")
    except Exception as e:
        error_msg = f"无法读取输入文件 {in_csv}: {e}"
        logging.error(error_msg)
        sys.exit(error_msg)
    
    # 检查必要的列
    required_columns = ['pmid', 'PubMed_ID']
    has_required = any(col in df_in.columns for col in required_columns)
    if not has_required:
        logging.warning(f"CSV中缺少pmid或PubMed_ID列，处理可能不完整")
    
    # 准备数据
    rows = df_in.to_dict("records")
    processed_results = []
    
    # 处理每一行
    for row in tqdm(rows, desc="分析文献"):
        try:
            result = process_article(row, force_llm)
            processed_results.append(result)
        except Exception as exc:
            pmid = row.get("pmid", "") or row.get("PubMed_ID", "unknown")
            error_msg = f"处理PMID:{pmid}时发生错误: {exc}"
            logging.error(error_msg)
            processed_results.append({
                "Structure_Prediction_Method": "Error",
                "LLM_Response": error_msg,
                "Decision_Path": "error",
                "References_Count": 0,
                "References": "",
                "Diagnostic_Info": str(exc)
            })
    
    # 将结果添加到原始DataFrame
    if inplace:
        # 添加结果列到原始DataFrame
        for i, result in enumerate(processed_results):
            for key, value in result.items():
                df_in.loc[i, key] = value
        
        # 保存结果
        df_final = df_in
    else:
        # 创建新的DataFrame
        df_results = pd.DataFrame(processed_results)
        
        # 如果两个DataFrame长度相同，则合并
        if len(df_in) == len(df_results):
            df_final = pd.concat([df_in, df_results], axis=1)
        else:
            logging.warning("处理结果数量与输入不匹配，创建单独的结果文件")
            df_final = df_results
    
    # 保存最终结果
    df_final.to_csv(out_csv, index=False)
    
    # 统计结果
    structure_methods = sum(1 for r in processed_results if r.get("Structure_Prediction_Method") == "Yes")
    references_count = sum(r.get("References_Count", 0) for r in processed_results)
    errors = sum(1 for r in processed_results if r.get("Structure_Prediction_Method") == "Error")
    
    # 输出统计信息
    logging.info(f"处理完成，结果已保存到 {out_csv}")
    logging.info(f"共处理 {len(processed_results)} 条记录")
    logging.info(f"  RNA结构预测方法: {structure_methods} ({structure_methods/len(processed_results)*100:.1f}%)")
    logging.info(f"  总参考文献数量: {references_count}")
    logging.info(f"  平均每篇文章参考文献: {references_count/max(1, structure_methods):.1f}")
    logging.info(f"  处理错误: {errors}")
    
    logging.info("分析完成！")


if __name__ == "__main__":
    # 命令行参数解析
    ap = argparse.ArgumentParser(description="筛选RNA结构预测相关文献并获取引文")
    ap.add_argument("in_csv", help="输入CSV文件")
    ap.add_argument("out_csv", help="输出CSV文件")
    ap.add_argument("--force-llm", action="store_true", help="对所有文章强制使用LLM判断")
    ap.add_argument("--new-file", action="store_true", help="创建新文件而不是在原始文件中添加列")
    ap.add_argument("--debug", action="store_true", help="启用详细调试信息")
    ap.add_argument("--use-web-search", action="store_true", help="启用网络搜索增强分析")
    ap.add_argument("--use-local-search", action="store_true", help="使用本地化搜索方案而非SerpAPI")
    ap.add_argument("--serpapi-key", help="设置SerpAPI密钥用于网络搜索（如不使用本地搜索）")
    args = ap.parse_args()
    
    # 设置调试模式
    if args.debug:
        os.environ["DEBUG"] = "true"
        logging.getLogger().setLevel(logging.DEBUG)
        logging.debug("调试模式已启用")
    
    # 设置网络搜索选项
    if args.use_web_search:
        os.environ["USE_WEB_SEARCH"] = "true"
        logging.info("网络搜索增强已启用")
    
    # 设置本地搜索选项
    if args.use_local_search:
        os.environ["USE_LOCAL_SEARCH"] = "true"
        logging.info("本地搜索方案已启用")
        
    # 安装依赖（如果缺少）
    try:
        import bs4
    except ImportError:
        logging.warning("缺少BeautifulSoup4库，尝试自动安装...")
        import subprocess
        subprocess.check_call([sys.executable, "-m", "pip", "install", "beautifulsoup4"])
        logging.info("BeautifulSoup4安装成功")
        
    # 设置SerpAPI密钥（如果提供且未使用本地搜索）
    if args.serpapi_key and not args.use_local_search:
        os.environ["SERPAPI_API_KEY"] = args.serpapi_key
        logging.info("已设置SerpAPI密钥")
        
    # 检查网络搜索依赖性
    if os.environ.get("USE_WEB_SEARCH", "").lower() == "true":
        if not os.environ.get("USE_LOCAL_SEARCH", "").lower() == "true" and not os.environ.get("SERPAPI_API_KEY"):
            logging.warning("启用了网络搜索但未设置SERPAPI_API_KEY且未启用本地搜索，搜索功能可能不可用")
    
    main(args.in_csv, args.out_csv, args.force_llm, not args.new_file)
