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
5. **Web search enhancement** - uses local searxng instance for context gathering.
6. **LangChain agent support** - uses LangChain for more advanced search capabilities.

Run:
    python rna_structure_workflow.py articles.csv output.csv

Env vars (new):
    ENTREZ_EMAIL     – for NCBI Entrez
    OPENAI_API_KEY   – OpenAI key
    OPENAI_BASE_URL  - 自定义API基础URL，如DeepSeek API
    USE_LOCAL_SEARCH - 使用本地搜索方案
    SEARXNG_URL      - searxng服务URL，默认为http://localhost:8080
    USE_LANGCHAIN_AGENT - 是否使用langchain agent进行搜索增强
    LANGCHAIN_MODEL  - langchain使用的模型，默认为deepseek-reasoner
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

# langchain相关导入
try:
    from langchain.tools import Tool
    from langchain_community.chat_models import ChatOpenAI
    from langchain.agents import initialize_agent
    LANGCHAIN_AVAILABLE = True
except ImportError:
    LANGCHAIN_AVAILABLE = False

# ---------------- configuration ---------------- #
MODEL_NAME           = os.getenv("OPENAI_MODEL", "deepseek-reasoner")  # 使用deepseek-reasoner模型
OPENAI_TEMP          = float(os.getenv("OPENAI_TEMP", 0))  # 设置温度为0以保证结果一致性
RATE_LIMIT_PER_SEC   = 3          # 对NCBI请求的礼貌延迟，每秒最多3个请求
OPENAI_BASE_URL      = os.getenv("OPENAI_BASE_URL", "https://api.deepseek.com")  # DeepSeek API基础URL
USE_WEB_SEARCH       = os.getenv("USE_WEB_SEARCH", "true").lower() == "true"  # 是否启用网络搜索，默认启用
USE_LOCAL_SEARCH     = os.getenv("USE_LOCAL_SEARCH", "true").lower() == "true"  # 是否使用本地搜索方案
LOCAL_EMBEDDINGS     = os.getenv("USE_LOCAL_EMBEDDINGS", "false").lower() == "true"  # 是否使用本地嵌入
DEBUG                = os.getenv("DEBUG", "false").lower() == "true"  # 是否启用详细调试模式
SEARXNG_URL          = os.getenv("SEARXNG_URL", "http://localhost:8080")  # searxng服务URL
USE_LANGCHAIN_AGENT  = os.getenv("USE_LANGCHAIN_AGENT", "false").lower() == "true"  # 是否使用langchain agent
LANGCHAIN_MODEL      = os.getenv("LANGCHAIN_MODEL", "deepseek-reasoner")  # langchain使用的模型
MAX_SEARCH_ITERATIONS = int(os.getenv("MAX_SEARCH_ITERATIONS", "2"))  # 最大搜索迭代次数

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

def search_with_searxng(query: str, num_results: int = 3) -> str:
    """
    使用本地部署的searxng服务进行搜索
    
    参数:
        query: 搜索查询词
        num_results: 返回结果数量
        
    返回:
        str: 搜索结果摘要
    """
    try:
        logging.info(f"使用searxng搜索: {query}")
        
        # 对查询进行URL编码
        encoded_query = urllib.parse.quote(query)
        
        # 构建searxng的JSON API URL
        url = f"{SEARXNG_URL}/search?q={encoded_query}&language=all&format=json"
        
        # 设置headers
        headers = {
            "Accept": "application/json",
            "X-Requested-With": "XMLHttpRequest"
        }
        
        # 发送请求
        response = requests.get(url, headers=headers, timeout=20)
        
        if not response.ok:
            logging.warning(f"searxng搜索返回错误: {response.status_code}")
            return ""
        
        # 解析结果
        data = response.json()
        search_results = []
        
        # 提取搜索结果
        results = data.get("results", [])
        if not results:
            logging.warning(f"searxng搜索未返回任何结果")
            return ""
            
        for i, result in enumerate(results[:num_results]):
            title = result.get("title", "")
            content = result.get("content", "")
            url = result.get("url", "")
            
            search_results.append(f"标题: {title}\n摘要: {content}\n链接: {url}\n")
        
        # 将所有结果合并为一个字符串
        combined_results = "\n".join(search_results)
        logging.info(f"searxng搜索返回{len(search_results)}条结果")
        
        return combined_results
            
    except Exception as e:
        logging.error(f"searxng搜索失败: {e}")
        return ""


def setup_langchain_agent():
    """
    设置langchain agent用于网络搜索
    
    返回:
        agent或None: 如果成功则返回初始化的agent，否则返回None
    """
    if not LANGCHAIN_AVAILABLE:
        logging.warning("langchain未安装，无法使用agent功能")
        return None
    
    try:
        # 定义搜索函数
        def search_func(query):
            return search_with_searxng(query)
        
        # 创建搜索工具
        search_tool = Tool(
            name="search_web",
            func=search_func,
            description="搜索网络获取当前信息",
        )
        
        # 创建LLM模型
        llm = ChatOpenAI(
            model_name=LANGCHAIN_MODEL, 
            temperature=0
        )
        
        # 初始化agent
        agent = initialize_agent(
            [search_tool], 
            llm, 
            agent="openai-functions", 
            verbose=DEBUG
        )
        
        logging.info(f"成功初始化langchain agent，使用模型: {LANGCHAIN_MODEL}")
        return agent
    except Exception as e:
        logging.warning(f"初始化langchain agent失败: {e}")
        return None


# 如果使用langchain agent，初始化它
langchain_agent = None
if USE_LANGCHAIN_AGENT and LANGCHAIN_AVAILABLE:
    langchain_agent = setup_langchain_agent()

def local_web_search(query: str, num_results: int = 3) -> str:
    """
    使用本地化方案进行网络搜索
    
    参数:
        query: 搜索查询词
        num_results: 返回结果数量
        
    返回:
        str: 搜索结果摘要
    """
    # 优先使用searxng搜索
    return search_with_searxng(query, num_results)


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
    使用本地搜索服务搜索相关信息
    
    参数:
        query: 搜索查询词
        num_results: 返回结果数量
        
    返回:
        str: 搜索结果摘要
    """
    # 如果使用langchain agent，优先使用它进行搜索
    if USE_LANGCHAIN_AGENT and langchain_agent:
        try:
            result = langchain_agent.run(f"搜索以下查询并返回最相关的信息: {query}")
            return result
        except Exception as e:
            logging.warning(f"langchain agent搜索失败: {e}")
    
    # 首先尝试PubMed搜索，因为我们主要处理生物医学文献
    pubmed_results = search_pubmed_articles(query, num_results)
    
    # 如果PubMed搜索成功，则返回结果
    if pubmed_results:
        return pubmed_results
        
    # 否则使用searxng搜索
    return search_with_searxng(query, num_results)


def iterative_search(query: str, context: str = "", max_iterations: int = None) -> str:
    """
    进行迭代搜索，让大语言模型根据搜索结果不断优化搜索查询
    
    参数:
        query: 初始搜索查询词
        context: 提供给模型的上下文（如文章标题、摘要等）
        max_iterations: 最大搜索迭代次数，若为None则使用全局设置
        
    返回:
        str: 合并后的搜索结果摘要和分析
    """
    if not USE_WEB_SEARCH:
        logging.warning("网络搜索功能未启用，无法进行迭代搜索")
        return ""
    
    # 使用全局配置的最大迭代次数    
    if max_iterations is None:
        max_iterations = MAX_SEARCH_ITERATIONS
        
    all_results = []
    search_queries = []
    current_query = query
    
    logging.info(f"开始迭代搜索，初始查询: {current_query}，最大迭代次数: {max_iterations}")
    
    try:
        # 首次搜索
        logging.info(f"执行初始搜索...")
        search_result = search_with_searxng(current_query)
        if not search_result:
            logging.warning("初始搜索未返回结果，尝试PubMed搜索...")
            search_result = search_pubmed_articles(current_query)
            
        if not search_result:
            logging.warning("所有搜索方法均未返回结果，终止迭代搜索")
            return ""
            
        all_results.append(f"初始搜索 [{current_query}]:\n{search_result}")
        search_queries.append(current_query)
        logging.info(f"初始搜索完成，结果长度: {len(search_result)}字符")
        
        # 如果启用了LangChain，使用agent进行后续迭代搜索
        if USE_LANGCHAIN_AGENT and langchain_agent:
            try:
                logging.info("使用LangChain Agent进行迭代搜索...")
                search_prompt = f"""
                我正在分析一篇可能是关于RNA结构预测的文章。
                
                文章上下文: {context}
                
                已有的搜索结果: {search_result}
                
                请生成1-2个更精确的搜索查询，以获取更相关的信息。每个查询应该专注于不同的方面：
                1. 这篇文章是否描述了RNA结构预测算法或方法
                2. 这类方法的技术细节或创新点
                
                返回格式:
                搜索查询1: [查询内容]
                搜索查询2: [查询内容]
                """
                
                logging.info("请求LangChain Agent生成后续搜索查询...")
                agent_response = langchain_agent.run(search_prompt)
                logging.info(f"LangChain Agent响应: {agent_response[:100]}...")
                
                # 解析agent返回的搜索查询
                new_queries = []
                for line in agent_response.split('\n'):
                    if "搜索查询" in line and ":" in line:
                        query_text = line.split(":", 1)[1].strip()
                        if query_text and query_text not in search_queries:
                            new_queries.append(query_text)
                
                logging.info(f"解析出{len(new_queries)}个后续搜索查询")
                
                # 执行额外的搜索迭代
                for i, new_query in enumerate(new_queries[:max_iterations-1]):
                    logging.info(f"执行迭代搜索 {i+1}: {new_query}")
                    iter_result = search_with_searxng(new_query)
                    if not iter_result:
                        logging.warning(f"迭代搜索 {i+1} 未返回结果，尝试PubMed搜索...")
                        iter_result = search_pubmed_articles(new_query)
                    
                    if iter_result:
                        search_queries.append(new_query)
                        all_results.append(f"迭代搜索 {i+1} [{new_query}]:\n{iter_result}")
                        logging.info(f"迭代搜索 {i+1} 完成，结果长度: {len(iter_result)}字符")
                    else:
                        logging.warning(f"迭代搜索 {i+1} 未能获取到任何结果")
                    
                # 最后让agent总结所有搜索结果
                if len(all_results) > 1:
                    logging.info("请求LangChain Agent总结所有搜索结果...")
                    summary_prompt = f"""
                    我正在分析一篇可能是关于RNA结构预测的文章。
                    
                    文章上下文: {context}
                    
                    搜索结果:
                    {''.join(all_results)}
                    
                    请分析以上搜索结果，并总结与该文章是否为RNA结构预测方法相关的重要信息。
                    """
                    
                    summary = langchain_agent.run(summary_prompt)
                    logging.info(f"LangChain Agent总结完成: {summary[:100]}...")
                    all_results.append(f"\n搜索结果分析:\n{summary}")
            
            except Exception as e:
                logging.warning(f"LangChain迭代搜索失败: {e}")
                logging.info("回退到普通LLM迭代搜索...")
                # 如果LangChain失败，回退到普通迭代搜索
        
        # 如果没有启用LangChain或LangChain失败，使用普通的LLM进行搜索迭代
        if not USE_LANGCHAIN_AGENT or not langchain_agent or len(all_results) <= 1:
            # 准备消息
            logging.info("使用标准LLM进行迭代搜索...")
            system_prompt = """你是一个专业的科学文献研究助手。
            你的任务是基于已有搜索结果生成新的搜索查询，以帮助确定一篇文章是否描述了RNA结构预测方法。
            生成的查询应该专注于获取更多相关信息，尤其是关于该文章是否提出或使用了预测RNA结构的算法。
            返回1-2个搜索查询，每行一个，格式为"搜索查询: [查询内容]"
            """
            
            user_content = f"""
            文章上下文: {context}
            
            已有的搜索结果: {search_result}
            
            请生成1-2个更精确的搜索查询，以获取更相关的信息。
            """
            
            messages = [
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": user_content},
            ]
            
            # 调用LLM生成新的搜索查询
            logging.info("请求LLM生成后续搜索查询...")
            response = client.chat.completions.create(
                model=MODEL_NAME,
                messages=messages,
                temperature=0.5,
                max_tokens=200,
                stream=False
            )
            
            llm_response = response.choices[0].message.content.strip()
            logging.info(f"LLM响应: {llm_response[:100]}...")
            
            # 解析LLM返回的搜索查询
            new_queries = []
            for line in llm_response.split('\n'):
                if "搜索查询" in line and ":" in line:
                    query_text = line.split(":", 1)[1].strip()
                    if query_text and query_text not in search_queries:
                        new_queries.append(query_text)
            
            logging.info(f"解析出{len(new_queries)}个后续搜索查询")
            
            # 执行额外的搜索迭代
            for i, new_query in enumerate(new_queries[:max_iterations-1]):
                logging.info(f"执行迭代搜索 {i+1}: {new_query}")
                iter_result = search_with_searxng(new_query)
                if not iter_result:
                    logging.warning(f"迭代搜索 {i+1} 未返回结果，尝试PubMed搜索...")
                    iter_result = search_pubmed_articles(new_query)
                
                if iter_result:
                    search_queries.append(new_query)
                    all_results.append(f"迭代搜索 {i+1} [{new_query}]:\n{iter_result}")
                    logging.info(f"迭代搜索 {i+1} 完成，结果长度: {len(iter_result)}字符")
                else:
                    logging.warning(f"迭代搜索 {i+1} 未能获取到任何结果")
                
            # 最后让LLM总结所有搜索结果
            if len(all_results) > 1:
                logging.info("请求LLM总结所有搜索结果...")
                summary_system_prompt = """你是一个专业的科学文献分析助手。
                你的任务是分析多轮搜索的结果，并总结与目标文章相关的重要信息。
                特别关注文章是否描述了RNA结构预测方法或算法。
                """
                
                summary_user_content = f"""
                文章上下文: {context}
                
                搜索结果:
                {''.join(all_results)}
                
                请分析以上搜索结果，并总结与该文章是否为RNA结构预测方法相关的重要信息。
                """
                
                summary_messages = [
                    {"role": "system", "content": summary_system_prompt},
                    {"role": "user", "content": summary_user_content},
                ]
                
                summary_response = client.chat.completions.create(
                    model=MODEL_NAME,
                    messages=summary_messages,
                    temperature=0,
                    max_tokens=400,
                    stream=False
                )
                
                summary = summary_response.choices[0].message.content.strip()
                logging.info(f"LLM总结完成: {summary[:100]}...")
                all_results.append(f"\n搜索结果分析:\n{summary}")
        
        # 返回所有搜索结果和分析
        logging.info(f"迭代搜索完成，共执行了{len(search_queries)}次搜索")
        return "\n\n".join(all_results)
        
    except Exception as e:
        logging.error(f"迭代搜索失败: {e}")
        if all_results:
            return "\n\n".join(all_results)
        return search_with_searxng(query)  # 回退到单次搜索


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
    
    # 构建上下文信息
    context = f"标题: {title}\n摘要: {abstract}"
    if doi:
        context += f"\nDOI: {doi}"
    if pmid:
        context += f"\nPMID: {pmid}"
        
    # 使用迭代搜索获取增强信息
    search_results = iterative_search(query, context)
    
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


def get_article_references(pmid: str, doi: Optional[str] = None) -> List[Dict[str, str]]:
    """
    获取文章的参考文献，使用PubMed API
    
    参数:
        pmid: PubMed ID
        doi: 文章DOI（未使用，保留参数兼容性）
        
    返回:
        List[Dict[str, str]]: 参考文献列表
    """
    # 使用PubMed获取参考文献
    if pmid:
        return get_pubmed_references(pmid)
    
    return []


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
    " Given a TITLE, ABSTRACT, and ITERATIVE WEB SEARCH RESULTS, analyze if the paper describes an RNA structure prediction method."
    " The search results include multiple iterations where I asked follow-up questions based on initial results."
    " Pay special attention to the search result analysis which summarizes key findings."
    " Use the additional context from iterative search to make a more informed decision."
    " Answer in JSON format: {\"decision\": \"yes/no\", \"reasoning\": \"your detailed explanation\", \"search_insights\": \"how search results influenced your decision\"}"
    " Answer \"yes\" **only** if the paper presents an ALGORITHM, METHOD, or computational PIPELINE whose goals include PREDICTING RNA secondary or tertiary (3D) STRUCTURE."
    " Reviews, surveys, databases, experimental protocols, purely protein‑structure methods, or papers without a predictive component => \"no\"."
    " Be specific in your reasoning, mentioning key terms or concepts that influenced your decision."
    " Your reasoning should be comprehensive and explain why this is or is not an RNA structure prediction method."
    " If iterative search results provide additional relevant information that helped your decision, mention it in the search_insights field."
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
            user_content = f"TITLE: {title}\nABSTRACT: {abstract}\nITERATIVE WEB SEARCH RESULTS: {search_results}"
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
            max_tokens=1000,  # 增加token上限以获取更详细的解释
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
    
    # 如果满足以下条件之一，则进行深入分析：
    # 1. 启发式规则判断为可能是结构预测方法
    # 2. 用户强制使用LLM分析
    if heuris_flag or force_llm:
        # 如果启用了网络搜索，先进行搜索获取额外上下文
        if USE_WEB_SEARCH:
            logging.info(f"为PMID:{pmid}启动网络搜索增强分析...")
            search_results = get_enhanced_context(title, abstract, doi, pmid)
            
            if search_results:
                diagnostic_info["web_search_performed"] = True
                diagnostic_info["search_results_length"] = len(search_results)
                logging.info(f"成功获取网络搜索结果，长度: {len(search_results)}字符")
            else:
                logging.warning(f"网络搜索未返回结果")
        
        # 调用LLM进行判断，如果有搜索结果则一并使用
        logging.info(f"开始LLM分析...")
        is_structure, llm_response = llm_is_structure_method(title, abstract, search_results)
        
        # 设置决策路径
        if search_results:
            decision_path = "llm_with_search" if heuris_flag else "force_llm_with_search"
        else:
            decision_path = "llm_after_heuristic" if heuris_flag else "force_llm"
            
        diagnostic_info["llm_called"] = True
        diagnostic_info["llm_result"] = is_structure
        
        logging.info(f"LLM分析结果: {'是' if is_structure else '否'} RNA结构预测方法")
    else:
        logging.info(f"跳过LLM调用，仅使用启发式结果: {heuris_flag}")
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
    ap.add_argument("--use-web-search", action="store_true", help="启用网络搜索增强分析（包括迭代搜索功能）")
    ap.add_argument("--searxng-url", help="设置searxng服务URL，默认为http://localhost:8080")
    ap.add_argument("--use-langchain", action="store_true", help="使用langchain agent进行搜索增强")
    ap.add_argument("--langchain-model", help="设置langchain使用的模型，默认为deepseek-reasoner") 
    ap.add_argument("--max-iterations", type=int, help="设置迭代搜索的最大迭代次数，默认为2")
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
    
    # 设置searxng URL
    if args.searxng_url:
        os.environ["SEARXNG_URL"] = args.searxng_url
        logging.info(f"已设置searxng服务URL: {args.searxng_url}")
    
    # 设置langchain选项
    if args.use_langchain:
        os.environ["USE_LANGCHAIN_AGENT"] = "true"
        logging.info("langchain agent已启用")
        
        # 安装依赖（如果缺少）
        try:
            from langchain.tools import Tool
            from langchain_community.chat_models import ChatOpenAI
            from langchain.agents import initialize_agent
        except ImportError:
            logging.warning("缺少langchain相关库，尝试自动安装...")
            import subprocess
            subprocess.check_call([sys.executable, "-m", "pip", "install", "langchain langchain-community"])
            logging.info("langchain安装成功")
            
    # 设置langchain模型
    if args.langchain_model:
        os.environ["LANGCHAIN_MODEL"] = args.langchain_model
        logging.info(f"已设置langchain模型: {args.langchain_model}")
    
    # 设置迭代搜索最大迭代次数
    if args.max_iterations:
        os.environ["MAX_SEARCH_ITERATIONS"] = str(args.max_iterations)
        logging.info(f"已设置最大迭代搜索次数: {args.max_iterations}")
        
    # 安装依赖（如果缺少）
    try:
        import bs4
    except ImportError:
        logging.warning("缺少BeautifulSoup4库，尝试自动安装...")
        import subprocess
        subprocess.check_call([sys.executable, "-m", "pip", "install", "beautifulsoup4"])
        logging.info("BeautifulSoup4安装成功")
    
    # 如果启用了网络搜索但未配置搜索引擎，提醒用户
    if os.environ.get("USE_WEB_SEARCH", "").lower() == "true":
        logging.info("将使用searxng进行网络搜索")
    
    main(args.in_csv, args.out_csv, args.force_llm, not args.new_file)
