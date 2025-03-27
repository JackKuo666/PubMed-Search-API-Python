from Bio import Entrez
import json

# 设置电子邮件（NCBI 要求）
Entrez.email = "email@gmail.com"

# 1. 搜索数据库（ESearch）
def search_database(term, db="pubmed", retmax=10):
    handle = Entrez.esearch(db=db, term=term, usehistory="y", retmax=retmax)
    results = Entrez.read(handle)
    handle.close()
    id_list = results["IdList"]
    records = fetch_details(id_list, db="pubmed", rettype="abstract")
    return records

# 2. 获取记录详情（EFetch）
def fetch_details(id_list, db="pubmed", rettype="abstract"):
    ids = ",".join(id_list)
    handle = Entrez.efetch(db=db, id=ids, rettype=rettype, retmode="xml")
    records = handle.read()
    records = parse_pubmed_xml(records)
    handle.close()
    return records


import xml.etree.ElementTree as ET

def parse_pubmed_xml(xml_data):
    # 解析 XML 数据
    tree = ET.ElementTree(ET.fromstring(xml_data))
    root = tree.getroot()

    articles = []

    # 遍历每个 PubmedArticle 元素
    for article in root.findall(".//PubmedArticle"):
        # 提取文章信息
        article_title = article.find(".//ArticleTitle").text if article.find(".//ArticleTitle") is not None else ""
        doi = article.find(".//ELocationID[@EIdType='doi']").text if article.find(".//ELocationID[@EIdType='doi']") is not None else ""
        pmid = article.find(".//ArticleId[@IdType='pubmed']").text if article.find(".//ArticleId[@IdType='pubmed']") is not None else ""
        abstract_texts = article.findall(".//AbstractText")
        abstract_text = " ".join([abstract.text if abstract.text is not None else "" for abstract in abstract_texts]) if abstract_texts else ""

        # 提取作者信息
        authors = []
        for author in article.findall(".//Author"):
            last_name = author.find("LastName").text if author.find("LastName") is not None else ""
            fore_name = author.find("ForeName").text if author.find("ForeName") is not None else ""
            authors.append(f"{fore_name} {last_name}")

        # 将每篇文章的信息添加到列表中
        articles.append({
            "title": article_title,
            "doi": doi,
            "pmid": pmid,
            "abstract": abstract_text,
            "authors": ", ".join(authors)
        })

    return articles

# 示例：搜索并获取记录
if __name__ == "__main__":
    # 搜索包含 "breast cancer" 的文章
    search_results = search_database("DNA sequencing", db="pubmed", retmax=5)
    print("search_results", search_results)
