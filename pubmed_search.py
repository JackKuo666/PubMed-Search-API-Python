from Bio import Entrez
import json
import xml.etree.ElementTree as ET


# 设置电子邮件（NCBI 要求）
Entrez.email = "email@gmail.com"

# 1. 搜索数据库（ESearch）
def search_database(term, retmax=10, search_type="keyword", db="pubmed"):
    if search_type == "keyword":
        art_type_list=['Address', 'Bibliography','Biography','Books and Documents','Clinical Conference','Clinical Study','Collected Works','Comment','Congress','Consensus Development Conference','Consensus Development Conference, NIH','Dictionary','Directory','Duplicate Publication','Editorial','Festschrift','Government Document','Guideline','Interactive Tutorial','Interview','Lecture','Legal Case','Legislation','Letter','News','Newspaper Article','Patient Education Handout','Periodical Index','Personal Narrative','Practice Guideline','Published Erratum','Technical Report','Video-Audio Media','Webcast']
        art_type = "("+" OR ".join(f'"{j}"[Filter]' for j in art_type_list)+")"
        term = "( "+ term + ")" 
        term += ' AND (fha[Filter]) NOT ' +art_type
        handle = Entrez.esearch(db=db, term=term, usehistory="y", retmax=retmax)
    elif search_type == "advanced":
        handle = Entrez.esearch(db=db, term=term, usehistory="y", retmax=retmax)
    else:
        raise ValueError("search_type must be either 'keyword' or 'advanced'")
    
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


def parse_pubmed_xml(xml_data):
    # with open("test.xml", "w", encoding='utf-8') as f:
    #     f.write(xml_data.decode('utf-8'))
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
        journal_title = article.find(".//Journal/Title").text if article.find(".//Journal/Title") is not None else ""
        issn = article.find(".//Journal/ISSN").text if article.find(".//Journal/ISSN") is not None else ""
        publication_date = article.find(".//ArticleDate[@DateType='Electronic']/Year").text if article.find(".//ArticleDate[@DateType='Electronic']/Year") is not None else ""
        language = article.find(".//Language").text if article.find(".//Language") is not None else ""
        country = article.find(".//MedlineJournalInfo/Country").text if article.find(".//MedlineJournalInfo/Country") is not None else ""
        publisher = article.find(".//MedlineJournalInfo/MedlineTA").text if article.find(".//MedlineJournalInfo/MedlineTA") is not None else ""

        # 提取作者信息
        authors = []
        for author in article.findall(".//Author"):
            last_name = author.find("LastName").text if author.find("LastName") is not None else ""
            fore_name = author.find("ForeName").text if author.find("ForeName") is not None else ""
            authors.append(f"{fore_name} {last_name}")

        # 提取关键词信息
        keywords = []
        keyword_list = article.findall(".//Keyword")
        for keyword in keyword_list:
            keywords.append(keyword.text if keyword.text is not None else "")

        # 提取文章类型信息
        article_types = []
        article_type_list = article.findall(".//PublicationType")
        for article_type in article_type_list:
            article_types.append(article_type.text if article_type.text is not None else "")

        # 将每篇文章的信息添加到列表中
        articles.append({
            "title": article_title,
            "doi": doi,
            "pmid": pmid,
            "abstract": abstract_text,
            "authors": ", ".join(authors),
            "journal_title": journal_title,
            "issn": issn,
            "publication_date": publication_date,
            "keywords": ", ".join(keywords),
            "article_types": ", ".join(article_types),
            "language": language,
            "country": country,
            "publisher": publisher
        })

    return articles

# 示例：搜索并获取记录
if __name__ == "__main__":
    # 搜索包含 "breast cancer" 的文章
    search_results = search_database("DNA sequencing", retmax=5, search_type="keyword", db="pubmed")
    print("search_results", search_results)

    query = '''
("Alzheimer Disease"[MeSH Terms] OR "Dementia"[MeSH Terms] OR "alzheimer*"[Title/Abstract] OR "dementia*"[Title/Abstract] OR "cognitive*"[Title/Abstract] OR "cognition*"[Title/Abstract]) 
AND 
(("Built Environment"[MeSH] OR "Built Environment"[Title/Abstract] OR landscape*[Title/Abstract] OR greenspace[Title/Abstract] OR "green space"[Title/Abstract] OR greenness[Title/Abstract] OR forest*[Title/Abstract] OR park*[Title/Abstract] OR greenway[Title/Abstract] OR outdoor*[Title/Abstract] OR walkabilit*[Title/Abstract] OR street[Title/Abstract] OR roadway[Title/Abstract] OR residen*[Title/Abstract] OR neighborhood*[Title/Abstract] OR neighbourhood*[Title/Abstract])) 
AND 
(("street view" OR streetview OR satellite* OR GIS OR "geographic information system*" OR "remote sensing" OR ArcGIS OR postal-code* OR zipcode* OR "postal code*" OR "zip code*" OR zip-code*))
'''
    search_results = search_database(query, retmax=5, search_type="advanced", db="pubmed")
    print("search_results", search_results)
