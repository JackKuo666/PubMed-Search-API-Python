from Bio import Entrez

# 设置电子邮件（NCBI 要求）
Entrez.email = "your_email@example.com"

# 1. 搜索数据库（ESearch）
def search_database(term, db="pubmed"):
    handle = Entrez.esearch(db=db, term=term, usehistory="y")
    results = Entrez.read(handle)
    handle.close()
    return results

# 2. 获取记录详情（EFetch）
def fetch_details(id_list, db="pubmed", rettype="abstract"):
    ids = ",".join(id_list)
    handle = Entrez.efetch(db=db, id=ids, rettype=rettype, retmode="text")
    records = handle.read()
    handle.close()
    return records

# 示例：搜索并获取记录
if __name__ == "__main__":
    # 搜索包含 "breast cancer" 的文章
    search_results = search_database("breast cancer", db="pubmed")
    id_list = search_results["IdList"]
    # 获取前 5 篇文章的摘要
    records = fetch_details(id_list[:5], db="pubmed", rettype="abstract")
    print(records)
