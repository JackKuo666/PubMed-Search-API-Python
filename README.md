<<<<<<< HEAD
# PubMed Search API Python Client

A Python client for interacting with the PubMed Search API to retrieve and manage biomedical literature.

## Project Structure

```
PubMed_search_api_python/
├── client.py          # Core client implementation
├── main.py            # Main application entry point
├── pubmed_search.py   # Search functionality and data handling
├── README.md          # This documentation file
└── requirements.txt   # Project dependencies
```


## 项目依赖
确保你已经安装了以下依赖项：
- `biopython`
- `fastapi`
- `uvicorn`
- `requests`

你可以通过以下命令安装这些依赖项：
```bash
pip install biopython fastapi uvicorn requests
=======
# PubMed Search API Python

This repository provides a simple Python interface to search and retrieve articles from PubMed, a database of biomedical literature. It utilizes the PubMed API (Entrez Programming Utilities, E-utilities) to fetch data, which can then be processed and analyzed.

## Features

- Search PubMed using keywords or query terms.
- Fetch metadata for articles, including title, authors, publication date, and more.
- Retrieve full article abstracts and other details.
- Easy-to-use and extensible Python API.
  
## Requirements

- Python 3.10+
- `requests` library

You can install the required dependencies using `pip`:

```bash
pip install requests
>>>>>>> f23d49a646f61cc6b6b8a79327b288f1c6ca5ca8
```

## Installation

<<<<<<< HEAD
1. Clone the repository
2. Install dependencies:
```bash
pip install -r requirements.txt
```


## 功能
### 1. PubMed 搜索
pubmed_search.py 文件中包含 search_database 函数，用于在 PubMed 数据库中搜索文章并获取摘要。

### 2. FastAPI 客户端
client.py 文件中包含 test_client 函数，用于测试 FastAPI 服务的搜索功能。

### 3. FastAPI 主应用
main.py 文件中定义了 FastAPI 应用的基本路由。

## 使用方法
### 1. 启动 FastAPI 服务
在终端中运行以下命令以启动 FastAPI 服务：

```bash
uvicorn main:app --reload
```
这将启动一个本地服务器，默认地址为 http://127.0.0.1:8000。

### 2. 使用客户端测试搜索功能
在 client.py 文件中，test_client 函数用于测试搜索功能。你可以根据需要修改 term、db 和 retmax 参数。

运行 client.py 文件：

```bash
python client.py
```
## 贡献
欢迎对该项目进行贡献！如果你发现任何问题或有改进建议，请提交一个 Issue 或 Pull Request。

## 许可证
本项目采用 MIT 许可证。详细信息请参阅 LICENSE 文件。
=======
To use this API in your project, you can either clone the repository or install it directly:

### Clone the repository:

```bash
git clone https://github.com/JackKuo666/PubMed_search_api_python.git
cd PubMed_search_api_python
```

### Install via pip (if available in PyPI):

```bash
pip install PubMed-search-api-python
```

## Usage

1. **Import the Library**  
   Start by importing the `PubMedSearch` class from the module:

   ```python
   from pubmed_search_api import PubMedSearch
   ```

2. **Search PubMed**  
   To search for articles, simply create an instance of `PubMedSearch` and use the `search` method. You can pass in keywords or a complex query string.

   ```python
   pubmed = PubMedSearch()
   results = pubmed.search("cancer immunotherapy")
   ```

3. **Fetch Detailed Information**  
   After performing a search, you can retrieve detailed information (such as abstracts and authors) for each article in the result.

   ```python
   for article in results:
       print(f"Title: {article['title']}")
       print(f"Authors: {', '.join(article['authors'])}")
       print(f"Abstract: {article['abstract']}")
   ```

## Example

Here's an example of a full search query to fetch articles on "neurodegenerative diseases":

```python
from pubmed_search_api import PubMedSearch

pubmed = PubMedSearch()
results = pubmed.search("neurodegenerative diseases")
for article in results:
    print(f"Title: {article['title']}")
    print(f"Authors: {', '.join(article['authors'])}")
    print(f"Journal: {article['journal']}")
    print(f"Published: {article['pub_date']}")
    print(f"Abstract: {article['abstract']}")
    print("="*50)
```

## Functions

### `search(query: str, max_results: int = 20) -> list`

Search PubMed with the provided query string and return a list of articles (up to `max_results`).

#### Parameters:
- `query` (str): The search query string (e.g., "cancer immunotherapy").
- `max_results` (int): The maximum number of results to return (default is 20).

#### Returns:
- A list of dictionaries, each containing article details such as `title`, `authors`, `journal`, `pub_date`, and `abstract`.

### `fetch_article_details(pubmed_id: str) -> dict`

Fetch detailed information for a specific article by its PubMed ID.

#### Parameters:
- `pubmed_id` (str): The PubMed ID of the article.

#### Returns:
- A dictionary with detailed article information, including the title, abstract, authors, and publication date.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- PubMed E-utilities for providing the API.
- `requests` library for handling HTTP requests.
>>>>>>> f23d49a646f61cc6b6b8a79327b288f1c6ca5ca8
