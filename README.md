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
```

## Installation

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