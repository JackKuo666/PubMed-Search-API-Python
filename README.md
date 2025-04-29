# PubMed Search API Python Client 🔍

A Python client for interacting with the PubMed Search API to retrieve and manage biomedical literature. 📚

## Project Structure 📁

```
PubMed_search_api_python/
├── client.py          # Core client implementation
├── main.py            # Main application entry point
├── pubmed_search.py   # Search functionality and data handling
├── README.md          # This documentation file
└── requirements.txt   # Project dependencies
```

## Dependencies 📦
Make sure you have installed the following dependencies:
- `biopython` 🧬
- `fastapi` ⚡
- `uvicorn` 🦄
- `requests` 🌐

You can install these dependencies using:
```bash
pip install biopython fastapi uvicorn requests
```

## Installation 💻

1. Clone the repository
```bash
https://github.com/JackKuo666/PubMed-Search-API-Python.git
```
2. Install dependencies:
```bash
pip install -r requirements.txt
```

## Features ✨
### 1. PubMed Search 🔎
The pubmed_search.py file contains the search_database function for searching articles and retrieving abstracts from the PubMed database.

### 2. FastAPI Client 🌟
The client.py file includes the test_client function for testing the search functionality of the FastAPI service.

### 3. FastAPI Main Application 🚀
The main.py file defines the basic routes for the FastAPI application.

## Usage 📖
### 1. Start FastAPI Service 🎯
Run the following command in the terminal to start the FastAPI service:

```bash
uvicorn main:app --reload
```
This will start a local server at http://127.0.0.1:8001 by default.

### 2. Test Search Function Using Client 🧪
In client.py, the test_client function is used to test the search functionality. You can modify the term, db, and retmax parameters as needed.

Run the client.py file:
```bash
python client.py
```

## Contributing 🤝
Contributions to this project are welcome! If you find any issues or have suggestions for improvements, please submit an Issue or Pull Request.

## License 📄
This project is licensed under the MIT License. See the LICENSE file for details.