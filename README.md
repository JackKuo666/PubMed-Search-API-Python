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
```

## Installation

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
