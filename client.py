import requests
import json

def test_client(term: str, db: str = "pubmed", retmax: int = 10):
    """
    Test client for the FastAPI PubMed search service
    """
    try:
        # Set up the request parameters
        params = {
            "term": term,
            "db": db,
            "retmax": retmax
        }
        
        # Make the GET request
        response = requests.get(
            "http://localhost:8001/search/",
            params=params
        )
        
        # Check if the request was successful
        if response.status_code == 200:
            result = response.json()
            print("Search successful!")
            print(json.dumps(result, indent=2))
        else:
            print(f"Request failed with status code {response.status_code}")
            print(response.text)
            
    except Exception as e:
        print(f"An error occurred: {str(e)}")

if __name__ == "__main__":
    # Get user input
    term = "DNA sequencing"
    db = "pubmed"
    retmax = 3
    
    # Run the test client
    test_client(term, db, retmax)
