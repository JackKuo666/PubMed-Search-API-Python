import requests
import json

def test_client(query: str, num_to_show: int = 10, search_type: str = "keyword", db: str = "pubmed"):
    """
    Test client for the FastAPI PubMed search service
    """
    try:
        # Set up the request parameters
        params = {
            "query": query,
            "num_to_show": num_to_show,
            "search_type": search_type,
            "db": db
        }
        
        # Make the GET request
        response = requests.get(
            # "http://localhost:8001/search/",
            "http://10.15.56.148:8001/search/",
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
    # Get keyword 
    query = "DNA sequencing"
    num_to_show = 3
    
    # Run the test client
    test_client(query, num_to_show, search_type="keyword", db="pubmed")

#     # advanced query
#     query2 = '''
# ("Alzheimer Disease"[MeSH Terms] OR "Dementia"[MeSH Terms] OR "alzheimer*"[Title/Abstract] OR "dementia*"[Title/Abstract] OR "cognitive*"[Title/Abstract] OR "cognition*"[Title/Abstract]) 
# AND 
# (("Built Environment"[MeSH] OR "Built Environment"[Title/Abstract] OR landscape*[Title/Abstract] OR greenspace[Title/Abstract] OR "green space"[Title/Abstract] OR greenness[Title/Abstract] OR forest*[Title/Abstract] OR park*[Title/Abstract] OR greenway[Title/Abstract] OR outdoor*[Title/Abstract] OR walkabilit*[Title/Abstract] OR street[Title/Abstract] OR roadway[Title/Abstract] OR residen*[Title/Abstract] OR neighborhood*[Title/Abstract] OR neighbourhood*[Title/Abstract])) 
# AND 
# (("street view" OR streetview OR satellite* OR GIS OR "geographic information system*" OR "remote sensing" OR ArcGIS OR postal-code* OR zipcode* OR "postal code*" OR "zip code*" OR zip-code*))
# '''
#     test_client(query2, num_to_show, search_type="advanced", db="pubmed")