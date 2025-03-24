from fastapi import FastAPI, Query
from pubmed_search import search_database
from typing import List

app = FastAPI()

@app.get("/")
async def root():
    return {"message": "Welcome to the PubMed Search API"}

@app.get("/search/")
async def search(term: str = Query(..., description="Search term"), 
                 db: str = Query("pubmed", description="Database to search"),
                 retmax: int = Query(20, description="Maximum number of results to return")):
    results = search_database(term, db, retmax)
    return {"results": results}

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8001)
