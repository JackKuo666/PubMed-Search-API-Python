from fastapi import FastAPI, Query
from pubmed_search import search_database
from typing import List

app = FastAPI()

@app.get("/")
async def root():
    return {"message": "Welcome to the PubMed Search API"}

@app.get("/search/")
async def search(query: str = Query(..., description="Search term"), 
                 num_to_show: int = Query(20, description="Maximum number of results to return"),
                 db: str = Query("pubmed", description="Database to search")):
    results = search_database(term=query, db=db, retmax=num_to_show)
    return {"results": results}

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8001)
