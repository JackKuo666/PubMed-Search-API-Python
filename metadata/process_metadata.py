import os
import gzip
import xml.etree.ElementTree as ET
import pandas as pd
import tqdm
import glob

import xml.etree.ElementTree as ET

def extract_meta_info(xml_content):
    root = ET.fromstring(xml_content)
    meta_info_list = []  # List to hold metadata of each article

    # Loop over each article in the XML (assuming it's in a root <PubmedArticle> list)
    articles = root.findall(".//PubmedArticle")  # Or adjust the XPath based on your XML structure
    
    for article in articles:
        meta_info = {}
        
        # Extract PMID
        pmid = article.find(".//PMID")
        meta_info['PMID'] = pmid.text if pmid is not None else None
        
        # Extract DateCompleted
        date_completed = article.find(".//DateCompleted")
        if date_completed is not None:
            year = date_completed.find(".//Year")
            month = date_completed.find(".//Month")
            day = date_completed.find(".//Day")
            meta_info['DateCompleted'] = f"{year.text}-{month.text}-{day.text}" if year is not None and month is not None and day is not None else None
        
        # Extract DateRevised
        date_revised = article.find(".//DateRevised")
        if date_revised is not None:
            year = date_revised.find(".//Year")
            month = date_revised.find(".//Month")
            day = date_revised.find(".//Day")
            meta_info['DateRevised'] = f"{year.text}-{month.text}-{day.text}" if year is not None and month is not None and day is not None else None
        
        # Extract ISSN
        issn = article.find(".//ISSN")
        meta_info['ISSN'] = issn.text if issn is not None else None
        
        # Extract Journal Title
        journal_title = article.find(".//Journal/Title")
        meta_info['JournalTitle'] = journal_title.text if journal_title is not None else None
        
        # Extract Article Title
        article_title = article.find(".//ArticleTitle")
        meta_info['ArticleTitle'] = article_title.text if article_title is not None else None
        
        # Extract Authors
        authors = article.findall(".//AuthorList/Author")
        author_names = []
        for author in authors:
            last_name = author.find(".//LastName")
            fore_name = author.find(".//ForeName")
            if last_name is not None and fore_name is not None:
                author_names.append(f"{last_name.text} {fore_name.text}")
        meta_info['Authors'] = ', '.join(author_names) if author_names else None
        
        # Extract Language
        language = article.find(".//Language")
        meta_info['Language'] = language.text if language is not None else None
        
        # Extract Grants
        grants = article.findall(".//GrantList/Grant")
        grant_info = []
        for grant in grants:
            grant_id = grant.find(".//GrantID")
            agency = grant.find(".//Agency")
            country = grant.find(".//Country")
            if grant_id is not None and agency is not None and country is not None:
                grant_info.append(f"{grant_id.text} ({agency.text}, {country.text})")
        meta_info['Grants'] = '; '.join(grant_info) if grant_info else None
        
        # Extract Publication Types
        publication_types = article.findall(".//PublicationTypeList/PublicationType")
        pub_types = []
        for pub_type in publication_types:
            pub_types.append(pub_type.text)
        meta_info['PublicationTypes'] = ', '.join(pub_types) if pub_types else None
        
        # Extract Chemicals
        chemicals = article.findall(".//ChemicalList/Chemical")
        chemical_info = []
        for chemical in chemicals:
            substance_name = chemical.find(".//NameOfSubstance")
            if substance_name is not None:
                chemical_info.append(substance_name.text)
        meta_info['Chemicals'] = ', '.join(chemical_info) if chemical_info else None
        
        # Extract CitationSubset
        citation_subset = article.find(".//CitationSubset")
        meta_info['CitationSubset'] = citation_subset.text if citation_subset is not None else None
        
        # Extract Article IDs (DOI, etc.)
        article_ids = article.findall(".//ArticleIdList/ArticleId")
        article_id_info = []
        for article_id in article_ids:
            article_id_info.append(article_id.text)
        meta_info['ArticleIds'] = ', '.join(filter(None, article_id_info)) if article_id_info else None
        
        # Extract Abstract
        abstract = article.find(".//Abstract/AbstractText")
        meta_info['Abstract'] = abstract.text if abstract is not None else None
        
        # Extract Mesh Terms
        mesh_terms = article.findall(".//MeshHeadingList/MeshHeading")
        mesh_terms_info = []
        for mesh_term in mesh_terms:
            descriptor_name = mesh_term.find(".//DescriptorName")
            if descriptor_name is not None:
                mesh_terms_info.append(descriptor_name.text)
        meta_info['MeshTerms'] = ', '.join(filter(None, mesh_terms_info)) if mesh_terms_info else None
        
        # Extract Keywords
        keywords = article.findall(".//KeywordList/Keyword")
        keyword_info = []
        for keyword in keywords:
            keyword_info.append(keyword.text)
        meta_info['Keywords'] = ', '.join(filter(None, keyword_info)) if keyword_info else None
        
        # Append the metadata for this article to the list
        meta_info_list.append(meta_info)
    
    return meta_info_list

def extract(input_dir, output_csv):
    # Create a temporary directory to store individual CSVs
    temp_dir = os.path.join(os.path.dirname(output_csv), 'temp')
    os.makedirs(temp_dir, exist_ok=True)
    
    # Iterate over all .gz files in the directory
    for filename in tqdm.tqdm(os.listdir(input_dir)):
        if filename.endswith('.xml.gz'):
            file_path = os.path.join(input_dir, filename)
            
            # Decompress and read the XML content
            with gzip.open(file_path, 'rb') as f:
                xml_content = f.read()
            
            # Extract meta information
            meta_info_list = extract_meta_info(xml_content)

            # Save meta information to a temporary CSV file
            temp_csv_path = os.path.join(temp_dir, f"{os.path.splitext(filename)[0]}.csv")

            # Create a DataFrame from the list of dictionaries (each dict represents an article's metadata)
            df = pd.DataFrame(meta_info_list)

            # Save the DataFrame to a CSV file
            df.to_csv(temp_csv_path, index=False)
    
    # Combine all temporary CSVs into a single large CSV file
    all_csv_files = glob.glob(os.path.join(temp_dir, '*.csv'))
    combined_df = pd.concat((pd.read_csv(f) for f in all_csv_files), ignore_index=True)
    combined_df.to_csv(output_csv, index=False)
    
    # Optionally, delete the temporary files
    # for f in all_csv_files:
    #     os.remove(f)
    # os.rmdir(temp_dir)
    
    print(f"Meta information extracted and saved to {output_csv}")

if __name__ == "__main__":
    # Define the input directory
    input_dir = '/home/sunshanxin/pubmed/pubmed_data'  # Replace with actual path
    output_csv = '/home/zjlab/code/PubMed_search_api_python/metadata/meta_info_2025_0327.csv'  # Output CSV file path
    extract(input_dir=input_dir, output_csv=output_csv)
