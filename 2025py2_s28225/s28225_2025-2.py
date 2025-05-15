#!/usr/bin/env python3
"""
NCBI GenBank Data Retriever
Pobiera rekordy z GenBanku dla podanego taxid, umożliwia filtrowanie po długości sekwencji (opcjonalnie).
"""

from Bio import Entrez
import os

class NCBIRetriever:
    def __init__(self, email, api_key):
        """Initialize with NCBI credentials."""
        self.email = email
        self.api_key = api_key

        Entrez.email = email
        Entrez.api_key = api_key
        Entrez.tool = 'BioScriptEx10'

    def search_taxid(self, taxid, min_len=None, max_len=None):
        """Search for all records associated with a taxonomic ID, with optional sequence length filtering."""
        print(f"Searching for records with taxID: {taxid}")

        try:
            # Fetch taxonomy info
            handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
            records = Entrez.read(handle)
            organism_name = records[0]["ScientificName"]
            print(f"Organism: {organism_name} (TaxID: {taxid})")

            # Build search term with optional sequence length filter
            search_term = f"txid{taxid}[Organism]"

            if min_len is not None and max_len is not None:
                search_term += f" AND {min_len}:{max_len}[SLEN]"

            # Search the nucleotide database with the created search term
            handle = Entrez.esearch(db="nucleotide", term=search_term, usehistory="y")
            search_results = Entrez.read(handle)
            count = int(search_results["Count"])

            if count == 0:
                print(f"No records found for {organism_name} in the specified length range.")
                return None

            print(f"Found {count} records")

            # Save search history info
            self.webenv = search_results["WebEnv"]
            self.query_key = search_results["QueryKey"]
            self.count = count

            return count

        except Exception as e:
            print(f"Error searching TaxID {taxid}: {e}")
            return None

    def fetch_records(self, start=0, max_records=10):
        """Fetch a batch of records using the stored search results."""
        if not hasattr(self, 'webenv') or not hasattr(self, 'query_key'):
            print("No search results to fetch. Run search_taxid() first.")
            return []

        try:
            batch_size = min(max_records, 50)
            handle = Entrez.efetch(
                db="nucleotide",
                rettype="gb",
                retmode="text",
                retstart=start,
                retmax=batch_size,
                webenv=self.webenv,
                query_key=self.query_key
            )
            records_text = handle.read()
            return records_text

        except Exception as e:
            print(f"Error fetching records: {e}")
            return ""

def main():
    print("=== NCBI GenBank Data Retriever ===")
    # User input
    email = input("Enter your email address for NCBI: ")
    api_key = input("Enter your NCBI API key: ")
    taxid = input("Enter taxonomic ID (taxid) of the organism: ")

    # Ask for length filters (optional)
    min_len = input("Enter minimum sequence length (e.g., 200) [Press Enter to skip]: ")
    max_len = input("Enter maximum sequence length (e.g., 3000) [Press Enter to skip]: ")

    # Convert to integers if provided
    min_len = int(min_len) if min_len else None
    max_len = int(max_len) if max_len else None

    # Create retriever object
    retriever = NCBIRetriever(email, api_key)

    # Search records with length filter if provided
    count = retriever.search_taxid(taxid, min_len=min_len, max_len=max_len)

    if not count:
        print("No records found. Exiting.")
        return

    # Fetch sample records
    print("\nFetching sample records...")
    sample_records = retriever.fetch_records(start=0, max_records=5)

    # Save to file
    output_file = f"taxid_{taxid}_sample.gb"
    with open(output_file, "w") as f:
        f.write(sample_records)

    print(f"Saved sample records to {output_file}")
    print("\nNote: This is just a basic retriever. You can now extend it with visualization.")

if __name__ == "__main__":
    main()
