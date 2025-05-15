#!/usr/bin/env python3
"""
NCBI GenBank Data Retriever
Rozszerzony skrypt do łączenia się z NCBI i pobierania rekordów sekwencji genetycznych dla danego identyfikatora taksonomicznego.
Dodano filtry długości sekwencji, generowanie raportu CSV oraz wizualizację danych.
"""

from Bio import Entrez
import time
import os
import pandas as pd
import matplotlib.pyplot as plt


class NCBIRetriever:
    def __init__(self, email, api_key):
        """Initialize with NCBI credentials."""
        self.email = email
        self.api_key = api_key

        # Ustawienia Entrez
        Entrez.email = email
        Entrez.api_key = api_key
        Entrez.tool = 'BioScriptEx10'

    def search_taxid(self, taxid, min_length=None, max_length=None):
        """Search for all records associated with a taxonomic ID and filter by sequence length."""
        print(f"Searching for records with taxID: {taxid}")
        try:
            # Najpierw pobierz informacje taksonomiczne
            handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
            records = Entrez.read(handle)
            organism_name = records[0]["ScientificName"]
            print(f"Organism: {organism_name} (TaxID: {taxid})")

            # Szukaj rekordów
            search_term = f"txid{taxid}[Organism]"
            handle = Entrez.esearch(db="nucleotide", term=search_term, usehistory="y")
            search_results = Entrez.read(handle)
            count = int(search_results["Count"])

            if count == 0:
                print(f"No records found for {organism_name}")
                return None

            print(f"Found {count} records")

            # Zapisz wyniki wyszukiwania do późniejszego wykorzystania
            self.webenv = search_results["WebEnv"]
            self.query_key = search_results["QueryKey"]
            self.count = count

            # Filtruj rekordy według długości, jeśli podano zakres
            if min_length or max_length:
                self.filtered_records = self.filter_records_by_length(min_length, max_length)
            else:
                self.filtered_records = None

            return count

        except Exception as e:
            print(f"Error searching TaxID {taxid}: {e}")
            return None

    def filter_records_by_length(self, min_length, max_length):
        """Filter the records by sequence length."""
        filtered_records = []
        # Załaduj rekordy, iterując po nich
        start = 0
        batch_size = 500
        while True:
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
            if not records_text:
                break

            # Przetwórz każdy rekord
            records = records_text.split("\n//\n")
            for record in records:
                if not record:
                    continue

                # Sprawdź długość sekwencji
                lines = record.splitlines()
                length_line = next((line for line in lines if line.startswith("     length")), None)
                if length_line:
                    length = int(length_line.split()[1])

                    if (min_length is None or length >= min_length) and (max_length is None or length <= max_length):
                        filtered_records.append(record)

            start += batch_size

        return filtered_records

    def fetch_records(self, start=0, max_records=10):
        """Fetch a batch of records using the stored search results."""
        if not hasattr(self, 'webenv') or not hasattr(self, 'query_key'):
            print("No search results to fetch. Run search_taxid() first.")
            return []

        try:
            # Limit, aby zapobiec przeciążeniu serwera
            batch_size = min(max_records, 500)

            handle = Entrez.efetch(
                db="nucleotide",
                rettype="gb",
                retmode="text",
                retstart=start,
                retmax=batch_size,
                webenv=self.webenv,
                query_key=self.query_key
            )

            # Surowy rekord GenBank
            records_text = handle.read()

            return records_text

        except Exception as e:
            print(f"Error fetching records: {e}")
            return ""

    def generate_csv_report(self, filename="genbank_report.csv"):
        """Generate a CSV report with accession numbers, sequence lengths, and descriptions."""
        if not hasattr(self, 'filtered_records'):
            print("No records to generate a report. Fetch some records first.")
            return

        data = []
        for record in self.filtered_records:
            lines = record.splitlines()
            accession = next((line.split()[1] for line in lines if line.startswith("ACCESSION")), None)
            description = next((line.split("  ")[1] for line in lines if line.startswith("  DEFINITION")), None)
            length = int(next((line.split()[1] for line in lines if line.startswith("     length")), None))

            if accession and description:
                data.append({"accession": accession, "length": length, "description": description})

        # Tworzymy DataFrame z listy danych
        df = pd.DataFrame(data)

        # Zapisz do pliku CSV
        df.to_csv(filename, index=False)
        print(f"CSV report saved to {filename}")

    def generate_length_plot(self, filename="sequence_lengths.png"):
        """Generate a plot of sequence lengths."""
        if not hasattr(self, 'filtered_records'):
            print("No records to generate a plot. Fetch some records first.")
            return

        lengths = []
        accession_numbers = []
        for record in self.filtered_records:
            lines = record.splitlines()
            accession = next((line.split()[1] for line in lines if line.startswith("ACCESSION")), None)
            length = int(next((line.split()[1] for line in lines if line.startswith("     length")), None))

            if accession:
                lengths.append(length)
                accession_numbers.append(accession)

        # Posortuj po długości
        data = list(zip(accession_numbers, lengths))
        data.sort(key=lambda x: x[1], reverse=True)

        accession_numbers_sorted, lengths_sorted = zip(*data)

        # Stwórz wykres
        plt.figure(figsize=(10, 6))
        plt.plot(accession_numbers_sorted, lengths_sorted, marker='o', linestyle='-', color='b')

        plt.xlabel("Accession Number")
        plt.ylabel("Sequence Length")
        plt.title("GenBank Records Sorted by Sequence Length")
        plt.xticks(rotation=90)
        plt.tight_layout()

        # Zapisz wykres jako plik PNG
        plt.savefig(filename)
        plt.close()
        print(f"Plot saved to {filename}")


def main():
    # Uzyskaj dane uwierzytelniające
    email = input("Enter your email address for NCBI: ")
    api_key = input("Enter your NCBI API key: ")

    # Utwórz obiekt retriever
    retriever = NCBIRetriever(email, api_key)

    # Uzyskaj taxid od użytkownika
    taxid = input("Enter taxonomic ID (taxid) of the organism: ")

    # Uzyskaj minimalną i maksymalną długość sekwencji
    min_length = input("Enter minimum sequence length (or press Enter to skip): ")
    min_length = int(min_length) if min_length else None

    max_length = input("Enter maximum sequence length (or press Enter to skip): ")
    max_length = int(max_length) if max_length else None

    # Szukaj rekordów
    count = retriever.search_taxid(taxid, min_length, max_length)

    if not count:
        print("No records found. Exiting.")
        return

    # Pobierz kilka pierwszych rekordów jako próbkę
    print("\nFetching sample records...")
    sample_records = retriever.fetch_records(start=0, max_records=5)

    # Zapisz do pliku
    output_file = f"taxid_{taxid}_sample.gb"
    with open(output_file, "w") as f:
        f.write(sample_records)

    print(f"Saved sample records to {output_file}")

    # Generowanie raportu CSV
    retriever.generate_csv_report()

    # Generowanie wykresu długości
    retriever.generate_length_plot()


if __name__ == "__main__":
    main()
