import argparse


from Bio import Entrez
Entrez.email = "<ghvprr@gmail.com>"
def search_pubmed(term, max_results=10):
    handle = Entrez.esearch(db="pubmed", term=term, retmax=max_results)
    record = Entrez.read(handle)
    return record["IdList"]
def main():
    parser = argparse.ArgumentParser(description="one or more pubmed ids+")
    parser.add_argument("terms", nargs="+", help="One or more search terms")
    args = parser.parse_args()

    for term in args.terms:
        ids = search_pubmed(term)
        print(f"\nResults for {term}:")
        print(ids)
if __name__ == "__main__":
    main()


    


