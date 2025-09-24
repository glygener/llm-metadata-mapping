import argparse
import csv
import os
import sys
import time
from typing import List, Tuple, Dict, Optional


from Bio import Entrez


Entrez.email = "<ghvprr@gmail.com>"

def esearch_tax_ids(term: str, max_return: int=5) -> List[str]:
    with Entrez.esearch(db= "taxonomy", term=term, retmax=max_return) as handle:
        records = Entrez.read(handle)
    return list(records.get["IdList"])

def fetch_tax_ids(ids: List[str]) -> List[Dict]:
    if not ids:
        return []
    with Entrez.efetch(db="taxonomy", id=",".join(ids)) as handle:
        records = Entrez.read(handle)
        return records.get["IdList"]

def best_taxon_match(name: str,  return_max: int=5, pause: float= 0.34) -> Tuple[str, str, str]:
    ids = esearch_tax_ids(f'"{name}"[Scientific Name]', max_return=return_max)
    time.sleep(pause)
    if not ids:
        ids = esearch_tax_ids(name)
        time.sleep(pause)
        if not ids:
            raise ValueError("no_match")

    records= fetch_tax_ids(ids)
    time.sleep(pause)

    name_cf = name.casefold()
    exact = None
    for r in records:
        sci = r.get("ScientificName", "") or ""
        if sci.casefold() == name_cf:
            exact = r
            break

    chosen = exact or (records[0] if records else None)
    if not chosen:
        raise ValueError("no_match")

    return chosen.get("TaxId", ""), chosen.get("ScientificName", ""), chosen.get("Rank", "")

def read_names_from_csv(csv_path: str, name_col: Optional[str])-> Tuple[List[str], List[Dict[str, str]], str]:
    with open(csv_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames or []
        if not fieldnames:
            raise ValueError("CSV has no header row.")
        chosen = name_col or fieldnames[0]
        if chosen not in fieldnames:
            raise ValueError(f"Column '{chosen}' not found. Available: {fieldnames}")
        rows = list(reader)
    names = [r.get(chosen, "").strip() for r in rows]
    return names, rows, chosen


def write_out_csv(csv_path: str, rows: List[Dict[str, str]], extra_cols: List[str]) -> str:
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        fieldnames = list(rows[0].keys()) if rows else []
        for c in extra_cols:
            if c not in fieldnames:
                fieldnames.append(c)
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return csv_path

def build_arg_parser() -> argparse:
    p = argparse.ArgumentParser(
        description="Match species names (from CSV or args) to NCBI Taxonomy IDs."
    )
    p.add_argument("--csv", dest="csv_in", help="Input CSV with species names.")
    p.add_argument("--name-col", help="Column name containing species (default: first column).")
    p.add_argument("--out", dest="csv_out", help="Output CSV path (default: <input>.tax_ids.csv).")
    p.add_argument("--return_max", type=int, default=5, help="Max taxonomy hits to consider (default 5).")
    p.add_argument("--email", help="Your email (or set NCBI_EMAIL).")
    p.add_argument("--api-key", help="NCBI API key (or set NCBI_API_KEY).")
    p.add_argument("terms", nargs="*", help="Optional: species names as positional args (if no CSV).")
def main() -> None:
    args = build_arg_parser().parse_args()
    cache: Dict[str, Tuple[str, str, str]] = {}
    results: List[Dict[str, str]] = []

    if args.csv_in:
        names, rows, chosen_col = read_names_from_csv(args.csv_in, args.name_col)
        for name, row in zip(names, rows):
            out_row = dict(row)
            if not name:
                out_row.update({"tax_id": "", "matched_name": "", "rank": "", "status": "empty_name"})
                results.append(out_row)
                continue

            key = name.casefold()
            try:
                if key not in cache:
                    cache[key] = best_taxon_match(name, return_max=args.retmax)
                tax_id, matched, rank = cache[key]
                out_row.update({"tax_id": tax_id, "matched_name": matched, "rank": rank, "status": "ok"})
            except Exception as e:
                out_row.update({"tax_id": "", "matched_name": "", "rank": "", "status": str(e)})
            results.append(out_row)

        out_path = args.csv_out or (os.path.splitext(args.csv_in)[0] + ".tax_ids.csv")
        write_out_csv(out_path, results, extra_cols=["tax_id", "matched_name", "rank", "status"])
        print(f"Wrote: {out_path}")
    else:
        if not args.terms:
            sys.stderr.write("Provide --csv <file.csv> or positional species names.\n")
            sys.exit(2)

        for name in args.terms:
            try:
                tax_id, matched, rank = best_taxon_match(name, return_max=args.retmax)
                print(f"{name}\t{matched}\t{rank}\t{tax_id}")
            except Exception as e:
                print(f"{name}\t\t\tstatus={e}")
if __name__ == "__main__":
    main()