import argparse
from pandas import read_csv

def build_arg_parser():
    parser = argparse.ArgumentParser(description="Read a CSV file and print preview.")
    parser.add_argument("--file", required=True, help="Path to the CSV file")
    parser.add_argument("--rows", type=int, default=5, help="How many rows to preview")
    return parser

def main():

    args = build_arg_parser().parse_args()


    df = read_csv(args.file)


    print(f"Rows: {len(df)}")
    print(f"Columns: {list(df.columns)}")
    print(df.head(args.rows))

if __name__ == "__main__":
    main()

