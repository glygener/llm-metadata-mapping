import pandas as pd
from pathlib import Path

FOLDER = Path("/Users/harivinaygujjula/Documents/GitHub/llm-metadata-mapping/data")
FILES = [FOLDER / f"updated__BS_GS-mapped{i}.csv" for i in range(11, 21)]

MASTER_OUT = FOLDER / "updated_BS_GS_master_2.csv"

frames = []

print(f"Found {len(FILES)} CSV files to merge")

for f in FILES:
    print(f"Reading: {f.name}")
    df = pd.read_csv(f)
    frames.append(df)

if frames:
    master_df = pd.concat(frames, ignore_index=True)
    master_df.to_csv(MASTER_OUT, index=False)
    print(f"\nMerged {len(FILES)} files -> {MASTER_OUT}")
    print(f"Total rows in master file: {master_df}")
else:
    print("No new files to merge")
