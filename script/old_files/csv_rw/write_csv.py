import pandas as pd



df = pd.read_csv("/Users/harivinaygujjula/Documents/GitHub/llm-metadata-mapping/data/mapping_BS_GS-mapped.csv")

if "chatgpt_name" not in df.columns:
    df["chatgpt_name"]= None

for i, row in df[df["chatgpt_name"].isna()]. iterrows():
    print(f"\nname: {row['name']}")
    ai_name = input("Enter ai_name (or press Enter to skip): ")

    if ai_name:
        df.loc[i, "chatgpt_name"] = ai_name

    df.to_csv("/Users/harivinaygujjula/Documents/Github/llm-metadata-mapping/data/updated_mapping_BS_GS-mapped.csv", index=False)
    print("\n Updates saved to mapping_BS_GS-mapped_updated.csv")
