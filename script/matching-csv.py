import pandas as pd
from pathlib import Path

data_path = Path("/Users/harivinaygujjula/Documents/GitHub/llm-metadata-mapping/data/updated_BS_GS_master2.csv")
df = pd.read_csv(data_path)
df = pd.read_csv("updated_BS_GS_master_2.csv")
df["match"] = df ["namespacename"] == df ["chatgpt_name"]
df["match"] = df["match"].map({True: "Y", False: "N"})

df.to_csv("updated_BS_GS_master_2match.csv", index=False)
