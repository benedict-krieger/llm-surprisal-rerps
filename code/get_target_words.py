import pandas as pd

studies = ["adsbc21","dbc19","adbc23"]

for s in studies:
    print(s)
    study_df = pd.read_csv(f"../data/{s}/{s}.csv", sep = ";")
    target_list = list(study_df["Target"].unique())
    print(len(target_list))
    with open(f"../data/{s}/{s}_target_words.txt", mode = "w", encoding = "utf-8") as tw_file:
        tw_file.write("\n".join(target_list))