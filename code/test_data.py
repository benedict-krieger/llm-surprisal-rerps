import pandas as pd

df = read_csv("../data/adsbc21/adsbc21_surp_erp.csv")

elec = ["Fz","Cz","Pz"]

df = df[df["Timestamp"].isin(range(600,801))] # Python: exclusive upper bound

df.shape # (408141, 56)

len(df.TrialNum.unique()) # 4041

df[df["TrialNum"]==18]
# 101 rows -> one for each time stamp
# timestamps and voltages differ, rest is constant

df_bytrial =df.groupby(["TrialNum","Subject","Condition"])[elec].mean() # shape: 4041,3
df_byitem =df.groupby(["Item","Subject","Condition"])[elec].mean() # shape: 4041,3

df_bytrial.summary()
df_byitem.summary()
# They are the same