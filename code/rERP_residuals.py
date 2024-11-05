import pandas as pd

study_ids = ["adsbc21","dbc19","dbc19_corrected","adbc23"]
surp_ids = ["leo13b_surp","gerpt2large_surp","gerpt2_surp"]


def get_mean_resids():

    study_ids = ["adsbc21","dbc19","dbc19_corrected","adbc23"]
    surp_ids = ["leo13b_surp","gerpt2large_surp","gerpt2_surp"]
    time_windows = {
        "adsbc21" : {
            "N400" : (350,450),
            "P600" : (600,800)
        },
        "dbc19" : {
            "N400" : (300,500),
            "P600" : (800,1000)
        },
        "dbc19_corrected" : {
            "N400" : (300,500),
            "P600" : (800,1000)
        },
        "adbc23" : {
            "N400" : (300,500),
            "P600" : (600,1000)
        }    
    }

    cond_dict = {1.0:"A", 2.0:"B", 3.0:"C", 4.0:"D"}

    tws = ["N400","P600"]

    # initialize empty dataframe
    cols = ["study_id","time_window","surp_id","condition","Pz","Pz_across_cond"]
    rows = []

    for study_id in study_ids:
        print(study_id)
        for surp_id in surp_ids:
            print(surp_id)
            df = pd.read_csv(f"../data/{study_id}/{study_id}_{surp_id}_rERP_data.csv")

            for tw in tws:
                interval = time_windows[study_id][tw]
                print(f"{tw}: {interval}")

                df_window = df[df["Timestamp"].isin(range(interval[0],interval[1]+1))] # +1 because of exclusive upper bound
                df_window = df_window[((df_window["Type"] == "res") & (df_window["Spec"] == f"[:Intercept, :{surp_id}]"))] # we only want residuals here
                df_window = df_window[["Condition","Pz"]]
                df_window["Pz"] = df_window["Pz"].abs() # absolute residuals
                df_window["Condition"] = df_window["Condition"].replace(cond_dict) # re-code conditions
                # This data should already be averaged across subjects
                # Avg by condition, across time samples
                df_window = df_window.groupby(["Condition"], as_index=False).mean()
                for c,m in zip(df_window["Condition"],df_window["Pz"]):
                    rows.append({"study_id":study_id, "time_window":tw, "surp_id":surp_id, "condition":c, "Pz":round(m,2), "Pz_across_cond":round(df_window["Pz"].mean(),2)}) 
                #df_window.to_csv(f"../data/{study_id}/{study_id}_{tw}.csv", index=False)
    
    out_df = pd.DataFrame(rows, columns = cols)
    return out_df 

df = get_mean_resids()
df.to_csv("../results/mean_resids.csv",index=False)