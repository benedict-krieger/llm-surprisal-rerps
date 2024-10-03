using DataFrames
using CSV
using StatsBase
using MixedModels

# Loading data frames
adsbc21_df = CSV.read("../data/adsbc21/adsbc21_surp_erp.csv", DataFrame)
dbc19_df = CSV.read("../data/dbc19/dbc19_surp_erp.csv", DataFrame)
dbc19_corrected_df = CSV.read("../data/dbc19_corrected/dbc19_corrected_surp_erp.csv", DataFrame)
adbc23_df = CSV.read("../data/adbc23/adbc23_surp_erp.csv", DataFrame)

studies = Dict(
    "adsbc21" => adsbc21_df,
    "dbc19" => dbc19_df,
    "dbc19_corrected" => dbc19_corrected_df,
    "adbc23" => adbc23_df
    )

time_windows = Dict(
    "adsbc21" => Dict("N400" => (350,450), "P600" => (600,800)),
    "dbc19" => Dict("N400" => (300,500), "P600" => (800,1000)),
    "dbc19_corrected" => Dict("N400" => (300,500), "P600" => (800,1000)),
    "adbc23" => Dict("N400" => (300,500), "P600" => (600,1000))    
)

elec = [:Fp1, :Fp2, :F7, :F3, :Fz, :F4, :F8, :FC5, :FC1, :FC2, :FC6, :C3,
        :Cz, :C4, :CP5, :CP1, :CP2, :CP6, :P7, :P3, :Pz, :P4, :P8, :O1, :Oz, :O2]
cols = [:Item, :Condition, :Subject, :Timestamp, :z_Zipf_freq,
        :z_Tw_position, :z_leo13b_surp, :z_gerpt2_surp, :z_gerpt2large_surp]


aic_diffs = Dict{String, Float64}() # Initialize to store AICs

for (study_id, df) in studies

    println(study_id)

    # Z-transform predictors
    df.z_leo13b_surp = zscore(df.leo13b_surp)
    df.z_gerpt2large_surp = zscore(df.gerpt2large_surp)
    df.z_gerpt2_surp = zscore(df.gerpt2_surp)
    df.z_Zipf_freq = zscore(df.Zipf_freq)
    df.z_Tw_position = zscore(df.Tw_position)

    for (window, interval) in time_windows[study_id]

    df_window = filter(row -> row.Timestamp >= interval[1] && row.Timestamp <= interval[2], df) # filter for time window
    df_window = stack(df, elec, cols; variable_name="Electrode", value_name=window) # melt data into long format, i.e. one col for elec

        if window == "N400"
            null_f = @formula N400 ~ z_Zipf_freq + z_Tw_position + (1|Item) + (1|Subject) + (1|Electrode)
            cond_f = @formula N400 ~ Condition + z_Zipf_freq + z_Tw_position + (1|Item) + (1|Subject) + (1|Electrode)
            leo13b_f = @formula N400 ~ z_leo13b_surp + z_Zipf_freq + z_Tw_position + (1|Item) + (1|Subject) + (1|Electrode)
            gerpt2large_f = @formula N400 ~ z_gerpt2large_surp + z_Zipf_freq + z_Tw_position + (1|Item) + (1|Subject) + (1|Electrode)
            gerpt2_f = @formula N400 ~ z_gerpt2_surp + z_Zipf_freq + z_Tw_position + (1|Item) + (1|Subject) + (1|Electrode)
            
        else
            null_f = @formula P600 ~ z_Zipf_freq + z_Tw_position + (1|Item) + (1|Subject) + (1|Electrode)
            cond_f = @formula P600 ~ Condition + z_Zipf_freq + z_Tw_position + (1|Item) + (1|Subject) + (1|Electrode)
            leo13b_f = @formula P600 ~ z_leo13b_surp + z_Zipf_freq + z_Tw_position + (1|Item) + (1|Subject) + (1|Electrode)
            gerpt2large_f = @formula P600 ~ z_gerpt2large_surp + z_Zipf_freq + z_Tw_position + (1|Item) + (1|Subject) + (1|Electrode)
            gerpt2_f = @formula P600 ~ z_gerpt2_surp + z_Zipf_freq + z_Tw_position + (1|Item) + (1|Subject) + (1|Electrode)
        end

    formulas = [("Condition",cond_f), ("Leo-13b",leo13b_f), ("GerPT-2 large",gerpt2large_f), ("GerPT-2",gerpt2_f)]
    
    # Null model
    null_model = fit(MixedModel, null_f, df_window)
    null_aic = aic(null_model)
    
        for f in formulas
            model_name = f[1] 
            model = fit(MixedModel, f[2], df_window)
            model_aic = aic(model)

            aic_diff = model_aic - null_aic # normalize by null model
            id_key = "$(study_id)_$(window)_$(model_name)"
            aic_diffs[id_key] = aic_diff
        end

    end
end

println(aic_diffs)


##############
#### N400 ####
##############

#df_n4 = filter(row -> row.Timestamp >= 350 && row.Timestamp <= 450, df) # filter for N400 window
#df_n4 = stack(df_n4, elec, cols; variable_name="Electrode", value_name="N400") # melt data into long format, i.e. one col for elec

# Null regression
#n4_null_formula = @formula N400 ~ z_Zipf_freq + z_Tw_position + (1|Item) + (1|Subject) + (1|Electrode)
#n4_null_model = fit(MixedModel, n4_null_formula, df_n4)
#n4_null_aic = aic(n4_null_model)

# Condition
#n4_cond_formula = @formula N400 ~ Condition + z_Zipf_freq + z_Tw_position + (1|Item) + (1|Subject) + (1|Electrode)
#n4_cond_model = fit(MixedModel, n4_cond_formula, df_n4)
#n4_cond_aic = aic(n4_cond_model)
#n4_cond_aic_norm = n4_cond_aic - n4_null_aic
#println(n4_cond_aic_norm)

# Leo-13b
#n4_leo13b_formula = @formula N400 ~ z_leo13b_surp + z_Zipf_freq + z_Tw_position + (1|Item) + (1|Subject) + (1|Electrode)
#n4_leo13b_model = fit(MixedModel, n4_leo13b_formula, df_n4)
#n4_leo13b_aic = aic(n4_leo13b_model)
#n4_leo13b_aic_norm = n4_leo13b_aic - n4_null_aic
#println(n4_leo13b_aic_norm)

# GerPT-2 large
#n4_gerpt2large_formula = @formula N400 ~ z_gerpt2large_surp + z_Zipf_freq + z_Tw_position + (1|Item) + (1|Subject) + (1|Electrode)
#n4_gerpt2large_model = fit(MixedModel, n4_gerpt2large_formula, df_n4)
#n4_gerpt2large_aic = aic(n4_gerpt2large_model)
#n4_gerpt2large_aic_norm = n4_gerpt2large_aic - n4_null_aic
#println(n4_gerpt2large_aic_norm)

# GerPT-2
#n4_gerpt2_formula = @formula N400 ~ z_gerpt2_surp + z_Zipf_freq + z_Tw_position + (1|Item) + (1|Subject) + (1|Electrode)
#n4_gerpt2_model = fit(MixedModel, n4_gerpt2_formula, df_n4)
#n4_gerpt2_aic = aic(n4_gerpt2_model)
#n4_gerpt2_aic_norm = n4_gerpt2_aic - n4_null_aic
#println(n4_gerpt2_aic_norm)


#function ()
#foreach(println, names(df_n400)) # show column names