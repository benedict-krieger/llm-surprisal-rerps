using CSV
using DataFrames
using MixedModels
using StatsBase
#using StatsModels

function load_data()

    println("Loading data")
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
    
    return studies

end



function predict_tws(studies::Dict{String, DataFrame})
    
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

    # Define LME formulas for time windows
    # N400
    n4_null_f = @formula N400 ~ z_Zipf_freq + z_Tw_position + (1|Item) + (1|Subject) + (1|Electrode)
    n4_cond_f = @formula N400 ~ Condition + z_Zipf_freq + z_Tw_position + (1|Item) + (1|Subject) + (1|Electrode)
    n4_leo13b_f = @formula N400 ~ z_leo13b_surp + z_Zipf_freq + z_Tw_position + (1|Item) + (1|Subject) + (1|Electrode)
    n4_gerpt2large_f = @formula N400 ~ z_gerpt2large_surp + z_Zipf_freq + z_Tw_position + (1|Item) + (1|Subject) + (1|Electrode)
    n4_gerpt2_f = @formula N400 ~ z_gerpt2_surp + z_Zipf_freq + z_Tw_position + (1|Item) + (1|Subject) + (1|Electrode)
    n4_formulas = [("null", n4_null_f), ("condition", n4_cond_f), ("leo13b", n4_leo13b_f), ("gerpt2large", n4_gerpt2large_f), ("gerpt2", n4_gerpt2_f)]

    # P600
    p6_null_f = @formula P600 ~ z_Zipf_freq + z_Tw_position + (1|Item) + (1|Subject) + (1|Electrode)
    p6_cond_f = @formula P600 ~ Condition + z_Zipf_freq + z_Tw_position + (1|Item) + (1|Subject) + (1|Electrode)
    p6_leo13b_f = @formula P600 ~ z_leo13b_surp + z_Zipf_freq + z_Tw_position + (1|Item) + (1|Subject) + (1|Electrode)
    p6_gerpt2large_f = @formula P600 ~ z_gerpt2large_surp + z_Zipf_freq + z_Tw_position + (1|Item) + (1|Subject) + (1|Electrode)
    p6_gerpt2_f = @formula P600 ~ z_gerpt2_surp + z_Zipf_freq + z_Tw_position + (1|Item) + (1|Subject) + (1|Electrode)
    p6_formulas = [("null", p6_null_f), ("condition", p6_cond_f), ("leo13b", p6_leo13b_f), ("gerpt2large", p6_gerpt2large_f), ("gerpt2", p6_gerpt2_f)]

    # Initialize data frame to store AICs
    aic_df = DataFrame(study = String[], time_window = String[], lme = String[], norm_aic = Float64[]) 

    for (study_id, df) in studies

        println("######## $(study_id) ########\n\n")

        # Z-transform predictors
        df.z_leo13b_surp = zscore(df.leo13b_surp)
        df.z_gerpt2large_surp = zscore(df.gerpt2large_surp)
        df.z_gerpt2_surp = zscore(df.gerpt2_surp)
        df.z_Zipf_freq = zscore(df.Zipf_freq)
        df.z_Tw_position = zscore(df.Tw_position)

        for (window, interval) in time_windows[study_id]

            println("#### $(window) ####\n")

            df_window = filter(row -> row.Timestamp >= interval[1] && row.Timestamp <= interval[2], df) # filter for time window
            df_window = stack(df_window, elec, cols; variable_name="Electrode", value_name=window) # melt data into long format, i.e. one col for elec

            window == "N400" ? formulas = n4_formulas : formulas = p6_formulas
            
            # Null model
            println("## Fitting $(formulas[1][1]) model ##\n")
            #println(formulas[1][2])
            null_model = fit(MixedModel, formulas[1][2], df_window)
            println(null_model)
            null_aic = aic(null_model)
            
            for f in formulas[2:end]
                model_name = f[1] 
                println("## Fitting $(model_name) model ## \n")
                model = fit(MixedModel, f[2], df_window)
                println("$(model)\n\n")
                model_aic = aic(model)

                aic_diff = model_aic - null_aic # normalize by null model
                push!(aic_df, (study_id,window,model_name,aic_diff))
            end

        end

    end

    return (aic_df)
end

studies = @time load_data()
aic_df = @time predict_tws(studies)
# Filter so that we have AICs for dbc19 only for the N4 and for dbc19_corrected only for the P6 windows
aic_df = filter(row -> !((row.study == "dbc19_corrected" && row.time_window == "N400") || (row.study == "dbc19" && row.time_window == "P600")), aic_df)
aic_df.study = replace(aic_df.study, "dbc19_corrected" => "dbc19")
aic_df[9:16, :] = vcat(aic_df[13:16, :], aic_df[9:12, :]) # Swapping rows to maintain N4/P6 order

CSV.write("../results/erp_aic/aic_diffs.csv", aic_df)