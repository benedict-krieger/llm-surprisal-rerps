using CSV
using DataFrames
using MixedModels
using StatsBase

function load_data()

    println("Loading data")
    adsbc21_n4_df = CSV.read("../data/adsbc21/adsbc21_N400.csv", DataFrame)
    adsbc21_p6_df = CSV.read("../data/adsbc21/adsbc21_P600.csv", DataFrame)
    
    dbc19_n4_df = CSV.read("../data/dbc19/dbc19_N400.csv", DataFrame)
    dbc19_p6_df = CSV.read("../data/dbc19/dbc19_P600.csv", DataFrame)
    
    dbc19_corrected_n4_df = CSV.read("../data/dbc19_corrected/dbc19_corrected_N400.csv", DataFrame)
    dbc19_corrected_p6_df = CSV.read("../data/dbc19_corrected/dbc19_corrected_P600.csv", DataFrame)
    
    adbc23_n4_df = CSV.read("../data/adbc23/adbc23_N400.csv", DataFrame)
    adbc23_p6_df = CSV.read("../data/adbc23/adbc23_P600.csv", DataFrame)

    studies = Dict(
        "adsbc21" => Dict(
            "N400" => adsbc21_n4_df,
            "P600" => adsbc21_p6_df
        ),
        "dbc19" => Dict(
            "N400" => dbc19_n4_df,
            "P600" => dbc19_p6_df
        ),
        "dbc19_corrected" => Dict(
            "N400" => dbc19_corrected_n4_df,
            "P600" => dbc19_corrected_p6_df
        ),
        "adbc23" => Dict(
            "N400" => adbc23_n4_df,
            "P600" => adbc23_p6_df
        )
        )
    
    return studies

end



function predict_tws(studies)
    
    # Define LME formulas for time windows
    # N400
    n4_null_f = @formula N400_mean ~ z_Zipf_freq + z_Tw_position + (1|Item) + (1|Subject) + (1|Electrode)
    n4_cond_f = @formula N400_mean ~ Condition + z_Zipf_freq + z_Tw_position + (1|Item) + (1|Subject) + (1|Electrode)
    n4_leo13b_f = @formula N400_mean ~ z_leo13b_surp + z_Zipf_freq + z_Tw_position + (1|Item) + (1|Subject) + (1|Electrode)
    n4_gerpt2large_f = @formula N400_mean ~ z_gerpt2large_surp + z_Zipf_freq + z_Tw_position + (1|Item) + (1|Subject) + (1|Electrode)
    n4_gerpt2_f = @formula N400_mean ~ z_gerpt2_surp + z_Zipf_freq + z_Tw_position + (1|Item) + (1|Subject) + (1|Electrode)
    n4_formulas = [("null", n4_null_f), ("condition", n4_cond_f), ("leo13b", n4_leo13b_f), ("gerpt2large", n4_gerpt2large_f), ("gerpt2", n4_gerpt2_f)]

    # P600
    p6_null_f = @formula P600_mean ~ z_Zipf_freq + z_Tw_position + (1|Item) + (1|Subject) + (1|Electrode)
    p6_cond_f = @formula P600_mean ~ Condition + z_Zipf_freq + z_Tw_position + (1|Item) + (1|Subject) + (1|Electrode)
    p6_leo13b_f = @formula P600_mean ~ z_leo13b_surp + z_Zipf_freq + z_Tw_position + (1|Item) + (1|Subject) + (1|Electrode)
    p6_gerpt2large_f = @formula P600_mean ~ z_gerpt2large_surp + z_Zipf_freq + z_Tw_position + (1|Item) + (1|Subject) + (1|Electrode)
    p6_gerpt2_f = @formula P600_mean ~ z_gerpt2_surp + z_Zipf_freq + z_Tw_position + (1|Item) + (1|Subject) + (1|Electrode)
    p6_formulas = [("null", p6_null_f), ("condition", p6_cond_f), ("leo13b", p6_leo13b_f), ("gerpt2large", p6_gerpt2large_f), ("gerpt2", p6_gerpt2_f)]

    # Initialize data frame to store AICs
    aic_df = DataFrame(study = String[], time_window = String[], lme = String[], norm_aic = Float64[]) 
    

    for (study_id, tw_dict) in pairs(studies)

        println("\n######## $(study_id) ########\n")

        for (tw,df) in pairs(tw_dict)

            println("#### $(tw) ####")
            
            # Z-transform continuous predictors
            df.z_leo13b_surp = zscore(df.leo13b_surp)
            df.z_gerpt2large_surp = zscore(df.gerpt2large_surp)
            df.z_gerpt2_surp = zscore(df.gerpt2_surp)
            df.z_Zipf_freq = zscore(df.Zipf_freq)
            df.z_Tw_position = zscore(df.Tw_position)

            tw == "N400" ? formulas = n4_formulas : formulas = p6_formulas
        
            # Store fitted models
            fitted_models = Any[]

            # Null model
            println("## Fitting $(formulas[1][1]) model ##")
            null_model = fit(MixedModel, formulas[1][2], df)
            null_aic = aic(null_model)
            push!(fitted_models, null_model)

            for f in formulas[2:end]
                model_name = f[1] 
                println("## Fitting $(model_name) model ## ")
                model = fit(MixedModel, f[2], df)
                model_aic = aic(model)
                push!(fitted_models, model)

                aic_diff = model_aic - null_aic # normalize by null model
                push!(aic_df, (study_id,tw,model_name,aic_diff))
            end
            
            open("../results/erp_aic/$(study_id)_$(tw)_lme.txt", "w") do file
                for m in fitted_models
                    write(file,string(m))
                    write(file,"\n\n\n\n#####################################\n\n\n\n")
                end
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