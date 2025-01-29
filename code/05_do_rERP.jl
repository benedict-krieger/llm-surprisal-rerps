# Christoph Aurnhammer, 2022
# rERPs have initially been described by Smith and Kutas (2015a, 2015b).

# Load the functions in rERPs.jl
# Packages are loaded from within rERPs.jl
include("rERP.jl");

function run_analysis(study_id, surp_id)
        println("Study id: $(study_id)")
        println("Surprisal id: $(surp_id)")
        # Define an array of electrodes on which to fit the models.
        elec = [:Fp1, :Fp2, :F7, :F3, :Fz, :F4, :F8, :FC5, :FC1, :FC2, :FC6, :C3,
                :Cz, :C4, :CP5, :CP1, :CP2, :CP6, :P7, :P3, :Pz, :P4, :P8, :O1, :Oz, :O2];

        # Define a models Structure. All arguments are arrays of column name Symbols, e.g. [:Item, :Subject]
        # Arguments:
        # - Descriptors: Fit separate models for each of these. You DO NOT have to add electrodes here!
        # - NonDescriptors: Columns that we don't use but want to carry over to our output data.
        # - Electrodes: Electrodes on which to fit separate models.
        # - Predictors: Predictors to use in the model.
        models = make_models([:Subject, :Timestamp], [:Item, :Condition], elec, [:Intercept, Symbol(surp_id)]);

        # Pre-process the data, using the following arguments:
        # - infile: Path to input file. Necessary.
        # - outfile: Path to output file. set to false to return DataFrame directly instead of writing. Necessary.
        # - models: Models structure specified above. Necessary.
        # - baseline_corr: Apply baseline correction, if you haven't done so already. Default false.
        # - sampling_rate: Downsample to a lower sampling rate. Default false.
        # - invert_preds: Takes an array of predictor Symbols (e.g., [:Cloze, :Association]) to invert. Default false.
        # - conds: Takes an array of condition labels to subset to (e.g. ["A", "B"]). Default false.
        # - components: Collect the average N400 and Segment amplitude (Special sauce for one of my projects). Default false.
        @time process_data("../data/$study_id/$(study_id)_surp_erp.csv", "../data/$study_id/$(study_id)_$(surp_id)_rERP.csv", models);


        #### Within subject

        # Read in processed data.
        @time dt = read_data("../data/$study_id/$(study_id)_$(surp_id)_rERP.csv", models);

        # Fit the rERP models, using the three arguments:
        # - data: data as processed by process_data()
        # - models: models Structure
        # - file: path to output file. Filenames will be automatically extended to *_data.csv and *_models.csv.
        @time fit_models(dt, models, "../data/$study_id/$(study_id)_$(surp_id)_rERP");


        #### Across subjects

        dt.Subject = ones(nrow(dt));
        @time fit_models(dt, models, "../data/$study_id/$(study_id)_$(surp_id)_across_subj_rERP");
end

#############################################################################################################################################

function dbc19_assocplaus()

        # Define an array of electrodes on which to fit the models.
        elec = [:Fp1, :Fp2, :F7, :F3, :Fz, :F4, :F8, :FC5, :FC1, :FC2, :FC6, :C3,
                :Cz, :C4, :CP5, :CP1, :CP2, :CP6, :P7, :P3, :Pz, :P4, :P8, :O1, :Oz, :O2];

        # Define a models Structure. All arguments are arrays of column name Symbols, e.g. [:Item, :Subject]
        models = make_models([:Subject, :Timestamp], [:Item, :Condition, :TrialNum], elec, [:Intercept, :Assoc, :Plaus]);

        # Pre-process the data, using the following arguments:
        # dbc19_erp.csv -> original dbc19 erp data
        @time process_data("../data/dbc19/dbc19_erp.csv", "../data/dbc19_corrected/dbc19_assocplaus_rERP.csv", models, invert_preds=[:Assoc, :Plaus]);

        # Read in processed data.
        @time dt = read_data("../data/dbc19_corrected/dbc19_assocplaus_rERP.csv", models);

        # Fit the rERP models _without_ averaging per time sample, condition and across subjects:
        @time fit_models_wo_avg(dt, models, "../data/dbc19_corrected/dbc19_assocplaus_rERP"); # y_hat = beta0 + beta1*assoc + beta2*plaus (both assoc and plaus standardized+inverted)

        # This will produce a very large _data.csv file (6.5 GB)


end

function dbc19_plausdata()

        @time dt = DataFrame(CSV.File("../data/dbc19_corrected/dbc19_assocplaus_rERP_data.csv"))
        dt_plaus = filter(row -> row.Spec == "[:Intercept, :Plaus]" && row.Type == "est", dt)
        dt_plaus = sort(dt_plaus, [:Item, :Condition]) # sort for item and condition in ascending order
        dt_plaus.Condition = convert(Vector{Any}, dt_plaus.Condition) # enable conversion from float to str
        replace!(dt_plaus.Condition, 1.0 => "A", 2.0 => "B", 3.0 => "C") # recode Condition column to A,B,C
        select!(dt_plaus, Not([:Spec, :Type]))# drop Spec & Type columns
        CSV.write("../data/dbc19_corrected/dbc19_corrected_erp.csv", dt_plaus)


end

#############################################################################################################################################
#############################################################################################################################################

study_ids = ["adsbc21","dbc19","adbc23","dbc19_corrected"]
surprisal_ids = ["leo13b_surp", "secretgpt2_surp", "gerpt2_surp", "gerpt2large_surp"]

for (study_id, surprisal_id) in Iterators.product(study_ids, surprisal_ids)
        run_analysis(study_id, surprisal_id)
end

#dbc19_assocplaus()
#dbc19_plausdata()