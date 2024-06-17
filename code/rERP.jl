##
# Christoph Aurnhammer <aurnhammer@coli.uni-saarland.de>
#
# Functions implementing rERP analysis
##

using DataFrames
using Combinatorics: combinations
using CSV: File, write
using Distributions: cdf, TDist
using StatsBase: mean, zscore, std
using LinearAlgebra: cholesky, diag
using CategoricalArrays: cut, levelcode

struct Models
    Descriptors::Array
    NonDescriptors::Array
    Electrodes::Array
    Predictors::Array
    Sets::Array
    Quantiles::Bool
end

function make_models(desc, nondesc, elec, pred; quant = false)
    models = Models(desc, nondesc, elec, pred, [], quant);
    Models(models.Descriptors, models.NonDescriptors, models.Electrodes, models.Predictors, pred_sets(models), models.Quantiles);
end

struct Ind
    numpred::Int
    pred_dict::Dict
    n::Int
    mn::Int
    s::Int
    e::Int
    m::Int
end

function make_Ind(data, models, m_indices, s_ind, e_ind, m_ind) 
    Ind(length(models.Predictors),  Dict([x => i-1 for (i, x) in enumerate(models.Predictors)]), nrow(data), length(m_indices), s_ind, e_ind, m_ind)
end

function process_data(infile, outfile, models; baseline_corr = false, sampling_rate = false, invert_preds = false, conds = false, components = false, keep_conds = false, time_windows = false)
    # Load Data from disk
    data = DataFrame(File(infile))
    
    # Make ItemNum vs Item coherent across datasets
    if "ItemNum" in names(data)
        rename!(data, :ItemNum => :Item);
    end

    # Downsample
    if sampling_rate != false
        data = downsample(data, sampling_rate)
    end

    # Take condition subsets at this point (i.e. before any z-scoring takes place)
    if conds != false
        data = data[subset_inds(data, conds),:]
    end

    # Select columns
    data.Intercept = ones(nrow(data));
    data = data[:,vcat(models.Descriptors, models.NonDescriptors, models.Predictors, models.Electrodes)];

    # Baseline correction. All our in house data has already been baseline corrected before export.
    if baseline_corr == true
        data = baseline_correction(data, models);
    end

    # Collect component predictor
    if components != false
        # Add an electrode specific component predictors
        data = collect_component(data, "N400", models ; tws=300, twe=500);
        data = collect_component(data, "Segment", models ; tws=0, twe=1200);
    end
    
    # Z-standardise predictors
    data = standardise(data, components, models);
    
    # Invert predictors
    if ((invert_preds != false) | (components != false))
        data = invert(data, components, models, invert_preds)
    end

    # Turn condition labels to numbers. Set Verbose to show them.
    if keep_conds == false
        data = transform_conds(data, verbose=true);
    end

    # Write data to file or return as DataFrame
    if typeof(outfile) == String
        write(outfile, data)
    else
        sort!(data, [x for x in reverse(models.Descriptors)])
    end
end

function downsample(data, sampling_rate)
    factor = Int8(1000 / sampling_rate)
    data = data[ [x in range(-200, stop=1199, step=factor) for x in data[!,:Timestamp]] , :]
    
    data
end

function subset_inds(data, conds)
    ind = Int.(zeros(nrow(data)))
    for x in conds
        ind_curr = (data.Condition .== x)
        ind = ind .| ind_curr
    end

    Bool.(ind)
end

function baseline_correction(data, models)
    # Collect Baseline amplitudes
    base = @view data[(data.Timestamp .< 0), vcat([:Subject, :Item, :Timestamp], models.Electrodes)];
    base = combine(groupby(base, [:Subject, :Item]), [x => mean => Symbol(x, "base") for x in models.Electrodes])
    data = innerjoin(data, base, on = [:Subject, :Item])
    
    # Baseline correct EEG data
    data[:,models.Electrodes] = data[:,models.Electrodes] .- Array(data[:,[Symbol(x, "base") for x in models.Electrodes]])

    data
end

function collect_component(data, name, models ; tws = 300, twe = 500)
        # Collect component amplitudes
        comp = @view data[((data.Timestamp .>= tws) .& (data.Timestamp .<= twe)), vcat([:Subject, :Item, :Timestamp], models.Electrodes)];
        comp = combine(groupby(comp, [:Subject, :Item]), [x => mean => Symbol(x, name) for x in models.Electrodes]);  
        data = innerjoin(data, comp, on = [:Subject, :Item]);
        data
end

function standardise(data, components, models)
    if length(models.Predictors) > 1
        for x in models.Predictors[2:end]
            data[!,x] = zscore(data[!,x])
        end
    elseif components != false
        for x in models.Electrodes
            for comp in components
                data[:,Symbol(x, comp)] = zscore(data[:,Symbol(x, comp)])
            end
        end
    end

    data
end

function invert(data, components, models, invert_preds)
    if components != false
        for x in models.Electrodes
            for comp in components
                data[:,Symbol(x, comp)] = data[:,Symbol(x, comp)] .* -1
            end
        end
    else     
        for x in invert_preds
            data[!,x] = (data[!,x]) .* -1
        end
    end
    
    data
end

function process_spr_data(infile, outfile, models; invert_preds = false, conds = false, keep_conds = false, logRT = true) 
    # Load data from disk
    data = DataFrame(File(infile));

    # Make ItemNum vs Item coherent across datasets
    if "ItemNum" in names(data)
        rename!(data, :ItemNum => :Item);
    end
    
    # Log-transform RT
    if logRT == true
        data.logRT = log.(data.ReadingTime)
    end

    # Pretend that RT regions are ERP timestamps
    rename!(data, :Region => :Timestamp)

    # Take condition subsets at this point (i.e. before any z-scoring takes place)
    if conds != false
        data = data[subset_inds(data, conds),:]
    end

    # Select columns
    data.Intercept = ones(nrow(data));
    data = data[:,vcat(models.Descriptors, models.NonDescriptors, models.Predictors, models.Electrodes)];
    
    # Z-standardise predictors
    data = standardise(data, false, models);
    
    # Invert predictors
    if ((invert_preds != false))
        data = invert(data, false, models, invert_preds)
    end

    # Turn condition labels to numbers. Set Verbose to show them.
    if keep_conds == false
        data = transform_conds(data, verbose=true);
    end
    # Write data to file or return as DataFrame
    if typeof(outfile) == String
        write(outfile, data)
    else
        sort!(data, [x for x in reverse(models.Descriptors)])
    end
end

function exclude_trial(df, lower, upper, lower_rc, upper_rc)
    df_out = df
    for s in unique(df[!,:Subject])
        # subset for current subject
        dts = df[(df.Subject .== s),:];
        for i in unique(dts[!,:Item])
            # subset current item
            dti = dts[(dts.Item .== i),:];
            
            if sum(dti.ReactionTime .=== missing) > 1
                exc = sum((dti.ReadingTime .> upper) .| (dti.ReadingTime .< lower)) > 0;
                if exc == true
                    #print(dti[:,[:ReadingTime]], "\n")
                    df_out = df_out[.!((df_out.Subject .== s) .& (df_out.Item .== i)),:];
                end
            elseif sum(dti.ReactionTime .=== missing) == 0
                exc = sum((dti.ReadingTime .> upper) .| (dti.ReadingTime .< lower) .| (dti.ReactionTime .> upper_rc) .| (dti.ReactionTime .< lower_rc)) > 0;
                if exc == true
                    #print(dti[:,[:ReactionTime, :ReadingTime]], "\n")
                    df_out = df_out[.!((df_out.Subject .== s) .& (df_out.Item .== i)),:];
                end
            end
        end
    end

    println(nrow(df)-nrow(df_out), " / ", nrow(df), " excluded (", (nrow(df)-nrow(df_out))*100/nrow(df) ,"%).")

    df_out
end

function read_data(infile, models)
    data = DataFrame(File(infile))

    sort!(data, [x for x in reverse(models.Descriptors)])
end

function fit_models(data, models, file)
    # Get subset indices
    Index = get_index(data, models)
    e_indices = findall(Index .!= 0);
    s_indices = flatten_ar([1, e_indices[1:end-1].+1]);
    m_indices = 1:length(e_indices);
    ind = make_Ind(data, models, m_indices, 0, 0, 0);

    # allocate output data frames
    out_data = allocate_data(data, models);
    out_models = allocate_models(data, models, ind);

    # Get number of models, for showing off.
    num = num_mod(data, models)
    print("Fitting ", num, " models using ", Threads.nthreads(), " threads. \n")   
    Threads.@threads for i in 1:length(s_indices)
	    local ind = make_Ind(data, models, m_indices, s_indices[i], e_indices[i], m_indices[i]);
        
        # Take subset
        df = @view data[ind.s:ind.e,:];

        # Insert coefficients.
        out_models = coef(out_models, df, models, ind);

        # Insert estimates.
        out_data = estimates(out_models, out_data, df, models, ind);
        
        # Insert residuals.
        out_data = residual(out_data, models, ind);
        
        # Insert SE on coefficients.
        out_models = standarderror(out_data, out_models, df, models, ind);
    end
    
    # compute t-values.
    out_models = tvalue(out_models, models, ind);

    # compute p-values.
    out_models = pvalue(out_models, models, ind);

    # Addition of intercept to coefs
    out_models = coef_addition(out_models, models, ind)

    if typeof(file) == String
        out_models = write_models(out_models, models, file)
        out_data = write_data(out_data, models, file)
    end

    [out_data, out_models]
end

function fit_models_components(dt, models, file)
    out_data = []
    out_models = []
    for (i, e) in enumerate(models.Electrodes)
        println("Electrode $e ")
        models_e = make_models(models.Descriptors, models.NonDescriptors, [e], [:Intercept, Symbol(e, "N400"), Symbol(e, "Segment")]; quant = models.Quantiles);
        output = fit_models(dt, models_e, "none")
        if i .== 1
            out_data = output[1]
            out_models = output[2]
        else
            out_data[!,e] = output[1][:,e]
            out_data[!,Symbol(e, "_CI")] = output[1][:,Symbol(e, "_CI")]
            out_models[!,e] = output[2][:,e]
            out_models[!,Symbol(e, "_CI")] = output[2][:,Symbol(e, "_CI")]     
        end
    end
    
    write(string("../data/", file, "_data.csv"), out_data)
    write(string("../data/", file, "_models.csv"), out_models)    
end

function transform_conds(data ; verbose=false, column = :Condition)
    # Turn conditions into numbers
    cond_labels = unique(data[:, column]);
    cond_dict = Dict(zip(cond_labels, [x for x in 1:length(cond_labels)]));
    data[!, column] = [cond_dict[x] for x in data[:, column]];

    if verbose
        println(cond_dict)
    end

    data
end

function get_index(data, models)
    if length(unique(data.Subject)) > 1
        Index = flatten_ar([data[[x for x in 1:nrow(data)-1] .+ 1, models.Descriptors[1]] - data[[x for x in 1:nrow(data)-1], models.Descriptors[1]], 1]);
    elseif length(unique(data.Subject)) == 1
        Index = flatten_ar([data[[x for x in 1:nrow(data)-1] .+ 1, models.Descriptors[2]] - data[[x for x in 1:nrow(data)-1], models.Descriptors[2]], 1]);
    end

    Index
end

function allocate_data(data, models)
    nperms = length(models.Sets);

    # Allocate out data.
    out_data_names = vcat(models.Descriptors, models.NonDescriptors, models.Predictors, [:Type, :Spec], models.Electrodes);
    
    out_nrow = nrow(data) + 2 * nperms * nrow(data); # (res + est) * nperms * nrow(data)
    
    out_data = DataFrame(zeros(out_nrow, length(out_data_names)), :auto);
    rename!(out_data, out_data_names);

    # add original data
    ind = (length(models.Descriptors)+length(models.NonDescriptors) + length(models.Predictors));
    out_data[1:nrow(data),vcat(out_data_names[1:ind], out_data_names[ind+3:end])] = @view data[!,vcat(out_data_names[1:ind], out_data_names[ind+3:end])];

    # add value 42 for original eeg data
    out_data[1:nrow(data),[:Type, :Spec]] = hcat(repeat([1], nrow(data)), repeat([42], nrow(data)))

    out_data
end

function allocate_models(data, models, ind)
    # Allocate out models
    out_models_names = vcat(models.Descriptors, [:Type, :Spec, :N], models.Electrodes);
    out_rows = length(unique(data.Subject)) * length(unique(data.Timestamp)) * (ind.numpred * 4)
    out_models = DataFrame(zeros(out_rows, length(out_models_names)), :auto);

    rename!(out_models, out_models_names);
end

function pred_sets(models)
    flatten_ar([[[:Intercept]], [flatten_ar([[:Intercept], x]) for x in combinations(models.Predictors[2:end])]])
end

function num_mod(data, models)
    num_mod = length(models.Electrodes)
    for desc in models.Descriptors
        num_mod = num_mod * length(unique(data[:,desc]))
    end

    num_mod
end

function coef(out_models, df, models, ind)
    # Fit model
    coefs = Array(df[:,models.Predictors]) \ Array(df[:,models.Electrodes]);

    # Insert coefficients. out_models.Type 1
    for pred in models.Predictors
        out_models[ind.m+ind.mn*ind.pred_dict[pred],:] = vcat(Array(df[1,models.Descriptors]), [1, ind.pred_dict[pred], nrow(df)], coefs[ind.pred_dict[pred]+1,:]);
    end

    out_models
end

function estimates(out_models, out_data, pred_df, models, ind)
    for (s_num, s) in enumerate(models.Sets)
        for pred in s
            out_data[ind.s+ind.n*s_num:ind.e+ind.n*s_num,models.Electrodes] = out_data[ind.s+ind.n*s_num:ind.e+ind.n*s_num,models.Electrodes] .+ DataFrame(out_models[ind.m+ind.mn*ind.pred_dict[pred],models.Electrodes]) .* pred_df[:,pred]
        end
        out_data[ind.s+ind.n*s_num:ind.e+ind.n*s_num,[:Type, :Spec]] = hcat(repeat([2], nrow(pred_df)), repeat([s_num], nrow(pred_df)));
        out_data[ind.s+ind.n*s_num:ind.e+ind.n*s_num,flatten_ar([models.Descriptors, models.NonDescriptors, models.Predictors])] = out_data[ind.s:ind.e,flatten_ar([models.Descriptors, models.NonDescriptors, models.Predictors])]
    end

    out_data
end

function residual(out_data, models, ind)
    for s_num in 1:length(models.Sets)
        d_start_est = ind.s+ind.n*s_num;
        d_end_est = ind.e+ind.n*s_num;
        d_start_res = d_start_est + ind.n * length(models.Sets)
        d_end_res = d_end_est + ind.n * length(models.Sets);
        out_data[d_start_res:d_end_res,models.Electrodes] = out_data[ind.s:ind.e,models.Electrodes] .- out_data[d_start_est:d_end_est,models.Electrodes]
        out_data[d_start_res:d_end_res,[:Type, :Spec]] = hcat(repeat([3], length(ind.s:ind.e)), repeat([s_num], length(ind.s:ind.e)))
        out_data[d_start_res:d_end_res, flatten_ar([models.Descriptors, models.NonDescriptors, models.Predictors])] = out_data[ind.s:ind.e,flatten_ar([models.Descriptors, models.NonDescriptors, models.Predictors])]
    end

    out_data
end

function standarderror(out_data, out_models, df, models, ind)
    res = out_data[ind.s + ind.n * 2 * length(models.Sets): ind.e + ind.n * 2 * length(models.Sets),models.Electrodes]
    preds = Array(@view df[:,models.Predictors])
    
    std_error = sqrt.(transpose(sum.(eachcol(res.^2)) ./ (nrow(df) - ind.numpred)) .* diag(inv(transpose(preds) * preds)))

    offset = 0
    for (p_num, p) in enumerate(models.Predictors)
        out_models[offset+ind.m+ind.mn*ind.numpred+p_num-1,:] = flatten_ar([Array(df[1,models.Descriptors]), [2,ind.pred_dict[p],nrow(df)], std_error[p_num,:]])
        offset += ind.mn-1
    end

    out_models
end

function coef_addition(out_models, models, ind)
    for (p_num, p) in enumerate(models.Predictors[2:end])
        out_models[p_num*ind.mn+1:(p_num+1)*ind.mn,models.Electrodes] = out_models[p_num*ind.mn+1:(p_num+1)*ind.mn,models.Electrodes] .+ out_models[1:ind.mn,models.Electrodes]
    end

    out_models
end

function tvalue(out_models, models, ind)
    m_start = ind.numpred * ind.mn * 2 + 1
    m_end = ind.numpred * ind.mn * 3 

    out_models[m_start:m_end,models.Electrodes] = (@view out_models[(out_models.Type .== 1),models.Electrodes]) ./ (@view out_models[(out_models.Type .== 2),models.Electrodes]);
    out_models[m_start:m_end,[:Subject, :Timestamp, :N]] = @view out_models[(out_models.Type .== 2),[:Subject, :Timestamp, :N]]
    out_models[m_start:m_end,[:Type, :Spec]] = hcat(repeat([3], m_end-m_start+1), @view out_models[(out_models.Type .== 2),:Spec])

    out_models
end

function pvalue(out_models, models, ind)
    n = @view out_models[(out_models.Type .== 3),:N];
    m_start = ind.numpred * ind.mn * 3 + 1
    m_end = ind.numpred * ind.mn * 4 

    out_models[m_start:m_end,models.Electrodes] = 2 .* (1 .- cdf.(TDist.(n .- ind.numpred), abs.(@view out_models[(out_models.Type .== 3),models.Electrodes])));
    out_models[m_start:m_end,[:Subject, :Timestamp, :N]] = @view out_models[(out_models.Type .== 2), [:Subject, :Timestamp, :N]];
    out_models[m_start:m_end,[:Type, :Spec]] = hcat(repeat([4], m_end-m_start+1), @view out_models[(out_models.Type .== 2),:Spec]);

    out_models
end

function flatten_ar(arrays)
    collect(Iterators.flatten(arrays))
end

function count_sig(pvals; alpha = 0.05)
    sum(pvals .< alpha) / length(pvals)
end

function se(x)
    std(x) / sqrt(length(x))
end

# The sample mean plus or minus 1.96 times its standard error gives ci
function ci(x)
    1.96 * se(x)
end

function write_models(out_models, models, file)
    out_models_cp = out_models[:,:]
    out_models = combine(groupby(out_models_cp, [:Timestamp, :Type, :Spec]), [x => mean => x for x in models.Electrodes]);

    for x in models.Electrodes
        out_models[:,Symbol(x,"_CI")] = zeros(nrow(out_models))
        out_models[out_models.Type .== 1,Symbol(x,"_CI")] = out_models[out_models.Type .== 2,x];
        out_models[out_models.Type .== 3,Symbol(x, "_CI")] = combine(groupby(out_models_cp[(out_models_cp.Type .== 3),:], [:Timestamp, :Type, :Spec]), [x => ci => x])[!,x]
        out_models[out_models.Type .== 4,Symbol(x, "_CI")] = combine(groupby(out_models_cp[(out_models_cp.Type .== 4),:], [:Timestamp, :Type, :Spec]), [x => count_sig => x])[!,x]
    end

    # delete SE rows (they are now cols next to their coef)
    out_models = @view out_models[out_models.Type .!= 2,:];
    
    mtype_dict = Dict(1 => "Coefficient", 2 => "SE", 3 => "t-value", 4 => "p-value");
    mspec_dict = Dict([i-1 => x for (i, x) in enumerate(models.Predictors)])
    out_models[!,:Type] = [mtype_dict[x] for x in out_models[:,:Type]];
    out_models[!,:Spec] = [mspec_dict[x] for x in out_models[:,:Spec]];

    if file != "none"
        write(string("../data/", file, "_models.csv"), out_models)
    end

    out_models
end

function assign_estimate_quantiles(data, models)
    est_n4 = data[((data.Type .== 2.0) .& (data.Spec .== 2.0)),:]
    est_n4 = collect_component(est_n4, "N400_est", models; tws=300, twe=500);
    est_n4.Quantile = levelcode.(cut(est_n4.PzN400_est, 3))
    data.Condition = repeat(est_n4.Quantile, (2 * length(unique(models.Sets)) + 1))

    data
end

function write_data(out_data, models, file)
    if models.Quantiles
        out_data = assign_estimate_quantiles(out_data, models)
    end

    out_data_cp = out_data[:,:]
    out_data = combine(groupby(out_data_cp, [:Timestamp, :Type, :Spec, :Condition]), [x => mean => x for x in models.Electrodes])
    
    for x in models.Electrodes
        out_data[!,Symbol(x, "_CI")] = zeros(nrow(out_data)); 
        out_data_subj = combine(groupby(out_data_cp, [:Timestamp, :Type, :Spec, :Condition, :Subject]), [x => mean => x])
        error_ribbons = combine(groupby(out_data_subj, [:Timestamp, :Type, :Spec, :Condition]), [x => ci => x])
        # Mysteriously, the output df of the previous line is not always sorted the same way. Hence, sort.
        sort!(error_ribbons, [:Timestamp, :Type, :Spec, :Condition])
        out_data[!,Symbol(x, "_CI")] = error_ribbons[!,x]
    end

    dtype_dict = Dict(1 => "EEG", 2 => "est", 3 => "res");
    dspec_dict = Dict([i => x for (i, x) in enumerate(models.Sets)]);
    dspec_dict[42] = [:EEG];
    out_data[!,:Type] = [dtype_dict[x] for x in out_data[:,:Type]];
    out_data[!,:Spec] = [dspec_dict[x] for x in out_data[:,:Spec]];

    if file != "none"
        write(string("../data/", file, "_data.csv"), out_data)
    end

    out_data
end
