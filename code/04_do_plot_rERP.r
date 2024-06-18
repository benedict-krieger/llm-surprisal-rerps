# This file originates from adbc23
library(here)
library(glue)
source("../code/plot_rERP_v2.r")
source("../code/benjamini-hochberg.r")

make_plots <- function(
    file,
    elec = c("F3", "Fz", "F4", "C3", "Cz", "C4", "P3", "Pz", "P4"),
    predictor = "Intercept",
    inferential = FALSE,
    model_labs,
    study_id
) {
    # make dirs
    plots_dir = glue("../results/{study_id}/plots/{file}/")
    system(paste0("mkdir -p ", plots_dir))
    system(paste0("mkdir -p ", plots_dir, "/Waveforms"))
    system(paste0("mkdir -p ", plots_dir, "/Topos"))

    ##################
    # Study-specific #
    ##################

    if (study_id == 'adsbc21') {
        time_windows <- list(c(350, 450), c(600, 800))
        data_labs <- c("A: A+E+", "B: A-E+", "C: A+E-", "D: A-E-")
        data_vals <- c("#000000", "#BB5566", "#004488", "#DDAA33")
        }
    else if (study_id == 'dbc19') {
        time_windows <- list(c(300, 500), c(800, 1000))
        data_labs <- c("A: Baseline",
                    "B: Event related violation",
                    "C: Event unrelated violation")
        data_vals <- c("#000000", "red", "blue")
        }
    else if (study_id == 'adbc23') {
            time_windows <- list(c(300, 500), c(600, 1000))
            data_labs <- c("A: Plausible",
                    "B: Less plausible, attraction",
                    "C: Implausible, no attraction")
            data_vals <- c("#000000", "red", "blue") 
        }
    
    
    if (grepl("across", file)) {
        ci = FALSE
    } else {
        ci = TRUE
    }

    ##########
    # MODELS #
    ##########
    mod <- fread(paste0("../data/",study_id, "/" , file, "_models.csv"))
    mod$Spec <- factor(mod$Spec, levels = predictor)

    # Models: coefficent
    coef <- mod[Type == "Coefficient", ]
    coef$Condition <- coef$Spec
    model_vals <- c("black", "#E349F6", "#00FFFF")

    plot_single_elec(
        data = coef,
        e = c("Pz"),
        file = paste0("../plots/", file, "/Waveforms/Coefficients.pdf"),
        title = "Coefficients",
        modus = "Coefficient",
        ylims = c(10.5, -7),
        leg_labs = model_labs,
        leg_vals = model_vals)

    plot_full_elec(
        data = coef,
        e = elec,
        file = paste0("../plots/",
        file,
        "/Waveforms/Coefficients_Full.pdf"),
        title = "Coefficients",
        modus = "Coefficient",
        ylims = c(7, -5),
        leg_labs = model_labs,
        leg_vals = model_vals)

    # Models: t-value
    if (inferential == TRUE) {
        # Specify subsets of time-windows and 
        # electrodes within which to correct
        time_windows <- time_windows
        elec_corr <- c("F3", "Fz", "F4", "C3", "Cz", "C4", "P3", "Pz", "P4")
        cols <- c("Timestamp", "Type", "Spec",
                elec_corr, paste0(elec_corr, "_CI"))
        mod <- mod[, ..cols]
        tval <- mod[Type == "t-value" & Spec != "Intercept", ]
        sig <- mod[Type == "p-value" & Spec != "Intercept", ]
        colnames(sig) <- gsub("_CI", "_sig", colnames(sig))

        # Apply correction
        sig_corr <- bh_apply_wide(
                        sig,
                        elec_corr,
                        alpha = 0.05,
                        tws = time_windows)
        sigcols <- grepl("_sig", colnames(sig_corr))
        tval <- cbind(tval, sig_corr[, ..sigcols])
        tval$Condition <- tval$Spec
        plot_nine_elec(
            data = tval,
            e = elec_corr,
            file = paste0("../plots/", file, "/Waveforms/t-values.pdf"),
            title = "Inferential Statistics",
            modus = "t-value",
            ylims = c(8, -9),
            tws = time_windows,
            leg_labs = model_labs[2:length(model_labs)],
            leg_vals = model_vals[2:length(model_vals)])
    }

    ########
    # DATA #
    ########
    eeg <- fread(paste0("../data/", study_id, file, "_data.csv"))
    eeg$Condition <- factor(plyr::mapvalues(eeg$Condition, c(2, 1, 3),
                        c("B", "A", "C")), levels = c("A", "B", "C"))
    data_labs <- data_labs
    data_vals <- data_vals

    # Data: Observed
    obs <- eeg[Type == "EEG", ]

    plot_single_elec(
        data = obs,
        e = c("Pz"),
        file = paste0("../plots/", file,  "/Waveforms/Observed.pdf"),
        title = "Observed",
        ylims = c(10.5, -7),
        modus = "Condition",
        ci = ci,
        leg_labs = data_labs,
        leg_vals = data_vals,
        omit_legend = TRUE,
        save_legend = TRUE,
        omit_x = TRUE,
        omit_y = TRUE,
        annotate = TRUE,
        text = "7")

    plot_full_elec(
        data = obs,
        e = elec_all,
        file = paste0("../plots/", file, "/Waveforms/Observed_Full.pdf"),
        title = "Observed",
        modus = "Condition",
        ci = ci,
        ylims = c(10.5, -7),
        leg_labs = data_labs,
        leg_vals = data_vals)

    plot_topo(
        data = obs,
        file = paste0("../plots/", file, "/Topos/Observed"),
        tw = c(600, 1000),
        cond_man = "B",
        cond_base = "A",
        omit_legend = TRUE,
        save_legend = TRUE)
    
    # Adding plots for midline
    plot_midline(
        data = obs,
        e = c("Fz","Cz","Pz"),
        file = paste0("../plots/", file, "/Waveforms/ObservedMidline.pdf"),
        title = "Observed",
        modus = "Condition",
        ci = ci,
        ylims = c(10.5, -7),
        leg_labs = data_labs,
        leg_vals = data_vals)

    # Data: Estimated
    combo <- c("Intercept", "Intercept + Leo-13b surprisal",
                "Intercept + Noun Association",
                "Intercept + Cloze + Noun Association")

    pred <- c("Intercept", "Leo-13b surprisal")
    est <- eeg[Type == "est", ]
    for (i in seq(1, length(unique(est$Spec)))) {
        spec <- unique(est$Spec)[i]
        print(spec)
        est_set <- est[Spec == spec, ]
        spec <- unique(est_set$Spec)
        name <- gsub("\\[|\\]|:|,| ", "", spec)
        plot_single_elec(
            data = est_set,
            e = c("Pz"),
            file = paste0("../plots/", file, 
                    "/Waveforms/Estimated_", name, ".pdf"),
            title = paste("Estimates", combo[i]),
            modus = "Condition",
            ci = ci,
            ylims = c(10.5, -7),
            leg_labs = data_labs,
            leg_vals = data_vals,
            omit_x = TRUE,
            omit_y=TRUE,
            omit_legend = TRUE,
            annotate = TRUE,
            text = "8")
        plot_topo(
            data = est_set,
            file = paste0("../plots/", file, "/Topos/Estimated_", name),
            tw = c(600, 1000),
            cond_man = "B",
            cond_base = "A",
            add_title = paste("\nEstimate", pred[i]),
            omit_legend = TRUE)
    
        # Adding plots for midline
        plot_midline(
            data = est_set,
            e = c("Fz","Cz","Pz"),
            file = paste0("../plots/", file, 
                    "/Waveforms/EstimatedMidline_", name, ".pdf"),
            title = paste("Estimates", combo[i]),
            modus = "Condition",
            ci = ci,
            ylims = c(10.5, -7),
            leg_labs = data_labs,
            leg_vals = data_vals)


    }
    
    # Data: Residual
    res <- eeg[Type == "res", ]
    for (i in seq(1, length(unique(res$Spec)))) {
        spec <- unique(res$Spec)[i]
        res_set <- res[Spec == spec, ]
        spec <- unique(res_set$Spec)
        name <- gsub("\\[|\\]|:|,| ", "", spec)
        plot_single_elec(
            data = res_set,
            e = c("Pz"),
            file = paste0("../plots/", file,
                    "/Waveforms/Residual_", name, ".pdf"),
            title = paste("Residuals", combo[i]),
            modus = "Condition",
            ci = ci,
            ylims = c(6, -6),
            leg_labs = data_labs,
            leg_vals = data_vals)
        plot_topo(
            data = res_set,
            file = paste0("../plots/", file, "/Topos/Residual_", name),
            tw = c(600, 1000),
            cond_man = "B",
            cond_base = "A",
            add_title = paste("\nResidual", pred[i]),
            omit_legend = TRUE)
    }
}

elec_all <- c("Fp1", "Fp2", "F7", "F3", "Fz", "F4", "F8", "FC5",
                "FC1", "FC2", "FC6", "C3", "Cz", "C4", "CP5", "CP1",
                "CP2", "CP6", "P7", "P3", "Pz", "P4", "P8", "O1", "Oz", "O2")

#make_plots("adbc23_rERP_leo13b", elec_all,
#    predictor = c("Intercept", "leo13b_s"), model_labs = c("Intercept", "Leo-13b surprisal"))

#make_plots("adbc23_rERP_leo13b_across_s", elec_all,
#    predictor = c("Intercept", "leo13b_s"), model_labs = c("Intercept", "Leo-13b surprisal"), inferential = TRUE)

study_ids = list("adsbc21", "dbc19", "adbc23")
surp_ids = list("leo13b_surp", "secretgpt2_surp")
infer_options = list(TRUE, FALSE)
surp_labs = c("leo13b_surp" = "Leo-13b surprisal", "secretgpt2_surp" = "GPT-2 surprisal")


for (o in infer_options) {
    for (st in study_ids) {
        for (su in surp_ids) {
            file = if (o) paste(st, su, "across_subj", "rERP", sep='_') else paste(st, su, "rERP", sep='_')
            print(file)
            s_lab = surp_labs[su]
            print(s_lab)
            make_plots(file, elec_all, predictor = c("Intercept", su), model_labs c("Intercept", s_lab), inferential = o, study_id = st)
        }
    }
}