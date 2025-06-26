## Original code:
# Christoph Aurnhammer <aurnhammer@coli.uni-saarland.de>
# EEG plotting options for (lme)rERPs
# Including my adaption of topography plotting from craddm
#
## Modified by:
# Benedict Krieger <bkrieger@lst.uni-saarland.de>
##

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
    study_id,
    surp_id
) {

    data_path = glue("../data/{study_id}/{file}")
    plots_path = glue("../results/{study_id}/plots/{file}")
    
    # make dirs
    system(glue("mkdir -p {plots_path}"))
    system(glue("mkdir -p {plots_path}/Waveforms"))
    #system(glue("mkdir -p {plots_path}/Topos"))

    ##################
    # Study-specific #
    ##################
    if (study_id == 'adsbc21') {
        time_windows <- list(c(300, 500), c(600, 1000))
        #data_labs <- c("A: A+E+", "B: A-E+", "C: A+E-", "D: A-E-")
        data_labs <- c("A: Assoc+Exp+", "B: Assoc-Exp+", "C: Assoc+Exp-", "D: Assoc-Exp-")
        data_vals <- c("#000000", "#BB5566", "#004488", "#DDAA33")
        observed_title <- "Observed"
        topo = FALSE
        }
    else if (study_id == 'dbc19') {
        time_windows <- list(c(300, 500), c(600, 1000))
        #data_labs <- c("A: Baseline",
        #            "B: Event related violation",
        #            "C: Event unrelated violation")
        data_labs <- c("Assoc+Exp+",
                    "Assoc+Exp-",
                    "Assoc-Exp-")
        data_vals <- c("#000000", "red", "blue")
        observed_title <- "Observed"
        topo = FALSE
        }
    else if (study_id == 'dbc19_corrected') {
        time_windows <- list(c(300, 500), c(600, 1000))
        #data_labs <- c("A: Baseline",
        #            "B: Event related violation",
        #            "C: Event unrelated violation")
        data_labs <- c("Assoc+Exp+",
                    "Assoc+Exp-",
                    "Assoc-Exp-")
        data_vals <- c("#000000", "red", "blue")
        observed_title <- "Re-estimated"
        topo = FALSE
        }
    else if (study_id == 'adbc23') {
            time_windows <- list(c(300, 500), c(600, 1000))
            #data_labs <- c("A: Plausible",
            #        "B: Less plausible, attraction",
            #        "C: Implausible, no attraction")
            data_labs <- c("Exp+",
                    "Exp-",
                    "Exp--")
            data_vals <- c("#000000", "red", "blue")
        observed_title <- "Observed"
            topo = FALSE 
        }
    

    if (grepl("across", file)) {
        ci = FALSE
    } else {
        ci = TRUE
    }

    ##########
    # MODELS #
    ##########
    mod <- fread(glue("{data_path}_models.csv"))
    mod$Spec <- factor(mod$Spec, levels = predictor)

    # Models: coefficent
    coef <- mod[Type == "Coefficient", ]
    coef$Condition <- coef$Spec
    model_vals <- c("black", "#E349F6", "#00FFFF")

    plot_single_elec(
        data = coef,
        e = c("Pz"),
        file = glue("{plots_path}/Waveforms/Coefficients.pdf"),
        title = "Coefficients",
        modus = "Coefficient",
        ylims = c(10.5, -7),
        leg_labs = model_labs,
        leg_vals = model_vals,
        omit_legend = TRUE,
        ci = TRUE
        )

    plot_full_elec(
        data = coef,
        e = elec,
        file = glue("{plots_path}/Waveforms/Coefficients_Full.pdf"),
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
        #elec_corr <- c("Pz")
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
        #plot_nine_elec(
        plot_single_elec(
            data = tval,
            #e = elec_corr,
            e = "Pz",
            file = glue("{plots_path}/Waveforms/t-values.pdf"),
            title = glue("Inferential Statistics ",tail(model_labs,n=1)),
            modus = "t-value",
            ylims = c(8, -9),
            omit_legend = TRUE,
            tws = time_windows,
            leg_labs = model_labs[2:length(model_labs)],
            leg_vals = model_vals[2:length(model_vals)])
    }

    ########
    # DATA #
    ########
    eeg <- fread(glue("{data_path}_data.csv"))
    if (study_id == 'dbc19' | study_id == 'dbc19_corrected' | study_id == 'adbc23') {
    eeg$Condition <- factor(plyr::mapvalues(eeg$Condition, c(2, 1, 3),
                        c("B", "A", "C")), levels = c("A", "B", "C"))
    } else if (study_id == 'adsbc21') {
    eeg$Condition <- factor(plyr::mapvalues(eeg$Condition, c(2, 1, 3, 4),
                        c("B", "A", "C", "D")), levels = c("A", "B", "C", "D"))
    }
    # data_labs <- data_labs
    # data_vals <- data_vals
    #print(head(eeg$Condition))

    # Data: Observed
    obs <- eeg[Type == "EEG", ]

    plot_single_elec(
        data = obs,
        e = c("Pz"),
        file = glue("{plots_path}/Waveforms/Observed.pdf"),
        title = observed_title,   #"Observed"
        ylims = c(10.5, -7),
        modus = "Condition",
        tws = time_windows,
        ci = TRUE,
        leg_labs = data_labs,
        leg_vals = data_vals,
        omit_legend = TRUE,
        save_legend = TRUE,
        omit_x = FALSE,
        omit_y = FALSE,
        annotate = FALSE,
        highlight_time_windows = TRUE)

    plot_full_elec(
        data = obs,
        e = elec_all,
        file = glue("{plots_path}/Waveforms/Observed_Full.pdf"),
        title = "Observed",
        modus = "Condition",
        ci = TRUE,
        ylims = c(10.5, -7),
        leg_labs = data_labs,
        leg_vals = data_vals)

    # if (topo == TRUE) {
    #     plot_topo(
    #     data = obs,
    #     file = glue("{plots_path}/Topos/Observed"),
    #     tw = c(600, 1000),
    #     cond_man = "B",
    #     cond_base = "A",
    #     omit_legend = TRUE,
    #     save_legend = TRUE)
    # }

    # Adding plots for midline
    plot_midline(
        data = obs,
        e = c("Fz","Cz","Pz"),
        file = glue("{plots_path}/Waveforms/ObservedMidline.pdf"),
        title = "Observed",
        modus = "Condition",
        ci = ci,
        ylims = c(10.5, -7),
        leg_labs = data_labs,
        leg_vals = data_vals)

    # Data: Estimated
    # combo <- c("Intercept", "Intercept + Leo-13b surprisal",
    #             "Intercept + Noun Association",
    #             "Intercept + Cloze + Noun Association")

    surp_labs = c("leo13b_surp" = "Leo-13b surprisal", "gerpt2_surp" = "GerPT-2 surprisal", "gerpt2large_surp" = "GerPT-2 large surprisal")
    s_lab = surp_labs[surp_id]
    combo <- c("Intercept", glue("Intercept + {s_lab}"))

    # pred <- c("Intercept", "Leo-13b surprisal")
    pred <- c("Intercept", glue("{s_lab}"))
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
            file = glue("{plots_path}/Waveforms/Estimated_{name}.pdf"),
            title = paste("Estimates", combo[i]),
            modus = "Condition",
            tws = time_windows, # overwrite default tws with study-specific ones
            ci = TRUE,
            ylims = c(10.5, -7),
            leg_labs = data_labs,
            leg_vals = data_vals,
            omit_x = TRUE,
            omit_y = FALSE,
            omit_legend = TRUE,
            annotate = FALSE,
            highlight_time_windows = TRUE)
    #     plot_topo(
    #         data = est_set,
    #         file = paste0("../plots/", file, "/Topos/Estimated_", name),
    #         tw = c(600, 1000),
    #         cond_man = "B",
    #         cond_base = "A",
    #         add_title = paste("\nEstimate", pred[i]),
    #         omit_legend = TRUE)
    
        # Adding plots for midline
        plot_midline(
            data = est_set,
            e = c("Fz","Cz","Pz"),
            file = glue("{plots_path}/Waveforms/EstimatedMidline_{name}.pdf"),
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
            file = glue("{plots_path}/Waveforms/Residual_{name}.pdf"),
            title = paste("Residuals", combo[i]),
            modus = "Condition",
            tws = time_windows,
            ci = TRUE,
            omit_legend = TRUE,
            omit_x = TRUE,
            ylims = c(6, -6),
            leg_labs = data_labs,
            leg_vals = data_vals,
            highlight_time_windows = TRUE)
    #     plot_topo(
    #         data = res_set,
    #         file = paste0("../plots/", file, "/Topos/Residual_", name),
    #         tw = c(600, 1000),
    #         cond_man = "B",
    #         cond_base = "A",
    #         add_title = paste("\nResidual", pred[i]),
    #         omit_legend = TRUE)
    }
}

######################################################################################################################################################################################
######################################################################################################################################################################################

elec_all <- c("Fp1", "Fp2", "F7", "F3", "Fz", "F4", "F8", "FC5",
                "FC1", "FC2", "FC6", "C3", "Cz", "C4", "CP5", "CP1",
                "CP2", "CP6", "P7", "P3", "Pz", "P4", "P8", "O1", "Oz", "O2")


study_ids = list("adsbc21", "dbc19", "adbc23", "dbc19_corrected")
surp_ids = list("leo13b_surp", "gerpt2_surp", "gerpt2large_surp")
infer_options = list(TRUE, FALSE)
surp_labs = c("leo13b_surp" = "Leo-13b surprisal", "gerpt2_surp" = "GerPT-2 surprisal", "gerpt2large_surp" = "GerPT-2 large surprisal")

for (o in infer_options) {
    for (st in study_ids) {
        for (su in surp_ids) {
            file = if (o) glue("{st}_{su}_across_subj_rERP") else glue("{st}_{su}_rERP")
            print(file)
            s_lab = surp_labs[su]
            print(s_lab)
            make_plots(file, elec_all, predictor = c("Intercept", su), model_labs = c("Intercept", s_lab), inferential = o, study_id = st, surp_id = su)
        }
    }
}

# make_plots("dbc19_leo13b_surp_rERP", elec_all, predictor = c("Intercept", "leo13b_surp"), model_labs = c("Intercept", "Leo-13b surprisal"), inferential = FALSE, study_id = "dbc19", surp_id = "leo13b_surp")
#make_plots("adsbc21_leo13b_surp_rERP", elec_all, predictor = c("Intercept", "leo13b_surp"), model_labs = c("Intercept", "Leo-13b surprisal"), inferential = TRUE, study_id = "adsbc21", surp_id = "leo13b_surp")
#make_plots("adbc23_leo13b_surp_rERP", elec_all, predictor = c("Intercept", "leo13b_surp"), model_labs = c("Intercept", "Leo-13b surprisal"), inferential = TRUE, study_id = "adbc23", surp_id = "leo13b_surp")
