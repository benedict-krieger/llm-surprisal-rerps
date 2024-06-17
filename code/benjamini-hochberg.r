# Benjamini-Hochberg procedure
# Both functions apply correction to all electrodes they receive.
# If correction should be applied to a subset of electrodes
# this subset has the be taken before the function call

# Apply bh correction to data in wide format (electrodes in columns)
bh_apply_wide <- function(
        data,
        elec_corr,
        alpha = 0.05,
        tws = list(c(300, 500), c(600, 1000))
) {
        preds <- unique(data$Spec)
        end_elec <- 4 + length(elec_corr) - 1
        elec_indices <- c(4:end_elec)
        keep_ts <- c()
        for (tw in tws) {
                keep_ts <- c(keep_ts, seq(tw[1], tw[2]))
                for (p in preds) {
                        # Select current predictor, and timewindow
                        df <- data[Spec == p & Timestamp >= tw[1] &
                                        Timestamp <= tw[2], ..elec_indices]
                        uncorrected <- unlist(df)
                        # Apply correction to vector containing
                        # all current p-values
                        corrected <- p.adjust(uncorrected, method = 'fdr')
                        corrected_matrix <- matrix(corrected,
                                                nrow = nrow(df),
                                                ncol = length(elec_corr))

                        # Replace uncorrected with corrected pvalues
                        data[Spec == p & Timestamp >= tw[1] &
                                Timestamp <= tw[2], elec_indices] <- data.table(
                                                        corrected_matrix)
                        # Store binary significance
                        data[Spec == p & Timestamp >= tw[1] &
                            Timestamp <= tw[2],
                            paste0(colnames(data)[elec_indices],
                                "_sig")] <- data.table(corrected_matrix < alpha)
                }
        }
        # Set everything outside of time-windows to nonsignificant
        sigcols <- grep("_sig", colnames(data))
        data[!(Timestamp %in% keep_ts), c(elec_indices, sigcols)] <- 0
        data
}

# Apply bh correction to data in long format (electrodes in rows)
bh_apply <- function(
        data,
        alpha = 0.05,
        time_windows = list(c(300, 500), c(600, 1000))
) {
        tws <- time_windows
        dt <- data[, lapply(.SD, mean), by = list(Timestamp, Electrode),
                .SDcols = colnames(data)[grep("pval", colnames(data))]]
        for (i in colnames(dt)[grep("pval", colnames(dt))]) {
                dt[, paste0("sig_", i)] <- rep(FALSE, nrow(dt))
                for (j in tws) {
                        signi <- (p.adjust(dt[Timestamp >= j[1] &
                                Timestamp <= j[2], ..i][[1]],
                                method = 'fdr') < alpha)
                        indi <- which(dt$Timestamp >= j[1] &
                                dt$Timestamp <= j[2])
                        dt[indi, paste0("sig_", i) := signi]
                }
        }

        cols <- c("Timestamp", "Electrode",
                colnames(dt)[grep("sig", colnames(dt))])
        merge(data, dt[,..cols], by = c("Timestamp", "Electrode"))
}