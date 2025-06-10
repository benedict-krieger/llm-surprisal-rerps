## Original code:
# Christoph Aurnhammer <aurnhammer@coli.uni-saarland.de>
# EEG plotting options for (lme)rERPs
# Including my adaption of topography plotting from craddm
#
## Modified by:
# Benedict Krieger <bkrieger@lst.uni-saarland.de>
##

library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(stringr) # Added for string wrap in legends

# compute standard error
se <- function(
    x,
    na_rm = FALSE
) {
    if (na_rm == TRUE) {
        sd(x, na.rm = TRUE) / sqrt(length(x[!is.na(x)]))
    } else {
        sd(x) / sqrt(length(x))
    }
}

# Return only the legend of an ggplot object
get_legend <- function(
    a_gplot
) { 
    tmp <- ggplot_gtable(ggplot_build(a_gplot +
            theme(legend.box = "vertical",
                  legend.spacing.y = unit(0.005, "inch"))))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}

# For a single Electrode, plot per-grouping mean.
plot_grandavg_ci <- function(
    df,
    ttl,
    yunit = "Amplitude (µVolt)",
    ylims = NULL,
    modus = "Condition",
    tws = list(c(300, 500), c(600, 1000)),
    ci = TRUE,
    leg_labs,
    leg_vals
) {
    ###### DATA PROC
    if (modus %in% c("Condition", "Quantile")) {
        colnames(df)[c(4, 5)] <- c("V", "V_CI")
    } else if (modus == "Coefficient") {
        colnames(df)[c(3, 4)] <- c("V", "V_CI")
    } else if (modus == "t-value") {
        colnames(df)[c(3, 4, 5)] <- c("V", "V_CI", "V_sig")
        sig_dt <- df[, c("Spec", "Timestamp", "V_sig")]
        sig_dt$posit <- rep(seq(ylims[1], ylims[1] - 2,
            length = length(unique(sig_dt$Spec))),
            length = nrow(sig_dt))
        sig_dt$sig <- factor(sig_dt$V_sig, levels = c("0", "1"),
            labels = c("insign", "sign"))
    }
    ##### PLOTTING
    # Initial plot call
    if (modus %in% c("Coefficient", "t-value")) {
        p <- ggplot(df, aes(x = Timestamp, y = V,
            color = Spec, fill = Spec)) + geom_line()
    }
    else if (modus %in% c("Condition", "Quantile")) {
        p <- ggplot(df, aes(x = Timestamp, y = V,
            color = df[[modus]], fill = df[[modus]])) + geom_line()
    }

    # For all plots
    if (ci == TRUE) {
    p <- p + geom_ribbon(aes(x = Timestamp,
        ymax = V + V_CI, ymin = V - V_CI), alpha = 0.15, color = NA)
    }
    p <- p + geom_hline(yintercept = 0, linetype = "dashed")
    p <- p + geom_vline(xintercept = 0, linetype = "dashed")
    p <- p + theme(panel.background = element_rect(fill = "#FFFFFF",
                    color = "#000000", linewidth = 0.1, linetype = "solid"),
                panel.grid.major = element_line(linewidth = 0.3,
                    linetype = "solid", color = "#A9A9A9"),
                panel.grid.minor = element_line(linewidth = 0.15,
                    linetype = "solid", color = "#A9A9A9"),
            legend.position = "top")

    # Conditional modifications
    if (is.vector(ylims) == TRUE) {
        p <- p + scale_y_reverse(limits = c(ylims[1], ylims[2]))
    } else {
        p <- p + scale_y_reverse()
    }
    if (modus == "Quantile") {
        p <- p + labs(y = yunit, x = "Time (ms)", title = ttl)
        p <- p + scale_color_manual(
                name = "N400 Quantile",
                labels = leg_labs,
                values = leg_vals)
        p <- p + scale_fill_manual(
                name = "N400 Quantile",
                labels = leg_labs,
                values = leg_vals)
    } else if (modus == "Condition") {
        p <- p + labs(y = yunit, x = "Time (ms)", title = ttl)
        p <- p + scale_color_manual(name = "Condition",
                labels = leg_labs,
                values = leg_vals)
        p <- p + scale_fill_manual(
                name = "Condition",
                labels = leg_labs,
                values = leg_vals)
    } else if (modus == "Coefficient") {
        p <- p + labs(  y = "Intercept + Coefficient", 
                        x = "Time (ms)",
                        title = ttl)
        p <- p + scale_color_manual(name = modus,
                labels = leg_labs, values = leg_vals)
        p <- p + scale_fill_manual(name = modus,
                labels = leg_labs, values = leg_vals)
    } else if (modus == "t-value") {
        p <- p + labs(y = "T-value", x = "Time (ms)", title = ttl)
        p <- p + scale_color_manual(name = "Predictor",
                labels = leg_labs, values = leg_vals)
        p <- p + scale_fill_manual(name = "Predictor",
                labels = leg_labs, values = leg_vals)
        p <- p + geom_point(data = sig_dt,
                    aes(x = Timestamp, y = posit, shape = sig),
                    size = 1)
        p <- p + scale_shape_manual(values = c(32, 124),
                    name = "Corrected p-values",
                    labels = c("Nonsignificant", "Significant"))
        p <- p + annotate("rect", xmin = tws[1][[1]][1], xmax = tws[1][[1]][2],
                    ymin = ylims[1], ymax = -Inf, alpha = .15)
        p <- p + annotate("rect", xmin = tws[2][[1]][1], xmax = tws[2][[1]][2],
                    ymin = ylims[1], ymax = -Inf, alpha = .15)
    }

    p
}

plot_single_elec <- function(
    data,
    e,
    file = FALSE,
    title = "ERPs",
    yunit = "Amplitude (µVolt)",
    ylims = NULL,
    modus = "Condition",
    tws = list(c(250, 400), c(600, 1000)),
    ci = TRUE,
    leg_labs,
    leg_vals,
    omit_legend = FALSE,
    save_legend = FALSE,
    omit_x = FALSE,
    omit_y = FALSE,
    annotate = FALSE,
    text = "",
    highlight_background = FALSE,
    highlight_time_windows = FALSE
) { 
    if (modus %in% c("Quantile", "Condition")) {
        cols <- c("Spec", "Timestamp", modus)
    } else if (modus %in% c("Coefficient", "t-value")) {
        data[,"Spec"] <- as.factor(data$Spec)
        cols <- c("Spec", "Timestamp")
    }

    # Make individual plots
    plotlist <- vector(mode = "list", length = length(e))
    if (modus == "t-value") {
        varforward <- c(e, paste0(e, "_CI"), paste0(e, "_sig"))
        p <- plot_grandavg_ci(
                df = cbind(data[, ..cols], data[, ..varforward]),
                ttl = e,
                yunit = yunit,
                ylims = ylims,
                modus = modus,
                tws = tws,
                ci = FALSE,
                leg_labs = leg_labs,
                leg_vals = leg_vals)
    } else if (modus %in% c("Coefficient")) {
        varforward <- c(e, paste0(e, "_CI"))
        p <- plot_grandavg_ci(
                df = cbind(data[, ..cols], data[, ..varforward]),
                ttl = e,
                yunit = yunit,
                ylims = ylims,
                modus = modus,
                ci = ci,
                leg_labs = leg_labs,
                leg_vals = leg_vals)
    } else if (modus %in% c("Quantile", "Condition")) {
        varforward <- c(e, paste0(e, "_CI"))
        p <- plot_grandavg_ci(
                df = cbind(data[, ..cols], data[, ..varforward]),
                ttl = e,
                yunit = yunit,
                ylims = ylims,
                modus = modus,
                ci = ci,
                leg_labs = str_wrap(leg_labs,20),
                leg_vals = leg_vals)
    }

    if (file != FALSE) {
        gg <- p
        gg <- gg + theme(legend.key.size = unit(0.5, 'cm'), # lgnd key size
            legend.key.height = unit(0.5, 'cm'),            # lgnd key height
            legend.key.width = unit(0.5, 'cm'),             # lgnd key width
            legend.title = element_text(size = 7),          # lgnd ttl font size; original: 7, prev: 11
            legend.text = element_text(size = 6))           # lgnd el font size; original: 6, prev: 11
        gg <- gg + theme(plot.title = element_text(size = 11.5)) # original: 7.5, prev: 11.5
        gg <- gg + theme(axis.text.x = element_text(size= 11), # prev: 14
                         axis.text.y = element_text(size = 11)) # prev 14

        if (omit_x != FALSE) {
            gg <- gg + theme(axis.title.x = element_blank(),
                             axis.text.x = element_blank(),
                             axis.ticks.x = element_blank() )
        }

        if (omit_y != FALSE) {
            gg <- gg + theme(axis.title.y = element_blank(),
                             axis.text.y = element_blank(),
                             axis.ticks.y = element_blank() )
        }

        if (annotate != FALSE) {
            text_label = grobTree(textGrob(text,x=0.01, y=0.95, hjust = 0, gp=gpar(fontsize=20, fontface="bold")))
            gg <- gg + annotation_custom(text_label) 
        }

        if (highlight_background != FALSE) {
            gg <- gg + theme(panel.background = element_rect(fill = '#f0f0f0')) # bit darker: #e6e6e6
        }

        if (highlight_time_windows != FALSE) {
            
        #gg <- gg + annotate("rect", xmin = tws[1][[1]][1], xmax = tws[1][[1]][2],
        #            ymin = ylims[1], ymax = ylims[2], alpha = .15)
        #gg <- gg + annotate("rect", xmin = tws[2][[1]][1], xmax = tws[2][[1]][2],
        #            ymin = ylims[1], ymax = ylims[2], alpha = .15)
        
        gg <- gg + annotate("rect", xmin = tws[1][[1]][1], xmax = tws[1][[1]][2],
                    ymin = -Inf, ymax = Inf, alpha = .15)
        gg <- gg + annotate("rect", xmin = tws[2][[1]][1], xmax = tws[2][[1]][2],
                    ymin = -Inf, ymax = Inf, alpha = .15)
        
        }        

        # Save legend to file
        if (save_legend) {
            lgnd <- get_legend(gg)
            file_trimmed <- strtrim(file, nchar(file) - 4)
            ggsave(paste0(file_trimmed, "_wavelegend.pdf"),
                    lgnd, device = cairo_pdf,
                    width = 5, height = 0.5) # original width=4.5, height=0.5
        }

        # Add legend or not
        if (omit_legend) {
            gg <- gg + theme(legend.position = "none")
            gg <- arrangeGrob(gg + ggtitle(paste0(e, ": ", title)),
                    heights = c(10, 0.25))
        } else {
            legend <- get_legend(gg)
            nl <- theme(legend.position = "none")
            gg <- arrangeGrob(gg + nl + ggtitle(paste0(e, ": ", title)),
            legend, heights = c(10, 2))
        }

        ggsave(file, gg, device = cairo_pdf, width = 4, height = 4)
    } else {
        p
    }
}


# Plot three midline electrodes
plot_midline <- function(
    data,
    e,
    file = FALSE,
    title = "ERPs",
    yunit = paste0("Amplitude (", "\u03BC", "Volt\u29"),
    ylims = NULL,
    modus = "Condition",
    tws = list(c(250, 400), c(600, 1000)),
    ci = TRUE,
    leg_labs,
    leg_vals
) {
    if (modus %in% c("Tertile", "Quantile", "Condition")) {
        cols <- c("Spec", "Timestamp", modus)
    } else if (modus %in% c("Coefficient", "t-value")) {
        data[,"Spec"] <- as.factor(data$Spec)
        cols <- c("Spec", "Timestamp")
    }

    # Make individual plots
    plotlist <- vector(mode = "list", length = length(e))
    if (modus == "t-value") {
        for (i in 1:length(e)) {
            varforward <- c(e[i], paste0(e[i], "_CI"), paste0(e[i], "_sig"))
            plotlist[[i]] <- plot_grandavg_ci(cbind(data[, ..cols],
                    data[, ..varforward]), e[i], yunit,
                    ylims, modus, tws, ci = FALSE,
                    leg_labs = leg_labs, leg_vals = leg_vals)
        }
    } else if (modus %in% c("Coefficient", "Tertile")) {
        for (i in 1:length(e)) {
            varforward <- c(e[i], paste0(e[i], "_CI"))
            plotlist[[i]] <- plot_grandavg_ci(cbind(data[, ..cols],
                    data[, ..varforward]), e[i], yunit = yunit,
                    ylims = ylims, modus = modus, ci = ci,
                    leg_labs = leg_labs, leg_vals = leg_vals)
        }
    } else if (modus %in% c("Tertile", "Quantile", "Condition")) {
        for (i in 1:length(e)) {
            varforward <- c(e[i], paste0(e[i], "_CI"))
            plotlist[[i]] <- plot_grandavg_ci(cbind(data[, ..cols],
                    data[, ..varforward]), e[i], yunit = yunit,
                    ylims = ylims, modus = modus, ci = ci,
                    leg_labs = leg_labs, leg_vals = leg_vals)
        }
    }
    # Get the legend
    legend <- get_legend(plotlist[[1]])
    # Arrange
    nl <- theme(legend.position = "none")
    gg <- arrangeGrob(arrangeGrob(
        plotlist[[1]] + nl,
        plotlist[[2]] + nl,
        plotlist[[3]] + nl,
        layout_matrix = matrix(1:3, ncol = 1, byrow = TRUE)), legend,
            heights = c(10, 1), top = textGrob(title))
    if (file != FALSE) {
       ggsave(file, gg, device = cairo_pdf, width = 4, height = 7)
    } else {
       gg
    }
}

# Plot nine electrode grid
plot_nine_elec <- function(
    data,
    e,
    file = FALSE,
    title = "ERPs",
    yunit = "Amplitude (µVolt)",
    ylims = NULL,
    modus = "Condition",
    tws = list(c(250, 400), c(600, 1000)),
    ci = TRUE,
    leg_labs,
    leg_vals
) {
    if (modus %in% c("Quantile", "Condition")) {
        cols <- c("Spec", "Timestamp", modus)
    } else if (modus %in% c("Coefficient", "t-value")) {
        data[,"Spec"] <- as.factor(data$Spec)
        cols <- c("Spec", "Timestamp")
    }

    # Make individual plots
    plotlist <- vector(mode = "list", length = length(e))
    if (modus == "t-value") {
        for (i in 1:length(e)) {
            varforward <- c(e[i], paste0(e[i], "_CI"), paste0(e[i], "_sig"))
            plotlist[[i]] <- plot_grandavg_ci(cbind(data[, ..cols],
                    data[, ..varforward]), e[i], yunit,
                    ylims, modus, tws, ci = FALSE,
                    leg_labs = leg_labs, leg_vals = leg_vals)
        }
    } else if (modus %in% c("Coefficient")) {
        for (i in 1:length(e)) {
            varforward <- c(e[i], paste0(e[i], "_CI"))
            plotlist[[i]] <- plot_grandavg_ci(cbind(data[, ..cols],
                    data[, ..varforward]), e[i], yunit = yunit,
                    ylims = ylims, modus = modus, ci = ci,
                    leg_labs = leg_labs, leg_vals = leg_vals)
        }
    } else if (modus %in% c("Quantile", "Condition")) {
        for (i in 1:length(e)) {
            varforward <- c(e[i], paste0(e[i], "_CI"))
            plotlist[[i]] <- plot_grandavg_ci(cbind(data[, ..cols],
                    data[, ..varforward]), e[i], yunit = yunit,
                    ylims = ylims, modus = modus, ci = ci,
                    leg_labs = leg_labs, leg_vals = leg_vals)
        }
    }
    # Get the legend
    legend <- get_legend(plotlist[[1]])
    # Arrange
    nl <- theme(legend.position = "none")
    gg <- arrangeGrob(arrangeGrob(
        plotlist[[1]] + nl
            + labs(x = ""),
        plotlist[[2]] + nl
            + labs(x = "", y = ""),
        plotlist[[3]] + nl
            + labs(x = "", y = ""),
        plotlist[[4]] + nl
            + labs(x = ""),
        plotlist[[5]] + nl
            + labs(x = "", y = ""),
        plotlist[[6]] + nl
            + labs(x = "", y = ""),
        plotlist[[7]] + nl,
        plotlist[[8]] + nl
            + labs(y = ""),
        plotlist[[9]] + nl
            + labs(y = ""),
        layout_matrix = matrix(1:9, ncol = 3, byrow = TRUE)), legend,
            heights = c(10, 1), top = ggtitle(title))
    if (file != FALSE) {
       ggsave(file, gg, device = cairo_pdf, width = 7, height = 7)
    } else {
       gg
    }
}

plot_full_elec <- function(
    data,
    e,
    file = FALSE,
    title = "ERPs",
    yunit = "Amplitude (µVolt)",
    ylims = NULL,
    modus = "Condition",
    tws = list(c(300, 500), c(600, 1000)),
    ci = FALSE,
    leg_labs,
    leg_vals
) {
    if (modus %in% c("Quantile", "Condition")) {
        cols <- c("Spec", "Timestamp", modus)
    } else if (modus %in% c("Coefficient", "t-value")) {
        data[,"Spec"] <- as.factor(data$Spec)
        cols <- c("Spec", "Timestamp")
    }

    # Make individual plots
    plotlist <- vector(mode = "list", length = length(e))
    if (modus == "t-value") {
        for (i in 1:length(e)) {
            varforward <- c(e[i], paste0(e[i], "_CI"), paste0(e[i], "_sig"))
            plotlist[[i]] <- plot_grandavg_ci(cbind(data[, ..cols],
                            data[, ..varforward]), e[i], yunit,
                            ylims, modus, tws, ci = FALSE,
                            leg_labs = leg_labs, leg_vals = leg_vals)
        }
    } else if (modus %in% c("Coefficient")) {
        for (i in 1:length(e)) {
            varforward <- c(e[i], paste0(e[i], "_CI"))
            plotlist[[i]] <- plot_grandavg_ci(cbind(data[, ..cols],
                            data[, ..varforward]), e[i], yunit = yunit,
                            ylims = ylims, modus = modus, ci = ci,
                            leg_labs = leg_labs, leg_vals = leg_vals)
        }
    } else if (modus %in% c("Quantile", "Condition")) {
        for (i in 1:length(e)) {
            varforward <- c(e[i], paste0(e[i], "_CI"))
            plotlist[[i]] <- plot_grandavg_ci(cbind(data[, ..cols],
                            data[, ..varforward]), e[i], yunit = yunit,
                            ylims = ylims, modus = modus, ci = ci,
                            leg_labs = leg_labs, leg_vals = leg_vals)
        }
    }

    no <- ggplot(data.frame()) + theme_minimal()
    legend_theme <- theme(legend.key.size = unit(4, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.key.width = unit(1, 'cm'),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18))
    legend <- get_legend(plotlist[[1]] + legend_theme) # Get the legend

    nl <- theme(legend.position = "none") # No legend
    nla <- labs(x = "", y = "")
    ngrid <- theme(panel.background = element_rect(fill=NA, color = NA),
                   panel.grid.major.x = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   panel.grid.major.y = element_blank(),
                   panel.grid.minor.y = element_blank())
    axes <- theme(axis.line.x = element_line(color = "black", linewidth = 0.2),
                  axis.line.y = element_line(color = "black", linewidth = 0.2))
    for (i in 1:length(e)) {
        plotlist[[i]] <- plotlist[[i]] + nl + nla + ngrid + axes
    }

    gg <- arrangeGrob(arrangeGrob(
        no,
        plotlist[[1]], #Fp1
        no,
        plotlist[[2]], #Fp2
        no,
        plotlist[[3]], #F7
        plotlist[[4]], #F3
        plotlist[[5]], #Fz
        plotlist[[6]], #F4
        plotlist[[7]], #F8
        plotlist[[8]], #FC5
        plotlist[[9]], #FC1
        no,
        plotlist[[10]], #FC2
        plotlist[[11]], #FC6
        no,
        plotlist[[12]], #C3
        plotlist[[13]], #Cz
        plotlist[[14]], #C4
        no,
        plotlist[[15]], #CP5
        plotlist[[16]], #CP1
        no,
        plotlist[[17]], #CP2
        plotlist[[18]], #CP6
        plotlist[[19]], #P7
        plotlist[[20]], #P3
        plotlist[[21]], #Pz
        plotlist[[22]], #P4
        plotlist[[23]], #P8
        no,
        plotlist[[24]], #O1
        plotlist[[25]], #Oz
        plotlist[[26]], #O2
        no,
        layout_matrix = matrix(1:35, ncol = 5, byrow = TRUE)),
            legend,
            heights = c(15, 1),
            top = ggtitle(title) + theme(plot.title = element_text(size = 20))
        )

    if (file != FALSE) {
       ggsave(file, gg, device = cairo_pdf, width = 12, height = 12)
    } else {
       gg
    }
}

# Plot topographic map
plot_topo <- function(
    data,
    file,
    tw = c(600, 1000),
    cond_man,
    cond_base,
    label = "Amplitude (µVolt)",
    title = "",
    add_title = "",
    omit_legend = FALSE,
    save_legend = FALSE
) {
    # get elec vector
    first_elec_ind <- grep("Fp1$", colnames(data))
    electrodes <- colnames(data)[c(first_elec_ind:(ncol(data) - 26))]

    # Subset data
    CI_elecs <- grep("_CI", colnames(data))
    data <- data[, !..CI_elecs]
    data <- data[, !c("Type", "Spec")]

    # melt to long format
    data_m <- melt(data, id.vars = c("Timestamp", "Condition"),
                    variable.name = "Electrode", value.name = "EEG")

    # select time and descriptors
    data_m_tw <- data_m[(Timestamp >= tw[1] & Timestamp <= tw[2]), ]

    if (title == "") {
        ttl <- paste0(cond_man, " minus ", cond_base, add_title)
    } else {
        ttl <- title
    }

    generate_topo(
                data_m_tw,
                file,
                tw,
                cond_man,
                cond_base,
                amplim = 2.8,
                elec = electrodes,
                title = ttl,
                label = label,
                omit_legend,
                save_legend)
}

generate_topo <- function(
    data,
    file,
    tw,
    cond_man,
    cond_base,
    amplim,
    elec,
    title,
    label,
    omit_legend,
    save_legend
) {
    # Average, compute difference, add topo sites
    data_agg <- data[, lapply(.SD, mean), by = list(Electrode, Condition),
                        .SDcols = c("EEG")]
    if (!grepl("t-values", file)) {
        data_diff <- compute_difference(data_agg, cond_man, cond_base, title)
    } else {
        data_diff <- data_agg[Condition == cond_man,
                              c("Condition", "Electrode", "EEG")]
    }

    # Load electrode location data
    if (file.exists("../../data/biosemi70elecs.loc")) {
        elec_locs <- fread("../../data/biosemi70elecs.loc", sep = "\t",
                    col.names = c("ElecNum", "Theta", "Radius", "Electrode"))
    } else if (file.exists("../../../data/biosemi70elecs.loc")) {
        elec_locs <- fread("../../../data/biosemi70elecs.loc", sep = "\t",
                    col.names = c("ElecNum", "Theta", "Radius", "Electrode"))
    } else {
        elec_locs <- fread("../data/biosemi70elecs.loc", sep = "\t",
                    col.names = c("ElecNum", "Theta", "Radius", "Electrode"))
    }
    elec_locs <- elec_locs[Electrode %in% elec]
    elec_locs$RadianTheta <- pi / 180 * elec_locs$Theta
    elec_locs$x <- elec_locs$Radius * sin(elec_locs$RadianTheta)
    elec_locs$y <- elec_locs$Radius * cos(elec_locs$RadianTheta)
    data_diff_locs <- merge(elec_locs, data_diff, by = "Electrode")

    # Create data frame to be used for interpolation
    # - the function needs columns labelled x, y, and z
    interprep <- data.table(x = data_diff_locs$x,
                            y = data_diff_locs$y,
                            z = data_diff_locs$EEG)

    # create the interpolation grid
    rmax <- .75  # specify a maximum boundary for the grid
    grid_res <- 134 # specify the interpolation grid resolution
    xo <- seq(min(-rmax, interprep$x), max(rmax, interprep$x),
                length = grid_res)
    yo <- seq(max(rmax, interprep$y), min(-rmax, interprep$y),
                length = grid_res)

    # TODO clean up v4_interp
    interpol <- data.table(v4_interp(interprep, xo, yo, 
                            rmax = rmax, gridRes = grid_res))
    interpol_m <- melt(interpol, id.vars = "x", variable.name = "y",
                       value.name = "EEG")
    interpol_m$y <- as.numeric(as.character(interpol_m$y))

    # Select data in circle
    interpol_m$incircle <- (interpol_m$x)^2 + (interpol_m$y)^2 < 0.7 ^ 2
    interpol_m <- interpol_m[interpol_m$incircle, ]

    # TODO find out why there is an NA
    interpol_m <- interpol_m[!is.na(interpol_m$EEG), ]

    ###### Prepare v4 Plot
    amplim_plusmin <- c(-amplim, amplim) # compute limits

    # Define color spectrum
    # my_spectrum <- colorRampPalette(c("#00007F", "blue", "#0080ff", "white",
    #                         "white", "white", "#ff6969", "red", "#a90101"))
    my_spectrum <- colorRampPalette(c(
                "#01015e", "#3434ff", "#0080ff",
                "#7fbfff", "white", "#fcd4d4",
                "#ff6969", "#ca2c2c", "#750202"))
    # my_spectrum <- colorRampPalette(c(
    #                 "#202245", "#45428A", "#645FD1", "#8980FF",
    #                 "#A29EFF", "#C3C4FF", "#E4E7FF",
    #                 "white",
    #                 "#FBE4E6", "#FDBDC3", "#FF989F",
    #                 "#FF7681", "#CF5861", "#9E2626", "#651313"))
    # my_spectrum <- colorRampPalette(c(
    #                 "#132f4b", "#1c4269", "#4877a6", "#5285b9",
    #                 "#6ea3d8", "#8ebce9", "#9ccbfa",
    #                 "white",
    #                 "#FBE4E6", "#FDBDC3", "#FF989F",
    #                 "#FF7681", "#CF5861", "#9E2626", "#651313"
    # ))

    # Create head and nose shape
    head_shape <- circle_fun(c(0, 0), round(max(data_diff_locs$x) * 2.35, 4),
                             npoints = 100)
    nose <- data.table(x = c(-0.075, 0, .075), y = c(.495, .575, .495))

    # create a circle round the outside of the plotting area
    # to mask the jagged edges of the interpolation
    mask_ring <- circle_fun(diameter = 1.42)

    ###### Plot
    p <- ggplot(interpol_m, aes(x = x, y = y, fill = EEG))
    p <- p + geom_raster() + stat_contour(aes(z = EEG), binwidth = 0.5)
    p <- p + theme_topo()
    p <- p + labs(title = title, subtitle = paste0(tw[1], "-", tw[2], " ms"))
    p <- p + geom_path(data = mask_ring, aes(x, y, z = NULL, fill = NULL),
              colour = "white", linewidth = 6)
    p <- p + scale_fill_gradientn(colours = my_spectrum(10),
                    limits = amplim_plusmin, guide = "colourbar",
                    oob = scales::squish, name = label)
    # Add electrode positions
    p <- p + geom_point(data = data_diff_locs, aes(x, y), size = 1)
    p <- p + geom_path(data = head_shape, aes(x, y, z = NULL, fill = NULL),
                        linewidth = 1.5)
    p <- p + geom_path(data = nose, aes(x, y, z = NULL, fill = NULL),
                        linewidth = 1.5)
    p <- p + coord_equal()
    # Optional legend processing
    if (omit_legend) {
        if (save_legend) {
            suppressWarnings(lgnd <- get_legend(p))
            ggsave(paste0(file, "_topolegend.pdf"), lgnd, device = cairo_pdf,
                    width = 1.5, height = 4)
        }
        p <- p + theme(legend.position = "none")
    }

    suppressWarnings(ggsave(paste0(file, "_", cond_man, "_", tw[1], "-", tw[2], ".pdf"),
            p, device = cairo_pdf, width = 4, height = 4))
}

compute_difference <- function(
        data, cond_man, cond_base, name
) {
    data_diff <- data[Condition == cond_man, ]$EEG - data[Condition == cond_base, ]$EEG

    data.table(Timestamp = unique(data$Timestamp), Condition = factor(name),
             Electrode = unique(data$Electrode), EEG = data_diff)
}

circle_fun <- function(
    center = c(0, 0),
    diameter = 1,
    npoints = 100
) {
  r <- diameter / 2
  tt <- seq(0, 2 * pi, length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.table(x = xx, y = yy))
}

v4_interp <- function(df, xo, yo, rmax = .75, gridRes = gridRes) {
  xo <- matrix(rep(xo, length(yo)), nrow = length(xo), ncol = length(yo))
  yo <- t(matrix(rep(yo, length(xo)), nrow = length(yo), ncol = length(xo)))
  xy <- df$x + df$y * sqrt(as.complex(-1))
  d <- matrix(rep(xy, length(xy)), nrow = length(xy), ncol = length(xy))
  d <- abs(d - t(d))
  diag(d) <- 1
  g <- (d^2) * (log(d) - 1) # Green's function.
  diag(g) <- 0
  weights <- qr.solve(g, df$z)
  xy <- t(xy)
  outmat <- matrix(nrow = gridRes,ncol = gridRes)
  for (i in 1:gridRes) {
    for (j in 1:gridRes) {
      test4 <- abs((xo[i, j] + sqrt(as.complex(-1)) * yo[i, j]) - xy)
      g <- (test4^2) * (log(test4) - 1)
      outmat[i, j] <- g %*% weights
    }
  }
  outDf <- data.table(x = xo[, 1], outmat)
  names(outDf)[1:length(yo[1, ]) + 1] <- yo[1, ]

  outDf
}

theme_topo <- function(
    base_size = 12
) {
  theme_bw(base_size = base_size) %+replace%
    theme(rect = element_blank(),
          line = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank()
    )
}

items_and_means <- function(data, Predictor) {
    data_items <- data[, lapply(.SD, mean), by = list(Item, Condition),
                    .SDcols = c(Predictor)]
    data_means <- data[, lapply(.SD, mean), by = list(Condition),
                    .SDcols = c(Predictor)]
    data_means <- aggregate(as.formula(paste0(Predictor, "~ Condition")),
                    data, FUN = mean)

    list(data_items, data_means)
}

plot_density <- function(data, data_means, ylab, xlab, predictor,
                        leg_labs, leg_vals, ylimits, xbreaks) {
    p <- ggplot(data, aes(x = data[[predictor]], color = Condition,
                fill = Condition))
    p <- p + geom_density(alpha = 0.4) + theme_minimal()
    p <- p + ylim(ylimits)
    p <- p + geom_vline(data = data_means, aes(xintercept = .data[[predictor]],
                color = Condition), linetype = "dashed")
    p <- p + scale_color_manual(labels = leg_labs, values = leg_vals)
    p <- p + scale_fill_manual(labels = leg_labs, values = leg_vals)
    p <- p + scale_x_continuous(name = xlab, breaks = xbreaks)
    p <- p + labs(y = ylab, x = xlab)
    p <- p + theme(legend.position = "bottom")
    p
}

plot_rSPR <- function(
    data,
    file,
    yunit,
    title,
    ylims = NULL,
    modus = "Condition",
    leg_labs,
    leg_vals,
    omit_legend = FALSE,
    save_legend = FALSE
) {
    # always exclude the word two words before
    data <- data[Region != "Pre-critical-2",]
    data$Region <- factor(data$Region,
            levels = c("Pre-critical", "Critical",
                       "Spillover", "Post-spillover"))

    if (modus == "t-value"){
        sig_dt <- data[, c("Region", "Spec", "sig")]
        sig_dt$posit <- rep(seq(ylims[1] + 2, ylims[1] + 4,
                length = length(unique(data$Spec))),
                length = nrow(sig_dt))
        sig_dt$sig <- factor(sig_dt$sig, levels = c(1, 0),
            labels = c("sign", "insign"))
    } else if (modus == "coefficients") { # pass
    } else {
        data$Spec <- data$Condition
    }

    # For all plots
    p <- ggplot(data, aes(x = Region, y = logRT, color = Spec, group = Spec)) +
            geom_point(size = 2.5, shape = "cross") + geom_line(linewidth = 0.5)
    p <- p + theme_minimal()
    p <- p + geom_errorbar(aes(ymin = logRT - logRT_CI,
                ymax = logRT + logRT_CI), width = .1, linewidth = 0.3)
    # Conditional modifications 
    if (modus == "coefficients") { # coefficients
        p <- p + scale_color_manual(name = "Coefficients",
                    labels = leg_labs, values = leg_vals)
    } else if (modus == "t-value") { # t-values
        p <- p + geom_hline(yintercept = 0, linetype = "dashed")
        p <- p + scale_color_manual(name = "Z-value",
            labels = leg_labs, values = leg_vals)
        p <- p + geom_point(data = sig_dt,
                aes(x = Region, y = posit, shape = sig), size = 2.5)
        p <- p + scale_shape_manual(values = c(20, 32), name = "P-values",
            labels = c("Significant", "Nonsignificant"))
    } else { # RTs, Residuals
        p <- p + scale_color_manual(name = "Condition",
                    labels = leg_labs, values = leg_vals)
    }
    # limit y-axis
    if (is.vector(ylims) == TRUE) {
        p <- p + ylim(ylims[1], ylims[2])
    }
    # Apply theme
    p <- p + theme(plot.title = element_text(size = 8),
                    axis.text.x = element_text(size = 7),
                    legend.position = "bottom",
                    legend.text = element_text(size = 5),
                    legend.title = element_text(size = 6),
                    legend.box = "vertical",
                    legend.spacing.y = unit(-0.2, "cm"),
                    legend.margin = margin(0, 0, 0, 0),
                    legend.box.margin = margin(-10, -10, -10, -50))
    p <- p + labs(x = "Region", y = yunit, title = title)

    if (save_legend) {
        lgnd <- get_legend(p)
        file_trimmed <- strtrim(file, nchar(file) - 4)
        ggsave(paste0(file_trimmed, "_rtlegend.pdf"),
                lgnd,
                device = cairo_pdf,
                width = 3.5,
                height = 0.5)
    }
    if (omit_legend) {
        p <- p + theme(legend.position = "none")
    }

    if (file != FALSE) {
       ggsave(
            file,
            p, 
            device = cairo_pdf, 
            width = 3, 
            height = 3)
    } else {
       p
    }
}