library(cowplot)
library(dplyr)
library(ggplot2)
library(PopGenome)
library(tidyr)
library(vcfR)

###############################################
######## Scan differentiation on scaffold 1 ###
###############################################

# This script compares A2G and A2S along scaffold 1.
# The idea is simple:
# 1. Calculate window statistics for WGS.
# 2. Calculate the same statistics for ddRAD.
# 3. Plot the two datasets.

theme_set(
  theme_minimal(base_size = 11) +
    theme(
      plot.title.position = "plot",
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
)

output_dir <- "outputs/03-scaffold1-window-scan"

# Window size and jump define the genomic resolution of the scan.
# Larger windows smooth the signal more.
# Smaller jumps make neighboring windows overlap more strongly.
window_width <- 10000
window_jump <- 1000

# The full tables are kept for export, but plotting every window makes the
# figure rendering unnecessarily slow for teaching use.
plot_step <- 25

keep_finite_rows <- function(data, columns) {
  valid_rows <- Reduce(
    `&`,
    lapply(columns, function(column_name) {
      is.finite(data[[column_name]])
    })
  )

  data[valid_rows, , drop = FALSE]
}

pick_interesting_windows <- function(data, n_windows = 3, min_spacing_mb = 0.3) {
  ranked_windows <- data[is.finite(data$fst) & is.finite(data$midpoint_mb), , drop = FALSE]
  ranked_windows <- ranked_windows[order(-ranked_windows$fst, ranked_windows$midpoint_mb), , drop = FALSE]

  if (nrow(ranked_windows) == 0) {
    return(data.frame())
  }

  chosen_rows <- integer()

  for (row_index in seq_len(nrow(ranked_windows))) {
    if (length(chosen_rows) == 0 || all(abs(ranked_windows$midpoint_mb[row_index] - ranked_windows$midpoint_mb[chosen_rows]) >= min_spacing_mb)) {
      chosen_rows <- c(chosen_rows, row_index)
    }

    if (length(chosen_rows) >= n_windows) {
      break
    }
  }

  chosen_windows <- ranked_windows[chosen_rows, , drop = FALSE]
  chosen_windows$start_mb <- chosen_windows$start_pos / 1e6
  chosen_windows$end_mb <- chosen_windows$end_pos / 1e6
  offset <- max(0.015, 0.08 * max(chosen_windows$fst, na.rm = TRUE))
  chosen_windows$label_y <- chosen_windows$fst + seq(offset, by = offset * 0.7, length.out = nrow(chosen_windows))
  chosen_windows$label <- paste0("Peak ", seq_len(nrow(chosen_windows)), "\n", sprintf("%.2f-%.2f Mb", chosen_windows$start_mb, chosen_windows$end_mb))
  chosen_windows
}

############################################
############ 1. Analyse WGS data ###########
############################################

# First open the VCF just to read the sample names.
# We need the sample names so we can split the individuals into the two focal
# populations, A2G and A2S.
wgs_vcf <- read.vcfR("dataset-practical3/03-scaffold1-window-scan/wgs_vcf/scaffold1/A2_wgs_scaffold1.vcf", verbose = FALSE)
wgs_names <- vcfR2genlight(wgs_vcf)$ind.names
wgs_populations <- list(
  A2G = grep("A2G", wgs_names, value = TRUE, fixed = TRUE),
  A2S = grep("A2S", wgs_names, value = TRUE, fixed = TRUE)
)

# PopGenome now slides a window across the scaffold and calculates summary statistics.
# Each window is one small genomic segment.
# By moving window by window, we can see where diversity or differentiation
# changes along the scaffold.
wgs_data <- readData("dataset-practical3/03-scaffold1-window-scan/wgs_vcf/scaffold1", format = "VCF", include.unknown = TRUE)
wgs_windows <- sliding.window.transform(wgs_data, width = window_width, jump = window_jump, type = 2)
wgs_windows <- set.populations(wgs_windows, wgs_populations, diploid = TRUE)
wgs_fst_data <- F_ST.stats(wgs_windows, mode = "nucleotide")
wgs_fst_table <- as.data.frame(get.F_ST(wgs_fst_data, mode = FALSE, pairwise = TRUE)[1])
names(wgs_fst_table)[1] <- "fst"

# Negative Fst values are not interpreted as negative differentiation, so we
# truncate them to zero.
wgs_fst_table$fst[wgs_fst_table$fst < 0] <- 0

# PopGenome stores the window coordinates as text.
# Here we split that text into a start and an end position.
wgs_positions <- data.frame(
  start_pos = as.numeric(sub(" -.*", "", wgs_windows@region.names)),
  end_pos = as.numeric(sub(" :", "", sub(".*- ", "", wgs_windows@region.names)))
)

# Nucleotide diversity is pi.
# Dxy measures absolute divergence between populations.
# Tajima's D summarizes the allele frequency spectrum.
# We calculate all three because they answer different biological questions.
wgs_pi_data <- diversity.stats(wgs_windows, pi = TRUE)
wgs_between_data <- diversity.stats.between(wgs_pi_data)
wgs_neutrality_data <- suppressWarnings(neutrality.stats(wgs_between_data, FAST = TRUE))
# If PopGenome returns a table with an unexpected number of rows, we replace it
# with NA values so the final output table still has one row per window.
wgs_pi_table <- as.data.frame(wgs_pi_data@Pi)
if (nrow(wgs_pi_table) != nrow(wgs_fst_table)) {
  wgs_pi_table <- data.frame(pi_A2G = rep(NA_real_, nrow(wgs_fst_table)), pi_A2S = rep(NA_real_, nrow(wgs_fst_table)))
} else {
  names(wgs_pi_table) <- c("pi_A2G", "pi_A2S")
}

wgs_dxy_table <- as.data.frame(wgs_between_data@nuc.diversity.between)
if (nrow(wgs_dxy_table) != nrow(wgs_fst_table)) {
  wgs_dxy_table <- data.frame(dxy_A2G_A2S = rep(NA_real_, nrow(wgs_fst_table)))
} else {
  names(wgs_dxy_table) <- "dxy_A2G_A2S"
}

wgs_tajima_table <- as.data.frame(wgs_neutrality_data@Tajima.D)
if (nrow(wgs_tajima_table) != nrow(wgs_fst_table)) {
  wgs_tajima_table <- data.frame(tajima_A2G = rep(NA_real_, nrow(wgs_fst_table)), tajima_A2S = rep(NA_real_, nrow(wgs_fst_table)))
} else {
  names(wgs_tajima_table) <- c("tajima_A2G", "tajima_A2S")
}

# We also calculate overall pi and overall Tajima's D without splitting into two populations.
# This gives a whole-dataset reference in addition to the two population-
# specific summaries.
wgs_all_data <- readData("dataset-practical3/03-scaffold1-window-scan/wgs_vcf/scaffold1", format = "VCF", include.unknown = TRUE)
wgs_all_windows <- sliding.window.transform(wgs_all_data, width = window_width, jump = window_jump, type = 2)
wgs_all_windows <- diversity.stats(wgs_all_windows, pi = TRUE)
wgs_all_windows <- suppressWarnings(neutrality.stats(wgs_all_windows, FAST = TRUE))
wgs_overall_pi <- as.data.frame(wgs_all_windows@Pi)
if (nrow(wgs_overall_pi) != nrow(wgs_fst_table)) {
  wgs_overall_pi <- data.frame(pi_overall = rep(NA_real_, nrow(wgs_fst_table)))
} else {
  names(wgs_overall_pi) <- "pi_overall"
}

wgs_overall_tajima <- as.data.frame(wgs_all_windows@Tajima.D)
if (nrow(wgs_overall_tajima) != nrow(wgs_fst_table)) {
  wgs_overall_tajima <- data.frame(tajima_overall = rep(NA_real_, nrow(wgs_fst_table)))
} else {
  names(wgs_overall_tajima) <- "tajima_overall"
}

# Put everything into one table.
# Each row now represents one genomic window with all statistics attached.
wgs_stats <- cbind(wgs_positions, wgs_fst_table, wgs_pi_table, wgs_dxy_table, wgs_tajima_table, wgs_overall_pi, wgs_overall_tajima)
wgs_stats$midpoint <- (wgs_stats$start_pos + wgs_stats$end_pos) / 2
wgs_stats$midpoint_mb <- wgs_stats$midpoint / 1e6

############################################
########### 2. Analyse ddRAD data ##########
############################################

# Now calculate the same statistics for ddRAD.
# The code is intentionally parallel to the WGS section so students can compare
# the two marker systems directly.
ddrad_vcf <- read.vcfR("dataset-practical3/03-scaffold1-window-scan/ddRAD_vcf/scaffold_1/A2_ddRAD_scaffold_1.vcf", verbose = FALSE)
ddrad_names <- vcfR2genlight(ddrad_vcf)$ind.names
ddrad_populations <- list(
  A2G = grep("A2G", ddrad_names, value = TRUE, fixed = TRUE),
  A2S = grep("A2S", ddrad_names, value = TRUE, fixed = TRUE)
)
ddrad_data <- readData("dataset-practical3/03-scaffold1-window-scan/ddRAD_vcf/scaffold_1", format = "VCF", include.unknown = TRUE)
ddrad_windows <- sliding.window.transform(ddrad_data, width = window_width, jump = window_jump, type = 2)
ddrad_windows <- set.populations(ddrad_windows, ddrad_populations, diploid = TRUE)
ddrad_fst_data <- F_ST.stats(ddrad_windows, mode = "nucleotide")
ddrad_fst_table <- as.data.frame(get.F_ST(ddrad_fst_data, mode = FALSE, pairwise = TRUE)[1])
names(ddrad_fst_table)[1] <- "fst"
ddrad_fst_table$fst[ddrad_fst_table$fst < 0] <- 0
# Extract genomic window coordinates from the PopGenome object.
ddrad_positions <- data.frame(
  start_pos = as.numeric(sub(" -.*", "", ddrad_windows@region.names)),
  end_pos = as.numeric(sub(" :", "", sub(".*- ", "", ddrad_windows@region.names)))
)
ddrad_pi_data <- diversity.stats(ddrad_windows, pi = TRUE)
ddrad_between_data <- diversity.stats.between(ddrad_pi_data)
ddrad_neutrality_data <- suppressWarnings(neutrality.stats(ddrad_between_data, FAST = TRUE))
ddrad_pi_table <- as.data.frame(ddrad_pi_data@Pi)
if (nrow(ddrad_pi_table) != nrow(ddrad_fst_table)) {
  ddrad_pi_table <- data.frame(pi_A2G = rep(NA_real_, nrow(ddrad_fst_table)), pi_A2S = rep(NA_real_, nrow(ddrad_fst_table)))
} else {
  names(ddrad_pi_table) <- c("pi_A2G", "pi_A2S")
}

ddrad_dxy_table <- as.data.frame(ddrad_between_data@nuc.diversity.between)
if (nrow(ddrad_dxy_table) != nrow(ddrad_fst_table)) {
  ddrad_dxy_table <- data.frame(dxy_A2G_A2S = rep(NA_real_, nrow(ddrad_fst_table)))
} else {
  names(ddrad_dxy_table) <- "dxy_A2G_A2S"
}

ddrad_tajima_table <- as.data.frame(ddrad_neutrality_data@Tajima.D)
if (nrow(ddrad_tajima_table) != nrow(ddrad_fst_table)) {
  ddrad_tajima_table <- data.frame(tajima_A2G = rep(NA_real_, nrow(ddrad_fst_table)), tajima_A2S = rep(NA_real_, nrow(ddrad_fst_table)))
} else {
  names(ddrad_tajima_table) <- c("tajima_A2G", "tajima_A2S")
}

ddrad_all_data <- readData("dataset-practical3/03-scaffold1-window-scan/ddRAD_vcf/scaffold_1", format = "VCF", include.unknown = TRUE)
ddrad_all_windows <- sliding.window.transform(ddrad_all_data, width = window_width, jump = window_jump, type = 2)
ddrad_all_windows <- diversity.stats(ddrad_all_windows, pi = TRUE)
ddrad_all_windows <- suppressWarnings(neutrality.stats(ddrad_all_windows, FAST = TRUE))
ddrad_overall_pi <- as.data.frame(ddrad_all_windows@Pi)
if (nrow(ddrad_overall_pi) != nrow(ddrad_fst_table)) {
  ddrad_overall_pi <- data.frame(pi_overall = rep(NA_real_, nrow(ddrad_fst_table)))
} else {
  names(ddrad_overall_pi) <- "pi_overall"
}

ddrad_overall_tajima <- as.data.frame(ddrad_all_windows@Tajima.D)
if (nrow(ddrad_overall_tajima) != nrow(ddrad_fst_table)) {
  ddrad_overall_tajima <- data.frame(tajima_overall = rep(NA_real_, nrow(ddrad_fst_table)))
} else {
  names(ddrad_overall_tajima) <- "tajima_overall"
}

ddrad_stats <- cbind(ddrad_positions, ddrad_fst_table, ddrad_pi_table, ddrad_dxy_table, ddrad_tajima_table, ddrad_overall_pi, ddrad_overall_tajima)
ddrad_stats$midpoint <- (ddrad_stats$start_pos + ddrad_stats$end_pos) / 2
ddrad_stats$midpoint_mb <- ddrad_stats$midpoint / 1e6

write.table(wgs_stats, file = "outputs/03-scaffold1-window-scan/wgs_scaffold1_window_stats.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ddrad_stats, file = "outputs/03-scaffold1-window-scan/ddrad_scaffold1_window_stats.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

wgs_plot_stats <- wgs_stats[seq(1, nrow(wgs_stats), by = plot_step), , drop = FALSE]
ddrad_plot_stats <- ddrad_stats[seq(1, nrow(ddrad_stats), by = plot_step), , drop = FALSE]
wgs_fst_plot_stats <- keep_finite_rows(wgs_plot_stats, c("midpoint_mb", "fst"))
wgs_dxy_plot_stats <- keep_finite_rows(wgs_plot_stats, c("midpoint_mb", "dxy_A2G_A2S"))
ddrad_fst_plot_stats <- keep_finite_rows(ddrad_plot_stats, c("midpoint_mb", "fst"))
ddrad_dxy_plot_stats <- keep_finite_rows(ddrad_plot_stats, c("midpoint_mb", "dxy_A2G_A2S"))
wgs_highlight_windows <- pick_interesting_windows(wgs_stats)
ddrad_highlight_windows <- pick_interesting_windows(ddrad_stats)

############################################
############### 3. Make plots ##############
############################################

# Convert the tables to long format for plotting.
# Long format is easier for ggplot when we want one line per series.
wgs_pi_long <- wgs_plot_stats |>
  dplyr::select(midpoint_mb, pi_A2G, pi_A2S) |>
  dplyr::rename(A2G = pi_A2G, A2S = pi_A2S) |>
  tidyr::pivot_longer(cols = -midpoint_mb, names_to = "series", values_to = "value") |>
  dplyr::filter(is.finite(midpoint_mb), is.finite(value))

wgs_tajima_long <- wgs_plot_stats |>
  dplyr::select(midpoint_mb, tajima_A2G, tajima_A2S) |>
  dplyr::rename(A2G = tajima_A2G, A2S = tajima_A2S) |>
  tidyr::pivot_longer(cols = -midpoint_mb, names_to = "series", values_to = "value") |>
  dplyr::filter(is.finite(midpoint_mb), is.finite(value))

ddrad_pi_long <- ddrad_plot_stats |>
  dplyr::select(midpoint_mb, pi_A2G, pi_A2S) |>
  dplyr::rename(A2G = pi_A2G, A2S = pi_A2S) |>
  tidyr::pivot_longer(cols = -midpoint_mb, names_to = "series", values_to = "value") |>
  dplyr::filter(is.finite(midpoint_mb), is.finite(value))

ddrad_tajima_long <- ddrad_plot_stats |>
  dplyr::select(midpoint_mb, tajima_A2G, tajima_A2S) |>
  dplyr::rename(A2G = tajima_A2G, A2S = tajima_A2S) |>
  tidyr::pivot_longer(cols = -midpoint_mb, names_to = "series", values_to = "value") |>
  dplyr::filter(is.finite(midpoint_mb), is.finite(value))

wgs_pi_long$series <- factor(wgs_pi_long$series, levels = c("A2G", "A2S"))
wgs_tajima_long$series <- factor(wgs_tajima_long$series, levels = c("A2G", "A2S"))
ddrad_pi_long$series <- factor(ddrad_pi_long$series, levels = c("A2G", "A2S"))
ddrad_tajima_long$series <- factor(ddrad_tajima_long$series, levels = c("A2G", "A2S"))

# These four WGS panels show differentiation, diversity, divergence, and
# Tajima's D along the same scaffold.
wgs_fst_plot <- ggplot(wgs_fst_plot_stats, aes(x = midpoint_mb, y = fst)) +
  geom_rect(
    data = wgs_highlight_windows,
    aes(xmin = start_mb, xmax = end_mb, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = "#4C78A8",
    alpha = 0.08
  ) +
  geom_point(size = 0.8, alpha = 0.45, colour = "#7A3E00") +
  geom_smooth(method = "loess", se = FALSE, span = 0.18, linewidth = 0.55, colour = "#7A3E00") +
  geom_point(
    data = wgs_highlight_windows,
    aes(x = midpoint_mb, y = fst),
    inherit.aes = FALSE,
    shape = 21,
    size = 2.2,
    fill = "white",
    colour = "#7A3E00",
    stroke = 0.55
  ) +
  geom_label(
    data = wgs_highlight_windows,
    aes(x = midpoint_mb, y = label_y, label = label),
    inherit.aes = FALSE,
    size = 2.5,
    label.size = 0.15,
    fill = "white"
  ) +
  labs(title = "WGS scaffold 1 differentiation", x = NULL, y = "Fst")

wgs_pi_plot <- ggplot(wgs_pi_long, aes(x = midpoint_mb, y = value, colour = series)) +
  geom_rect(
    data = wgs_highlight_windows,
    aes(xmin = start_mb, xmax = end_mb, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = "#4C78A8",
    alpha = 0.08
  ) +
  geom_point(size = 0.75, alpha = 0.35) +
  geom_smooth(method = "loess", se = FALSE, span = 0.2, linewidth = 0.55) +
  facet_wrap(~series, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = c(A2G = "#1B9E77", A2S = "#D95F02")) +
  labs(title = "Within-population nucleotide diversity", x = NULL, y = expression(pi), colour = NULL) +
  theme(legend.position = "none")

wgs_dxy_plot <- ggplot(wgs_dxy_plot_stats, aes(x = midpoint_mb, y = dxy_A2G_A2S)) +
  geom_rect(
    data = wgs_highlight_windows,
    aes(xmin = start_mb, xmax = end_mb, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = "#4C78A8",
    alpha = 0.08
  ) +
  geom_point(size = 0.8, alpha = 0.45, colour = "#B22222") +
  geom_smooth(method = "loess", se = FALSE, span = 0.18, linewidth = 0.55, colour = "#B22222") +
  labs(title = "Between-population divergence", x = NULL, y = "Dxy")

wgs_tajima_plot <- ggplot(wgs_tajima_long, aes(x = midpoint_mb, y = value, colour = series)) +
  geom_rect(
    data = wgs_highlight_windows,
    aes(xmin = start_mb, xmax = end_mb, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = "#4C78A8",
    alpha = 0.08
  ) +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = "grey55") +
  geom_point(size = 0.75, alpha = 0.35) +
  geom_smooth(method = "loess", se = FALSE, span = 0.2, linewidth = 0.55) +
  facet_wrap(~series, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = c(A2G = "#1B9E77", A2S = "#D95F02")) +
  labs(title = "Tajima's D", x = "Scaffold 1 position (Mb)", y = "Tajima's D", colour = NULL) +
  theme(legend.position = "none")

wgs_summary_plot <- cowplot::plot_grid(wgs_fst_plot, wgs_pi_plot, wgs_dxy_plot, wgs_tajima_plot, ncol = 1, align = "v")

# The ddRAD panels show the same statistics for the reduced-representation
# marker set.
ddrad_fst_plot <- ggplot(ddrad_fst_plot_stats, aes(x = midpoint_mb, y = fst)) +
  geom_rect(
    data = ddrad_highlight_windows,
    aes(xmin = start_mb, xmax = end_mb, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = "#F4A261",
    alpha = 0.1
  ) +
  geom_point(size = 0.8, alpha = 0.45, colour = "#7A3E00") +
  geom_smooth(method = "loess", se = FALSE, span = 0.18, linewidth = 0.55, colour = "#7A3E00") +
  geom_point(
    data = ddrad_highlight_windows,
    aes(x = midpoint_mb, y = fst),
    inherit.aes = FALSE,
    shape = 21,
    size = 2.2,
    fill = "white",
    colour = "#7A3E00",
    stroke = 0.55
  ) +
  geom_label(
    data = ddrad_highlight_windows,
    aes(x = midpoint_mb, y = label_y, label = label),
    inherit.aes = FALSE,
    size = 2.5,
    label.size = 0.15,
    fill = "white"
  ) +
  labs(title = "ddRAD scaffold 1 differentiation", x = NULL, y = "Fst")

ddrad_pi_plot <- ggplot(ddrad_pi_long, aes(x = midpoint_mb, y = value, colour = series)) +
  geom_rect(
    data = ddrad_highlight_windows,
    aes(xmin = start_mb, xmax = end_mb, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = "#F4A261",
    alpha = 0.1
  ) +
  geom_point(size = 0.75, alpha = 0.35) +
  geom_smooth(method = "loess", se = FALSE, span = 0.2, linewidth = 0.55) +
  facet_wrap(~series, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = c(A2G = "#1B9E77", A2S = "#D95F02")) +
  labs(title = "Within-population nucleotide diversity", x = NULL, y = expression(pi), colour = NULL) +
  theme(legend.position = "none")

ddrad_dxy_plot <- ggplot(ddrad_dxy_plot_stats, aes(x = midpoint_mb, y = dxy_A2G_A2S)) +
  geom_rect(
    data = ddrad_highlight_windows,
    aes(xmin = start_mb, xmax = end_mb, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = "#F4A261",
    alpha = 0.1
  ) +
  geom_point(size = 0.8, alpha = 0.45, colour = "#B22222") +
  geom_smooth(method = "loess", se = FALSE, span = 0.18, linewidth = 0.55, colour = "#B22222") +
  labs(title = "Between-population divergence", x = NULL, y = "Dxy")

ddrad_tajima_plot <- ggplot(ddrad_tajima_long, aes(x = midpoint_mb, y = value, colour = series)) +
  geom_rect(
    data = ddrad_highlight_windows,
    aes(xmin = start_mb, xmax = end_mb, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = "#F4A261",
    alpha = 0.1
  ) +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = "grey55") +
  geom_point(size = 0.75, alpha = 0.35) +
  geom_smooth(method = "loess", se = FALSE, span = 0.2, linewidth = 0.55) +
  facet_wrap(~series, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = c(A2G = "#1B9E77", A2S = "#D95F02")) +
  labs(title = "Tajima's D", x = "Scaffold 1 position (Mb)", y = "Tajima's D", colour = NULL) +
  theme(legend.position = "none")

ddrad_summary_plot <- cowplot::plot_grid(ddrad_fst_plot, ddrad_pi_plot, ddrad_dxy_plot, ddrad_tajima_plot, ncol = 1, align = "v")

# This last figure isolates only Fst so students can compare ddRAD and WGS
# directly without the other statistics.
comparison_ddrad_plot <- ggplot(ddrad_fst_plot_stats, aes(x = midpoint_mb, y = fst)) +
  geom_rect(
    data = ddrad_highlight_windows,
    aes(xmin = start_mb, xmax = end_mb, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = "#F4A261",
    alpha = 0.1
  ) +
  geom_point(size = 0.8, alpha = 0.45, colour = "#7A3E00") +
  geom_smooth(method = "loess", se = FALSE, span = 0.18, linewidth = 0.55, colour = "#7A3E00") +
  labs(title = "ddRAD scaffold scan", x = NULL, y = "Fst")

comparison_wgs_plot <- ggplot(wgs_fst_plot_stats, aes(x = midpoint_mb, y = fst)) +
  geom_rect(
    data = wgs_highlight_windows,
    aes(xmin = start_mb, xmax = end_mb, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = "#4C78A8",
    alpha = 0.08
  ) +
  geom_point(size = 0.8, alpha = 0.45, colour = "#4C78A8") +
  geom_smooth(method = "loess", se = FALSE, span = 0.18, linewidth = 0.55, colour = "#4C78A8") +
  labs(title = "WGS scaffold scan", x = "Scaffold 1 position (Mb)", y = "Fst")

comparison_plot <- cowplot::plot_grid(comparison_ddrad_plot, comparison_wgs_plot, ncol = 1, align = "v")

ggsave(
  filename = "outputs/03-scaffold1-window-scan/wgs_scaffold1_summary.png",
  plot = wgs_summary_plot,
  width = 10,
  height = 14,
  dpi = 320
)

ggsave(
  filename = "outputs/03-scaffold1-window-scan/ddrad_scaffold1_summary.png",
  plot = ddrad_summary_plot,
  width = 10,
  height = 14,
  dpi = 320
)

ggsave(
  filename = "outputs/03-scaffold1-window-scan/ddrad_vs_wgs_fst.png",
  plot = comparison_plot,
  width = 9,
  height = 7,
  dpi = 320
)

message("Finished sliding-window differentiation analyses. Results were written to: ", output_dir)











