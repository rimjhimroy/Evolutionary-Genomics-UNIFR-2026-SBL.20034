library(cowplot)
library(dplyr)
library(ggplot2)
library(tidyr)

##########################################
####### Plot the STRUCTURE results #######
##########################################

# This script is written like a worksheet.
# Read it from top to bottom.
#
# 1. Read the STRUCTURE result files.
# 2. Summarise support across values of K.
# 3. Plot the K = 4 ancestry results.

theme_set(
  theme_minimal(base_size = 11) +
    theme(
      plot.title.position = "plot",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position = "bottom"
    )
)

output_dir <- "outputs/01-structure-ddrad"

get_switzerland_outline <- function() {
  swiss_outline <- tryCatch(
    ggplot2::map_data("world", region = "Switzerland"),
    error = function(...) data.frame()
  )

  if (nrow(swiss_outline) > 0) {
    return(swiss_outline[, c("long", "lat", "group")])
  }

  data.frame(
    long = c(5.96, 6.08, 6.35, 6.73, 7.05, 7.55, 8.04, 8.54, 8.98, 9.55, 10.23, 10.45, 10.08, 9.48, 8.96, 8.47, 7.88, 7.36, 6.92, 6.48, 6.03, 5.96),
    lat = c(46.14, 46.40, 46.55, 46.52, 45.93, 45.93, 46.11, 46.14, 46.03, 46.35, 46.84, 47.55, 47.74, 47.78, 47.80, 47.73, 47.69, 47.61, 47.67, 47.63, 47.38, 46.14),
    group = 1
  )
}

build_pie_polygons <- function(location_table, cluster_columns, radius = 0.18, vertices_per_slice = 30) {
  polygon_rows <- vector("list", nrow(location_table) * length(cluster_columns))
  polygon_index <- 1L

  for (row_index in seq_len(nrow(location_table))) {
    one_population <- location_table[row_index, ]
    memberships <- as.numeric(one_population[, cluster_columns])
    membership_total <- sum(memberships, na.rm = TRUE)

    if (!is.finite(membership_total) || membership_total <= 0) {
      next
    }

    memberships <- memberships / membership_total
    start_angle <- 0

    for (cluster_index in seq_along(cluster_columns)) {
      membership_value <- memberships[cluster_index]

      if (!is.finite(membership_value) || membership_value <= 0) {
        next
      }

      end_angle <- start_angle + (2 * pi * membership_value)
      slice_angles <- seq(start_angle, end_angle, length.out = vertices_per_slice)
      polygon_rows[[polygon_index]] <- data.frame(
        long = c(one_population$X, one_population$X + radius * cos(slice_angles), one_population$X),
        lat = c(one_population$Y, one_population$Y + radius * sin(slice_angles), one_population$Y),
        pop = one_population$pop,
        cluster = cluster_columns[cluster_index],
        slice_id = paste(one_population$pop, cluster_columns[cluster_index], sep = "_"),
        stringsAsFactors = FALSE
      )
      polygon_index <- polygon_index + 1L
      start_angle <- end_angle
    }
  }

  do.call(rbind, polygon_rows[seq_len(polygon_index - 1L)])
}

# We keep K = 4 as the main teaching example because it gives a clear ancestry
# plot for discussion in class.
selected_k <- 4L
cluster_colours <- c("#264653", "#2A9D8F", "#E9C46A", "#E76F51", "#8AB17D", "#6D597A")
population_order <- c("A2G", "A2S", "V2B", "B2MV", "P2", "K2S2")

structure_zip <- "dataset-practical3/01-structure-ddrad/02-precomputed-structure-runs/ddRAD_results.zip"
structure_dir <- "dataset-practical3/01-structure-ddrad/02-precomputed-structure-runs/ddRAD_results"
samples_file <- "dataset-practical3/01-structure-ddrad/04-sample-metadata/samples_xy.txt"
elevation_file <- "dataset-practical3/02-ddrad-wgs-population-comparison/elevation.txt"

run_files <- list.files(structure_dir, pattern = "_f$", full.names = TRUE)
if (length(run_files) == 0) {
  stop("No STRUCTURE result files were found in: ", structure_dir, call. = FALSE)
}

############################################
####### 1. Summarise all K values ##########
############################################

# Each file is one STRUCTURE replicate.
# A replicate means STRUCTURE was run once for one value of K.
# From each file we want three things:
# 1. the value of K
# 2. the replicate number
# 3. the log-likelihood reported by STRUCTURE
# That log-likelihood is later used to compare support across K values.
run_summary <- data.frame(
  result_file = character(),
  file_name = character(),
  k = integer(),
  replicate = integer(),
  ln_prob = numeric(),
  stringsAsFactors = FALSE
)

for (result_file in run_files) {
  result_lines <- readLines(result_file, warn = FALSE)
  file_name <- basename(result_file)
  short_name <- sub("_f", "", file_name, fixed = TRUE)
  file_parts <- strsplit(short_name, "_", fixed = TRUE)[[1]]

  # Example file name: K4_rep_3
  # This means K = 4 and replicate = 3.
  k_value <- as.integer(sub("K", "", file_parts[1], fixed = TRUE))
  replicate_value <- as.integer(file_parts[3])

  # STRUCTURE writes a line called "Estimated Ln Prob of Data".
  # We keep that value because it is used in the model-choice summary.
  ln_prob_line <- result_lines[grepl("Estimated Ln Prob of Data", result_lines, fixed = TRUE)][1]
  ln_prob_text <- strsplit(ln_prob_line, "=", fixed = TRUE)[[1]][2]
  ln_prob_value <- as.numeric(trimws(ln_prob_text))

  one_row <- data.frame(
    result_file = result_file,
    file_name = short_name,
    k = k_value,
    replicate = replicate_value,
    ln_prob = ln_prob_value,
    stringsAsFactors = FALSE
  )

  run_summary <- rbind(run_summary, one_row)
}

# Now summarise each value of K across replicates.
# For each K we calculate:
# - the number of replicates
# - the mean log-likelihood
# - the standard deviation across replicates
evanno_table <- summarise(
  group_by(run_summary, k),
  replicates = n(),
  mean_ln_prob = mean(ln_prob, na.rm = TRUE),
  sd_ln_prob = sd(ln_prob, na.rm = TRUE)
)

evanno_table <- arrange(evanno_table, k)

evanno_table$ln_prime <- c(NA_real_, diff(evanno_table$mean_ln_prob))
evanno_table$ln_double_prime <- NA_real_
evanno_table$delta_k <- NA_real_

if (nrow(evanno_table) >= 3) {
  for (index in 2:(nrow(evanno_table) - 1)) {
    # Delta K is based on how quickly support changes between neighboring K
    # values, relative to the variation among replicates.
    second_difference <- abs(
      evanno_table$mean_ln_prob[index + 1] -
        (2 * evanno_table$mean_ln_prob[index]) +
        evanno_table$mean_ln_prob[index - 1]
    )
    evanno_table$ln_double_prime[index] <- second_difference
    if (!is.na(evanno_table$sd_ln_prob[index]) && evanno_table$sd_ln_prob[index] > 0) {
      evanno_table$delta_k[index] <- second_difference / evanno_table$sd_ln_prob[index]
    }
  }
}

best_k <- evanno_table$k[which.max(evanno_table$delta_k)]
best_k <- best_k[!is.na(best_k)][1]

write.table(
  evanno_table,
  file = "outputs/01-structure-ddrad/structure_evanno_summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

############################################
####### 2. Read the K = 4 results ##########
############################################

# Read the sample coordinates used for the map.
# These coordinates tell us where each sample was collected.
sampling_coordinates <- read.table(samples_file, header = TRUE, stringsAsFactors = FALSE)

# Read one small table giving the elevation of each population.
# We use this to classify populations as low or high elevation.
elevation_table <- read.table(elevation_file, header = FALSE, stringsAsFactors = FALSE)
names(elevation_table) <- c("pop", "elevation_m")

# Join coordinates and elevation information into one sample table.
sampling_coordinates <- merge(sampling_coordinates, elevation_table, by = "pop", all.x = TRUE, sort = FALSE)
sampling_coordinates$pop <- factor(sampling_coordinates$pop, levels = population_order)
sampling_coordinates$elevation_band <- ifelse(sampling_coordinates$elevation_m >= 1000, "High", "Low")

# Choose the single K = 4 replicate with the highest likelihood.
# The ancestry plot must come from one specific run, so we pick the best one.
k4_summary <- run_summary[run_summary$k == selected_k, ]
k4_summary <- k4_summary[order(-k4_summary$ln_prob, k4_summary$replicate), ]
representative_file <- k4_summary$result_file[1]

# Read ancestry coefficients from the chosen K = 4 output file.
# These coefficients describe the fraction of ancestry assigned to each cluster
# for each individual.
representative_lines <- readLines(representative_file, warn = FALSE)
start_line <- which(representative_lines == "Inferred ancestry of individuals:")
end_line <- which(representative_lines == "Estimated Allele Frequencies in each cluster")

if (length(start_line) != 1) {
  stop("Could not find the ancestry section in: ", representative_file, call. = FALSE)
}

if (length(end_line) != 1 || end_line <= start_line) {
  stop("Could not find the end of the ancestry section in: ", representative_file, call. = FALSE)
}

ancestry_lines <- representative_lines[(start_line + 2):(end_line - 2)]

# Keep only the lines that begin with an individual number.
# The STRUCTURE ancestry table is stored as plain text inside the result file.
# We therefore scan one line at a time and keep only the rows for samples.
sample_ids <- character()
missing_percent <- numeric()
q_matrix <- matrix(nrow = 0, ncol = selected_k)

for (one_line in ancestry_lines) {
  one_line_split <- scan(text = one_line, what = character(), quiet = TRUE)

  if (length(one_line_split) == 0) {
    next
  }

  if (is.na(suppressWarnings(as.numeric(one_line_split[1])))) {
    next
  }

  # Column 2 is the sample name.
  sample_ids <- c(sample_ids, one_line_split[2])

  # Column 3 stores the missing-data percentage for that individual.
  one_missing <- gsub("(", "", one_line_split[3], fixed = TRUE)
  one_missing <- gsub(")", "", one_missing, fixed = TRUE)
  missing_percent <- c(missing_percent, as.numeric(one_missing))

  # All values after the colon are the ancestry coefficients.
  colon_position <- which(one_line_split == ":")
  if (length(colon_position) != 1) {
    next
  }

  q_values <- as.numeric(one_line_split[(colon_position + 1):length(one_line_split)])
  q_matrix <- rbind(q_matrix, q_values)
}

if (ncol(q_matrix) != selected_k) {
  stop("Expected ", selected_k, " ancestry columns but found ", ncol(q_matrix), call. = FALSE)
}

membership_table <- data.frame(sample_id = sample_ids, missing_percent = missing_percent, q_matrix)
names(membership_table)[3:ncol(membership_table)] <- paste0("cluster_", seq_len(selected_k))

# Join the ancestry table to the sampling coordinates.
# After this join, one table contains geography, elevation, and ancestry.
plot_table <- merge(sampling_coordinates, membership_table, by.x = "Ind", by.y = "sample_id", all.x = TRUE, sort = FALSE)

cluster_columns <- paste0("cluster_", seq_len(selected_k))
if (any(!stats::complete.cases(plot_table[, cluster_columns]))) {
  stop("Some samples in samples_xy.txt were not matched to the STRUCTURE result file.", call. = FALSE)
}

write.table(
  plot_table,
  file = "outputs/01-structure-ddrad/structure_k4_membership.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

############################################
####### 3. Prepare the plot tables #########
############################################

# Make one row per population for the map.
# We average coordinates only to place one symbol per population.
population_summary <- summarise(
  group_by(plot_table, pop, elevation_band),
  X = mean(X),
  Y = mean(Y),
  n_samples = n(),
  .groups = "drop"
)

# Order individuals inside each population for a cleaner barplot.
# This makes the ancestry blocks easier to see by eye.
ordered_table <- plot_table[order(plot_table$pop, -plot_table$cluster_1, -plot_table$cluster_2, -plot_table$cluster_3, -plot_table$cluster_4), ]
ordered_table$sample_id <- factor(ordered_table$Ind, levels = ordered_table$Ind)

# Turn the ancestry table into a long table for plotting.
# ggplot draws stacked bars most easily when the data are in long format.
ancestry_wide <- ordered_table[, c("Ind", "pop", cluster_columns)]
ancestry_long <- pivot_longer(ancestry_wide, cols = all_of(cluster_columns), names_to = "cluster", values_to = "membership")
ancestry_long$sample_id <- factor(ancestry_long$Ind, levels = levels(ordered_table$sample_id))
ancestry_long$pop <- factor(ancestry_long$pop, levels = population_order)

# Calculate mean ancestry for each population.
# This gives a simpler population-level summary alongside the individual plot.
population_means <- summarise(
  group_by(plot_table, pop, elevation_band),
  cluster_1 = mean(cluster_1),
  cluster_2 = mean(cluster_2),
  cluster_3 = mean(cluster_3),
  cluster_4 = mean(cluster_4),
  .groups = "drop"
)
population_means$pop <- factor(population_means$pop, levels = population_order)
population_map_table <- merge(
  population_summary,
  population_means[, c("pop", cluster_columns)],
  by = "pop",
  all.x = TRUE,
  sort = FALSE
)
population_pies <- build_pie_polygons(population_map_table, cluster_columns)
switzerland_outline <- get_switzerland_outline()

population_means <- pivot_longer(population_means, cols = all_of(cluster_columns), names_to = "cluster", values_to = "membership")

############################################
############ 4. Make the figures ###########
############################################

# Figure 1. Population coordinates over a Switzerland reference map.
# Each pie chart shows the mean ancestry proportion for one population.
map_x_limits <- range(c(population_summary$X, switzerland_outline$long), na.rm = TRUE) + c(-0.5, 0.5)
map_y_limits <- range(c(population_summary$Y, switzerland_outline$lat), na.rm = TRUE) + c(-0.35, 0.35)

sampling_map <- ggplot() +
  geom_polygon(
    data = switzerland_outline,
    aes(x = long, y = lat, group = group),
    fill = "#EEF1EA",
    colour = "#7A8A6E",
    linewidth = 0.35
  ) +
  geom_point(
    data = population_summary,
    aes(x = X, y = Y, size = n_samples),
    inherit.aes = FALSE,
    shape = 21,
    fill = "white",
    colour = "#2F3E46",
    stroke = 0.3,
    alpha = 0.45
  ) +
  geom_polygon(
    data = population_pies,
    aes(x = long, y = lat, group = slice_id, fill = cluster),
    inherit.aes = FALSE,
    colour = "#1F1F1F",
    linewidth = 0.15
  ) +
  geom_text(
    data = population_summary,
    aes(x = X, y = Y, label = pop),
    inherit.aes = FALSE,
    nudge_y = 0.28,
    size = 3.1,
    fontface = "bold"
  ) +
  scale_fill_manual(values = cluster_colours[seq_along(cluster_columns)], labels = paste("Cluster", seq_along(cluster_columns))) +
  scale_size(range = c(5.5, 8.5)) +
  coord_equal(xlim = map_x_limits, ylim = map_y_limits, expand = FALSE) +
  labs(
    title = "Population coordinates and mean admixture proportions",
    subtitle = "Switzerland is shown for geographic reference",
    x = "Longitude",
    y = "Latitude",
    fill = "Mean ancestry",
    size = "Samples"
  )

# Figure 2A. How likelihood changes across K.
# This shows whether model support increases when more clusters are allowed.
ln_prob_plot <- ggplot(evanno_table, aes(x = k, y = mean_ln_prob)) +
  geom_line(linewidth = 0.6, colour = "#264653") +
  geom_point(size = 2.3, colour = "#264653") +
  geom_errorbar(
    aes(ymin = mean_ln_prob - sd_ln_prob, ymax = mean_ln_prob + sd_ln_prob),
    width = 0.15,
    linewidth = 0.45,
    colour = "#264653"
  ) +
  scale_x_continuous(breaks = evanno_table$k) +
  labs(title = "How support changes across K", x = "K", y = "Mean ln P(D)")

# Figure 2B. Delta K, used to choose a good value of K.
# The Evanno method highlights where the increase in support is strongest.
delta_plot <- ggplot(evanno_table, aes(x = k, y = delta_k)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
  geom_col(width = 0.75, fill = "#E76F51", na.rm = TRUE) +
  geom_vline(xintercept = best_k, linetype = "dashed", linewidth = 0.45, colour = "#6D597A") +
  scale_x_continuous(breaks = evanno_table$k) +
  labs(title = paste("Best support at K =", best_k), x = "K", y = "Delta K")

model_choice_plot <- cowplot::plot_grid(ln_prob_plot, delta_plot, ncol = 1, align = "v", rel_heights = c(1.1, 1))

# Figure 3. Ancestry coefficients for each individual.
# Each vertical bar is one sample.
# The colored fractions show estimated ancestry membership in each cluster.
structure_barplot <- ggplot(ancestry_long, aes(x = sample_id, y = membership, fill = cluster)) +
  geom_col(width = 0.95) +
  facet_grid(. ~ pop, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = cluster_colours[seq_along(cluster_columns)]) +
  labs(title = "Individual ancestry proportions", x = NULL, y = "Membership coefficient", fill = "Cluster") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing.x = grid::unit(0.15, "lines"),
    strip.background = element_rect(fill = "grey95", colour = NA),
    strip.text = element_text(face = "bold")
  )

# Figure 4. Average ancestry for each population.
# This is easier to compare across populations than the individual-level plot.
population_barplot <- ggplot(population_means, aes(x = pop, y = membership, fill = cluster)) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = cluster_colours[seq_along(cluster_columns)]) +
  labs(title = "Mean ancestry per population", x = NULL, y = "Mean membership coefficient", fill = "Cluster")

ggsave(
  filename = "outputs/01-structure-ddrad/population_sampling_map.png",
  plot = sampling_map,
  width = 8,
  height = 5,
  dpi = 320
)

ggsave(
  filename = "outputs/01-structure-ddrad/structure_model_choice.png",
  plot = model_choice_plot,
  width = 8,
  height = 7,
  dpi = 320
)

ggsave(
  filename = "outputs/01-structure-ddrad/structure_k4_barplot.png",
  plot = structure_barplot,
  width = 12,
  height = 5,
  dpi = 320
)

ggsave(
  filename = "outputs/01-structure-ddrad/structure_k4_population_summary.png",
  plot = population_barplot,
  width = 8,
  height = 4.5,
  dpi = 320
)

message("Finished STRUCTURE plotting workflow. Results were written to: ", output_dir)




