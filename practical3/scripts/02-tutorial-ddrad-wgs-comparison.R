library(adegenet)
library(cowplot)
library(dplyr)
library(ggplot2)
library(hierfstat)
library(tidyr)
library(vcfR)
library(vegan)

###############################################
####### Compare ddRAD and WGS results #########
###############################################

# This script compares the same six populations in two datasets.
# The idea is simple:
# 1. Make a PCA for ddRAD.
# 2. Make a PCA for WGS.
# 3. Compare pairwise Fst, geography, and elevation.

theme_set(
        theme_minimal(base_size = 11) +
                theme(
                        plot.title.position = "plot",
                        legend.position = "bottom",
                        panel.grid.minor = element_blank()
                )
)

output_dir <- "outputs/02-ddrad-wgs-population-comparison"

# We keep the same population order in every figure and table so students can
# compare panels directly.
population_order <- c("A2G", "A2S", "V2B", "B2MV", "P2", "K2S2")
population_colours <- c(
        A2G = "#1B9E77",
        A2S = "#D95F02",
        V2B = "#7570B3",
        B2MV = "#E7298A",
        P2 = "#66A61E",
        K2S2 = "#E6AB02"
)

eligible_ellipse_rows <- function(scores_table) {
        scores_table |>
                dplyr::filter(is.finite(PC1), is.finite(PC2), !is.na(pop)) |>
                dplyr::group_by(pop) |>
                dplyr::filter(dplyr::n() >= 3, dplyr::n_distinct(PC1) > 1, dplyr::n_distinct(PC2) > 1) |>
                dplyr::ungroup()
}

############################################
########## 1. Read the input data ##########
############################################

# The metadata tables are tiny text files.
# The VCF files contain the SNP genotypes.
# We load both marker systems at the start so the rest of the script can ask
# the same questions of ddRAD and WGS.
elevation_table <- read.table("dataset-practical3/02-ddrad-wgs-population-comparison/elevation.txt", header = FALSE, stringsAsFactors = FALSE)
names(elevation_table) <- c("pop", "elevation_m")

geo_table <- read.table("dataset-practical3/02-ddrad-wgs-population-comparison/geo_dist.txt", header = TRUE, check.names = FALSE)
geo_matrix <- as.matrix(geo_table)

wgs_vcf <- read.vcfR(
        "dataset-practical3/02-ddrad-wgs-population-comparison/diploids_bisnps_filt1rem_minDP4nc_maxDP1009rm_filtMAC4MAF01_MD02_pruned.vcf.gz",
        verbose = FALSE
)

ddrad_vcf <- read.vcfR("dataset-practical3/02-ddrad-wgs-population-comparison/ddRAD_filtered_pruned.vcf.gz", verbose = FALSE)

############################################
########### 2. Analyse ddRAD data ##########
############################################

# First convert the VCF into a genlight object used by adegenet.
# A genlight object is useful for multivariate analyses such as PCA.
ddrad_genlight <- vcfR2genlight(ddrad_vcf)

# Sample names start with the population name.
# The only exception is K2S2, which starts with s2.
# We extract the population label directly from the sample names.
ddrad_sample_ids <- indNames(ddrad_genlight)
ddrad_name_parts <- strsplit(ddrad_sample_ids, "-", fixed = TRUE)
ddrad_population_names <- sapply(ddrad_name_parts, function(one_sample) one_sample[1])
ddrad_population_names[ddrad_population_names == "s2"] <- "K2S2"

ddrad_metadata <- data.frame(
        sample_id = ddrad_sample_ids,
        pop = ddrad_population_names,
        stringsAsFactors = FALSE
)

# Add the elevation of each population.
# This lets us compare low- and high-elevation groups in the figures.
ddrad_metadata <- merge(ddrad_metadata, elevation_table, by = "pop", all.x = TRUE, sort = FALSE)
ddrad_metadata$pop <- factor(ddrad_metadata$pop, levels = population_order)
ddrad_metadata$elevation_band <- ifelse(ddrad_metadata$elevation_m >= 1000, "High", "Low")

# Add marker names and population labels to the ddRAD genlight object.
# The marker names combine scaffold and genomic position.
locNames(ddrad_genlight) <- paste(ddrad_vcf@fix[, 1], ddrad_vcf@fix[, 2], sep = "_")
pop(ddrad_genlight) <- ddrad_metadata$pop

# Run PCA.
# PCA reduces thousands of SNPs to a few axes that summarize the main pattern
# of genetic structure.
ddrad_pca <- glPca(ddrad_genlight, nf = min(10L, max(3L, nInd(ddrad_genlight) - 1L)))

# Calculate how much variation each PCA axis explains.
# Students can use these percentages to judge how much of the structure is
# captured by PC1 and PC2.
ddrad_explained_variance <- data.frame(
        principal_component = paste0("PC", seq_along(ddrad_pca$eig)),
        explained_variance_percent = round(100 * ddrad_pca$eig / sum(ddrad_pca$eig), 2)
)

# Keep the first four PCA axes together with the sample information.
ddrad_scores <- as.data.frame(ddrad_pca$scores[, 1:4, drop = FALSE])
ddrad_scores$sample_id <- ddrad_metadata$sample_id
ddrad_scores$pop <- ddrad_metadata$pop
ddrad_scores$elevation_band <- ddrad_metadata$elevation_band

write.table(
        ddrad_explained_variance,
        file = "outputs/02-ddrad-wgs-population-comparison/ddrad_pca_explained_variance.tsv",
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
)

# Convert the same VCF to a genind object for population genetic distances.
ddrad_genind <- vcfR2genind(ddrad_vcf)
pop(ddrad_genind) <- ddrad_metadata$pop

# Pairwise Fst summarizes relative genetic differentiation between populations.
ddrad_fst <- as.matrix(genet.dist(ddrad_genind, method = "Fst"))

# Observed heterozygosity summarizes how much variation is present within each
# population. Here we average heterozygosity across loci.
ddrad_basic_stats <- hierfstat::basic.stats(ddrad_genind)
ddrad_heterozygosity_values <- colMeans(ddrad_basic_stats$Ho, na.rm = TRUE)
ddrad_heterozygosity <- data.frame(
        pop = names(ddrad_heterozygosity_values),
        observed_heterozygosity = as.numeric(ddrad_heterozygosity_values),
        dataset = "ddRAD",
        stringsAsFactors = FALSE
)

############################################
############ 3. Analyse WGS data ###########
############################################

# Now do exactly the same thing for the WGS dataset.
# The parallel structure of the script is intentional because students should
# be able to compare ddRAD and WGS step by step.
wgs_genlight <- vcfR2genlight(wgs_vcf)
wgs_sample_ids <- indNames(wgs_genlight)
wgs_name_parts <- strsplit(wgs_sample_ids, "-", fixed = TRUE)
wgs_population_names <- sapply(wgs_name_parts, function(one_sample) one_sample[1])
wgs_population_names[wgs_population_names == "s2"] <- "K2S2"

wgs_metadata <- data.frame(
        sample_id = wgs_sample_ids,
        pop = wgs_population_names,
        stringsAsFactors = FALSE
)

wgs_metadata <- merge(wgs_metadata, elevation_table, by = "pop", all.x = TRUE, sort = FALSE)
wgs_metadata$pop <- factor(wgs_metadata$pop, levels = population_order)
wgs_metadata$elevation_band <- ifelse(wgs_metadata$elevation_m >= 1000, "High", "Low")

locNames(wgs_genlight) <- paste(wgs_vcf@fix[, 1], wgs_vcf@fix[, 2], sep = "_")
pop(wgs_genlight) <- wgs_metadata$pop
wgs_pca <- glPca(wgs_genlight, nf = min(10L, max(3L, nInd(wgs_genlight) - 1L)))

wgs_explained_variance <- data.frame(
        principal_component = paste0("PC", seq_along(wgs_pca$eig)),
        explained_variance_percent = round(100 * wgs_pca$eig / sum(wgs_pca$eig), 2)
)

wgs_scores <- as.data.frame(wgs_pca$scores[, 1:4, drop = FALSE])
wgs_scores$sample_id <- wgs_metadata$sample_id
wgs_scores$pop <- wgs_metadata$pop
wgs_scores$elevation_band <- wgs_metadata$elevation_band

write.table(
        wgs_explained_variance,
        file = "outputs/02-ddrad-wgs-population-comparison/wgs_pca_explained_variance.tsv",
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
)

wgs_genind <- vcfR2genind(wgs_vcf)
pop(wgs_genind) <- wgs_metadata$pop
# Again, Fst summarizes how different the populations are from one another.
wgs_fst <- as.matrix(genet.dist(wgs_genind, method = "Fst"))

wgs_basic_stats <- hierfstat::basic.stats(wgs_genind)
wgs_heterozygosity_values <- colMeans(wgs_basic_stats$Ho, na.rm = TRUE)
wgs_heterozygosity <- data.frame(
        pop = names(wgs_heterozygosity_values),
        observed_heterozygosity = as.numeric(wgs_heterozygosity_values),
        dataset = "WGS",
        stringsAsFactors = FALSE
)

############################################
######## 4. Compare the two datasets #######
############################################

# Put all distance matrices in the same population order.
# This is essential because each row and column must refer to the same
# population in every matrix before we compare them.
reference_order <- population_order[population_order %in% rownames(ddrad_fst)]
ddrad_fst <- ddrad_fst[reference_order, reference_order, drop = FALSE]
wgs_fst <- wgs_fst[reference_order, reference_order, drop = FALSE]
geo_matrix <- geo_matrix[reference_order, reference_order, drop = FALSE]

# Convert the population elevation vector into a distance matrix.
# This matrix stores the elevation difference between every pair of populations.
elevation_matrix <- as.matrix(dist(elevation_table$elevation_m, method = "euclidean", diag = TRUE))
rownames(elevation_matrix) <- elevation_table$pop
colnames(elevation_matrix) <- elevation_table$pop
elevation_matrix <- elevation_matrix[reference_order, reference_order, drop = FALSE]

ddrad_heterozygosity$pop <- factor(ddrad_heterozygosity$pop, levels = population_order)
wgs_heterozygosity$pop <- factor(wgs_heterozygosity$pop, levels = population_order)
heterozygosity_table <- rbind(ddrad_heterozygosity, wgs_heterozygosity)

# Keep only one copy of each population pair for the ddRAD versus WGS Fst
# comparison scatterplot.
pair_index <- which(lower.tri(ddrad_fst), arr.ind = TRUE)
fst_correlation_table <- data.frame(
        population_1 = rownames(ddrad_fst)[pair_index[, 1]],
        population_2 = colnames(ddrad_fst)[pair_index[, 2]],
        ddrad_fst = ddrad_fst[lower.tri(ddrad_fst)],
        wgs_fst = wgs_fst[lower.tri(wgs_fst)],
        stringsAsFactors = FALSE
)

write.table(ddrad_fst, file = "outputs/02-ddrad-wgs-population-comparison/ddrad_pairwise_fst.tsv", sep = "\t", quote = FALSE, col.names = NA)
write.table(wgs_fst, file = "outputs/02-ddrad-wgs-population-comparison/wgs_pairwise_fst.tsv", sep = "\t", quote = FALSE, col.names = NA)
write.table(geo_matrix, file = "outputs/02-ddrad-wgs-population-comparison/geographic_distance_matrix.tsv", sep = "\t", quote = FALSE, col.names = NA)
write.table(elevation_matrix, file = "outputs/02-ddrad-wgs-population-comparison/elevation_distance_matrix.tsv", sep = "\t", quote = FALSE, col.names = NA)
write.table(heterozygosity_table, file = "outputs/02-ddrad-wgs-population-comparison/population_heterozygosity.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(fst_correlation_table, file = "outputs/02-ddrad-wgs-population-comparison/ddrad_vs_wgs_fst_pairs.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# IBD asks whether more distant populations are also more differentiated.
# IBE asks whether elevation differences matter after accounting for geography.
# Mantel tests compare whole distance matrices, not just single values.
ddrad_ibd <- mantel(geo_matrix, ddrad_fst, permutations = 1000)
wgs_ibd <- mantel(geo_matrix, wgs_fst, permutations = 1000)
ddrad_ibe <- mantel.partial(ddrad_fst, elevation_matrix, geo_matrix, permutations = 1000)
wgs_ibe <- mantel.partial(wgs_fst, elevation_matrix, geo_matrix, permutations = 1000)

test_summary <- data.frame(
        analysis = c(
                "ddRAD IBD",
                "WGS IBD",
                "ddRAD IBE controlling for geography",
                "WGS IBE controlling for geography"
        ),
        statistic = c(
                unname(ddrad_ibd$statistic),
                unname(wgs_ibd$statistic),
                unname(ddrad_ibe$statistic),
                unname(wgs_ibe$statistic)
        ),
        p_value = c(ddrad_ibd$signif, wgs_ibd$signif, ddrad_ibe$signif, wgs_ibe$signif)
)

write.table(
        test_summary,
        file = "outputs/02-ddrad-wgs-population-comparison/mantel_test_summary.tsv",
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
)

############################################
############### 5. Make plots ##############
############################################

# PCA plot: each point is one individual.
# If ddRAD and WGS tell a similar story, the broad clustering pattern should be
# similar in both panels.
ddrad_centroids <- aggregate(cbind(PC1, PC2) ~ pop, data = ddrad_scores, FUN = mean)
wgs_centroids <- aggregate(cbind(PC1, PC2) ~ pop, data = wgs_scores, FUN = mean)
ddrad_ellipse_scores <- eligible_ellipse_rows(ddrad_scores)
wgs_ellipse_scores <- eligible_ellipse_rows(wgs_scores)

ddrad_pca_plot <- ggplot(ddrad_scores, aes(x = PC1, y = PC2, colour = pop, shape = elevation_band)) +
        stat_ellipse(data = ddrad_ellipse_scores, aes(group = pop, colour = pop), type = "norm", linewidth = 0.4, alpha = 0.15, show.legend = FALSE) +
        geom_point(size = 2.4, alpha = 0.9) +
        scale_colour_manual(values = population_colours, drop = FALSE) +
        labs(
                title = "ddRAD population structure",
                x = paste0("PC1 (", ddrad_explained_variance$explained_variance_percent[1], "%)"),
                y = paste0("PC2 (", ddrad_explained_variance$explained_variance_percent[2], "%)"),
                colour = "Population",
                shape = "Elevation"
        )

ddrad_pca_centroid_plot <- ggplot(ddrad_scores, aes(x = PC1, y = PC2, colour = pop, shape = elevation_band)) +
        stat_ellipse(data = ddrad_ellipse_scores, aes(group = pop, colour = pop), type = "norm", linewidth = 0.4, alpha = 0.15, show.legend = FALSE) +
        geom_point(size = 2, alpha = 0.55) +
        geom_point(data = ddrad_centroids, aes(x = PC1, y = PC2, colour = pop), inherit.aes = FALSE, shape = 4, size = 4.5, stroke = 1.2, show.legend = FALSE) +
        geom_text(data = ddrad_centroids, aes(x = PC1, y = PC2, label = pop, colour = pop), inherit.aes = FALSE, nudge_y = 0.02, fontface = "bold", show.legend = FALSE) +
        scale_colour_manual(values = population_colours, drop = FALSE) +
        labs(
                title = "ddRAD PCA with population centroids",
                x = paste0("PC1 (", ddrad_explained_variance$explained_variance_percent[1], "%)"),
                y = paste0("PC2 (", ddrad_explained_variance$explained_variance_percent[2], "%)"),
                colour = "Population",
                shape = "Elevation"
        )

wgs_pca_plot <- ggplot(wgs_scores, aes(x = PC1, y = PC2, colour = pop, shape = elevation_band)) +
        stat_ellipse(data = wgs_ellipse_scores, aes(group = pop, colour = pop), type = "norm", linewidth = 0.4, alpha = 0.15, show.legend = FALSE) +
        geom_point(size = 2.4, alpha = 0.9) +
        scale_colour_manual(values = population_colours, drop = FALSE) +
        labs(
                title = "WGS population structure",
                x = paste0("PC1 (", wgs_explained_variance$explained_variance_percent[1], "%)"),
                y = paste0("PC2 (", wgs_explained_variance$explained_variance_percent[2], "%)"),
                colour = "Population",
                shape = "Elevation"
        )

wgs_pca_centroid_plot <- ggplot(wgs_scores, aes(x = PC1, y = PC2, colour = pop, shape = elevation_band)) +
        stat_ellipse(data = wgs_ellipse_scores, aes(group = pop, colour = pop), type = "norm", linewidth = 0.4, alpha = 0.15, show.legend = FALSE) +
        geom_point(size = 2, alpha = 0.55) +
        geom_point(data = wgs_centroids, aes(x = PC1, y = PC2, colour = pop), inherit.aes = FALSE, shape = 4, size = 4.5, stroke = 1.2, show.legend = FALSE) +
        geom_text(data = wgs_centroids, aes(x = PC1, y = PC2, label = pop, colour = pop), inherit.aes = FALSE, nudge_y = 0.02, fontface = "bold", show.legend = FALSE) +
        scale_colour_manual(values = population_colours, drop = FALSE) +
        labs(
                title = "WGS PCA with population centroids",
                x = paste0("PC1 (", wgs_explained_variance$explained_variance_percent[1], "%)"),
                y = paste0("PC2 (", wgs_explained_variance$explained_variance_percent[2], "%)"),
                colour = "Population",
                shape = "Elevation"
        )

pca_panel <- cowplot::plot_grid(
        ddrad_pca_plot,
        wgs_pca_plot,
        ncol = 2,
        labels = c("A", "B")
)

ggsave(
        filename = "outputs/02-ddrad-wgs-population-comparison/pca_overview.png",
        plot = pca_panel,
        width = 11,
        height = 5.5,
        dpi = 320
)

pca_centroid_panel <- cowplot::plot_grid(
        ddrad_pca_centroid_plot,
        wgs_pca_centroid_plot,
        ncol = 2,
        labels = c("A", "B")
)

ggsave(
        filename = "outputs/02-ddrad-wgs-population-comparison/pca_population_centroids.png",
        plot = pca_centroid_panel,
        width = 11,
        height = 5.5,
        dpi = 320
)

# Heterozygosity by population.
# This summarizes within-population diversity and makes it easy to compare the
# two marker systems with a simple barplot.
heterozygosity_plot <- ggplot(heterozygosity_table, aes(x = pop, y = observed_heterozygosity, fill = dataset)) +
        geom_col(position = "dodge", width = 0.75) +
        labs(title = "Observed heterozygosity by population", x = NULL, y = "Observed heterozygosity", fill = "Dataset") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
        filename = "outputs/02-ddrad-wgs-population-comparison/population_heterozygosity.png",
        plot = heterozygosity_plot,
        width = 8,
        height = 5,
        dpi = 320
)

# Pairwise Fst heatmaps.
# These plots make it easy to see which pairs of populations are weakly or
# strongly differentiated.
ddrad_fst_long <- as.data.frame(as.table(ddrad_fst), stringsAsFactors = FALSE)
names(ddrad_fst_long) <- c("population_1", "population_2", "value")
ddrad_fst_long$population_1 <- factor(ddrad_fst_long$population_1, levels = population_order)
ddrad_fst_long$population_2 <- factor(ddrad_fst_long$population_2, levels = rev(population_order))

wgs_fst_long <- as.data.frame(as.table(wgs_fst), stringsAsFactors = FALSE)
names(wgs_fst_long) <- c("population_1", "population_2", "value")
wgs_fst_long$population_1 <- factor(wgs_fst_long$population_1, levels = population_order)
wgs_fst_long$population_2 <- factor(wgs_fst_long$population_2, levels = rev(population_order))

ddrad_heatmap <- ggplot(ddrad_fst_long, aes(x = population_1, y = population_2, fill = value)) +
        geom_tile(colour = "white", linewidth = 0.35) +
        geom_text(aes(label = sprintf("%.3f", value)), size = 3) +
        scale_fill_gradient(low = "#F7FBFF", high = "#08306B") +
        labs(title = "ddRAD pairwise differentiation", x = NULL, y = NULL, fill = "Fst") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())

wgs_heatmap <- ggplot(wgs_fst_long, aes(x = population_1, y = population_2, fill = value)) +
        geom_tile(colour = "white", linewidth = 0.35) +
        geom_text(aes(label = sprintf("%.3f", value)), size = 3) +
        scale_fill_gradient(low = "#F7FBFF", high = "#08306B") +
        labs(title = "WGS pairwise differentiation", x = NULL, y = NULL, fill = "Fst") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())

fst_panel <- cowplot::plot_grid(
        ddrad_heatmap,
        wgs_heatmap,
        ncol = 2,
        labels = c("A", "B")
)

ggsave(
        filename = "outputs/02-ddrad-wgs-population-comparison/pairwise_fst_heatmaps.png",
        plot = fst_panel,
        width = 11,
        height = 5.5,
        dpi = 320
)

# ddRAD versus WGS Fst comparison.
# If the two marker systems agree, the population pairs should fall close to a
# positive diagonal relationship.
fst_correlation_plot <- ggplot(fst_correlation_table, aes(x = ddrad_fst, y = wgs_fst)) +
        geom_point(size = 2.6, alpha = 0.85, colour = "#264653") +
        geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, colour = "black") +
        geom_text(aes(label = paste(population_1, population_2, sep = "-")), nudge_y = 0.002, size = 3, check_overlap = TRUE) +
        labs(title = "ddRAD versus WGS pairwise Fst", x = "ddRAD pairwise Fst", y = "WGS pairwise Fst")

ggsave(
        filename = "outputs/02-ddrad-wgs-population-comparison/ddrad_vs_wgs_fst_correlation.png",
        plot = fst_correlation_plot,
        width = 7,
        height = 5.5,
        dpi = 320
)

# Mantel plots: each point is one pair of populations.
# The fitted line is only a visual guide. The actual result is the Mantel
# statistic and p-value shown in the panel label.
ddrad_geo_pairs <- data.frame(predictor = geo_matrix[upper.tri(geo_matrix)], fst = ddrad_fst[upper.tri(ddrad_fst)])
wgs_geo_pairs <- data.frame(predictor = geo_matrix[upper.tri(geo_matrix)], fst = wgs_fst[upper.tri(wgs_fst)])
ddrad_elev_pairs <- data.frame(predictor = elevation_matrix[upper.tri(elevation_matrix)], fst = ddrad_fst[upper.tri(ddrad_fst)])
wgs_elev_pairs <- data.frame(predictor = elevation_matrix[upper.tri(elevation_matrix)], fst = wgs_fst[upper.tri(wgs_fst)])

ddrad_geo_plot <- ggplot(ddrad_geo_pairs, aes(x = predictor, y = fst)) +
        geom_point(size = 2.6, alpha = 0.85, colour = "#1B9E77") +
        geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, colour = "black") +
        annotate("label", x = Inf, y = Inf, hjust = 1.02, vjust = 1.3,
                 label = sprintf("Mantel r = %.3f\np = %.3f", unname(ddrad_ibd$statistic), ddrad_ibd$signif),
                 size = 3, linewidth = 0.15, fill = "white") +
        labs(title = "ddRAD: isolation by distance", x = "Geographic distance (km)", y = "Pairwise Fst")

wgs_geo_plot <- ggplot(wgs_geo_pairs, aes(x = predictor, y = fst)) +
        geom_point(size = 2.6, alpha = 0.85, colour = "#7570B3") +
        geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, colour = "black") +
        annotate("label", x = Inf, y = Inf, hjust = 1.02, vjust = 1.3,
                 label = sprintf("Mantel r = %.3f\np = %.3f", unname(wgs_ibd$statistic), wgs_ibd$signif),
                 size = 3, linewidth = 0.15, fill = "white") +
        labs(title = "WGS: isolation by distance", x = "Geographic distance (km)", y = "Pairwise Fst")

ddrad_elev_plot <- ggplot(ddrad_elev_pairs, aes(x = predictor, y = fst)) +
        geom_point(size = 2.6, alpha = 0.85, colour = "#D95F02") +
        geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, colour = "black") +
        annotate("label", x = Inf, y = Inf, hjust = 1.02, vjust = 1.3,
                 label = sprintf("Partial Mantel r = %.3f\np = %.3f", unname(ddrad_ibe$statistic), ddrad_ibe$signif),
                 size = 3, linewidth = 0.15, fill = "white") +
        labs(title = "ddRAD: elevation contrast", x = "Elevation difference (m)", y = "Pairwise Fst")

wgs_elev_plot <- ggplot(wgs_elev_pairs, aes(x = predictor, y = fst)) +
        geom_point(size = 2.6, alpha = 0.85, colour = "#E7298A") +
        geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, colour = "black") +
        annotate("label", x = Inf, y = Inf, hjust = 1.02, vjust = 1.3,
                 label = sprintf("Partial Mantel r = %.3f\np = %.3f", unname(wgs_ibe$statistic), wgs_ibe$signif),
                 size = 3, linewidth = 0.15, fill = "white") +
        labs(title = "WGS: elevation contrast", x = "Elevation difference (m)", y = "Pairwise Fst")

mantel_panel <- cowplot::plot_grid(
        ddrad_geo_plot,
        wgs_geo_plot,
        ddrad_elev_plot,
        wgs_elev_plot,
        ncol = 2,
        labels = c("A", "B", "C", "D")
)

ggsave(
        filename = "outputs/02-ddrad-wgs-population-comparison/mantel_relationships.png",
        plot = mantel_panel,
        width = 11,
        height = 8,
        dpi = 320
)

message("Finished PCA, Fst, IBD, and IBE analyses. Results were written to: ", output_dir)


