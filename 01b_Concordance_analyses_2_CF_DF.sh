#!/bin/bash
#------------------------------------------------------------------------------#
#                             PhyloReconcile                                   #
#                                                                              #
#                   Script_01b - Concordance Analysis 2                        #
#              Analysis of concordance and discordance factors                 #
#------------------------------------------------------------------------------#

# Set paths
BASE_DIR="$(pwd)/PhyloReconcile"
DATA_DIR="${BASE_DIR}/02_phylogenies/01b_30AX_ML_ASTRAL"
CONCORDANCE_DIR="${BASE_DIR}/03_Concordance_analyses"
SUPERTRI_OUT="${CONCORDANCE_DIR}/SuperTRI"
CF_DF_DIR="${CONCORDANCE_DIR}/CF_DF_analysis"
CONCAT_ALIGNMENT="${BASE_DIR}/01_initial_data/30AX_MSAs/30AX_concatenated/30AX_concatenated.fasta"

# Ensure Conda is initialized properly
source "$(conda info --base)/etc/profile.d/conda.sh"
# Activate the Conda environment
conda activate PhyloReconcile

# Testing assumption of ILS using concordance factors (CF) and discordance factors (DF) values
# Run CF/EF analyses 
# prepare data
mkdir -p "$CF_DF_DIR"
cd "$CF_DF_DIR"
mkdir -p ./CF_plots
cp "${DATA_DIR}/30AX_ASTRAL_ML_species_tree.cf.stat" "${CF_DF_DIR}"
cp "${SUPERTRI_OUT}/30AX_Nreps_tree.cf.stat" "${CF_DF_DIR}"
cp "${SUPERTRI_OUT}/30AX_MPP_tree.cf.stat" "${CF_DF_DIR}"
cp "${SUPERTRI_OUT}/30AX_Bootstrap_tree.cf.stat" "${CF_DF_DIR}"

# Prepare R script for CF/DF analysis
cat > "${CF_DF_DIR}/analyze_CF_DF.R" << EOF
library(viridis)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(entropy)
library(patchwork)

# Read the data
d <- read.delim("${CF_DF_DIR}/30AX_ASTRAL_ML_species_tree.cf.stat", header = TRUE, comment.char='#')
names(d)[20] <- "bootstrap"
names(d)[21] <- "branchlength"
d <- d[-1, ]

# Calculate correlations and plot relationships

# 1. Check correlation gCF/Bootstrap and sCF/bootstrap
correlation_value_gCF <- cor(d\$gCF, d\$bootstrap)
p_value_gCF <- cor.test(d\$gCF, d\$bootstrap)\$p.value

correlation_value_sCF <- cor(d\$sCF, d\$bootstrap)
p_value_sCF <- cor.test(d\$sCF, d\$bootstrap)\$p.value

# Plot
g1 <- ggplot(d, aes(x = gCF, y = sCF, label = ID)) + 
    geom_point(aes(colour = bootstrap)) + 
    scale_colour_viridis(direction = -1) + 
    xlim(0, 100) +
    ylim(0, 100) +
    geom_text_repel(max.overlaps = Inf) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggtitle("a")
    

ggsave("${CF_DF_DIR}/CF_plots/Relationship_boot_CF_nuDNA.png", plot = g1)

# 2. Plot relationship of branch length and bootstrap
correlation_value_brlength <- cor(d\$branchlength, d\$bootstrap)
p_value_brlength <- cor.test(d\$branchlength, d\$bootstrap)\$p.value

g2 <- ggplot(d, aes(x = branchlength, y = bootstrap, label = ID)) + 
    geom_point(aes(colour = bootstrap)) + 
    scale_colour_viridis(direction = -1) + 
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_text_repel(max.overlaps = Inf) +
    ggtitle("b") +
    annotate(
        "text",
        x = 2.5, 
        y = 55, 
        label = paste("Correlation (R) =", round(correlation_value_brlength, 2))
    ) +
    annotate(
        "text",
        x = 2.5, 
        y = 51, 
        label = paste("p-value =", formatC(p_value_brlength, format = "e"))
    )

ggsave("${CF_DF_DIR}/CF_plots/Correlation_brlength_bootstr.png", plot = g2)

# 3. Plot relationship of branch length and gCF/sCF
correlation_value_brlength_gCF <- cor(d\$branchlength, d\$gCF)
p_value_brlength_gCF <- cor.test(d\$branchlength, d\$gCF)\$p.value

g3 <- ggplot(d, aes(x = branchlength, y = gCF, label = ID)) + 
    geom_point(aes(colour = bootstrap)) + 
    scale_colour_viridis(direction = -1) + 
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
	geom_text_repel(max.overlaps = Inf) +
    ggtitle("c") +
    annotate(
        "text",
        x = 2.5, 
        y = 30, 
        label = paste("Correlation (R) =", round(correlation_value_brlength_gCF, 2))
    ) +
    annotate(
        "text",
        x = 2.5, 
        y = 23, 
        label = paste("p-value =", formatC(p_value_brlength_gCF, format = "e"))
    )

ggsave("${CF_DF_DIR}/CF_plots/Correlation_brlength_gCF.png", plot = g3)

# 4. Plot relationship of branch length and sCF
correlation_value_brlength_sCF <- cor(d\$branchlength, d\$sCF)
p_value_brlength_sCF <- cor.test(d\$branchlength, d\$sCF)\$p.value

g4 <- ggplot(d, aes(x = branchlength, y = sCF, label = ID)) + 
    geom_point(aes(colour = bootstrap)) + 
    scale_colour_viridis(direction = -1) + 
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_text_repel(max.overlaps = Inf) +
    ggtitle("d") +
    annotate(
        "text",
        x = 2.5, 
        y = 35, 
        label = paste("Correlation (R) =", round(correlation_value_brlength_sCF, 2))
    ) +
    annotate(
        "text",
        x = 2.5, 
        y = 30, 
        label = paste("p-value =", formatC(p_value_brlength_sCF, format = "e"))
    )

ggsave("${CF_DF_DIR}/CF_plots/Correlation_brlength_sCF.png", plot = g4)

# Plot all together
combined_plot <- (g1) + (g2) + 
    (g3) + (g4) +
    plot_layout(ncol = 2, nrow = 2)

ggsave("${CF_DF_DIR}/CF_plots/CF_analysis.png", plot = combined_plot, width = 17, height = 10, units = "in")


# 5. Calculate gEF_p and sEF_p
chisq <- function(DF1, DF2, N) {
    tryCatch({
        chisq.test(c(round(DF1*N)/100, round(DF2*N)/100))\$p.value
    },
    error = function(err) {
        return(1.0)
    })
}

e <- d %>%
    group_by(ID) %>%
    mutate(gEF_p = chisq(gDF1, gDF2, gN),
           sEF_p = chisq(sDF1, sDF2, sN))

# Filter significant results
f <- subset(data.frame(e), (gEF_p < 0.05 | sEF_p < 0.05))

# Save results
write.table(f, file = "${CF_DF_DIR}/significant_EF_p_values.txt", sep = "\t", quote = FALSE, row.names = FALSE)

EOF

# Run R script
Rscript "${CF_DF_DIR}/analyze_CF_DF.R"

echo "Analysis of Concordance and Discordance Factors (CF/DF) completed."

