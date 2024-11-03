# Required Libraries
library(edgeR)
library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(pheatmap)
library(grid)

# Read Raw Count Data
file_path <- "~/Downloads/normalized_sorted_Amerge1.csv"
data <- read_delim(file_path, delim = ";", col_names = FALSE, skip = 1)
colnames(data) <- c("chr", "start", "end", paste0("A", 1:9))

# Create Unique circRNA ID and Filter Data
data <- data %>%
  mutate(circRNA_ID = paste(chr, start, end, sep = "_")) %>%
  select(circRNA_ID, A1:A9) %>%
  column_to_rownames(var = "circRNA_ID")

# Data Filtering Based on Thresholds
filtered_data <- data %>%
  rowwise() %>%
  mutate(max_diff = max(c_across(A1:A9)) - min(c_across(A1:A9))) %>%
  filter(max_diff > 50 & max(c_across(A1:A9)) < 9000) %>%
  ungroup() %>%
  select(-max_diff)

# Save Filtered Data
write.csv(filtered_data, "~/Downloads/filtered_data_cleaned.csv", row.names = TRUE)

# Sample Metadata and Normalization
sample_metadata <- data.frame(
  row.names = colnames(filtered_data),
  condition = rep(c("control", "cmv", "zika"), times = 3),
  time = rep(c("24h", "48h", "72h"), each = 3)
)
sample_metadata$condition <- relevel(as.factor(sample_metadata$condition), ref = "control")
sample_metadata$time <- relevel(as.factor(sample_metadata$time), ref = "24h")

# Create DGEList and Estimate Dispersion
group <- factor(sample_metadata$condition)
y <- DGEList(counts = filtered_data, group = group)
y <- estimateDisp(y)

# Design Matrix and Fit Model with Interaction Terms
design_interaction <- model.matrix(~ condition * time, data = sample_metadata)
fit_interaction <- glmFit(y, design_interaction)

# Likelihood Ratio Tests for Each Comparison
lrt_list <- list(
  cmv_vs_control_24h = glmLRT(fit_interaction, coef = 2),
  zika_vs_control_24h = glmLRT(fit_interaction, coef = 3),
  cmv_vs_control_48h = glmLRT(fit_interaction, coef = 5),
  zika_vs_control_48h = glmLRT(fit_interaction, coef = 6),
  cmv_vs_control_72h = glmLRT(fit_interaction, coef = 7),
  zika_vs_control_72h = glmLRT(fit_interaction, coef = 8)
)

# Function to Generate Volcano Plot
generate_volcano_plot <- function(lrt_result, title, output_path) {
  results <- topTags(lrt_result, n = nrow(y))$table %>%
    mutate(circRNA_ID = rownames(.), Significant = ifelse(FDR < 0.05 & abs(logFC) > 1, "Yes", "No"))
  
  top_circrnas <- results %>%
    filter(Significant == "Yes")
  
  p <- ggplot(results, aes(x = logFC, y = -log10(PValue), color = Significant)) +
    geom_point() +
    geom_text_repel(data = top_circrnas, aes(label = circRNA_ID), size = 3, max.overlaps = Inf, box.padding = 0.4) +
    labs(title = title, x = "Log2 Fold Change", y = "-Log10 P-value") +
    scale_color_manual(values = c("No" = "red", "Yes" = "blue")) +
    theme_minimal()
  
  ggsave(output_path, plot = p, width = 14, height = 8, dpi = 300)
}

# Generate and Save Volcano Plots for Each Comparison
for (name in names(lrt_list)) {
  generate_volcano_plot(
    lrt_list[[name]],
    paste("Volcano Plot:", gsub("_", " ", name)),
    paste0("/Users/martina/Documents/Template_tesi_UniPD_DEI/Immagini/Filtered_Volcano_", name, ".png")
  )
}

# Function to Save Significant circRNAs as Table Image
save_table_as_image <- function(results, file_path, top_n = 15) {
  significant_circRNAs <- results %>%
    filter(FDR < 0.05 & abs(logFC) > 1) %>%
    arrange(FDR, desc(abs(logFC))) %>%
    head(top_n)
  
  table_grob <- tableGrob(
    significant_circRNAs,
    theme = ttheme_minimal(core = list(fg_params = list(cex = 0.6)), colhead = list(fg_params = list(cex = 0.7)))
  )
  
  ggsave(file_path, plot = table_grob, width = 10, height = 6, dpi = 300)
}

# Save Tables for Each Comparison
for (name in names(lrt_list)) {
  results <- topTags(lrt_list[[name]], n = nrow(y))$table
  save_table_as_image(
    results,
    paste0("/Users/martina/Documents/Template_tesi_UniPD_DEI/Immagini/Filtered_Top15_Table_", name, ".png")
  )
}

# Comparisons for Time Only
design_time <- model.matrix(~ time, data = sample_metadata)
fit_time <- glmFit(y, design_time)

lrt_time_list <- list(
  time_48h_vs_24h = glmLRT(fit_time, coef = 2),
  time_72h_vs_24h = glmLRT(fit_time, coef = 3)
)

# Save Volcano Plots and Tables for Time Comparisons
for (name in names(lrt_time_list)) {
  generate_volcano_plot(
    lrt_time_list[[name]],
    paste("Volcano Plot:", gsub("_", " ", name)),
    paste0("/Users/martina/Documents/Template_tesi_UniPD_DEI/Immagini/Filtered_Volcano_", name, ".png")
  )
  results <- topTags(lrt_time_list[[name]], n = nrow(y))$table
  save_table_as_image(
    results,
    paste0("/Users/martina/Documents/Template_tesi_UniPD_DEI/Immagini/Filtered_Top15_Table_", name, ".png")
  )
}

# Comparisons for Condition Only
design_condition <- model.matrix(~ condition, data = sample_metadata)
fit_condition <- glmFit(y, design_condition)

lrt_condition_list <- list(
  cmv_vs_control = glmLRT(fit_condition, coef = 2),
  zika_vs_control = glmLRT(fit_condition, coef = 3)
)

# Save Volcano Plots and Tables for Condition Comparisons
for (name in names(lrt_condition_list)) {
  generate_volcano_plot(
    lrt_condition_list[[name]],
    paste("Volcano Plot:", gsub("_", " ", name)),
    paste0("/Users/martina/Documents/Template_tesi_UniPD_DEI/Immagini/Filtered_Volcano_", name, ".png")
  )
  results <- topTags(lrt_condition_list[[name]], n = nrow(y))$table
  save_table_as_image(
    results,
    paste0("/Users/martina/Documents/Template_tesi_UniPD_DEI/Immagini/Filtered_Top15_Table_", name, ".png")
  )
}

# Find Common circRNAs Across Comparisons (Flexible Matching)
split_circRNA_ID <- function(circRNA_ID) {
  strsplit(circRNA_ID, "_")[[1]]
}

is_similar_circRNA <- function(id1, id2, tolerance = 2) {
  parts1 <- split_circRNA_ID(id1)
  parts2 <- split_circRNA_ID(id2)
  
  if (parts1[1] != parts2[1]) return(FALSE)
  start_diff <- abs(as.numeric(parts1[2]) - as.numeric(parts2[2]))
  end_diff <- abs(as.numeric(parts1[3]) - as.numeric(parts2[3]))
  
  return(start_diff <= tolerance & end_diff <= tolerance)
}

find_similar_circRNAs <- function(outliers, significant_ids, tolerance = 2) {
  matches <- c()
  for (outlier in outliers) {
    for (sig_id in significant_ids) {
      if (is_similar_circRNA(outlier, sig_id, tolerance)) {
        matches <- c(matches, outlier)
        break
      }
    }
  }
  return(matches)
}

# Example Usage for Finding Common circRNAs
common_circRNAs <- find_similar_circRNAs(
  rownames(filtered_data),
  rownames(filtered_data),
  tolerance = 2
)
print(common_circRNAs)
