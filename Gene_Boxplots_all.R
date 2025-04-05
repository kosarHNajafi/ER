# ---------------------------
# 1. Load Required Libraries
# ---------------------------
# ggplot2 is used for plotting and gridExtra for arranging multiple plots on one page.
library(ggplot2)
library(gridExtra)

# Set a random seed for reproducibility (e.g., for jitter positions)
#set.seed(123456)

# ---------------------------
# 2. Set Working Directory and Read Data
# ---------------------------

# Read the discovery and validation datasets from text files.
disc <- read.delim("Disc.Boxplot.txt", header = TRUE, stringsAsFactors = FALSE)
valid <- read.delim("Valid.Boxplot.txt", header = TRUE, stringsAsFactors = FALSE)

# ---------------------------
# 3. Pre-process Data
# ---------------------------
# Reorder the 'ER.status' factor so that "Positive" appears before "Negative".
disc$`ER.status` <- factor(disc$`ER.status`, levels = c("Positive", "Negative"))
valid$`ER.status` <- factor(valid$`ER.status`, levels = c("Positive", "Negative"))

# Define a custom color mapping for the two groups.
custom_colors <- c("Positive" = "#2C7BB6", "Negative" = "#D7191C")

# Identify the columns that contain gene expression data.
# (Here, it is assumed that the first two columns are not gene data.)
gene_cols <- colnames(disc)[3:ncol(disc)]

# ---------------------------
# 4. Define Helper Function for Significance Stars
# ---------------------------
# This function returns a string of asterisks based on the p-value.
add_significance <- function(p_value) {
  if (p_value < 0.0001) return("****")
  else if (p_value < 0.001) return("***")
  else if (p_value < 0.01) return("**")
  else if (p_value < 0.05) return("*")
  else return("ns")
}

# ---------------------------
# 5. Initialize Results Data Frame
# ---------------------------
# This data frame will store the t-test results for each dataset:
# Discovery, Validation, and the METABRIC dataset.
results_df <- data.frame(
  gene = character(),
  dataset = character(),
  p_value = numeric(),
  conf_lower = numeric(),
  conf_upper = numeric(),
  mean_diff = numeric(),
  t_statistic = numeric(),
  df = numeric(),
  significance = character(),
  stringsAsFactors = FALSE
)

# ---------------------------
# 6. Open PDF Device for Plot Output
# ---------------------------
# Increase the PDF width to accommodate three plots side by side.
pdf(paste0("Gene_Boxplots_", Sys.Date(), ".pdf"), width = 18, height = 6)

# ---------------------------
# 7. Loop Through Each Gene to Create Plots
# ---------------------------
for(gene in gene_cols) {
  
  # ---------------------------
  # a) Perform t-tests for Discovery and Validation Separately
  # ---------------------------
  # T-test for Discovery data:
  t_test_disc <- t.test(disc[[gene]] ~ disc$`ER.status`)
  disc_label <- add_significance(t_test_disc$p.value)
  
  # T-test for Validation data:
  t_test_valid <- t.test(valid[[gene]] ~ valid$`ER.status`)
  valid_label <- add_significance(t_test_valid$p.value)
  
  # ---------------------------
  # b) Create METABRIC Data and Perform t-test
  # ---------------------------
  # Combine Discovery and Validation data for the current gene.
  disc_combined <- data.frame(Expression = disc[[gene]], ER_status = disc$`ER.status`, Dataset = "Discovery")
  valid_combined <- data.frame(Expression = valid[[gene]], ER_status = valid$`ER.status`, Dataset = "Validation")
  combined_data <- rbind(disc_combined, valid_combined)
  
  # Perform t-test on the combined METABRIC data.
  t_test_combined <- t.test(combined_data$Expression ~ combined_data$ER_status)
  combined_label <- add_significance(t_test_combined$p.value)
  
  # ---------------------------
  # c) Calculate Mean Differences (Positive minus Negative)
  # ---------------------------
  mean_diff_disc <- as.numeric(t_test_disc$estimate[1] - t_test_disc$estimate[2])
  mean_diff_valid <- as.numeric(t_test_valid$estimate[1] - t_test_valid$estimate[2])
  mean_diff_combined <- as.numeric(t_test_combined$estimate[1] - t_test_combined$estimate[2])
  
  # ---------------------------
  # d) Append t-test Results to the Results Data Frame
  # ---------------------------
  results_df <- rbind(results_df, data.frame(
    gene = gene,
    dataset = "Discovery",
    p_value = t_test_disc$p.value,
    conf_lower = t_test_disc$conf.int[1],
    conf_upper = t_test_disc$conf.int[2],
    mean_diff = mean_diff_disc,
    t_statistic = t_test_disc$statistic,
    df = t_test_disc$parameter,
    significance = disc_label,
    stringsAsFactors = FALSE
  ))
  
  results_df <- rbind(results_df, data.frame(
    gene = gene,
    dataset = "Validation",
    p_value = t_test_valid$p.value,
    conf_lower = t_test_valid$conf.int[1],
    conf_upper = t_test_valid$conf.int[2],
    mean_diff = mean_diff_valid,
    t_statistic = t_test_valid$statistic,
    df = t_test_valid$parameter,
    significance = valid_label,
    stringsAsFactors = FALSE
  ))
  
  results_df <- rbind(results_df, data.frame(
    gene = gene,
    dataset = "METABRIC",
    p_value = t_test_combined$p.value,
    conf_lower = t_test_combined$conf.int[1],
    conf_upper = t_test_combined$conf.int[2],
    mean_diff = mean_diff_combined,
    t_statistic = t_test_combined$statistic,
    df = t_test_combined$parameter,
    significance = combined_label,
    stringsAsFactors = FALSE
  ))
  
  # ---------------------------
  # e) Calculate Positions for Significance Stars and Comparison Bars
  # ---------------------------
  # These calculations determine where to place the stars and bars above the boxplots.
  # For Discovery:
  range_disc <- diff(range(disc[[gene]], na.rm = TRUE))
  star_pos_disc <- max(disc[[gene]], na.rm = TRUE) + 0.05 * range_disc
  bar_pos_disc <- star_pos_disc - 0.03 * range_disc
  
  # For Validation:
  range_valid <- diff(range(valid[[gene]], na.rm = TRUE))
  star_pos_valid <- max(valid[[gene]], na.rm = TRUE) + 0.05 * range_valid
  bar_pos_valid <- star_pos_valid - 0.03 * range_valid
  
  # For METABRIC:
  range_combined <- diff(range(combined_data$Expression, na.rm = TRUE))
  star_pos_combined <- max(combined_data$Expression, na.rm = TRUE) + 0.05 * range_combined
  bar_pos_combined <- star_pos_combined - 0.03 * range_combined
  
  # ---------------------------
  # f) Create Boxplots for Each Dataset (with Jitter)
  # ---------------------------
  
  ## Discovery Plot (with jitter points)
  p_disc <- ggplot(disc, aes(x = `ER.status`, y = .data[[gene]], fill = `ER.status`)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    #geom_jitter(width = 0.2, shape = 21, color = "black") +  # Add individual data points
    labs(title = paste("Discovery_" ,gene, "Expression")) +
    scale_fill_manual(values = custom_colors, limits = c("Positive", "Negative")) +
    coord_cartesian(clip = "off") +
    geom_text(x = 1.5, y = star_pos_disc, label = disc_label, size = 5, fontface = "plain") +
    geom_segment(aes(x = 1, xend = 2, y = bar_pos_disc, yend = bar_pos_disc), linewidth = 0.6) +
    geom_segment(aes(x = 1, xend = 1, y = bar_pos_disc, yend = bar_pos_disc - 0.005 * range_disc), linewidth = 0.6) +
    geom_segment(aes(x = 2, xend = 2, y = bar_pos_disc, yend = bar_pos_disc - 0.005 * range_disc), linewidth = 0.6) +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 18),   # Change x axis labels font (e.g., Arial, bold)
          axis.text.y = element_text(size = 18),   # Change y axis tick label font size
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.position = "right",
          panel.grid = element_blank(),      # Remove gridlines
          axis.line = element_line(linewidth = 0.6, color = "black")  # Keep axis lines
    )
  
  ## Validation Plot (with jitter points)
  p_valid <- ggplot(valid, aes(x = `ER.status`, y = .data[[gene]], fill = `ER.status`)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    #geom_jitter(width = 0.2, shape = 21, color = "black") +  # Add individual data points
    labs(title = paste("Validation_",gene, "Expression")) +
    scale_fill_manual(values = custom_colors, limits = c("Positive", "Negative")) +
    coord_cartesian(clip = "off") +
    geom_text(x = 1.5, y = star_pos_valid, label = valid_label, size = 5, fontface = "plain") +
    geom_segment(aes(x = 1, xend = 2, y = bar_pos_valid, yend = bar_pos_valid), linewidth = 0.6) +
    geom_segment(aes(x = 1, xend = 1, y = bar_pos_valid, yend = bar_pos_valid - 0.005 * range_valid), linewidth = 0.6) +
    geom_segment(aes(x = 2, xend = 2, y = bar_pos_valid, yend = bar_pos_valid - 0.005 * range_valid), linewidth = 0.6) +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 18),   # Change x axis labels font (e.g., Arial, bold)
          axis.text.y = element_text(size = 18),   # Change y axis tick label font size
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.position = "right",
          panel.grid = element_blank(),      # Remove gridlines
          axis.line = element_line(linewidth = 0.6, color = "black") ) # Keep axis lines
  
  ## METABRIC Plot (with jitter points)
  p_combined <- ggplot(combined_data, aes(x = ER_status, y = Expression, fill = ER_status)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    #geom_jitter(width = 0.2, shape = 21, color = "black") +  # Add individual data points
    labs(title = paste(gene, "Expression")) +
    scale_fill_manual(values = custom_colors, limits = c("Positive", "Negative")) +
    coord_cartesian(clip = "off") +
    geom_text(x = 1.5, y = star_pos_combined, label = combined_label, size = 5, fontface = "plain") +
    geom_segment(aes(x = 1, xend = 2, y = bar_pos_combined, yend = bar_pos_combined), linewidth = 0.6) +
    geom_segment(aes(x = 1, xend = 1, y = bar_pos_combined, yend = bar_pos_combined - 0.005 * range_combined), linewidth = 0.6) +
    geom_segment(aes(x = 2, xend = 2, y = bar_pos_combined, yend = bar_pos_combined - 0.005 * range_combined), linewidth = 0.6) +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 18),   # Change x axis labels font (e.g., Arial, bold)
          axis.text.y = element_text(size = 18),   # Change y axis tick label font size
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.position = "right",
          panel.grid = element_blank(),      # Remove gridlines
          axis.line = element_line(linewidth = 0.6, color = "black")  # Keep axis lines
    )
  # ---------------------------
  # g) Arrange the Three Plots Side by Side on One Page
  # ---------------------------
  grid.arrange(p_disc, p_valid, p_combined, ncol = 3)
}

# ---------------------------
# 8. Close PDF Device and Save Results
# ---------------------------
# Close the PDF device so that the file is written.
dev.off()

# Write the t-test results to a text file.
write.table(results_df, file = paste0("Gene_ttest_results_", Sys.Date(), ".txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Save session information and all objects for future reference.
se <- sessionInfo()
save(list = ls(), file = paste0("Gene_Exp_Boxplot_", Sys.Date(), ".RData"))