########################################################################
# Section 1: Load Libraries
########################################################################
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)    # for combining plots
library(ggforce)    # for custom annotations if needed
library(ggtext)     # for rich text x-axis labels
library(grid)       # for grobs

########################################################################
# Section 2: Define Helper Functions
########################################################################

# Function to load and filter HR data from a given file (p-value < 0.05)
load_and_filter_HR_data <- function(file_path) {
  data <- read.delim(file_path, header = TRUE, stringsAsFactors = FALSE)
  filtered_data <- data %>%
    filter(Pvalue < 0.05) %>%
    mutate(HR_status = ifelse(HR >= 1, "HR>=1", "HR<1"))
  return(filtered_data)
}

# Function to process one survival type for a given ER group.
# It uses:
#   - survival: "OS", "DSS", "RFS", or "MFS"
#   - er_group: "ERpos" or "ERneg"
#   - discovery_HR and validation_HR: HR data for discovery and validation (already filtered)
#   - fdr_common: the common FDR pathways (FDR < 0.01) for this ER group (unused in this version)
#   - discovery_top10 and validation_top10: top10 FDR data for discovery and validation for this ER group
# It returns a data frame with:
#   - Pathway name, Survival type, ER group, averaged HR, maximum p-value,
#     a flag for whether both HR values are >=1, and top10 membership flags.
process_survival_ER <- function(survival, er_group,
                                discovery_HR, validation_HR,
                                fdr_common, discovery_top10, validation_top10) {
  # Find common pathways in HR data (after filtering) present in both discovery and validation
  common_hr <- intersect(discovery_HR$Gene, validation_HR$Gene)
  # Further restrict to those pathways that are also in the FDR data (from discovery and validation)
  D.pathways <- intersect(common_hr, discovery$pathways)
  V.pathways <- intersect(common_hr, validation$pathways)
  common_pathways <- intersect(D.pathways, V.pathways)
  
  results <- data.frame()
  for(p in common_pathways) {
    disc_row <- discovery_HR %>% filter(Gene == p)
    val_row  <- validation_HR %>% filter(Gene == p)
    if(nrow(disc_row) > 0 & nrow(val_row) > 0) {
      avg_HR <- mean(c(disc_row$HR, val_row$HR))
      # Bubble size based on the maximum p-value between discovery and validation
      max_p <- max(c(disc_row$Pvalue, val_row$Pvalue))
      
      # Flag for drawing a grey square outline: TRUE if both HR statuses are "HR>=1"
      if(disc_row$HR_status == "HR<1" & val_row$HR_status == "HR<1") {
        square_flag <- TRUE
      } else if(disc_row$HR_status == "HR>=1" & val_row$HR_status == "HR>=1") {
        square_flag <- FALSE
      } else {
        square_flag <- NA
      }
      
      
      # Determine top10 membership for annotation: if the pathway is in the top10 list
      top10_disc <- p %in% discovery_top10$pathways
      top10_val  <- p %in% validation_top10$pathways
      
      results <- rbind(results, data.frame(
        Pathway = p,
        Survival = survival,
        ER_Status = er_group,
        Avg_HR = avg_HR,
        Max_Pvalue = max_p,
        Grey_Square = square_flag,
        Top10_Discovery = top10_disc,
        Top10_Validation = top10_val,
        stringsAsFactors = FALSE
      ))
    }
  }
  return(results)
}

########################################################################
# Section 3: Load FDR Data and Split into ER Groups
########################################################################
# Load discovery and validation FDR datasets and rename the pathways column.
discovery <- read.delim("D.FDR.txt", check.names = FALSE)
colnames(discovery)[colnames(discovery) == "DISCOVERY; FDR<0.01; Metabolic pathways"] <- "pathways"

validation <- read.delim("V.FDR.txt", check.names = FALSE)
colnames(validation)[colnames(validation) == "VALIDATION; FDR<0.01; Metabolic pathways"] <- "pathways"

# Split the FDR datasets based on MD value:
# ERpos: MD >= 0; ERneg: MD < 0.
discovery_ERpos <- discovery %>% filter(MD >= 0)
discovery_ERneg <- discovery %>% filter(MD < 0)
validation_ERpos <- validation %>% filter(MD >= 0)
validation_ERneg <- validation %>% filter(MD < 0)

# For each ER group, get the common FDR set (i.e. pathways available in both discovery and validation FDR)
fdr_common_ERpos <- intersect(discovery_ERpos$pathways, validation_ERpos$pathways)
fdr_common_ERneg <- intersect(discovery_ERneg$pathways, validation_ERneg$pathways)

########################################################################
# Section 4: Define Top10 Pathways for ERpos and ERneg (Discovery and Validation)
########################################################################
# ERpos top10 pathways
top10_ERpos_names_discovery <- c("Methionine metabolism",
                                 "Bile acid biosynthesis",
                                 "Valine leucine and isoleucine metabolism",
                                 "Nucleotide metabolism",
                                 "O Glycan biosynthesis",
                                 "GABA shunt",
                                 "Ascorbate uptake",
                                 "Phosphoinositide Signalling",
                                 "Inositol phosphate metabolism",
                                 "Keratan sulfate biosynthesis")
top10_ERpos_names_validation <- c("Methionine metabolism",
                                  "Bile acid biosynthesis",
                                  "Valine leucine and isoleucine metabolism",
                                  "Nucleotide metabolism",
                                  "O Glycan biosynthesis",
                                  "Carnitine shuttle",
                                  "GABA shunt",
                                  "Phosphoinositide Signalling",
                                  "Keratan sulfate biosynthesis",
                                  "Ascorbate uptake")
discovery_top10_ERpos <- discovery_ERpos %>% filter(pathways %in% top10_ERpos_names_discovery)
validation_top10_ERpos <- validation_ERpos %>% filter(pathways %in% top10_ERpos_names_validation)

# ERneg top10 pathways
top10_ERneg_names_discovery <- c("Biotin metabolism",
                                 "Polyamines metabolism",
                                 "Folate metabolism",
                                 "Pyrimidine metabolism",
                                 "Fatty acid metabolism",
                                 "Cysteine metabolism",
                                 "Steroid metabolism",
                                 "NAD metabolism",
                                 "Tryptophan metabolism",
                                 "Ketone bodies metabolism")
top10_ERneg_names_validation <- c("Biotin metabolism",
                                  "Tryptophan metabolism",
                                  "Pyrimidine metabolism",
                                  "Cysteine metabolism",
                                  "NAD metabolism",
                                  "Steroid metabolism",
                                  "Fatty acid metabolism",
                                  "Polyamines metabolism",
                                  "Folate metabolism",
                                  "Ketone bodies metabolism")
discovery_top10_ERneg <- discovery_ERneg %>% filter(pathways %in% top10_ERneg_names_discovery)
validation_top10_ERneg <- validation_ERneg %>% filter(pathways %in% top10_ERneg_names_validation)

########################################################################
# Section 5: Load HR Data Files for Each Survival Type and ER Group
########################################################################
# ERpos HR files (Discovery and Validation) for OS, DSS, RFS, MFS
discovery_OS_HR_ERpos <- load_and_filter_HR_data("Discovery_OS_HR_ERP_08102024.txt")
validation_OS_HR_ERpos  <- load_and_filter_HR_data("Validation_OS_HR_ERP_08102024.txt")
discovery_DSS_HR_ERpos <- load_and_filter_HR_data("Discovery_DSS_HR_ERP_08102024.txt")
validation_DSS_HR_ERpos  <- load_and_filter_HR_data("Validation_DSS_HR_ERP_08102024.txt")
discovery_RFS_HR_ERpos <- load_and_filter_HR_data("Discovery_RFS_HR_ERP_08102024.txt")
validation_RFS_HR_ERpos  <- load_and_filter_HR_data("Validation_RFS_HR_ERP_08102024.txt")
discovery_MFS_HR_ERpos <- load_and_filter_HR_data("Discovery_MFS_HR_ERP_08102024.txt")
validation_MFS_HR_ERpos  <- load_and_filter_HR_data("Validation_MFS_HR_ERP_08102024.txt")

# ERneg HR files (Discovery and Validation) for OS, DSS, RFS, MFS
discovery_OS_HR_ERneg <- load_and_filter_HR_data("Discovery_OS_HR_ERN_08102024.txt")
validation_OS_HR_ERneg  <- load_and_filter_HR_data("Validation_OS_HR_ERN_08102024.txt")
discovery_DSS_HR_ERneg <- load_and_filter_HR_data("Discovery_DSS_HR_ERN_08102024.txt")
validation_DSS_HR_ERneg  <- load_and_filter_HR_data("Validation_DSS_HR_ERN_08102024.txt")
discovery_RFS_HR_ERneg <- load_and_filter_HR_data("Discovery_RFS_HR_ERN_08102024.txt")
validation_RFS_HR_ERneg  <- load_and_filter_HR_data("Validation_RFS_HR_ERN_08102024.txt")
discovery_MFS_HR_ERneg <- load_and_filter_HR_data("Discovery_MFS_HR_ERN_08102024.txt")
validation_MFS_HR_ERneg  <- load_and_filter_HR_data("Validation_MFS_HR_ERN_08102024.txt")

########################################################################
# Section 6: Process Each Survival Type for Both ER Groups
########################################################################
# Process each survival type separately for ERpos and ERneg.
res_OS_ERpos  <- process_survival_ER("OS", "ERpos",
                                     discovery_OS_HR_ERpos, validation_OS_HR_ERpos,
                                     fdr_common_ERpos, discovery_top10_ERpos, validation_top10_ERpos)
res_DSS_ERpos <- process_survival_ER("DSS", "ERpos",
                                     discovery_DSS_HR_ERpos, validation_DSS_HR_ERpos,
                                     fdr_common_ERpos, discovery_top10_ERpos, validation_top10_ERpos)
res_RFS_ERpos <- process_survival_ER("RFS", "ERpos",
                                     discovery_RFS_HR_ERpos, validation_RFS_HR_ERpos,
                                     fdr_common_ERpos, discovery_top10_ERpos, validation_top10_ERpos)
res_MFS_ERpos <- process_survival_ER("MFS", "ERpos",
                                     discovery_MFS_HR_ERpos, validation_MFS_HR_ERpos,
                                     fdr_common_ERpos, discovery_top10_ERpos, validation_top10_ERpos)

res_OS_ERneg  <- process_survival_ER("OS", "ERneg",
                                     discovery_OS_HR_ERneg, validation_OS_HR_ERneg,
                                     fdr_common_ERneg, discovery_top10_ERneg, validation_top10_ERneg)
res_DSS_ERneg <- process_survival_ER("DSS", "ERneg",
                                     discovery_DSS_HR_ERneg, validation_DSS_HR_ERneg,
                                     fdr_common_ERneg, discovery_top10_ERneg, validation_top10_ERneg)
res_RFS_ERneg <- process_survival_ER("RFS", "ERneg",
                                     discovery_RFS_HR_ERneg, validation_RFS_HR_ERneg,
                                     fdr_common_ERneg, discovery_top10_ERneg, validation_top10_ERneg)
res_MFS_ERneg <- process_survival_ER("MFS", "ERneg",
                                     discovery_MFS_HR_ERneg, validation_MFS_HR_ERneg,
                                     fdr_common_ERneg, discovery_top10_ERneg, validation_top10_ERneg)

# Combine results for all survival types and both ER groups.
common_results <- bind_rows(res_OS_ERpos, res_DSS_ERpos, res_RFS_ERpos, res_MFS_ERpos,
                            res_OS_ERneg, res_DSS_ERneg, res_RFS_ERneg, res_MFS_ERneg)

# --- New Step: Merge Duplicate Pathways from Both ER Groups ---
# If a pathway appears in both ERpos and ERneg, label it as "Common"
common_results_combined <- common_results %>%
  group_by(Pathway, Survival) %>%
    summarise(
    ER_Status = if(length(unique(ER_Status)) > 1) "Common" else first(ER_Status),
    Avg_HR = mean(Avg_HR),
    Max_Pvalue = max(Max_Pvalue),
    Grey_Square = if(all(!is.na(Grey_Square)) & all(Grey_Square == TRUE)) TRUE else NA,
    Top10_ERpos = first(Pathway) %in% c(top10_ERpos_names_discovery, top10_ERpos_names_validation),
    Top10_ERneg = first(Pathway) %in% c(top10_ERneg_names_discovery, top10_ERneg_names_validation),
    .groups = "drop"
  )


########################################################################
# Section 7: Create X-axis Label Mapping Based on Top10 Membership
#INACTIVATE SECTION 7 IF YOU DON'T WANT FONTCOLOR IN YOUR PLOT
########################################################################
label_df <- common_results_combined %>%
  group_by(Pathway) %>%
  summarise(
    ER_Status = if(length(unique(ER_Status)) > 1) "Common" else first(ER_Status),
    Top10 = any(Top10_ERpos) | any(Top10_ERneg),
    .groups = "drop"
  )

labels_map <- setNames(
  ifelse(label_df$Pathway %in% intersect(top10_ERpos_names_discovery, top10_ERpos_names_validation),
         paste0("<span style='color:steelblue'>", label_df$Pathway, "</span>"),
         ifelse(label_df$Pathway %in% intersect(top10_ERneg_names_discovery, top10_ERneg_names_validation),
                paste0("<span style='color:red'>", label_df$Pathway, "</span>"),
                label_df$Pathway)
  ),
  label_df$Pathway
)


########################################################################
# Section 8: Prepare Plot Data and Create the Bubble Plot
########################################################################
# Ensure that all pathways appear on the x-axis.
all_pathways <- sort(unique(common_results_combined$Pathway))
common_results_combined$Pathway <- factor(common_results_combined$Pathway, levels = all_pathways)

# Define manual colors for ER groups, including common
er_colors <- c("ERpos" = "steelblue", "ERneg" = "violet", "Common" = "yellow")

# To make higher p-values produce smaller circles, you can invert the size scale.
# One approach is to use -log(Max_Pvalue) if p-values are small.
# Here we assume that lower p-values are more significant, so we use -log10.
common_results_combined <- common_results_combined %>%
  mutate(SizeMetric = -log10(Max_Pvalue))

common_results_combined <- common_results_combined %>%
  mutate(SizeMetric2 = SizeMetric * 1.5)

# Build the base bubble plot.
p <- ggplot(common_results_combined, aes(x = Pathway, y = Survival)) +
  geom_point(aes(size = SizeMetric, color = ER_Status), shape = 16, fill = NA) +
  scale_color_manual(values = er_colors) +
  scale_size_continuous(name = "-log (P-value)", range = c(1, 10)) +
  labs(x = "Metabolic Pathway", y = "Survival Type", color = "ER Status") +
  theme_minimal() +
  theme(
    axis.text.x = element_markdown(size = 12, angle = 90, hjust = 1),  # Rotate x-axis labels
    axis.title.x = element_text(hjust = 0.5, vjust = -3.5),  # Adjust "Metabolic Pathway" title
    #vjust = -1 : it will be very close to y-axis
    axis.title.y = element_text(hjust = 0.5, vjust = 5, angle = 90),  # Rotate "Survival Type" title to vertical
    #vjust = 0.5 will be very close to x-axis legends
    plot.margin = margin(t = 10, r = 10, b = 10, l = 30)  # Adjust left margin to provide space for y-axis title
  ) +
  guides(color = guide_legend(order = 1), size = guide_legend(order = 2))

# Overlay grey square outlines for pathways where both HR values are >= 1.
# Multiply the size metric by 1.5 to make the square larger.
p <- p + geom_point(data = subset(common_results_combined, !is.na(Grey_Square) & Grey_Square == TRUE),
                    aes(x = Pathway, y = Survival, size = SizeMetric2),
                    shape = 22, fill = NA, color = "grey", stroke = 2,
                    show.legend = FALSE)

########################################################################
# Section 9: Apply the Custom X-axis Labels and Display the Plot
########################################################################
p <- p + scale_x_discrete(labels = labels_map)

#to remove the background grid
#p <- p + theme(
#  panel.grid.major = element_blank(),
#  panel.grid.minor = element_blank(),
#  axis.line = element_line(color = "black")
#)

# Create a custom legend for the square border.
# Here, we use a thicker border (lwd = 3) in the custom legend.
legend_plot <- ggdraw() + 
  draw_label("HR<1",  x = 0.94, y = 6.7, size = 10, fontface = "italic") +   #size to fit US letter pdf: x = 0.95, y = 6.2 
  draw_grob(rectGrob(gp = gpar(col = "grey", fill = NA, lwd = 3)),        
            x = 0.89, y = 6.5, width = 0.02, height = 0.3)                #size to fit US letter pdf:  x = 0.905, y = 6.1 
           
# Combine the main plot and custom legend.
final_plot <- plot_grid(p, legend_plot, ncol = 1, rel_heights = c(1, 0.1))
final_plot

# Optionally, save the plot to a file:
ggsave(filename = paste0("Final_HR_ER_", Sys.Date(), ".TIFF"), plot = final_plot, width = 10, height = 8)
se <- sessionInfo()
save(list = ls(),file = paste0("Final_HR_ER_",Sys.Date(),".RData"))
