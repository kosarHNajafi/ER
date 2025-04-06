########################################
# 0. Load libraries
########################################
library(dplyr)
library(ggplot2)
library(patchwork)  # for arranging multiple plots
library(ggtext)     # to allow HTML formatting in legend and axis text
setwd("E:/2.ER/7.NCA.ER.Final/4.TF gene expression and PDS/3.PDS.Cor.Barplot/")
########################################
# 1. Define helper function to pick correlation
########################################
pick_larger_pval <- function(corD, pD, corV, pV, alpha = 0.05) {
  # Decide which correlation to keep based on significance & p-value.
  dSig <- (pD < alpha)
  vSig <- (pV < alpha)
  
  if (dSig && vSig) {
    # Both are significant; pick the one with the larger p-value.
    if (pD > pV) {
      return(c(corD, pD))
    } else {
      return(c(corV, pV))
    }
  } else if (dSig && !vSig) {
    return(c(corD, pD))
  } else if (!dSig && vSig) {
    return(c(corV, pV))
  } else {
    return(c(NA, NA))
  }
}

########################################
# 2. List TFs and ER status
########################################
tf_names  <- c("GATA3", "ESR1", "BCL11A", "CBX2", "YBX1")
er_status <- c("ERpos", "ERpos", "ERneg", "ERneg", "ERneg")
stopifnot(length(tf_names) == length(er_status))

########################################
# 3. Loop over TFs, read and process data
########################################
all_data_list <- list()

for (i in seq_along(tf_names)) {
  tf     <- tf_names[i]
  status <- er_status[i]
  
  # Construct file name, e.g., "GATA3_ERpos.txt"
  fileName <- paste0(tf, "_", status, ".txt")
  
  # Read data (assuming tab-delimited with header)
  df <- read.table(fileName, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Compute final correlation and corresponding p-value row by row
  df$Cor_final  <- NA
  df$Pval_final <- NA
  for (j in seq_len(nrow(df))) {
    out <- pick_larger_pval(
      corD = df$Cor_D.ER[j],
      pD   = df$P_value_D.ER[j],
      corV = df$Cor_V.ER[j],
      pV   = df$P_value_V.ER[j]
    )
    df$Cor_final[j]  <- as.numeric(out[1])
    df$Pval_final[j] <- as.numeric(out[2])
  }
  
  # Keep only the Pathway and computed Cor_final; add a TF column.
  df_final <- df %>%
    select(Pathway, Cor_final) %>%
    mutate(TF = tf)
  
  # Store the data frame for later combining.
  all_data_list[[tf]] <- df_final
}

########################################
# 4. Combine all TF data frames
########################################
all_data <- bind_rows(all_data_list)

########################################
# 5. Prepare data for bar plots
########################################
# (A) Top bar plot data: For each (Pathway, TF), sum the absolute correlations.
top_data <- all_data %>%
  group_by(Pathway, TF) %>%
  summarize(AbsCor = sum(abs(Cor_final), na.rm = TRUE), .groups = "drop")

# (B) Side bar plot data: Count the number of pathways for each TF.
side_data <- all_data %>%
  group_by(TF) %>%
  summarize(nPath = n(), .groups = "drop")

########################################
# 6. Factor ordering for meaningful layout
########################################
# Order pathways by total impact (sum of absolute correlations, descending)
pathway_order <- top_data %>%
  group_by(Pathway) %>%
  summarize(TotalImpact = sum(AbsCor), .groups = "drop") %>%
  arrange(desc(TotalImpact)) %>%
  pull(Pathway)

# Order TFs by how many pathways they appear in (descending)
tf_order <- side_data %>%
  arrange(desc(nPath)) %>%
  pull(TF)

# Define the custom order for TFs (GATA3, ESR1, YBX1, CBX2, BCL11A)
custom_tf_order <- c("BCL11A", "CBX2", "YBX1", "ESR1","GATA3")

# Apply this custom order to the TF column in all_data and top_data
all_data$TF <- factor(all_data$TF, levels = custom_tf_order)
top_data$TF <- factor(top_data$TF, levels = custom_tf_order)

# Convert character variables to factors with defined ordering.
all_data$Pathway <- factor(all_data$Pathway, levels = pathway_order)
top_data$Pathway <- factor(top_data$Pathway, levels = pathway_order)
side_data$TF     <- factor(side_data$TF, levels = custom_tf_order)

########################################
# 7. Define TF colors and HTML legend labels
########################################
tf_colors <- c("GATA3"  = "#230066",
               "ESR1"   = "deepskyblue3",
               "BCL11A" = "pink1",
               "CBX2"   = "red3",
               "YBX1"   = "violetred4" )


# HTML-colored legend labels for TFs.
tf_labels_legend <- c(
  "GATA3"  = "<span style='color:#230066'>GATA3</span>",
  "ESR1"   = "<span style='color:deepskyblue3'>ESR1</span>",
  "BCL11A" = "<span style='color:pink1'>BCL11A</span>",
  "CBX2"   = "<span style='color:red3'>CBX2</span>",
  "YBX1"   = "<span style='color:violetred4'>YBX1</span>"
)

########################################
# 9. Create the main matrix plot (Triangles) with grid background and adjusted clipping
########################################
# Create a sign indicator for triangle shape
all_data <- all_data %>%
  mutate(Sign = ifelse(Cor_final >= 0, "Positive", "Negative"))

# Now create the main plot again with colored y-axis labels
myplot <- ggplot(all_data, aes(x = Pathway, y = TF)) +
  geom_point(aes(fill = TF, shape = Sign, size = abs(Cor_final)), color = "black") +
  scale_shape_manual(values = c("Positive" = 24, "Negative" = 25)) +
  # Set the fill color to black for both Positive and Negative
  scale_fill_manual(values = tf_colors, labels = tf_labels_legend, guide = "none") +
  scale_size_area(name = "Correlation Magnitude", max_size = 10) + #, breaks = c(0.1, 1.25,0.5)
  scale_x_discrete(expand = c(0,0)) +
  theme_minimal(base_size = 20) +
  scale_y_discrete(
    expand = c(0, 0), 
    labels = function(x) {
      # Apply color to TF names in y-axis
      sapply(x, function(tf) {
        paste0("<span style='color:", tf_colors[tf], "'>", tf, "</span>")
      })
    }
  ) +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    #margin(top, right, bottom, left)
    plot.margin = margin(0, 30, -10, 30),
    axis.text.y = element_markdown(size = 18,face= "bold", margin = margin(r = 18)),  # Apply markdown formatting to y-axis labels
    axis.text.x = element_text(angle = 90, vjust = 0,hjust = 1, size = 20 ,margin = margin(t = 20)),  # Adjust x-axis labels (angle and size)
    axis.title.x = element_text(vjust = 0, size = 40),  # Adjust the position of the x-axis title without moving the plot
    legend.title = element_markdown(size = 25),
    legend.text = element_markdown(size = 28),  # Increase text size in legend
    legend.justification = c(0, 0.5),  # Adjust the position of the legend
    legend.margin = margin(r = 10,l = 9),  # Optional: add a bit of space to the left of the legend if needed
    legend.spacing = unit(3, "lines"),  # Add space between legend items 
    aspect.ratio = 1/4  # Space between TF rows
  ) +
  coord_cartesian(clip = "off") +  # So that it does not clip the top and bottom of the plot
  labs(size = "Correlation Magnitude", 
       shape = "Correlation Direction",
       x = "Correlation between gene expression and metabolic pathways") +
  guides(
    shape = guide_legend(
      keyheight = unit(5, "lines"),    # Increase vertical space for each key
      keywidth  = unit(3, "lines"),    # Increase horizontal space for each key
      override.aes = list(size = 10)   # Increase the triangle (shape) size
    )
)

# Save the plot
ggsave(paste0("TFs and PDS COR_", Sys.Date(), ".pdf"), myplot, width = 37, height = 17)

se <- sessionInfo()
# Save the workspace data
save(list = ls(), file = paste0("TFs and PDS COR_", Sys.Date(), ".RData"))
