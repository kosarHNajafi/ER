#Read data
Barplot <- read.delim("G:/My Drive/%ER/Barplot.txt")
View(Barplot)

# Load necessary libraries
library(ggplot2)
library(dplyr)

# Group by 'pathway' and count the unique 'gene' entries for each pathway
pathway_counts <- Barplot %>%
  group_by(Pathways) %>%
  summarise(gene_count = n_distinct(Genes)) # Count unique genes per pathway
View(pathway_counts)

#Remove the pathways not mentioned in 90PDS
barplot_90 <- pathway_counts[-c(24,60,67,69,94),]

#Add Vitamin_B as the avg of Vit_B6&12 
vitamin_B <- data.frame(
  Pathways= "Vitamin_B",
  gene_count = 4
)
barplot_90 <- rbind(barplot_90, vitamin_B)
barplot_90 <- barplot_90[-c(88,89),]

# Ensure data columns are correctly named and counted
View(barplot_90)

# Assuming 'pathway_counts' data frame already exists with 'Pathways' and 'gene_count' columns
ggplot(barplot_90, aes(x = reorder(Pathways, -gene_count), y = gene_count, fill = Pathways)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(
    x = "Metabolic Pathways", # Custom title for x-axis
    y = "Number of Genes", # Custom title for y-axis
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 5), # Font size for x-axis labels
    axis.text.y = element_text(size = 3),                       # Font size for y-axis labels
    axis.title.x = element_text(size = 10),                     # Font size for x-axis title
    axis.title.y = element_text(size = 10),                     # Font size for y-axis title
    plot.title = element_text(size = 14, face = "bold"),        # Font size and style for plot title
    legend.position = "none"                                    # Hide legend if not needed
  )

