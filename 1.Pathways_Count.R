# Load required libraries
library(dplyr)

#Pathway_GAUDE sheet from GAUDE geneset.excel
# Step 1: Convert Pathway_GAUDE to a Data Frame
# Filter out pathways with no genes and convert to long format
MP <- lapply(names(Pathway_GAUDE), function(pathway) {
  genes <- Pathway_GAUDE[[pathway]]
  if (!is.null(genes) && length(genes) > 0 && any(!is.na(genes))) {
    data.frame(
      Metabolic_Pathways = pathway,
      Genes = genes
    )
  } else {
    NULL  # Exclude pathways with no genes
  }
}) %>%
  bind_rows() %>%
  filter(!is.na(Genes), Genes != "")  # Remove remaining invalid entries

# Save the cleaned pathways to a file
write.table(
  MP, 
  file = "Pathway_Counts_PathwaySheet.txt", 
  sep = "\t", 
  row.names = FALSE
)

# Step 2: Identify Unique Genes
MP.Unique <- MP %>%
  group_by(Genes) %>%
  filter(n() == 1) %>%  # Keep only genes that appear once
  ungroup()

# Step 3: Count Unique Genes per Pathway
pathway_counts_MP.Unique <- MP.Unique %>%
  group_by(Metabolic_Pathways) %>%
  summarise(gene_count = n_distinct(Genes))  # Count unique genes per pathway

# View and save the unique gene counts
View(pathway_counts_MP.Unique)
sum(pathway_counts_MP.Unique$gene_count)

write.table(
  pathway_counts_MP.Unique, 
  file = file.path("G:/My Drive/ER/2.Metabolomics(GoogleDrive)/Metabolic Pathways/90 Metabolic Pathways.Patients", 
                   "Pathways_Unique_Count.txt"),
  sep = "\t", 
  row.names = FALSE
)

# Step 4: Count Total Genes per Pathway
pathway_counts_MP <- MP %>%
  group_by(Metabolic_Pathways) %>%
  summarise(gene_count = n_distinct(Genes))  # Count total unique genes per pathway

# View and save the total gene counts
View(pathway_counts_MP)
sum(pathway_counts_MP$gene_count)

write.table(
  pathway_counts_MP, 
  file = file.path("G:/My Drive/ER/2.Metabolomics(GoogleDrive)/Metabolic Pathways/90 Metabolic Pathways.Patients", 
                   "Pathways_Total_Count.txt"),
  sep = "\t", 
  row.names = FALSE
)

#Step5(optional): melt the regulon to see the regulon names in one column and target genes in another
library(tidyverse)

# Convert the list into a data frame with only TF and Target columns
regulon_df <- tibble(
  TF = names(regulon.NA), 
  TargetValues = regulon.NA
) %>%
  unnest_longer(TargetValues, indices_to = "Target") %>%
  select(TF, Target)

# View the result
head(regulon_df)
write.table(regulon_df, "reguln.df.txt", row.names = FALSE)
