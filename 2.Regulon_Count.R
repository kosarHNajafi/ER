# Load necessary libraries
library(stats)

# Load data
load("~/NCA.ER/NCA.METABRIC/NCA.Disc.All.Genes/NCA.Disc.All.Gene.Rproject/Disc.allgenes.Regulon.NA.RData")
pathway <- read.delim("C:/NCA.ER/NCA.METABRIC/NCA.Disc.All.Genes/Disc.Counts/Pathway_Counts_PathwaySheet.txt", header = TRUE)

# Inspect data structure
str(regulon.NA)
#List of 1508
#$ ADNP     : Named num [1:127] -0.0359 -0.0819 -0.0794 -0.0917 -0.0165 ...
#..- attr(*, "names")= chr [1:127] "ABHD6" "ACVRL1" "ADAMTS9" "ADCY4" ...
#$ ADNP2    : Named num [1:92] -0.08 -0.0452 -0.0806 -0.0333 0.029 ...

str(pathway)
#'data.frame':	1932 obs. of  2 variables:
#$ Metabolic_Pathways: chr  "Alanine.and.Aspartate.Metabolism" "Alanine.and.Aspartate.Metabolism" "Alanine.and.Aspartate.Metabolism" "Alanine.and.Aspartate.Metabolism" ...
#$ Genes             : chr  "GPT" "GPT2" "ACY3" "AGXT" ...

# Extract the list of genes from the metabolic pathway
metabolic_pathway_genes <- pathway$Genes

# Compare genes in regulons with metabolic pathway genes
LH <- lapply(regulon.NA, function(regulon) {
  intersect(names(regulon), metabolic_pathway_genes)
})

# Count how many target genes from each regulon are in the metabolic pathway list
LH_counts <- sapply(LH, length)

# Filter out regulons with matches
LH_filtered <- LH[LH_counts > 0]

# Calculate LH ratios with additional details
LH_ratio <- sapply(names(LH_filtered), function(regulon) {
  num_LH <- length(LH_filtered[[regulon]])  # Number of LH genes
  total_genes <- length(names(regulon.NA[[regulon]]))  # LT in the regulon
  ratio <- num_LH / total_genes  # Calculate ratio
  
  # Return details as a list
  return(c(Regulon = regulon, 
           LH = num_LH, 
           Total_Genes = total_genes, 
           Ratio = ratio))
})

# Convert ratios to a data frame
ratio_df <- as.data.frame(t(LH_ratio), stringsAsFactors = FALSE)
ratio_df$LH <- as.numeric(ratio_df$LH)
ratio_df$Total_Genes <- as.numeric(ratio_df$Total_Genes)
ratio_df$Ratio <- as.numeric(ratio_df$Ratio)

# Define global variables
PH <- 1420  # Total number of recon metabolic genes (Population Hits)
PT <- 30000 # Total number of genes in human genome (Population Total)

# Initialize a list to store contingency tables
contingency_tables <- list()

# Perform Fisher's exact test and save contingency tables
fisher_p_values <- sapply(names(LH_filtered), function(regulon) {
  num_LH <- length(LH_filtered[[regulon]])  # Number of LH genes (Hits)
  total_genes <- length(names(regulon.NA[[regulon]]))  # LT in the regulon (Total Genes)
  num_not_LH <- total_genes - num_LH  # Not hits in regulon
  num_in_metabolic <- PH - num_LH  # Remaining metabolic genes not in regulon
  num_outside_all <- PT - total_genes - num_in_metabolic  # Remaining genes not in regulon or metabolic pathways
  
  # Construct contingency table
  contingency_table <- matrix(c(
    num_LH, num_in_metabolic,         # In Metabolic Pathways
    num_not_LH, num_outside_all       # Not In Metabolic Pathways
  ), nrow = 2, byrow = TRUE,
  dimnames = list(c("In_Metabolic_Pathways", "Not_In_Metabolic_Pathways"),
                  c("In_Regulon", "Not_In_Regulon")))
  
  # Save the contingency table for later inspection
  contingency_tables[[regulon]] <<- contingency_table
  
  # Perform Fisher's exact test
  fisher_result <- fisher.test(contingency_table, alternative = "greater")
  return(fisher_result$p.value)
})

# Adjust p-values using Benjamini-Hochberg correction
fisher_adjusted_p_values <- p.adjust(fisher_p_values, method = "BH")

# Consolidate all results into a data frame
fisher_result_df <- data.frame(
  Regulon = names(LH_filtered),
  LH = sapply(LH_filtered, length),
  Total_Genes = sapply(names(LH_filtered), function(regulon) length(names(regulon.NA[[regulon]]))),
  Ratio = sapply(names(LH_filtered), function(regulon) {
    num_LH <- length(LH_filtered[[regulon]])
    total_genes <- length(names(regulon.NA[[regulon]]))
    return(num_LH / total_genes)
  }),
  Fisher_P_Value = fisher_p_values,
  Fisher_FDR = fisher_adjusted_p_values,
  stringsAsFactors = FALSE
)

# Sort by FDR for better interpretation
fisher_result_df <- fisher_result_df[order(fisher_result_df$Fisher_FDR), ]

# Save results to files
write.table(fisher_result_df, "disc.regulon_fisher_enrichment_results.txt", sep = "\t", row.names = FALSE)
save(contingency_tables, file = "contingency_tables.RData")

# View a specific table, e.g., for "ADNP"
contingency_tables[["ADNP"]]

disc.session.info <- sessionInfo()
sessionInfo()
#R version 4.4.1 (2024-06-14 ucrt)
#Platform: x86_64-w64-mingw32/x64
#Running under: Windows 11 x64 (build 26100)

#Matrix products: default


#locale:
#  [1] LC_COLLATE=English_United Arab Emirates.utf8 
#[2] LC_CTYPE=English_United Arab Emirates.utf8   
#[3] LC_MONETARY=English_United Arab Emirates.utf8
#[4] LC_NUMERIC=C                                 
#[5] LC_TIME=English_United Arab Emirates.utf8    

#time zone: Asia/Dubai
#tzcode source: internal

#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] dplyr_1.1.4    tidyr_1.3.1    reshape2_1.4.4

#loaded via a namespace (and not attached):
#  [1] utf8_1.2.4        R6_2.5.1          tidyselect_1.2.1  magrittr_2.0.3   
#[5] glue_1.8.0        stringr_1.5.1     tibble_3.2.1      pkgconfig_2.0.3  
#[9] generics_0.1.3    lifecycle_1.0.4   cli_3.6.3         fansi_1.0.6      
#[13] vctrs_0.6.5       compiler_4.4.1    plyr_1.8.9        purrr_1.0.2      
#[17] rstudioapi_0.17.1 tools_4.4.1       pillar_1.9.0      Rcpp_1.0.13-1    
#[21] rlang_1.1.4       stringi_1.8.4 

save(list = ls(),file = "1,2.Disc.TF.Choosing.RData")
save.image(file = "1,2.Disc.TF.chooseing.RData")
contingency_tables1 <- stack(contingency_tables)

write.table(contingency_tables,file = "contingency.txt",sep = "\t")
