

#getwd()  #[1] "C:/NCA.ER/NCA.METABRIC/Disc.all.final"
#---Load libraries-----------------
library(biomaRt)
library(RTN)
#---1.Load Gene Expression Data of 90 metabolic pathways and TFs = a mtrix for counts/assays---------------------------------
#load txt. Format for available genes in the discovery dataset
counts.Disc.genes <- read.delim("C:\\Users\\mohdb\\OneDrive - Towheed International School\\Documents\\NCA.ER\\data/Disc.Raw.Genes.all.txt",header = TRUE, as.is = TRUE)

sum(is.na(counts.Disc.genes[, 1]))  # Count the number of NA values
#[1] 1

which(is.na(counts.Disc.genes[, 1]))  # Identify their positions
#[1] 11028 in the main file not sorted, only metabolic clusters row was removed.
# [1] 10390 in the sorted data based on gene names without 1 and 0 metabolic clusters.

counts.Disc.genes <- counts.Disc.genes[-11028,]  # Remove the row with 'NA' as gene name
rownames(counts.Disc.genes) <- counts.Disc.genes[,1]
counts.Disc.genes <- counts.Disc.genes[,-1]
dim(counts.Disc.genes)
#[1] 19219  993

str(counts.Disc.genes)
#'data.frame':	19219 obs. of  993 variables:
#$ MB.0005: chr  "6.24393549" "5.228923847" "6.781157168" "5.583909882" ...
#$ MB.0006: num  6.22 5.91 6.75 5.36 6.08 ...
#$ MB.0008: num  6.06 5.3 7.66 5.3 5.69 ...
#$ MB.0014: num  6.16 5.4 6.68 5.85 5.65 ...

View(counts.Disc.genes)

colnames(counts.Disc.genes) <- gsub("\\.", "_", colnames(counts.Disc.genes))

#Sort data by columnnames(Sample IDs) and rownames(Gene IDs)
counts.Disc.genes <- counts.Disc.genes[order(rownames(counts.Disc.genes)),order(colnames(counts.Disc.genes))]
View(counts.Disc.genes)
head(counts.Disc.genes)
str(counts.Disc.genes)
#'data.frame':	19219 obs. of  993 variables:
#$ MB_0005: chr  "5.059348923" "9.924912584" "5.421198528" "6.24393549" ...
#$ MB_0006: num  5.37 11.31 5.14 6.22 5.91 ...
#$ MB_0008: num  5.28 9.83 5.56 6.06 5.3 ...
#$ MB_0014: num  5.53 10.68 5.36 6.16 5.4 ...
#$ MB_0020: num  5.28 9.31 5.69 6.16 5.53 ...

#check for NAs
sum(is.na(counts.Disc.genes))
#[1] 10

which(is.na((rownames(counts.Disc.genes))))
#integer(0)

#check for non-numeric elements
counts.Disc.genes <- as.matrix(counts.Disc.genes)

# Find elements that will turn into NA when coerced to numeric
non_numeric_elements1 <- counts.Disc.genes[is.na(as.numeric(counts.Disc.genes))]
#Warning message:
#NAs introduced by coercion 

length(non_numeric_elements1)
#[1] 204

write.table(non_numeric_elements1,"non_numeric_all.txt",sep = "\t",row.names = FALSE)

# Remove everything after the first space, including trailing digits and an optional dot
counts.Disc.genes <- gsub("\\s+\\d+\\.?$", "", counts.Disc.genes)

# Find elements that will turn into NA when coerced to numeric
non_numeric_elements2 <- counts.Disc.genes[is.na(as.numeric(counts.Disc.genes))]
#Warning message:
#NAs introduced by coercion 

length(non_numeric_elements2)
#[1] 13

print(non_numeric_elements2)
#[1] "11.60540448 7.88" "5.455753203 5.36" "5.3912438332 5.4" NA                
#[5] NA                 NA                 NA                 NA                
#[9] NA                 NA                 NA                 NA                
#[13] NA  

write.table(counts.Disc.genes,"pre.curated.counts.Disc.Genes.txt",sep = "\t",row.names = TRUE)

#Three non_numeric values the end digits were removed manually
counts.Disc.genes2 <- read.delim("C:\\Users\\mohdb\\OneDrive - Towheed International School\\Documents\\NCA.ER\\Finalyzing\\1.Curating Data\\Disc.Curating.Data\\2.pre.counts.Disc.all.Genes.txt",header = TRUE, as.is = TRUE, row.names = 1)
str(counts.Disc.genes2)
#'data.frame':	19219 obs. of  993 variables:
#$ MB_0005: num  5.06 9.92 5.42 6.24 5.23 ...
#$ MB_0006: num  5.37 11.31 5.14 6.22 5.91 ...
#$ MB_0008: num  5.28 9.83 5.56 6.06 5.3 ...
#$ MB_0014: num  5.53 10.68 5.36 6.16 5.4 ...
#$ MB_0020: num  5.28 9.31 5.69 6.16 5.53 ...

counts.Disc.genes <- as.matrix(counts.Disc.genes2)
non_numeric_elements3 <- counts.Disc.genes[is.na(as.numeric(counts.Disc.genes))]
#Warning message:
#NAs introduced by coercion 

length(non_numeric_elements3)
#[1] 10

print(non_numeric_elements3)
#[1] NA NA NA NA NA NA NA NA NA NA

# Convert the matrix to numeric values (after applying gsub)
counts.Disc.genes <- apply(counts.Disc.genes, c(1, 2), as.numeric)

# Check the structure of the matrix after conversion
str(counts.Disc.genes)
#num [1:19219, 1:993] 5.06 9.92 5.42 6.24 5.23 ...
#- attr(*, "dimnames")=List of 2
#..$ : chr [1:19219] "A1CF" "A2M" "A2ML1" "A4GALT" ...
#..$ : chr [1:993] "MB_0005" "MB_0006" "MB_0008" "MB_0014" ...

counts.Disc.genes <- t(apply(counts.Disc.genes, 1, function(row) {
  # Replace NA in the row with the mean of non-NA values in that row
  row[is.na(row)] <- mean(row, na.rm = TRUE)
  return(row)
}))

View(counts.Disc.genes)

non_numeric_elements4 <- counts.Disc.genes[is.na(as.numeric(counts.Disc.genes))]
length(non_numeric_elements4)
#[1] 0

print(non_numeric_elements4)
#numeric(0)

write.table(counts.Disc.genes,file = "3.Disc.Genes.Curated.txt",sep = "\t")
#---Above is done in 1.Curating Data; here it will be loaded to continue----
Disc.Genes <- read.delim("C:\\Users\\mohdb\\OneDrive - Towheed International School\\Documents\\NCA.ER\\Finalyzing\\1.Curating Data\\Disc.Curating.Data/3.Disc.Genes.Curated.txt",header = TRUE, as.is = TRUE, row.names = 1)
counts.Disc.genes <- as.matrix(Disc.Genes)
non_numeric_elements <- counts.Disc.genes[is.na(as.numeric(counts.Disc.genes))]
print(non_numeric_elements) #numeric(0)
#---2.Sample Annotation Data---------------------------------
Metabric_Manual_disc <- read.delim("C:\\Users\\mohdb\\OneDrive - Towheed International School\\Documents\\NCA.ER\\data/METABRIC Clinical.txt", row.names = 1)
View(Metabric_Manual_disc)

Metabric_Manual_disc <- Metabric_Manual_disc[order(rownames(Metabric_Manual_disc)),]
View(Metabric_Manual_disc)

# Replace dots with underscores in column names if needed
rownames(Metabric_Manual_disc) <- gsub("\\-", "_", rownames(Metabric_Manual_disc))
View(Metabric_Manual_disc)

# Load required library
library(dplyr)

# Assuming your dataset is named Metabric_Manual_disc
# Create a new binary dataframe based on transformations
Metabric_Manual_disc <- Metabric_Manual_disc %>%
  transmute(
    IDs = rownames(Metabric_Manual_disc),
    Cohort = Cohort,
    
    OS.time = Overall.Survival..Months.,
    OS.event = ifelse(Overall.Survival.Status == "1:DECEASED", 1, 0),
    DSS.event = ifelse(Patient.s.Vital.Status == "Died of Disease", 1, 0),
    RFS.event = ifelse(Relapse.Free.Status == "1:Recurred", 1, 0),
    RFS.time = Relapse.Free.Status..Months.,
    
    Grade = Neoplasm.Histologic.Grade,
    Size = Tumor.Size,
    LN = Lymph.nodes.examined.positive,
    Age = Age.at.Diagnosis,
    LN = Lymph.nodes.examined.positive,
    Age = Age.at.Diagnosis,
    
    
    # Subtypes for LumA, LumB, Basal, Her2, Normal based on Pam50 subtype
    LumA = ifelse(Pam50...Claudin.low.subtype == "LumA", 1, 0),
    LumB = ifelse(Pam50...Claudin.low.subtype == "LumB", 1, 0),
    Basal = ifelse(Pam50...Claudin.low.subtype == "Basal", 1, 0),
    Her2 = ifelse(Pam50...Claudin.low.subtype == "Her2", 1, 0),
    Normal = ifelse(Pam50...Claudin.low.subtype == "Normal", 1, 0),
    
    # ER and PR status (positive and negative)
    `ER+` = ifelse(ER.Status == "Positive", 1, 0),
    `ER-` = ifelse(ER.Status == "Negative", 1, 0),
    
    
    # Histologic Grade categories G1, G2, G3
    G1 = ifelse(Neoplasm.Histologic.Grade == 1, 1, 0),
    G2 = ifelse(Neoplasm.Histologic.Grade == 2, 1, 0),
    G3 = ifelse(Neoplasm.Histologic.Grade == 3, 1, 0),
    
    # Hormone Therapy (HT)
    HT = ifelse(Hormone.Therapy == "YES", 1, 0),
    
  )
View(Metabric_Manual_disc)

# Print head of the colAnnotation for verification
head(Metabric_Manual_disc)

#---Converting Metabric Clinical data into Sample Annotation---
Disc_sample_all_ids <- colnames(counts.Disc.genes)
View(Disc_sample_all_ids)
length(Disc_sample_all_ids)
#[1] 993
sort(Disc_sample_all_ids)
View(Disc_sample_all_ids)
length(Disc_sample_all_ids)
#[1] 993

all.equal(Disc_sample_all_ids,colnames(counts.Disc.genes))
#[1] TRUE

# Find the common sample IDs between the two datasets
common_samples_Disc <- intersect(rownames(Metabric_Manual_disc), Disc_sample_all_ids)
common_samples_Disc <- sort(common_samples_Disc)
View(common_samples_Disc)
length(common_samples_Disc)
#[1] 988

# Subset and reorder both datasets to only include common samples
counts.Disc.genes <- counts.Disc.genes[, common_samples_Disc]
all.equal(colnames(counts.Disc.genes), common_samples_Disc)
#[1] TRUE
View(counts.Disc.genes)
dim(counts.Disc.genes)
#[1] 19219  988

sampleAnnotation.Disc <- Metabric_Manual_disc[common_samples_Disc, ]
dim(sampleAnnotation.Disc)
#[1] 988   22

# Verify that they are aligned
all.equal(colnames(counts.Disc.genes), rownames(sampleAnnotation.Disc))  # Should return TRUE
#[1] TRUE

# If the first is not TRUE, you can match up the samples/columns in
# counts with the samples/rows in sampleAnnotation.Disc like this (which is fine
# to run even if the first was TRUE):

#tempindex <- match(colnames(counts.Disc.genes), rownames(sample_metadata))
#sampleAnnotation.Disc <- sample_metadata[tempindex, ]

#Check again

#all.equal(colnames(counts.Disc.genes), rownames(sampleAnnotation.Disc))

#--- 3. Gene Annotation Data----
#Duplicates Removed Manually
Disc.Gene.Annot.Curated <- read.delim("C:\\Users\\mohdb\\OneDrive - Towheed International School\\Documents\\NCA.ER\\Annotation/Gene_Annot_All_Used_ENS113.txt")
table(duplicated(Disc.Gene.Annot.Curated$external_gene_name))
#FALSE 
#18515 

dim(Disc.Gene.Annot.Curated)
#[1] 18515     7
length(Disc.Gene.Annot.Curated)
View(Disc.Gene.Annot.Curated)
dim(counts.Disc.genes)
#[1] 19219  988

#Only Validated IDs Taken
# Extract gene IDs from both datasets
ensemble_ids <- Disc.Gene.Annot.Curated$external_gene_name
View(ensemble_ids)
ensemble_ids <- sort(ensemble_ids)
length(ensemble_ids) #[1] 18515

Disc_all_gene_ids <- rownames(counts.Disc.genes )
View(Disc_all_gene_ids)

common_genes <- intersect(Disc.Gene.Annot.Curated$external_gene_name, Disc_all_gene_ids) # Should match or be close to 1420
common_genes <- sort(common_genes)
all.equal(ensemble_ids, common_genes) #[1] TRUE
length(common_genes)
#[1] 18515

# Subset and reorder both datasets to only include common samples
counts.Disc.genes <- counts.Disc.genes[common_genes, ]
View(counts.Disc.genes)

all.equal(ensemble_ids,rownames(counts.Disc.genes))
all.equal(common_genes,rownames(counts.Disc.genes))


#Set rownames for gene annotation, %n% doesn't let your data go NA 
Disc.Gene.Annot.Curated <- Disc.Gene.Annot.Curated[Disc.Gene.Annot.Curated$external_gene_name %in% common_genes, ]
rownames(Disc.Gene.Annot.Curated) <- Disc.Gene.Annot.Curated$external_gene_name
Disc.Gene.Annot.Curated <- Disc.Gene.Annot.Curated[order(rownames(Disc.Gene.Annot.Curated)),]
View(Disc.Gene.Annot.Curated)
all.equal(Disc.Gene.Annot.Curated$external_gene_name, rownames(counts.Disc.genes))
all.equal(Disc.Gene.Annot.Curated$external_gene_name,rownames(Disc.Gene.Annot.Curated))
all.equal(Disc.Gene.Annot.Curated$external_gene_name,common_genes)
all.equal(common_genes,ensemble_ids)

#Sort counts.Disc.genes and gene_annot_all_Disc based on rownames which is gene names
counts.Disc.genes <- counts.Disc.genes[order(rownames(counts.Disc.genes)), ]
Disc.Gene.Annot.Curated <- Disc.Gene.Annot.Curated[order(rownames(Disc.Gene.Annot.Curated)), ]

all.equal(rownames(counts.Disc.genes), rownames(Disc.Gene.Annot.Curated))  #[1] TRUE 
dim(counts.Disc.genes) #[1] 18515   988
dim(Disc.Gene.Annot.Curated) #[1] 18515     7
View(Disc.Gene.Annot.Curated)

#In rtni analysis it needs "SYMBOL" in rowAnnotation
colnames(Disc.Gene.Annot.Curated)
colnames(Disc.Gene.Annot.Curated)[colnames(Disc.Gene.Annot.Curated) == "external_gene_name"] <- "SYMBOL"
colnames(Disc.Gene.Annot.Curated)[colnames(Disc.Gene.Annot.Curated) == "ensembl_gene_id"] <- "ENSEMBL"

# Verify the change
colnames(Disc.Gene.Annot.Curated)


# One final check:
stopifnot(rownames(Disc.Gene.Annot.Curated) == rownames(counts.Disc.genes), # features
          rownames(sampleAnnotation.Disc) == colnames(counts.Disc.genes)) # samples

save(list = ls(), file = "1.Data before List.RData")

#---4. Prepare for rtni---
Disc_all <- list(
  expData = counts.Disc.genes,                       # expData as an assay
  rowAnnotation = Disc.Gene.Annot.Curated,             # Gene annotations (row metadata)
  colAnnotation = sampleAnnotation.Disc           # Sample annotations (column metadata)
)

View(Disc_all)

#---Load RegulatoryElements---
library(RTN)
# Load TF annotation
data("tfsData")

# Check TF annotation:
# Intersect TFs from Lambert et al. (2018) with gene annotation
# from the gene expression of 90 metabolic pathway cohort
regulatoryElements <- intersect(tfsData$Lambert2018$SYMBOL, Disc_all$rowAnnotation$SYMBOL)
View(regulatoryElements)
regulatoryElements <- sort(regulatoryElements)
View(regulatoryElements)

#---constructing rtni---
#This dataset consists of a list with 3 objects:
##a named gene expression matrix (tniData$expData),
##a data frame with gene annotations (tniData$rowAnnotation),
##and a data frame with sample annotations (tniData$colAnnotation).
##alternatively, 'expData' can be a 'SummarizedExperiment' object
rtni_disc_all <- tni.constructor(expData = Disc_all$expData,
                                   regulatoryElements = regulatoryElements,
                                   rowAnnotation = Disc_all$rowAnnotation,
                                   colAnnotation = Disc_all$colAnnotation)

#-Preprocessing for input data...
#--Mapping 'expData' to 'rowAnnotation'...
#--Mapping 'expData' to 'colAnnotation'...
#--Checking 'regulatoryElements' in 'rowAnnotation'...
#-Checking 'expData'...
#-Preprocessing complete!
  
save(rtni_disc_all, file = "2.rtni_disc_all.constructed.RData")
save.image("2.Disc.Before.Run.RData")
#--- 4. Check Normal Distribution-----
hist(counts.Disc.genes, breaks = 30, main = "Histogram of Data", xlab = "Values")

save(list = ls(), file = "3.rtni_hist_before_permutation.RData")
#---5.Building Regulons-----
rtni_Disc_all <- tni.permutation(rtni_disc_all, 
                                 nPermutations = 1000) 
#-Performing permutation analysis...
#--For 1508 regulons...
#|=================================================================================================| 100%
#-Permutation analysis complete! 
save(rtni_Disc_all,file = "5.rtni_Disc_all.RData")

#Unstable interactions are subsequently removed by bootstrap analysis,
##creates a consensus bootstrap network, referred here as refnet (reference network).
rtni_Disc_all <- tni.bootstrap(rtni_Disc_all)
#-Performing bootstrap analysis...
#--For 1508 regulons...
#|=================================================================================================| 100%
#-Bootstrap analysis complete! 
save(rtni_Disc_all,file = "6.rtni_Disc_all.Afterbootstrap.RData")

# Compute the DPI-filtered regulatory network
rtni_Disc_all.NA <- tni.dpi.filter(rtni_Disc_all, eps = NA)
#-Applying dpi filter...
#-DPI filter complete! 
tni.regulon.summary(rtni_Disc_all.NA)
#Regulatory network comprised of 1508 regulons. 
#-- DPI-filtered network: 
#  regulatoryElements            Targets              Edges 
#1508              17372             132674 
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0      33      67      88     117     802 
#-- Reference network: 
#  regulatoryElements            Targets              Edges 
#1508              17372            4280823 
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0     788    2655    2839    4624    7974 
#---

tni.regulon.summary(rtni_Disc_all.NA, regulatoryElements = "GATA3")
#The GATA3 regulon has 338 targets, it's a large and unbalanced regulon. 
#-- DPI filtered network targets:
#  Total Positive Negative 
#338       41      297 
#-- Reference network targets:
#  Total Positive Negative 
#7453     3560     3893 
#-- Regulators with mutual information:
#  AKAP8, AKAP8L, AKNA, AR, ARID2, ARID3B, ARID5B, ARNT2, ASCL1, ASCL2...[540 more]
#
#---
save(list = ls(), file = "6.rtni_Disc_all.NA.RData")

regulon.NA <- tni.get(rtni_Disc_all.NA, what = 'regulons.and.mode', idkey = "SYMBOL")
View(regulon.NA)
head(regulon.NA)
#---Extract Regulons----
# Find the maximum number of genes across all regulons
max_genes <- max(sapply(regulon.NA, function(x) if (is.null(x)) 0 else length(x)))

# Create a list of data frames with equal row lengths
regulon_list <- lapply(names(regulon.NA), function(regulon_name) {
  regulon_data <- regulon.NA[[regulon_name]]
  
  # Handle empty or NULL regulons
  if (is.null(regulon_data) || length(regulon_data) == 0) {
    df <- data.frame(
      Gene = rep(NA, max_genes),
      Value = rep(NA, max_genes),
      stringsAsFactors = FALSE
    )
  } else if (is.vector(regulon_data)) {
    # Ensure matching lengths for Gene and Value
    genes <- names(regulon_data)
    values <- regulon_data
    if (length(genes) == 0) genes <- rep(NA, length(values))
    df <- data.frame(
      Gene = genes,
      Value = values,
      stringsAsFactors = FALSE
    )
  } else if (is.matrix(regulon_data) || is.data.frame(regulon_data)) {
    df <- data.frame(
      Gene = rownames(regulon_data),
      Value = regulon_data[, 1], # Assuming values are in the first column
      stringsAsFactors = FALSE
    )
  } else {
    df <- data.frame(
      Gene = rep(NA, max_genes),
      Value = rep(NA, max_genes),
      stringsAsFactors = FALSE
    )
  }
  
  # Extend to max_genes rows if needed
  if (nrow(df) < max_genes) {
    df <- rbind(df, data.frame(
      Gene = rep(NA, max_genes - nrow(df)),
      Value = rep(NA, max_genes - nrow(df))
    ))
  }
  
  # Rename columns with regulon names
  colnames(df) <- c(paste0(regulon_name, "_Gene"), paste0(regulon_name, "_Value"))
  return(df)
})

# Combine all into one data frame
regulon_df <- do.call(cbind, regulon_list)

# Write to file
write.table(regulon_df, file = ("Disc.All.regulon.NA.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
save(regulon.NA,file = "Disc.All.Genes.Regulon.NA.RData")

#---9. Regulon activity profiles----------

library(Fletcher2013b)
library(pheatmap)
library(grid)
library(gridExtra)

# Load 'rtni1st' data object, which includes regulons and expression profiles
#here I will use: rtni_Disc_all.NA.NUM

# A list of transcription factors of interest ( #9 )
# Compute regulon activity for individual samples
rtni_Disc_all.NA <- tni.gsea2(rtni_Disc_all.NA, regulatoryElements = rtni_Disc_all.NA@regulatoryElements)
metabric_regact_disc <- tni.get(rtni_Disc_all.NA, what = "regulonActivity")
View(metabric_regact_disc)

# Get sample attributes from the 'rtni_Disc_all.NA.NUM' dataset
metabric_annot_disc <- tni.get(rtni_Disc_all.NA, "colAnnotation")

# Get ER+/- and PAM50 attributes for pheatmap
attribs_disc <- c("ER+","ER-")
metabric_annot_disc <- metabric_annot_disc[,attribs_disc]

save(rtni_Disc_all.NA, file = "rtni_Disc_all.NA.gsea.RData")

# Step 1: Identify samples with ER+ = 1
ER_positive_samples <- rownames(metabric_annot_disc)[metabric_annot_disc$`ER+` == 1]
ER_negative_samples <- rownames(metabric_annot_disc)[metabric_annot_disc$`ER+` == 0]

# Step 2: Order columns in metabric_regact_disc with ER+ samples on the left
ordered_sample_names <- c(ER_positive_samples, ER_negative_samples)
metabric_regact_disc_ordered <- metabric_regact_disc$differential[ordered_sample_names,]

# Define custom colors for each category in annotation_col
ER_annotation <- metabric_annot_disc[,c("ER+","ER-")]
ER_annotation_colors <- list(
  "ER+" = c("0" = "lightgrey", "1" = "blue"),
  "ER-" = c("0" = "lightgrey", "1" = "red")
)

pdf("Disc.all.Heatmap.pdf", width = 10, height = 10)

# Plot regulon activity profiles
disc.heatmap <- pheatmap(t(metabric_regact_disc_ordered), 
                         main="Discovery Set (n=988 samples)",
                         annotation_col = ER_annotation,
                         show_colnames = FALSE, annotation_legend = FALSE, 
                         clustering_method = "ward.D2", fontsize_row = 3,
                         clustering_distance_rows = "correlation",
                         cluster_cols = FALSE,
                         legend = TRUE,
                         annotation_colors = ER_annotation_colors,
                         fontsize_col = 3, fontsize = 6,
                         border_color = NA,
)

grid.text("Regulons", x= 0.97 , y=0.3, rot=270)

dev.off()

disc.session.info <- sessionInfo()
save(list = ls(),file = "Disc.all.RData")
