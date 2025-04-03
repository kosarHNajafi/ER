#Conditions-----------------
#Error: NOTE: all names in 'phenotype' should be available in col1 of 'phenoIDs'!
library(dplyr)
gene_annot_mp_disc_tna_p <- gene_annot_mp.tf_disc %>%
  select(SYMBOL, everything())

gene_annot_mp_disc_tna_common_p <- intersect(gene_annot_mp_disc_tna_p$SYMBOL,rownames(phenotype_mp_pos))
gene_annot_mp_disc_tna_p <- gene_annot_mp_disc_tna_p[gene_annot_mp_disc_tna_common_p,]
View(gene_annot_mp_disc_tna_p)

phenotype_mp_pos <- phenotype_mp_pos[gene_annot_mp_disc_tna_common_p,]
View(phenotype_mp_pos)

all(gene_annot_mp_disc_tna_p$SYMBOL == gene_annot_mp_disc_tna_common_p)
all(rownames(phenotype_mp_pos) == rownames(gene_annot_mp_disc_tna_p))

write.table(phenotype_mp_pos,file = "C3.Disc.Phenotype.txt",sep = "\t")

# Extract 'logFC' as a named numeric vector
logFC_disc_mp_p <- setNames(phenotype_mp_pos$logFC, rownames(phenotype_mp_pos))
View(logFC_disc_mp_p)

all.equal(names(logFC_disc_mp_p),rownames(phenotype_mp_pos))
all.equal(as.vector(logFC_disc_mp_p),phenotype_mp_pos$logFC)

tna.disc_mp_p <- list(
  phenotype = logFC_disc_mp_p,
  phenoID = gene_annot_mp_disc_tna_p,
  hits = rownames(phenotype_mp_pos)
)

View(tna.disc_mp_p)
save(phenotype_mp_pos,logFC_disc_mp_p,file = "C3.Differential.Exp.Disc.MP.RData")
#---9.RTNA----
# Input 1: 'object', a TNI object with regulons
# Input 2: 'phenotype', a named numeric vector, usually log2 differential expression levels
# Input 3: 'hits', a character vector, usually a set of differentially expressed genes
# Input 4: 'phenoIDs', an optional data frame with gene anottation mapped to the phenotype

#rtni_disc_mp.tf.dpi after dpi.filter
#CHECK "rtnaData"
rtna_disc.mp_p <- tni2tna.preprocess(object = rtni_disc_mp.tf.NA, 
                                     phenotype = tna.disc_mp_p$phenotype, 
                                     hits = tna.disc_mp_p$hits, 
                                     phenoIDs = tna.disc_mp_p$phenoID)
#-Preprocessing for input data...
#--Mapping 'phenotype' to 'phenoIDs'...
#--Mapping 'hits' to 'phenoIDs'...
#-Mapping 'transcriptionalNetwork' annotation to 'phenotype'...
#--Checking agreement between 'transcriptionalNetwork' and 'phenotype'... 82.9% ! 
#--Extracting regulons...
#-Preprocessing complete!
#Warning message:
#NOTE: 17.1% of 'transcriptionalNetwork' targets not represented in the 'phenotype'! 

# Run the MRA method
rtna_disc.mp_p <- tna.mra(rtna_disc.mp_p)
#-Performing master regulatory analysis...
#--For 254 regulons...
#|==============================================================================| 100%
#Master regulatory analysis complete

# Get MRA results;
#..setting 'ntop = -1' will return all results, regardless of a threshold
mra_disc.mp_p <- tna.get(rtna_disc.mp_p, what="mra", ntop = -1)
View(mra_disc.mp_p)
write.table(mra_disc.mp_p,file = "C3.DISC.MRA.MP.txt",sep = "\t")

# Run the GSEA method
# Please set nPermutations >= 1000
rtna_disc.mp_p <- tna.gsea1(rtna_disc.mp_p, nPermutations=1000)
#-Performing gene set enrichment analysis...
#--For 182 regulons...
#|==============================================================================| 100%
#-Gene set enrichment analysis complete 

# Get GSEA results
gsea1_disc.mp_p <- tna.get(rtna_disc.mp_p, what="gsea1", ntop = -1)
head(gsea1_disc.mp_p)
#Regulon Regulon.Size Observed.Score     Pvalue Adjusted.Pvalue
#ENSG00000065978    YBX1           52           0.75 7.5264e-07      0.00013322
#ENSG00000091831    ESR1           44           0.74 3.4307e-06      0.00030362
#ENSG00000184221   OLIG1           30           0.79 2.7111e-05      0.00137300
#ENSG00000131668   BARX1           21           0.86 3.1029e-05      0.00137300
#ENSG00000107485   GATA3           45           0.71 5.9439e-05      0.00210410
#ENSG00000173894    CBX2           41           0.70 2.2908e-04      0.00675770

View(gsea1_disc.mp_p)
write.table(gsea1_disc.mp_p,file = "C3.DISC.tna.gsea1.mp.txt",sep = "\t",row.names = TRUE, col.names = TRUE)

# Filter for significant TFs
gsea1_disc.mp.sig_p <- gsea1_disc.mp_p[gsea1_disc.mp_p$Adjusted.Pvalue <= 0.05, ]
View(gsea1_disc.mp.sig_p)
write.table(gsea1_disc.mp.sig_p,file = "C3.DISC.sig.tna.gsea1.mp.txt",sep = "\t")

# Plot GSEA results
tna.plot.gsea1(rtna_disc.mp_p,
               labPheno="abs(log2 fold changes)", 
               #ntop = 5,
               tfs = c("GATA3","ESR1","YBX1","CBX2","BCL11A"),
               file = "Final.C3.disc.gsea1_top5_phenotype", 
               ylimPanels = c(0.0,3.5,0.0,2),
               regulon.order = "score"
)


# Run the GSEA-2T method
# Please set nPermutations >= 1000
rtna_disc.mp_p <- tna.gsea2(rtna_disc.mp_p, nPermutations = 1000)
#-Performing two-tailed GSEA analysis...
#--For 182 regulons...
#|==============================================================================| 100%
#|==============================================================================| 100%
#GSEA2 analysis complete 
View(rtna_disc.mp_p)

# Get GSEA-2T results
gsea2_disc.mp_p <- tna.get(rtna_disc.mp_p, what = "gsea2", ntop = -1)
head(gsea2_disc.mp_p$differential)
#Regulon Regulon.Size Observed.Score   Pvalue Adjusted.Pvalue
#ENSG00000131668   BARX1           21          -1.87 0.000999       0.0028987
#ENSG00000119866  BCL11A           41          -1.55 0.000999       0.0028987
#ENSG00000173894    CBX2           41          -1.68 0.000999       0.0028987
#ENSG00000172216   CEBPB           25          -1.47 0.000999       0.0028987
#ENSG00000115163   CENPA           50          -1.57 0.000999       0.0028987
#ENSG00000257923    CUX1           20          -1.63 0.000999       0.0028987

write.table(gsea2_disc.mp_p,file = "C3.DISC.tna.gsea2.mp.txt",sep = "\t",row.names = TRUE, col.names = TRUE)

# Plot GSEA-2T results
tna.plot.gsea2(rtna_disc.mp_p, labPheno="log2 fold changes", tfs="CBX2", file = "C3.Disc.gsea2_mp_phenotype")
save(list = ls(),file = "C3.Disc.MP.295tf.RData")
