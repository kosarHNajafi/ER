#4.3Regulon activity profiles
#https://bioconductor.org/packages/devel/bioc/vignettes/RTN/inst/doc/RTN.html#ref-Margolin2006b

#---9. Regulon activity profiles----------
library(pheatmap)
library(grid)
library(gridExtra)

# Compute regulon activity for individual samples
rtni1st <- tni.gsea2(rtni_disc_mp.tf.NA, regulatoryElements = rtni_disc_mp.tf.NA@regulatoryElements)
#-Checking log space... OK!
#  -Performing two-tailed GSEA...
#--For 177 regulon(s) and 988 sample(s)...
#|==============================================================================| 100%
#-GSEA2 complete! 
  
metabric_regact <- tni.get(rtni1st, what = "regulonActivity")
View(metabric_regact)

# Get sample attributes from the 'rtni1st' dataset
metabric_annot <- tni.get(rtni1st, "colAnnotation")
# Create a new annotation data frame for ER status.
# Here we assume that the "ER+" column has a value of 1 for ER-positive samples.
er_status <- ifelse(metabric_annot$`ER+` == 1, "ERpos", "ERneg")
# Ensure that the factor levels are ordered so that ERpos comes before ERneg.
metabric_er_annot <- data.frame(ER = factor(er_status, levels = c("ERpos", "ERneg")))
rownames(metabric_er_annot) <- rownames(metabric_annot)

# Transpose the differential matrix so that samples become columns.
# In the transposed matrix, column names should correspond to sample names.
View(metabric_regact$differential)
mat <- t(metabric_regact$differential)
View(mat)
mat <- mat[order(rownames(mat)),order(colnames(mat))]

# Order the samples so that ERpos samples are on the left and ERneg on the right.
# We assume that the rownames of metabric_annot match the column names of mat.
stopifnot(all.equal(rownames(metabric_er_annot), colnames(mat)))

# Define custom annotation colors: blue for ERpos and red for ERneg.
annotation_colors <- list(ER = c("ERpos" = "steelblue", "ERneg" = "tomato"))

# Order metabric_er_annot to ensure ERpos samples are on the left and ERneg on the right
metabric_er_annot$ER <- factor(metabric_er_annot$ER, levels = c("ERpos", "ERneg"))

# Ensure that ERpos samples come first, followed by ERneg
ordered_samples <- rownames(metabric_er_annot)[order(metabric_er_annot$ER)]

stopifnot(all.equal(rownames(metabric_er_annot), colnames(mat)))

# Reorder the columns of mat based on the ER status from metabric_annot
mat <- mat[, ordered_samples]  # order samples based on ER status

# Open PDF device to save the heatmap.
pdf(paste0("Disc.MP.288tf_",Sys.Date(),".pdf"), width = 10, height = 10)

# Plot the heatmap.
# Set cluster_cols = FALSE to preserve the manual sample order.
disc.heatmap <- pheatmap(mat, 
                         main = "Discovery Set (n=988 samples)",
                         cluster_cols = FALSE,          # preserve the manually set order
                         clustering_method = "ward.D2",   # cluster rows (regulons) if desired
                         clustering_distance_rows = "correlation",
                         show_colnames = FALSE,
                         annotation_col = metabric_er_annot,
                         annotation_legend = FALSE,
                         annotation_colors = annotation_colors,
                         border_color = NA,               # remove cell borders (no black lines)
                         fontsize_row = 5)

# Optionally add a text label for "Regulons"
grid.text("Regulons", x = 0.97, y = 0.6, rot = 270)

# Close the PDF device.
dev.off()

write.table(mat,file = paste0("Disc.Regulon.for.Heatmap,",Sys.Date(),".txt"),sep = "\t")
se <- sessionInfo()
