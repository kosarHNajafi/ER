# Function to process dataset
process_data <- function(file_path, output_file) {
  # Read the data
  data <- read.delim(file_path, check.names = FALSE, stringsAsFactors = FALSE)
  
  # Renaming columns as per the specified mapping
  colnames(data) <- gsub("Age at Diagnosis", "Age (years)", colnames(data))
  colnames(data) <- gsub("Lymph nodes examined positive", "Lymph Node Status", colnames(data))
  colnames(data) <- gsub("Tumor Size", "Tumor Size (mm)", colnames(data))
  colnames(data) <- gsub("Radio Therapy", "Radiotherapy", colnames(data))
  colnames(data) <- gsub("Oncotree Code", "Histology", colnames(data))
  colnames(data) <- gsub("Pam50 \\+ Claudin-low subtype", "PAM50 status", colnames(data))
  
  # Create Receptor Status column
  data$`Receptor Status` <- ifelse(
    is.na(data$`ER Status`) | is.na(data$`PR Status`) | is.na(data$`HER2 Status`), 
    "NA",  # If any of the values are missing, assign "NA"
    ifelse(
      data$`ER Status` == "Negative" & data$`PR Status` == "Negative" & data$`HER2 Status` == "Negative", 
      "TNBC",  # If all three are negative, assign "TNBC"
      "non-TNBC"  # Otherwise, assign "non-TNBC"
    )
  )
  
  # Reclassify Cancer Type Detailed
  data$`Cancer Type Detailed` <- ifelse(
    data$`Cancer Type Detailed` %in% c("Breast Invasive Lobular Carcinoma", "Invasive Breast Carcinoma"),
    data$`Cancer Type Detailed`, "Others"
  )
  
  # Reclassify Histology
  data$Histology <- ifelse(
    is.na(data$Histology), "NA",
    ifelse(data$Histology %in% c("IDC", "ILC"), data$Histology, "Others")
  )
  
  # Reclassify Lymph Node Status
  data$`Lymph Node Status` <- ifelse(
    is.na(data$`Lymph Node Status`), "NA",
    ifelse(data$`Lymph Node Status` >= 1, "≥1", "0")
  )
  
  # Reclassify Tumor Size
  data$`Tumor Size (mm)` <- ifelse(
    is.na(data$`Tumor Size (mm)`), "NA",
    ifelse(data$`Tumor Size (mm)` < 20, "<20", "≥20")
  )
  
  # Reclassify Tumor Stage
  data$`Tumor Stage` <- ifelse(
    is.na(data$`Tumor Stage`), "NA", 
    ifelse(data$`Tumor Stage` %in% c("3", "4"), "3,4", data$`Tumor Stage`)
  )
  
  data$`3-Gene classifier subtype` <- ifelse(
    is.na(data$`3-Gene classifier subtype`), "NA", 
    data$`3-Gene classifier subtype`
  )
  
  data$`Cellularity` <- ifelse(
    is.na(data$`Cellularity`), "NA",
    ifelse(data$`Cellularity` %in% c("High", "Moderate"), data$`Cellularity`, "Low")
  )
  
  data$`Hormone Therapy` <- ifelse(
    is.na(data$`Hormone Therapy`), "NA",
    ifelse(data$`Hormone Therapy` %in% c("Yes", "No"), data$`Hormone Therapy`, "NA")
  )
  
  data$`Chemotherapy` <- ifelse(
    is.na(data$`Chemotherapy`), "NA",
    ifelse(data$`Chemotherapy` %in% c("YES", "NO"), data$`Chemotherapy`, "NA")
  )
  data$`PAM50 status` <- ifelse(
    is.na(data$`PAM50 status`), "NA",
    ifelse(data$`PAM50 status` %in% c("Basal","Her2","LumA", "LumB", "Normal", "claudin-low"), data$`PAM50 status`, "NA")
  )
  
  data$`PR Status` <- ifelse(
    is.na(data$`PR Status`), "NA",
    ifelse(data$`PR Status` %in% c("Positive", "Negative"), data$`PR Status`, "NA")
  )
  
  data$`HER2 Status` <- ifelse(
    is.na(data$`HER2 Status`), "NA",
    ifelse(data$`HER2 Status` %in% c("Positive", "Negative"), data$`HER2 Status`, "NA")
  )
  selected_columns <- c(
    "Age (years)", "Cancer Type Detailed", "Cellularity", "Chemotherapy", "ER Status", 
    "HER2 Status", "Histology", "Hormone Therapy","Lymph Node Status", 
    "PAM50 status", "PR Status", "Receptor Status", "Radiotherapy", 
    "Neoplasm Histologic Grade", "Tumor Size (mm)", "Tumor Stage", "3-Gene classifier subtype"
  )
  
  data <- data[,selected_columns]
  
  # Define the grouping column
  group_column <- "ER Status"  # Ensure this matches exactly with your dataset
  
  # Convert the grouping column to a factor
  data[[group_column]] <- as.factor(data[[group_column]])
  
  # Group "Age (years)" into <50 and ≥50
  data$`Age (years)` <- factor(ifelse(data$`Age (years)` < 50, "<50", "≥50"))
  
  # Calculate total sample sizes
  total_n <- nrow(data)
  er_counts <- table(data[[group_column]])
  
  # Extract counts for ER Positive & ER Negative
  er_positive_label <- names(er_counts)[1]
  er_negative_label <- names(er_counts)[2]
  er_positive_n <- er_counts[er_positive_label]
  er_negative_n <- er_counts[er_negative_label]
  
  # Function to calculate percentages
  calculate_percentages_matrix <- function(ct) {
    cpct <- matrix("", nrow = nrow(ct), ncol = ncol(ct))
    rownames(cpct) <- rownames(ct)
    colnames(cpct) <- colnames(ct)
    for (j in seq_len(ncol(ct))) {
      col_sum <- sum(ct[, j])
      for (i in seq_len(nrow(ct))) {
        cpct[i, j] <- paste0(ct[i, j], " (", round((ct[i, j] / col_sum) * 100, 1), ")")
      }
    }
    return(cpct)
  }
  
  # Function to apply tests and format output
  apply_tests <- function(df, group_col) {
    results_list <- list()
    
    for (col in colnames(df)) {
      if (col != group_col && col != "Patient ID") {
        df[[col]] <- as.factor(df[[col]])
        df_filtered <- droplevels(df[!is.na(df[[col]]), ])
        
        if (length(unique(df_filtered[[col]])) < 2) next
        
        # Create contingency table
        ct <- as.matrix(table(df_filtered[[col]], df_filtered[[group_col]]))
        cpct <- calculate_percentages_matrix(ct)
        
        # Capture p-value
        p_value <- NA
        error_message <- ""
        
        test_result <- withCallingHandlers(
          tryCatch({
            if (nrow(ct) == 2 && ncol(ct) == 2) {
              fisher.test(ct, workspace = 2e8)
            } else {
              chisq.test(ct)
            }
          }, error = function(e) {
            error_message <<- paste("Error:", e$message)
            return(NULL)
          }),
          warning = function(w) {
            error_message <<- paste("Warning:", w$message)
            invokeRestart("muffleWarning")
          }
        )
        
        if (!is.null(test_result)) {
          p_value <- format.pval(test_result$p.value, digits = 4, eps = 0.0001)
        }
        
        for (row_name in rownames(ct)) {
          total_count <- sum(ct[row_name, ])
          total_pct <- round(total_count / sum(ct) * 100, 1)
          total_str <- paste0(total_count, " (", total_pct, ")")
          
          # Extract values correctly
          er_positive_val <- ifelse(er_positive_label %in% colnames(ct), cpct[row_name, er_positive_label], "0 (0)")
          er_negative_val <- ifelse(er_negative_label %in% colnames(ct), cpct[row_name, er_negative_label], "0 (0)")
          
          results_list[[length(results_list) + 1]] <- data.frame(
            Variable = col,
            Category = row_name,
            Total = total_str,
            `ER.Positive` = er_positive_val,
            `ER.Negative` = er_negative_val,
            P_Value = ifelse(row_name == rownames(ct)[1], p_value, ""),
            Error_Message = ifelse(row_name == rownames(ct)[1], error_message, ""),
            stringsAsFactors = FALSE
          )
        }
      }
    }
    
    results_df <- do.call(rbind, results_list)
    
    # Formatting: Remove duplicate variable names
    results_df$Variable <- ifelse(duplicated(results_df$Variable), "", results_df$Variable)
    
    return(results_df)
  }
  
  # Apply the tests and save the results
  test_results <- apply_tests(data, group_column)
  
  # Insert the total sample size header row with new format
  header_row <- data.frame(
    Variable = "Total",
    Category = "",
    Total = paste0("Total (n = ", total_n, ")"),
    `ER.Positive` = paste0("Total (n = ", er_positive_n, ")"),
    `ER.Negative` = paste0("Total (n = ", er_negative_n, ")"),
    P_Value = "",
    Error_Message = "",
    stringsAsFactors = FALSE
  )
  
  # Merge the header row with the main table
  final_results <- rbind(header_row, test_results)
  
  # Rename the column headers
  colnames(final_results) <- c(
    "Variable", "Category", "Total;n(%)", "ER.Positive;n(%)", "ER.Negative;n(%)", "P_Value", "Error_Message"
  )
  
  # Save the formatted results
  write.table(final_results, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
  
  se <- sessionInfo()
  
  # Save session info and objects
  save(list = ls(), file = paste0(sub("\\.txt$", "", output_file), "_", Sys.Date(), ".RData"))
}

# Process the dataset
process_data("Validation.txt", "Final_Valid.ER.Clinicopathological.txt")
process_data("Discovery.txt", "Final_Disc.ER.Clinicopathological.txt")
