# ------------------------------------------------------------
# Load Required Packages
# ------------------------------------------------------------
library(survival)      # For survival analysis functions
library(survminer)     # For creating survival plots
library(cowplot)       # For combining multiple plots
library(gridExtra)     # For advanced grid layouts
library(grid)          # For low-level grid functionalities
library(openxlsx)      # For reading and writing Excel files

# ------------------------------------------------------------
# Read in the Survival Data
# ------------------------------------------------------------
# This file contains the survival information for the patients.
# NotD: Make sure the file path is correct.
  survival_data <- read.delim("E:/2.ER/5.Survival Analysis/PDS_KM/ER.Status.Survival/Disc.ER.Surv.txt", 
                            check.names = FALSE)
colnames(survival_data) <- gsub(" ", "_", colnames(survival_data))

# ------------------------------------------------------------
# Prepare Excel Workbook for Saving Results
# ------------------------------------------------------------
wb_all <- createWorkbook()
sheet_counter <- 1

# ------------------------------------------------------------
# Helper Function: Adjust Data for 5-Year Survival Analysis
# ------------------------------------------------------------
# For the 5-year analysis, if a patient's time is greater than 60 months,
# we consider them as censored (set event to 0)
adjust_to_5_year <- function(data, columns) {
  data_5yr <- data
  data_5yr[[columns$event]] <- ifelse(data_5yr[[columns$time]] > 60, 0, data_5yr[[columns$event]])
  return(data_5yr)
}

# Helper function to format p-value in a more readable way
format_p_value <- function(p) {
  if (p < 0.0001) {
    return("p < 0.0001")
  } else {
    return(paste0("p = ", round(p, 3)))
  }
}

# ------------------------------------------------------------
# Define Endpoints and Their Corresponding Columns
# ------------------------------------------------------------
endpoints <- list(
  OS  = list(time = 'OS_(Months)',  event = 'OS_Status'),
  DSS = list(time = 'OS_(Months)',  event = "DSS"),
  RFS = list(time = 'RFS_(Months)', event = "RFS"),
  MFS = list(time = 'OS_(Months)',  event = "MFS")
)

# ------------------------------------------------------------
# Define Endpoints and Their Full Names
# ------------------------------------------------------------
endpoints_full_name <- list(
  OS  = "Overall survival",
  DSS = "Disease-specific survival",
  RFS = "Relapse-free survival",
  MFS = "Metastasis-free survival"
)

# ------------------------------------------------------------
# Open PDF Device to Save All Plots in One File
# ------------------------------------------------------------
pdf(paste0("Disc.ER_", Sys.Date(), ".pdf"))
plot_created <- FALSE
# ------------------------------------------------------------
# Loop Through Each Endpoint to Create and Save Survival Curves
# ------------------------------------------------------------
for (endpoint in names(endpoints)) {
  cols <- endpoints[[endpoint]]
  
  # Filter rows with complete information in the time, event, and grouping ("ER_Status") columns
  current_data <- survival_data[complete.cases(survival_data[[cols$time]], 
                                               survival_data[[cols$event]],
                                               survival_data[["ER_Status"]]), ]
  
  # Create a patient identifier; here we assume "DISCOVERY" contains the patient IDs
  current_data$PatientID <- current_data[["DISCOVERY"]]
  
  # Convert the time column to numeric and round it to 4 decimals
  current_data[[cols$time]] <- round(as.numeric(current_data[[cols$time]]), 4)
  
  # Create the survival object and fit the Kaplan–Meier survival curve using ER_Status as grouping
  surv_obj <- Surv(time = current_data[[cols$time]], event = as.numeric(current_data[[cols$event]]))
  surv_fit <- survfit(surv_obj ~ `ER_Status`, data = current_data, na.action = na.exclude)
  
  # Calculate HR for the regular analysis using Cox regression
  #cox_model <- coxph(surv_obj ~ `ER_Status`, data = current_data)
  #hr_val <- summary(cox_model)$coefficients[1, "exp(coef)"]
  #p_val <- summary(cox_model)$coefficients[1, "Pr(>|z|)"]
  #hr_text <- paste0("HR: ", round(hr_val, 2), " (", format_p_value(p_val), ")")
  
  # Generate the survival plot for the regular analysis
  plot_regular <- ggsurvplot(surv_fit, data = current_data, pval = TRUE,
                             pval.coord = c(0, 0.15), pval.size = 4.8,
                             pval.method = TRUE, pval.method.size = 4.8,
                             pval.method.coord = c(0, 0.22),
                             xlab = "Time (months)", 
                             ylab = endpoints_full_name[[endpoint]],  # Use the full name here
                             legend.labs = c("Negative", "Positive"),  #by inactivating first is Negative then Positive
                             risk.table = TRUE, conf.int = FALSE,
                             risk.table.y.text = FALSE,
                             ggtheme = theme_minimal() +
                               theme(panel.grid = element_blank(),
                                     axis.line = element_line(linewidth = 0.5, color = "black"),
                                     axis.ticks = element_line(linewidth = 0.5, color = "black"),
                                     text         = element_text(size = 16),
                                     legend.text  = element_text(size = 14),
                                     
                                     # Increase the size of the months (X-axis) and probability numbers (Y-axis)
                                     axis.text.x = element_text(size = 14),  # Bigger months
                                     axis.text.y = element_text(size = 14),  # Bigger probability numbers
                                     legend.title = element_text(size = 14),
                                     # -------------------------------------------------------------
                                     # Add a bounding box around the entire plotting region:
                                     # -------------------------------------------------------------
                                     panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)
                                     )) 
  # Modify the ggplot to update the legend title
  plot_regular$plot <- plot_regular$plot + labs(color = "ER status") 
  # Post-process the risk table object
  plot_regular$table <- plot_regular$table +
    labs(y = NULL) +                             # Removes the y-axis label
    theme(axis.title.y = element_blank(),         # Also ensure the title is blank
          plot.title = element_text(size = 16)
    )
  # Annotate the plot with the HR information and add a title
  #plot_regular$plot <- plot_regular$plot +
   # annotate("text", x = 0, y = 0.08, label = hr_text, size = 4.8, hjust = 0) +
   # ggtitle(paste(endpoint,"\n", "Survival by ER_Status"))
  
  # Print the regular analysis plot to the PDF
  print(plot_regular)
  
  # --------------------------
  # 5-Year Analysis
  # --------------------------
  # Adjust the data for 5-year survival analysis
  current_data_5yr <- adjust_to_5_year(current_data, cols)
  current_data_5yr[[cols$time]] <- round(as.numeric(current_data_5yr[[cols$time]]), 4)
  
  # Create the survival object and fit the Kaplan–Meier curve for the 5-year analysis
  surv_obj_5yr <- Surv(time = current_data_5yr[[cols$time]], event = as.numeric(current_data_5yr[[cols$event]]))
  surv_fit_5yr <- survfit(surv_obj_5yr ~ `ER_Status`, data = current_data_5yr, na.action = na.exclude)
  
  # Calculate the HR for the 5-year analysis using Cox regression
 # cox_model_5yr <- coxph(surv_obj_5yr ~ `ER_Status`, data = current_data_5yr)
 # hr_val_5yr <- summary(cox_model_5yr)$coefficients[1, "exp(coef)"]
 # p_val_5yr <- summary(cox_model_5yr)$coefficients[1, "Pr(>|z|)"]
 # hr_text_5yr <- paste0("HR: ", round(hr_val_5yr, 2), " (", format_p_value(p_val_5yr), ")")
  
  # Generate the survival plot for the 5-year analysis
  plot_5yr <- ggsurvplot(surv_fit_5yr, data = current_data_5yr, pval = TRUE,
                         pval.coord = c(0, 0.15), pval.size = 4.8,
                         pval.method = TRUE, pval.method.size = 4.8,
                         pval.method.coord = c(0, 0.22),
                         xlab = "Time (months)", 
                         ylab =paste0(endpoints_full_name[[endpoint]]," (5-year)"),  # Use the full name here
                         legend.labs = c("Negative", "Positive"),  #by inactivating first is Negative then Positive
                         risk.table = TRUE, conf.int = FALSE,
                         risk.table.y.text = FALSE,
                         ggtheme = theme_minimal() +
                           theme(panel.grid = element_blank(),
                                 axis.line = element_line(linewidth = 0.5, color = "black"),
                                 axis.ticks = element_line(linewidth = 0.5, color = "black"),
                                 text         = element_text(size = 16),
                                 legend.text  = element_text(size = 14),
                                 
                                 # Increase the size of the months (X-axis) and probability numbers (Y-axis)
                                 axis.text.x = element_text(size = 14),  # Bigger months
                                 axis.text.y = element_text(size = 14),  # Bigger probability numbers
                                 legend.title = element_text(size = 14),
                                 # -------------------------------------------------------------
                                 # Add a bounding box around the entire plotting region:
                                 # -------------------------------------------------------------
                                 panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
                           ))
  # Modify the ggplot to update the legend title
  plot_5yr$plot <- plot_5yr$plot + labs(color = "ER status")
  # Post-process the risk table object
  plot_5yr$table <- plot_5yr$table +
    labs(y = NULL) +                             # Removes the y-axis label
    theme(axis.title.y = element_blank(),         # Also ensure the title is blank
          plot.title = element_text(size = 16))
          #plot_5yr$plot <- plot_5yr$plot +
  #  annotate("text", x = 0, y = 0.08, label = hr_text_5yr, size = 4.8, hjust = 0) +
   # ggtitle(paste(endpoint,"\n", "5-Year Survival by ER_Status"))
  
  # Print the 5-year analysis plot to the PDF
  print(plot_5yr)
  
  # ------------------------------------------------------------
  # Save Data Used for Each Plot into the Excel Workbook
  # ------------------------------------------------------------
  # Define the columns that were used for the plot
  columns_used <- c(cols$time, cols$event, "ER_Status", "PatientID")
  
  # Subset the data to only include these columns
  current_data_subset <- current_data[, columns_used, drop = FALSE]
  
  # For regular analysis
  sheet_name_regular <- endpoint
  addWorksheet(wb_all, sheetName = sheet_name_regular)
  writeData(wb_all, sheet = sheet_name_regular, x = current_data_subset)
  
  # Subset the data to only include these columns
  current_data_5yr_subset <- current_data_5yr[, columns_used, drop = FALSE]
  
  # For 5-year analysis
  sheet_name_5yr <- paste(endpoint, "5yr", sep = "_")
  addWorksheet(wb_all, sheetName = sheet_name_5yr)
  writeData(wb_all, sheet = sheet_name_5yr, x = current_data_5yr_subset)
  
  sheet_counter <- sheet_counter + 2
}

# Close the PDF device
dev.off()

# ------------------------------------------------------------
# Save the Excel Workbook with All Data Sheets
# ------------------------------------------------------------
saveWorkbook(wb_all, file = paste0("Disc.ER_Survival_Data_", Sys.Date(), ".xlsx"), overwrite = TRUE)
se <- sessionInfo()
save(list = ls(), file = paste0('Disc.ER_',Sys.Date(),'.RData'))

rm(list = ls())

# ------------------------------------------------------------
# Read in the Survival Data
# ------------------------------------------------------------
# This file contains the survival information for the patients.
# NotD: Make sure the file path is correct.
survival_data <- read.delim("E:/2.ER/5.Survival Analysis/PDS_KM/ER.Status.Survival/Valid.ER.Surv.txt", 
                            check.names = FALSE)
colnames(survival_data) <- gsub(" ", "_", colnames(survival_data))

# ------------------------------------------------------------
# Prepare Excel Workbook for Saving Results
# ------------------------------------------------------------
wb_all <- createWorkbook()
sheet_counter <- 1

# ------------------------------------------------------------
# Helper Function: Adjust Data for 5-Year Survival Analysis
# ------------------------------------------------------------
# For the 5-year analysis, if a patient's time is greater than 60 months,
# we consider them as censored (set event to 0)
adjust_to_5_year <- function(data, columns) {
  data_5yr <- data
  data_5yr[[columns$event]] <- ifelse(data_5yr[[columns$time]] > 60, 0, data_5yr[[columns$event]])
  return(data_5yr)
}

# Helper function to format p-value in a more readable way
format_p_value <- function(p) {
  if (p < 0.001) {
    return("p < 0.001")
  } else {
    return(paste0("p = ", round(p, 3)))
  }
}

# ------------------------------------------------------------
# Define Endpoints and Their Corresponding Columns
# ------------------------------------------------------------
endpoints <- list(
  OS  = list(time = 'OS_(Months)',  event = 'OS_Status'),
  DSS = list(time = 'OS_(Months)',  event = "DSS"),
  RFS = list(time = 'RFS_(Months)', event = "RFS"),
  MFS = list(time = 'OS_(Months)',  event = "MFS")
)

# ------------------------------------------------------------
# Define Endpoints and Their Full Names
# ------------------------------------------------------------
endpoints_full_name <- list(
  OS  = "Overall survival",
  DSS = "Disease-specific survival",
  RFS = "Relapse-free survival",
  MFS = "Metastasis-free survival"
)

# ------------------------------------------------------------
# Open PDF Device to Save All Plots in One File
# ------------------------------------------------------------
pdf(paste0("Valid.ER_", Sys.Date(), ".pdf"))
plot_created <- FALSE

# ------------------------------------------------------------
# Loop Through Each Endpoint to Create and Save Survival Curves
# ------------------------------------------------------------
for (endpoint in names(endpoints)) {
  cols <- endpoints[[endpoint]]
  
  # Filter rows with complete information in the time, event, and grouping ("ER_Status") columns
  current_data <- survival_data[complete.cases(survival_data[[cols$time]], 
                                               survival_data[[cols$event]],
                                               survival_data[["ER_Status"]]), ]
  
  # Create a patient identifier; here we assume "Validation" contains the patient IDs
  current_data$PatientID <- current_data[["Validation"]]
  
  # Convert the time column to numeric and round it to 4 decimals
  current_data[[cols$time]] <- round(as.numeric(current_data[[cols$time]]), 4)
  
  # Create the survival object and fit the Kaplan–Meier survival curve using ER_Status as grouping
  surv_obj <- Surv(time = current_data[[cols$time]], event = as.numeric(current_data[[cols$event]]))
  surv_fit <- survfit(surv_obj ~ `ER_Status`, data = current_data, na.action = na.exclude)
  
  # Calculate HR for the regular analysis using Cox regression
  #cox_model <- coxph(surv_obj ~ `ER_Status`, data = current_data)
  #hr_val <- summary(cox_model)$coefficients[1, "exp(coef)"]
  #p_val <- summary(cox_model)$coefficients[1, "Pr(>|z|)"]
  #hr_text <- paste0("HR: ", round(hr_val, 2), " (", format_p_value(p_val), ")")
  
  # Generate the survival plot for the regular analysis
  plot_regular <- ggsurvplot(surv_fit, data = current_data, pval = TRUE,
                             pval.coord = c(0, 0.15), pval.size = 4.8,
                             pval.method = TRUE, pval.method.size = 4.8,
                             pval.method.coord = c(0, 0.22),
                             xlab = "Time (months)", 
                             ylab =endpoints_full_name[[endpoint]],  # Use the full name here
                             legend.labs = c("Negative", "Positive"),  #by inactivating first is Negative then Positive
                             risk.table = TRUE, conf.int = FALSE,
                             risk.table.y.text = FALSE,
                             ggtheme = theme_minimal() +
                               theme(panel.grid = element_blank(),
                                     axis.line = element_line(linewidth = 0.5, color = "black"),
                                     axis.ticks = element_line(linewidth = 0.5, color = "black"),
                                     text         = element_text(size = 16),
                                     legend.text  = element_text(size = 14),
                                     
                                     # Increase the size of the months (X-axis) and probability numbers (Y-axis)
                                     axis.text.x = element_text(size = 14),  # Bigger months
                                     axis.text.y = element_text(size = 14),  # Bigger probability numbers
                                     
                                     # -------------------------------------------------------------
                                     # Add a bounding box around the entire plotting region:
                                     # -------------------------------------------------------------
                                     panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
                               ))
  # Modify the ggplot to update the legend title
  plot_regular$plot <- plot_regular$plot + labs(color = "ER status")
  # Post-process the risk table object
  plot_regular$table <- plot_regular$table +
    labs(y = NULL) +                             # Removes the y-axis label
    theme(axis.title.y = element_blank(),         # Also ensure the title is blank
          plot.title = element_text(size = 16))
  # Annotate the plot with the HR information and add a title
  #plot_regular$plot <- plot_regular$plot +
  # annotate("text", x = 0, y = 0.08, label = hr_text, size = 4.8, hjust = 0) +
  # ggtitle(paste(endpoint,"\n", "Survival by ER_Status"))
  
  # Print the regular analysis plot to the PDF
  print(plot_regular)
  
  # --------------------------
  # 5-Year Analysis
  # --------------------------
  # Adjust the data for 5-year survival analysis
  current_data_5yr <- adjust_to_5_year(current_data, cols)
  current_data_5yr[[cols$time]] <- round(as.numeric(current_data_5yr[[cols$time]]), 4)
  
  # Create the survival object and fit the Kaplan–Meier curve for the 5-year analysis
  surv_obj_5yr <- Surv(time = current_data_5yr[[cols$time]], event = as.numeric(current_data_5yr[[cols$event]]))
  surv_fit_5yr <- survfit(surv_obj_5yr ~ `ER_Status`, data = current_data_5yr, na.action = na.exclude)
  
  # Calculate the HR for the 5-year analysis using Cox regression
  # cox_model_5yr <- coxph(surv_obj_5yr ~ `ER_Status`, data = current_data_5yr)
  # hr_val_5yr <- summary(cox_model_5yr)$coefficients[1, "exp(coef)"]
  # p_val_5yr <- summary(cox_model_5yr)$coefficients[1, "Pr(>|z|)"]
  # hr_text_5yr <- paste0("HR: ", round(hr_val_5yr, 2), " (", format_p_value(p_val_5yr), ")")
  
  # Generate the survival plot for the 5-year analysis
  plot_5yr <- ggsurvplot(surv_fit_5yr, data = current_data_5yr, pval = TRUE,
                         pval.coord = c(0, 0.15), pval.size = 4.8,
                         pval.method = TRUE, pval.method.size = 4.8,
                         pval.method.coord = c(0, 0.22),
                         xlab = "Time (months)", 
                         ylab =paste0(endpoints_full_name[[endpoint]]," (5-year)"),  # Use the full name here
                         legend.labs = c("Negative", "Positive"),  #by inactivating first is Negative then Positive
                         risk.table = TRUE, conf.int = FALSE,
                         risk.table.y.text = FALSE,
                         ggtheme = theme_minimal() +
                           theme(panel.grid = element_blank(),
                                 axis.line = element_line(linewidth = 0.5, color = "black"),
                                 axis.ticks = element_line(linewidth = 0.5, color = "black"),
                                 text         = element_text(size = 16),
                                 legend.text  = element_text(size = 14),
                                 
                                 # Increase the size of the months (X-axis) and probability numbers (Y-axis)
                                 axis.text.x = element_text(size = 14),  # Bigger months
                                 axis.text.y = element_text(size = 14),  # Bigger probability numbers
                                 
                                 # -------------------------------------------------------------
                                 # Add a bounding box around the entire plotting region:
                                 # -------------------------------------------------------------
                                 panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
                           ))
  # Modify the ggplot to update the legend title
  plot_5yr$plot <- plot_5yr$plot + labs(color = "ER status")
  # Post-process the risk table object
  plot_5yr$table <- plot_5yr$table +
    labs(y = NULL) +                             # Removes the y-axis label
    theme(axis.title.y = element_blank(),         # Also ensure the title is blank
          plot.title = element_text(size = 16))
  #plot_5yr$plot <- plot_5yr$plot +
  #  annotate("text", x = 0, y = 0.08, label = hr_text_5yr, size = 4.8, hjust = 0) +
  # ggtitle(paste(endpoint,"\n", "5-Year Survival by ER_Status"))
  
  # Print the 5-year analysis plot to the PDF
  print(plot_5yr)
  
  # ------------------------------------------------------------
  # Save Data Used for Each Plot into the Excel Workbook
  # ------------------------------------------------------------
  # Define the columns that were used for the plot
  columns_used <- c(cols$time, cols$event, "ER_Status", "Patient_ID")
  
  # Subset the data to only include these columns
  current_data_subset <- current_data[, columns_used, drop = FALSE]
  
  # For regular analysis
  sheet_name_regular <- endpoint
  addWorksheet(wb_all, sheetName = sheet_name_regular)
  writeData(wb_all, sheet = sheet_name_regular, x = current_data_subset)
  
  # Subset the data to only include these columns
  current_data_5yr_subset <- current_data_5yr[, columns_used, drop = FALSE]
  
  # For 5-year analysis
  sheet_name_5yr <- paste(endpoint, "5yr", sep = "_")
  addWorksheet(wb_all, sheetName = sheet_name_5yr)
  writeData(wb_all, sheet = sheet_name_5yr, x = current_data_5yr_subset)
  
  sheet_counter <- sheet_counter + 2
}


# Close the PDF device
dev.off()

# ------------------------------------------------------------
# Save the Excel Workbook with All Data Sheets
# ------------------------------------------------------------
saveWorkbook(wb_all, file = paste0("Valid.ER_Survival_Data_", Sys.Date(), ".xlsx"), overwrite = TRUE)
se <- sessionInfo()
save(list = ls(),file = paste0('Valid.ER_',Sys.Date(),'.RData'))
