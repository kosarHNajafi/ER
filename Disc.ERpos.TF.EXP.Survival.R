##############################################################################
# 1) Load Required Packages
##############################################################################
# These libraries provide survival analysis, plotting, and Excel-writing functions.
library(survival)      # Core survival analysis routines
library(survminer)     # Convenient functions for Kaplan-Meier plots, etc.
library(cowplot)       # Useful for combining multiple plots (not mandatory here)
library(gridExtra)     # Advanced grid layout utilities
library(grid)          # Low-level grid capabilities
library(openxlsx)      # For reading/writing Excel files

##############################################################################
# 2) Load set working directory
##############################################################################
setwd("E:/2.ER/7.NCA.ER.Final/3.Gene Expression/3.Survival Analysis/Disc.ERpos/")

##############################################################################
# 3) Read in the Survival Data
##############################################################################
# This file contains survival times, event indicators, and TF columns.
survival_data_full <- read.delim(
  "E:/2.ER/7.NCA.ER.Final/3.Gene Expression/3.Survival Analysis/Disc.ERpos/Disc.ERpos.Genes.Survival.txt",
  check.names = FALSE
)

##############################################################################
# 4) Prepare Excel Workbooks for Saving Results
##############################################################################
# We'll have separate Excel workbooks for each endpoint, both in a "regular" version
# and a "5-year" version (censored at 60 months).
wb_OS      <- createWorkbook()
wb_DSS     <- createWorkbook()
wb_RFS     <- createWorkbook()
wb_MFS     <- createWorkbook()
wb_OS_5yr  <- createWorkbook()
wb_DSS_5yr <- createWorkbook()
wb_RFS_5yr <- createWorkbook()
wb_MFS_5yr <- createWorkbook()

# We'll track which sheet we are on for each workbook with counters.
sheet_counter_regular <- 1
sheet_counter_5year   <- 1

##############################################################################
# 5) Helper Function: Adjust Data for 5-Year Survival Analysis
##############################################################################
# If a patient's time > 60 months, set the event to 0 (censored).
adjust_to_5_year <- function(data, columns) {
  data_5yr <- data
  # columns$time is the time column (e.g. OS(Months)), columns$event is the event column
  data_5yr[[columns$event]] <- ifelse(
    data_5yr[[columns$time]] > 60,
    0,
    data_5yr[[columns$event]]
  )
  data_5yr
}


endpoints_full_name <- list(
  OS  = "Overall survival",
  DSS = "Disease-specific survival",
  RFS = "Relapse-free survival",
  MFS = "Metastasis-free survival"
)


##############################################################################
# 6) Block 1: Create Survival Plots for All Data (No p-value Filtering)
##############################################################################
# We'll save these plots in a PDF file.
pdf(paste0("disc.ERpos_", Sys.Date(), "_TF_EXP_survival.pdf"))

# We assume that columns 11 to 16 in the survival_data_full are "TF" columns.
for (col_index in 11:16) {
  
  # Get the TF name (column header) and its values.
  TF_name   <- colnames(survival_data_full)[col_index]
  TF_values <- survival_data_full[[col_index]]
  
  # Only proceed if the column is not all NA.
  if (!all(is.na(TF_values))) {
    
    # Compute the average (mean) for this TF, ignoring NA values.
    disc_values   <- na.omit(TF_values)
    average_value <- mean(disc_values, na.rm = TRUE)
    cat("Average for", TF_name, ":", average_value, "\n")
    
    # Create a new grouping column: "High" if above average, "Low" if below or equal.
    survival_data_full$PDS_group <- ifelse(
      TF_values > average_value,
      paste("High", TF_name),
      paste("Low",  TF_name)
    )
    
    # Define a list of endpoints and the corresponding time/event columns.
    endpoints <- list(
      OS  = list(time = "OS(Months)",   event = "OS status"),
      DSS = list(time = "OS(Months)",   event = "DSS"),
      RFS = list(time = "RFS (Months)", event = "RFS"),
      MFS = list(time = "OS(Months)",   event = "MFS")
    )
    
    # Loop over each endpoint to create Kaplan-Meier plots.
    for (endpoint in names(endpoints)) {
      
      # Retrieve the time/event column names for the current endpoint.
      cols <- endpoints[[endpoint]]
      
      # Subset data to rows with valid (non-NA) time and event.
      current_data <- survival_data_full[
        complete.cases(
          survival_data_full[[cols$time]],
          survival_data_full[[cols$event]]
        ),
      ]
      
      # Create a unique ID for each patient (optional).
      current_data$SampleID <- current_data$`Patient ID`
      
      # Convert time to numeric and round to 4 decimals for clarity.
      current_data[[cols$time]] <- round(
        as.numeric(current_data[[cols$time]]), 4
      )
      
      # Create a Surv object for the regular (full) analysis.
      surv_obj <- Surv(
        time  = current_data[[cols$time]],
        event = as.numeric(current_data[[cols$event]])
      )
      # Fit the KM survival curve by the grouping variable PDS_group (High vs Low).
      surv_fit <- survfit(
        surv_obj ~ PDS_group,
        data = current_data,
        na.action = na.exclude
      )
      
      # Prepare the data for a 5-year analysis (censor if > 60 months).
      current_data_5yr <- adjust_to_5_year(current_data, cols)
      current_data_5yr[[cols$time]] <- round(
        as.numeric(current_data_5yr[[cols$time]]), 4
      )
      surv_obj_5yr <- Surv(
        time  = current_data_5yr[[cols$time]],
        event = as.numeric(current_data_5yr[[cols$event]])
      )
      surv_fit_5yr <- survfit(
        surv_obj_5yr ~ PDS_group,
        data = current_data_5yr,
        na.action = na.exclude
      )
      
      ############################################################################
      # Plot for Regular (Non-5-Year) Analysis
      ############################################################################
      # We'll use ggsurvplot from survminer. We add a bounding box around the plot
      # by specifying panel.border in the theme. This ensures the entire plot is enclosed.
      plot_regular <- ggsurvplot(
        fit        = surv_fit,
        data       = current_data,
        pval       = TRUE,
        pval.coord = c(0, 0.15),
        pval.size  = 4.8,
        pval.method       = TRUE,
        pval.method.size  = 4.8,
        pval.method.coord = c(0, 0.22),
        xlab       = "Time (months)",
        ylab = endpoints_full_name[[endpoint]],  # Use the full name here
        legend.labs= c("High", "Low"), #check correct in this experiment
        risk.table = TRUE,
        conf.int   = FALSE,
        risk.table.y.text = FALSE,
        ggtheme    = theme_minimal() +
          theme(
            panel.grid   = element_blank(),                  # Remove grid lines
            axis.line    = element_line(linewidth = 0.5, color = "black"),
            axis.ticks   = element_line(linewidth = 0.5, color = "black"),
            text         = element_text(size = 16), #endpoint of the y-axis size
            legend.text  = element_text(size = 14), #size of "high" and "low"
            
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
      plot_regular$plot <- plot_regular$plot + labs(color = paste0(TF_name,"-expression")) 
      
      # Post-process the risk table object
      plot_regular$table <- plot_regular$table +
        labs(y = NULL) +                             # Removes the y-axis label
        theme(axis.title.y = element_blank(),         # Reduce legend title size) # Also ensure the title is blank
              plot.title = element_text(size = 16))   #Number at risk size   
      
      # Print the plot to the currently open PDF device
      print(plot_regular)
      
      # Save the data used in the regular analysis to the appropriate Excel sheet.
      data_export_regular <- current_data[, c("SampleID", cols$time, cols$event, "PDS_group")]
      if (endpoint == "OS") {
        addWorksheet(wb_OS, paste0("Sheet", sheet_counter_regular))
        writeData(wb_OS, paste0("Sheet", sheet_counter_regular), data_export_regular)
      } else if (endpoint == "DSS") {
        addWorksheet(wb_DSS, paste0("Sheet", sheet_counter_regular))
        writeData(wb_DSS, paste0("Sheet", sheet_counter_regular), data_export_regular)
      } else if (endpoint == "RFS") {
        addWorksheet(wb_RFS, paste0("Sheet", sheet_counter_regular))
        writeData(wb_RFS, paste0("Sheet", sheet_counter_regular), data_export_regular)
      } else if (endpoint == "MFS") {
        addWorksheet(wb_MFS, paste0("Sheet", sheet_counter_regular))
        writeData(wb_MFS, paste0("Sheet", sheet_counter_regular), data_export_regular)
      }
      
      ############################################################################
      # Plot for 5-Year Analysis
      ############################################################################
      # Similar approach, but using the 5-year–censored data.
      plot_5year <- ggsurvplot(
        fit        = surv_fit_5yr,
        data       = current_data_5yr,
        pval       = TRUE,
        pval.coord = c(0, 0.15),
        pval.size  = 4.8,
        pval.method       = TRUE,
        pval.method.size  = 4.8,
        pval.method.coord = c(0, 0.22),
        xlab       = "Time (months)",
        ylab = paste0(endpoints_full_name[[endpoint]]," (5-year)"),  # Use the full name here
        legend.labs= c("High", "Low"),
        risk.table = TRUE,
        conf.int   = FALSE,
        risk.table.y.text = FALSE,
        ggtheme    = theme_minimal() +
          theme(
            panel.grid   = element_blank(),
            axis.line    = element_line(linewidth = 0.5, color = "black"),
            axis.ticks   = element_line(linewidth = 0.5, color = "black"),
            text         = element_text(size = 16), #endpoint of the y-axis size
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
      plot_5year$plot <- plot_5year$plot + labs(color = paste0(TF_name,"-expression")) 
      
      # Post-process the risk table object
      plot_5year$table <- plot_5year$table +
        labs(y = NULL) +                             # Removes the y-axis label
        theme(axis.title.y = element_blank(),         # Reduce legend title size) # Also ensure the title is blank
              plot.title = element_text(size = 16))         
      
      print(plot_5year)
      
      # Save the 5-year data in the appropriate Excel workbook.
      data_export_5year <- current_data_5yr[, c("SampleID", cols$time, cols$event, "PDS_group")]
      if (endpoint == "OS") {
        addWorksheet(wb_OS_5yr, paste0("Sheet", sheet_counter_5year))
        writeData(wb_OS_5yr, paste0("Sheet", sheet_counter_5year), data_export_5year)
      } else if (endpoint == "DSS") {
        addWorksheet(wb_DSS_5yr, paste0("Sheet", sheet_counter_5year))
        writeData(wb_DSS_5yr, paste0("Sheet", sheet_counter_5year), data_export_5year)
      } else if (endpoint == "RFS") {
        addWorksheet(wb_RFS_5yr, paste0("Sheet", sheet_counter_5year))
        writeData(wb_RFS_5yr, paste0("Sheet", sheet_counter_5year), data_export_5year)
      } else if (endpoint == "MFS") {
        addWorksheet(wb_MFS_5yr, paste0("Sheet", sheet_counter_5year))
        writeData(wb_MFS_5yr, paste0("Sheet", sheet_counter_5year), data_export_5year)
      }
      
      # Increment counters for the next iteration
      sheet_counter_regular <- sheet_counter_regular + 1
      sheet_counter_5year   <- sheet_counter_5year + 1
    }
  }
}

# Close the PDF device after finishing Block 1
dev.off()

##############################################################################
# 7) Block 2: Create Survival Plots for Significant Results (p-value < 0.05)
##############################################################################
pdf(paste0("disc.ERpos_", Sys.Date(), "_sig_TF_EXP_survival.pdf"))

# Create Excel workbooks for significant results
wb_OS_sig      <- createWorkbook()
wb_DSS_sig     <- createWorkbook()
wb_RFS_sig     <- createWorkbook()
wb_MFS_sig     <- createWorkbook()
wb_OS_5yr_sig  <- createWorkbook()
wb_DSS_5yr_sig <- createWorkbook()
wb_RFS_5yr_sig <- createWorkbook()
wb_MFS_5yr_sig <- createWorkbook()

sheet_counter_sig_regular <- 1
sheet_counter_sig_5year   <- 1

# Repeat a similar loop for columns 11 to 16, but only create plots if p < 0.05
for (col_index in 11:16) {
  
  TF_name   <- colnames(survival_data_full)[col_index]
  TF_values <- survival_data_full[[col_index]]
  
  if (!all(is.na(TF_values))) {
    disc_values    <- na.omit(TF_values)
    average_value  <- mean(disc_values, na.rm = TRUE)
    cat("Average for", TF_name, ":", average_value, "\n")
    
    survival_data_full$PDS_group <- ifelse(
      TF_values > average_value,
      paste("High", TF_name),
      paste("Low",  TF_name)
    )
    
    endpoints <- list(
      OS  = list(time = "OS(Months)",   event = "OS status"),
      DSS = list(time = "OS(Months)",   event = "DSS"),
      RFS = list(time = "RFS (Months)", event = "RFS"),
      MFS = list(time = "OS(Months)",   event = "MFS")
    )
    
    for (endpoint in names(endpoints)) {
      cols <- endpoints[[endpoint]]
      
      current_data <- survival_data_full[
        complete.cases(
          survival_data_full[[cols$time]],
          survival_data_full[[cols$event]]
        ),
      ]
      current_data$SampleID <- current_data$`Patient ID`
      current_data[[cols$time]] <- round(
        as.numeric(current_data[[cols$time]]), 4
      )
      
      surv_obj <- Surv(
        time  = current_data[[cols$time]],
        event = as.numeric(current_data[[cols$event]])
      )
      surv_fit <- survfit(
        surv_obj ~ PDS_group,
        data = current_data,
        na.action = na.exclude
      )
      
      # Calculate log-rank test p-value for the regular analysis
      logrank_test   <- survdiff(surv_obj ~ PDS_group, data = current_data)
      p_val_regular  <- 1 - pchisq(
        logrank_test$chisq,
        df = length(unique(current_data$PDS_group)) - 1
      )
      
      # Prepare 5-year data
      current_data_5yr <- adjust_to_5_year(current_data, cols)
      current_data_5yr[[cols$time]] <- round(
        as.numeric(current_data_5yr[[cols$time]]), 4
      )
      surv_obj_5yr <- Surv(
        time  = current_data_5yr[[cols$time]],
        event = as.numeric(current_data_5yr[[cols$event]])
      )
      surv_fit_5yr <- survfit(
        surv_obj_5yr ~ PDS_group,
        data = current_data_5yr,
        na.action = na.exclude
      )
      
      logrank_test_5yr <- survdiff(surv_obj_5yr ~ PDS_group, data = current_data_5yr)
      p_val_5yr        <- 1 - pchisq(
        logrank_test_5yr$chisq,
        df = length(unique(current_data_5yr$PDS_group)) - 1
      )
      
      ########################################################################
      # If p_val_regular < 0.05, create the regular significant plot
      ########################################################################
      if (p_val_regular < 0.05) {
        plot_sig_regular <- ggsurvplot(
          fit        = surv_fit,
          data       = current_data,
          pval       = TRUE,
          pval.coord = c(0, 0.15),
          pval.size  = 4.8,
          pval.method       = TRUE,
          pval.method.size  = 4.8,
          pval.method.coord = c(0, 0.22),
          xlab       = "Time (months)",
          ylab = endpoints_full_name[[endpoint]],  # Use the full name here
          legend.labs= c("High", "Low"),
          risk.table = TRUE,
          conf.int   = FALSE,
          risk.table.y.text = FALSE,
          ggtheme    = theme_minimal() +
            theme(
              panel.grid   = element_blank(),
              axis.line    = element_line(linewidth = 0.5, color = "black"),
              axis.ticks   = element_line(linewidth = 0.5, color = "black"),
              text         = element_text(size = 16), #endpoint of the y-axis size
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
        plot_sig_regular$plot <- plot_sig_regular$plot + labs(color = paste0(TF_name,"-expression")) 
        
        # Post-process the risk table object
        plot_sig_regular$table <- plot_sig_regular$table +
          labs(y = NULL) +                             # Removes the y-axis label
          theme(axis.title.y = element_blank(),         # Reduce legend title size) # Also ensure the title is blank
                plot.title = element_text(size = 16))         
        
        print(plot_sig_regular)
        
        # Save the data used for this significant regular analysis
        data_export_sig_regular <- current_data[, c("SampleID", cols$time, cols$event, "PDS_group")]
        if (endpoint == "OS") {
          addWorksheet(wb_OS_sig, paste0("Sheet", sheet_counter_sig_regular))
          writeData(wb_OS_sig, paste0("Sheet", sheet_counter_sig_regular), data_export_sig_regular)
        } else if (endpoint == "DSS") {
          addWorksheet(wb_DSS_sig, paste0("Sheet", sheet_counter_sig_regular))
          writeData(wb_DSS_sig, paste0("Sheet", sheet_counter_sig_regular), data_export_sig_regular)
        } else if (endpoint == "RFS") {
          addWorksheet(wb_RFS_sig, paste0("Sheet", sheet_counter_sig_regular))
          writeData(wb_RFS_sig, paste0("Sheet", sheet_counter_sig_regular), data_export_sig_regular)
        } else if (endpoint == "MFS") {
          addWorksheet(wb_MFS_sig, paste0("Sheet", sheet_counter_sig_regular))
          writeData(wb_MFS_sig, paste0("Sheet", sheet_counter_sig_regular), data_export_sig_regular)
        }
        sheet_counter_sig_regular <- sheet_counter_sig_regular + 1
      }
      
      ########################################################################
      # If p_val_5yr < 0.05, create the 5-year significant plot
      ########################################################################
      if (p_val_5yr < 0.05) {
        plot_sig_5year <- ggsurvplot(
          fit        = surv_fit_5yr,
          data       = current_data_5yr,
          pval       = TRUE,
          pval.coord = c(0, 0.15),
          pval.size  = 4.8,
          pval.method       = TRUE,
          pval.method.size  = 4.8,
          pval.method.coord = c(0, 0.22),
          xlab       = "Time (months)",
          ylab = paste0(endpoints_full_name[[endpoint]]," (5-year)"),  # Use the full name here
          legend.labs= c("High", "Low"),
          risk.table = TRUE,
          conf.int   = FALSE,
          risk.table.y.text = FALSE,
          ggtheme    = theme_minimal() +
            theme(
              panel.grid   = element_blank(),
              axis.line    = element_line(linewidth = 0.5, color = "black"),
              axis.ticks   = element_line(linewidth = 0.5, color = "black"),
              text         = element_text(size = 16), #endpoint of the y-axis size
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
        plot_sig_5year$plot <- plot_sig_5year$plot + labs(color = paste0(TF_name,"-expression")) 
        
        # Post-process the risk table object
        plot_sig_5year$table <- plot_sig_5year$table +
          labs(y = NULL) +                             # Removes the y-axis label
          theme(axis.title.y = element_blank(),         # Reduce legend title size) # Also ensure the title is blank
                plot.title = element_text(size = 16))         
        
        print(plot_sig_5year)
        
        # Save the data used for this significant 5-year analysis
        data_export_sig_5year <- current_data_5yr[, c("SampleID", cols$time, cols$event, "PDS_group")]
        if (endpoint == "OS") {
          addWorksheet(wb_OS_5yr_sig, paste0("Sheet", sheet_counter_sig_5year))
          writeData(wb_OS_5yr_sig, paste0("Sheet", sheet_counter_sig_5year), data_export_sig_5year)
        } else if (endpoint == "DSS") {
          addWorksheet(wb_DSS_5yr_sig, paste0("Sheet", sheet_counter_sig_5year))
          writeData(wb_DSS_5yr_sig, paste0("Sheet", sheet_counter_sig_5year), data_export_sig_5year)
        } else if (endpoint == "RFS") {
          addWorksheet(wb_RFS_5yr_sig, paste0("Sheet", sheet_counter_sig_5year))
          writeData(wb_RFS_5yr_sig, paste0("Sheet", sheet_counter_sig_5year), data_export_sig_5year)
        } else if (endpoint == "MFS") {
          addWorksheet(wb_MFS_5yr_sig, paste0("Sheet", sheet_counter_sig_5year))
          writeData(wb_MFS_5yr_sig, paste0("Sheet", sheet_counter_sig_5year), data_export_sig_5year)
        }
        sheet_counter_sig_5year <- sheet_counter_sig_5year + 1
      }
    }
  }
}

# Close the PDF device for the significant results
dev.off()

##############################################################################
# 8) Save the Excel Workbooks for Regular and 5-Year Analyses
##############################################################################
saveWorkbook(wb_OS,     paste0("disc.ERpos_OS_",        Sys.Date(), ".xlsx"), overwrite = TRUE)
saveWorkbook(wb_DSS,    paste0("disc.ERpos_DSS_",       Sys.Date(), ".xlsx"), overwrite = TRUE)
saveWorkbook(wb_RFS,    paste0("disc.ERpos_RFS_",       Sys.Date(), ".xlsx"), overwrite = TRUE)
saveWorkbook(wb_MFS,    paste0("disc.ERpos_MFS_",       Sys.Date(), ".xlsx"), overwrite = TRUE)
saveWorkbook(wb_OS_5yr, paste0("disc.ERpos_OS_5YEAR_",  Sys.Date(), ".xlsx"), overwrite = TRUE)
saveWorkbook(wb_DSS_5yr,paste0("disc.ERpos_DSS_5YEAR_", Sys.Date(), ".xlsx"), overwrite = TRUE)
saveWorkbook(wb_RFS_5yr,paste0("disc.ERpos_RFS_5YEAR_", Sys.Date(), ".xlsx"), overwrite = TRUE)
saveWorkbook(wb_MFS_5yr,paste0("disc.ERpos_MFS_5YEAR_", Sys.Date(), ".xlsx"), overwrite = TRUE)

##############################################################################
# Save the Excel Workbooks for Significant Results
##############################################################################
saveWorkbook(wb_OS_sig,     paste0("disc.ERpos_OS_sig_",        Sys.Date(), ".xlsx"), overwrite = TRUE)
saveWorkbook(wb_DSS_sig,    paste0("disc.ERpos_DSS_sig_",       Sys.Date(), ".xlsx"), overwrite = TRUE)
saveWorkbook(wb_RFS_sig,    paste0("disc.ERpos_RFS_sig_",       Sys.Date(), ".xlsx"), overwrite = TRUE)
saveWorkbook(wb_MFS_sig,    paste0("disc.ERpos_MFS_sig_",       Sys.Date(), ".xlsx"), overwrite = TRUE)
saveWorkbook(wb_OS_5yr_sig, paste0("disc.ERpos_OS_5YEAR_sig_",  Sys.Date(), ".xlsx"), overwrite = TRUE)
saveWorkbook(wb_DSS_5yr_sig,paste0("disc.ERpos_DSS_5YEAR_sig_", Sys.Date(), ".xlsx"), overwrite = TRUE)
saveWorkbook(wb_RFS_5yr_sig,paste0("disc.ERpos_RFS_5YEAR_sig_", Sys.Date(), ".xlsx"), overwrite = TRUE)
saveWorkbook(wb_MFS_5yr_sig,paste0("disc.ERpos_MFS_5YEAR_sig_", Sys.Date(), ".xlsx"), overwrite = TRUE)

##############################################################################
# 9) Save the R Environment
##############################################################################
session_info <- sessionInfo()  # Capture session info for reproducibility
save(list = ls(), file = paste0("Disc.ERpos_", Sys.Date(), ".RData"))
save.image(paste0("Image.disc.ERpos_", Sys.Date(), ".RData"))
