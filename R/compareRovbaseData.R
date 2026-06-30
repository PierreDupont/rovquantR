#' @title Data set comparison.
#'
#' @description
#' \code{compareRovbaseData} compares two versions of the same dataframe and prints out an \code{.html} report summarizing the differences.
#' It can be used to look for changes in columns names, identify new data, flag data with differences between versions.
#' 
#' @name compareRovbaseData
#' 
#' @param data.dir the \code{path} pointing to the directory containing the raw data from Rovbase.
#' @param working.dir the \code{path} pointing to the working directory. By default, the cleaned data will be stored in a subfolder of this working directory called 'data'.
#' @param species A \code{character} string with the name of the focal species
#'  ("bear", "wolf", or "wolverine").
#' @param years A \code{numeric} vector containing the years of interest. 
#' Only data for those years will be cleaned and returned.
#' @param two.sex A \code{logical} determining whether the analysis will be done by sex (two.sex = T) or both together (two.sex = F).
#' @param sampling.months (Optional) A \code{list} containing the sampling period months. If the sampling period overlaps two calendar years, the list should contain one element per year (e.g. samplingMonths <- list(c(11,12), c(1,2,3,4))) for a sampling period extending from November to April of the following year.
#' @param rename.list (Optional) A named \code{character} vector used to rename columns in the raw Rovbase files.
#' @param legal.dead (Optional) A named \code{character} vector used to identify legal death causes.
#' @param print.report A \code{logical} denoting whether to print out a \code{.html} report summarizing the cleaning process or not.
#' @param Rmd.template (Optional) The \code{path} to the \code{.rmd} template to be used for cleaning the data. By default, the \code{.rmd} template provided with the \code{rovquantR} package is used.  
#' @param overwrite A \code{logical} (default = FALSE) to force overwriting of previously existing clean data.
#'  If FALSE, the function checks for any pre-existing clean data files and ask whether to overwrite it or not.
#' @param output.dir (Optional) the \code{path} pointing to the directory where the \code{.html} report will be printed.
#' By default, the \code{.html} report describing the content of the clean data will 
#' be placed in a subfolder of the working directory (\code{working.dir}) called 'reports'.
#' 
#' @return This function returns:
#' \enumerate{
#' \item A \code{.RData} file with the clean NGS and dead recovery data objects 
#' for the species and period specified. The clean data file is saved as an \code{.RData} 
#' file named using the species name and the date of extraction of the raw Rovbase data 
#' to facilitate replicability (e.g. 'CleanData_bear_2024-08-10.RData').
#' \item A \code{.html} report summarizing the data cleaning process. 
#' The \code{.RData} report is using the same naming convention as the clean \code{.RData} (e.g. 'CleanData_bear_2024-08-10.html').
#' \item Additional \code{.png} images and summary \code{.csv} tables that can be reused somewhere else.
#' }
#'
#' @author Pierre Dupont
#' 
#' @importFrom readxl render
#' @importFrom dplyr png
#' @importFrom purrr mtext 
#' @import writexl 
#' @import tibble
#' 
#' @rdname compareRovbaseData
#' @export
compareRovbaseData <- function(
    df1,
    df2,
    keys,
    compare_cols = NULL,
    
    ##-- miscellanious
    print.report = FALSE,
    Rmd.template = NULL,
    output.dir = NULL,
    verbose = TRUE,
    overwrite = FALSE) {
  
  ## ----- 0. Helper functions ----- 
  
  values_equal <- function(a, b) {
    if (length(a) == 0 || length(b) == 0) {
      return(FALSE)
    }
    if (is.na(a) && is.na(b)) {
      return(TRUE)
    }
    if (is.na(a) || is.na(b)) {
      return(FALSE)
    }
    identical(a, b)
  }
  
  make_key_id <- function(df, keys) {
    if (length(keys) == 1) {
      return(as.character(df[[keys]]))
    }
    apply(df[, keys, drop = FALSE], 1, function(x) paste(x, collapse = " | "))
  }
  
  
  ## ----- 1. Load files ----- 

  if (verbose) {
    cat("Dataset 1 rows:", nrow(df1), "| columns:", ncol(df1), "\n")
    cat("Dataset 2 rows:", nrow(df2), "| columns:", ncol(df2), "\n\n")
  }
  
  
  ## ----- 2. Compare column names ----- 
  
  cols1 <- names(df1)
  cols2 <- names(df2)
  
  column_comparison <- tibble(
    column = union(cols1, cols2),
    in_dataset1 = column %in% cols1,
    in_dataset2 = column %in% cols2
  )
  
  exact_match <- identical(cols1, cols2)
  same_names_ignore_order <- setequal(cols1, cols2)
  
  if (verbose) {
    cat("Column name comparison\n")
    cat("----------------------\n")
    cat("Exact match:", exact_match, "\n")
    cat("Same names ignoring order:", same_names_ignore_order, "\n\n")
  }
  
  
  ## ----- 3. Check keys exist -----
  
  missing_keys_df1 <- setdiff(keys, names(df1))
  missing_keys_df2 <- setdiff(keys, names(df2))
  
  if (length(missing_keys_df1) > 0) {
    stop(paste("Missing key(s) in dataset 1:", paste(missing_keys_df1, collapse = ", ")))
  }
  
  if (length(missing_keys_df2) > 0) {
    stop(paste("Missing key(s) in dataset 2:", paste(missing_keys_df2, collapse = ", ")))
  }
  
  if (verbose) {
    cat("All key columns exist in both datasets:\n")
    cat(paste(keys, collapse = ", "), "\n\n")
  }
  
  
  ## ----- 4. Check key uniqueness -----
  
  dup_df1 <- df1 %>%
    count(across(all_of(keys)), name = "n") %>%
    filter(n > 1)
  
  dup_df2 <- df2 %>%
    count(across(all_of(keys)), name = "n") %>%
    filter(n > 1)
  
  if (verbose) {
    cat("Duplicate key check\n")
    cat("-------------------\n")
    cat("Duplicated key combinations in dataset 1:", nrow(dup_df1), "\n")
    cat("Duplicated key combinations in dataset 2:", nrow(dup_df2), "\n\n")
  }
  
  if (nrow(dup_df1) > 0) {
    stop("Dataset 1 contains duplicated key values/combinations. Comparison stopped.")
  }
  
  if (nrow(dup_df2) > 0) {
    stop("Dataset 2 contains duplicated key values/combinations. Comparison stopped.")
  }
  
  
  ## ----- 5. Determine columns to compare ----- 
  
  common_cols <- intersect(names(df1), names(df2))
  non_key_common_cols <- setdiff(common_cols, keys)
  
  if (is.null(compare_cols)) {
    cols_to_check <- non_key_common_cols
  } else {
    missing_compare_cols_df1 <- setdiff(compare_cols, names(df1))
    missing_compare_cols_df2 <- setdiff(compare_cols, names(df2))
    
    if (length(missing_compare_cols_df1) > 0) {
      stop(paste("Selected compare_cols missing from dataset 1:",
                 paste(missing_compare_cols_df1, collapse = ", ")))
    }
    
    if (length(missing_compare_cols_df2) > 0) {
      stop(paste("Selected compare_cols missing from dataset 2:",
                 paste(missing_compare_cols_df2, collapse = ", ")))
    }
    
    cols_to_check <- setdiff(compare_cols, keys)
  }
  
  if (length(cols_to_check) == 0) {
    stop("No non-key columns available for comparison.")
  }
  
  if (verbose) {
    cat("Columns selected for value comparison:\n")
    cat(paste(cols_to_check, collapse = ", "), "\n\n")
  }
  
  
  
  ## ----- 6. Find rows only in one dataset -----
  
  keys_only_in_df1 <- anti_join(df1, df2, by = keys) %>%
    select(all_of(keys))
  
  keys_only_in_df2 <- anti_join(df2, df1, by = keys) %>%
    select(all_of(keys))
  
  rows_only_in_df1 <- anti_join(df1, df2, by = keys)
  rows_only_in_df2 <- anti_join(df2, df1, by = keys)
  
  
  
  ## ----- 7. Compare shared rows -----
  
  common <- inner_join(df1, df2, by = keys, suffix = c(".old", ".new"))
  
  diff_table <- map_dfr(
    cols_to_check,
    function(col) {
      old_col <- paste0(col, ".old")
      new_col <- paste0(col, ".new")
      
      keep <- mapply(
        function(a, b) !values_equal(a, b),
        common[[old_col]],
        common[[new_col]])
      
      common[keep, , drop = FALSE] %>%
        transmute(across(all_of(keys)),
                  column = col,
                  old_value = as.character(.data[[old_col]]),
                  new_value = as.character(.data[[new_col]])) 
    }) %>%
    arrange(across(all_of(keys)), column)
  
  
  
  ## ----- 8. Row-level summary -----
  
  if (nrow(diff_table) > 0) {
    changed_rows <- diff_table %>%
      group_by(across(all_of(keys))) %>%
      summarise(
        diff_cols = paste(column, collapse = ", "),
        n_diff_cols = n(),
        .groups = "drop") %>%
      arrange(across(all_of(keys)))
  } else {
    changed_rows <- tibble()
  }
  
  
  
  ## ----- 9: Summary table -----
  
  summary_table <- tibble(
    metric = c( "dataset1_rows",
                "dataset2_rows",
                "dataset1_columns",
                "dataset2_columns",
                "exact_column_name_match",
                "same_column_names_ignore_order",
                "duplicate_key_combinations_dataset1",
                "duplicate_key_combinations_dataset2",
                "rows_only_in_dataset1",
                "rows_only_in_dataset2",
                "shared_rows",
                "cell_level_differences",
                "changed_rows"),
    value = c( nrow(df1),
               nrow(df2),
               ncol(df1),
               ncol(df2),
               exact_match,
               same_names_ignore_order,
               nrow(dup_df1),
               nrow(dup_df2),
               nrow(rows_only_in_df1),
               nrow(rows_only_in_df2),
               nrow(common),
               nrow(diff_table),
               ifelse(nrow(diff_table) > 0,
                      nrow(distinct(diff_table, across(all_of(keys)))),
                      0)))
  
  if (verbose) {
    cat("Summary\n")
    cat("-------\n")
    print(summary_table)
    cat("\n")
  }
  
  
  ## ----- 10: Optional export to xlsx -----
  
  # if (export_xlsx) {
  #   export_list <- list(
  #     summary = summary_table,
  #     column_comparison = column_comparison,
  #     duplicate_keys_dataset1 = dup_df1,
  #     duplicate_keys_dataset2 = dup_df2,
  #     keys_only_in_dataset1 = keys_only_in_df1,
  #     keys_only_in_dataset2 = keys_only_in_df2,
  #     rows_only_in_dataset1 = rows_only_in_df1,
  #     rows_only_in_dataset2 = rows_only_in_df2,
  #     changed_rows_summary = changed_rows,
  #     cell_level_differences = diff_table
  #   )
  #   
  #   writexl::write_xlsx(export_list, path = output_file)
  #   
  #   if (verbose) {
  #     cat("Results exported to:\n")
  #     cat(output_file, "\n\n")
  #   }
  # }
  
  
  
  info.ls <- list(
    data1 = df1,
    data2 = df2,
    summary = summary_table,
    column_comparison = column_comparison,
    duplicate_keys_dataset1 = dup_df1,
    duplicate_keys_dataset2 = dup_df2,
    keys_only_in_dataset1 = keys_only_in_df1,
    keys_only_in_dataset2 = keys_only_in_df2,
    rows_only_in_dataset1 = rows_only_in_df1,
    rows_only_in_dataset2 = rows_only_in_df2,
    changed_rows_summary = changed_rows,
    cell_level_differences = diff_table,
    compared_columns = cols_to_check,
    keys = keys)
  
  
  
  
  
  ## ----- 10. PRINT REPORT -----

  
  if (print.report) {

    ##-- List of rendering parameters for the .Rmd report
    if (engSpecies == "bear") {
      info.ls$remove.alive <- remove.alive
      info.ls$remove.dead <- remove.dead
    }

    if (engSpecies == "wolverine") {
      info.ls$youngDeads <- youngDeads
      info.ls$lowWeightDeads <- lowWeightDeads
      info.ls$zeroWeightDeads <- zeroWeightDeads
    }

    if (engSpecies == "wolf") {}

    ##-- Find the .rmd template for the report
    if(is.null(Rmd.template)) {
      Rmd.template <- system.file("rmd", "RovQuant_CompareReport.Rmd", package = "rovquantR")
      if(!file.exists(Rmd.template)) {
        stop("Can not find the Rmarkdown document to use for cleaning Rovbase.3.0 data.\n You must provide the path to the Rmarkdown template through the \"Rmd_template\" argument.")
      }
    }

    ##-- Find the directory to print the report
    if(is.null(output.dir)){ output.dir <- file.path(working.dir, "reports") }

    ##-- Render .Rmd report
    rmarkdown::render(
      input = Rmd.template,
      params = list( species = SPECIES,
                     years = years,
                     sampling.months = sampling.months,
                     data.dir = data.dir,
                     working.dir = working.dir,
                     date = DATE,
                     info.ls = info.ls),
      output_dir = output.dir,
      output_file = paste0("CompareData_", engSpecies, "_", DATE,".html"))
  }
  
  
  

  # Return results
  return(info.ls)
}



##------------------------------------------------------------------------------