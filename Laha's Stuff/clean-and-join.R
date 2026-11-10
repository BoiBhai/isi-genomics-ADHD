#!/usr/bin/env Rscript

library(tidyverse)
library(glue) # For string interpolation
library(fs) # For file/path manipulation

folder_path <- "./Data Cleaning"
file_names <- c("genotype.csv", "expressions.csv", "concentrations.csv", "cprs-r.csv", "pharma.csv")

# Create the full file paths
file_paths <- fs::path(folder_path, file_names)

# Get the "clean" names for prefixes (e.g., "concentrations")
# This removes the folder path and the ".csv" extension
sheet_names <- file_names %>%
  fs::path_file() %>%
  fs::path_ext_remove()

load_and_process <- function(file_path, sheet_name) {
  df <- read_csv(file_path, show_col_types = FALSE)
  df_processed <- df %>%
    mutate(
      ID = ID %>%
        stringr::str_remove_all("P") %>%
        stringr::str_remove("^0+")
    )

  # Pivot is necessary if there is a column ending in "Of"
  col_names <- names(df_processed)
  pivot_names_col <- col_names[stringr::str_ends(col_names, "Of")]
  # Check if we found exactly one such column
  if (length(pivot_names_col) == 1) {
    # Automatically derive the value column name
    pivot_values_col <- sub("Of$", "", pivot_names_col)

    # Check if this derived value column *actually exists*
    if (pivot_values_col %in% col_names) {
      # Perform the pivot
      df_final <- df_processed %>%
        tidyr::pivot_wider(
          id_cols = ID, # Keep the processed ID as the identifier
          names_from = all_of(pivot_names_col),
          values_from = all_of(pivot_values_col)
        )
    } else {
      # The derived value column does not exist; skip pivoting
      warning(glue(
        "In file '{sheet_name}', found pivot column '{pivot_names_col}' ",
        "but derived value column '{pivot_values_col}' does not exist. ",
        "SKIPPING PIVOT for this file."
      ))
      df_final <- df_processed
    }
  } else {
    # No column ending in "Of" was found. No pivot needed.
    df_final <- df_processed
  }

  # Prefix all columns (except ID) with the sheet name
  df_final %>%
    rename_with(
      .fn = ~ glue("{sheet_name}.{.x}"), # e.g., "concentrations.Age"
      .cols = -ID # Apply to all columns *except* ID
    )
}

# Use purrr::map2() to apply the function to each file path and sheet name
all_data_list <- purrr::map2(
  .x = file_paths,
  .y = sheet_names,
  .f = ~ load_and_process(file_path = .x, sheet_name = .y)
)

# Use purrr::reduce() to "fold" the list, joining one data frame at a time
# This performs a full (outer) join on all elements by the "ID" column
final_joined_df <- all_data_list %>%
  reduce(full_join, by = "ID") %>%
  mutate(
    Type = case_when(
      str_detect(ID, "^C") ~ "Control",
      str_detect(ID, "^KC") ~ "Control",
      str_detect(ID, "^UC") ~ "Control",
      str_detect(ID, "^RR") ~ "Control",
      str_detect(ID, "^AS") ~ "Control",
      str_detect(ID, "^H") ~ "Case",
      str_detect(ID, "^SB") ~ "Case",
      str_detect(ID, "^RG") ~ "Case",
      TRUE ~ NA
    )
  ) %>%
  mutate(
    ID = {
      m <- stringr::str_match(ID, "^([A-Z]+)([0-9]+)$")
      prefix <- m[, 2]
      num <- m[, 3]
      valid <- !is.na(prefix) & !is.na(num)

      num_width <- 5 - nchar(prefix)

      new_id <- ifelse(
        valid,
        ifelse(
          num_width >= 0,
          paste0(prefix, stringr::str_pad(num, width = num_width, side = "left", pad = "0")),
          stringr::str_sub(prefix, 1, 5)
        ),
        ID
      )

      new_id
    }
  )

final_joined_df %>%
  arrange(Type, ID) %>%
  write_csv("cleaned_joined_data.csv", na = "")
