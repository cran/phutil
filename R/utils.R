check_2column_matrix <- function(x, warn = TRUE) {
  if (!is.matrix(x)) {
    cli::cli_abort(
      c(
        "{.arg x} must be a matrix.",
        "i" = "{.arg x} is of class {.cls {class(x)}}."
      )
    )
  }

  if (ncol(x) != 2L) {
    cli::cli_abort(
      c(
        "{.arg x} must be a matrix with 2 columns.",
        "i" = "{.arg x} has {ncol(x)} columns."
      )
    )
  }

  if (!is.numeric(x)) {
    cli::cli_abort("Input must be a matrix with numeric entries.")
  }

  if (anyNA(x)) {
    cli::cli_abort("Input must be a matrix with no missing values.")
  }

  if (any(x[, 1L] > x[, 2L])) {
    if (warn) {
      cli::cli_alert_warning(
        "Birth values are expected to be smaller than death values."
      )
    }
    return(FALSE)
  }

  TRUE
}

capitalize <- function(x) {
  gsub("(?<=\\b)([a-z])", "\\U\\1", tolower(x), perl = TRUE)
}
