#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(jsonlite)
    library(purrr)
    library(tibble)
})

base_dir <- file.path("output", "simulated_instances")
if (!dir.exists(base_dir)) {
    stop(sprintf("Directory not found: %s", base_dir))
}

json_files <- list.files(base_dir, pattern = "\\.json$", recursive = TRUE, full.names = TRUE)
json_files <- json_files[!grepl("/detail_", json_files)]
json_files <- json_files[!grepl("(_inv|_sym)/", json_files)]
if (!length(json_files)) {
    stop("No JSON files found under ", normalizePath(base_dir))
}

as_numeric_matrix <- function(x) {
    if (is.null(x)) return(NULL)
    if (is.array(x) || is.matrix(x)) {
        return(apply(x, c(1, 2), as.numeric))
    }
    if (is.list(x) && !is.data.frame(x)) {
        rows <- lapply(x, function(r) as.numeric(unlist(r)))
        mat <- do.call(rbind, rows)
        return(mat)
    }
    as.matrix(as.data.frame(x))
}

compute_mae_p2 <- function(js) {
    p_true <- as_numeric_matrix(js$p_true)
    p_est <- as_numeric_matrix(js$p_est)
    if (is.null(p_true) || is.null(p_est)) return(NA_real_)
    if (!all(dim(p_true) == dim(p_est))) return(NA_real_)
    G <- nrow(p_true)
    if (!length(G) || !is.finite(G) || G == 0) return(NA_real_)
    total_abs <- sum(abs(p_est - p_true))
    total_abs / (2 * G)
}

update_file <- function(path) {
    js <- tryCatch(fromJSON(path, simplifyVector = FALSE), error = function(e) NULL)
    if (is.null(js)) return(FALSE)
    mae_p2 <- compute_mae_p2(js)
    if (!is.finite(mae_p2)) return(FALSE)
    js$MAE_p2 <- mae_p2
    tmp <- paste0(path, ".tmp")
    write_json(js, tmp, pretty = TRUE, auto_unbox = TRUE)
    file.rename(tmp, path)
    TRUE
}

results <- map_lgl(json_files, update_file)
message(sprintf("Updated %d/%d simulated-instance JSON files with MAE_p2.", sum(results), length(results)))
