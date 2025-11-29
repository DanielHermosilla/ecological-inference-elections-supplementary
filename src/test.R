#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(purrr)
    library(stringr)
    library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
outfile <- NULL
if (length(args) && !startsWith(args[1], "--")) {
    outfile <- args[1]
    args <- args[-1]
}

inst_like <- {
    m <- grep("^--inst-like=", args, value = TRUE)
    if (length(m)) sub("^--inst-like=", "", m[1]) else "ei_"
}
limit_inst <- {
    m <- grep("^--limit-inst=", args, value = TRUE)
    if (length(m)) suppressWarnings(as.integer(sub("^--limit-inst=", "", m[1]))) else NA_integer_
}

base_dir <- file.path("output", "ei_instances")
if (!dir.exists(base_dir)) {
    stop(sprintf("'%s' was not found. Run this script from the project root.", base_dir))
}

all_instances <- sort(list.dirs(base_dir, recursive = FALSE, full.names = FALSE))
instances <- all_instances[str_detect(all_instances, paste0("^", inst_like))]
if (length(instances) == 0) {
    stop(sprintf("No instance folders starting with '%s' under '%s'.", inst_like, base_dir))
}
if (is.finite(limit_inst) && limit_inst > 0) {
    instances <- head(instances, limit_inst)
}

method_pairs <- tibble::tibble(
    method_key = c(
        "mult_project_lp_FALSE",
        "mvn_cdf_project_lp_FALSE",
        "mvn_pdf_project_lp_FALSE"
    ),
    standard_method = method_key,
    inv_method = paste0(method_key, "_inv")
) %>%
    pivot_longer(
        cols = c(standard_method, inv_method),
        names_to = "method_type_src",
        values_to = "method"
    ) %>%
    mutate(
        method_type = if_else(method_type_src == "inv_method", "inv", "standard")
    ) %>%
    select(method_key, method_type, method)

target_methods <- unique(method_pairs$method)

first_non_na <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) NA_real_ else x[1]
}

extract_value <- function(name, pattern) {
    m <- str_match(name, pattern)
    ifelse(is.na(m[, 2]), NA_real_, as.numeric(m[, 2]))
}

read_ei_v <- function(district_path, method, field = "EI_V") {
    method_path <- file.path(district_path, method)
    if (!dir.exists(method_path)) {
        return(NA_real_)
    }
    files <- Sys.glob(file.path(method_path, "*.json"))
    files <- files[!grepl("^detail_", basename(files))]
    if (length(files) == 0) {
        return(NA_real_)
    }
    vals <- vapply(
        files,
        function(f) {
            tryCatch(
                {
                    js <- fromJSON(f, simplifyVector = TRUE)
                    v <- js[[field]]
                    if (is.null(v)) NA_real_ else as.numeric(v)
                },
                error = function(e) NA_real_
            )
        },
        numeric(1)
    )
    vals <- vals[is.finite(vals)]
    if (length(vals) == 0) {
        return(NA_real_)
    }
    mean(vals)
}

strip_suffix <- function(name) {
    name <- sub("(_G[0-9]+_C[0-9]+|_C[0-9]+_G[0-9]+)$", "", name)
    name
}

collect_folder_records <- function(inst_name, district_folder) {
    district_path <- file.path(base_dir, inst_name, district_folder)
    base_district <- strip_suffix(district_folder)
    vals <- vapply(
        target_methods,
        function(m) read_ei_v(district_path, m),
        numeric(1)
    )
    tibble(
        inst = inst_name,
        district_folder = district_folder,
        district = base_district,
        G = extract_value(district_folder, "_G(\\d+)"),
        C = extract_value(district_folder, "_C(\\d+)"),
        method = target_methods,
        EI_V = vals
    )
}

collect_instance_records <- function(inst_name) {
    district_folders <- sort(list.dirs(file.path(base_dir, inst_name), recursive = FALSE, full.names = FALSE))
    if (length(district_folders) == 0) {
        return(tibble(
            inst = character(),
            district_folder = character(),
            district = character(),
            G = numeric(),
            C = numeric(),
            method = character(),
            EI_V = numeric()
        ))
    }
    map_dfr(district_folders, ~ collect_folder_records(inst_name, .x))
}

records <- map_dfr(instances, collect_instance_records)
if (nrow(records) == 0) {
    stop("No method results were found for the requested instances.")
}

records <- records %>%
    inner_join(method_pairs, by = "method") %>%
    filter(!is.na(EI_V))

orientation <- records %>%
    filter(method_type == "standard") %>%
    group_by(inst, district) %>%
    summarise(
        G = first_non_na(G),
        C = first_non_na(C),
        .groups = "drop"
    )

folders_by_type <- records %>%
    group_by(inst, district, method_type) %>%
    summarise(folder = first(district_folder), .groups = "drop") %>%
    pivot_wider(
        names_from = method_type,
        values_from = folder,
        names_prefix = "folder_"
    )

comparison <- records %>%
    select(inst, district, method_key, method_type, EI_V) %>%
    group_by(inst, district, method_key, method_type) %>%
    summarise(EI_V = mean(EI_V), .groups = "drop") %>%
    pivot_wider(
        names_from = method_type,
        values_from = EI_V,
        names_prefix = "EI_V_"
    ) %>%
    left_join(orientation, by = c("inst", "district")) %>%
    left_join(folders_by_type, by = c("inst", "district")) %>%
    mutate(diff_EI_V = EI_V_inv - EI_V_standard) %>%
    arrange(inst, district, method_key)

if (!is.null(outfile)) {
    utils::write.csv(comparison, file = outfile, row.names = FALSE)
    message(sprintf("Saved comparison table to '%s'.", outfile))
}

print(comparison, n = min(nrow(comparison), 50))

write.csv(comparison, file = "comparacion.csv")
