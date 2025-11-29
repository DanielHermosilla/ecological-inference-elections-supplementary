#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(stringr)
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(future.apply)
    library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)

parse_opt_int <- function(pos, default = NA_integer_) {
    if (length(args) >= pos) suppressWarnings(as.integer(args[pos])) else default
}
parse_opt_chr <- function(pos, default = NULL) {
    if (length(args) >= pos) args[pos] else default
}

field_to_average <- {
    ft <- parse_opt_chr(1, NULL)
    if (is.null(ft) || is.na(ft)) {
        stop("Missing <field_to_average>. Example: Rscript src/table_simulated.R MAE_p2 output/table_simulated.tex")
    }
    ft
}

outfile <- {
    tmp <- parse_opt_chr(2, NA_character_)
    if (!is.na(tmp)) tmp else NULL
}
req_workers <- parse_opt_int(3, NA_integer_)
scale_100 <- {
    m <- grep("^--percent$", args, value = TRUE)
    length(m) > 0
}

inst_like <- {
    m <- grep("^--inst-like=", args, value = TRUE)
    if (length(m)) sub("^--inst-like=", "", m[1]) else "I"
}
limit_inst <- {
    m <- grep("^--limit-inst=", args, value = TRUE)
    if (length(m)) suppressWarnings(as.integer(sub("^--limit-inst=", "", m[1]))) else NA_integer_
}

base_sim <- file.path("output", "simulated_instances")
if (!dir.exists(base_sim)) {
    stop(sprintf("Directory not found: %s", base_sim))
}

combos_interest <- expand.grid(
    G = c(2L, 3L, 4L, 6L, 8L),
    C = c(2L, 3L, 5L, 10L),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
) %>%
    mutate(combo = paste0("G", G, "_C", C))
combo_levels <- combos_interest$combo

base_methods <- c(
    "ei.MD.bayes",
    "nslphom_dual_w",
    "mvn_cdf_project_lp_FALSE",
    "mvn_pdf_project_lp_FALSE",
    "mult_project_lp_FALSE",
    "exact_project_lp_FALSE"
)

instances_all <- sort(list.dirs(base_sim, recursive = FALSE, full.names = FALSE))
instances_all <- instances_all[str_detect(instances_all, paste0("^", inst_like))]

parse_inst <- function(name) {
    m <- regexec("^I(\\d+)_B(\\d+)_G(\\d+)_C(\\d+)_lambda(\\d+)$", name)
    g <- regmatches(name, m)[[1]]
    if (length(g) != 6) {
        return(NULL)
    }
    tibble(
        inst = name,
        I = as.integer(g[2]),
        B = as.integer(g[3]),
        G = as.integer(g[4]),
        C = as.integer(g[5]),
        lambda = as.integer(g[6]),
        combo = paste0("G", g[4], "_C", g[5])
    )
}

inst_info <- bind_rows(lapply(instances_all, parse_inst))
inst_info <- inst_info %>% filter(combo %in% combo_levels)
if (!nrow(inst_info)) {
    stop("No simulated instances match the requested G/C combinations.")
}

instances <- inst_info$inst
if (is.finite(limit_inst) && limit_inst > 0) {
    instances <- head(instances, limit_inst)
    inst_info <- inst_info[match(instances, inst_info$inst), , drop = FALSE]
}

detect_methods <- function(instances) {
    method_sets <- lapply(instances, function(inst) {
        path <- file.path(base_sim, inst)
        if (!dir.exists(path)) character(0) else sort(list.dirs(path, recursive = FALSE, full.names = FALSE))
    })
    sort(unique(unlist(method_sets)))
}
available_methods <- detect_methods(instances)
methods <- base_methods[base_methods %in% available_methods]
if (!length(methods)) {
    stop("None of the expected methods were found in simulated instances.")
}

if (is.finite(req_workers) && req_workers > 0) {
    plan(multisession, workers = req_workers)
} else {
    plan(multisession, workers = future::availableCores())
}

read_field_from_json <- function(path, field) {
    tryCatch(
        {
            js <- jsonlite::fromJSON(path, simplifyVector = TRUE)
            val <- js[[field]]
            if (is.null(val)) NA_real_ else as.numeric(val)
        },
        error = function(e) NA_real_
    )
}

cache_env <- new.env(parent = emptyenv())
get_cached_field <- function(path, field) {
    info <- suppressWarnings(file.info(path))
    key <- paste0(path, "||", info$mtime, "||", field)
    cached <- cache_env[[key]]
    if (!is.null(cached)) {
        return(cached)
    }
    val <- read_field_from_json(path, field)
    cache_env[[key]] <- val
    val
}

by_instance_method <- function(inst, method, field) {
    dir_path <- file.path(base_sim, inst, method)
    if (!dir.exists(dir_path)) {
        return(tibble(inst = inst, method = method, mean_val = NA_real_, n_files = 0L))
    }
    files <- Sys.glob(file.path(dir_path, "*.json"))
    files <- files[!grepl("/detail_", files)]
    if (!length(files)) {
        return(tibble(inst = inst, method = method, mean_val = NA_real_, n_files = 0L))
    }
    vals <- vapply(files, get_cached_field, numeric(1), field = field)
    vals <- vals[is.finite(vals)]
    n <- length(vals)
    tibble(
        inst = inst,
        method = method,
        mean_val = if (n == 0) NA_real_ else mean(vals),
        n_files = n
    )
}

grid <- expand.grid(inst = instances, method = methods, stringsAsFactors = FALSE)

message(sprintf("Computing %s over %d instances x %d methods...", field_to_average, length(instances), length(methods)))
detail_list <- future_lapply(
    seq_len(nrow(grid)),
    function(i) by_instance_method(grid$inst[i], grid$method[i], field_to_average),
    future.seed = TRUE
)
detail_df <- bind_rows(detail_list)

mean_df <- detail_df %>%
    filter(is.finite(mean_val), n_files > 0) %>%
    right_join(grid, by = c("inst", "method")) %>%
    left_join(inst_info, by = "inst") %>%
    filter(combo %in% combo_levels)

combo_summary <- mean_df %>%
    filter(is.finite(mean_val)) %>%
    group_by(combo, method) %>%
    summarise(mean_val = mean(mean_val), .groups = "drop")

if (!nrow(combo_summary)) {
    stop("No valid data found for the requested field/method combination.")
}

combo_summary <- combo_summary %>%
    mutate(
        combo = factor(combo, levels = combo_levels, ordered = TRUE),
        method = factor(method, levels = methods, ordered = TRUE)
    )

table_mat <- combos_interest %>%
    select(combo) %>%
    left_join(
        combo_summary %>%
            pivot_wider(
                id_cols = combo,
                names_from = method,
                values_from = mean_val
            ),
        by = "combo"
    ) %>%
    arrange(match(combo, combo_levels)) %>%
    as.data.frame()

numeric_table <- table_mat
original_table <- table_mat

decimals <- if (scale_100) 2 else 3
tol_tie <- 10^(-decimals - 6)
if (scale_100) {
    numeric_table[, -1] <- numeric_table[, -1] * 100
}
latex_numeric <- numeric_table
colnames(latex_numeric)[1] <- "G/C"
display_df <- as.data.frame(numeric_table, stringsAsFactors = FALSE)
colnames(display_df)[1] <- "G/C"

for (i in seq_len(nrow(display_df))) {
    row_vals <- display_df[i, -1, drop = FALSE]
    numeric_vals <- suppressWarnings(as.numeric(row_vals))
    rounded_vals <- round(numeric_vals, decimals)
    finite_idx <- which(is.finite(rounded_vals))
    if (!length(finite_idx)) {
        idx_min <- integer(0)
    } else {
        min_val <- min(rounded_vals[finite_idx], na.rm = TRUE)
        idx_min <- which(is.finite(rounded_vals) & abs(rounded_vals - min_val) < tol_tie)
    }
    for (j in seq_along(row_vals)) {
        val <- numeric_vals[j]
        if (is.na(val)) {
            display_df[i, j + 1] <- "--"
        } else {
            val_str <- format(round(val, decimals), nsmall = decimals)
            if (j %in% idx_min) {
                display_df[i, j + 1] <- paste0("\\textbf{", val_str, "}")
            } else {
                display_df[i, j + 1] <- val_str
            }
        }
    }
}

make_latex_table <- function(df, caption, label, methods) {
    tidy <- df %>%
        mutate(
            G = as.integer(sub("^G(\\d+)_C.*$", "\\1", `G/C`)),
            C = as.integer(sub("^G\\d+_C(\\d+)$", "\\1", `G/C`))
        ) %>%
        select(-`G/C`) %>%
        pivot_longer(
            cols = -c(G, C),
            names_to = "method",
            values_to = "value"
        )

    method_labels <- c(
        "ei.MD.bayes" = "ei.md.bayes",
        "nslphom_dual_w" = "nslphom\\_dual\\_w",
        "mvn_cdf_project_lp_FALSE" = "mvn\\_cdf",
        "mvn_pdf_project_lp_FALSE" = "mvn\\_pdf",
        "mult_project_lp_FALSE" = "mult",
        "exact_project_lp_FALSE" = "exact"
    )
    active_methods <- methods
    display_methods <- method_labels[active_methods]
    display_methods[is.na(display_methods)] <- active_methods[is.na(display_methods)]
    groups <- c(2, 3, 4, 6, 8)

    format_value <- function(val, is_min) {
        if (is.na(val)) {
            return("--")
        }
        str <- format(round(val, decimals), nsmall = decimals)
        if (is_min) paste0("\\textbf{", str, "}") else str
    }

    block_lines <- function(cands) {
        block_data <- tidy %>% filter(C %in% cands, method %in% active_methods)
        mins <- block_data %>%
            filter(is.finite(value)) %>%
            group_by(C, G) %>%
            summarise(min_val = min(round(value, decimals)), .groups = "drop")

        cell_val <- function(cand, grp, method) {
            val <- block_data %>%
                filter(C == cand, G == grp, method == !!method) %>%
                pull(value)
            val <- if (length(val)) val[1] else NA_real_
            min_val <- mins %>%
                filter(C == cand, G == grp) %>%
                pull(min_val)
            min_val <- if (length(min_val)) min_val[1] else NA_real_
            val_round <- round(val, decimals)
            format_value(val, is.finite(min_val) && isTRUE(abs(val_round - min_val) < tol_tie))
        }

        rows <- vapply(active_methods, function(m) {
            vals <- c(
                vapply(groups, function(g) cell_val(cands[1], g, m), character(1)),
                vapply(groups, function(g) cell_val(cands[2], g, m), character(1))
            )
            method_label <- if (!is.null(method_labels[[m]])) method_labels[[m]] else m
            paste0(method_label, " & ", paste(vals, collapse = " & "), " \\\\")
        }, character(1))

        c(
            "\\begin{tabularx}{\\textwidth}{l *{10}{>{\\centering\\arraybackslash}X}}",
            "\\toprule",
            paste0("Candidates & \\multicolumn{5}{c}{", cands[1], "} "),
            paste0("& \\multicolumn{5}{c}{", cands[2], "} \\\\"),
            "\\cmidrule(lr){2-6} \\cmidrule(lr){7-11}",
            "Groups & 2 & 3 & 4 & 6 & 8 & 2 & 3 & 4 & 6 & 8 \\\\",
            "\\midrule",
            rows,
            "\\bottomrule",
            "\\end{tabularx}"
        )
    }

    lines <- c(
        "\\begin{table}[H]",
        "\\centering",
        paste0("\\caption{", caption, "}"),
        paste0("\\label{", label, "}"),
        "",
        "\\begingroup",
        "\\setlength{\\tabcolsep}{4pt}",
        "",
        "% ======================================================",
        "% ===================== C = 2 and C = 3 =================",
        "% ======================================================",
        "",
        block_lines(c(2, 3)),
        "",
        "\\vspace{0.9em}",
        "",
        "% ======================================================",
        "% ===================== C = 5 and C = 10 ================",
        "% ======================================================",
        "",
        block_lines(c(5, 10)),
        "",
        "\\endgroup",
        "\\end{table}"
    )

    paste(lines, collapse = "\n")
}

unit_label <- if (scale_100) " (\\%)" else ""
caption_txt <- sprintf("Average \\texttt{%s}%s by G/C combination and method", field_to_average, unit_label)
label_txt <- sprintf("tab:sim_%s_gc", tolower(field_to_average))
latex_tbl <- make_latex_table(latex_numeric, caption_txt, label_txt, methods)

ascii_table <- (function(df) {
    cols <- colnames(df)
    df_chr <- as.data.frame(lapply(df, as.character), stringsAsFactors = FALSE)
    widths <- mapply(function(h, v) max(nchar(c(h, v), type = "width"), na.rm = TRUE), cols, df_chr)
    pad <- function(x, w) sprintf(paste0("%-", w, "s"), x)
    sep_line <- paste0("+", paste0(vapply(widths, function(w) strrep("-", w + 2), ""), collapse = "+"), "+")
    row_to_line <- function(r) paste0("| ", paste(mapply(pad, r, widths), collapse = " | "), " |")
    header <- row_to_line(cols)
    body <- apply(df_chr, 1, row_to_line)
    paste(c(sep_line, header, sep_line, body, sep_line), collapse = "\n")
})

cat(latex_tbl, sep = "\n")

if (!is.null(outfile)) {
    dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
    writeLines(latex_tbl, outfile)
    message(sprintf("[OK] LaTeX table saved to: %s", normalizePath(outfile, mustWork = FALSE)))
}

cat("\n========== TABLE (ASCII pretty print) ==========\n")
cat(ascii_table(display_df), sep = "\n")
cat("\n================================================\n")

xlsx_path <- if (!is.null(outfile)) {
    sub("\\.tex$", ".xlsx", outfile)
} else {
    paste0("table_simulated_", tolower(field_to_average), ".xlsx")
}

save_xlsx <- function(df_matrix, path_xlsx) {
    ok <- requireNamespace("openxlsx", quietly = TRUE)
    if (ok) {
        wb <- openxlsx::createWorkbook()
        openxlsx::addWorksheet(wb, "Summary")
        openxlsx::writeData(wb, "Summary", df_matrix)
        openxlsx::freezePane(wb, "Summary", firstRow = TRUE, firstCol = TRUE)
        dir.create(dirname(path_xlsx), recursive = TRUE, showWarnings = FALSE)
        openxlsx::saveWorkbook(wb, path_xlsx, overwrite = TRUE)
    } else if (requireNamespace("writexl", quietly = TRUE)) {
        dir.create(dirname(path_xlsx), recursive = TRUE, showWarnings = FALSE)
        writexl::write_xlsx(list(Summary = df_matrix), path = path_xlsx)
    } else {
        warning("Install 'openxlsx' or 'writexl' to export XLSX tables.")
    }
}

numeric_summary <- original_table
colnames(numeric_summary)[1] <- "G/C"
save_xlsx(numeric_summary, xlsx_path)
message(sprintf("[OK] XLSX written to: %s", normalizePath(xlsx_path, mustWork = FALSE)))
