#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(stringr)
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(future.apply)
    library(jsonlite)
})

# =========================================
# Parse arguments
# =========================================
args <- commandArgs(trailingOnly = TRUE)

parse_opt_int <- function(pos, default = NA_integer_) {
    if (length(args) >= pos) suppressWarnings(as.integer(args[pos])) else default
}
parse_opt_chr <- function(pos, default = NULL) {
    if (length(args) >= pos) args[pos] else default
}

# Usage (positional):
#   Rscript src/table.R <field_to_average> [outfile.tex] [workers] [--inst-like=ei_] [--limit-inst=9999]
# Examples:
#   Rscript src/table.R EI_Z output/table1.tex
#   Rscript src/table.R runtime output/table_runtime.tex 8 --inst-like=ei_NZ_ --limit-inst=100

field_to_average <- {
    ft <- parse_opt_chr(1, NULL)
    if (is.null(ft) || is.na(ft)) {
        stop("Missing <field_to_average> argument. Example: Rscript src/table.R EI_Z output/table1.tex")
    }
    ft
}

outfile <- if (!is.na(tmp <- parse_opt_chr(2))) tmp else NULL
req_workers <- parse_opt_int(3, NA_integer_)

inst_like <- {
    m <- grep("^--inst-like=", args, value = TRUE)
    if (length(m)) sub("^--inst-like=", "", m[1]) else "ei_"
}
limit_inst <- {
    m <- grep("^--limit-inst=", args, value = TRUE)
    if (length(m)) suppressWarnings(as.integer(sub("^--limit-inst=", "", m[1]))) else NA_integer_
}

# =========================================
# Config
# =========================================
base_ei <- file.path("output", "ei_instances")

methods <- c(
    # "lclphom",
    # "lphom",
    # "mult_project_lp_FALSE",
    "mvn_pdf_project_lp_FALSE_sym"
    # "mvn_pdf_project_lp_FALSE",
    # "mult_project_lp_FALSE_sym"
    # "mvn_pdf_project_lp_FALSE_sym"
    # "mvn_pdf_project_lp_FALSE_sym"
    # "ei.md.bayes",
    # "nslphom_joint",
    # "lphom_joint",
    # "nslphom_dual_w",
    # "nslphom_dual_a",
    # "mult_project_lp_FALSE_inv",
    # "mvn_cdf_project_lp_FALSE_inv",
    # "mvn_pdf_project_lp_FALSE_inv"
    # "nslphom_dual_w",
    # "ecolRxC"
)

if (!dir.exists(base_ei)) {
    stop(sprintf("'%s' doesn't exist. Run this script from the project root.", base_ei))
}

# --- NEW: sólo errar si existe 'ei_instances' dentro de base_ei ---
all_instances <- sort(list.dirs(base_ei, recursive = FALSE, full.names = FALSE))
has_ei_instances <- dir.exists(file.path(base_ei, "ei_instances"))

instances <- all_instances[str_detect(all_instances, paste0("^", inst_like))]
if (length(instances) == 0) {
    if (has_ei_instances) {
        stop(sprintf("No subfolders found starting with '%s' in '%s'.", inst_like, base_ei))
    } else {
        message(sprintf(
            "No subfolders starting with '%s', and '%s' has no 'ei_instances' folder; continuing without error.",
            inst_like, base_ei
        ))
    }
}
if (is.finite(limit_inst) && limit_inst > 0) {
    instances <- head(instances, limit_inst)
}

# =========================================
# Parallelization
# =========================================
if (is.finite(req_workers) && req_workers > 0) {
    plan(multisession, workers = req_workers)
} else {
    plan(multisession, workers = future::availableCores())
}

# =========================================
# Generic JSON field reader + cache
# =========================================
read_field_from_json <- function(f, field) {
    # Reads JSON and returns the specified field as numeric or NA
    tryCatch(
        {
            js <- jsonlite::fromJSON(f, simplifyVector = TRUE)
            v <- js[[field]]
            if (is.null(v)) NA_real_ else as.numeric(v)
        },
        error = function(e) NA_real_
    )
}

cache_env <- new.env(parent = emptyenv())
get_cached_field <- function(f, field) {
    fi <- suppressWarnings(file.info(f))
    key <- paste0(f, "||", fi$mtime, "||", field)
    got <- cache_env[[key]]
    if (!is.null(got)) {
        return(got)
    }
    v <- read_field_from_json(f, field)
    cache_env[[key]] <- v
    v
}

# =========================================
# Helpers
# =========================================
# Average by district (for xlsx sheets)
by_district_for <- function(inst, method, field) {
    base_inst <- file.path(base_ei, inst)
    districts <- sort(list.dirs(base_inst, recursive = FALSE, full.names = FALSE))

    if (length(districts) == 0) {
        return(tibble(
            inst = character(0),
            method = character(0),
            district = character(0),
            mean_val = double(0),
            n_files = integer(0)
        ))
    }

    rows <- lapply(districts, function(d) {
        files <- Sys.glob(file.path(base_inst, d, method, "*.json"))
        files <- files[!grepl("/detail_", files)] # exclude detail jsons
        if (length(files) == 0) {
            return(tibble(
                inst = inst, method = method, district = d,
                mean_val = NA_real_, n_files = 0L
            ))
        }
        vals <- vapply(files, get_cached_field, numeric(1), field = field)
        vals <- vals[is.finite(vals)]
        n <- length(vals)
        mean_val <- if (n == 0) NA_real_ else mean(vals)
        tibble(
            inst = inst, method = method, district = d,
            mean_val = mean_val, n_files = as.integer(n)
        )
    })

    bind_rows(rows)
}

# =========================================
# Grid (instance x method)
# =========================================
grid <- expand.grid(inst = instances, method = methods, stringsAsFactors = FALSE)
# =========================================
# Parallel calculation (detail by district + global mean)
# =========================================
detail_list <- future_lapply(
    seq_len(nrow(grid)),
    function(i) by_district_for(grid$inst[i], grid$method[i], field_to_average),
    future.seed = TRUE
)
detail_df <- bind_rows(detail_list)

# Global mean per cell with equal weight per district
mean_df <- detail_df %>%
    filter(is.finite(mean_val), n_files > 0) %>%
    group_by(inst, method) %>%
    summarise(mean_val = mean(mean_val), .groups = "drop") %>%
    right_join(grid, by = c("inst", "method")) %>%
    select(inst, method, mean_val)

# =========================================
# Matrix of means (rows = method, columns = instance)
# =========================================

# --- Definición del orden deseado ---
# order_methods <- c(
#     "lclphom",
#     "lphom",
#     "mult_project_lp_FALSE",
#     "mvn_cdf_project_lp_FALSE",
#     "mvn_pdf_project_lp_FALSE",
#     "nslphom_dual_a",
#     "nslphom_dual_w",
#     "nslphom_joint",
#     "ei.md.bayes",
#     "ecolRxC"
# )
order_methods <- c(
    # "mvn_pdf_project_lp_FALSE_sym",
    # "mvn_pdf_project_lp_FALSE"
    "mvn_pdf_project_lp_FALSE_sym"
    # "mult_project_lp_FALSE_sym"
)



# --- Aplicar el orden antes de pivotear ---
mean_df <- mean_df %>%
    mutate(method = factor(method, levels = order_methods, ordered = TRUE))

mat <- mean_df %>%
    tidyr::pivot_wider(
        id_cols = method,
        names_from = inst,
        values_from = mean_val
    ) %>%
    arrange(method) %>%
    mutate(method = as.character(method)) %>%
    as.data.frame()

# ====== pesos por instancia: # de distritos con datos ======
w_by_inst <- detail_df %>%
    dplyr::filter(is.finite(mean_val), n_files > 0) %>%
    dplyr::group_by(inst) %>%
    dplyr::summarise(n_districts = dplyr::n_distinct(district), .groups = "drop")

inst_cols <- setdiff(colnames(mat), "method")

# alinear pesos al orden de columnas de la tabla
w <- w_by_inst$n_districts[match(inst_cols, w_by_inst$inst)]
# si alguna instancia no aparece en w_by_inst, ponemos 0 para que no pese
w[is.na(w)] <- 0

# ====== promedio ponderado por método (fila) ======
weighted_mean_row <- function(x, w) {
    mask <- !is.na(x) & (w > 0)
    if (!any(mask)) {
        return(NA_real_)
    }
    sum(x[mask] * w[mask]) / sum(w[mask])
}

mat$mean <- apply(as.matrix(mat[inst_cols]), 1, weighted_mean_row, w = w)

# =========================================
# LaTeX with bold MINIMUMS per column
# (Assumes "lower is better": works for EI_* and runtime)
# =========================================
decimals <- 3
df_fmt <- mat
colnames(df_fmt)[1] <- "Method"

for (j in 2:ncol(df_fmt)) {
    col_vals <- df_fmt[[j]]
    if (all(is.na(col_vals))) {
        df_fmt[[j]] <- ""
        next
    }
    idx_min <- which(col_vals == min(col_vals, na.rm = TRUE))
    df_fmt[[j]] <- vapply(seq_along(col_vals), function(i) {
        if (is.na(col_vals[i])) {
            return("")
        }
        val_str <- format(round(col_vals[i], decimals), nsmall = decimals)
        if (i %in% idx_min) {
            paste0("\\textbf{", val_str, "}")
        } else {
            val_str
        }
    }, character(1))
}
make_latex_table <- function(df, caption = "Average by method and instance",
                             label = "tab:mean_table") {
    cols <- colnames(df)
    align <- paste0("l", paste(rep("c", length(cols) - 1), collapse = ""))
    header <- paste(cols, collapse = " & ")
    lines <- apply(df, 1, function(r) paste(r, collapse = " & "))
    body <- paste0("    ", lines, " \\\\")
    paste(
        "\\begin{table}[H]",
        "  \\centering",
        paste0("  \\caption{", caption, "}"),
        paste0("  \\label{", label, "}"),
        paste0("  \\begin{tabular}{", align, "}"),
        "    \\hline",
        paste0("    ", header, " \\\\"),
        "    \\hline",
        paste(body, collapse = "\n"),
        "    \\hline",
        "  \\end{tabular}",
        "\\end{table}",
        sep = "\n"
    )
}

caption_txt <- sprintf("Average \\texttt{%s} by method and instance", field_to_average)
label_txt <- sprintf("tab:%s_mean", tolower(field_to_average))
latex_tbl <- make_latex_table(df_fmt, caption = caption_txt, label = label_txt)

# =========================================
# Pretty print
# =========================================
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

# =========================================
# LaTeX + ASCII output
# =========================================
cat(latex_tbl, sep = "\n")

if (!is.null(outfile)) {
    dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
    writeLines(latex_tbl, outfile)
    message(sprintf("\n[OK] LaTeX table saved to: %s", normalizePath(outfile, mustWork = FALSE)))
}

cat("\n========== TABLE (ASCII pretty print) ==========\n")
cat(ascii_table(df_fmt), sep = "\n")
cat("\n================================================\n")

# =========================================
# Export XLSX (one sheet per instance + summary)
# =========================================
xlsx_path <- if (!is.null(outfile)) {
    sub("\\.tex$", ".xlsx", outfile)
} else {
    paste0("mean_", tolower(field_to_average), ".xlsx")
}

save_xlsx <- function(detail_df, mat_means, path_xlsx) {
    ok <- requireNamespace("openxlsx", quietly = TRUE)
    if (ok) {
        wb <- openxlsx::createWorkbook()
        for (inst in unique(detail_df$inst)) {
            df_inst <- detail_df %>%
                filter(inst == !!inst) %>%
                select(district, method, mean_val) %>%
                mutate(method = factor(method, levels = methods)) %>%
                arrange(district, method) %>%
                tidyr::pivot_wider(names_from = method, values_from = mean_val)
            openxlsx::addWorksheet(wb, sheetName = substr(inst, 1, 31))
            openxlsx::writeData(wb, sheet = substr(inst, 1, 31), x = df_inst)
            openxlsx::freezePane(wb, substr(inst, 1, 31), firstRow = TRUE, firstCol = TRUE)
        }
        df_means <- mat_means
        colnames(df_means)[1] <- "Method"
        openxlsx::addWorksheet(wb, sheetName = "Summary_means")
        openxlsx::writeData(wb, "Summary_means", df_means)
        openxlsx::freezePane(wb, "Summary_means", firstRow = TRUE, firstCol = TRUE)
        dir.create(dirname(path_xlsx), recursive = TRUE, showWarnings = FALSE)
        openxlsx::saveWorkbook(wb, file = path_xlsx, overwrite = TRUE)
    } else {
        if (!requireNamespace("writexl", quietly = TRUE)) {
            stop("You need 'openxlsx' or 'writexl' to export to XLSX. Please install one of them.")
        }
        sheets <- list()
        for (inst in unique(detail_df$inst)) {
            df_inst <- detail_df %>%
                filter(inst == !!inst) %>%
                select(district, method, mean_val) %>%
                mutate(method = factor(method, levels = methods)) %>%
                arrange(district, method) %>%
                tidyr::pivot_wider(names_from = method, values_from = mean_val)
            sheets[[inst]] <- df_inst
        }
        df_means <- mat_means
        colnames(df_means)[1] <- "Method"
        sheets[["Summary_means"]] <- df_means
        dir.create(dirname(path_xlsx), recursive = TRUE, showWarnings = FALSE)
        writexl::write_xlsx(sheets, path = path_xlsx)
    }
}

means_numeric <- mat
colnames(means_numeric)[1] <- "Method"

save_xlsx(detail_df, means_numeric, xlsx_path)
message(sprintf("[OK] XLSX written to: %s", normalizePath(xlsx_path, mustWork = FALSE)))
