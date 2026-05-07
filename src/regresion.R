#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(dplyr)
    library(jsonlite)
    library(purrr)
    library(stargazer)
    library(stringr)
    library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
Sys.setenv(PATH = paste(
    unique(c("/bin", "/usr/bin", "/usr/local/bin", "/opt/homebrew/bin", Sys.getenv("PATH"))),
    collapse = ":"
))

parse_arg <- function(flag, default = NULL) {
    hit <- grep(paste0("^--", flag, "="), args, value = TRUE)
    if (!length(hit)) {
        return(default)
    }
    sub(paste0("^--", flag, "="), "", hit[1])
}

parse_flag <- function(flag, default = FALSE) {
    hit <- grep(paste0("^--", flag, "(=|$)"), args, value = TRUE)
    if (!length(hit)) {
        return(default)
    }
    val <- sub(paste0("^--", flag, "=?"), "", hit[1])
    if (!nzchar(val)) {
        return(TRUE)
    }
    tolower(val) %in% c("1", "true", "t", "yes", "y")
}

out_tex <- parse_arg("out", file.path("output", "regresion_table.tex"))
use_common_sample <- parse_flag("common-sample", TRUE)
exclude_maori <- parse_flag("exclude-maori", FALSE)
include_nslphom <- parse_flag("include-nslphom", FALSE)

base_root <- file.path("output", "ei_instances")
if (!dir.exists(base_root)) {
    stop(sprintf("Directory '%s' not found. Run the script from the project root.", base_root))
}

maori_districts <- c(
    "Hauraki_Waikato",
    "Ikaroa_Rawhiti",
    "Tamaki_Makaurau",
    "Te_Tai_Hauauru",
    "Te_Tai_Tokerau",
    "Te_Tai_Tonga",
    "Waiariki"
)

method_candidates <- list(
    mult = c("mult_project_lp_FALSE_pt0p001_llminf"),
    mult_sym = c("mult_lp_FALSE_sym_EM_weight_pt0p001_llminf"),
    mvn_pdf = c("mvn_pdf_project_lp_FALSE_pt0p001_llminf"),
    mvn_pdf_sym = c("mvn_pdf_lp_FALSE_sym_EM_weight_pt0p001_llminf"),
    mvn_cdf = c("mvn_cdf_project_lp_FALSE_pt0p001_llminf"),
    mvn_cdf_sym = c("mvn_cdf_lp_FALSE_sym_EM_weight_pt0p001_llminf"),
    nslphom_dual_w = c("nslphom_dual_w"),
    nslphom_joint = c("nslphom_joint")
)

base_predictor_names <- c(
    "ratio_ijk",
    "log_sum_voters",
    "hete",
    "chi2e",
    "dwr_dwc_sum",
    # "dwr_dwc_log_gap",
    "var_vac_sum"
    # "var_vac_log_gap"
)

predictor_names <- c(
    "ratio_ijk",
    "ratio_ijk_sq",
    "log_sum_voters",
    "hete",
    "chi2e",
    "chi2e_sq",
    "dwr_dwc_sum",
    # "dwr_dwc_log_gap",
    "var_vac_sum"
    # "var_vac_log_gap"
    # "log_sum_voters_sq"
)

quadratic_predictor_map <- c(
    ratio_ijk_sq = "ratio_ijk",
    chi2e_sq = "chi2e"
)

table_method_order <- if (include_nslphom) {
    c("mvn_pdf_sym", "mvn_cdf_sym", "mult_sym", "nslphom_dual_w", "nslphom_joint")
} else {
    c("mvn_pdf_sym", "mvn_cdf_sym", "mult_sym")
}

z_score <- function(x) {
    mu <- mean(x, na.rm = TRUE)
    sig <- stats::sd(x, na.rm = TRUE)
    if (!is.finite(sig) || sig <= 0) {
        return(rep(0, length(x)))
    }
    (x - mu) / sig
}

canonical_district_name <- function(x) {
    sub("_[Gg][0-9]+_C[0-9]+$", "", basename(x))
}

pick_method_dir <- function(case_dir, candidates) {
    direct_hits <- file.path(case_dir, candidates)
    direct_hits <- direct_hits[dir.exists(direct_hits)]
    if (length(direct_hits)) {
        return(direct_hits[[1]])
    }

    subdirs <- list.dirs(case_dir, recursive = FALSE, full.names = TRUE)
    if (!length(subdirs)) {
        return(NULL)
    }

    idx <- match(candidates, basename(subdirs), nomatch = 0L)
    idx <- idx[idx > 0L]
    if (!length(idx)) {
        return(NULL)
    }
    subdirs[[idx[[1]]]]
}

pick_json_file <- function(method_dir, detail = FALSE) {
    files <- list.files(method_dir, pattern = "\\.json$", full.names = TRUE)
    if (!length(files)) {
        return(NULL)
    }

    if (detail) {
        files <- files[grepl("(^|/)detail_.*\\.json$", files)]
        preferred <- c("detail_1.json", "detail_1_pinit123.json")
    } else {
        files <- files[!grepl("(^|/)detail_.*\\.json$", files)]
        files <- files[!grepl("(^|/)metrics_.*\\.json$", files)]
        preferred <- c("1.json", "1_pinit123.json")
    }

    if (!length(files)) {
        return(NULL)
    }

    hits <- files[match(preferred, basename(files), nomatch = 0L)]
    hits <- hits[nzchar(hits)]
    if (length(hits)) {
        return(hits[[1]])
    }

    info <- file.info(files)
    files[[order(info$mtime, decreasing = TRUE)[[1]]]]
}

pick_case_detail_file <- function(case_dir) {
    method_dirs <- list.dirs(case_dir, recursive = FALSE, full.names = TRUE)
    if (!length(method_dirs)) {
        return(NULL)
    }

    detail_files <- purrr::map_chr(method_dirs, function(method_dir) {
        out <- pick_json_file(method_dir, detail = TRUE)
        if (is.null(out)) "" else out
    })
    detail_files <- detail_files[nzchar(detail_files)]
    if (!length(detail_files)) {
        return(NULL)
    }

    preferred <- c("detail_1.json", "detail_1_pinit123.json")
    hits <- detail_files[match(preferred, basename(detail_files), nomatch = 0L)]
    hits <- hits[nzchar(hits)]
    if (length(hits)) {
        return(hits[[1]])
    }

    info <- file.info(detail_files)
    detail_files[[order(info$mtime, decreasing = TRUE)[[1]]]]
}

# Normalize JSON-derived objects into a plain numeric matrix.
# The same conceptual object can be parsed by jsonlite as a matrix, a data
# frame, or a nested list, depending on file structure. All regression metrics
# below are defined in matrix notation, so this helper makes their inputs
# uniform before any calculations happen.
#
# In this script, matrix roles follow the fastei naming convention:
#   X = precinct-by-candidate matrix
#   W = precinct-by-demographic-group matrix
#   V_est = estimated aggregate group-by-candidate table
as_num_matrix <- function(x) {
    if (is.null(x) || !length(x)) {
        return(matrix(numeric(), nrow = 0, ncol = 0))
    }

    if (is.matrix(x)) {
        storage.mode(x) <- "double"
        return(x)
    }

    if (is.data.frame(x)) {
        return(as.matrix(data.frame(lapply(x, as.numeric), check.names = FALSE)))
    }

    if (!is.list(x)) {
        stop("Could not convert the object to a numeric matrix.")
    }

    rows <- lapply(x, function(row) as.numeric(unlist(row, recursive = TRUE, use.names = FALSE)))
    do.call(rbind, rows)
}

# Rebuild the estimated internal cell array Z_est and orient it to the shape
# expected by the metric definitions.
# Conceptually, Z_est[j, k, i] is the estimated number of voters in precinct i
# who belong to demographic group j and vote for candidate k. Because the JSON
# nesting order is not guaranteed to match that convention, we try all six
# permutations of the raw dimensions and keep the one that matches c(J, K, I).
# This is especially important for HETe, which is defined from precinct-level
# jk tables rather than only from margins.
as_num_array3 <- function(x, target_dims) {
    if (is.null(x) || !length(x)) {
        stop("Z_est is empty or missing.")
    }

    raw_dims <- c(length(x), length(x[[1]]), length(x[[1]][[1]]))
    arr <- array(
        as.numeric(unlist(x, recursive = TRUE, use.names = FALSE)),
        dim = raw_dims
    )

    permutations <- list(
        c(1L, 2L, 3L), c(1L, 3L, 2L), c(2L, 1L, 3L),
        c(2L, 3L, 1L), c(3L, 1L, 2L), c(3L, 2L, 1L)
    )

    for (perm in permutations) {
        if (identical(raw_dims[perm], target_dims)) {
            return(aperm(arr, perm))
        }
    }

    stop(sprintf(
        "Could not orient Z_est. Found dimensions: [%s]; expected: [%s].",
        paste(raw_dims, collapse = ", "),
        paste(target_dims, collapse = ", ")
    ))
}

# Compute average component variance of precinct-level compositions.
# If mat = W, this returns VAR: heterogeneity of demographic compositions across
# precincts. If mat = X, this returns VAC: heterogeneity of candidate vote
# compositions across precincts.
#
# Each precinct row is first normalized to shares that sum to one. Then, for
# each component/column, we compute the variance of that share across precincts.
# The final statistic is the average of those component-wise variances.
# Under this definition, two precincts with different sizes but identical
# compositions contribute zero dispersion.
composition_total_variance <- function(mat) {
    mat <- as_num_matrix(mat)
    if (!nrow(mat) || ncol(mat) < 2) {
        return(NA_real_)
    }

    row_totals <- rowSums(mat, na.rm = TRUE)
    valid_rows <- is.finite(row_totals) & row_totals > 0
    if (sum(valid_rows) < 2) {
        return(NA_real_)
    }

    I <- sum(mat)
    # print(I)

    shares <- mat[valid_rows, , drop = FALSE] / row_totals[valid_rows]
    component_vars <- apply(shares, 2L, stats::sd, na.rm = TRUE)

    weighted <- sum(component_vars * colSums(shares, na.rm = TRUE), na.rm = TRUE)
    # weighted <- sum(component_vars, na.rm = TRUE)
    # sum(component_vars, na.rm = TRUE)
    # return(weighted)
}

# DWR: standard deviation of aggregate row shares.
# Using W, we sum each demographic group over all precincts, convert the totals
# into shares of the full population, and take the standard deviation across
# groups. This is a simple imbalance measure for row margins in the aggregate
# table: it increases when some groups are much larger than others.
compute_dwr <- function(w_mat) {
    w_mat <- as_num_matrix(w_mat)
    totals <- colSums(w_mat, na.rm = TRUE)
    stats::sd(totals / sum(totals), na.rm = TRUE)
}

# DWC: standard deviation of aggregate column shares.
# This is the column-side analogue of DWR. Using X, we aggregate candidate vote
# totals across precincts, convert them into vote shares, and compute the
# standard deviation across candidates. It is larger when the election is more
# concentrated on a subset of candidates and smaller when candidate totals are
# more evenly balanced.
compute_dwc <- function(x_mat) {
    x_mat <- as_num_matrix(x_mat)
    totals <- colSums(x_mat, na.rm = TRUE)
    stats::sd(totals / sum(totals), na.rm = TRUE)
}

max_log_ratio_gap <- function(a, b, pseudocount = 0.5) {
    if (!is.finite(a) || !is.finite(b)) {
        return(NA_real_)
    }

    max(
        log((a + pseudocount) / (b + pseudocount)),
        log((b + pseudocount) / (a + pseudocount))
    )
}

# chi2e: Pearson chi-square statistic for the estimated aggregate table, scaled
# by its degrees of freedom.
# Starting from V_est, we form the independence benchmark implied by row and
# column margins:
#   E_jk = row_total_j * col_total_k / grand_total.
# We then sum ((V_est_jk - E_jk)^2 / E_jk) over valid cells and divide by
# (J - 1)(K - 1). The metric is near zero when the estimated aggregate table is
# close to independence and grows as the estimated association between
# demographic groups and candidates becomes stronger.
compute_chi2e <- function(v_est) {
    v_est <- as_num_matrix(v_est)
    j <- nrow(v_est)
    k <- ncol(v_est)
    if (j < 2 || k < 2) {
        return(NA_real_)
    }

    expected <- outer(rowSums(v_est), colSums(v_est)) / sum(v_est)
    ok <- is.finite(expected) & expected > 0
    sum((v_est[ok] - expected[ok])^2 / expected[ok], na.rm = TRUE) / ((j - 1) * (k - 1))
}

# HETe: precinct-level heterogeneity of the estimated internal tables.
# The idea is to compare each local estimated jk table against two benchmark
# structures implied by the aggregate estimate V_est.
#
# First, p_hat stores candidate shares within each demographic group, and q_hat
# stores demographic shares within each candidate, both derived from V_est.
# For each precinct i:
#   1. expected_from_rows rescales p_hat using that precinct's demographic
#      totals W[i, ].
#   2. expected_from_cols rescales q_hat using that precinct's candidate totals
#      X[i, ].
#   3. We measure the absolute deviation of the local Z_est table from both
#      benchmark tables.
#
# After summing those deviations over all precincts and all jk cells, we divide
# by 4 times the total number of votes. A larger HETe means the local precinct
# tables differ more from the common association pattern suggested by the
# aggregate table, i.e. greater within-instance heterogeneity.
compute_hete <- function(z_est, v_est, w_mat, x_mat) {
    v_est <- as_num_matrix(v_est)
    w_mat <- as_num_matrix(w_mat)
    x_mat <- as_num_matrix(x_mat)

    j <- nrow(v_est)
    k <- ncol(v_est)
    i_units <- nrow(w_mat)

    z_arr <- as_num_array3(z_est, target_dims = c(j, k, i_units))

    p_hat <- sweep(v_est, 1, rowSums(v_est), "/")
    q_hat <- t(sweep(v_est, 2, colSums(v_est), "/"))

    acc <- 0
    for (i in seq_len(i_units)) {
        local_votes <- z_arr[, , i, drop = FALSE][, , 1]
        expected_from_rows <- p_hat * w_mat[i, ]
        expected_from_cols <- t(q_hat * x_mat[i, ])
        acc <- acc +
            sum(abs(local_votes - expected_from_rows), na.rm = TRUE) +
            sum(abs(local_votes - expected_from_cols), na.rm = TRUE)
    }

    acc / (4 * sum(v_est))
}

compute_logvoters <- function(w_mat, total_units) {
    total_voters <- sum(as_num_matrix(w_mat), na.rm = TRUE)
    if (!is.finite(total_voters) || total_voters <= 0) {
        return(NA_real_)
    }
    log(total_voters / total_units)
}

safe_read_json <- function(path) {
    tryCatch(
        jsonlite::fromJSON(path, simplifyVector = FALSE),
        error = function(e) NULL
    )
}

collect_one_case <- function(base_dir, district_dir, method_label, candidates) {
    method_dir <- pick_method_dir(district_dir, candidates)
    if (is.null(method_dir)) {
        return(NULL)
    }

    summary_file <- pick_json_file(method_dir, detail = FALSE)
    detail_file <- pick_json_file(method_dir, detail = TRUE)
    if (is.null(detail_file)) {
        detail_file <- pick_case_detail_file(district_dir)
    }
    if (is.null(summary_file) || is.null(detail_file)) {
        return(NULL)
    }

    summary_js <- safe_read_json(summary_file)
    detail_js <- safe_read_json(detail_file)
    if (is.null(summary_js) || is.null(detail_js)) {
        return(NULL)
    }

    if (is.null(summary_js$EI_V) || is.null(summary_js$V_est) || is.null(summary_js$Z_est) ||
        is.null(detail_js$X) || is.null(detail_js$W)) {
        return(NULL)
    }

    x_mat <- as_num_matrix(detail_js$X)
    w_mat <- as_num_matrix(detail_js$W)
    v_est <- as_num_matrix(summary_js$V_est)

    if (!nrow(x_mat) || !nrow(w_mat) || !nrow(v_est)) {
        return(NULL)
    }

    j <- ncol(w_mat)
    k <- ncol(x_mat)
    ratio_ijk <- nrow(x_mat) / (j * k) # Here we obtain B / CG
    chi2e <- compute_chi2e(v_est)
    var_w <- composition_total_variance(w_mat)
    var_x <- composition_total_variance(x_mat)
    dwr <- compute_dwr(w_mat)
    dwc <- compute_dwc(x_mat)

    tibble(
        method = method_label,
        base = basename(base_dir),
        district_folder = basename(district_dir),
        district = canonical_district_name(district_dir),
        case_id = paste(basename(base_dir), canonical_district_name(district_dir), sep = "::"),
        is_maori = canonical_district_name(district_dir) %in% maori_districts,
        EI = as.numeric(summary_js$EI_V[[1]]),
        ratio_ijk = ratio_ijk,
        hete = compute_hete(summary_js$Z_est, v_est, w_mat = w_mat, x_mat = x_mat),
        chi2e = chi2e,
        var_w = var_w,
        var_x = var_x,
        dwr = dwr,
        dwc = dwc,
        log_sum_voters = compute_logvoters(w_mat, nrow(w_mat)),
        # log_sum_voters_sq = compute_logvoters(w_mat, nrow(w_mat) * j * k)^2,
        dwr_dwc_sum = dwr + dwc,
        # dwr_dwc_log_gap = max_log_ratio_gap(dwr, dwc),
        var_vac_sum = var_w + var_x
        # var_vac_log_gap = max_log_ratio_gap(var_w, var_x)
    )
}

base_dirs <- list.dirs(base_root, recursive = FALSE, full.names = TRUE)
district_dirs <- unlist(lapply(base_dirs, function(base_dir) {
    list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
}), use.names = FALSE)

raw_rows <- purrr::imap_dfr(method_candidates, function(candidates, method_label) {
    purrr::map_dfr(district_dirs, function(district_dir) {
        base_dir <- dirname(district_dir)
        collect_one_case(base_dir, district_dir, method_label, candidates)
    })
})

if (!nrow(raw_rows)) {
    stop("Could not build rows for any method.")
}

if (exclude_maori) {
    raw_rows <- raw_rows %>% filter(!is_maori)
}

analysis_rows <- raw_rows %>%
    filter(if_all(all_of(c("EI", base_predictor_names)), ~ is.finite(.x)))

if (!nrow(analysis_rows)) {
    stop("No observations remained after filtering missing values.")
}

if (use_common_sample) {
    sample_rows <- analysis_rows %>%
        filter(method %in% table_method_order)
    if (!nrow(sample_rows)) {
        stop("No observations for the required methods in the final table.")
    }
    ids_by_method <- split(sample_rows$case_id, sample_rows$method)
    common_ids <- Reduce(intersect, ids_by_method)
    if (length(common_ids)) {
        analysis_rows <- analysis_rows %>% filter(case_id %in% common_ids)
    } else {
        warning(sprintf(
            "No common sample across the required methods (%s); using each method's available sample.",
            paste(table_method_order, collapse = ", ")
        ))
    }
}

prepare_predictors <- function(df_method) {
    out <- df_method %>%
        mutate(across(all_of(base_predictor_names), z_score))

    for (sq_name in names(quadratic_predictor_map)) {
        base_name <- quadratic_predictor_map[[sq_name]]
        out[[sq_name]] <- out[[base_name]]^2
    }

    out
}

table_rows <- analysis_rows %>%
    filter(method %in% table_method_order)

if (length(unique(table_rows$method)) != length(table_method_order)) {
    stop(sprintf(
        "Not all required methods are available for the final table. Required: %s",
        paste(table_method_order, collapse = ", ")
    ))
}

predictor_base <- table_rows %>%
    filter(method == table_method_order[[1]]) %>%
    select(case_id, base, all_of(base_predictor_names)) %>%
    prepare_predictors() %>%
    select(case_id, base, all_of(predictor_names))

ei_wide <- reshape(
    as.data.frame(table_rows[, c("case_id", "method", "EI")]),
    idvar = "case_id",
    timevar = "method",
    direction = "wide"
)

required_ei_cols <- paste0("EI.", table_method_order)
if (!all(required_ei_cols %in% names(ei_wide))) {
    stop(sprintf(
        "Could not build the expected wide EI columns. Available columns: %s",
        paste(names(ei_wide), collapse = ", ")
    ))
}

wide_data <- predictor_base %>%
    left_join(as_tibble(ei_wide), by = "case_id")

for (method_name in table_method_order) {
    wide_data[[method_name]] <- wide_data[[paste0("EI.", method_name)]]
    wide_data[[paste0("EI.", method_name)]] <- NULL
}

if (include_nslphom) {
    response_order <- c(
        "mvn_pdf_sym",
        "mvn_cdf_sym",
        "mult_sym",
        "nslphom_dual_w",
        "nslphom_joint"
    )
    response_labels <- response_order
} else {
    wide_data <- wide_data %>%
        mutate(
            diff_mult_pdf = mult_sym - mvn_pdf_sym,
            diff_cdf_pdf = mvn_cdf_sym - mvn_pdf_sym,
            diff_mult_cdf = mult_sym - mvn_cdf_sym
        )

    response_order <- c(
        "mvn_pdf_sym",
        "mvn_cdf_sym",
        "mult_sym",
        "diff_mult_pdf",
        "diff_cdf_pdf",
        "diff_mult_cdf"
    )

    response_labels <- c(
        "mvn_pdf_sym",
        "mvn_cdf_sym",
        "mult_sym",
        "mult_sym - mvn_pdf_sym",
        "mvn_cdf_sym - mvn_pdf_sym",
        "mult_sym - mvn_cdf_sym"
    )
}

build_model_frame <- function(df, response_name) {
    out <- df[, c(response_name, predictor_names, "base")]
    names(out)[[1]] <- "EI"
    out
}

render_correlation_matrix <- function(df) {
    corr_mat <- stats::cor(df[, predictor_names], use = "pairwise.complete.obs")
    corr_df <- data.frame(variable = rownames(corr_mat), corr_mat, check.names = FALSE)
    rownames(corr_df) <- NULL
    capture.output(print(corr_df, digits = 3, row.names = FALSE))
}

# formula_table5 <- EI ~ ratio_ijk + ratio_ijk_sq + log_sum_voters + hete + chi2e +
#    chi2e_sq + dwr_dwc_sum + dwr_dwc_log_gap + var_vac_sum + var_vac_log_gap

# formula_table5 <- EI ~ ratio_ijk + ratio_ijk_sq + log_sum_voters + hete + chi2e +
#    chi2e_sq + dwr_dwc_sum + dwr_dwc_log_gap + var_vac_sum + var_vac_log_gap

formula_table5 <- EI ~ log_sum_voters + hete + chi2e +
    chi2e_sq + dwr_dwc_sum + var_vac_sum

# formula_table5 <- EI ~ log_sum_voters + hete + chi2e +
# chi2e_sq + dwr_dwc_sum + var_vac_sum

# INSTANCE DUMMIES: set to FALSE, or comment the next two lines, to remove
# election/instance fixed effects such as SCOTLAND, NZ2002, etc.
use_instance_dummies <- FALSE
if (use_instance_dummies) formula_table5 <- update(formula_table5, . ~ . + factor(base))

model_frames <- setNames(
    lapply(response_order, function(response_name) build_model_frame(wide_data, response_name)),
    response_order
)

models <- lapply(model_frames, function(df_model) lm(formula_table5, data = df_model))
model_list <- unname(models)

dir.create(dirname(out_tex), recursive = TRUE, showWarnings = FALSE)

table_note <- paste(
    "* p<0.1; ** p<0.05; *** p<0.01.",
    "Linear predictors are standardized before fitting; quadratic terms are squares of those standardized predictors and are not re-standardized.",
    if (use_common_sample) {
        "All columns use the common sample across methods."
    } else {
        "Each column uses its available sample."
    }
)

render_stargazer <- function(type, ...) {
    m1 <- model_list[[1]]
    m2 <- if (length(model_list) >= 2) model_list[[2]] else NULL
    m3 <- if (length(model_list) >= 3) model_list[[3]] else NULL
    m4 <- if (length(model_list) >= 4) model_list[[4]] else NULL
    m5 <- if (length(model_list) >= 5) model_list[[5]] else NULL
    m6 <- if (length(model_list) >= 6) model_list[[6]] else NULL

    if (length(model_list) == 5) {
        return(stargazer(m1, m2, m3, m4, m5, type = type, ...))
    }

    if (length(model_list) == 6) {
        return(stargazer(m1, m2, m3, m4, m5, m6, type = type, ...))
    }

    stop(sprintf("Unexpected number of models for stargazer: %d", length(model_list)))
}

stargazer_label_map <- c(
    ratio_ijk = "R",
    ratio_ijk_sq = "R2",
    log_sum_voters = "logVoters",
    log_sum_voters_sq = "logVoters**2",
    hete = "HETe",
    chi2e = "chi2e",
    chi2e_sq = "chi2e2",
    dwr_dwc_sum = "DWR + DWC",
    dwr_dwc_log_gap = "max log-ratio(DWR, DWC)",
    var_vac_sum = "VAR + VAC",
    var_vac_log_gap = "max log-ratio(VAR, VAC)"
)

latex_label_map <- c(
    ratio_ijk = "$B/(CG)$",
    ratio_ijk_sq = "$(B/(CG))^2$",
    log_sum_voters = "$\\log (I/BCG)$",
    log_sum_voters_sq = "$(\\log (I/BCG))^2$",
    hete = "HETe",
    chi2e = "chi2e",
    chi2e_sq = "(chi2e)$^2$",
    dwr_dwc_sum = "TD",
    dwr_dwc_log_gap = "max log-ratio(DWR, DWC)",
    var_vac_sum = "TVA",
    var_vac_log_gap = "max log-ratio(VAR, VAC)"
)

active_formula_terms <- function(formula_obj) {
    terms_raw <- attr(stats::terms(formula_obj), "term.labels")
    terms_raw[terms_raw != "factor(base)"]
}

active_stargazer_labels <- function(formula_obj) {
    terms_raw <- active_formula_terms(formula_obj)
    unname(vapply(terms_raw, function(term) {
        label <- stargazer_label_map[term]
        if (length(label) == 1L && !is.na(label)) label else term
    }, character(1)))
}

active_latex_terms <- function(formula_obj) {
    terms_raw <- active_formula_terms(formula_obj)
    lapply(terms_raw, function(term) {
        label <- latex_label_map[term]
        if (!(length(label) == 1L && !is.na(label))) label <- term
        list(label = unname(label), term = term)
    })
}

sink_file <- tempfile(fileext = ".log")
sink(sink_file)
text_table <- render_stargazer(
    type = "text",
    title = "Performance across instances covariates",
    dep.var.caption = "Dependent variable: EI",
    dep.var.labels.include = FALSE,
    covariate.labels = active_stargazer_labels(formula_table5),
    keep.stat = c("rsq", "ser"),
    digits = 3,
    notes = table_note,
    notes.align = "l",
    notes.append = FALSE,
    omit.table.layout = "n"
)
sink()
unlink(sink_file)
cat(paste(text_table, collapse = "\n"), "\n")

format_latex_number <- function(x, digits = 3) {
    if (!is.finite(x)) {
        return("")
    }
    out <- sprintf(paste0("%.", digits, "f"), abs(x))
    if (x < 0) {
        paste0("$-$", out)
    } else {
        out
    }
}

format_stars <- function(p_value) {
    if (!is.finite(p_value)) {
        return("")
    }
    if (p_value < 0.01) {
        return("$^{***}$")
    }
    if (p_value < 0.05) {
        return("$^{**}$")
    }
    if (p_value < 0.1) {
        return("$^{*}$")
    }
    ""
}

format_coef_cell <- function(model, term) {
    coef_mat <- summary(model)$coefficients
    if (!term %in% rownames(coef_mat)) {
        return("")
    }
    paste0(format_latex_number(coef_mat[term, "Estimate"]), format_stars(coef_mat[term, "Pr(>|t|)"]))
}

format_se_cell <- function(model, term) {
    coef_mat <- summary(model)$coefficients
    if (!term %in% rownames(coef_mat)) {
        return("")
    }
    paste0("(", sprintf("%.3f", coef_mat[term, "Std. Error"]), ")")
}

format_r2_cell <- function(model) {
    sprintf("%.3f", summary(model)$r.squared)
}

latex_row <- function(label, term) {
    coef_cells <- vapply(model_list, format_coef_cell, character(1), term = term)
    se_cells <- vapply(model_list, format_se_cell, character(1), term = term)
    c(
        paste0(" ", label, " & ", paste(coef_cells, collapse = " & "), " \\\\"),
        paste0("  & ", paste(se_cells, collapse = " & "), " \\\\"),
        "  & & & & & & \\\\"
    )
}

if (length(model_list) != 6) {
    stop("The exact LaTeX format requires exactly 6 models.")
}

latex_terms <- c(
    active_latex_terms(formula_table5),
    list(list(label = "Constant", term = "(Intercept)"))
)

latex_body <- unlist(lapply(latex_terms, function(x) latex_row(x$label, x$term)), use.names = FALSE)
latex_body <- latex_body[-length(latex_body)]

latex_table <- c(
    "\\begin{table}[!htbp] \\centering",
    "  \\caption{Performance across instances covariates}",
    "  \\label{tab:regresion}",
    "\\begin{tabular}{@{\\extracolsep{5pt}}lcccccc}",
    "\\\\[-1.8ex]\\hline",
    "\\hline \\\\[-1.8ex]",
    "%& \\multicolumn{6}{c}{Dependent variable: $\\text{EI}_V$} \\\\",
    "%\\cmidrule(lr){2-7}",
    " & \\multicolumn{3}{c}{Error Index ($\\text{EI}_V$)} & \\multicolumn{3}{c}{Error Index Differences} \\\\",
    "\\cmidrule(lr){2-4} \\cmidrule(lr){5-7}",
    " & (1) & (2) & (3) & (4) & (5) & (6) \\\\",
    "\\midrule",
    latex_body,
    "\\midrule",
    paste0("R$^{2}$ & ", paste(vapply(model_list, format_r2_cell, character(1)), collapse = " & "), " \\\\"),
    "\\hline",
    "\\hline \\\\[-1.8ex]",
    "\\end{tabular}",
    "",
    "\\vspace{0.3em}",
    "",
    "\\begin{minipage}{\\textwidth}",
    "\\footnotesize \\textit{Note:} * p$<$0.1; ** p$<$0.05; *** p$<$0.01. Each predictor is standardized before fitting the regression. The used models are \\texttt{mvn\\_pdf\\_sym} (1), \\texttt{mvn\\_cdf\\_sym} (2), \\texttt{mult\\_sym} (3), \\texttt{mult\\_sym} - \\texttt{mvn\\_pdf\\_sym} (4), \\texttt{mvn\\_cdf\\_sym} - \\texttt{mvn\\_pdf\\_sym} (5) and \\texttt{mult\\_sym} - \\texttt{mvn\\_cdf\\_sym} (6)",
    "\\end{minipage}",
    "\\end{table}"
)
writeLines(latex_table, out_tex)
cat("\nCorrelation matrix of predictors:\n")
cat(paste(render_correlation_matrix(predictor_base), collapse = "\n"), "\n")

message(sprintf("Analysis rows: %d", nrow(analysis_rows)))
message(sprintf(
    "Stargazer column order: %s",
    paste(sprintf("(%d) %s", seq_along(response_labels), response_labels), collapse = " | ")
))
message(sprintf("LaTeX table written to: %s", normalizePath(out_tex, mustWork = FALSE)))
