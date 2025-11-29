#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(dplyr)
    library(stringr)
    library(tidyr)
    library(tibble)
    library(jsonlite)
    library(foreach)
    library(doParallel)
    library(stargazer)
})



# ============ CLI ============
args <- commandArgs(trailingOnly = TRUE)
stop_if <- function(cond, msg) {
    if (isTRUE(cond)) stop(msg, call. = FALSE)
}
arg_chr <- function(pos, default = NULL) {
    if (length(args) >= pos) args[pos] else default
}
get_flag_val <- function(flag, default = NULL) {
    m <- grep(paste0("^--", flag, "="), args, value = TRUE)
    if (length(m)) sub(paste0("^--", flag, "="), "", m[1]) else default
}

# Posicionales:
# 1) target_method (método objetivo a comparar)
# 2) field_to_test (campo del JSON, por defecto "EI_V")
target_method <- arg_chr(1, NULL)
field_to_test <- arg_chr(2, "EI_V")

# Flags
inst_like <- get_flag_val("inst-like", "ei_")
limit_inst <- suppressWarnings(as.integer(get_flag_val("limit-inst", NA_integer_)))
workers_arg <- suppressWarnings(as.integer(get_flag_val("workers", NA_integer_)))

stop_if(
    is.null(target_method) || is.na(target_method),
    "Usage: Rscript src/pval_matrix_by_instance.R <target_method> [field] [--inst-like=ei_] [--limit-inst=N] [--workers=M]"
)

# ============ Descubrimiento ============
base_ei <- file.path("output", "ei_instances")
stop_if(!dir.exists(base_ei), sprintf("'%s' doesn't exist. Run from project root.", base_ei))

all_instances <- sort(list.dirs(base_ei, recursive = FALSE, full.names = FALSE))
instances <- all_instances[str_detect(all_instances, paste0("^", inst_like))]
stop_if(length(instances) == 0, sprintf("No instances starting with '%s' found in '%s'.", inst_like, base_ei))
if (is.finite(limit_inst) && limit_inst > 0) instances <- head(instances, limit_inst)

# Baselines
baselines <- c("mult_project_lp_FALSE_sym", "mvn_cdf_project_lp_FALSE_sym", "mvn_pdf_project_lp_FALSE_sym")

discover_methods <- function(inst) {
    di <- file.path(base_ei, inst)
    dists <- sort(list.dirs(di, recursive = FALSE, full.names = FALSE))
    if (!length(dists)) {
        return(character(0))
    }
    sample_dists <- head(dists, 10)
    unique(unlist(lapply(sample_dists, function(d) {
        p <- file.path(di, d)
        basename(list.dirs(p, recursive = FALSE, full.names = FALSE))
    })))
}
methods_found <- sort(unique(unlist(lapply(head(instances, 5), discover_methods))))
stop_if(!length(methods_found), "No method folders detected under the given instances.")

row_methods <- setdiff(intersect(baselines, methods_found), target_method)
stop_if(!length(row_methods), "None of the specified baselines are present in the detected instances.")

# ============ Lectura JSON con caché ============
read_field_from_json <- function(f, field) {
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

# Promedio por distrito (inst, method, district)
by_district_for <- function(inst, method, field) {
    base_inst <- file.path(base_ei, inst)
    districts <- sort(list.dirs(base_inst, recursive = FALSE, full.names = FALSE))
    if (!length(districts)) {
        return(tibble(inst = character(0), method = character(0), district = character(0), mean_val = double(0)))
    }
    rows <- lapply(districts, function(d) {
        files <- Sys.glob(file.path(base_inst, d, method, "*.json"))
        files <- files[!grepl("/detail_", files)]
        if (!length(files)) {
            return(tibble(inst = inst, method = method, district = d, mean_val = NA_real_))
        }
        vals <- vapply(files, get_cached_field, numeric(1), field = field)
        vals <- vals[is.finite(vals)]
        tibble(
            inst = inst, method = method, district = d,
            mean_val = if (!length(vals)) NA_real_ else mean(vals)
        )
    })
    bind_rows(rows)
}

# ============ Dif y p-valor por instancia ============
# Para cada instancia y método fila:
#   diffs por distrito = (target - metodo)
#   mean_diff = mean(diffs)
#   pval = t.test(diffs == 0)
diff_and_stats_for_instance <- function(inst, target_method, row_methods, field) {
    xt <- by_district_for(inst, target_method, field) %>%
        select(district, mean_val) %>%
        rename(val_t = mean_val)

    rows <- lapply(row_methods, function(meth) {
        xm <- by_district_for(inst, meth, field) %>%
            select(district, mean_val) %>%
            rename(val_m = mean_val)
        j <- inner_join(xt, xm, by = "district") %>%
            filter(is.finite(val_t), is.finite(val_m))

        if (nrow(j) < 2) {
            list(
                per_inst_row = {
                    t <- tibble(
                        Method = meth, inst = inst, mean_diff = NA_real_, pval = NA_real_,
                        n_pairs = 0L, sum_diff = 0.0, sumsq_diff = 0.0
                    )
                    t
                }
            )
        } else {
            dif <- j$val_t - j$val_m # <-- target - metodo (como pediste)
            mean_d <- mean(dif)
            # p-valor (t-test a 1 muestra sobre difs)
            pv <- tryCatch(stats::t.test(dif, mu = 0), error = function(e) NULL)
            pv <- if (is.null(pv)) NA_real_ else pv$p.value
            n <- length(dif)
            s1 <- sum(dif)
            s2 <- sum(dif^2)

            list(
                per_inst_row = tibble(
                    Method = meth, inst = inst, mean_diff = mean_d, pval = pv,
                    n_pairs = as.integer(n), sum_diff = s1, sumsq_diff = s2
                )
            )
        }
    })

    bind_rows(lapply(rows, `[[`, "per_inst_row"))
}

# ============ Paralelización por instancia ============
n_cores <- if (is.finite(workers_arg) && workers_arg > 0) workers_arg else max(1L, parallel::detectCores(logical = TRUE) - 1L)
cl <- parallel::makeCluster(n_cores)
on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
doParallel::registerDoParallel(cl)

parallel::clusterExport(
    cl,
    varlist = c(
        "base_ei", "instances", "row_methods", "target_method",
        "field_to_test", "cache_env",
        "read_field_from_json", "get_cached_field",
        "by_district_for", "diff_and_stats_for_instance"
    ),
    envir = environment()
)
parallel::clusterEvalQ(cl, {
    library(dplyr)
    library(tibble)
    library(stringr)
    NULL
})

parts <- foreach(inst = instances, .packages = c("dplyr", "tibble")) %dopar% {
    diff_and_stats_for_instance(inst, target_method, row_methods, field_to_test)
}

# Long con todas las filas/instancias
stats_long <- bind_rows(parts) %>%
    arrange(factor(Method, levels = row_methods), inst)

# ============ Agregar TOTAL (global) ============
global_stats <- stats_long %>%
    group_by(Method) %>%
    summarise(
        n_total = sum(n_pairs, na.rm = TRUE),
        sum1 = sum(sum_diff, na.rm = TRUE),
        sum2 = sum(sumsq_diff, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    mutate(
        mean_diff = ifelse(n_total > 0, sum1 / n_total, NA_real_),
        var_d = ifelse(n_total > 1, (sum2 - n_total * mean_diff^2) / (n_total - 1), NA_real_),
        se = sqrt(var_d / pmax(n_total, 1L)),
        tstat = mean_diff / se,
        pval = ifelse(is.finite(tstat), 2 * stats::pt(-abs(tstat), df = pmax(n_total - 1, 1L)), NA_real_)
    ) %>%
    select(Method, mean_diff, pval) %>%
    arrange(factor(Method, levels = row_methods))

# ============ Wide: una columna por instancia (dos filas por método) ============
# Formateadores
fmt_dif <- function(x) ifelse(is.na(x), "", sprintf("%.4f", as.numeric(x)))
fmt_p <- function(p) ifelse(is.na(p), "", sprintf("(%.3g)", as.numeric(p)))

# 1) Tablas separadas de diferencias y p-valores por instancia (wide)
by_inst_diff <- stats_long %>%
    select(Method, inst, mean_diff) %>%
    tidyr::pivot_wider(id_cols = Method, names_from = inst, values_from = mean_diff)

by_inst_p <- stats_long %>%
    select(Method, inst, pval) %>%
    tidyr::pivot_wider(id_cols = Method, names_from = inst, values_from = pval)

# Asegurar mismo orden de métodos
by_inst_diff <- by_inst_diff %>% arrange(factor(Method, levels = row_methods))
by_inst_p <- by_inst_p %>% arrange(factor(Method, levels = row_methods))

# 2) Columna Total (global): diferencias y p-valores por separado
total_diff <- tibble(Total = global_stats$mean_diff) %>% mutate(Total = fmt_dif(Total))
total_p <- tibble(Total = global_stats$pval) %>% mutate(Total = fmt_p(Total))

# 3) Formateo (string) de cada celda
inst_names <- setdiff(colnames(by_inst_diff), "Method")

diff_wide_str <- by_inst_diff
p_wide_str <- by_inst_p

if (length(inst_names)) {
    diff_wide_str[, inst_names] <- lapply(by_inst_diff[, inst_names], fmt_dif)
    p_wide_str[, inst_names] <- lapply(by_inst_p[, inst_names], fmt_p)
}

# Añadir Total formateado
diff_wide_str <- diff_wide_str %>%
    bind_cols(total_diff)

p_wide_str <- p_wide_str %>%
    bind_cols(total_p)

# 4) Intercalar filas: por cada método, primero la fila con ∆, después la fila con (p)
interleave_rows <- function(diff_df, p_df) {
    stopifnot(identical(diff_df$Method, p_df$Method))
    out_list <- vector("list", length = nrow(diff_df) * 2)
    k <- 1L
    for (i in seq_len(nrow(diff_df))) {
        # Fila 1: diferencias
        out_list[[k]] <- diff_df[i, , drop = FALSE]
        # Fila 2: p-valores (entre paréntesis); vaciamos la etiqueta del método
        row_p <- p_df[i, , drop = FALSE]
        row_p$Method[1] <- ""
        out_list[[k + 1]] <- row_p
        k <- k + 2L
    }
    dplyr::bind_rows(out_list)
}

table_wide_2rows <- interleave_rows(diff_wide_str, p_wide_str)

# ============ Salidas ============
# ASCII bonitillo
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

cat("\n====== MEAN DIFFERENCE (target - row method) with p-values on the next row ======\n")
cat(ascii_table(table_wide_2rows), sep = "\n")
cat("\nTarget method:", target_method, " | Field:", field_to_test, " | Workers:", n_cores, "\n")
cat("Rows: for each method, the first row shows mean difference across districts; the second row shows the paired t-test p-value in parentheses.\n")
cat("The 'Total' column pools ALL districts across ALL instances.\n")
cat("=================================================================================\n")

# ============ LaTeX tipo STARGAZER: dos filas por método ============
latex_df <- as.data.frame(table_wide_2rows, stringsAsFactors = FALSE)

cat("\n% ====== STARGAZER-STYLE TABLE (coef row + p-value row) ======\n")
stargazer(
    latex_df,
    summary = FALSE,
    rownames = FALSE,
    type = "latex",
    title = paste0("Mean Difference (", target_method, " − row method) by instance"),
    label = "tab:mean_diff_pvals_twolines",
    header = FALSE,
    digits = 4,
    font.size = "small",
    notes = c(
        "Top row: mean difference (target − method) across districts.",
        "Bottom row: paired t-test p-value in parentheses.",
        "Total pools all districts across all instances."
    ),
    notes.align = "l",
    notes.append = FALSE
)
cat("% ======================================================================\n")
