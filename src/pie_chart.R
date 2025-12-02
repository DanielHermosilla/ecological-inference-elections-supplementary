#!/usr/bin/env Rscript
# Builds pairwise win/tie pie charts comparing EI methods across instances.

suppressPackageStartupMessages({
    library(dplyr)
    library(stringr)
    library(tidyr)
    library(tibble)
    library(jsonlite)
    library(ggplot2)
    library(foreach)
    library(doParallel)
})

# ================= CLI =================
args <- commandArgs(trailingOnly = TRUE)
arg_chr <- function(pos, default = NULL) if (length(args) >= pos) args[pos] else default
get_flag <- function(flag, default = NULL) {
    m <- grep(paste0("^--", flag, "="), args, value = TRUE)
    if (length(m)) sub(paste0("^--", flag, "="), "", m[1]) else default
}
get_flag_lgl <- function(flag, default = FALSE) {
    v <- tolower(get_flag(flag, NA_character_))
    if (is.na(v)) {
        return(default)
    }
    v %in% c("1", "true", "t", "yes", "y")
}
stop_if <- function(cond, msg) {
    if (isTRUE(cond)) stop(msg, call. = FALSE)
}

# Positional: output file
outfile <- arg_chr(1, "figures/pairwise_pie_grid.pdf")

# Flags
field_to_use <- get_flag("field", "EI_V")
inst_like <- get_flag("inst-like", "ei_")
limit_inst <- suppressWarnings(as.integer(get_flag("limit-inst", NA_integer_)))
workers_arg <- suppressWarnings(as.integer(get_flag("workers", NA_integer_)))
higher_better <- get_flag_lgl("higher-better", FALSE) # default: lower is better
use_parallel <- get_flag_lgl("parallel", FALSE)

# ================= Paths =================
base_ei <- file.path("output", "ei_instances")
stop_if(!dir.exists(base_ei), sprintf("'%s' doesn't exist. Run from project root.", base_ei))

# ================= Instances =================
all_instances <- sort(list.dirs(base_ei, recursive = FALSE, full.names = FALSE))
instances <- all_instances[str_detect(all_instances, paste0("^", inst_like))]
stop_if(length(instances) == 0, sprintf("No instances starting with '%s' found in '%s'.", inst_like, base_ei))
if (is.finite(limit_inst) && limit_inst > 0) instances <- head(instances, limit_inst)

# ================= Methods (ORDERED) =================
# Row facets in EXACT order you want
row_methods <- c(
    "mvn_cdf_project_lp_FALSE_sym",
    "mvn_pdf_project_lp_FALSE_sym",
    "mult_project_lp_FALSE_sym"
)

# Master order for ALL methods; column facets will follow this order (minus row_methods)
methods <- c(
    "lphom",
    "lclphom",
    "mult_project_lp_FALSE_sym",
    "mvn_cdf_project_lp_FALSE_sym",
    "mvn_pdf_project_lp_FALSE_sym",
    "nslphom_dual_a",
    "nslphom_dual_w",
    "nslphom_joint",
    "ei.md.bayes",
    "ecolRxC"
)

# Exclude any row method from the columns (PRESERVE the order in `methods`)
col_methods <- methods[!methods %in% row_methods]

all_needed_methods <- sort(unique(c(row_methods, col_methods)))

# ================= JSON reading + cache =================
# Read one numeric field from a JSON result file, safely returning NA on error.
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
# Cache JSON reads to avoid repeated disk I/O across workers.
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

# District-level means per (inst, method, district)
# Average a JSON field across all files for one (instance, method, district).
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

# ================= Parallel collect (inst × method) =================
grid <- expand.grid(inst = instances, method = all_needed_methods, stringsAsFactors = FALSE)

n_cores <- 1L
district_means_list <- NULL

if (use_parallel) {
    n_cores <- if (is.finite(workers_arg) && workers_arg > 0) {
        workers_arg
    } else {
        max(1L, parallel::detectCores(logical = TRUE) - 1L)
    }
    cl <- parallel::makeCluster(n_cores)
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    doParallel::registerDoParallel(cl)
    parallel::clusterExport(cl,
        varlist = c(
            "base_ei", "instances", "all_needed_methods", "field_to_use", "cache_env",
            "read_field_from_json", "get_cached_field", "by_district_for"
        ),
        envir = environment()
    )
    parallel::clusterEvalQ(cl, {
        library(dplyr)
        library(tibble)
        NULL
    })

    district_means_list <- foreach(i = seq_len(nrow(grid)), .packages = c("dplyr", "tibble")) %dopar% {
        by_district_for(grid$inst[i], grid$method[i], field_to_use)
    }
} else {
    message("Running sequentially; use --parallel=true to enable parallel workers.")
    district_means_list <- lapply(seq_len(nrow(grid)), function(i) {
        by_district_for(grid$inst[i], grid$method[i], field_to_use)
    })
}

district_means <- bind_rows(district_means_list)

# ================= Wide panel (inst × district) =================
wide_panel <- district_means %>%
    filter(method %in% all_needed_methods, is.finite(mean_val)) %>%
    select(inst, district, method, mean_val) %>%
    tidyr::pivot_wider(id_cols = c(inst, district), names_from = method, values_from = mean_val)

# ================= Pairwise comparisons =================
wins_row_over_col <- function(x_row, x_col, higher_is_better = FALSE, eps = 0) {
    if (is.na(x_row) || is.na(x_col)) {
        return(NA_character_)
    }
    if (!higher_is_better) {
        if (x_row + eps < x_col - eps) {
            return("row")
        }
        if (x_col + eps < x_row - eps) {
            return("col")
        }
        return("tie")
    } else {
        if (x_row - eps > x_col + eps) {
            return("row")
        }
        if (x_col - eps > x_row + eps) {
            return("col")
        }
        return("tie")
    }
}

pair_results <- lapply(row_methods, function(rm) {
    lapply(col_methods, function(cm) {
        cols_needed <- c(rm, cm)
        sub <- wide_panel %>% select(any_of(cols_needed))
        ok <- stats::complete.cases(sub)
        sub <- sub[ok, , drop = FALSE]
        if (!nrow(sub)) {
            tibble(
                row_method = rm, col_method = cm,
                row_wins = NA_real_, col_wins = NA_real_, ties = NA_real_, n_comp = 0L
            )
        } else {
            res <- mapply(
                wins_row_over_col,
                x_row = sub[[rm]], x_col = sub[[cm]],
                MoreArgs = list(higher_is_better = higher_better, eps = 0),
                SIMPLIFY = TRUE, USE.NAMES = FALSE
            )
            row_w <- mean(res == "row")
            col_w <- mean(res == "col")
            tie_w <- mean(res == "tie")
            tibble(
                row_method = rm, col_method = cm,
                row_wins = row_w, col_wins = col_w, ties = tie_w, n_comp = as.integer(length(res))
            )
        }
    }) %>% bind_rows()
}) %>%
    bind_rows() %>%
    mutate(diag_flag = (row_method == col_method)) # robust flag (won't happen now but safe)

# ================= Pie data =================
pie_long <- pair_results %>%
    tidyr::pivot_longer(cols = c(row_wins, col_wins, ties), names_to = "who", values_to = "prop") %>%
    mutate(
        who = dplyr::recode(who,
            col_wins = "Benchmark method",
            row_wins = "Proposed method",
            ties = "Ties"
        ),
        prop = ifelse(is.na(prop), 0, prop),
        row_lab_id = row_method,
        col_lab_id = col_method
    )

# ================= Labels (safe recode) =================
row_labels_map <- c(
    "mvn_cdf_project_lp_FALSE_sym" = "mvn_cdf_sym",
    "mvn_pdf_project_lp_FALSE_sym" = "mvn_pdf_sym",
    "mult_project_lp_FALSE_sym"    = "mult_sym"
)
col_labels_map <- c(
    "ei.md.bayes" = "ei.MD.bayes",
    "lclphom" = "lclphom",
    "lphom" = "lphom",
    "mult_project_lp_FALSE_sym" = "mult_sym",
    "mvn_cdf_project_lp_FALSE_sym" = "mvn_cdf_sym",
    "mvn_pdf_project_lp_FALSE_sym" = "mvn_pdf_sym",
    "nslphom_joint" = "nslphom_joint",
    "nslphom_dual_a" = "nslphom_dual_a",
    "nslphom_dual_w" = "nslphom_dual_w"
)

pie_long <- pie_long %>%
    mutate(
        row_lab = dplyr::recode(row_lab_id, !!!row_labels_map, .default = row_lab_id, .missing = row_lab_id),
        col_lab = dplyr::recode(col_lab_id, !!!col_labels_map, .default = col_lab_id, .missing = col_lab_id)
    )

# (Diagonal safeguard: if there ever were overlaps, force 100% ties for diag)
pie_long <- pie_long %>%
    mutate(prop = ifelse(diag_flag & who == "Ties", 1,
        ifelse(diag_flag, 0, prop)
    ))

# ================= Center label: % row wins =================
label_df <- pair_results %>%
    mutate(
        row_lab = dplyr::recode(row_method, !!!row_labels_map, .default = row_method, .missing = row_method),
        col_lab = dplyr::recode(col_method, !!!col_labels_map, .default = col_method, .missing = col_method),
        pct_row = ifelse(is.na(row_wins), NA_real_, 100 * row_wins)
    )

# ================= Facet order from mapped labels =================
map_levels <- function(keys, map) {
    mapped <- unname(map[keys])
    mapped[is.na(mapped)] <- keys[is.na(mapped)]
    mapped
}
row_levels <- map_levels(row_methods, row_labels_map)
col_levels <- map_levels(col_methods, col_labels_map)

pie_long$row_lab <- factor(pie_long$row_lab, levels = row_levels)
pie_long$col_lab <- factor(pie_long$col_lab, levels = col_levels)
label_df$row_lab <- factor(label_df$row_lab, levels = row_levels)
label_df$col_lab <- factor(label_df$col_lab, levels = col_levels)

# ================= Style =================
if (requireNamespace("showtext", quietly = TRUE) &&
    requireNamespace("sysfonts", quietly = TRUE)) {
    if (!"Fira Sans" %in% sysfonts::font_families()) {
        sysfonts::font_add_google("Fira Sans", "Fira Sans")
    }
    showtext::showtext_auto()
    base_family <- "Fira Sans"
} else {
    base_family <- "Fira Sans"
}

oi <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#F0E442", "#56B4E9", "#999999")
fill_map <- c(
    "Proposed method" = "#B22222", # red B22222 for the paper, blue 0000FF for presentation
    "Benchmark method" = "#1f77b400", # blue with alpha 0 (transparent)A
    "Ties" = "grey80"
)

# Flip the filling direction
who_levels <- c("Proposed method", "Benchmark method", "Ties")
pie_long <- pie_long %>%
    mutate(
        who = recode(who,
            "Proposed method" = "Proposed method",
            "Competing method" = "Benchmark method", # rename
            "Ties" = "Ties"
        ),
        who = factor(who, levels = who_levels)
    )

# ================= Plot =================
p <- ggplot(pie_long, aes(x = "1", y = prop, fill = who)) +
    geom_col(width = 1, color = "white") +
    coord_polar(theta = "y") +
    facet_grid(rows = vars(row_lab), cols = vars(col_lab), switch = "y") +
    scale_fill_manual(
        values = fill_map,
        breaks = c("Proposed method", "Benchmark method") # omit "Ties" from legend
    ) +
    theme_bw(base_family = base_family, base_size = 8) +
    labs(fill = NULL, x = NULL, y = NULL) +
    theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "none",
        panel.spacing.x = unit(0.2, "lines"),
        panel.spacing.y = unit(0.2, "lines")
    )

# Center % label
p <- p + geom_text(
    data = label_df,
    aes(
        x = 0, y = 0,
        label = dplyr::case_when(
            is.na(pct_row) ~ "",
            TRUE ~ sprintf("%.0f%%", pct_row)
        )
    ),
    inherit.aes = FALSE,
    size = 2,
    fontface = "bold"
)

# ================= Save =================
dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)
ext <- tolower(tools::file_ext(outfile))
if (ext %in% c("png", "jpg", "jpeg", "tiff", "bmp")) {
    # ggsave(outfile, plot = p, width = 18, height = 10, dpi = 300)
    ggsave(outfile, plot = p, width = 6, height = 3.3)
} else if (ext %in% c("pdf")) {
    ggsave(outfile, plot = p, width = 6, height = 3.3, device = grDevices::cairo_pdf)
} else {
    out_pdf <- sub("\\.[A-Za-z0-9]+$", ".pdf", outfile)
    ggsave(out_pdf, plot = p, width = 18, height = 10, device = grDevices::cairo_pdf)
    message(sprintf("Unrecognized extension; saved as PDF: %s", out_pdf))
}

message(sprintf("[OK] Figure saved to: %s", normalizePath(outfile, mustWork = FALSE)))
