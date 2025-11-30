#!/usr/bin/env Rscript
# Visualises per-instance method wins across districts for EI datasets.

suppressPackageStartupMessages({
    library(stringr)
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(purrr)
    library(future.apply)
    library(jsonlite)
    library(ggplot2)
    library(forcats)
})

# ============== Font (optional Fira Sans via showtext) ==============
# Optional Fira Sans font loader (no-op if unavailable).
.use_fira <- function() {
    if (requireNamespace("showtext", quietly = TRUE) &&
        requireNamespace("sysfonts", quietly = TRUE)) {
        if (!"Fira Sans" %in% sysfonts::font_families()) {
            sysfonts::font_add_google("Fira Sans", "Fira Sans")
        }
        showtext::showtext_auto()
        return("Fira Sans")
    }
    "Fira Sans"
}

# ======================== CLI args helpers ===========================
args <- commandArgs(trailingOnly = TRUE)

# Simple CLI parsing helpers for flags/ints/bools.
parse_opt_chr <- function(flag, default = NULL) {
    m <- grep(paste0("^", flag, "="), args, value = TRUE)
    if (length(m)) sub(paste0("^", flag, "="), "", m[1]) else default
}
parse_opt_int <- function(flag, default = NA_integer_) {
    v <- parse_opt_chr(flag, NULL)
    if (is.null(v)) default else suppressWarnings(as.integer(v))
}
parse_opt_lgl <- function(flag, default = FALSE) {
    v <- tolower(parse_opt_chr(flag, NA_character_))
    if (is.na(v)) {
        return(default)
    }
    if (v == "") {
        return(TRUE)
    }
    v %in% c("1", "true", "t", "yes", "y")
}

# Flags (keep workers/limit/out; inst-like is not used):
#   --workers=8
#   --limit-inst=20
#   --out=figures/nz_relative_wins.pdf
req_workers <- parse_opt_int("--workers", NA_integer_)
limit_inst <- parse_opt_int("--limit-inst", NA_integer_)
outfile <- parse_opt_chr("--out", file.path("figures", "nz_relative_wins.pdf"))
use_parallel <- parse_opt_lgl("--parallel", FALSE)

# ============================ Config ================================
base_ei <- file.path("output", "ei_instances")

# Canonical method order + labels
methods <- c(
    "lclphom",
    "lphom",
    "nslphom_dual_a",
    "nslphom_dual_w",
    "nslphom_joint",
    "ei.MD.bayes",
    "ecolRxC",
    "mvn_cdf_project_lp_FALSE",
    "mvn_pdf_project_lp_FALSE",
    "mult_project_lp_FALSE"
    # "exact_project_lp_FALSE"
)
method_labels <- c(
    "lclphom",
    "lphom",
    "nslphom_dual_a",
    "nslphom_dual_w",
    "nslphom_joint",
    "ei.md.bayes",
    "ecolRxC",
    "mvn_cdf",
    "mvn_pdf",
    "mult"
    # "exact"
)
stopifnot(length(methods) == length(method_labels))

if (!dir.exists(base_ei)) {
    stop(sprintf("'%s' doesn't exist. Run this script from the project root.", base_ei))
}

# ========= Manual instance order (preserved on the x-axis) ========
instances <- c(
    "ei_NZ_2002",
    "ei_NZ_2005",
    "ei_SCO_2007",
    "ei_NZ_2008",
    "ei_NZ_2011",
    "ei_NZ_2014",
    "ei_NZ_2017",
    "ei_NZ_2020"
)

# Optional limit via flag
if (is.finite(limit_inst) && limit_inst > 0) {
    instances <- head(instances, limit_inst)
}

# Validate presence on disk; warn on missing but continue
exist_mask <- file.exists(file.path(base_ei, instances))
if (!all(exist_mask)) {
    missing <- instances[!exist_mask]
    warning(sprintf(
        "These instances do not exist in '%s' and will be skipped: %s",
        base_ei, paste(missing, collapse = ", ")
    ))
}
instances <- instances[exist_mask]
if (length(instances) == 0) {
    stop("No valid instances present in 'output/ei_instances'.")
}

# ========================= Parallel setup ===========================
if (use_parallel) {
    if (is.finite(req_workers) && req_workers > 0) {
        future::plan(future::multisession, workers = req_workers)
    } else {
        future::plan(future::multisession, workers = future::availableCores())
    }
    message("Running with parallel workers; disable via --parallel=false.")
} else {
    future::plan(future::sequential)
    message("Running sequentially; pass --parallel=true to enable workers.")
}

# ===================== JSON field cache helpers =====================
# Safe JSON reader that extracts one numeric field.
read_field_from_json <- function(f, field) {
    tryCatch(
        {
            js <- jsonlite::fromJSON(f, simplifyVector = TRUE)
            v <- js[[field]]

            if (is.null(v)) {
                return(NA_real_)
            }

            v <- suppressWarnings(as.numeric(v))

            if (!length(v)) {
                return(NA_real_)
            }

            v[[1]]
        },
        error = function(e) NA_real_
    )
}

.cache_env <- new.env(parent = emptyenv())
# Cache JSON reads by path/mtime to reduce repeated IO.
get_cached_field <- function(f, field) {
    fi <- suppressWarnings(file.info(f))
    key <- paste0(f, "||", fi$mtime, "||", field)
    got <- .cache_env[[key]]
    if (!is.null(got)) {
        return(got)
    }
    v <- read_field_from_json(f, field)
    .cache_env[[key]] <- v
    v
}

# ===================== Core: per-instance tallies ====================
# For one instance, tally per-district wins (lower is better by default).
ei_wins_for_instance <- function(inst, field = "EI_V") {
    base_inst <- file.path(base_ei, inst)
    districts <- sort(list.dirs(base_inst, recursive = FALSE, full.names = FALSE))

    if (length(districts) == 0) {
        return(tibble(
            inst = character(0), method = character(0),
            wins = double(0), districts = integer(0), prop = double(0)
        ))
    }

    rows <- lapply(districts, function(d) {
        per_method <- lapply(methods, function(m) {
            files <- Sys.glob(file.path(base_inst, d, m, "*.json"))
            files <- files[!grepl("/detail_", files)]
            if (length(files) == 0) {
                return(tibble(district = d, method = m, mean_val = NA_real_, n_files = 0L))
            }
            vals <- vapply(files, get_cached_field, numeric(1), field = field)
            vals <- vals[is.finite(vals)]
            if (!length(vals)) {
                tibble(district = d, method = m, mean_val = NA_real_, n_files = 0L)
            } else {
                tibble(district = d, method = m, mean_val = mean(vals), n_files = length(vals))
            }
        })
        bind_rows(per_method)
    })
    df <- bind_rows(rows)

    districts_valid <- df %>%
        group_by(district) %>%
        summarise(has_data = any(is.finite(mean_val)), .groups = "drop") %>%
        filter(has_data) %>%
        pull(district)

    if (!length(districts_valid)) {
        return(tibble(inst = inst, method = methods, wins = 0, districts = 0L, prop = 0))
    }

    df_valid <- df %>% filter(district %in% districts_valid)

    winners <- df_valid %>%
        group_by(district) %>%
        mutate(min_val = suppressWarnings(min(mean_val, na.rm = TRUE))) %>%
        mutate(is_winner = is.finite(mean_val) & (mean_val == min_val)) %>%
        mutate(win_share = ifelse(is_winner, 1 / sum(is_winner), 0)) %>%
        ungroup()

    out <- winners %>%
        group_by(method) %>%
        summarise(wins = sum(win_share), .groups = "drop") %>%
        mutate(inst = inst, districts = length(unique(districts_valid))) %>%
        mutate(prop = ifelse(districts > 0, wins / districts, 0)) %>%
        select(inst, method, wins, districts, prop)

    out <- right_join(out, tibble(method = methods), by = "method") %>%
        mutate(
            inst = inst,
            wins = replace_na(wins, 0),
            districts = replace_na(districts, length(unique(districts_valid))),
            prop = replace_na(prop, 0)
        ) %>%
        arrange(match(method, methods))

    out
}

# ================= Aggregate across instances (parallel) =============
message("Scanning instances (custom order):")
message(paste0(" - ", paste(instances, collapse = "\n - ")))

res_list <- future_lapply(instances, ei_wins_for_instance, field = "EI_V", future.seed = TRUE)
props_df <- bind_rows(res_list)

# === Literal instance labels to display ===
instance_labels <- c(
    "NZ_2002",
    "NZ_2005",
    "SC_2007",
    "NZ_2008",
    "NZ_2011",
    "NZ_2014",
    "NZ_2017",
    "NZ_2020"
)

# === Apply manual labels ===
props_df <- props_df %>%
    mutate(
        method = factor(method, levels = methods, labels = method_labels),
        inst   = factor(inst, levels = instances, labels = instance_labels, ordered = TRUE)
    )

# ============================== Plot =================================
base_family <- .use_fira()

oi <- c(
    "#0072B2", "#E69F00", "#999933", "#D55E00", "#CC79A7",
    "#009E73", "#56B4E9", "#999999", "#9999CC", "#663399"
)
oi <- oi[seq_len(nlevels(props_df$method))]

p <- ggplot(props_df, aes(x = inst, y = prop, fill = method)) +
    geom_bar(stat = "identity", width = 0.9, color = "white", linewidth = 0.2) +
    scale_fill_manual(values = oi, name = "Method") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
    labs(
        x = "Election",
        y = "Share of districts won"
    ) +
    theme_bw(base_size = 8, base_family = base_family) +
    theme(
        legend.position = "right",
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.15, "cm"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 7),
        axis.text.y = element_text(size = 7)
    )

# ============================== Save =================================
dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
ggsave(outfile, p, width = 6.2, height = 3, dpi = 300, device = grDevices::cairo_pdf)
message(sprintf("[OK] Saved relative stacked bar chart to: %s", normalizePath(outfile, mustWork = FALSE)))

print(p)
