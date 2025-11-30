# Build MAE_p2 error boxplots for simulated instances.
suppressPackageStartupMessages({
    library(jsonlite)
    library(dplyr)
    library(purrr)
    library(stringr)
    library(forcats)
    library(ggplot2)
})
source("src/R_functions.R")

# ------------------ PLOT ------------------
# Boxplot of MAE_p2 across methods, faceting by groups and candidates.
plot_EM_error_boxplot_R <- function(
    df,
    save_dir = "figures",
    save_name = "fig2-p-error.pdf",
    # --- Additional parameters ---
    groups_keep = NULL, # groups G to keep, e.g. c("1","2","3")
    cand_keep = NULL, # candidates C to keep (overrides 'candidates')
    candidates = c("all", "uar"), # "all" = all candidates; "uar" = uniform-at-random sampling
    n_cand = NULL, # number of candidates to sample if candidates = "uar"
    seed = 123, # seed for UAR
    letter_size = 6) {
    candidates <- match.arg(candidates)

    # ---- Optional Fira Sans font ----
    base_family <- "sans"
    if (requireNamespace("showtext", quietly = TRUE) &&
        requireNamespace("sysfonts", quietly = TRUE)) {
        have_font <- "Fira Sans" %in% sysfonts::font_families()
        if (!have_font) {
            have_font <- tryCatch({
                sysfonts::font_add_google("Fira Sans", "Fira Sans")
                TRUE
            }, error = function(e) FALSE)
        }
        if (have_font) {
            showtext::showtext_auto()
            base_family <- "Fira Sans"
        }
    }

    # ---- Local copy ----
    d <- df

    # 1) Methods to display (and labels)
    full_order <- c(
        "ei.MD.bayes", "nslphom_dual_w",
        "mvn_cdf_project_lp_FALSE", "mvn_pdf_project_lp_FALSE",
        "mult_project_lp_FALSE", "exact_project_lp_FALSE"
    )
    full_labels <- c(
        "ei.md.bayes", "nslphom_dual_w",
        "mvn_cdf", "mvn_pdf", "mult", "exact"
    )

    # 2) Methods actually present in d
    present <- intersect(full_order, unique(d$method))

    # 3) Keep only present methods and drop problematic NAs
    d <- d %>%
        dplyr::filter(method %in% present) %>%
        dplyr::filter(!is.na(G), !is.na(C)) %>%
        dplyr::filter(!is.na(MAE_p2))

    # ---- Filters by group and candidate ----
    # Normalize to character for robust filtering
    d <- d %>%
        mutate(
            G = as.character(G),
            C = as.character(C)
        )

    # 3a) Filter by groups if provided
    if (!is.null(groups_keep)) {
        groups_keep <- as.character(groups_keep)
        d <- d %>% filter(G %in% groups_keep)
        if (nrow(d) == 0) stop("Filtering by 'groups_keep' removed all rows.")
    }

    # 3b) Filter by candidates:
    #     - If 'cand_keep' is passed, use it as-is.
    #     - Otherwise, if candidates == "uar", sample 'n_cand' candidates.
    #     - If candidates == "all", no candidate filtering is applied.
    cand_levels_all <- sort(unique(d$C))

    if (!is.null(cand_keep)) {
        cand_keep <- as.character(cand_keep)
        missing_c <- setdiff(cand_keep, cand_levels_all)
        if (length(missing_c)) {
            warning(
                "Some 'cand_keep' candidates are missing in the data and will be ignored: ",
                paste(missing_c, collapse = ", ")
            )
        }
        cand_keep_in <- intersect(cand_keep, cand_levels_all)
        if (!length(cand_keep_in)) stop("None of the 'cand_keep' candidates exist in the data.")
        d <- d %>% filter(C %in% cand_keep_in)
    } else if (candidates == "uar") {
        if (is.null(n_cand) || !is.numeric(n_cand) || n_cand < 1) {
            stop("For candidates = 'uar' you must provide an 'n_cand' >= 1.")
        }
        if (n_cand > length(cand_levels_all)) {
            warning("n_cand > available candidates; using n_cand = ", length(cand_levels_all))
            n_cand <- length(cand_levels_all)
        }
        set.seed(seed)
        sampled_c <- sample(cand_levels_all, n_cand, replace = FALSE)
        d <- d %>% filter(C %in% sampled_c)
    } # if candidates == "all": no filtering on C

    # Re-factor for ordered, clean facets
    # --- Ordered C levels ---
    if (!is.null(cand_keep)) {
        c_levels <- cand_keep_in
    } else {
        c_levels <- sort(unique(as.numeric(d$C)))
        c_levels <- as.character(c_levels)
    }

    d <- d %>%
        mutate(
            G = factor(G, levels = sort(unique(G))),
            C = factor(C, levels = c_levels)
        )

    # 4) Recode method labels
    method_lut <- stats::setNames(full_labels[match(present, full_order)], present)
    d <- d %>%
        mutate(method = factor(method, levels = present, labels = unname(method_lut)))

    # ---- Palette (Okabe–Ito) ----
    oi <- c(
        "#0072B2", "#E69F00", "#009E73", "#D55E00",
        "#CC79A7", "#F0E442", "#56B4E9", "#999999"
    )
    lvls <- levels(d$method)
    pal <- stats::setNames(oi[seq_along(lvls)], lvls)

    # ---- Plot ----
    p <- ggplot(d, aes(x = method, y = MAE_p2, fill = method)) +
        geom_boxplot(outlier.shape = 21, outlier.size = 0.4, linewidth = 0.2) +
        facet_grid(rows = vars(G), cols = vars(C), labeller = label_both) +
        scale_fill_manual(values = pal) +
        labs(x = "Method", y = "MAE_p2") +
        theme_bw(base_size = letter_size, base_family = base_family) +
        theme(
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            strip.text = element_text(size = letter_size, face = "bold"),
            legend.position = "none",
            # panel.grid.major = element_line(color = "grey90", size = 0.2),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.spacing.x = grid::unit(2, "pt"),
            panel.spacing.y = grid::unit(2, "pt"),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank()
        )

    width <- 6.3
    height <- 6.3
    # height <- 5


    dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
    ggsave(file.path(save_dir, save_name), p,
        width = width, height = height, dpi = 300,
        device = grDevices::cairo_pdf
    )
    print(p)

    # Return filtered data and plot invisibly for reuse
    invisible(list(data = d, plot = p))
}

# Coerce lists/arrays to a numeric matrix for metric computation.
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

# Compute MAE_p2 from raw JSON entry (fallback when field missing).
compute_mae_p2 <- function(entry) {
    js <- entry$json
    p_true <- as_numeric_matrix(js$p_true)
    p_est <- as_numeric_matrix(js$p_est)
    if (is.null(p_true) || is.null(p_est)) return(NA_real_)
    if (!all(dim(p_true) == dim(p_est))) return(NA_real_)
    G <- nrow(p_true)
    if (!is.finite(G) || G <= 0) return(NA_real_)
    total_abs <- sum(abs(p_est - p_true))
    total_abs / (2 * G)
}

main <- function() {
    INSTANCE_DIR <- "output/simulated_instances" # relative to repo root
    FIGURES_DIR <- "figures"

    message("Reading instances…")
    result <- read_simulated_instances(INSTANCE_DIR)
    if (!"MAE_p2" %in% names(result$df)) {
        mae_vec <- vapply(result$raw, compute_mae_p2, numeric(1))
        result$df$MAE_p2 <- mae_vec
    }

    # Example usages:
    # 1) Keep groups "1","2","3" and sample 3 candidates uniformly at random
    # plot_EM_error_boxplot_R(
    #     result$df,
    #     save_dir = FIGURES_DIR,
    #     save_name = "sim_prob_error_g123_cUAR3.pdf",
    #     groups_keep = c("1", "2", "3"),
    #     candidates = "uar",
    #     n_cand = 3,
    #     seed = 2025
    # )

    # 2) Keep groups "1","3" and specific candidates "A","B","C" (if present)
    plot_EM_error_boxplot_R(
        result$df,
        save_dir = FIGURES_DIR,
        save_name = "sim_prob_error_3.pdf",
        groups_keep = c("2", "3", "4", "6", "8"),
        cand_keep = c("2", "3", "5", "10"),
        letter_size = 6
    )

    # 3) Only filter groups and keep all candidates
    # plot_EM_error_boxplot_R(
    #     result$df,
    #     save_dir  = FIGURES_DIR,
    #     save_name = "sim_prob_error_g12_allC.pdf",
    #     groups_keep = c("1", "2"),
    #     candidates  = "all"
    # )
}

main()
