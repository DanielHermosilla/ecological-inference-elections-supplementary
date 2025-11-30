# Plot runtime vs number of candidates, faceted by groups.
suppressPackageStartupMessages({
    library(jsonlite)
    library(dplyr)
    library(purrr)
    library(stringr)
    library(forcats)
    library(ggplot2)
})
source("src/R_functions.R")

# ------------------ PLOT (runtime vs C, facet por G) ------------------
plot_runtime_vsC_byG <- function(
    df,
    letter_size = 8,
    save_dir = "figures",
    save_name = "fig-runtime-vsC-byG.pdf") {
    # Optional font setup (use Fira Sans if available via showtext)
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

    # Canonical method order and display labels
    full_order <- c(
        "ei.MD.bayes", "nslphom_dual_w",
        "mvn_cdf_project_lp_FALSE", "mvn_pdf_project_lp_FALSE",
        "mult_project_lp_FALSE", "exact_project_lp_FALSE"
    )
    full_labels <- c("ei.md.bayes", "nslphom_dual_w", "mvn_cdf", "mvn_pdf", "mult", "exact")

    # Keep valid rows and coerce G to factor
    d <- df %>%
        dplyr::filter(method %in% full_order, !is.na(G), !is.na(C), !is.na(runtime)) %>%
        dplyr::mutate(G = as.factor(G))

    # Map internal method ids to display labels with a stable order
    present <- intersect(full_order, unique(d$method))
    method_lut <- stats::setNames(full_labels[match(present, full_order)], present)
    d <- d %>%
        dplyr::mutate(method = factor(method, levels = present, labels = unname(method_lut)))

    # Build a stable, explicit discrete x scale:
    # - convert C to integer then character to avoid lexical sort issues
    # - fix levels so all facets share the same x positions
    all_C_values <- sort(unique(as.integer(d$C)))
    all_C_levels <- as.character(all_C_values)
    d <- d %>%
        dplyr::mutate(C = factor(as.character(C), levels = all_C_levels))

    # Aggregate: mean runtime per (G, C, method)
    d_mean <- d %>%
        dplyr::group_by(G, C, method) %>%
        dplyr::summarise(runtime = mean(runtime, na.rm = TRUE), .groups = "drop")

    # Color palette and shapes per method
    oi <- c(
        "#0072B2", "#E69F00", "#009E73", "#D55E00",
        "#CC79A7", "#F0E442", "#56B4E9", "#999999"
    )
    lvls <- levels(d$method)
    pal <- stats::setNames(oi[seq_along(lvls)], lvls)
    shape_vals <- c(21, 24, 22, 25, 23, 21)[seq_along(lvls)]

    # y-axis labels: minimal decimals where significant
    fmt_min_decimals <- function(vals) {
        vapply(vals, function(x) {
            if (is.na(x)) {
                return(NA_character_)
            }
            if (x >= 1) {
                return(formatC(x, format = "f", digits = 0, big.mark = ""))
            }
            s <- formatC(x, format = "f", digits = 6, drop0trailing = TRUE)
            if (substr(s, 1, 1) == ".") s <- paste0("0", s)
            s
        }, character(1))
    }

    # Global y-limits; keep strictly positive for log scale
    min_pos <- min(d_mean$runtime[d_mean$runtime > 0], na.rm = TRUE)
    max_pos <- max(d_mean$runtime, na.rm = TRUE)
    ylims <- c(min_pos * 0.9, max_pos * 1.1)

    # Breaks: powers of ten only within observed range
    pow_seq <- seq(floor(log10(ylims[1])), ceiling(log10(ylims[2])), by = 1)
    candidate_breaks <- 10^pow_seq
    y_breaks <- candidate_breaks[candidate_breaks >= ylims[1] & candidate_breaks <= ylims[2]]

    # Plot
    p <- ggplot2::ggplot(
        d_mean,
        ggplot2::aes(x = C, y = runtime, fill = method, shape = method)
    ) +
        # Center points exactly on category ticks (no dodge)
        ggplot2::geom_point(
            size = 2, # stroke = 0.9,
            position = ggplot2::position_identity(),
            color = "black"
        ) +
        ggplot2::scale_shape_manual(values = shape_vals) +
        ggplot2::scale_fill_manual(values = pal) +
        # Discrete x scale:
        # - fix limits/breaks/labels so every tick (1..10..) appears
        # - small additive padding so points at the ends (e.g., C=1, C=10) are not clipped
        ggplot2::scale_x_discrete(
            limits = all_C_levels,
            breaks = all_C_levels,
            labels = all_C_levels,
            drop   = FALSE,
            expand = ggplot2::expansion(add = 0.23, mult = 0)
        ) +
        # Log y scale with powers-of-ten only; no minor breaks/lines
        ggplot2::scale_y_log10(
            limits = ylims,
            breaks = y_breaks,
            labels = fmt_min_decimals,
            minor_breaks = NULL
        ) +
        ggplot2::facet_grid(
            cols = ggplot2::vars(G),
            labeller = ggplot2::label_both,
            scales = "fixed"
        ) +
        ggplot2::labs(
            x = "Candidates (C)",
            y = "Average Runtime [s]",
            fill = "Method",
            shape = "Method"
        ) +
        ggplot2::theme_bw(base_size = letter_size, base_family = base_family) +
        ggplot2::theme(
            legend.position = "bottom",
            panel.spacing.x = grid::unit(12, "pt"),
            strip.text = ggplot2::element_text(face = "bold"),
            panel.grid.minor.y = ggplot2::element_blank(), # no minor grid lines
            panel.grid.major.y = ggplot2::element_line(size = 0.4), # keep only major (powers of ten)
            axis.text.x = element_text(size = letter_size - 2),
            axis.text.y = element_text(size = letter_size - 2),
            axis.title.x = element_text(size = letter_size - 2), # ← X-axis name
            axis.title.y = element_text(size = letter_size - 2) # ← Y-axis name
        ) +
        # Prevent clipping at panel edges after adding x padding
        ggplot2::coord_cartesian(clip = "off") +
        # Ensure x guide does not drop labels due to potential overlap
        ggplot2::guides(x = ggplot2::guide_axis(n.dodge = 1))

    # Save and print
    dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
    ggsave(file.path(save_dir, save_name), p,
        width = 6.3, height = 4, dpi = 300, device = grDevices::cairo_pdf
    )
    # paper -> width = 6.3, height = 4
    # presentation -> width = 7.3, height = 5
    print(p)
}
# ------------------ MAIN ------------------
main <- function() {
    INSTANCE_DIR <- "output/simulated_instances"
    FIGURES_DIR <- "figures"

    message("Reading instances…")
    result <- read_simulated_instances(INSTANCE_DIR)

    plot_runtime_vsC_byG(
        result$df,
        save_dir = FIGURES_DIR,
        save_name = "sim_runtime_3.pdf"
        # letter_size = 13
    )
}

main()
