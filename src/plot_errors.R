suppressPackageStartupMessages({
    library(jsonlite)
    library(dplyr)
    library(purrr)
    library(stringr)
    library(forcats)
    library(ggplot2)
})
source("src/utilities.R")

# ------------------ PLOT ------------------
plot_EM_error_boxplot_R <- function(
    df,
    save_dir = "figures",
    save_name = "fig2-p-error.pdf",
    # --- NUEVOS PARÁMETROS ---
    groups_keep = NULL, # vector con grupos G a mantener, p.ej. c("1","2","3")
    cand_keep = NULL, # vector con candidatos C a mantener (override a 'candidates')
    candidates = c("all", "uar"), # "all" = todos los candidatos; "uar" = muestreo uniforme al azar
    n_cand = NULL, # número de candidatos a muestrear si candidates = "uar"
    seed = 123, # semilla para UAR
    letter_size = 6) {
    candidates <- match.arg(candidates)

    # ---- Fuente Fira Sans (opcional) ----
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

    # ---- Copia local ----
    d <- df

    # 1) Métodos a mostrar (y etiquetas)
    full_order <- c(
        "ei.MD.bayes", "nslphom_dual_w",
        "mvn_cdf_project_lp_FALSE", "mvn_pdf_project_lp_FALSE",
        "mult_project_lp_FALSE", "exact_project_lp_FALSE"
    )
    full_labels <- c(
        "ei.md.bayes", "nslphom_dual_w",
        "mvn_cdf", "mvn_pdf", "mult", "exact"
    )

    # 2) Qué métodos están realmente en d
    present <- intersect(full_order, unique(d$method))

    # 3) Filtrar SOLO métodos presentes y quitar NAs problemáticos
    d <- d %>%
        dplyr::filter(method %in% present) %>%
        dplyr::filter(!is.na(G), !is.na(C)) %>%
        dplyr::filter(!is.na(MAE_p2))

    # ---- Filtros por grupos y candidatos (NUEVO) ----
    # Normalizar a character para filtrar de forma robusta
    d <- d %>%
        mutate(
            G = as.character(G),
            C = as.character(C)
        )

    # 3a) Filtrar por grupos si se indicó
    if (!is.null(groups_keep)) {
        groups_keep <- as.character(groups_keep)
        d <- d %>% filter(G %in% groups_keep)
        if (nrow(d) == 0) stop("Tras filtrar por 'groups_keep' no quedan filas.")
    }

    # 3b) Filtrar por candidatos:
    #     - Si se pasa 'cand_keep', se usa tal cual.
    #     - Si no, y candidates == "uar", se muestrean 'n_cand' candidatos UAR.
    #     - Si candidates == "all", no se filtra por C.
    cand_levels_all <- sort(unique(d$C))

    if (!is.null(cand_keep)) {
        cand_keep <- as.character(cand_keep)
        missing_c <- setdiff(cand_keep, cand_levels_all)
        if (length(missing_c)) {
            warning(
                "Algunos candidatos de 'cand_keep' no existen en los datos y se ignorarán: ",
                paste(missing_c, collapse = ", ")
            )
        }
        cand_keep_in <- intersect(cand_keep, cand_levels_all)
        if (!length(cand_keep_in)) stop("Ninguno de los candidatos de 'cand_keep' existe en los datos.")
        d <- d %>% filter(C %in% cand_keep_in)
    } else if (candidates == "uar") {
        if (is.null(n_cand) || !is.numeric(n_cand) || n_cand < 1) {
            stop("Para candidates = 'uar' debes especificar un 'n_cand' >= 1.")
        }
        if (n_cand > length(cand_levels_all)) {
            warning("n_cand > número de candidatos disponibles; se usará n_cand = ", length(cand_levels_all))
            n_cand <- length(cand_levels_all)
        }
        set.seed(seed)
        sampled_c <- sample(cand_levels_all, n_cand, replace = FALSE)
        d <- d %>% filter(C %in% sampled_c)
    } # si candidates == "all": no filtramos por C

    # Volver a factorizar para facet ordenado/limpio
    # --- NIVELES DE C ORDENADOS ---
    if (!is.null(cand_keep)) {
        # Usa el orden en que pasaste cand_keep
        c_levels <- cand_keep_in
    } else {
        # Si no hay cand_keep, ordena numéricamente
        c_levels <- sort(unique(as.numeric(d$C)))
        c_levels <- as.character(c_levels)
    }

    d <- d %>%
        mutate(
            G = factor(G, levels = sort(unique(G))),
            C = factor(C, levels = c_levels)
        )

    # 4) Recodificar etiquetas de métodos
    method_lut <- stats::setNames(full_labels[match(present, full_order)], present)
    d <- d %>%
        mutate(method = factor(method, levels = present, labels = unname(method_lut)))

    # ---- Paleta (Okabe–Ito) ----
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

    # Retornar invisiblemente lo filtrado y el plot por si quieres reusar
    invisible(list(data = d, plot = p))
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
    INSTANCE_DIR <- "output/simulated_instances" # relativo a code/
    FIGURES_DIR <- "figures"

    message("Leyendo instancias…")
    result <- read_simulated_instances(INSTANCE_DIR)
    if (!"MAE_p2" %in% names(result$df)) {
        mae_vec <- vapply(result$raw, compute_mae_p2, numeric(1))
        result$df$MAE_p2 <- mae_vec
    }

    # EJEMPLOS DE USO:
    # 1) Mantener grupos "1","2","3" y tomar 3 candidatos UAR
    # plot_EM_error_boxplot_R(
    #     result$df,
    #     save_dir = FIGURES_DIR,
    #     save_name = "sim_prob_error_g123_cUAR3.pdf",
    #     groups_keep = c("1", "2", "3"),
    #     candidates = "uar",
    #     n_cand = 3,
    #     seed = 2025
    # )

    # 2) Mantener grupos "1","3" y candidatos específicos "A","B","C" (si existen)
    plot_EM_error_boxplot_R(
        result$df,
        save_dir = FIGURES_DIR,
        save_name = "sim_prob_error_3.pdf",
        groups_keep = c("2", "3", "4", "6", "8"),
        cand_keep = c("2", "3", "5", "10"),
        letter_size = 6
    )

    # 3) Sólo filtrar grupos y dejar "todos" los candidatos presentes
    # plot_EM_error_boxplot_R(
    #     result$df,
    #     save_dir  = FIGURES_DIR,
    #     save_name = "sim_prob_error_g12_allC.pdf",
    #     groups_keep = c("1", "2"),
    #     candidates  = "all"
    # )
}

main()
