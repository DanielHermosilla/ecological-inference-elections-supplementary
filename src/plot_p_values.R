#!/usr/bin/env Rscript

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
arg_chr <- function(pos, default = NULL) {
    if (length(args) >= pos) args[pos] else default
}
get_flag <- function(flag, default = NULL) {
    m <- grep(paste0("^--", flag, "="), args, value = TRUE)
    if (length(m)) sub(paste0("^--", flag, "="), "", m[1]) else default
}
stop_if <- function(cond, msg) {
    if (isTRUE(cond)) stop(msg, call. = FALSE)
}

# Posicional:
# 1) outfile (ej: figures/eiv_facets_methods3.pdf)
outfile <- arg_chr(1, "figures/eiv_facets_methods3.pdf")

# Flags:
field_to_plot <- get_flag("field", "EI_V")
inst_like <- get_flag("inst-like", "ei_")
limit_inst <- suppressWarnings(as.integer(get_flag("limit-inst", 8)))
workers_arg <- suppressWarnings(as.integer(get_flag("workers", NA_integer_)))

# ================= Paths =================
base_ei <- file.path("output", "ei_instances")
stop_if(
    !dir.exists(base_ei),
    sprintf("'%s' doesn't exist. Run from project root.", base_ei)
)

# ================= Instancias (máximo 8) =================
all_instances <- sort(list.dirs(base_ei, recursive = FALSE, full.names = FALSE))
instances <- all_instances[str_detect(all_instances, paste0("^", inst_like))]
stop_if(
    length(instances) == 0,
    sprintf("No instances starting with '%s' found in '%s'.", inst_like, base_ei)
)
if (is.finite(limit_inst) && limit_inst > 0) instances <- head(instances, limit_inst)
instances <- head(instances, 8)

# ================= Métodos (sólo los 3 pedidos) =================
method_ids <- c("mvn_cdf_project_lp_FALSE_sym", "mvn_pdf_project_lp_FALSE_sym", "mult_project_lp_FALSE_sym")
method_labs <- c("mvn_cdf", "mvn_pdf", "mult")
names(method_labs) <- method_ids
method_order_display <- method_labs[method_ids] # c("mvn_cdf","mvn_pdf","mult")

# ================= Lectura JSON + caché =================
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

# Promedio por distrito (para IC tomaremos estos promedios como muestras)
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

# ================= Paralelización por (inst, method) =================
grid <- expand.grid(inst = instances, method = method_ids, stringsAsFactors = FALSE)

n_cores <- if (is.finite(workers_arg) && workers_arg > 0) {
    workers_arg
} else {
    max(1L, parallel::detectCores(logical = TRUE) - 1L)
}
cl <- parallel::makeCluster(n_cores)
on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
doParallel::registerDoParallel(cl)
parallel::clusterExport(
    cl,
    varlist = c(
        "base_ei", "instances", "method_ids", "field_to_plot",
        "cache_env", "read_field_from_json", "get_cached_field", "by_district_for"
    ),
    envir = environment()
)
parallel::clusterEvalQ(cl, {
    library(dplyr)
    library(tibble)
    NULL
})

district_means_list <- foreach(i = seq_len(nrow(grid)), .packages = c("dplyr", "tibble")) %dopar% {
    by_district_for(grid$inst[i], grid$method[i], field_to_plot)
}
district_means <- bind_rows(district_means_list)

# ================= Resumen (media + IC 95%) por (inst, method) =================
# IC 95%: mean ± t_{0.975, n-1} * sd/sqrt(n)
summ <- district_means %>%
    filter(method %in% method_ids, is.finite(mean_val)) %>%
    mutate(method = factor(method, levels = method_ids, labels = method_labs[method_ids])) %>%
    group_by(inst, method) %>%
    summarise(
        n = n(),
        mean = mean(mean_val),
        sd = sd(mean_val),
        se = sd / sqrt(n),
        tcrit = ifelse(n > 1, qt(0.975, df = n - 1), NA_real_),
        ci_low = ifelse(n > 1, mean - tcrit * se, NA_real_),
        ci_high = ifelse(n > 1, mean + tcrit * se, NA_real_),
        .groups = "drop"
    ) %>%
    # Mantener orden de métodos solicitado
    mutate(method = factor(method, levels = method_order_display))

# ================= Estilo tipo tu script =================
# Fuente Fira Sans si está
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

# Paleta y shapes similares
oi <- c(
    "#0072B2", "#E69F00", "#009E73", "#D55E00",
    "#CC79A7", "#F0E442", "#56B4E9", "#999999"
)
lvls <- levels(summ$method)
pal <- stats::setNames(oi[seq_along(lvls)], lvls)
shape_vals <- c(21, 24, 22, 25, 23, 21)[seq_along(lvls)]

# ================= Plot =================
p <- ggplot(summ, aes(x = method, y = mean, fill = method, shape = method)) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.25, alpha = 0.7) +
    geom_point(size = 3.8, stroke = 0.9, color = "black") +
    scale_fill_manual(values = pal) +
    scale_shape_manual(values = shape_vals) +
    facet_wrap(~inst, ncol = 4, nrow = 2, scales = "free_y") +
    labs(x = "Método", y = field_to_plot, fill = "Método", shape = "Método") +
    theme_bw(base_size = 18, base_family = base_family) +
    theme(
        legend.position = "right",
        strip.text = element_text(face = "bold"),
        panel.grid.minor.y = element_blank()
    )

# ================= Guardar =================
dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)
ext <- tolower(tools::file_ext(outfile))
if (ext %in% c("png", "jpg", "jpeg", "tiff", "bmp")) {
    ggsave(outfile, plot = p, width = 16, height = 8.5, dpi = 300)
} else if (ext %in% c("pdf")) {
    ggsave(outfile, plot = p, width = 16, height = 8.5, device = grDevices::cairo_pdf)
} else {
    out_pdf <- sub("\\.[A-Za-z0-9]+$", ".pdf", outfile)
    ggsave(out_pdf, plot = p, width = 16, height = 8.5, device = grDevices::cairo_pdf)
    message(sprintf("Extensión no reconocida; guardado como PDF: %s", out_pdf))
}

message(sprintf("[OK] Figura guardada en: %s", normalizePath(outfile, mustWork = FALSE)))
