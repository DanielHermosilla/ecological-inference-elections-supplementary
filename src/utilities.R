suppressPackageStartupMessages({
    library(jsonlite)
    library(dplyr)
})

# Reads simulated instances, returns a data.frame with summary metrics and a list with raw JSON content.
read_simulated_instances <- function(
    base_dir = "output/simulated_instances",
    include_details = FALSE, # if TRUE, keeps detail_*.json
    verbose = TRUE) {
    if (!dir.exists(base_dir)) {
        stop("The directory doesn't exist: ", normalizePath(base_dir))
    }

    # List JSON files
    files <- list.files(base_dir, pattern = "\\.json$", recursive = TRUE, full.names = TRUE)
    if (!include_details) {
        files <- files[!grepl("/detail_\\d+\\.json$", files)]
    }
    if (length(files) == 0) stop("A json file wasn't found in: ", normalizePath(base_dir))

    # Methods (folder names)
    methods <- c(
        "mvn_cdf_project_lp_FALSE", "mvn_pdf_project_lp_FALSE", "mult_project_lp_FALSE", "exact_project_lp_FALSE",
        "lphom", "lclphom", "nslphom_dual_w", "lphom_joint", "ei.MD.bayes"
    )

    # Keep files like .../I*_B*_G*_C*_lambda*/<method>/<n>.json
    parent_dir <- basename(dirname(files))
    keep <- parent_dir %in% methods & grepl("/[0-9]+\\.json$", files)
    files <- files[keep]
    parent_dir <- parent_dir[keep]
    if (length(files) == 0) {
        stop(
            "There aren't any JSON that match '.../<method>/<n>.json' under: ",
            normalizePath(base_dir)
        )
    }

    # Parse params from folder name I*_B*_G*_C*_lambda*
    parse_params <- function(path) {
        prm_dir <- basename(dirname(dirname(path)))
        m <- regexec("^I(\\d+)_B(\\d+)_G(\\d+)_C(\\d+)_lambda(\\d+)$", prm_dir)
        g <- regmatches(prm_dir, m)[[1]]
        if (length(g) == 6) {
            list(
                I = as.integer(g[2]),
                B = as.integer(g[3]),
                G = as.integer(g[4]),
                C = as.integer(g[5]),
                lambda = as.integer(g[6])
            )
        } else {
            list(
                I = NA_integer_, B = NA_integer_, G = NA_integer_,
                C = NA_integer_, lambda = NA_integer_
            )
        }
    }

    # Helpers (silent)
    scalarize_num <- function(val) {
        if (is.null(val) || length(val) == 0L) {
            return(NA_real_)
        }
        if (is.list(val)) val <- unlist(val, recursive = TRUE, use.names = FALSE)
        as.numeric(val[[1]])
    }
    probs_to_vec <- function(z) {
        if (is.null(z)) {
            return(NULL)
        }
        if (is.list(z) && !is.matrix(z)) z <- unlist(z, recursive = TRUE, use.names = FALSE)
        as.numeric(z)
    }

    rows <- vector("list", length(files))
    raw <- vector("list", length(files))
    names(raw) <- files

    for (i in seq_along(files)) {
        f <- files[i]
        meth <- parent_dir[i]
        prm <- parse_params(f)

        # lectura normal (si el JSON está malformado, dejará error — sin prints)
        txt <- readLines(f, warn = FALSE)
        x <- jsonlite::fromJSON(paste(txt, collapse = "\n"), simplifyVector = TRUE)

        # Guardar crudo
        raw[[i]] <- list(file = f, method = meth, params = prm, json = x)

        # Tomar escalar o NA en cada métrica
        EI_Z <- scalarize_num(x$EI_Z)
        runtime <- scalarize_num(x$runtime)
        EI_V <- scalarize_num(x$EI_V)
        MAE_Z <- scalarize_num(x$MAE_Z)
        MAE_V <- scalarize_num(x$MAE_V)
        MAE_p <- scalarize_num(x$MAE_p)

        # Si no viene MAE_p, lo calculamos desde p_true/p_est (si son conformables)
        if (is.na(MAE_p)) {
            p_true <- probs_to_vec(x$p_true)
            p_est <- probs_to_vec(x$p_est)
            if (!is.null(p_true) && !is.null(p_est) && length(p_true) == length(p_est)) {
                MAE_p <- mean(abs(p_est - p_true))
            }
        }

        rows[[i]] <- data.frame(
            file = f,
            method = meth,
            I = prm$I,
            B = prm$B,
            G = prm$G,
            C = prm$C,
            lambda = prm$lambda,
            EI_Z = EI_Z,
            EI_V = EI_V,
            MAE_Z = MAE_Z,
            MAE_V = MAE_V,
            MAE_p = MAE_p,
            runtime = runtime,
            status = 0L,
            stringsAsFactors = FALSE
        )
    }

    df <- dplyr::bind_rows(rows)

    if (verbose) {
        message(sprintf(
            "Read: %d rows | methods: %s",
            nrow(df),
            paste(sort(unique(df$method)), collapse = ", ")
        ))
    }

    list(df = df, raw = raw)
}
