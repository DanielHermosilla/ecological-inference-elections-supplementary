# Core utilities for generating instances, running EI methods, and summarising results.
cat("\014") # clear console

# ======== Paquetes ========
suppressPackageStartupMessages({
    # devtools::install_github("DanielHermosilla/ecological-inference-elections")
    # remove.packages("fastei"); install.packages("fastei")
    library(fastei)
    library(jsonlite) # to save a named list as a json file
    library(dplyr)
    library(eiPack)
    library(ei.Datasets) # load voting datasets
    library(lphom) # benchmark library
    library(ecolRxC) # benchmark library (si lo usas)
    library(microbenchmark)
})

# ======== Paths ========
PATH_SIM_INSTANCES <- file.path("output", "simulated_instances")
PATH_EI_INSTANCES <- file.path("output", "ei_instances")

# ======== Utils ========

# Create folder recursively when missing (no-op otherwise).
create_folder <- function(folder_name, verbose = FALSE) {
    if (!dir.exists(folder_name)) {
        dir.create(folder_name, recursive = TRUE, showWarnings = FALSE)
        if (verbose) cat("Folder created.\n")
    } else if (verbose) {
        cat("Folder already exists.\n")
    }
}

# Sanitize strings for use in file/folder names.
sanitize <- function(x) gsub("[^A-Za-z0-9._-]", "_", as.character(x))

# Build simulated-instance folder name.
folder_name_simulated_instance <- function(I, B, G, C, lambda = 0.5) {
    paste0("I", I, "_B", B, "_G", G, "_C", C, "_lambda", round(100 * lambda))
}

# Build method-specific folder name (encoding options).
folder_name_method <- function(method, mcmc_samples = NULL,
                               adjust_prob_cond_method = NULL,
                               adjust_prob_cond_every = NULL,
                               invert = NULL, symmetric = NULL) {
    paste0(
        method,
        ifelse(method == "mcmc", paste0("_", mcmc_samples), ""),
        ifelse(!is.null(adjust_prob_cond_method), paste0("_", adjust_prob_cond_method), ""),
        ifelse(!is.null(adjust_prob_cond_every), paste0("_", adjust_prob_cond_every), ""),
        ifelse(!is.null(invert) && invert, "_inv", ""), ifelse(!is.null(symmetric) && symmetric, "_sym", "")
    )
}

# Build per-seed JSON filename.
file_name <- function(seed, seed_initial = -123, file_type = ".json") {
    paste0(seed, ifelse(seed_initial == -123, "", paste0("_pinit", seed_initial)), file_type)
}

# Skip heavy exact cases.
fun_skip_cases <- function(method, G, C, skip_cases = TRUE) {
    skip_cases & (method == "exact") & (((G > 2) & (C == 5)) | (C > 5))
}

# EI and MAE helpers (return NULL if missing).
compute_EI <- function(x_est, x_true) {
    if (is.null(x_est)) {
        return(NULL)
    }
    100 * 0.5 * sum(abs(x_est - x_true)) / sum(x_true)
}
compute_MAE <- function(x_est, x_true) {
    if (is.null(x_est)) {
        return(NULL)
    }
    mean(abs(x_est - x_true))
}

# Runtime in seconds (portable wall-clock).
nano_time <- function(expr) {
    t0 <- proc.time()[["elapsed"]]
    force(expr)
    t1 <- proc.time()[["elapsed"]]
    t1 - t0
}

# Reads simulated instance JSON outputs and returns summary metrics plus raw content.
read_simulated_instances <- function(
    base_dir = "output/simulated_instances",
    include_details = FALSE,
    verbose = TRUE) {
    if (!dir.exists(base_dir)) {
        stop("The directory doesn't exist: ", normalizePath(base_dir))
    }

    files <- list.files(base_dir, pattern = "\\.json$", recursive = TRUE, full.names = TRUE)
    if (!include_details) {
        files <- files[!grepl("/detail_\\d+\\.json$", files)]
    }
    if (length(files) == 0) {
        stop("A json file wasn't found in: ", normalizePath(base_dir))
    }

    methods <- c(
        "mvn_cdf_project_lp_FALSE", "mvn_pdf_project_lp_FALSE",
        "mult_project_lp_FALSE", "exact_project_lp_FALSE",
        "lphom", "lclphom", "nslphom_dual_w", "lphom_joint", "ei.MD.bayes"
    )

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
            list(I = NA_integer_, B = NA_integer_, G = NA_integer_, C = NA_integer_, lambda = NA_integer_)
        }
    }

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

        txt <- readLines(f, warn = FALSE)
        x <- jsonlite::fromJSON(paste(txt, collapse = "\n"), simplifyVector = TRUE)

        raw[[i]] <- list(file = f, method = meth, params = prm, json = x)

        EI_Z <- scalarize_num(x$EI_Z)
        runtime <- scalarize_num(x$runtime)
        EI_V <- scalarize_num(x$EI_V)
        MAE_Z <- scalarize_num(x$MAE_Z)
        MAE_V <- scalarize_num(x$MAE_V)
        MAE_p <- scalarize_num(x$MAE_p)

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

# Validate method argument early for clearer errors.
validate_method <- function(method, method_variable_name = "method", upper_function_name = "function") {
    valid_method <- c(
        "mcmc", "mvn_cdf", "mvn_pdf", "mult", "exact",
        "nslphom_dual_w", "lphom", "lclphom", "ecolRxC",
        "nslphom_joint", "nslphom_dual_a", "ei.MD.bayes"
    )
    if (!(method %in% valid_method)) {
        stop(sprintf(
            "The variable '%s' of the function '%s' takes value '%s', but should be one of these: %s",
            method_variable_name, upper_function_name, method, paste(valid_method, collapse = ", ")
        ))
    }
}

# Convert ei.MD.bayes output to lphom-style format.
ei.MD.bayes2lphom <- function(ei.object, votes_election1 = NULL, votes_election2 = NULL) {
    X <- ei.object$draws$Cell.counts
    if (length(dim(X)) != 3) stop("Recompute ei.MD.bayes output with the option ret.mcmc = FALSE")
    X <- apply(X, c(1, 2), mean)
    pjk <- X / rowSums(X)
    pkj <- round(X / colSums(X) * 100, 2)
    pjk <- round(pjk * 100, 2)
    list("VTM" = pjk, "OTM" = pkj)
}

# compute_loglik_pdf <- function(X, W, P, ridge = 1e-8) {
#     X <- as.matrix(X)
#     W <- as.matrix(W)
#     P <- as.matrix(P)
#
#     B <- nrow(X)
#     C <- ncol(X)
#     G <- ncol(W)
#     n <- C - 1
#
#     P_red <- P[, 1:n, drop = FALSE]
#
#     getParams_R <- function(b, W, P_red) {
#         w_b <- W[b, ]
#         mu <- as.numeric(crossprod(w_b, P_red))
#         diagW <- diag(w_b, nrow = G)
#         temp <- t(P_red) %*% diagW
#         Sigma <- temp %*% P_red
#
#         for (j in seq_len(n)) {
#             for (i in seq_len(n)) {
#                 if (i == j) {
#                     Sigma[i, j] <- mu[i] - Sigma[i, j]
#                 } else {
#                     Sigma[i, j] <- -Sigma[i, j]
#                 }
# ======== Core: run one EI method (real or simulated) ========
run_ecological_inference_algorithm <- function(
    W = NULL, X = NULL, Z_true = NULL, V_true = NULL, p_true = NULL,
    method,
    adjust_prob_cond_method = NULL,
    adjust_prob_cond_every = NULL,
    maxtime,
    invert = FALSE,
    symmetric = FALSE,
    initial_prob = NULL,
    seed_initial,
    mcmc_samples,
    save_json = TRUE,
    data_source_folder,
    instance_folder,
    instance_seed = NULL) {
    validate_method(method, "method", "run_ecological_inference_algorithm")

    if (!is.null(Z_true) && is.null(V_true)) V_true <- apply(Z_true, c(1, 2), sum)
    if (!is.null(Z_true) && (is.null(W) || is.null(X))) {
        X <- t(apply(Z_true, c(2, 3), sum)) # B x C
        W <- t(apply(Z_true, c(1, 3), sum)) # B x G
    }
    if (is.null(p_true) && !is.null(V_true)) p_true <- V_true / rowSums(V_true)

    if (method %in% c("mcmc", "mvn_cdf", "mvn_pdf", "mult", "exact")) {
        initial_prob <- ifelse(seed_initial == -123, "group_proportional", "random")
        runtime <- nano_time({
            output <- run_em(
                X = X, W = W, method = method, maxtime = maxtime,
                initial_prob = initial_prob, seed = seed_initial,
                mcmc_samples = mcmc_samples,
                compute_ll = if (method == "mvn_cdf") TRUE else FALSE,
                adjust_prob_cond_method = adjust_prob_cond_method,
                verbose = FALSE,
                mvncdf_error = 1e-3,
                adjust_prob_cond_every = if (!is.null(adjust_prob_cond_every)) adjust_prob_cond_every else FALSE, symmetric = symmetric
            )
        })

        Z_est <- output$expected_outcome
        V_est <- apply(Z_est, c(1, 2), sum)
        p_est <- output$prob
        cat(sprintf("Total iterations: %d\n", output$iterations))
    } else if (method == "nslphom_dual_w") {
        runtime <- nano_time({
            output <- nslphom_dual(W, X)
        })
        V_est <- output$VTM.votes.w
        Z_est <- output$VTM.votes.units.w
        p_est <- output$VTM12.w
    } else if (method == "lphom") {
        runtime <- nano_time({
            output <- lphom(W, X)
        })
        V_est <- output$VTM.votes
        Z_est <- output$VTM.votes.units
        p_est <- output$VTM.complete
    } else if (method == "lclphom") {
        runtime <- nano_time({
            output <- lclphom(W, X)
        })
        V_est <- output$VTM.votes
        Z_est <- output$VTM.votes.units
        p_est <- output$VTM.complete
    } else if (method == "ecolRxC") {
        runtime <- nano_time({
            output <- ecolRxC(W, X)
        })
        V_est <- output$VTM.votes
        Z_est <- output$VTM.votes.units
        p_est <- output$VTM.global
    } else if (method == "nslphom_joint") {
        runtime <- nano_time({
            output <- nslphom_joint(W, X)
        })
        V_est <- output$VTM.votes
        Z_est <- output$VTM.votes.units
        p_est <- output$VTM12
    } else if (method == "nslphom_dual_a") {
        runtime <- nano_time({
            output <- nslphom_dual(W, X)
        })
        V_est <- output$VTM.votes.a
        Z_est <- output$VTM.votes.units.a
        p_est <- output$VTM12.a
    } else if (method == "ei.MD.bayes") {
        print(colnames(X))
        colnames(X) <- paste0("cand_", seq_len(ncol(X)))
        print(colnames(X))
        colnames(W) <- paste0("grp_", seq_len(ncol(W)))
        var.y <- paste(paste0("`", colnames(X), "`"), collapse = ",")
        var.x <- paste(paste0("`", colnames(W), "`"), collapse = ",")
        formula <- paste("cbind(", var.y, ") ~ cbind(", var.x, ")")
        X_df <- as.data.frame(X, stringsAsFactors = FALSE)
        W_df <- as.data.frame(W, stringsAsFactors = FALSE)
        base <- data.frame(W_df, X_df, check.names = FALSE)
        runtime <- nano_time({
            tune.nocov <- tuneMD(formula = eval(formula), data = base, ntunes = 10, totaldraws = 100000)
            output <- ei.MD.bayes(
                formula = eval(formula), data = base, tune.list = tune.nocov,
                sample = 1000, thin = 100, burnin = 100000, ret.mcmc = F
            )
        })
        output <- ei.MD.bayes2lphom(output, votes_election1 = X, votes_election2 = W)
        G <- ncol(W)
        B <- nrow(W)
        C <- ncol(X)
        p_est <- output$VTM / 100
        V_est <- sweep(p_est, 1, colSums(W), `*`)
        Z_est <- array(NA_real_, dim = c(G, C, B))
        for (b in seq_len(B)) Z_est[, , b] <- sweep(p_est, 1, W[b, ], `*`)
    } else if (method == "ecolRxC") {
        runtime <- nano_time({
            output <- ecolRxC(W, X)
        })
        V_est <- output$VTM.votes
        Z_est <- output$VTM.votes.units
        p_est <- output$VTM.global
    } else {
        stop("Method not recognized")
    }

    EI_Z <- if (!is.null(Z_true)) compute_EI(Z_est, Z_true) else NULL
    EI_V <- if (!is.null(V_true)) compute_EI(V_est, V_true) else NULL
    MAE_Z <- if (!is.null(Z_true)) compute_MAE(Z_est, Z_true) else NULL
    MAE_V <- if (!is.null(V_true)) compute_MAE(V_est, V_true) else NULL
    MAE_p <- if (!is.null(p_true)) compute_MAE(p_est, p_true) else NULL

    cat(sprintf(
        "\n- EI_Z: %s, EI_V: %s, MAE_Z: %s, MAE_V: %s, MAE_p: %s, runtime: %.2f seconds\n\n",
        if (is.null(EI_Z)) "NA" else format(EI_Z, digits = 4),
        if (is.null(EI_V)) "NA" else format(EI_V, digits = 4),
        if (is.null(MAE_Z)) "NA" else format(MAE_Z, digits = 4),
        if (is.null(MAE_V)) "NA" else format(MAE_V, digits = 4),
        if (is.null(MAE_p)) "NA" else format(MAE_p, digits = 4),
        runtime
    ))

    if (isTRUE(save_json)) {
        root <- if (identical(data_source_folder, "simulated_instances")) {
            PATH_SIM_INSTANCES
        } else if (identical(data_source_folder, "ei_instances")) {
            PATH_EI_INSTANCES
        } else {
            stop("data_source_folder must be 'simulated_instances' or 'ei_instances'")
        }

        method_folder <- folder_name_method(method, mcmc_samples, adjust_prob_cond_method, adjust_prob_cond_every, invert, symmetric)
        out_dir <- file.path(root, instance_folder, method_folder)
        if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

        seed_for_name <- if (is.null(instance_seed)) "1" else instance_seed
        pj_name <- file_name(seed = seed_for_name, seed_initial = seed_initial, file_type = ".json")

        jsonlite::write_json(
            list(
                Z_true = Z_true, V_true = V_true, p_true = p_true,
                Z_est = Z_est, V_est = V_est, p_est = p_est,
                EI_Z = EI_Z, EI_V = EI_V, MAE_Z = MAE_Z, MAE_V = MAE_V, MAE_p = MAE_p,
                runtime = runtime
            ),
            path = file.path(out_dir, pj_name),
            pretty = TRUE, auto_unbox = TRUE
        )

        if (method %in% c("mcmc", "mvn_cdf", "mvn_pdf", "mult", "exact")) {
            base_name <- tools::file_path_sans_ext(pj_name)
            save_eim(output, filename = file.path(out_dir, paste0("detail_", base_name, ".json")))
        }
    }

    invisible(list(
        Z_true = Z_true, V_true = V_true, p_true = p_true,
        Z_est = Z_est, V_est = V_est, p_est = p_est,
        EI_Z = EI_Z, EI_V = EI_V, MAE_Z = MAE_Z, MAE_V = MAE_V, MAE_p = MAE_p
    ))
}

# ======== Worker top-level (no closure) ========
run_one_ei_task <- function(tk, base_map, invert, symmetric, maxtime, save_json) {
    method <- tk$method
    adjust_prob_cond_method <- tk$adjust_prob_cond_method
    adjust_prob_cond_every <- tk$adjust_prob_cond_every
    mcmc_samples <- tk$mcmc_samples
    base_name <- tk$base_name
    i <- tk$i
    seed_initial <- tk$seed_initial

    df <- base_map[[base_name]]
    cruzada <- df$District_cross_votes[[i]]
    candidatos <- df$Votes_to_candidates[[i]]
    names(candidatos) <- paste0(names(candidatos), ".c")
    partidos <- df$Votes_to_parties[[i]]
    names(partidos) <- paste0(names(partidos), ".p")
    district0 <- df$District[i]

    V_true <- as.matrix(cruzada[, -1, drop = FALSE]) # G x C
    p_true <- sweep(V_true, 1, rowSums(V_true), "/")

    if (!invert) {
        X <- as.matrix(candidatos[, -c(1, 2), drop = FALSE]) # B x C
        W <- as.matrix(partidos[, -c(1, 2), drop = FALSE]) # B x G
    } else {
        W <- as.matrix(candidatos[, -c(1, 2), drop = FALSE]) # B x C
        X <- as.matrix(partidos[, -c(1, 2), drop = FALSE]) # B x G
        V_true <- t(V_true)
        p_true <- t(p_true)
    }

    G <- nrow(V_true)
    C <- ncol(V_true)
    cat(sprintf("- Base=%s\tDist=%s (G=%d, C=%d)\n", base_name, as.character(district0), G, C))

    run_ecological_inference_algorithm(
        W = W, X = X,
        V_true = V_true, p_true = p_true,
        method = method, maxtime = maxtime,
        adjust_prob_cond_method = adjust_prob_cond_method,
        adjust_prob_cond_every = adjust_prob_cond_every,
        seed_initial = seed_initial, mcmc_samples = mcmc_samples,
        save_json = save_json,
        invert = invert,
        symmetric = symmetric,
        data_source_folder = "ei_instances",
        instance_folder = file.path(
            paste0(sanitize(base_name)),
            paste0(sanitize(district0), "_G", G, "_C", C)
        )
    )
    invisible(NULL)
}

# ======== EI datasets (parallel with doParallel) ========
run_eidatasets_instances <- function(
    method_arr,
    adjust_prob_cond_arr,
    adjust_prob_cond_every_arr,
    mcmc_samples_arr,
    bases_arr,
    seed_initial_arr = -123,
    maxtime = 60,
    invert = FALSE,
    symmetric = FALSE,
    save_json = TRUE,
    parallel = FALSE,
    workers = NULL,
    reproducible_seed = NULL) {
    create_folder("output")
    create_folder(PATH_EI_INSTANCES)
    if (is.null(seed_initial_arr)) seed_initial_arr <- -123

    adjustable_methods <- c("mvn_cdf", "mvn_pdf", "mult", "mcmc", "exact")

    # Prepare bases once
    base_list <- lapply(bases_arr, function(base_name) {
        df <- merge_small_options(get(base_name), 3, 3) # defined in your base code
        list(name = base_name, df = df)
    })
    base_map <- setNames(
        lapply(base_list, `[[`, "df"),
        sapply(base_list, `[[`, "name")
    )

    # Tasks
    tasks <- list()
    for (method in method_arr) {
        samples_arr <- if (identical(method, "mcmc")) mcmc_samples_arr else 1000L
        adj_methods_vec <- if (method %in% adjustable_methods) adjust_prob_cond_arr else list(NULL)

        for (adjust_prob_cond_method in adj_methods_vec) {
            adj_every_vec_curr <- if (method %in% adjustable_methods && !identical(adjust_prob_cond_method, "")) {
                adjust_prob_cond_every_arr
            } else {
                list(NULL)
            }

            for (adjust_prob_cond_every in adj_every_vec_curr) {
                for (mcmc_samples in samples_arr) {
                    for (base_name in bases_arr) {
                        df <- base_map[[base_name]]
                        for (i in seq_len(nrow(df))) {
                            for (seed_initial in seed_initial_arr) {
                                tasks[[length(tasks) + 1L]] <- list(
                                    method = method,
                                    adjust_prob_cond_method = adjust_prob_cond_method,
                                    adjust_prob_cond_every = adjust_prob_cond_every,
                                    mcmc_samples = mcmc_samples,
                                    base_name = base_name,
                                    i = i,
                                    seed_initial = seed_initial
                                )
                            }
                        }
                    }
                }
            }
        }
    }
    if (!length(tasks)) {
        return(invisible(NULL))
    }

    if (parallel) {
        # Attach packages for %dopar%
        if (!"package:foreach" %in% search()) library(foreach)
        if (!"package:doParallel" %in% search()) library(doParallel)

        if (is.null(workers)) workers <- max(1L, parallel::detectCores(logical = TRUE) - 1L)
        cl <- parallel::makeCluster(workers)
        on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
        doParallel::registerDoParallel(cl)
        if (!is.null(reproducible_seed) && requireNamespace("doRNG", quietly = TRUE)) {
            doRNG::registerDoRNG(reproducible_seed)
        }

        # Export GLOBAL helpers (defined above in this script)
        parallel::clusterExport(
            cl,
            varlist = c(
                # workers/functions
                "run_one_ei_task",
                "run_ecological_inference_algorithm",
                # helpers llamados por run_ecological_inference_algorithm
                "validate_method", "compute_EI", "compute_MAE", "nano_time",
                "folder_name_method", "file_name", "ei.MD.bayes2lphom", "ecolRxC",
                # paths used when saving
                "PATH_SIM_INSTANCES", "PATH_EI_INSTANCES",
                # util
                "sanitize"
            ),
            envir = .GlobalEnv
        )
        # Export local objects from this function
        parallel::clusterExport(
            cl,
            varlist = c("base_map", "invert", "maxtime", "save_json"),
            envir = environment()
        )
        # Load packages on workers
        parallel::clusterEvalQ(cl, {
            library(fastei)
            library(jsonlite)
            library(lphom)
            library(eiPack)
            library(ei.Datasets)
            library(ecolRxC)
            NULL
        })

        foreach(
            idx = seq_along(tasks),
            .export = c("run_one_ei_task", "run_ecological_inference_algorithm", "sanitize"),
            .packages = character(0)
        ) %dopar% {
            run_one_ei_task(tasks[[idx]],
                base_map = base_map,
                invert = invert,
                maxtime = maxtime,
                symmetric = symmetric,
                save_json = save_json
            )
            NULL
        }
    } else {
        lapply(tasks, run_one_ei_task,
            base_map = base_map, invert = invert, symmetric = symmetric,
            maxtime = maxtime, save_json = save_json
        )
    }

    invisible(NULL)
}

# ======== Simulated runs (sequential; parallelize similarly if desired) ========
run_simulated_instances <- function(
    method_arr, mcmc_samples_arr,
    adjust_prob_cond_arr, adjust_prob_cond_every_arr,
    I_arr, B_arr, C_arr, G_arr, lambda_arr,
    seed_arr, seed_initial_arr,
    maxtime = 3600, skip_cases = TRUE, save_json = TRUE) {
    create_folder("output")
    create_folder(PATH_SIM_INSTANCES)
    if (is.null(seed_initial_arr) || (length(seed_initial_arr) == 1L && is.na(seed_initial_arr))) {
        seed_initial_arr <- -123L
    }
    seed_initial_arr <- vapply(seed_initial_arr, function(z) if (is.null(z)) -123L else as.integer(z), integer(1))

    adjustable_methods <- c("mvn_cdf", "mvn_pdf", "mult", "mcmc", "exact")

    for (method in method_arr) {
        samples_arr <- if (identical(method, "mcmc")) mcmc_samples_arr else 1000L
        adj_methods_vec <- if (method %in% adjustable_methods) adjust_prob_cond_arr else list(NULL)
        adj_every_vec <- if (method %in% adjustable_methods) adjust_prob_cond_every_arr else list(NULL)

        for (adjust_prob_cond_method in adj_methods_vec) {
            for (adjust_prob_cond_every in adj_every_vec) {
                cat(sprintf(
                    "----------------\nRunning method: %s\n- Adjust conditional probability method: %s\n- Adjust conditional probability every iteration: %s\n----------------\n",
                    method, adjust_prob_cond_method, adjust_prob_cond_every
                ))
                for (mcmc_samples in samples_arr) {
                    for (I in I_arr) {
                        for (B in B_arr) {
                            for (C in C_arr) {
                                for (G in G_arr) {
                                    for (lambda in lambda_arr) {
                                        if (fun_skip_cases(method, G, C, skip_cases)) next
                                        for (seed in seed_arr) {
                                            instance <- simulate_election(
                                                num_ballots = B, num_candidates = C, num_groups = G,
                                                ballot_voters = rep(I, B), seed = seed, lambda = lambda
                                            )
                                            Z_true <- instance$outcome
                                            p_true <- instance$real_prob
                                            instance_folder <- folder_name_simulated_instance(I, B, G, C, lambda)
                                            for (seed_initial in seed_initial_arr) {
                                                cat(sprintf(
                                                    "\n- B=%d C=%d G=%d lambda=%.2f seed=%d seed_init=%s",
                                                    B, C, G, lambda, seed, if (is.na(seed_initial)) "NA" else as.character(seed_initial)
                                                ))
                                                run_ecological_inference_algorithm(
                                                    X = instance$X, W = instance$W,
                                                    adjust_prob_cond_method = adjust_prob_cond_method,
                                                    adjust_prob_cond_every = adjust_prob_cond_every,
                                                    Z_true = Z_true, p_true = p_true,
                                                    method = method, maxtime = maxtime,
                                                    seed_initial = seed_initial, mcmc_samples = mcmc_samples,
                                                    save_json = save_json,
                                                    data_source_folder = "simulated_instances",
                                                    instance_folder = instance_folder,
                                                    instance_seed = seed
                                                )
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    invisible(NULL)
}

# ===================== Ejemplo de uso (comentar/ajustar) =====================
# run_eidatasets_instances(
#   method_arr = c("mvn_cdf"),
#   adjust_prob_cond_arr = list(""),
#   adjust_prob_cond_every_arr = list(NULL),
#   mcmc_samples_arr = 500L,
#   bases_arr = c("NombreObjetoBase"),
#   parallel = TRUE,
#   workers = 4,
#   reproducible_seed = 12345
# )
