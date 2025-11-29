cat("\014") # limpia consola

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

# crea carpeta si no existe
create_folder <- function(folder_name, verbose = FALSE) {
    if (!dir.exists(folder_name)) {
        dir.create(folder_name, recursive = TRUE, showWarnings = FALSE)
        if (verbose) cat("Folder created.\n")
    } else if (verbose) {
        cat("Folder already exists.\n")
    }
}

# cross-platform
sanitize <- function(x) gsub("[^A-Za-z0-9._-]", "_", as.character(x))

# nombre carpeta instancia simulada
folder_name_simulated_instance <- function(I, B, G, C, lambda = 0.5) {
    paste0("I", I, "_B", B, "_G", G, "_C", C, "_lambda", round(100 * lambda))
}

# nombre carpeta def método
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

# nombre archivo json
file_name <- function(seed, seed_initial = -123, file_type = ".json") {
    paste0(seed, ifelse(seed_initial == -123, "", paste0("_pinit", seed_initial)), file_type)
}

# saltar casos costosos para exact
fun_skip_cases <- function(method, G, C, skip_cases = TRUE) {
    skip_cases & (method == "exact") & (((G > 2) & (C == 5)) | (C > 5))
}

# EI y MAE
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

# runtime en segundos (portable)
nano_time <- function(expr) {
    t0 <- proc.time()[["elapsed"]]
    force(expr)
    t1 <- proc.time()[["elapsed"]]
    t1 - t0
}

# validar método
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

# convierte ei.MD.bayes a formato Pavia
ei.MD.bayes2lphom <- function(ei.object, votes_election1 = NULL, votes_election2 = NULL) {
    X <- ei.object$draws$Cell.counts
    if (length(dim(X)) != 3) stop("Vuelva a calcular la salida de ei.MD.bayes con la opcion ret.mcmc= F")
    X <- apply(X, c(1, 2), mean)
    pjk <- X / rowSums(X)
    pkj <- round(X / colSums(X) * 100, 2)
    pjk <- round(pjk * 100, 2)
    list("VTM" = pjk, "OTM" = pkj)
}

compute_loglik_mult <- function(X, W, P) {
    B <- nrow(X) # ballot boxes
    C <- ncol(X) # candidates
    G <- ncol(W) # groups

    # WP[b,c] = sum_g W[b,g] * P[g,c]
    WP <- W %*% P # (b x g) %*% (g x c) = (b x c)

    # totalWP[b]
    totalWP <- rowSums(WP)

    # pi[b,c]
    pi_bc <- WP / totalWP

    # Totals in each ballot box
    N_b <- rowSums(X)

    # Compute log-likelihood
    ll <- 0
    for (b in 1:B) {
        # multinomial coefficient part
        ll_b <- lgamma(N_b[b] + 1) - sum(lgamma(X[b, ] + 1))

        # sum_{c} X_{bc} * log(pi_{bc})
        ll_b <- ll_b + sum(X[b, ] * log(pi_bc[b, ]))

        ll <- ll + ll_b
    }

    return(ll)
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
#             }
#         }
#         list(mu = mu, Sigma = Sigma)
#     }
#
#     getAverageConditional_R <- function(P_red, mu, Sigma) {
#         cond_mu <- matrix(0, nrow = G, ncol = n)
#         cond_sigma <- vector("list", G)
#         for (g in seq_len(G)) {
#             p_g <- P_red[g, ]
#             cond_mu[g, ] <- mu - p_g
#             cond_sigma[[g]] <- Sigma + tcrossprod(p_g) - diag(p_g, n)
#         }
#         list(mu = cond_mu, Sigma = cond_sigma)
#     }
#
#     M_2_PI <- 2 / pi
#     ll <- 0.0
#
#     for (b in seq_len(B)) {
#         pars <- getParams_R(b, W, P_red)
#         mu <- pars$mu
#         Sigma0 <- pars$Sigma
#
#         cond <- getAverageConditional_R(P_red, mu, Sigma0)
#         mu_cond <- cond$mu
#         sigma_cond <- cond$Sigma
#
#         if (C == 2) {
#             normalizeConstant <- 1.0
#         } else {
#             diag_prod <- prod(diag(sigma_cond[[1]]))
#             det <- 1 / (diag_prod * diag_prod)
#             normalizeConstant <- sqrt((M_2_PI^n) * det)
#         }
#
#         feature <- as.numeric(X[b, ])
#         maha <- matrix(0, nrow = G, ncol = C)
#
#         for (g in seq_len(G)) {
#             mu_g <- mu_cond[g, ]
#             Sigma_g <- sigma_cond[[g]]
#
#             # ---- AQUI EL AJUSTE PARA SINGULARIDAD ----
#             # si Σ_g es singular o casi singular, añadir ridge
#             attempt_inverse <- try(solve(Sigma_g), silent = TRUE)
#
#             if (inherits(attempt_inverse, "try-error")) {
#                 Sigma_g2 <- Sigma_g + ridge * diag(n)
#                 attempt_inverse <- solve(Sigma_g2)
#             }
#             # ------------------------------------------
#
#             Sigma_inv <- attempt_inverse
#
#             x_trunc <- feature[1:n]
#             diff <- x_trunc - mu_g
#             z <- as.numeric(Sigma_inv %*% diff)
#             baseline <- as.numeric(diff %*% z)
#             diag_Sinv <- diag(Sigma_inv)
#
#             maha[g, C] <- baseline
#             for (c in seq_len(n)) {
#                 maha[g, c] <- baseline - 2 * z[c] + diag_Sinv[c]
#             }
#         }
#
#         g0 <- 1
#         logw <- rep(-Inf, C)
#         logw_max <- -Inf
#
#         for (c in seq_len(C)) {
#             prior <- P[g0, c]
#             logP <- if (prior > 0) log(prior) else -Inf
#             logw[c] <- -0.5 * maha[g0, c] + logP
#             if (is.finite(logw[c]) && logw[c] > logw_max) logw_max <- logw[c]
#         }
#
#         den <- sum(exp(logw - logw_max))
#         if (den > 0 && is.finite(logw_max)) {
#             logden <- logw_max + log(den)
#             ll <- ll + logden * log(normalizeConstant)
#         }
#     }
#
#     if (!is.finite(ll)) ll <- 0
#     ll
# }


# ======== Núcleo: ejecuta un método de EI (real o simulado) ========
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
        # if (symmetric) {
        #     em_sym <- function(X, W) {
        #         # 1) Correr EM en esta dirección
        #         output <- run_em(
        #             X = X, W = W, method = method, maxtime = maxtime,
        #             initial_prob = initial_prob, seed = seed_initial,
        #             mcmc_samples = mcmc_samples,
        #             compute_ll = FALSE, # no necesitamos el logLik interno
        #             adjust_prob_cond_method = adjust_prob_cond_method,
        #             verbose = FALSE,
        #             mvncdf_error = 1e-3,
        #             adjust_prob_cond_every = if (!is.null(adjust_prob_cond_every)) {
        #                 adjust_prob_cond_every
        #             } else {
        #                 FALSE
        #             }
        #         )
        #
        #         P <- output$prob # p_{gc} directo
        #         E <- output$expected_outcome # z_{g c b}, dim = (G x C x B)
        #
        #         dims <- dim(E)
        #         G <- dims[1]
        #         C <- dims[2]
        #         B <- dims[3]
        #
        #         # 2) Log-likelihood multinomial "directo" en esta dirección
        #         ll_mult <- compute_loglik_pdf(X, W, P)
        #
        #         # 3) Construir q_induced_rev en "dirección reversa":
        #         #    q_induced_rev[g,c,b] = z_{gcb} / x_{bc}
        #         #    (reemplazando w_{bg} por x_{bc} como denominador)
        #
        #         X_expanded <- array(0, dim = dims) # g x c x b
        #         for (b in seq_len(B)) {
        #             for (c in seq_len(C)) {
        #                 # repetir X[b,c] sobre la dimensión g
        #                 X_expanded[, c, b] <- X[b, c]
        #             }
        #         }
        #
        #         q_induced_rev <- array(0, dim = dims)
        #         mask <- X_expanded > 0
        #         q_induced_rev[mask] <- E[mask] / X_expanded[mask]
        #         # donde X == 0, dejamos q_induced_rev = 0 (no aporta al M-step)
        #
        #         # 4) M-step inducido usando x_{bc} y q_induced_rev:
        #         #    p_induced_rev[g,c] = sum_b x_{bc} q_{b,g,c} / sum_b x_{bc}
        #
        #         q_bgc <- aperm(q_induced_rev, c(3, 1, 2)) # b x g x c
        #
        #         numerator <- matrix(0, nrow = G, ncol = C)
        #         denom <- colSums(X) # sum_b x_{bc}
        #
        #         for (c in seq_len(C)) {
        #             if (denom[c] > 0) {
        #                 numerator[, c] <- t(X[, c]) %*% q_bgc[, , c] # 1xB %*% BxG = 1xG
        #             }
        #         }
        #
        #         p_induced_rev <- matrix(0, nrow = G, ncol = C)
        #
        #         nonzero_cols <- denom > 0
        #         if (any(nonzero_cols)) {
        #             p_induced_rev[, nonzero_cols] <-
        #                 sweep(
        #                     numerator[, nonzero_cols, drop = FALSE],
        #                     2, denom[nonzero_cols], "/"
        #                 )
        #         }
        #         # Para candidatos con denom == 0, podemos reutilizar P (no afectan al LL)
        #         if (any(!nonzero_cols)) {
        #             p_induced_rev[, !nonzero_cols] <- P[, !nonzero_cols, drop = FALSE]
        #         }
        #
        #         # 5) Log-likelihood multinomial inducido
        #         ll_induced_mult <- compute_loglik_mult(X, W, p_induced_rev)
        #
        #         invisible(list(
        #             cond_prob            = output$cond_prob,
        #             q_induced_rev        = q_induced_rev,
        #             prob                 = P,
        #             p_induced_rev        = p_induced_rev,
        #             logLik_mult          = ll_mult,
        #             logLik_induced_mult  = ll_induced_mult,
        #             expected_outcome     = E
        #         ))
        #     }
        #
        #     # Dirección normal y reversa
        #     normal_direction <- em_sym(X = X, W = W)
        #     inv_direction <- em_sym(X = W, W = X)
        #
        #     # tau  y tau_rev con SIEMPRE loglik multinomial
        #     tau <- (normal_direction$logLik_mult - inv_direction$logLik_induced_mult) /
        #         abs(normal_direction$logLik_mult)
        #
        #     tau_rev <- (inv_direction$logLik_mult - normal_direction$logLik_induced_mult) /
        #         abs(inv_direction$logLik_mult)
        #
        #     A1 <- normal_direction$expected_outcome # G x C x B
        #     A2 <- aperm(inv_direction$expected_outcome, c(2, 1, 3)) # G x C x B
        #
        #     if (tau > 0 && tau_rev > 0) {
        #         pond1 <- tau / (tau + tau_rev)
        #         pond2 <- tau_rev / (tau + tau_rev)
        #         Z_est <- pond1 * A1 + pond2 * A2
        #         print("Los ponderadores son:")
        #         print(pond1)
        #         print(pond2)
        #     } else if (tau < 0 && tau_rev > 0) {
        #         Z_est <- A2
        #     } else if (tau > 0 && tau_rev < 0) {
        #         Z_est <- A1
        #     } else {
        #         pond1 <- abs(tau_rev) / (abs(tau) + abs(tau_rev))
        #         pond2 <- abs(tau) / (abs(tau) + abs(tau_rev))
        #         Z_est <- pond1 * A1 + pond2 * A2
        #     }
        #     Z_est_final <- Z_est
        #     output$expected_outcome <- Z_est_final
        # }

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

# ======== EI datasets (paralelo con doParallel) ========
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

    # Prepara bases una sola vez
    base_list <- lapply(bases_arr, function(base_name) {
        df <- merge_small_options(get(base_name), 3, 3) # <- definida en tu código base
        list(name = base_name, df = df)
    })
    base_map <- setNames(
        lapply(base_list, `[[`, "df"),
        sapply(base_list, `[[`, "name")
    )

    # Tareas
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
        # Adjunta paquetes para %dopar%
        if (!"package:foreach" %in% search()) library(foreach)
        if (!"package:doParallel" %in% search()) library(doParallel)

        if (is.null(workers)) workers <- max(1L, parallel::detectCores(logical = TRUE) - 1L)
        cl <- parallel::makeCluster(workers)
        on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
        doParallel::registerDoParallel(cl)
        if (!is.null(reproducible_seed) && requireNamespace("doRNG", quietly = TRUE)) {
            doRNG::registerDoRNG(reproducible_seed)
        }

        # Exporta helpers GLOBALS (definidos arriba en este script)
        parallel::clusterExport(
            cl,
            varlist = c(
                # workers/funciones
                "run_one_ei_task",
                "run_ecological_inference_algorithm",
                # helpers llamados por run_ecological_inference_algorithm
                "validate_method", "compute_EI", "compute_MAE", "nano_time",
                "folder_name_method", "file_name", "ei.MD.bayes2lphom", "ecolRxC",
                # rutas usadas al guardar
                "PATH_SIM_INSTANCES", "PATH_EI_INSTANCES",
                # util
                "sanitize"
            ),
            envir = .GlobalEnv
        )
        # Exporta OBJETOS LOCALES de esta función
        parallel::clusterExport(
            cl,
            varlist = c("base_map", "invert", "maxtime", "save_json"),
            envir = environment()
        )
        # Carga paquetes en workers
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

# ======== Simuladas (versión secuencial; puedes paralelizar igual si quieres) ========
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
