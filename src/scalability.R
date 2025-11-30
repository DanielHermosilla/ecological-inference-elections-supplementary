# Benchmarks runtime feasibility over large ballot/voter combinations.
# ---- Paquetes ----
library(fastei)
library(lphom)
library(ecolRxC)
library(microbenchmark)
library(eiPack)
library(dplyr)
library(R.utils) # for withTimeout
library(future.apply) # for future_lapply

args <- commandArgs(trailingOnly = TRUE)
run_parallel <- any(args %in% c("--parallel", "--par", "-p"))

if (run_parallel) {
    future::plan(future::multisession, workers = max(1, parallel::detectCores() - 1))
    message("Running in parallel (omit --parallel to stay sequential).")
} else {
    future::plan(future::sequential)
    message("Running sequentially; enable workers with --parallel.")
}

# ballot_sizes <- c(1000, 5000, 10000)
ballot_sizes <- c(50000)
voters_per_ballot <- c(100, 1000, 10000)
# voters_per_ballot <- c(10000)
TIME_LIMIT <- 3600 * 10

# Combination grids
grid <- expand.grid(
    ballots = ballot_sizes,
    voters = voters_per_ballot,
    stringsAsFactors = FALSE
)

# Wrap expression with timeout handling.
safe_timeout <- function(expr, limit) {
    tryCatch(
        R.utils::withTimeout(expr, timeout = limit, onTimeout = "silent"),
        TimeoutException = function(e) {
            message("⏱️ Timeout: ", conditionMessage(e))
            NULL
        },
        error = function(e) {
            message("⛔ Error dentro de safe_timeout: ", conditionMessage(e))
            NULL
        }
    )
}

# Catch and log non-timeout errors without stopping the benchmark.
safe_run <- function(expr) {
    tryCatch(expr, error = function(e) {
        message("Error interno: ", conditionMessage(e))
        NULL
    })
}

# Runs one combinations of (ballots, voters), so it can be parallelized
# Runs one (ballots, voters) combination; intended for parallel map.
run_one <- function(ballots, voters, TIME_LIMIT) {
    suppressPackageStartupMessages({
        library(fastei)
        library(lphom)
        library(ecolRxC)
        library(eiPack)
        library(dplyr)
        library(R.utils)
    })

    instance <- simulate_election(
        num_groups = 6,
        num_candidates = 6,
        num_ballots = ballots,
        ballot_voters = rep(voters, ballots)
    )

    # fastei offers a timeout
    pdf <- run_em(instance, method = "mvn_pdf", maxtime = TIME_LIMIT)
    cdf <- run_em(instance, method = "mvn_cdf", maxtime = TIME_LIMIT)
    mult <- run_em(instance, method = "mult", maxtime = TIME_LIMIT)
    state_pdf <- !is.null(pdf$prob)
    state_cdf <- !is.null(cdf$prob)
    state_mult <- !is.null(mult$prob)

    # bayes
    colnames(instance$X) <- paste0("cand_", seq_len(ncol(instance$X)))
    colnames(instance$W) <- paste0("grp_", seq_len(ncol(instance$W)))
    var.y <- paste(paste0("`", colnames(instance$X), "`"), collapse = ",")
    var.x <- paste(paste0("`", colnames(instance$W), "`"), collapse = ",")
    formula_str <- paste("cbind(", var.y, ") ~ cbind(", var.x, ")")
    X_df <- as.data.frame(instance$X, stringsAsFactors = FALSE)
    W_df <- as.data.frame(instance$W, stringsAsFactors = FALSE)
    base <- data.frame(W_df, X_df, check.names = FALSE)

    # bayes
    bayes <- safe_timeout(
        ei.MD.bayes(
            formula = as.formula(formula_str),
            data = base,
            sample = 1000,
            thin = 100,
            burnin = 100000,
            ret.mcmc = FALSE
        ),
        limit = TIME_LIMIT
    )
    state_bayes <- !is.null(bayes) && !is.null(bayes$draws$Cell.counts)


    # nslphom
    nslphom <- safe_timeout(
        safe_run(
            nslphom_dual(instance$W, instance$X)
        ),
        limit = TIME_LIMIT
    )
    state_nslphom <- !is.null(nslphom) && !is.null(nslphom$VTM.votes.a)

    # ecolRxC
    ecol <- safe_timeout(safe_run(ecolRxC(instance$W, instance$X)), limit = TIME_LIMIT)
    state_ecol <- !is.null(ecol) && !is.null(ecol$VTM.votes)

    # joint
    joint <- safe_timeout(safe_run(nslphom_joint(instance$W, instance$X)), limit = TIME_LIMIT)
    state_joint <- !is.null(joint) && !is.null(joint$VTM.votes)

    # lphom
    lphom_ <- safe_timeout(safe_run(lphom(instance$W, instance$X)), limit = TIME_LIMIT)
    state_lphom <- !is.null(lphom_) && !is.null(lphom_$VTM.votes)

    # lclphom
    lclphom_ <- safe_timeout(safe_run(lclphom(instance$W, instance$X)), limit = TIME_LIMIT)
    state_lclphom <- !is.null(lclphom_) && !is.null(lclphom_$VTM.votes)


    tib <- tibble(
        ballots = ballots,
        voters = voters,
        mvn_pdf = state_pdf,
        mvn_cdf = state_cdf,
        mult = state_mult,
        bayes = state_bayes,
        nslphom = state_nslphom,
        ecolRxC = state_ecol,
        nslphom_joint = state_joint,
        lphom = state_lphom,
        lclphom = state_lclphom
    )

    # print when done
    message(sprintf("✅ Finished: ballots=%d, voters=%d", ballots, voters))
    print(tib)

    return(tib)
}

results_list <- future_lapply(
    seq_len(nrow(grid)),
    function(i) run_one(grid$ballots[i], grid$voters[i], TIME_LIMIT),
    future.seed = TRUE # reproducibilidad de RNG entre workers
)

results <- bind_rows(results_list)

write.csv(results, "output/scalability_results50000.csv", row.names = FALSE)
print(results, n = Inf)
