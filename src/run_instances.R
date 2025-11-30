# Generate the data for simulated instances / EI datasets.
source("src/R_functions.R")

# Runs EI datasets (real data) with an optional symmetry flag.
fun_run_eidatasets_instances <- function(symmetric = FALSE) {
    run_eidatasets_instances(
        method_arr = c(
            # "mcmc",
            "mvn_cdf"
            # "nslphom_dual_a",
            # "ei.MD.bayes"
            # "nslphom_joint"
            # "ecolRxC"
            # "mvn_pdf"
            # "mult"
            # "nslphom_dual_w",
            # "lphom",
            # "lclphom"
        ),
        adjust_prob_cond_arr = c("project_lp"),
        adjust_prob_cond_every_arr = c(FALSE),
        mcmc_samples_arr = c(1000),
        bases_arr = c(
            "ei_NZ_2002",
            "ei_NZ_2005",
            "ei_SCO_2007",
            "ei_NZ_2008",
            "ei_NZ_2011",
            "ei_NZ_2014",
            "ei_NZ_2017",
            "ei_NZ_2020"
        ),
        workers = 12,
        parallel = FALSE,
        seed_initial_arr = c(-123), # same convention: NULL -> -123 ("group_proportional")
        maxtime = 3600,
        invert = FALSE, # whether to invert the X and W matrices
        symmetric = symmetric,
        save_json = TRUE
    )
}

# Default run executed when sourcing this file
fun_run_eidatasets_instances(symmetric = TRUE)
