# Generate the data for simulated instances / EI datasets.
source("src/R_functions.R")

# Runs the simulated instances grid and stores estimated parameters and timing.
fun_run_simulated_instances <- function() {
    run_simulated_instances(
        # method_arr = c("mcmc","mvn_cdf","mvn_pdf","mult", "exact", "nslphom_dual_w", "lphom", "lclphom"),
        method_arr = c(
            # "mvn_cdf"
            # "mvn_pdf",
            # "mult"
            # "ei.MD.bayes"
            "exact"
            # "nslphom_dual_w",
            # "lphom",
            # "lphom_joint",
            # "lclphom",
            # # "nslphom_dual_w",
            # "nslphom_dual"
        ), # todos sin exact
        mcmc_samples_arr = c(100, 1000),
        I_arr = c(100),
        B_arr = c(50),
        C_arr = c(4),
        # C_arr = c(2, 3),
        G_arr = c(8),
        # G_arr = c(3, 4, 6, 8),
        adjust_prob_cond_arr = c("project_lp"),
        adjust_prob_cond_every_arr = c(FALSE),
        lambda_arr = c(0.5),
        # seed_arr = 1:20, # seed for the generated instance
        seed_arr = 9:20,
        seed_initial_arr = c(NULL), # seed for the initial p of the EM
        maxtime = 3600,
        skip_cases = FALSE
    )
}

# Convenience runner for a smaller simulated grid when testing.
test_fun_run_simulated_instances <- function() {
    run_simulated_instances(
        # method_arr = c("mcmc","mvn_cdf","mvn_pdf","mult", "exact", "nslphom_dual_w", "lphom", "lclphom"),
        # method_arr = c("mcmc","mvn_cdf","mvn_pdf","mult", "nslphom_dual_w", "lphom", "lclphom"),
        adjust_prob_cond_arr = c("project_lp", "lp", ""),
        adjust_prob_cond_every_arr = c(TRUE, FALSE),
        method_arr = c("mult", "lphom", "lclphom", "nslphom_dual_w"),
        mcmc_samples_arr = c(100, 1000),
        I_arr = c(100),
        B_arr = c(50),
        C_arr = c(2),
        G_arr = c(3),
        lambda_arr = c(0.5),
        seed_arr = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), # seed for the generated instance
        seed_initial_arr = c(NULL), # seed for the initial p of the EM
        maxtime = 3600,
        skip_cases = TRUE
    )
}


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
        ), # todos sin exact
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

# Narrow test runner for EI datasets.
test_fun_run_eidatasets_instances <- function() {
    run_eidatasets_instances(
        method_arr = c( # "mcmc",
            "lphom_joint"
            # "mvn_cdf"
            # "mvn_pdf"
            # "mult",
            # "nslphom_dual_w",
            # "lphom",
            # "lclphom"
        ), # todos sin exact
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
        seed_initial_arr = c(123), # same convention: NULL -> -123 ("group_proportional")
        maxtime = 3600
    )
}

# Default run executed when sourcing this file
fun_run_eidatasets_instances(symmetric = TRUE)
