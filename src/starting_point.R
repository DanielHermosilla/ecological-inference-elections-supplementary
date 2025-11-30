# Explore EM sensitivity by computing Frobenius distances across random starts.
library(fastei)

SEEDS <- c(1:20)
CANDIDATES <- c(2:4)
GROUPS <- c(2:4)

df <- data.frame()
for (c in CANDIDATES) {
    for (g in GROUPS) {
        problem <- simulate_election(num_groups = g, num_candidates = c, num_ballots = 50, seed = 42)
        P_list <- list()

        for (i in seq_along(SEEDS)) {
            solution <- run_em(problem, initial_prob = "random", seed = i)
            P_list[[i]] <- solution$prob
        }
        # Obtain the mean, sd and max frobenius norm for combination c, g. Bind it to the dataframe
        # compute pairwise Frobenius norms between all EM solutions
        pairs <- combn(length(P_list), 2)
        frob_values <- apply(pairs, 2, function(idx) {
            D <- P_list[[idx[1]]] - P_list[[idx[2]]]
            sqrt(sum(D^2)) # Frobenius norm
        })

        mean_frob <- mean(frob_values)
        sd_frob <- sd(frob_values)
        max_frob <- max(frob_values)

        # Bind results to dataframe
        df <- rbind(df, data.frame(
            candidates = c,
            groups = g,
            mean_frobenius = mean_frob,
            sd_frobenius = sd_frob,
            max_frobenius = max_frob
        ))
    }
}
