simulate_multi <- function(
    N = 100, # Number of sites,
    J = 1, # Number of visits per site,
    n_obs = 5, # Number of observers,
    beta0 = -2, # Intercept for eco proc
    beta1 = 0.3, # Beta for salinity
    beta2 = 0.1, # beta for temperature
    alpha0 = 1.5, # baseline detection intercept
    sigma_obs = 0.5, # SD of random observer effects
    observer_eff = NULL # random effect for each observer
) {
    # Site-level covariates
    salinity <- runif(N, 0, 10)
    temperature <- runif(N, 15, 30)

    # --- ECOLOGICAL PROCESS: occupancy ---
    psi <- plogis(beta0 + beta1 * salinity + beta2 * temperature)
    Z <- rbinom(N, 1, psi)

    # --- OBSERVATION PROCESS: detection with observer random effect ---
    if (is.null(observer_eff)) {
        observer_eff <- rnorm(n_obs, 0, sigma_obs)
    } else if (length(observer_eff) != n_obs) {
        stop("observer_eff should be equal to n_obs")
    }

    # Assign observers to visits
    observer_id <- sample(seq_len(n_obs), N * J, replace = TRUE)

    # Calculate detection probability for each observation
    p_ij <- plogis(alpha0 + observer_eff[observer_id])

    # Simulate detection histories
    Y <- matrix(NA, nrow = N, ncol = J)
    for (i in 1:N) {
        for (j in 1:J) {
            idx <- (i - 1) * J + j
            Y[i, j] <- rbinom(1, 1, Z[i] * p_ij[idx])
        }
    }

    # --- Combine into data frame ---
    data <- data.frame(
        site = rep(seq_len(N), each = J),
        visit = rep(seq_len(J), N),
        observer = factor(observer_id),
        salinity = rep(salinity, each = J),
        temperature = rep(temperature, each = J),
        Z = rep(Z, each = J),
        Y = c(t(Y)),#as.vector(Y),
        observer_eff = observer_eff[observer_id]
    )
}
