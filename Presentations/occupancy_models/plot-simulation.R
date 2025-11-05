set.seed(2025) 
source("simulate-multi.R")

N <- 500 # Number of sites
n_obs <- 5 # Number of observers

data <- simulate_multi(
    N = N, # Number of sites,
    J = 3, # Number of visits per site,
    n_obs = 5, # Number of observers,
    beta0 = -2, # Intercept for eco proc
    beta1 = 0.3, # Beta for salinity
    beta2 = 0.1, # beta for temperature
    alpha0 = 0.6, # baseline detection intercept
    sigma_obs = 0.5, # SD of random observer effects
    observer_eff = NULL # random effect for each observer
)

site_z <- aggregate(data$Z, list(data$site), max)
site_y <- aggregate(data$Y, list(data$site), max)
sum(site_z[,2])
sum(site_y[,2])
total_occupancy <- sum(site_z$x)
total_visits <- N * 1 # For plotting, consider multiple visits as a single site

plot_points <- TRUE
skip_large <- TRUE # To avoid many steps when simulation number is large
plot_speed <- 0.1

if (plot_points) {
    plot(NULL, xlim = c(0, total_visits), ylim = c(0, total_occupancy+5),
     xlab = "Site visited", ylab = c("Sites occupied"))
    observer_names <- c("Silas", "Lotte", "Damian", "Johannes", "Clyde") #paste("Observer", 1:5)
    observer_color <- RColorBrewer::brewer.pal(5, "Dark2")
    legend("topleft", 
        legend = c("Occupancy", "Observed", observer_names), 
        col = c("black", "grey", observer_color), 
        lty = c(1, 2, rep(1, 5)),
        bty = "n", 
        pt.cex = 2, 
        cex = 1.2, 
        text.col = "black", 
        horiz = F
    )
}

user_obs <- vector("list", 5)

for (i in seq_len(N)) {
    sel_data <- data[data$site == i, ]
    if (i == 1) {
        total_observation <- data.frame(
                site = i,
                Z = max(sel_data$Z),
                Y = max(sel_data$Y)
            )
    } else {
        total_observation <- rbind(
            total_observation[,1:3],
            data.frame(
                site = i,
                Z = max(sel_data$Z),
                Y = max(sel_data$Y)
            )
        )
    }
    total_observation$cum_Z <- cumsum(
        total_observation$Z
    )
    total_observation$cum_Y <- cumsum(
        total_observation$Y
    )

    for (o in unique(sel_data$observer)) {
        sel_data_o <- sel_data[sel_data$observer == o,]
        o <- as.integer(o)
        user_obs[[o]] <- rbind(
        user_obs[[o]][,1:4],
            data.frame(
                site = i, observer = sel_data_o$observer[1],
                Z = max(sel_data_o$Z), Y = max(sel_data_o$Y)
            )
        )
        user_obs[[o]]$cum_Z <- cumsum(
            user_obs[[o]]$Z
        )
        user_obs[[o]]$cum_Y <- cumsum(
            user_obs[[o]]$Y
        )
    }
    
   if (plot_points) {
        if (N >= 200) {
            if (!i %in% seq(0, N, by = 10)) next
        }
        points(x = total_observation$site, y = total_observation$cum_Z, type ="l", col = "black")
        points(x = total_observation$site, y = total_observation$cum_Y, lty = 2, type ="l", col = "grey")

        for (k in seq_along(user_obs)) {
            points(x = user_obs[[k]]$site, y = user_obs[[k]]$cum_Z, type ="l", col = observer_color[k])
            points(x = user_obs[[k]]$site, y = user_obs[[k]]$cum_Y, lty = 2, type ="l", col = observer_color[k])
        }
        Sys.sleep(plot_speed)
    }
}


# Model:
library(spOccupancy)

occ_formula <- ~ salinity + temperature
det_formula <- ~ (1 | observer)

n.samples <- 15000
n.report <- 1000

# Create data object
N <- length(unique(data$site))
J <- length(unique(data$visit))
y_mat <- matrix(data$Y, nrow = N, ncol = J, byrow = TRUE)
salinity_vec <- as.numeric(tapply(data$salinity, data$site, function(x) x[1]))
temp_vec      <- as.numeric(tapply(data$temperature, data$site, function(x) x[1]))
occ_covs_mat  <- cbind(salinity = salinity_vec, temperature = temp_vec)
observer_int <- as.integer(factor(data$observer))
observer_mat <- matrix(observer_int, nrow = N, ncol = J, byrow = TRUE)

stopifnot(all(dim(y_mat) == c(N, J)))
stopifnot(all(dim(occ_covs_mat)[1] == N))
stopifnot(all(dim(observer_mat) == c(N, J)))

# Build the data list for PGOcc
data_obj <- list(
  y = y_mat,
  occ.covs = occ_covs_mat,
  det.covs = list(observer = observer_mat)
)
str(data_obj)
# y should be 500 rows (1 each site) x 3 columns (1 each visit)
# occ.covs should be 500 rows (1 each site) for 2 variables
# det.covs should be 500 rows (1 each site) x 3 columns (1 each visit)
inits <- list(alpha = 0, beta = 0, sigma.sq.psi = 0.5, 
              sigma.sq.p = 0.5, z = apply(data_obj$y, 1, max, na.rm = TRUE), 
              sigma.sq = 1, phi = 3 / 0.5)

out <- PGOcc(
  occ.formula = occ_formula,
  det.formula = det_formula,
  data = data_obj,
  n.samples = n.samples,
  n.burn = 4000,
  n.thin = 1,
  n.chains = 4,
  n.omp.threads = 1,
  verbose = TRUE,
  n.report = n.report,
  inits = inits
)

summary(out)
#reference values
# beta0 = -2, # Intercept for eco proc
# beta1 = 0.3, # Beta for salinity
# beta2 = 0.1, # beta for temperature
# alpha0 = 0.6, # baseline detection intercept
# sigma_obs = 0.5, # SD of random observer effects

occ <- predict(out, as.matrix(cbind(intercept = 1, data_obj$occ.covs)), type = "occupancy")
occ_quants <- apply(occ$psi.0.samples, 2, quantile, c(0.025, 0.5, 0.975))
occ_med <- occ_quants[2,]

plot(plogis(-2 + 0.3 * data_obj$occ.covs[,1] + 0.1 * data_obj$occ.covs[,2]),
     occ_med, pch = 21, bg = "lightblue",
     xlab = "Simulated", ylab = "Estimated")

det <- predict(out, as.matrix(cbind(intercept = 1, observer = as.integer(data$observer))), type = "detection")
det_quants <- apply(det$p.0.samples, 2, quantile, c(0.025, 0.5, 0.975))
det_med <- det_quants[2,]

true_p <- plogis(0.6 + data$observer_eff)
plot(true_p, det_med, pch = 21, bg = "orange",
     xlab = "True detection probability", ylab = "Estimated detection probability")
