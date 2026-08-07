# ============================================================================
# Spatial Analysis of Lithology Data - Idaho National Engineering Laboratory
# Using GAM (Generalized Additive Model) for Smooth Surface Interpolation
# WITH 3D PERSPECTIVE PLOTS
# ============================================================================

# 1. SETUP AND PACKAGES ------------------------------------------------

# Set working directory (modify this path as needed)
setwd("D:/queensland health/Hierachical Modeling and Analysis for spatial data/chapter_6/exercise")
getwd()

# Load required packages
library(tidyverse)
library(mgcv)         # For GAM smoothing
library(viridisLite)  # For color palettes
library(fields)       # For image.plot with legend

# 2. IMPORT AND CLEAN DATA ---------------------------------------------

# Import the lithology dataset
litho <- read.table(
  "lithology.dat",
  header = FALSE,
  skip = 1,
  col.names = c(
    "Well_Name",
    "Northing",
    "Easting",
    "Surf_Elevation",
    "Thickness_AB",
    "Elevation_AB",
    "Thickness_BC",
    "Elevation_BC",
    "Thickness_CD",
    "Elevation_CD",
    "Thickness_DE"
  ),
  na.strings = c("NA", "N/A", "", " "),
  fill = TRUE,
  stringsAsFactors = FALSE
)

# Select only the variables needed and clean the data
litho2 <- litho |>
  select(
    Northing,
    Easting,
    Surf_Elevation,
    Thickness = Thickness_AB,
    AB_Elevation = Elevation_AB
  ) |>
  mutate(
    across(everything(), ~ as.numeric(gsub(",", "", as.character(.x))))
  ) |>
  drop_na()

# Check the cleaned data
cat("\n=== DATA SUMMARY ===\n")
cat("Number of complete records:", nrow(litho2), "\n\n")
head(litho2)
summary(litho2)

# 3. GAM-BASED PLOTTING FUNCTION ----------------------------------------

plot_lithology_gam <- function(data, variable, plot_title, 
                               show_points = TRUE, 
                               n_contours = 10,
                               grid_size = 150,
                               color_palette = viridis,
                               k_smooth = 20) {
  
  # Extract the variable to model
  z <- data[[variable]]
  
  # Fit GAM with thin plate regression splines
  # k controls the maximum degrees of freedom (smoothness)
  formula <- as.formula(paste(variable, "~ s(Northing, Easting, k =", k_smooth, ")"))
  gam_model <- gam(formula, data = data)
  
  # Create prediction grid
  northing_range <- seq(min(data$Northing), max(data$Northing), length.out = grid_size)
  easting_range <- seq(min(data$Easting), max(data$Easting), length.out = grid_size)
  
  # Expand grid for prediction
  pred_grid <- expand.grid(
    Northing = northing_range,
    Easting = easting_range
  )
  
  # Predict on the grid
  pred_grid$value <- predict(gam_model, newdata = pred_grid, type = "response")
  
  # Convert to matrix format for image/contour functions
  # Note: image() expects x = columns, y = rows, so we need to reshape
  value_matrix <- matrix(pred_grid$value, 
                         nrow = length(easting_range), 
                         ncol = length(northing_range),
                         byrow = FALSE)
  
  # Calculate proper aspect ratio
  asp_ratio <- diff(range(data$Northing)) / diff(range(data$Easting))
  
  # Determine z range for consistent color scaling
  zlim <- range(pred_grid$value, na.rm = TRUE)
  
  # Create the image plot with legend using fields package
  if (requireNamespace("fields", quietly = TRUE)) {
    # Plot with color legend
    fields::image.plot(
      x = easting_range,
      y = northing_range,
      z = value_matrix,
      col = color_palette(100),
      xlab = "Easting (ft)",
      ylab = "Northing (ft)",
      main = plot_title,
      asp = 1/asp_ratio,
      zlim = zlim,
      legend.width = 1.2,
      legend.mar = 3.5,
      legend.args = list(text = variable, side = 4, line = 2.5, cex = 0.8)
    )
  } else {
    # Fallback: plot without legend
    image(
      x = easting_range,
      y = northing_range,
      z = value_matrix,
      col = color_palette(100),
      xlab = "Easting (ft)",
      ylab = "Northing (ft)",
      main = plot_title,
      asp = 1/asp_ratio,
      zlim = zlim
    )
  }
  
  # Add contour lines
  contour(
    x = easting_range,
    y = northing_range,
    z = value_matrix,
    add = TRUE,
    drawlabels = TRUE,
    col = "black",
    lwd = 0.8,
    nlevels = n_contours,
    labcex = 0.6
  )
  
  # Add observed borehole locations
  if (show_points) {
    points(
      data$Easting,
      data$Northing,
      pch = 16,
      cex = 0.4,
      col = adjustcolor("white", alpha.f = 0.6)
    )
  }
  
  box()
  
  # Return the model and predictions for further analysis
  invisible(list(
    model = gam_model,
    grid = pred_grid,
    matrix = value_matrix,
    x = easting_range,
    y = northing_range,
    zlim = zlim
  ))
}

# 4. FUNCTION FOR 3D PERSPECTIVE PLOT ----------------------------------

plot_3d_perspective <- function(data, variable, plot_title, 
                                theta = 30, phi = 30, 
                                color_palette = viridis,
                                k_smooth = 20,
                                grid_size = 50) {
  
  # Fit GAM
  formula <- as.formula(paste(variable, "~ s(Northing, Easting, k =", k_smooth, ")"))
  gam_model <- gam(formula, data = data)
  
  # Create prediction grid
  northing_range <- seq(min(data$Northing), max(data$Northing), length.out = grid_size)
  easting_range <- seq(min(data$Easting), max(data$Easting), length.out = grid_size)
  
  # Expand grid for prediction
  pred_grid <- expand.grid(
    Northing = northing_range,
    Easting = easting_range
  )
  
  # Predict on the grid
  pred_grid$value <- predict(gam_model, newdata = pred_grid, type = "response")
  
  # Convert to matrix
  value_matrix <- matrix(pred_grid$value, 
                         nrow = length(easting_range), 
                         ncol = length(northing_range),
                         byrow = FALSE)
  
  # Create color palette for the surface
  n_colors <- 50
  colors <- color_palette(n_colors)
  
  # Create perspective plot
  persp(
    x = easting_range,
    y = northing_range,
    z = value_matrix,
    theta = theta,
    phi = phi,
    col = colors[cut(value_matrix, n_colors)],
    shade = 0.3,
    border = NA,
    xlab = "Easting (ft)",
    ylab = "Northing (ft)",
    zlab = variable,
    main = plot_title,
    ticktype = "detailed",
    cex.axis = 0.8,
    cex.lab = 0.9,
    expand = 0.8
  )
  
  # Return the matrix for potential further use
  invisible(list(
    matrix = value_matrix,
    x = easting_range,
    y = northing_range,
    model = gam_model
  ))
}

# 5. GENERATE 2D PLOTS USING GAM -------------------------------------

cat("\n=== GENERATING 2D IMAGE PLOTS ===\n")

# Set up plotting area for three side-by-side plots
par(
  mfrow = c(1, 3),
  mar = c(4.5, 4.5, 3.5, 4.5),  # Increased right margin for legend
  oma = c(1, 1, 3, 1)           # Outer margins for overall title
)

# Plot 1: A-B Thickness
cat("\nFitting GAM for Thickness...\n")
result_thickness <- plot_lithology_gam(
  data = litho2,
  variable = "Thickness",
  plot_title = "A-B Thickness",
  n_contours = 10,
  grid_size = 150,
  color_palette = viridis,
  k_smooth = 20
)

# Plot 2: Surface Elevation
cat("\nFitting GAM for Surface Elevation...\n")
result_surface <- plot_lithology_gam(
  data = litho2,
  variable = "Surf_Elevation",
  plot_title = "Surface Elevation",
  n_contours = 10,
  grid_size = 150,
  color_palette = viridis,
  k_smooth = 20
)

# Plot 3: A-B Elevation
cat("\nFitting GAM for A-B Elevation...\n")
result_ab <- plot_lithology_gam(
  data = litho2,
  variable = "AB_Elevation",
  plot_title = "A-B Elevation",
  n_contours = 10,
  grid_size = 150,
  color_palette = viridis,
  k_smooth = 20
)

# Add overall title for 2D plots
mtext(
  "Spatial Distribution of Subsurface Variables (GAM Smoothing)",
  outer = TRUE,
  cex = 1.4,
  font = 2,
  line = 1
)

# Add subtitle with number of observations
mtext(
  paste("Based on", nrow(litho2), "borehole locations | Smoothing parameter k = 20"),
  outer = TRUE,
  cex = 0.9,
  font = 3,
  line = 0
)

# Reset plotting parameters
par(mfrow = c(1, 1))

# 6. GENERATE 3D PERSPECTIVE PLOTS -------------------------------------

cat("\n=== GENERATING 3D PERSPECTIVE PLOTS ===\n")

# Option 1: Three 3D plots in a row
par(
  mfrow = c(1, 3),
  mar = c(1, 1, 3, 1),
  oma = c(2, 2, 3, 2)
)

# 3D Plot 1: A-B Thickness
plot_3d_perspective(
  data = litho2,
  variable = "Thickness",
  plot_title = "A-B Thickness (3D)",
  theta = 30,
  phi = 30,
  color_palette = viridis,
  k_smooth = 20,
  grid_size = 50
)

# 3D Plot 2: Surface Elevation
plot_3d_perspective(
  data = litho2,
  variable = "Surf_Elevation",
  plot_title = "Surface Elevation (3D)",
  theta = 30,
  phi = 30,
  color_palette = viridis,
  k_smooth = 20,
  grid_size = 50
)

# 3D Plot 3: A-B Elevation
plot_3d_perspective(
  data = litho2,
  variable = "AB_Elevation",
  plot_title = "A-B Elevation (3D)",
  theta = 30,
  phi = 30,
  color_palette = viridis,
  k_smooth = 20,
  grid_size = 50
)

# Add overall title for 3D plots
mtext(
  "3D Perspective Views of Subsurface Variables",
  outer = TRUE,
  cex = 1.4,
  font = 2,
  line = 1
)

# Reset plotting parameters
par(mfrow = c(1, 1))



# ============================================================
# Exercise 5(b): Bayesian univariate Gaussian spatial models
# Exponential and Matern covariance functions
# ============================================================

library(tidyverse)
library(coda)

set.seed(2026)

# ============================================================
# 1. Prepare the data
# ============================================================

# This assumes litho2 was created in part (a) and contains:
# Northing, Easting, Thickness, Surf_Elevation, AB_Elevation

spatial_data <- litho2 |>
  filter(
    is.finite(Northing),
    is.finite(Easting),
    is.finite(Thickness),
    is.finite(Surf_Elevation),
    is.finite(AB_Elevation),
    Thickness > 0
  ) |>
  mutate(
    log_Thickness = log(Thickness),
    
    # Convert coordinates to units of 1,000 original coordinate units
    Easting_scaled = (Easting - mean(Easting)) / 1000,
    Northing_scaled = (Northing - mean(Northing)) / 1000,
    
    # Standardise the covariates for numerical stability
    Surf_Elevation_scaled =
      as.numeric(scale(Surf_Elevation)),
    
    AB_Elevation_scaled =
      as.numeric(scale(AB_Elevation))
  )

nrow(spatial_data)
summary(spatial_data)
head(spatial_data)

# Response
y <- spatial_data$log_Thickness

# Design matrix
X <- model.matrix(
  ~ Surf_Elevation_scaled + AB_Elevation_scaled,
  data = spatial_data
)

# Spatial coordinates
coords <- as.matrix(
  spatial_data[, c("Easting_scaled", "Northing_scaled")]
)

n <- length(y)
p <- ncol(X)

# Pairwise distance matrix
distance_matrix <- as.matrix(dist(coords))

# ============================================================
# 2. Covariance functions
# ============================================================

spatial_correlation <- function(
    distance_matrix,
    phi,
    covariance = c("exponential", "matern"),
    nu = 1.5
) {
  
  covariance <- match.arg(covariance)
  
  if (!is.finite(phi) || phi <= 0) {
    stop("phi must be positive.")
  }
  
  if (covariance == "exponential") {
    
    # Exponential correlation:
    # R[j,k] = exp(-d[j,k] / phi)
    
    R <- exp(-distance_matrix / phi)
    
  } else {
    
    # Matern correlation:
    #
    # R(d) =
    # 2^(1-nu) / Gamma(nu) *
    # (sqrt(2*nu)*d/phi)^nu *
    # K_nu(sqrt(2*nu)*d/phi)
    
    x <- sqrt(2 * nu) * distance_matrix / phi
    
    R <- matrix(0, nrow(x), ncol(x))
    
    positive <- x > 0
    
    R[positive] <-
      (2^(1 - nu) / gamma(nu)) *
      x[positive]^nu *
      besselK(x[positive], nu = nu)
    
    diag(R) <- 1
  }
  
  # Small numerical stabilisation term
  diag(R) <- diag(R) + 1e-8
  
  R
}

# Check the functions
R_exp_test <- spatial_correlation(
  distance_matrix,
  phi = 1,
  covariance = "exponential"
)

R_mat_test <- spatial_correlation(
  distance_matrix,
  phi = 1,
  covariance = "matern",
  nu = 1.5
)

range(R_exp_test)
range(R_mat_test)

# ============================================================
# 3. Inverse-Gamma random-number generator
# ============================================================

# Parameterisation:
#
# p(x) proportional to
# x^(-a-1) exp(-b/x)

rinvgamma <- function(n, shape, scale) {
  1 / rgamma(
    n,
    shape = shape,
    rate = scale
  )
}

# ============================================================
# 4. Log-target for the range parameter phi
# ============================================================

log_phi_target <- function(
    log_phi,
    w,
    sigma2,
    distance_matrix,
    covariance,
    nu,
    phi_shape,
    phi_rate
) {
  
  phi <- exp(log_phi)
  
  R <- tryCatch(
    spatial_correlation(
      distance_matrix = distance_matrix,
      phi = phi,
      covariance = covariance,
      nu = nu
    ),
    error = function(e) NULL
  )
  
  if (is.null(R)) {
    return(-Inf)
  }
  
  chol_R <- tryCatch(
    chol(R),
    error = function(e) NULL
  )
  
  if (is.null(chol_R)) {
    return(-Inf)
  }
  
  # log determinant of R
  log_det_R <- 2 * sum(log(diag(chol_R)))
  
  # w' R^(-1) w
  R_inv_w <- backsolve(
    chol_R,
    forwardsolve(t(chol_R), w)
  )
  
  quadratic <- drop(crossprod(w, R_inv_w))
  
  # Log density of w | sigma2, phi
  log_w_density <-
    -0.5 * log_det_R -
    0.5 * quadratic / sigma2
  
  # Gamma prior for phi
  log_phi_prior <- dgamma(
    phi,
    shape = phi_shape,
    rate = phi_rate,
    log = TRUE
  )
  
  # Jacobian because the MH proposal is on log(phi)
  log_jacobian <- log_phi
  
  log_w_density +
    log_phi_prior +
    log_jacobian
}

# ============================================================
# 5. Spatial MCMC function
# ============================================================

run_spatial_mcmc <- function(
    y,
    X,
    distance_matrix,
    covariance = c("exponential", "matern"),
    nu = 1.5,
    n_iter = 30000,
    burn_in = 10000,
    thin = 10,
    prior = list(
      sigma_shape = 0.1,
      sigma_scale = 0.1,
      tau_shape = 0.1,
      tau_scale = 0.1,
      phi_shape = 0.1,
      phi_rate = 0.1
    ),
    proposal_sd = 0.12,
    seed = 2026
) {
  
  set.seed(seed)
  
  covariance <- match.arg(covariance)
  
  n <- length(y)
  p <- ncol(X)
  
  if (burn_in >= n_iter) {
    stop("burn_in must be smaller than n_iter.")
  }
  
  # Number of posterior samples to retain
  keep_iterations <- seq(
    burn_in + thin,
    n_iter,
    by = thin
  )
  
  n_keep <- length(keep_iterations)
  
  # ----------------------------------------------------------
  # Initial values
  # ----------------------------------------------------------
  
  beta <- drop(solve(crossprod(X), crossprod(X, y)))
  
  residuals_initial <- y - drop(X %*% beta)
  
  total_variance <- var(residuals_initial)
  
  sigma2 <- total_variance / 2
  tau2 <- total_variance / 2
  
  # Initial spatial range
  nonzero_distances <- distance_matrix[distance_matrix > 0]
  phi <- median(nonzero_distances)
  
  # Initial spatial random effects
  w <- rep(0, n)
  
  # ----------------------------------------------------------
  # Storage
  # ----------------------------------------------------------
  
  samples <- matrix(
    NA_real_,
    nrow = n_keep,
    ncol = p + 3
  )
  
  colnames(samples) <- c(
    colnames(X),
    "sigma2",
    "tau2",
    "phi"
  )
  
  w_sum <- rep(0, n)
  
  accepted_phi <- 0
  stored <- 0
  
  XtX <- crossprod(X)
  XtX_inv <- solve(XtX)
  
  # ----------------------------------------------------------
  # MCMC iterations
  # ----------------------------------------------------------
  
  for (iteration in seq_len(n_iter)) {
    
    # ========================================================
    # Step 1: Update beta
    #
    # beta | rest ~ Normal(
    #   (X'X)^(-1) X'(y-w),
    #   tau2 (X'X)^(-1)
    # )
    # ========================================================
    
    beta_mean <- drop(
      XtX_inv %*% crossprod(X, y - w)
    )
    
    beta_covariance <- tau2 * XtX_inv
    
    beta <- drop(
      beta_mean +
        t(chol(beta_covariance)) %*% rnorm(p)
    )
    
    # ========================================================
    # Step 2: Construct the spatial correlation matrix
    # ========================================================
    
    R <- spatial_correlation(
      distance_matrix = distance_matrix,
      phi = phi,
      covariance = covariance,
      nu = nu
    )
    
    chol_R <- chol(R)
    R_inverse <- chol2inv(chol_R)
    
    # ========================================================
    # Step 3: Update spatial random effects w
    #
    # Precision:
    # Qw = I/tau2 + R^(-1)/sigma2
    #
    # Mean:
    # muw = Qw^(-1)(y-X beta)/tau2
    # ========================================================
    
    Q_w <- diag(n) / tau2 +
      R_inverse / sigma2
    
    b_w <- (y - drop(X %*% beta)) / tau2
    
    w_mean <- drop(solve(Q_w, b_w))
    
    chol_Q_w <- chol(Q_w)
    
    w <- drop(
      w_mean +
        backsolve(chol_Q_w, rnorm(n))
    )
    
    # ========================================================
    # Step 4: Update sigma2
    #
    # sigma2 | rest ~ IG(
    #   a_sigma + n/2,
    #   b_sigma + w'R^(-1)w/2
    # )
    # ========================================================
    
    spatial_quadratic <- drop(
      crossprod(w, R_inverse %*% w)
    )
    
    sigma2 <- rinvgamma(
      n = 1,
      shape = prior$sigma_shape + n / 2,
      scale = prior$sigma_scale +
        spatial_quadratic / 2
    )
    
    # ========================================================
    # Step 5: Update tau2
    #
    # tau2 | rest ~ IG(
    #   a_tau + n/2,
    #   b_tau + e'e/2
    # )
    # ========================================================
    
    residual <- y - drop(X %*% beta) - w
    
    tau2 <- rinvgamma(
      n = 1,
      shape = prior$tau_shape + n / 2,
      scale = prior$tau_scale +
        drop(crossprod(residual)) / 2
    )
    
    # ========================================================
    # Step 6: Metropolis-Hastings update for phi
    # ========================================================
    
    current_log_phi <- log(phi)
    
    proposed_log_phi <- rnorm(
      1,
      mean = current_log_phi,
      sd = proposal_sd
    )
    
    current_target <- log_phi_target(
      log_phi = current_log_phi,
      w = w,
      sigma2 = sigma2,
      distance_matrix = distance_matrix,
      covariance = covariance,
      nu = nu,
      phi_shape = prior$phi_shape,
      phi_rate = prior$phi_rate
    )
    
    proposed_target <- log_phi_target(
      log_phi = proposed_log_phi,
      w = w,
      sigma2 = sigma2,
      distance_matrix = distance_matrix,
      covariance = covariance,
      nu = nu,
      phi_shape = prior$phi_shape,
      phi_rate = prior$phi_rate
    )
    
    log_acceptance_ratio <-
      proposed_target - current_target
    
    # This is the actual rejection/acceptance step
    if (
      is.finite(log_acceptance_ratio) &&
      log(runif(1)) < log_acceptance_ratio
    ) {
      phi <- exp(proposed_log_phi)
      accepted_phi <- accepted_phi + 1
    }
    
    # If the condition is false, phi is unchanged:
    # the proposed value is rejected.
    
    # ========================================================
    # Store posterior samples
    # ========================================================
    
    if (iteration %in% keep_iterations) {
      
      stored <- stored + 1
      
      samples[stored, ] <- c(
        beta,
        sigma2,
        tau2,
        phi
      )
      
      w_sum <- w_sum + w
    }
    
    if (iteration %% 5000 == 0) {
      message(
        covariance,
        " model: iteration ",
        iteration,
        " of ",
        n_iter
      )
    }
  }
  
  list(
    samples = as.mcmc(samples),
    posterior_mean_w = w_sum / n_keep,
    acceptance_phi = accepted_phi / n_iter,
    covariance = covariance,
    nu = nu,
    prior = prior
  )
}

# ============================================================
# 6. Baseline priors
# ============================================================

baseline_prior <- list(
  sigma_shape = 0.1,
  sigma_scale = 0.1,
  tau_shape = 0.1,
  tau_scale = 0.1,
  phi_shape = 0.1,
  phi_rate = 0.1
)

# Priors:
#
# sigma2 ~ IG(0.1, 0.1)
# tau2   ~ IG(0.1, 0.1)
# phi    ~ Gamma(0.1, 0.1)
#
# beta has an improper flat prior.

# ============================================================
# 7. Fit the exponential model
# ============================================================

fit_exponential <- run_spatial_mcmc(
  y = y,
  X = X,
  distance_matrix = distance_matrix,
  covariance = "exponential",
  n_iter = 30000,
  burn_in = 10000,
  thin = 10,
  prior = baseline_prior,
  proposal_sd = 0.12,
  seed = 1001
)

fit_exponential$acceptance_phi

summary(fit_exponential$samples)

# ============================================================
# 8. Fit the Matern model
# ============================================================

# nu = 1.5 is treated as fixed.
# Larger nu produces a smoother spatial process.

fit_matern <- run_spatial_mcmc(
  y = y,
  X = X,
  distance_matrix = distance_matrix,
  covariance = "matern",
  nu = 1.5,
  n_iter = 30000,
  burn_in = 10000,
  thin = 10,
  prior = baseline_prior,
  proposal_sd = 0.12,
  seed = 2001
)

fit_matern$acceptance_phi

summary(fit_matern$samples)

# ============================================================
# 9. Posterior summaries
# ============================================================

posterior_summary <- function(fit) {
  
  samples <- as.matrix(fit$samples)
  
  data.frame(
    Parameter = colnames(samples),
    Mean = apply(samples, 2, mean),
    SD = apply(samples, 2, sd),
    Lower_2.5 = apply(
      samples,
      2,
      quantile,
      probs = 0.025
    ),
    Median = apply(
      samples,
      2,
      median
    ),
    Upper_97.5 = apply(
      samples,
      2,
      quantile,
      probs = 0.975
    ),
    row.names = NULL
  )
}

exponential_summary <- posterior_summary(fit_exponential)
matern_summary <- posterior_summary(fit_matern)

exponential_summary
matern_summary

# Combine summaries
model_comparison <- bind_rows(
  Exponential = exponential_summary,
  Matern = matern_summary,
  .id = "Model"
)

model_comparison

# ============================================================
# 10. MCMC diagnostic plots
# ============================================================

# Trace plots
plot(fit_exponential$samples)
plot(fit_matern$samples)

# Autocorrelation plots
autocorr.plot(fit_exponential$samples)
autocorr.plot(fit_matern$samples)

# Effective sample sizes
effectiveSize(fit_exponential$samples)
effectiveSize(fit_matern$samples)

# Geweke convergence diagnostics
geweke.diag(fit_exponential$samples)
geweke.diag(fit_matern$samples)

# Range-parameter MH acceptance rates
fit_exponential$acceptance_phi
fit_matern$acceptance_phi







