#' Generate MLM IDA data with known parameters
#'
#' Generates data from a two-level integrative data analysis multilevel model,
#' \deqn{d_{ij} = \beta_{0j} + \beta_{1k} (X_{ij} - \bar{X}_j) + \epsilon_{ij}}
#' \deqn{\beta_{0j} = \gamma_{00} + \gamma_{B} \bar{X}_k + u_{0j}}
#' \deqn{\beta_{1j} = \gamma_{W} + u_{1j}}
#'
#' @param J The number of studies (level-2 groups).
#' @param nj A vector of sample sizes for each study with length J.
#' @param gamma_00 Fixed effect coefficient for mean intercept.
#' @param beta_w Within-study fixed effect coefficient for
#' \eqn{(X_i - \bar{X}_k)}.
#' @param beta_b Between-study fixed effect coefficient for \eqn{\bar{X}_k}.
#' @param sigma2_x Level-2 variance of covariate X (i.e., variance of the study
#'   means of X). Setting to 0 yields intraclass correlation \eqn{ICC(X) = 0}.
#' @param sigma2_r Level-1 variance of covariate X (i.e., conditional residual
#'   variance of X given study k.
#' @param sigma2_u0 Level-2 variance of random intercept coefficients.
#' @param sigma2_u1 Level-2 variance of random slope coefficients for
#'   \eqn{(X_i - \bar{X}_k)}.
#' @param sigma_u01 Level-2 covariance of random intercept and slope
#'   coefficients.
#' @param sigma2_e Level-1 conditional residual variance. Default value: 0.25.
#' @return A list containing (1) \code{data}: a tibble of generated data and (2)
#'   \code{init_params}: a list of initialization parameters \code{J, nj,
#'   gamma_00, beta_w, beta_b, sigma2_x, sigma2_r, sigma2_u0, sigma2_u1,
#'   sigma2_u01, sigma2_e}.
#'
#' @importFrom stats rnorm
#' @importFrom mvtnorm rmvnorm
#' @import dplyr
#' @importFrom readr read_csv write_csv
#'
#' @export
sim_ida <- function(
  J, nj, gamma_00 = 0.5, beta_w = 0.0, beta_b = 0.0,
  sigma2_x = 1.0, sigma2_r = 1.0, sigma2_u0 = 0.1,
  sigma2_u1 = 0.1, sigma_u01 = 0.05, sigma2_e = 0.25) {

  ## Set simulation parameters ##
  # Number of studies
  if (!is.numeric(J)) stop("J is not numeric")
  if (length(J) == 1L & !is.na(J) & is.finite(J) & J > 0 &
      identical(J %% 1, 0)) {
      J_ <- as.integer(J)
  } else {
    stop("J is not a finite positive integer.")
  }
  # Study sizes
  if (!is.numeric(nj)) stop("nj is not numeric")
  if (identical(length(nj), J_) & identical(sum(!is.finite(nj)), 0L) &
      identical(sum(nj %% 1 != 0), 0L) & identical(sum(nj <= 0), 0L)) {
    nj_ <- nj
  } else {
    stop("nj is not a vector of finite positive integers.")
  }
  # Check gamma_00
  if (!is.numeric(gamma_00)) stop("gamma_00 is not numeric")
  if (identical(length(gamma_00), 1L) & is.finite(gamma_00)) {
    gamma_00_ <- gamma_00
  } else {
    stop("gamma_00 is not a finite number.")
  }
  # Check beta_w
  if (!is.numeric(beta_w)) stop("beta_w is not numeric")
  if (identical(length(beta_w), 1L) & is.finite(beta_w)) {
      beta_w_ <- beta_w
  } else {
    stop("beta_w is not a finite number.")
  }
  # Check beta_b
  if (!is.numeric(beta_b)) stop("beta_b is not numeric")
  if (identical(length(beta_b), 1L) & is.finite(beta_b)) {
    beta_b_ <- beta_b
  } else {
    stop("beta_b is not a finite number.")
  }
  # Check sigma2_x
  if (!is.numeric(sigma2_x)) stop("sigma2_x is not numeric")
  if (identical(length(sigma2_x), 1L) & !is.na(sigma2_x) &
      is.finite(sigma2_x) & sigma2_x >= 0.0) {
    sigma2_x_ <- sigma2_x
    if (near(sigma2_x, 0.0)) warning("sigma2_x is effectively 0.0")
  } else {
    stop("sigma2_x is not a finite non-negative number.")
  }
  # Check sigma_2r
  if (!is.numeric(sigma2_r)) stop("sigma2_r is not numeric")
  if (identical(length(sigma2_r), 1L) & is.finite(sigma2_r) & sigma2_r >= 0.0) {
    sigma2_r_ <- sigma2_r
    if (near(sigma2_r, 0.0)) warning("sigma2_r is effectively 0.0")
  } else {
    stop("sigma2_r is not a finite non-negative number.")
  }
  # Check sigma2_u0
  if (!is.numeric(sigma2_u0)) stop("sigma2_u0 is not numeric")
  if (identical(length(sigma2_u0), 1L) & is.finite(sigma2_u0) &
      sigma2_u0 >= 0.0) {
    sigma2_u0_ <- sigma2_u0
    if (near(sigma2_u0, 0.0)) warning("sigma2_u0 is effectively 0.0")
  } else {
    stop("sigma2_u0 is not a non-negative number.")
  }
  # Check sigma2_u1
  if (!is.numeric(sigma2_u1)) stop("sigma2_u1 is not numeric")
  if (identical(length(sigma2_u1), 1L) & is.finite(sigma2_u1) &
      sigma2_u1 >= 0.0) {
    sigma2_u1_ <- sigma2_u1
    if (near(sigma2_u1, 0.0)) warning("sigma2_u1 is effectively 0.0")
  } else {
    stop("sigma2_u1 is not a finite non-negative number.")
  }
  # Check sigma_u01
  if (!is.numeric(sigma_u01)) stop("sigma_u01 is not numeric")
  if (identical(length(sigma_u01), 1L) & is.finite(sigma_u01)) {
    sigma_u01_ <- sigma_u01
  } else {
    stop("sigma2_u01 is not finite.")
  }
  # Check sigma2_e
  if (!is.numeric(sigma2_e)) stop("sigma2_e is not numeric")
  if (identical(length(sigma2_e), 1L) & is.finite(sigma2_e) &
      sigma2_e >= 0.0) {
    sigma2_e_ <- sigma2_e
    if (near(sigma2_e, 0.0)) warning("sigma2_e is effectively 0.0")
  } else {
    stop("sigma2_e is not a finite positive number.")
  }
  # Check for zero L1 or L2 random effect variance parameters
  if (near(sigma2_u0_, 0) | near(sigma2_u1_, 0)) {
    if (!near(sigma_u01, 0.0)) {
      warning("sigma_u01 was set to 0 since either sigma2_u0 or sigma2_u1 was 0")
    }
    sigma_u01_ <- sigma_u01 <- 0.0
  }
  # Ensure level-2 covariance matrix is positive definite
  sigl2 = matrix(c(sigma2_u0, sigma_u01,
                   sigma_u01, sigma2_u1), nrow = 2)
  if (sum(eigen(sigl2)$values > 0.0) != 2) {
    warning("Level-2 covariance matrix is not positive definite")
  }

  ## Draw covariate X
  # Draw study population means for X
  mu_x = draw_covmeans(J = J_, sigma2_x = sigma2_x_)
  # Draw n_k observations of covariate for each study
  dat_df = draw_l1cov(J = J_, nj = nj_, mu_x = mu_x, sigma2_r = sigma2_r_)

  ## Compute Xbar and X - Xbar for each study
  dat_df <- dat_df %>%
    group_by(k) %>%
    mutate(xbar = mean(x),
           grp_centered = x - xbar) %>%
    ungroup()
  xbars <- dat_df %>%
    arrange(k) %>%
    distinct(k, xbar) %>%
    select(xbar) %>%
    unlist()

  ## Draw random intercept and slope betas
  betas <- tibble(k = seq(J_), beta_0k = NaN, beta_1k = NaN)
  for (k in seq(J)) {
    betas[k, c("beta_0k", "beta_1k")] <- draw_rancoef(
      gamma_00 = gamma_00_, beta_b = beta_b_, beta_w = beta_w_,
      l2cov = sigl2, xbar = xbars[k])
  }

  ## Draw n_k observations of d
  dat_df <- dat_df %>%
    left_join(betas, by = "k") %>%
    rowwise() %>%
    mutate(delta_ik = rnorm(
      1, beta_0k + beta_1k * grp_centered, sqrt(sigma2_e_))) %>%
    ungroup()

  return(
    list(
      data = dat_df,
      init_params = list(
        "J" = J, "nj" = nj,
        "gamma_00" = gamma_00, "beta_w" = beta_w, "beta_b" = beta_b,
        "sigma2_x" = sigma2_x, "sigma2_r" = sigma2_r, "sigma2_u0" = sigma2_u0,
        "sigma2_u1" = sigma2_u1, "sigma2_u01" = sigma_u01,
        "sigma2_e" =  sigma2_e)
    )
  )
}

#' Draw continous predictor means
#'
#' Draw normal-distributed level-1 covariate means for J studies
#'
#' @inheritParams sim_ida_r2
#' @return A vector of study means for the covariate.
#'
#' @importFrom stats rnorm
draw_covmeans = function(J, sigma2_x) {
  if (sigma2_x > 0.0) {
    mu_x <- rnorm(J, 0, sqrt(sigma2_x))
  } else {
    mu_x <- rep(0.0, J)
  }
  return(mu_x)
}

#' Draw continuous covariate
#'
#' Draw normal-distributed level-1 covariate for J studies
#'
#' @inheritParams sim_ida_r2
#' @param mu_x Vector of study means for the covariate.
#' @return A data frame containing the level-1 covariate \code{x} and study
#'   index \code{k}.
#'
#' @importFrom tibble tibble
draw_l1cov <- function(J, nj, mu_x, sigma2_r) {
  N = sum(nj)
  x <- numeric(N)
  dat_df <- tibble(x = NaN, k = NaN)
  obs <- 1
  for (k in seq(J)) {
    x[obs:(obs + nj[k] - 1)] <- rnorm(nj[k], mu_x[k], sqrt(sigma2_r))
    dat_df[obs:(obs + nj[k] - 1), "x"] <- x[obs:(obs + nj[k] - 1)]
    dat_df[obs:(obs + nj[k] - 1), "k"] <- k
    obs <- obs + nj[k]
  }
  return(dat_df)
}

#' Draw balanced dichotomous predictor
#'
#' Draw dichotomous balanced level-1 covariate for J studies
#'
#' @inheritParams sim_ida_r2_dichotx
#' @return A data frame containing the level-1 dichotomous predictor \code{x}
#'   and study index \code{j}.
#'
#' @importFrom tibble tibble
draw_l1cov_dichotx <- function(J, nj) {
  N = sum(nj)
  x <- numeric(N) # Initialize to 0
  dat_df <- tibble(x = NaN, j = NaN)
  obs <- 1
  for (j in seq(J)) {
    ind_one <- sample(seq_len(nj[j]), size = floor(nj[j] / 2), replace = FALSE)
    x[obs:(obs + nj[j] - 1)][ind_one] <- 1L
    dat_df[obs:(obs + nj[j] - 1), "x"] <- x[obs:(obs + nj[j] - 1)]
    dat_df[obs:(obs + nj[j] - 1), "j"] <- j
    obs <- obs + nj[j]
  }
  return(dat_df)
}

#' Draw random coefficients for a level-2 group
#'
#' @inheritParams sim_ida
#' @param l2cov A 2 x 2 symmetric matrix of level-2 covariance parameters.
#' @param xbar A level-1-covariate mean for one study.
#' @return Two-element vector of random intercept and random slope.
#'
#' @importFrom mvtnorm rmvnorm
#' @importFrom dplyr near
draw_rancoef <- function(gamma_00, beta_b, beta_w, l2cov, xbar) {
  if (
    (near(l2cov[1, 1], 0.0) & l2cov[1, 2] > 0.0) |
    (near(l2cov[1, 1], 0.0) & l2cov[2, 1] > 0.0) |
    (near(l2cov[2, 2], 0.0) & l2cov[1, 2] > 0.0) |
    (near(l2cov[2, 2], 0.0) & l2cov[2, 1] > 0.0) ) {
    stop("Variances near 0 with covariances != 0 detected")
  }
  betas <- numeric(2)
  betas <- rmvnorm(
    1,
    c(gamma_00 + beta_b * xbar,
      beta_w),
    l2cov)
  return(betas)
}

#' Draw random coefficients for a level-2 group
#'
#' @inheritParams sim_ida_r2_dichotx
#' @param gamma_10 Fixed effect coefficient for \eqn{X_{ij}}.
#' @param l2cov A 2 x 2 symmetric matrix of level-2 covariance parameters.
#' @return Two-element vector of random intercept and random slope.
#'
#' @importFrom mvtnorm rmvnorm
#' @importFrom dplyr near
draw_rancoef_dichotx <- function(gamma_00, gamma_10, l2cov) {
  if (
    (near(l2cov[1, 1], 0.0) & l2cov[1, 2] > 0.0) |
    (near(l2cov[1, 1], 0.0) & l2cov[2, 1] > 0.0) |
    (near(l2cov[2, 2], 0.0) & l2cov[1, 2] > 0.0) |
    (near(l2cov[2, 2], 0.0) & l2cov[2, 1] > 0.0) ) {
    stop("Variances near 0 with covariances != 0 detected")
  }
  betas <- numeric(2)
  betas <- rmvnorm(
    n = 1L,
    mean = c(gamma_00, gamma_10),
    sigma = l2cov)
  return(betas)
}

#' Generate MLM IDA data based on R^2 effect size framework (Rights & Sterba, 2019)
#'
#' \code{sim_ida_r2} generates data from a two-level integrative data analysis
#' multilevel model,
#' \deqn{d_{ij} = \beta_{0j} + \beta_{1k} (X_{ij} - \bar{X}_j) + \epsilon_{ij}}
#' \deqn{\beta_{0j} = \gamma_{00} + \gamma_{B} \bar{X}_k + u_{0j}}
#' \deqn{\beta_{1j} = \gamma_{W} + u_{1j}}
#'
#' @param J The number of studies (level-2 groups).
#' @param nj A vector of sample sizes for each study with length J.
#' @param gamma_00 Fixed effect coefficient for mean intercept.
#' @param gamma_w_es Effect size of within-study fixed effect coefficient for
#' \eqn{(X_{ij} - \bar{X}_j)}.
#' @param gamma_b_es Effect size of Between-study fixed effect coefficient for \eqn{\bar{X}_j}.
#' @param sigma2_x Level-2 variance of covariate X (i.e., variance of the study
#'   means of X). Setting to 0 yields intraclass correlation \eqn{ICC(X) = 0}.
#' @param sigma2_r Level-1 variance of covariate X (i.e., conditional residual
#'   variance of X given study k.
#' @param sigma2_u0_es Effect size of Level-2 variance of random intercept coefficients.
#' @param sigma2_u1_es Effect size of Level-2 variance of random slope coefficients for
#'   \eqn{(X_{ij} - \bar{X}_j)}.
#' @param sigma_u01_es Level-2 correlation of random intercept and slope
#'   coefficients.
#' @param total_vary Total variance of outcome \eqn{y_{ij}}. Default value: 1.
#' @return A list containing (1) \code{data}: a tibble of generated data and (2)
#'   \code{init_params}: a list of initialization parameters \code{J, nj,
#'   gamma_00, beta_w, beta_b, sigma2_x, sigma2_r, sigma2_u0, sigma2_u1,
#'   sigma2_u01, sigma2_e}.
#'
#' @importFrom stats rnorm
#' @importFrom mvtnorm rmvnorm
#' @import dplyr
#' @importFrom readr read_csv write_csv
#'
#' @export
sim_ida_r2 <- function(
  J, nj, gamma_00 = 0, gamma_w_es = 0.0, gamma_b_es = 0.0,
  sigma2_x = 1.0, sigma2_r = 1.0, sigma2_u0_es = 0.0,
  sigma2_u1_es = 0.0, sigma_u01_es = 0.0, total_vary = 1.0) {

  ## Set simulation parameters ##
  # Number of studies
  if (!is.numeric(J)) stop("J is not numeric")
  if (length(J) == 1L & !is.na(J) & is.finite(J) & J > 0 &
      identical(J %% 1, 0)) {
    J_ <- as.integer(J)
  } else {
    stop("J is not a finite positive integer.")
  }
  # Study sizes
  if (!is.numeric(nj)) stop("nj is not numeric")
  if (identical(length(nj), J_) & identical(sum(!is.finite(nj)), 0L) &
      identical(sum(nj %% 1 != 0), 0L) & identical(sum(nj <= 0), 0L)) {
    nj_ <- nj
  } else {
    stop("nj is not a vector of finite positive integers.")
  }
  # Check gamma_00
  if (!is.numeric(gamma_00)) stop("gamma_00 is not numeric")
  if (identical(length(gamma_00), 1L) & is.finite(gamma_00)) {
    gamma_00_ <- gamma_00
  } else {
    stop("gamma_00 is not a finite number.")
  }
  # Check beta_w_es
  if (!is.numeric(gamma_w_es)) stop("beta_w is not numeric")
  if (identical(length(gamma_w_es), 1L) & is.finite(gamma_w_es)) {
    gamma_w_es_ <- gamma_w_es
  } else {
    stop("gamma_w_es is not a finite number.")
  }
  # Check beta_b_es
  if (!is.numeric(gamma_b_es)) stop("gamma_b is not numeric")
  if (identical(length(gamma_b_es), 1L) & is.finite(gamma_b_es)) {
    gamma_b_es_ <- gamma_b_es
  } else {
    stop("gamma_b_es is not a finite number.")
  }
  # Check sigma2_x
  if (!is.numeric(sigma2_x)) stop("sigma2_x is not numeric")
  if (identical(length(sigma2_x), 1L) & !is.na(sigma2_x) &
      is.finite(sigma2_x) & sigma2_x >= 0.0) {
    sigma2_x_ <- sigma2_x
    if (near(sigma2_x, 0.0)) warning("sigma2_x is effectively 0.0")
  } else {
    stop("sigma2_x is not a finite non-negative number.")
  }
  # Check sigma_2r
  if (!is.numeric(sigma2_r)) stop("sigma2_r is not numeric")
  if (identical(length(sigma2_r), 1L) & is.finite(sigma2_r) & sigma2_r >= 0.0) {
    sigma2_r_ <- sigma2_r
    if (near(sigma2_r, 0.0)) warning("sigma2_r is effectively 0.0")
  } else {
    stop("sigma2_r is not a finite non-negative number.")
  }
  # Check sigma2_u0_es
  if (!is.numeric(sigma2_u0_es)) stop("sigma2_u0 is not numeric")
  if (identical(length(sigma2_u0_es), 1L) & is.finite(sigma2_u0_es) &
      sigma2_u0_es >= 0.0) {
    sigma2_u0_es_ <- sigma2_u0_es
    if (near(sigma2_u0_es, 0.0)) warning("sigma2_u0_es is effectively 0.0")
  } else {
    stop("sigma2_u0_es is not a non-negative number.")
  }
  # Check sigma2_u1_es
  if (!is.numeric(sigma2_u1_es)) stop("sigma2_u1 is not numeric")
  if (identical(length(sigma2_u1_es), 1L) & is.finite(sigma2_u1_es) &
      sigma2_u1_es >= 0.0) {
    sigma2_u1_es_ <- sigma2_u1_es
    if (near(sigma2_u1_es, 0.0)) warning("sigma2_u1_es is effectively 0.0")
  } else {
    stop("sigma2_u1_es is not a finite non-negative number.")
  }
  # Check sigma_u01_es
  if (!is.numeric(sigma_u01_es)) stop("sigma_u01 is not numeric")
  if (identical(length(sigma_u01_es), 1L) & is.finite(sigma_u01_es)) {
    sigma_u01_es_ <- sigma_u01_es
  } else {
    stop("sigma2_u01_es is not finite.")
  }

  # Check for zero L1 or L2 random effect variance parameters
  if (near(sigma2_u0_es_, 0) | near(sigma2_u1_es_, 0)) {
    if (!near(sigma_u01_es, 0.0)) {
      warning("sigma_u01 was set to 0 since either sigma2_u0_es or sigma2_u1_es was 0")
    }
    sigma_u01_es_ <- sigma_u01_es <- 0.0
  }

  if (sum(gamma_w_es_ + gamma_b_es_ + sigma2_u0_es_ + sigma2_u1_es_) > 1)
    stop("Total R2 effect size cannot exceed 1.")

  # Calculate model parameters based on R2 effect size
  sigma2_u0 <- sigma2_u0_es_ * total_vary
  sigma2_u1 <- sigma2_u1_es_ * total_vary / (sigma2_x + sigma2_r)
  sigma_u01 <- sigma_u01_es_ * sqrt(sigma2_u0) * sqrt(sigma2_u1)
  gamma_w <- sqrt(gamma_w_es_ * total_vary / (sigma2_x + sigma2_r))
  gamma_b <- sqrt(gamma_b_es_ * total_vary / sigma2_x)
  sigma2_e <- total_vary - sigma2_u0 -
    (sigma2_u1 + gamma_w ^ 2) * (sigma2_x + sigma2_r) - gamma_b ^ 2 * sigma2_x

  # Ensure level-2 covariance matrix is positive definite
  sigl2 = matrix(c(sigma2_u0, sigma_u01,
                   sigma_u01, sigma2_u1), nrow = 2)
  if (sum(eigen(sigl2)$values > 0.0) != 2) {
    warning("Level-2 covariance matrix is not positive definite")
  }

  ## Draw covariate X
  # Draw study population means for X
  mu_x = draw_covmeans(J = J_, sigma2_x = sigma2_x_)
  # Draw n_k observations of covariate for each study
  dat_df = draw_l1cov(J = J_, nj = nj_, mu_x = mu_x, sigma2_r = sigma2_r_)

  ## Compute Xbar and X - Xbar for each study
  dat_df <- dat_df %>%
    group_by(k) %>%
    mutate(xbar = mean(x),
           grp_centered = x - xbar) %>%
    ungroup()
  xbars <- dat_df %>%
    arrange(k) %>%
    distinct(k, xbar) %>%
    select(xbar) %>%
    unlist()

  ## Draw random intercept and slope betas
  betas <- tibble(k = seq(J_), beta_0j = NaN, beta_1j = NaN)
  for (j in seq(J)) {
    betas[j, c("beta_0j", "beta_1j")] <- draw_rancoef(
      gamma_00 = gamma_00_, beta_b = gamma_b, beta_w = gamma_w,
      l2cov = sigl2, xbar = xbars[j])
  }

  ## Draw n_k observations of d
  dat_df <- dat_df %>%
    left_join(betas, by = "k") %>%
    rowwise() %>%
    mutate(delta_ij = rnorm(
      1, beta_0j + beta_1j * grp_centered, sqrt(sigma2_e))) %>%
    ungroup()

  return(
    list(
      data = dat_df,
      init_params = list(
        "J" = J, "nj" = nj,
        "gamma_00" = gamma_00, "gamma_w" = gamma_w, "gamma_b" = gamma_b,
        "sigma2_x" = sigma2_x, "sigma2_r" = sigma2_r, "sigma2_u0" = sigma2_u0,
        "sigma2_u1" = sigma2_u1, "sigma2_u01" = sigma_u01,
        "sigma2_e" =  sigma2_e)
    )
  )
}

#' Generate MLM IDA data with a dichotomous predictor based on R^2 effect sizes framework (Rights & Sterba, 2019)
#'
#' \code{sim_ida_r2_dichotx} generates data from a two-level multilevel model
#'   with a single dichotomous predictor,
#' \deqn{y_{ij} = \beta_{0j} + \beta_{1j} x_{ij} + \epsilon_{ij}}
#' \deqn{\beta_{0j} = \gamma_{00} + u_{0j}}
#' \deqn{\beta_{1j} = \gamma_{10} + u_{1j}}
#'
#' @param J The number of studies (level-2 groups).
#' @param nj A vector of sample sizes for each study with length J.
#' @param gamma_00 Fixed effect coefficient for mean intercept.
#' @param gamma_10_es Effect size of fixed effect coefficient for \eqn{X_{ij}}.
#' @param sigma2_x Variance of dichotomous predictor. Default value: 0.25
#'   (corresponding to balanced classes).
#' @param sigma2_u0_es Effect size of Level-2 variance of random intercept coefficients.
#' @param sigma2_u1_es Effect size of Level-2 variance of random slope coefficients for
#'   \eqn{(X_{ij} - \bar{X}_j)}.
#' @param total_vary Total variance of outcome \eqn{y_{ij}}. Default value: 1.
#' @param rho_01 Level-2 correlation of random intercept and slope
#'   coefficients.
#' @return A list containing (1) \code{data}: a tibble of generated data and (2)
#'   \code{init_params}: a list of initialization parameters \code{J, nj,
#'   gamma_00, gamma_10, sigma2_u0, sigma2_u1, rho_01, sigma2_e}.
#'
#' @importFrom stats rnorm
#' @importFrom mvtnorm rmvnorm
#' @import dplyr
#' @importFrom readr read_csv write_csv
#'
#' @export
sim_ida_r2_dichotx <- function(
  J, nj, gamma_00 = 0, gamma_10_es = 0.0, sigma2_x = 0.25,
  sigma2_u0_es = 0.0, sigma2_u1_es = 0.0, rho_01 = 0.0, total_vary = 1.0) {

  ## Set simulation parameters ##
  # Number of studies
  if (!is.numeric(J)) stop("J is not numeric")
  if (length(J) == 1L & !is.na(J) & is.finite(J) & J > 0 &
      identical(J %% 1, 0)) {
    J_ <- as.integer(J)
  } else {
    stop("J is not a finite positive integer.")
  }
  # Study sizes
  if (!is.numeric(nj)) stop("nj is not numeric")
  if (identical(length(nj), J_) & identical(sum(!is.finite(nj)), 0L) &
      identical(sum(nj %% 1 != 0), 0L) & identical(sum(nj <= 0), 0L)) {
    nj_ <- nj
  } else {
    stop("nj is not a vector of finite positive integers.")
  }
  # Check gamma_00
  if (!is.numeric(gamma_00)) stop("gamma_00 is not numeric")
  if (identical(length(gamma_00), 1L) & is.finite(gamma_00)) {
    gamma_00_ <- gamma_00
  } else {
    stop("gamma_00 is not a finite number.")
  }
  # Check gamma_10_es
  if (!is.numeric(gamma_10_es)) stop("beta_w is not numeric")
  if (identical(length(gamma_10_es), 1L) & is.finite(gamma_10_es)) {
    gamma_10_es_ <- gamma_10_es
  } else {
    stop("gamma_10_es is not a finite number.")
  }
  # Check sigma2_x
  if (!is.numeric(sigma2_x)) stop("sigma2_x is not numeric")
  if (identical(length(sigma2_x), 1L) & !is.na(sigma2_x) &
      is.finite(sigma2_x) & sigma2_x >= 0.0) {
    sigma2_x_ <- sigma2_x
    if (near(sigma2_x, 0.0)) warning("sigma2_x is effectively 0.0")
  } else {
    stop("sigma2_x is not a finite non-negative number.")
  }
  # Check sigma2_u0_es
  if (!is.numeric(sigma2_u0_es)) stop("sigma2_u0 is not numeric")
  if (identical(length(sigma2_u0_es), 1L) & is.finite(sigma2_u0_es) &
      sigma2_u0_es >= 0.0) {
    sigma2_u0_es_ <- sigma2_u0_es
    if (near(sigma2_u0_es, 0.0)) warning("sigma2_u0_es is effectively 0.0")
  } else {
    stop("sigma2_u0_es is not a non-negative number.")
  }
  # Check sigma2_u1_es
  if (!is.numeric(sigma2_u1_es)) stop("sigma2_u1 is not numeric")
  if (identical(length(sigma2_u1_es), 1L) & is.finite(sigma2_u1_es) &
      sigma2_u1_es >= 0.0) {
    sigma2_u1_es_ <- sigma2_u1_es
    if (near(sigma2_u1_es, 0.0)) warning("sigma2_u1_es is effectively 0.0")
  } else {
    stop("sigma2_u1_es is not a finite non-negative number.")
  }
  # Check rho_01
  if (!is.numeric(rho_01)) stop("rho_01 is not numeric")
  if (identical(length(rho_01), 1L) & is.finite(rho_01) & between(rho_01, -1.0, 1.0)) {
    rho_01_ <- rho_01
  } else {
    stop("rho_01 is not finite.")
  }

  # Check for zero L1 or L2 random effect variance parameters
  if (near(sigma2_u0_es_, 0) | near(sigma2_u1_es_, 0)) {
    if (!near(rho_01, 0.0)) {
      warning("rho_01 was set to 0 since either sigma2_u0_es or sigma2_u1_es was 0")
    }
    rho_01_ <- rho_01 <- 0.0
  }

  if (sum(gamma_10_es_ + sigma2_u0_es_ + sigma2_u1_es_) > 1)
    stop("Total R2 effect size cannot exceed 1.")

  # Calculate model parameters based on R2 effect size
  sigma2_u0 <- sigma2_u0_es_ * total_vary
  sigma2_u1 <- sigma2_u1_es_ * total_vary / sigma2_x
  sigma_u01 <- rho_01_ * sqrt(sigma2_u0) * sqrt(sigma2_u1)
  gamma_10 <- sqrt(gamma_10_es_ * total_vary / sigma2_x)
  sigma2_e <- total_vary - sigma2_u0 - (sigma2_u1 + gamma_10 ^ 2) * sigma2_x

  # Ensure level-2 covariance matrix is positive definite
  sigl2 = matrix(c(sigma2_u0, sigma_u01,
                   sigma_u01, sigma2_u1), nrow = 2)
  if (sum(eigen(sigl2)$values > 0.0) != 2) {
    warning("Level-2 covariance matrix is not positive definite")
  }

  ## Draw covariate X
  dat_df = draw_l1cov_dichotx(J = J_, nj = nj_)

  ## Draw random intercept and slope betas
  betas <- tibble(j = seq_len(J_), beta_0j = NaN, beta_1j = NaN)
  for (j in seq(J)) {
    betas[j, c("beta_0j", "beta_1j")] <- draw_rancoef_dichotx(
      gamma_00 = gamma_00_, gamma_10 = gamma_10,
      l2cov = sigl2)
  }

  ## Draw n_k observations of y
  dat_df <- dat_df %>%
    left_join(betas, by = "j") %>%
    rowwise() %>%
    mutate(y_ij = rnorm(
      n = 1, mean = beta_0j + beta_1j * x, sd = sqrt(sigma2_e))) %>%
    ungroup()

  return(
    list(
      data = dat_df,
      init_params = list(
        "J" = J, "nj" = nj,
        "gamma_00" = gamma_00, "gamma_10" = gamma_10,
        "sigma2_u0" = sigma2_u0, "sigma2_u1" = sigma2_u1, "rho_01" = rho_01,
        "sigma2_e" =  sigma2_e)
    )
  )
}
