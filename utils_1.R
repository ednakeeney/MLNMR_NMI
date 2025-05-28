# Utility functions for survival feasibility assessment

# Generate all pairs of powers of fractional polynomials
generate_fp_so_powers <- function(fp_powers) {
  fp_so_powers <- as.data.frame(matrix(nrow = sum(1:length(fp_powers)),
                         ncol = 2))
  colnames(fp_so_powers) <- c("P1", "P2")
  fp_so_powers[, "P1"] <- rep(fp_powers, times = seq(from = length(fp_powers), to = 1, by = -1))
  for(i in 1:length(fp_powers)) {
    fp_so_powers[fp_so_powers[, "P1"] == fp_powers[i], "P2"] <- 
      fp_powers[i:length(fp_powers)]
  }
  return(fp_so_powers)
}

# Time transformations for second order FP
fp_transform_times <- function(times, powers = c(-1, 0.5)) 
{
  nobs <- length(times)
  npoly <- length(powers)
  timetrans <- matrix(0, nrow = nobs, ncol = npoly)
  timetrans1 <- ifelse(powers[1] != rep(0, nobs), times^powers[1], log(times))
  timetrans[, 1] <- timetrans1
  if (npoly >= 2) {
    for (i in 2:npoly) {
      if (powers[i] == powers[(i - 1)])
        timetrans2 <- log(times) * timetrans1
      else timetrans2 <- ifelse(powers[i] != rep(0, nobs), times^powers[i],
                                log(times))
      timetrans[, i] <- timetrans2
      timetrans1 <- timetrans2
    }
  }
  return(timetrans)
}

# Function to generate initial values for FP models
fp_inits <- function(t, aux){
  mh <- muhaz::muhaz(t, max.time = max(t))
  lhdat <- data.frame(loghaz = log(mh$haz.est),
                      time = mh$est.grid)
  lhdat <- na.omit(lhdat)
  lhdat <- lhdat[is.finite(lhdat$loghaz),]
  X <- flexsurv:::bfp(lhdat$time, powers=aux$powers)
  coef(lm(loghaz ~ X, data=lhdat))
}

# Format a vector of mean, 2.5% CI and 97.5% CI
format_results <- function(x, ndigits = 2) {
  paste0(format(x[1], digits = ndigits, nsmall = ndigits), " (",
         format(x[2], digits = ndigits, nsmall = ndigits), ", ",
         format(x[3], digits = ndigits, nsmall = ndigits), ")")
}

format_results_se <- function(x, ndigits = 2) {
  paste0(format(x[1], digits = ndigits, nsmall = ndigits), " (",
         format(x[2], digits = ndigits, nsmall = ndigits), ")")
}

# Function to summarise relative effects from multinma
# Designed for survival outcomes
format_relative_effects <- function(multinma_rel, hazard_scale = TRUE,
                                    invert = TRUE,
                                    treatment_names = NULL) {
  multinma_rel <- as.data.frame(multinma_rel$summary)
  formatted_results <- matrix(NA, nrow = nrow(multinma_rel))
  for(i in 1:nrow(multinma_rel)) {
    if(invert){
      if(hazard_scale) {
        # Exponentiate the log hazard ratios
        formatted_results[i, ] <- 
          format_results(exp(-multinma_rel[i, c("mean", "97.5%", "2.5%")])  )
      } else {
        formatted_results[i, ] <- 
          format_results(-multinma_rel[i, c("mean", "97.5%", "2.5%")])
      }
    } else {
      if(hazard_scale) {
        # Exponentiate the log hazard ratios
        formatted_results[i, ] <- 
          format_results(exp( multinma_rel[i, c("mean", "2.5%", "97.5%")])  )
      } else {
        formatted_results[i, ] <- 
          format_results(multinma_rel[i, c("mean", "2.5%", "97.5%")])
      }
    }
  }
  colnames(formatted_results) <- c("Mean (95% CrI)")
  # Use treatment names if provided.
  if(!is.null(treatment_names)) {
    rownames(formatted_results) <- treatment_names
  } else {
    rownames(formatted_results) <- multinma_rel[, "parameter"]
  }
  return(formatted_results)
}

