# ML-NMR using multinma package
# Howard Thom 14-January-2024

# Updated from v4 to a 'new' analysis
# This assumes SNA+ADT!=ADT and incorporates new baseline characteristics and KM data

# If we need the developer version of multinma best to take from r-universe
# install.packages("devtools")
# install.packages("multinma", repos = c("https://dmphillippo.r-universe.dev", getOption("repos")))

library(multinma)
library(survival)
library(dplyr)
library(ggplot2)
library(loo)
library(readxl)

# Set to use more than one core to speed up calculation
max_cores <- parallel::detectCores()
options(mc.cores = max_cores)

# Directory from which to load ARASENS data and save any objects that may embed the IPD
arasens_directory <- "C:/Users/thomh/Thoms Statistical Consultants Limited/Clifton Insight SharePoint - Documents/Clients Projects/Bayer/Prostate Cancer - Darolutamide/ARASENS Trial/ARASENS IPD"

#arasens_directory <- "C:/Users/DaisyGaunt/OneDrive - Thoms Statistical Consultants Limited/Prostate Cancer - Darolutamide/ARASENS Trial/ARASENS IPD"

# If using actual data
load(file = paste0(arasens_directory, "/arasens_ipd_maic.R"))

# If using dummy data
# load(file = "data/arasens_dummy_ipd.rda")
# arasens_ipd <- arasens_dummy_ipd

source("code/utils_1.R")

# Rename CROD to PFS and remove outcomes other than OS and PFS
arasens_ipd$PFS <- arasens_ipd$CROD
arasens_ipd$CRPC <- NULL
arasens_ipd$CROD <- NULL

outcome_names <- c("OS", "PFS")

agd_filename_mapping <- readxl::read_excel(
  path = "data/IPD reconstructed new/mlnmr_file_mapping.xlsx", 
  sheet = "filename_mapping")

agd_study_names <- unique(agd_filename_mapping$mapped_study)

treatment_names <- c(unique(agd_filename_mapping$mapped_trt), "DAR+ADT+Doc")

covariate_names <- c("AGE", "ECOG_0", "DOCELI", "GLEASON_M8", 
                     "VOLUME_HIGH", "M1", "RACE_WHITE", "VISC_MET" )



#################################################################################
# ML-NMR data preparation
#################################################################################


# List of aggregate data 
mcspc_agd <- list()
mcspc_agd$OS <- mcspc_agd$PFS <- 
  matrix(NA, nrow = 0, ncol = 6,
         dimnames = list(NULL, 
                         c("study", "studyf", "trt", "trtf", 
                           "eventtime", "status")))

for(i_file in 1:dim(agd_filename_mapping)[1]) {
  # Load the reconstructed IPD
  agd_temp <- read.csv(paste0("data/IPD reconstructed new/", agd_filename_mapping$filename[i_file]))
  outcome_name <- agd_filename_mapping$outcome[i_file]
  treatment_name <- agd_filename_mapping$mapped_trt[i_file]
  study_name <- agd_filename_mapping$mapped_study[i_file]    
  
  # status is 1 if observed and 0 if censored
  colnames(agd_temp)[is.element(colnames(agd_temp), c("time", "event"))] <-
    c("eventtime", "status")
  # Drop columns not needed for ML-NMR
  agd_temp <-
    agd_temp[, c("eventtime", "status")]
  
  n_patients <- dim(agd_temp)[1]
  
  # Create study and treatment indicator dataset
  temp_study_treatment_indicator <- matrix(c(
    rep(study_name, n_patients),
    rep(study_name, n_patients),
    rep(treatment_name, n_patients),
    rep(treatment_name, n_patients)),
    nrow = n_patients, ncol = 4, 
    dimnames = list(NULL, c("study", "studyf", "trt", "trtf"))
  )
  
  agd_temp <- cbind(temp_study_treatment_indicator, agd_temp)
  
  # Combine with rest of aggregate data studies
  mcspc_agd[[outcome_name]] <- rbind(mcspc_agd[[outcome_name]],
                                     agd_temp)
}

# Ensure the aggregate datasets are dataframes and use the correct data type
for(outcome_name in outcome_names) {
  mcspc_agd[[outcome_name]] <- as.data.frame(mcspc_agd[[outcome_name]])
  rownames(mcspc_agd[[outcome_name]]) <- c(1:(dim(mcspc_agd[[outcome_name]])[1]))
}

# Store the aggregate data for sensitivity analyses as well
mcspc_sens_agd <- mcspc_agd

# Aggregate data baseline characteristics
mcspc_agd_covs <- readxl::read_excel(path = "data/agd_baseline_characteristics_new.xlsx", sheet = "baseline_clean")
# Convert to data frame
mcspc_agd_covs <- as.data.frame(mcspc_agd_covs)
for(covariate_name in covariate_names) {
  mcspc_agd_covs[, covariate_name] <- as.numeric(mcspc_agd_covs[, covariate_name])
}

# Convert from percentage (0 to 100) to proportion (0 to 1)
mcspc_agd_covs[, c("ECOG_0", "DOCELI", "GLEASON_M8", 
                   "VOLUME_HIGH", "M1", "RACE_WHITE", "VISC_MET" )] <-
  mcspc_agd_covs[, c("ECOG_0", "DOCELI", "GLEASON_M8", 
                     "VOLUME_HIGH", "M1", "RACE_WHITE", "VISC_MET" )] / 100

# Set missing covariates to mean of other studies (including ARASENS)
for(covariate_name in covariate_names) {
  mcspc_agd_covs[is.na(mcspc_agd_covs[, covariate_name]), covariate_name] <-
    mean(mcspc_agd_covs[, covariate_name], na.rm = TRUE)
}

# Remove ARASENS from aggregate data covariates
# But keep a full dataset for use in population specification
mcspc_agd_covs_full <- mcspc_agd_covs
mcspc_agd_covs <- mcspc_agd_covs[-which(mcspc_agd_covs$study == "ARASENS"), ]



# Define the IPD data 
arasens_ipd$OS$ECOG_0 <- as.numeric(arasens_ipd$OS$ECOGBL == 1)
arasens_ipd$PFS$ECOG_0 <- as.numeric(arasens_ipd$PFS$ECOGBL == 1)
arasens_ipd$OS$DOCELI <- 1
arasens_ipd$PFS$DOCELI <- 1
arasens_ipd$OS$GLEASON_M8 <- as.numeric(arasens_ipd$OS$GLEASEN >= 8)
arasens_ipd$PFS$GLEASON_M8 <- as.numeric(arasens_ipd$PFS$GLEASEN >= 8)
arasens_ipd$OS$VOLUME_HIGH <- as.numeric(arasens_ipd$OS$volume == "High")
arasens_ipd$PFS$VOLUME_HIGH <- as.numeric(arasens_ipd$PFS$volume == "High")
arasens_ipd$OS$M1 <- as.numeric(arasens_ipd$OS$meta == "De-novo")
arasens_ipd$PFS$M1 <- as.numeric(arasens_ipd$PFS$meta == "De-novo")
arasens_ipd$OS$RACE_WHITE <- as.numeric(arasens_ipd$OS$RACE == "WHITE")
arasens_ipd$PFS$RACE_WHITE <- as.numeric(arasens_ipd$PFS$RACE == "WHITE")
arasens_ipd$OS$VISC_MET <- 1 # FIX THIS
arasens_ipd$PFS$VISC_MET <- 1 # FIX THIS

mcspc_ipd <- list()
for(outcome_name in outcome_names) {
  mcspc_ipd[[outcome_name]] <-
    arasens_ipd[[outcome_name]][, c("treatment", "time", "event", "AGE", 
                                    "ECOG_0", "DOCELI", "GLEASON_M8", "VOLUME_HIGH",
                                    "M1", "RACE_WHITE", "VISC_MET")]
  
  colnames(mcspc_ipd[[outcome_name]])[
    is.element(colnames(mcspc_ipd[[outcome_name]]), c("treatment", "time", "event"))] <-
    c("trt", "eventtime", "status")
  
  mcspc_ipd[[outcome_name]] <- as.data.frame(mcspc_ipd[[outcome_name]])
  mcspc_ipd[[outcome_name]]$study <- 
    mcspc_ipd[[outcome_name]]$studyf <- "ARASENS"
  mcspc_ipd[[outcome_name]]$trt[mcspc_ipd[[outcome_name]]$trt == "Darolutamide+docetaxel arm"] <-
    "DAR+ADT+Doc"
  mcspc_ipd[[outcome_name]]$trt[mcspc_ipd[[outcome_name]]$trt == "Placebo+docetaxel arm"] <-
    "ADT+Doc"
  mcspc_ipd[[outcome_name]]$trtf <- mcspc_ipd[[outcome_name]]$trt
  
  # Reorder columns
  mcspc_ipd[[outcome_name]] <- mcspc_ipd[[outcome_name]][, c("study", "studyf", "trt", "trtf", "eventtime", "status",
                                                             "AGE", "ECOG_0", "DOCELI", "GLEASON_M8", "VOLUME_HIGH", "M1", "RACE_WHITE", "VISC_MET" )]
  
  # Remove the missing rows
  mcspc_ipd[[outcome_name]] <- 
    mcspc_ipd[[outcome_name]][mcspc_ipd[[outcome_name]]$trt != "", ]
  
  # Remove any with NA
  mcspc_ipd[[outcome_name]] <- 
    mcspc_ipd[[outcome_name]][!is.na(mcspc_ipd[[outcome_name]]$eventtime), ]
  # Remove any with NA
  mcspc_ipd[[outcome_name]] <- 
    mcspc_ipd[[outcome_name]][!is.na(mcspc_ipd[[outcome_name]]$status), ]
  
  # Set any missing covariates to the average value
  for(covariate_name in covariate_names) {
    mcspc_ipd[[outcome_name]][is.na(mcspc_ipd[[outcome_name]][, covariate_name]),
                              covariate_name] <-
      mean(mcspc_ipd[[outcome_name]][, covariate_name], na.rm = TRUE)
    
  }
}

mcspc_sens_ipd <- mcspc_ipd

#################################################################################
# ML-NMR analysis
#################################################################################
mcspc_net <- list()

# Network layout is specified manually
# In order of treatments object in the mcspc[[outcome_name]] object
net_x <- net_y <- list()
net_x$OS <- net_x$PFS <- 4 * c(-0.33, -0.33, 1, 0.33, -1, 0.33, -0.33, 1)
net_y$OS <- net_y$PFS <- c(0, -1, 0, 0, -1, 1, 1, 1)

# Times at which to predict hazards
timepoints <- list()
timepoints$OS <- seq(0, 60, length.out = 50) #c(12, 24, 36, 48, 60)
timepoints$PFS <- seq(0, 48, length.out = 50) #c(12, 24, 36, 48)

# Object to store the estimated models and summaries
mcspc_fit <- list()
mcspc_hazards <- list()
mcspc_ave_survival <- list()
arasens_survival <- list()
mcspc_arasens_effects <- mcspc_average_effects <- mcspc_aranote_effects <- list()

# Options for survival modelling
survival_models <- c("exponential", "weibull", "gompertz", "exponential-aft", "weibull-aft", "lognormal", "loglogistic", "gamma", "mspline")
do_survival_selection <- FALSE
survival_models_fit <- list()
survival_models_looic <- list()
survival_comparison_results <- matrix(NA, ncol = 8, nrow = length(survival_models) + 1)
rownames(survival_comparison_results) <- c(survival_models, "NPH")
colnames(survival_comparison_results) <- paste0(rep(c("OS", "PFS"), each = 4),
                                                c(" elpd_loo", " p_loo", " looic", " looic_raw"))
# List to store non-proportional hazards models
mcspc_fit_nph <- list()
# If we don't do survival model selection we'll use the below
selected_survival_model <- list()
selected_survival_model[["OS"]] <- selected_survival_model[["PFS"]] <- "mspline"
#names(selected_survival_model) <- outcome_names

# Option to run covariate assessment
# This runs the selected survival model with each covariate one-at-a-time
# Reports treatment effects and looic
do_covariate_assessment <- list()
do_covariate_assessment[["OS"]] <- do_covariate_assessment[["PFS"]] <- FALSE
covariate_fit <- list()
covariates_to_explore <- c("AGE", "ECOG_0", "GLEASON_M8", "VOLUME_HIGH")
covariate_assessment <- list()

for(outcome_name in outcome_names) {
  
  
  # Making the shared effect modifier assumption 
  # We group non-darolutamide triplets into a single class, and effects between
  # these treatments don't change.
  mcspc_ipd[[outcome_name]]$trtclass <- case_match(mcspc_ipd[[outcome_name]]$trtf,
                                                   "DAR+ADT+Doc" ~ "Darolutamide triplet",
                                                   treatment_names[treatment_names != "DAR+ADT+Doc"] ~ "Non-darolutamide")
  
  mcspc_agd[[outcome_name]]$trtclass <- case_match(mcspc_agd[[outcome_name]]$trtf,
                                                   "DAR+ADT+Doc" ~ "Darolutamide triplet",
                                                   treatment_names[treatment_names != "DAR+ADT+Doc"] ~ "Non-darolutamide")
  
  # Set up the IPD and AgD network using the Surv function for survival data
  
  mcspc_net[[outcome_name]] <- combine_network(
    set_ipd(mcspc_ipd[[outcome_name]],
            study = studyf,
            trt = trtf,
            trt_class = trtclass,
            Surv = Surv(eventtime, status)),
    set_agd_surv(mcspc_agd[[outcome_name]],
                 study = studyf,
                 trt = trtf,
                 trt_class = trtclass,
                 Surv = Surv(eventtime, status),
                 covariates = mcspc_agd_covs)
  )
  
  
  
  if(do_survival_selection) {
    survival_models_fit[[outcome_name]] <- list()
    survival_models_looic[[outcome_name]] <- list()
    
    # Only using n_int 16 rather than 64 for model selection/exploration
    mcspc_net[[outcome_name]] <- add_integration(mcspc_net[[outcome_name]],
                                                 AGE = distr(qgamma, 
                                                             mean = mean(mcspc_agd_covs_full[, "AGE"]), 
                                                             sd = sd(mcspc_agd_covs_full[, "AGE"])),
                                                 ECOG_0 = distr(qbern, ECOG_0),
                                                 GLEASON_M8 = distr(qbern, GLEASON_M8),
                                                 VOLUME_HIGH = distr(qbern, VOLUME_HIGH),
                                                 M1 = distr(qbern, M1),
                                                 RACE_WHITE = distr(qbern, RACE_WHITE),
                                                 n_int = 8)
    # error in gamma, gengamma
    system.time({
    for(survival_model in survival_models) {
      print(paste0(outcome_name, " : ", survival_model, " (", which(survival_model == survival_models), "/", length(survival_models), ")"))
      if(survival_model == "gengamma") {
        prior_aux_list <- list(sigma = half_normal(10), k = half_normal(10))
      } else {
        prior_aux_list <- half_normal(1)
      }
      survival_models_fit[[outcome_name]][[survival_model]] <- 
        nma(mcspc_net[[outcome_name]],
            regression = ~(AGE + ECOG_0 + GLEASON_M8 + VOLUME_HIGH)*.trt,
            likelihood = survival_model,
            prior_intercept = normal(0, 100),
            prior_trt = normal(0, 100),
            prior_reg = normal(0, 100),
            prior_aux = half_normal(1),
            chains = 2, iter = 1000,
            init = 0,
            QR = TRUE,
            # Reduce max_treedepth to speed up warm-up (default = 10)
            # Increase if max_treedepth warnings occur
            control = list(max_treedepth = 7))
      
      
      survival_models_looic[[outcome_name]][[survival_model]] <- 
        loo(survival_models_fit[[outcome_name]][[survival_model]])
      
      survival_comparison_results[survival_model, paste(outcome_name, "looic_raw")] <-
        survival_models_looic[[outcome_name]][[survival_model]]$estimates["looic", "Estimate"]
      for(quantity in c("elpd_loo", "p_loo", "looic")) {
        survival_comparison_results[survival_model, paste(outcome_name, quantity)] <- 
          format_results_se(survival_models_looic[[outcome_name]][[survival_model]]$estimates[quantity, ])
      }
    }
    })
    
    selected_survival_model[[outcome_name]] <- 
      survival_models[which.min(survival_comparison_results[survival_models, paste(outcome_name, "looic_raw")])]
    
    # Assess proportional hazards on the selected model
    mcspc_fit_nph[[outcome_name]] <- 
      nma(mcspc_net[[outcome_name]],
          regression = ~(AGE + ECOG_0 + GLEASON_M8 + VOLUME_HIGH)*.trt,
          likelihood = selected_survival_model[[outcome_name]],
          prior_intercept = normal(0, 100),
          prior_trt = normal(0, 100),
          prior_reg = normal(0, 100),
          prior_aux = half_normal(1),
          aux_by = c(.study, .trt),
          chains = 2, iter = 1000,
          QR = TRUE,
          # Reduce max_treedepth to speed up warm-up (default = 10)
          # Increase if max_treedepth warnings occur
          control = list(max_treedepth = 7))
    
    loo_nph_temp <- loo(mcspc_fit_nph[[outcome_name]])
    
    survival_comparison_results["NPH", paste(outcome_name, "looic_raw")] <-
      loo_nph_temp$estimates["looic", "Estimate"]
    for(quantity in c("elpd_loo", "p_loo", "looic")) {
      survival_comparison_results["NPH", paste(outcome_name, quantity)] <- 
        format_results_se(loo_nph_temp$estimates[quantity, ])
    }
    
    write.csv(survival_comparison_results, file = "results/ml-nmr_new/survival_comparison.csv")
  } # End of survival model comparison
  
  if(do_covariate_assessment[[outcome_name]]) {
    covariate_fit[[outcome_name]] <- list()
    
    covariate_assessment[[outcome_name]] <- matrix(NA, ncol = length(covariates_to_explore),
                                                   nrow = (length(treatment_names) + 4),
                                                   dimnames = list(c(treatment_names[treatment_names != "DAR+ADT+Doc"], 
                                                                     "coefficient", "coefficient:trt", "elpd_loo", "p_loo", "looic"),
                                                                   covariates_to_explore))
    
    
    # Only using n_int 16 rather than 64 for model selection/exploration
    mcspc_net[[outcome_name]] <- add_integration(mcspc_net[[outcome_name]],
                                                 AGE = distr(qgamma, 
                                                             mean = mean(mcspc_agd_covs_full[, "AGE"]), 
                                                             sd = sd(mcspc_agd_covs_full[, "AGE"])),
                                                 ECOG_0 = distr(qbern, ECOG_0),
                                                 GLEASON_M8 = distr(qbern, GLEASON_M8),
                                                 VOLUME_HIGH = distr(qbern, VOLUME_HIGH),
                                                 M1 = distr(qbern, M1),
                                                 RACE_WHITE = distr(qbern, RACE_WHITE),
                                                 n_int = 16)
    
    # Fit the regression one covariate at a time
    for(covariate_name in covariates_to_explore) {
      regression_formula <- as.formula(paste0("~(", covariate_name, ")*.trt"))
      
      covariate_fit[[outcome_name]][[covariate_name]] <- 
        nma(mcspc_net[[outcome_name]],
            regression = regression_formula,
            likelihood = selected_survival_model[[outcome_name]],
            prior_intercept = normal(0, 100),
            prior_trt = normal(0, 100),
            prior_reg = normal(0, 100),
            prior_aux = half_normal(1),
            chains = 2, iter = 1000,
            QR = TRUE,
            # Reduce max_treedepth to speed up warm-up (default = 10)
            # Increase if max_treedepth warnings occur
            control = list(max_treedepth = 7))
      
      # Calculate LOOIC
      loo_cov_temp <- loo(covariate_fit[[outcome_name]][[covariate_name]])
      
      
      # Calculate effects in average population
      average_effects_temp <- 
        relative_effects(covariate_fit[[outcome_name]][[covariate_name]], trt_ref = "DAR+ADT+Doc",
                         newdata = data.frame(
                           AGE = mean(mcspc_agd_covs_full$AGE),
                           ECOG_0 = mean(mcspc_agd_covs_full$ECOG_0),
                           GLEASON_M8 = mean(mcspc_agd_covs_full$GLEASON_M8),
                           VOLUME_HIGH = mean(mcspc_agd_covs_full$VOLUME_HIGH)),
                         probs = c(0.025, 0.975),
                         predictive_distribution = FALSE,
                         summary = TRUE)
      
      # Format results
      formatted_results_temp <-
        format_relative_effects(average_effects_temp,
                                hazard_scale = TRUE,
                                invert = TRUE)
      # Put the relative treatment effect in the correct place
      for(treatment_name in treatment_names[-which(treatment_names =="DAR+ADT+Doc")]) {
        covariate_assessment[[outcome_name]][treatment_name, covariate_name] <- 
          formatted_results_temp[paste0("d[New 1: ", treatment_name, "]"), ]
      }
      
      # Extract the effects relative to darolutamide triplet
      # And the looic
      for(quantity in c("elpd_loo", "p_loo", "looic")) {
        covariate_assessment[[outcome_name]][quantity, covariate_name] <-  
          format_results_se(loo_cov_temp$estimates[quantity, ])
      }
      
      

      # Extract the coefficients
        beta_temp <- as.data.frame(summary(covariate_fit[[outcome_name]][[covariate_name]], pars = "beta")) 
        
      covariate_assessment[[outcome_name]]["coefficient", covariate_name] <-
        format_results(beta_temp[1, c("mean", "2.5%", "97.5%")])
      covariate_assessment[[outcome_name]]["coefficient:trt", covariate_name] <-
        format_results(beta_temp[2, c("mean", "2.5%", "97.5%")])

      
      
      write.csv(covariate_assessment[[outcome_name]], file = paste0("results/ml-nmr_new/", outcome_name, "_covariate_assessment.csv"))
      save(covariate_fit, file = "covariate_fit.rda")
      
    }
    

    
    
  } # End of covariate assessment
  
  # Binary variables given Bernoulli distributions and age given gamma as it's skewed
  # Correlation can be specified but is otherwise estimated from the IPD
  # Omitted DOCELI and VISC_MET as IPD distribution can't be fit with ARASENS
  mcspc_net[[outcome_name]] <- add_integration(mcspc_net[[outcome_name]],
                                               AGE = distr(qgamma, 
                                                           mean = mean(mcspc_agd_covs_full[, "AGE"]), 
                                                           sd = sd(mcspc_agd_covs_full[, "AGE"])),
                                               ECOG_0 = distr(qbern, ECOG_0),
                                               #DOCELI = distr(qbern, DOCELI))
                                               GLEASON_M8 = distr(qbern, GLEASON_M8),
                                               VOLUME_HIGH = distr(qbern, VOLUME_HIGH),
                                               M1 = distr(qbern, M1),
                                               RACE_WHITE = distr(qbern, RACE_WHITE),
                                               n_int = 64)
  
  # Plot the network diagram
  plot(mcspc_net[[outcome_name]],
       weight_nodes = TRUE,
       weight_edges = TRUE,
       # Nudge treatment labels away from nodes
       nudge = 0.4, 
       # Manual layout
       layout = data.frame(x = net_x[[outcome_name]],
                           y = net_y[[outcome_name]])) +
    guides(edge_colour = guide_legend(override.aes = list(edge_width = 2))) +
    theme(legend.position = "bottom", legend.direction = "vertical", 
          plot.margin = margin(1, 3, 1, 2, unit = "cm"))
  
  ggsave(filename = paste0("results/ml-nmr_new/", outcome_name, "_evidence_network.png"),
         device = ragg::agg_png, # For antialiased output
         width = 6, height = 3.5, units = "in", dpi = 300, 
         # Rescale output
         scale = 1.2)
  
  # Plot the KM data without model fit
  ggplot() +
    geom_km(mcspc_net[[outcome_name]]) +
    facet_wrap(~.study) +
    labs(y = "Survival probability", x = "Time") +
    coord_cartesian(ylim = c(0, 1)) +
    theme_multinma() +
    theme(legend.position = "top", legend.box.spacing = unit(0, "lines"))
  
  ggsave(filename = paste0("results/ml-nmr_new/", outcome_name, "_raw_km.png"),
         device = ragg::agg_png, # For antialiased output
         width = 6, height = 4, units = "in", dpi = 300, 
         # Rescale output
         scale = 1.5)
  
  # Run the ML-NMR (assuming proportional hazards other than individual level
  # Could set iter= 1000 for model selection/exploration
  # variation)
  # With 1 core, n_int=64, 4 chains, took 27,458 seconds
  # With 1 core, n_int =4 and chains = 2 this took 2728 seconds
  # With 6 cores, n_int=16, chains=4, this took 6063
  mcspc_fit[[outcome_name]] <- nma(mcspc_net[[outcome_name]],
                                   regression = ~(AGE + ECOG_0 + GLEASON_M8 + VOLUME_HIGH)*.trt,
                                   likelihood = selected_survival_model[[outcome_name]],
                                   prior_intercept = normal(0, 100),
                                   prior_trt = normal(0, 100),
                                   prior_reg = normal(0, 100),
                                   prior_aux = half_normal(1),
                                   chains = 4, iter = 2000,
                                   QR = TRUE,
                                   # Reduce max_treedepth to speed up warm-up (default = 10)
                                   # Increase if max_treedepth warnings occur
                                   control = list(max_treedepth = 7))
  
  # Timings
  rstan::get_elapsed_time(as.stanfit(mcspc_fit[[outcome_name]])) %>% 
    cbind(total = rowSums(.))
  
  # Set up target populations for outputs
  pred_pop <- bind_rows(
    filter(mcspc_net[[outcome_name]]$ipd, study == "ARASENS") %>% select(study, AGE:VISC_MET),
    # For AgD studies, combine integration points across arms.
    # NOTE: We can do this because the sample sizes are balanced here. Otherwise
    # we should explicitly define the summary statistics for the combined
    # population and generate a new set of integration points for that.
    filter(mcspc_net[[outcome_name]]$agd_arm, study %in% c("ENZAMET", "ENZAMET_Doc")) %>% unnest_integration() %>% select(study, DOCELI:RACE_WHITE)
  )
  par_names <- list(ARASENS = "ARASENS", ENZAMET = "ENZAMET", `ENZAMET_Doc` = "ENZAMET_Doc")
  
  # Omit some outputs for now
  if(0) {
    mcspc_hazards[[outcome_name]] <- 
      predict(mcspc_fit[[outcome_name]], type = "hazard", level = "aggregate",
              times = timepoints[[outcome_name]])
    
    # Plot the hazards at population average (i.e., level is aggregate)
    plot(mcspc_hazards[[outcome_name]]) + facet_wrap("Study")
    ggsave(filename = paste0("results/ml-nmr_new/", outcome_name, "_hazards.pdf"),
           width = 6, height = 6, units = "in", scale = 1.2)
    
    # # Create dataset representing population at which we want to see hazards
    # refdat <- tibble(study = mcspc_net[[outcome_name]]$studies,
    #                  AGE = mcspc_agd_covs_full[mcspc_agd_covs_full$study == "ARASENS" &
    #                                              mcspc_agd_covs_full$trt == "DAR+ADT+Doc", "AGE"],
    #                  ECOG_0 = mcspc_agd_covs_full[mcspc_agd_covs_full$study == "ARASENS" &
    #                                                 mcspc_agd_covs_full$trt == "DAR+ADT+Doc", "ECOG_0"],
    #                  GLEASON_M8 = mcspc_agd_covs_full[mcspc_agd_covs_full$study == "ARASENS" &
    #                                                     mcspc_agd_covs_full$trt == "DAR+ADT+Doc", "ECOG_0"],
    #                  VOLUME_HIGH = mcspc_agd_covs_full[mcspc_agd_covs_full$study == "ARASENS" &
    #                                                      mcspc_agd_covs_full$trt == "DAR+ADT+Doc", "VOLUME_HIGH"])
    # 
    # # Specify distributions of baseline (intercept) and 
    # # auxiliary (spline coefficient) parameters
    # # Predict the hazards at evenly spaced times between the boundary knots
    # tdat <- purrr::imap_dfr(mcspc_net[[outcome_name]]$basis,
    #                         ~tibble(study = factor(.y, levels = mcspc_net[[outcome_name]]$studies),
    #                                 lower = attr(.x, "Boundary.knots")[1],
    #                                 upper = attr(.x, "Boundary.knots")[2],
    #                                 times = seq(lower, upper, length = 50)))
    # refdat <- left_join(refdat, tdat, by = "study")
    # studies <- as.list(setNames(nm = levels(mcspc_net[[outcome_name]]$studies)))
    
    plot(predict(mcspc_fit[[outcome_name]], type = "hazard", level = "aggregate",
                 times = timepoints[[outcome_name]], 
                 newdata = pred_pop, study = "study", baseline = par_names, aux = par_names))
    # 
  
  
  # Population average survival curves 
  mcspc_ave_survival[[outcome_name]] <- 
    predict(mcspc_fit[[outcome_name]], type = "survival",
            times = timepoints[[outcome_name]])
  
  # Plot with raw KM data included using geom_km
  plot(mcspc_ave_survival[[outcome_name]]) + 
    geom_km(mcspc_net[[outcome_name]]) +
    facet_wrap("Study") +
    theme(legend.position = "top", legend.box.spacing = unit(0, "lines"))
  
  ggsave(filename = paste0("results/ml-nmr_new/", outcome_name, "_ave_survival.pdf"),
         width = 6, height = 6, units = "in", scale = 1.4)
  

  # Survival curves at ARASENS population characteristics
    arasens_survival[[outcome_name]] <- predict(mcspc_fit[[outcome_name]], type = "hazard", level = "individual",
                                                newdata = refdat, study = study, times = times,
                                                baseline = studies, aux = studies)
    plot(arasens_survival[[outcome_name]])
  
  
  # Survival curves for target populations
  plot(predict(mcspc_fit[[outcome_name]], type = "survival", level = "aggregate",
               times = timepoints[[outcome_name]], 
               newdata = pred_pop, study = "study", baseline = par_names, aux = par_names))
  }
  
  results_matrix <- matrix(NA, nrow = length(treatment_names) - 1, 
                           ncol = 3)
  colnames(results_matrix) <- c("Average HR (95% CrI)", 
                                "ARASENS HR (95% CrI)",
                                "ARANOTE HR (95% CrI)")
  
  # Summarise the log hazard ratios in average population
  mcspc_average_effects[[outcome_name]] <- 
    relative_effects(mcspc_fit[[outcome_name]], trt_ref = "DAR+ADT+Doc",
                     newdata = data.frame(
                       AGE = mean(mcspc_agd_covs_full$AGE),
                       ECOG_0 = mean(mcspc_agd_covs_full$ECOG_0),
                       GLEASON_M8 = mean(mcspc_agd_covs_full$GLEASON_M8),
                       VOLUME_HIGH = mean(mcspc_agd_covs_full$VOLUME_HIGH)),
                     probs = c(0.025, 0.975),
                     predictive_distribution = FALSE,
                     summary = TRUE)
  
  # Format results. 
  results_matrix[, "Average HR (95% CrI)"] <-
    format_relative_effects(mcspc_average_effects[[outcome_name]],
                            hazard_scale = TRUE,
                            invert = TRUE)
  
  # Give names to the results matrix
  rownames(results_matrix) <- mcspc_average_effects[[outcome_name]]$summary$.trtb       
  
  # Log hazard ratios in ARASENS population
  mcspc_arasens_effects[[outcome_name]] <- 
    relative_effects(mcspc_fit[[outcome_name]], trt_ref = "DAR+ADT+Doc",
                     newdata = data.frame(
                       AGE = mcspc_agd_covs_full[mcspc_agd_covs_full$study == "ARASENS" &
                                                   mcspc_agd_covs_full$trt == "DAR+ADT+Doc", "AGE"],
                       ECOG_0 = mcspc_agd_covs_full[mcspc_agd_covs_full$study == "ARASENS" &
                                                      mcspc_agd_covs_full$trt == "DAR+ADT+Doc", "ECOG_0"],
                       GLEASON_M8 = mcspc_agd_covs_full[mcspc_agd_covs_full$study == "ARASENS" &
                                                          mcspc_agd_covs_full$trt == "DAR+ADT+Doc", "ECOG_0"],
                       VOLUME_HIGH = mcspc_agd_covs_full[mcspc_agd_covs_full$study == "ARASENS" &
                                                           mcspc_agd_covs_full$trt == "DAR+ADT+Doc", "VOLUME_HIGH"]),
                     probs = c(0.025, 0.975),
                     predictive_distribution = FALSE,
                     summary = TRUE)
  
  results_matrix[, "ARASENS HR (95% CrI)"] <-
    format_relative_effects(mcspc_arasens_effects[[outcome_name]],
                            hazard_scale = TRUE,
                            invert = TRUE)
  
  # Summarise the log hazard ratios in ARANOTE population
  mcspc_aranote_effects[[outcome_name]] <- 
    relative_effects(mcspc_fit[[outcome_name]], trt_ref = "DAR+ADT+Doc",
                     newdata = data.frame(
                       AGE = 69.5,
                       ECOG_0 = 0.502,
                       GLEASON_M8 = 0.682,
                       VOLUME_HIGH = 0.706),
                     probs = c(0.025, 0.975),
                     predictive_distribution = FALSE,
                     summary = TRUE)
  
  results_matrix[, "ARANOTE HR (95% CrI)"] <-
    format_relative_effects(mcspc_aranote_effects[[outcome_name]],
                            hazard_scale = TRUE,
                            invert = TRUE)
  
  write.csv(results_matrix, file = paste0("results/ml-nmr_new/", outcome_name, "_results_matrix.csv"))
  
  save(file = paste0(arasens_directory, "/Clifton Insight Analyses/mlnmr_new_objects.rda"),
       mcspc_net, mcspc_fit, mcspc_ave_survival)
  #mcspc_hazards,  mcspc_relative_effects
}

# Run the sensitivity analysis with all non-ADT treatments in a separate class
mcspc_sens_net <- list()
mcspc_sens_fit <- list()
mcspc_sens_arasens_effects <- list()
mcspc_sens_ave_effects <- list()
mcspc_sens_aranote_effects <- list()

for(outcome_name in outcome_names) {
  
  
  # Making the shared effect modifier assumption 
  # We group non-darolutamide triplets into a single class, and effects between
  # these treatments don't change.
  mcspc_sens_ipd[[outcome_name]]$trtclass <- case_match(mcspc_sens_ipd[[outcome_name]]$trtf,
                                                        "ADT" ~ "ADT",
                                                        treatment_names[treatment_names != "ADT"] ~ "Non-ADT")
  
  mcspc_sens_agd[[outcome_name]]$trtclass <- case_match(mcspc_sens_agd[[outcome_name]]$trtf,
                                                        "ADT" ~ "ADT",
                                                        treatment_names[treatment_names != "ADT"] ~ "Non-ADT")
  
  # Set up the IPD and AgD network using the Surv function for survival data
  mcspc_sens_net[[outcome_name]] <- combine_network(
    set_ipd(mcspc_sens_ipd[[outcome_name]],
            study = studyf,
            trt = trtf,
            trt_class = trtclass,
            Surv = Surv(eventtime, status)),
    set_agd_surv(mcspc_sens_agd[[outcome_name]],
                 study = studyf,
                 trt = trtf,
                 trt_class = trtclass,
                 Surv = Surv(eventtime, status),
                 covariates = mcspc_agd_covs)
  )
  
  
  # Binary variables given Bernoulli distributions and age given gamma as it's skewed
  # Correlation can be specified but is otherwise estimated from the IPD
  # Omitted DOCELI and VISC_MET as IPD distribution can't be fit with ARASENS
  mcspc_sens_net[[outcome_name]] <- add_integration(mcspc_sens_net[[outcome_name]],
                                                    AGE = distr(qgamma, 
                                                                mean = mean(mcspc_agd_covs_full[, "AGE"]), 
                                                                sd = sd(mcspc_agd_covs_full[, "AGE"])),
                                                    ECOG_0 = distr(qbern, ECOG_0),
                                                    #DOCELI = distr(qbern, DOCELI))
                                                    GLEASON_M8 = distr(qbern, GLEASON_M8),
                                                    VOLUME_HIGH = distr(qbern, VOLUME_HIGH),
                                                    M1 = distr(qbern, M1),
                                                    RACE_WHITE = distr(qbern, RACE_WHITE),
                                                    n_int = 32)
  
  mcspc_sens_fit[[outcome_name]] <- nma(mcspc_sens_net[[outcome_name]],
                                        regression = ~(AGE + ECOG_0 + GLEASON_M8 + VOLUME_HIGH)*.trt,
                                        likelihood = selected_survival_model[[outcome_name]],
                                        prior_intercept = normal(0, 100),
                                        prior_trt = normal(0, 100),
                                        prior_reg = normal(0, 100),
                                        prior_aux = half_normal(1),
                                        chains = 4, iter = 2000,
                                        QR = TRUE,
                                        # Reduce max_treedepth to speed up warm-up (default = 10)
                                        # Increase if max_treedepth warnings occur
                                        control = list(max_treedepth = 7))
  
  # Set up target populations for outputs
  pred_pop <- bind_rows(
    filter(mcspc_sens_net[[outcome_name]]$ipd, study == "ARASENS") %>% select(study, AGE:VISC_MET),
    # For AgD studies, combine integration points across arms.
    # NOTE: We can do this because the sample sizes are balanced here. Otherwise
    # we should explicitly define the summary statistics for the combined
    # population and generate a new set of integration points for that.
    filter(mcspc_sens_net[[outcome_name]]$agd_arm, study %in% c("ENZAMET", "ENZAMET_Doc")) %>% unnest_integration() %>% select(study, DOCELI:RACE_WHITE)
  )
  par_names <- list(ARASENS = "ARASENS", ENZAMET = "ENZAMET", `ENZAMET_Doc` = "ENZAMET_Doc")
  
  results_matrix_sens <- matrix(NA, nrow = length(treatment_names) - 1, 
                                ncol = 3)
  colnames(results_matrix_sens) <- c("Average HR (95% CrI)", 
                                     "ARASENS HR (95% CrI)",
                                     "ARANOTE HR (95% CrI)")
  
  # Summarise the log hazard ratios in average population
  mcspc_sens_average_effects[[outcome_name]] <- 
    relative_effects(mcspc_sens_fit[[outcome_name]], trt_ref = "DAR+ADT+Doc",
                     newdata = data.frame(
                       AGE = mean(mcspc_agd_covs_full$AGE),
                       ECOG_0 = mean(mcspc_agd_covs_full$ECOG_0),
                       GLEASON_M8 = mean(mcspc_agd_covs_full$GLEASON_M8),
                       VOLUME_HIGH = mean(mcspc_agd_covs_full$VOLUME_HIGH)),
                     probs = c(0.025, 0.975),
                     predictive_distribution = FALSE,
                     summary = TRUE)
  
  # Format results. 
  results_matrix_sens[, "Average HR (95% CrI)"] <-
    format_relative_effects(mcspc_sens_average_effects[[outcome_name]],
                            hazard_scale = TRUE,
                            invert = TRUE)
  
  # Give names to the results matrix
  rownames(results_matrix_sens) <- mcspc_sens_average_effects[[outcome_name]]$summary$.trtb       
  
  # Log hazard ratios in ARASENS population
  mcspc_sens_arasens_effects[[outcome_name]] <- 
    relative_effects(mcspc_sens_fit[[outcome_name]], trt_ref = "DAR+ADT+Doc",
                     newdata = data.frame(
                       AGE = mcspc_agd_covs_full[mcspc_agd_covs_full$study == "ARASENS" &
                                                   mcspc_agd_covs_full$trt == "DAR+ADT+Doc", "AGE"],
                       ECOG_0 = mcspc_agd_covs_full[mcspc_agd_covs_full$study == "ARASENS" &
                                                      mcspc_agd_covs_full$trt == "DAR+ADT+Doc", "ECOG_0"],
                       GLEASON_M8 = mcspc_agd_covs_full[mcspc_agd_covs_full$study == "ARASENS" &
                                                          mcspc_agd_covs_full$trt == "DAR+ADT+Doc", "ECOG_0"],
                       VOLUME_HIGH = mcspc_agd_covs_full[mcspc_agd_covs_full$study == "ARASENS" &
                                                           mcspc_agd_covs_full$trt == "DAR+ADT+Doc", "VOLUME_HIGH"]),
                     probs = c(0.025, 0.975),
                     predictive_distribution = FALSE,
                     summary = TRUE)
  
  results_matrix_sens[, "ARASENS HR (95% CrI)"] <-
    format_relative_effects(mcspc_sens_arasens_effects[[outcome_name]],
                            hazard_scale = TRUE,
                            invert = TRUE)
  
  # Summarise the log hazard ratios in ARANOTE population
  mcspc_sens_aranote_effects[[outcome_name]] <- 
    relative_effects(mcspc_sens_fit[[outcome_name]], trt_ref = "DAR+ADT+Doc",
                     newdata = data.frame(
                       AGE = 69.5,
                       ECOG_0 = 0.502,
                       GLEASON_M8 = 0.682,
                       VOLUME_HIGH = 0.706),
                     probs = c(0.025, 0.975),
                     predictive_distribution = FALSE,
                     summary = TRUE)
  
  results_matrix_sens[, "ARANOTE HR (95% CrI)"] <-
    format_relative_effects(mcspc_sens_aranote_effects[[outcome_name]],
                            hazard_scale = TRUE,
                            invert = TRUE)
  
  write.csv(results_matrix_sens, file = paste0("results/ml-nmr_new/", outcome_name, "_results_matrix_sens.csv"))
  
  save(file = paste0(arasens_directory, "/Clifton Insight Analyses/mlnmr_new_objects_sens.rda"),
       mcspc_sens_net, mcspc_sens_fit)
  
}