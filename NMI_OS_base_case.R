# Darolutamide in mCSPC NMI #
# Outcome: OS
# Target population: ARANOTE (base case)

# 1. Packages used --------------------------------------------------------
pkgs <- c("tidyverse", "readxl", "here", "haven", "multinma", "R2OpenBUGS",
          "igraph", "doParallel", "doSNOW", "tictoc", "tibble", "broom.mixed", 
          "reshape2", "lemon", "kableExtra", "plotly", "meta") # Do we need all these packages?
lapply(pkgs, library, character.only = T)
set.seed(18022015)
source(here("NMI", "code", "NMI_functions.R")) 
source(here("NMI", "code", "utils.R")) 

# 2. Import data ---------------------------------------------------------
#2.1 Individual-level data
aranote_ipd <- read.csv(here("data", "ARANOTE_OS.csv")) #IPD from ARANOTE created in 'ARANOTE data preparation.R'
IPD_EM_cols <- AgD_EM_cols <- c('ECOG_0', 'GLEASON_L7', 'AGE_L65', 'VOLUME_HIGH') # Effect modifiers columns

aranote_ipd$ECOG_0 <- if_else(aranote_ipd$ecog == 0, 1, 0) # Create new columns with dummy variables for subgroups
aranote_ipd$GLEASON_L7 <- if_else(aranote_ipd$gleason == "<8", 1, 0)
aranote_ipd$AGE_L65 <- if_else(aranote_ipd$age65 == "<65", 1, 0)
aranote_ipd$VOLUME_HIGH <- if_else(aranote_ipd$volume == 1, 1, 0)
aranote_ipd <- select(aranote_ipd, time, event, treatment, IPD_EM_cols) # Select only column of interest to merge with arasens

# Treatment effect from ARANOTE 0.813 (0.591, 1.118) based on NUBEQA PT ARANOTE July 15 2024 TLR _ FINAL_v1.2.pdf
aranote_te <- log(0.813)
aranote_se <- (1.118-0.591)/3.92

arasens_ipd <- read.csv(here("NMI", "arasens_ipd_os.csv")) #IPD from ARASENS to improve correlation estimation
arasens_ipd <- select(arasens_ipd, time, event, treatment, IPD_EM_cols) # Select only column of interest to merge with aranote

load(file = here("data", "CHAARTED Trial", "chaarted_ipd_maic.R")) # load CHAARTED IPD with covariates created in prostate_prepare_chaarted.R
chaarted_ipd <- chaarted_ipd$OS
chaarted_ipd <- select(chaarted_ipd, time, event, treatment, ecog, gleason, age, dz_extent)
chaarted_ipd$ECOG_0 <- if_else(chaarted_ipd$ecog == 0, 1, 0) # Create new columns with dummy variables for subgroups
chaarted_ipd$GLEASON_L7 <- if_else(chaarted_ipd$gleason == "<8", 1, 0)
chaarted_ipd$AGE_L65 <- if_else(chaarted_ipd$age == "<70", 1, 0)
chaarted_ipd$VOLUME_HIGH <- if_else(chaarted_ipd$dz_extent == "High", 1, 0)
chaarted_ipd <- select(chaarted_ipd, time, event, treatment, IPD_EM_cols) # Select only column of interest to merge with aranote

# ECOG data obtained from https://pubmed.ncbi.nlm.nih.gov/29384722/ (supplement)
ecog_0 <- c(0.75,	0.58,	0.97) # HR in patients with ECOG 0
ecog_1 <- c(0.58,	0.41,	0.83) # HR in patients with ECOG ≥1
ecog_0_n <- 549/790 # Sample size in patients with ECOG 0
ecog_1_n <- 241/790 # Sample size in patients with ECOG ≥1
log(ecog_0[1])      # log of treatment effect
(ecog_0[3]-ecog_0[2])/3.92 # standard error
log(ecog_1[1])
(ecog_1[3]-ecog_1[2])/3.92

# Gleason data obtained from https://pubmed.ncbi.nlm.nih.gov/29384722/ (supplement)
gleason_7 <- c(0.66,	0.42,	1.03) # HR in patients with gleason <8
gleason_8 <- c(0.68,	0.53,	0.87) # HR in patients with gleason ≥8
gleason_7_n <- 221/790 # Sample size in patients with gleason 0
gleason_8_n <- 48/790 # Sample size in patients with gleason ≥1
log(gleason_7[1])      # log of treatment effect
(gleason_7[3]-gleason_7[2])/3.92 # standard error
log(gleason_8[1])
(gleason_8[3]-gleason_8[2])/3.92

IPD <- bind_rows(aranote_ipd, arasens_ipd, chaarted_ipd) # ARASENS and ARANOTE merged, including effect modifiers
IPD <- drop_na(IPD) # Delete missing values, otherwise imputation doesn't work
cor_ipd <- as_tibble(cor(IPD[, c(4:7)]))
write_csv(cor_ipd, here("NMI", "Results", "OS_correlation_ipd_with_chaarted.csv"))

#2.1 Aggregate-level data
AgD <- read.csv(here('NMI', 'Data', 'OS_data_updated_with_arasens.csv')) #AgD for NMI (updated with new subgroups by David Aceituno)
Study_col <- 'Study_n' #study column name
samp_sizes <- AgD[match(unique(na.omit(AgD$n)), AgD$n),]$n #sample sizes for AgD studies
Trt_cols <- AgD_Trt_cols <- c('Trt1', 'Trt2') #AgD treatment column names
TE_col <- 'TE' #AgD treatment effect estimate column name
SE_col <- 'SE' #AgD standard error column name
IPD_Trt_col <- 'Tr' #IPD treatment column name

# 3. Network meta-interpolation--------------------------
#3.1 Step 1: data enrichement using BLUP
imputed <- BLUP_impute(IPD, AgD, AgD_EM_cols, IPD_EM_cols, Study_col, 
            samp_sizes, AgD_Trt_cols, TE_col, IPD_Trt_col, SE_col)

#3.2 Step 2: Estimation of the treatment effect and its associated variance at specified EM values)
studies_to_impute <- read.csv(here("NMI","Data", "OS_studies_to_impute.csv")) # Use ARANOTE as target population

#value of effect modifiers to be used in ITC (calculated directly from ARANOTE IPD)
x1 <- sum(aranote_ipd$ECOG_0 == 1)/nrow(aranote_ipd)
x2 <- sum(aranote_ipd$GLEASON_L7 == 1, na.rm = TRUE)/nrow(aranote_ipd)
x3 <- sum(aranote_ipd$AGE_L65 == 1, na.rm = TRUE)/nrow(aranote_ipd)
x4 <- sum(aranote_ipd$VOLUME_HIGH == 1)/nrow(aranote_ipd)

x_vect = c(x1, x2, x3, x4)

# This function include the previous BLUP function inside. It also estimate the treatment effect at specified EM values
NMI_object <- NMI_interpolation(IPD, AgD, x_vect, AgD_EM_cols, IPD_EM_cols, 
                  Study_col, samp_sizes, AgD_Trt_cols, TE_col, 
                  SE_col, IPD_Trt_col)

# The output of this function includes an imputed dataset, treatment effects by study at desired EM values and diagnostics 
NMI_object$Imputed
NMI_object$Final
NMI_object$Diagnostics

original_te <- filter(AgD, !is.na(n))
predicted_te <- select(NMI_object$Final, TE, se)
comparison_te <- bind_cols(original_te$Study, round(original_te[,c("TE", "SE")], 4), round(predicted_te, 4))
names(comparison_te) <- c("Study","Original LHR","Original SE","Predicted LHR","Predicted SE")
write_csv(comparison_te, here("NMI", "Results", "OS_original_vs_predicted_TE_v2.csv"))

#3.3 Step 3: Application of standard NMA on the interpolated treatment effect and variance.

#Estimate Vaishampayan TE from ARCHES
#value of effect modifiers in Vaishampayan # This was added to borrow information from another trial

x1 = studies_to_impute[studies_to_impute$Study == 'Vaishampayan', 'ECOG_0']
x2 = studies_to_impute[studies_to_impute$Study == 'Vaishampayan', 'GLEASON_L7']
x3 = studies_to_impute[studies_to_impute$Study == 'Vaishampayan', 'AGE_L65']
x4 = studies_to_impute[studies_to_impute$Study == 'Vaishampayan', 'VOLUME_HIGH']

x_vect = c(x1, x2, x3, x4)

imputed = BLUP_impute(IPD, AgD, AgD_EM_cols, IPD_EM_cols, Study_col, 
                      samp_sizes, AgD_Trt_cols, TE_col, IPD_Trt_col, SE_col)
imputed$Study = AgD$Study
dat = imputed %>% filter(Study == 'ARCHES') 
X = apply(as.matrix(dat[, AgD_EM_cols]), 2, as.numeric)
x_orig = X[1,]
m = ncol(X)
M2 = as.matrix(cbind(1, X^2, 2*X, 
                     apply(combn(1:m, 2), 2, 
                           function(u){2*X[,u[1]]*X[,u[2]]})))

beta_hat = lm(as.numeric(dat$TE) ~ X)$coef
sigma_hat = c(t(M2)%*%solve(M2%*%t(M2), (dat[,SE_col])^2))

sigma_hat_vec_to_mat = function(sigma_hat){
  M = length(sigma_hat)
  K = (sqrt(1 + 8*M) - 1)/2
  
  C = matrix(0, nrow = K, ncol = K)
  t = K + 1
  for(i in 1:(K-1)){
    C[i,(i+1):K] = sigma_hat[t:(t + K - i - 1)]
    t = t + K - i
  }
  
  C = C + t(C)
  diag(C) = sigma_hat[1:K]
  
  return(C)
}

#C = sigma_hat_vec_to_mat(sigma_hat)
#eigen(C)$values

#  x_vect_star = c(1, x_vect)
#  SE = c(sqrt(t(x_vect_star)%*%C%*%x_vect_star))
#  SE


C = sigma_hat_vec_to_mat(sigma_hat)
eigen_vals = eigen(C)$values
lambda_min = min(eigen_vals)

if(lambda_min <= 0){
  diag(C) = diag(C) - lambda_min + 1e-6
}

x_vect_star = c(1, x_vect)

TE = beta_hat%*%x_vect_star
se = sqrt(t(x_vect_star) %*% C %*% x_vect_star)
# TE = beta_hat%*%c(1, x_vect)
# se = sqrt(t(sigma_hat)%*%u)


#Difference between two study baselines

delta_hat <- studies_to_impute[studies_to_impute$Study == 'Vaishampayan', 'TE'] - TE

#TE for Vaishampayan at the conditions of the eventual NMA

x1 <- sum(aranote_ipd$ECOG_0 == 1)/nrow(aranote_ipd)
x2 <- sum(aranote_ipd$GLEASON_L7 == 1, na.rm = TRUE)/nrow(aranote_ipd)
x3 <- sum(aranote_ipd$AGE_L65 == 1, na.rm = TRUE)/nrow(aranote_ipd)
x4 <- sum(aranote_ipd$VOLUME_HIGH == 1)/nrow(aranote_ipd)


x_vect = c(x1, x2, x3, x4)
beta_hat["(Intercept)"] <- beta_hat["(Intercept)"] + delta_hat
TE_vaish = beta_hat %*%c(1, x_vect)

x_vect_star = c(1, x_vect)


se_nma = sqrt(t(x_vect_star) %*% C %*% x_vect_star)

SE_vaish = (se_nma/se) * studies_to_impute[1, 'SE']
new_row <- nrow(NMI_object$Final)+1
NMI_object$Final[new_row, "Study"] <- nrow(NMI_object$Final)+1
NMI_object$Final[new_row, "Trt1"] <- studies_to_impute[studies_to_impute$Study == 'Vaishampayan', 'Trt1']
NMI_object$Final[new_row, "Trt2"] <- studies_to_impute[studies_to_impute$Study == 'Vaishampayan', 'Trt2']
NMI_object$Final[new_row, "TE"] <- TE_vaish
NMI_object$Final[new_row, "se"] <- SE_vaish

NMI_object$Final$study_name <- c(unique(AgD$Study), "Vaishampayan")

studies_not_changed <- c("STAMPEDE C vs G", "ENZAMET", "ARANOTE") # deleted PEACE-1 because now we have subgroup data and ENZ+DOC+ADT because isn't connected anymore
OS_agg_data <- read.csv(here("NMI", "Data", "mCSPC OS updated.csv")) # Previous network meta-regression results (updated with ARANOTE)
treatment_numbers <- read.csv(here("NMI", "Data", "treatment numbers updated.csv"))
new_row <- nrow(NMI_object$Final)

for (i in 1:length(studies_not_changed)) {
  
  NMI_object$Final[new_row+i, "study_name"] = studies_not_changed[i]
  NMI_object$Final[new_row+i, "Trt1"] = OS_agg_data[OS_agg_data$X.ID == studies_not_changed[i], "t1"]
  NMI_object$Final[new_row+i, "Trt2"] = OS_agg_data[OS_agg_data$X.ID == studies_not_changed[i], "t2"]
  NMI_object$Final[new_row+i, "TE"] = OS_agg_data[OS_agg_data$X.ID == studies_not_changed[i], "y"]
  NMI_object$Final[new_row+i, "se"] = OS_agg_data[OS_agg_data$X.ID == studies_not_changed[i], "se"]
}

NMI_object$Final <- merge(NMI_object$Final, treatment_numbers, by="study_name")

write.csv(NMI_object$Final, "NMI/Imputed Data/OS_data_imputed_ARANOTE.csv")

# 4. Use imputed data in the NMA models------------------------------------------
OS_imputed <- read.csv(here("NMI", "Imputed Data", "OS_data_imputed_ARANOTE.csv"))
OS_imputed$t1 <- c(2, 1, 3, 7, 3, 7, 4, 6, 7, 4, 4, 5, 3) # change treatment indices so DAR triplet is reference

# 5. NMA models --------------------------------------------------------
# Normal likelihood, identity link, fixed effects
# With adjustment for Multiarm multistage (MAMS) design
# Assumes only a single MAMS trial with na_mams arms

# Number of MCMC chains and samples
n_chains <- 3
num_sims <- 10000 * n_chains 
burn_in <- 10000 * n_chains	

# Define the bugs data 
# also, to get the correct number of dimensions is good to use a "comparator" arm with 0 for the lhr and the se
ns <- nrow(OS_imputed)
t  <- array(c(OS_imputed$t1, OS_imputed$t2), dim = c(ns, 2)) 
nt <- max(t) 
y  <- array(c(rep(0, ns), OS_imputed$TE), dim = c(ns, 2))
se <- array(c(rep(0, ns), OS_imputed$se), dim = c(ns, 2))
study_names <- gsub("#", "", OS_imputed$study_name)
rownames(t) <- rownames(y) <- rownames(se) <- study_names

# Separate out the STAMPEDE/MAMS studies
mams_indices <- grep("STAMPEDE", rownames(y))

# Data for the MAMS trials
y_mams <- y[mams_indices, 2]
se_mams <- se[mams_indices, ] # Not used as SE may not be correct and need covariance
t_mams <- t[mams_indices, ]
na_mams <- length(y_mams)

# Data for the non-MAMS trials (analysed as usual)
y_non_mams <- y[-mams_indices, ]
se_non_mams <- se[-mams_indices, ]
t_non_mams <- t[-mams_indices, ] 
ns_non_mams <- dim(y_non_mams)[1]

# Covariance matrix for the MAMS trial
# From Section 4.3 of the SAP (FFS/PFS can be found there as well)
var_mams <- matrix(c(0.01179151, 0.00203, 0.01113,
                     0.00203, 0.00665433, 0.00085,
                     0.01113, 0.00085, 0.0318173),
                   nrow = 3)

# Inverse of covariance matrix is precision
prec_mams <- solve(var_mams)

# Bugs data for unadjusted model
bugs_data <- list(
  y = y,
  se = se,
  t = t,
  ns = ns, 
  nt = nt)

# Bugs data for adjusted model
bugs_data_mams <- list(
  y = y_non_mams,
  se = se_non_mams,
  t = t_non_mams,
  ns = ns_non_mams, 
  y_mams = y_mams,
  prec_mams = prec_mams,
  t_mams = t_mams,
  na_mams = na_mams,
  nt = nt)

# Create initial values for MCMC simulation 
# initial values according to the number of parameters
# These are the same for both adjusted and unadjusted models
inits1 <- list(d=c( NA, rep(0, nt - 1)))
inits2 <- list(d=c( NA, rep(-1, nt - 1)))
inits3 <- list(d=c( NA, rep(2, nt - 1)))
bugs_inits <- list(inits1, inits2, inits3)

# Use wine to call OpenBUGS from a Mac
WINE="/usr/local/bin/wine"
WINEPATH="/usr/local/bin/winepath"
OpenBUGS.pgm="/Users/daceituno/.wine/drive_c/OpenBUGS323/OpenBUGS.exe"

# 6. Call OpenBUGS --------------------------------------------------------
bugs_object_fe_mams <- bugs(data = bugs_data_mams, inits = bugs_inits,
                            parameters.to.save = c("d", "totresdev", "rk", "best", "prob"),
                            model = model_mams_adjust_fe, clearWD = TRUE, 
                            summary.only = FALSE,
                            OpenBUGS.pgm=OpenBUGS.pgm, # This line to run on a Mac
                            WINE=WINE,                 # This line to run on a Mac
                            WINEPATH=WINEPATH,         # This line to run on a Mac
                            useWINE=TRUE,              # This line to run on a Mac
                            n.iter = (num_sims + burn_in), n.burnin = burn_in,
                            n.chains = n_chains, bugs.seed = 1, debug = FALSE, save.history = TRUE)


bugs_object_re_mams <- bugs(data = bugs_data_mams, inits = bugs_inits,
                            parameters.to.save = c("d", "totresdev", "rk", "best", "prob"),
                            model = model_mams_adjust_re, clearWD = TRUE, 
                            summary.only = FALSE,
                            OpenBUGS.pgm=OpenBUGS.pgm,
                            WINE=WINE,
                            WINEPATH=WINEPATH,
                            useWINE=TRUE,
                            n.iter = (num_sims + burn_in), n.burnin = burn_in,
                            n.chains = n_chains, bugs.seed = 1, debug = FALSE)

# Models fit (uses extract_fit function from utils.R)
extract_fit(bugs_object_fe_mams)
extract_fit(bugs_object_re_mams)

# 7. Results (fixed effects)-----------------------------------------------------
#7.1 Treatment effects against DAR+ADT
releff_FE <- bugs_object_fe_mams$summary[grep("d", rownames(bugs_object_fe_mams$summary)), c("mean", "2.5%", "97.5%")] #EK added
releff_FE <- releff_FE[1:8,]

# Transforming to hazard ratios
LHR_array <- as.array(releff_FE)
HR_array <- exp(LHR_array)
HR_FE <- as_tibble(apply(HR_array, 1, format_results)) # format results in CI style

# Include a column with treatment names
HR_FE$treatment <-  c("DAR+ADT", "ENZ+ADT", "ABI+ADT", "APA+ADT", "ABI+DOC+ADT", "DOC+ADT", "ADT", "SNA+ADT")
write_csv(HR_FE, here("NMI", "Results", "OS_relative_effects_FE_v3.csv"))

#7.2 Forest plot
HR_tibble <- as_tibble(HR_array)
jpeg(here("NMI", "Results", "OS_forest_FE_v3.jpg"), width = 900, height = 700, res = 120)
forest_FE <- forest(x=(HR_tibble$mean),          #change to x=(data$mean) if you would like to plot posterior means
                    ci.lb=(HR_tibble$`2.5%`),
                    ci.ub=(HR_tibble$`97.5%`),
                    digits=2,
                    showweights=FALSE,
                    slab=HR_FE$treatment,
                    psize=1,                      #size of observed effect
                    xlab="HR",                   #change to HR, MD, etc. if appropriate
                    cex=1.2,
                    refline=1)
dev.off()

#7.3 Cross-tables
treatments <- c("DAR+DOC+ADT", HR_FE$treatment)
cross_meandiff_fe <- cross_effect(bugs_object = bugs_object_fe_mams, t_names = treatments, med = TRUE, exp = TRUE)
write.csv(x = cross_meandiff_fe, file = here("NMI", "Results", "OS_cross_table_FE_ARANOTE_v3.csv"))

#7.4 SUCRA and ranks
sucra_FE <- calculate_sucra(bugs_object = bugs_object_fe_mams, bugs_data = bugs_data_mams, t_names = treatments)

#Ranks
FE_ranks <- bugs_object_fe_mams$summary[grep("rk", rownames(bugs_object_fe_mams$summary)), c("mean", "2.5%", "97.5%")]
FE_ranks <- as_tibble(apply(FE_ranks, 1, format_results))

#Table of ranks (NMI)
FE_rank_table <- tibble(treatments, FE_ranks, round(sucra_FE, digits = 2))
colnames(FE_rank_table) <- c("Treatment", "Mean rank (95% CI) - NMI",  "SUCRA - NMI")
FE_rank_export <- FE_rank_table[order(FE_rank_table$`SUCRA - NMI`, decreasing = TRUE), ] 
write.csv(FE_rank_export, here("NMI", "Results", "OS_rank_table_FE_ARANOTE_v3.csv"))

# 8. Results (Random effects)-----------------------------------------------------
#8.1 Treatment effects against DAR+ADT
releff_RE <- bugs_object_re_mams$summary[grep("d", rownames(bugs_object_re_mams$summary)), c("mean", "2.5%", "97.5%")] #EK added
releff_RE <- releff_RE[1:8,]

# Transforming to hazard ratios
LHR_array <- as.array(releff_RE)
HR_array <- exp(LHR_array)
HR_RE <- as_tibble(apply(HR_array, 1, format_results)) # format results in CI style

# Include a column with treatment names
HR_RE$treatment <-  c("DAR+DOC+ADT", "ENZ+ADT", "ABI+ADT", "APA+ADT", "ABI+DOC+ADT", "DOC+ADT", "ADT", "SNA+ADT")
write_csv(HR_RE, here("NMI", "Results", "OS_relative_effects_RE.csv"))

#8.2 Forest plot
HR_tibble_ <- as_tibble(HR_array)

jpeg(here("NMI", "Results", "OS_forest_RE.jpg"), width = 900, height = 700, res = 120)
forest_RE <- forest(x=(HR_tibble$mean),          #change to x=(data$mean) if you would like to plot posterior means
                    ci.lb=(HR_tibble$`2.5%`),
                    ci.ub=(HR_tibble$`97.5%`),
                    digits=2,
                    showweights=FALSE,
                    slab=HR_FE$treatment,
                    psize=1,                      #size of observed effect
                    xlab="HR",                   #change to HR, MD, etc. if appropriate
                    cex=1.2,
                    refline=1)
dev.off()

#8.3 Cross-tables
cross_meandiff_re <- cross_effect(bugs_object = bugs_object_re_mams, t_names = treatments, med = TRUE, exp = TRUE)
write.csv(x = cross_meandiff_re, file = here("NMI", "Results", "OS_cross_meandiff_RE_ARANOTE.csv"))

#8.4 SUCRA and ranks
sucra_RE <- calculate_sucra(bugs_object = bugs_object_re_mams, bugs_data = bugs_data_mams, t_names = treatments)

#Ranks
RE_ranks <- bugs_object_re_mams$summary[grep("rk", rownames(bugs_object_re_mams$summary)), c("mean", "2.5%", "97.5%")]
RE_ranks <- as_tibble(apply(RE_ranks, 1, format_results))

#Table of ranks (NMI)
RE_rank_table <- tibble(treatments, RE_ranks, round(sucra_RE, digits = 2))
colnames(RE_rank_table) <- c("Treatment", "Mean rank (95% CI) - NMI",  "SUCRA - NMI")
RE_rank_export <- RE_rank_table[order(RE_rank_table$`SUCRA - NMI`, decreasing = TRUE), ] 
write.csv(RE_rank_export, here("NMI", "Results", "OS_rank_table_RE_ARANOTE.csv"))
