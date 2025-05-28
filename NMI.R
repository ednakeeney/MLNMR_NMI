
#NMI
# Code associated with manuscript titled: "Population-Adjusted Network Meta-Analyses (NMA) Provide New Insights 
#into the Efficacy of Treatment Alternatives for Metastatic Castration-Sensitive Prostate Cancer (mCSPC)"

library(R2OpenBUGS)
library(dplyr)
library(doParallel)
library(doSNOW)
library(tictoc)
library(tibble)
library(multinma)
library(tidyr)
library(ggplot2)
library(broom)
library(broom.mixed)
library(reshape2)
library(lemon)
library(kableExtra)
library(purrr)
library(plotly)


IPD = read.csv("C:/Users/EdnaKeeney/OneDrive - Thoms Statistical Consultants Limited/jrsm1608-sup-0001-supinfo/jrsm1608-sup-0001-supinfo/Data/arasens_ipd_os.csv") #IPD 
AgD = read.csv('NMI/Data/OS_data.csv') #AgD for NMI
IPD_EM_cols = AgD_EM_cols = c('ECOG_0', 'GLEASON_L7', 'AGE_L65', 'VOLUME_HIGH')
Study_col = 'Study_n' #study column name
samp_sizes = AgD[match(unique(na.omit(AgD$n)), AgD$n),]$n #sample sizes for AgD studies
Trt_cols = AgD_Trt_cols = c('Trt1', 'Trt2') #AgD treatment column names
TE_col = 'TE' #AgD treatment effect estimate column name
SE_col = 'SE' #AgD standard error column name
IPD_Trt_col = 'Tr' #IPD treatment column name


#Function to impute missing proportions of EMs (Step 1)
BLUP_impute = function(IPD, AgD, AgD_EM_cols, IPD_EM_cols, Study_col, 
                       samp_sizes, AgD_Trt_cols, TE_col, IPD_Trt_col, SE_col){
  rho = cor(IPD[, IPD_EM_cols]) #estimate correlation between effect modifiers
  n_studs = length(unique(AgD[,Study_col]))
  studies = 1:n_studs
  p = length(AgD_EM_cols)
  
  single_BLUP_impute = function(study){
    n = samp_sizes[study]
    Xbar = AgD[AgD[, Study_col] == study, AgD_EM_cols][1,] #Means of effect modifiers
    Sbar = sqrt(Xbar*(1 - Xbar)/n) #Variance of effect modifiers
    
    missing_mat = function(i){
      a = rep(NA, 2*p + 1)
      a[1] = Xbar[i]
      a[2*i] = 1
      a[2*i + 1] = 0
      
      return(a)
    }
    
    X = sapply(1:p, missing_mat)
    Y = X
    for(i in 1:(2*p)){
      ind = ceiling(i/2)
      for(j in 1:p){
       # Y[-1,][i,j] = Sbar[j]/Sbar[ind]*rho[ind,j]*(X[-1,][i,ind] - Xbar[ind]) + Xbar[j]  #Equation 13 in paper
        temp = Sbar[j]/Sbar[ind]*rho[ind,j]*(X[-1,][i,ind] - Xbar[ind]) + Xbar[j]
        temp = unlist(ifelse(Xbar[j] == 1, 1, ifelse(Xbar[j] == 0, 0, temp)))
        Y[-1,][i,j] = max(min(temp, 1), 0)
      }
    }
    
    return(Y)
  }
  
  
  imputed = do.call(rbind, lapply(as.list(studies), single_BLUP_impute))
  out = cbind(AgD[, c(Study_col, AgD_Trt_cols)], 
              imputed, AgD[, c(TE_col, SE_col)])
  names(out)[p:(p+3)] = AgD_EM_cols 
  names(out)[1] = 'Study'

  return(out)
}

BLUP_impute(IPD, AgD, AgD_EM_cols, IPD_EM_cols, Study_col, 
                       samp_sizes, AgD_Trt_cols, TE_col, IPD_Trt_col, SE_col)


#' NMI interpolation pre-NMA
#'
#' Interpolating treatment effect estimates and standard errors at new 
#' effect modifier values (Step 2)
#'

studies_to_impute <- read.csv("NMI/Data/OS_studies_to_impute.csv")

#value of effect modifiers to be used in ITC
x1 = studies_to_impute[studies_to_impute$Study == 'ARASENS', 'ECOG_0']
x2 = studies_to_impute[studies_to_impute$Study == 'ARASENS', 'GLEASON_L7']
x3 = studies_to_impute[studies_to_impute$Study == 'ARASENS', 'AGE_L65']
x4 = studies_to_impute[studies_to_impute$Study == 'ARASENS', 'VOLUME_HIGH']

x_vect = c(x1, x2, x3, x4)

NMI_interpolation = function(IPD, AgD, x_vect, AgD_EM_cols, IPD_EM_cols, 
                             Study_col, samp_sizes, AgD_Trt_cols, TE_col, 
                             SE_col, IPD_Trt_col){
  imputed = BLUP_impute(IPD, AgD, AgD_EM_cols, IPD_EM_cols, Study_col, 
                        samp_sizes, AgD_Trt_cols, TE_col, IPD_Trt_col, SE_col) #From step 1
  studies = unique(imputed$Study)
  
  single_study_interpolation = function(study){
    dat = imputed %>% filter(Study == study)
    X = apply(as.matrix(dat[, AgD_EM_cols]), 2, as.numeric)
    x_orig = X[1,]
    m = ncol(X)
    M2 = as.matrix(cbind(1, X^2, 2*X, 
                         apply(combn(1:m, 2), 2, 
                               function(u){2*X[,u[1]]*X[,u[2]]})))  #Equation 17 in paper
    
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
    
    
    C = sigma_hat_vec_to_mat(sigma_hat)
    eigen_vals = eigen(C)$values
    lambda_min = min(eigen_vals)
    
    if(lambda_min <= 0){
      diag(C) = diag(C) - lambda_min + 1e-6
    }
    
    x_vect_star = c(1, x_vect)
    
    TE = beta_hat%*%x_vect_star
    se = sqrt(t(x_vect_star) %*% C %*% x_vect_star)
    
    TE_orig = dat[, TE_col]
    TE_pred = cbind(1, X)%*%beta_hat
    
    
    NMI_out = data.frame(Study = study, 
                         Trt1 = unique(dat[,AgD_Trt_cols[1]]),
                         Trt2 = unique(dat[,AgD_Trt_cols[2]]),
                         x = rbind(x_vect),
                         TE = TE, se = se) 
    
    colnames(NMI_out) = gsub('\\.', '', colnames(NMI_out))
    
    Diag_out = data.frame(Study = study,
                          X,
                          TE_orig = TE_orig,
                          TE_pred = TE_pred)
    
    list(NMI_out = NMI_out, Diag_out = Diag_out)
  }
  
  out = lapply(studies, single_study_interpolation)
  
  Final = do.call(rbind, lapply(out, `[[`, 1))
  Diagnostics = do.call(rbind, lapply(out, `[[`, 2))
  
  return(list(Imputed = imputed, Final = Final, Diagnostics = Diagnostics))
}

NMI_interpolation(IPD, AgD, x_vect, AgD_EM_cols, IPD_EM_cols, 
                  Study_col, samp_sizes, AgD_Trt_cols, TE_col, 
                  SE_col, IPD_Trt_col)




#' Goodness of interpolation interactive plot
#'
#' @param NMI_object The output of \code{\link{NMI_interpolation}}
#' 
#' @return An interactive plot displaying the predicted TE estimates from the 
#' NMI algorithm vs. the observed TE estimates.
#'
#' @export
NMI_diagnostic_plotly = function(NMI_object){
  df = NMI_object$Diagnostics
  df$Text = apply(df, 1,
                  function(u){
                    v = u[2:(length(u)-2)]
                    paste0('Study ', u[1], '\n',
                           paste0(names(v), ' = ', 
                                  sprintf('%.1f', 100*v), '%', '\n', 
                                  collapse="")
                    )
                  })
  
  p = ggplot(df, aes(TE_orig, TE_pred, text = Text,
                     fill = as.factor(Study))) + 
    geom_segment(x = min(min(df$TE_orig), min(df$TE_pred)), 
                 y = min(min(df$TE_orig), min(df$TE_pred)),
                 xend = max(max(df$TE_orig), max(df$TE_pred)),
                 yend = max(max(df$TE_orig), max(df$TE_pred)),
                 linetype = "dashed", color = "red") + 
    geom_point(size = 3, shape = 21, col = 'royal blue') + 
    theme(panel.background = element_blank(),
          legend.position = 'none',
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_line(colour = "grey92"),
          axis.title.x = element_text(size=13, face = "bold"),
          axis.title.y = element_text(size=13, face = "bold"),
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10),
          plot.title = element_text(size=14, face = "bold", hjust = .5)) + 
    xlab('Observed treatment effect estimate') + 
    ylab('Predicted treatment effect estimate') + 
    ggtitle('NMI goodness of fit') + 
    scale_fill_brewer(palette = "Set1")
  
  ggplotly(p, tooltip = "text")
}




#*******************************************************************************
#** NMI run
#*******************************************************************************
#imputing AgD
NMI_object = NMI_interpolation(IPD, AgD, x_vect, AgD_EM_cols, IPD_EM_cols, 
                               Study_col, samp_sizes, AgD_Trt_cols, TE_col, 
                               SE_col, IPD_Trt_col)
#Have a look at the imputed dataset
NMI_object$Imputed 
#The data submitted to NMA for ITC
NMI_object$Final

#Interpolation goodness of fit
NMI_diagnostic_plotly(NMI_object)


#*******************************************************************************
     
     #Estimate Vaishampayan TE from ARCHES
     #value of effect modifiers in Vaishampayan
     		
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
     
     x1 = studies_to_impute[studies_to_impute$Study == 'ARASENS', 'ECOG_0']
     x2 = studies_to_impute[studies_to_impute$Study == 'ARASENS', 'GLEASON_L7']
     x3 = studies_to_impute[studies_to_impute$Study == 'ARASENS', 'AGE_L65']
     x4 = studies_to_impute[studies_to_impute$Study == 'ARASENS', 'VOLUME_HIGH']
     
     
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
       
       studies_not_changed <- c("ARASENS", "PEACE-1", "STAMPEDE C vs G", "ENZAMET (docetaxel planned)", "ENZAMET (non-docetaxel)" )
       mHSPC_OS_data_PO_NMR <- read.csv("NMI/Data/mHSPC OS north america.csv")
       treatment_numbers <- read.csv("NMI/Data/treatment numbers.csv")
       new_row <- nrow(NMI_object$Final)
       
       for (i in 1:length(studies_not_changed)) {
         
         NMI_object$Final[new_row+i, "study_name"] = studies_not_changed[i]
         NMI_object$Final[new_row+i, "Trt1"] = mHSPC_OS_data_PO_NMR[mHSPC_OS_data_PO_NMR$X.ID == studies_not_changed[i], "t1"]
         NMI_object$Final[new_row+i, "Trt2"] = mHSPC_OS_data_PO_NMR[mHSPC_OS_data_PO_NMR$X.ID == studies_not_changed[i], "t2"]
         NMI_object$Final[new_row+i, "TE"] = mHSPC_OS_data_PO_NMR[mHSPC_OS_data_PO_NMR$X.ID == studies_not_changed[i], "y"]
         NMI_object$Final[new_row+i, "se"] = mHSPC_OS_data_PO_NMR[mHSPC_OS_data_PO_NMR$X.ID == studies_not_changed[i], "se"]
          }
       
       NMI_object$Final <- merge(NMI_object$Final, treatment_numbers, by="study_name")
       
     write.csv(NMI_object$Final, "NMI/Imputed Data/OS_data_imputed_arasens.csv")
     
     
