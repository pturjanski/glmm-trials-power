#------------------------------------------------------------
# Script Name:    build_trials_dataframe_GLMM_for_power_analysis.R
# Description:    This script performs the data processing and statistical
#                 analysis for the manuscript titled:
#                 "Why accounting for data clustering matters in parasitology: simulation-based insights
#                  from hierarchical infection models"
#
#                 In particular, it focuses on the generation of the dataframe 
#                 containing the trials and the corresponding GLMM models.
# Authors:        
#                 Helman, Elisa* (CONICET / GBA-FCEN-UBA / LAINPA-FCV-UNLP, ehelman@fcv.unlp.edu.ar)
#                 Perez, Adriana (GBA-FCEN-UBA, aaperez@ege.fcen.uba.ar)
#                 Turjanski, Pablo (GBA-FCEN-UBA // ICC-UBA-CONICET, pturjanski@dc.uba.ar)
#                 Unzaga, Juan Manuel (LAINPA-FCV-UNLP, junzaga2003@yahoo.es)
#                 Fernández, Maria Soledad  (GBA-FCEN-UBA // IC-UBA-CONICET, sfernandez@ic.fcen.uba.ar)
#
#                 *Corresponding Author 
#
# Date Created:   2025-03-31
# Last Modified:  2026-05-19
# Version:        1.1
#------------------------------------------------------------
# Notes:
# - This script is intended to be run in R version 4.6.0
# - Ensure all required packages are installed before execution.
#
# =========================================================
# Interactive execution notes
# =========================================================
#
# This script uses readline() to allow users to provide
# simulation parameters interactively.
#
# To run the script interactively in RStudio:
#   - Use the "Source" button located in the top-right corner
#     of the script editor, or
#   - Use the keyboard shortcut Ctrl + Shift + S.
#
# Another way to run the script interactively is:
#   1. Start an R session from the terminal.
#   2. Execute:
#
#      source("build_trials_dataframe_GLMM_for_power_analysis.R")
#
# In non-interactive execution modes, readline() may not
# pause execution to wait for user input, and default values
# may be used automatically.
#
# Default values:
#   input_bp             = 0.12
#   input_OR1            = 2.117
#   input_OR2            = 1.492
#   input_sigmaFarms     = 1
#   input_numberOfFarms  = c(10, 20, 50, 80, 120, 160)
#   input_numberOfPigs   = c(10, 20, 50, 80)
#
# =========================================================


#------------------------------------------------------------
# Load required packages
#------------------------------------------------------------
library(methods) 
library(boot)
library(glmmTMB)


#------------------------------------------------------------
# Define auxiliary functions
#------------------------------------------------------------


#-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
# functions for interactive console input
#-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
# Reads and validates a numeric value entered from the console
read_numeric <- function(prompt, default, validator, error_msg) {
  
  input <- readline(
    paste0(prompt, " [default: ", default, "]: ")
  )
  
  value <- if (nchar(trimws(input)) == 0) {
    default
  } else {
    as.numeric(input)
  }
  
  if (is.na(value) || !validator(value)) {
    stop(error_msg)
  }
  
  return(value)
}

# Reads and validates a vector of numeric values entered from the console
read_numeric_vector <- function(prompt, default) {
  
  default_text <- paste(default, collapse = ", ")
  
  input <- readline(
    paste0(prompt, " [default: ", default_text, "]: ")
  )
  
  values <- if (nchar(trimws(input)) == 0) {
    
    default
    
  } else {
    
    vals <- as.numeric(
      trimws(
        strsplit(input, ",")[[1]]
      )
    )
    
    # Numeric validation
    if (any(is.na(vals))) {
      stop(paste(prompt, ": all values must be numeric."))
    }
    
    # Values must be greater than or equal to 2
    if (any(vals < 2)) {
      stop(paste(prompt, ": all values must be greater than or equal to 2."))
    }
    
    # No duplicated values
    if (length(unique(vals)) != length(vals)) {
      stop(paste(prompt, ": duplicated values are not allowed."))
    }
    
    # Strictly increasing order
    if (is.unsorted(vals, strictly = TRUE)) {
      stop(paste(prompt, ": values must be sorted in increasing order."))
    }
    
    vals
  }
  
  return(values)
}


#-.-.-.-.-.-.-.-.-.-.-.-.-
# Original data generation
#-.-.-.-.-.-.-.-.-.-.-.-.-

# Generation of a data frame with predictor variables
predictorVars <- function(nFarms, nMinPigs, nMaxPigs, sigmaFarms, initialSeed){
  # Set a seed since random functions will be used, ensuring reproducibility
  set.seed(initialSeed)
  
  # Resulting data frame. Initialized as empty
  dataFrameResult <- c()

  for (iFarm in 1:nFarms) {
    # Farm's type (Intensive System(IS)/Semi-Extensive System (SES))
    farmType <- ifelse(iFarm <= (nFarms/2), "IS","SES") 
    # Random effect for each farm
    randomEffect <- rnorm(1,0,sigmaFarms)
    # Number of pigs in the farms (between nMinPigs and nMaxPigs)
    nPigs <- as.integer(runif(1, min = nMinPigs, max = nMaxPigs))
    for (iPig in 1:nPigs) {
      # Pig's age (Piglet (P)/ Weaned Pig (WP))
      pigAge <-  ifelse(iPig <= (nPigs/2), "P","WP")       
      # We generate the row with all the data
      dataRow <- c(iFarm, (iFarm*(nPigs+1))+iPig, farmType, pigAge, randomEffect)
      # We add the row to the resulting dataframe
      dataFrameResult <- cbind.data.frame(rbind(dataFrameResult, dataRow))
    }
  }
  
  # Assign names to the columns
  colnames(dataFrameResult)<-c("idFarm", "iPig", "farmType", "pigAge", "randomEffect")

  # Assign the class to each column in the dataframe
  dataFrameResult$idFarm <- as.factor(dataFrameResult$idFarm)
  dataFrameResult$iPig <- as.factor(dataFrameResult$iPig)
  dataFrameResult$farmType <- as.factor(dataFrameResult$farmType)
  dataFrameResult$pigAge <- as.factor(dataFrameResult$pigAge)
  dataFrameResult$randomEffect <- as.numeric(dataFrameResult$randomEffect)
  
  
  return(dataFrameResult)
}

# Add the response variable to the predictor variables dataframe
responseVar <- function(dfPredictorVars, initialSeed, beta_0, beta_1, beta_2){

  # Set a seed since random functions will be used, ensuring reproducibility
  set.seed(initialSeed)
  
  # Calculate the linear predictor using given coefficients and predictor variables
  dfPredictorVars$linear <- beta_0 + 
                             beta_1 * (dfPredictorVars$farmType == "SES") +  
                             beta_2 * (dfPredictorVars$pigAge == "WP")    + 
                             dfPredictorVars$randomEffect
 
  # Apply the inverse logit function to get probabilities (pp) from the linear predictor
  dfPredictorVars$pp <- inv.logit(dfPredictorVars$linear)
  
  # Get the total number of rows (observations) in the data frame
  nTotal <- nrow(dfPredictorVars)
  
  # Simulate binary response variable 'y' from a binomial distribution with probability 'pp'
  dfPredictorVars$y <- rbinom(nTotal,1, dfPredictorVars$pp)
  
  # Return the modified data frame with new variables added
  return(dfPredictorVars)
}


#-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
# Functions to extract p-values from predictors 
# using LTR tests with drop1 function ()
#-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.


# Get p-value from b1 ('farmType')
drop1_farmType <- function(model){     
   return(drop1(model, test = "Chisq")$`Pr(>Chi)`[2]) 
}

# Get p-value from b2 ('pigAge')
drop1_pigAge <- function(model){
   return(drop1(model, test = "Chisq")$`Pr(>Chi)`[3])
}


#-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
# Retrieving model data
#-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-


# Generate the row with model data
getModelDataAnalysis <- function(xtrial){
  return(
    c( 
      xtrial$nFarms,
      xtrial$nMinPigs,
      xtrial$nMaxPigs,
      xtrial$sigmaFarms,
      xtrial$initialPredictorSeed,
      xtrial$beta_0,
      xtrial$beta_1,
      xtrial$beta_2,
      xtrial$initialResponseSeed,
      drop1_farmType(xtrial$GLMM),
      drop1_pigAge(xtrial$GLMM)
      )                          
  )
}

# Transform trialAnalysis into a data frame (and apply the appropriate format)
setTrialAnalysisFormat <- function(xtrialAnalysis){
  # Convert to data frame
  xtrialAnalysis <- as.data.frame(xtrialAnalysis, row.names = FALSE)
  
  # Assign column names
  colnames(xtrialAnalysis) <- c(
    "nFarms",
    "nMinPigs",
    "nMaxPigs",
    "sigmaFarms",
    "initialPredictorSeed",
    "beta_0",
    "beta_1",
    "beta_2",
    "initialResponseSeed",
    "drop1_farmType",
    "drop1_pigAge"
  )  
  
  # Assign appropriate data types to each variable
  xtrialAnalysis$nFarms               <- as.numeric(xtrialAnalysis$nFarms)
  xtrialAnalysis$nMinPigs             <- as.numeric(xtrialAnalysis$nMinPigs)
  xtrialAnalysis$nMaxPigs             <- as.numeric(xtrialAnalysis$nMaxPigs)
  xtrialAnalysis$sigmaFarms           <- as.numeric(xtrialAnalysis$sigmaFarms)
  xtrialAnalysis$initialPredictorSeed <- as.numeric(xtrialAnalysis$initialPredictorSeed)
  xtrialAnalysis$beta_0               <- as.numeric(xtrialAnalysis$beta_0)
  xtrialAnalysis$beta_1               <- as.numeric(xtrialAnalysis$beta_1)
  xtrialAnalysis$beta_2               <- as.numeric(xtrialAnalysis$beta_2)
  xtrialAnalysis$initialResponseSeed  <- as.numeric(xtrialAnalysis$initialResponseSeed)
  xtrialAnalysis$drop1_farmType       <- as.numeric(xtrialAnalysis$drop1_farmType)
  xtrialAnalysis$drop1_pigAge         <- as.numeric(xtrialAnalysis$drop1_pigAge)

  return(xtrialAnalysis)
}

#-.-.-.-.-.-.-.-.-.-.-.-
# Data class definitions
#-.-.-.-.-.-.-.-.-.-.-.-


# 'trial' class. Allows generating a dataframe with the trial and GLMM
trial <- setRefClass("trial", fields = list(nFarms                = "numeric",
                                            nMinPigs              = "numeric",
                                            nMaxPigs              = "numeric",
                                            sigmaFarms            = "numeric",
                                            initialPredictorSeed  = "numeric",
                                            beta_0                = "numeric",
                                            beta_1                = "numeric",
                                            beta_2                = "numeric",
                                            initialResponseSeed   = "numeric",
                                            df                    = "data.frame",
                                            GLMM                  = "ANY"
                                            )
                     , methods = list( 
                        initialize = function(xnFarms, 
                                              xnMinPigs,
                                              xnMaxPigs,
                                              xsigmaFarms,
                                              xinitialPredictorSeed,
                                              xbeta_0,
                                              xbeta_1,
                                              xbeta_2,
                                              xinitialResponseSeed
                                              ) {

                          # Assign values to the internal attributes of the object
                          nFarms                <<- xnFarms
                          nMinPigs              <<- xnMinPigs
                          nMaxPigs              <<- xnMaxPigs
                          sigmaFarms            <<- xsigmaFarms
                          initialPredictorSeed  <<- xinitialPredictorSeed
                          beta_0                <<- xbeta_0
                          beta_1                <<- xbeta_1
                          beta_2                <<- xbeta_2
                          initialResponseSeed   <<- xinitialResponseSeed
                          
                          # Generate data frame with predictor variables
                          df <<- predictorVars(xnFarms, xnMinPigs, xnMaxPigs, xsigmaFarms, xinitialPredictorSeed)
                          # Add the response variable to the data frame
                          df  <<- responseVar(df, xinitialResponseSeed, xbeta_0, xbeta_1, xbeta_2)
                          # Generate GLMM 
                          GLMM <<- glmmTMB(y ~ farmType + pigAge + (1|idFarm), data = df, family = binomial)
                          
                        }
                     )
) 



#-----------------------------------------------------------
# Generation of multiple simulated trials and their analysis
#-----------------------------------------------------------


# Initialize containers to store trials and their corresponding analysis
trials = c()
trialAnalysis = c()


# Interactive inputs
# ..................
# Define the input parameters used in the simulation:
# bp: Baseline prevalence, representing the expected proportion of positive cases
#     in the absence of risk factors.
# OR1: Odds ratio associated with the farm-level predictor.
# OR2: Odds ratio associated with the pig-level predictor.
# sigmaFarms: Between-farm variability, expressed as the standard deviation
#              of the random farm effect.
# numberOfFarms: Numbers of farms to be evaluated in the power curves.
# numberOfPigs: Numbers of pigs per farm to be evaluated in the power curves.


# Baseline prevalence
input_bp <- read_numeric(
  prompt     = "Baseline prevalence (bp)",
  default    = 0.12,
  validator  = function(x) x > 0 && x < 1,
  error_msg  = "bp must be between 0 and 1."
)

# Group-level predictor odds ratio
input_OR1 <- read_numeric(
  prompt     = "Group-level predictor odds ratio (OR1)",
  default    = 2.117,
  validator  = function(x) x > 0,
  error_msg  = "OR1 must be greater than 0."
)

# Individual-level predictor odds ratio
input_OR2 <- read_numeric(
  prompt     = "Individual-level predictor odds ratio (OR2)",
  default    = 1.492,
  validator  = function(x) x > 0,
  error_msg  = "OR2 must be greater than 0."
)

# Between-group variability
input_sigmaFarms <- read_numeric(
  prompt     = "Between-group variability (sigmaFarms)",
  default    = 1,
  validator  = function(x) x >= 0,
  error_msg  = "sigmaFarms must be greater than or equal to 0."
)

# Number of farms
input_numberOfFarms <- read_numeric_vector(
  prompt  = "numberOfFarms (comma-separated values)",
  default = c(10, 20, 50, 80, 120, 160)
)

# Number of pigs
input_numberOfPigs <- read_numeric_vector(
  prompt  = "numberOfPigs (comma-separated values)",
  default = c(10, 20, 50, 80)
)

#-.-.-.-.-.-.-.-.-.-.-.-.-
# Display selected values
#-.-.-.-.-.-.-.-.-.-.-.-.-

cat("\nSelected values:\n")

cat("bp =", input_bp, "\n")
cat("OR1 =", input_OR1, "\n")
cat("OR2 =", input_OR2, "\n")
cat("sigmaFarms =", input_sigmaFarms, "\n")

cat(
  "numberOfFarms =",
  paste(input_numberOfFarms, collapse = ", "),
  "\n"
)

cat(
  "numberOfPigs =",
  paste(input_numberOfPigs, collapse = ", "),
  "\n"
)


# Assign interactive input values to the model variables used in the simulation:
bp <-  input_bp
OR1 <- input_OR1
OR2 <- input_OR2


# Loop over different number of farms
for(iNFarms in input_numberOfFarms){	  
  print(paste("Running number of farms (iNFarms): ", iNFarms))		
  
  # Loop over different number of pigs per farm
  for(iNPigs in input_numberOfPigs){	      	
    print(paste("Running number of pigs per farm (iNPigs): ", iNPigs))		
    
    # Repeat each scenario a fixed number of times (number of simulations)
    for(iSim in 1:500){	              

      # Set simulation parameters
      nFarms               <- iNFarms          # Number of farms
      nMinPigs             <- iNPigs           # Minimum number of pigs per farm
      nMaxPigs             <- iNPigs           # Maximum number of pigs per farm
      sigmaFarms           <- input_sigmaFarms # Inter-farm variability (Standard deviation)
      initialPredictorSeed <- iSim             # Seed for random effect in the predictor generation function
      beta_0               <- logit(bp)        # Intercept
      beta_1               <- log(OR1)         # Coefficient for predictor 1
      beta_2               <- log(OR2)         # Coefficient for predictor 2
      initialResponseSeed  <- iSim             # Seed for response variable
      
      
      # Generate a single trial based on the parameters
      t1 <- trial(nFarms, nMinPigs, nMaxPigs, sigmaFarms, initialPredictorSeed, beta_0, beta_1, beta_2, initialResponseSeed)
      
      # Store the generated trial
      trials = cbind(rbind(trials, c(t1) ))
      
      # Store the analysis result of the trial 
      trialAnalysis = cbind(rbind(trialAnalysis, getModelDataAnalysis(t1) ))
      
     }
  }

}
  

# Format 'trialAnalysis' as a structured data frame for further analysis
trialAnalysis <- setTrialAnalysisFormat(trialAnalysis)

# ------------------------------------------------------------
# Save results
# ------------------------------------------------------------

write.csv(trialAnalysis, "simulated_trial_analysis.csv", row.names = FALSE)

#------------------------------------------------------------
# End of script 
#------------------------------------------------------------


