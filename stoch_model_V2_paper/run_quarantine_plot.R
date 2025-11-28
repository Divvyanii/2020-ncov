# Standalone script to generate quarantine effect plot
# This script loads everything needed and generates the graph

# Set working directory
setwd("/Users/divyanisingh/2020-ncov/stoch_model_V2_paper")

# Load required libraries
library(foreach)
library(doMC)
library(lubridate)
library(magrittr)
library(coda)
library(tidyverse)
library(rootSolve)
library(mgcv)

# Register cores
registerDoMC(4)

# Load datasets
dropbox_path <- ""
travel_data_mobs <- read_csv(paste0(dropbox_path,"data/connectivity_data_mobs.csv"))
international_conf_data_in <- read_csv(paste0(dropbox_path,"data/international_case_data.csv"))
international_onset_data_in <- read_csv(paste0(dropbox_path,"data/time_series_WHO_report.csv"))
china_onset_data_in <- read_csv(paste0(dropbox_path,"data/time_series_data_bioRvix_Liu_et_al.csv"))
wuhan_onset_data_in <- read_csv(paste0(dropbox_path,"data/time_series_data_lancet_huang_et_al.csv"))
wuhan_onset_2020_01_30 <- read_csv(paste0(dropbox_path,"data/time_series_data_qui_li_nejm_wuhan.csv"))
wuhan_conf_data_in <- read_csv(paste0(dropbox_path,"data/time_series_HKU_Wuhan.csv"))
data_hubei_Feb <- read_csv(paste0(dropbox_path,"data/hubei_confirmed_cases.csv"))

case_data_in <- international_conf_data_in
travel_data <- travel_data_mobs
t_step <- 0.25

# Load model and plotting functions
source("R/model_functions.R")
source("R/plotting_functions.R")

# Load model parameters
thetaR_IC <- read_csv("inputs/theta_initial_conditions.csv")
theta <- c( r0=as.numeric(thetaR_IC[thetaR_IC$param=="r0","value"]),
            beta=NA,
            betavol=as.numeric(thetaR_IC[thetaR_IC$param=="betavol","value"]),
            gentime=as.numeric(thetaR_IC[thetaR_IC$param=="gentime","value"]),
            incubation = 1/as.numeric(thetaR_IC[thetaR_IC$param=="incubation","value"]),
            report = 1/as.numeric(thetaR_IC[thetaR_IC$param=="report","value"]),
            report_local = 1/as.numeric(thetaR_IC[thetaR_IC$param=="report_local","value"]),
            recover = 1/as.numeric(thetaR_IC[thetaR_IC$param=="recover","value"]),
            init_cases=as.numeric(thetaR_IC[thetaR_IC$param=="init_cases","value"]),
            passengers=as.numeric(thetaR_IC[thetaR_IC$param=="outbound_travel","value"]),
            pop_travel=as.numeric(thetaR_IC[thetaR_IC$param=="population_travel","value"]),
            local_rep_prop=as.numeric(thetaR_IC[thetaR_IC$param=="local_rep_prop","value"]),
            onset_prop=as.numeric(thetaR_IC[thetaR_IC$param=="onset_prop","value"]),
            onset_prop_int=as.numeric(thetaR_IC[thetaR_IC$param=="onset_prop_int","value"]),
            confirmed_prop=as.numeric(thetaR_IC[thetaR_IC$param=="confirmed_prop","value"]),
            travel_frac=NA,
            r0_decline =as.numeric(thetaR_IC[thetaR_IC$param=="r0_decline","value"]),
            rep_local_var =as.numeric(thetaR_IC[thetaR_IC$param=="rep_local_var","value"]),
            pre_symp =as.numeric(thetaR_IC[thetaR_IC$param=="pre_symp","value"]),
            quarantine_effectiveness =as.numeric(thetaR_IC[thetaR_IC$param=="quarantine_effectiveness","value"]),
            quarantine_start_date =as.character(thetaR_IC[thetaR_IC$param=="quarantine_start_date","value"])
)

theta[["travel_frac"]] <- theta[["passengers"]]/theta[["pop_travel"]]
theta[["beta"]] <- theta[["r0"]]*(theta[["recover"]])
theta_initNames <- c("sus","tr_exp1","tr_exp2","exp1","exp2","inf1","inf2","tr_waiting","cases","reports","waiting_local","cases_local","reports_local")

# Load timeseries data
source("R/load_timeseries_data.R")

# Update quarantine start time
if(!is.null(theta[["quarantine_start_date"]])){
  quarantine_start_date <- as.Date(theta[["quarantine_start_date"]])
  quarantine_start_time <- as.numeric(quarantine_start_date - start_date + 1)
}

cat("Setup complete. Generating quarantine effect plot...\n")

# Load and run quarantine plotting function
source("R/plot_quarantine_effect.R")

# Ensure plots directory exists
if(!dir.exists("plots")){
  dir.create("plots")
}

# Generate quick comparison plot
# This will save as PDF automatically
plot_quarantine_quick(
  quarantine_levels = c(0, 0.3, 0.5, 0.7),
  nn = 1e3,
  dt = 0.25,
  filename = "quarantine_effect_quick"  # Saves as plots/quarantine_effect_quick.pdf
)

cat("\n✓ Plot generated successfully!\n")
cat("  Saved to: plots/quarantine_effect_quick.pdf\n")
cat("  File size:", file.info("plots/quarantine_effect_quick.pdf")$size, "bytes\n")

