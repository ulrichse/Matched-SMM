
## This code calculates cumulative lag relative risk for a heatwave period (lag0 to lag3, lag0 to lag7)
## and plots the cumulative lag for different outcomes and subgroups using the forestploter package
## For subgroup analyses, you have to restrict the delivery data to the subgroup before the matching procedure. 

#### DATA PREP ####

library(dplyr)
library(arrow)
library(lubridate)
library(data.table)
library(haven)
library(sjPlot)
library(MASS)
library(grid)
library(forestploter)
library(dlnm)
library(lme4)

setwd("~/RStudio/Matched SMM")
getwd()

#### SET UP FUNCTIONS ####

# Function to exclude ZCTAs with no heatwave days during the study period from matching
zctas_without_hw <- function(data){
  
  no_hw <- data %>%
    group_by(zip)%>%
    summarize(hw=sum(heatwave))%>% # Create list of ZCTAs without heatwaves
    filter(hw==0)%>%
    dplyr::select(zip)
  
  data <- data %>% # Remove ZCTAs without heatwaves from the data
    dplyr::filter(!zip %in% no_hw)
  
}

# Create matched dataframe with lags -2 to 7
matched_df_lag7 <- function(data, n_controls = 3, lag_range = -2:7, control_doy_range = -3:3) {
  
  data <- zctas_without_hw(data)
  dat <- data
  zip_list <- unique(dat$zip)
  setorder(dat, zip, date)
  
  for (i in 1:length(zip_list)) {
    df <- subset(dat, zip == zip_list[i])
    
    # Exclude the 3 days within any other heatwave
    df$time <- 1:nrow(df)
    cand_control <- unique(c(which(df$heatwave == 1), which(df$heatwave == 1) + 1, which(df$heatwave == 1) - 1))
    df$cand_control <- TRUE
    df$cand_control[cand_control[cand_control <= nrow(df)]] <- FALSE
    
    case_dates <- subset(df, heatwave == 1)
    control_dates <- subset(df, heatwave == 0)
    
    for (j in 1:nrow(case_dates)) {
      # Choose lags (lagged 0 to 3)
      lag_dates <- case_dates[j, ]$date + lag_range
      lag_case <- subset(df, date %in% lag_dates)
      
      # Choose 10 comparable unexposed days for each heatwave-exposed day
      control_range <- case_dates[j, ]$doy + control_doy_range
      control_subset <- subset(control_dates,
                               control_dates$year != case_dates[j, ]$year &
                                 doy %in% control_range &
                                 cand_control)
      controls <- dplyr::sample_n(control_subset, n_controls)
      
      # Choose lagged days for selected unexposed days
      la_con <- c(-2:-1, 1:7)
      for (p in 1:length(la_con)) {
        lag_control_dates <- controls$date + la_con[p]
        lag_control_each <- subset(df, date %in% lag_control_dates)
        
        if (p == 1) {
          lag_control <- lag_control_each
        } else {
          lag_control <- rbind(lag_control, lag_control_each)
        }
      }
      j_stratum <- rbind(lag_case, controls, lag_control)
      stratum <- paste("stratum", j, sep = ".")
      j_stratum$stratum <- stratum
      status <- c(rep("case", nrow(lag_case)), rep("control", nrow(controls)), rep("control", nrow(lag_control)))
      j_stratum$status <- status
      lag <- c(rep(-2:7, length.out = nrow(lag_case)), rep(0, length.out = nrow(controls)), rep(c(-2:-1, 1:7), length.out = nrow(lag_control)))
      j_stratum$lag <- lag
      
      if (j == 1) {
        new_df <- j_stratum
      } else {
        new_df <- rbind(new_df, j_stratum)
      }
    }
    if (i == 1) {
      matched_df <- new_df
    } else {
      matched_df <- rbind(matched_df, new_df)
    }
  }
  
  return(matched_df)
  
  gc()
  
}

# Create crossbasis for lags -2 to 7
matched_cb_lagn2_lag7 <- function(data, n_controls = 3, lag_range = -2:7, control_doy_range = -3:3){
  
  data <- zctas_without_hw(data)
  dat <- data
  zip_list <- unique(dat$zip)
  setorder(dat, zip, date)
  
  # Use "dlnm" package to generate the distributed lag function for "heatwave"
  for (i in 1:length(zip_list)) {
    orig_dat <- subset(dat, zip == zip_list[i])
    match_dat <- subset(matched_df, zip == zip_list[i])
    
    orig_cb <- dlnm::crossbasis(orig_dat$heatwave, lag = c(-2, 7),
                                argvar = list(fun = "lin"),
                                arglag = list(fun = "integer"))
    obs_n <- nrow(orig_dat)
    orig_cb_matr <- as.data.frame(subset(orig_cb, nrow = obs_n))
    orig_cb_matr$date <- orig_dat$date
    matched_date <- match_dat %>% dplyr::select(date)
    matched_cb_matr <- orig_cb_matr %>%
      dplyr::right_join(matched_date, by = "date") %>%
      dplyr::select(-date) %>% as.matrix()
    
    if (i == 1) {
      matched_cb_matrix <- matched_cb_matr
    } else {
      matched_cb_matrix <- rbind(matched_cb_matrix, matched_cb_matr)
    }
    # Add attributes
    matched_dim <- dim(matched_cb_matrix)
    attr <- attributes(orig_cb)
    attr$dim <- matched_dim
    matched_cb <- matched_cb_matrix
    attributes(matched_cb) <- attr
  }
  
  return(matched_cb)
  
  gc()
  
}

# Function to generate individual lags from -2 to 7
get_daily_lags <- function(matched_df, matched_cb){
  
  datasets <- list()
  combined_data_list <- list()
  
  for (outcome in outcomes) {
    formula <- as.formula(paste0(outcome, " ~ matched_cb + factor(dow) + year * zip"))
    
    fit_am5 <- gnm::gnm(formula,
                        eliminate = factor(zip), family = quasipoisson(link = "log"),
                        data = matched_df)
    pred_am5 <- dlnm::crosspred(matched_cb, fit_am5, at = 1)
    
    # Print matRRfit, matRRlow, and matRRhigh
    print(pred_am5$matRRfit)
    print(pred_am5$matRRlow)
    print(pred_am5$matRRhigh)
    
    datasets[[paste0(outcome, "_pred_am5")]] <- pred_am5
    
    combined_data_name <- paste0(outcome, "_combined_data")
    
    # Create an empty data frame with dynamic column names
    max_columns <- max(sapply(datasets, function(x) ncol(x$matRRfit), simplify = TRUE),
                       sapply(datasets, function(x) ncol(x$matRRlow), simplify = TRUE),
                       sapply(datasets, function(x) ncol(x$matRRhigh), simplify = TRUE))
    
    combined_data <- data.frame(matrix(NA, nrow = 0, ncol = max_columns + 2))
    colnames(combined_data) <- c(names(datasets[[1]]$matRRfit), "Type", "Pred_Model")
    
    # Populate the combined data frame
    for (name in names(datasets)) {
      matRRfit <- datasets[[name]]$matRRfit
      matRRlow <- datasets[[name]]$matRRlow
      matRRhigh <- datasets[[name]]$matRRhigh
      
      combined_data <- rbind(combined_data, cbind(matRRfit, Type = "matRRfit", Pred_Model = name))
      combined_data <- rbind(combined_data, cbind(matRRlow, Type = "matRRlow", Pred_Model = name))
      combined_data <- rbind(combined_data, cbind(matRRhigh, Type = "matRRhigh", Pred_Model = name))
    }
    
    # Assign the combined data frame to the list with dynamic name
    combined_data_list[[combined_data_name]] <- combined_data
  }
  
  return(combined_data)
  
  gc()
  
}

# Function to create crossbasis for lags 0 to 3
matched_cb_lag3 <- function(data, n_controls = 3, lag_range = 0:3, control_doy_range = -3:3){
  
  data <- zctas_without_hw(data)
  dat <- data
  zip_list <- unique(dat$zip)
  setorder(dat, zip, date)
  
  # Use "dlnm" package to generate the distributed lag function for "heatwave"
  for (i in 1:length(zip_list)) {
    orig_dat <- subset(dat, zip == zip_list[i])
    match_dat <- subset(matched_df, zip == zip_list[i])
    
    orig_cb <- dlnm::crossbasis(orig_dat$heatwave, lag = c(0, 3),
                                argvar = list(fun = "lin"),
                                arglag = list(fun = "integer"))
    obs_n <- nrow(orig_dat)
    orig_cb_matr <- as.data.frame(subset(orig_cb, nrow = obs_n))
    orig_cb_matr$date <- orig_dat$date
    matched_date <- match_dat %>% dplyr::select(date)
    matched_cb_matr <- orig_cb_matr %>%
      dplyr::right_join(matched_date, by = "date") %>%
      dplyr::select(-date) %>% as.matrix()
    
    if (i == 1) {
      matched_cb_matrix <- matched_cb_matr
    } else {
      matched_cb_matrix <- rbind(matched_cb_matrix, matched_cb_matr)
    }
    # Add attributes
    matched_dim <- dim(matched_cb_matrix)
    attr <- attributes(orig_cb)
    attr$dim <- matched_dim
    matched_cb <- matched_cb_matrix
    attributes(matched_cb) <- attr
  }
  
  return(matched_cb)
  
  gc()
  
}

# Function to create crossbasis for lags 0 to 7
matched_cb_lag7 <- function(data, n_controls = 3, lag_range = 0:7, control_doy_range = -3:3){
  
  data <- zctas_without_hw(data)
  dat <- data
  zip_list <- unique(dat$zip)
  setorder(dat, zip, date)
  
  # Use "dlnm" package to generate the distributed lag function for "heatwave"
  for (i in 1:length(zip_list)) {
    orig_dat <- subset(dat, zip == zip_list[i])
    match_dat <- subset(matched_df, zip == zip_list[i])
    
    orig_cb <- dlnm::crossbasis(orig_dat$heatwave, lag = c(0, 7),
                                argvar = list(fun = "lin"),
                                arglag = list(fun = "integer"))
    obs_n <- nrow(orig_dat)
    orig_cb_matr <- as.data.frame(subset(orig_cb, nrow = obs_n))
    orig_cb_matr$date <- orig_dat$date
    matched_date <- match_dat %>% dplyr::select(date)
    matched_cb_matr <- orig_cb_matr %>%
      dplyr::right_join(matched_date, by = "date") %>%
      dplyr::select(-date) %>% as.matrix()
    
    if (i == 1) {
      matched_cb_matrix <- matched_cb_matr
    } else {
      matched_cb_matrix <- rbind(matched_cb_matrix, matched_cb_matr)
    }
    # Add attributes
    matched_dim <- dim(matched_cb_matrix)
    attr <- attributes(orig_cb)
    attr$dim <- matched_dim
    matched_cb <- matched_cb_matrix
    attributes(matched_cb) <- attr
  }
  
  return(matched_cb)
  
  gc()
  
}

# Function to generate cumulative RR from lag0 to lag7
cumulative_lag0_lag7 <- function(matched_df, matched_cb){
  
  outcome_data <- list()
  
  for (outcome in outcomes) {
    formula <- as.formula(paste0(outcome, " ~ matched_cb + factor(dow) + year * zip"))
    
    fit <- gnm::gnm(formula,
                    eliminate = factor(zip), family = quasipoisson(link = "log"),
                    data = matched_df)
    
    pred <- dlnm::crosspred(matched_cb, fit, at = 1)
    
    over_rr <- sum(pred$matRRfit) / 8
    
    library(msm)
    estvar <- pred$vcov
    estmean <- c(pred$coefficients)
    
    over_rr_se <- msm::deltamethod(~ (exp(x1) + exp(x2) + exp(x3) + exp(x4) + exp(x5) + exp(x6) + exp(x7) + exp(x8)) / 8, 
                                   estmean, estvar)
    over_rr_low <- over_rr / exp(1.96 * over_rr_se)
    over_rr_high <- over_rr * exp(1.96 * over_rr_se)
    
    # Create a data frame for the current outcome
    outcome_df <- data.frame(outcome = outcome,
                             over_rr = over_rr,
                             over_rr_se = over_rr_se,
                             over_rr_low = over_rr_low,
                             over_rr_high = over_rr_high,
                             estvar = estvar,
                             estmean = estmean)
    
    outcome_data[[paste0(outcome, "_cumlag")]] <- outcome_df
  }
  
  data_frames <- unname(outcome_data)
  merged_df <- do.call(rbind, data_frames)
  return(merged_df)
  
  gc()
  
}

# Cumulative RR for lag0 to lag3
cumulative_lag0_lag3 <- function(matched_df, matched_cb){
  
  outcome_data <- list()
  
  for (outcome in outcomes) {
    formula <- as.formula(paste0(outcome, " ~ matched_cb + factor(dow) + year * zip"))
    
    fit <- gnm::gnm(formula,
                    eliminate = factor(zip), family = quasipoisson(link = "log"),
                    data = matched_df)
    
    pred <- dlnm::crosspred(matched_cb, fit, at = 1)
    
    over_rr <- sum(pred$matRRfit) / 4
    
    library(msm)
    estvar <- pred$vcov
    estmean <- c(pred$coefficients)
    
    over_rr_se <- msm::deltamethod(~ (exp(x1) + exp(x2) + exp(x3) + exp(x4)) / 4, 
                                   estmean, estvar)
    over_rr_low <- over_rr / exp(1.96 * over_rr_se)
    over_rr_high <- over_rr * exp(1.96 * over_rr_se)
    
    # Create a data frame for the current outcome
    outcome_df <- data.frame(outcome = outcome,
                             over_rr = over_rr,
                             over_rr_se = over_rr_se,
                             over_rr_low = over_rr_low,
                             over_rr_high = over_rr_high,
                             estvar = estvar,
                             estmean = estmean)
    
    outcome_data[[paste0(outcome, "_cumlag")]] <- outcome_df
  }
  
  data_frames <- unname(outcome_data)
  merged_df <- do.call(rbind, data_frames)
  return(merged_df)
  
  gc()
  
}

# Cumulative RR for lag-2 to lag7
cumulative_lagneg2_lag7 <- function(matched_df, matched_cb){
  
  outcome_data <- list()
  
  for (outcome in outcomes) {
    formula <- as.formula(paste0(outcome, " ~ matched_cb + factor(dow) + year * zip"))
    
    fit <- gnm::gnm(formula,
                    eliminate = factor(zip), family = quasipoisson(link = "log"),
                    data = matched_df)
    
    pred <- dlnm::crosspred(matched_cb, fit, at = 1)
    
    over_rr <- sum(pred$matRRfit) / 10
    
    library(msm)
    estvar <- pred$vcov
    estmean <- c(pred$coefficients)
    
    over_rr_se <- msm::deltamethod(~ (exp(x1) + exp(x2) + exp(x3) + exp(x4) + exp(x5) + exp(x6) + exp(x7) + exp(x8) + exp(x9) + exp(x10)) / 10, 
                                   estmean, estvar)
    over_rr_low <- over_rr / exp(1.96 * over_rr_se)
    over_rr_high <- over_rr * exp(1.96 * over_rr_se)
    
    # Create a data frame for the current outcome
    outcome_df <- data.frame(outcome = outcome,
                             over_rr = over_rr,
                             over_rr_se = over_rr_se,
                             over_rr_low = over_rr_low,
                             over_rr_high = over_rr_high,
                             estvar = estvar,
                             estmean = estmean)
    
    outcome_data[[paste0(outcome, "_cumlag")]] <- outcome_df
  }
  
  data_frames <- unname(outcome_data)
  merged_df <- do.call(rbind, data_frames)
  return(merged_df)
  
  gc()
  
}


##### FULL #####

delivzip <- read_parquet("Data/Outcomes_SMM_v2.parquet")%>%
  group_by(zip, date)%>%
  summarize(n=n(), 
            SMM=sum(SMM21))

data <- read_parquet("Data/Heatwave_Metrics_v2.parquet")%>%
  left_join(delivzip, by=c('zip', 'date'))
data[is.na(data)]<-0

# Create matched dataframe

matched_df <- matched_df_lag7(data)
table(matched_df$lag)
#matched_filter <- matched_df %>%
  #dplyr::select(-SMM)
#write_parquet(matched_filter, "matched_df_dates.parquet")

# Create matched crossbasis

matched_cb <- matched_cb_lagn2_lag7(data)

outcomes <- c("SMM")

combined_data <- get_daily_lags(matched_df, matched_cb)
write.csv(combined_data, file = "Results/Heatwaves/Daily Lags/combined_matrices_full.csv", row.names = FALSE)
gc()

merged_df <- cumulative_lagneg2_lag7(matched_df, matched_cb)
merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="All")
write.csv(merged_df, "Results/Heatwaves/Cumulative_Lagn2_Lag7_Full.csv")
gc()

matched_df <- matched_df %>%
  filter(lag <= 7 & lag >= 0)
table(matched_df$lag)

matched_cb <- matched_cb_lag7(data)
merged_df <- cumulative_lag0_lag7(matched_df, matched_cb)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="All ED Visits")
write.csv(merged_df, "Results/Heatwaves/Lag7/Cumulative_Lag7_Full.csv")
gc()

matched_df <- matched_df %>%
  filter(lag <= 3 & lag >= 0)
table(matched_df$lag)

matched_cb <- matched_cb_lag3(data)
merged_df <- cumulative_lag0_lag3(matched_df, matched_cb)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="All ED Visits")
write.csv(merged_df, "Results/Heatwaves/Lag3/Cumulative_Lag3_Full.csv")
gc()


##### BLACK #####

delivzip <- read_parquet("Data/Outcomes_SMM_v2.parquet")%>%
  filter(Race=="Black") %>%
  group_by(zip, date) %>%
  summarize(n=n(),
            SMM=sum(SMM21))

data <- read_parquet("Data/Heatwave_Metrics_v2.parquet")%>%
  left_join(delivzip, by=c('zip', 'date'))

matched_df <- read_parquet("Data/matched_df_dates.parquet")%>%
  left_join(delivzip, by=c('zip', 'date'))
matched_df[is.na(matched_df)]<-0

matched_cb <- matched_cb_lagn2_lag7(data)

outcomes <- c("SMM")

combined_data <- get_daily_lags(matched_df, matched_cb)
write.csv(combined_data, file = "Results/Heatwaves/Daily Lags/combined_matrices_Black.csv", row.names = FALSE)
gc()

merged_df <- cumulative_lagneg2_lag7(matched_df, matched_cb)
merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="Black")
write.csv(merged_df, "Results/Heatwaves/Lagn2/Cumulative_Lagn2_Lag7_Black.csv")
gc()

matched_df <- matched_df %>%
  filter(lag <= 7 & lag >= 0)
table(matched_df$lag)

matched_cb <- matched_cb_lag7(data)
merged_df <- cumulative_lag0_lag7(matched_df, matched_cb)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="Black")
write.csv(merged_df, "Results/Heatwaves/Lag7/Cumulative_Lag7_Black.csv")
gc()

matched_df <- matched_df %>%
  filter(lag <= 3 & lag >= 0)
table(matched_df$lag)

matched_cb <- matched_cb_lag3(data)
merged_df <- cumulative_lag0_lag3(matched_df, matched_cb)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="Black")
write.csv(merged_df, "Results/Heatwaves/Lag3/Cumulative_Lag3_Black.csv")
gc()

#### HISPANIC ####

delivzip <- read_parquet("Data/Outcomes_SMM_v2.parquet")%>%
  filter(Ethnicity=="Hispanic") %>%
  group_by(zip, date) %>%
  summarize(n=n(),
            SMM=sum(SMM21))

data <- read_parquet("Data/Heatwave_Metrics_v2.parquet")%>%
  left_join(delivzip, by=c('zip', 'date'))

matched_df <- read_parquet("Data/matched_df_dates.parquet")%>%
  left_join(delivzip, by=c('zip', 'date'))
matched_df[is.na(matched_df)]<-0

data[is.na(data)] <- 0

matched_cb <- matched_cb_lagn2_lag7(data)

outcomes <- c("SMM")

combined_data <- get_daily_lags(matched_df, matched_cb)
write.csv(combined_data, file = "Results/Heatwaves/Daily Lags/combined_matrices_Hispanic.csv", row.names = FALSE)
gc()

merged_df <- cumulative_lagneg2_lag7(matched_df, matched_cb)
merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="Hispanic")
write.csv(merged_df, "Results/Heatwaves/Lagn2/Cumulative_Lagn2_Lag7_Hispanic.csv")
gc()

matched_df <- matched_df %>%
  filter(lag <= 7 & lag >= 0)
table(matched_df$lag)

matched_cb <- matched_cb_lag7(data)
merged_df <- cumulative_lag0_lag7(matched_df, matched_cb)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="Hispanic")
write.csv(merged_df, "Results/Heatwaves/Lag7/Cumulative_Lag7_Hispanic.csv")
gc()

matched_df <- matched_df %>%
  filter(lag <= 3 & lag >= 0)
table(matched_df$lag)

matched_cb <- matched_cb_lag3(data)
merged_df <- cumulative_lag0_lag3(matched_df, matched_cb)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="Hispanic")
write.csv(merged_df, "Results/Heatwaves/Lag3/Cumulative_Lag3_Hispanic.csv")
gc()

#### WHITE ####

delivzip <- read_parquet("Data/Outcomes_SMM_v2.parquet")%>%
  filter(Race=="White") %>%
  group_by(zip, date) %>%
  summarize(n=n(),
            SMM=sum(SMM21))

matched_df <- read_parquet("Data/matched_df_dates.parquet")%>%
  left_join(delivzip, by=c('zip', 'date'))
matched_df[is.na(matched_df)] <- 0

data <- read_parquet("Data/Heatwave_Metrics_v2.parquet")%>%
  left_join(delivzip, by=c('zip', 'date'))

data[is.na(data)] <- 0

matched_cb <- matched_cb_lagn2_lag7(data)

outcomes <- c("SMM")

combined_data <- get_daily_lags(matched_df, matched_cb)
write.csv(combined_data, file = "Results/Heatwaves/Daily Lags/combined_matrices_White.csv", row.names = FALSE)
gc()

merged_df <- cumulative_lagneg2_lag7(matched_df, matched_cb)
merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="White")
write.csv(merged_df, "Results/Heatwaves/Lagn2/Cumulative_Lagn2_Lag7_White.csv")
gc()

matched_df <- matched_df %>%
  filter(lag <= 7 & lag >= 0)
table(matched_df$lag)

matched_cb <- matched_cb_lag7(data)
merged_df <- cumulative_lag0_lag7(matched_df, matched_cb)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="White")
write.csv(merged_df, "Results/Heatwaves/Lag7/Cumulative_Lag7_White.csv")
gc()

matched_df <- matched_df %>%
  filter(lag <= 3 & lag >= 0)
table(matched_df$lag)

matched_cb <- matched_cb_lag3(data)
merged_df <- cumulative_lag0_lag3(matched_df, matched_cb)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="White")
write.csv(merged_df, "Results/Heatwaves/Lag3/Cumulative_Lag3_White.csv")
gc()


#### AGE ####

delivzip <- read_parquet("Data/Outcomes_SMM_v2.parquet")%>%
  filter(Age=="35-39" | Age=="40+") %>%
  group_by(zip, date) %>%
  summarize(n=n(),
            SMM=sum(SMM21))

data <- read_parquet("Data/Heatwave_Metrics_v2.parquet")%>%
  left_join(delivzip, by=c('zip', 'date'))

matched_df <- read_parquet("Data/matched_df_dates.parquet")%>%
  left_join(delivzip, by=c('zip', 'date'))
matched_df[is.na(matched_df)]<-0

data[is.na(data)] <- 0

matched_cb <- matched_cb_lagn2_lag7(data)

outcomes <- c("SMM")

combined_data <- get_daily_lags(matched_df, matched_cb)
write.csv(combined_data, file = "Results/Heatwaves/Daily Lags/combined_matrices_Age.csv", row.names = FALSE)
gc()

merged_df <- cumulative_lagneg2_lag7(matched_df, matched_cb)
merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="Age")
write.csv(merged_df, "Results/Heatwaves/Lagn2/Cumulative_Lagn2_Lag7_Age.csv")
gc()

matched_df <- matched_df %>%
  filter(lag <= 7 & lag >= 0)
table(matched_df$lag)

matched_cb <- matched_cb_lag7(data)
merged_df <- cumulative_lag0_lag7(matched_df, matched_cb)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="Age")
write.csv(merged_df, "Results/Heatwaves/Lag7/Cumulative_Lag7_Age.csv")
gc()

matched_df <- matched_df %>%
  filter(lag <= 3 & lag >= 0)
table(matched_df$lag)

matched_cb <- matched_cb_lag3(data)
merged_df <- cumulative_lag0_lag3(matched_df, matched_cb)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="Age")
write.csv(merged_df, "Results/Heatwaves/Lag3/Cumulative_Lag3_Age.csv")
gc()

#### MEDICAID ####

delivzip <- read_parquet("Data/Outcomes_SMM_v2.parquet")%>%
  filter(Insurance=="Medicaid") %>%
  group_by(zip, date) %>%
  summarize(n=n(),
            SMM=sum(SMM21))

matched_df <- read_parquet("Data/matched_df_dates.parquet")%>%
  left_join(delivzip, by=c('zip', 'date'))
matched_df[is.na(matched_df)] <- 0

data <- read_parquet("Data/Heatwave_Metrics_v2.parquet")%>%
  left_join(delivzip, by=c('zip', 'date'))

data[is.na(data)] <- 0

matched_cb <- matched_cb_lagn2_lag7(data)

outcomes <- c("SMM")

combined_data <- get_daily_lags(matched_df, matched_cb)
write.csv(combined_data, file = "Results/Heatwaves/Daily Lags/combined_matrices_Medicaid.csv", row.names = FALSE)
gc()

merged_df <- cumulative_lagneg2_lag7(matched_df, matched_cb)
merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="Medicaid")
write.csv(merged_df, "Results/Heatwaves/Lagn2/Cumulative_Lagn2_Lag7_Medicaid.csv")
gc()

matched_df <- matched_df %>%
  filter(lag <= 7 & lag >= 0)
table(matched_df$lag)

matched_cb <- matched_cb_lag7(data)
merged_df <- cumulative_lag0_lag7(matched_df, matched_cb)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="Medicaid")
write.csv(merged_df, "Results/Heatwaves/Lag7/Cumulative_Lag7_Medicaid.csv")
gc()

matched_df <- matched_df %>%
  filter(lag <= 3 & lag >= 0)
table(matched_df$lag)

matched_cb <- matched_cb_lag3(data)
merged_df <- cumulative_lag0_lag3(matched_df, matched_cb)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="Medicaid")
write.csv(merged_df, "Results/Heatwaves/Lag3/Cumulative_Lag3_Medicaid.csv")
gc()

#### COMMERCIAL/PRIVATE INS ####

delivzip <- read_parquet("Data/Outcomes_SMM_v2.parquet")%>%
  filter(Insurance=="Commercial") %>%
  group_by(zip, date) %>%
  summarize(n=n(),
            SMM=sum(SMM21))

matched_df <- read_parquet("Data/matched_df_dates.parquet")%>%
  left_join(delivzip, by=c('zip', 'date'))
matched_df[is.na(matched_df)] <- 0

data <- read_parquet("Data/Heatwave_Metrics_v2.parquet")%>%
  left_join(delivzip, by=c('zip', 'date'))

data[is.na(data)] <- 0

matched_cb <- matched_cb_lagn2_lag7(data)

outcomes <- c("SMM")

combined_data <- get_daily_lags(matched_df, matched_cb)
write.csv(combined_data, file = "Results/Heatwaves/Daily Lags/combined_matrices_PrivateIns.csv", row.names = FALSE)
gc()

merged_df <- cumulative_lagneg2_lag7(matched_df, matched_cb)
merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="PrivateIns")
write.csv(merged_df, "Results/Heatwaves/Lagn2/Cumulative_Lagn2_Lag7_PrivateIns.csv")
gc()

matched_df <- matched_df %>%
  filter(lag <= 7 & lag >= 0)
table(matched_df$lag)

matched_cb <- matched_cb_lag7(data)
merged_df <- cumulative_lag0_lag7(matched_df, matched_cb)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="PrivateIns")
write.csv(merged_df, "Results/Heatwaves/Lag7/Cumulative_Lag7_PrivateIns.csv")
gc()

matched_df <- matched_df %>%
  filter(lag <= 3 & lag >= 0)
table(matched_df$lag)

matched_cb <- matched_cb_lag3(data)
merged_df <- cumulative_lag0_lag3(matched_df, matched_cb)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="PrivateIns")
write.csv(merged_df, "Results/Heatwaves/Lag3/Cumulative_Lag3_PrivateIns.csv")
gc()


#### OTHER RACE ####

delivzip <- read_parquet("Data/Outcomes_SMM_v2.parquet")%>%
  filter(Race=="Other") %>%
  group_by(zip, date) %>%
  summarize(n=n(),
            SMM=sum(SMM21))

matched_df <- read_parquet("Data/matched_df_dates.parquet")%>%
  left_join(delivzip, by=c('zip', 'date'))
matched_df[is.na(matched_df)] <- 0

data <- read_parquet("Data/Heatwave_Metrics_v2.parquet")%>%
  left_join(delivzip, by=c('zip', 'date'))

data[is.na(data)] <- 0

matched_cb <- matched_cb_lagn2_lag7(data)

outcomes <- c("SMM")

combined_data <- get_daily_lags(matched_df, matched_cb)
write.csv(combined_data, file = "Results/Heatwaves/Daily Lags/combined_matrices_OtherRace.csv", row.names = FALSE)
gc()

merged_df <- cumulative_lagneg2_lag7(matched_df, matched_cb)
merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="OtherRace")
write.csv(merged_df, "Results/Heatwaves/Lagn2/Cumulative_Lagn2_Lag7_OtherRace.csv")
gc()

matched_df <- matched_df %>%
  filter(lag <= 7 & lag >= 0)
table(matched_df$lag)

matched_cb <- matched_cb_lag7(data)
merged_df <- cumulative_lag0_lag7(matched_df, matched_cb)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="OtherRace")
write.csv(merged_df, "Results/Heatwaves/Lag7/Cumulative_Lag7_OtherRace.csv")
gc()

matched_df <- matched_df %>%
  filter(lag <= 3 & lag >= 0)
table(matched_df$lag)

matched_cb <- matched_cb_lag3(data)
merged_df <- cumulative_lag0_lag3(matched_df, matched_cb)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="OtherRace")
write.csv(merged_df, "Results/Heatwaves/Lag3/Cumulative_Lag3_OtherRace.csv")
gc()

#### SELF-PAY #####

delivzip <- read_parquet("Data/Outcomes_SMM_v2.parquet")%>%
  filter(Insurance=="Self-pay") %>%
  group_by(zip, date) %>%
  summarize(n=n(),
            SMM=sum(SMM21))

matched_df <- read_parquet("Data/matched_df_dates.parquet")%>%
  left_join(delivzip, by=c('zip', 'date'))
matched_df[is.na(matched_df)] <- 0

data <- read_parquet("Data/Heatwave_Metrics_v2.parquet")%>%
  left_join(delivzip, by=c('zip', 'date'))

data[is.na(data)] <- 0

matched_cb <- matched_cb_lagn2_lag7(data)

outcomes <- c("SMM")

combined_data <- get_daily_lags(matched_df, matched_cb)
write.csv(combined_data, file = "Results/Heatwaves/Daily Lags/combined_matrices_SelfPay.csv", row.names = FALSE)
gc()

merged_df <- cumulative_lagneg2_lag7(matched_df, matched_cb)
merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="SelfPay")
write.csv(merged_df, "Results/Heatwaves/Lagn2/Cumulative_Lagn2_Lag7_SelfPay.csv")
gc()

matched_df <- matched_df %>%
  filter(lag <= 7 & lag >= 0)
table(matched_df$lag)

matched_cb <- matched_cb_lag7(data)
merged_df <- cumulative_lag0_lag7(matched_df, matched_cb)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="SelfPay")
write.csv(merged_df, "Results/Heatwaves/Lag7/Cumulative_Lag7_SelfPay.csv")
gc()

matched_df <- matched_df %>%
  filter(lag <= 3 & lag >= 0)
table(matched_df$lag)

matched_cb <- matched_cb_lag3(data)
merged_df <- cumulative_lag0_lag3(matched_df, matched_cb)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="SelfPay")
write.csv(merged_df, "Results/Heatwaves/Lag3/Cumulative_Lag3_SelfPay.csv")
gc()
#### OTHER INS ####

delivzip <- read_parquet("Data/Outcomes_SMM_v2.parquet")%>%
  filter(Insurance=="Other gov't" | Insurance=="Other/Unknown" ) %>%
  group_by(zip, date) %>%
  summarize(n=n(),
            SMM=sum(SMM21))

matched_df <- read_parquet("Data/matched_df_dates.parquet")%>%
  left_join(delivzip, by=c('zip', 'date'))
matched_df[is.na(matched_df)] <- 0

data <- read_parquet("Data/Heatwave_Metrics_v2.parquet")%>%
  left_join(delivzip, by=c('zip', 'date'))

data[is.na(data)] <- 0

matched_cb <- matched_cb_lagn2_lag7(data)

outcomes <- c("SMM")

combined_data <- get_daily_lags(matched_df, matched_cb)
write.csv(combined_data, file = "Results/Heatwaves/Daily Lags/combined_matrices_OtherIns.csv", row.names = FALSE)
gc()

merged_df <- cumulative_lagneg2_lag7(matched_df, matched_cb)
merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="OtherIns")
write.csv(merged_df, "Results/Heatwaves/Lagn2/Cumulative_Lagn2_Lag7_OtherIns.csv")
gc()

matched_df <- matched_df %>%
  filter(lag <= 7 & lag >= 0)
table(matched_df$lag)

matched_cb <- matched_cb_lag7(data)
merged_df <- cumulative_lag0_lag7(matched_df, matched_cb)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="OtherIns")
write.csv(merged_df, "Results/Heatwaves/Lag7/Cumulative_Lag7_OtherIns.csv")
gc()

matched_df <- matched_df %>%
  filter(lag <= 3 & lag >= 0)
table(matched_df$lag)

matched_cb <- matched_cb_lag3(data)
merged_df <- cumulative_lag0_lag3(matched_df, matched_cb)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="OtherIns")
write.csv(merged_df, "Results/Heatwaves/Lag3/Cumulative_Lag3_OtherIns.csv")
gc()

#### AGE < 35 ####

delivzip <- read_parquet("Data/Outcomes_SMM_v2.parquet")%>%
  filter(Age != "35-39" & Age != "40+") %>%
  group_by(zip, date) %>%
  summarize(n=n(),
            SMM=sum(SMM21))

data <- read_parquet("Data/Heatwave_Metrics_v2.parquet")%>%
  left_join(delivzip, by=c('zip', 'date'))

matched_df <- read_parquet("Data/matched_df_dates.parquet")%>%
  left_join(delivzip, by=c('zip', 'date'))
matched_df[is.na(matched_df)]<-0

data[is.na(data)] <- 0

matched_cb <- matched_cb_lagn2_lag7(data)

outcomes <- c("SMM")

combined_data <- get_daily_lags(matched_df, matched_cb)
write.csv(combined_data, file = "Results/Heatwaves/Daily Lags/combined_matrices_AgeUnder35.csv", row.names = FALSE)
gc()

merged_df <- cumulative_lagneg2_lag7(matched_df, matched_cb)
merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="AgeUnder35")
write.csv(merged_df, "Results/Heatwaves/Lagn2/Cumulative_Lagn2_Lag7_AgeUnder35.csv")
gc()

matched_df <- matched_df %>%
  filter(lag <= 7 & lag >= 0)
table(matched_df$lag)

matched_cb <- matched_cb_lag7(data)
merged_df <- cumulative_lag0_lag7(matched_df, matched_cb)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="AgeUnder35")
write.csv(merged_df, "Results/Heatwaves/Lag7/Cumulative_Lag7_AgeUnder35.csv")
gc()

matched_df <- matched_df %>%
  filter(lag <= 3 & lag >= 0)
table(matched_df$lag)

matched_cb <- matched_cb_lag3(data)
merged_df <- cumulative_lag0_lag3(matched_df, matched_cb)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="AgeUnder35")
write.csv(merged_df, "Results/Heatwaves/Lag3/Cumulative_Lag3_AgeUnder35.csv")
gc()


#### NON-HISPANIC ####

delivzip <- read_parquet("Data/Outcomes_SMM_v2.parquet")%>%
  filter(Ethnicity=="Not Hispanic") %>%
  group_by(zip, date) %>%
  summarize(n=n(),
            SMM=sum(SMM21))


data <- read_parquet("Data/Heatwave_Metrics_v2.parquet")%>%
  left_join(delivzip, by=c('zip', 'date'))

matched_df <- read_parquet("Data/matched_df_dates.parquet")%>%
  left_join(delivzip, by=c('zip', 'date'))
matched_df[is.na(matched_df)]<-0

data[is.na(data)] <- 0

matched_cb <- matched_cb_lagn2_lag7(data)

outcomes <- c("SMM")

combined_data <- get_daily_lags(matched_df, matched_cb)
write.csv(combined_data, file = "Results/Heatwaves/Daily Lags/combined_matrices_NonHispanic.csv", row.names = FALSE)
gc()

merged_df <- cumulative_lagneg2_lag7(matched_df, matched_cb)
merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="NonHispanic")
write.csv(merged_df, "Results/Heatwaves/Lagn2/Cumulative_Lagn2_Lag7_NonHispanic.csv")
gc()

matched_df <- matched_df %>%
  filter(lag <= 7 & lag >= 0)
table(matched_df$lag)

matched_cb <- matched_cb_lag7(data)
merged_df <- cumulative_lag0_lag7(matched_df, matched_cb)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="NonHispanic")
write.csv(merged_df, "Results/Heatwaves/Lag7/Cumulative_Lag7_NonHispanic.csv")
gc()

matched_df <- matched_df %>%
  filter(lag <= 3 & lag >= 0)
table(matched_df$lag)

matched_cb <- matched_cb_lag3(data)
merged_df <- cumulative_lag0_lag3(matched_df, matched_cb)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="NonHispanic")
write.csv(merged_df, "Results/Heatwaves/Lag3/Cumulative_Lag3_NonHispanic.csv")
gc()


#### ICE INCOME LOW ####

delivzip <- read_parquet("Data/Outcomes_SMM_v2.parquet")%>%
  group_by(zip, date) %>%
  summarize(n=n(),
            SMM=sum(SMM21))

ice <- read.csv("Data/NC_ZCTA_ICE.csv") %>%
  mutate(zip = as.numeric(ZCTA),
         ICE_Income = as.numeric(AB_ICE)) %>%
  dplyr::select(zip, ICE_Income) %>%
  filter(!is.na(ICE_Income)) %>%
  mutate(Tertile = ntile(ICE_Income, 3)) %>%
  mutate(ICE_Income_Tertile = case_when(
    Tertile == 1 ~ "Low",
    Tertile == 2 ~ "Mid",
    Tertile == 3 ~ "High"
  )) %>%
  dplyr::select(-Tertile)  # Optional: Remove the numeric tertile indicator if not needed

data <- read_parquet("Data/Heatwave_Metrics_v2.parquet")%>%
  left_join(delivzip, by=c('zip', 'date'))%>%
  left_join(ice, by=c('zip'))%>%
  filter(ICE_Income_Tertile=="Low")

matched_df <- read_parquet("Data/matched_df_dates.parquet") %>%
  mutate(date=as.Date(date))%>%
  left_join(ice, by=c('zip'))%>%
  filter(ICE_Income_Tertile=="Low")%>%
  dplyr::select(-n)%>%
  left_join(delivzip, by=c('zip', 'date'))
matched_df[is.na(matched_df)] <- 0

data[is.na(data)] <- 0

matched_cb <- matched_cb_lagn2_lag7(data)

outcomes <- c("SMM")

combined_data <- get_daily_lags(matched_df, matched_cb)
write.csv(combined_data, file = "Results/Heatwaves/Daily Lags/combined_matrices_ICEIncome_Low.csv", row.names = FALSE)
gc()

merged_df <- cumulative_lagneg2_lag7(matched_df, matched_cb)
merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="ICE Income: Low")
write.csv(merged_df, "Results/Heatwaves/Lagn2/Cumulative_Lagn2_Lag7_ICEIncome_Low.csv")
gc()

matched_df <- matched_df %>%
  filter(lag <= 7 & lag >= 0)
table(matched_df$lag)

matched_cb <- matched_cb_lag7(data)
merged_df <- cumulative_lag0_lag7(matched_df, matched_cb)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="ICEIncome: Low")
write.csv(merged_df, "Results/Heatwaves/Lag7/Cumulative_Lag7_ICEIncome_Low.csv")
gc()

matched_df <- matched_df %>%
  filter(lag <= 3 & lag >= 0)
table(matched_df$lag)

matched_cb <- matched_cb_lag3(data)
merged_df <- cumulative_lag0_lag3(matched_df, matched_cb)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="ICEIncome: Low")
write.csv(merged_df, "Results/Heatwaves/Lag3/Cumulative_Lag3_ICEIncome_Low.csv")
gc()

#### ICE INCOME MID ####

delivzip <- read_parquet("Data/Outcomes_SMM_v2.parquet")%>%
  group_by(zip, date) %>%
  summarize(n=n(),
            SMM=sum(SMM21))

ice <- read.csv("Data/NC_ZCTA_ICE.csv") %>%
  mutate(zip = as.numeric(ZCTA),
         ICE_Income = as.numeric(AB_ICE)) %>%
  dplyr::select(zip, ICE_Income) %>%
  filter(!is.na(ICE_Income)) %>%
  mutate(Tertile = ntile(ICE_Income, 3)) %>%
  mutate(ICE_Income_Tertile = case_when(
    Tertile == 1 ~ "Low",
    Tertile == 2 ~ "Mid",
    Tertile == 3 ~ "High"
  )) %>%
  dplyr::select(-Tertile)  # Optional: Remove the numeric tertile indicator if not needed

data <- read_parquet("Data/Heatwave_Metrics_v2.parquet")%>%
  left_join(delivzip, by=c('zip', 'date'))%>%
  left_join(ice, by=c('zip'))%>%
  filter(ICE_Income_Tertile=="Mid")

matched_df <- read_parquet("Data/matched_df_dates.parquet") %>%
  mutate(date=as.Date(date))%>%
  left_join(ice, by=c('zip'))%>%
  filter(ICE_Income_Tertile=="Mid")%>%
  dplyr::select(-n)%>%
  left_join(delivzip, by=c('zip', 'date'))
matched_df[is.na(matched_df)] <- 0

data[is.na(data)] <- 0

matched_cb <- matched_cb_lagn2_lag7(data)

outcomes <- c("SMM")

combined_data <- get_daily_lags(matched_df, matched_cb)
write.csv(combined_data, file = "Results/Heatwaves/Daily Lags/combined_matrices_ICEIncome_Mid.csv", row.names = FALSE)
gc()

merged_df <- cumulative_lagneg2_lag7(matched_df, matched_cb)
merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="ICE Income: Mid")
write.csv(merged_df, "Results/Heatwaves/Lagn2/Cumulative_Lagn2_Lag7_ICEIncome_Mid.csv")
gc()

matched_df <- matched_df %>%
  filter(lag <= 7 & lag >= 0)
table(matched_df$lag)

matched_cb <- matched_cb_lag7(data)
merged_df <- cumulative_lag0_lag7(matched_df, matched_cb)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="ICEIncome: Mid")
write.csv(merged_df, "Results/Heatwaves/Lag7/Cumulative_Lag7_ICEIncome_Mid.csv")
gc()

matched_df <- matched_df %>%
  filter(lag <= 3 & lag >= 0)
table(matched_df$lag)

matched_cb <- matched_cb_lag3(data)
merged_df <- cumulative_lag0_lag3(matched_df, matched_cb)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="ICEIncome: Mid")
write.csv(merged_df, "Results/Heatwaves/Lag3/Cumulative_Lag3_ICEIncome_Mid.csv")
gc()

#### ICE INCOME HIGH ####

delivzip <- read_parquet("Data/Outcomes_SMM_v2.parquet")%>%
  group_by(zip, date) %>%
  summarize(n=n(),
            SMM=sum(SMM21))

ice <- read.csv("Data/NC_ZCTA_ICE.csv") %>%
  mutate(zip = as.numeric(ZCTA),
         ICE_Income = as.numeric(AB_ICE)) %>%
  dplyr::select(zip, ICE_Income) %>%
  filter(!is.na(ICE_Income)) %>%
  mutate(Tertile = ntile(ICE_Income, 3)) %>%
  mutate(ICE_Income_Tertile = case_when(
    Tertile == 1 ~ "Low",
    Tertile == 2 ~ "Mid",
    Tertile == 3 ~ "High"
  )) %>%
  dplyr::select(-Tertile)  # Optional: Remove the numeric tertile indicator if not needed

data <- read_parquet("Data/Heatwave_Metrics_v2.parquet")%>%
  left_join(delivzip, by=c('zip', 'date'))%>%
  left_join(ice, by=c('zip'))%>%
  filter(ICE_Income_Tertile=="High")

matched_df <- read_parquet("Data/matched_df_dates.parquet") %>%
  mutate(date=as.Date(date))%>%
  left_join(ice, by=c('zip'))%>%
  filter(ICE_Income_Tertile=="High")%>%
  dplyr::select(-n)%>%
  left_join(delivzip, by=c('zip', 'date'))
matched_df[is.na(matched_df)] <- 0

data[is.na(data)] <- 0

matched_cb <- matched_cb_lagn2_lag7(data)

outcomes <- c("SMM")

combined_data <- get_daily_lags(matched_df, matched_cb)
write.csv(combined_data, file = "Results/Heatwaves/Daily Lags/combined_matrices_ICEIncome_High.csv", row.names = FALSE)
gc()

merged_df <- cumulative_lagneg2_lag7(matched_df, matched_cb)
merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="ICE Income: High")
write.csv(merged_df, "Results/Heatwaves/Lagn2/Cumulative_Lagn2_Lag7_ICEIncome_High.csv")
gc()

matched_df <- matched_df %>%
  filter(lag <= 7 & lag >= 0)
table(matched_df$lag)

matched_cb <- matched_cb_lag7(data)
merged_df <- cumulative_lag0_lag7(matched_df, matched_cb)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="ICE Income: High")
write.csv(merged_df, "Results/Heatwaves/Lag7/Cumulative_Lag7_ICEIncome_High.csv")
gc()

matched_df <- matched_df %>%
  filter(lag <= 3 & lag >= 0)
table(matched_df$lag)

matched_cb <- matched_cb_lag3(data)
merged_df <- cumulative_lag0_lag3(matched_df, matched_cb)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="ICE Income: High")
write.csv(merged_df, "Results/Heatwaves/Lag3/Cumulative_Lag3_ICEIncome_High.csv")
gc()

#### ICE RACE LOW ####

delivzip <- read_parquet("Data/Outcomes_SMM_v2.parquet")%>%
  group_by(zip, date) %>%
  summarize(n=n(),
            SMM=sum(SMM21))

ice <- read.csv("Data/NC_ZCTA_ICE.csv") %>%
  mutate(zip = as.numeric(ZCTA),
         ICE_Race = as.numeric(CD_ICE)) %>%
  dplyr::select(zip, ICE_Race) %>%
  filter(!is.na(ICE_Race)) %>%
  mutate(Tertile = ntile(ICE_Race, 3)) %>%
  mutate(ICE_Race_Tertile = case_when(
    Tertile == 1 ~ "Low",
    Tertile == 2 ~ "Mid",
    Tertile == 3 ~ "High"
  )) %>%
  dplyr::select(-Tertile)  # Optional: Remove the numeric tertile indicator if not needed

data <- read_parquet("Data/Heatwave_Metrics_v2.parquet")%>%
  left_join(delivzip, by=c('zip', 'date'))%>%
  left_join(ice, by=c('zip'))%>%
  filter(ICE_Race_Tertile=="Low")

matched_df <- read_parquet("Data/matched_df_dates.parquet") %>%
  mutate(date=as.Date(date))%>%
  left_join(ice, by=c('zip'))%>%
  filter(ICE_Race_Tertile=="Low")%>%
  dplyr::select(-n)%>%
  left_join(delivzip, by=c('zip', 'date'))
matched_df[is.na(matched_df)] <- 0

data[is.na(data)] <- 0

matched_cb <- matched_cb_lagn2_lag7(data)

outcomes <- c("SMM")

combined_data <- get_daily_lags(matched_df, matched_cb)
write.csv(combined_data, file = "Results/Heatwaves/Daily Lags/combined_matrices_ICERace_Low.csv", row.names = FALSE)
gc()

merged_df <- cumulative_lagneg2_lag7(matched_df, matched_cb)
merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="ICE Race: Low")
write.csv(merged_df, "Results/Heatwaves/Lagn2/Cumulative_Lagn2_Lag7_ICERace_Low.csv")
gc()

matched_df <- matched_df %>%
  filter(lag <= 7 & lag >= 0)
table(matched_df$lag)

matched_cb <- matched_cb_lag7(data)
merged_df <- cumulative_lag0_lag7(matched_df, matched_cb)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="ICERace: Low")
write.csv(merged_df, "Results/Heatwaves/Lag7/Cumulative_Lag7_ICERace_Low.csv")
gc()

matched_df <- matched_df %>%
  filter(lag <= 3 & lag >= 0)
table(matched_df$lag)

matched_cb <- matched_cb_lag3(data)
merged_df <- cumulative_lag0_lag3(matched_df, matched_cb)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="ICERace: Low")
write.csv(merged_df, "Results/Heatwaves/Lag3/Cumulative_Lag3_ICERace_Low.csv")
gc()

#### ICE RACE MID ####

delivzip <- read_parquet("Data/Outcomes_SMM_v2.parquet")%>%
  group_by(zip, date) %>%
  summarize(n=n(),
            SMM=sum(SMM21))

ice <- read.csv("Data/NC_ZCTA_ICE.csv") %>%
  mutate(zip = as.numeric(ZCTA),
         ICE_Race = as.numeric(CD_ICE)) %>%
  dplyr::select(zip, ICE_Race) %>%
  filter(!is.na(ICE_Race)) %>%
  mutate(Tertile = ntile(ICE_Race, 3)) %>%
  mutate(ICE_Race_Tertile = case_when(
    Tertile == 1 ~ "Low",
    Tertile == 2 ~ "Mid",
    Tertile == 3 ~ "High"
  )) %>%
  dplyr::select(-Tertile)  # Optional: Remove the numeric tertile indicator if not needed

data <- read_parquet("Data/Heatwave_Metrics_v2.parquet")%>%
  left_join(delivzip, by=c('zip', 'date'))%>%
  left_join(ice, by=c('zip'))%>%
  filter(ICE_Race_Tertile=="Mid")

matched_df <- read_parquet("Data/matched_df_dates.parquet") %>%
  mutate(date=as.Date(date))%>%
  left_join(ice, by=c('zip'))%>%
  filter(ICE_Race_Tertile=="Mid")%>%
  dplyr::select(-n)%>%
  left_join(delivzip, by=c('zip', 'date'))
matched_df[is.na(matched_df)] <- 0

data[is.na(data)] <- 0

matched_cb <- matched_cb_lagn2_lag7(data)

outcomes <- c("SMM")

combined_data <- get_daily_lags(matched_df, matched_cb)
write.csv(combined_data, file = "Results/Heatwaves/Daily Lags/combined_matrices_ICERace_Mid.csv", row.names = FALSE)
gc()

merged_df <- cumulative_lagneg2_lag7(matched_df, matched_cb)
merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="ICE Race: Mid")
write.csv(merged_df, "Results/Heatwaves/Lagn2/Cumulative_Lagn2_Lag7_ICERace_Mid.csv")
gc()

matched_df <- matched_df %>%
  filter(lag <= 7 & lag >= 0)
table(matched_df$lag)

matched_cb <- matched_cb_lag7(data)
merged_df <- cumulative_lag0_lag7(matched_df, matched_cb)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="ICERace: Mid")
write.csv(merged_df, "Results/Heatwaves/Lag7/Cumulative_Lag7_ICERace_Mid.csv")
gc()

matched_df <- matched_df %>%
  filter(lag <= 3 & lag >= 0)
table(matched_df$lag)

matched_cb <- matched_cb_lag3(data)
merged_df <- cumulative_lag0_lag3(matched_df, matched_cb)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="ICERace: Mid")
write.csv(merged_df, "Results/Heatwaves/Lag3/Cumulative_Lag3_ICERace_Mid.csv")
gc()

#### ICE RACE HIGH ####

delivzip <- read_parquet("Data/Outcomes_SMM_v2.parquet")%>%
  group_by(zip, date) %>%
  summarize(n=n(),
            SMM=sum(SMM21))

ice <- read.csv("Data/NC_ZCTA_ICE.csv") %>%
  mutate(zip = as.numeric(ZCTA),
         ICE_Race = as.numeric(CD_ICE)) %>%
  dplyr::select(zip, ICE_Race) %>%
  filter(!is.na(ICE_Race)) %>%
  mutate(Tertile = ntile(ICE_Race, 3)) %>%
  mutate(ICE_Race_Tertile = case_when(
    Tertile == 1 ~ "Low",
    Tertile == 2 ~ "Mid",
    Tertile == 3 ~ "High"
  )) %>%
  dplyr::select(-Tertile)  # Optional: Remove the numeric tertile indicator if not needed

data <- read_parquet("Data/Heatwave_Metrics_v2.parquet")%>%
  left_join(delivzip, by=c('zip', 'date'))%>%
  left_join(ice, by=c('zip'))%>%
  filter(ICE_Race_Tertile=="High")

matched_df <- read_parquet("Data/matched_df_dates.parquet") %>%
  mutate(date=as.Date(date))%>%
  left_join(ice, by=c('zip'))%>%
  filter(ICE_Race_Tertile=="High")%>%
  dplyr::select(-n)%>%
  left_join(delivzip, by=c('zip', 'date'))
matched_df[is.na(matched_df)] <- 0

data[is.na(data)] <- 0

matched_cb <- matched_cb_lagn2_lag7(data)

outcomes <- c("SMM")

combined_data <- get_daily_lags(matched_df, matched_cb)
write.csv(combined_data, file = "Results/Heatwaves/Daily Lags/combined_matrices_ICERace_High.csv", row.names = FALSE)
gc()

merged_df <- cumulative_lagneg2_lag7(matched_df, matched_cb)
merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="ICE Race: High")
write.csv(merged_df, "Results/Heatwaves/Lagn2/Cumulative_Lagn2_Lag7_ICERace_High.csv")
gc()

matched_df <- matched_df %>%
  filter(lag <= 7 & lag >= 0)
table(matched_df$lag)

matched_cb <- matched_cb_lag7(data)
merged_df <- cumulative_lag0_lag7(matched_df, matched_cb)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="ICE Race: High")
write.csv(merged_df, "Results/Heatwaves/Lag7/Cumulative_Lag7_ICERace_High.csv")
gc()

matched_df <- matched_df %>%
  filter(lag <= 3 & lag >= 0)
table(matched_df$lag)

matched_cb <- matched_cb_lag3(data)
merged_df <- cumulative_lag0_lag3(matched_df, matched_cb)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="ICE Race: High")
write.csv(merged_df, "Results/Heatwaves/Lag3/Cumulative_Lag3_ICERace_High.csv")
gc()

#### URBAN ####

delivzip <- read_parquet("Data/Outcomes_SMM_v2.parquet")%>%
  group_by(zip, date) %>%
  summarize(n=n(),
            SMM=sum(SMM21))

ruca <- read.csv("Data/RUCA.csv")%>%
  rename(zip=ZCTA)

data <- read_parquet("Data/Heatwave_Metrics_v2.parquet")%>%
  left_join(delivzip, by=c('zip', 'date'))%>%
  left_join(ruca, by=c('zip'))%>%
  filter(RUCA10 <= 3)

matched_df <- read_parquet("Data/matched_df_dates.parquet") %>%
  mutate(date=as.Date(date))%>%
  left_join(ruca, by=c('zip'))%>%
  filter(RUCA10 <= 3)%>%
  dplyr::select(-n)%>%
  left_join(delivzip, by=c('zip', 'date'))
matched_df[is.na(matched_df)] <- 0

data[is.na(data)] <- 0

matched_cb <- matched_cb_lagn2_lag7(data)

outcomes <- c("SMM")

combined_data <- get_daily_lags(matched_df, matched_cb)
write.csv(combined_data, file = "Results/Heatwaves/Daily Lags/combined_matrices_Urban.csv", row.names = FALSE)
gc()

merged_df <- cumulative_lagneg2_lag7(matched_df, matched_cb)
merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="Urban")
write.csv(merged_df, "Results/Heatwaves/Lagn2/Cumulative_Lagn2_Lag7_Urban.csv")
gc()

matched_df <- matched_df %>%
  filter(lag <= 7 & lag >= 0)
table(matched_df$lag)

matched_cb <- matched_cb_lag7(data)
merged_df <- cumulative_lag0_lag7(matched_df, matched_cb)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="Urban")
write.csv(merged_df, "Results/Heatwaves/Lag7/Cumulative_Lag7_Urban.csv")
gc()

matched_df <- matched_df %>%
  filter(lag <= 3 & lag >= 0)
table(matched_df$lag)

matched_cb <- matched_cb_lag3(data)
merged_df <- cumulative_lag0_lag3(matched_df, matched_cb)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="Urban")
write.csv(merged_df, "Results/Heatwaves/Lag3/Cumulative_Lag3_Urban.csv")
gc()

#### RURAL ####

delivzip <- read_parquet("Data/Outcomes_SMM_v2.parquet")%>%
  group_by(zip, date) %>%
  summarize(n=n(),
            SMM=sum(SMM21))

ruca <- read.csv("Data/RUCA.csv")%>%
  rename(zip=ZCTA)

data <- read_parquet("Data/Heatwave_Metrics_v2.parquet")%>%
  left_join(delivzip, by=c('zip', 'date'))%>%
  left_join(ruca, by=c('zip'))%>%
  filter(RUCA10 >= 7)%>%
  dplyr::select(-n)%>%
  left_join(delivzip, by=c('zip', 'date'))
matched_df[is.na(matched_df)] <- 0

data[is.na(data)]<-0

matched_df <- read_parquet("Data/matched_df_dates.parquet") %>%
  mutate(date=as.Date(date))%>%
  left_join(ruca, by=c('zip'))%>%
  filter(RUCA10 >= 7)

data[is.na(data)] <- 0

matched_cb <- matched_cb_lagn2_lag7(data)

outcomes <- c("SMM")

combined_data <- get_daily_lags(matched_df, matched_cb)
write.csv(combined_data, file = "Results/Heatwaves/Daily Lags/combined_matrices_Rural.csv", row.names = FALSE)
gc()

merged_df <- cumulative_lagneg2_lag7(matched_df, matched_cb)
merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="Rural")
write.csv(merged_df, "Results/Heatwaves/Lagn2/Cumulative_Lagn2_Lag7_Rural.csv")
gc()

matched_df <- matched_df %>%
  filter(lag <= 7 & lag >= 0)
table(matched_df$lag)

matched_cb <- matched_cb_lag7(data)
merged_df <- cumulative_lag0_lag7(matched_df, matched_cb)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="Rural")
write.csv(merged_df, "Results/Heatwaves/Lag7/Cumulative_Lag7_Rural.csv")
gc()

matched_df <- matched_df %>%
  filter(lag <= 3 & lag >= 0)
table(matched_df$lag)

matched_cb <- matched_cb_lag3(data)
merged_df <- cumulative_lag0_lag3(matched_df, matched_cb)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="Rural")
write.csv(merged_df, "Results/Heatwaves/Lag3/Cumulative_Lag3_Rural.csv")
gc()

#### SUBURBAN ####

delivzip <- read_parquet("Data/Outcomes_SMM_v2.parquet")%>%
  group_by(zip, date) %>%
  summarize(n=n(),
            SMM=sum(SMM21))

ruca <- read.csv("Data/RUCA.csv")%>%
  rename(zip=ZCTA)

data <- read_parquet("Data/Heatwave_Metrics_v2.parquet")%>%
  left_join(delivzip, by=c('zip', 'date'))%>%
  left_join(ruca, by=c('zip'))%>%
  filter(RUCA10 > 3 & RUCA10 <7)

data <- read_parquet("Data/Heatwave_Metrics_v2.parquet")%>%
  left_join(delivzip, by=c('zip', 'date'))%>%
  left_join(ruca, by=c('zip'))%>%
  filter(RUCA10 > 3 & RUCA10 <7)%>%
  dplyr::select(-n)%>%
  left_join(delivzip, by=c('zip', 'date'))
matched_df[is.na(matched_df)] <- 0

data[is.na(data)] <- 0

matched_cb <- matched_cb_lagn2_lag7(data)

outcomes <- c("SMM")

combined_data <- get_daily_lags(matched_df, matched_cb)
write.csv(combined_data, file = "Results/Heatwaves/Daily Lags/combined_matrices_Suburban.csv", row.names = FALSE)
gc()

merged_df <- cumulative_lagneg2_lag7(matched_df, matched_cb)
merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="Suburban")
write.csv(merged_df, "Results/Heatwaves/Lagn2/Cumulative_Lagn2_Lag7_Suburban.csv")
gc()

matched_df <- matched_df %>%
  filter(lag <= 7 & lag >= 0)
table(matched_df$lag)

matched_cb <- matched_cb_lag7(data)
merged_df <- cumulative_lag0_lag7(matched_df, matched_cb)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="Suburban")
write.csv(merged_df, "Results/Heatwaves/Lag7/Cumulative_Lag7_Suburban.csv")
gc()

matched_df <- matched_df %>%
  filter(lag <= 3 & lag >= 0)
table(matched_df$lag)

matched_cb <- matched_cb_lag3(data)
merged_df <- cumulative_lag0_lag3(matched_df, matched_cb)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="Suburban")
write.csv(merged_df, "Results/Heatwaves/Lag3/Cumulative_Lag3_Suburban.csv")
gc()

#### WESTERN ####

delivzip <- read_parquet("Data/Outcomes_SMM_v2.parquet")%>%
  group_by(zip, date) %>%
  summarize(n=n(),
            SMM=sum(SMM21))

geo <- read.csv("Data/nc_zcta_geo_regions.csv")%>%
  dplyr::select(zip, geo_region)

data <- read_parquet("Data/Heatwave_Metrics_v2.parquet")%>%
  left_join(delivzip, by=c('zip', 'date'))%>%
  left_join(geo, by=c('zip'))%>%
  filter(geo_region=="Western")

matched_df <- read_parquet("Data/matched_df_dates.parquet") %>%
  mutate(date=as.Date(date))%>%
  left_join(geo, by=c('zip'))%>%
  filter(geo_region=="Western")%>%
  dplyr::select(-n)%>%
  left_join(delivzip, by=c('zip', 'date'))
matched_df[is.na(matched_df)] <- 0

data[is.na(data)] <- 0

matched_cb <- matched_cb_lagn2_lag7(data)

outcomes <- c("SMM")

combined_data <- get_daily_lags(matched_df, matched_cb)
write.csv(combined_data, file = "Results/Heatwaves/Daily Lags/combined_matrices_Western.csv", row.names = FALSE)
gc()

merged_df <- cumulative_lagneg2_lag7(matched_df, matched_cb)
merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="Western")
write.csv(merged_df, "Results/Heatwaves/Lagn2/Cumulative_Lagn2_Lag7_Western.csv")
gc()

matched_df <- matched_df %>%
  filter(lag <= 7 & lag >= 0)
table(matched_df$lag)

matched_cb <- matched_cb_lag7(data)
merged_df <- cumulative_lag0_lag7(matched_df, matched_cb)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="Western")
write.csv(merged_df, "Results/Heatwaves/Lag7/Cumulative_Lag7_Western.csv")
gc()

matched_df <- matched_df %>%
  filter(lag <= 3 & lag >= 0)
table(matched_df$lag)

matched_cb <- matched_cb_lag3(data)
merged_df <- cumulative_lag0_lag3(matched_df, matched_cb)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="Western")
write.csv(merged_df, "Results/Heatwaves/Lag3/Cumulative_Lag3_Western.csv")
gc()

#### COASTAL ####

delivzip <- read_parquet("Data/Outcomes_SMM_v2.parquet")%>%
  group_by(zip, date) %>%
  summarize(n=n(),
            SMM=sum(SMM21))

geo <- read.csv("Data/nc_zcta_geo_regions.csv")%>%
  dplyr::select(zip, geo_region)

data <- read_parquet("Data/Heatwave_Metrics_v2.parquet")%>%
  left_join(delivzip, by=c('zip', 'date'))%>%
  left_join(geo, by=c('zip'))%>%
  filter(geo_region=="Eastern")

matched_df <- read_parquet("Data/matched_df_dates.parquet") %>%
  mutate(date=as.Date(date))%>%
  left_join(geo, by=c('zip'))%>%
  filter(geo_region=="Eastern")%>%
  dplyr::select(-n)%>%
  left_join(delivzip, by=c('zip', 'date'))
matched_df[is.na(matched_df)] <- 0

data[is.na(data)] <- 0

matched_cb <- matched_cb_lagn2_lag7(data)

outcomes <- c("SMM")

combined_data <- get_daily_lags(matched_df, matched_cb)
write.csv(combined_data, file = "Results/Heatwaves/Daily Lags/combined_matrices_Coastal.csv", row.names = FALSE)
gc()

merged_df <- cumulative_lagneg2_lag7(matched_df, matched_cb)
merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="Coastal")
write.csv(merged_df, "Results/Heatwaves/Lagn2/Cumulative_Lagn2_Lag7_Coastal.csv")
gc()

matched_df <- matched_df %>%
  filter(lag <= 7 & lag >= 0)
table(matched_df$lag)

matched_cb <- matched_cb_lag7(data)
merged_df <- cumulative_lag0_lag7(matched_df, matched_cb)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="Coastal")
write.csv(merged_df, "Results/Heatwaves/Lag7/Cumulative_Lag7_Coastal.csv")
gc()

matched_df <- matched_df %>%
  filter(lag <= 3 & lag >= 0)
table(matched_df$lag)

matched_cb <- matched_cb_lag3(data)
merged_df <- cumulative_lag0_lag3(matched_df, matched_cb)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="Coastal")
write.csv(merged_df, "Results/Heatwaves/Lag3/Cumulative_Lag3_Coastal.csv")
gc()

#### PIEDMONT ####


delivzip <- read_parquet("Data/Outcomes_SMM_v2.parquet")%>%
  group_by(zip, date) %>%
  summarize(n=n(),
            SMM=sum(SMM21))

geo <- read.csv("Data/nc_zcta_geo_regions.csv")%>%
  dplyr::select(zip, geo_region)

data <- read_parquet("Data/Heatwave_Metrics_v2.parquet")%>%
  left_join(delivzip, by=c('zip', 'date'))%>%
  left_join(geo, by=c('zip'))%>%
  filter(geo_region=="Piedmont")

matched_df <- read_parquet("Data/matched_df_dates.parquet") %>%
  mutate(date=as.Date(date))%>%
  left_join(geo, by=c('zip'))%>%
  filter(geo_region=="Piedmont")%>%
dplyr::select(-n)%>%
  left_join(delivzip, by=c('zip', 'date'))
matched_df[is.na(matched_df)] <- 0

data[is.na(data)] <- 0

matched_cb <- matched_cb_lagn2_lag7(data)

outcomes <- c("SMM")

combined_data <- get_daily_lags(matched_df, matched_cb)
write.csv(combined_data, file = "Results/Heatwaves/Daily Lags/combined_matrices_Piedmont.csv", row.names = FALSE)
gc()

merged_df <- cumulative_lagneg2_lag7(matched_df, matched_cb)
merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="Piedmont")
write.csv(merged_df, "Results/Heatwaves/Lagn2/Cumulative_Lagn2_Lag7_Piedmont.csv")
gc()

matched_df <- matched_df %>%
  filter(lag <= 7 & lag >= 0)
table(matched_df$lag)

matched_cb <- matched_cb_lag7(data)
merged_df <- cumulative_lag0_lag7(matched_df, matched_cb)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="Piedmont")
write.csv(merged_df, "Results/Heatwaves/Lag7/Cumulative_Lag7_Piedmont.csv")
gc()

matched_df <- matched_df %>%
  filter(lag <= 3 & lag >= 0)
table(matched_df$lag)

matched_cb <- matched_cb_lag3(data)
merged_df <- cumulative_lag0_lag3(matched_df, matched_cb)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="Piedmont")
write.csv(merged_df, "Results/Heatwaves/Lag3/Cumulative_Lag3_Piedmont.csv")
gc()















