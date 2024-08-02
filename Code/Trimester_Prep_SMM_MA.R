
# Code to create trimesters for SMM cases

library(tableone)
library(gmodels)
library(dplyr)
library(lubridate)
library(arrow)
library(data.table)
library(dplyr)
library(geepack)
library(haven)
library(lubridate)
library(lme4)
library(tidyverse)
setwd("~/RStudio/Matched SMM/")

## Read in data files

# Hospital delivery data for North Carolina

df_wrong <- read_sas("Data/smm.sas7bdat")%>%
  mutate(Date=as.Date(admitdt))%>%
  filter(ptzip >= 27006 & ptzip <= 28909)%>%
  filter(year(Date) >= 2011 & year(Date) <= 2019)%>%
  mutate(Year=year(Date))
table(df_wrong$Year)
table(df_wrong$fyear)

# Code for age, race, ethnicity, and insurance categories

#### AGE ####

df_wrong$Age <- as.factor(ifelse(df_wrong$agey < 20, '18-19', 
                                 ifelse(df_wrong$agey<25, '20-24',
                                        ifelse(df_wrong$agey <30, '25-29',
                                               ifelse(df_wrong$agey<35,'30-34',
                                                      ifelse(df_wrong$agey<40, '35-39',
                                                             ifelse(df_wrong$agey <= 44, '40+', 0)))))))
df_wrong$Age <- factor(df_wrong$Age)

df_wrong$Race <- fifelse(df_wrong$Race == "5", "White",
                         fifelse(df_wrong$Race == "3", "Black",
                                 fifelse(df_wrong$Race %in% c("2", "1", "4", "6"), "Other",
                                         "Unknown")))
df_wrong$Race <- as.factor(df_wrong$Race)
table(df_wrong$Race)

df_wrong$Ethnicity <- fifelse(df_wrong$Ethnicity == "1", "Not Hispanic",
                              fifelse(df_wrong$Ethnicity == "2", "Hispanic",
                                      "Unknown"))
df_wrong$Ethnicity <- factor(df_wrong$Ethnicity)
table(df_wrong$Ethnicity)

df_wrong$Age <- as.factor(ifelse(df_wrong$agey < 20, '18-19', 
                                 ifelse(df_wrong$agey<25, '20-24',
                                        ifelse(df_wrong$agey <30, '25-29',
                                               ifelse(df_wrong$agey<35,'30-34',
                                                      ifelse(df_wrong$agey<40, '35-39',
                                                             ifelse(df_wrong$agey <= 44, '40+', 0)))))))
df_wrong$Age <- factor(df_wrong$Age)
table(df_wrong$Age)

df_wrong$Insurance <- as.factor(ifelse(df_wrong$payer1 == "09", "Self-pay",
                                       ifelse(df_wrong$payer1 == "MC", "Medicaid",
                                              ifelse(df_wrong$payer1 %in% c("MA", "MB", "OF", "VA", "TV", "11", "16", "CH", "MM"), "Other gov't", 
                                                     ifelse(df_wrong$payer1 %in% c("LM", "LI", "HM", "DS", "CI", "BL", "AM", "12", "13", "14", "15"),
                                                            "Commercial", "Other/Unknown")))))

df_wrong$Insurance <- factor(df_wrong$Insurance)
table(df_wrong$Insurance)

#### Get week of gestation using diag codes ####

# Use ICD codes to identify the weeks of gestation for each delivery; this will
# cut down sample sizes but there is no other indicator
# ICD-10 codes for gestation go from 39 to 8 weeks; anything under 8 is coded as 
# 8 weeks or less
# ICD-10 code for <8 weeks is "Z3A01"

for (week in 39:8) {
  
  code_to_check <- paste0("Z3A", week)
  
  check_code <- function(code) {
    any(grepl(code_to_check, gsub("\\.", "", code)))
  }
  
  # Create column name dynamically
  week_column <- paste0("Week_", week)
  
  # Apply the function across rows and store the result in the correct column
  df_wrong[[week_column]] <- apply(df_wrong[, paste0("diag", 1:25)], 1, function(x) as.numeric(any(check_code(x))))
}

# Create weeks_gest column to indicate weeks of gestation for each delivery

df_wrong <- df_wrong %>%
  mutate(
    weeks_gest = case_when(
      Week_39 == 1 ~ 39,
      Week_38 == 1 ~ 38,
      Week_37 == 1 ~ 37,
      Week_36 == 1 ~ 36,
      Week_35 == 1 ~ 35,
      Week_34 == 1 ~ 34,
      Week_33 == 1 ~ 33,
      Week_32 == 1 ~ 32,
      Week_31 == 1 ~ 31,
      Week_30 == 1 ~ 30,
      Week_29 == 1 ~ 29,
      Week_28 == 1 ~ 28,
      Week_27 == 1 ~ 27,
      Week_26 == 1 ~ 26,
      Week_25 == 1 ~ 25,
      Week_24 == 1 ~ 24,
      Week_23 == 1 ~ 23,
      Week_22 == 1 ~ 22,
      Week_21 == 1 ~ 21,
      Week_20 == 1 ~ 20,
      Week_19 == 1 ~ 19,
      Week_18 == 1 ~ 18,
      Week_17 == 1 ~ 17,
      Week_16 == 1 ~ 16,
      Week_15 == 1 ~ 15,
      Week_14 == 1 ~ 14,
      Week_13 == 1 ~ 13,
      Week_12 == 1 ~ 12,
      Week_11 == 1 ~ 11,
      Week_10 == 1 ~ 10,
      Week_9 == 1 ~ 9,
      Week_8 == 1 ~ 8,
     # Add more conditions for Week_35 through Week_0 as needed
      TRUE ~ NA_real_  # Default case (if none of the conditions are met)
    )
  )

# Check distribution and see how many NA values there are
# Usually a little more than half the sample is still NA which is fine

table(df_wrong$weeks_gest)
table(is.na(df_wrong$weeks_gest))

# Use gestation weeks to calculate trimester dates 

birth1 <- df_wrong %>%
  mutate(zip=as.numeric(ptzip))%>%
  filter(weeks_gest > 0)
birth1$GEST <- as.numeric(birth1$weeks_gest)
birth1$Pre <- as.Date(birth1$Date) - ((birth1$GEST * 7) + (13*7))
birth1$Tr1 <- as.Date(birth1$Date) - (birth1$GEST * 7)
birth1$Tr2 <- as.Date(birth1$Date) - ((birth1$GEST * 7) - (13*7))
birth1$Tr3 <- as.Date(birth1$Date) - ((birth1$GEST * 7) - (26*7))
birth1$Last_Gest_Week <- as.Date(birth1$Date) - 7
birth1$Tr3_end <- birth1$Date
birth1 <- birth1 %>%
  dplyr::select(zip, Year, Date, Shepsid, GEST, Date, Pre, Tr1, Tr2, Tr3, Last_Gest_Week, Tr3_end)

# Read in heatwave file

heat <- read_parquet("Data/Heatwave_Metrics_v2.parquet")%>%
  rename(Date=date)%>%
  mutate(zip=as.numeric(zip))

# Read in coldwave file

cold <- read_parquet("Data/Coldwave_Metrics_v2.parquet")%>%
  rename(Date=date)%>%
  mutate(zip=as.numeric(zip))

create_days <- function(birth_data){
  
  birth_data <- birth_data %>%
    data.frame()
  
  # Calculate number of days between start and end date
  birth_data$duration <- as.numeric(birth_data$Tr3_end - birth_data$Pre) + 1
  
  start_day <- birth_data$Pre
  end_day <- birth_data$Tr3_end
  num_days <- birth_data$duration
  column_names <- paste0("Day_", 1:max(num_days))
  
  # Create a list of date ranges, padding shorter ranges with NA values
  date_ranges <- mapply(function(start, end, days) {
    seq(start, end, by = "day")[1:days]
  }, start_day, end_day, num_days)
  
  # Pad shorter date ranges with NA values
  max_length <- max(num_days)
  date_ranges <- lapply(date_ranges, function(x) {
    if (length(x) < max_length) {
      c(x, rep(NA, max_length - length(x)))
    } else {
      x
    }
  })
  
  birth_data[column_names] <- as.data.frame(do.call(rbind, date_ranges))
  
  # Find columns starting with "Day_"
  day_columns <- grep("^Day_", names(birth_data), value = TRUE)
  
  # Convert numeric columns to date format
  birth_data[, day_columns] <- lapply(birth_data[, day_columns], function(x) as.Date(as.numeric(x), format="%Y-%m-%d"))
  
  birth_data <- birth_data %>%
    dplyr::select(-duration)
  
  return(birth_data)
}

birth2 <- create_days(birth1)

birth_Tr3 <- birth2[,-c(12:294)] # Trimester 3
birth_Tr2 <- birth2[,-c(12:194, 284:376)] # Trimester 2
birth_Tr1 <- birth2[,-c(12: 102, 194:376)] # Trimester 1

# Number of heatwaves experienced during trimester 3

birth3 <- birth_Tr3 %>%
  mutate(row = row_number())%>%
  dplyr::select(starts_with("Day"), zip, row, Shepsid) %>%
  tidyr::pivot_longer(c(-"zip", -"row", -"Shepsid")) %>%
  group_by(zip)%>%
  mutate(Date=as.Date(value))%>%
  filter(month(Date) >= 5 & month(Date) <= 9)%>%
  left_join(heat, by = c('zip', 'Date'))%>% 
  pivot_wider(names_from = name)%>%
  data.frame() 

Tr3_grp_hw <- birth3 %>%
  group_by(Shepsid)%>%
  summarize(Tr3_hw_count=sum(heatwave))

# Number of heatwaves experienced during trimester 2

birth3 <- birth_Tr2 %>%
  mutate(row = row_number())%>%
  select(starts_with("Day"), zip, row, Shepsid) %>%
  tidyr::pivot_longer(c(-"zip", -"row", -"Shepsid")) %>%
  group_by(zip)%>%
  mutate(Date=as.Date(value))%>%
  filter(month(Date) >= 5 & month(Date) <= 9)%>%
  left_join(heat, by = c('zip', 'Date'))%>% 
  pivot_wider(names_from = name)%>%
  data.frame() 

Tr2_grp_hw <- birth3 %>%
  group_by(Shepsid)%>%
  summarize(Tr2_hw_count=sum(heatwave))

# Number of heatwaves experienced during trimester 1

birth3 <- birth_Tr1 %>%
  mutate(row = row_number())%>%
  select(starts_with("Day"), zip, row, Shepsid) %>%
  tidyr::pivot_longer(c(-"zip", -"row", -"Shepsid")) %>%
  group_by(zip)%>%
  mutate(Date=as.Date(value))%>%
  filter(month(Date) >= 5 & month(Date) <= 9)%>%
  left_join(heat, by = c('zip', 'Date'))%>% 
  pivot_wider(names_from = name)%>%
  data.frame() 

Tr1_grp_hw <- birth3 %>%
  group_by(Shepsid)%>%
  summarize(Tr1_hw_count=sum(heatwave))

# Number of coldwaves experienced during trimester 3

birth3 <- birth_Tr3 %>%
  mutate(row = row_number())%>%
  dplyr::select(starts_with("Day"), zip, row, Shepsid) %>%
  tidyr::pivot_longer(c(-"zip", -"row", -"Shepsid")) %>%
  group_by(zip)%>%
  mutate(Date=as.Date(value))%>%
  filter(month(Date) >= 10 | month(Date) <= 4)%>%
  left_join(cold, by = c('zip', 'Date'))%>% 
  pivot_wider(names_from = name)%>%
  data.frame() 

Tr3_grp_cw <- birth3 %>%
  group_by(Shepsid)%>%
  summarize(Tr3_cw_count=sum(coldwave))

# Number of coldwaves experienced during trimester 2

birth3 <- birth_Tr2 %>%
  mutate(row = row_number())%>%
  select(starts_with("Day"), zip, row, Shepsid) %>%
  tidyr::pivot_longer(c(-"zip", -"row", -"Shepsid")) %>%
  group_by(zip)%>%
  mutate(Date=as.Date(value))%>%
  filter(month(Date) >= 10 | month(Date) <= 4)%>%
  left_join(cold, by = c('zip', 'Date'))%>% 
  pivot_wider(names_from = name)%>%
  data.frame() 

Tr2_grp_cw <- birth3 %>%
  group_by(Shepsid)%>%
  summarize(Tr2_cw_count=sum(coldwave))

# Number of coldwaves experienced during trimester 1

birth3 <- birth_Tr1 %>%
  mutate(row = row_number())%>%
  select(starts_with("Day"), zip, row, Shepsid) %>%
  tidyr::pivot_longer(c(-"zip", -"row", -"Shepsid")) %>%
  group_by(zip)%>%
  mutate(Date=as.Date(value))%>%
  filter(month(Date) >= 10 | month(Date) <= 4)%>%
  left_join(cold, by = c('zip', 'Date'))%>% 
  pivot_wider(names_from = name)%>%
  data.frame() 

Tr1_grp_cw <- birth3 %>%
  group_by(Shepsid)%>%
  summarize(Tr1_cw_count=sum(coldwave))

# Merge trimester data frames back to sample

df_wrong$Tr1 <- as.Date(df_wrong$Date) - (df_wrong$weeks_gest * 7)
df_wrong$Tr2 <- as.Date(df_wrong$Date) - ((df_wrong$weeks_gest * 7) - (13*7))
df_wrong$Tr3 <- as.Date(df_wrong$Date) - ((df_wrong$weeks_gest * 7) - (26*7))
df_wrong$Last_Gest_Week <- as.Date(df_wrong$Date) - 7
df_wrong$Tr3_end <- df_wrong$Date

births_hw <- df_wrong %>%
  filter(weeks_gest > 0)%>%
  left_join(Tr1_grp_hw, by=c('Shepsid'))%>%
  left_join(Tr2_grp_hw, by=c('Shepsid'))%>%
  left_join(Tr3_grp_hw, by=c('Shepsid'))%>%
  left_join(Tr1_grp_cw, by=c('Shepsid'))%>%
  left_join(Tr2_grp_cw, by=c('Shepsid'))%>%
  left_join(Tr3_grp_cw, by=c('Shepsid'))%>%
  dplyr::select(Shepsid, SMM21, ptzip, Date, Tr1, Tr2, Tr3, Last_Gest_Week, Year, Age, Race, Insurance, Ethnicity,
                weeks_gest, Tr1_hw_count, Tr2_hw_count, Tr3_hw_count, Tr1_cw_count, Tr2_cw_count, Tr3_cw_count)

write_parquet(births_hw, "Data/Trimester_Exposure.parquet")

















