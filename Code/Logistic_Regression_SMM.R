
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
library(sjPlot)
setwd("~/RStudio/Matched SMM/")

# This data file includes SMM deliveries with and without gestation information. 
# To filter for only deliveries with gestation data, set weeks_gest > 0. 
# Around 6,200 observations have full gestation information
# Deliveries without gestation information do not have counts of heatwave or coldwaves per trimester, 
# So Tr1_hw_count:Tr3_cw_count will all be NA, as well as the dates for Tr1, Tr2, Tr3, and Last_Gest_Week

data <- read_parquet("Data/Trimester_Exposure_DELIVERIES.parquet")%>%
  rename(date=Date)%>%
  mutate(zip=as.numeric(ptzip))
data$SMM21[is.na(data$SMM21)] <- 0

regionvars <- read.csv("Data/Regional_Vars.csv")

data <- data %>%
  left_join(regionvars, by=c('zip'))

data <- data %>%
  mutate(Race_Group=ifelse(Race=="Black", "Black", "White or Other"), 
         Age_Group=ifelse((Age=="35-39" | Age=="40+"), "Age >= 35", "Age < 35"),
         Ethnicity_Group=ifelse(Ethnicity=="Hispanic", "Hispanic", "Not Hispanic/Unknown"),
         Insurance_Group=ifelse(Insurance=="Medicaid", "Medicaid", "Non-Medicaid"))

# Restrict to warm season for heatwave models
hw_data <- data%>%
  filter(month(date) >= 5 & month(date) <= 9)%>%
  filter(!is.na(LGW_hw_count))%>%
  mutate(heatwave=ifelse(LGW_hw_count > 0, 1, 0))

# Restrict to cold season for coldwave models
cw_data <- data%>%
  filter(month(date) <= 4 | month(date) >= 10)%>%
  filter(!is.na(LGW_cw_count))%>%
  mutate(coldwave=ifelse(LGW_cw_count > 0, 1, 0))

# Set reference categories
hw_data$Race_Group <- as.factor(hw_data$Race_Group)
hw_data$Insurance_Group <- as.factor(hw_data$Insurance_Group)
hw_data$Ethnicity_Group <- as.factor(hw_data$Ethnicity_Group)
hw_data$Age_Group <- as.factor(hw_data$Age_Group)
hw_data$RUCA_Cat <- as.factor(hw_data$RUCA_Cat)
hw_data$geo_region <- as.factor(hw_data$geo_region)
hw_data$ICE_Income_Tertile <- as.factor(hw_data$ICE_Income_Tertile)
hw_data$ICE_Race_Tertile <- as.factor(hw_data$ICE_Race_Tertile)

hw_data$Race_Group <- relevel(hw_data$Race_Group, ref = "White or Other")
hw_data$Insurance_Group <- relevel(hw_data$Insurance_Group, ref = "Non-Medicaid")
hw_data$Ethnicity_Group <- relevel(hw_data$Ethnicity_Group, ref = "Not Hispanic/Unknown")
hw_data$Age_Group <- relevel(hw_data$Age_Group, ref = "Age < 35")
hw_data$RUCA_Cat <- relevel(hw_data$RUCA_Cat, ref = "Rural")
hw_data$geo_region <- relevel(hw_data$geo_region, ref = "Piedmont")
hw_data$ICE_Income_Tertile <- relevel(hw_data$ICE_Income_Tertile, ref = "High")
hw_data$ICE_Race_Tertile <- relevel(hw_data$ICE_Race_Tertile, ref = "High")

cw_data$Race_Group <- as.factor(cw_data$Race_Group)
cw_data$Insurance_Group <- as.factor(cw_data$Insurance_Group)
cw_data$Ethnicity_Group <- as.factor(cw_data$Ethnicity_Group)
cw_data$Age_Group <- as.factor(cw_data$Age_Group)
cw_data$RUCA_Cat <- as.factor(cw_data$RUCA_Cat)
cw_data$geo_region <- as.factor(cw_data$geo_region)
cw_data$ICE_Income_Tertile <- as.factor(cw_data$ICE_Income_Tertile)
cw_data$ICE_Race_Tertile <- as.factor(cw_data$ICE_Race_Tertile)

cw_data$Race_Group <- relevel(cw_data$Race_Group, ref = "White or Other")
cw_data$Insurance_Group <- relevel(cw_data$Insurance_Group, ref = "Non-Medicaid")
cw_data$Ethnicity_Group <- relevel(cw_data$Ethnicity_Group, ref = "Not Hispanic/Unknown")
cw_data$Age_Group <- relevel(cw_data$Age_Group, ref = "Age < 35")
cw_data$RUCA_Cat <- relevel(cw_data$RUCA_Cat, ref = "Rural")
cw_data$geo_region <- relevel(cw_data$geo_region, ref = "Piedmont")
cw_data$ICE_Income_Tertile <- relevel(cw_data$ICE_Income_Tertile, ref = "High")
cw_data$ICE_Race_Tertile <- relevel(cw_data$ICE_Race_Tertile, ref = "High")

gc()

heatwave_indiv <- glmer(SMM21 ~ heatwave + Race_Group +  Ethnicity_Group + Age_Group + Insurance_Group +  (1 | zip), 
                      data = hw_data, 
                      family = poisson(link = "log"),
                      control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
gc()

coldwave_indiv <- glmer(SMM21 ~ coldwave + Race_Group +  Ethnicity_Group + Age_Group + Insurance_Group + (1 | zip), 
                        data = cw_data, 
                        family = poisson(link = "log"),
                        control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
gc()

heatwave_comm <- glmer(SMM21 ~ heatwave + RUCA_Cat + geo_region + ICE_Income_Tertile + ICE_Race_Tertile + (1 | zip), 
                        data = hw_data, 
                        family = poisson(link = "log"),
                        control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
gc()

coldwave_comm <- glmer(SMM21 ~ coldwave + RUCA_Cat + geo_region + ICE_Income_Tertile + ICE_Race_Tertile + (1 | zip), 
                       data = cw_data, 
                       family = poisson(link = "log"),
                       control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
gc()

heatwave_full <- glmer(SMM21 ~ heatwave + Race_Group +  Ethnicity_Group + Age_Group + Insurance_Group + RUCA_Cat + geo_region + ICE_Income_Tertile + ICE_Race_Tertile + (1 | zip), 
                       data = hw_data, 
                       family = poisson(link = "log"),
                       control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
gc()

coldwave_full <- glmer(SMM21 ~ coldwave + Race_Group +  Ethnicity_Group + Age_Group + Insurance_Group + RUCA_Cat + geo_region + ICE_Income_Tertile + ICE_Race_Tertile + (1 | zip), 
                       data = cw_data, 
                       family = poisson(link = "log"),
                       control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
gc()

tab_model(heatwave_full, heatwave_indiv, heatwave_comm)
tab_model(coldwave_full, coldwave_indiv, coldwave_comm)

write.csv(cw_data, "Data/coldwave_SMM_2011_2019.csv")
write.csv(hw_data, "Data/heatwave_SMM_2011_2019.csv")


