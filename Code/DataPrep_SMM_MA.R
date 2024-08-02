
## This code pulls selected maternal mental health outcomes from ED delivery data ("deliveries.sas7bdat") using ICD-10 codes. 

## Call packages and set working directory 

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
setwd("~/RStudio/Matched SMM/")

## Read in data files

# SMM delivery data for North Carolina

df_wrong <- read_sas("Data/smm.sas7bdat")%>%
  mutate(Date=as.Date(admitdt))%>%
  filter(ptzip >= 27006 & ptzip <= 28909)%>%
  filter(year(Date) >= 2011 & year(Date) <= 2019)%>%
  mutate(Year=year(Date))
table(df_wrong$Year)
table(df_wrong$fyear)

#### AGE ####
df_wrong$Age <- as.factor(ifelse(df_wrong$agey < 20, '18-19', 
                                 ifelse(df_wrong$agey<25, '20-24',
                                        ifelse(df_wrong$agey <30, '25-29',
                                               ifelse(df_wrong$agey<35,'30-34',
                                                      ifelse(df_wrong$agey<40, '35-39',
                                                             ifelse(df_wrong$agey <= 44, '40+', 0)))))))
df_wrong$Age <- factor(df_wrong$Age)
table(df_wrong$Age)

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

df_wrong$Insurance <- as.factor(ifelse(df_wrong$payer1 == "09", "Self-pay",
                                     ifelse(df_wrong$payer1 == "MC", "Medicaid",
                                            ifelse(df_wrong$payer1 %in% c("MA", "MB", "OF", "VA", "TV", "11", "16", "CH", "MM"), "Other gov't", 
                                                   ifelse(df_wrong$payer1 %in% c("LM", "LI", "HM", "DS", "CI", "BL", "AM", "12", "13", "14", "15"),
                                                          "Commercial", "Other/Unknown")))))

df_wrong$Insurance <- factor(df_wrong$Insurance)
table(df_wrong$Insurance)

table(df_wrong$SMM20)
table(df_wrong$SMM21)

tabdat <- df_wrong %>%
  dplyr::select(Shepsid,
                Date,
                Year,
                ptzip,
                Age, 
                Race, 
                Ethnicity, 
                Insurance,
                SMM20,
                SMM21
  )

#### Combine with regional data ####

regionvars <- read.csv("Data/Regional_Vars.csv")

cw <- read.csv("Data/Crosswalk_NC.csv")

tabdat <- tabdat %>%
  mutate(ZIP=as.numeric(ptzip))

tabdat <- tabdat %>%
  left_join(cw, by=c('ZIP'))

tabdat <- tabdat %>%
  rename(zip=ZIP)

tabdat <- tabdat %>%
  left_join(regionvars, by=c('zip'))

table(tabdat$geo_region)

#### Create table one ####

library(dplyr)
library(tableone)

# Define categorical variables
catVars <- c("Year", "Age", "Race", "Ethnicity", "Insurance", "ICE_Income_Tertile", "ICE_Race_Tertile", "RUCA_Cat", "geo_region")

# Create table one 
smm_tab <- CreateTableOne(vars = catVars, data = tabdat, factorVars = catVars)
smm_tabmat <- print(smm_tab, showAllLevels = TRUE, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE, smd = FALSE)
write.csv(smm_tabmat, file = "Results/TableOne_SMM_All.csv")

# Filter for warm season
tabwarm <- tabdat %>%
  filter(month(Date) <= 9 & month(Date) >= 5)
smm_tab <- CreateTableOne(vars = catVars, data = tabwarm, factorVars = catVars)
smm_tabmat <- print(smm_tab, showAllLevels = TRUE, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE, smd = FALSE)
write.csv(smm_tabmat, file = "Results/TableOne_SMM_Warm_Season.csv")

tabcold <- tabdat %>%
  filter(month(Date) <= 4 | month(Date) >= 10)
smm_tab <- CreateTableOne(vars = catVars, data = tabcold, factorVars = catVars)
smm_tabmat <- print(smm_tab, showAllLevels = TRUE, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE, smd = FALSE)
write.csv(smm_tabmat, file = "Results/TableOne_SMM_Cold_Season.csv")

#### Create SMM file for use in matching ####

deliv <- tabdat%>% # Read in deliveries file
  mutate(zip=as.numeric(ptzip),
         date=as.Date(Date))%>% # This is not in ZCTA form so we will have to crosswalk it below
  dplyr::select(Shepsid,
                zip, 
                date, 
                SMM21, 
                Race, 
                Age, 
                Insurance, 
                Ethnicity) # Filter for variables of interest 

# Read in crosswalk file to convert from zip code to ZCTA
cw <- read.csv("Data/Crosswalk_NC.csv")%>%
  setDT()%>%
  rename(zip=ZIP)

# Crosswalk the delivery data 
deliv <- deliv %>%
  left_join(cw, by=c('zip'))%>%
  dplyr::select(-zip, ZIP_TYPE)%>%
  rename(zip=ZCTA) #The 'zip' variable in the deliv dataframe will now represent ZCTAs, not zip codes. There should be ~804 unique ZCTAs. 

write_parquet(deliv,"Data/Outcomes_SMM_v2.parquet")

#### Create shapefile of SMM counts for mapping ####

data <- read_parquet("Data/Outcomes_SMM_v2.parquet")%>%
  group_by(zip)%>%
  summarize(SMM=sum(SMM21))

shp <- read_sf("Data/shp/tl_2020_us_zcta510.shp")%>%
  rename(zip=ZCTA5CE10)%>%
  filter(zip >= 27006 & zip <= 28909)%>%
  mutate(zip=as.numeric(zip))

shp <- shp %>%
  left_join(data, by=c('zip'))

write_sf(shp, "Data/shp/SMM.shp")

#### HEATWAVE PREP ####

data <- read_parquet("Data/Heatwave_Metrics.parquet")%>%
  setDT()%>%
  filter(Zip >= 27006 & Zip <= 28909)%>%  # Filter for NC zip codes
  rename(date=Date, # Rename date column to lowercase
         heatwave=Heatwave, # Select Heatwave column as our heatwave variable
         zip=Zip)%>% # Rename zip code column to lowercase
  mutate(year=year(date),
         doy=yday(date),
         dow=wday(date))%>%
  dplyr::select(date, heatwave, year, doy, dow, zip)%>% # Select variables of interest
  mutate(heatwave=as.integer(heatwave)) %>% # Set as an integer just in case
  filter(year(date) >= 2011 & year(date) <= 2019)%>% # Filter for years of interest
  filter(month(date) >= 5 & month(date) <= 9) # Filter for months of interest

#Filter out heatwaves with lag periods outside of May-September (heatwaves occurring during the last week of September)
data <- data %>% # Create sept_heatwave variable that =1 if heatwave occurs during September 23-30
  mutate(sept_heatwave=ifelse(heatwave==1 & date %in% c("2016-09-23","2016-09-24","2016-09-25", "2016-09-26", "2016-09-27", "2016-09-28", "2016-09-29", "2016-09-30", 
                                                        "2017-09-23","2017-09-24","2017-09-25", "2017-09-26", "2017-09-27", "2017-09-28", "2017-09-29", "2017-09-30", 
                                                        "2018-09-23","2018-09-24","2018-09-25", "2018-09-26", "2018-09-27", "2018-09-28", "2018-09-29", "2018-09-30", 
                                                        "2019-09-23","2019-09-24","2019-09-25", "2019-09-26", "2019-09-27", "2019-09-28", "2019-09-29", "2019-09-30"), 1, 0))

data <- data %>%
  mutate(heatwave=ifelse(heatwave==1 & sept_heatwave==0, 1, 0)) %>% # Remove late September heatwaves from the dataset
  dplyr::select(-sept_heatwave)

write_parquet(data, "Data/Heatwave_Metrics_v2.parquet")

heat <- read_parquet("Data/Heatwave_Metrics_v2.parquet")%>%
  group_by(zip)%>%
  summarize(cwcount=sum(coldwave))

shp <- shp %>%
  left_join(heat, by=c('zip'))

write_sf(shp, "Data/shp/Heatwave_Count.shp")

#### COLDWAVE PREP ####

data <- read_parquet("Data/NC_Coldwave.parquet")%>%
  setDT()%>%
  filter(Zip >= 27006 & Zip <= 28909)%>%  # Filter for NC zip codes
  rename(date=Date, # Rename date column to lowercase
         zip=Zip)%>% # Rename zip code column to lowercase
  mutate(year=year(date),
         doy=yday(date),
         dow=wday(date))%>%
  mutate(coldwave=ifelse(ECF > 1, 1, 0))%>%
  dplyr::select(date, coldwave, year, doy, dow, zip)%>% # Select variables of interest
  mutate(coldwave=as.integer(coldwave)) %>% # Set as an integer just in case
  filter(year(date) >= 2011 & year(date) <= 2019)%>% # Filter for years of interest
  filter(month(date) >= 10 | month(date) <= 4) # Filter for months of interest
table(data$coldwave)

#Filter out coldwaves with lag periods outside of October-April (coldwaves occurring during the last week of April)
data <- data %>% # Create sept_heatwave variable that =1 if heatwave occurs during September 23-30
  mutate(april_coldwave=ifelse(coldwave==1 & date %in% c("2016-04-23","2016-04-24","2016-04-25", "2016-04-26", "2016-04-27", "2016-04-28", "2016-04-29", "2016-04-30", 
                                                         "2017-04-23","2017-04-24","2017-04-25", "2017-04-26", "2017-04-27", "2017-04-28", "2017-04-29", "2017-04-30", 
                                                         "2018-04-23","2018-04-24","2018-04-25", "2018-04-26", "2018-04-27", "2018-04-28", "2018-04-29", "2018-04-30", 
                                                         "2019-04-23","2019-04-24","2019-04-25", "2019-04-26", "2019-04-27", "2019-04-28", "2019-04-29", "2019-04-30"), 1, 0))

data <- data %>%
  mutate(coldwave=ifelse(coldwave==1 & april_coldwave==0, 1, 0)) %>% # Remove late September coldwaves from the dataset
  dplyr::select(-april_coldwave)

write_parquet(data, "Data/Coldwave_Metrics_v2.parquet")

# Create coldwave file for mapping 

cold <- read_parquet("Data/Coldwave_Metrics_v2.parquet")%>%
  group_by(zip)%>%
  summarize(cwcount=sum(coldwave))

shp <- shp %>%
  left_join(cold, by=c('zip'))

write_sf(shp, "Data/shp/Coldwave_Count.shp")











