

library(grid)
library(forestploter)
library(data.table)
library(tidyr)

# Plot of daily relative risk for lag-2 to lag7

forest_data <- read.csv("Results/Daily Lags/combined_matrices_full.csv")%>%
  rename("lag-2"=lag.2,
         "lag-1"=lag.1)%>%
  dplyr::select(-Pred_Model)

forest_data2 <- pivot_longer(forest_data, cols = c("lag-2", "lag-1", "lag0", "lag1", "lag2", "lag3", "lag4", "lag5", "lag6", "lag7"), names_to = "Lag")
forest_data2 <- pivot_wider(forest_data2, names_from = c("Type"), values_from = value)
dt <- forest_data2

dt$`SMM` <- paste(rep(" ", 20), collapse = " ")

dt$'95% CI' <- paste(sprintf("%.2f (%.2f, %.2f)", dt$'matRRfit', dt$'matRRlow', dt$'matRRhigh'))

est = list(dt$'matRRfit')
lower =list(dt$'matRRlow')
upper = list(dt$'matRRhigh')

dt <- dt %>%
  dplyr::select(Lag, 'SMM', '95% CI')

p <- forest(dt,
            est = est,
            lower = lower, 
            upper = upper,
            ci_column = c(2),
            ref_line = 1,
            ticks_at = c(1),
            title = "Daily RR (May-Sept 2011-2019): Heatwaves"
)

plot(p)

# Plots of subgroups 

full_cumlag <- read.csv("Results/Heatwaves/Lag7/Cumulative_Lag7_Full.csv")%>%
  mutate(Category="All Deliveries")%>%
  mutate(Subgroup=" ")
other_cumlag <- read.csv("Results/Heatwaves/Lag7/Cumulative_Lag7_OtherRace.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="Other Race")
medicaid_cumlag <- read.csv("Results/Heatwaves/Lag7/Cumulative_Lag7_Medicaid.csv")%>%
  mutate(Category="Insurance")%>%
  mutate(Subgroup="Medicaid")
commercial_cumlag <- read.csv("Results/Heatwaves/Lag7/Cumulative_Lag7_PrivateIns.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="Private Insurance")
age_cumlag <- read.csv("Results/Heatwaves/Lag7/Cumulative_Lag7_Age.csv")%>%
  mutate(Category="Age")%>%
  mutate(Subgroup=">= 35")
black_cumlag <- read.csv("Results/Heatwaves/Lag7/Cumulative_Lag7_Black.csv")%>%
  mutate(Category="Race")%>%
  mutate(Subgroup="Black")
white_cumlag <- read.csv("Results/Heatwaves/Lag7/Cumulative_Lag7_White.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="White")
hispanic_cumlag <- read.csv("Results/Heatwaves/Lag7/Cumulative_Lag7_Hispanic.csv")%>%
  mutate(Category="Ethnicity")%>%
  mutate(Subgroup="Hispanic")
nonhisp_cumlag <- read.csv("Results/Heatwaves/Lag7/Cumulative_Lag7_NonHispanic.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="Not Hispanic")
selfpay_cumlag <- read.csv("Results/Heatwaves/Lag7/Cumulative_Lag7_SelfPay.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="Self Pay")
agelow_cumlag <- read.csv("Results/Heatwaves/Lag7/Cumulative_Lag7_AgeUnder35.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="< 35")
otherins_cumlag <- read.csv("Results/Heatwaves/Lag7/Cumulative_Lag7_OtherIns.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="Other Insurance")

urban <- read.csv("Results/Heatwaves/Lag7/Cumulative_Lag7_Urban.csv")%>%
  mutate(Category="RUCA")%>%
  mutate(Subgroup="Urban")
suburban <- read.csv("Results/Heatwaves/Lag7/Cumulative_Lag7_Suburban.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="Suburban")
rural <- read.csv("Results/Heatwaves/Lag7/Cumulative_Lag7_Rural.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="Rural")
western <- read.csv("Results/Heatwaves/Lag7/Cumulative_Lag7_Western.csv")%>%
  mutate(Category="Region")%>%
  mutate(Subgroup="Western")
coastal <- read.csv("Results/Heatwaves/Lag7/Cumulative_Lag7_Coastal.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="Coastal")
piedmont <- read.csv("Results/Heatwaves/Lag7/Cumulative_Lag7_Piedmont.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="Piedmont")
iceracelow <- read.csv("Results/Heatwaves/Lag7/Cumulative_Lag7_ICERace_Low.csv")%>%
  mutate(Category="ICE Race")%>%
  mutate(Subgroup="Mostly Nonwhite")
iceracemid <- read.csv("Results/Heatwaves/Lag7/Cumulative_Lag7_ICERace_Mid.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="Mixed-Race")
iceracehigh <- read.csv("Results/Heatwaves/Lag7/Cumulative_Lag7_ICERace_High.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="Mostly White")
iceincomelow <- read.csv("Results/Heatwaves/Lag7/Cumulative_Lag7_ICEIncome_Low.csv")%>%
  mutate(Category="ICE Income")%>%
  mutate(Subgroup="Low Income")
iceincomemid <- read.csv("Results/Heatwaves/Lag7/Cumulative_Lag7_ICEIncome_Mid.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="Mixed Income")
iceincomehigh <- read.csv("Results/Heatwaves/Lag7/Cumulative_Lag7_ICEIncome_High.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="High Income")

forest_data_lag7_hw <- rbind(full_cumlag, 
                     age_cumlag, 
                     agelow_cumlag, 
                     black_cumlag, 
                     white_cumlag,  
                     other_cumlag,  
                     medicaid_cumlag, 
                     commercial_cumlag, 
                     selfpay_cumlag, 
                     otherins_cumlag,
                     urban, 
                     suburban, 
                     rural, 
                     western,
                     coastal,
                     piedmont,
                     iceracelow, 
                     iceracemid, 
                     iceracehigh,
                     iceincomelow,
                     iceincomemid,
                     iceincomehigh
                     )

dt <- forest_data_lag7_hw %>%
  dplyr::select(-X, -estmean, -outcome)%>%
  distinct(over_rr, over_rr_se, over_rr_high, over_rr_low, Subgroup, Category)

# Add two blank columns for CI
dt$`SMM Lag7` <- paste(rep(" ", 20), collapse = " ")

# Generate point estimation and 95% CI. Paste two CIs together and separate by line break.
dt$`SMM Lag7 RR (95% CI)` <- ifelse(is.na(dt$over_rr_se), "",
                                sprintf("%.2f (%.2f to %.2f)",
                                        dt$over_rr, dt$over_rr_low, dt$over_rr_high))

est = list(dt$over_rr)
lower = list(dt$over_rr_low) 
upper = list(dt$over_rr_high)

dt <- dt %>%
  dplyr::select(Category, 
                Subgroup, 
                `SMM Lag7 RR (95% CI)`, 
                'SMM Lag7')

tm <- forest_theme(base_size = 10)
p <- forest(dt,
            est = est,
            lower = lower, 
            upper = upper,
            ci_column = c(4),
            ticks_at = c(1),
            ref_line = 1,
            title = "Cumulative Lag0 to Lag7, SMM",
            theme = tm)
plot(p)

# Plot subgroups for lag3

# Plots of subgroups 

full_cumlag <- read.csv("Results/Heatwaves/Lag3/Cumulative_Lag3_Full.csv")%>%
  mutate(Category="All Deliveries")%>%
  mutate(Subgroup=" ")
other_cumlag <- read.csv("Results/Heatwaves/Lag3/Cumulative_Lag3_OtherRace.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="Other Race")
medicaid_cumlag <- read.csv("Results/Heatwaves/Lag3/Cumulative_Lag3_Medicaid.csv")%>%
  mutate(Category="Insurance")%>%
  mutate(Subgroup="Medicaid")
commercial_cumlag <- read.csv("Results/Heatwaves/Lag3/Cumulative_Lag3_PrivateIns.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="Private Insurance")
age_cumlag <- read.csv("Results/Heatwaves/Lag3/Cumulative_Lag3_Age.csv")%>%
  mutate(Category="Age")%>%
  mutate(Subgroup=">= 35")
black_cumlag <- read.csv("Results/Heatwaves/Lag3/Cumulative_Lag3_Black.csv")%>%
  mutate(Category="Race")%>%
  mutate(Subgroup="Black")
white_cumlag <- read.csv("Results/Heatwaves/Lag3/Cumulative_Lag3_White.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="White")
hispanic_cumlag <- read.csv("Results/Heatwaves/Lag3/Cumulative_Lag3_Hispanic.csv")%>%
  mutate(Category="Ethnicity")%>%
  mutate(Subgroup="Hispanic")
nonhisp_cumlag <- read.csv("Results/Heatwaves/Lag3/Cumulative_Lag3_NonHispanic.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="Not Hispanic")
selfpay_cumlag <- read.csv("Results/Heatwaves/Lag3/Cumulative_Lag3_SelfPay.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="Self Pay")
agelow_cumlag <- read.csv("Results/Heatwaves/Lag3/Cumulative_Lag3_AgeUnder35.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="< 35")
otherins_cumlag <- read.csv("Results/Heatwaves/Lag3/Cumulative_Lag3_OtherIns.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="Other Insurance")

urban <- read.csv("Results/Heatwaves/Lag3/Cumulative_Lag3_Urban.csv")%>%
  mutate(Category="RUCA")%>%
  mutate(Subgroup="Urban")
suburban <- read.csv("Results/Heatwaves/Lag3/Cumulative_Lag3_Suburban.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="Suburban")
rural <- read.csv("Results/Heatwaves/Lag3/Cumulative_Lag3_Rural.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="Rural")
western <- read.csv("Results/Heatwaves/Lag3/Cumulative_Lag3_Western.csv")%>%
  mutate(Category="Region")%>%
  mutate(Subgroup="Western")
coastal <- read.csv("Results/Heatwaves/Lag3/Cumulative_Lag3_Coastal.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="Coastal")
piedmont <- read.csv("Results/Heatwaves/Lag3/Cumulative_Lag3_Piedmont.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="Piedmont")
iceracelow <- read.csv("Results/Heatwaves/Lag3/Cumulative_Lag3_ICERace_Low.csv")%>%
  mutate(Category="ICE Race")%>%
  mutate(Subgroup="Mostly Nonwhite")
iceracemid <- read.csv("Results/Heatwaves/Lag3/Cumulative_Lag3_ICERace_Mid.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="Mixed-Race")
iceracehigh <- read.csv("Results/Heatwaves/Lag3/Cumulative_Lag3_ICERace_High.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="Mostly White")
iceincomelow <- read.csv("Results/Heatwaves/Lag3/Cumulative_Lag3_ICEIncome_Low.csv")%>%
  mutate(Category="ICE Income")%>%
  mutate(Subgroup="Low Income")
iceincomemid <- read.csv("Results/Heatwaves/Lag3/Cumulative_Lag3_ICEIncome_Mid.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="Mixed Income")
iceincomehigh <- read.csv("Results/Heatwaves/Lag3/Cumulative_Lag3_ICEIncome_High.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="High Income")

forest_data_lag3_hw <- rbind(full_cumlag, 
                          age_cumlag, 
                          agelow_cumlag, 
                          black_cumlag, 
                          white_cumlag,  
                          other_cumlag,  
                          medicaid_cumlag, 
                          commercial_cumlag, 
                          selfpay_cumlag, 
                          otherins_cumlag,
                          urban, 
                          suburban, 
                          rural, 
                          western,
                          coastal,
                          piedmont,
                          iceracelow, 
                          iceracemid, 
                          iceracehigh,
                          iceincomelow,
                          iceincomemid,
                          iceincomehigh)

dt <- forest_data_lag3_hw %>%
  dplyr::select(-estmean, -outcome)%>%
  distinct(over_rr, over_rr_se, over_rr_high, over_rr_low, Subgroup, Category)

# Add two blank columns for CI
dt$`SMM Lag3` <- paste(rep(" ", 20), collapse = " ")

# Generate point estimation and 95% CI. Paste two CIs together and separate by line break.
dt$`SMM Lag3 RR (95% CI)` <- ifelse(is.na(dt$over_rr_se), "",
                               sprintf("%.2f (%.2f to %.2f)",
                                       dt$over_rr, dt$over_rr_low, dt$over_rr_high))

est = list(dt$over_rr)
lower = list(dt$over_rr_low) 
upper = list(dt$over_rr_high)

dt <- dt %>%
  dplyr::select(Category, 
                Subgroup, 
                `SMM Lag3 RR (95% CI)`, 
                'SMM Lag3')

tm <- forest_theme(base_size = 10)
p <- forest(dt,
            est = est,
            lower = lower, 
            upper = upper,
            ci_column = c(4),
            ticks_at = c(1),
            ref_line = 1,
            theme = tm)
plot(p)

### Lag3 and lag7 side by side

dt_lag7 <- forest_data_lag7 %>%
  dplyr::select(-X, -estmean, -outcome)%>%
  distinct(over_rr, over_rr_se, over_rr_high, over_rr_low, Subgroup, Category)

dt_lag3 <- forest_data_lag3 %>%
  dplyr::select(-X, -estmean, -outcome)%>%
  distinct(over_rr, over_rr_se, over_rr_high, over_rr_low, Subgroup, Category)

dt_lag3$`SMM Lag3` <- paste(rep(" ", 20), collapse = " ")

dt_lag3$`SMM Lag3 RR (95% CI)` <- ifelse(is.na(dt_lag3$over_rr_se), "",
                                     sprintf("%.2f (%.2f to %.2f)",
                                             dt_lag3$over_rr, dt_lag3$over_rr_low, dt_lag3$over_rr_high))

dt_lag7$`SMM Lag7` <- paste(rep(" ", 20), collapse = " ")

dt_lag7$`SMM Lag7 RR (95% CI)` <- ifelse(is.na(dt_lag7$over_rr_se), "",
                                          sprintf("%.2f (%.2f to %.2f)",
                                                  dt_lag7$over_rr, dt_lag7$over_rr_low, dt_lag7$over_rr_high))

est = list(dt_lag7$over_rr,
           dt_lag3$over_rr)
lower = list(dt_lag7$over_rr_low,
             dt_lag3$over_rr_low) 
upper = list(dt_lag7$over_rr_high,
             dt_lag3$over_rr_high)

dt <- dt_lag3%>%
  left_join(dt_lag7, by=c("Category", "Subgroup"))%>%
  dplyr::select(Category, 
                Subgroup, 
                `SMM Lag3 RR (95% CI)`, 
                'SMM Lag3', 
                `SMM Lag7 RR (95% CI)`, 
                'SMM Lag7')

tm <- forest_theme(base_size = 10)
p <- forest(dt,
            est = est,
            lower = lower, 
            upper = upper,
            ci_column = c(4, 6),
            ticks_at = c(1),
            ref_line = 1,
            title = "Heatwaves Cumulative RR",
            theme = tm)
plot(p)

# Lag3 and Lag 7 for cold waves and heat waves

heat_lag3 <- read.csv("Results/Lag3/Cumulative_Lag3_Full.csv")%>%
  mutate(Lag="Lag0 to Lag3")
heat_lag7 <- read.csv("Results/Lag7/Cumulative_Lag7_Full.csv")%>%
  mutate(Lag="Lag0 to Lag7")
cold_lag3 <- read.csv("Results/Cumulative_Lag3_Full_Cold.csv")%>%
  mutate(Lag="Lag0 to Lag3")
cold_lag7 <- read.csv("Results/Cumulative_Lag7_Full_Cold.csv")%>%
  mutate(Lag="Lag0 to Lag7")

heat_data <- rbind(heat_lag3, heat_lag7)
cold_data <- rbind(cold_lag3, cold_lag7)

heat_dt <- heat_data %>%
  dplyr::select(-X, -estmean, -outcome)%>%
  distinct(over_rr, over_rr_se, over_rr_high, over_rr_low, Lag)

cold_dt <- cold_data %>%
  dplyr::select(-X, -estmean, -outcome)%>%
  distinct(over_rr, over_rr_se, over_rr_high, over_rr_low, Lag)

# Add two blank columns for CI
heat_dt$`Heatwave` <- paste(rep(" ", 20), collapse = " ")

# Generate point estimation and 95% CI. Paste two CIs together and separate by line break.
heat_dt$`Heatwave RR (95% CI)` <- ifelse(is.na(heat_dt$over_rr_se), "",
                                     sprintf("%.2f (%.2f to %.2f)",
                                             heat_dt$over_rr, heat_dt$over_rr_low, heat_dt$over_rr_high))

# Add two blank columns for CI
cold_dt$`Coldwave` <- paste(rep(" ", 20), collapse = " ")

# Generate point estimation and 95% CI. Paste two CIs together and separate by line break.
cold_dt$`Coldwave RR (95% CI)` <- ifelse(is.na(cold_dt$over_rr_se), "",
                                         sprintf("%.2f (%.2f to %.2f)",
                                                 cold_dt$over_rr, cold_dt$over_rr_low, cold_dt$over_rr_high))

est = list(cold_dt$over_rr,
           heat_dt$over_rr)
lower = list(cold_dt$over_rr_low,
             heat_dt$over_rr_low) 
upper = list(cold_dt$over_rr_high,
             heat_dt$over_rr_high)

dt <- cold_dt%>%
  left_join(heat_dt, by=c("Lag"))%>%
  dplyr::select(Lag, 
                `Coldwave RR (95% CI)`, 
                'Coldwave', 
                `Heatwave RR (95% CI)`, 
                'Heatwave')

tm <- forest_theme(base_size = 10)
p <- forest(dt,
            est = est,
            lower = lower, 
            upper = upper,
            ci_column = c(3,5),
            ticks_at = c(1),
            ref_line = 1,
            theme = tm)
plot(p)

## Trimester model forest plots

dt <- read.csv("Data/trimester_models_fp.csv")%>%
  dplyr::filter(Trimester > 0)

# Add two blank columns for CI
dt$`Heatwave` <- paste(rep(" ", 20), collapse = " ")

# Generate point estimation and 95% CI. Paste two CIs together and separate by line break.
dt$`HW RR (95% CI)` <- paste(sprintf("%.2f (%.2f to %.2f)",
                                                 dt$over_rr_heat, dt$over_rr_low_heat, dt$over_rr_high_heat))

# Add two blank columns for CI
dt$`Coldwave` <- paste(rep(" ", 20), collapse = " ")

# Generate point estimation and 95% CI. Paste two CIs together and separate by line break.
dt$`CW RR (95% CI)` <- paste(sprintf("%.2f (%.2f to %.2f)",
                                                 dt$over_rr_cold, dt$over_rr_low_cold, dt$over_rr_high_cold))

est = list(dt$over_rr_heat,
           dt$over_rr_cold)
lower = list(dt$over_rr_low_heat,
             dt$over_rr_low_cold) 
upper = list(dt$over_rr_high_heat,
             dt$over_rr_high_cold)

dt <- dt%>%
  dplyr::select(Trimester, 
                `HW RR (95% CI)`, 
                'Heatwave',
                `CW RR (95% CI)`, 
                'Coldwave'
                )

tm <- forest_theme(base_size = 10)
p <- forest(dt,
            est = est,
            lower = lower, 
            upper = upper,
            ci_column = c(3,5),
            ticks_at = c(1),
            ref_line = 1,
            theme = tm,
            title = "SMM Risk by Trimester of Exposure")
plot(p)


#### Plot heatwave and coldwave subgroups ####

full_cumlag <- read.csv("Results/Coldwaves/Lag3/Cumulative_Lag3_Full_Cold.csv")%>%
  mutate(Category="All Deliveries")%>%
  mutate(Subgroup=" ")
other_cumlag <- read.csv("Results/Coldwaves/Lag3/Cumulative_Lag3_OtherRace_Cold.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="Other Race")
medicaid_cumlag <- read.csv("Results/Coldwaves/Lag3/Cumulative_Lag3_Medicaid_Cold.csv")%>%
  mutate(Category="Insurance")%>%
  mutate(Subgroup="Medicaid")
commercial_cumlag <- read.csv("Results/Coldwaves/Lag3/Cumulative_Lag3_PrivateIns_Cold.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="Private Insurance")
age_cumlag <- read.csv("Results/Coldwaves/Lag3/Cumulative_Lag3_Age_Cold.csv")%>%
  mutate(Category="Age")%>%
  mutate(Subgroup=">= 35")
black_cumlag <- read.csv("Results/Coldwaves/Lag3/Cumulative_Lag3_Black_Cold.csv")%>%
  mutate(Category="Race")%>%
  mutate(Subgroup="Black")
white_cumlag <- read.csv("Results/Coldwaves/Lag3/Cumulative_Lag3_White_Cold.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="White")
hispanic_cumlag <- read.csv("Results/Coldwaves/Lag3/Cumulative_Lag3_Hispanic_Cold.csv")%>%
  mutate(Category="Ethnicity")%>%
  mutate(Subgroup="Hispanic")
nonhisp_cumlag <- read.csv("Results/Coldwaves/Lag3/Cumulative_Lag3_NonHispanic_Cold.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="Not Hispanic")
selfpay_cumlag <- read.csv("Results/Coldwaves/Lag3/Cumulative_Lag3_SelfPay_Cold.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="Self Pay")
agelow_cumlag <- read.csv("Results/Coldwaves/Lag3/Cumulative_Lag3_AgeUnder35_Cold.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="< 35")
otherins_cumlag <- read.csv("Results/Coldwaves/Lag3/Cumulative_Lag3_OtherIns_Cold.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="Other Insurance")

urban <- read.csv("Results/Coldwaves/Lag3/Cumulative_Lag3_Urban_Cold.csv")%>%
  mutate(Category="RUCA")%>%
  mutate(Subgroup="Urban")
suburban <- read.csv("Results/Coldwaves/Lag3/Cumulative_Lag3_Suburban_Cold.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="Suburban")
rural <- read.csv("Results/Coldwaves/Lag3/Cumulative_Lag3_Rural_Cold.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="Rural")
western <- read.csv("Results/Coldwaves/Lag3/Cumulative_Lag3_Western.csv")%>%
  mutate(Category="Region")%>%
  mutate(Subgroup="Western")
coastal <- read.csv("Results/Coldwaves/Lag3/Cumulative_Lag3_Coastal.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="Coastal")
piedmont <- read.csv("Results/Coldwaves/Lag3/Cumulative_Lag3_Piedmont.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="Piedmont")
iceracelow <- read.csv("Results/Coldwaves/Lag3/Cumulative_Lag3_ICERace_Low.csv")%>%
  mutate(Category="ICE Race")%>%
  mutate(Subgroup="Mostly Nonwhite")
iceracemid <- read.csv("Results/Coldwaves/Lag3/Cumulative_Lag3_ICERace_Mid.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="Mixed-Race")
iceracehigh <- read.csv("Results/Coldwaves/Lag3/Cumulative_Lag3_ICERace_High.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="Mostly White")
iceincomelow <- read.csv("Results/Coldwaves/Lag3/Cumulative_Lag3_ICEIncome_Low.csv")%>%
  mutate(Category="ICE Income")%>%
  mutate(Subgroup="Low Income")
iceincomemid <- read.csv("Results/Coldwaves/Lag3/Cumulative_Lag3_ICEIncome_Mid.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="Mixed Income")
iceincomehigh <- read.csv("Results/Coldwaves/Lag3/Cumulative_Lag3_ICEIncome_High.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="High Income")

forest_data_lag3_cw <- rbind(full_cumlag, 
                             age_cumlag, 
                             agelow_cumlag, 
                             black_cumlag, 
                             white_cumlag,  
                             other_cumlag,  
                             medicaid_cumlag, 
                             commercial_cumlag, 
                             selfpay_cumlag, 
                             otherins_cumlag,
                             urban, 
                             suburban, 
                             rural, 
                             western,
                             coastal,
                             piedmont,
                             iceracelow, 
                             iceracemid, 
                             iceracehigh,
                             iceincomelow,
                             iceincomemid,
                             iceincomehigh)


full_cumlag <- read.csv("Results/Coldwaves/Lag7/Cumulative_Lag7_Full_Cold.csv")%>%
  mutate(Category="All Deliveries")%>%
  mutate(Subgroup=" ")
other_cumlag <- read.csv("Results/Coldwaves/Lag7/Cumulative_Lag7_OtherRace_Cold.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="Other Race")
medicaid_cumlag <- read.csv("Results/Coldwaves/Lag7/Cumulative_Lag7_Medicaid_Cold.csv")%>%
  mutate(Category="Insurance")%>%
  mutate(Subgroup="Medicaid")
commercial_cumlag <- read.csv("Results/Coldwaves/Lag7/Cumulative_Lag7_PrivateIns_Cold.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="Private Insurance")
age_cumlag <- read.csv("Results/Coldwaves/Lag7/Cumulative_Lag7_Age_Cold.csv")%>%
  mutate(Category="Age")%>%
  mutate(Subgroup=">= 35")
black_cumlag <- read.csv("Results/Coldwaves/Lag7/Cumulative_Lag7_Black_Cold.csv")%>%
  mutate(Category="Race")%>%
  mutate(Subgroup="Black")
white_cumlag <- read.csv("Results/Coldwaves/Lag7/Cumulative_Lag7_White_Cold.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="White")
hispanic_cumlag <- read.csv("Results/Coldwaves/Lag7/Cumulative_Lag7_Hispanic_Cold.csv")%>%
  mutate(Category="Ethnicity")%>%
  mutate(Subgroup="Hispanic")
nonhisp_cumlag <- read.csv("Results/Coldwaves/Lag7/Cumulative_Lag7_NonHispanic_Cold.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="Not Hispanic")
selfpay_cumlag <- read.csv("Results/Coldwaves/Lag7/Cumulative_Lag7_SelfPay_Cold.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="Self Pay")
agelow_cumlag <- read.csv("Results/Coldwaves/Lag7/Cumulative_Lag7_AgeUnder35_Cold.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="< 35")
otherins_cumlag <- read.csv("Results/Coldwaves/Lag7/Cumulative_Lag7_OtherIns_Cold.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="Other Insurance")

urban <- read.csv("Results/Coldwaves/Lag7/Cumulative_Lag7_Urban_Cold.csv")%>%
  mutate(Category="RUCA")%>%
  mutate(Subgroup="Urban")
suburban <- read.csv("Results/Coldwaves/Lag7/Cumulative_Lag7_Suburban_Cold.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="Suburban")
rural <- read.csv("Results/Coldwaves/Lag7/Cumulative_Lag7_Rural_Cold.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="Rural")
western <- read.csv("Results/Coldwaves/Lag7/Cumulative_Lag7_Western.csv")%>%
  mutate(Category="Region")%>%
  mutate(Subgroup="Western")
coastal <- read.csv("Results/Coldwaves/Lag7/Cumulative_Lag7_Coastal.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="Coastal")
piedmont <- read.csv("Results/Coldwaves/Lag7/Cumulative_Lag7_Piedmont.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="Piedmont")
iceracelow <- read.csv("Results/Coldwaves/Lag7/Cumulative_Lag7_ICERace_Low.csv")%>%
  mutate(Category="ICE Race")%>%
  mutate(Subgroup="Mostly Nonwhite")
iceracemid <- read.csv("Results/Coldwaves/Lag7/Cumulative_Lag7_ICERace_Mid.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="Mixed-Race")
iceracehigh <- read.csv("Results/Coldwaves/Lag7/Cumulative_Lag7_ICERace_High.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="Mostly White")
iceincomelow <- read.csv("Results/Coldwaves/Lag7/Cumulative_Lag7_ICEIncome_Low.csv")%>%
  mutate(Category="ICE Income")%>%
  mutate(Subgroup="Low Income")
iceincomemid <- read.csv("Results/Coldwaves/Lag7/Cumulative_Lag7_ICEIncome_Mid.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="Mixed Income")
iceincomehigh <- read.csv("Results/Coldwaves/Lag7/Cumulative_Lag7_ICEIncome_High.csv")%>%
  mutate(Category=" ")%>%
  mutate(Subgroup="High Income")

forest_data_lag7_cw <- rbind(full_cumlag, 
                             age_cumlag, 
                             agelow_cumlag, 
                             black_cumlag, 
                             white_cumlag,  
                             other_cumlag,  
                             medicaid_cumlag, 
                             commercial_cumlag, 
                             selfpay_cumlag, 
                             otherins_cumlag,
                             urban, 
                             suburban, 
                             rural, 
                             western,
                             coastal,
                             piedmont,
                             iceracelow, 
                             iceracemid, 
                             iceracehigh,
                             iceincomelow,
                             iceincomemid,
                             iceincomehigh)

dt_lag7_hw <- forest_data_lag7_hw %>%
  dplyr::select(-X, -estmean, -outcome)%>%
  distinct(over_rr, over_rr_se, over_rr_high, over_rr_low, Subgroup, Category)

dt_lag3_hw <- forest_data_lag3_hw %>%
  dplyr::select(-X, -estmean, -outcome)%>%
  distinct(over_rr, over_rr_se, over_rr_high, over_rr_low, Subgroup, Category)

dt_lag7_cw <- forest_data_lag7_cw %>%
  dplyr::select(-X, -estmean, -outcome)%>%
  distinct(over_rr, over_rr_se, over_rr_high, over_rr_low, Subgroup, Category)

dt_lag3_cw <- forest_data_lag3_cw %>%
  dplyr::select(-X, -estmean, -outcome)%>%
  distinct(over_rr, over_rr_se, over_rr_high, over_rr_low, Subgroup, Category)


dt_lag3_hw$`Heatwave Lag0-Lag3` <- paste(rep(" ", 20), collapse = " ")

dt_lag3_hw$`HW Lag0-Lag3 RR (95% CI)` <- ifelse(is.na(dt_lag3_hw$over_rr_se), "",
                                         sprintf("%.2f (%.2f to %.2f)",
                                                 dt_lag3_hw$over_rr, dt_lag3_hw$over_rr_low, dt_lag3_hw$over_rr_high))

dt_lag7_hw$`Heatwave Lag0-Lag7` <- paste(rep(" ", 20), collapse = " ")

dt_lag7_hw$`HW Lag0-Lag7 RR (95% CI)` <- ifelse(is.na(dt_lag7_hw$over_rr_se), "",
                                         sprintf("%.2f (%.2f to %.2f)",
                                                 dt_lag7_hw$over_rr, dt_lag7_hw$over_rr_low, dt_lag7_hw$over_rr_high))

dt_lag3_cw$`Coldwave Lag0-Lag3` <- paste(rep(" ", 20), collapse = " ")

dt_lag3_cw$`CW Lag0-Lag3 RR (95% CI)` <- ifelse(is.na(dt_lag3_cw$over_rr_se), "",
                                           sprintf("%.2f (%.2f to %.2f)",
                                                   dt_lag3_cw$over_rr, dt_lag3_cw$over_rr_low, dt_lag3_cw$over_rr_high))

dt_lag7_cw$`Coldwave Lag0-Lag7` <- paste(rep(" ", 20), collapse = " ")

dt_lag7_cw$`CW Lag0-Lag7 RR (95% CI)` <- ifelse(is.na(dt_lag7_cw$over_rr_se), "",
                                           sprintf("%.2f (%.2f to %.2f)",
                                                   dt_lag7_cw$over_rr, dt_lag7_cw$over_rr_low, dt_lag7_cw$over_rr_high))


est = list(dt_lag3_hw$over_rr,
           dt_lag7_hw$over_rr,
           dt_lag3_cw$over_rr,
           dt_lag7_cw$over_rr)
lower = list(dt_lag3_hw$over_rr_low,
             dt_lag7_hw$over_rr_low,
             dt_lag3_cw$over_rr_low,
             dt_lag7_cw$over_rr_low) 
upper = list(dt_lag3_hw$over_rr_high,
             dt_lag7_hw$over_rr_high,
             dt_lag3_cw$over_rr_high,
             dt_lag7_cw$over_rr_high)

dt <- dt_lag3_hw%>%
  left_join(dt_lag7_hw, by=c("Category", "Subgroup"))%>%
  left_join(dt_lag3_cw, by=c("Category", "Subgroup"))%>%
  left_join(dt_lag7_cw, by=c("Category", "Subgroup"))%>%
  dplyr::select(Category, 
                Subgroup, 
                `HW Lag0-Lag3 RR (95% CI)`, 
                'Heatwave Lag0-Lag3', 
                `HW Lag0-Lag7 RR (95% CI)`, 
                'Heatwave Lag0-Lag7',
                `CW Lag0-Lag3 RR (95% CI)`, 
                'Coldwave Lag0-Lag3', 
                `CW Lag0-Lag7 RR (95% CI)`, 
                'Coldwave Lag0-Lag7')

tm <- forest_theme(base_size = 10)
p <- forest(dt,
            est = est,
            lower = lower, 
            upper = upper,
            ci_column = c(4, 6, 8, 10),
            ticks_at = c(1),
            ref_line = 1,
            title = "Cumulative SMM Relative Risk (RR)",
            theme = tm)
plot(p)



