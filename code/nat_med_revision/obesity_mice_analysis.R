######################################################################
## Code author: Steven Kerr steven.kerr@ed.ac.uk
## Description: Load in all data and create a cohort dataframe
######################################################################

library(tidyverse)
library(mice)
library(finalfit)
library(survival)
library(broom)

df = readRDS('/conf/EAVE/GPanalysis/progs/UA/second_booster_dose_failures/data/xdf_full_obesity_covid_hosp_death.RDS')  %>%
  select(-bmi_impute, -bmi_gp)

df <- filter(df, pv_period_f!="uv")
df <- filter(df, pv_period_f!="v1_0:3")
df <- filter(df, pv_period_f!="v1_4+")
df <- filter(df, pv_period_f!="v2_0:1")

df$pv_period_f <- as.factor(df$pv_period_f)
df$pv_period_f <- relevel(df$pv_period_f, ref = "v2_2:9")

df$simd2020_sc_quintile <- as.factor(df$simd2020_sc_quintile)
df$simd2020_sc_quintile <- relevel(df$simd2020_sc_quintile, ref = "5-Low")

df$vacc_gap <- as.factor(df$vacc_gap)
df$vacc_gap <- relevel(df$vacc_gap, ref = "<7wk")

df$ur_combined <- as.factor(df$ur_combined)
df$ur_combined <- relevel(df$ur_combined, ref = "2")

df$age_gp <- as.factor(df$age_gp)
df$age_gp <- relevel(df$age_gp, ref = "18-49")

# Get non-imputed BMI
q_bmi =  readRDS('/conf/EAVE/GPanalysis/data/cleaned_data/QCOVID_feb22.rds') %>%
  select(EAVE_LINKNO, Q_BMI) %>%
  mutate(Q_BMI = as.numeric(Q_BMI)) %>%
  # Set unrealistic values of BMI to NA
  mutate(Q_BMI = case_when(Q_BMI < 10 | Q_BMI > 100 ~ NA_real_,
                           TRUE ~ Q_BMI))

df = left_join(df, q_bmi)

# If they had an event, select that row; else select random
df_first = arrange(df, EAVE_LINKNO, desc(event)) %>%
  filter(!duplicated(EAVE_LINKNO))

# missing = missing_glimpse(df_first)

# Variables to be used in imputation
q_names = grep( 'Q', names(df), value = TRUE)
vars = c('Sex', 'ageYear', 'event', 'simd2020_sc_quintile',  setdiff(q_names, 'Q_ETHNICITY'))

# Carry out imputation
set.seed(123)

# For testing
#df_first = sample_n(df_first, 1000)

df_base = select(df, -Q_BMI)

imputation = mice(df_first[vars], m=10)

analysis = mice::complete(imputation, "all") 

# Replace elements of analysis with aggregated count data
for (i in 1:10){
  
  complete_i = analysis[[i]] %>%
    mutate(EAVE_LINKNO = df_first$EAVE_LINKNO,
           bmi_gp = cut(Q_BMI, breaks = c(0, 18.5, 25, 30, 40, Inf), right = FALSE,
                        labels = c('<18.5', '18.5-24.9', '25-29.9', '30-39.9', '40+' ))) %>% 
    select(EAVE_LINKNO, Q_BMI, bmi_gp) %>%
    right_join(df_base)
  
  complete_i$bmi_gp <- as.factor(complete_i$bmi_gp)
  complete_i$bmi_gp <- relevel(complete_i$bmi_gp, ref = "18.5-24.9")
  
  z.agg <- pyears(Surv(as.numeric(tstart), as.numeric(tstop), endpt) ~
                    pv_period_f + period_f + Sex + age_gp + n_risk_gps +
                    bmi_gp + simd2020_sc_quintile + pos_before_start_u + vacc_gap +
                    n_tests_gp + ur_combined,
                  data = complete_i, 
                  weight = ew, 
                  scale = 365.25, 
                  data.frame = TRUE)
  
  analysis[[i]] = z.agg$data
}

# Pool results
multi_imp_results = analysis %>%
  map(glm, formula = event ~ offset(log(pyears)) + pv_period_f + period_f + Sex + age_gp +
    n_risk_gps + bmi_gp + simd2020_sc_quintile + pos_before_start_u + vacc_gap +
    n_tests_gp + ur_combined,
    family=poisson) %>%
    pool()

# Create results table
results_table = tidy(multi_imp_results) %>%
  mutate(ucl = estimate + 1.96 * std.error,
         lcl = estimate - 1.96 * std.error) %>%
  select(term, estimate, lcl, ucl) %>%
  mutate_if(is.numeric, ~formatC(round(exp(.), 3), format = "f", big.mark = ",", drop0trailing = TRUE)) %>%
  mutate(estimate = paste0(estimate, ' (', lcl, ' - ', ucl, ')' )) %>%
  select(-lcl, -ucl)

write.csv(results_table, '/conf/EAVE/GPanalysis/progs/UA/second_booster_dose_failures/obesity_codes/nat_med_revision/obesity_multi_imp.csv')
