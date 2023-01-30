##########################################################
# Name of file: 01a_Vaccinations_Input.R
# Original author(s): Utkarsh Agrawal
# Description of content: Analysis of vaccination failure cohort 
# Approximate run time: Unknown
##########################################################

# 01 Setup ####
#Libraries
library(plyr)
library(tidyverse)
library(survival)
#library(dplyr)
#library(mgcv)
library(tidyr)
library(ggplot2)

Location <- "/conf/"
project_path <- paste0(Location,"EAVE/GPanalysis/progs/UA/second_booster_dose_failures")
project_path_vaccine <- paste0(Location,"EAVE/GPanalysis/progs/CR/Vaccine")
df_cohort_vacc1 <- readRDS(paste0(project_path,"/data/xdf_full_obesity.RDS"))
xdf_full_covid_hosp_death <- readRDS(paste0(project_path,"/data/xdf_full_obesity_covid_hosp_death.RDS"))
##################
# data prep for the analysis

rg <- readRDS(paste0(project_path_vaccine,"/output/temp/Qcovid_all.rds"))
rg <- rg %>% dplyr::select(-(Sex:ur6_2016_name))
rg <- filter(rg, !duplicated(EAVE_LINKNO))
rg <- select(rg, EAVE_LINKNO, Q_BMI)
df_cohort_vacc1 <- left_join(df_cohort_vacc1,rg, by="EAVE_LINKNO")

df_cohort_vacc1 <- df_cohort_vacc1 %>% 
  mutate(bmi_gp_r = cut(Q_BMI, breaks = c(-1, 18.5, 24.9, 29.9,39.9, 51),
                      labels=c("<18.5", "18.5-24.9","25-29.9","30-39.9","40+")))
###############