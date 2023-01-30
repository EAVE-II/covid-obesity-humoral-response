##########################################################
# Name of file: 01a_cohort_prep.R
# Data release (if applicable):
# Original author(s): Utkarsh Agrawal
# Original date: 7 Sep 2021
# Type of script: Descriptive stats
# Version of R that the script was most recently run on: R 3.5.1
# Description of content: reads in the cohort and merges in vaccination data 
# Approximate run time: Unknown
##########################################################

# 01 Setup ####
#Libraries
library(plyr)
library(tidyverse)
library(survival)
library(lubridate)
#Load data

Location <- "/conf/"  # Server
#Location <- "//isdsf00d03/"  # Desktop

project_path <- paste0(Location,"EAVE/GPanalysis/progs/UA/second_booster_dose_failures")
project_path_vaccine <- paste0(Location,"EAVE/GPanalysis/progs/CR/Vaccine")
#a_begin <- as.Date("2021-04-05")  #second dose
a_begin <- as.Date("2021-09-14")  #beginning
a_end <- as.Date("2022-03-19")  # event end date
#get covid hospitalisations in Study period  - 14 days of a positive test
#first positive test in study period

smr01 <- readRDS(paste0(Location,"EAVE/GPanalysis/data/SMR01_allstays.rds"))
z <- smr01 %>%  
 filter(covid_main_diag_admit==1)
z <- select(z, EAVE_LINKNO, ADMISSION_DATE)%>%
  filter(ADMISSION_DATE >= a_begin & ADMISSION_DATE <= a_end )
all_hosp_covid_cert <- z

all_deaths  <- readRDS(paste0(Location,"EAVE/GPanalysis/data/all_deaths.rds"))
summary(all_deaths)
#get covid death certificate deaths
z <- all_deaths %>%  
  mutate(across(UNDERLYING_CAUSE_OF_DEATH:CAUSE_OF_DEATH_CODE_9, ~if_else(. %in% c("U071","U072"), 1,0)))
z <- z %>% rowwise() %>% mutate(rowsum = sum(c_across(UNDERLYING_CAUSE_OF_DEATH:CAUSE_OF_DEATH_CODE_9))) %>% 
  mutate(covid_death_cert = if_else(rowsum>=1,1,0)) %>% 
  dplyr::select(EAVE_LINKNO, NRS.Date.Death, covid_death_cert)
all_deaths_covid_dth_cert <- z
z <- all_deaths_covid_dth_cert %>% filter(covid_death_cert==1) %>% 
  filter(NRS.Date.Death >= a_begin & NRS.Date.Death <= a_end )
all_deaths_covid_dth_cert <- z
all_deaths_covid_dth_cert$covid_death_cert<-NULL
colnames(all_deaths_covid_dth_cert)[2] <- "ADMISSION_DATE"

# make sure all the events are in the time period

#combined hospital death endpoint
z1 <- all_hosp_covid_cert %>% mutate(hosp=1)
z2 <- all_deaths_covid_dth_cert %>% mutate(hosp=2)
z <- bind_rows(z1,z2) %>% 
  arrange(EAVE_LINKNO, hosp) %>%  #take hositalisation first
  filter(!duplicated(EAVE_LINKNO))
covid_hosp_death <- z %>% dplyr::select(-hosp)
covid_hosp_death$ADMISSION_DATE <- NULL
covid_hosp_death <- mutate(covid_hosp_death, certi_pre = 1)


xdf_full_covid_hosp_death <- readRDS(paste0(project_path,"/data/xdf_full_obesity_covid_hosp_death.RDS"))
z <- xdf_full_covid_hosp_death %>%
  select(EAVE_LINKNO, event)
z <- left_join(z, covid_hosp_death, by = "EAVE_LINKNO")
colnames(z)[3] <- "event_n"
z$event<-NULL
z$EAVE_LINKNO<-NULL
xdf_full_covid_hosp_death <- cbind(xdf_full_covid_hosp_death, z)

z <- xdf_full_covid_hosp_death %>%
  select(EAVE_LINKNO, endpt) %>%
  filter(endpt==1) %>%
  filter(!duplicated(EAVE_LINKNO))
z <- left_join(z, covid_hosp_death, by = "EAVE_LINKNO")
z <- filter(z, is.na(certi_pre))
z$certi_pre[is.na(z$certi_pre)] <- 0
z$endpt<-NULL
xdf_full_covid_hosp_death <- left_join(xdf_full_covid_hosp_death, z,
                                       by = "EAVE_LINKNO")

xdf_full_covid_hosp_death$event_n[is.na(xdf_full_covid_hosp_death$event_n)] <- 0
xdf_full_covid_hosp_death$certi_pre <- coalesce(xdf_full_covid_hosp_death$certi_pre,
                                                xdf_full_covid_hosp_death$endpt)
colnames(xdf_full_covid_hosp_death)[82] <- "endpt_n"


saveRDS(xdf_full_covid_hosp_death, paste0(project_path,"/data/xdf_full_obesity_covid_hosp_death_certi.RDS"))
## now colace the two speciment dates for covid-19

