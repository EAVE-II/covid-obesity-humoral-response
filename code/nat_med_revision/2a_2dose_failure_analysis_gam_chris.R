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
library(survminer)
#library(dplyr)
#library(mgcv)
library(tidyr)
library(ggplot2)
library(scales)

Location <- "/conf/"
project_path <- paste0(Location,"EAVE/GPanalysis/progs/UA/second_booster_dose_failures")
xdf_full_covid_hosp_death <- readRDS(paste0(project_path,"/data/xdf_full_obesity_covid_hosp_death_certi.RDS"))
xdf_full_covid_hosp_death <- readRDS(paste0(project_path,"/data/xdf_full_obesity_covid_hosp_death.RDS"))
#df_cohort_vacc <- readRDS(paste0(project_path,"/data/df_cohort_vacc_15-03-2022_gam.rds"))
#df_cohort_vacc_g <- readRDS(paste0(project_path,"/data/df_cohort_vacc_g_15-03-2022_gam.rds"))
##################
# data prep for the analysis
# method 1

df_cohort_vacc_g <- filter(xdf_full_covid_hosp_death, pv_period_f!="uv")
df_cohort_vacc_g <- filter(df_cohort_vacc_g, pv_period_f!="v1_0:3")
df_cohort_vacc_g <- filter(df_cohort_vacc_g, pv_period_f!="v1_4+")
df_cohort_vacc_g <- filter(df_cohort_vacc_g, pv_period_f!="v2_0:1")

## run this for certi dataset
df_cohort_vacc_g$event<-NULL
colnames(df_cohort_vacc_g)[80]<-"event"
##

df_cohort_vacc_g$pv_period_f <- as.factor(df_cohort_vacc_g$pv_period_f)
df_cohort_vacc_g$pv_period_f <- relevel(df_cohort_vacc_g$pv_period_f, ref = "v2_2:9")

df_cohort_vacc_g$simd2020_sc_quintile <- as.factor(df_cohort_vacc_g$simd2020_sc_quintile)
df_cohort_vacc_g$simd2020_sc_quintile <- relevel(df_cohort_vacc_g$simd2020_sc_quintile, ref = "5-Low")

df_cohort_vacc_g$bmi_gp <- as.factor(df_cohort_vacc_g$bmi_gp)
df_cohort_vacc_g$bmi_gp <- relevel(df_cohort_vacc_g$bmi_gp, ref = "18.5-24.9")

df_cohort_vacc_g$vacc_gap <- as.factor(df_cohort_vacc_g$vacc_gap)
df_cohort_vacc_g$vacc_gap <- relevel(df_cohort_vacc_g$vacc_gap, ref = "<7wk")

df_cohort_vacc_g$ur_combined <- as.factor(df_cohort_vacc_g$ur_combined)
df_cohort_vacc_g$ur_combined <- relevel(df_cohort_vacc_g$ur_combined, ref = "2")

df_cohort_vacc_g$age_gp <- as.factor(df_cohort_vacc_g$age_gp)
df_cohort_vacc_g$age_gp <- relevel(df_cohort_vacc_g$age_gp, ref = "18-49")

# Method 1
z.agg <- pyears(Surv(as.numeric(tstart),as.numeric(tstop),endpt_n) ~
                  pv_period_f+period_f + Sex + age_gp+n_risk_gps+
                  bmi_gp+simd2020_sc_quintile+pos_before_start_u+vacc_gap+
                  n_tests_gp+ur_combined,
                data=df_cohort_vacc_g, weight=ew , scale=365.25, data.frame=TRUE)
glm_pos1 <- glm(event ~ offset(log(pyears)) +pv_period_f+period_f + Sex + age_gp+
                  n_risk_gps+bmi_gp+simd2020_sc_quintile+pos_before_start_u+vacc_gap+
                  n_tests_gp+ur_combined,
                family=poisson, data=z.agg$data)

z<-as.data.frame(round(exp(glm_pos1$coefficients),3))
z<-cbind(z,round(exp(glm_pos1$coefficients-coef(summary(glm_pos1))[,2]*1.96),3))
z<-cbind(z,round(exp(glm_pos1$coefficients+coef(summary(glm_pos1))[,2]*1.96),3))
z<-cbind(z, paste(z[,1]," (",z[,2],"-",z[,3],")", sep = ""))
colnames(z)[1]<-"Adjusted Rate Ratio"
colnames(z)[2]<-"LCI"
colnames(z)[3]<-"UCI"
colnames(z)[4]<-"ARR (LCI-UCI)"
View(z)
z1<-z

write.csv(z,paste0(project_path,"/obesity_codes/output/ARR_obesity_certi.csv"))

########
# interaction analysis

glm_pos2 <- glm(event ~ offset(log(pyears)) +pv_period_f+period_f + Sex + age_gp+
                 n_risk_gps+bmi_gp+simd2020_sc_quintile+pos_before_start_u+vacc_gap+
                 n_tests_gp+ur_combined+bmi_gp*pv_period_f, 
               family=poisson, data=z.agg$data)

z<-as.data.frame(round(exp(glm_pos2$coefficients),3))
z<-cbind(z,round(exp(glm_pos2$coefficients-coef(summary(glm_pos2))[,2]*1.96),3))
z<-cbind(z,round(exp(glm_pos2$coefficients+coef(summary(glm_pos2))[,2]*1.96),3))
z<-cbind(z, paste(z[,1]," (",z[,2],"-",z[,3],")", sep = ""))
colnames(z)[1]<-"Adjusted Rate Ratio"
colnames(z)[2]<-"LCI"
colnames(z)[3]<-"UCI"
colnames(z)[4]<-"ARR (LCI-UCI)"
View(z)

anova(glm_pos1,glm_pos2,test = "Chisq")
########
# ARR by vacc type
df_cohort_vacc_g_full1 <- df_cohort_vacc_g
df_cohort_vacc_g <- filter(df_cohort_vacc_g, vacc_type=="Mo")
df_cohort_vacc_g <- df_cohort_vacc_g_full1
df_cohort_vacc_g <- filter(df_cohort_vacc_g, vacc_type=="AZ")
df_cohort_vacc_g <- df_cohort_vacc_g_full1
df_cohort_vacc_g <- filter(df_cohort_vacc_g, vacc_type=="PB")

# ARR for individuals BMI groups
df_cohort_vacc_g_full <- df_cohort_vacc_g
df_cohort_vacc_g <- filter(df_cohort_vacc_g, bmi_gp=="18.5-24.9")
df_cohort_vacc_g <- df_cohort_vacc_g_full
df_cohort_vacc_g <- filter(df_cohort_vacc_g, bmi_gp=="<18.5")
df_cohort_vacc_g <- df_cohort_vacc_g_full
df_cohort_vacc_g <- filter(df_cohort_vacc_g, bmi_gp=="25-29.9")
df_cohort_vacc_g <- df_cohort_vacc_g_full
df_cohort_vacc_g <- filter(df_cohort_vacc_g, bmi_gp=="30-39.9")
df_cohort_vacc_g <- df_cohort_vacc_g_full
df_cohort_vacc_g <- filter(df_cohort_vacc_g, bmi_gp=="40+")

z.agg <- pyears(Surv(as.numeric(tstart),as.numeric(tstop),endpt) ~ 
                  pv_period_f+period_f + Sex + age_gp+n_risk_gps+
                  simd2020_sc_quintile+pos_before_start_u+vacc_gap+
                  n_tests_gp+ur_combined, 
                data=df_cohort_vacc_g, weight=ew , scale=365.25, data.frame=TRUE)
glm_pos1 <- glm(event ~ offset(log(pyears)) +pv_period_f+period_f + Sex + age_gp+
                  n_risk_gps+simd2020_sc_quintile+pos_before_start_u+vacc_gap+
                  n_tests_gp+ur_combined, 
                family=poisson, data=z.agg$data)

z<-as.data.frame(round(exp(glm_pos1$coefficients),3))
z<-cbind(z,round(exp(glm_pos1$coefficients-coef(summary(glm_pos1))[,2]*1.96),3))
z<-cbind(z,round(exp(glm_pos1$coefficients+coef(summary(glm_pos1))[,2]*1.96),3))
z<-cbind(z, paste(z[,1]," (",z[,2],"-",z[,3],")", sep = ""))
colnames(z)[1]<-"Adjusted Rate Ratio"
colnames(z)[2]<-"LCI"
colnames(z)[3]<-"UCI"
colnames(z)[4]<-"ARR (LCI-UCI)"
View(z)

# write.csv(z,paste0(project_path,"/obesity_codes/output/ARR_obesity_18-24.csv"))
# write.csv(z,paste0(project_path,"/obesity_codes/output/ARR_obesity_<18.csv"))
# write.csv(z,paste0(project_path,"/obesity_codes/output/ARR_obesity_25-29.csv"))
# write.csv(z,paste0(project_path,"/obesity_codes/output/ARR_obesity_30-39.csv"))
# write.csv(z,paste0(project_path,"/obesity_codes/output/ARR_obesity_>40.csv"))
# 
# write.csv(z,paste0(project_path,"/obesity_codes/output/ARR_obesity_18-24_AZ.csv"))
# write.csv(z,paste0(project_path,"/obesity_codes/output/ARR_obesity_<18_AZ.csv"))
# write.csv(z,paste0(project_path,"/obesity_codes/output/ARR_obesity_25-29_AZ.csv"))
# write.csv(z,paste0(project_path,"/obesity_codes/output/ARR_obesity_30-39_AZ.csv"))
# write.csv(z,paste0(project_path,"/obesity_codes/output/ARR_obesity_>40_AZ.csv"))
# 
# write.csv(z,paste0(project_path,"/obesity_codes/output/ARR_obesity_18-24_PB.csv"))
# write.csv(z,paste0(project_path,"/obesity_codes/output/ARR_obesity_<18_PB.csv"))
# write.csv(z,paste0(project_path,"/obesity_codes/output/ARR_obesity_25-29_PB.csv"))
# write.csv(z,paste0(project_path,"/obesity_codes/output/ARR_obesity_30-39_PB.csv"))
# write.csv(z,paste0(project_path,"/obesity_codes/output/ARR_obesity_>40_PB.csv"))
#################

df_cohort_vacc_g <- mutate(df_cohort_vacc_g,Q_DIAG_cardio=if_else(
  (Q_DIAG_CHD==1|Q_DIAG_PVD==1|Q_DIAG_STROKE==1),1,0))

conditions <- c("Q_DIAG_ASTHMA", "Q_DIAG_DIABETES_2",
                "Q_DIAG_CCF","Q_DIAG_CKD","Q_DIAG_cardio")
ARR <- rbind()
for (i in 1:length(conditions)){#length(conditions)
  temp <- mutate(df_cohort_vacc_g, risk_gp = n_risk_gps)
  temp$risk_gp <- as.numeric(temp$risk_gp)
  temp <- mutate(temp, risk_gp=risk_gp-1)
  temp <- mutate(temp, risk_gp=risk_gp-get(conditions[i]))
  temp$risk_gp<-as.factor(temp$risk_gp)
  tt<-temp[conditions[i]]
  colnames(tt)<-"cond"
  temp<-cbind(temp,tt)
  temp$cond <- as.factor(temp$cond)
  temp$cond <- relevel(temp$cond, ref = "0")
  z.agg <- pyears(Surv(as.numeric(tstart),as.numeric(tstop),endpt) ~ +
                    cond+ period_f+Sex + age_gp+n_risk_gps+
                    simd2020_sc_quintile+pos_before_start_u+vacc_gap+
                    n_tests_gp+ur_combined, 
                  data=temp, weight=ew , scale=365.25, data.frame=TRUE)
  glm_pos <- glm(event ~ offset(log(pyears)) +
                   cond+ period_f+Sex + age_gp+n_risk_gps+
                   simd2020_sc_quintile+pos_before_start_u+vacc_gap+
                   n_tests_gp+ur_combined, 
                  family=poisson, data=z.agg$data)
  #summary(glm_pos)
  z<-as.data.frame(round(exp(glm_pos$coefficients),3))
  z<-cbind(z,round(exp(glm_pos$coefficients-coef(summary(glm_pos))[,2]*1.96),2))
  z<-cbind(z,round(exp(glm_pos$coefficients+coef(summary(glm_pos))[,2]*1.96),2))
  colnames(z)[1]<-"Adjusted Rate Ratio"
  colnames(z)[2]<-"LCI"
  colnames(z)[3]<-"UCI"
  rownames(z)[2]<-conditions[i]
  ARR <- rbind(ARR,z[2,])
  #print(z)
}
ARR<-mutate(ARR, ARR_LCI_UCI = paste(ARR[,1]," (",ARR[,2],"-",ARR[,3],")"))
colnames(ARR)[4]<-"ARR (LCI-UCI)"

write.csv(ARR,paste0(project_path,"/obesity_codes/output/ARR_conds_14-28_days_as_ref_ob.csv"))
write.csv(ARR,paste0(project_path,"/obesity_codes/output/ARR_conds_14-28_days_as_ref_ob_30-39.csv"))
write.csv(ARR,paste0(project_path,"/obesity_codes/output/ARR_conds_14-28_days_as_ref_ob_40+.csv"))

df_cohort_vacc_g_full <- df_cohort_vacc_g
df_cohort_vacc_g <- filter(df_cohort_vacc_g, bmi_gp=="30-39.9")
df_cohort_vacc_g <- df_cohort_vacc_g_full
df_cohort_vacc_g <- filter(df_cohort_vacc_g, bmi_gp=="40+")
#################

p_years<-pyears(Surv(as.numeric(tstart),as.numeric(tstop),endpt)~Sex, scale = 365.25 ,
                data = df_cohort_vacc_g, data.frame = TRUE, weights = ew)
sum(p_years$data$event)*1000/sum(p_years$data$pyears)

#################
########
# interaction analysis age and bmi for nature reviewer
z.agg <- pyears(Surv(as.numeric(tstart),as.numeric(tstop),endpt) ~
                  pv_period_f+period_f + Sex + age_gp+n_risk_gps+
                  bmi_gp+simd2020_sc_quintile+pos_before_start_u+vacc_gap+
                  n_tests_gp+ur_combined,
                data=df_cohort_vacc_g, weight=ew , scale=365.25, data.frame=TRUE)
glm_pos2 <- glm(event ~ offset(log(pyears)) +pv_period_f+period_f + Sex + age_gp+
                  n_risk_gps+bmi_gp+simd2020_sc_quintile+pos_before_start_u+vacc_gap+
                  n_tests_gp+ur_combined+bmi_gp*age_gp, 
                family=poisson, data=z.agg$data)

z<-as.data.frame(round(exp(glm_pos2$coefficients),3))
z<-cbind(z,round(exp(glm_pos2$coefficients-coef(summary(glm_pos2))[,2]*1.96),3))
z<-cbind(z,round(exp(glm_pos2$coefficients+coef(summary(glm_pos2))[,2]*1.96),3))
z<-cbind(z, paste(z[,1]," (",z[,2],"-",z[,3],")", sep = ""))
colnames(z)[1]<-"Adjusted Rate Ratio"
colnames(z)[2]<-"LCI"
colnames(z)[3]<-"UCI"
colnames(z)[4]<-"ARR (LCI-UCI)"
View(z)

anova(glm_pos1,glm_pos2,test = "Chisq")

# model only for 80+ for nature reviewer
df_cohort_vacc_g_full <- df_cohort_vacc_g
df_cohort_vacc_g <- filter(df_cohort_vacc_g, age_gp=="80+")

z.agg <- pyears(Surv(as.numeric(tstart),as.numeric(tstop),endpt) ~ 
                  pv_period_f+period_f + Sex +n_risk_gps+bmi_gp+
                  simd2020_sc_quintile+pos_before_start_u+vacc_gap+
                  n_tests_gp+ur_combined, 
                data=df_cohort_vacc_g, weight=ew , scale=365.25, data.frame=TRUE)
glm_pos1 <- glm(event ~ offset(log(pyears)) +pv_period_f+period_f + Sex +
                  n_risk_gps+bmi_gp+simd2020_sc_quintile+pos_before_start_u+vacc_gap+
                  n_tests_gp+ur_combined, 
                family=poisson, data=z.agg$data)

z<-as.data.frame(round(exp(glm_pos1$coefficients),3))
z<-cbind(z,round(exp(glm_pos1$coefficients-coef(summary(glm_pos1))[,2]*1.96),3))
z<-cbind(z,round(exp(glm_pos1$coefficients+coef(summary(glm_pos1))[,2]*1.96),3))
z<-cbind(z, paste(z[,1]," (",z[,2],"-",z[,3],")", sep = ""))
colnames(z)[1]<-"Adjusted Rate Ratio"
colnames(z)[2]<-"LCI"
colnames(z)[3]<-"UCI"
colnames(z)[4]<-"ARR (LCI-UCI)"
View(z)

# spline recorded bmi - method 2
# run a couple of lines from recorded bmi r file
project_path_vaccine <- paste0(Location,"EAVE/GPanalysis/progs/CR/Vaccine")
rg <- readRDS(paste0(project_path_vaccine,"/output/temp/Qcovid_all.rds"))
rg <- rg %>% dplyr::select(-(Sex:ur6_2016_name))
rg <- filter(rg, !duplicated(EAVE_LINKNO))
rg <- select(rg, EAVE_LINKNO, Q_BMI)
df_cohort_vacc_g <- left_join(df_cohort_vacc_g,rg, by="EAVE_LINKNO")
df_cohort_vacc_g <- mutate(df_cohort_vacc_g,Q_BMI=as.integer(bmi_impute))

z.agg <- pyears(Surv(as.numeric(tstart),as.numeric(tstop),endpt) ~
                  pv_period_f+period_f + Sex + age_gp+n_risk_gps+
                  Q_BMI+simd2020_sc_quintile+pos_before_start_u+vacc_gap+
                  n_tests_gp+ur_combined,
                data=df_cohort_vacc_g, weight=ew , scale=365.25, data.frame=TRUE)
z.agg$data$Q_BMI <- as.numeric(as.character(z.agg$data$Q_BMI))
glm_pos11 <- glm(event ~ offset(log(pyears)) +pv_period_f+period_f + Sex + age_gp+
                  n_risk_gps+pspline(Q_BMI)+simd2020_sc_quintile+pos_before_start_u+vacc_gap+
                  n_tests_gp+ur_combined,
                family=poisson, data=z.agg$data)

model_fit <- glm_pos11
term <- "Q_BMI"
plot_HR <- function(model_fit, term){ 
  # plots hazard ratios for a single term in a fitted model
  
  # hr <- termplot(model_fit, terms = "Q_BMI", se = T, plot = F)
  # var <- names(hr)
  # hr <- hr[[var]]
  
  hr <- termplot(model_fit, se = T, plot = F)
  hr <- hr$Q_BMI
  
  hr <- mutate(hr, ucl = y + 1.96*se,
               lcl = y - 1.96*se) %>%
    mutate_at(c('y', 'ucl', 'lcl'), exp)
  
  hr <- do.call(data.frame,lapply(hr, function(x) replace(x, is.infinite(x),NA)))
  
  # output <- ggplot(data=hr, aes(x=x, y=y)) + geom_line() +
  #   geom_ribbon(aes(ymin=lcl, ymax=ucl), linetype=2, alpha=0.1, fill = 'steelblue')  + 
  #   ylab("Hazard Ratio")
  # 
  # output
  
  lkp_colour_qcovid <- "#fb8c61"
  lkp_colour <- c(
    "reference" = "#000000",
    "estimate"  = "#0a6aa6"
  )
  p_coef_qcovid <-
    hr %>% 
    ggplot(aes(
      x = x,
      ymin = lcl,
      ymax = ucl,
      y = y
    )) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    geom_pointrange(colour = lkp_colour_qcovid) +
    scale_y_continuous(
      name = "Adjusted rate ratio (95% CI)",
      breaks = pretty_breaks()
    ) +
    scale_x_continuous(
      name = "BMI",
      breaks = pretty_breaks()
    ) +
    scale_colour_manual(values = lkp_colour) +
    theme(
      #axis.title.y      = element_blank(),
      strip.placement   = "outside",
      strip.background  = element_blank(),
      strip.text.y      = element_text(face = "bold"),
      strip.text.x = element_text(face = "bold"),
      legend.position   = "none",
      axis.text = element_text(face="bold"),
      axis.text.x = element_text(angle=0,vjust = 0.5),
      plot.title = element_text(size = 11)
    ) +
    coord_cartesian(
      ylim = c(0, 4)
    ) +
    labs(
      #title = "ARRs for "
    )
  p_coef_qcovid
}


####################