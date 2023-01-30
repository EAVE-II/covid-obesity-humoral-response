library(tidyverse)
library(scales)
library(cowplot)

# plot QCovid items ==========================
ARR_data <- read.csv("/conf/EAVE/GPanalysis/progs/UA/second_booster_dose_failures/obesity_codes/output/ARR_conds_14-28_days_as_ref_ob.csv")
ARR_data <- read.csv("/conf/EAVE/GPanalysis/progs/UA/second_booster_dose_failures/obesity_codes/output/ARR_conds_14-28_days_as_ref_ob_30-39.csv")
ARR_data <- read.csv("/conf/EAVE/GPanalysis/progs/UA/second_booster_dose_failures/obesity_codes/output/ARR_conds_14-28_days_as_ref_ob_40+.csv")

expr_coef <- "([\\.0-9]+) \\(([\\.0-9]+), ([\\.0-9]+)\\)"

# conditions
z <- ARR_data %>% 
  mutate(
    term=X,
    doseb_est       = as.numeric(str_replace(Adjusted.Rate.Ratio, expr_coef, "\\1")),
    doseb_conf_low  = as.numeric(str_replace(LCI, expr_coef, "\\2")), 
    doseb_conf_high = as.numeric(str_replace(UCI, expr_coef, "\\3"))
  ) 
t_coef_qcovid <- z[,c(6:9)]

lkp_qcovid <- c(
  "Kidney disease"        = "Q_DIAG_CKD",
  "Diabetes type 2" = "Q_DIAG_DIABETES_2",
  "Heart failure" = "Q_DIAG_CCF",
  "Asthma" = "Q_DIAG_ASTHMA",
  "Cardiovascular diseases" = "Q_DIAG_cardio"
)

lkp_qcovid <- c(
  "Diabetes type 2" = "Q_DIAG_DIABETES_2",
  "Heart failure" = "Q_DIAG_CCF",
  "Kidney disease"        = "Q_DIAG_CKD",
  "Asthma" = "Q_DIAG_ASTHMA",
  "Cardiovascular diseases" = "Q_DIAG_cardio"
)

lkp_colour_qcovid <- "#fb8c61"

p_coef_qcovid <-
  t_coef_qcovid %>% 
  filter(term != "Q_DIAG_CEREBRALPALSY") %>% 
  filter(term %in% lkp_qcovid) %>% 
  mutate(
    term = factor(term, lkp_qcovid, names(lkp_qcovid)),
    term = fct_rev(term)
  ) %>% 
  ggplot(aes(
    x = doseb_est,
    xmin = doseb_conf_low,
    xmax = doseb_conf_high,
    y = term
  )) +
  geom_vline(xintercept = 1) +
  geom_pointrange(colour = lkp_colour_qcovid) +
  scale_x_continuous(
    name = "Adjusted risk ratio (95% CI)",
    breaks = pretty_breaks()
  ) +
  theme(
    axis.title.y      = element_blank(),
    strip.background  = element_blank(),
    axis.text = element_text(face="bold")
  ) +
  coord_cartesian(
    xlim = c(0.5, 2.0)
  ) +
  labs(
    #title=expression(Obese~(BMI~30.0-39.9~kg/m^2~", "~N~"="~"568,420"))
    title = expression(Severly~obese~(BMI~40.0+~kg/m^2~", "~N~"="~"98,706"))
  )

p_coef_qcovid

#####