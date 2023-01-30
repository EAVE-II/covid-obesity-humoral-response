library(tidyverse)
library(scales)
library(cowplot)

ARR_data <- read.csv("/conf/EAVE/GPanalysis/progs/UA/second_booster_dose_failures/output/ARR_obesity_int.csv")
expr_coef <- "([\\.0-9]+) \\(([\\.0-9]+), ([\\.0-9]+)\\)"

### BMI
ARR_data1 <- ARR_data[c(2,3,6,7,9,10,12,13),]
z <- ARR_data1 %>% 
  mutate(
    xvar="vacc_dose_interval",
    xlbl=X,
    doseb_est       = as.numeric(str_replace(Adjusted.Rate.Ratio, expr_coef, "\\1")),
    doseb_conf_low  = as.numeric(str_replace(LCI, expr_coef, "\\2")), 
    doseb_conf_high = as.numeric(str_replace(UCI, expr_coef, "\\3"))
  ) 
t_coef_overall <- z[,c(6:10)]

# tidy plot of demographics and health characteristics =======


t_coef_ref <- tribble(
  ~xvar,                ~xlbl, ~doseb_est, ~doseb_conf_low, ~doesb_conf_high,
  "vacc_dose_interval", "pv_period_fv2_2:9",  1.000, 1.000, 1.000
)

lkp_xvar <- c(
  "Dose\ninterval"           = "vacc_dose_interval"
)

lkp_xlbl <- c(
  "Day 14-69" = "pv_period_fv2_2:9",
  "Day 70-139" = "pv_period_fv2_10:19",
  "Day 140+" = "pv_period_fv2_20+",
  "Day 14-34 Mo" = "pv_period_fv3_2:4_Mo",
  "Day 35-55 Mo" = "pv_period_fv3_5:8_Mo",
  "Day 56+ Mo" = "pv_period_fv3_9+_Mo",
  "Day 35-55 PB" = "pv_period_fv3_5:8_PB",
  "Day 14-34 PB" = "pv_period_fv3_2:4_PB",
  "Day 56+ PB"   = "pv_period_fv3_9+_PB"
)

lkp_colour <- c(
  "reference" = "#000000",
  "estimate"  = "#0a6aa6"
)

p_coef_overall <-
  bind_rows(
    t_coef_overall %>% mutate(est_type = "estimate"),
    t_coef_ref     %>% mutate(est_type = "reference")
  ) 
p_coef_overall$doseb_conf_high<-coalesce(p_coef_overall$doseb_conf_high,
                                         p_coef_overall$doesb_conf_high)
p_coef_overall$doesb_conf_high<-NULL
p_coef_overall_o <-p_coef_overall%>% 
  filter(
    xvar %in% lkp_xvar,
    xlbl %in% lkp_xlbl
  ) %>% 
  mutate(
    xvar = factor(xvar, lkp_xvar, names(lkp_xvar)),
    xlbl = factor(xlbl, lkp_xlbl, names(lkp_xlbl))
  ) %>% 
  ggplot(aes(
    x = xlbl,
    ymin = doseb_conf_low,
    ymax = doseb_conf_high,
    y = doseb_est,
    colour = est_type
  )) +
  geom_hline(yintercept = 1,
             linetype = "dashed") +
  geom_pointrange() +
  scale_y_continuous(
    name = "Adjusted risk ratio (95% CI)",
    breaks = pretty_breaks()
  ) +
  scale_colour_manual(values = lkp_colour) +
  theme(
    axis.title.x      = element_blank(),
    strip.placement   = "outside",
    strip.background  = element_blank(),
    strip.text.y      = element_text(face = "bold"),
    strip.text.x = element_text(face = "bold"),
    legend.position   = "none",
    axis.text = element_text(face="bold"),
    axis.text.x = element_text(angle=90,vjust = 0.5)
  ) +
  coord_cartesian(
    ylim = c(0,9)
  ) +
  labs(
    title = "BMI all categories"
  )

p_coef_overall_o
######

# <18.5 BMI
ARR_data1 <- ARR_data[c(73,74,77,78,80,81,83,84),]
z <- ARR_data1 %>% 
  mutate(
    xvar="vacc_dose_interval",
    xlbl=X,
    doseb_est       = as.numeric(str_replace(Adjusted.Rate.Ratio, expr_coef, "\\1")),
    doseb_conf_low  = as.numeric(str_replace(LCI, expr_coef, "\\2")), 
    doseb_conf_high = as.numeric(str_replace(UCI, expr_coef, "\\3"))
  ) 
z[,7] <- t_coef_overall[,2]
t_coef_overall <- z[,c(6:10)]

# tidy plot of demographics and health characteristics =======


t_coef_ref <- tribble(
  ~xvar,                ~xlbl, ~doseb_est, ~doseb_conf_low, ~doesb_conf_high,
  "vacc_dose_interval", "pv_period_fv2_2:9",  1.000, 1.000, 1.000
)


p_coef_overall <-
  bind_rows(
    t_coef_overall %>% mutate(est_type = "estimate"),
    t_coef_ref     %>% mutate(est_type = "reference")
  ) 
p_coef_overall$doseb_conf_high<-coalesce(p_coef_overall$doseb_conf_high,
                                         p_coef_overall$doesb_conf_high)
p_coef_overall$doesb_conf_high<-NULL
p_coef_overall_o1 <-p_coef_overall%>% 
  filter(
    xvar %in% lkp_xvar,
    xlbl %in% lkp_xlbl
  ) %>% 
  mutate(
    xvar = factor(xvar, lkp_xvar, names(lkp_xvar)),
    xlbl = factor(xlbl, lkp_xlbl, names(lkp_xlbl))
  ) %>% 
  ggplot(aes(
    x = xlbl,
    ymin = doseb_conf_low,
    ymax = doseb_conf_high,
    y = doseb_est,
    colour = est_type
  )) +
  geom_hline(yintercept = 1,
             linetype = "dashed") +
  geom_pointrange() +
  scale_y_continuous(
    name = "Adjusted risk ratio (95% CI)",
    breaks = pretty_breaks()
  ) +
  scale_colour_manual(values = lkp_colour) +
  theme(
    axis.title.x      = element_blank(),
    strip.placement   = "outside",
    strip.background  = element_blank(),
    strip.text.y      = element_text(face = "bold"),
    strip.text.x = element_text(face = "bold"),
    legend.position   = "none",
    axis.text = element_text(face="bold"),
    axis.text.x = element_text(angle=90,vjust = 0.5)
  ) +
  coord_cartesian(
    ylim = c(0,9)
  ) +
  labs(
    title = "BMI <18.5"
  )

p_coef_overall_o1
####

# 25-29 BMI
ARR_data1 <- ARR_data[c(85,86,89,90,92,93,95,96),]
z <- ARR_data1 %>% 
  mutate(
    xvar="vacc_dose_interval",
    xlbl=X,
    doseb_est       = as.numeric(str_replace(Adjusted.Rate.Ratio, expr_coef, "\\1")),
    doseb_conf_low  = as.numeric(str_replace(LCI, expr_coef, "\\2")), 
    doseb_conf_high = as.numeric(str_replace(UCI, expr_coef, "\\3"))
  ) 
z[,7] <- t_coef_overall[,2]
t_coef_overall <- z[,c(6:10)]

# tidy plot of demographics and health characteristics =======


t_coef_ref <- tribble(
  ~xvar,                ~xlbl, ~doseb_est, ~doseb_conf_low, ~doesb_conf_high,
  "vacc_dose_interval", "pv_period_fv2_2:9",  1.000, 1.000, 1.000
)

p_coef_overall <-
  bind_rows(
    t_coef_overall %>% mutate(est_type = "estimate"),
    t_coef_ref     %>% mutate(est_type = "reference")
  ) 
p_coef_overall$doseb_conf_high<-coalesce(p_coef_overall$doseb_conf_high,
                                         p_coef_overall$doesb_conf_high)
p_coef_overall$doesb_conf_high<-NULL
p_coef_overall_o3 <-p_coef_overall%>% 
  filter(
    xvar %in% lkp_xvar,
    xlbl %in% lkp_xlbl
  ) %>% 
  mutate(
    xvar = factor(xvar, lkp_xvar, names(lkp_xvar)),
    xlbl = factor(xlbl, lkp_xlbl, names(lkp_xlbl))
  ) %>% 
  ggplot(aes(
    x = xlbl,
    ymin = doseb_conf_low,
    ymax = doseb_conf_high,
    y = doseb_est,
    colour = est_type
  )) +
  geom_hline(yintercept = 1,
             linetype = "dashed") +
  geom_pointrange() +
  scale_y_continuous(
    name = "Adjusted risk ratio (95% CI)",
    breaks = pretty_breaks()
  ) +
  scale_colour_manual(values = lkp_colour) +
  theme(
    axis.title.x      = element_blank(),
    strip.placement   = "outside",
    strip.background  = element_blank(),
    strip.text.y      = element_text(face = "bold"),
    strip.text.x = element_text(face = "bold"),
    legend.position   = "none",
    axis.text = element_text(face="bold"),
    axis.text.x = element_text(angle=90,vjust = 0.5)
  ) +
  coord_cartesian(
    ylim = c(0,9)
  ) +
  labs(
    title = "BMI 25.0-29.9"
  )

p_coef_overall_o3
####

# 30-39 BMI
ARR_data1 <- ARR_data[c(97,98,101,102,104,105,107,108),]
z <- ARR_data1 %>% 
  mutate(
    xvar="vacc_dose_interval",
    xlbl=X,
    doseb_est       = as.numeric(str_replace(Adjusted.Rate.Ratio, expr_coef, "\\1")),
    doseb_conf_low  = as.numeric(str_replace(LCI, expr_coef, "\\2")), 
    doseb_conf_high = as.numeric(str_replace(UCI, expr_coef, "\\3"))
  ) 
z[,7] <- t_coef_overall[,2]
t_coef_overall <- z[,c(6:10)]

# tidy plot of demographics and health characteristics =======


t_coef_ref <- tribble(
  ~xvar,                ~xlbl, ~doseb_est, ~doseb_conf_low, ~doesb_conf_high,
  "vacc_dose_interval", "pv_period_fv2_2:9",  1.000, 1.000, 1.000
)

p_coef_overall <-
  bind_rows(
    t_coef_overall %>% mutate(est_type = "estimate"),
    t_coef_ref     %>% mutate(est_type = "reference")
  ) 
p_coef_overall$doseb_conf_high<-coalesce(p_coef_overall$doseb_conf_high,
                                         p_coef_overall$doesb_conf_high)
p_coef_overall$doesb_conf_high<-NULL
p_coef_overall_o4 <-p_coef_overall%>% 
  filter(
    xvar %in% lkp_xvar,
    xlbl %in% lkp_xlbl
  ) %>% 
  mutate(
    xvar = factor(xvar, lkp_xvar, names(lkp_xvar)),
    xlbl = factor(xlbl, lkp_xlbl, names(lkp_xlbl))
  ) %>% 
  ggplot(aes(
    x = xlbl,
    ymin = doseb_conf_low,
    ymax = doseb_conf_high,
    y = doseb_est,
    colour = est_type
  )) +
  geom_hline(yintercept = 1,
             linetype = "dashed") +
  geom_pointrange() +
  scale_y_continuous(
    name = "Adjusted risk ratio (95% CI)",
    breaks = pretty_breaks()
  ) +
  scale_colour_manual(values = lkp_colour) +
  theme(
    axis.title.x      = element_blank(),
    strip.placement   = "outside",
    strip.background  = element_blank(),
    strip.text.y      = element_text(face = "bold"),
    strip.text.x = element_text(face = "bold"),
    legend.position   = "none",
    axis.text = element_text(face="bold"),
    axis.text.x = element_text(angle=90,vjust = 0.5)
  ) +
  coord_cartesian(
    ylim = c(0,9)
  ) +
  labs(
    title = "BMI 30.0-39.9"
  )

p_coef_overall_o4
####

# >40 BMI
ARR_data1 <- ARR_data[c(109,110,113,114,116,117,119,120),]
z <- ARR_data1 %>% 
  mutate(
    xvar="vacc_dose_interval",
    xlbl=X,
    doseb_est       = as.numeric(str_replace(Adjusted.Rate.Ratio, expr_coef, "\\1")),
    doseb_conf_low  = as.numeric(str_replace(LCI, expr_coef, "\\2")), 
    doseb_conf_high = as.numeric(str_replace(UCI, expr_coef, "\\3"))
  ) 
z[,7] <- t_coef_overall[,2]
t_coef_overall <- z[,c(6:10)]

# tidy plot of demographics and health characteristics =======


t_coef_ref <- tribble(
  ~xvar,                ~xlbl, ~doseb_est, ~doseb_conf_low, ~doesb_conf_high,
  "vacc_dose_interval", "pv_period_fv2_2:9",  1.000, 1.000, 1.000
)


p_coef_overall <-
  bind_rows(
    t_coef_overall %>% mutate(est_type = "estimate"),
    t_coef_ref     %>% mutate(est_type = "reference")
  ) 
p_coef_overall$doseb_conf_high<-coalesce(p_coef_overall$doseb_conf_high,
                                         p_coef_overall$doesb_conf_high)
p_coef_overall$doesb_conf_high<-NULL
p_coef_overall_o5 <-p_coef_overall%>% 
  filter(
    xvar %in% lkp_xvar,
    xlbl %in% lkp_xlbl
  ) %>% 
  mutate(
    xvar = factor(xvar, lkp_xvar, names(lkp_xvar)),
    xlbl = factor(xlbl, lkp_xlbl, names(lkp_xlbl))
  ) %>% 
  ggplot(aes(
    x = xlbl,
    ymin = doseb_conf_low,
    ymax = doseb_conf_high,
    y = doseb_est,
    colour = est_type
  )) +
  geom_hline(yintercept = 1,
             linetype = "dashed") +
  geom_pointrange() +
  scale_y_continuous(
    name = "Adjusted risk ratio (95% CI)",
    breaks = pretty_breaks()
  ) +
  scale_colour_manual(values = lkp_colour) +
  theme(
    axis.title.x      = element_blank(),
    strip.placement   = "outside",
    strip.background  = element_blank(),
    strip.text.y      = element_text(face = "bold"),
    strip.text.x = element_text(face = "bold"),
    legend.position   = "none",
    axis.text = element_text(face="bold"),
    axis.text.x = element_text(angle=90,vjust = 0.5)
  ) +
  coord_cartesian(
    ylim = c(0,9)
  ) +
  labs(
    title = "BMI 40+"
  )

p_coef_overall_o5
####
# combine =====

p_coef <- plot_grid(p_coef_overall_o,p_coef_overall_o1,
                    p_coef_overall_o3,p_coef_overall_o4,p_coef_overall_o5)
print(p_coef)
ggsave(
  plot = p_coef,
  filename = "/conf/EAVE/GPanalysis/progs/UA/p_vacc_breakthrough_coef.png",
  width = 35,
  height = 20,
  units = "cm"
)


