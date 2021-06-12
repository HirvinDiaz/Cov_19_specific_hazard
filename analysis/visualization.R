#### Visualization ####
library(writexl)


#### Efficient frontier Not Intubated ####
gg_eff_hosp <- ggplot(data = filter(df_cea_psa_Hosp, 
                     Status == "ND"), 
       aes(x = Effect, y = Cost/100)) +
  geom_point(size = 2, shape = 15) +
  geom_text(aes(label = Strategy, 
                vjust = -1, hjust =1), size = 3.1)+
  geom_line(size = 1) +
  geom_point(data = filter(df_cea_psa_Hosp, 
                           Status == "ED"), size = 2, shape = 15) +
  geom_text(data = filter(df_cea_psa_Hosp, 
                           Status == "ED"), aes(label = Strategy, 
                 vjust = -1, hjust =1), size = 3.1)+
  theme(plot.title = element_text(face = "bold", 
                                  size = 12,
                                  family =, 
                                  hjust = 0.5),
        axis.line.x = element_line(color="#757373", size = 1),
        axis.line.y = element_line(color="#757373", size = 1),
        plot.caption = element_text(hjust = 0,
                                    colour = "#777777",
                                    size = 10),
        panel.background = element_rect(fill = "white", 
                                        colour = "white", 
                                        size = 0.15, 
                                        linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, 
                                        linetype = 'solid',
                                        colour = "gray"),
        legend.position = "bottom") +
  scale_y_continuous(labels=dollar_format(prefix="$"),breaks = number_ticks(4), limits = c(2500,4100)) +
  scale_x_continuous(breaks = number_ticks(4), limits = c(5.0,6.8)) +
  labs(title = "Hospitalized, not intubated",
       x = "Effect (QALYs)",
       y = "Cost (in $1,000 Mexican Pesos)")

#### Efficient frontier Intubated ####
gg_eff_ICU <- ggplot(data = df_cea_psa_ICU, 
       aes(x = Effect, y = Cost/100)) +
  geom_point(size = 4, shape = 15) +
  geom_text(aes(label = Strategy, 
                vjust = -1, hjust =1), size = 3.1)+
  geom_line(size = 1) +
  theme(plot.title = element_text(face = "bold", 
                                  size = 12,
                                  family =, hjust = 0.5),
        axis.line.x = element_line(color="#757373", size = 1),
        axis.line.y = element_line(color="#757373", size = 1),
        plot.caption = element_text(hjust = 0,
                                    colour = "#777777",
                                    size = 10),
        panel.background = element_rect(fill = "white", 
                                        colour = "white", 
                                        size = 0.15, 
                                        linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, 
                                        linetype = 'solid',
                                        colour = "gray"),
        legend.position = "bottom") +
  scale_y_continuous(labels=dollar_format(prefix="$"),breaks = number_ticks(4), limits = c(6800,7250)) +
  scale_x_continuous(breaks = number_ticks(4), limits = c(1.0,2.8)) +
  labs(title = "Hospitalized, Intubated",
       x = "Effect (QALYs)",
       y = "Cost (in $1,000 Mexican Pesos)")

figure <- ggarrange(gg_eff_hosp + rremove("ylab") + rremove("xlab"), 
                    gg_eff_ICU + rremove("ylab") + rremove("xlab"), 
                    ncol = 2, 
                    common.legend = TRUE)

annotate_figure(figure, left = text_grob("Costs (in $1,000 Mexican Pesos)", rot = 90, vjust = 1),
                bottom = text_grob("Effect (QALYs)"))


ggsave(paste0("figs/Efficient frontier",
              format(Sys.Date(), "%F"), ".png"), 
       width =10, height = 5)


#### PSA Cost-effectiveness Intubated ####
gg_hosp_CEAs <-  ggplot(data = df_CEAs_hosp_long, 
       aes(x = Effect, 
           y = Costs/100, 
           color = Strategy))+
  geom_point(alpha = 0.35)+ 
  theme(plot.title = element_text(face = "bold", 
                                              size = 12,
                                              family =, hjust = 0.5),
                    axis.line.x = element_line(color="#757373", size = 1),
                    axis.line.y = element_line(color="#757373", size = 1),
                    plot.caption = element_text(hjust = 0,
                                                colour = "#777777",
                                                size = 10),
                    panel.background = element_rect(fill = "white", 
                                                    colour = "white", 
                                                    size = 0.15, 
                                                    linetype = "solid"),
                    panel.grid.major = element_line(size = 0.15, 
                                                    linetype = 'solid',
                                                    colour = "gray"),
                    legend.position = "top",
                    legend.text = element_text(size = 8)) +
  scale_y_continuous(labels=dollar_format(prefix="$"),breaks = number_ticks(4)) +
  scale_x_continuous(breaks = number_ticks(4)) +
  scale_color_manual(values=c ("#474140", "#eb542f", "#265796")) +
  labs(title = "Hospitalized, Not intubated",
       x = "Effect (QALYs)",
       y = "Cost (in $1,000 Mexican Pesos)")


gg_hosp <- gg_hosp_CEAs + 
  geom_label(x = 6.7, 
             y = 4900, 
             label = "Mean Remdesivir & Baricitinib\n6.70, $401,300",
             label.padding = unit(0.15, "lines"), # Rectangle size around label
             label.size = 0.17,
             color = "#4d8ee3",
             fill="#f7fafc",
             size = 3.5)+
  geom_point(x = 6.7, 
             y = 4012,
             alpha = 1.5,
             color = "#265796",
             size = 5,
             shape = 17,
             fill = "#265796")+
  geom_label(x = 5.97, 
             y = 4180, 
             label = "Mean Remdesivir\n5.97, $332,000",
             label.padding = unit(0.15, "lines"), # Rectangle size around label
             label.size = 0.17,
             color = "#eb542f",
             fill="#f7fafc",
             size = 3.5)+
  geom_point(x = 5.97, 
             y = 3320,
             alpha = 1.5,
             color = "#eb542f",
             size = 5,
             shape = 17,
             fill = "#eb542f")+
  geom_label(x = 5.27, 
             y = 3400, 
             label = "Mean No treatment\n5.27, $250,000",
             label.padding = unit(0.15, "lines"), # Rectangle size around label
             label.size = 0.17,
             color = "#474140",
             fill="#f7fafc",
             size = 3.5)+
  geom_point(x = 5.27, 
             y = 2500,
             alpha = 1.5,
             color = "#474140",
             size = 5,
             shape = 17,
             fill = "#474140")
  
ggsave(paste0("figs/CEAs_PSA_hosp",
            format(Sys.Date(), "%F"), ".png"), 
     width = 7, height = 5)


gg_icu_CEAs <- ggplot(data = df_CEAs_icu_long, 
       aes(x = Effect, 
           y = Costs/100, 
           color = Strategy))+
  geom_point(alpha = 0.3)+ 
  theme(plot.title = element_text(face = "bold", 
                                  size = 12,
                                  family =, hjust = 0.5),
        axis.line.x = element_line(color="#757373", size = 1),
        axis.line.y = element_line(color="#757373", size = 1),
        plot.caption = element_text(hjust = 0,
                                    colour = "#777777",
                                    size = 10),
        panel.background = element_rect(fill = "white", 
                                        colour = "white", 
                                        size = 0.15, 
                                        linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, 
                                        linetype = 'solid',
                                        colour = "gray"),
        legend.position = "top",
        legend.text = element_text(size = 8)) +
  scale_y_continuous(labels=dollar_format(prefix="$"), limits = c(0,50000)) +
  scale_x_continuous(breaks = number_ticks(4)) +
  scale_color_manual(values=c ("#474140", "#265796")) +
  labs(title = "Hospitalized, Intubated",
       x = "Effect (QALYs)",
       y = "Cost (in $1,000 Mexican Pesos)")

gg_icu <- gg_icu_CEAs + 
  geom_label(x = 1.57, 
             y = 11700, 
             label = "Mean No treatment\n1.45, $684,900",
             label.padding = unit(0.15, "lines"), # Rectangle size around label
             label.size = 0.17,
             color = "#474140",
             fill="#f7fafc",
             size = 3.5)+
  geom_point(x = 1.45, 
             y = 6800,
             alpha = 1.5,
             color = "#474140",
             size = 5,
             shape = 17,
             fill = "#474140")+
  geom_label(x = 2.73, 
             y = 12100, 
             label = "Mean Dexamethasone\n2.73, $719,500",
             label.padding = unit(0.15, "lines"), # Rectangle size around label
             label.size = 0.17,
             color = "#265796",
             fill="#f7fafc",
             size = 3.5)+
  geom_point(x = 2.73, 
             y = 7200,
             alpha = 1.5,
             color = "#265796",
             size = 5,
             shape = 17,
             fill = "#265796")

ggarrange(gg_hosp, gg_icu, ncol = 1)

ggsave(paste0("figs/CEAs_labels",
              format(Sys.Date(), "%F"), ".png"), 
       width = 7, height = 10)

# figure <- ggarrange(gg_hosp_CEAs + rremove("ylab") + rremove("xlab"), 
#                     gg_icu_CEAs + rremove("ylab") + rremove("xlab"), 
#                     ncol = 2, 
#                     common.legend = FALSE)
# 
# annotate_figure(figure, left = text_grob("Costs (in $1,000 Mexican Pesos)", rot = 90, vjust = 1),
#                 bottom = text_grob("Effect (QALYs)"))

ggsave(paste0("figs/CEAs_PSA_icu",
              format(Sys.Date(), "%F"), ".png"), 
       width = 7, height = 5)


#### Cost-Effectiveness Acceptability Curve ####

df_frontier <- ceac_obj_hosp %>% 
  filter(On_Frontier == "TRUE") %>% 
  select(WTP, Proportion) %>% 
  mutate(Strategy = "Frontier")

# Not Intubated
CEAC_hosp <- ggplot()+ 
  geom_line(data = ceac_obj_hosp, 
            aes(x = WTP, 
                y = Proportion, 
                linetype = Strategy, 
                color = Strategy), 
            size = 1) + 
  geom_point(data = df_frontier, 
                           aes_(x = as.name("WTP"), 
                               y = as.name("Proportion"),
                               shape = as.name("Strategy")),
             size = 3)+
  # geom_vline(xintercept = PIB_pc, linetype = "dotted") +
  theme(plot.title = element_text(face = "bold", 
                                  size = 12,
                                  family =, 
                                  hjust = 0.5),
        axis.line.x = element_line(color="#ebe4e4", size = 1),
        axis.line.y = element_line(color="#ebe4e4", size = 1),
        plot.caption = element_text(hjust = 0,
                                    colour = "#777777",
                                    size = 10),
        panel.background = element_rect(fill = "white", 
                                        colour = "white", 
                                        size = 0.15, 
                                        linetype = "dashed"),
        panel.grid.major = element_line(size = 0.15, 
                                        linetype = 'solid',
                                        colour = "#ebe4e4"),
        legend.position = c(0.75,0.48), 
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8.5))+
  scale_y_continuous(breaks = number_ticks(4)) +
  scale_x_continuous(breaks = number_ticks(4), 
                     labels=dollar_format(prefix="$")) +
  scale_color_manual(values=c ("#eb542f", "#2cc9de", "#575c59")) +
  scale_shape_manual(name = NULL, values = 0, labels = "Frontier") +
  labs(title = "Hospitalized, Not intubated",
       x = "Willingness-to-Pay Threshold, $ Mexican pesos/QALY",
       y = "Probability of being Cost-Effective")

ggsave(paste0("figs/CEAC_not_int_new",
              format(Sys.Date(), "%F"), ".png"), 
       width = 7, height = 5)

# Intubated
df_frontier_icu <- ceac_obj_icu %>% 
  filter(On_Frontier == "TRUE") %>% 
  select(WTP, Proportion) %>% 
  mutate(Strategy = "Frontier")

ceac_obj_icu <- ceac_obj_icu %>% 
  mutate(Strategy = fct_relevel(Strategy, 
                                "No treatment", "Dexamethasone"))

CEAC_icu <- ggplot()+ 
  geom_line(data = ceac_obj_icu, 
            aes(x = WTP, 
                y = Proportion, 
                linetype = Strategy, 
                color = Strategy), 
            size = 1) + 
  geom_point(data = df_frontier_icu, 
             aes_(x = as.name("WTP"), 
                  y = as.name("Proportion"),
                  shape = as.name("Strategy")),
             size = 3)+
  theme(plot.title = element_text(face = "bold", 
                                  size = 12,
                                  family =, 
                                  hjust = 0.5),
        axis.line.x = element_line(color="#ebe4e4", size = 1),
        axis.line.y = element_line(color="#ebe4e4", size = 1),
        plot.caption = element_text(hjust = 0,
                                    colour = "#777777",
                                    size = 10),
        panel.background = element_rect(fill = "white", 
                                        colour = "white", 
                                        size = 0.15, 
                                        linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, 
                                        linetype = 'solid',
                                        colour = "#ebe4e4"),
        legend.position = c(0.8,0.61), 
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8.5))+
  scale_y_continuous(breaks = number_ticks(4)) +
  scale_x_continuous(breaks = number_ticks(4), labels=dollar_format(prefix="$")) +
  scale_color_manual(values=c ("#eb542f", "#575c59")) +
  scale_shape_manual(name = NULL, values = 0, labels = "Frontier") +
  labs(title = "Hospitalized, Intubated",
       x = "Willingness-to-Pay Threshold, $ Mexican pesos/QALY",
       y = "Probability of being Cost-Effective")


ggarrange(CEAC_hosp, CEAC_icu, ncol = 1)

ggsave(paste0("figs/CEACs_label",
              format(Sys.Date(), "%F"), ".png"), 
       width = 7, height = 10)

#### Parameter distribution ####
# Hospitalization time
ggplot(data = df_hosp_time_long, 
       aes(x = Days)) + 
  geom_histogram(color = "#575c59", binwidth=7) +
  facet_wrap(~Cohort, scales = "free") +
  theme(plot.title = element_text(face = "bold", 
                                  size = 12,
                                  family =, hjust = 0),
        axis.line.x = element_line(color="gray", size = 1),
        axis.line.y = element_line(color="gray", size = 1),
        plot.caption = element_text(hjust = 0,
                                    colour = "#777777",
                                    size = 10),
        panel.background = element_rect(fill = "white", 
                                        colour = "white", 
                                        size = 0.15, 
                                        linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, 
                                        linetype = 'solid',
                                        colour = "gray"))+
  scale_y_continuous(breaks = number_ticks(5)) +
  scale_x_continuous(breaks = number_ticks(5)) +
  #scale_color_manual(values=c ("#eb542f", "#575c59")) +
  labs(title = "Hospitalization time",
       x = "Days",
       y = " ")

ggsave(paste0("figs/params_Hosp_time",
              format(Sys.Date(), "%F"), ".png"), 
       width = 7, height = 5)

# Treatment Effect
ggplot(data = df_effects_params_long, 
       aes(x = HR)) + 
  geom_histogram(color = "#575c59", binwidth=0.1) +
  facet_wrap(~Cohort, scales = "free") +
  theme(plot.title = element_text(face = "bold", 
                                  size = 12,
                                  family =, hjust = 0),
        axis.line.x = element_line(color="gray", size = 1),
        axis.line.y = element_line(color="gray", size = 1),
        plot.caption = element_text(hjust = 0,
                                    colour = "#777777",
                                    size = 10),
        panel.background = element_rect(fill = "white", 
                                        colour = "white", 
                                        size = 0.15, 
                                        linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, 
                                        linetype = 'solid',
                                        colour = "gray"))+
  #scale_y_continuous(breaks = number_ticks(5)) +
  #scale_x_continuous(breaks = number_ticks(5)) +
  #scale_color_manual(values=c ("#eb542f", "#575c59")) +
  labs(title = "Treatment effects",
       x = "Hazard Ratio",
       y = " ")

ggsave(paste0("figs/params_treat_effts",
              format(Sys.Date(), "%F"), ".png"), 
       width = 7, height = 5)


# Hospitalization Costs
ggplot(data = df_costs_params_long, 
       aes(x = Costs)) + 
  geom_histogram(color = "#575c59") +
  facet_wrap(~Cohort, scales = "free") +
  theme(plot.title = element_text(face = "bold", 
                                  size = 12,
                                  family =, hjust = 0),
        axis.line.x = element_line(color="gray", size = 1),
        axis.line.y = element_line(color="gray", size = 1),
        plot.caption = element_text(hjust = 0,
                                    colour = "#777777",
                                    size = 10),
        panel.background = element_rect(fill = "white", 
                                        colour = "white", 
                                        size = 0.15, 
                                        linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, 
                                        linetype = 'solid',
                                        colour = "gray"))+
  #scale_y_continuous(breaks = number_ticks(5)) +
  scale_x_continuous(breaks = number_ticks(3), labels=dollar_format(prefix="$")) +
  #scale_color_manual(values=c ("#eb542f", "#575c59")) +
  labs(title = "Hospitalization Costs",
       x = "Costs $",
       y = " ")

ggsave(paste0("figs/params_hosp_costs",
              format(Sys.Date(), "%F"), ".png"), 
       width = 7, height = 5)


write_csv(df_QALEs, path = "data/df_QALEs_csv.csv")

#### EVPI ####

# Not Intubated
GG_EVPI_hosp <- ggplot()+ 
  geom_line(data = filter(elc_obj_hosp_long,
                          Strategy != "Frontier & EVPI"), 
            aes(x = WTP, 
                y = Value, 
                linetype = Strategy, 
                color = Strategy), 
            size = 1) + 
  geom_point(data = filter(elc_obj_hosp_long,
                           Strategy == "Frontier & EVPI"), 
             aes_(x = as.name("WTP"), 
                  y = as.name("Value"),
                  shape = as.name("Strategy")),
             size = 3)+
  # geom_vline(xintercept = PIB_pc, linetype = "dotted") +
  theme(plot.title = element_text(face = "bold", 
                                  size = 12,
                                  family =, 
                                  hjust = 0.5),
        axis.line.x = element_line(color="#ebe4e4", size = 1),
        axis.line.y = element_line(color="#ebe4e4", size = 1),
        plot.caption = element_text(hjust = 0,
                                    colour = "#777777",
                                    size = 10),
        panel.background = element_rect(fill = "white", 
                                        colour = "white", 
                                        size = 0.15, 
                                        linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, 
                                        linetype = 'solid',
                                        colour = "#ebe4e4"),
        legend.position = c(0.25,0.77), 
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8.5))+
  scale_y_continuous(breaks = number_ticks(4), labels=dollar_format(prefix="$")) +
  scale_x_continuous(breaks = number_ticks(4), 
                     labels=dollar_format(prefix="$")) +
  scale_color_manual(values=c ("#575c59", "#2cc9de", "#eb542f")) +
  scale_shape_manual(name = NULL, values = 0, labels = "Frontier & EVPI") +
  labs(title = "Hospitalized, Not intubated",
       x = "Willingness-to-Pay Threshold, $ Mexican pesos/QALY",
       y = "Expected Loss ($)")

ggsave(paste0("figs/EVPI_not_int",
              format(Sys.Date(), "%F"), ".png"), 
       width = 7, height = 5)

# Intubated
GG_EVPI_icu <- ggplot()+ 
  geom_line(data = filter(elc_obj_icu_long,
                          Strategy != "Frontier & EVPI"), 
            aes(x = WTP, 
                y = Value, 
                linetype = Strategy, 
                color = Strategy), 
            size = 1) + 
  geom_point(data = filter(elc_obj_icu_long,
                           Strategy == "Frontier & EVPI"), 
             aes_(x = as.name("WTP"), 
                  y = as.name("Value"),
                  shape = as.name("Strategy")),
             size = 3)+
  # geom_vline(xintercept = PIB_pc, linetype = "dotted") +
  theme(plot.title = element_text(face = "bold", 
                                  size = 12,
                                  family =, 
                                  hjust = 0.5),
        axis.line.x = element_line(color="#ebe4e4", size = 1),
        axis.line.y = element_line(color="#ebe4e4", size = 1),
        plot.caption = element_text(hjust = 0,
                                    colour = "#777777",
                                    size = 10),
        panel.background = element_rect(fill = "white", 
                                        colour = "white", 
                                        size = 0.15, 
                                        linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, 
                                        linetype = 'solid',
                                        colour = "#ebe4e4"),
        legend.position = c(0.25,0.77), 
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8.5))+
  scale_y_continuous(breaks = number_ticks(4), labels=dollar_format(prefix="$")) +
  scale_x_continuous(breaks = number_ticks(4), 
                     labels=dollar_format(prefix="$")) +
  scale_color_manual(values=c ("#eb542f", "#575c59", "#eb542f")) +
  scale_shape_manual(name = NULL, values = 0, labels = "Frontier & EVPI") +
  labs(title = "Hospitalized, Intubated",
       x = "Willingness-to-Pay Threshold, $ Mexican pesos/QALY",
       y = "Expected Loss ($)")

ggarrange(GG_EVPI_hosp, GG_EVPI_icu, ncol = 1)

ggsave(paste0("figs/EVPI",
              format(Sys.Date(), "%F"), ".png"), 
       width = 7, height = 10)

#### Deterministic Sensitivity Analysis ####
o_hosp$parameter[o_hosp$parameter == "hr_Remd"] <- "HR Remdesivir" 
o_hosp$parameter[o_hosp$parameter == "hr_RemdBa_vs_Remd"] <- "HR Baricitinib" 
o_hosp$parameter[o_hosp$parameter == "c_hosp"] <- "Hospitalization Costs" 
o_hosp$parameter[o_hosp$parameter == "mean_days"] <- "Mean Days" 

o_hosp$strategy[o_hosp$strategy == "No.Treatment"] <- "No Treatment" 
o_hosp$strategy[o_hosp$strategy == "Remdesivir.Baricitinib"] <- "Remdesivir & Baricitinib" 


unique(o_hosp$parameter)



GG_ow <- ggplot(data = filter(o_hosp,
                     parameter == "Hospitalization Costs" | 
                     parameter == "HR Baricitinib" | 
                     parameter == "Mean Days" | 
                     parameter == "HR Remdesivir"), 
       aes(x = param_val, 
           y = outcome_val/100, 
           color = strategy))+ 
  geom_line(size = 1) + 
  facet_wrap(~parameter, ncol = 2, scales = "free") +
  # geom_vline(xintercept = PIB_pc, linetype = "dotted") +
  theme(plot.title = element_text(face = "bold", 
                                  size = 12,
                                  family =, 
                                  hjust = 0.5),
        axis.line.x = element_line(color="#ebe4e4", size = 1),
        axis.line.y = element_line(color="#ebe4e4", size = 1),
        plot.caption = element_text(hjust = 0,
                                    colour = "#777777",
                                    size = 10),
        panel.background = element_rect(fill = "white", 
                                        colour = "white", 
                                        size = 0.15, 
                                        linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, 
                                        linetype = 'solid',
                                        colour = "#ebe4e4"),
        legend.position = "bottom", 
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8.5))+
  scale_y_continuous(breaks = number_ticks(4), 
                     labels=dollar_format(prefix="$")) +
  #scale_x_continuous(breaks = number_ticks(4), 
  #                   labels=dollar_format(prefix="$")) +
  scale_color_manual(values=c ("#575c59", "#2cc9de", "#eb542f")) +
  labs(title = "One-way Sensitivity Analysis",
       x = "Parameter Values",
       y = "Net Monetary Benefit (in $1,000 Mexican Pesos)",
       color = "Strategy")

GG_tw <- plot(tw_hosp) + 
  theme(plot.title = element_text(face = "bold", 
                                  size = 12,
                                  family =, 
                                  hjust = 0.5),
        axis.line.x = element_line(color="#ebe4e4", size = 1),
        axis.line.y = element_line(color="#ebe4e4", size = 1),
        plot.caption = element_text(hjust = 0,
                                    colour = "#777777",
                                    size = 10),
        panel.background = element_rect(fill = "white", 
                                        colour = "white", 
                                        size = 0.15, 
                                        linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, 
                                        linetype = 'solid',
                                        colour = "#ebe4e4"),
        legend.position = "bottom", 
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8.5)) +
  labs(title = "Two-way Sensitivity Analysis",
       fill = "Strategy") +
  scale_fill_manual(values=c ("#575c59", "#2cc9de", "#eb542f"))

ggarrange(GG_ow, GG_tw, ncol = 1)

ggsave(paste0("figs/Sensitivity Analysis",
              format(Sys.Date(), "%F"), ".png"), 
       width = 7, height = 10)
