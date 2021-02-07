install.packages("Stats")
library(stats)

rs.surv.exp <- function (formula = formula(data), data = parent.frame(), ratetable = relsurv::slopop, 
          na.action, fin.date, method = "pohar-perme", conf.type = "log", 
          conf.int = 0.95, type = "kaplan-meier", add.times, precision = 1, 
          rmap) 
{
  call <- match.call()
  if (!missing(rmap)) {
    rmap <- substitute(rmap)
  }
  rform <- rformulate(formula, data, ratetable, na.action, 
                      rmap)
  data <- rform$data
  type <- match.arg(type, c("kaplan-meier", "fleming-harrington"))
  type <- match(type, c("kaplan-meier", "fleming-harrington"))
  method <- match.arg(method, c("pohar-perme", "ederer2", 
                                "hakulinen", "ederer1"))
  method <- match(method, c("pohar-perme", "ederer2", "hakulinen", 
                            "ederer1"))
  conf.type <- match.arg(conf.type, c("plain", "log", "log-log"))
  if (method == 3) {
    R <- rform$R
    coll <- match("year", attributes(ratetable)$dimid)
    year <- R[, coll]
    if (missing(fin.date)) 
      fin.date <- max(rform$Y + year)
    Y2 <- rform$Y
    if (length(fin.date) == 1) 
      Y2[rform$status == 1] <- fin.date - year[rform$status == 
                                                 1]
    else if (length(fin.date) == nrow(rform$R)) 
      Y2[rform$status == 1] <- fin.date[rform$status == 
                                          1] - year[rform$status == 1]
    else stop("fin.date must be either one value or a vector of the same length as the data")
    status2 <- rep(0, nrow(rform$X))
  }
  p <- rform$m
  if (p > 0) 
    data$Xs <- strata(rform$X[, , drop = FALSE])
  else data$Xs <- rep(1, nrow(data))
  se.fac <- sqrt(qchisq(conf.int, 1))
  out <- NULL
  out$n <- table(data$Xs)
  out$time <- out$n.risk <- out$n.event <- out$n.censor <- out$surv <- out$std.err <- out$strata <- NULL
  for (kt in 1:length(out$n)) {
    inx <- which(data$Xs == names(out$n)[kt])
    tis <- sort(unique(rform$Y[inx]))
    if (method == 1 & !missing(add.times)) {
      add.times <- pmin(as.numeric(add.times), max(rform$Y[inx]))
      tis <- sort(union(rform$Y[inx], as.numeric(add.times)))
    }
    if (method == 3) 
      tis <- sort(unique(pmin(max(tis), c(tis, Y2[inx]))))
    temp <- exp.prep(rform$R[inx, , drop = FALSE], rform$Y[inx], 
                     rform$ratetable, rform$status[inx], times = tis, 
                     fast = (method < 3), prec = precision)
    out$time <- c(out$time, tis)
    out$n.risk <- c(out$n.risk, temp$yi)
    out$n.event <- c(out$n.event, temp$dni)
    out$n.censor <- c(out$n.censor, c(-diff(temp$yi), temp$yi[length(temp$yi)]) - 
                        temp$dni)
    if (method == 1) {
      approximate <- temp$yidlisiw
      haz <- temp$dnisi/temp$yisi - approximate
      out$hazard_exc <- haz
      out$hazard_obs <- temp$dnisi/temp$yisi
      out$hazard_pop <- approximate 
      out$std.err <- c(out$std.err, sqrt(cumsum(temp$dnisisq/(temp$yisi)^2)))
    }
    else if (method == 2) {
      haz <- temp$dni/temp$yi - temp$yidli/temp$yi
      out$hazard_exc <- haz
      out$hazard_obs <- temp$dni/temp$yi
      out$hazard_pop <- temp$yidli/temp$yi 
      out$std.err <- c(out$std.err, sqrt(cumsum(temp$dni/(temp$yi)^2)))
    }
    else if (method == 3) {
      temp2 <- exp.prep(rform$R[inx, , drop = FALSE], 
                        Y2[inx], ratetable, status2[inx], times = tis)
      popsur <- exp(-cumsum(temp2$yisidli/temp2$yisis))
      haz <- temp$dni/temp$yi
      out$hazard_exc <- haz
      out$hazard_obs <- haz
      out$hazard_pop <- haz 
      out$std.err <- c(out$std.err, sqrt(cumsum(temp$dni/(temp$yi)^2)))
    }
    else if (method == 4) {
      popsur <- temp$sis/length(inx)
      haz <- temp$dni/temp$yi
      out$hazard_exc <- haz
      out$hazard_obs <- haz
      out$hazard_pop <- haz 
      out$std.err <- c(out$std.err, sqrt(cumsum(temp$dni/(temp$yi)^2)))
    }
    if (type == 2) 
      survtemp <- exp(-cumsum(haz))
    else survtemp <- cumprod(1 - haz)
    if (method > 2) {
      survtemp <- survtemp/popsur
    }
    out$surv <- c(out$surv, survtemp)
    out$strata <- c(out$strata, length(tis))
  }
  if (conf.type == "plain") {
    out$lower <- as.vector(out$surv - out$std.err * se.fac * 
                             out$surv)
    out$upper <- as.vector(out$surv + out$std.err * se.fac * 
                             out$surv)
  }
  else if (conf.type == "log") {
    out$lower <- exp(as.vector(log(out$surv) - out$std.err * 
                                 se.fac))
    out$upper <- exp(as.vector(log(out$surv) + out$std.err * 
                                 se.fac))
  }
  else if (conf.type == "log-log") {
    out$lower <- exp(-exp(as.vector(log(-log(out$surv)) - 
                                      out$std.err * se.fac/log(out$surv))))
    out$upper <- exp(-exp(as.vector(log(-log(out$surv)) + 
                                      out$std.err * se.fac/log(out$surv))))
  }
  names(out$strata) <- names(out$n)
  if (p == 0) {
    out$strata <- NULL
  }
  out$n <- as.vector(out$n)
  out$conf.type <- conf.type
  out$conf.int <- conf.int
  out$method <- method
  out$call <- call
  out$type <- "right"
  class(out) <- c("survfit", "rs.surv")
  out
}


KM_surv_exp <- rs.surv.exp(Surv(time, state) ~ 1,
                   data = Covid_rr,
                   ratetable = rate_exp_mx_2020, 
                   type="kaplan-meier",
                   conf.type="log",
                   conf.int=0.95,
                   rmap = list(age = age, sex = sex, year = diag))

Covid_19 <- head(KM_surv_exp$hazard_exc, 61)
Observed <- head(KM_surv_exp$hazard_obs, 61)
Population <- head(KM_surv_exp$hazard_pop, 61)
time <- seq(0,60, by = 1)

df_hazard_exp <- data.frame(time, Covid_19, Observed, Population)

df_hazard_exp_long <- gather(data = df_hazard_exp, 
                           key = "Survival", 
                           value = "Hazard",
                           -time)


ggplot(df_hazard_exp_long, 
       aes(x = time, 
           y = Hazard,
           color = Survival)) +
  geom_line() 
  
#### Hazard by groups ####
KM_surv_sex <- rs.surv(Surv(time, state) ~ sex,
                           data = Covid_rr,
                           ratetable = rate_exp_mx_2020, 
                           type="kaplan-meier",
                           conf.type="log",
                           conf.int=0.95,
                           rmap = list(age = age, sex = sex, year = diag))

Covid_rr_muj <- Covid_rr %>% 
  filter(sex == 2)

dim(Covid_rr_muj)

KM_surv_sex_age <- rs.surv(Surv(time, state) ~ sex + age_range,
                       data = Covid_rr,
                       ratetable = rate_exp_mx_2020, 
                       type="kaplan-meier",
                       conf.type="log",
                       conf.int=0.95,
                       rmap = list(age = age, sex = sex, year = diag))
#### No salió bien, haré 8 data frames de Covid_rr

Covid_rr_male_45 <- Covid_rr %>% 
  filter(sex == 1 & age_range == "45 - 54")

dim(Covid_rr_male_45)

Covid_rr_male_55 <- Covid_rr %>% 
  filter(sex == 1 & age_range == "55 - 64")

dim(Covid_rr_male_55)

Covid_rr_male_70 <- Covid_rr %>% 
  filter(sex == 1 & age_range == "70 +")

dim(Covid_rr_male_70)

Covid_rr_male_65 <- Covid_rr %>% 
  filter(sex == 1 & age_range == "65 - 69")

dim(Covid_rr_male_65)

### Mujeres ####
Covid_rr_female_45 <- Covid_rr %>% 
  filter(sex == 2 & age_range == "45 - 54")

dim(Covid_rr_female_45)

Covid_rr_female_55 <- Covid_rr %>% 
  filter(sex == 2 & age_range == "55 - 64")

dim(Covid_rr_female_55)

Covid_rr_female_70 <- Covid_rr %>% 
  filter(sex == 2 & age_range == "70 +")

dim(Covid_rr_female_70)

Covid_rr_female_65 <- Covid_rr %>% 
  filter(sex == 2 & age_range == "65 - 69")

dim(Covid_rr_female_65)


KM_surv_sex_male_45 <- rs.surv.exp(Surv(time, state) ~ 1,
                           data = Covid_rr_male_45,
                           ratetable = rate_exp_mx_2020, 
                           type="kaplan-meier",
                           conf.type="log",
                           conf.int=0.95,
                           rmap = list(age = age, sex = sex, year = diag))

KM_surv_sex_male_55 <- rs.surv.exp(Surv(time, state) ~ 1,
                                   data = Covid_rr_male_55,
                                   ratetable = rate_exp_mx_2020, 
                                   type="kaplan-meier",
                                   conf.type="log",
                                   conf.int=0.95,
                                   rmap = list(age = age, sex = sex, year = diag))

KM_surv_sex_male_65 <- rs.surv.exp(Surv(time, state) ~ 1,
                                   data = Covid_rr_male_65,
                                   ratetable = rate_exp_mx_2020, 
                                   type="kaplan-meier",
                                   conf.type="log",
                                   conf.int=0.95,
                                   rmap = list(age = age, sex = sex, year = diag))

KM_surv_sex_male_70 <- rs.surv.exp(Surv(time, state) ~ 1,
                                   data = Covid_rr_male_70,
                                   ratetable = rate_exp_mx_2020, 
                                   type="kaplan-meier",
                                   conf.type="log",
                                   conf.int=0.95,
                                   rmap = list(age = age, sex = sex, year = diag))

#### Female ####
KM_surv_sex_female_45 <- rs.surv.exp(Surv(time, state) ~ 1,
                                   data = Covid_rr_female_45,
                                   ratetable = rate_exp_mx_2020, 
                                   type="kaplan-meier",
                                   conf.type="log",
                                   conf.int=0.95,
                                   rmap = list(age = age, sex = sex, year = diag))

KM_surv_sex_female_55 <- rs.surv.exp(Surv(time, state) ~ 1,
                                   data = Covid_rr_female_55,
                                   ratetable = rate_exp_mx_2020, 
                                   type="kaplan-meier",
                                   conf.type="log",
                                   conf.int=0.95,
                                   rmap = list(age = age, sex = sex, year = diag))

KM_surv_sex_female_65 <- rs.surv.exp(Surv(time, state) ~ 1,
                                   data = Covid_rr_female_65,
                                   ratetable = rate_exp_mx_2020, 
                                   type="kaplan-meier",
                                   conf.type="log",
                                   conf.int=0.95,
                                   rmap = list(age = age, sex = sex, year = diag))

KM_surv_sex_female_70 <- rs.surv.exp(Surv(time, state) ~ 1,
                                   data = Covid_rr_female_70,
                                   ratetable = rate_exp_mx_2020, 
                                   type="kaplan-meier",
                                   conf.type="log",
                                   conf.int=0.95,
                                   rmap = list(age = age, sex = sex, year = diag))

#### Data Frame with hazards by each group ####

# Males

# 45
Covid_19 <- head(KM_surv_sex_male_45$hazard_exc, 60)
Observed <- head(KM_surv_sex_male_45$hazard_obs, 60)
Population <- head(KM_surv_sex_male_45$hazard_pop, 60)
time <- seq(1,60, by = 1)

df_hazard_exp_male_45 <- data.frame(time, 
                                    Covid_19, 
                                    Observed, 
                                    Population)

df_hazard_exp_long_45 <- gather(data = df_hazard_exp_male_45, 
                             key = "Type", 
                             value = "Hazard",
                             -time)

df_hazard_exp_long_45 <- df_hazard_exp_long_45 %>% 
  mutate(group = "45 - 54")

# 55
Covid_19 <- head(KM_surv_sex_male_55$hazard_exc, 60)
Observed <- head(KM_surv_sex_male_55$hazard_obs, 60)
Population <- head(KM_surv_sex_male_55$hazard_pop, 60)
time <- seq(1,60, by = 1)

df_hazard_exp_male_55 <- data.frame(time, 
                                    Covid_19, 
                                    Observed, 
                                    Population)

df_hazard_exp_long_55 <- gather(data = df_hazard_exp_male_55, 
                                key = "Type", 
                                value = "Hazard",
                                -time)

df_hazard_exp_long_55 <- df_hazard_exp_long_55 %>% 
  mutate(group = "55 - 64")

# 65
Covid_19 <- head(KM_surv_sex_male_65$hazard_exc, 60)
Observed <- head(KM_surv_sex_male_65$hazard_obs, 60)
Population <- head(KM_surv_sex_male_65$hazard_pop, 60)
time <- seq(1,60, by = 1)

df_hazard_exp_male_65 <- data.frame(time, 
                                    Covid_19, 
                                    Observed, 
                                    Population)

df_hazard_exp_long_65 <- gather(data = df_hazard_exp_male_65, 
                                key = "Type", 
                                value = "Hazard",
                                -time)

df_hazard_exp_long_65 <- df_hazard_exp_long_65 %>% 
  mutate(group = "65 - 69")

# 70
Covid_19 <- head(KM_surv_sex_male_70$hazard_exc, 60)
Observed <- head(KM_surv_sex_male_70$hazard_obs, 60)
Population <- head(KM_surv_sex_male_70$hazard_pop, 60)
time <- seq(1,60, by = 1)

df_hazard_exp_male_70 <- data.frame(time, 
                                    Covid_19, 
                                    Observed, 
                                    Population)

df_hazard_exp_long_70 <- gather(data = df_hazard_exp_male_70, 
                                key = "Type", 
                                value = "Hazard",
                                -time)

df_hazard_exp_long_70 <- df_hazard_exp_long_70 %>% 
  mutate(group = "70 +")

df_hazard_long_male <- rbind(df_hazard_exp_long_45, df_hazard_exp_long_55,
                             df_hazard_exp_long_65, df_hazard_exp_long_70)

df_hazard_long_male <- df_hazard_long_male %>% 
  mutate(sex = "male")

gg_males_hazard <- ggplot(df_hazard_long_male, 
       aes(x = time, 
           y = Hazard,
           color = Type)) +
  geom_line(size = 1)+
  facet_wrap(~group, scales = "free") +
  theme(plot.title = element_text(face = "bold", 
                                  size = 16,
                                  family =),
        plot.caption = element_text(hjust = 0,
                                    colour = "#777777",
                                    size = 10),
        panel.background = element_rect(fill = "white", 
                                        colour = "gray", 
                                        size = 0.15, 
                                        linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, 
                                        linetype = 'solid',
                                        colour = "gray"))+ 
       # axis.text.x = element_text(angle = 90, hjust = 0))+
    scale_x_continuous(breaks = number_ticks(4))+
    scale_y_continuous(breaks = number_ticks(3))+
    scale_color_manual(values=c ("#2e7bff", "#e6ad02", "#02ebdb")) +
    labs(title = "Male age groups",
         x = " Days ",
         y = " ",
         color = "Hazard")


# Females

# 45
Covid_19 <- head(KM_surv_sex_female_45$hazard_exc, 60)
Observed <- head(KM_surv_sex_female_45$hazard_obs, 60)
Population <- head(KM_surv_sex_female_45$hazard_pop, 60)
time <- seq(1,60, by = 1)

df_hazard_exp_female_45 <- data.frame(time, 
                                    Covid_19, 
                                    Observed, 
                                    Population)

df_hazard_exp_long_45 <- gather(data = df_hazard_exp_female_45, 
                                key = "Type", 
                                value = "Hazard",
                                -time)

df_hazard_exp_long_45 <- df_hazard_exp_long_45 %>% 
  mutate(group = "45 - 54")

# 55
Covid_19 <- head(KM_surv_sex_female_55$hazard_exc, 60)
Observed <- head(KM_surv_sex_female_55$hazard_obs, 60)
Population <- head(KM_surv_sex_female_55$hazard_pop, 60)
time <- seq(1,60, by = 1)

df_hazard_exp_female_55 <- data.frame(time, 
                                    Covid_19, 
                                    Observed, 
                                    Population)

df_hazard_exp_long_55 <- gather(data = df_hazard_exp_female_55, 
                                key = "Type", 
                                value = "Hazard",
                                -time)

df_hazard_exp_long_55 <- df_hazard_exp_long_55 %>% 
  mutate(group = "55 - 64")

# 65
Covid_19 <- head(KM_surv_sex_female_65$hazard_exc, 60)
Observed <- head(KM_surv_sex_female_65$hazard_obs, 60)
Population <- head(KM_surv_sex_female_65$hazard_pop, 60)
time <- seq(1,60, by = 1)

df_hazard_exp_female_65 <- data.frame(time, 
                                    Covid_19, 
                                    Observed, 
                                    Population)

df_hazard_exp_long_65 <- gather(data = df_hazard_exp_female_65, 
                                key = "Type", 
                                value = "Hazard",
                                -time)

df_hazard_exp_long_65 <- df_hazard_exp_long_65 %>% 
  mutate(group = "65 - 69")

# 70
Covid_19 <- head(KM_surv_sex_female_70$hazard_exc, 60)
Observed <- head(KM_surv_sex_female_70$hazard_obs, 60)
Population <- head(KM_surv_sex_female_70$hazard_pop, 60)
time <- seq(1,60, by = 1)

df_hazard_exp_female_70 <- data.frame(time, 
                                    Covid_19, 
                                    Observed, 
                                    Population)

df_hazard_exp_long_70 <- gather(data = df_hazard_exp_female_70, 
                                key = "Type", 
                                value = "Hazard",
                                -time)

df_hazard_exp_long_70 <- df_hazard_exp_long_70 %>% 
  mutate(group = "70 +")

df_hazard_long_female <- rbind(df_hazard_exp_long_45, df_hazard_exp_long_55,
                             df_hazard_exp_long_65, df_hazard_exp_long_70)

df_hazard_long_female <- df_hazard_long_female %>% 
  mutate(sex = "female")

gg_females_hazard <- ggplot(df_hazard_long_female, 
                          aes(x = time, 
                              y = Hazard,
                              color = Type)) +
  geom_line(size = 1)+
  facet_wrap(~group, scales = "free") +
  theme(plot.title = element_text(face = "bold", 
                                  size = 16,
                                  family =),
        plot.caption = element_text(hjust = 0,
                                    colour = "#777777",
                                    size = 10),
        panel.background = element_rect(fill = "white", 
                                        colour = "gray", 
                                        size = 0.15, 
                                        linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, 
                                        linetype = 'solid',
                                        colour = "gray"))+
  scale_x_continuous(breaks = number_ticks(4))+
  scale_y_continuous(breaks = number_ticks(3))+
  # axis.text.x = element_text(angle = 90, hjust = 0))+
  scale_color_manual(values=c ("#2e7bff", "#e6ad02", "#02ebdb")) +
  labs(title = "Female age groups",
       x = " ",
       y = " ",
       color = "Hazard")

ggarrange(gg_females_hazard, gg_males_hazard,
          #labels        = c("Females", "Males"),
          ncol          = 1, 
          common.legend = TRUE,
          label.y = c(1,0),
          legend        = "bottom")

ggsave(paste0("figs/Hazards",
              format(Sys.Date(), "%F"), ".pdf"), 
       width = 7, height = 5)


df_hazards <- rbind(df_hazard_long_female, df_hazard_long_male)

save(df_hazards, file = "data/df_hazards.Rdata")

#### Rsadd experiment - Not good ####

fit <- rsadd(Surv(time,state)~sex,
             ratetable=rate_exp_mx_2020,
             data=Covid_rr,
             rmap = list(age = age, sex = sex, year = diag))

d <- rs.surv.rsadd(fit,
                   newdata=data.frame(sex=1))

plot(d)

             