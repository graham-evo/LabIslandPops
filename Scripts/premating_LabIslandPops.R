# Graham C. McLaughlin
# Created: 09-01-2023
# Modified: 09-01-2023
# Analysis of August 26th experiment characterizing RI between our four Lab Island Pops
# ----

# Libraries:
library(lme4)
library(emmeans)
library(ggplot2)
library(RColorBrewer)
library(afex)
library(survival)
library(survminer)
library(olsrr)
library(rstatix)
library(dplyr)
library(ggpubr)

## Data Exploration/Management:

ncsi <- read.table("Data/LabIslandPops_premating_R.txt",
                   header = TRUE,
                   comment.char = "",
                   na.strings = "#N/A")
head(ncsi)

# Let's first define our factors:
ncsi$female <- factor(ncsi$female,
                      levels = c("iv","lhm","dhm","val"),
                      labels = c("IV", "LHm", "Dhm", "Val"))
ncsi$male <- factor(ncsi$male,
                    levels = c("iv","lhm","dhm","val"),
                    labels = c("IV","LHm","Dhm","Val"))
head(ncsi) # Factors defined and relabelled

# Let's add a cross column to make comparisons easier:
ncsi$cross <- paste0(ncsi$female, "_", ncsi$male)
head(ncsi) # We now have a cross column

# Trimming out the vials lost due to experimental failure that won't be used for analysis of any measures:
dim(ncsi) # 400 observations before trimming
ncsi <- subset(ncsi, failure == 0)
dim(ncsi) # Lost 52 to experimental failure (most of those were due to lack of virgins)

## Looking first at Willingness To Mate: -----
# Starting with just the simplest model assessing the effect of cross on probability of mating
# We will add in the effect of male and female and any interaction terms later.
hist(ncsi$mated)

num_mated <- ncsi %>%
  group_by(cross) %>%
  summarise(count = sum(mated))

glm_willingness <- glm(mated ~ male*female,
                       family = binomial,
                       data = ncsi)

glm_willingness
summary(glm_willingness)
anova(glm_willingness)
emmeans(glm_willingness, ~cross, type = "response", level = 0.69)
pairs(emmeans(glm_willingness, ~cross, type = "response"), adjust = "none")
# We aren't interested in all of the above pairwise comparisons so this is not too helpful.

plot_frame <- as.data.frame(emmeans(glm_willingness, ~cross, type = "response", level = 0.69))

ggplot(data = plot_frame, aes(x = cross, y = prob, color = male)) +
  geom_point(shape = 21) +
  geom_errorbar(aes(x = cross, ymin = prob-SE, ymax = prob + SE)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Heat Map of Willingness to Mate:

plot_frame %>%
  arrange(male) %>%
  mutate(female = factor(female, levels=c("Val","Dhm","LHm","IV"))) %>%
  ggplot(aes(x = male, y = female)) +
  geom_raster(aes(fill = prob)) +
  scale_fill_gradient(low = "grey100", high = "red") +
  labs(x = "Male Genotype", y = "Female Genotype", fill = "Probability\nof\nMating") +
  theme_bw() +
  scale_x_discrete(position = "top")

# It's easier to visualize and model if we only look at plots of the within-pop cross versus a single across population cross. We can then run a model with male/female and interaction effects:

## Starting with Dahomey versus IV:
IV_Dhm <- subset(ncsi, cross %in% c("IV_IV", "Dhm_Dhm", "IV_Dhm", "Dhm_IV"))

glm_IV_Dhm_willingness <- glm(mated ~ cross,
                              family = binomial,
                              data = IV_Dhm)
summary(glm_IV_Dhm_willingness)
emmeans(glm_IV_Dhm_willingness, ~cross, type = "response", level = 0.69)
pairs(emmeans(glm_IV_Dhm_willingness, ~cross, type = "response"), adjust = "none")

plot_frame <- as.data.frame(emmeans(glm_IV_Dhm_willingness,
                                    ~cross,
                                    type = "response"))

ggplot(aes(x = cross, y = prob), data = plot_frame) +
  geom_point() +
  geom_errorbar(aes(x = cross, ymin = prob-SE, ymax = prob + SE))

# Adding male/female and interaction term to IV versus Dahomey comparison:
glm_IV_Dhm_willingness_interaction <- glm(mated ~ male*female,
                              family = binomial,
                              data = IV_Dhm)
summary(glm_IV_Dhm_willingness_interaction)

## Dahomey versus LHm-GFP:
Dhm_LHm <- subset(ncsi, cross %in% c("LHm-GFP_LHm-GFP", "Dhm_Dhm", "LHm-GFP_Dhm", "Dhm_LHm-GFP"))

glm_Dhm_LHm_willingness <- glm(mated ~ cross,
                              family = binomial,
                              data = Dhm_LHm)
summary(glm_Dhm_LHm_willingness)
emmeans(glm_Dhm_LHm_willingness, ~cross, type = "response")

plot_frame <- as.data.frame(emmeans(glm_Dhm_LHm_willingness,
                                    ~cross,
                                    type = "response"))

ggplot(aes(x = cross, y = prob), data = plot_frame) +
  geom_point() +
  geom_errorbar(aes(x = cross, ymin = prob-SE, ymax = prob + SE))

# Adding male/female and interaction term to IV versus Dahomey comparison:
glm_Dhm_LHm_willingness_interaction <- glm(mated ~ male*female,
                                          family = binomial,
                                          data = Dhm_LHm)
summary(glm_Dhm_LHm_willingness_interaction)

## Dahomey versus Valais:
Dhm_Val <- subset(ncsi, cross %in% c("Val_Val", "Dhm_Dhm", "Val_Dhm", "Dhm_Val"))

glm_Val_Dhm_willingness <- glm(mated ~ cross,
                              family = binomial,
                              data = Dhm_Val)
summary(glm_Val_Dhm_willingness)
emmeans(glm_Val_Dhm_willingness, ~cross, type = "response")
pairs(emmeans(glm_Val_Dhm_willingness, ~cross, type = "response"), adjust = "none")

plot_frame <- as.data.frame(emmeans(glm_Val_Dhm_willingness,
                                    ~cross,
                                    type = "response"))

ggplot(aes(x = cross, y = prob), data = plot_frame) +
  geom_point() +
  geom_errorbar(aes(x = cross, ymin = prob-SE, ymax = prob + SE))

# Adding male/female and interaction term to IV versus Dahomey comparison:
glm_IV_Dhm_willingness_interaction <- glm(mated ~ male*female,
                                          family = binomial,
                                          data = Dhm_Val)
summary(glm_IV_Dhm_willingness_interaction)

## IV Versus LHm-GFP:
iv_lhm <- subset(ncsi, cross %in% c("IV_IV", "LHm-GFP_LHm-GFP", "IV_LHm-GFP", "LHm-GFP_IV"))

glm_iv_lhm_willingness <- glm(mated ~ cross,
                               family = binomial,
                               data = iv_lhm)
summary(glm_iv_lhm_willingness)
emmeans(glm_iv_lhm_willingness, ~cross, type = "response")
pairs(emmeans(glm_iv_lhm_willingness, ~cross, type = "response"), adjust = "none")

plot_frame <- as.data.frame(emmeans(glm_iv_lhm_willingness,
                                    ~cross,
                                    type = "response"))

ggplot(aes(x = cross, y = prob), data = plot_frame) +
  geom_point() +
  geom_errorbar(aes(x = cross, ymin = prob-SE, ymax = prob + SE))

# Adding male/female and interaction term to IV versus Dahomey comparison:
glm_iv_lhm_willingness_interaction <- glm(mated ~ male*female,
                                          family = binomial,
                                          data = iv_lhm)
summary(glm_iv_lhm_willingness_interaction)

## IV versus Valais:
iv_val <- subset(ncsi, cross %in% c("IV_IV", "Val_Val", "IV_Val", "Val_IV"))

glm_iv_val_willingness <- glm(mated ~ cross,
                              family = binomial,
                              data = iv_val)
summary(glm_iv_val_willingness)
emmeans(glm_iv_val_willingness, ~cross, type = "response")
pairs(emmeans(glm_iv_val_willingness, ~cross, type = "response"), adjust = "none")

plot_frame <- as.data.frame(emmeans(glm_iv_val_willingness,
                                    ~cross,
                                    type = "response"))

ggplot(aes(x = cross, y = prob), data = plot_frame) +
  geom_point() +
  geom_errorbar(aes(x = cross, ymin = prob-SE, ymax = prob + SE))

# Adding male/female and interaction term to IV versus Dahomey comparison:
glm_iv_val_willingness_interaction <- glm(mated ~ male*female,
                                          family = binomial,
                                          data = iv_val)
summary(glm_iv_val_willingness_interaction)

## Val versus Dahomey:
val_dhm <- subset(ncsi, cross %in% c("Dhm_Dhm", "Val_Val", "Dhm_Val", "Val_Dhm"))

glm_dhm_val_willingness <- glm(mated ~ cross,
                              family = binomial,
                              data = val_dhm)
summary(glm_dhm_val_willingness)
emmeans(glm_dhm_val_willingness, ~cross, type = "response")
pairs(emmeans(glm_dhm_val_willingness, ~cross, type = "response"))

plot_frame <- as.data.frame(emmeans(glm_dhm_val_willingness,
                                    ~cross,
                                    type = "response"))

ggplot(aes(x = cross, y = prob), data = plot_frame) +
  geom_point() +
  geom_errorbar(aes(x = cross, ymin = prob-SE, ymax = prob + SE))

# Adding male/female and interaction term to IV versus Dahomey comparison:
glm_dhm_val_willingness_interaction <- glm(mated ~ male*female,
                                          family = binomial,
                                          data = val_dhm)
summary(glm_dhm_val_willingness_interaction)

## Val versus LHm-GFP:
val_LHm <- subset(ncsi, cross %in% c("Val_Val", "LHm-GFP_LHm-GFP", "Val_LHm-GFP", "LHm-GFP_Val"))

glm_val_lhm_willingness <- glm(mated ~ cross,
                               family = binomial,
                               data = val_LHm)
summary(glm_val_lhm_willingness)
emmeans(glm_val_lhm_willingness, ~cross, type = "response")
pairs(emmeans(glm_val_lhm_willingness, ~cross, type = "response"), adjust = "none")

plot_frame <- as.data.frame(emmeans(glm_val_lhm_willingness,
                                    ~cross,
                                    type = "response"))

ggplot(aes(x = cross, y = prob), data = plot_frame) +
  geom_point() +
  geom_errorbar(aes(x = cross, ymin = prob-SE, ymax = prob + SE))

# Adding male/female and interaction term to IV versus Dahomey comparison:
glm_lhm_val_willingness_interaction <- glm(mated ~ male*female,
                                           family = binomial,
                                           data = val_LHm)
summary(glm_lhm_val_willingness_interaction)

# Plotting Joint Isolation Indices of each genetic background comparison (Merrell 1950)
i_index <- c(1.2333, 0.88095, 0.40625, 0.94117, 1.090, 1.03125)
cross <- c("IV\n&\nDhm","IV\n&\nVal","IV\n&\nLHm","LHm\n&\nVal","LHm\n&\nDhm","Val\n&\nDhm")

si <- data.frame(i_index, cross)

ggplot(data = si, aes(x = cross, y = i_index)) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_hline(yintercept = 1, linetype = 2) +
  scale_y_continuous(name = "Joint Sexual Isolation Index (I)", limits = c(0,1.25), breaks = seq(0,1.25,0.2)) +
  scale_x_discrete(name = "Population Pair") +
  theme_classic()

## Moving on to Latency To Mate:
head(ncsi, n = 100)
ncsi$failure == TRUE # Double checking that all failures are removed:

# Setting all did not mate (dnm) vials to 180:
ncsi$latency[ncsi$dnm==1] <- 180
head(ncsi, n =100)
ncsi$latency < 0
View(ncsi) #Vials 113 and 248 have negative latencies even though they are coded
# as having mated, so lets just remove them for now
ncsi <- subset(ncsi, id != c(113, 248))
ncsi$latency<0 # negative latency observations removed

# Running a model assessing differences in latency based on cross:
coxph_latency <- coxph(formula = Surv(latency, mated) ~ male*female, data = ncsi)
coxph_latency
summary(coxph_latency)
anova(coxph_latency)
emmeans(coxph_latency, ~cross, type = "response")
pairs(emmeans(coxph_latency, ~cross, type = "response"), adjust = "none")

# Plot of all cross latencies:
ggplot(data = subset(ncsi, latency != 180), aes(x = cross, y = latency, fill = female)) +
  geom_boxplot()

# HeatMap of all median latencies by male and female:
library(gplots)
library(reshape2)
library(cowplot)

test <- as.data.frame(aggregate(latency~male+female, data = ncsi, FUN = median))

test %>%
  arrange(male) %>%
  mutate(female = factor(female, levels=c("Val","Dhm","LHm","IV"))) %>%
  ggplot(aes(x = male, y = female)) +
  geom_raster(aes(fill = latency)) +
  scale_fill_gradient(low = "grey100", high = "red") +
  labs(x = "Male Genotype", y = "Female Genotype", fill = "Latency\n(min)") +
  theme_bw() +
  scale_x_discrete(position = "top")



# dhm versus iv latency
IV_Dhm <- subset(ncsi, cross %in% c("IV_IV", "Dhm_Dhm", "IV_Dhm", "Dhm_IV"))

coxph_latency <- coxph(formula = Surv(latency, mated) ~ cross + male + female, data = IV_Dhm)
coxph_latency

ggplot(data = IV_Dhm, aes(x = cross, y = latency)) +
  geom_boxplot()

# dahomey versus lhm latency
Dhm_LHm <- subset(ncsi, cross %in% c("LHm-GFP_LHm-GFP", "Dhm_Dhm", "LHm-GFP_Dhm", "Dhm_LHm-GFP"))

coxph_latency <- coxph(formula = Surv(latency, mated) ~ cross + male + female, data = Dhm_LHm)
coxph_latency

ggplot(data = subset(Dhm_LHm), aes(x = cross, y = latency)) +
  geom_boxplot()

# dahomey versus val
Dhm_Val <- subset(ncsi, cross %in% c("Val_Val", "Dhm_Dhm", "Val_Dhm", "Dhm_Val"))

coxph_latency <- coxph(formula = Surv(latency, mated) ~ cross, data = Dhm_Val)
coxph_latency
pairs(emmeans(coxph_latency, ~cross, type = "response"), adjust = "none")

ggplot(data =Dhm_Val, aes(x = cross, y = latency)) +
  geom_boxplot()

# iv versus lhm
iv_lhm <- subset(ncsi, cross %in% c("IV_IV", "LHm-GFP_LHm-GFP", "IV_LHm-GFP", "LHm-GFP_IV"))

coxph_latency <- coxph(formula = Surv(latency, mated) ~ cross, data = iv_lhm)
coxph_latency
pairs(emmeans(coxph_latency, ~cross, type = "response"), adjust = "none")

ggplot(data = iv_lhm, aes(x = cross, y = latency)) +
  geom_boxplot()

# iv versus val
iv_val <- subset(ncsi, cross %in% c("IV_IV", "Val_Val", "IV_Val", "Val_IV"))

coxph_latency <- coxph(formula = Surv(latency, mated) ~ cross, data = iv_val)
coxph_latency
pairs(emmeans(coxph_latency, ~cross, type = "response"), adjust = "none")

ggplot(data = iv_val, aes(x = cross, y = latency)) +
  geom_boxplot()

# val versus lhm
val_lhm <- subset(ncsi, cross %in% c("Val_LHm-GFP", "Val_Val", "LHm-GFP_Val", "LHm-GFP_LHm-GFP"))

coxph_latency <- coxph(formula = Surv(latency, mated) ~ cross, data = val_lhm)
coxph_latency
pairs(emmeans(coxph_latency, ~cross, type = "response"), adjust = "none")

ggplot(data = val_lhm, aes(x = cross, y = latency)) +
  geom_boxplot()

# Plotting each within-population cross against only meaningful/relevant cross comparisons:
ggplot(data = subset(ncsi, male == "Dhm" | female == "Dhm"), aes(x = cross, y = latency)) +
  geom_boxplot()

ggplot(data = subset(ncsi, male == "IV" | female == "IV"), aes(x = cross, y = latency)) +
  geom_boxplot()

ggplot(data = subset(ncsi, male == "Val" | female == "Val"), aes(x = cross, y = latency)) +
  geom_boxplot()

ggplot(data = subset(ncsi, male == "LHm-GFP" | female == "LHm-GFP"), aes(x = cross, y = latency)) +
  geom_boxplot()

# Mating Duration:

# Removing the vials for which they didn't mate:
ncsi <- subset(ncsi, dnm == 0)
ncsi$mating_duration > 0
View(ncsi)
# Vial 384 had a start time but not stop time so let's remove it:
ncsi <- subset(ncsi, id != 384)

# Running the model:

lm_mating_duration <- lm(mating_duration ~ male*female, data = ncsi)
summary(lm_mating_duration)
anova(lm_mating_duration)

ggplot(data = ncsi, aes(x = cross, y = mating_duration, fill = male)) +
  geom_boxplot() +
  theme_classic()

hist(ncsi$mating_duration, breaks = 100)

# individual comparisons:
# Plot of all cross mating durations:
test <- as.data.frame(aggregate(mating_duration~male+female, data = ncsi, FUN = mean))

test %>%
  arrange(male) %>%
  mutate(female = factor(female, levels=c("Val","Dhm","LHm","IV"))) %>%
  ggplot(aes(x = male, y = female)) +
  geom_raster(aes(fill = mating_duration)) +
  scale_fill_gradient(low = "grey100", high = "red") +
  labs(x = "Male Genotype", y = "Female Genotype", fill = "Mating\nDuration\n(min)") +
  theme_bw() +
  scale_x_discrete(position = "top")

# dhm versus iv
IV_Dhm <- subset(ncsi, cross %in% c("IV_IV", "Dhm_Dhm", "IV_Dhm", "Dhm_IV"))

lm_mating_duration <- lm(mating_duration ~ cross, data = IV_Dhm)
summary(lm_mating_duration)
pairs(emmeans(lm_mating_duration, ~cross, type = "response"), adjust = "none")

ggplot(data = subset(IV_Dhm, mating_duration > 0), aes(x = cross, y = mating_duration)) +
  geom_boxplot()

# dahomey versus lhm 
Dhm_LHm <- subset(ncsi, cross %in% c("LHm-GFP_LHm-GFP", "Dhm_Dhm", "LHm-GFP_Dhm", "Dhm_LHm-GFP"))

lm_mating_duration <- lm(mating_duration ~ cross, data = Dhm_LHm)
summary(lm_mating_duration)
pairs(emmeans(lm_mating_duration, ~cross, type = "response"), adjust = "none")

ggplot(data = subset(Dhm_LHm, mating_duration > 0), aes(x = cross, y = mating_duration)) +
  geom_boxplot()

# dahomey versus val
Dhm_Val <- subset(ncsi, cross %in% c("Val_Val", "Dhm_Dhm", "Val_Dhm", "Dhm_Val"))

lm_mating_duration <- lm(mating_duration ~ cross, data = Dhm_Val)
summary(lm_mating_duration)
pairs(emmeans(lm_mating_duration, ~cross, type = "response"), adjust = "none")

ggplot(data = subset(Dhm_Val, mating_duration > 0), aes(x = cross, y = mating_duration)) +
  geom_boxplot()

# iv versus lhm
iv_lhm <- subset(ncsi, cross %in% c("IV_IV", "LHm-GFP_LHm-GFP", "IV_LHm-GFP", "LHm-GFP_IV"))

lm_mating_duration <- lm(mating_duration ~ cross, data = iv_lhm)
summary(lm_mating_duration)
pairs(emmeans(lm_mating_duration, ~cross, type = "response"), adjust = "none")

ggplot(data = subset(iv_lhm, mating_duration > 0), aes(x = cross, y = mating_duration)) +
  geom_boxplot()

# iv versus val
iv_val <- subset(ncsi, cross %in% c("IV_IV", "Val_Val", "IV_Val", "Val_IV"))

lm_mating_duration <- lm(mating_duration ~ cross, data = iv_val)
summary(lm_mating_duration)
pairs(emmeans(lm_mating_duration, ~cross, type = "response"), adjust = "none")

ggplot(data = subset(iv_val, mating_duration > 0), aes(x = cross, y = mating_duration)) +
  geom_boxplot()

# val versus lhm
val_lhm <- subset(ncsi, cross %in% c("Val_LHm-GFP", "Val_Val", "LHm-GFP_Val", "LHm-GFP_LHm-GFP"))

lm_mating_duration <- lm(mating_duration ~ cross, data = val_lhm)
summary(lm_mating_duration)
pairs(emmeans(lm_mating_duration, ~cross, type = "response"), adjust = "none")

ggplot(data = subset(val_lhm, mating_duration > 0), aes(x = cross, y = mating_duration)) +
  geom_boxplot()


