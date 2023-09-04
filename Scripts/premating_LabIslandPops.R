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

ncsi <- read.table("LabIslandPops_premating_R.txt",
                   header = TRUE,
                   comment.char = "",
                   na.strings = "#N/A")
head(ncsi)

# Let's first define our factors:
ncsi$female <- factor(ncsi$female,
                      levels = c("iv","lhm","dhm","val"),
                      labels = c("IV", "LHm-GFP", "Dhm", "Val"))
ncsi$male <- factor(ncsi$male,
                    levels = c("iv","lhm","dhm","val"),
                    labels = c("IV","LHm-GFP","Dhm","Val"))
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

glm_willingness <- glm(mated ~ cross + male + female,
                       family = binomial,
                       data = ncsi)
glm_willingness
summary(glm_willingness)
emmeans(glm_willingness, ~cross + male + female, type = "response", level = 0.69)
pairs(emmeans(glm_willingness, ~cross, type = "response"), adjust = "none")
# We aren't interested in all of the above pairwise comparisons so this is not too helpful.

plot_frame <- as.data.frame(emmeans(glm_willingness, ~cross, type = "response", level = 0.69))

ggplot(data = plot_frame, aes(x = cross, y = prob)) +
  geom_point() +
  geom_errorbar(aes(x = cross, ymin = prob-SE, ymax = prob + SE))

# It's easier to visualize and model if we only look at plots of the within-pop cross versus a single across population cross. We can then run a model with male/female and interaction effects:

## Starting with Dahomey versus IV:
IV_Dhm <- subset(ncsi, cross %in% c("IV_IV", "Dhm_Dhm", "IV_Dhm", "Dhm_IV"))

glm_IV_Dhm_willingness <- glm(mated ~ cross,
                              family = binomial,
                              data = IV_Dhm)
summary(glm_IV_Dhm_willingness)
emmeans(glm_IV_Dhm_willingness, ~cross, type = "response")

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
pairs(emmeans(glm_Val_Dhm_willingness, ~cross, type = "response"))

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
pairs(emmeans(glm_iv_lhm_willingness, ~cross, type = "response"))

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
pairs(emmeans(glm_iv_val_willingness, ~cross, type = "response"))

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
pairs(emmeans(glm_val_lhm_willingness, ~cross, type = "response"))

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
coxph_latency <- coxph(formula = Surv(latency, mated) ~ cross + male + female, data = ncsi)
coxph_latency
summary(coxph_latency)
anova(coxph_latency)
emmeans(coxph_latency, ~cross, type = "response")
pairs(emmeans(coxph_latency, ~cross, type = "response"))

# Plot of all cross latencies:
ggplot(data = subset(ncsi, latency != 180), aes(x = cross, y = latency)) +
  geom_boxplot()

# dhm versus iv latency
ggplot(data = subset(IV_Dhm, latency > 0), aes(x = cross, y = latency)) +
  geom_boxplot()

# dahomey versus lhm latency
ggplot(data = subset(Dhm_LHm, latency > 0), aes(x = cross, y = latency)) +
  geom_boxplot()

# dahomey versus val
ggplot(data = subset(val_dhm, latency > 0), aes(x = cross, y = latency)) +
  geom_boxplot()

# iv versus lhm
ggplot(data = subset(iv_lhm, latency > 0), aes(x = cross, y = latency)) +
  geom_boxplot()

# iv versus val
ggplot(data = subset(iv_val, latency > 0), aes(x = cross, y = latency)) +
  geom_boxplot()

# val versus lhm
ggplot(data = subset(val_LHm, latency > 0), aes(x = cross, y = latency)) +
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

lm_mating_duration <- lm(mating_duration ~ cross, data = ncsi)
summary(lm_mating_duration)

ggplot(data = ncsi, aes(x = cross, y = mating_duration)) +
  geom_boxplot()

# individual comparisons:
# Plot of all cross mating durations:

# dhm versus iv
ggplot(data = subset(IV_Dhm, mating_duration > 0), aes(x = cross, y = mating_duration)) +
  geom_boxplot()

# dahomey versus lhm 
ggplot(data = subset(Dhm_LHm, mating_duration > 0), aes(x = cross, y = mating_duration)) +
  geom_boxplot()

# dahomey versus val
ggplot(data = subset(val_dhm, mating_duration > 0), aes(x = cross, y = mating_duration)) +
  geom_boxplot()

# iv versus lhm
ggplot(data = subset(iv_lhm, mating_duration > 0), aes(x = cross, y = mating_duration)) +
  geom_boxplot()

# iv versus val
ggplot(data = subset(iv_val, mating_duration > 0), aes(x = cross, y = mating_duration)) +
  geom_boxplot()

# val versus lhm
ggplot(data = subset(val_LHm, mating_duration > 0), aes(x = cross, y = mating_duration)) +
  geom_boxplot()


