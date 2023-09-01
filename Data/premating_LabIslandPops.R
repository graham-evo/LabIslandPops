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
# Starting with just the simplest model assessing the effect of cross on probabiltiy of mating
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

# It's easier to visualize if we only look at plots of the within-pop cross versus
# all crosses that contain the within-pop genotype:

# Starting with Dahomey:
ggplot(data = subset(plot_frame, male == "Dhm" | female == "Dhm"), aes(x = cross, y = prob)) +
  geom_point() +
  geom_errorbar(aes(x = cross, ymin = prob-SE, ymax = prob + SE))

# IV:
ggplot(data = subset(plot_frame, male == "IV" | female == "IV"), aes(x = cross, y = prob)) +
  geom_point() +
  geom_errorbar(aes(x = cross, ymin = prob-SE, ymax = prob + SE))

# Val: 
ggplot(data = subset(plot_frame, male == "Val" | female == "Val"), aes(x = cross, y = prob)) +
  geom_point() +
  geom_errorbar(aes(x = cross, ymin = prob-SE, ymax = prob + SE))

# LHm-GFP:
ggplot(data = subset(plot_frame, male == "LHm-GFP" | female == "LHm-GFP"), aes(x = cross, y = prob)) +
  geom_point() +
  geom_errorbar(aes(x = cross, ymin = prob-SE, ymax = prob + SE))

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
ggplot(data = ncsi, aes(x = cross, y = latency)) +
  geom_boxplot()

# Plotting each within-population cross against only meaningful/relevent cross comparisons:
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

