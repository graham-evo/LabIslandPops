# Graham C. McLaughlin
# Created: 09-01-2023
# Modified: 12-12-2023
# Analysis of August 26th and October 23rd experiment characterizing RI between our four Lab Island Pops
# ----

# Libraries:
##---------------
library(lme4)
library(coxme)
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
library(gplots)
library(reshape2)
library(cowplot)
##---------------

# DATA CLEANING/MANAGEMENT:
##----------------------------
ncsi_b1 <- read.table("LabIslandPops_premating_R.txt", # August Experiment
                   header = TRUE,
                   comment.char = "",
                   na.strings = "#N/A")
ncsi_b1$block <- 1
tail(ncsi_b1)
head(ncsi_b1)


ncsi_b2 <- read.table("oct23_block/premating_oct23.txt", # October Experiment
                      header = TRUE,
                      comment.char = "",
                      na.strings = "#N/A")
ncsi_b2$block <- 2
head(ncsi_b2)
tail(ncsi_b2)

# Combine blocks:
ncsi <- bind_rows(ncsi_b1,ncsi_b2)
ncsi$id <- seq(1,800,1)
head(ncsi)
tail(ncsi)

# Let's first define our factors:
ncsi$female <- factor(ncsi$female,
                      levels = c("iv","lhm","dhm","val"),
                      labels = c("IV", "LHm", "Dhm", "Val"))
ncsi$male <- factor(ncsi$male,
                    levels = c("iv","lhm","dhm","val"),
                    labels = c("IV","LHm","Dhm","Val"))
head(ncsi) # Factors leveled and re-labelled

# Let's add a cross column to make comparisons easier:
ncsi$cross <- paste0(ncsi$female, "_", ncsi$male)
head(ncsi) # We now have a cross column

# Trimming out the vials lost due to experimental failure that won't be used for analysis of any measures:
dim(ncsi[ncsi$block == 1,]) # 400 observations before trimming in b1
dim(ncsi[ncsi$block == 2,]) # 400 observations before trimming in b2
ncsi <- subset(ncsi, failure != 1)
dim(ncsi) # Lost 56 to experimental failure, most of them from block 1

# PHENOTYPE MEASUREMENTS:
## ----------------------

# Looking first at WILLINGNESS TO MATE:
# Starting with just the simplest model assessing the effect of cross on probability of mating
# We will add in the effect of male and female and any interaction terms later.
# WILLINGNESS BLOCK 1:
hist(ncsi$mated[ncsi$block == 1]) # Binomial distribution of course
hist(ncsi$mated[ncsi$block == 2]) # Binomial distribution of course

glm_willingness_b1 <- glm(mated ~ cross + male + female, #These male and female terms are only here so that we can plot the emmeans df by male and female later. They don't alter the fitting of the cross variable
                       family = binomial,
                       data = ncsi[ncsi$block == 1,])

glm_willingness_b1
summary(glm_willingness_b1) # moderate effects but no obvious patterns
anova(glm_willingness_b1, test = "Chisq")
emmeans(glm_willingness_b1, ~cross, type = "response", level = 0.69)

# Modeling male/female effects:
glm_willingness_interaction <- glm(mated ~ male*female, 
                                   family = binomial(),
                                   data = ncsi[ncsi$block == 1,])
glm_willingness_interaction # Strong effect of female and male LHm plus interactions
summary(glm_willingness_interaction)
anova(glm_willingness_interaction, test = "Chisq")

# Let's plot block 1:
plot_frame_b1 <- as.data.frame(emmeans(glm_willingness_b1, ~cross + male + female, type = "response", level = 0.69)) 
plot_frame_b1 %>%
  arrange(male) %>%
  mutate(female = factor(female, levels=c("Val","Dhm","LHm","IV"))) %>%
  ggplot(aes(x = male, y = female)) +
  geom_raster(aes(fill = prob)) +
  scale_fill_gradient(low = "grey100", high = "red") +
  labs(x = "Male Genotype", y = "Female Genotype", fill = "Probability\nof\nMating") +
  theme_bw() +
  scale_x_discrete(position = "top")

# Different way of visualizing differences:
ggplot() +
  geom_errorbar(aes(x = cross, ymin = prob-SE, ymax = prob + SE), data = plot_frame_b1) +
  geom_point(data = plot_frame_b1, aes(x = cross, y = prob, color = male)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# WILLINGESS BLOCK 2:
# Modelling the effect of cross
glm_willingness_b2 <- glm(mated ~ cross + male + female,
                          family = binomial,
                          data = ncsi[ncsi$block == 2,])

glm_willingness_b2
summary(glm_willingness_b2) # Not as many differences between crosses as block 1
anova(glm_willingness_b2, test = "Chisq")
emmeans(glm_willingness_b2, ~cross, type = "response", level = 0.69)
pairs(emmeans(glm_willingness_b2, ~cross, type = "response"), adjust = "none")

# Modelling male and female effects
glm_willingness_interaction <- glm(mated ~ male*female, 
                                   family = binomial(),
                                   data = ncsi[ncsi$block == 2,])
glm_willingness_interaction
summary(glm_willingness_interaction) # LHm male and female effects disappear in block 2. Male DHm and Val drive the differences here
anova(glm_willingness_interaction, test = "Chisq")

# Plotting willingness to mate:
plot_frame_b2 <- as.data.frame(emmeans(glm_willingness_b2, ~cross + male + female, type = "response", level = 0.69))
plot_frame_b2 %>%
  arrange(male) %>%
  mutate(female = factor(female, levels=c("Val","Dhm","LHm","IV"))) %>%
  ggplot(aes(x = male, y = female)) +
  geom_raster(aes(fill = prob)) +
  scale_fill_gradient(low = "grey100", high = "red") +
  labs(x = "Male Genotype", y = "Female Genotype", fill = "Probability\nof\nMating") +
  theme_bw() +
  scale_x_discrete(position = "top")

# Dot plot:
ggplot() +
  geom_errorbar(aes(x = cross, ymin = prob-SE, ymax = prob + SE), data = plot_frame_b2) +
  geom_point(data = plot_frame_b2, aes(x = cross, y = prob, color = male)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# WILLINGNESS COMBINED BLOCK ANALYSIS
glmer_willingness <- mixed(mated ~ cross + male + female + (1|block),
                       method = "LRT",
                       data = ncsi)
glmer_willingness_null <- mixed(mated ~ cross + (1|block),
                           method = "LRT",
                           data = ncsi)
anova(glmer_willingness, glmer_willingness_null) # same models, just checking to make sure male and female terms don't mess up model fit

# Looking at effect of cross:
glmer_willingness
summary(glmer_willingness)
anova(glmer_willingness)
emmeans(glmer_willingness, ~cross, type = "response", level = 0.69)
pairs(emmeans(glmer_willingness, ~cross, type = "response"), adjust = "none")

# Looking at effect of male and female:
glmer_willingness_interaction_null <- glm(mated ~ male*female,
                                       family = binomial(),
                                       data = ncsi)
glmer_willingness_interaction_alt <- mixed(mated ~ male*female + (1|block),
                                            method = "LRT",
                                            data = ncsi)
anova(glmer_willingness_interaction_null, glmer_willingness_interaction_alt,
      test = "Chisq") # Blocks differ by male female and interactions???

#Comparing blocks
plot_frame_combined_blocks <- as.data.frame(emmeans(glmer_willingness, ~cross + male + female, type = "response", level = 0.69))

plot_frame_combined_blocks %>%
  arrange(male) %>%
  mutate(female = factor(female, levels=c("Val","Dhm","LHm","IV"))) %>%
  ggplot(aes(x = male, y = female)) +
  geom_raster(aes(fill = emmean)) +
  scale_fill_gradient(low = "grey100", high = "red") +
  labs(x = "Male Genotype", y = "Female Genotype", fill = "Probability\nof\nMating") +
  theme_bw() +
  scale_x_discrete(position = "top")

# Comparing blocks:
plot(plot_frame_b1$prob, plot_frame_b2$prob, type = "p",
     xlab = "Block 1 Proportion Mated", ylab = "Block 2 Proportion Mated",
     xlim = c(0,1), ylim = c(0,1))
cor(plot_frame_b1$prob,plot_frame_b2$prob) # Super inconsistent results between blocks

lm_block <- lm(mated ~ cross*block, data = ncsi)
summary(lm_block) # Which crosses differed in willingness to mate between blocks?

# Different way of visualizing inconsistencies between blocks
ggplot() +
  geom_point(data = plot_frame_b1, aes(x = cross, y = prob), color = "red") +
  geom_point(data = plot_frame_b2, aes(x = cross, y = prob), color = "blue") +
  geom_errorbar(aes(x = cross, ymin = prob-SE, ymax = prob + SE), data = plot_frame_b1) +
  geom_errorbar(aes(x = cross, ymin = prob-SE, ymax = prob + SE), data = plot_frame_b2, color = "red") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# It's easier to visualize and model if we only look at plots of the within-pop
##cross versus a single across population cross.
# We can then run a model with male/female and interaction effects
# To make things simpler I'm going to first create subsetted data frames with only the
# Population cross comparisons of interest:
iv_dhm <- subset(ncsi, cross %in% c("IV_IV", "Dhm_Dhm", "IV_Dhm", "Dhm_IV"))
dhm_lhm <- subset(ncsi, cross %in% c("LHm-GFP_LHm-GFP", "Dhm_Dhm", "LHm-GFP_Dhm", "Dhm_LHm-GFP"))
dhm_val <- subset(ncsi, cross %in% c("Val_Val", "Dhm_Dhm", "Val_Dhm", "Dhm_Val"))
iv_lhm <- subset(ncsi, cross %in% c("IV_IV", "LHm-GFP_LHm-GFP", "IV_LHm-GFP", "LHm-GFP_IV"))
iv_val <- subset(ncsi, cross %in% c("IV_IV", "Val_Val", "IV_Val", "Val_IV"))
val_lhm <- subset(ncsi, cross %in% c("Val_Val", "LHm-GFP_LHm-GFP", "Val_LHm-GFP", "LHm-GFP_Val"))
# For each subdataset above we can plug the comparison of interest into the models below:

# Combined block glmer models:
glmer_willingness <- glm(mated ~ cross + (1|block),
                              family = binomial,
                              data = iv_dhm)
summary(glmer_willingness)
anova(glmer_willingness, test = "Chisq")
emmeans(glmer_willingness, ~cross, type = "response", level = 0.69, data = iv_dhm)
pairs(emmeans(glmer_willingness, ~cross, type = "response", data = iv_dhm), adjust = "none")

plot_frame <- as.data.frame(emmeans(glmer_willingness,
                                    ~cross,
                                    type = "response",
                                    data = iv_dhm))

ggplot(aes(x = cross, y = prob), data = plot_frame) +
  geom_point() +
  geom_errorbar(aes(x = cross, ymin = prob-SE, ymax = prob + SE))

# Adding male/female and interaction term to cross comparison model:
glm_willingness_interaction <- glm(mated ~ male*female,
                              family = binomial,
                              data = iv_dhm)
summary(glm_willingness_interaction)

# Plotting Joint Isolation Indices of each genetic background comparison (Merrell 1950)
# This is just from the first block of data, should redo with the second block and combined.
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

# LATENCY TO MATE:
##----------------------------
## Moving on to LATENCY TO MATE:
## Let's read in our data frames again to avoid any mistakes

ncsi_b1 <- read.table("LabIslandPops_premating_R.txt", # August Experiment
                      header = TRUE,
                      comment.char = "",
                      na.strings = "#N/A")
ncsi_b1$block <- 1
tail(ncsi_b1)


ncsi_b2 <- read.table("oct23_block/premating_oct23.txt", # October Experiment
                      header = TRUE,
                      comment.char = "",
                      na.strings = "#N/A")
ncsi_b2$block <- 2
tail(ncsi_b2)

# Combine blocks:
ncsi <- bind_rows(ncsi_b1,ncsi_b2)
ncsi$id <- seq(1,800,1)
tail(ncsi)

# Let's first define our factors:
ncsi$female <- factor(ncsi$female,
                      levels = c("iv","lhm","dhm","val"),
                      labels = c("IV", "LHm", "Dhm", "Val"))
ncsi$male <- factor(ncsi$male,
                    levels = c("iv","lhm","dhm","val"),
                    labels = c("IV","LHm","Dhm","Val"))
head(ncsi) # Factors levelled and relabelled

# Let's add a cross column to make comparisons easier:
ncsi$cross <- paste0(ncsi$female, "_", ncsi$male)
head(ncsi) # We now have a cross column
dim(ncsi) # 800 observations before trimming

ncsi <- subset(ncsi, failure == 0) # Get rid of experimental failures
head(ncsi)
tail(ncsi)
dim(ncsi)

ncsi$failure == TRUE # Double checking that all failures are removed:

# Setting all did not mate (dnm) vials to 180:
ncsi$latency[ncsi$dnm==1] <- 180
head(ncsi, n =100)
ncsi$latency < 0
#View(ncsi) #Vials 113 and 248 have negative latencies even though they are coded
# as having mated, so lets just remove them for now
ncsi <- subset(ncsi, latency > 0)
ncsi$latency<0 # negative latency observations removed

# LATENCY TO MATE BLOCK 1:
# Effect of cross:
coxph_latency_b1 <- 
  coxph(formula = Surv(latency, mated)
                          ~ cross + male + female,
                          data = ncsi[ncsi$block == 1,])
coxph_latency_b1
summary(coxph_latency_b1)
anova(coxph_latency_b1)
emmeans(coxph_latency_b1, ~cross, type = "response")

# Effect of male and female:
coxph_latency_b1_interaction <- coxph(formula = Surv(latency, mated)
        ~ male*female,
        data = ncsi[ncsi$block == 1,])
coxph_latency_b1_interaction # Strong effects of LHm males and females
summary(coxph_latency_b1_interaction)
anova(coxph_latency_b1_interaction)

# BLOCK 2 LATENCY TO MATE:
# Effect of cross:
coxph_latency_b2 <- 
  coxph(formula = Surv(latency, mated)
        ~ cross + male + female,
        data = ncsi[ncsi$block == 2,])
coxph_latency_b2
summary(coxph_latency_b2)
anova(coxph_latency_b2)

# Effect of male and female:
coxph_latency_b2_interaction <- coxph(formula = Surv(latency, mated)
                                      ~ male*female,
                                      data = ncsi[ncsi$block == 2,])
coxph_latency_b2_interaction # Strong effects of male DHm and Val
summary(coxph_latency_b2_interaction)
anova(coxph_latency_b2_interaction)

# Running a combined model:
coxph_latency <- coxme(formula = Surv(latency, mated) ~ cross + (1|block), data = ncsi)
coxph_latency
emmeans(coxph_latency, ~cross, type = "response")
pairs(emmeans(coxph_latency, ~cross, type = "response"), adjust = "none")

# Block by block comparison plot:
ggplot(data = ncsi, aes(x = cross, y = latency, fill = female)) +
  geom_boxplot() +
  facet_wrap(~block) +
  theme_bw()

# Combined blocks
ggplot(data = ncsi, aes(x = cross, y = latency, fill = female)) +
  geom_boxplot()

# Correlation between blocks:

latencyb1<-as.data.frame(emmeans(coxph_latency_b2, ~cross, type = "response"))
latencyb2<-as.data.frame(emmeans(coxph_latency_b1, ~cross, type = "response"))
ncsi_test <- inner_join(ncsi[ncsi$block==1,],ncsi[ncsi$block==2,], by = "id")
plot(latencyb1$response,latencyb2$response[1:15], type = "p",
     xlab = "Block 1 Latency", ylab = "Block 2 Latency")
cor(latencyb1$response,latencyb2$response[1:15]) # Super inconsistent results between blocks

# Heatmap
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

#-----------------
# Mating Duration:
## Let's read in our data frames again to avoid any mistakes

ncsi_b1 <- read.table("LabIslandPops_premating_R.txt", # August Experiment
                      header = TRUE,
                      comment.char = "",
                      na.strings = "#N/A")
ncsi_b1$block <- 1
tail(ncsi_b1)


ncsi_b2 <- read.table("oct23_block/premating_oct23.txt", # October Experiment
                      header = TRUE,
                      comment.char = "",
                      na.strings = "#N/A")
ncsi_b2$block <- 2
tail(ncsi_b2)

# Combine blocks:
ncsi <- bind_rows(ncsi_b1,ncsi_b2)
ncsi$id <- seq(1,800,1)
tail(ncsi)

# Let's first define our factors:
ncsi$female <- factor(ncsi$female,
                      levels = c("iv","lhm","dhm","val"),
                      labels = c("IV", "LHm", "Dhm", "Val"))
ncsi$male <- factor(ncsi$male,
                    levels = c("iv","lhm","dhm","val"),
                    labels = c("IV","LHm","Dhm","Val"))
head(ncsi) # Factors levelled and relabelled

# Let's add a cross column to make comparisons easier:
ncsi$cross <- paste0(ncsi$female, "_", ncsi$male)
head(ncsi) # We now have a cross column
dim(ncsi) # 800 observations before trimming

ncsi <- subset(ncsi, failure == 0) # Get rid of experimental failures
head(ncsi)
tail(ncsi)
dim(ncsi)

ncsi$failure == 1 # Double checking that all failures are removed:

# Removing the vials for which they didn't mate:
ncsi <- subset(ncsi, mating_duration > 0)
ncsi$mating_duration > 0
View(ncsi)
head(ncsi, n = 400)
dim(ncsi) # 503 mated

hist(ncsi$mating_duration, breaks = 100) # Clear outliers here
# Looking for outliers based on QQ residuals:
lm_mating_duration <- lm(mating_duration ~ cross, data = ncsi)
plot(lm_mating_duration, 2) # Q-Q Plot
boxplot(ncsi$mating_duration)

# Removing Outliers for mating duration:
# Vial 269 = 2 minutes
# Vial 302 = 48 minutes
# Vial 436 = 120 minutes
ncsi <- subset(ncsi, id != c(269,302,436))
# Run this again make sure they are gone and now look at distribution
hist(ncsi$mating_duration, breaks = 100) # Clear outliers here
lm_mating_duration <- lm(mating_duration ~ cross, data = ncsi)
plot(lm_mating_duration, 2) # Q-Q Plot
boxplot(ncsi$mating_duration)

# MATING DURATION BLOCK 1:
# Effect of cross:
lm_mating_duration_b1 <- lm(mating_duration ~ cross,
                            data = ncsi[ncsi$block==1,])
summary(lm_mating_duration_b1)

# Effect of male and female:
lm_mating_duration_b1_interaction <- lm(mating_duration ~ male*female,
                            data = ncsi[ncsi$block==1,])
summary(lm_mating_duration_b1_interaction) # Really strong male effectss of mating duration and strong LHm female effects of mating duration

# MATING DURATION BLOCK 2:
# Effect of cross:
lm_mating_duration_b2 <- lm(mating_duration ~ cross,
                            data = ncsi[ncsi$block==2,])
summary(lm_mating_duration_b2)

# Effect of male and female:
lm_mating_duration_b2_interaction <- lm(mating_duration ~ male*female,
                                        data = ncsi[ncsi$block==2,])
summary(lm_mating_duration_b2_interaction) # Really strong male effects of mating duration. Effect of LHm female disappeared from block 2

# Comparing blocks visually with boxplots:
ggplot(data = ncsi, aes(x = cross, y = mating_duration, fill = male)) +
  geom_boxplot() +
  theme_classic() +
  facet_wrap(~block)

# Comparing blocks visually with Heatmaps:
mating_duration_plotframe_b1 <- as.data.frame(aggregate(mating_duration~male+female, data = ncsi[ncsi$block==1,], FUN = mean))
mating_duration_plotframe_b2 <- as.data.frame(aggregate(mating_duration~male+female, data = ncsi[ncsi$block==2,], FUN = mean))

mating_duration_plotframe_b2 %>%
  arrange(male) %>%
  mutate(female = factor(female, levels=c("Val","Dhm","LHm","IV"))) %>%
  ggplot(aes(x = male, y = female)) +
  geom_raster(aes(fill = mating_duration)) +
  scale_fill_gradient(low = "grey100", high = "red") +
  labs(x = "Male Genotype", y = "Female Genotype", fill = "Mating\nDuration\n(min)") +
  theme_bw() +
  scale_x_discrete(position = "top")
# Combined block analysis
lmer_mating_duration <- lmer(mating_duration ~ male*female + (1|block),
     data = ncsi)
summary(lmer_mating_duration)
anova(lmer_mating_duration)

# Finally lets compare repeatability of mating duration measures between blocks:
plot(ncsi_b1$mating_duration, ncsi_b2$mating_duration,
     xlim = c(0,40), ylim = c(0,40)) # correlation if i remove all the failures and non-maters it looks like
anova(lm_mating_duration_b1_interaction,lm_mating_duration_b2_interaction)

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


