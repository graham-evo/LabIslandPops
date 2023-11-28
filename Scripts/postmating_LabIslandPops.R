# Graham C. McLaughlin
# Date Created: 09-04-2023
# Date Modified: 09-04-2023
# Analysis of postmating egg counts and egg-adult-viability
## ----------

library(emmeans)
library(tidyverse)

# Read in egg count data
pmpz <- read.table(file = "Data/LabIslandCrosses_Postmating_R.txt",
                   header = TRUE,
                   comment.char = "",
                   na.strings = "#N/A")

# Let's first define our factors:
pmpz$female <- factor(pmpz$female,
                      levels = c("iv","lhm","dhm","val"),
                      labels = c("IV", "LHm", "Dhm", "Val"))
pmpz$male <- factor(pmpz$male,
                    levels = c("iv","lhm","dhm","val"),
                    labels = c("IV","LHm","Dhm","Val"))
head(ncsi) # Factors defined and relabelled

# Let's add a cross column to make comparisons easier:
pmpz$cross <- paste0(pmpz$female, "_", pmpz$male)
head(pmpz) # We now have a cross column

# Let's remove all of the rows where failure == 1, since these should be removed
# from all downstream analyses.
dim(pmpz) # 584 observations before trimming
pmpz <- subset(pmpz, failure != 1)
dim(pmpz) # Lost 24 due to experimental failure

str(pmpz) # eggs need to be a number
pmpz$eggs <- as.integer(pmpz$eggs)
str(pmpz)

# Let's start by making a new data frame with total egg counts per vial across the 4 days
cumulative_pmpz <- pmpz %>%
  group_by(vial, cross, male, female) %>%
  summarise(eggs = sum(eggs)) %>%
  ungroup()

# And also remove any females for which less than 10 eggs were laid (assuming a non-successful mating)
dim(cumulative_pmpz)
cumulative_pmpz <- subset(cumulative_pmpz, eggs > 10)
dim(cumulative_pmpz)

# Let's look at the egg count distributions of each of our crosses:
egg_dist <- ggplot(data = cumulative_pmpz, aes(eggs, fill = cross)) +
  geom_density() +
  facet_wrap(~cross)

# Running a model looking at total eggs:
total_eggs_lm <- lm(eggs ~ male*female, data = cumulative_pmpz)
summary(total_eggs_lm)
anova(total_eggs_lm)

emmeans(total_eggs_lm, ~cross, type="response")
pairs(emmeans(total_eggs_lm, ~cross, type="response"))

plot_frame <- as.data.frame(emmeans(total_eggs_lm, ~cross + female + male, type="response"))

ggplot(data = plot_frame, aes(x = cross,  y = emmean, color = female)) +
  geom_errorbar(aes(x = cross, ymax = emmean+SE, ymin=emmean-SE)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Making Heat Map of Total Eggs Laid:
heatmap.df <- plot_frame %>%
  select(male, female, emmean)

heatmap.df %>%
  arrange(male) %>%
  mutate(female = factor(female, levels=c("Val","Dhm","LHm","IV"))) %>%
  ggplot(aes(x = male, y = female)) +
  geom_raster(aes(fill = emmean)) +
  scale_fill_gradient(low = "grey90", high = "red") +
  labs(x = "Male Genotype", y = "Female Genotype", fill = "Eggs\nLaid") +
  theme_bw() +
  scale_x_discrete(position = "top")

# Day by day eggs:

days_plotFrame <- as.data.frame (aggregate(eggs ~ cross+day, data = pmpz, FUN = summary)) #Produce our initial summary statistics of eggs laid per day across treatments
days_plotFrame$n <- aggregate(eggs ~ cross*day, data = pmpz, FUN = "length")[,3] #Get the sample size of each treatment per day
days_plotFrame$var <- aggregate(eggs ~ cross*day, data = pmpz, FUN = "var")[,3] #Calculate the variance
days_plotFrame$sd <- sqrt(days_plotFrame$var) #Standard deviation
days_plotFrame$se <- days_plotFrame$sd/sqrt(days_plotFrame$n) #Standard error or SD divided by variance

#Subset by population comparison:
days_plotFrame <- subset(days_plotFrame, cross %in% c("IV_IV", "Dhm_Dhm", "IV_Dhm", "Dhm_IV"))
days_plotFrame <- subset(days_plotFrame, cross %in% c("LHm-GFP_LHm-GFP", "Dhm_Dhm", "LHm-GFP_Dhm", "Dhm_LHm-GFP"))
days_plotFrame <- subset(days_plotFrame, cross %in% c("Val_Val", "Dhm_Dhm", "Val_Dhm", "Dhm_Val"))
days_plotFrame <- subset(days_plotFrame, cross %in% c("IV_IV", "Val_Val", "IV_Val", "Val_IV"))



eggCountTimeSeriesPlot <- ggplot(data=days_plotFrame, aes(x=day, y=eggs[,'Mean'], fill=cross, color=cross, group = cross)) +
  geom_line(size=0.5, data=days_plotFrame, aes(x=day, y=eggs[,'Mean'], color=cross)) +
  geom_errorbar(color="black", size=0.5, width=0.1, data=days_plotFrame, aes(ymax=eggs[,'Mean']+se, ymin=eggs[,'Mean']-se)) +
  geom_point(color="black", size=2, stroke=1, shape=21, aes(fill=cross)) +
  scale_y_continuous(name="Eggs laid", limits=c(0,32), breaks=seq(0,30,5), expand=c(0,0)) +
  scale_x_continuous(name="Day", breaks=seq(1,6,1)) +
  labs(fill = "Cross") +
  theme_bw() +
  theme(legend.position = "top") +
  theme(axis.line.x = element_line(size=1), axis.line.y = element_line(size=1)) +
  theme(axis.ticks.x = element_line(size=1), axis.ticks.y = element_line(size=1)) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(colour="black", size=1)) +
  theme(axis.title.x=element_text(size=16, color="black"), axis.title.y=element_text(size=16, color="black", margin=margin(1,1,1,1,"pt"))) +
  theme(axis.text.x = element_text(size=14, margin=margin(2,2,2,2,"pt"))) +
  theme(axis.text.y = element_text(size=14, margin=margin(2,2,2,2,"pt"))) +
  theme(plot.margin = unit(c(5, 5, 5, 5), "points"))
eggCountTimeSeriesPlot 

# Looking at egg to adult viability of days 1 and 2
pzi <- read.table(file = "data/LabIslandCrosses_EggToAdultViability_R.txt",
                  header = TRUE,
                  comment.char = "",
                  na.strings = "#N/A", fill = TRUE)
pzi$female <- factor(pzi$female,
                      levels = c("iv","lhm","dhm","val"),
                      labels = c("IV", "LHm", "Dhm", "Val"))
pzi$male <- factor(pzi$male,
                    levels = c("iv","lhm","dhm","val"),
                    labels = c("IV","LHm","Dhm","Val"))
pzi$cross <- paste0(pzi$female,"_",pzi$male)

pzi <- subset(pzi, failure == 0)
pzi <- subset(pzi, eggs > 10)
pzi$prop_eggs_emerged <- as.numeric(pzi$prop_eggs_emerged)

hist(pzi$prop_eggs_emerged, breaks = 100)

ggplot(data = pzi) +
  geom_density(aes(x = prop_eggs_emerged, color = female)) +
  facet_wrap(~cross)

lm_viability <- lm(prop_eggs_emerged ~ female*male, data = pzi)
summary(lm_viability)
anova(lm_viability)

emmeans(lm_viability, ~ male + female, type = "response")

plot_frame <- as.data.frame(emmeans(lm_viability, ~ male + female, type = "response"))

viability_plot <- ggplot(data = plot_frame, aes(y = emmean, x = female, group = male, fill = male)) +
  geom_errorbar(size = 1, width = 0.1, aes(ymin = emmean - SE, ymax = emmean + SE))+
  geom_point(size = 3, stroke = 1, shape = 21, color = "black")

# Heatmap of egg to adult viability:

heatmap.df <- plot_frame %>%
  select(male, female, emmean)

heatmap.df %>%
  arrange(male) %>%
  mutate(female = factor(female, levels=c("Val","Dhm","LHm","IV"))) %>%
  ggplot(aes(x = male, y = female)) +
  geom_raster(aes(fill = emmean)) +
  scale_fill_gradient(low = "grey90", high = "red") +
  labs(x = "Male Genotype", y = "Female Genotype", fill = "Egg-to-Adult\nViability") +
  theme_bw() +
  scale_x_discrete(position = "top")

# Development Time
dev_time <- read.table(file = "data/LabIslandCrosses_DevelopmentTime_R.txt",
                  header = TRUE,
                  comment.char = "",
                  na.strings = "#N/A", fill = TRUE)
dev_time$female <- factor(dev_time$female,
                     levels = c("iv","lhm","dhm","val"),
                     labels = c("IV", "LHm", "Dhm", "Val"))
dev_time$male <- factor(dev_time$male,
                   levels = c("iv","lhm","dhm","val"),
                   labels = c("IV","LHm","Dhm","Val"))
dev_time$cross <- paste0(dev_time$female,"_",dev_time$male)


dev_time <- subset(dev_time, vial != 326)

dev_time$emerged <- as.numeric(dev_time$emerged)
dev_time$day <- as.integer(dev_time$day)

# Calculate development time
dev_time$weightedemerged <- dev_time$emerged*dev_time$day
for (i in dev_time$vial){
  dev_time$finaldeveloptime[dev_time$vial == i] <-
    sum(dev_time[dev_time$vial==i,]$weightedemerged)/
    sum(dev_time[dev_time$vial==i,]$emerged)
}

dev_time$vial <- dev_time$vial[dev_time$finaldeveloptime != 0]

# Development time of each brood for each day development time ~ cross + eggs

lm_devTime <- lm(finaldeveloptime ~ cross + eggs +male+female, data = dev_time)
summary(lm_devTime)
lm_devTime_sex <- lm(finaldeveloptime ~ male*female, data = dev_time)
summary(lm_devTime_sex)
anova(lm_devTime_sex)

emmeans(lm_devTime, ~cross + male + female, type = "response")

heatmap.df <- as.data.frame(emmeans(lm_devTime, ~cross + male + female, type = "response"))

heatmap.df %>%
  arrange(male) %>%
  mutate(female = factor(female, levels=c("Val","Dhm","LHm","IV"))) %>%
  ggplot(aes(x = male, y = female)) +
  geom_raster(aes(fill = emmean)) +
  scale_fill_gradient(low = "grey90", high = "red") +
  labs(x = "Male Genotype", y = "Female Genotype", fill = "Mean\nDevelopment\nTime (Days)") +
  theme_bw() +
  scale_x_discrete(position = "top")
# Egg to adult viability
# lmer(viability ~ cross + eggs + day + (1|female))
# Add column of female

dev_time$percent_eclosed <- as.numeric(dev_time$percent_eclosed)

dev_time$cross <- paste0(dev_time$female,"_",dev_time$male)

plot_frame <- as.data.frame(aggregate(percent_eclosed ~ cross + day + male + female, data = dev_time, FUN = mean))

ggplot(data = plot_frame, aes(x = day, y = percent_eclosed, color = female)) +
  geom_point(size = 0) +
  geom_path(size = 1) +
  facet_grid(~male) +
  theme_classic()

lm_dev_time <- lm(percent_eclosed ~ male*female, data = dev_time)
summary(lm_dev_time)

