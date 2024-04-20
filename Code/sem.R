### SEM for all SPRUCE variables ###

spruce <- read.csv("./Data/Clean/complete_combined_spruce_data.csv")
str(spruce)

# Expect bacteria and archaea copy num to be influenced by measured temp, depth,
# DOC, DN, gravimetric water content GWC
# Expect measured temp to be influenced by temp_experimental, co2_treatment, sample_date
# MBC, MBN influenced by bacteria/archae copy num, GWC

# First, investigate correlations
GGally::ggpairs(spruce[,c("Bacteria_copy_dry","Archaea_copy_dry","MBN","MBC")],
                lower=list(continuous="smooth_lm"))

# Remove outlier from MBN
GGally::ggpairs(spruce[-which(spruce$MBN>100),c("Bacteria_copy_dry","Archaea_copy_dry","MBN","MBC")],
                lower=list(continuous="smooth_lm"))

# for app: allow user to choose response variables, predictor variables -- 
# construct SEM from there. but DEFAULT to what we think the most reasonable 
# SEM is

library("piecewiseSEM")
library("lme4")

# Since there's only one plot per combination of Temp_experimental and CO2_treatment,
# can we ignore the plot number altogether and assume it's implicitly taken care
# of by the inclusion of these two variables? Or do we need to include it as a random
# effect, since each plot was sampled on two dates? If we do this, can we still include
# date as a fixed effect?

#### Testing model assumptions ####

## Dissolved organic carbon
DOC_lm <- lm(DOC_unfumigated_soil ~ depth2, data=spruce, na.action=na.omit)
par(mfrow=c(2,2))
plot(DOC_lm)
# Mostly fine, little bit of a right skew but I think it's good enough
# Check if the DOC is right-skewed
par(mfrow=c(1,2))
hist(spruce$DOC_unfumigated_soil)

# Does it look better log-transformed?
hist(log(spruce$DOC_unfumigated_soil))
# It does. Does transformation improve the Q-Q plot?
DOC_log_lm <- lm(log(DOC_unfumigated_soil) ~ depth2 + Sample_date, data=spruce, na.action=na.omit)
par(mfrow=c(2,2))
plot(DOC_log_lm)
# Somewhat, but I don't know that it's enough to be worth it.

## Dissolved nitrogen
DN_lm <- lm(DN_unfumigated_soil ~ depth2, data=spruce, na.action=na.omit)
par(mfrow=c(2,2))
plot(DN_lm)
# Same as DOC. Check histograms
par(mfrow=c(1,2))
hist(spruce$DN_unfumigated_soil)
hist(log(spruce$DN_unfumigated_soil))
# Here, log-transforming actually leads to more of a left-skew than otherwise

## Temperature
temp_lm <- lm(temp ~ Temp_experimental + CO2_treatment + depth2 + Sample_date,
              data=spruce, na.action=na.omit)
par(mfrow=c(2,2))
plot(temp_lm)
# Also fine, check histograms for completeness' sake
par(mfrow=c(1,2))
hist(spruce$temp)
hist(log(spruce$temp))
# Untransformed definitely better

## Gravimetric water content
GWC_lm <- lm(GWC ~ depth2 + Sample_date, data=spruce, na.action=na.omit)
par(mfrow=c(2,2))
plot(GWC_lm)
# Mostly fine but definitely a skew in the Q-Q plot
par(mfrow=c(1,2))
hist(spruce$GWC)
hist(log(spruce$GWC))
# I'm not sure that log-transformming is the move here. Does it change the Q-Q plot?
GWC_log_lm <- lm(log(GWC) ~ depth2 + Sample_date, data=spruce, na.action=na.omit)
par(mfrow=c(2,2))
plot(GWC_log_lm)
# It actually skews it worse. Stick to the untransformed data.

## Bacteria copy number
BCN_lm <- lm(Bacteria_copy_dry ~ DOC_unfumigated_soil + DN_unfumigated_soil + temp,
             data=spruce, na.action=na.omit)
par(mfrow=c(2,2))
plot(BCN_lm)
# Heavy skew
par(mfrow=c(1,2))
hist(spruce$Bacteria_copy_dry)
hist(log(spruce$Bacteria_copy_dry))
# log-transformation seems to improve things considerably
BCN_log_lm <- lm(log(Bacteria_copy_dry) ~ DOC_unfumigated_soil + DN_unfumigated_soil + temp,
                 data=spruce, na.action=na.omit)
par(mfrow=c(2,2))
plot(BCN_log_lm)
# This is much more respectable. We may have to remove a couple outliers, though.
# Check points 135 and 108: 
spruce[c(108,135),]
# Keep these in mind while going through the rest.

## Archaea copy number
ACN_lm <- lm(Archaea_copy_dry ~ DOC_unfumigated_soil + DN_unfumigated_soil + temp,
             data=spruce, na.action=na.omit)
par(mfrow=c(2,2))
plot(ACN_lm)
# Same skewed behavior as bacteria
par(mfrow=c(1,2))
hist(spruce$Archaea_copy_dry)
hist(log(spruce$Archaea_copy_dry))
# log-transformation seems to improve things considerably
ACN_log_lm <- lm(log(Archaea_copy_dry) ~ DOC_unfumigated_soil + DN_unfumigated_soil + temp,
                 data=spruce, na.action=na.omit)
par(mfrow=c(2,2))
plot(ACN_log_lm)
# This is also better, and the points that were outliers for bacteria aren't here.

## Microbial biomass nitrogen
MBN_lm <- lm(MBN ~ DN_unfumigated_soil + temp + GWC,
             data=spruce, na.action=na.omit)
plot(MBN_lm)
# Point 139 is a huge outlier, we noticed it in the covariance testing as well.
# Remove it by setting it to NA but keeping the row intact (is this legit?)
spruce_no_out <- spruce
spruce_no_out[139,"MBN"] <- NA
MBN_lm <- lm(MBN ~ DN_unfumigated_soil + temp + GWC,
             data=spruce_no_out, na.action=na.omit)
plot(MBN_lm)
# Looking better, but very skewed.
par(mfrow=c(1,2))
hist(spruce_no_out$MBN)
hist(log(spruce_no_out$MBN))
# Looks nicer, let's test a new model.
MBN_log_lm <- lm(log(MBN) ~ DN_unfumigated_soil + temp + GWC,
                 data=spruce_no_out, na.action=na.omit)
par(mfrow=c(2,2))
plot(MBN_log_lm)
# Very nice!

## Microbial biomass carbon
MBC_lm <- lm(MBC ~ DOC_unfumigated_soil + temp + GWC,
             data=spruce, na.action=na.omit)
plot(MBC_lm)
# Also somewhat skewed, not as bad as the others but probably worth log-transforming
par(mfrow=c(1,2))
hist(spruce$MBC)
hist(log(spruce$MBC))
# It just kind of shifts to left-skewed... still, see if this helps the model
MBC_log_lm <- lm(log(MBC) ~ DOC_unfumigated_soil + temp + GWC,
                     data=spruce, na.action=na.omit)
par(mfrow=c(2,2))
plot(MBC_log_lm)
# Honestly not sure which is better! Maybe go with log-transformed for consistency. 
# Ask for feedback.

#### Piecewise structural Equation Modeling ####
library(piecewiseSEM)

# Take the preferred models determined above, but use the outlier-removed
# dataset for all so that they're all on the same dataset.
# Package recommends removing all rows with NAs; consider doing this...
spruce_psem <- psem(
  
  # Intermediate layer: DOC, DN, temperature, GWC
  lm(DOC_unfumigated_soil ~ depth2, data=spruce_no_out, na.action=na.omit),
  lm(DN_unfumigated_soil ~ depth2, data=spruce_no_out, na.action=na.omit),
  lm(temp ~ Temp_experimental + CO2_treatment + depth2 + Sample_date,
     data=spruce_no_out, na.action=na.omit),
  lm(GWC ~ depth2 + Sample_date, data=spruce_no_out, na.action=na.omit),
  
  # Predicted variables layer: Bacteria and Archaea copy numbers, MBN, MBC
  lm(log(Bacteria_copy_dry) ~ DOC_unfumigated_soil + DN_unfumigated_soil + temp,
     data=spruce_no_out, na.action=na.omit),
  lm(log(Archaea_copy_dry) ~ DOC_unfumigated_soil + DN_unfumigated_soil + temp,
     data=spruce_no_out, na.action=na.omit),
  lm(log(MBN) ~ DN_unfumigated_soil + temp + GWC,
     data=spruce_no_out, na.action=na.omit),
  lm(log(MBC) ~ DOC_unfumigated_soil + temp + GWC,
     data=spruce, na.action=na.omit)
  
)

summary(spruce_psem, .progressBar = FALSE)
