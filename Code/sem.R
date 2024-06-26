### SEM for all SPRUCE variables ###

spruce_orig <- read.csv("./Data/Clean/complete_combined_spruce_data_noprecip.csv")
# str(spruce_orig)
# 
# # Expect bacteria and archaea copy num to be influenced by measured temp, depth,
# # DOC, DN, gravimetric water content GWC
# # Expect measured temp to be influenced by temp_experimental, co2_treatment, sample_date
# # MBC, MBN influenced by bacteria/archae copy num, GWC
# 
# ### Add precipitation back into this dataframe
# env <- read.csv("./Data/Clean/dailymean_environmental_data_2021.csv")
# 
# # Convert day of year to date
# env$Sample_date <- as.character(as.Date(env$Day_of_Year,origin="2020-12-31"))
# 
# # Grab precip data
# spruce <- merge(spruce_orig,env[,c("Sample_date","Plot","PREC_6")],
#                  by=c("Sample_date","Plot"),all.x=T,all.y=F)
# 
# # Due to "heterogeneity between the two replicate samples" used to calculate MBN
# # and MBC, there are a few negative values that naturally become NAs when log-transformed.
# # Since negative values are nonsensical, we'll change these values to 0 
# # to reflect that they have a very small microbial biomass. We'll then add 1
# # to DN and DOC for the log transformation.
# spruce[which(spruce$MBN<0),"MBN"] <- 0
# spruce[which(spruce$MBC<0),"MBC"] <- 0
# 
# # Set MBN outlier to NA
# spruce[which(spruce$MBN>100),"MBN"] <- NA
# 
# # Remove extra column
# spruce$X <- NULL
# 
# # Save this new dataset
# write.csv(spruce,"./Data/Clean/complete_combined_spruce_data.csv",row.names=F)

# Read in this dataset
spruce <- read.csv("./Data/Clean/complete_combined_spruce_data.csv")

#### Covariance matrix ####
# First, investigate correlations
GGally::ggpairs(spruce[,c("Bacteria_copy_dry","Archaea_copy_dry","MBN","MBC")],
                lower=list(continuous="smooth_lm"))

# Full covariance matrix, excluding non-numeric treatments/depth
GGally::ggpairs(subset(spruce, MBN < 100, 
                       select = -c(Sample_date,Plot,datesitedepth,Temp_experimental,CO2_treatment,depth2)),
                lower=list(continuous="smooth_lm")) + theme_light()
ggsave("C:/Users/linne/Documents/School/Cornell/MultivariateAnalysis/FinalProject/spruce_analysis/Plots/cov_matrix_orig.png",
       width=5000,height=3000,units="px")

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
### For all, use the dataset without the ambient temperature treatment

## Dissolved organic carbon
DOC_lm <- lm(DOC_unfumigated_soil ~ depth2 + CO2_treatment + depth2, data=spruce_noAmb, na.action=na.omit)
par(mfrow=c(2,2))
plot(DOC_lm)
# Mostly fine, little bit of a right skew but I think it's good enough
# Check if the DOC is right-skewed
par(mfrow=c(1,2))
hist(spruce_noAmb$DOC_unfumigated_soil)

# Does it look better log-transformed?
hist(log(spruce_noAmb$DOC_unfumigated_soil))
# It does. Does transformation improve the Q-Q plot?
DOC_log_lm <- lm(log(DOC_unfumigated_soil) ~ depth2 + CO2_treatment + depth2, data=spruce_noAmb, na.action=na.omit)
par(mfrow=c(2,2))
plot(DOC_log_lm)
# Somewhat, but I don't know that it's enough to be worth it.

## Dissolved nitrogen
DN_lm <- lm(DN_unfumigated_soil ~ depth2 + Temp_experimental, data=spruce_noAmb, na.action=na.omit)
par(mfrow=c(2,2))
plot(DN_lm)
# Same as DOC. Check histograms
par(mfrow=c(1,2))
hist(spruce_noAmb$DN_unfumigated_soil)
hist(log(spruce_noAmb$DN_unfumigated_soil))
# Here, log-transforming actually leads to more of a left-skew than otherwise

## Temperature
temp_lm <- lm(temp ~ Temp_experimental + depth2,
              data=spruce_noAmb, na.action=na.omit)
par(mfrow=c(2,2))
plot(temp_lm)
# Also fine, check histograms for completeness' sake
par(mfrow=c(1,2))
hist(spruce_noAmb$temp)
hist(log(spruce_noAmb$temp))
# Untransformed definitely better

## Gravimetric water content
GWC_lm <- lm(GWC ~ depth2 + PREC_6, data=spruce_noAmb, na.action=na.omit)
par(mfrow=c(2,2))
plot(GWC_lm)
# Mostly fine but definitely a skew in the Q-Q plot
par(mfrow=c(1,2))
hist(spruce_noAmb$GWC)
hist(log(spruce_noAmb$GWC))
# I'm not sure that log-transformming is the move here. Does it change the Q-Q plot?
GWC_log_lm <- lm(log(GWC) ~ depth2 + PREC_6, data=spruce_noAmb, na.action=na.omit)
par(mfrow=c(2,2))
plot(GWC_log_lm)
# It actually skews it worse. Stick to the untransformed data.

## Bacteria copy number
BCN_lm <- lm(Bacteria_copy_dry ~ DOC_unfumigated_soil + DN_unfumigated_soil + temp + GWC,
             data=spruce_noAmb, na.action=na.omit)
par(mfrow=c(2,2))
plot(BCN_lm)
# Heavy skew
par(mfrow=c(1,2))
hist(spruce_noAmb$Bacteria_copy_dry)
hist(log(spruce_noAmb$Bacteria_copy_dry))
# log-transformation seems to improve things considerably
BCN_log_lm <- lm(log(Bacteria_copy_dry) ~ DOC_unfumigated_soil + DN_unfumigated_soil + temp + GWC,
                 data=spruce_noAmb, na.action=na.omit)
par(mfrow=c(2,2))
plot(BCN_log_lm)
# This is much more respectable. We may have to remove a couple outliers, though.
# Check points 135 and 108: 
spruce_noAmb[c(108,135),]
# Keep these in mind while going through the rest.

## Archaea copy number
ACN_lm <- lm(Archaea_copy_dry ~ DOC_unfumigated_soil + DN_unfumigated_soil + temp + GWC,
             data=spruce_noAmb, na.action=na.omit)
par(mfrow=c(2,2))
plot(ACN_lm)
# Same skewed behavior as bacteria
par(mfrow=c(1,2))
hist(spruce_noAmb$Archaea_copy_dry)
hist(log(spruce_noAmb$Archaea_copy_dry))
# log-transformation seems to improve things considerably
ACN_log_lm <- lm(log(Archaea_copy_dry) ~ DOC_unfumigated_soil + DN_unfumigated_soil + temp + GWC,
                 data=spruce_noAmb, na.action=na.omit)
par(mfrow=c(2,2))
plot(ACN_log_lm)
# This is also better, and the points that were outliers for bacteria aren't here.

## Microbial biomass nitrogen
MBN_lm <- lm(MBN ~ DN_unfumigated_soil + temp + GWC,
             data=spruce_noAmb, na.action=na.omit)
plot(MBN_lm)

# Skewed, log-transform (add 1 because there are 0s in the dataset)
par(mfrow=c(1,2))
hist(spruce_noAmb$MBN)
hist(log(spruce_noAmb$MBN+1))
# This actually doesn't help matters all that much...
MBN_log_lm <- lm(log(MBN+1) ~ DN_unfumigated_soil + temp + GWC,
                 data=spruce_noAmb, na.action=na.omit)
par(mfrow=c(2,2))
plot(MBN_log_lm)
# A little nicer

## Microbial biomass carbon
MBC_lm <- lm(MBC ~ DOC_unfumigated_soil + temp + GWC,
             data=spruce_noAmb, na.action=na.omit)
plot(MBC_lm)
# Also somewhat skewed, not as bad as the others but probably worth log-transforming
par(mfrow=c(1,2))
hist(spruce_noAmb$MBC)
hist(log(spruce_noAmb$MBC+1))
# It just kind of shifts to left-skewed... still, see if this helps the model
MBC_log_lm <- lm(log(MBC+1) ~ DOC_unfumigated_soil + temp + GWC,
                 data=spruce_noAmb, na.action=na.omit)
par(mfrow=c(2,2))
plot(MBC_log_lm)
# Honestly not sure which is better! Maybe go with log-transformed for consistency. 
# Ask for feedback.

#### Scale, log-transform data ####
# Log-transform the variables that need log-transformation as determined above.
# Also scale the entire dataset.

# The ambient temperature treatment is very different to the 0.00 treatment. Since
# the ambient plots were treated differently to the rest of the plots and may
# not be comparable, so we'll remove them.
spruce_noAmb <- spruce[-which(spruce$Temp_experimental=="Amb"),]

# Make Temp_experimental numeric now that we've removed ambient
spruce_noAmb$Temp_experimental <- as.numeric(spruce_noAmb$Temp_experimental)

# Remove all NA rows, to make sure that the SEM is fitting the same dataset for each model
spruce_noAmb <- na.omit(spruce_noAmb)

# Initialize dataframe for transformation/scaling
spruce_log_scale <- spruce_noAmb

### Log-transformation

# Variables to transform
to_transform <- c("Bacteria_copy_dry","Archaea_copy_dry","MBN","MBC")

# Check for 0s
apply(spruce_log_scale[,to_transform],2,FUN=min,na.rm=T)

# Add 1 to MBC and MBN because they include 0
spruce_log_scale$MBN <- spruce_log_scale$MBN + 1
spruce_log_scale$MBC <- spruce_log_scale$MBC + 1

# Transform
spruce_log_scale[,to_transform] <- apply(spruce_log_scale[,to_transform],2,log)

### Scale all numerical values
# Vector of numeric columns
to_scale <- unlist(lapply(spruce_log_scale, is.numeric), use.names = FALSE)

# Double check these are only columns we want
colnames(spruce_log_scale)[to_scale]

# Make plot into a character
spruce_log_scale$Plot <- as.character(spruce_log_scale$Plot)

# Redo the vector of numeric columns
to_scale <- unlist(lapply(spruce_log_scale, is.numeric), use.names = FALSE)
colnames(spruce_log_scale)[to_scale]

# Scale columns
spruce_log_scale[,to_scale] <- scale(spruce_log_scale[,to_scale])

### Covariance matrix with log-transformed scaled data

# Column names for better legibility
vars <- c("Exp.\ntemp.","Bacteria","Archaea","GWC","DN","MBN","DOC","MBC",
          "Soil\ntemp.","Precip.")

# Covariance matrix
library(GGally)
log_cov_plot <- GGally::ggpairs(subset(spruce_log_scale,
                       select = -c(Sample_date,Plot,datesitedepth,CO2_treatment,depth2)),
                lower=list(continuous="smooth_lm"),
                columnLabels = vars,
                upper = list(continuous = wrap("cor", size = 8.5))) + 
  theme_light(base_size = 18) +
  theme(strip.text.x = element_text(size = 20),
        strip.text.y = element_text(size = 20))
ggsave("C:/Users/linne/Documents/School/Cornell/MultivariateAnalysis/FinalProject/spruce_analysis/Plots/cov_matrix_logscale.png",
       width=6300,height=3300,units="px")

### Redo covariance matrix with comparable original data
cov_plot <- GGally::ggpairs(subset(spruce_noAmb,
                                   select = -c(Sample_date,Plot,datesitedepth,CO2_treatment,depth2)),
                            lower=list(continuous="smooth_lm"),
                            columnLabels = vars,
                            upper = list(continuous = wrap("cor", size = 8.5))) + 
  theme_light(base_size = 18) +
  theme(strip.text.x = element_text(size = 20),
        strip.text.y = element_text(size = 20))
ggsave("C:/Users/linne/Documents/School/Cornell/MultivariateAnalysis/FinalProject/spruce_analysis/Plots/cov_matrix_noAmb.png",
       width=6300,height=3300,units="px")

# Combine both into one figure for paper
library(cowplot)
fig_covplot <- cowplot::plot_grid(ggmatrix_gtable(cov_plot),
                   ggmatrix_gtable(log_cov_plot),
                   nrow=2,labels=c("(A)","(B)"),label_size=18)
ggsave("C:/Users/linne/Documents/School/Cornell/MultivariateAnalysis/FinalProject/spruce_analysis/Plots/figure_covariance_matrix.png",
       fig_covplot,width=7000,height=7200,units="px")


#### Piecewise structural Equation Modeling ####
library(piecewiseSEM)

# Take the preferred models determined above, but use the outlier-removed
# dataset for all so that they're all on the same dataset.
spruce_psem <- psem(
  
  # Intermediate layer: DOC, DN, temperature, GWC
  lm(DOC_unfumigated_soil ~ depth2 + Temp_experimental + CO2_treatment, data=spruce_log_scale, na.action=na.omit),
  lm(DN_unfumigated_soil ~ depth2 + Temp_experimental, data=spruce_log_scale, na.action=na.omit),
  lm(temp ~ Temp_experimental + depth2,
     data=spruce_log_scale, na.action=na.omit),
  lm(GWC ~ depth2 + PREC_6, data=spruce_log_scale, na.action=na.omit),
  
  # Predicted variables layer: Bacteria and Archaea copy numbers, MBN, MBC
  lm(Bacteria_copy_dry ~ DOC_unfumigated_soil + DN_unfumigated_soil + temp + GWC,
     data=spruce_log_scale, na.action=na.omit),
  lm(Archaea_copy_dry ~ DOC_unfumigated_soil + DN_unfumigated_soil + temp + GWC,
     data=spruce_log_scale, na.action=na.omit),
  lm(MBN ~ DN_unfumigated_soil + temp + GWC,
     data=spruce_log_scale, na.action=na.omit),
  lm(MBC ~ DOC_unfumigated_soil + temp + GWC,
     data=spruce_log_scale, na.action=na.omit)
  
)

summary(spruce_psem, .progressBar = FALSE)

## Add some parameters in based on d-sep tests (currently can afford to add 12 relationships)

# Add depth, experimental temperature to MBN, MBC (8 relationships left to add)
spruce_psem <- update(spruce_psem,
       MBN ~ DN_unfumigated_soil + temp + GWC + depth2 + Temp_experimental,
       MBC ~ DOC_unfumigated_soil + temp + GWC + depth2 + Temp_experimental)
summary(spruce_psem)
# R-squared for MBN and MBC much higher now!

# Add temperature to DOC, DN
spruce_psem <- update(spruce_psem,
                      DOC_unfumigated_soil ~ depth2 + Temp_experimental + CO2_treatment + temp,
                      DN_unfumigated_soil ~ depth2 + Temp_experimental + temp)
summary(spruce_psem)

#### Plot PSEM
source("./Code/plot_psem.R")
library(DiagrammeR)
plot.psem(spruce_psem)

## What happens if we make all variables numeric? We need to to use the plotting
## function anyway! (Can't deal with ANOVAs)

#### Make categorical variables numeric ####
spruce_num <- spruce_noAmb
### Temperature treatment
unique(spruce_num$Temp_experimental)
spruce_num$Temp_experimental <- as.numeric(spruce_num$Temp_experimental)

### CO2 treatment
unique(spruce_num$CO2_treatment)

# A is ambient, and E is elevated 500 ppm above ambient. Change these to 0 and 1.
spruce_num[which(spruce_num$CO2_treatment=="A"),"CO2_treatment"] <- 0
spruce_num[which(spruce_num$CO2_treatment=="E"),"CO2_treatment"] <- 1
spruce_num$CO2_treatment <- as.numeric(spruce_num$CO2_treatment)

### Depth
unique(spruce_num$depth2)

# They sampled from the top of each range, not throughout the whole range, so we can
# take this as the smaller number for each range given.
spruce_num$depth2 <- stringr::str_extract(spruce_num$depth2,"(\\d+)-",group=1)
spruce_num$depth2 <- as.numeric(spruce_num$depth2)

### Date
unique(spruce_num$Sample_date)

# We only have two sampling dates, so this makes sense to do as a binary.
# 0 = 6/23/21, 1 = 8/24/21
spruce_num[which(spruce_num$Sample_date=="2021-06-23"),"Sample_date"] <- 0
spruce_num[which(spruce_num$Sample_date=="2021-08-24"),"Sample_date"] <- 1
spruce_num$Sample_date <- as.numeric(spruce_num$Sample_date)

#### Log-transform and scale numeric data ####
spruce_num_log_scale <- spruce_num

# Variables to transform
to_transform <- c("Bacteria_copy_dry","Archaea_copy_dry","MBN","MBC")

# Check for 0s
apply(spruce_num_log_scale[,to_transform],2,FUN=min,na.rm=T)

# Add 1 to MBC and MBN because they include 0
spruce_num_log_scale$MBN <- spruce_num_log_scale$MBN + 1
spruce_num_log_scale$MBC <- spruce_num_log_scale$MBC + 1

# Log-transformation
spruce_num_log_scale[,to_transform] <- log(spruce_num_log_scale[,to_transform])

### Scale all numerical values
# Vector of numeric columns
to_scale <- unlist(lapply(spruce_num_log_scale, is.numeric), use.names = FALSE)

# Double check these are only columns we want
colnames(spruce_num_log_scale)[to_scale]

# Make plot into a character
spruce_num_log_scale$Plot <- as.character(spruce_num_log_scale$Plot)

# Redo the vector of numeric columns
to_scale <- unlist(lapply(spruce_num_log_scale, is.numeric), use.names = FALSE)
colnames(spruce_num_log_scale)[to_scale]

# Scale columns
spruce_num_log_scale[,to_scale] <- scale(spruce_num_log_scale[,to_scale])

# Remove NA rows
spruce_num_log_scale_noNA <- na.omit(spruce_num_log_scale)

## Try fitting the same model with numerical data
spruce_num_psem <- psem(
  
  # Intermediate layer: DOC, DN, temperature, GWC
  lm(DOC_unfumigated_soil ~ depth2 + Temp_experimental + CO2_treatment + temp, data=spruce_num_log_scale, na.action=na.omit),
  lm(DN_unfumigated_soil ~ depth2 + Temp_experimental + temp, data=spruce_num_log_scale, na.action=na.omit),
  lm(temp ~ Temp_experimental + depth2,
     data=spruce_num_log_scale, na.action=na.omit),
  lm(GWC ~ depth2 + PREC_6, data=spruce_num_log_scale, na.action=na.omit),
  
  # Predicted variables layer: Bacteria and Archaea copy numbers, MBN, MBC
  lm(Bacteria_copy_dry ~ DOC_unfumigated_soil + DN_unfumigated_soil + temp + GWC,
     data=spruce_num_log_scale, na.action=na.omit),
  lm(Archaea_copy_dry ~ DOC_unfumigated_soil + DN_unfumigated_soil + temp + GWC,
     data=spruce_num_log_scale, na.action=na.omit),
  lm(MBN ~ DN_unfumigated_soil + temp + GWC + depth2 + Temp_experimental,
     data=spruce_num_log_scale, na.action=na.omit),
  lm(MBC ~ DOC_unfumigated_soil + temp + GWC + depth2 + Temp_experimental,
     data=spruce_num_log_scale, na.action=na.omit)
  
)

summary(spruce_num_psem, .progressBar = FALSE)
summary(spruce_num_psem, .progressBar = FALSE)$R2
summary(spruce_psem, .progressBar = F)$R2

# This is worse than the categorical SEM for explaining all endogenous variables.
# Compare AICs:
AIC_psem(spruce_psem)
AIC_psem(spruce_num_psem)
# Categorical is better than numeric. Still, try using it to plot.

diagram <- plot.psem(spruce_num_psem,digits=1, return=T, layout="tree")
plot.psem(spruce_num_psem,layout="tree")


##### Plotting #####
#### Hypothesis path diagram ###
# Original model run with numeric data
spruce_psem <- psem(
  
  # Intermediate layer: DOC, DN, temperature, GWC
  lm(DOC_unfumigated_soil ~ depth2 + Temp_experimental + CO2_treatment, data=spruce_num_log_scale, na.action=na.omit),
  lm(DN_unfumigated_soil ~ depth2 + Temp_experimental, data=spruce_num_log_scale, na.action=na.omit),
  lm(temp ~ Temp_experimental + depth2,
     data=spruce_num_log_scale, na.action=na.omit),
  lm(GWC ~ depth2 + PREC_6, data=spruce_num_log_scale, na.action=na.omit),
  
  # Predicted variables layer: Bacteria and Archaea copy numbers, MBN, MBC
  lm(Bacteria_copy_dry ~ DOC_unfumigated_soil + DN_unfumigated_soil + temp + GWC,
     data=spruce_num_log_scale, na.action=na.omit),
  lm(Archaea_copy_dry ~ DOC_unfumigated_soil + DN_unfumigated_soil + temp + GWC,
     data=spruce_num_log_scale, na.action=na.omit),
  lm(MBN ~ DN_unfumigated_soil + temp + GWC,
     data=spruce_num_log_scale, na.action=na.omit),
  lm(MBC ~ DOC_unfumigated_soil + temp + GWC,
     data=spruce_num_log_scale, na.action=na.omit)
  
)


diagram <- plot.psem(spruce_num_psem,digits=0, return=T, layout="tree")
plot.psem(spruce_psem,layout="tree",digits=0)

## Modify the nodes DF
# Readable labels
diagram$nodes_df$label <- c("DOC","DN","Soil\ntemp.","GWC","Bacteria",
                            "Archaea","MBN","MBC","Depth","Exp.\ntemp.",
                            "CO2\nlevel","Precip.")

# Different shapes for exogenous and endogenous variables
diagram$nodes_df$shape <- c(rep("oval",times=8),rep("rectangle",times=4))

# Change layout
diagram$global_attrs[1,2] <- "dot"

# Make nodes a little wider
diagram$global_attrs[8,2] <- 0.6

# Make all arrows solid
diagram$edges_df$style <- "solid"

# Change colors of edges
diagram$edges_df$color <- RColorBrewer::brewer.pal(n = 8, name = "Dark2")[diagram$edges_df$to]

# Change colors of boxes
diagram$nodes_df$color <- c(RColorBrewer::brewer.pal(n = 8, name = "Dark2"),rep("black",times=4))

# Remove edge numbers
diagram$edges_df$label <- NA

# Check how it looks
render_graph(diagram)

# Save diagram
export_graph(diagram,"./Plots/hypothesis_sem.png",width = 1200, height = 800)

#### Final path diagram ####
## Modify the nodes DF
# Readable labels
diagram$nodes_df$label <- c("DOC","DN","Soil\ntemp.","GWC","Bacteria",
                            "Archaea","MBN","MBC","Depth","Exp.\ntemp.",
                            "CO2\nlevel","Precip.")

# Different shapes for exogenous and endogenous variables
diagram$nodes_df$shape <- c(rep("oval",times=8),rep("rectangle",times=4))

# Change layout
diagram$global_attrs[1,2] <- "dot"

# Make nodes a little wider
diagram$global_attrs[8,2] <- 0.6

# Convert node ID numbers to descriptions so the mouseover will be descriptive
node_ids <- diagram$nodes_df$label
names(node_ids) <-  diagram$nodes_df$id
diagram$nodes_df$id <-diagram$nodes_df$label

diagram$edges_df$from <- node_ids[diagram$edges_df$from]
diagram$edges_df$to <- node_ids[diagram$edges_df$to]

# See how it all looks
render_graph(diagram)

# Convert to graphViz/DOT format to put into shiny
dot <- generate_dot(diagram)

# Save DOT output for use in graphViz/Shiny
writeLines(dot,"./Plots/dot_psem.txt")

# Also save the PSEM summary
saveRDS(summary(spruce_num_psem),file = "./Analyses/psem_numeric_summary.RDS")


## Model identification
# Model still doesn't fit. Check the t-rule. The model currently includes
# 12 variables; according to the t-rule, t <= [n(n+1)]/2, where n is the number
# of unique cells in the variance-covariance matrix (the diagonal and the lower
# off-diagonal; see: https://jslefche.github.io/sem_book/global-estimation.html#model-identifiability)
# With 12 variables, n=78. 
sum(1:12)
(78*(78+1))/2

# We have 34 relationships and 12 parameters to estimate, so t = 34+12 = 46. This
# is valid, assuming I'm using the t-rule right. We move on to checking sample size.

# "The basic rule-of-thumb is that the level of replication should be at least 5 times the 
# number of estimated coefficients (not error variances or other correlations). So in our 
# previous path model, we are estimating two relationships, so we require at least n=10 to
# fit that model. However, this value is a lower limit: ideally, replication is 5-20x 
# the number of estimated parameters. "
# We're estimating 34 relationships with our current model structure.
34*5
# This is 170, which is fewer than our 177 samples.
177/5
# At the most, we can include 35 relationships; ideally, we include fewer than this.


## What if we do one SEM per response variable to bring down the total number of relationships?

# Bacteria copy number
bcn_psem <- psem(
  
  # Intermediate layer: DOC, DN, temperature, GWC
  lm(DOC_unfumigated_soil ~ depth2 + CO2_treatment + GWC, data=spruce_log_scale, na.action=na.omit),
  lm(DN_unfumigated_soil ~ depth2 + CO2_treatment + GWC, data=spruce_log_scale, na.action=na.omit),
  lm(temp ~ Temp_experimental + CO2_treatment + depth2,
     data=spruce_log_scale, na.action=na.omit),
  lm(GWC ~ depth2 + temp + PREC_6, data=spruce_log_scale, na.action=na.omit),
  
  # Predicted variables layer: Bacteria and Archaea copy numbers, MBN, MBC
  lm(Bacteria_copy_dry ~ DOC_unfumigated_soil + DN_unfumigated_soil + temp + GWC + depth2 + Temp_experimental,
     data=spruce_log_scale, na.action=na.omit)#,
  # lm(Archaea_copy_dry ~ DOC_unfumigated_soil + DN_unfumigated_soil + temp + GWC + depth2 + Temp_experimental,
  #    data=spruce_log_scale, na.action=na.omit),
  # lm(MBN ~ DN_unfumigated_soil + temp + GWC + depth2,
  #    data=spruce_log_scale, na.action=na.omit),
  # lm(MBC ~ DOC_unfumigated_soil + temp + GWC + depth2,
  #    data=spruce_log_scale, na.action=na.omit)
  
)

summary(bcn_psem, .progressBar = FALSE)
# Still not an identified model but WAY better than before

# Archaea copy number
acn_psem <- psem(
  
  # Intermediate layer: DOC, DN, temperature, GWC
  lm(DOC_unfumigated_soil ~ depth2 + CO2_treatment + GWC + temp, data=spruce_log_scale_noNA, na.action=na.omit),
  lm(DN_unfumigated_soil ~ depth2 + CO2_treatment + GWC + temp, data=spruce_log_scale_noNA, na.action=na.omit),
  lm(temp ~ Temp_experimental + CO2_treatment + depth2,
     data=spruce_log_scale_noNA, na.action=na.omit),
  lm(GWC ~ depth2 + temp + PREC_6, data=spruce_log_scale_noNA, na.action=na.omit),
  
  # Predicted variables layer: Bacteria and Archaea copy numbers, MBN, MBC
  # lm(Bacteria_copy_dry ~ DOC_unfumigated_soil + DN_unfumigated_soil + temp + GWC + depth2 + Temp_experimental,
  #    data=spruce_log_scale_noNA, na.action=na.omit)#,
  lm(Archaea_copy_dry ~ DOC_unfumigated_soil + DN_unfumigated_soil + temp + GWC + depth2 + Temp_experimental + CO2_treatment,
     data=spruce_log_scale_noNA, na.action=na.omit)#,
  # lm(MBN ~ DN_unfumigated_soil + temp + GWC + depth2,
  #    data=spruce_log_scale_noNA, na.action=na.omit),
  # lm(MBC ~ DOC_unfumigated_soil + temp + GWC + depth2,
  #    data=spruce_log_scale_noNA, na.action=na.omit)
  
)
summary(acn_psem, .progressBar = FALSE)

#### Compare numeric vs categorical model fits ####

