########################
### Spruce Analysis ###
########################

# Housekeeping: R project defaults to Git working directory, go up one
# level to get to the whole project's working directory (including data etc.)
setwd("../")

##### Data combination #####

#### Environmental data ####
# These data were taken every half hour, every day, in each plot.
# Goal: Combine all environmental data for all plots into a single dataset,
# limit it to our year of interest (2021)

# Vector of all files
env_filenames <- list.files("./Data/Raw/WEW_Complete_Environ_20220518/",full.names = T)

# Read in all as single list object
env_data_list <- lapply(env_filenames,read.csv)

# Extract plot numbers from file names
plot_nums <- stringr::str_extract(env_filenames,"(PLOT_\\d\\d).",group=1)
names(env_data_list) <- plot_nums

# Add column for plot number to each DF (will need for later combination into single DF)
for(i in 1:length(env_data_list)){
  env_data_list[[i]]$plot_number <- stringr::str_extract(names(env_data_list)[i],"\\d+")
}

## Check that all datasets have the same column names
# Column names of first dataset
plot4_env_cols <- colnames(env_data_list$PLOT_04)

# Compare to the rest
lapply(env_data_list[-1],FUN=function(x){
  if(all(colnames(x)==plot4_env_cols)){
    print("All columns equal")
  } else{
    print(which(colnames(x)!=plot4_env_cols))
  }
})
# All datasets have the same columns!

## Limit to just the year 2021 for all data (do this before combining into one DF
## to avoid having a ridiculously large DF)
env_data_list_2021 <- lapply(env_data_list,FUN=function(x){
  return(x[x$Year==2021,])
})

# Check that the data was really reduced to only 2021
lapply(env_data_list_2021,FUN=function(x)unique(x$Year))

# What's the dimensions of each dataset?
lapply(env_data_list_2021,FUN=dim)
# Each plot has 17,520 observations for the year 2021

## Combine into one dataframe
env_data <- data.table::rbindlist(env_data_list_2021)

# Check for irregularities in data
summary(env_data)
# Several columns have a few hundred NAs, but in a dataset with 210,240 observations
# this is neither surprising nor concerning. The deeper soil temperature measurements
# (TS_.100_A8, TS_.200_A9,TS_.100__C8, TS_.200__C9) have a couple thousand NAs -- maybe there were more
# equipment failures at lower depths, e.g. due to the higher pressure required to
# insert and remove equipment. It's interesting that this doesn't seem to have been
# the case in zone B, only in zones A and C. 
# PREC_6 has 43,401 NAs, and seems to consist primarily of 0s. PREC_6 is the total
# mm of rainfall at 6 m on the central tower, at 30 minute intervals, so it
# makes sense why it's mostly 0s (it's not raining more often than it rains).
# I'm not sure why there are so many NAs, though.

## Are the NAs mostly the same rows?

# Get the total number of NAs in each row
row_nas <- rowSums(is.na(env_data))
summary(row_nas)
hist(row_nas)
# Impossible to get a good look from this histogram. Narrow down to 0-20
hist(row_nas,breaks=c(seq(0,20,1),length(row_nas)),xlim=c(0,20))

# The majority of rows have 0 NAs. Remove all of these to get a better look
row_nas <- row_nas[row_nas>0]
summary(row_nas)
hist(row_nas)
hist(row_nas,breaks=c(seq(0,91,1)))
# Looks like most rows that have NAs only have 1, and almost all have fewer than 5.
sum(row_nas == 91); sum(row_nas>50); sum(row_nas>20)
# 269 rows have 91 NAs; 471 more than 50; 491 more than 20.
sum(row_nas < 20); sum(row_nas==1)
# By contrast, 46,993 rows have fewer than 20 NAs; of these, 40,672 rows have
# only 1 NA. So no, the NAs are not mostly the same rows. However, this is probably
# because PREC_6 has so many more NAs than any other column. What happens if we
# do the same thing with PREC_6 removed?

row_nas_nop <- rowSums(is.na(as.data.frame(env_data)[,c(which(colnames(env_data)!="PREC_6"))]))
summary(row_nas_nop)

row_nas_nop <- row_nas_nop[row_nas_nop>0]
summary(row_nas_nop)
hist(row_nas_nop,seq(1,max(row_nas_nop),1))
# No, the vast majority are still singlet NAs. No further action needed - there are
# more than enough datapoints in each variable to make scattered individual NAs
# non-important. Rows with more NAs are so rare as to be unimportant as well.

## Year.Fraction, WS_10, PAR_2, and plot_number should all be numeric.
# Year.Fraction
head(env_data$Year.Fraction) # no immediately visible reason to not be numeric
test <- as.numeric(env_data$Year.Fraction)
sum(is.na(test)) # no NAs ==> no incorrect conversions from character to numeric
env_data$Year.Fraction <- as.numeric(env_data$Year.Fraction)

# WS_10
head(env_data$WS_10,50)
test <- as.numeric(env_data$WS_10)
sum(is.na(test))==sum(is.na(env_data$WS_10))
# Which columns became NA when converting to numeric?
bad <- which(is.na(test) & !is.na(env_data$WS_10))
unique(env_data[bad,"WS_10"])
# They're all just empty. It's correct to convert them to NA so we can make the conversion.
env_data$WS_10 <- as.numeric(env_data$WS_10)

# PAR_2
head(env_data$PAR_2,50)
test <- as.numeric(env_data$PAR_2)
sum(is.na(test))==sum(is.na(env_data$PAR_2))
bad <- which(is.na(test) & !is.na(env_data$PAR_2))
unique(env_data[bad,"PAR_2"])
# It looks like entries in the 1000s were sometimes written with a comma, causing
# them to register as character and not convert to numeric. Remove all commas
# in the column and then convert to numeric
test<-env_data$PAR_2
test2<-gsub(",","",test)
test3 <- as.numeric(test2)
sum(is.na(test3))==sum(is.na(test)) # Confirm the commas were the only problem

env_data$PAR_2 <- gsub(",","",env_data$PAR_2)
env_data$PAR_2 <- as.numeric(env_data$PAR_2)

## Save this
write.csv(env_data,"./Data/Clean/environmental_data_2021.csv",row.names=F)

#### Soil temperature ####
soiltemp_raw <- read.csv("./Data/Figure3_Soil_Temperature_20170223.csv")
head(soiltemp_raw)
# Is it actually necessary to read this in separately? We can calculate it from the
# data we've already got, if we need it. Does it maybe make sense to calculate
# day-long averages for each variable regardless? Choose a depth of interest for
# some of them?

#### Phenology ####
# Read in raw data
phen_raw <- read.csv("./Data/Raw/SPRUCE_Ground_phenology_observations_2022.csv")
