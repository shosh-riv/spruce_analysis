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
                lower=list(continuous="smooth_loess"))
