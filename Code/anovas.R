#########################################
### Calculate one ANOVA and Tukey HSD ###
### for each response variable        ###
#########################################

spruce <- read.csv("./Data/Clean/complete_combined_spruce_data.csv")
colnames(spruce)

# Replace dash with underscore in depth ranges so the names don't confuse the modeling
spruce$depth2 <- gsub("-","_",spruce$depth2)

# Response variables
response_vars <- c("Bacteria Copy Number"="Bacteria_copy_dry",
                   "Archaeal Copy Number"="Archaea_copy_dry",
                   "Gravimetric Water Content"="GWC",
                   "Dissolved Nitrogen"="DN_unfumigated_soil",
                   "Microbial Biomass Nitrogen"="MBN",
                   "Dissolved Organic Carbon"="DOC_unfumigated_soil",
                   "Microbial Biomass Carbon"="MBC",
                   "Soil Temperature"="temp")

# Explanatory variables
explanatory_vars <- c("Experimental Temperature"="Temp_experimental",
                      "CO2 Treatment"="CO2_treatment","Depth"="depth2")

# How many possible tests can be run?
length(response_vars) * 7 # 7 because each response can be run against one of the
# 3 explanatory variables, one of 3 possible combinations of 2 of them, or all 3
# explanatory variables at once. 3 + 3 + 1 possible combinations = 7 possible
# ANOVAs for each response variable.

# Select variable for ANOVA
lhs <- "Bacteria_copy_dry"

# Select explantory variable(s)
x_axis <- "depth2"
fill_var <- "Temp_experimental"

# Create formula
aov_formula <- paste0(lhs,"~",paste(x_axis,fill_var,sep ="*"))

# Run ANOVA
aov_object <- aov(formula(aov_formula), data=spruce)

# Run TukeyHSD
tukey <- TukeyHSD(aov_object)

# Extract p values
tukey_p <- lapply(tukey,FUN=function(x){
  y <- x[,4]
  names(y) <- rownames(x)
  return(y)
})

# Get letters for plotting
library(multcompView)
if(!is.null(tukey[[3]])){
  # If the interaction was included in the model, check if there were any significant interactions
  # Which interactions were significant?
  int <- which(tukey_p[[3]]<0.05)
}

# If there were any significant interactions, print them as a table
if(length(int) > 0){
  interaction_table <- tukey[[3]][int,]
  interaction_table <- knitr::kable(interaction_table,
                                    colnames=c("Comparison","Difference","Lower","Upper","Adjusted p-value"))
}

# Get letters for boxplot

x_letters <- multcompLetters(tukey_p[[x_axis]])$Letters
x_letters <- data.frame(level = names(x_letters),letter = x_letters)

fill_letters <- multcompLetters(tukey_p[[fill_var]])

# Boxplot
library(ggplot2)
ggplot(data=spruce)+
  geom_boxplot(aes(x=.data[[x_axis]],y=.data[[lhs]],fill=.data[[fill_var]]))+
  theme_classic() +
  geom_text(data=x_letters,aes(x=level,y=(max(spruce[[lhs]])+max(spruce[[lhs]])/10)),label=x_letters$letter)

