#####################################
### SPRUCE data visualization app ###
### Autumn Pereira, Linnea Smith  ###
### Cornell U, Spring 2024        ###
#####################################
library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(ggplot2)
library(ggfortify)
library(plotly)
library(multcompView)
library(DiagrammeR)

# Read in original (non-transformed, non-scaled) data
d_orig <- read.csv("./Data/complete_combined_spruce_data.csv")
# Replace dash with underscore in depth ranges so the names don't confuse the modeling
d_orig$depth2 <- gsub("-","_",d_orig$depth2)
# Relevel depth variable
d_orig$depth2 <- factor(d_orig$depth2,levels=c("0_10","10_20","20_30","30_40",
                                               "40_50","50_75","75_100","100_125","150_175"))
# Remove NA
d_orig <- na.omit(d_orig)

# Read in scaled, log-transformed data
d_log_scale <- read.csv("../Data/Clean/complete_spruce_logged_scaled.csv")
# Replace dash with underscore in depth ranges so the names don't confuse the modeling
d_log_scale$depth2 <- gsub("-","_",d_log_scale$depth2)
# Relevel depth variable
d_log_scale$depth2 <- factor(d_log_scale$depth2,levels=c("0_10","10_20","20_30","30_40",
                                               "40_50","50_75","75_100","100_125","150_175"))
# Make experimental temperature a factor so that its levels can be distinguished in the PCA
d_log_scale$Temp_experimental <- factor(d_log_scale$Temp_experimental)


# Read in DOT plot for SEM
dot_diagram <-readLines("./Data/dot_psem.txt")
# Read in PSEM summary object
spruce_num_psem_summary <- readRDS("./Data/psem_numeric_summary.RDS")

ui <- dashboardPage(
  title = "SPRUCE",
  
  # Title page
  dashboardHeader(
    title = "SPRUCE"
  ),
  
  # Sidebar
  dashboardSidebar(
    sidebarMenu(
      #menuItem("Introduction",tabName="intro"),
      menuItem("Exploration",tabName="exploration"),
      menuItem("Principal Component Analysis",tabName="PCA"),
      menuItem("SEM",tabName="sem")
    )
  ),
  
  # Body
  dashboardBody(
    tabItems(
      #### Introduction tab ####
      # tabItem(
      #   tabName = "intro",
      #   # Upload introduction (HTML file) here
      # ),
      #### Data exploration tab ####
      tabItem(
        tabName = "exploration",
        fluidRow(column(width=4,
                        box(width=NULL,
                            radioButtons(inputId = "y_choice",
                                               "Y-Axis Selection",
                                               choices = c("Bacteria Copy Number"="Bacteria_copy_dry","Archaeal Copy Number"="Archaea_copy_dry","Gravimetric Water Content"="GWC", "Dissolved Nitrogen"="DN_unfumigated_soil", "Microbial Biomass Nitrogen"="MBN", "Dissolved Organic Carbon"="DOC_unfumigated_soil", "Microbial Biomass Carbon"="MBC", "Soil Temperature"="temp")
                            )),
                        box(width=NULL,
                            radioButtons(inputId = "x_choice",
                                         "X-Axis Selection",
                                         choices = c("Experimental Temperature"="Temp_experimental","CO2 Treatment"="CO2_treatment","Depth"="depth2")
                                         
                            )),
                        box(width=NULL,
                            radioButtons(inputId = "colour_choice",
                                         "Colour By",
                                         choices = c("Experimental Temperature"="Temp_experimental","CO2 Treatment"="CO2_treatment","Depth"="depth2")
                                         
                            )),
                        box(width=NULL,
                            "Include ambient temperature treatment?",
                             switchInput(inputId = "include_ambient",
                                            label = "Include\nambient",
                                         onLabel = "Yes",
                                         offLabel = "No")
                            ),
                        box(width=NULL,
                            "Log-transform response variable?",
                            switchInput(inputId = "log_transform_anova",
                                           label = "Log-transform",
                                        onLabel = "Yes",
                                        offLabel = "No")
                        )
        ),
        column(width=8,
               box(width=NULL,
                   plotOutput("explore_plot")),
               box(width=NULL,
                   title = "Interactions",
                   tableOutput("interaction_table"),
                   # tableOutput("test_xletters"), # this and below: for testing
                   # tableOutput("test_forlegend"),
                   # textOutput("tukey_test"),
                   # tableOutput("tukeyp_test") 
               )
        ))
      ),
      #### PCA tab ####
      tabItem(
        tabName = "PCA",
        fluidRow(
          column(width=4,
                box(width=NULL,
                    checkboxGroupInput(inputId = "selected",
                                       "Explanatory Variables",
                                       choices = c("Bacteria Copy Number"="Bacteria_copy_dry","Archaeal Copy Number"="Archaea_copy_dry","Gravimetric Water Content"="GWC", "Dissolved Nitrogen"="DN_unfumigated_soil", "Microbial Biomass Nitrogen"="MBN", "Dissolved Organic Carbon"="DOC_unfumigated_soil", "Microbial Biomass Carbon"="MBC", "Soil Temperature"="temp"),
                                       selected = c("Bacteria_copy_dry","Archaea_copy_dry", "GWC", "DN_unfumigated_soil", "MBN", "DOC_unfumigated_soil", "MBC", "temp")
                    )),
                box(width=NULL,
                    radioButtons(inputId = "selection",
                                       "Colour By",
                                       choices = c("Experimental Temperature"="Temp_experimental","CO2 Treatment"="CO2_treatment","Depth"="depth2")
                                       
                    ))
                ),
          column(width=8,
                 box(width=NULL,
                     plotOutput("pca_plot"))
                 )
        )
      ),
      #### SEM tab ####
      tabItem(
        tabName = "sem",
        fluidRow(
          column(width=12,
                 box(width=NULL,
                     title = "Piecewise Structural Equation Model",
                     grVizOutput("psem_plot")),
                 box(width=6,
                     title = "Goodness of fit",
                     tableOutput("psem_chisq"),
                     tableOutput("psem_c"),
                     tableOutput("psem_aic")),
                 box(width=6,
                     title = "R-squared",
                     tableOutput("psem_r2"))
                 )
        )
      )
    )
  )
)

###### Server ###### 

server <- function(input, output){
  #### PCA plot ####
  output$pca_plot <- renderPlot({
    autoplot(prcomp(d_log_scale[,input$selected], scale = T), 
             data = d_log_scale, 
             colour= input$selection,
             frame = T,
             frame.type = "norm") + 
      theme_bw()
    })
  
  #### ANOVAs ####
  
  # Get data with or without ambient temperature
  d_amb <- reactive({
    if(input$include_ambient){
      d_orig
    } else{
      subset(d_orig,Temp_experimental != "Amb")
    }
  })
  

  
  # Select dependent variable for ANOVA
  lhs <- reactive(input$y_choice)
  
  # Log-transform if requested
  d <- reactive({
    dl <- d_amb()
    if(input$log_transform_anova){
      # If it includes 0, add 1
      if(any(dl[,lhs()]==0)){
        dl[,lhs()] <- dl[,lhs()] + 1
      }
      dl[,lhs()] <- log(dl[,lhs()])
    }
    return(dl)
  })
  
  # Select explantory variable(s)
  x_axis <- reactive(input$x_choice)
  fill_var <- reactive(input$colour_choice)
  
  # Create formula
  aov_formula <- reactive(paste0(lhs(),"~",paste(x_axis(),fill_var(),sep ="*")))
  
  # Run ANOVA
  aov_object <- reactive(aov(formula(aov_formula()), data=d()))
  
  # Run TukeyHSD
  tukey <- reactive(TukeyHSD(aov_object()))
  
  # Extract p values
  tukey_p <- reactive({
      lapply(tukey(),FUN=function(x){
      y <- x[,4]
      names(y) <- rownames(x)
      return(y)
    })
  })
  
  # Get letters for plotting
  int <- reactive({
    int_exists <- tryCatch(expr = (length(tukey()[[3]]) > -1), error = function(e){
      return(F)
    }
    )
    if(int_exists){
      # If the interaction was included in the model, check if there were any significant interactions
      # Which interactions were significant?
      which(tukey_p()[[3]]<0.05)
    } else{
      0
    }
  })
  
  # If there were any significant interactions, print them as a table
  interaction_table <- reactive({
    if(any(int() != 0)){
      # If you only have one interaction, it won't format it as a dataframe - deal with this case specially
      if(length(int())==1){
        temp_inttable <- t(as.data.frame(tukey()[[3]][int(),]))
        rownames(temp_inttable) <- rownames(tukey()[[3]])[int()]
      } else{
        temp_inttable <- tukey()[[3]][int(),]
      }
      colnames(temp_inttable) <- c("Difference","Lower","Upper","Adjusted p-value")
    } else{
      temp_inttable <- "No interaction effects"
    }
    return(temp_inttable)
  })
  
  # Return table
  output$interaction_table <- renderTable(interaction_table(),rownames = T)
  
  # Get letters for boxplot
  x_letters <- reactive({
    temp_xletters <- multcompLetters(tukey_p()[[x_axis()]])$Letters
    temp_xletters <- data.frame(level = names(temp_xletters),letter = temp_xletters)
    return(temp_xletters)
  })
  
  fill_letters <- reactive(multcompLetters(tukey_p()[[fill_var()]])$Letters)
  
  for_legend <- reactive({
    temp_legend <- paste(names(fill_letters())," (",fill_letters(),")",sep="")
    names(temp_legend) <- names(fill_letters())
    return(temp_legend)
  })
  
  # test output
 #  output$test_xletters <- renderTable(x_letters())
 #  output$test_forlegend <- renderTable(for_legend())
 #  output$tukey_test <- renderPrint(tukey())
 #  output$tukey_test <- renderPrint({    # Lettering DF
 #    x_let <- x_letters()
 #    for_x <- input$x_choice
 #    
 #    # Merge label with full dataframe
 #    forplot <- merge(d,x_let,by.x=for_x,by.y=colnames(x_let)[2])
 #    
 #    return(colnames(forplot))})
 #  output$tukeyp_test <- renderTable(tukey()[[3]])
 #  output$tukeyp_test <- renderTable({
 #    # Lettering DF
 #    x_let <- x_letters()
 #    for_x <- input$x_choice
 #    
 #    # Merge label with full dataframe
 #    forplot <- merge(d,x_let,by.x=for_x,by.y="level")
 #    
 #    return(head(forplot[,10:ncol(forplot)]))
 #  })
  
  #### Boxplots ####

  
  ## Make plot
  exploratory_plot <- reactive({
    # Get y-coordinate of letters
    letter_y <- max(d()[[lhs()]])+max(d()[[lhs()]])/10
    
    # Lettering DF
    x_let <- x_letters()
    
    # Get choices
    for_x <- input$x_choice
    for_y <- input$y_choice
    for_fill <- input$colour_choice
    for_letter <- "letter"
    
    # Merge label with full dataframe
    forplot <- merge(d(),x_let,by.x=for_x,by.y="level")
    
    # Make plot
    p <- ggplot(data=forplot, aes_string(x=for_x, y=for_y, fill = for_fill))+
      geom_boxplot()+
      geom_text(aes_string(x=for_x,y=letter_y,label=for_letter))+
      scale_fill_discrete(labels=for_legend())+
      theme_bw()

    return(p)
  })
  
  ## Show plot
  output$explore_plot <- renderPlot(exploratory_plot())
  
  #### SEM ####
  # Interactive plot
  output$psem_plot <- renderGrViz(grViz(dot_diagram))
  
  # Model fit
  output$psem_chisq <- renderTable(spruce_num_psem_summary$ChiSq)
  output$psem_c <- renderTable(spruce_num_psem_summary$Cstat)
  output$psem_aic <- renderTable(spruce_num_psem_summary$AIC)
  
  # R-squared
  output$psem_r2 <- renderTable({
    response_vars <- c("Dissolved organic carbon (DOC)","Dissolved nitrogen (DN)",
                       "Soil temperature","Gravimetric water content (GWC)",
                       "Bacteria Abundance","Archaea Abundance",
                       "Microbial biomass nitrogen (MBN)",
                       "Microbial biomass carbon (MBC)")
    out <- spruce_num_psem_summary$R2[,c("Response","R.squared")]
    out$Response <- response_vars
    return(out)
    })
}

shinyApp(ui = ui, server = server)
