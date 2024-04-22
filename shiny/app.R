#####################################
### SPRUCE data visualization app ###
### Autumn Pereira, Linnea Smith  ###
### Cornell U, Spring 2024        ###
#####################################
library(shiny)
library(shinydashboard)
library(ggplot2)
library(ggfortify)
library(plotly)
library(multcompView)

# Read in data (note: once dataset finalized, move to shiny folder)
d <- read.csv("../Data/Clean/complete_combined_spruce_data.csv")

ui <- dashboardPage(
  title = "SPRUCE",
  
  # Title page
  dashboardHeader(
    title = "SPRUCE"
  ),
  
  # Sidebar
  dashboardSidebar(
    sidebarMenu(
      menuItem("Introduction",tabName="intro"),
      menuItem("Exploration",tabName="exploration"),
      menuItem("Principal Component Analysis",tabName="PCA"),
      menuItem("SEM",tabName="sem")
    )
  ),
  
  # Body
  dashboardBody(
    tabItems(
      #### Introduction tab ####
      tabItem(
        tabName = "intro",
        # Upload introduction (HTML file) here
      ),
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
                                         
                            ))
        ),
        column(width=8,
               box(width=NULL,
                   tableOutput("interaction_table"),
                   #tableOutput("test_xletters"),
                   #textOutput("test_forlegend") # for testing
               ),
               box(width=NULL,
                   plotOutput("explore_plot"))
        ))
      ),
      #### Regression tab ####
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
        fluidRow()
      )
    )
  )
)

###### Server ###### 

server <- function(input, output){
  #### PCA plot ####
  output$pca_plot <- renderPlot({
    autoplot(prcomp(d[,input$selected], scale = T), 
             data = d, 
             colour= input$selection) + 
      theme_bw()
    })
  
  #### ANOVAs ####
  
  # Select dependent variable for ANOVA
  lhs <- reactive(input$y_choice)
  
  # Select explantory variable(s)
  x_axis <- reactive(input$x_choice)
  fill_var <- reactive(input$colour_choice)
  
  # Create formula
  aov_formula <- reactive(paste0(lhs(),"~",paste(x_axis(),fill_var(),sep ="*")))
  
  # Run ANOVA
  aov_object <- reactive(aov(formula(aov_formula()), data=d))
  
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
    if(!is.null(tukey()[[3]])){
      # If the interaction was included in the model, check if there were any significant interactions
      # Which interactions were significant?
      which(tukey_p()[[2]]<0.05) # change 2 back to 3!
    }
  })
  
  # If there were any significant interactions, print them as a table
  interaction_table <- reactive({
    if(length(int()) > 0){
      temp <- tukey()[[2]][int,] # change 2 back to 3!
      colnames(temp) <- c("Comparison","Difference","Lower","Upper","Adjusted p-value")
    }
    return(temp)
  })
  
  # Return table
  output$interaction_table <- renderTable(interaction_table())
  
  # Get letters for boxplot
  x_letters <- reactive({
    temp <- multcompLetters(tukey_p()[[x_axis()]])$Letters
    temp <- data.frame(level = names(temp),letter = temp)
    return(temp)
  })
  
  fill_letters <- reactive(multcompLetters(tukey_p()[[fill_var()]])$Letters)
  
  for_legend <- reactive({
    temp <- paste(names(fill_letters())," (",fill_letters(),")",sep="")
    names(temp) <- names(fill_letters())
    return(temp)
  })
  
  # test output
  # output$test_xletters <- renderTable(x_letters())
  # output$test_forlegend <- renderText(for_legend())
  
  #### Boxplots ####
  output$explore_plot <- renderPlot({
    ggplot(data=d, aes_string(x=input$x_choice, y=input$y_choice, fill = input$colour_choice))+
      geom_boxplot()+
      geom_text(data=x_letters(),aes(x=level,y=(max(d[[lhs()]])+max(d[[lhs()]])/10)),label=x_letters()$letter)+
      scale_fill_discrete(labels=for_legend())
      theme_bw()
    })
}

shinyApp(ui = ui, server = server)
