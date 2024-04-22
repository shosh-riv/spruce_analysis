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
  
  
  #### Boxplots ####
  output$explore_plot <- renderPlot({
    ggplot(data=d, aes_string(x=input$x_choice, y=input$y_choice, fill = input$colour_choice))+
      geom_boxplot()+
      theme_bw()
    })
}

shinyApp(ui = ui, server = server)
