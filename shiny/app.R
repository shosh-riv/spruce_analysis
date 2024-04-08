#####################################
### SPRUCE data visualization app ###
### Autumn Pereira, Linnea Smith  ###
### Cornell U, Spring 2024        ###
#####################################
library(shiny)
library(shinydashboard)

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
      menuItem("Regression",tabName="regression"),
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
        fluidRow()
      ),
      #### Regression tab ####
      tabItem(
        tabName = "regression",
        fluidRow(
          column(width=4,
                box(width=NULL,
                    title="PCA Select",
                    checkboxGroupInput(inputId = "pca_incl_vars",
                                       "Variables to include",
                                       choices = c("Var1"="Var1","Var2"="Var2","Var3"="Var3"),
                                       selected = "Var1"
                          
                    ))
                ),
          column(width=8,
                 box(width=NULL,
                     title="PCA",
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
  output$pca_plot <- renderPlot(plot(1:10,1:10,ylab=input$pca_incl_vars[1]))
}

shinyApp(ui = ui, server = server)
