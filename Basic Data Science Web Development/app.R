#BASIC DATA SCIENCE WEB DEVELOPMENT
#-------------------------------------------------------------------------------

#libraries
library(shiny)
library(plotly)
library(tibble)
library(tidyverse)
library(tidyr)
library(readr)
library(dplyr)
library(ggplot2)
library(gapminder)



#read the data
gm <- read_csv("gapminder_clean.csv")

#separate dataframe to remove categorical data from 
gm2 <- gm %>%
  select(-`Country Name`, -continent, -Year,-...1)


# Define UI ----
ui <- fluidPage(
  column(3,offset = 4, titlePanel("Explore Gapminder Data with Shiny")),
  headerPanel('Graphs'),
  mainPanel(
    plotlyOutput('plot')
  ),
  sidebarPanel(
    #variable selection for x-axis
    selectInput(inputId ='xvrbl', #The input slot that will be used to access the value.
                label = 'X-Axis Variable', #Display label for the control, or NULL for no label.
                choices = colnames(gm2), #List of values to select from
                selected = 'CO2 emissions (metric tons per capita)'
    ),
    
    checkboxInput(inputId = "LogX", 
                  label = "Log Transform", 
                  value = FALSE),
    
    #variable selection for y-axis
    selectInput(inputId ='yvrbl', #The input slot that will be used to access the value.
                label = 'Y-Axis Variable', #Display label for the control, or NULL for no label.
                choices = colnames(gm2), #List of values to select from
                selected = 'gdpPercap'
    ),
    
    checkboxInput(inputId = "LogY", 
                  label = "Log Transform", 
                  value = FALSE),
    
    # date range - slider
    sliderInput(
      inputId = "time",
      label = "Years",
      min = min(gm$Year),
      max = max(gm$Year),
      step = 5,
      value = range(gm$Year)
    )
  )
)

server <- function(input, output) {
  x <- reactive({
    x <- dat()[[input$xvrbl]]
    if (input$LogX) x <- log(x)
    return(x)
  })
  
  y <- reactive({
    y <- dat()[[input$yvrbl]]
    if (input$LogY) y <- log(y)
    return(y)
  })
  
  dat <- reactive({
    subset(gm, Year >= input$time[[1]], Year <= input$time[[2]])
  })
  
  output$plot <- renderPlotly({
    plot_ly(
      x = x(),
      y = y(),
      type = "scatter",
      mode = "markers",
      color = dat()$continent
    ) %>% 
      layout( xaxis = list( title=input$xvrbl), 
              yaxis = list( title=input$yvrbl ) )
  })
}

# Run the app
shinyApp(ui = ui, server = server)
