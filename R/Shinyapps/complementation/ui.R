library(shiny)

# Define UI
shinyUI(fluidPage(
  
  # Application title
  titlePanel("Mutation Rate Calculator"),
  
  # Sidebar
  sidebarLayout(
    sidebarPanel(
       numericInput("N",
                   "Average Cells Per Culture",
                   min = 0,
                   max = 1E128,
                   value = 1),
       fileInput(inputId = "mutants",
                    label = "Upload Mutants")
    ),
    
    # Draw plot
    mainPanel(
       #tableOutput("out")
       textOutput("debug")
    )
  )
))
