library(shiny)

# Define UI
shinyUI(fluidPage(
  
  # Application title
  titlePanel("Complementation Table Analyzer"),
  
  # Sidebar
  sidebarLayout(
    sidebarPanel(
      checkboxInput("checkMap", label = "Use Mapping", value = TRUE),
      fileInput(inputId = "complementation", label = "Upload Complementation"),
      fileInput(inputId = "mapping", label = "Upload Mapping")           
    ),
    
    # Draw plot
    mainPanel(
      plotOutput("out")
    )
  )
))
