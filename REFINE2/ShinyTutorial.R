library(shiny)
ui <- fluidPage(
  sliderInput(inputId = "num",
              label = "Choose a number",
              value = 25, min = 1, max = 100),
  textInput(inputId = "title",
            label = "Write a title", 
            value = "Histogram of Random Normal Values"), 
  plotOutput("hist") #,
  # verbatimTextOutput("stats") 
)

server <- function(input, output) {
  output$hist <- renderPlot({
      hist(rnorm(input$num), main = isolate({input$title}))
    })
  # data <- reactive({     
  #     rnorm(input$num)
  #   })
  # output$hist <- renderPlot({    
  #     hist(data(), main = input$title)
  #   })
  # output$stats <- renderPrint({     
  #     summary(data())   
  #   }) 
}
shinyApp(ui = ui, server = server)


library(shiny)
runApp()
library(rsconnect)
deployApp()
