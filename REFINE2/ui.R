# 
# fluidPage(
#   titlePanel("Realistic Evaluations of Finite sample INference using Efficient Estimators (REFINE2)"),
#   # SMILE4: Simulating Multiple ImpLementations to Evaluate Efficient Effect Estimators
#   # "REVUE: Realistic Evaluations of Variability Using Efficient Estimators" 
#   # titlePanel("T-PEE: Tool for Pratical Evaluation of Effect Estimation"),
#   numericInput("obs", "Number of Simulation Samples:", 3, min = 2, max = 500),
#   textInput("path", "Path of the data","./plas_data.csv"),
#   sidebarLayout(
#   sidebarPanel(h4("Models for simulation"),
#                textInput(inputId = "expForm",
#                          label = "PS Model", 
#                          value = "A ~ 1+VAR_1 + VAR_2"), 
#                textInput(inputId = "outForm",
#                          label = "Outcome Model", 
#                          value = "Y ~ 1 +A + VAR_1 + VAR_2"),
#                h4("Models for estimation"),
#     textInput(inputId = "expForm.est",
#               label = "PS Model", 
#               value = "A ~ 1+VAR_1 + VAR_2"), 
#     textInput(inputId = "outForm.est",
#               label = "Outcome Model", 
#               value = "Y ~ 1 +A + VAR_1 + VAR_2"), 
#     h4("Estimation Method"),
#     selectInput('mtd', 'Estimation Method', c("IPW","GComp","TMLE","AIPW","DCTMLE","DCAIPW"),
#                 selected = "IPW"),
#     selectInput('is.par', 'Library (not for IPW and GComp)', c("Smooth","Non-Smooth"),
#                 selected = "Smooth"),
#     actionButton("run", "Run")
#   ),
#   
#   mainPanel(
#     textOutput("exp.Form"),textOutput("out.Form"),textOutput("exp.Form.est"),textOutput("out.Form.est"),
#     plotOutput("plot"),
#     htmlOutput("res.text.1"),htmlOutput("res.text.2"),htmlOutput("res.text.3")
#   )
#   ),
#   fluidRow(
#     column(12,
#            tableOutput('table')
#     )
#   )
# )

fluidPage(
  titlePanel("Realistic Evaluations of Finite sample INference using Efficient Estimators (REFINE2)"),
  # SMILE4: Simulating Multiple ImpLementations to Evaluate Efficient Effect Estimators
  # "REVUE: Realistic Evaluations of Variability Using Efficient Estimators" 
  # titlePanel("T-PEE: Tool for Pratical Evaluation of Effect Estimation"),
  numericInput("obs", "Number of Simulation Samples:", 3, min = 2, max = 500),
  textInput("path", "Path of the data","./plas_data.csv"),
  sidebarLayout(
    sidebarPanel(h4("Models for simulation"),
                 textInput(inputId = "expForm",
                           label = "PS Model", 
                           value = "A ~ 1+VAR_1 + VAR_2"), 
                 textInput(inputId = "outForm",
                           label = "Outcome Model", 
                           value = "Y ~ 1 +A + VAR_1 + VAR_2"),
                 h4("Models for estimation"),
                 textInput(inputId = "expForm.est",
                           label = "PS Model", 
                           value = "A ~ 1+VAR_1 + VAR_2"), 
                 textInput(inputId = "outForm.est",
                           label = "Outcome Model", 
                           value = "Y ~ 1 +A + VAR_1 + VAR_2"), 
                 h4("Estimation Method"),
                 selectInput('mtd', 'Estimation Method', c("IPW","GComp","TMLE","AIPW","DCTMLE","DCAIPW"),
                             selected = "IPW"),
                 selectInput('is.par', 'Library (not for IPW and GComp)', c("Smooth","Non-Smooth"),
                             selected = "Smooth"),
                 actionButton("run", "Run")
    ),
    
    mainPanel(
      textOutput("exp.Form"),textOutput("out.Form"),textOutput("exp.Form.est"),textOutput("out.Form.est"),
      htmlOutput("res.text.1"),htmlOutput("res.text.2"),htmlOutput("res.text.3"),
      plotOutput("plot")
    )
  ),
  fluidRow(
    column(12,
           tableOutput('table')
    )
  )
)


# shiny::runApp("~/Desktop/HuangGroup/cvtmle_plasmode/Code/REFINE2")
# rsconnect::deployApp()
# rsconnect::configureApp("shinyapp", size="xxlarge")