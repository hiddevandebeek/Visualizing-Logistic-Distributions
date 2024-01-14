
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
library(plotly)
library(MASS)
library(shinythemes)
library(matlib)

ui <- shinyUI(navbarPage(
  theme = shinytheme("flatly"),  # Applying a pre-designed theme
  "Distribution Visualizer", 
  # Home tab
  tabPanel("Graph", id = "navibar", value = "home",
           fluidRow(
             div(style = "padding-left: 20px; padding-top: 20px;", # Adding padding to the left and top of the title
                 h2("Visualizing Distribution Changes")
                 )
             ),
           # Row for sliders and plots
           fluidRow(
             # Column for sliders
             column(4,
                    sliderInput("auc",
                                "AUC Size:",
                                min = 0.5,
                                max = 0.99,
                                value = 0.75,
                                step = 0.01),
                    sliderInput("population",
                                "Population Size:",
                                min = 100,
                                max = 50000,
                                value = 10000,
                                step = 100),
                    sliderInput("delta_sigma",
                                "Delta Sigma:",
                                min = 0,
                                max = 0.99,
                                value = 0,
                                step = 0.01),
                    sliderInput("z", 
                                "Correlation between Predictors:",
                                min = 0,
                                max = 0.99,
                                value = 0,
                                step = 0.01)
             ),
             
             # Column for plots
             column(8,
                    tabsetPanel(
                      tabPanel("2D Plot", plotOutput("distPlot2d")),
                      tabPanel("3D Plot", plotlyOutput("distPlot3d")),
                      tabPanel("Model Matrices", 
                               verbatimTextOutput("matrixOutput"))
                    ),
             )
           )
  ),
  tabPanel("Explanation", id = "explanation",
           fluidPage(
             titlePanel("Welcome to the Distribution Visualizer"),
             fluidRow(
               div(
                 p("Welcome to the Distribution Visualizer! This web application provides you with a powerful tool to gain deeper insights into data distribution changes within a two-class dataset. By adjusting various parameters, you can visualize and understand how they influence the data distributions."),
                 h4("How to Use:"),
                 p("To make the most of this tool, follow these steps:"),
                 HTML("<ol>
              <li>Use the 'AUC Size' slider to control the area under the curve (AUC) of the distributions, influencing the shape and spread of data.</li>
              <li>Customize the 'Population Size' slider to specify the number of samples in each class, reflecting your dataset's characteristics.</li>
              <li>Experiment with the 'Delta Sigma' slider to change the difference in variances between classes, providing insights into variance-related effects.</li>
              <li>Fine-tune the 'Correlation between Predictors' slider to set the correlation between predictors, which plays a crucial role in data distribution.</li>
              <li>Explore the '2D Plot' and '3D Plot' tabs to gain visual perspectives on data distribution changes in different dimensions.</li>
              <li>Analyze the 'Model Matrices' tab for detailed information about covariance matrices and class means, aiding in statistical understanding.</li>
              <li>Feel free to interact with the controls, iterate, and interpret the visualizations to enhance your dataset understanding and exploratory data analysis.</li>
            </ol>"),
                 p("With its user-friendly interface and powerful visualization capabilities, the Distribution Visualizer empowers you to uncover hidden patterns, outliers, and trends in your data. It serves as an invaluable tool for data scientists, analysts, and researchers seeking to deepen their understanding of complex datasets.")
               )
             )
           )
  )
  ,
  selected = "home")
)


# Define server logic required to draw a histogram
server <- shinyServer(function(input, output) {
  data <- reactive({
  auc    <- input$auc                        # AUC
  y      <- input$delta_sigma                # delta sigma 
  z      <- input$z                              # z 
  p      <- 1                                # proportion of correlated predictors
  npred  <- 2                                # number of predictors
  n      <- input$population                 # Number of samples in each class
  
  get_delta_mu <- function(){
    inn = scenario()      # get covariance matrices
    s0  = inn[[1]]              
    s1  = inn[[2]]                       
    A   = inv(s0 + s1)         # calculate matrix A
    Y   = sum(colSums(A))                  
    dmu = qnorm(auc) / sqrt(Y)
    return(dmu)
  }
  scenario <- function(){
    # set up correlations
    corr0  <-  matrix(0, npred, npred)         # matrix: set up for cov matrix, 0 on diagonals
    corr0[1:npred*p, 1:npred*p] = z            # class 0 
    diag(corr0) = 0
    corr1  <-  matrix(0, npred, npred)         # matrix: set up for cov matrix, 0 on diagonals
    corr1[1:npred*p, 1:npred*p] = (1-y)*z      # class 1
    diag(corr1) = 0
    # covariance structures
    sigma0 <-  diag(npred)  + corr0            # matrix: cov matrix of class 0 
    dsig   <-  diag(c(rep(y, npred)))          # matrix: difference in variances between classes 
    sigma1 <-  diag(npred) - dsig + corr1      # matrix: cov matrix of class 1
    return(list("sigma0"      = sigma0, 
                "sigma1"      = sigma1))
  }
  
  mu0    <-  c(rep(0,npred))            # vector: class 0 means
  mu1    <-  mu0 + get_delta_mu()       # vector: class 1 means 
  
  inn = scenario()      # get covariance matrices
  sigma1 <- inn[[2]]         # class 1
  sigma0 <- inn[[1]]         # class 0
  
  diseased_data <- mvrnorm(n, mu = mu1, Sigma = sigma1)
  non_diseased_data <- mvrnorm(n, mu = mu0, Sigma = sigma0)
  
  combined_data <- rbind(diseased_data, non_diseased_data)
  data <- data.frame(combined_data, class = c(rep(1, n), rep(0, n)))
  colnames(data) <- c(paste0("predictor", 1:(ncol(combined_data))), "class")
  
  data
  })
  
  output$distPlot2d <- renderPlot({
    req(data())
    data <- data()
    # create plot
    ggplot(data, aes(x = predictor1, fill = as.factor(class))) +
      geom_histogram(alpha = 0.5, position = 'identity', bins = round(input$population*0.02)) +
      ggtitle("Distribution of Diseased and Non-Diseased Classes for Predictor 1")
  })
  
  output$distPlot3d <- renderPlotly({
    req(data())
    data <- data()
    density1 <- with(data[data$class == 1, ], kde2d(predictor1, predictor2, n = 30))
    density0 <- with(data[data$class == 0, ], kde2d(predictor1, predictor2, n = 30))
    
    p1 <- plot_ly(z = ~density1$z) %>% 
      add_surface(
        x = ~density1$x,
        y = ~density1$y,
        colorscale = list(c(0, "red"), list(1, "pink")),
        opacity = 0.8
      )
    p2 <- plot_ly(z = ~density0$z) %>% 
      add_surface(
        x = ~density0$x,
        y = ~density0$y,
        colorscale = list(c(0, "blue"), list(1, "lightblue")),
        opacity = 0.8
      )
    plotly::subplot(p1, p2, nrows = 1, shareX = TRUE, shareY = TRUE)
  })
  
  output$matrixOutput <- renderPrint({
    # Reactive block to get input values and calculate matrices
    auc <- input$auc
    y <- input$delta_sigma
    z <- input$z
    npred <- 2
    p <- 1
    
    # Function to compute delta mu based on AUC
    get_delta_mu <- function(auc) {
      qnorm(auc) / sqrt(sum(colSums(inv(scenario()[[1]] + scenario()[[2]]))))
    }
    
    # Compute scenario covariance matrices
    scenario <- function() {
      corr0 <- matrix(0, npred, npred)
      corr0[1:npred*p, 1:npred*p] <- z
      diag(corr0) <- 0
      corr1 <- matrix(0, npred, npred)
      corr1[1:npred*p, 1:npred*p] <- (1-y)*z
      diag(corr1) <- 0
      sigma0 <- diag(npred) + corr0
      dsig <- diag(c(rep(y, npred)))
      sigma1 <- diag(npred) - dsig + corr1
      list(sigma0 = sigma0, sigma1 = sigma1)
    }
    
    # Call scenario function to get covariance matrices
    matrices <- scenario()
    sigma0 <- matrices$sigma0
    sigma1 <- matrices$sigma1
    mu0 <- rep(0, npred)
    mu1 <- mu0 + get_delta_mu(auc)
    
    # Combine mu and sigma into a list for output
    list("mu0" = mu0, "mu1" = mu1, "sigma0" = sigma0, "sigma1" = sigma1)
  })
  
})
# Run the application
shinyApp(ui = ui, server = server)