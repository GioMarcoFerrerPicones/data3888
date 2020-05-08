library(shiny)
shinyUI(
    pageWithSidebar(
        headerPanel("KNN Simulator"),


        sidebarPanel(


            checkboxGroupInput("hyper_checkbox", "Hyperparameter Accuracy Measures", choices = list("Accuracy" = 1, "F1 Score" = 2, "False Negative Rate" = 3), selected = c(1,2,3)),

            textInput("n_genes", "Please select number of genes", 4),
            textInput("k_fold", "Please select number of folds", 5),
            textInput("n_sim", "Please select number of simulations", 10),

            sliderInput("range_k", "Range of K Neighbours", min = 1, max = 30, value = c(1, 10))


        ),


        mainPanel(

            plotOutput("hyperPlotAccLine"),

            plotOutput("hyperPlotAccBox"),

            conditionalPanel(
                condition = "hyper_checkbox.items.includes('2')",
                plotOutput("hyperPlotF1Line")
            ),

            conditionalPanel(
                condition = "hyper_checkbox.items.includes('2')",
                plotOutput("hyperPlotF1Box")
            ),

            conditionalPanel(
                condition = "hyper_checkbox.items.includes('3')",
                plotOutput("hyperPlotFNLine")
            ),

            conditionalPanel(
                condition = "hyper_checkbox.items.includes('3')",
                plotOutput("hyperPlotFNBox")
            )

            #plotOutput("hyperPlotFPBox")

            #plotOutput("myPlot")
        )

    )
)
