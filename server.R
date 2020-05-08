shinyServer(function(input, output){

    result <- reactive({
        result <- get_acc_knn_hyper(as.numeric(input$n_genes), as.numeric(input$k_n), as.numeric(input$k_fold), as.numeric(input$range_k), as.numeric(input$n_sim))
    })


    result_f1 <- reactive({
        result <- get_f1_knn_hyper(as.numeric(input$n_genes), as.numeric(input$k_n), as.numeric(input$k_fold), as.numeric(input$range_k))
    })


    output$hyperPlotAccLine <- renderPlot({
        result = result()
        plot(result$k_vector, result$acc_vector, type="l", panel.first=grid(), xlab="K", ylab="Accuracy", main="Mean Accuracy")

        #df = data.frame(result$k_vector, result$acc_vector)
        #colnames(df) = c('k_vector','acc_vector')
        #ggplot(data = df, aes(x=k_vector, y=acc_vector)) +
        #    geom_line()
    })

    output$hyperPlotAccBox <- renderPlot({
        result = result()
        boxplot(result$acc_vector_2d, names=result$k_vector, main = "Accuracy", xlab="K", ylab="Accuracy Spread")
    })


    output$hyperPlotFNLine <- renderPlot({
        result = result()
        #plot(result$k_vector, result$fn_vector, type="l",panel.first=grid(), main="Mean False Negative Rate", xlab="K", ylab="Rate")
        # na.rm to remove weird ylim bug!
        plot(result$k_vector, colMeans(result$fn_vector_2d, na.rm=TRUE), type="l",panel.first=grid(), main="Mean False Negative Rate", xlab="K", ylab="Rate")
    })

    output$hyperPlotFNBox <- renderPlot({
        result = result()
        boxplot(result$fn_vector_2d, names=result$k_vector, main="False Negative Rate Spread", xlab="K", ylab="Spread")
        
    })


    output$hyperPlotF1Line <- renderPlot({
        #result = result_f1()
        result = result()
        plot(result$k_vector, result$f1_vector, type="l",panel.first=grid(), main="Mean F1 Score", xlab="K", ylab="F1 Score")
    })

    output$hyperPlotF1Box <- renderPlot({
        #result = result_f1()
        result = result()
        boxplot(result$f1_vector_2d, names=result$k_vector, main="F1 Score Spread", xlab="K", ylab="F1 Score Spread")
    })


    output$hyperPlotFPBox <- renderPlot({
        result = result()
        boxplot(result$fp_vector_2d, names=result$k_vector, main="False Positive Rate Spread", xlab="K", ylab="Spread")
        
    })


})
