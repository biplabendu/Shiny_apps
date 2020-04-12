

# Example -----------------------------------------------------------------

function(input, output, session) {
  
  # # Return the requested dataset ----
  # datasetInput <- reactive({
  #   switch(input$dataset,
  #          "foragers" = cflo.zscores.for,
  #          "nurses" = cflo.zscores.nur)
  # })
  # 
  # 
  # output$gene_desc <- renderText({
  #   paste('You selected', if (input$github == '') 'no' else input$gene,
  #         'gene.')
  # })  
  
  # Generate a summary of the dataset ----
  output$summary <- renderTable({
    df <- annot(input$gene)
    df
  })
  
  output$blast.annot <- renderTable({
    df <- annot(input$gene)
    df <- df %>%
      dplyr::select("Blast annotation" = Blast)
    df
  })
  
  output$DEG <- renderTable({
    df <- cuffdiff.results(input$gene) %>% 
      filter(treatment_1 == "Biting" & treatment_2 == "Control") %>% 
      select(treatment_1, treatment_2,
             logFC,
             significant,
             fpkm_1 = value_1, 
             fpkm_2 = value_2)
  })
  
  output$rhy.summary <- renderTable({
    df <- rhy.results(input$gene)
    df
  })
  
  output$de.summary <- renderTable({
    df <- cuffdiff.results(input$gene)
    df
  })
  
  
  output$gene <- renderPrint({
    g <- input$gene
    g
  })
  
  
  
  # # Show the first "n" observations ----
  # output$view <- renderTable({
  #   head(datasetInput(), n = input$obs)
  # })
  
  output$myplot <- renderPlotly({
    if(input$plots == "Zscores") {
      g <- z.plot(input$gene)
    }
    
    else if (input$plots == "log10(exp)"){
      g <- log.plot(input$gene, log=T)
    }
    
    else {
      g <- log.plot(input$gene, log=F)
    }
    
    ggplotly(g) %>%
      layout(legend = list(x = 0.5, y = 1.1, orientation = "h", xanchor = "center"),
             hovermode = "x")
  })
}



# Simon’s code ------------------------------------------------------------

# function(input, output) {
# 
#   res <- reactive({
#     parms <- list(iniN = input$N,
#                   iniS = rep(0,2),
#                   iniQ = rep(0,2),
#                   l = c(input$l1, input$l2),
#                   q = c(input$q1, input$q2),
#                   qe = 0,
#                   alpha = 0.01,
#                   beta = 5,
#                   gamma = 1200,
#                   eta = 2,
#                   psi = 0.05,
#                   rho = 0.01,
#                   k = c(100 - input$k1 + 1, 100 - input$k2 + 1),
#                   n = 2)
# 
#     out <- dede(c(N = parms$iniN, S = parms$iniS, Q = parms$iniQ),
#                 seq(0, 3600, length.out = 360),
#                 ode_sys,
#                 parms)
# 
#     data.frame(time = rep(out[, 1], 2),
#                S = c(rep("Source 1", nrow(out)), rep("Source 2", nrow(out))),
#                val = c(out[, 3], out[, 4]))
#   })
# 
#   output$the_display <- renderPlotly({
#     g <- ggplot(data = res(), aes(x = time, y = val, color = S)) +
#       geom_path(size = 0.75) +
#       theme_minimal(base_size = 16) +
#       theme(legend.title = element_blank()) +
#       xlab("Time") + ylab("Number of ants") +
#       scale_color_manual(values = c("#2678B2", "#FD7F28"))
# 
#     ggplotly(g) %>%
#       layout(legend = list(x = 0.5, y = 1.1, orientation = "h", xanchor = "center"),
#              hovermode = "x")
#   })
# 
# }
# 


# # My COVID_19 code ------------------------------------------------------------
# 
# # Code borrowed from: https://towardsdatascience.com/create-a-coronavirus-app-using-r-shiny-and-plotly-6a6abf66091d
# # Data from:
# # World: Johns Hopkins University Center for System Science and Engineering (JHU CCSE)
# # India: https://www.kaggle.com/sudalairajkumar/covid19-in-india#covid_19_india.csv
# 
# library(dplyr)
# library(tidyr)
# f1 = list(family="Courier New, monospace", size=12, color="rgb(30,30,30)")
# 
# # Add the server logic ----------------------------------------------------
# function(input, output, session) {
#   
#   data = reactive({
#     d = allData %>%
#       filter(Country == input$country)
#     if(input$state != "<all>") {
#       d = d %>% 
#         filter(State == input$state) 
#     } else {
#       d = d %>% 
#         group_by(date) %>% 
#         summarise_if(is.numeric, sum, na.rm=TRUE)
#     }
#     
#     d %>%
#       mutate(
#         dateStr = format(date, format="%b %d, %Y"),    # Jan 20, 2020
#         NewConfirmed=CumConfirmed - lag(CumConfirmed, default=0),
#         NewRecovered=CumRecovered - lag(CumRecovered, default=0),
#         NewDeaths=CumDeaths - lag(CumDeaths, default=0)
#       )
#   })
#   
#   observeEvent(input$country, {
#     states = allData %>%
#       filter(Country == input$country) %>% 
#       pull(State)
#     states = c("<all>", sort(unique(states)))
#     updateSelectInput(session, "state", choices=states, selected=states[1])
#   })
#   
#   countries = sort(unique(allData$Country))
#   
#   updateSelectInput(session, "country", choices=countries, selected="India")
#   
#   # renderBarPlot = function(varPrefix, legendPrefix, yaxisTitle) {
#   #   renderPlotly({
#   #     data = data()
#   #     plt = data %>% 
#   #       plot_ly() %>%
#   #       config(displayModeBar=FALSE) %>%
#   #       layout(
#   #         barmode='group', 
#   #         xaxis=list(
#   #           title="", tickangle=-90, type='category', ticktext=as.list(data$dateStr), 
#   #           tickvals=as.list(data$date), gridwidth=1), 
#   #         yaxis=list(
#   #           title=yaxisTitle
#   #         ),
#   #         legend=list(x=0.05, y=0.95, font=list(size=15), bgcolor='rgba(240,240,240,0.5)'),
#   #         font=f1
#   #       )
#   #     for(metric in input$metrics) 
#   #       plt = plt %>%
#   #       add_trace(
#   #         x= ~date, y=data[[paste0(varPrefix, metric)]], type='bar', 
#   #         name=paste(legendPrefix, metric, "Cases"),
#   #         marker=list(
#   #           color=switch(metric, Deaths='rgb(200,30,30)', Recovered='rgb(30,200,30)', Confirmed='rgb(100,140,240)'),
#   #           line=list(color='rgb(8,48,107)', width=1.0)
#   #         )
#   #       )
#   #     plt
#   #   })
#   # }
#   
#   # output$cumulatedMetrics = renderBarPlot("Cum", legendPrefix="Cumulated", yaxisTitle="Cumulated Cases")
#   # output$dailyMetrics = renderBarPlot("New", legendPrefix="New", yaxisTitle="New Cases per Day")
#   output$cumulatedMetrics <- renderPlotly({
#     g <- ggplot(data = data, aes(x = CumConfirmed, y = CumDeaths)) +
#       # geom_path(size = 0.75) +
#       theme_minimal(base_size = 16) +
#       geom_point() +
#       theme(legend.title = element_blank()) +
#       xlab("X") + ylab("Y") 
#       # scale_color_manual(values = c("#2678B2", "#FD7F28"))
#     
#     ggplotly(g) %>%
#       layout(legend = list(x = 0.5, y = 1.1, orientation = "h", xanchor = "center"),
#              hovermode = "x")
#   })
# }
