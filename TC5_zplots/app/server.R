

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
  output$summary <- DT::renderDataTable({
    df <- annot(input$gene)
    df
  }, options = list(scrollX = TRUE))
  
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
             DEG_biting,
             fpkm_1 = value_1, 
             fpkm_2 = value_2) 
  })
  
  output$rhy.summary <- renderTable({
    df <- rhy.results(input$gene)
    df
  })
  
  output$de.summary <- DT::renderDataTable({
    df <- cuffdiff.results(input$gene) %>% 
      select(-DEG_biting)
    df
  }, options = list(scrollX = TRUE))
  
  
  output$gene <- renderPrint({
    g <- input$gene
    g
  })
  
  output$blastp.summary <- DT::renderDataTable({
    g <- as.character(input$gene)
    df <- tbl(db, "blastp_rhy_summary") %>% filter(cflo_gene == g) %>% collect() %>% as.data.frame()
    df
  }, options = list(scrollX = TRUE, 
                    columnDefs = list(list(width = '400px', targets = c(3,4)))
                    ))
  
  output$cflo.summary <- DT::renderDataTable({
    
    rhy.for <- for.ejtk %>% filter(rhy=="yes") %>% pull(gene_name)
    rhy.nur <- nur.ejtk %>% filter(rhy=="yes") %>% pull(gene_name)
    # add degs for-nur
    # add degs control-biting
    
    
    df <- 
      cflo.annots.exp %>% 
      select(gene_name, blast_annotation = old_annotation,
             GOs, pfams, signalP, TMHMM) %>% 
      filter(gene_name %in% union(rhy.for, rhy.nur)) %>% 
      mutate(rhy_for = ifelse(gene_name %in% rhy.for, "yes", "no")) %>% 
      mutate(rhy_nur = ifelse(gene_name %in% rhy.nur, "yes", "no")) %>% 
      select(gene_name, blast_annotation, rhy_for, rhy_nur, everything()) %>% 
      left_join(tbl(blast_db,"dmel") %>%
                    select(cflo_gene, dmel_name) %>%
                    collect() %>%
                    rename(gene_name = cflo_gene) %>% 
                    distinct(gene_name, dmel_name, .keep_all = T),
                by=c("gene_name"))
    
    df
  }, options = list(scrollX = TRUE, 
                    columnDefs = list(list(width = '400px', targets = c(3,4)))
  ))
  
  
  
  # # Show the first "n" observations ----
  # output$view <- renderTable({
  #   head(datasetInput(), n = input$obs)
  # })
  
  output$myplot <- renderPlotly({
    if(input$plots == "Zscores") {
      g <- z.plot(input$gene)
    }
    
    else if (input$plots == "log2(exp)"){
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



