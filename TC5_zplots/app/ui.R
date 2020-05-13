shinyUI(
  navbarPage(
    title = "TimeCourse 5 - Cflo",
    theme = shinytheme("flatly"),
    fluid = T,
    collapsible = TRUE,

    tabPanel("Z-plots",

             includeCSS("www/custom.css"),

             # Sidebar
             sidebarLayout(
               sidebarPanel(
                 width = 3,

                 h4(tags$b("Select Cflo gene")),
                 # sliderInput('l1', 'Distance:', ticks = FALSE,
                 #             min = 0, max = 100, value = 20, step = 1, width = "100%"),
                 # sliderInput('q1', 'Quality:', ticks = FALSE,
                 #             min = 0, max = 2, value = 1, step = 0.01, width = "100%"),
                 # sliderInput('k1', 'Light intensity:', ticks = FALSE,
                 #             min = 0, max = 100, value = 20, step = 1, width = "100%"),
                 # Drop-down menu to select a gene ----
                 selectizeInput(
                     'gene', label = "", choices = cflo.annots.exp %>% filter(gene_name %in% bg.genes) %>% pull(gene_name),
                     options = list(create = TRUE)),

               # Input: Selector for choosing dataset ----
               selectInput(inputId = "plots",
                           label = "Choose plot type:",
                           choices = c("Zscores", "log2(exp)", "Raw FPKMs"))),
               # 
               # # Input: Numeric entry for number of obs to view ----
               # numericInput(inputId = "obs",
               #              label = "Number of observations to view:",
               #              value = 10)),


               #   hr(),
               #
               #   h4(tags$b("Source #2")),
               #   sliderInput('l2', 'Distance:', ticks = FALSE,
               #               min = 0, max = 100, value = 25, step = 1, width = "100%"),
               #   sliderInput('q2', 'Quality:', ticks = FALSE,
               #               min = 0, max = 2, value = 1, step = 0.01, width = "100%"),
               #   sliderInput('k2', 'Light intensity:', ticks = FALSE,
               #               min = 0, max = 100, value = 20, step = 1, width = "100%")
               # ),

               # Main panel for displaying outputs ----
               mainPanel(
                 
                 # # Output: Verbatim text for data summary ----
                 # verbatimTextOutput("summary"),
                 
                 fluidRow(
                   plotlyOutput("myplot")
                 ),
                 
                 # # Output: HTML table with requested number of observations ----
                 # tableOutput("summary"),
                 
                 # Output: HTML table with requested number of observations ----
                 tableOutput("rhy.summary"),
                 
                 # # Output: HTML table with requested number of observations ----
                 tableOutput("blast.annot"),
                 
                 tableOutput("DEG")
                 
               )
              )
             ),
    
    tabPanel("Gene Info",
             
             includeCSS("www/custom.css"),
             
             # Sidebar
             sidebarLayout(
               sidebarPanel(
                 width = 1),
               
               # Main panel for displaying outputs ----
               mainPanel(
                 
                 # # Output: HTML table with requested number of observations ----
                 DT::dataTableOutput("summary", width = "120%")
               )
             )
    ),
    
    tabPanel("Manipulation?",
             
             includeCSS("www/custom.css"),
             
             # Sidebar
             sidebarLayout(
               sidebarPanel(
                 width = 1),
               
               # Main panel for displaying outputs ----
               mainPanel(
                 
                 # # Output: HTML table with requested number of observations ----
                 DT::dataTableOutput("de.summary", width = "120%")
                  )
               )
      ),
    
      tabPanel("Orthologs",
             
             includeCSS("www/custom.css"),
             
             # Sidebar
             sidebarLayout(
               sidebarPanel(
                 width = 1),
               
               # Main panel for displaying outputs ----
               mainPanel(
                 
                 # # Output: HTML table with requested number of observations ----
                 DT::dataTableOutput("blastp.summary", width = "120%")
                 
               )
            )
      ),
    
    tabPanel("Cflo_annots",
             
             includeCSS("www/custom.css"),
             
             # Sidebar
             sidebarLayout(
               sidebarPanel(
                 width = 1),
               
               # Main panel for displaying outputs ----
               mainPanel(
                 
                 # # Output: HTML table with requested number of observations ----
                 DT::dataTableOutput("cflo.summary", width = "120%")
                 
               )
             )
      )
  )
)



