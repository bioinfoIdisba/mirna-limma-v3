#library(plotly)
#library(rintrojs)
# BiocManager::install("enrichplot")
library(shinymanager)
library(enrichplot)



library(DiagrammeR)
library(shinydashboard)
library(shinydashboardPlus)
#library(shinyhelper)
library(shiny)
library(shinyFiles)
library(shinybusy)



# dataframe that holds usernames, passwords and other user data

source("/usr/local/lib/R/site-library/AnalisisArrays2/templates/llibreries_Informes.R")

add_class <- function(x, class) {
  x$attribs <- append(x$attribs, list(class = class))
  x
}


header<- dashboardHeader(title = span(icon("fire-burner"),"miRNA with Limma"),
                         titleWidth = 400, 
                         dropdownMenuOutput("messageMenu"))


sidebar<- dashboardSidebar(  
  hr(),
  (" UPLOAD DATA AND SELECT ANALYSIS"),
  hr(),
  
  disable = FALSE,
  width = 400,
  collapsed = F,
  
  sidebarMenu(

#menu####    
menuItem("Introduce data to analyze",
     menuItem("Upload celfiles file",
              
              fileInput("file_input_celfiles","Select file for celfiles",accept = ".CEL",
                                               multiple = TRUE))%>% add_class("upload"),
 
     ###INPUT FILE FENOTIP_IN###  
     menuItem("Select file of fenotips",fileInput("file_input_fenotip",
                                                      "Select file for fenotip",
                                                      accept=".csv"
                                                      ),
              
              fluidRow(
                column(3,radioButtons("sep", "Separator",choices = c(Comma = ",",Semicolon = ";",Tab = "\t"),
                                      selected = ";")),
                column(3,radioButtons("dec", "Decimal",choices = c(Comma = ",",Dot = "."),selected = ","))
                # column(3,radioButtons("quote", "Quote",choices = c(Comma = "\"'",Dot = ""),selected = ","))
              )
              
              
              ),
     div(style="display:inline;vertical-align:left;",
         actionButton("execute_pre","Prepare data and normalize",
                      style = "color: #fff; background-color: #367fa9; border-color: #fff;padding: 5px 14px 5px 14px;margin: 50px 5px 5px 5px; ")
     )
     ),
menuItem("Select analysis",status = "primary", solidHeader = TRUE,collapsible = TRUE,
         checkboxGroupInput("analisis", "",inline =F,
                            choices = c(
                              "Quality Control" = "qc",
                              "Fold Change"="fc",
                              "Find Targets of miRNA"="target",
                              "Biological significance"="enrichment"
                            ),
                            selected=""),
         # width = 250,
         collapsed=TRUE),
     menuItem("Define comparations",
    selectInput("batch_effect"," Batch Effect",list("empty")),          
    selectInput("gruped_variable"," Grouping variable",list("empty")),
    selectInput("variable","Select variable",list("empty")),
     
    radioButtons("group1", "Select group 1:",list("empty")),
    radioButtons("group2", "Select groups 2",list("empty")),
    
    numericInput("p_val","p value",value = 0.05,min = 0,max = 1,step = 0.01),
    numericInput("fc_val","Fold change",value = 1.5,step = 0.01)
     ),
     


 add_busy_spinner(spin = "semipolar",position="full-page" ),

hr(),

tags$head(tags$style(HTML( ".dataTables_scroll      {overflow-x:scroll;}"))),

tags$head(tags$style(HTML( ".dataTables_scrollBody  {overflow: unset !important;}"))),
tags$head(tags$style(HTML( ".dataTables_scrollHead  {
                                       overflow: unset !important;
                                       z-index: 10;
                             }"))),

div(style="display:inline;vertical-align:left;",
    actionButton("do_analysis","Execute analysis",
                 style = "color: #fff; background-color: #d98609; border-color: #fff;
                 padding: 5px 14px 5px 14px;margin: 5px 5px 5px 5px; ")
),
# div(style="display:inline;vertical-align:left;",
#     actionButton("target_analysis","Find targets",
#                  style = "color: #fff; background-color: #d98609; border-color: #fff;
#                  padding: 5px 14px 5px 14px;margin: 5px 5px 5px 5px; ")
# ),
# div(style="display:inline;vertical-align:left;",
#     actionButton("target_enrich","Get enrichments",
#                  style = "color: #fff; background-color: #d98609; border-color: #fff;
#                  padding: 5px 14px 5px 14px;margin: 5px 5px 5px 5px; ")
# ),
tags$style(type="text/css", "#report1 {background-color:#367fa9}"),
downloadButton("report1", class = "butt1",
               style = " border-color: #367fa9; padding: 5px 14px 5px 14px;margin: 5px 5px 5px 5px; ",
               label = "Generate report html",
               tags$head(tags$style(".butt1{background-color:white;}.butt1{color: black;}")))%>%
  add_class("report"),
downloadButton("report2", class = "butt1",
               style = " border-color: #367fa9; padding: 5px 14px 5px 14px;margin: 5px 5px 5px 5px; ",
               label = "Generate report doc",
               tags$head(tags$style(".butt1{background-color:white;}.butt1{color: black;}")))%>%
  add_class("report")

# )

,br(),
hr()
 
  ))
 


#body####

body<-  dashboardBody(
  tags$head(tags$style(HTML('
         .skin-blue .left-side, .skin-blue .wrapper {
                        background-color: #ecf0f5;
                        }
         '))),
  mainPanel(
    navbarPage(
      # position =  c("fixed-bottom"), 
                 "RESULTS", id = "tabs",
               ## Data introduced####
               tabPanel("Data introduced",
                        tabsetPanel(
                          
               tabPanel("Cel files RAW",
                        tabItem(tabName = "0",verticalLayout(
                                  box(title = "Table of actual celfiles introduced",
                                      "This table show the cel files upload",
                                      
                                      
                                      hr(),
                                      DTOutput("cel_in"),collapsible = TRUE,collapsed = F,solidHeader = TRUE,width = 12)))),
               tabPanel("Phenotype",
                        tabItem(tabName = "0",verticalLayout(
                          box(title = "Table with phenotype",
                              "This table show the phenotype data upload. Filter the subjects that you want to analyze",
                              hr(),
                              DTOutput("feno"),
                              collapsible = TRUE,collapsed = F,solidHeader = TRUE,width = 12)
                        ))),
               
               tabPanel("Your data to analyze",
                        tabItem(tabName = "0",verticalLayout(
                          box(title = "Cel files annotated ",
                              "This is table that is going to use to perform the  analysis",
                              hr(),
                              DTOutput("cel_files_feno"),
                              collapsible = TRUE,collapsed = F,solidHeader = TRUE,width = 12)
                        ))))),
               ## RAW data####
               tabPanel("Processed data",
                        tabsetPanel(
               tabPanel("RAW DATA",
                        tabItem(tabName = "0",verticalLayout(
                        box(title = "Raw Data ",
                   hr(h4("Dimension of arrays")),
                   "Raw data contains all probsets from cel files unproccessed",
                   hr(),
                   DTOutput("dim_raw"),
                   DTOutput("raw_data"),
                   
                   
                   hr(h4("Quality control of raw data")),
                   "To see the quality contol results, first choosse the analysis",
                   hr(),
                   DTOutput("raw_qc"),
                   hr(h4("Gene annotation detected")),
                   textOutput("annotation"),
                   # DTOutput("annotate_file"),
                   plotOutput("annotation_plot"),
                   collapsible = TRUE,collapsed = F,solidHeader = TRUE,width = 12),
    textAreaInput("input", "", "Escribe aqui tu comentario...", width = "1600px",resize = "vertical"),
  verbatimTextOutput("data_raw_text", placeholder = FALSE)))),
  ## NORM DATA data####
  tabPanel("Normalize DATA",
           tabItem(tabName = "0",verticalLayout(
             box(title = "Norm Data ",
                 hr(h4("Dimension of arrays")),
                 DTOutput("dim_norm"),
                 DTOutput("exprs_data"),
                 hr(h4("Quality control of norm data")),
                 
                 
                 DTOutput("data_norm_qc"),
                
                 
                 
                 collapsible = TRUE,collapsed = F,solidHeader = TRUE,width = 12),
             textAreaInput("input", "", "Escribe aqui tu comentario...", width = "1600px",resize = "vertical"),
             verbatimTextOutput("data_rma_text", placeholder = FALSE)))))),
  ## FOLD CHANGE ####
  tabPanel("Fold change",tabItem(tabName = "3",verticalLayout(
    # actionButton("execute_fc","Calculate Fold Change"),
    textOutput("a"),
    "This table show differentialy expressed miRNAs. Choose one to see the validated targets",
    hr(),
    DTOutput("data_norm_filt"),
    plotOutput("volvano"),
    
    
    plotOutput("heatmap"),
    box(DTOutput("tab_all"), solidHeader = T, collapsible = T, 
        title = "Table for Volcanoplot", status = "primary",width = 12),
    box(DTOutput("heatmap_df"), solidHeader = T, collapsible = T, 
        title = "Table for Heatmap", status = "primary",width = 12)
  ))),
  ## Mirna ####
  tabPanel("MiRNA Targets",
           tabItem(tabName = "",verticalLayout(
             # actionButton("targets","Get targets of miRNA"),
             "This table show validated targets for the miRNA DE selected. Select one to see the biological functions ",
             hr(),
             DTOutput("top_miRNAs")))),
  ## ENRICHMENT ####
  tabPanel("Biological significance",
           "This tables show the significative biological terms and pathways that are altered for the targets of selected miRNA ",
           hr(),
           # actionButton("enrich","Enrichment"),
           tabsetPanel( 
  tabPanel("GO ORA",
           tabItem(tabName = "",verticalLayout(
             textOutput("mirna_1"),
              DTOutput("GO_ORA_table"),
             plotOutput("GO_ORA_dotplot", click = "plot_click")
              
             )
             )
           ),
  tabPanel("WIKI ORA",
           tabItem(tabName = "",verticalLayout(
         
             DTOutput("wiki_ora_test"),
             plotOutput("WIKI_ORA_dotplot", click = "plot_click")
             ))),
  
  tabPanel("KEGG ORA",
           tabItem(tabName = "",verticalLayout(
             
             DTOutput("kegg_ora_test"),
             plotOutput("KEGG_ORA_dotplot", click = "plot_click")
             
             )))))
    )),
  tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "CSS.css")
  )
)

a<-mermaid("
graph TB
  A-->B
  A-->C
  C-->E
  B-->D
  C-->D
  D-->F
  E-->F
")

graph <-
  create_graph() %>%
  add_node(label = "Cel files",node_aes = node_aes(color = "darkgreen")) %>%
  add_node(label = "Phenotype",node_aes = node_aes(color = "darkgreen")) %>%
  add_node(label = "DATA") %>%
  add_node(label = "Normalize") %>%
  add_node(label = "Filter") %>% 
  add_node(label = "Annotation") %>% 
  add_node(label = "Fold Change") %>% 
  add_node(label = "Targets") %>% 
  add_node(label = "GO") %>%
  add_node(label = "KEGG") %>%
  add_node(label = "WIKI") %>% 

  add_edge(from = "Cel files",to = "DATA") %>% 
  add_edge(from = "Phenotype",to = "DATA") %>% 
  add_edge(from = "DATA",to = "Normalize") %>% 
  
  # add_edge(from = "Normalize",to = "Annotation") %>%
add_edge(from = "Annotation",to = "Normalize") %>%
# add_edge(from = "Normalize",to = "Annotation") %>% 
  add_edge(from = "Normalize",to = "Filter",edge_aes = edge_aes(color = "red",arrowhead = "vee",tooltip = "Red Arrow")) %>%
  add_edge(from = "Filter",to = "Fold Change",edge_aes = edge_aes(color = "red",arrowhead = "vee",tooltip = "Red Arrow")) %>% 
  add_edge(from = "Normalize",to = "Fold Change") %>%
  
  add_edge(from = "Fold Change",to = "Targets") %>%
  add_edge(from = "Targets",to = "GO") %>%
  add_edge(from = "Targets",to = "KEGG") %>%
  add_edge(from = "Targets",to = "WIKI")  %>% 


 render_graph(layout = "tree",width = 300,height = 700)
graph



ui <- dashboardPage (header,sidebar,  
                     controlbar = dashboardControlbar(width = 300,skin = "light", 
                                                      
                                                      graph
                        
                     ),body)

shinyUI(ui)