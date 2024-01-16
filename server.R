library(shinymanager)
library(multiMiR)
library(cowplot)
library(oligo)
library(shinyFiles)
library(data.table)

# library(arrayQualityMetrics)
options(shiny.maxRequestSize = 3000 * 1024^2)
dir_u<-"RESULTATS_USUARIS"
server<- (function(input, output, session) {



  
  
  logname<-reactive({
    logname<-Sys.getenv()
    logname<-data.frame(logname)
    logname<-logname["RSTUDIO_USER_IDENTITY",][[1]]
    if(dir.exists(logname)){
      print("Delete previous results?")
    }
    logname 
  })
  dir_usuari <-reactive({
    dir_u<-"RESULTATS_USUARIS"
    paste0(dir_u,"/resultats_",logname())})
  observe({
    
    print(session$clientData$url_hostname)
    url_pathname<-session$clientData$url_pathname
    
    print(url_pathname)
    dir_u<-"RESULTATS_USUARIS"
    dir.create(dir_u)
    dir <-paste0(dir_u,"/resultats_",logname())
    
    if (!dir.exists(dir)) {
      dir.create(dir)
    }
    
    dirplots <- file.path(dir, "Plots")
    if (!dir.exists(dirplots)) {
      dir.create(dirplots)
    }
    
    
    dirtaules <- file.path(dir, "Taules")
    if (!dir.exists(dirtaules)) {
      dir.create(dirtaules)
    }
    
    
  })  
  source("./functions/06_fun_tab.R", local = TRUE)
  source("./functions/01_fun_celfiles_in.R", local = TRUE)
  
  output$cel_in <- renderDT({
    if (is.null(input$file_input_celfiles)) {
      datatable_mod4(cel_files_in())
    } else {
      datatable_mod4(cel_files_in()[, c(1:2)])
    }
  })
  
  source("./functions/02_fun_feno_in.R", local = TRUE)
  
  observe({
    
    s <- input$cel_files_feno_rows_selected
    if (length(s)) {
      cel_files_feno <- cel_files_feno()[s, ]
    } else {
      cel_files_feno <- cel_files_feno()
    }
    
    x <- colnames(cel_files_feno)
    updateSelectInput(session, "batch_effect",
                      label = "Batch Effect",
                      choices = c("no", x),
                      selected = "no")
    
    updateSelectInput(session, "gruped_variable",
                      label = "Grouping variable",
                      choices = c("no", x),
                      selected = "no")
    updateSelectInput(session, "variable",
                      label = "Select variable",
                      choices = c("empty", x),
                      selected = "")
  })
  
  observe({
    if (input$variable != "empty") {
      s <- input$cel_files_feno_rows_selected
      if (length(s)) {
        cel_files_feno <- cel_files_feno()[s, ]
      } else {
        cel_files_feno <- cel_files_feno()
      }
      
      levels_var <- make.names(levels(factor(cel_files_feno[, input$variable])))
      updateRadioButtons(session, "group1",
                         label = "Select group 1",
                         choices = levels_var)
      updateRadioButtons(session, "group2",
                         label = "Select group 2",
                         choices = levels_var)
    }
  })
  
  output$feno <- renderDT(datatable_mod4(feno()))
  
  source("./functions/03_fun_celfiles_pheno.R", local = TRUE)
  output$cel_files_feno <- renderDataTable({
    datatable_mod4(cel_files_feno()[, -c(1:4)])
  })
  unlink("./Resultats/Plots/*")
  observeEvent(input$execute_pre, {
    
    source("./functions/05_fun_data_raw.R", local = TRUE)
    source("./functions/04_fun_annotation.R", local = TRUE)
    
    output$dim_raw <- renderDT(datatable(data.frame(samples=dim(data_raw()@phenoData@data)[1],
                                                         probeSets=dim(data_raw()@featureData@data)[1])))
   
     output$raw_data <- renderDT(datatable_mod4(data.frame(data_raw()@phenoData@data)))
    
    output$annotation <- renderText(data_raw()@annotation)
    
   output$annotation_plot <- renderPlot(annotation_plot())
  output$dim_norm <- renderDT(datatable(data.frame(samples=dim(data_norm()@phenoData@data)[1],
                                                          probeSets=dim(data_norm()@featureData@data)[1]
    )))
    outputOptions(output, "dim_raw", suspendWhenHidden = TRUE)
    outputOptions(output, "dim_norm", suspendWhenHidden = TRUE)
    output$exprs_data <- renderDT(datatable_mod4(exprs_data()))
    outputOptions(output, "annotation", suspendWhenHidden = FALSE)
    
    
    observeEvent(input$do_analysis, {
      
      unlink("./Resultats/Plots/*")
      if ("qc" %in% input$analisis) {
       source("./functions/07_fun_qc.R", local = TRUE)
      output$raw_qc <- renderDataTable(datatable_mod4(data_raw_qc(), "QC RAW"))
      output$data_norm <- renderDataTable(datatable_mod4(data_norm_qc(), "QC RAW"))
      outputOptions(output, "raw_qc", suspendWhenHidden = FALSE)
      outputOptions(output, "data_norm", suspendWhenHidden = FALSE)
      }
      if ("fc" %in% input$analisis) {
        source("./functions/08_fun_fc.R", local = TRUE)
        output$data_norm_filt <- renderDataTable({
          cols <- c("logFC", "AveExpr")
          cols1 <- c("P.Value", "adj.P.Val")
          FC_tab <- data_norm_filt() %>%
            filter(P.Value <= input$p_val, abs(logFC) >= log2(input$fc_val)) %>%
            mutate(across(cols, round, 3)) %>%
            mutate(across(all_of(cols1), format.pval))
          datatable_mod4(FC_tab, "FC")
        })
        output$volvano <- renderPlot(volcano())
        output$heatmap <- renderPlot(heatmap())
        outputOptions(output, "data_norm_filt", suspendWhenHidden = FALSE)
      }
      
      if ("target" %in% input$analisis) {
        source("./functions/09_fun_targets.R", local = TRUE)
        output$top_miRNAs <- renderDataTable({
          datatable_mod4(top_miRNAs())})
      } 
      # else {
      #   data <- data.frame(Error = "Seleccione análisis")
      #   output$top_miRNAs <- renderDataTable(data.frame(data, row.names = ""))
      #   unlink(file.path(dir, "Taules", "go_class.csv"), recursive = TRUE)
      # }
      
      if ("enrichment" %in% input$analisis) {
        source("./functions/10_fun_go_ora.R", local = TRUE)
        output$mirna_1 <- renderText(mirna_1())
        output$GO_ORA_table <- renderDataTable(datatable_mod4(go_ora_test(), col_inv = "geneID"))
        source("./functions/11_fun_wiki_ora.R", local = TRUE)
        output$wiki_ora_test <- renderDataTable(datatable_mod4(wiki_ora_test(), col_inv = "geneID"))
        source("./functions/12_fun_kegg_ora.R", local = TRUE)
        output$kegg_ora_test <- renderDataTable(datatable_mod4(kegg_ora_test(), col_inv = "geneID"))
       }
      # else {
      #   data <- data.frame(Error = "Seleccione análisis")
      #   output$top_miRNAs <- renderDataTable(data.frame(data, row.names = ""))
      #   unlink(file.path(dir, "Taules", "go_class.csv"), recursive = TRUE)
      # }
    })
  })
  
  output$report1 <- downloadHandler(
    filename = "report.html",
    content = function(file) {
      withProgress(
        message = "Report in progress",
        detail = "This may take a while...",
        value = 0,
        {
          for (i in 1:6) {
            incProgress(1 / 6)
          }
          fold_result<-paste0(dir_u,"/resultats_",logname())
          # tempReport <- file.path(tempdir(), "biagra.Rmd")
          # file.copy("biagra.Rmd", tempReport, overwrite = TRUE)
          rmarkdown::render("./mirna_limma.Rmd", output_file = file, envir = new.env(parent = globalenv()),params = list(logname=paste0(dir_u,"/resultats_",logname())))
           rmarkdown::render("./mirna_limma.Rmd", output_dir = paste0(dir_u,"/resultats_",logname()),envir = new.env(parent = globalenv()),params = list(logname=paste0(dir_u,"/resultats_",logname())))
        }
      )
    }
  )
  
  output$report2 <- downloadHandler(
    filename = "report.doc",
    content = function(file) {
      withProgress(
        message = "Report in progress",
        detail = "This may take a while...",
        value = 0,
        {
          for (i in 1:6) {
            incProgress(1 / 6)
          }
          # tempReport <- file.path(tempdir(), "biagra2.Rmd")
          # file.copy("biagra2.Rmd", tempReport, overwrite = TRUE)
          rmarkdown::render("./mirna_limma2.Rmd", output_file = file, envir = new.env(parent = globalenv()),params = list(logname=paste0(dir_u,"/resultats_",logname())))
          rmarkdown::render("./mirna_limma2.Rmd",output_dir = paste0(dir_u,"/resultats_",logname()), envir = new.env(parent = globalenv()),params = list(logname=paste0(dir_u,"/resultats_",logname())))
        }
      )
    }
  )
  

})
shinyServer(server)
