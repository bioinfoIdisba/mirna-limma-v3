

top_miRNAs<-
      # eventReactive(input$target_analysis,{   
   reactive({
  # observeEvent({
  # observe({
   
    withProgress(message = 'Targets',
                 detail = 'This may take a while...', value = 0, {
                   for (i in 1:6) {incProgress(1/6)}
                   # tempReport <- file.path(tempdir(), "biagra.Rmd")
                   # file.copy("biagra.Rmd", tempReport, overwrite = TRUE)
                   
  s<-input$data_norm_filt_rows_selected
  if(length(s)){
    FC_tab<-data_norm_filt() %>% filter(P.Value<=input$p_val,abs(logFC)>=log2(input$fc_val))
    example1 <- get_multimir(mirna = FC_tab$ID[s], summary = F,table = "validated")
    
  }
  write.csv(data.frame(example1@data %>%distinct(target_entrez,mature_mirna_id,.keep_all = T) ),"./Resultats/Taules/targets.csv")
  data.frame(example1@data %>%distinct(target_entrez,mature_mirna_id,.keep_all = T) )
  })
    # })
   })