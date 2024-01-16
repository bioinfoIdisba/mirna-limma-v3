mirna_1<-renderText({
  s<-input$data_norm_filt_rows_selected
  FC_tab<-data_norm_filt() %>% filter(P.Value<=input$p_val,abs(logFC)>=log2(input$fc_val))
  mirna = FC_tab$ID[s]
  return(mirna)
  
})
go_ora_test<-
  # eventReactive(input$do_analysis,{
    
  reactive({
    
    
  if("enrichment"%in%input$analisis){
  # OrgDb         = input$org
  # ont           = input$ont
  # pAdjustMethod = input$pAdjustMethod
  # pvalueCutoff  = input$pvalueCutoff
  # qvalueCutoff  = input$qvalueCutof
  # readable      = input$readable_GroupGO 
  # 
  OrgDb         = "org.Hs.eg.db"
  ont           = "BP"
  pAdjustMethod = "BH"
  pvalueCutoff  = 0.05
  qvalueCutoff  = 0.2
  readable      = T
  
  withProgress(message = 'GO ORA in progress',
               detail = 'This may take a while...', value = 0, {
                 for (i in 1:6) {
                   incProgress(1/6)
                   
                 }
                 
                 # gens<-gene_convert()
                 gens<-top_miRNAs()$target_entrez
                 # gens<-gens %>% 
                 #   filter(get(input$pvalue_column)<=input$pvalue_ora) %>% 
                 #   filter(abs(get(input$FC_column))>=input$fc_ora)
                 go_ora_test<-enrichGO(gene      = as.character(gens),
                                       # universe      = genes_ID_FC()[,1],
                                       # keyType       = input$keytype,
                                       OrgDb         = OrgDb,
                                       ont           = ont,
                                       pAdjustMethod = pAdjustMethod,
                                       pvalueCutoff  = pvalueCutoff,
                                       qvalueCutoff  = qvalueCutoff,
                                       readable      = readable)
                 if(dim(go_ora_test)[1]>0){
                   
                 go_ora_test_df<-data.frame(go_ora_test)
                 }else{
                   go_ora_test_df<-data.frame(go_ora_test)
                   # output$GO_ORA_dotplot<-renderPlot(plot.new())
                 }
                 data.frame(go_ora_test_df)
                 return(data.frame(go_ora_test_df))
                 
                 })
  }
})


  GO_ORA_dotplot<-reactive({
    go_ora_test<-data.frame(go_ora_test())
    go_ora_test$geneRatio<-unlist(lapply(go_ora_test$GeneRatio, function(x){ as.numeric(strsplit(x,"/")[[1]][1])/as.numeric(strsplit(x,"/")[[1]][2])} ))
    
    
      s_GO <- input$GO_ORA_table_rows_selected
      if(dim(go_ora_test)[1]>1){
      if(length (s_GO)<1){
        sel_GO<-1}
      else{
        sel_GO<-s_GO}
      }else{sel_GO<-1}
    ggplot(go_ora_test[sel_GO,], aes(x = geneRatio, y = reorder(Description, geneRatio),fill=p.adjust)) +
      geom_point(size = 5,aes(color=p.adjust)) +  # Use a larger dot
      geom_segment(aes(color=p.adjust,yend = reorder(Description, geneRatio)), xend = 0,linetype = "dashed") +
      theme_bw() +
      theme(axis.text.y = element_text(hjust=1,size = 15),
            axis.text.x = element_text(angle=45, hjust=1,size = 15))+ 
      
      theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_line( linetype = "dashed"))+
      scale_y_discrete(labels = function(x) str_trunc(x,width = 50))
    })                             
  output$GO_ORA_dotplot<-renderPlot({
    print(GO_ORA_dotplot())
    ggsave("./Resultats/Plots/go_ora_dotplot.png")
  })
  
  

  
  



