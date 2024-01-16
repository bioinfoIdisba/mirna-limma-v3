mirna_1<-renderText({
  s<-input$data_norm_filt_rows_selected
  FC_tab<-data_norm_filt() %>% filter(P.Value<=input$p_val,abs(logFC)>=log2(input$fc_val))
  mirna = FC_tab$ID[s]
  return(mirna)
  
})
wiki_ora_test<-
  # eventReactive(input$do_analysis,{
  
  reactive({
    
    
    if("enrichment"%in%input$analisis){
      # OrgDb         = input$org
      # ont           = input$ont
      # pAdjustMethod = input$pAdjustMethod
      # pvalueCutoff  = input$pvalueCutoff
      # qvalueCutoff  = input$qvalueCutof
      # readable      = input$readable_GroupWIKI 
      organism         = "Homo sapiens"
      pAdjustMethod = "BH"
      pvalueCutoff  = 0.05
      qvalueCutoff  = 0.2
      readable      = T
      
      
      
      withProgress(message = 'WIKI ORA in progress',
                   detail = 'This may take a while...', value = 0, {
                     for (i in 1:6) {
                       incProgress(1/6)
                       
                     }
                     
                     # gens<-gene_convert()
                     gens<-top_miRNAs()$target_entrez
                     # gens<-gens %>% 
                     #   filter(get(input$pvalue_column)<=input$pvalue_ora) %>% 
                     #   filter(abs(get(input$FC_column))>=input$fc_ora)
          
                     
                     wiki_ora_test<-enrichWP( gene =as.character(gens),
                                              pAdjustMethod = pAdjustMethod,
                                              # minGSSize    =  input$minGSSize,
                                              # maxGSSize    = input$maxGSSize,
                                              # qvalueCutoff  = input$qvalueCutoff,
                                              pvalueCutoff = pvalueCutoff,
                                              organism = organism)
                     if(dim(wiki_ora_test)[1]>0){
                       
                       wiki_ora_test_df<-data.frame(wiki_ora_test)
                     }else{
                       wiki_ora_test_df<-data.frame(wiki_ora_test)
                       # output$WIKI_ORA_dotplot<-renderPlot(plot.new())
                     }
                     data.frame(wiki_ora_test_df)
                     return(data.frame(wiki_ora_test_df))
                     
                   })
    }
  })


WIKI_ORA_dotplot<-reactive({
  wiki_ora_test<-data.frame(wiki_ora_test())
  wiki_ora_test$geneRatio<-unlist(lapply(wiki_ora_test$GeneRatio, function(x){ as.numeric(strsplit(x,"/")[[1]][1])/as.numeric(strsplit(x,"/")[[1]][2])} ))
  
  
  s_WIKI <- input$wiki_ora_test_rows_selected
  if(length (s_WIKI)<1){
    sel_WIKI<-1:10}
  else{
    sel_WIKI<-s_WIKI}
  ggplot(wiki_ora_test[sel_WIKI,], aes(x = geneRatio, y = reorder(Description, geneRatio))) +
    geom_point(size = 5) +  # Use a larger dot
    geom_point(size = 5,aes(color=p.adjust)) +  # Use a larger dot
    geom_segment(aes(color=p.adjust,yend = reorder(Description, geneRatio)), xend = 0,linetype = "dashed") +
    theme_bw() +
    theme(axis.text.y = element_text(hjust=1,size = 15),
          axis.text.x = element_text(angle=45, hjust=1,size = 15))+ 
    
    theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_line( linetype = "dashed"))+
    scale_y_discrete(labels = function(x) str_trunc(x,width = 50))
})                                  
output$WIKI_ORA_dotplot<-renderPlot({
  print(WIKI_ORA_dotplot())
  ggsave("./Resultats/Plots/wiki_ora_dotplot.png")
})








