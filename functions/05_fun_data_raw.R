data_raw<-eventReactive(input$execute_pre,{
  # source("./functions/feno_in.R",local = T)
 
   s<-input$cel_files_feno_rows_selected
   cel_files_feno_all<-cel_files_feno()
  if(length(s)){
    cel_files_feno<-cel_files_feno()[s,]
  }else{
    cel_files_feno<-cel_files_feno()
  }
  
  pheno<-cel_files_feno_all
  
   # rownames(feno)<-feno[,1]
  # pheno_ano<-AnnotatedDataFrame(feno)
   rownames(pheno)<-pheno[,1]
   
   pheno_ano<-AnnotatedDataFrame(pheno)
  # 
    # cel_files_feno<-list.celfiles("./amanda/Amanda_miRNA_serum_cel files/",full.name=T)
   # data_raw<-read.celfiles(cel_files_feno,phenoData = pheno_ano)
  data_raw<-read.celfiles(cel_files_feno_all$datapath, phenoData = pheno_ano)
  
  data_raw[,colnames(data_raw)%in%cel_files_feno[,1]]
})
data_norm<-reactive({
  # source("./functions/feno_in.R",local = T)
  # 
    # data_norm<-rma(data_raw)
   # feat<-annotate_file[match(rownames(data_norm@featureData@data),annotate_file$Probe.Set.Name),]
  data_norm<-oligo::rma(data_raw())
  
  if(data_raw()@annotation=="pd.mirna.4.0"){
  feat<-annotate_file()[match(rownames(data_norm@featureData@data),annotate_file()$Probe.Set.Name),]
  }
  
  if(data_raw()@annotation=="pd.clariom.s.human"|data_raw()@annotation=="pd.hugene.1.0.st.v1"){
    feat<-annotate_file()[match(rownames(data_norm@featureData@data),annotate_file()$transcript_cluster_id),]
    # feat<-feat[!is.na(feat$transcript_cluster_id),]
    
  }
  
  
  feat<-AnnotatedDataFrame(feat)
  rownames(feat@data)<-
    rownames(data_norm@featureData)
  data_norm@featureData<-feat
  dup.ids <- feat@data$Accession[duplicated(feat@data$Accession)] %>% 
    unique %>%
    sort
  
  
  data_norm
  
})

exprs_data<-reactive({
  oligo::exprs(data_norm())
  
})
