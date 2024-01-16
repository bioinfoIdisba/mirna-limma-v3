cel_files_feno<-reactive({
  if(!is.null(input$file_input_fenotip)){
   cel_files_feno<-
    data.frame(cbind(cel_files_in(),feno()[match(cel_files_in()$name,feno()[,1]),]))
  
  }else{
    cel_files_feno<-data.frame("No pheno data upload")
  }
  # cel_files_feno<-    cel_files_feno[!is.na(cel_files_feno[,1]),]
  # cel_files_feno<-    data.frame(a=1)
  cel_files_feno
   
  return(cel_files_feno)
})