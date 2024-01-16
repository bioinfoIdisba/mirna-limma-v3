feno<-reactive({
 if(!is.null(input$file_input_fenotip)){
   sep<-input$sep
   dec<-input$dec
   
feno_file<-input$file_input_fenotip$datapath
sep<-input$sep
dec<-input$dec

feno<-read.csv(feno_file,sep = sep,dec=dec)
feno<-data.frame(feno)
dim(feno)


return(feno)
}else{
  feno<-data.frame("No pheno data upload")
}
})

