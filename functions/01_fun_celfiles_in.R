cel_files_in<-reactive({
  
  
  if(!is.null(input$file_input_celfiles)){
    cel_files_in<-input$file_input_celfiles
    
    cel_files_in<-data.frame(cel_files_in)
  }
  else{
    cel_files_in<-data.frame("No CEL files upload")
    
  }
  
  return(cel_files_in)
  
})
