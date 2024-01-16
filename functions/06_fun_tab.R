# observe(
#   if(!"qc"%in%input$analisis){
#     hideTab(inputId = "tabs", target = "Quality control")
#   }    
# )
# observe(
#   if("qc"%in%input$analisis){

#     showTab(inputId = "tabs", target = "Quality control")
#   }    
# )
observe(
  if(!"fc"%in%input$analisis){
    hideTab(inputId = "tabs", target = "Fold change")
  }    
)
observe(
  if("fc"%in%input$analisis){
    showTab(inputId = "tabs", target = "Fold change")
  }    
)

observe(
  if(!"target"%in%input$analisis){
    hideTab(inputId = "tabs", target = "MiRNA Targets")
  }    
)
observe(
  if("target"%in%input$analisis){
    showTab(inputId = "tabs", target = "MiRNA Targets")
  }    
)


observe(
  if(!"enrichment"%in%input$analisis){
    hideTab(inputId = "tabs", target = "Biological significance")
  }
)
observe(
  if("enrichment"%in%input$analisis){
    showTab(inputId = "tabs", target = "Biological significance")
  }
)



