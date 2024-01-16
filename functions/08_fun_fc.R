
data_norm_filt<-
   eventReactive(input$do_analysis,{
  # reactive({
  g1<-input$group1
  g2<-input$group2
  columna<-input$variable
  data_norm<-data_norm()
  s<-input$cel_files_feno_rows_selected
  # if(length(s)){
  # data_norm_filt<-data_norm[,s]
  # }else{
  #   data_norm_filt<-data_norm
  # }
  data_norm_filt<-data_norm
  data_norm_filt<-data_norm_filt[,pData(data_norm_filt)[,columna]%in%c(g1,g2)]
  
if(data_raw()@annotation=="pd.mirna.4.0"){
data_norm_filt<-
    data_norm_filt[fData(data_norm_filt)$Species.Scientific.Name=="Homo sapiens",]
data_norm_filt<-data_norm_filt[data_norm_filt@featureData@data$Sequence.Type=="miRNA",]
}
  
if(data_raw()@annotation=="pd.clariom.s.human"|data_raw()@annotation=="pd.hugene.1.0.st.v1"){
     # data_norm_filt<-data_norm_filt[!is.na(fData(data_norm_filt)$gene_entrez),]
     data_norm_filt<-data_norm_filt[!is.na(data_norm_filt@featureData@data$gene_entrez),]
     
     print(data_norm_filt@featureData@data$gene_entrez)
}
  
data_norm_filt<-data_norm_filt
exp_rma <- exprs(data_norm_filt)

# Non Paired ####
if(input$batch_effect=="no"){
if(input$gruped_variable=="no"){
variables = make.names(pData(data_norm_filt)[,columna])

# variables = (pData(data_norm_filt)[,columna])
f = factor(variables,levels = c(g2,g1))

design = model.matrix(~ 0 + f)

 data.fit = lmFit(exp_rma,design)
 a<-c(colnames(design)[1],"-",colnames(design)[2])
 b<-c(colnames(design)[1],"vs",colnames(design)[2])
 output$a<-renderText(b)
 astr=paste(a, collapse="")
 prestr="makeContrasts("
 poststr=",levels=design)"
  commandstr=paste(prestr,astr,poststr,sep="")
contrast.matrix=eval(parse(text=commandstr))

data.fit.con = contrasts.fit(data.fit,contrast.matrix)
data.fit.eb = eBayes(data.fit.con)



if(data_raw()@annotation=="pd.mirna.4.0"){
  tab = topTable(data.fit.eb,
                 # coef=length(colnames(design)),
                 lfc = log2(input$fc_val),
                 number=Inf,
                 adjust.method="BH",
                 genelist = data_norm_filt@featureData@data$Transcript.ID.Array.Design. )
  tab_all = topTable(data.fit.eb,
                     # coef=length(colnames(design)),
                     number=Inf,
                     adjust.method="BH",
                     genelist = data_norm_filt@featureData@data$Transcript.ID.Array.Design.)
}

if(data_raw()@annotation=="pd.clariom.s.human"|data_raw()@annotation=="pd.hugene.1.0.st.v1"){
  tab = topTable(data.fit.eb,
                 # coef=length(colnames(design)),
                 lfc = log2(input$fc_val),
                 number=Inf,
                 adjust.method="BH",
                 genelist = data.frame(symbol=data_norm_filt@featureData@data$gene_symbol,entrez=data_norm_filt@featureData@data$gene_entrez)  )
  tab_all = topTable(data.fit.eb,
                     # coef=length(colnames(design)),
                     number=Inf,
                     adjust.method="BH",
                     genelist = data.frame(symbol=data_norm_filt@featureData@data$gene_symbol,entrez=data_norm_filt@featureData@data$gene_entrez))
}
cols<-c("logFC","AveExpr")
cols1<-c("P.Value", "adj.P.Val")
#write.csv(tab_all %>%  filter(P.Value<=input$p_val,abs(logFC)>=log2(input$fc_val)),"./Resultats/Taules/tab_all.csv")
tab_all

}else{
# Paired ####
  # columna_grup<-"pacient"
  columna_grup<-input$gruped_variable
  
  
  variables = make.names(pData(data_norm_filt)[,columna])

  # variables = (pData(data_norm_filt)[,columna])
  variable_group = make.names(pData(data_norm_filt)[,columna_grup])
  f = factor(variables,levels = c(g2,g1))
  fg = factor(variable_group)
  paired.design = model.matrix(~  fg+f)
  
  # colnames(design) = c(g1,g2)
  
  # # 
  # # lev <- c(g1_mod, g2_mod)
  # 
  # 
  # 
  # # Parsing
  data.fit = lmFit(exp_rma,paired.design)
  data.fit$coefficients
  # # 
  # 
  # 
  # 
  # 
  data.fit.eb = eBayes(data.fit)
  
  if(data_raw()@annotation=="pd.mirna.4.0"){
  tab = topTable(data.fit.eb,coef=length(colnames(paired.design)),
                 lfc = log2(input$fc_val),number=Inf,adjust.method="BH",
                 genelist = data_norm_filt@featureData@data$Transcript.ID.Array.Design. )
  tab_all = topTable(data.fit.eb,coef=length(colnames(paired.design)),number=Inf,adjust.method="BH",genelist = data_norm_filt@featureData@data$Transcript.ID.Array.Design.)
  }
  
  if(data_raw()@annotation=="pd.clariom.s.human"|data_raw()@annotation=="pd.hugene.1.0.st.v1"){
    tab = topTable(data.fit.eb,coef=length(colnames(paired.design)),
                   lfc = log2(input$fc_val),number=Inf,adjust.method="BH",
                   genelist = data_norm_filt@featureData@data$gene_entrez )
    tab_all = topTable(data.fit.eb,coef=length(colnames(paired.design)),number=Inf,adjust.method="BH",genelist = data_norm_filt@featureData@data$gene_entrez)
  }
  cols<-c("logFC","AveExpr")
  cols1<-c("P.Value", "adj.P.Val")
  # 
  # 
  # 
  # FC_tab<-tab_all %>% filter(P.Value<=0.05,abs(logFC)>=log2(1.5)) %>%
    
    # dplyr::select(-c(GeneChip.Array,Annotation.Date,Sequence,Sequence.Source,Probe.Set.ID,B,t,Probe.Set.Name,Alignments,Clustered.miRNAs.within.10kb,Genome.Context,Target.Genes)) %>%
    # mutate(across(cols, round, 3)) %>%
    # mutate(across(all_of(cols1), format.pval))
    # 
  #write.csv(tab_all %>%  filter(P.Value<=input$p_val,abs(logFC)>=log2(input$fc_val)),"./Resultats/Taules/tab_all.csv")
  tab_all
}
}else{
  
  if(input$gruped_variable=="no"){
    variables = make.names(pData(data_norm_filt)[,columna])
    
    # variables = (pData(data_norm_filt)[,columna])
    f = factor(variables,levels = c(g2,g1))
    
    design = model.matrix(~ 0 + f+input$batch_effect)
    
    data.fit = lmFit(exp_rma,design)
    a<-c(colnames(design)[1],"-",colnames(design)[2])
    b<-c(colnames(design)[1],"vs",colnames(design)[2])
    output$a<-renderText(b)
    astr=paste(a, collapse="")
    prestr="makeContrasts("
    poststr=",levels=design)"
    commandstr=paste(prestr,astr,poststr,sep="")
    contrast.matrix=eval(parse(text=commandstr))
    
    data.fit.con = contrasts.fit(data.fit,contrast.matrix)
    data.fit.eb = eBayes(data.fit.con)
    
    
    
    if(data_raw()@annotation=="pd.mirna.4.0"){
      tab = topTable(data.fit.eb,
                     # coef=length(colnames(design)),
                     lfc = log2(input$fc_val),
                     number=Inf,
                     adjust.method="BH",
                     genelist = data_norm_filt@featureData@data$Transcript.ID.Array.Design. )
      tab_all = topTable(data.fit.eb,
                         # coef=length(colnames(design)),
                         number=Inf,
                         adjust.method="BH",
                         genelist = data_norm_filt@featureData@data$Transcript.ID.Array.Design.)
    }
    
    if(data_raw()@annotation=="pd.clariom.s.human"|data_raw()@annotation=="pd.hugene.1.0.st.v1"){
      tab = topTable(data.fit.eb,
                     # coef=length(colnames(design)),
                     lfc = log2(input$fc_val),
                     number=Inf,
                     adjust.method="BH",
                     genelist = data.frame(symbol=data_norm_filt@featureData@data$gene_symbol,entrez=data_norm_filt@featureData@data$gene_entrez)  )
      tab_all = topTable(data.fit.eb,
                         # coef=length(colnames(design)),
                         number=Inf,
                         adjust.method="BH",
                         genelist = data.frame(symbol=data_norm_filt@featureData@data$gene_symbol,entrez=data_norm_filt@featureData@data$gene_entrez))
    }
    cols<-c("logFC","AveExpr")
    cols1<-c("P.Value", "adj.P.Val")
    #write.csv(tab_all %>%  filter(P.Value<=input$p_val,abs(logFC)>=log2(input$fc_val)),"./Resultats/Taules/tab_all.csv")
    tab_all
    
  }else{
    # Paired ####
    # columna_grup<-"pacient"
    columna_grup<-input$gruped_variable
    
    
    variables = make.names(pData(data_norm_filt)[,columna])
    
    # variables = (pData(data_norm_filt)[,columna])
    variable_group = make.names(pData(data_norm_filt)[,columna_grup])
    f = factor(variables,levels = c(g2,g1))
    fg = factor(variable_group)
    paired.design = model.matrix(~  fg+f+input$batch_effect)
    
    # colnames(design) = c(g1,g2)
    
    # # 
    # # lev <- c(g1_mod, g2_mod)
    # 
    # 
    # 
    # # Parsing
    data.fit = lmFit(exp_rma,paired.design)
    data.fit$coefficients
    # # 
    # 
    # 
    # 
    # 
    data.fit.eb = eBayes(data.fit)
    
    if(data_raw()@annotation=="pd.mirna.4.0"){
      tab = topTable(data.fit.eb,coef=length(colnames(paired.design)),
                     lfc = log2(input$fc_val),number=Inf,adjust.method="BH",
                     genelist = data_norm_filt@featureData@data$Transcript.ID.Array.Design. )
      tab_all = topTable(data.fit.eb,coef=length(colnames(paired.design)),number=Inf,adjust.method="BH",genelist = data_norm_filt@featureData@data$Transcript.ID.Array.Design.)
    }
    
    if(data_raw()@annotation=="pd.clariom.s.human"|data_raw()@annotation=="pd.hugene.1.0.st.v1"){
      tab = topTable(data.fit.eb,coef=length(colnames(paired.design)),
                     lfc = log2(input$fc_val),number=Inf,adjust.method="BH",
                     genelist = data_norm_filt@featureData@data$gene_entrez )
      tab_all = topTable(data.fit.eb,coef=length(colnames(paired.design)),number=Inf,adjust.method="BH",genelist = data_norm_filt@featureData@data$gene_entrez)
    }
    cols<-c("logFC","AveExpr")
    cols1<-c("P.Value", "adj.P.Val")
    # 
    # 
    # 
    # FC_tab<-tab_all %>% filter(P.Value<=0.05,abs(logFC)>=log2(1.5)) %>%
    
    # dplyr::select(-c(GeneChip.Array,Annotation.Date,Sequence,Sequence.Source,Probe.Set.ID,B,t,Probe.Set.Name,Alignments,Clustered.miRNAs.within.10kb,Genome.Context,Target.Genes)) %>%
    # mutate(across(cols, round, 3)) %>%
    # mutate(across(all_of(cols1), format.pval))
    # 
    #write.csv(tab_all %>%  filter(P.Value<=input$p_val,abs(logFC)>=log2(input$fc_val)),"./Resultats/Taules/tab_all.csv")
    tab_all
  }
  
}
# data.frame(c(g1,g2))
# data.frame(variables)
})




# volcano<-eventReactive(input$do_analysis,{
volcano<-reactive({
  tab_all<-data_norm_filt()
  tab_all$diffexpressed <- "NO"
  # # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
  tab_all$diffexpressed[tab_all$logFC > log2(input$fc_val) & tab_all$P.Value < input$p_val] <- "UP"
  # # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  tab_all$diffexpressed[tab_all$logFC < (-log2(input$fc_val)) & tab_all$P.Value < input$p_val] <- "DOWN"


  # tab_all$delabel <- NA
  # tab_all$delabel[tab_all$diffexpressed != "NO"] <- tab_all$gene_symbol[tab_all$diffexpressed != "NO"]
# plot.new()
  # plot adding up all layers we have seen so far
  if(data_raw()@annotation=="pd.mirna.4.0"){
    lab<-tab_all$ID
  }  
  if(data_raw()@annotation=="pd.clariom.s.human"|data_raw()@annotation=="pd.hugene.1.0.st.v1"){
    lab<-tab_all$symbol
  } 
  output$tab_all<-renderDataTable(datatable_mod4(tab_all))
  
  p<-EnhancedVolcano(tab_all,
                  lab = lab,
                  # lab = rownames(tab_all),
                  x = 'logFC',
                  y = 'P.Value',
                  pCutoff = input$p_val,
                  FCcutoff = input$fc_val,
                  pointSize = 2.0,
                  labSize = 4.0,
                  labCol = 'black',
                  labFace = 'bold',

                  # boxedLabels = TRUE,
                  colAlpha = 4/5,
                  # legendPosition = 'right',
                  legendLabSize = 14,
                  legendIconSize = 4.0,
                  drawConnectors = TRUE,
                  widthConnectors = 1.0,
                  colConnectors = 'black',
                  ylim = c(0, max(-log10(tab_all$P.Value), na.rm = TRUE) + 0),
                  xlim = c(min(tab_all$logFC, na.rm = TRUE) - 0, max(tab_all$logFC, na.rm = TRUE)+0))
#                                    # encircle =rownames(FC_tab),
# encircleCol = 'black',
# encircleSize = 2.5,
# encircleFill = 'pink',
# encircleAlpha = 1/2
#   )
  
  png(paste0(dir_usuari(),"/Plots/volcano.png"))
  print(p)
  dev.off()
  p
})

# heatmap<-eventReactive(input$do_analysis,{
heatmap<-reactive({
  if(data_raw()@annotation=="pd.mirna.4.0"){
  tab_all<-data_norm_filt()
  data_norm<-data_norm()
  gens_SIG<-head(tab_all,25)$ID
  columna<-input$variable
  gens_SIG_1<-(tab_all[tab_all$diffexpressed!="NO",])$gene_symbol
  
  if(length(gens_SIG_1)>=length(gens_SIG)){# gens_SIG<-gens_SIG_1
    gens_SIG<-gens_SIG}
  
  data_norm_SIG<-data_norm[data_norm@featureData@data$Transcript.ID.Array.Design.%in%gens_SIG,]
  data_norm_SIG_df<-exprs(data_norm_SIG)
  
  rownames(data_norm_SIG_df)<-data_norm_SIG@featureData@data$Transcript.ID.Array.Design.
  # colnames(data_norm_SIG_df)<-data_norm_SIG@phenoData@data$pacient
  hmcol <- rev(colorRampPalette(RColorBrewer::brewer.pal(5, "YlOrRd"))(255))
  # colnames(dists) <- NULL
  # diag(dists) <- NA
  
  ann_colors <- list(
    Phenotype = c(g1 = "chartreuse4", g2 = "burlywood3"))
  }
  
  
  if(data_raw()@annotation=="pd.clariom.s.human"|data_raw()@annotation=="pd.hugene.1.0.st.v1"){
    
    tab_all<-data_norm_filt()
    data_norm<-data_norm()
    gens_SIG<-head(tab_all,25)$symbol
     columna<-input$variable
    # columna<-"visita"
    gens_SIG_1<-(tab_all[tab_all$diffexpressed!="NO",])$symbol
    
    if(length(gens_SIG_1)>=length(gens_SIG)){# gens_SIG<-gens_SIG_1
      gens_SIG<-gens_SIG}
    
    data_norm_SIG<-data_norm[data_norm@featureData@data$gene_symbol%in%gens_SIG,]
    data_norm_SIG_df<-exprs(data_norm_SIG)
    
    
    rownames(data_norm_SIG_df)<-data_norm_SIG@featureData@data$gene_symbol
    colnames(data_norm_SIG_df)<-data_norm_SIG@phenoData@data$pacient
    
    hmcol <- rev(colorRampPalette(RColorBrewer::brewer.pal(5, "YlOrRd"))(255))
    # colnames(dists) <- NULL
    # diag(dists) <- NA
    
    ann_colors <- list(
      Phenotype = c(g1 = "chartreuse4", g2 = "burlywood3"))
    
  }
  # names(ann_colors$Phenotype)<-c(g1_mod,g2_mod)
  phenotype_names <- pData(data_norm_SIG)[,columna]
  
  annotation_for_heatmap <- data.frame(Phenotype = phenotype_names)
  # annotation_for_heatmap<-as.factor(annotation_for_heatmap[,1])
  
  # row.names(annotation_for_heatmap) <- pData(data_norm)$pacient
  # rownames(annotation_for_heatmap) <- make.names(colnames(data_norm_SIG_df),unique = T)
  rownames(annotation_for_heatmap) <- colnames(data_norm_SIG_df)
  
  output$heatmap_df<-renderDataTable(datatable_mod4(data_norm_SIG_df))
  p1<-pheatmap::pheatmap(data_norm_SIG_df,
                         labels_col =   pData(data_norm)$pacient,
                            
               cutree_rows = 2,
               cutree_cols = 2,
               
                annotation_col = annotation_for_heatmap,
               # annotation_colors = ann_colors,
               scale = "row",
               legend = F,
               treeheight_row = 50,
               # legend_breaks = c(min(dists, na.rm = TRUE),
               #                   max(dists, na.rm = TRUE)),
               # legend_labels = (c("small distance", "large distance")),
               main = "top DE")
  
  png(paste0(dir_usuari(),"/Plots/heat.png"))
  # png("./Resultats/Plots/heat.png")
  print(p1)
  dev.off()
  
  p1
  
  
  
  
})
