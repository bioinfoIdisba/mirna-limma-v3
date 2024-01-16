
annotate_file<-eventReactive(input$execute_pre,{
  
   # anno<-"pd.mirna.4.0"
anno<-data_raw()@annotation

files_grep<-list.files("/home/SHARED/annotation_affymetrix/",full.names=T)
# files_grep<-list.files("./anotacions/",full.names=T)
annotate_file<-read.csv(files_grep[grepl(anno,files_grep)],sep=",",comment.char = "#")
annotate_file<-annotate_file[,!grepl("Target.Genes",colnames(annotate_file))]
annotate_file<-annotate_file[,!grepl("Genome.Context",colnames(annotate_file))]
annotate_file<-annotate_file[,!grepl("Clustered.miRNAs.within.10kb",colnames(annotate_file))]

if(anno=="pd.mirna.4.0"){
annotate_file<-
  annotate_file %>% 
  group_by(Accession) %>% 
  mutate(gene=paste(Accession, collapse="|"))
}

if(anno=="pd.clariom.s.human"|anno=="pd.hugene.1.0.st.v1"){
  
annotate_file$gene_symbol<-
  unlist(lapply(annotate_file$gene_assignment, function(x) gsub(" ","",strsplit(x,"//")[[1]][2])))
annotate_file$gene_entrez<-
  unlist(lapply(annotate_file$gene_assignment, function(x) gsub(" ","",strsplit(x,"//")[[1]][5])))

annotate_file$gene_cytoband<-
  unlist(lapply(annotate_file$gene_assignment, function(x) gsub(" ","",strsplit(x,"//")[[1]][4])))

}

print(annotate_file$gene_symbol)
return(data.frame(annotate_file))
})


annotation_plot<-eventReactive(input$execute_pre,{
 
annotate_file<-annotate_file()


if(data_raw()@annotation=="pd.clariom.s.human"|data_raw()@annotation=="pd.hugene.1.0.st.v1"){
  
 
  # data_rma<-oligo::rma(data_raw)
  ### Plot by category####
  data<-data.frame(table(annotate_file$category))
  colnames(data)<-c("category","count")
  data$category<-as.character(data$category)
  data$category[grepl("control",data$category)]<-"control"
  data<-
    data %>% group_by(category) %>%
    mutate(count=sum(count)) %>%
    ungroup() %>%
    distinct(category,.keep_all = T)
  data$fraction <- round(data$count / sum(data$count),3)*100
  data$ymax <- cumsum(data$fraction)
  data$ymin <- c(0, head(data$ymax, n=-1))
  data$labelPosition <- (data$ymax + data$ymin) / 2
  data$label <- paste0(data$category, ":\n ", data$count)
  p1<-
    ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
    geom_rect()+
    geom_text_repel( x=3.5, aes(y=labelPosition, label=paste0(fraction,"%")),
                     size=6,box.padding = 2,
                     point.padding = 0.2,
                     nudge_x = 0,
                     nudge_y = 0,
                     segment.curvature = -1e-20,
                     arrow = arrow(length = unit(0.015, "npc"))) +
    scale_fill_discrete() +
    coord_polar(theta="y") +
    xlim(c(2, 4)) +
    theme_void() +
    theme(legend.position = "left")+
    
    ggtitle(paste0("Probe sets in array by category (", sum(data$count)," probes)")) + guides(fill=guide_legend(title=""))
  
  
  data<-data.frame(table(annotate_file$seqname))
  colnames(data)<-c("category","count")
  data$category<-as.character(data$category)
  cromosomes<-(unique(data$category))
  cromosomes_1<-cromosomes[sapply(strsplit(cromosomes,""),length)==5]
  cromosomes_2<-cromosomes[sapply(strsplit(cromosomes,""),length)==4]
  cromosomes_3<-cromosomes[sapply(strsplit(cromosomes,""),length)==3]
  cromosomes<-c(cromosomes_1,cromosomes_2,cromosomes_3)
  data<-
    data %>% mutate(category=ifelse(category%in%cromosomes,category,"other"))
  data$category<-as.factor(data$category)
  data$category<-factor(data$category,levels = c("chr1",
                                                 "chr2",
                                                 "chr3",
                                                 "chr4",
                                                 "chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14",
                                                 "chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrM","chrUn","other","---"))
  data<-data[order(data$category),]
  data<-
    data %>% group_by(category) %>% 
    mutate(count=sum(count)) %>% 
    ungroup() %>% 
    distinct(category,.keep_all = T)
  data$fraction <- round(data$count / sum(data$count),3)*100
  data$ymax <- cumsum(data$fraction)
  data$ymin <- c(0, head(data$ymax, n=-1))
  data$labelPosition <- (data$ymax + data$ymin) / 2
  data$label <- paste0(data$category, ":\n ", data$count)
  p2<-
    data %>%
    # mutate(category = sort(category)) %>% 
    
    ggplot( aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=((category)))) +
    geom_rect()+
    geom_text_repel( x=3.5, aes(y=labelPosition, label=paste0(fraction,"%\n",category)), 
                     size=3,box.padding = 0.5,
                     point.padding = 2, 
                     nudge_x = .15,
                     nudge_y = 0,
                     segment.curvature = -1e-20,
                     arrow = arrow(length = unit(0.015, "npc"))) +
    scale_fill_discrete() +
    coord_polar(theta="y") +
    xlim(c(2, 4)) +
    theme_void() +
    theme(legend.position = "left")+
    ggtitle(paste0("Probe sets in array by locus type (", sum(data$count)," probes)")) +
    guides(fill=guide_legend(title=""))
  
  
  
  
  print( plot_grid(p1, p2))
}
if(data_raw()@annotation=="pd.mirna.4.0"){
# colnames(annotate)  
  annotate_file$chr<-unlist(lapply(annotate_file$Alignments, function(x){strsplit(x,":")[[1]][1]} ))
  data<-data.frame(table(annotate_file$Species.Scientific.Name))
  colnames(data)<-c("category","count")
  data$category<-as.character(data$category)
  
  data<-
    data %>% group_by(category) %>%
    mutate(count=sum(count)) %>%
    ungroup() %>%
    distinct(category,.keep_all = T)
  data<-data[order(data$count,decreasing = T),]
  data<-head(data,n=10)
  data$fraction <- round(data$count / sum(data$count),3)*100
  data$ymax <- cumsum(data$fraction)
  data$ymin <- c(0, head(data$ymax, n=-1))
  data$labelPosition <- (data$ymax + data$ymin) / 2
  data$label <- paste0(data$category, ":\n ", data$count)
  
  
  p1<-
    ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
    geom_rect()+
    geom_text_repel( x=3.5, aes(y=labelPosition, label=paste0(fraction,"%\n",category)),
                     size=3,box.padding = 1,
                     point.padding = 0.2,
                     nudge_x = 0,
                     nudge_y = 0,
                     segment.curvature = -1e-20,
                     arrow = arrow(length = unit(0.015, "npc"))) +
    scale_fill_discrete() +
    coord_polar(theta="y") +
    xlim(c(2, 4)) +
    theme_void() +
    theme(legend.position = "left")+
    
    ggtitle(paste0("Probe sets in array by category (", sum(data$count)," probes)")) + guides(fill=guide_legend(title=""))
  
  
  data<-data.frame(table(annotate_file$Sequence.Type))
  colnames(data)<-c("category","count")
  data$category<-as.character(data$category)
  
  data<-
    data %>% group_by(category) %>%
    mutate(count=sum(count)) %>%
    ungroup() %>%
    distinct(category,.keep_all = T)
  data<-data[order(data$count,decreasing = T),]
  data<-head(data,n=25)
  data$fraction <- round(data$count / sum(data$count),3)*100
  data$ymax <- cumsum(data$fraction)
  data$ymin <- c(0, head(data$ymax, n=-1))
  data$labelPosition <- (data$ymax + data$ymin) / 2
  data$label <- paste0(data$category, ":\n ", data$count)
  
  
  p2<-
    ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
    geom_rect()+
    geom_text_repel( x=3.5, aes(y=labelPosition, label=paste0(fraction,"%\n",category)),
                     size=3,box.padding = 1,
                     point.padding = 0.2,
                     nudge_x = 0,
                     nudge_y = 0,
                     segment.curvature = -1e-20,
                     arrow = arrow(length = unit(0.015, "npc"))) +
    scale_fill_discrete() +
    coord_polar(theta="y") +
    xlim(c(2, 4)) +
    theme_void() +
    theme(legend.position = "left")+
    
    ggtitle(paste0("Probe sets in array by category (", sum(data$count)," probes)")) + guides(fill=guide_legend(title=""))

  
  plot_grid(p1, p2)  
  
}
})

