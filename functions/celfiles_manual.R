feno<-read.csv("./amanda/Amanda_miRNA_serum_cel files/amanda_feno.csv",sep = ";",dec=",")
cel_files_in<-data.frame(datapath=list.celfiles("./amanda/Amanda_miRNA_serum_cel files/",full.name=T),
                         name=list.celfiles("./amanda/Amanda_miRNA_serum_cel files/",full.name=F))

cel_files_feno<-    data.frame(cbind(cel_files_in,feno[match(cel_files_in$name,feno[,1]),]))
cel_files_feno_all<-cel_files_feno

pheno<-cel_files_feno_all

# rownames(feno)<-feno[,1]
# pheno_ano<-AnnotatedDataFrame(feno)
rownames(pheno)<-pheno[,2]

pheno_ano<-AnnotatedDataFrame(pheno)
# 
# cel_files_feno<-list.celfiles("./amanda/Amanda_miRNA_serum_cel files/",full.name=T)
# data_raw<-read.celfiles(cel_files_feno,phenoData = pheno_ano)
data_raw<-read.celfiles(cel_files_feno_all$datapath, phenoData = pheno_ano)

data_norm<-oligo::rma(data_raw)

anno<-data_norm@annotation

# files_grep<-list.files("/home/josep/Documents/01_IDISBA_2020/22_CELL_FILES/anotacions/",full.names=T)
files_grep<-list.files("./anotacions/",full.names=T)
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

if(data_norm@annotation=="pd.mirna.4.0"){
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
dev.off()
feat<-annotate_file[match(rownames(data_norm@featureData@data),annotate_file$Probe.Set.Name),]

feat<-AnnotatedDataFrame(feat)
rownames(feat@data)<-
  rownames(data_norm@featureData)
data_norm@featureData<-feat
dup.ids <- feat@data$Accession[duplicated(feat@data$Accession)] %>% 
  unique %>%
  sort



exprs_data<-oligo::exprs(data_norm)
dim(exprs_data)

columna<-"visita"
phenotype_names <- pData(data_norm)[,columna]

annotation_for_heatmap <- data.frame(Phenotype = phenotype_names)
rownames(annotation_for_heatmap) <- colnames(data_norm_SIG_df)
# annotation_for_heatmap<-as.factor(annotation_for_heatmap[,1])

# row.names(annotation_for_heatmap) <- make.names(pData(data_norm)$pacient,unique = T)
# rownames(annotation_for_heatmap) <- make.names(colnames(data_norm_SIG_df),unique = T)
# rownames(annotation_for_heatmap) <- colnames(data_norm_SIG_df)


gens_SIG<-data_norm@featureData@data$Transcript.ID.Array.Design.[1:25]

data_norm_SIG<-data_norm[data_norm@featureData@data$Transcript.ID.Array.Design.%in%gens_SIG,]

data_norm_SIG_df<-exprs(data_norm_SIG)

p1<-pheatmap::pheatmap(data_norm_SIG_df,
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
