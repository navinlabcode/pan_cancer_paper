#Read summary data
summary<-readRDS("Data/summary.rds")

#Assign cancer type colors
cancer_type_colors<-c('hotpink','slategray3', 'lawngreen', 'gold', 'thistle', 'mediumorchid2', 'dodgerblue2')
names(cancer_type_colors)<-unique(summary$cancer_type)

#Plot MPD
ggplot(summary, aes(y=MPD, x=cancer_type))+geom_boxplot(aes(x=reorder(cancer_type, MPD, FUN=median)),
                                                                 show.legend = FALSE, outlier.shape = NA)+theme_bw()+
  xlab("Cancer type")+ylab("MPD")+ theme(axis.line = element_line(linetype = 1), axis.text.x=element_text(size=15), 
                                         axis.title.y =element_text(size=15), axis.title.x  =element_text(size=15),
                                         axis.text.y=element_text(size=15), panel.grid.minor = element_blank(), 
                                         panel.border = element_blank(),panel.grid.major = element_blank(),
                                         panel.background = element_blank() , legend.position = 'none')+ 
  geom_point(position = position_jitter(seed=17, width = 0.2, height = 0.1), aes(color=cancer_type), size=2, alpha=1)+
  scale_color_manual(values = cancer_type_colors)

#Correlate MPD with CNA burden
ggplot(summary)+geom_point(aes(x=CNA_burden, y=MPD, color=cancer_type))+
  geom_smooth(aes(x=CNA_burden, y=MPD),method='lm',color='black', formula= y~x)+
  ylab('MPD')+xlab('CNA burden')+scale_color_manual(values=cancer_type_colors)+theme_classic()

#group patients based on ploidy
summary<-cbind(summary, rep("<2.5N", times=nrow(summary)))
colnames(summary)[ncol(summary)]<-'ploidy.group'
summary$ploidy.group<-as.character(summary$ploidy.group)
for(i in 1: nrow(summary)){
  if(summary$mean_ploidy[i]>=2.5){
    summary$ploidy.group[i]<-'>=2.5N'
  }
}
ggplot(summary, aes(y=MPD, x=ploidy.group))+geom_boxplot(show.legend = FALSE, outlier.shape = NA)+theme_bw()+
  xlab("DNA Ploidy")+ylab("MPD")+ theme(axis.line = element_line(linetype = 1), axis.text.x=element_text(size=15), 
                                        axis.title.y =element_text(size=15), axis.title.x  =element_text(size=15),
                                        axis.text.y=element_text(size=15), panel.grid.minor = element_blank(), 
                                        panel.border = element_blank(),panel.grid.major = element_blank(),
                                        panel.background = element_blank() )+ 
  geom_point(position = position_jitter(seed=17, width = 0.2, height = 0.1), aes(color=cancer_type), size=2, alpha=1)+
  scale_color_manual(values = cancer_type_colors)

#Group patients based on TP53 mutation status
all_mutations<-read.table("Data/All_mutations.txt", header = T)
all_tp53<-all_mutations[all_mutations$Gene.refGene=='TP53',]
all_tp53<-all_tp53[all_tp53$ExonicFunc.refGene!='synonymous SNV', ]
tp53.mut.patients<-as.character(all_tp53$patient)
tp53.mut.patients<-tp53.mut.patients[!duplicated(tp53.mut.patients)]
all.patients<-unique(all_mutations$patient)
tp53.wt.patients<-all.patients[!all.patients %in% tp53.mut.patients]
summary_exome<-summary[summary$patient %in% c(tp53.mut.patients, tp53.wt.patients),]
TP53.status<-c()
for (i in 1:nrow(summary_exome)){
  if(summary_exome$patient[i] %in% tp53.mut.patients){
    TP53.status<-c(TP53.status, "MUT")
  }
  else{
    TP53.status<-c(TP53.status, "WT")
  }
}
summary_exome<-cbind(summary_exome, as.data.frame(TP53.status))
summary_exome$TP53.status<-factor(summary_exome$TP53.status, levels = c('WT', 'MUT'))
ggplot(summary_exome, aes(y=MPD, x=TP53.status))+geom_boxplot(show.legend = FALSE, outlier.shape = NA)+theme_bw()+
  xlab("TP53 status")+ylab("MPD")+ theme(axis.line = element_line(linetype = 1), axis.text.x=element_text(size=15),
                                         axis.title.y =element_text(size=15), axis.title.x  =element_text(size=15),
                                         axis.text.y=element_text(size=15), panel.grid.minor = element_blank(), 
                                         panel.border = element_blank(),panel.grid.major = element_blank(),
                                         panel.background = element_blank() )+ 
  geom_point(position = position_jitter(seed=17, width = 0.2, height = 0.1), aes(color=cancer_type), size=2, alpha=1)+
  scale_color_manual(values = cancer_type_colors)

#Plot the correlation between chromosome gains and MPD in clonal WGD samples
df<-readRDS("Data/Pseudobulk_MPD_cWGD_samples.rds")
#In this data frame, if rel_gain is 1 and p_value<0.05, tumors with gains in that bin have higher MPD, if rel_gain is -1 and p_value<0.05, tumors with gains in that bin have lower MPD

bins_in_cna_pipeline<-read.delim("Data/bins_in_cna_pipeline_bands.bed")
bins_in_cna_pipeline<-bins_in_cna_pipeline[1:12167,]
bins_in_cna_pipeline$chr<-gsub("chr", '', bins_in_cna_pipeline$chr)
bins_in_cna_pipeline$chr[bins_in_cna_pipeline$chr=='X']<-23
bins_in_cna_pipeline$chr<-as.numeric(bins_in_cna_pipeline$chr)
bins_in_cna_pipeline<-cbind(bins_in_cna_pipeline, df$abspos)
colnames(bins_in_cna_pipeline)[ncol(bins_in_cna_pipeline)]<-'abspos'
chr_lengths <- rle(bins_in_cna_pipeline$chr)$lengths
chr_ranges_start <- bins_in_cna_pipeline %>%
  dplyr::group_by(chr) %>%
  dplyr::arrange(chr, start) %>%
  dplyr::filter(dplyr::row_number() == 1) %>%
  dplyr::ungroup()
chr_ranges_end <- bins_in_cna_pipeline %>%
  dplyr::group_by(chr) %>%
  dplyr::arrange(chr, start) %>%
  dplyr::filter(dplyr::row_number() == dplyr::n()) %>%
  dplyr::ungroup()

# Creating data frame object for chromosome rectangles shadows
chrom_rects <- data.frame(
  chr = chr_ranges_start$chr,
  xstart = as.numeric(chr_ranges_start$abspos),
  xend = as.numeric(chr_ranges_end$abspos)
)
xbreaks <- rowMeans(chrom_rects %>%
                      dplyr::select(
                        xstart,
                        xend
                      ))

chrom_rects$colors <- rep_len(
  c("white", "gray"),
  nrow(chrom_rects))

# Creating the geom_rect object
ggchr_back <-
  list(
    geom_rect(
      data = chrom_rects,
      aes(
        xmin = xstart,
        xmax = xend,
        ymin = -Inf,
        ymax = Inf,
        fill = colors
      ),
      alpha = .2
    ),
    scale_fill_identity()
  )

sec_breaks <- c(0, 0.5e9, 1e9, 1.5e9, 2e9, 2.5e9, 3e9)
sec_labels <- c(0, 0.5, 1, 1.5, 2, 2.5, 3)

# theme
ggaes <- list(
  scale_x_continuous(
    breaks = xbreaks,
    labels = gsub("23", "X", chrom_rects$chr),
    position = "top",
    expand = c(0, 0),
    sec.axis = sec_axis(
      ~.,
      breaks = sec_breaks,
      labels = sec_labels,
      name = "genome position (Gb)"
    )
  ),
  theme_classic(),
  xlab(""),
  ylab("-log10(p.value)"),
  theme(
    axis.text.x = element_text(
      angle = 0,
      vjust = .5,
      size = 15
    ),
    axis.text.y = element_text(size = 15),
    legend.position = "none",
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 15),
    plot.title = element_text(size = 15),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )
)
color_bin<-  structure(
  c('blue', 'red', 'grey'),
  names = c(-1, 1, 0)  
)

#Plot final results
ggplot(df) + ggchr_back + ggaes +geom_point(
  aes(abspos, -log10(p_value), color=rel_gain)
)+scale_color_manual(values=color_bin)

#Correlate MPD with overall survival
survival<-read.csv('Data/survival_data.csv')
colnames(survival)<-c('cancer_type', 'patient', 'age', 'os_months', 'os_event')
survival$os_event[survival$os_event=='D']<-1
survival$os_event[survival$os_event!=1]<-0
survival$os_event<-as.numeric(survival$os_event)
summary_ordered<-summary[match(survival$patient, summary$patient,),]
ploidy.group<-rep('<2.5N', times=nrow(summary_ordered))
for(i in 1:nrow(summary_ordered)){
  if(summary_ordered$mean_ploidy[i]>2.5){
    ploidy.group[i]<-'>2.5N'
  }
}
mpd_group<-ifelse(summary_ordered$MPD>median(summary_ordered$MPD), 'High', 'Low')
survival<-cbind(survival, dplyr::select(summary_ordered, 'CNA_burden'), ploidy.group, mpd_group)

# Define the cut-off point
cutoff <- 120

# Modify the time and status variables
survival$modified_time <- ifelse(survival$os_months > cutoff, cutoff, survival$os_months)
survival$modified_status <- ifelse(survival$os_months > cutoff, 0, survival$os_event)

library(survminer)
library(survival)

#Adjust cancer type order
survival$cancer_type<-factor(survival$cancer_type, levels = c('Breast', 'Bladder', 'Colon', 'GBM', 'Kidney', 'Lung', 'Ovary'))
fit<-coxph(Surv(modified_time, modified_status) ~ age +  mpd_group  + cancer_type  + ploidy.group + CNA_burden, data = survival)
ggforest(fit)

#Plot gene clonality
pan_cancer_gene_clonality<-readRDS("Data/pan_cancer_gene_clonality.rds")
pan_cancer_genes<-c('TP53', 'KRAS', 'TERT', 'PIK3CA', 'APC', 'ARID1A', 'KMT2D', 'PTEN', 'KMT2C', 'EGFR', 'CDKN2A', 'CDKN2B',
                    'CCND1', 'MYC', 'ERBB2', 'FGF19', 'MDM2', 'FGF4', 'FGF3', 'FGFR1', 'CDK4', 'CDK12', 'MCL1', 
                    'RECQL4', 'RB1', 'CCNE1', 'FAT1', 'NF1', 'BRAF', 'ATM', 'ATRX', 'PTPRT', 'ZFHX3', 'SMAD4')
pan_cancer_genes_category<-c('TSG', 'ONC', 'ONC', 'ONC', 'TSG', 'TSG', 'TSG', 'TSG', 'TSG', 'ONC', 'TSG', 'TSG', 
                             'ONC', 'ONC', 'ONC', 'ONC', 'ONC', 'ONC', 'ONC', 'ONC', 'ONC', 'ONC', 'ONC',
                            'ONC', 'TSG', 'ONC', 'TSG', 'TSG', 'ONC', 'TSG', 'TSG', 'TSG', 'TSG', 'TSG')
# Make a barplot
gene<-c()
clonality<-rep(c('Clonal', 'Subclonal'), length(pan_cancer_genes))
count<-c()
for(i in 1:nrow(pan_cancer_gene_clonality)){
  gene<-c(gene, rep(pan_cancer_genes[i],2))
  count<-c(count, length(which(pan_cancer_gene_clonality[i,]=='Clonal')),
           length(which(pan_cancer_gene_clonality[i,]=='Subclonal')))
}
gene_clonality_count<-data.frame(gene,  clonality, count)
gene_clonality_count$gene <-factor(gene_clonality_count$gene , levels = unique(gene_clonality_count$gene)[order(gene_clonality_count$count[which(gene_clonality_count$clonality=='Subclonal')], decreasing = F)])

pan_cancer_genes_color<-pan_cancer_genes_category
pan_cancer_genes_color[pan_cancer_genes_color=='TSG']<-'blue'
pan_cancer_genes_color[pan_cancer_genes_color=='ONC']<-'red'
names(pan_cancer_genes_color)<-pan_cancer_genes
ggplot(gene_clonality_count, aes(fill=clonality, y=count, x=gene)) + 
  geom_bar(position="stack", stat="identity") + 
  scale_fill_manual(values = c('darkolivegreen1', 'darkorchid1')) + theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1, size = 10, 
                                                                                                     color = pan_cancer_genes_color[levels(gene_clonality_count$gene)]),
                                                                          legend.title = element_blank(), panel.spacing = unit(0.1, "lines"),  panel.grid.minor = element_blank(), 
                                                                          panel.border = element_blank(),panel.grid.major = element_blank(),
                                                                          panel.background = element_blank(),
                                                                          axis.title.y =element_text(size=10))+ xlab("")+ylab("Number of patients")+
  scale_y_continuous(breaks=seq(0, 75, by=5), expand = c(0,0,0,0)) 
