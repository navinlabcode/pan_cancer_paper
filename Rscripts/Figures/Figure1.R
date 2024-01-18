#Read summary data
summary<-readRDS("Data/summary.rds")

#assign colors
cancer_type_colors<-c('hotpink','slategray3', 'lawngreen', 'gold', 'thistle', 'mediumorchid2', 'dodgerblue2')
names(cancer_type_colors)<-unique(summary$cancer_type)

library(ggplot2)

#Plot CNA burden
ggplot(summary, aes(y=CNA_burden, x=cancer_type))+geom_boxplot(aes(x=reorder(cancer_type, CNA_burden, FUN=median)),
                                                              show.legend = FALSE, outlier.shape = NA)+theme_bw()+
  xlab("Cancer type")+ylab("CNA burden")+ theme(axis.line = element_line(linetype = 1), axis.text.x=element_text(size=15), 
                                                axis.title.y =element_text(size=15), axis.title.x  =element_text(size=15),
                                                axis.text.y=element_text(size=15), panel.grid.minor = element_blank(), 
                                                panel.border = element_blank(),panel.grid.major = element_blank(),
                                                panel.background = element_blank() , legend.position = 'none')+ 
  geom_point(position = position_jitter(seed=17, width = 0.2), aes(color=cancer_type), size=2, alpha=1)+
  scale_color_manual(values = cancer_type_colors)

#Plot mean ploidy
ggplot(summary, aes(y=mean_ploidy, x=cancer_type))+geom_boxplot(aes(x=reorder(cancer_type, mean_ploidy, FUN=median)),
                                                                show.legend = FALSE, outlier.shape = NA)+theme_bw()+
  xlab("Cancer type")+ylab("DNA Ploidy (N)")+ theme(axis.line = element_line(linetype = 1), axis.text.x=element_text(size=15), 
                                                    axis.title.y =element_text(size=15), axis.title.x  =element_text(size=15),
                                                    axis.text.y=element_text(size=15), panel.grid.minor = element_blank(), 
                                                    panel.border = element_blank(),panel.grid.major = element_blank(),
                                                    panel.background = element_blank() , legend.position = 'none')+ 
  geom_point(position = position_jitter(seed=17, width = 0.2), aes(color=cancer_type), size=2, alpha=1)+
  scale_color_manual(values = cancer_type_colors)

#PLot exonic mutation burden
all_snv<-readRDS("Data/All_SNV.rds")
all_snv_exonic<-all_snv[all_snv$Func.refGene=='exonic',]
count_burden<-table(all_snv_exonic$patient)
count_burden<-as.data.frame(count_burden)
count_burden<-cbind(distinct(all_snv_exonic, all_snv_exonic$patient, .keep_all = T)$cancer_type, count_burden)
colnames(count_burden)<-c('cancer_type', 'patient', 'count')
count_burden$log_count<-log(count_burden$count)
ggplot(count_burden, aes(y=log_count, x=cancer_type))+geom_boxplot(aes(x=reorder(cancer_type, count, FUN=median)),
                                                                   show.legend = FALSE, outlier.shape = NA)+theme_bw()+
  xlab("Cancer type")+ylab("Mutation burden (log)")+ theme(axis.line = element_line(linetype = 1), axis.text.x=element_text(size=15),
                                                           axis.title.y =element_text(size=15), axis.title.x  =element_text(size=15),
                                                           axis.text.y=element_text(size=15), panel.grid.minor = element_blank(), 
                                                           panel.border = element_blank(),panel.grid.major = element_blank(),
                                                           panel.background = element_blank() , legend.position = 'none')+ 
  geom_point(position = position_jitter(seed=17, width = 0.2), aes(color=cancer_type), size=2, alpha=1)+
  scale_color_manual(values = cancer_type_colors)+scale_y_continuous(breaks = seq(2, 10, by=1))

#Oncomap
library(maftools)
#Generate a maf file
pan.annovar.maf = annovarToMaf(annovar = 'Data/all_mutations.txt', refBuild = 'hg19',
                               tsbCol = 'patient', MAFobj = TRUE)

#Prepare the data for MutsigCV analysis
pan.annovar.maf.corrected<-prepareMutSig(maf=pan.annovar.maf, fn='Data/pan.annovar')

#Read top mutated genes
top.genes<-pan.annovar.maf@gene.summary$Hugo_Symbol[1:400]

#Read mutsig results
mutsig<-read.table("Data/pan.cancer.mutsig.sig_genes.txt", header=TRUE)
gene.symbols.match<-read.table("Data/pan.annovar.correctedSymbols.tsv", header = T)

#Mutsig significant genes
sig_genes<-mutsig$gene[mutsig$p<=0.001]
#Convert back
for(i in 1:length(sig_genes)){
  for (j in 1:nrow(gene.symbols.match)){
    if(sig_genes[i]==gene.symbols.match$MutSig_Synonym[j]){
      sig_genes[i]<-gene.symbols.match$Hugo_Symbol[j]
    }
  }
}

#Plot mutsig significant genes
#Add cancer type info
Cancer_type<-c()
for(i in 1:nrow(pan.annovar.maf@clinical.data)){
  Cancer_type<-c(Cancer_type, summary$cancer_type[which(summary$patient==pan.annovar.maf@clinical.data$Tumor_Sample_Barcode[i])])
}
pan.annovar.maf@clinical.data<-cbind(pan.annovar.maf@clinical.data, Cancer_type)
pdf("~/pan_cancer_project/Figures/Figure1_summary/oncomap.pdf", height = 8, width = 8)
oncoplot(pan.annovar.maf,showTumorSampleBarcodes = F,clinicalFeatures = 'Cancer_type',sortByAnnotation = T, 
         annotationColor =list(Cancer_type=cancer_type_colors), barcode_mar = 6,gene_mar = 6, 
         genes = sig_genes[sig_genes%in%top.genes])
dev.off()