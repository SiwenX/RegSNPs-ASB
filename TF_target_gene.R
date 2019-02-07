#Extract MCF7 data from all 4D-Genome database
setwd("/Users/xsw/bio/AS-TF_occupancy/data/4Dgenome")
IN<-read.csv("4DGenome_HomoSapiens_hg19.txt", fill = T, sep = "\t")
IN[which(IN$Cell.Tissue == 'MCF7'),]->IN_MCF7
write.table(IN_MCF7,"Interaction_MCF7",quote=FALSE,row.names=FALSE,col.names= T, sep = "\t")

#整理我们得到的RefSNPs
MCF7<-read.table("/Users/xsw/bio/AS-TF_occupancy/result/MCF7_ASB", header = T)
cbind(MCF7,MCF7[,2])->MCF7_bed
write.table(MCF7_bed,"MCF7.bed",quote=FALSE,row.names=FALSE,col.names= F, sep = '\t')
#
write.table(IN_MCF7,"IN_MCF7_AB.bed",quote=FALSE,row.names=FALSE,col.names= F, sep = '\t')
IN_MCF7[,4:6] -> IN_MCF7[,1:3]
write.table(IN_MCF7,"IN_MCF7_BA.bed",quote=FALSE,row.names=FALSE,col.names= F, sep = '\t')
#$ bedtools intersect -wa -wb -a IN_MCF7_BA.bed -b MCF7.bed > ASB_with_col7_interaction.bed
#$ bedtools intersect -wa -wb -a IN_MCF7_AB.bed -b MCF7.bed > ASB_with_col8_interaction.bed
#
ASB_with_col7_interaction<-read.table("ASB_with_col7_interaction.bed", header = F)
ASB_with_col8_interaction<-read.table("ASB_with_col8_interaction.bed", header = F)
ASB_with_col7_interaction[,c(16,17,7)] -> ASB_INTERACT_GENE_1
ASB_with_col8_interaction[,c(16,17,8)] -> ASB_INTERACT_GENE_2
as.data.frame(ASB_INTERACT_GENE_1) -> ASB_INTERACT_GENE_1
as.data.frame(ASB_INTERACT_GENE_2) -> ASB_INTERACT_GENE_2
colnames(ASB_INTERACT_GENE_1) <- c('chr','snp','tarGene')
colnames(ASB_INTERACT_GENE_2) <- c('chr','snp','tarGene')
rbind(ASB_INTERACT_GENE_1, ASB_INTERACT_GENE_2) -> ASB_INTERACT_GENE
ASB_INTERACT_GENE[which(!duplicated(ASB_INTERACT_GENE$tarGene) == T),] -> ASB_INTERACT_GENE
ASB_INTERACT_GENE[-which(is.na(ASB_INTERACT_GENE$tarGene) == T),] -> ASB_INTERACT_GENE
####### gene_region
write.table(ASB_INTERACT_GENE,"ASB_INTERACT_GENE.csv",quote=FALSE,row.names=FALSE,col.names= F, sep = '\t')

ASB_INTERACT_GENE[1,3]
ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")
gene <- read.table("/Users/xsw/bio/AS-TF_occupancy/data/human_MCF7/ASE/gene_ID1", header = F)
aa<-c()
for(i in 1:nrow(gene)){
  getBM(c("hgnc_symbol","chromosome_name","start_position","end_position","ensembl_gene_id"), "ensembl_gene_id", gene$V1[i], ensembl) -> a
  rbind(aa, a) -> aa
}
write.table(aa,"/Users/xsw/bio/AS-TF_occupancy/data/human_MCF7/ASE/gene_ID2_region",quote=FALSE,row.names=FALSE,col.names= F, sep = '\t')
#
gene <- read.table("/Users/xsw/bio/AS-TF_occupancy/data/human_MCF7/ASE/gene_ID2", header = F)
aa<-c()
for(i in 1:nrow(gene)){
  getBM(c("hgnc_symbol","chromosome_name","start_position","end_position","ensembl_transcript_id"), "ensembl_transcript_id", gene$V1[i], ensembl) -> a
  rbind(aa, a) -> aa
}
write.table(aa,"/Users/xsw/bio/AS-TF_occupancy/data/human_MCF7/ASE/transcript_ID2_region",quote=FALSE,row.names=FALSE,col.names= F, sep = '\t')
####### exon_region
gene1 <- read.table("/Users/xsw/bio/AS-TF_occupancy/data/human_MCF7/ASE/gene_ID1", header = F)
gene2 <- read.table("/Users/xsw/bio/AS-TF_occupancy/data/human_MCF7/ASE/gene_ID2", header = F)
rbind(gene1, gene2) -> gene
unique(gene)->gene
aa<-c()
for(i in 1:nrow(gene)){
  getBM(c("external_gene_name","chromosome_name","exon_chrom_start","exon_chrom_end","ensembl_gene_id"), "ensembl_gene_id", gene$V1[i], ensembl) -> a
  rbind(aa, a) -> aa
}
write.table(aa,"/Users/xsw/bio/AS-TF_occupancy/data/human_MCF7/ASE/gene_ID_exon_region",quote=FALSE,row.names=FALSE,col.names= F, sep = '\t')
####
gene_with_snp <- read.table("/Users/xsw/bio/AS-TF_occupancy/data/RNA-seq/target_gene/ASE.txt", header = F)
paste(gene_with_snp$V2, gene_with_snp$V3)->gene_with_snp[,8]
length(unique(gene_with_snp[,8]))
gene_with_snp[!duplicated(gene_with_snp[,8]),]->gene_with_snp
write.table(gene_with_snp,"/Users/xsw/bio/AS-TF_occupancy/data/human_MCF7/ASE/gene_with_snp",quote=FALSE,row.names=FALSE,col.names= F, sep = '\t')
##### promotor_region
gene <- read.table("/Users/xsw/bio/AS-TF_occupancy/data/human_MCF7/ASE/ASB_in_promotor", header = F)
aa<-c()
for(i in 1:nrow(gene)){
  getBM(c("hgnc_symbol","chromosome_name","start_position","end_position","ensembl_gene_id"), "hgnc_symbol", gene$V2[i], ensembl) -> a
  rbind(aa, a) -> aa
}
##### binomial test
gene <- read.table("/Users/xsw/bio/AS-TF_occupancy/data/human_MCF7/ASE/binomial test/target_gene", header = F)
aa<-c()
for(i in 1:nrow(gene)){
  binom.test(as.numeric(gene[i,]), p = 0.5)$p.value->a
  rbind(aa, a)->aa
}
library(fdrtool)
as.matrix(aa)->a; as.numeric(a)->a
fdrtool(a, statistic = "pvalue", plot = FALSE)->fdr
fdr$lfdr->fdr
write.table(fdr,"/Users/xsw/bio/AS-TF_occupancy/data/human_MCF7/ASE/binomial test/fdr.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

write.table(aa,"/Users/xsw/bio/AS-TF_occupancy/data/human_MCF7/ASE/binomial test/target_gene_binomialtest",quote=FALSE,row.names=FALSE,col.names= F, sep = '\t')
# 
gene <- read.table("/Users/xsw/bio/AS-TF_occupancy/data/human_MCF7/ASE/binomial test/ASB_impact_ASE", header = F)
dim(unique(gene))
