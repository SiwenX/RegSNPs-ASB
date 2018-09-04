setwd("/Users/xsw/bio/AS-TF_occupancy/data/human_MCF7")
hotplot<-read.table("MCF7_ana_input",header = T)  #MSC(2418)   #MCF7(923)
#hotplot<-read.table("/Users/xsw/bio/AS-TF_occupancy/result/MCF7_FDR0.05",header = T)  #(1029)

paste(hotplot$chr,hotplot$snp, sep = ":")->bindingsite
bindingsite->hotplot[,10]
hotplot$deta_PSSM * hotplot$beta -> hotplot[,11]
colnames(hotplot) <- c('TF','chr','sta','end','snp','deta_PSSM','beta','ref','alt','BS','*')
#Fiter by TF species
db<-read.csv("/Users/xsw/bio/database/JASPAR2018_CORE_vertebrates_non-redundant_pfms_transfac.txt",header=FALSE,stringsAsFactors = FALSE)
as.matrix.data.frame(db) -> db
#
results<-c()
for(i in 1:nrow(hotplot)){
  if(length(grep(pattern = hotplot[i,1], db)) > 0){
    hotplot[i,]->temp
  } else {
      next
    }
  rbind(results, temp)->results
}
#还剩MSC(1354),MCF7(524)
#Fiter by PSSM & BETA
results[which(results$`*`>0),]->results   #MSC(886);MCF7(300)
results->hotplot2
hotplot2[which(hotplot2$deta_PSSM > 5 | hotplot2$deta_PSSM < -5),]->hotplot3 #MSC(406),MCF7(122)
#
#cytoinput
hotplot3[,c(1,6,5)]->cytoinput
Names<-read.table("factorNames.txt",fill = T,col.names = paste("V", 1:9, sep = ""))
results<-matrix(0,nrow = nrow(cytoinput),ncol = 1)
for(i in 1:nrow(cytoinput)){
  which(Names[,1]==as.character(cytoinput[i,1]))->a
  if(length(a)==0){'NA'->results[i,]}
  else{as.character(Names[a,2])->results[i,]}
}
results->cytoinput[,1]

write.table(cytoinput,"cytoinput.txt",quote=FALSE,row.names=FALSE,col.names=T)

#TF frequency
unique(hotplot3$TF)->TF   #137
number<-c()
for(i in 1:length(TF)){
  length(which(hotplot3$TF == TF[i])) -> number[i]
}
sort(number,decreasing = T)
TF[order(number,decreasing = T)]->TF_order
TF_order[1:20]->TF_order_top20
TF_order->TF_order_top20

#MATLIGN INPUT
PFM_ALL<-c()
for(i in 1:length(TF_order_top20)){
grep(pattern = TF_order_top20[i], db)[2]+2->start
db[start:nrow(db),]->db2
grep(pattern = 'XX', db2)[1]-1->end
db2[1:end]->motif
PFM<-c()
for(j in 1:length(motif)){
  strsplit(motif[j],'\t')[[1]]->temp
  as.numeric(temp[2:length(temp)])->temp2
  rbind(PFM,temp2)->PFM
}
which(Names[,1]==as.character(TF_order_top20[i]))->a
as.character(Names[a,2])->TF
paste('>', TF, sep = '')->tf
rbind(tf,PFM)->a
as.matrix(a)->a
''->a[1,2:4]
rbind(PFM_ALL,a)->PFM_ALL
}
write.table(PFM_ALL,"PFM_ALL_TF_LPS.txt",quote=FALSE,row.names=FALSE,col.names=F,sep = '\t')



#EffSNP
unique(ASO2$V10)
ASO3[,c(2,5,5,8,9)]->SNP
write.table(SNP,"EffSNP_input",quote=FALSE,row.names=FALSE,col.names=F,sep = '\t')
#count.bed
hotplot3[!duplicated(hotplot3$BS),]->hotplot3_unique_SNP
hotplot3_unique_SNP[,c(2,5,5,8,9)]->SNP
write.table(SNP,"/Users/xsw/bio/AS-TF_occupancy/data/human_MCF7/Fig.3/count.bed",quote=FALSE,row.names=FALSE,col.names=F,sep = '\t')
hotplot3_unique_SNP[,c(2,3,4,7)]->SNP2
write.table(SNP2,"/Users/xsw/bio/AS-TF_occupancy/fig/Fig.4/heatmap_MCF7",quote=FALSE,row.names=FALSE,col.names=F,sep = '\t')

#
ASO<-read.csv("ASO_2",header = F,sep = '\t')  #(2905,732)
length(unique(ASO$V4))
ASO[!duplicated(ASO$V4),]->ASO_unique_SNP
write.table(ASO_unique_SNP,"ASO_unique_SNP",quote=FALSE,row.names=FALSE,col.names=F,sep = '\t')
