##################################################
###########Fiter by TF species
###########Fiter by PSSM & BETA
##################################################

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
