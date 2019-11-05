##################################################
###########Fiter by TF species
###########Fiter by PSSM & BETA
##################################################

setwd("")
hotplot <- read.table("MCF7_ana_input", header = T)  

paste(hotplot$chr, hotplot$snp, sep = ":") -> bindingsite
bindingsite -> hotplot[, 10]
hotplot$deta_PSSM * hotplot$beta -> hotplot[, 11]
colnames(hotplot) <- c('TF', 'chr', 'sta', 'end', 'snp', 'deta_PSSM', 'beta', 'ref', 'alt', 'BS', '*')
#Fiter by TF species
db <- read.csv("JASPAR2018_CORE_vertebrates_non-redundant_pfms_transfac.txt", header = FALSE, stringsAsFactors = FALSE)
as.matrix.data.frame(db) -> db
#
results <- c()
for(i in 1:nrow(hotplot)){
  if(length(grep(pattern = hotplot[i, 1], db)) > 0){
    hotplot[i, ] -> temp
  } else {
      next
    }
  rbind(results, temp) -> results
}
#Fiter by PSSM & BETA
results[which(results$`*` > 0), ] -> results   
results -> hotplot2
hotplot2[which(hotplot2$deta_PSSM > 5 | hotplot2$deta_PSSM < -5), ] -> hotplot3 
#
write.table(hotplot3, "MCF7_filterd", quote = FALSE, row.names = FALSE, col.names = FALSE)
