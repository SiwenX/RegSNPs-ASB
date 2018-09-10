arg <- commandArgs(T)
library(stringr)
setwd("/path/to/file")

snp<-read.table("SNP_in_TFBS.bed",header=FALSE)
tfbs<-read.table("TFBS_with_SNP.bed",header=FALSE)
sam<-read.table("fragments.bed",header=FALSE)

snp[arg[3],]->snp; tfbs[arg[3],]->tfbs
as.numeric(arg[2])->sam[,6]
colnames(snp)<-c("chr","start","end","ref","alt");
colnames(tfbs)<-c("ID1","ID2","chr","start","end");
colnames(sam)<-c("chr", "start", "end", "seq","cigar", "snp");

##N or n
myfunction1 <- function(x){
  ifelse(x[2] < tfbs$start && x[3] > tfbs$end, "N"->n, "n"->n)
}
apply(sam, 1, myfunction1)->aa

##o or x
###adjust the snp position by cigar
adjust_cigar<-function(start,pos,cigar){
  len<-as.numeric(pos)-as.numeric(start)+1;
  
  num<-as.numeric(str_split(cigar,"D|I|M")[[1]]);
  num<-num[!is.na(num)]
  label<-(str_split(cigar,"\\d")[[1]] );
  label<-label[nchar(label)>0]
  
  cur_len<-0;insert_len<-0;del_len<-0;
  for(i in 1:length(num)){
    
    if(label[i]=="I"){
      insert_len<-insert_len+num[i];
    }
    
    if(label[i]=="D"){
      del_len<-del_len+num[i];
    }
    
    if(label[i]!="I"){
      cur_len<-cur_len+num[i]
    }
    
    if(cur_len>len){
      break;
    }
    
    
  }
  
  return(as.numeric(pos)+as.numeric(insert_len)-as.numeric(del_len))
  
}
myfunction2 <- function(x){
  loci<-adjust_cigar(x[2],x[6],x[5])
  strsplit(as.character(x[4]), split = "", fixed = TRUE)[[1]]->seq
  seq[as.numeric(loci)-as.numeric(x[2])+1]->SNP 
}
apply(sam, 1, myfunction2)->SNP
SNP[which(SNP==str_split(snp$ref,"")[[1]][1])]<-"o"; SNP[which(SNP==str_split(snp$alt,"")[[1]][1])]<-"x"

##No, Nx, no, nx
count<-matrix(0,1,length(SNP))
for(i in 1:length(SNP)){
  paste(aa[i],SNP[i],sep = "")->count[,i]
}
length(which(count=="No"))->No
length(which(count=="Nx"))->Nx
length(which(count=="no"))->no
length(which(count=="nx"))->nx
x = c(No, no, Nx, nx)

as.character(tfbs[[2]])->ID ;as.character(tfbs[[3]])->Pos; as.numeric(tfbs[4])->Sta; as.numeric(tfbs[5])->End; as.numeric(snp[6])->Snp
c(ID,Pos,Sta,End,Snp,No,no,Nx,nx)->zz
cat(zz, file = "ASO.txt", append = TRUE)
cat("\n", file = "ASO.txt", append = TRUE)


