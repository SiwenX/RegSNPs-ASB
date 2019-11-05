arg <- commandArgs(T)
setwd('')
arg[1] -> gene_id
gene <- read.table("gene.bed", header=F, stringsAsFactors = FALSE)
gene[arg[2], ] -> gene
tryCatch({aa <- read.table("gene4.vcf", header=F, stringsAsFactors = FALSE)}, error = function(e) 
{  cat("No SNPs")   })
for(i in 1:nrow(aa)){
aa[i, ] -> aaa
aa[i, 8] -> a
strsplit(a, split = ";")[[1]][5] -> b
strsplit(b, split = "=")[[1]][2] -> c
as.numeric(strsplit(c, split = ",")[[1]]) -> d
d[1] + d[2] -> ref_allele
d[3] + d[4] -> alt_allele
if(is.na(ref_allele)){
  strsplit(a, split = ";")[[1]][4] -> b
  strsplit(b, split = "=")[[1]][2] -> c
  as.numeric(strsplit(c, split = ",")[[1]]) -> d
  d[1]+d[2] -> ref_allele
  d[3]+d[4] -> alt_allele
}
c(gene_id, aaa$V1, aaa$V2, aaa$V4, aaa$V5, ref_allele, alt_allele) -> zz
if(aaa$V2 >= gene[3] && aaa$V2 <= gene[4]){cat(zz, file = "ASE.txt", append = TRUE)
cat("\n", file = "ASE.txt", append = TRUE)
  }
}
