# RegSNPs-ASB
RegSNPs-ASB is a pipeline for extracting regulatory SNPs from ATAC-seq data. RegSNPs-ASB first call TFBS and heterozygote SNP using published tools. In each TFBS with heterozygote SNP, RegSNPs-ASB will take a generalized linear model to identify allele-specific TF binding events. These functional candidate variants can help us understand the molecular mechanism of complex diseases and provide an essential foundation for experimental follow-up analysis.
## Overall flow chart
![Image of flow chart](https://github.com/SiwenX/RegSNP-ASB/blob/master/Figures/Fig2.png)
## Installation
`git clone https://github.com/SiwenX/RegSNP-ASB.git`
## Requirements
  - Linux working environment 
  - R packages
      - GenomicRanges
      - BSgenome
      - MASS
      - biomaRt
      - fdrtools
  - Linux tools
      - samtools
      - plink2
      - vcffilter
## Usage
  - Step 1. Call heterozygote variants
    - `$samtools merge input.bam files` # if you have multiple treatments, you should merge them together to increase the read coverage
    - `$samtools mpileup -uf reference.fa input.bam | bcftools view -Nvcg - > SNP.vcf`
    - `$grep "0/1:" SNP.vcf > hete_SNP.vcf`
    - `$vcffilter -f "DP > 10 & MQ > 20" hete_SNP.vcf > hete_SNP_filtered.vcf` # filter SNP by depth and quality
  - Step 2. Call TFBS
    - prepare sequence file
      ```rscript
      # load required R packages
      library("GenomicRanges")
      library('BSgenome')
      library("BSgenome.Hsapiens.UCSC.hg19")
      
      # load peak file
      peak <- read.table("ATAC.narrowPeak", header = T, stringsAsFactors = FALSE)
      
      # get sequence within peaks
      Range <- with(peak, GRanges(seqnames = chr, ranges = IRanges(start=start, end = end), strand = "*"))
      as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, peak$seqnames, start = peak$start, end = peak$end, strand = "+"))->peak_seq
      for(i in 1 : nrow(peak)){
      cat(">", as.character(peak[i, 1]), file = "fasta.txt", sep = "", append = TRUE)
      cat("-", file = "fasta.txt", append = TRUE)
      cat(as.character(peak[i, 2]),"-",as.character(peak[i, 3]),file = "fasta.txt", sep = "", append = TRUE)
      cat("\n", file = "fasta.txt", append = TRUE)
      cat(peak_seq[i], file = "fasta.txt", sep = "", append = TRUE)
      cat("\n", file = "fasta.txt", append = TRUE)
       }
      ``` 
    - `$fimo <motif file> <fasta.txt> --o <output dir>` # motif file can be downloaded from http://jaspar.genereg.net/downloads/ 
  - Step 3. Call potential allele-specific TFBS
    ```rscript
    # load required R packages
    library("GenomicRanges") 
    
    # load input files
    peak1 <- read.table("TFBS.txt", header = FALSE, stringsAsFactors = FALSE)
    peak2 <- read.table("Hete_SNP.bed", header = FALSE, stringsAsFactors = FALSE)
    colnames(peak1) <- c("ID", "chr", "start", "end");
    colnames(peak2) <- c("chr", "start", "end", "REF", "ALT")
    
    # get potential AS-TFBS
    Range1 <- with(peak1, GRanges(seqnames = chr, ranges = IRanges(start = start,end = end), strand = "*"));
    Range2 <- with(peak2, GRanges(seqnames = chr, ranges = IRanges(start = start,end = end), strand = "*"));
    findOverlaps(Range1, Range2) -> a
    peak1[a@from, ] -> TFBS_with_SNP.bed
    peak2[a@to, ] -> SNP_in_TFBS.bed
    write.table(TFBS_with_SNP.bed, "TFBS_with_SNP.bed", quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(SNP_in_TFBS.bed, "SNP_in_TFBS.bed", quote = FALSE, row.names = FALSE, col.names = FALSE)
    ```
  - Step 4. Using GLM to identify AS-TFBS
    ```
    $Detect_ASB.sh -i SNP_in_TFBS.bed
    ```
  - Step 5. Filtering AS-TFBS
    ```
    $Filter_ASB.sh -i ASB.txt
    ```
