#R Script to annotate RBP Binding Sites given an RBP binding consensus pattern

#Read in the the raw data and filter as needed

mbnl12_ko <- fread("rmats_mbnl12_cdko_all.txt")
splicing_df <- mbnl12_ko[df$AStype=="SE", c(1,2,8,11,13)]

#Install and load mouse genomes
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(BSgenome.Mmusculus.UCSC.mm10)
library(org.Mm.eg.db)

##note: to annotate human genome, follow the same steps but install and load the human ref genome for annotation:
#library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#txdb = TxDb.Hsapiens.UCSC.hg38.knownGene

#describe the RBP binding pattern
mbnl_pattern <- "((C|T)GC(C|T))(([(A|T|C|G)]){0,5})?((C|T)GC(C|T))|((C|T)GC(C|T)GC(C|T))"


#note: annotating mouse genomes; change or rename variable to annotate human genome
txdb = TxDb.Mmusculus.UCSC.mm10.knownGene

#annotate the data frame
  splicing_df$sp_coord <- strsplit(splicing_df$AS_ID, split = ";|:")
  splicing_df$chr <- sapply((splicing_df$sp_coord),"[[",1)
  splicing_df$start <- sapply(strsplit(sapply((splicing_df$sp_coord),"[[",2), "-"),"[[",1)
  splicing_df$end <- sapply(strsplit(sapply((splicing_df$sp_coord),"[[",2), "-"),"[[",2)
  splicing_df$strand <- sapply((splicing_df$sp_coord),"[[",5)
  
  #make the start and end coordinates numeric so we can get the sequence strings from the mouse genome
  splicing_df$start <- as.numeric(splicing_df$start)
  splicing_df$end <- as.numeric(splicing_df$end)
  
  #remove alt chromosomes
  splicing_df <- splicing_df[which(splicing_df$chr %in% chromosomes),]

  #get 1000bp of sequence of upstream and downstream of exon junction
  splicing_df$upstream <- as.character(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38, splicing_df$chr, splicing_df$start-1000, splicing_df$start))
  splicing_df$downstream <- as.character(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38, splicing_df$chr, splicing_df$end, splicing_df$end+1000))
  
  #for genes on the reverse (-) strand, we need to flip the upstream and downstream sequences
  library(spgs)

  splicing_df$upstream[which(splicing_df$strand=="-")] <- 
    as.character(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38, 
                                    splicing_df$chr[which(splicing_df$strand=="-")],
                                    splicing_df$end[which(splicing_df$strand=="-")],
                                    splicing_df$end[which(splicing_df$strand=="-")]+1000))
  
  
  splicing_df$downstream[which(splicing_df$strand=="-")] <- 
    as.character(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38,
                                    splicing_df$chr[which(splicing_df$strand=="-")],
                                    splicing_df$start[which(splicing_df$strand=="-")]-1000,
                                    splicing_df$start[which(splicing_df$strand=="-")]))
  
  
  splicing_df$upstream[which(splicing_df$strand=="-")] <- spgs::reverseComplement(splicing_df$upstream[which(splicing_df$strand=="-")])
  splicing_df$downstream[which(splicing_df$strand=="-")] <- spgs::reverseComplement(splicing_df$downstream[which(splicing_df$strand=="-")])
  
  
  #determine the closest distance to the RBP binding pattern

  #annotate closest upstream RBP binding
  splicing_df$mbnl.up_dist_pattern <- 1000 - stri_locate_last_regex(splicing_df$upstream, mbnl_pattern)[,1]
  #annotate closest downstream RBP binding
  splicing_df$mbnl.down_dist_pattern <- stri_locate_first_regex(splicing_df$downstream, mbnl_pattern)[,2]
  
  #for instances where RBP pattern is not found, numerically categorize out of bounds
  splicing_df$mbnl.up_dist_pattern[is.na(splicing_df$mbnl.up_dist_pattern)] <- 1001
  splicing_df$mbnl.down_dist_pattern[is.na(splicing_df$mbnl.down_dist_pattern)] <- 1001
  
  #determine the closest RBP binding pattern, regardless of whether upstream or downstream
  splicing_df$mbnl.mindist_pattern <- apply(splicing_df[,c(13,14)], 1, FUN = min, na.rm = TRUE)
  
  #if both patterns are the same, R throws an error - simply label those cases as NA
  splicing_df$mbnl.mindist_pattern[is.infinite(splicing_df$mbnl.mindist_pattern)] <- NA 

#once nearest RBP binding site is annotated, statistics and analyses can be performed based on any distance parameters
