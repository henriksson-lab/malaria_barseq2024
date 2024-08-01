
library(ggplot2)
library(stringr)

#previous counter in /corgi/websites/malariascreenviewer/countcrispr.R 


# Counts/Total
# Counts/Reference
# Normalized growth rate


# Weight by parasitemia
# Shrinkage by average depth


################################################################################
####################### Alignment based approach for barseq ####################
################################################################################

#17S??
#ACTTCAAAATAGAATTCTACCACCCAGCTTTCTTGTACAAAGTGGTTGATATCTCTATAGTCTCAGTAGTCGG CCCGCTTTCAG CTGACGCACACGAATTACAGATCCGAAGAGAGCTTCAGCAGGAATGCCGAGACCGA

#### R1
#CCTTCAATTTCGATGGGTACCACCCAGCTTTCTTGTACAAAGTGGTTGATATCTCTATAGTCGCAGTAGGCGG TCGCGGCATTT CTGACGCGCACGAATTACAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTC AGTGATT
#CCTTCAATTTCGATGGGTACCACCCAGCTTTCTTGTACAAAGTGGTTGATATCTCTATAGTCGCAGTAGGCGG TAAAGAGAAGC CTGACGCGCACGAATTACAAATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTC AGTGATT
#CCTTCAATTTCGATGGGTACCACCCAGCTTTCTTGTACAAAGTGGTTGATATCTCTATAGTCGCAGTAGGCGG TCGCGGCATTT CTGACGCGCACGAATTACAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTC AGTGATT
#CCTTCAATTTCGATGGGTACCACCCAGCTTTCTTGTACAAAGTGGTTGATATCTCTATAGTCGCAGTAGGCGG ACACCGATGGG CTGACGCGCACGAATTACAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTC AGTGATT
#CCTTCAATTTCGATGGGTACCACCCAGCTTTCTTGTACAAAGTGGTTGATATCTCTATAGTCGCAGTAGGCGG TAAAGAGAAGC CTGACGCGCACGAATTACAAATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTC AGTGATT
#                                                                          CGTTCTTCGGC                                                               lib index?
#CCTTCAATTTCGATGGGTACCACCCAGCTTTCTTGTACAAAGTGGTTGATATCTCTATAGTCGCAGTAGGCGG GCTTTATTCTA CTGACGCGCACGAATTACAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTC GAATCCC
#CCTTCAATTTCGATGGGTACCACCCAGCTTTCTTGTACAAAGTGGTTGATATCTCTATAGTCGCAGTAGGCGG GAACATGCCAT CTGACGCGCACGAATTACAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTC GAATGCC
#CCTTCAATTTCGATGGGTACCACCCAGCTTTCTTGTACAAAGTGGTTGATATCTCTATAGTCGCAGTAGGCGG ACGACATCTAC CTGACGCGCACGAATTACAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTC GAATCCC
#CCTTCAATTTCGATGGGTACCACCCAGCTTTCTTGTACAAAGTGGTTGATATCTCTATAGTCGCAGTAGGCGG GACCATCTACT CTGACGCGCACGAATTACAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTC GAATCAC

#### R2
#GTAATTCGTGCGCGTCAG AGCGTCTGCTT CCGCCTACTGCGACTATAGAGATATCAACCACTTTGTACAAGAAAGCTGGGTGGTACCCATCGAAATTGAAGGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGC
#GTAATTCGTGCGCGTCAG GCCTATCCTAT CCGCCTACTGCGACTATAGAGATATCAACCACTTTGTACAAGAAAGCTGGGTGGTACCCATCGAAATTGAAGGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGC
#GTAATTCGTGCGCGTCAG CGAGTATGCTT CCGCCTACTGCGACTATAGAGATATCAACCACTTTGTACAAGAAAGCTGGGTGGTACCCATCGAAATTGAAGGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGC
#GTAATTCGTGCGCGTCAG GAGATGTTAGG CCGCCTACTGCGACTATAGAGATATCAACCACTTTGTACAAGAAAGCTGGGTGGTACCCATCGAAATTGAAGGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGC
#GTAATTCGTGCGCGTCAG CACGGCGCTCA CCGCCTACTGCGACTATAGAGATATCAACCACTTTGTACAAGAAAGCTGGGTGGTACCCATCGAAATTGAAGGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGC
#                   CGTTCTTCGGC

### could attempt using fastadapt to remove surrounding!


pwalign::pairwiseAlignment(
  "CCTTCAATTTCGATGGGTACCACCCAGCTTTCTTGTACAAAGTGGTTGATATCTCTATAGTCGCAGTAGGCGG NNNNNNNNNNN CTGACGCGCACGAATTACAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTC AGTGATT", 
  "CCTTCAATTTCGATGGGTACCACCCAGCTTTCTTGTACAAAGTGGTTGATATCTCTATAGTCGCAGTAGGCGG TCGCGGCATTT CTGACGCGCACGAATTACAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTC AGTGATT")

outp <- pwalign::pairwiseAlignment(
  "CCTTCAATTTCGATGGGTACCACCCAGCTTTCTTGTACAAAGTGGTTGATATCTCTATAGTCGCAGTAGGCGGNNNNNNNNNNNCTGACGCGCACGAATTACAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCAGTGATT", 
  "ACTTCAAAATAGAATTCTACCACCCAGCTTTCTTGTACAAAGTGGTTGATATCTCTATAGAATCTCAGTAGTCGGCCCGCTTTCAGCTGACGCACACGAATTACAGATCCGAAGAGAGCTTCAGCAGGAATGCCGAGACCGA")
pwalign::alignedPattern(outp)
pwalign::alignedSubject(outp)

ssubject <- as.character(pwalign::alignedSubject(outp))
spattern <- as.character(pwalign::alignedPattern(outp))
toget <- stringr::str_locate_all(spattern, "N")[[1]]#[,1]
str_sub(ssubject,toget[1,1], toget[nrow(toget),1])

if(FALSE){
  #Generate indices for bwa
  frank_bc <- read.csv("/corgi/otherdataset/ellenbushell/barcode_to_gene_210920_FRANK.csv")
  frank_bc$sgrna <- str_split_fixed(frank_bc$gene,"\\|",2)[,1]
  frank_bc$seq <- str_to_upper(frank_bc$barcode)
  frank_bc$seq_rc <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(frank_bc$seq))) 

  tostore <- Biostrings::DNAStringSet(paste0(
    "CCTTCAATTTCGATGGGTACCACCCAGCTTTCTTGTACAAAGTGGTTGATATCTCTATAGTCGCAGTAGGCGG",
    frank_bc$seq_rc,
    "CTGACGCGCACGAATTACAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTC"))
  names(tostore) <- frank_bc$seq
  Biostrings::writeXStringSet(tostore, "/corgi/otherdataset/ellenbushell/bwa_barseq/bc.fa")
}

################################################################################
####################### Generate count using alignment -- barseq ###############   
################################################################################


listpools <- c(
  "EB_deepseq_barseqpool3",
  "EB_priming_barseqpool3"
)

####### Alignment part
for(curpool in listpools){
  print(curpool)
  allpooldir <- "/corgi/otherdataset/ellenbushell/barseq_pools"
  pooldir <- file.path(allpooldir, curpool)
  fastqdir <- file.path(pooldir, "fastq")

  list_bclist <- list()
  for(onef in list.files(fastqdir)){
    if(str_ends(onef,"001.fastq.gz")){
      print(onef)

      #Align
      the_cmd <- paste0("bwa mem /corgi/otherdataset/ellenbushell/bwa_barseq/bc.fa ",file.path(fastqdir,onef)," -t 10 -o ",file.path(fastqdir,onef),".sam")
      print(the_cmd)
      system(the_cmd)
      
      #Sort, make bam
      the_cmd <- paste0("samtools sort ",file.path(fastqdir,onef),".sam"," -o ",file.path(fastqdir,onef),".bam")
      print(the_cmd)
      system(the_cmd)

      #Index the bam
      the_cmd <- paste0("samtools index ",file.path(fastqdir,onef),".bam")
      print(the_cmd)
      system(the_cmd)
    }
  }
}


curpool <- "EB_deepseq_barseqpool3"

####### Counting part
for(curpool in listpools){
  
  print(curpool)
  
  allpooldir <- "/corgi/otherdataset/ellenbushell/barseq_pools"
  pooldir <- file.path(allpooldir, curpool)
  fastqdir <- file.path(pooldir, "fastq")
  bcfile <- file.path(pooldir, "used_bc.csv")
  countfile <- file.path(pooldir,"counts_bwa.RDS")  #note, different!
  countfile_csv <- file.path(pooldir,"counts_bwa.csv")
  
  frank_bc <- read.csv("/corgi/otherdataset/ellenbushell/barcode_to_gene_210920_FRANK.csv")
  frank_bc$sgrna <- str_split_fixed(frank_bc$gene,"\\|",2)[,1]
  frank_bc$seq <- str_to_upper(frank_bc$barcode)
  
  list_bclist <- list()
  for(onef in list.files(fastqdir)){
    if(str_ends(onef,"R1_001.fastq.gz.bam")){
      print(onef)
      onef2 <- str_replace_all(onef, "R1","R2")  #Check also corresponding second read if present
      count_r1 <- read.table(pipe(paste0("samtools idxstats ",file.path(fastqdir,onef))))

#      bclist <- sqldf::sqldf("select count(bc) as cnt, bc from bclist group by bc order by cnt desc")
#      bclist$file <- onef
      
      onedf <- data.frame(bc=count_r1$V1, cnt=count_r1$V3)
      onedf$file <- onef
      list_bclist[[onef]] <- onedf
    }
  }
  counts <- do.call(rbind,list_bclist)
  counts <- reshape2::acast(counts, bc~file, value.var = "cnt", fill = 0)
  colSums(counts)
  rowSums(counts) #no *!
  
  #Which BC to expect here? Only keep these
  if(FALSE){
    usedbc <- read.csv(bcfile,sep="\t")
    usedbc <- frank_bc[frank_bc$sgrna %in% usedbc$gene,]
  } else {
    usedbc <- frank_bc
  }
  

  #Prepare column names (samples), and order by depth
  counts <- counts[order(rowSums(counts), decreasing = TRUE),]
  counts <- counts[rownames(counts) %in% usedbc$seq,]
  colnames(counts) <- str_split_fixed(colnames(counts),"_S",2)[,1]
  
  #Set rownames to gene name
  rownames(usedbc) <- usedbc$seq
  rownames(counts) <- usedbc[rownames(counts),]$sgrna
  
  if(TRUE){
    saveRDS(counts, countfile)
    write.csv(counts, countfile_csv)  #temp
  } 
}


################################################################################
####################### Test bwa count matrices ################################
################################################################################


############# Compare BWA counts vs simpler method
if(FALSE){
  
  dat_simple <- read.csv("/corgi/otherdataset/ellenbushell/barseq_pools/EB_priming_barseqpool3/counts.csv")
  dat_simple <- dat_simple[!str_starts(dat_simple$X,"_"),,drop=FALSE]
  dat_simple <- dat_simple[!duplicated(dat_simple$X),,drop=FALSE] #horrible
  rownames(dat_simple) <- dat_simple$X
  dat_simple <- dat_simple[,-1,drop=FALSE]
  
  dat_simple <- dat_simple[,str_detect(colnames(dat_simple),"R1")]
  colnames(dat_simple) <- str_remove(colnames(dat_simple),".R1")
  dat_simple <- dat_simple[,order(colnames(dat_simple)),drop=FALSE]
  
  dat_bwa <- readRDS("/corgi/otherdataset/ellenbushell/barseq_pools/EB_priming_barseqpool3/counts_bwa.RDS")
  dat_bwa <- dat_bwa[,order(colnames(dat_bwa)),drop=FALSE]
  
  
  compare_genes <- sort(intersect(rownames(dat_bwa), rownames(dat_simple)))
  dat_bwa <- dat_bwa[compare_genes,,drop=FALSE]
  dat_simple <- dat_simple[compare_genes,,drop=FALSE]
  
  plot(dat_bwa[,1], dat_simple[,1])
  plot(log10(1+dat_bwa[,1]), log10(1+dat_simple[,1]))
  
}


############# Compare BWA counts, old lib vs deep seq
if(FALSE){
  
  dat_shallow <- readRDS("/corgi/otherdataset/ellenbushell/barseq_pools/EB_priming_barseqpool3/counts_bwa.RDS")
  dat_shallow <- dat_shallow[,order(colnames(dat_shallow)),drop=FALSE]
  
  dat_deep <- readRDS("/corgi/otherdataset/ellenbushell/barseq_pools/EB_deepseq_barseqpool3/counts_bwa.RDS")
  dat_deep <- dat_deep[,order(colnames(dat_deep)),drop=FALSE]
  
  
  compare_genes <- sort(intersect(rownames(dat_shallow), rownames(dat_deep)))
  dat_shallow <- dat_shallow[compare_genes,,drop=FALSE]
  dat_deep <- dat_deep[compare_genes,,drop=FALSE]
  
  plot(dat_shallow[,1], dat_deep[,1])
  plot(dat_shallow[,2], dat_deep[,2])
  plot(dat_shallow[,3], dat_deep[,3])

  plot(log10(1+dat_shallow[,1]), log10(1+dat_deep[,1]))
  plot(log10(1+dat_shallow[,2]), log10(1+dat_deep[,2]))
  plot(log10(1+dat_shallow[,3]), log10(1+dat_deep[,3]))
  
  colnames(dat_shallow)[1]
  colnames(dat_deep)[1]
}








################################################################################
####################### Generate count using alignment -- barseq ###############   
################################################################################


listpools <- c(
  "EB_deepseq_barseqpool3",
  "EB_priming_barseqpool3"
)

####### Alignment part
for(curpool in listpools){
  print(curpool)
  allpooldir <- "/corgi/otherdataset/ellenbushell/barseq_pools"
  pooldir <- file.path(allpooldir, curpool)
  fastqdir <- file.path(pooldir, "fastq")
  bcfile <- file.path(pooldir, "used_bc.csv")
  bwa_file <- file.path(pooldir, "bwa_bc.fa")
  usedbc <- read.csv(bcfile,sep="\t")
  
  
  ################### Generate indices for bwa
  frank_bc <- read.csv("/corgi/otherdataset/ellenbushell/barcode_to_gene_210920_FRANK.csv")
  frank_bc$sgrna <- str_split_fixed(frank_bc$gene,"\\|",2)[,1]
  frank_bc$seq <- str_to_upper(frank_bc$barcode)
  frank_bc$seq_rc <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(frank_bc$seq))) 
  frank_bc <- frank_bc[frank_bc$sgrna %in% usedbc$gene,]
  
  tostore <- Biostrings::DNAStringSet(paste0(
    "CCTTCAATTTCGATGGGTACCACCCAGCTTTCTTGTACAAAGTGGTTGATATCTCTATAGTCGCAGTAGGCGG",
    frank_bc$seq_rc,
    "CTGACGCGCACGAATTACAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTC"))
  names(tostore) <- frank_bc$seq
  Biostrings::writeXStringSet(tostore, bwa_file)
  
  #Index
  the_cmd <- paste0("bwa index ",bwa_file)
  print(the_cmd)
  system(the_cmd)
  
  
  ################### Align with bwa
  list_bclist <- list()
  for(onef in list.files(fastqdir)){
    if(str_ends(onef,"001.fastq.gz") & str_detect(onef,"_S1_")){
      print(onef)
      
      #Align
      the_cmd <- paste0("bwa mem ",bwa_file," ",file.path(fastqdir,onef)," -t 10 -o ",file.path(fastqdir,onef),".sam")
      print(the_cmd)
      system(the_cmd)
      
      #Sort, make bam
      the_cmd <- paste0("samtools sort ",file.path(fastqdir,onef),".sam"," -o ",file.path(fastqdir,onef),".bam")
      print(the_cmd)
      system(the_cmd)
      
      #Index the bam
      the_cmd <- paste0("samtools index ",file.path(fastqdir,onef),".bam")
      print(the_cmd)
      system(the_cmd)
    }
  }
}


#curpool <- "EB_deepseq_barseqpool3"

####### Counting part
for(curpool in listpools){
  
  print(curpool)
  
  allpooldir <- "/corgi/otherdataset/ellenbushell/barseq_pools"
  pooldir <- file.path(allpooldir, curpool)
  fastqdir <- file.path(pooldir, "fastq")
  bcfile <- file.path(pooldir, "used_bc.csv")
  countfile <- file.path(pooldir,"counts_bwa.RDS")  #note, different!
  countfile_csv <- file.path(pooldir,"counts_bwa.csv")
  
  frank_bc <- read.csv("/corgi/otherdataset/ellenbushell/barcode_to_gene_210920_FRANK.csv")
  frank_bc$sgrna <- str_split_fixed(frank_bc$gene,"\\|",2)[,1]
  frank_bc$seq <- str_to_upper(frank_bc$barcode)
  
  list_bclist <- list()
  for(onef in list.files(fastqdir)){
    if(str_ends(onef,"R1_001.fastq.gz.bam")){
      print(onef)
      onef2 <- str_replace_all(onef, "R1","R2")  #Check also corresponding second read if present
      count_r1 <- read.table(pipe(paste0("samtools idxstats ",file.path(fastqdir,onef))))
      
      #      bclist <- sqldf::sqldf("select count(bc) as cnt, bc from bclist group by bc order by cnt desc")
      #      bclist$file <- onef
      
      onedf <- data.frame(bc=count_r1$V1, cnt=count_r1$V3)
      onedf$file <- onef
      list_bclist[[onef]] <- onedf
    }
  }
  counts <- do.call(rbind,list_bclist)
  counts <- reshape2::acast(counts, bc~file, value.var = "cnt", fill = 0)
  colSums(counts)
  rowSums(counts) #no *!
  
  #Which BC to expect here? Only keep these
  if(FALSE){
    usedbc <- read.csv(bcfile,sep="\t")
    usedbc <- frank_bc[frank_bc$sgrna %in% usedbc$gene,]
  } else {
    usedbc <- frank_bc
  }
  
  
  #Prepare column names (samples), and order by depth
  counts <- counts[order(rowSums(counts), decreasing = TRUE),,drop=FALSE]
  counts <- counts[rownames(counts) %in% usedbc$seq,,drop=FALSE]
  colnames(counts) <- str_split_fixed(colnames(counts),"_S",2)[,1]
  
  #Set rownames to gene name
  rownames(usedbc) <- usedbc$seq
  rownames(counts) <- usedbc[rownames(counts),,drop=FALSE]$sgrna
  
  if(TRUE){
    saveRDS(counts, countfile)
    write.csv(counts, countfile_csv)  #temp
  } 
}

























################################################################################
################ Generate count matrices -- barseq, exact matching #############
################################################################################

listpools <- c(
  "EB_priming_barseqpool3"
)

listpools <- c(
  "EB_organs1030_240707"  
)

listpools <- c(
  "EB_slowpool_staging_2024apr",
  "EB_slowpool_organs_2024apr"
)

listpools <- c(
  #  "barseq_minipool2",  #just a cnt file!
  "barseq_priming_Candidatepool1",
  "barseq_priming_Candidatepool2",
  "barseq_priming_barseqpool1",
  "EB_priming_barseqpool3",  #was renamed at some point
  "barseq_priming_barseqpool4",
  "barseq_slowhires_2023dec",
  "barseq_slowpool_1",
  "barseq_slowpool_2"
)
#listpools <- c("barseq_slowhires_2023dec")
for(curpool in listpools){
  
  print(curpool)
  
  # re.compile("GGCGG"+"(\w{8,16})"+"CTGAC")
  
  seqbefore <- "TAGTCGCAGTAGGCGG"
  
  allpooldir <- "/corgi/otherdataset/ellenbushell/barseq_pools"
  pooldir <- file.path(allpooldir, curpool)
  fastqdir <- file.path(pooldir, "fastq")
  bcfile <- file.path(pooldir, "used_bc.csv")
  countfile <- file.path(pooldir,"counts.RDS")
  countfile_csv <- file.path(pooldir,"counts.csv")
  
  frank_bc <- read.csv("/corgi/otherdataset/ellenbushell/barcode_to_gene_210920_FRANK.csv")
  frank_bc$sgrna <- str_split_fixed(frank_bc$gene,"\\|",2)[,1]
  frank_bc$seq <- str_to_upper(frank_bc$barcode)
  
  #Subset by the BCs expected here
  usedbc <- read.csv(bcfile,sep="\t")
  usedbc <- frank_bc[frank_bc$sgrna %in% usedbc$gene,]
  bclength <- str_length(usedbc$seq[1]) #  TCTTTTCCCAG
  
  #R1 needs RC
  usedbc$seq <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(usedbc$seq)))
  
  
  list_bclist <- list()
  for(onef in list.files(fastqdir)){
    if(str_ends(onef,"R1_001.fastq.gz")){
      print(onef)
      onep <- pipe(paste("zcat",file.path(fastqdir,onef),"| grep TAGTCGCAGTAGGCGG"))
      li <- readLines(onep)
      close(onep)

      onef2 <- str_replace_all(onef, "R1","R2")  #Corresponding second read
      
      ### TODO also use R2!
      
      #revcomp BC for barseq R1??
      
      bclist <- str_split_fixed(li, seqbefore,2)[,2]
      bclist <- data.frame(bc=str_sub(bclist,1,bclength))
      bclist <- sqldf::sqldf("select count(bc) as cnt, bc from bclist group by bc order by cnt desc")
      
      if(nrow(bclist)>0){
        bclist$file <- onef
        list_bclist[[onef]] <- bclist
      } else {
        print(paste("Skipping because no counts found in file",onef))
      }
    }
  }
  counts <- do.call(rbind,list_bclist)
  counts <- reshape2::acast(counts, bc~file, value.var = "cnt", fill = 0)
  colSums(counts)
  rowSums(counts)
  
  #
  counts <- counts[order(rowSums(counts), decreasing = TRUE),]
  counts <- counts[rownames(counts) %in% usedbc$seq,]
  #colnames(counts) <- str_sub(colnames(counts),1,11)  #dangerous. split by _S instead?
  colnames(counts) <- str_split_fixed(colnames(counts),"_S",2)[,1]
  
  rownames(usedbc) <- usedbc$seq
  rownames(counts) <- usedbc[rownames(counts),]$sgrna
 
  if(TRUE){
    saveRDS(counts, countfile)
    write.csv(counts, countfile_csv)  #temp
  } 
}

#write.csv(counts[rowSums(counts)>5,], "/home/mahogny/temp_cnt.csv") 



################################################################################
####################### Hack for CRISPR minipool2 ##############################
################################################################################

if(FALSE){
  counts <- read.csv("/corgi/otherdataset/ellenbushell/crispr_pools/barseq_minipool2/counts.csv", row.names = "X")
  counts <- counts[rownames(counts) != "_other" & rownames(counts) != "_err",]
  colnames(counts) <- str_split_fixed(colnames(counts),"\\.",2)[,1]
  saveRDS(counts, "/corgi/otherdataset/ellenbushell/crispr_pools/barseq_minipool2/counts.RDS")
}



################################################################################
####################### Generate count matrices -- CRISPR ######################
################################################################################

listpools <- c("cr_2024march_half1","cr_2024march_p1","cr_2024march_p12","cr_2024march_p2")

# "tags" = EB_barseq_slowpool_CRISPR


#listpools <- c("tags")
listpools <- c("2023march_screen_noD4")
listpools <- c("2023march_screen")
#listpools <- c("aug_p192","aug_p24","aug_p96")
listpools <- c("2023aug_p192","2023aug_p24","2023aug_p96","2023jan_tags", "2023march_screen", "2023march_screen_noD4")
listpools <- c("2023dec_ligpool")
for(curpool in listpools){
  
  print(curpool)
  
  seqbefore <- "CAATATTATT"
  
  allpooldir <- "/corgi/otherdataset/ellenbushell/crispr_pools"
  pooldir <- file.path(allpooldir, curpool)
  fastqdir <- file.path(pooldir, "fastq")
  bcfile <- file.path(pooldir, "bc.csv")
  countfile <- file.path(pooldir,"counts.RDS")
  
  usedbc <- read.csv(bcfile,sep="\t")
  bclength <- str_length(usedbc$seq[1])
  
  list_bclist <- list()
  for(onef in list.files(fastqdir)){
    print(onef)
    if(str_ends(onef,"R1_001.fastq.gz")){
      onep <- pipe(paste("zcat",file.path(fastqdir,onef),"| grep CAATATTATT"))
      li <- readLines(onep)
      close(onep)
      
      bclist <- str_split_fixed(li, seqbefore,2)[,2]
      bclist <- data.frame(bc=str_sub(bclist,1,bclength))
      bclist <- sqldf::sqldf("select count(bc) as cnt, bc from bclist group by bc order by cnt desc")
      bclist$file <- onef
      
      list_bclist[[onef]] <- bclist
    }
  }
  counts <- do.call(rbind,list_bclist)
  counts <- reshape2::acast(counts, bc~file, value.var = "cnt", fill = 0)

  #getting out the right barcodes  
  counts <- counts[order(rowSums(counts), decreasing = TRUE),]
  counts <- counts[rownames(counts) %in% usedbc$seq,]
  
  #Keep name up to _S##
  colnames(counts) <- str_split_fixed(colnames(counts),"_S",2)[,1]
  #colnames(counts) <- str_sub(colnames(counts),1,11) ##### only applies if
  
  rownames(usedbc) <- usedbc$seq
  rownames(counts) <- usedbc[rownames(counts),]$sgrna
  
  saveRDS(counts, countfile)
}



listpools <- c("cr_2024march_half1","cr_2024march_p1","cr_2024march_p12","cr_2024march_p2")
for(curpool in listpools){
  print(curpool)
  allpooldir <- "/corgi/otherdataset/ellenbushell/crispr_pools"
  pooldir <- file.path(allpooldir, curpool)
  countfile <- file.path(pooldir,"counts.RDS")
  counts <- readRDS(countfile)
  print(rowMeans(counts))
}


