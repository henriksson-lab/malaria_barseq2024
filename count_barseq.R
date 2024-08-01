library(ggplot2)
library(stringr)


################################################################################
####################### barseq, what reads look like ###########################
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


################################################################################
######## Generate count matrices -- barseq, improved approximate matching ######
################################################################################


listpools <- c(
  "EB_priming_Candidatepool1",
  "EB_priming_Candidatepool2",
  "EB_priming_barseqpool1",
  "EB_priming_barseqpool3", 
  "EB_priming_barseqpool4",
  "slowhires_2023dec",
  "EB_barseq_slowpool_1",
  "EB_barseq_slowpool_2",
  "EB_deepseq_barseqpool3"
)
for(curpool in listpools){
  
  print(paste("=========================================================",curpool))

  #align_to_R1 <- "CCTTCAATTT CGATGGGTAC CACCCAGCTT TCTTGTACAA AGTGGTTGAT ATCTCTATAG TCGCAGTAGG CGGNNNNNNN NNNNCTGACG CGCACGAATT ACAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCAGTGATT" #full seq
  full_R1 <- "CCTTCAATTTCGATGGGTACCACCCAGCTTTCTTGTACAAAGTGGTTGATATCTCTATAGTCGCAGTAGGCGGNNNNNNNNNNNCTGACGCGCACGAATTACAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCAGTGATT" #full seq
  #align_to_R1 <- "GGTTGATATCTCTATAGTCGCAGTAGGCGGNNNNNNNNNNNCTGACGCGCACGAATTACAGATCGGAAGAGC"
  align_to_R1 <- "TCGCAGTAGGCGGNNNNNNNNNNNCTGACGCGCACGAA"
  
  seqbefore <- "TAGTCGCAGTAGGCGG"
  
  allpooldir <- "/corgi/otherdataset/ellenbushell/barseq_pools"
  pooldir <- file.path(allpooldir, curpool)
  fastqdir <- file.path(pooldir, "fastq")
  bcfile <- file.path(pooldir, "used_bc.csv")
  countfile <- file.path(pooldir,"counts.v2.RDS")
  countfile_csv <- file.path(pooldir,"counts.v2.csv")
  
  frank_bc <- read.csv("/corgi/otherdataset/ellenbushell/barcode_to_gene_210920_FRANK.csv")
  frank_bc$sgrna <- str_split_fixed(frank_bc$gene,"\\|",2)[,1]
  frank_bc$seq <- str_to_upper(frank_bc$barcode)
  
  #Subset by the BCs expected here
  if(TRUE){
    usedbc <- read.csv(bcfile,sep="\t")  ### why not here??
    usedbc <- frank_bc[frank_bc$sgrna %in% usedbc$gene,]
  } else {
    usedbc <- frank_bc
  }
  bclength <- str_length(usedbc$seq[1]) #  TCTTTTCCCAG
  usedbc$seq <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(usedbc$seq))) ##R1 needs RC
  
  
  list_bclist <- list()
  for(onef in list.files(fastqdir)){
    if(str_ends(onef,"R1_001.fastq.gz")){ #& str_detect(onef,"_S1_"
      print(onef)
      onep <- pipe(paste("zcat",file.path(fastqdir,onef)))
      li <- readLines(onep)
      li <- li[seq(from=2,by=4,to=length(li))]
      close(onep)
      

      ######## Dirty approximate approach: See if correct position is anywhere near expected. If smallest score is unique then use this instead
      expected_position <- str_locate(full_R1,seqbefore)
      relpos <- -3:3
      distmat <- matrix(nrow=length(li), ncol=length(relpos))
      for(j in 1:length(relpos)){
        distmat[,j] <- adist(str_sub(li,expected_position[1]+relpos[j],expected_position[2]+relpos[j]),seqbefore)
      }
      number_of_mins <- rowSums(distmat == rowMins(distmat)) #see if there is a unique optimum
      shift_position <- relpos[apply(distmat,1,which.min)]
      
      
      #Only keep the ones with decent confidence; use other method for remaining
      keep_direct_test <- rowMins(distmat) <= 5 & number_of_mins==1  #62%
      mean(keep_direct_test)
      direct_detect <- li[keep_direct_test]
      no_detect <- li[!keep_direct_test]
      print(paste("Detected start for",sum(keep_direct_test),"but not",sum(!keep_direct_test),"reads"))
      
      bclist_direct <- data.frame(bc=str_sub(
        direct_detect, 
        shift_position[keep_direct_test] + expected_position[2] + 1,
        shift_position[keep_direct_test] + expected_position[2] + bclength
      ))
      
      if(FALSE){
        table(rowMins(distmat))
        table(shift_position)
      }
      
      #Very slow
      if(FALSE){
        ######## Pairwise alignment of each read to the pattern expected; extract then BC
        print(paste("Getting BCs by alignment",length(no_detect)))
        outp <- pwalign::pairwiseAlignment(
          #str_sub(no_detect,40,110),
          str_sub(no_detect,50,100),
          align_to_R1
        )
        trimmed_bc <- rep(NA, length(outp))
        for(i in 1:length(outp)){ #28s just get aligned vectorized, to 36s for it all in a loop
          one_outp <- outp[i]
          ssubject <- as.character(pwalign::alignedSubject(one_outp)) #significantly faster if trimming input
          spattern <- as.character(pwalign::alignedPattern(one_outp))
          toget <- stringr::str_locate_all(ssubject, "N")[[1]]
          trimmed_bc[i] <- str_sub(spattern,toget[1,1], toget[nrow(toget),1])
        }
        trimmed_bc <- str_remove_all(trimmed_bc,"-")
      }

      ###################### Put together the two list of BCs, and count   
      print("Counting BCs")
      bclist_aligned <- data.frame(bc=trimmed_bc)
      bclist <- bclist_direct
      #bclist <- rbind(bclist_direct, bclist_aligned)
      bclist <- sqldf::sqldf("select count(bc) as cnt, bc from bclist group by bc order by cnt desc")

      
      ###################### Correct the BCs that are not recognized
      print("Correcting BCs")
      bclist$in_list <- bclist$bc %in% usedbc$seq
      bclist$corrected_bc <- bclist$bc
      bclist$score <- NA
      print(paste("BC in list:",sum(bclist$in_list), " ----- BC not in list:",sum(!bclist$in_list)))
      for(i in which(!bclist$in_list)){
        #print(paste("correcting ",i))
        scored_bc_match <- adist(usedbc$seq, bclist$bc[i])[,1]
        wmin <- which.min(scored_bc_match)
        bclist$score[i] <- scored_bc_match[wmin]
        bclist$corrected_bc[i] <- usedbc$seq[wmin]
        if(sum(scored_bc_match==wmin)>1){
          bclist$score[i] <- 15 #if more than one is the most similar, put full penalty
        }
        if(i%%1000 == 0) print(i)
      }

      ## Correct BCs if good enough
      #sort(bclist$score)
      hist(bclist$score)
      table(bclist$score)
      bclist_keep <- bclist[is.na(bclist$score) | bclist$score <= 3,,drop=FALSE]  #almost all are <=3

      ## Map grna -> gene
      #Note: Some grnas for the same gene. But we can safely sum these up to avoid issues; only one grna dominates
      rownames(usedbc) <- usedbc$seq
      bclist_keep$gene <- usedbc[bclist_keep$corrected_bc,,drop=FALSE]$sgrna
      
      print("Recounting BCs after correction")
      bclist_keep <- sqldf::sqldf("select sum(cnt) as cnt, gene from bclist_keep group by gene order by cnt desc")
      if(FALSE){
        plot(bclist_keep$cnt)
        plot(bclist_keep$cnt)
        bclist_keep$gene
      }

      ## Only store samples having any counts
      if(nrow(bclist_keep)>0){
        bclist_keep$file <- onef
        list_bclist[[onef]] <- bclist_keep
      } else {
        print(paste("Skipping because no counts found in file",onef))
      }
    }
  }
  counts <- do.call(rbind,list_bclist)
  counts <- reshape2::acast(counts, gene~file, value.var = "cnt", fill = 0)

  counts <- counts[order(rowSums(counts), decreasing = TRUE),,drop=FALSE]
  #counts <- counts[rownames(counts) %in% usedbc$seq,,drop=FALSE]
  colnames(counts) <- str_split_fixed(colnames(counts),"_S",2)[,1]

  saveRDS(counts, countfile)
  write.csv(counts, countfile_csv)
}





