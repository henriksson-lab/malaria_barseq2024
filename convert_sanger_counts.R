 "/corgi/otherdataset/ellenbushell/barseq_pools/sanger_all"
 #/corgi/otherdataset/ellenbushell/barseq_pools/sanger_all

 #d.4, 5, 6, 7 -> WT vs Rag, P vs NP -> x2 mice
 #counts_25253.csv  PbPbSTM140 - Prime barseq PCR2
 #counts_25302.csv  PbPbSTM139 - Prime barseq PCR1  
 #counts_25792.csv  PbPbSTM145 - Prime barseq PCR2_repeat
 #counts_25800.csv  PbPbSTM144 - Prime barseq PCR1_repeat

 
 #d.4, 5, 6, 7 -> WT vs Rag, P vs NP -> x 4 mice
 #counts_26059.csv  PbSTM155   merge as PCR1
 #counts_26072.csv  PbSTM156
 
 #d.4, 5, 6, 7 -> WT vs Rag, P vs NP -> x 4 mice
 #counts_26073.csv  PbSTM157   merge as PCR2
 #counts_26080.csv  PbSTM158
 

library(stringr)

remove_outliers <- c(
 "PbPbSTM139_P_Jax_d7_m1_r1",
 "PbSTM140_NP_Rag_d5_m2_r2",
 "PbSTM158_NP_Jax_d4_m2_PCR2_2_r2",
 "PbSTM158_NP_Jax_d5_m2_PCR2_2_r2",
 "PbSTM158_NP_Jax_d6_m2_PCR2_2_r2",
 "PbSTM158_NP_Jax_d7_m2_PCR2_2_r2",
 "PbSTM158_NP_Jax_d7_m3_PCR2_2_r2"
)

for(f in list.files("/corgi/otherdataset/ellenbushell/barseq_pools/sanger_all/orig_counts")){
 dat <- read.csv(file.path("/corgi/otherdataset/ellenbushell/barseq_pools/sanger_all/orig_counts",f))
 print(f)
 
 for(cname in colnames(dat)[str_ends(colnames(dat),"2")]){
   cname2 <- cname
   cname1 <- str_replace(cname2,fixed(".2"),".1")
   dat[,cname1] <- dat[,cname1]+dat[,cname2]
 }
 dat <- dat[,!str_ends(colnames(dat),"2")]
 colnames(dat) <- str_remove_all(colnames(dat),fixed(".1"))
 
 #print(unique(str_split_fixed(colnames(dat),"_",2)[,1]))
 #print(nrow(dat))
 
 dat <- dat[dat$barcode!="no_match",]
 rownames(dat) <- dat$barcode  #note, column "gene" can be duplicated!
 dat <- dat[,-(1:2)]
 
 #Remove outliers
 to_keep <- !(colnames(dat) %in% remove_outliers)
 dat <- dat[,to_keep]
 
 
 #Fix odd naming
 colnames(dat) <- str_replace(colnames(dat),"_day","_d")
 colnames(dat) <- str_replace(colnames(dat),"_r1_","_R1_") #not the read but RAG

 #read1 ends in .1, and read2 in .2
 sampleinfo <- data.frame(NGI.ID=colnames(dat), User.ID=colnames(dat))
 sampleinfo$Mreads <- 666
 sampleinfo$X30 <- 666
 sampleinfo$mouse <- str_sub(str_split_fixed(colnames(dat),"_m",2)[,2],1,1)
 sampleinfo$day <- str_sub(str_split_fixed(colnames(dat),"_d",2)[,2],1,1)
 sampleinfo$primed <- "NP" #default  --- wrong???? 4 pools are primed?
 sampleinfo$primed[str_detect(colnames(dat),"_PP_")] <- "P"
 sampleinfo$primed[str_detect(colnames(dat),"_UP_")] <- "NP"
 sampleinfo$is_input <- FALSE
 sampleinfo$is_input[str_detect(colnames(dat),"_Input_")] <- TRUE
 #no sampleinfo$is_input[str_detect(colnames(dat),"_r1_")] <- TRUE  #yes, weird
 sampleinfo$genotype <- ""
 sampleinfo$genotype[str_detect(colnames(dat),"JACS")] <- "BL6"
 sampleinfo$genotype[str_detect(colnames(dat),"RAGG")] <- "RAG1KO"
 sampleinfo$genotype[str_detect(colnames(dat),"_J")] <- "BL6"
 sampleinfo$genotype[str_detect(colnames(dat),"_R")] <- "RAG1KO"
 
 sampleinfo$primed[str_detect(colnames(dat),"J1")] <- "NP"
 sampleinfo$primed[str_detect(colnames(dat),"J2")] <- "NP"
 sampleinfo$primed[str_detect(colnames(dat),"J3")] <- "P"
 sampleinfo$primed[str_detect(colnames(dat),"J4")] <- "P"
 
 sampleinfo$primed[str_detect(colnames(dat),"R1")] <- "NP"
 sampleinfo$primed[str_detect(colnames(dat),"R2")] <- "NP"
 sampleinfo$primed[str_detect(colnames(dat),"R3")] <- "P"
 sampleinfo$primed[str_detect(colnames(dat),"R4")] <- "P"
 
 for(i in 1:4){
   sampleinfo$mouse[str_detect(colnames(dat),paste0("J",i))] <- i
   sampleinfo$mouse[str_detect(colnames(dat),paste0("R",i))] <- i
 }
 
 #Default to unprimed?
 
 sampleinfo$User.ID <- paste0("poolsanger_",sampleinfo$primed,"_",sampleinfo$genotype,"_d",sampleinfo$day,"_m",sampleinfo$mouse)
 sampleinfo$User.ID[sampleinfo$is_input] <- paste0("poolsanger_input")
 colnames(dat) <- sampleinfo$User.ID
 
 write.csv(sampleinfo, file.path("/corgi/otherdataset/ellenbushell/barseq_pools/sanger_all/newmeta",f))
}



 
 
 
 
###################### alternative way; renaming table is incomplete
###################### alternative way; renaming table is incomplete
###################### alternative way; renaming table is incomplete

 
 
 
 
 
 
 
tab_replace <- read.csv("/corgi/otherdataset/ellenbushell/barseq_pools/sanger_all/renaming.csv",sep="\t")
tab_replace$to <- str_remove(tab_replace$to,"_r1")
tab_replace$to <- str_remove(tab_replace$to,"_r2")

for(f in list.files("/corgi/otherdataset/ellenbushell/barseq_pools/sanger_all/orig_counts")){
 dat <- read.csv(file.path("/corgi/otherdataset/ellenbushell/barseq_pools/sanger_all/orig_counts",f))
 print(f)
 
 #Sum up R1 and R2
 for(cname in colnames(dat)[str_ends(colnames(dat),"2")]){
   cname2 <- cname
   cname1 <- str_replace(cname2,fixed(".2"),".1")
   dat[,cname1] <- dat[,cname1]+dat[,cname2]
 }
 dat <- dat[,!str_ends(colnames(dat),"2")]
 colnames(dat) <- str_remove_all(colnames(dat),fixed(".1"))
 
 dat <- dat[dat$barcode!="no_match",]
 rownames(dat) <- dat$barcode  #note, column "gene" can be duplicated!
 dat <- dat[,-(1:2)]

 sampleinfo <- data.frame(NGI.ID=colnames(dat), User.ID=colnames(dat))
 sampleinfo$Mreads <- 666
 sampleinfo$X30 <- 666
 
 for(i in 1:nrow(sampleinfo)){
   cname <- sampleinfo$NGI.ID[i]
   for(j in 1:nrow(tab_replace)){
     if(str_starts(cname,tab_replace$from[j])){
       cname <- tab_replace$to[j]
     }
   sampleinfo$User.ID[i] <- cname
   }
   
 }
 
 #write.csv(sampleinfo, file.path("/corgi/otherdataset/ellenbushell/barseq_pools/sanger_all/newmeta",f))
}

















#############################################################
############################################################# final run
#############################################################





remove_outliers <- c(
  "PbPbSTM139_P_Jax_d7_m1_r1",
  "PbSTM140_NP_Rag_d5_m2_r2",
  "PbSTM158_NP_Jax_d4_m2_PCR2_2_r2",
  "PbSTM158_NP_Jax_d5_m2_PCR2_2_r2",
  "PbSTM158_NP_Jax_d6_m2_PCR2_2_r2",
  "PbSTM158_NP_Jax_d7_m2_PCR2_2_r2",
  "PbSTM158_NP_Jax_d7_m3_PCR2_2_r2"
)


#invalid_bc <- c(
#  "taagttcgat","ttcagctcat",
#  "aaaaggggtct","caccagcaccc",
#  "tcctcaatat",
#  "aatgaagagtc",
#  "ccgcaccgctt",
#  "tctcggttat")

invalid_bc <- c(
  "ttccaccttac" #this gene is actually present twice. removing low abundant one
)

read_count_file <- function(f){
  dat <- read.csv(file.path("/corgi/otherdataset/ellenbushell/barseq_pools/sanger_all/orig_counts",f))
  print(f)

  for(cname in colnames(dat)[str_ends(colnames(dat),"2")]){
    cname2 <- cname
    cname1 <- str_replace(cname2,fixed(".2"),".1")
    dat[,cname1] <- dat[,cname1]+dat[,cname2]
  }
  dat <- dat[,!str_ends(colnames(dat),"2")]
  colnames(dat) <- str_remove_all(colnames(dat),fixed(".1"))

  #print(unique(str_split_fixed(colnames(dat),"_",2)[,1]))
  #print(nrow(dat))
  
  #dat_other <- dat[dat$barcode!="no_match",,drop=FALSE]
  #dat <- dat[dat$barcode!="no_match",]
  #rownames(dat) <- dat$barcode  #note, column "gene" can be duplicated!
  
  ### Map BC to gene
  dat$gene[dat$barcode=="no_match"] <- "_other"
  dat <- dat[!(dat$barcode %in% invalid_bc),] #Remove a duplicate PBANKA_030600
  #rownames(dat) <- dat$gene
  #print("Missing genes:")
  #print(missing_genes)

  #Check if all genes in old table
  frank_bc <- read.csv("/corgi/otherdataset/ellenbushell/barcode_to_gene_210920_FRANK.csv")
  frank_bc$sgrna <- str_split_fixed(frank_bc$gene,"\\|",2)[,1]
  #frank_bc$seq <- str_to_upper(frank_bc$barcode)
  rownames(frank_bc) <- frank_bc$barcode
  missing_bc <- dat$barcode[!(dat$barcode %in% frank_bc$barcode)]
  if(length(missing_bc)>1){
    #print("Missing BCs")
    #print(missing_bc)
    #print(dat$gene[!(dat$barcode %in% frank_bc$barcode)])
  }
  #print(frank_bc[dat$barcode,]$sgrna==dat$gene)
  
  ### Look up gene name anew; remove ones we cannot tell
  dat$gene <- frank_bc[dat$barcode,]$sgrna
  dat$gene[dat$barcode=="no_match"] <- "_other"
  dat <- dat[!is.na(dat$gene),]
  #print(dat)
  #print(666)
  #print(777)
  #print(dat[dat$gene=="PBANKA_030600",])
  rownames(dat) <- dat$gene

  ### Properly format table  
  dat <- dat[,-(1:2)]

  #Remove outliers
  to_keep <- !(colnames(dat) %in% remove_outliers)
  dat <- dat[,to_keep]

  #Fix odd naming
  colnames(dat) <- str_replace(colnames(dat),"_day","_d")
  colnames(dat) <- str_replace(colnames(dat),"_r1_","_R1_") #not the read but RAG
  
  return(dat)
}


#read_count_file("counts_26059.csv")
#'PBANKA_093370', 'PBANKA_101330' 


store_count_file <- function(dat, tofile){
  
  orig_name <- colnames(dat)
  
  #read1 ends in .1, and read2 in .2
  sampleinfo <- data.frame(NGI.ID=colnames(dat), User.ID=colnames(dat))
  sampleinfo$Mreads <- 666
  sampleinfo$X30 <- 666
  sampleinfo$mouse <- str_sub(str_split_fixed(colnames(dat),"_m",2)[,2],1,1)
  sampleinfo$day <- str_sub(str_split_fixed(colnames(dat),"_d",2)[,2],1,1)
  sampleinfo$primed <- "" #default  --- wrong???? 4 pools are primed?
  sampleinfo$primed[str_detect(colnames(dat),"_PP_")] <- "P"
  sampleinfo$primed[str_detect(colnames(dat),"_UP_")] <- "NP"
  sampleinfo$is_input <- FALSE
  sampleinfo$is_input[str_detect(colnames(dat),"_Input_")] <- TRUE
  #no sampleinfo$is_input[str_detect(colnames(dat),"_r1_")] <- TRUE  #yes, weird
  sampleinfo$genotype <- ""
  sampleinfo$genotype[str_detect(colnames(dat),"JACS")] <- "BL6"
  sampleinfo$genotype[str_detect(colnames(dat),"RAGG")] <- "RAG1KO"
  sampleinfo$genotype[str_detect(colnames(dat),"_J")] <- "BL6"
  sampleinfo$genotype[str_detect(colnames(dat),"_R")] <- "RAG1KO"
  
  sampleinfo$primed[str_detect(colnames(dat),"_J1")] <- "NP"
  sampleinfo$primed[str_detect(colnames(dat),"_J2")] <- "NP"
  sampleinfo$primed[str_detect(colnames(dat),"_J3")] <- "P"
  sampleinfo$primed[str_detect(colnames(dat),"_J4")] <- "P"
  
  sampleinfo$primed[str_detect(colnames(dat),"_R1")] <- "NP"
  sampleinfo$primed[str_detect(colnames(dat),"_R2")] <- "NP"
  sampleinfo$primed[str_detect(colnames(dat),"_R3")] <- "P"
  sampleinfo$primed[str_detect(colnames(dat),"_R4")] <- "P"
  
  for(i in 1:4){
    sampleinfo$mouse[str_detect(colnames(dat),paste0("_J",i))] <- i
    sampleinfo$mouse[str_detect(colnames(dat),paste0("_R",i))] <- i
  }
  
  #Default to unprimed?
  
  sampleinfo$User.ID <- paste0("poolsanger_",sampleinfo$primed,"_",sampleinfo$genotype,"_d",sampleinfo$day,"_m",sampleinfo$mouse)
  sampleinfo$User.ID[sampleinfo$is_input] <- sprintf("poolsanger_input_%s",1:sum(sampleinfo$is_input))
  if(any(duplicated(sampleinfo$User.ID))){
    print("duplicated sample IDs")
    print(sampleinfo)
    error()
  }
  sampleinfo$NGI.ID <- sampleinfo$User.ID
  colnames(dat) <- sampleinfo$User.ID
  
  
  
  outdir <- file.path("/corgi/otherdataset/ellenbushell/barseq_pools/",tofile)
  if(!file.exists(outdir)){
    dir.create(outdir)
  }

  write.csv(dat, file.path(outdir,"counts.csv"), quote = FALSE)
  write.table(sampleinfo[,c("NGI.ID","User.ID","Mreads","X30")], file.path(outdir,"sampleinfo.txt"), quote = FALSE, row.names = FALSE, sep="\t")
  saveRDS(dat, file.path(outdir,"counts.RDS"))
  #note, has _other; keep?
  
  
  formaria <- data.frame(oldid=orig_name, newid=sampleinfo$User.ID)
  formaria$Mreads <- 666
  formaria$X30 <- 666
  write.table(formaria, file.path("/corgi/otherdataset/ellenbushell/barseq_pools/sanger_all/newmeta/",paste0(tofile,".csv")), quote = FALSE, row.names = FALSE, sep="\t")
  
  system(paste("cp /corgi/otherdataset/ellenbushell/barseq_pools/sanger_all/controls.csv ",file.path(outdir,"controls.csv")))
}



#dat1 <- read_count_file("counts_26059.csv")
#dat2 <- read_count_file("counts_26072.csv")
#rbind(reshape2::melt(dat1)

#read_count_file("counts_26059.csv")

#d.4, 5, 6, 7 -> WT vs Rag, P vs NP -> x2 mice
#counts_25253.csv  PbPbSTM140 - Prime barseq PCR2
#counts_25302.csv  PbPbSTM139 - Prime barseq PCR1  
#counts_25792.csv  PbPbSTM145 - Prime barseq PCR2_repeat
#counts_25800.csv  PbPbSTM144 - Prime barseq PCR1_repeat


#d.4, 5, 6, 7 -> WT vs Rag, P vs NP -> x 4 mice
#counts_26059.csv  PbSTM155   merge as PCR1
#counts_26072.csv  PbSTM156

#d.4, 5, 6, 7 -> WT vs Rag, P vs NP -> x 4 mice
#counts_26073.csv  PbSTM157   merge as PCR2
#counts_26080.csv  PbSTM158

store_count_file(read_count_file("counts_25253.csv"),"sanger_primed_barseq_PCR2")
store_count_file(read_count_file("counts_25302.csv"),"sanger_primed_barseq_PCR1")
store_count_file(read_count_file("counts_25792.csv"),"sanger_primed_barseq_PCR2_repeat")
store_count_file(read_count_file("counts_25800.csv"),"sanger_primed_barseq_PCR2_repeat")

store_count_file(read_count_file("counts_26059.csv"),"sanger_some_PCR1a")
store_count_file(read_count_file("counts_26072.csv"),"sanger_some_PCR1b")
store_count_file(read_count_file("counts_26073.csv"),"sanger_some_PCR2a")
store_count_file(read_count_file("counts_26080.csv"),"sanger_some_PCR2b")










################################################################################
################################################################################
################# plug in new metadata, maria corrected ########################
################################################################################
################################################################################
################################################################################






store_count_file_fixed <- function(dat, usemeta, tofile, rename_mice=FALSE){
  
  orig_name <- colnames(dat)
  
  renaming_table <- read.csv(file.path("/corgi/otherdataset/ellenbushell/barseq_pools/sanger_all/fixedmeta/",paste0(usemeta,".csv")), sep="\t")

  if(rename_mice){
    renaming_table$newid <- str_replace(renaming_table$newid,"_m1","_m3")
    renaming_table$newid <- str_replace(renaming_table$newid,"_m2","_m4")
  }
  
  #Remove some samples
  print(colnames(dat))
  print(renaming_table)
  dat <- dat[,colnames(dat) %in% renaming_table$oldid,drop=FALSE]
  print(colnames(dat))
    
  sampleinfo <- data.frame(NGI.ID=colnames(dat), User.ID=colnames(dat))

  rownames(renaming_table) <- renaming_table$oldid
  sampleinfo$NGI.ID <- renaming_table[sampleinfo$NGI.ID,]$newid
  sampleinfo$User.ID <- sampleinfo$NGI.ID
  colnames(dat) <- sampleinfo$User.ID
  
  sampleinfo$Mreads <- 666
  sampleinfo$X30 <- 666
  

  outdir <- file.path("/corgi/otherdataset/ellenbushell/barseq_pools/",tofile)
  if(!file.exists(outdir)){
    dir.create(outdir)
  }
  
  write.csv(dat, file.path(outdir,"counts.csv"), quote = FALSE)
  write.table(sampleinfo[,c("NGI.ID","User.ID","Mreads","X30")], file.path(outdir,"sampleinfo.txt"), quote = FALSE, row.names = FALSE, sep="\t")
  saveRDS(dat, file.path(outdir,"counts.RDS"))
  
  system(paste("cp /corgi/otherdataset/ellenbushell/barseq_pools/sanger_all/controls.csv ",file.path(outdir,"controls.csv")))
}


#dat <- read_count_file("counts_25253.csv")
#tofile <- "sanger_primed_barseq_PCR2"

store_count_file_fixed(read_count_file("counts_25302.csv"),"sanger_primed_barseq_PCR1",        "EB_priming_barseqpool2s_biorep1_PCR1a")
store_count_file_fixed(read_count_file("counts_25253.csv"),"sanger_primed_barseq_PCR2",        "EB_priming_barseqpool2s_biorep1_PCR1b", rename_mice=TRUE) #concatenate these; 4 mice total

store_count_file_fixed(read_count_file("counts_25792.csv"),"sanger_primed_barseq_PCR1_repeat", "EB_priming_barseqpool2s_biorep1_PCR2a") 
store_count_file_fixed(read_count_file("counts_25800.csv"),"sanger_primed_barseq_PCR2_repeat", "EB_priming_barseqpool2s_biorep1_PCR2b", rename_mice=TRUE) #all of these should be precisely the same. 2m each?


store_count_file_fixed(read_count_file("counts_26059.csv"),"sanger_some_PCR1a",                "EB_priming_barseqpool2s_biorep2_PCR1a")
store_count_file_fixed(read_count_file("counts_26072.csv"),"sanger_some_PCR1b",                "EB_priming_barseqpool2s_biorep2_PCR1b") #concatenate these; two sep cond. 4 mice total

store_count_file_fixed(read_count_file("counts_26073.csv"),"sanger_some_PCR2a",                "EB_priming_barseqpool2s_biorep2_PCR2a") 
store_count_file_fixed(read_count_file("counts_26080.csv"),"sanger_some_PCR2b",                "EB_priming_barseqpool2s_biorep2_PCR2b") #concatenate these; two sep cond. 4 mice total





################################
################################
################################ Merge pool PCR1 & PCR2 - EB_priming_barseqpool2s_PCR1+2
################################
################################


library(reshape2)


concat_straight <- function(f1, f2){
  count1 <- readRDS(file.path("/corgi/otherdataset/ellenbushell/barseq_pools",f1,"counts.RDS"))
  count2 <- readRDS(file.path("/corgi/otherdataset/ellenbushell/barseq_pools",f2,"counts.RDS"))
  totc <- rbind(melt(as.matrix(count1)), melt(as.matrix(count2)))
  dat <- acast(totc, Var1~Var2, fill = 0)
  dat
}

sum_straight <- function(f1, f2){
  count1 <- readRDS(file.path("/corgi/otherdataset/ellenbushell/barseq_pools",f1,"counts.RDS"))
  count2 <- readRDS(file.path("/corgi/otherdataset/ellenbushell/barseq_pools",f2,"counts.RDS"))
  totc <- rbind(melt(as.matrix(count1)), melt(as.matrix(count2)))
  totc <- sqldf::sqldf("select Var1, Var2, sum(value) as value from totc group by Var1, Var2") #sum them up
  dat <- acast(totc, Var1~Var2, fill = 0)
  dat
}

sum_straight_from_c <- function(count1, count2){
  totc <- rbind(melt(as.matrix(count1)), melt(as.matrix(count2)))
  totc <- sqldf::sqldf("select Var1, Var2, sum(value) as value from totc group by Var1, Var2") #sum them up
  dat <- acast(totc, Var1~Var2, fill = 0)
  dat
}


store_straight <- function(dat, tofile){
  outdir <- file.path("/corgi/otherdataset/ellenbushell/barseq_pools/",tofile)
  if(!file.exists(outdir)){
    dir.create(outdir)
  }
  write.csv(dat, file.path(outdir,"counts.csv"), quote = FALSE)
  sampleinfo <- data.frame(NGI.ID=colnames(dat),User.ID=colnames(dat))
  sampleinfo$Mreads <- 666
  sampleinfo$X30 <- 666
  write.table(sampleinfo, file.path(outdir,"sampleinfo.txt"), quote = FALSE, row.names = FALSE, sep="\t")
  saveRDS(dat, file.path(outdir,"counts.RDS"))
  
  system(paste("cp /corgi/otherdataset/ellenbushell/barseq_pools/sanger_all/controls.csv ",file.path(outdir,"controls.csv")))
}


######
###### biorep #1
######


###### These two PCRs of separate samples; make it into comparable libraries
store_straight(
  concat_straight("EB_priming_barseqpool2s_biorep1_PCR1a","EB_priming_barseqpool2s_biorep1_PCR1b"),
  "EB_priming_barseqpool2s_biorep1_PCR1")
store_straight(
  concat_straight("EB_priming_barseqpool2s_biorep1_PCR2a","EB_priming_barseqpool2s_biorep1_PCR2b"),
  "EB_priming_barseqpool2s_biorep1_PCR2")

###### PCR of the same thing, so can sum them up
store_straight(
  sum_straight("EB_priming_barseqpool2s_biorep1_PCR1","EB_priming_barseqpool2s_biorep1_PCR2"),
  "EB_priming_barseqpool2s_biorep1")


######
###### biorep #2
######

###### These two PCRs of separate samples; make it into comparable libraries
store_straight(
  concat_straight("EB_priming_barseqpool2s_biorep2_PCR1a","EB_priming_barseqpool2s_biorep2_PCR1b"),
  "EB_priming_barseqpool2s_biorep2_PCR1")
store_straight(
  concat_straight("EB_priming_barseqpool2s_biorep2_PCR2a","EB_priming_barseqpool2s_biorep2_PCR2b"),
  "EB_priming_barseqpool2s_biorep2_PCR2")

###### PCR of the same thing, so can sum them up
store_straight(
  sum_straight("EB_priming_barseqpool2s_biorep2_PCR1","EB_priming_barseqpool2s_biorep2_PCR2"),
  "EB_priming_barseqpool2s_biorep2")

######
###### biorep #1+#2 --  so 8 mice!!?
######

concat_mice <- function(f1, f2){
  count1 <- readRDS(file.path("/corgi/otherdataset/ellenbushell/barseq_pools",f1,"counts.RDS"))
  count2 <- readRDS(file.path("/corgi/otherdataset/ellenbushell/barseq_pools",f2,"counts.RDS"))
  
  colnames(count2) <- str_replace(colnames(count2),"_m1","_m5")
  colnames(count2) <- str_replace(colnames(count2),"_m2","_m6")
  colnames(count2) <- str_replace(colnames(count2),"_m3","_m7")
  colnames(count2) <- str_replace(colnames(count2),"_m4","_m8")
  
  totc <- rbind(melt(as.matrix(count1)), melt(as.matrix(count2)))
  dat <- acast(totc, Var1~Var2, fill = 0)
  dat
}

store_straight(
  concat_mice("EB_priming_barseqpool2s_biorep1","EB_priming_barseqpool2s_biorep2"),
  "EB_priming_barseqpool2s")
