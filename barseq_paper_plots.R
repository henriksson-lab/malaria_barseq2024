
library(stringr)
library(ggplot2)
library(patchwork)
library(sqldf)
library(reshape2)
library(cowplot)

make_white_plot <- function(p){
  p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
}


################################################################################
########### Comparison of shallow vs deep seq data #############################
################################################################################

outdir <- "/corgi/otherdataset/ellenbushell/newplots_barseq/comparison_depth"

############# 
############# Comparison on the level of counts -- look at a single sample (#1)
############# 

dat_shallow <- readRDS("/corgi/otherdataset/ellenbushell/barseq_pools/EB_priming_barseqpool3/counts.v2.RDS")
dat_shallow <- dat_shallow[,order(colnames(dat_shallow)),drop=FALSE]

dat_deep <- readRDS("/corgi/otherdataset/ellenbushell/barseq_pools/EB_deepseq_barseqpool3/counts.v2.RDS")
dat_deep <- dat_deep[,order(colnames(dat_deep)),drop=FALSE]


compare_genes <- sort(intersect(rownames(dat_shallow), rownames(dat_deep)))
dat_shallow <- dat_shallow[compare_genes,,drop=FALSE]
dat_deep <- dat_deep[compare_genes,,drop=FALSE]
ncol(dat_shallow)==ncol(dat_shallow)

df <- data.frame(
  shallow=(1+log10(dat_shallow[,1])),
  deep=(1+log10(dat_deep[,1]))
)
ptot <- make_white_plot(
  ggplot(df, aes(shallow, deep)) + 
  geom_point() + 
  xlab("log10 1+shallow") + 
  ylab("log10 1+deep"))
ggsave(plot=ptot, file.path(outdir,"sample1_counts.pdf"), width = 4, height = 4)

# + scale_x_log10() + scale_y_log10()
#plot(dat_shallow[,1], dat_deep[,1])
#plot(dat_shallow[,2], dat_deep[,2])
#plot(dat_shallow[,3], dat_deep[,3])


############# 
############# Comparison on the level of FC: P RAG1KO
############# 

dat_shallow <- all_grstats$EB_priming_barseqpool3$stats_per_grna$`P RAG1KO`
dat_deep    <- all_grstats$EB_deepseq_barseqpool3$stats_per_grna$`P RAG1KO`

df <- merge(
  data.frame(grna=dat_shallow$grna, fc_shallow=dat_shallow$fc, genedesc=dat_shallow$genecat),
  data.frame(grna=dat_deep$grna, fc_deep=dat_deep$fc), 
  all=TRUE #same coverage
)  
ptot <- make_white_plot(
  ggplot(df, aes(fc_shallow, fc_deep, label=grna, color=genedesc)) + 
  geom_point() + 
  xlab("Shallow: FC RAG1KO")+
  ylab("Deep: FC RAG1KO")) 
ggsave(plot=ptot, file.path(outdir,"sample1_fc_rag1ko.pdf"), width = 4, height = 4)

make_white_plot(ggplot(df, aes(fc_shallow, fc_deep, label=grna, color=genedesc)) + geom_point(color="gray") + geom_text())



############# 
############# Comparison on the level of FC: P BL6 / P RAG1KO
############# 

dat_shallow <- all_grstats$EB_priming_barseqpool3$scatterplot$`P BL6 / P RAG1KO`
dat_deep    <- all_grstats$EB_deepseq_barseqpool3$scatterplot$`P BL6 / P RAG1KO`
df <- merge(
  data.frame(grna=dat_shallow$gene, fc_shallow=dat_shallow$diff_fc, genedesc=dat_shallow$genedesc),
  data.frame(grna=dat_deep$gene, fc_deep=dat_deep$diff_fc), 
  all=TRUE #same coverage
)  
ptot <- make_white_plot(
  ggplot(df, aes(fc_shallow, fc_deep, label=grna, color=genedesc)) + 
    geom_point() + 
    xlab("Shallow: P BL6 / P RAG1KO")+
    ylab("Deep: P BL6 / P RAG1KO")) 
ggsave(plot=ptot, file.path(outdir,"sample1_fc_delta.pdf"), width = 4, height = 4)

#ggplot(df, aes(fc_shallow, fc_deep, label=grna, color=genedesc)) + geom_point()
ggplot(df, aes(fc_shallow, fc_deep, label=grna, color=genedesc)) + geom_point(color="gray") + geom_text()



################################################################################
########### Comparison of rounds ###############################################
################################################################################



list_round1 <- c(
  #Initial 4 pools
  "EB_priming_barseqpool1",
  "EB_priming_barseqpool2s",
  "EB_priming_barseqpool3", 
  "EB_priming_barseqpool4"
)

list_round2 <- c(
  #Two biological replicates, picked from round #1
  "EB_priming_Candidatepool1",
  "EB_priming_Candidatepool2"
)

list_round3 <- c(
  #The final pool, few mutants
  "EB_minipool2"
)


#different proj
list_not_included <- c(
  "EB_barseq_slowpool_1",
  "EB_barseq_slowpool_2",
  "slowhires_2023dec"  #subset of above
)


### Function to compare delta between two pools
compare_two_delta <- function(pool1, pool2, cond){
  dat1 <- all_grstats[[pool1]]$scatterplot[[cond]]
  dat2 <- all_grstats[[pool2]]$scatterplot[[cond]]
  
  df <- merge(
    data.frame(grna=dat1$gene, fc1=dat1$diff_fc, genedesc=dat1$genedesc),
    data.frame(grna=dat2$gene, fc2=dat2$diff_fc), 
    all=TRUE #same coverage
  )  
  ptot <- make_white_plot(
    ggplot(df, aes(fc1, fc2, label=grna, color=genedesc)) + 
      geom_point() + 
      xlab(paste0(pool1,": ", cond))+
      ylab(paste0(pool2,": ", cond)))
  ggsave(plot=ptot, file.path(outdir,paste0(pool1," VS ",pool2,".pdf")), width = 4, height = 4)
  ptot
}


outdir <- "/corgi/otherdataset/ellenbushell/newplots_barseq/comparison_rounds"
for(pool1 in list_round1){
  for(pool2 in list_round2){
    print(paste(pool1,pool2))
    compare_two_delta(pool1,pool2, "P BL6 / P RAG1KO")
  }
}
for(pool2 in list_round2){
  for(pool3 in list_round3){
    print(paste(pool2,pool3))
    compare_two_delta(pool2,pool3, "P BL6 / P RAG1KO")  #not present
  }
}



################################################################################
########### Comparison of biological replicates ################################
################################################################################

list_biorep <- c(
  "EB_priming_barseqpool2s_biorep1",
  "EB_priming_barseqpool2s_biorep2",
)

outdir <- "/corgi/otherdataset/ellenbushell/newplots_barseq/comparison_biorep"
compare_two_delta(
  "EB_priming_barseqpool2s_biorep1",
  "EB_priming_barseqpool2s_biorep2", 
  "P BL6 / P RAG1KO")

################################################################################
########### Comparison of PCR variability ######################################
################################################################################

list_pcr1 <- c(
  "EB_priming_barseqpool2s_biorep1_PCR1",
  "EB_priming_barseqpool2s_biorep1_PCR2"
)

list_pcr2 <- c(
  "EB_priming_barseqpool2s_biorep2_PCR1",
  "EB_priming_barseqpool2s_biorep2_PCR2"

)

outdir <- "/corgi/otherdataset/ellenbushell/newplots_barseq/comparison_pcr"
compare_two_delta(
  "EB_priming_barseqpool2s_biorep1_PCR1",
  "EB_priming_barseqpool2s_biorep1_PCR2", 
  "P BL6 / P RAG1KO")

compare_two_delta(
  "EB_priming_barseqpool2s_biorep2_PCR1",
  "EB_priming_barseqpool2s_biorep2_PCR2", 
  "P BL6 / P RAG1KO")


################################################################################
########### Plot of UMAPs ######################################################
################################################################################

all_samplemeta <- readRDS("samplemeta.rds")
all_grstats <- readRDS("grstats.rds")
all_timecourses <- readRDS("timecourses.rds")
all_coverage_stat <- readRDS("coverage_stat.rds")


outdir <- "/corgi/otherdataset/ellenbushell/newplots_barseq/umaps"


library(stringr)
library(ggplot2)
library(patchwork)
library(sqldf)
library(reshape2)
library(cowplot)

list_pools <- c(list_round1, list_round2)   #, list_round3)  #minipool is not happy
for(current_pool in list_pools){
  print(current_pool)
  samplemeta <- all_samplemeta[[current_pool]]
  #coverage_stat <- all_coverage_stat[[current_pool]]
  #print(samplemeta)
  
  samplemeta$day <- sprintf("d%s", samplemeta$day)
  samplemeta$day[samplemeta$day=="dNA"] <- "0"
  p1 <- make_white_plot(ggplot(samplemeta, aes(umap1,umap2,color=mouse_ref))+geom_point())
  p2 <- make_white_plot(ggplot(samplemeta, aes(umap1,umap2,color=day))+geom_point())
  p3 <- make_white_plot(ggplot(samplemeta, aes(umap1,umap2,color=is_input))+geom_point())
  p4 <- make_white_plot(ggplot(samplemeta, aes(umap1,umap2,color=genotype))+geom_point())
  p5 <- make_white_plot(ggplot(samplemeta, aes(umap1,umap2,color=primed))+geom_point())
  p6 <- make_white_plot(ggplot(samplemeta, aes(umap1,umap2,color=total_count))+geom_point())
  ptot <- p1/p2|p3/p4|p5/p6
 
  #UMAP of mice, color by background; symbol by treatment
  p7 <- make_white_plot(ggplot(samplemeta, aes(umap1,umap2,color=genotype, shape=primed))+geom_point())

  p8 <- make_white_plot(ggplot(samplemeta, aes(umap1,umap2,label=samplename))+geom_point()+geom_text())
  
  #TODO: why the odd separation in PCR? need plot with names

  ggsave(plot = ptot, file.path(outdir, paste0("many ",current_pool,".pdf")),               width = 12,  height = 6)
  ggsave(plot = p7,   file.path(outdir, paste0("genotype_treatment ",current_pool,".pdf")), width = 3,  height = 3)
  ggsave(plot = p8,   file.path(outdir, paste0("name ",current_pool,".pdf")),               width = 20, height = 20)
  
  
  # entropy: sum p log(p)
  
}




