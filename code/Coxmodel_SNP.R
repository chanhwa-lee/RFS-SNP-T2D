#*******************************************************************************
#**********           < RFS-SNP interaction on T2D survival >         **********
#**********           Cox Survival analysis (with SNP)                **********
#**********                                                           **********	
#**********           Written by:				                              **********
#**********           Chanhwa Lee - chanhwa@email.unc.edu             **********
#**********                                                           **********
#**********           Version: 1.0                                    **********
#**********           Created: Aug 5, 2019                            **********
#**********           Modified: Oct 16, 2020                          **********
#*******************************************************************************

install.packages("remotes")
remotes::install_github("suchestoncampbelllab/gwasurvivr")
library("gwasurvivr")

###########
datacox.o <- AS3.total
colnames(datacox.o)
datacox.o$SEX <- as.numeric(datacox.o$SEX)
datacox.o$AREA <- as.numeric(datacox.o$AREA)
datacox.o <- datacox.o[,c(1,2,3,4,5,32,33,35)]

plinkCoxSurv(bed.file="SNPdata/indfiltered.bed",
             covariate.file=datacox.o,
             id.column="NIHID",
             time.to.event="time",
             event="status",
             covariates=c("SEX", "AGE", "AREA", "BMI", "RFS.o"),
             inter.term="RFS.o",
             print.covs="some",
             out.file="CoxSNP_RFS_o",
             chunk.size=50,
             maf.filter=F,
             flip.dosage=TRUE,
             verbose=TRUE,
             clusterObj=NULL)  

coxresultRFS <- read.table("SNPdata/CoxSNP_RFS_o.coxph", header = T)

################################
# RFS and SNP interaction only #
################################

colnames(coxresultRFS)
result_final_inter<-cbind(as.vector(coxresultRFS[,c(2,1,8)]))
colnames(result_final_inter) <- c("CHR","SNP", "P")

result_final_inter2 <- cbind(as.vector(coxresultRFS[,c(2,1,11,8)]))
colnames(result_final_inter2) <- c("CHR","SNP", "BETA","P")

write.table(result_final_inter, "RFS_SNP_interaction_result_plot.txt", quote = F, row.names = F)
write.table(result_final_inter2, "RFS_SNP_interaction_result.txt", quote = F, row.names = F)
system("Rscript qqplot.R RFS_SNP_interaction_result_plot.txt RFS_SNP_interaction_result")
system("Rscript manhattan.R RFS_SNP_interaction_result_plot.txt RFS_SNP_interaction_result")

#########################
# SNP significance only #
#########################

result_final_SNP<-cbind(as.vector(coxresultRFS[,c(2,1,9)]))
colnames(result_final_SNP) <- c("CHR","SNP", "P")

result_final_SNP2 <- cbind(as.vector(coxresultRFS[,c(2,1,12,9)]))
colnames(result_final_SNP2) <- c("CHR","SNP", "BETA","P")

write.table(result_final_SNP, "RFS_SNP_result_plot.txt", quote = F, row.names = F)
write.table(result_final_SNP2, "RFS_SNP_result.txt", quote = F, row.names = F)
system("Rscript qqplot.R RFS_SNP_result_plot.txt RFS_SNP_result")
system("Rscript manhattan.R RFS_SNP_result_plot.txt RFS_SNP_result")
