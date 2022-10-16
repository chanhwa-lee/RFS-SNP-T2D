#*******************************************************************************
#**********           < RFS-SNP interaction on T2D survival >         **********
#**********           KEGG Pathway analysis from Cox survival         **********
#**********           analysis with RFS-SNP interaction               **********
#**********                                                           **********	
#**********           Written by:				                              **********
#**********           Chanhwa Lee - chanhwa@email.unc.edu             **********
#**********                                                           **********
#**********           Version: 1.0                                    **********
#**********           Created: Aug 10, 2019                           **********
#**********           Modified: Oct 16, 2020                          **********
#*******************************************************************************

sum(coxresultRFS$PVALUE_INTER.TERM < 10^-6) # 0
sum(coxresultRFS$PVALUE_INTER.TERM < 10^-5) # 4
sum(coxresultRFS$PVALUE_INTER.TERM < 10^-4) # 49
sum(coxresultRFS$PVALUE_INTER.TERM < 10^-3) # 489
sum(coxresultRFS$PVALUE_INTER.TERM < 10^-2) # 4334
sum(coxresultRFS$PVALUE_INTER.TERM < 5*10^-2) # 18705
sum(coxresultRFS$PVALUE_INTER.TERM < 10^-1) # 35388

write.csv(AS3.final, "AS3.final.csv")

# Pathway Database from KEGG
kegg <- read.csv("c2.cp.kegg.v7.0.symbols.gmt",header = F)
str(kegg)
path.df <- data.frame()
for(i in 1:nrow(kegg)){
  str <- unlist(strsplit(as.character(kegg[i,]), split = "\t"))[-2]
  geneIDs <- str[-1]
  pathway_name <- rep(str[1],length(geneIDs))
  path.df <- rbind(path.df, data.frame(pathway_name,geneIDs))
}
path.df

# Gene ID from NCBI37
genematch <- read.csv("NCBI37.3.gene.loc",header = F)
str(genematch)
genematch.df <- data.frame(geneID = c(),genenum = c())
for(i in 1:nrow(genematch)){
  str <- unlist(strsplit(as.character(genematch[i,]), split = "\t"))
  genematch.df <- rbind(genematch.df, data.frame(str[6], str[1]))
}
genematch.df

colnames(path.df)
colnames(genematch.df) <- c("geneIDs","genenum")
colnames(genematch.df)
pathnum <- inner_join(path.df,genematch.df, by = "geneIDs")[-2]
write.csv(pathnum, "pathgene.csv")
getwd()

# p-value < 0.01 significant
RSID_0.01 <- coxresultRFS[which(coxresultRFS$PVALUE_INTER.TERM < 10^-2),c(1,8)] # 4334
write.csv(RSID_0.01, "RSID_0.01.csv")
# next, do snp-level and gene-level analysis using tools such as MAGMA, 
# finally obtain gene-set analysis 
# for this, gene ID matching is needed.
# that is why pathgene.csv file needed.

# p-value < 0.05 significant
RSID_0.05 <- coxresultRFS[which(coxresultRFS$PVALUE_INTER.TERM < 5*10^-2),c(1,8)] # 4334
write.csv(RSID_0.05, "RSID_0.05.csv")

# p-value < 0.001 significant
RSID_0.001 <- coxresultRFS[which(coxresultRFS$PVALUE_INTER.TERM < 10^-3),c(1,8)] # 489
write.csv(RSID_0.001, "RSID_0.001.csv")

# all SNPs and p-value
write.csv(coxresultRFS[,c(1,8)], "SNPpval.csv")

# SNP probe ID 2 RSID
a <- coxresultRFS[,c(1,8)]
colnames(a) <- c("ID","p")
b <- read.csv("pathway/gsasnp2/gsasnp2-windows-gui/rsid_Affyid.csv")
colnames(b) <- c("ID","n","rs")
c <- inner_join(a,b)
c[,1] <- c[,4]
c <- c[,-c(3,4)]
c
write.csv(c, "pathway/gsasnp2/gsasnp2-windows-gui/pval all/data/pval.csv")
