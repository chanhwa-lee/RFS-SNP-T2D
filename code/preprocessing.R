#*******************************************************************************
#**********           < RFS-SNP interaction on T2D survival >         **********
#**********           Data preprocessing                              **********
#**********                                                           **********	
#**********           Written by:				                              **********
#**********           Chanhwa Lee - chanhwa@email.unc.edu             **********
#**********                                                           **********
#**********           Version: 1.0                                    **********
#**********           Created: July 11, 2019                          **********
#**********           Modified: Oct 16, 2020                          **********
#*******************************************************************************

Sys.setlocale("LC_ALL", "C")
path.mac <- getwd()
geno <- read.table("160719_3_8840_Affy5.fam", header = T,stringsAsFactors = F,sep = " ")
setwd(paste(path.mac,"CSV",sep="/"))

# Stage 3

AS3_Drug <- read.csv('AS3_DRUG.csv',header = T,stringsAsFactors = F,sep = ',')

AS3_Pastd <- read.csv('AS3_PASTD.csv',header = T,stringsAsFactors = F,sep = ',')
AS3_Pastd <- AS3_Pastd[match(AS3_Drug$NIHID,AS3_Pastd$NIHID),]

AS3_Treat <- read.csv('AS3_Treat.csv',header = T,stringsAsFactors = F,sep = ',')
AS3_Treat <- AS3_Treat[match(AS3_Drug$NIHID,AS3_Treat$NIHID),]

AS3_Biochem <- read.csv('AS3_BIOCHEM.csv',header = T,stringsAsFactors = F,sep = ',')
AS3_Biochem <- AS3_Biochem[match(AS3_Drug$NIHID,AS3_Biochem$NIHID),]

AS3_Anthro <- read.csv('AS3_Anthro.csv',header = T,stringsAsFactors = F,sep = ',')
AS3_Anthro <- AS3_Anthro[match(AS3_Drug$NIHID,AS3_Anthro$X...NIHID),]

AS3_BODYCOM <- read.csv('AS3_BODYCOM.csv',header = T,stringsAsFactors = F,sep = ',')
AS3_BODYCOM <- AS3_BODYCOM[match(AS3_Drug$NIHID,AS3_BODYCOM$X...NIHID),]

AS3_GEN <- read.csv('AS3_GEN.csv',header = T,stringsAsFactors = F,sep = ',')
AS3_GEN <- AS3_GEN[match(AS3_Drug$NIHID,AS3_GEN$NIHID),]

AS3_Drsm <- read.csv("AS3_Drsm.csv",h=T,stringsAsFactors = F,sep=',')
AS3_Drsm <- AS3_Drsm[match(AS3_Drug$NIHID,AS3_Drsm$NIHID),]

AS3_FFQ <- read.csv("AS3_19_FFQ.csv",h=T,stringsAsFactors = F,sep=',')
AS3_FFQ <- AS3_FFQ[match(AS3_Drug$NIHID,AS3_FFQ$DIST_ID),]

AS3_ACTIVE <- read.csv("AS3_ACTIVE.csv",h=T,stringsAsFactors = F,sep=',')
AS3_ACTIVE <- AS3_ACTIVE[match(AS3_Drug$NIHID,AS3_ACTIVE$NIHID),]


AS3 <- data.frame(NIHID=AS3_Drug$NIHID,
                  SEX=AS3_GEN$AS3_SEX,
                  AGE=AS3_GEN$AS3_AGE,
                  AREA=AS3_GEN$AS3_AREA,
                  BMI=AS3_BODYCOM$AS3_BMI,
                  BMI.cal=AS3_Anthro$AS3_WEIGHT/(AS3_Anthro$AS3_HEIGHT/100)^2,
                  
                  T2D=AS3_Pastd$AS3_PDFDM,
                  T2D.age=AS3_Pastd$AS3_PDFDMAG,
                  Treat=AS3_Treat$AS3_TRTDM,
                  Drugins.1=AS3_Drug$AS3_DRUGINS,
                  Drugdm.1=AS3_Drug$AS3_DRUGDM,
                  Druginscu.1=AS3_Drug$AS3_DRUGINSCU,
                  Drugdmcu.1=AS3_Drug$AS3_DRUGDMCU,
                  Drugins.2=AS3_Drug$AS3_DRUGFINS,
                  Drugdm.2=AS3_Drug$AS3_DRUGFDM,
                  Druginscu.2=AS3_Drug$AS3_DRUGFINSCU,
                  Drugdmcu.2=AS3_Drug$AS3_DRUGFDMCU,
                  Glu0=AS3_Biochem$AS3_GLU0,
                  Glu120=AS3_Biochem$AS3_GLU120,
                  
                  Edu=AS3_GEN$AS3_EDU,
                  income=AS3_GEN$AS3_INCOME,
                  Smoke=AS3_Drsm$AS3_SMOKE,
                  Drink=AS3_Drsm$AS3_DRINK,
                  HT=AS3_Pastd$AS3_PDFHT,
                  Coffee = AS3_FFQ$AS3_F089_FQ,
                  
                  Redmeat = AS3_FFQ$AS3_F062_FQ + AS3_FFQ$AS3_F063_FQ +
                    AS3_FFQ$AS3_F057_FQ + AS3_FFQ$AS3_F058_FQ + AS3_FFQ$AS3_F059_FQ,
                  
                  cancer1=AS3_Pastd$AS3_PDFLCA,cancer2=AS3_Pastd$AS3_PDFGCA,
                  cancer3=AS3_Pastd$AS3_PDFHCA,cancer4=AS3_Pastd$AS3_PDFCOLCA,
                  cancer5=AS3_Pastd$AS3_PDFPACA,cancer6=AS3_Pastd$AS3_PDFUTCA,
                  cancer7=AS3_Pastd$AS3_PDFBRCA,cancer8=AS3_Pastd$AS3_PDFCA1,
                  cancer9=AS3_Pastd$AS3_PDFCA2,
                  
                  cardio1=AS3_Pastd$AS3_PDFMI, cardio2=AS3_Pastd$AS3_PDFLP,
                  cardio3=AS3_Pastd$AS3_PDFCH, cardio4=AS3_Pastd$AS3_PDFCD,
                  METs = AS3_ACTIVE$AS3_ACT24_1 * 1.0 + AS3_ACTIVE$AS3_ACT24_5 * 1.5 + 
                    AS3_ACTIVE$AS3_ACT24_4 * 2.4 + AS3_ACTIVE$AS3_ACT24_3 * 5.0 + 
                    AS3_ACTIVE$AS3_ACT24_2 * 7.5,
                  # (1.0 for sedentary, 1.5 for very light, 2.4 for light, 5.0 for moderate, and 7.5 for intense)
                  Insdate = AS3_Drug$AS3_DRUGFINSYR1,
                  Dmdate = AS3_Drug$AS3_DRUGDMYR1,
                  Edate = AS3_GEN$AS3_EDATE3) 

str(AS3)
AS3 <- AS3[-which(AS3$METs > 1000),]
hist(AS3$METs)
AS3 <- na.omit(AS3)
dim(AS3) # 7078 

# Stage 4

AS4_01_EXAMINEE <- read.csv('AS4_01_EXAMINEE.csv',header = T,stringsAsFactors = F,sep = ',')

AS4_04_MEDIC <- read.csv('AS4_04_MEDIC.csv',header = T,stringsAsFactors = F,sep = ',')
AS4_04_MEDIC <- AS4_04_MEDIC[match(AS4_01_EXAMINEE$NIHID,AS4_04_MEDIC$NIHID),]

AS4_07_DISEASE <- read.csv('AS4_07_DISEASE.csv',header = T,stringsAsFactors = F,sep = ',')
AS4_07_DISEASE <- AS4_07_DISEASE[match(AS4_01_EXAMINEE$NIHID,AS4_07_DISEASE$NIHID),]

AS4_08_TREAT <- read.csv('AS4_08_TREAT.csv',header = T,stringsAsFactors = F,sep = ',')
AS4_08_TREAT <- AS4_08_TREAT[match(AS4_01_EXAMINEE$NIHID,AS4_08_TREAT$X...NIHID),]

AS4_09_DRUG <- read.csv('AS4_09_DRUG.csv',header = T,stringsAsFactors = F,sep = ',')
AS4_09_DRUG <- AS4_09_DRUG[match(AS4_01_EXAMINEE$NIHID,AS4_09_DRUG$NIHID),]

AS4_20_BIOCHEM <- read.csv('AS4_20_BIOCHEM.csv',header = T,stringsAsFactors = F,sep = ',')
AS4_20_BIOCHEM <- AS4_20_BIOCHEM[match(AS4_01_EXAMINEE$NIHID,AS4_20_BIOCHEM$NIHID),]

AS4_27_RESULT1 <- read.csv('AS4_27_RESULT1.csv',header = T,stringsAsFactors = F,sep = ',')
AS4_27_RESULT1 <- AS4_27_RESULT1[match(AS4_01_EXAMINEE$NIHID,AS4_27_RESULT1$X...NIHID),]

AS4_30_RESULT4_INBODY <- read.csv('AS4_30_RESULT4_INBODY.csv',header = T,stringsAsFactors = F,sep = ',')
AS4_30_RESULT4_INBODY <- AS4_30_RESULT4_INBODY[match(AS4_01_EXAMINEE$NIHID,AS4_30_RESULT4_INBODY$X...NIHID),]

AS4_GEN <- read.csv('AS4_GEN.csv',header = T,stringsAsFactors = F,sep = ',')
AS4_GEN <- AS4_GEN[match(AS4_01_EXAMINEE$NIHID,AS4_GEN$NIHID),]

AS4 <- data.frame(NIHID=AS4_09_DRUG$NIHID,SEX=AS4_01_EXAMINEE$AS4_SEX,AGE=AS4_01_EXAMINEE$AS4_AGE,AREA=AS4_01_EXAMINEE$AS4_DATA_CLASS,
                  BMI=AS4_30_RESULT4_INBODY$AS4_BMI,BMI.cal=AS4_27_RESULT1$AS4_WEIGHT/(AS4_27_RESULT1$AS4_HEIGHT/100)^2,
                  T2D=AS4_07_DISEASE$AS4_DM,T2D.age=AS4_07_DISEASE$AS4_DMAG,
                  Treat=AS4_08_TREAT$AS4_TREATD2,
                  Drugins.1=AS4_09_DRUG$AS4_D54A,Drugdm.1=AS4_09_DRUG$AS4_D62A,
                  Druginscu.1=AS4_09_DRUG$AS4_D54D,Drugdmcu.1=AS4_09_DRUG$AS4_D62D,
                  Glu0=AS4_20_BIOCHEM$AS4_GLU0,Glu120=AS4_20_BIOCHEM$AS4_GLU120,
                  Insdate = AS4_09_DRUG$AS4_D54C,
                  Dmdate = AS4_09_DRUG$AS4_D62C,
                  T2Dyr = AS4_04_MEDIC$AS4_DMDIAGYR,
                  T2Dmo = AS4_04_MEDIC$AS4_DMDIAGMO,
                  Edate = AS4_01_EXAMINEE$AS4_EDATE4)

dim(AS4)

# Stage 5

AS5_01_EXAMINEE <- read.csv('AS5_01_EXAMINEE.csv',header = T,stringsAsFactors = F,sep = ',')

AS5_04_MEDIC <- read.csv('AS5_04_MEDIC.csv',header = T,stringsAsFactors = F,sep = ',')
AS5_04_MEDIC <- AS5_04_MEDIC[match(AS5_01_EXAMINEE$NIHID,AS5_04_MEDIC$X...NIHID),]

AS5_06_DISEASE <- read.csv('AS5_06_DISEASE.csv',header = T,stringsAsFactors = F,sep = ',')
AS5_06_DISEASE <- AS5_06_DISEASE[match(AS5_01_EXAMINEE$NIHID,AS5_06_DISEASE$NIHID),]

AS5_05_Treat <- read.csv('AS5_05_Treat.csv',header = T,stringsAsFactors = F,sep = ',')
AS5_05_Treat <- AS5_05_Treat[match(AS5_01_EXAMINEE$NIHID,AS5_05_Treat$X...NIHID),]

AS5_08_DRUG <- read.csv('AS5_08_DRUG.csv',header = T,stringsAsFactors = F,sep = ',')
AS5_08_DRUG <- AS5_08_DRUG[match(AS5_01_EXAMINEE$NIHID,AS5_08_DRUG$NIHID),]

AS5_36_BIOCHEM <- read.csv('AS5_36_BIOCHEM.csv',header = T,stringsAsFactors = F,sep = ',')
AS5_36_BIOCHEM <- AS5_36_BIOCHEM[match(AS5_01_EXAMINEE$NIHID,AS5_36_BIOCHEM$NIHID),]

AS5_27_RESULT1 <- read.csv('AS5_27_RESULT1.csv',header = T,stringsAsFactors = F,sep = ',')
AS5_27_RESULT1 <- AS5_27_RESULT1[match(AS5_01_EXAMINEE$NIHID,AS5_27_RESULT1$X...NIHID),]

AS5_30_RESULT4_INBODY <- read.csv('AS5_30_RESULT4_INBODY.csv',header = T,stringsAsFactors = F,sep = ',')
AS5_30_RESULT4_INBODY <- AS5_30_RESULT4_INBODY[match(AS5_01_EXAMINEE$NIHID,AS5_30_RESULT4_INBODY$X...NIHID),]

AS5_07_Hospitalization <- read.csv('AS5_07_Hospitalization.csv',header = T,stringsAsFactors = F,sep = ',')
AS5_07_Hospitalization <- AS5_07_Hospitalization[match(AS5_01_EXAMINEE$NIHID,AS5_07_Hospitalization$X...NIHID),]

AS5 <- data.frame(NIHID=AS5_08_DRUG$NIHID,SEX=AS5_01_EXAMINEE$AS5_SEX,AGE=AS5_01_EXAMINEE$AS5_AGE,AREA=AS5_01_EXAMINEE$AS5_DATA_CLASS,
                  BMI=AS5_30_RESULT4_INBODY$AS5_BMI,BMI.cal=AS5_27_RESULT1$AS5_WEIGHT/(AS5_27_RESULT1$AS5_HEIGHT/100)^2,
                  T2D=AS5_06_DISEASE$AS5_DM,T2D.age=AS5_06_DISEASE$AS5_DMAG,
                  Treat=AS5_05_Treat$AS5_TREATD2,
                  Drugins.1=AS5_08_DRUG$AS5_D54A,Drugdm.1=AS5_08_DRUG$AS5_D62A,
                  Druginscu.1=AS5_08_DRUG$AS5_D54D,Drugdmcu.1=AS5_08_DRUG$AS5_D62D,
                  Glu0=AS5_36_BIOCHEM$AS5_GLU0,Glu120=AS5_36_BIOCHEM$AS5_GLU120,
                  Insdate = AS5_08_DRUG$AS5_D54C,
                  Dmdate = AS5_08_DRUG$AS5_D62C,
                  T2Dyr = AS5_04_MEDIC$AS5_DMDIAGYR,
                  T2Dmo = AS5_04_MEDIC$AS5_DMDIAGMO,
                  Edate = AS5_01_EXAMINEE$AS5_EDATE5)
dim(AS5)


# Stage 6

AS6_01_EXAMINEE <- read.csv('AS6_01_EXAMINEE.csv',header = T,stringsAsFactors = F,sep = ',')

AS6_06_MEDIC <- read.csv('AS6_06_MEDIC.csv',header = T,stringsAsFactors = F,sep = ',')
AS6_06_MEDIC <- AS6_06_MEDIC[match(AS6_01_EXAMINEE$NIHID,AS6_06_MEDIC$X...NIHID),]

AS6_08_DISEASE <- read.csv('AS6_08_DISEASE.csv',header = T,stringsAsFactors = F,sep = ',')
AS6_08_DISEASE <- AS6_08_DISEASE[match(AS6_01_EXAMINEE$NIHID,AS6_08_DISEASE$NIHID),]

AS6_07_TREAT <- read.csv('AS6_07_TREAT.csv',header = T,stringsAsFactors = F,sep = ',')
AS6_07_TREAT <- AS6_07_TREAT[match(AS6_01_EXAMINEE$NIHID,AS6_07_TREAT$X...NIHID),]

AS6_10_DRUG <- read.csv('AS6_10_DRUG.csv',header = T,stringsAsFactors = F,sep = ',')
AS6_10_DRUG <- AS6_10_DRUG[match(AS6_01_EXAMINEE$NIHID,AS6_10_DRUG$NIHID),]

AS6_37_BIOCHEM <- read.csv('AS6_37_BIOCHEM.csv',header = T,stringsAsFactors = F,sep = ',')
AS6_37_BIOCHEM <- AS6_37_BIOCHEM[match(AS6_01_EXAMINEE$NIHID,AS6_37_BIOCHEM$NIHID),]

AS6_29_RESULT1_BODY <- read.csv('AS6_29_RESULT1_BODY.csv',header = T,stringsAsFactors = F,sep = ',')
AS6_29_RESULT1_BODY <- AS6_29_RESULT1_BODY[match(AS6_01_EXAMINEE$NIHID,AS6_29_RESULT1_BODY$X...NIHID),]

AS6_32_RESULT4_INBODY <- read.csv('AS6_32_RESULT4_INBODY.csv',header = T,stringsAsFactors = F,sep = ',')
AS6_32_RESULT4_INBODY <- AS6_32_RESULT4_INBODY[match(AS6_01_EXAMINEE$NIHID,AS6_32_RESULT4_INBODY$X...NIHID),]

AS6_09_Hospitalization <- read.csv('AS6_09_HOSPITALIZATION.csv',header = T,stringsAsFactors = F,sep = ',')
AS6_09_Hospitalization <- AS6_09_Hospitalization[match(AS6_01_EXAMINEE$NIHID,AS6_09_Hospitalization$X...NIHID),]

AS6 <- data.frame(NIHID=AS6_10_DRUG$NIHID,SEX=AS6_01_EXAMINEE$AS6_SEX,AGE=AS6_01_EXAMINEE$AS6_AGE,AREA=AS6_01_EXAMINEE$AS6_DATA_CLASS,
                  BMI=AS6_32_RESULT4_INBODY$AS6_BMI,BMI.cal=AS6_29_RESULT1_BODY$AS6_WEIGHT/(AS6_29_RESULT1_BODY$AS6_HEIGHT/100)^2,
                  T2D=AS6_08_DISEASE$AS6_DM,T2D.age=AS6_08_DISEASE$AS6_DMAG,
                  Treat=AS6_07_TREAT$AS6_TREATD2,
                  Drugins.1=AS6_10_DRUG$AS6_D54A,Drugdm.1=AS6_10_DRUG$AS6_D62A,
                  Druginscu.1=AS6_10_DRUG$AS6_D54D,Drugdmcu.1=AS6_10_DRUG$AS6_D62D,
                  Glu0=AS6_37_BIOCHEM$AS6_GLU0_TR,Glu120=AS6_37_BIOCHEM$AS6_GLU120_TR,
                  Edate = AS6_01_EXAMINEE$AS6_EDATE6)

############################## 
#        Data filtering      #
##############################

# start with 2,580 samples

# STEP1 : no inforamtion of adjusting factors filtered out
AS3 <- AS3[which(AS3$income!=99999 &
                   AS3$Drink!=99999 &
                   AS3$Smoke!=99999 &
                   AS3$Coffee!=99999 &
                   AS3$Edu!=99999 &
                   AS3$METs<1000 &
                   AS3$BMI.cal > 10 &
                   AS3$Redmeat<99999),]

AS3$SEX <- as.factor(AS3$SEX)
AS3$AREA <- as.factor(AS3$AREA)

AS3$Smoke <- ifelse(AS3$Smoke==3,2,1)
AS3$Drink <- ifelse(AS3$Drink==3,2,1)
AS3$Smoke <- as.factor(AS3$Smoke)
AS3$Drink <- as.factor(AS3$Drink)

AS3 <- na.omit(AS3)
dim(AS3)

# STEP2 : no inforamtion or having T2D filtered out
noinfoind <- which(AS3$Drugins.1 >=66666 & AS3$Drugdm.1 >=66666 
                   & AS3$Glu0==99999 &  AS3$Glu120==99999) 
length(noinfoind)# length 12
AS3.info <- AS3[-noinfoind,]
dim(AS3.info) # 2568 by 40

# Exisiting Diabetes
T2Dind <- which(AS3.info$T2D==2 | AS3.info$Drugins.1==2 | AS3.info$Drugdm.1==2 | 
                  (AS3.info$Glu0 >= 126 & AS3.info$Glu0 < 10000) |
                  (AS3.info$Glu120 >= 200 & AS3.info$Glu120 < 10000) ) 
length(T2Dind) # 528
AS3.noT2D <- AS3.info[-T2Dind,]
dim(AS3.noT2D) # 2007 by 40


# STEP3 : having Cancer filtered out
cancerind <- unique(which(AS3.noT2D$cancer1==2 | AS3.noT2D$cancer2==2 | AS3.noT2D$cancer3==2 |
                            AS3.noT2D$cancer4==2 | AS3.noT2D$cancer5==2 | AS3.noT2D$cancer6==2 |
                            AS3.noT2D$cancer7==2 | AS3.noT2D$cancer8==2 | AS3.noT2D$cancer9==2))

length(cancerind) # 82
AS3.nocan <- AS3.noT2D[-cancerind,]
dim(AS3.nocan) # 4,244 by 34

# STEP4 : having cardiovascular diseases filtered out
CDind <- which(AS3.nocan$cardio1==2 | AS3.nocan$cardio2==2 | 
                 AS3.nocan$cardio3==2 | AS3.nocan$cardio4==2 )
length(CDind) # 70

AS3.nocardio <- AS3.nocan[-CDind,]
dim(AS3.nocardio) # 4,174 by 34

colnames(AS3.nocardio)
AS3.nocardio <- AS3.nocardio[,-match( c("cancer1", "cancer2", "cancer3", "cancer4","cancer5", 
                                        "cancer6", "cancer7", "cancer8", "cancer9", 
                                        "cardio1", "cardio2", "cardio3", "cardio4"), colnames(AS3.nocardio))]
dim(AS3.nocardio) # 4,174 by 21
