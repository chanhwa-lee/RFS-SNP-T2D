#*******************************************************************************
#**********           < RFS-SNP interaction on T2D survival >         **********
#**********           Cox Survival analysis (without SNP)             **********
#**********                                                           **********	
#**********           Written by:				                              **********
#**********           Chanhwa Lee - chanhwa@email.unc.edu             **********
#**********                                                           **********
#**********           Version: 1.0                                    **********
#**********           Created: July 20, 2019                          **********
#**********           Modified: Oct 16, 2020                          **********
#*******************************************************************************

library(dplyr)
library(survival)

# T2D diagnosis
# First, Modified RFS 
AS3.final <- AS3.nocardio ; dim(AS3.final)
AS4.final <- AS4[which(!is.na(match(AS4$NIHID,AS3.final$NIHID))),] ; dim(AS4.final)
AS5.final <- AS5[which(!is.na(match(AS5$NIHID,AS4.final$NIHID))),] ; dim(AS5.final)
AS6.final <- AS6[which(!is.na(match(AS6$NIHID,AS5.final$NIHID))),] ; dim(AS6.final)

# Transform to Date type
# Stage 3

colnames(AS3.final)

for(j in 1:length(AS3.final$NIHID)){
  if(AS3.final$Insdate[j] < 100000){
    AS3.final$Insdate[j] <- NA
  }
  if(AS3.final$Dmdate[j] < 100000){
    AS3.final$Dmdate[j] <- NA
  }
}

AS3.final$Insdate <- as.Date(as.character(AS3.final$Insdate), "%Y%m%d")  
AS3.final$Dmdate <- as.Date(as.character(AS3.final$Dmdate), "%Y%m%d")  
AS3.final$Diagdate <- pmin(AS3.final$Insdate, AS3.final$Dmdate, na.rm = T)

AS3.final$Edate <- as.Date(paste0(as.character(AS3.final$Edate),"01"), "%Y%m%d")  

str(AS3.final)

# Stage 4
colnames(AS4.final)
# T2D diagnostics date
AS4.final$T2Ddate <- 0
for(j in 1:length(AS4.final$NIHID)){
  yr <- AS4.final$T2Dyr[j]
  mo <- AS4.final$T2Dmo[j]
  if(yr < 10000 & mo < 10){
    AS4.final$T2Ddate[j] <- paste0(yr,"0",mo)
  }
  if(yr < 10000 & mo < 10000 & mo >= 10){
    AS4.final$T2Ddate[j] <- paste0(yr,mo)
  }
  if(yr < 10000 & mo >= 10000){
    AS4.final$T2Ddate[j] <- paste0(yr,"06")
  }
  if(yr >= 10000){
    AS4.final$T2Ddate[j] <- NA
  }
}
# Insuline & Drug onset date
colnames(AS4.final)

for(j in 1:length(AS4.final$NIHID)){
  if(AS4.final$Insdate[j] < 100000){
    AS4.final$Insdate[j] <- NA
  }
  if(AS4.final$Dmdate[j] < 100000){
    AS4.final$Dmdate[j] <- NA
  }
}

AS4.final$Insdate <- as.Date(as.character(AS4.final$Insdate), "%Y%m%d")  
AS4.final$Dmdate <- as.Date(as.character(AS4.final$Dmdate), "%Y%m%d")  
AS4.final$T2Ddate <- as.Date(paste0(as.character(AS4.final$T2Ddate),"01"), "%Y%m%d")
AS4.final$Diagdate <- pmin(AS4.final$Insdate, AS4.final$Dmdate, AS4.final$T2Ddate, na.rm = T)

AS4.final$Edate <- as.Date(paste0(as.character(AS4.final$Edate),"01"), "%Y%m%d")
AS4.final$priEdate <- AS3.final$Edate[match(AS4.final$NIHID, AS3.final$NIHID)] 
AS4.final$baseEdate <- AS3.final$Edate[match(AS4.final$NIHID, AS3.final$NIHID)] 

str(AS4.final)

# Error index eliminating
errID.4 <- as.character(AS4.final$NIHID[which(AS4.final$Diagdate < AS4.final$baseEdate)])
AS3.final <- AS3.final[-na.omit(match(errID.4,AS3.final$NIHID)),]; dim(AS3.final)
AS4.final <- AS4.final[-na.omit(match(errID.4,AS4.final$NIHID)),]; dim(AS4.final)
AS5.final <- AS5.final[-na.omit(match(errID.4,AS5.final$NIHID)),]; dim(AS5.final)
AS6.final <- AS6.final[-na.omit(match(errID.4,AS6.final$NIHID)),]; dim(AS6.final)

AS4.final[which(AS4.final$NIHID == "NIH16K7211996"),]

# Stage 5
colnames(AS5.final)
# T2D diagnostics date
AS5.final$T2Ddate <- 0
for(j in 1:length(AS5.final$NIHID)){
  yr <- AS5.final$T2Dyr[j]
  mo <- AS5.final$T2Dmo[j]
  if(yr < 10000 & mo < 10){
    AS5.final$T2Ddate[j] <- paste0(yr,"0",mo)
  }
  if(yr < 10000 & mo < 10000 & mo >= 10){
    AS5.final$T2Ddate[j] <- paste0(yr,mo)
  }
  if(yr < 10000 & mo >= 10000){
    AS5.final$T2Ddate[j] <- paste0(yr,"06")
  }
  if(yr >= 10000){
    AS5.final$T2Ddate[j] <- NA
  }
}
# Insuline & Drug onset date
for(j in 1:length(AS5.final$NIHID)){
  if(AS5.final$Insdate[j] < 100000){
    AS5.final$Insdate[j] <- NA
  }
  if(AS5.final$Dmdate[j] < 100000){
    AS5.final$Dmdate[j] <- NA
  }
}

AS5.final$Insdate <- as.Date(as.character(AS5.final$Insdate), "%Y%m%d")  
AS5.final$Dmdate <- as.Date(as.character(AS5.final$Dmdate), "%Y%m%d")  
AS5.final$T2Ddate <- as.Date(paste0(as.character(AS5.final$T2Ddate),"01"), "%Y%m%d")
AS5.final$Diagdate <- pmin(AS5.final$Insdate, AS5.final$Dmdate, AS5.final$T2Ddate, na.rm = T)

AS5.final$Edate <- as.Date(paste0(as.character(AS5.final$Edate),"01"), "%Y%m%d") 
AS5.final$priEdate <- AS4.final$Edate[match(AS5.final$NIHID, AS4.final$NIHID)] 
AS5.final$baseEdate <- AS3.final$Edate[match(AS5.final$NIHID, AS3.final$NIHID)] 

str(AS5.final)

# Error index eliminating
errID.5 <- as.character(AS5.final$NIHID[which(AS5.final$Diagdate < AS5.final$baseEdate)])
AS3.final <- AS3.final[-na.omit(match(errID.5,AS3.final$NIHID)),]; dim(AS3.final)
AS4.final <- AS4.final[-na.omit(match(errID.5,AS4.final$NIHID)),]; dim(AS4.final)
AS5.final <- AS5.final[-na.omit(match(errID.5,AS5.final$NIHID)),]; dim(AS5.final)
AS6.final <- AS6.final[-na.omit(match(errID.5,AS6.final$NIHID)),]; dim(AS6.final)

# Stage 6
colnames(AS6.final)
AS6.final$Edate <- as.Date(paste0(as.character(AS6.final$Edate),"01"), "%Y%m%d")  
AS6.final$priEdate <- AS5.final$Edate[match(AS6.final$NIHID, AS5.final$NIHID)] 
AS6.final$baseEdate <- AS3.final$Edate[match(AS6.final$NIHID, AS3.final$NIHID)] 

str(AS6.final)

# Error index eliminating - not needed
errID.6 <- as.character(AS6.final$NIHID[which(AS6.final$Diagdate < AS6.final$baseEdate)])
# AS3.final <- AS3.final[-na.omit(match(errID.4,AS3.final$NIHID)),]; dim(AS3.final)
# AS4.final <- AS4.final[-na.omit(match(errID.4,AS4.final$NIHID)),]; dim(AS4.final)
# AS5.final <- AS5.final[-na.omit(match(errID.4,AS5.final$NIHID)),]; dim(AS5.final)
# AS6.final <- AS6.final[-na.omit(match(errID.4,AS6.final$NIHID)),]; dim(AS6.final)

# Cox analysis
AS3.final$time <- 2250
AS3.final$status <- 0

## Stage 4
# loss of follow-up
lossind.4 <- which(is.na(match(AS3.final$NIHID,AS4.final$NIHID)))
for(i in lossind.4){
  AS3.final$time[i] <- 0
}

# Diabetes occurence at stage 4
self.rep.ind.4 <- which(!is.na(AS4.final$Diagdate))
length(self.rep.ind.4) # 3
no.self.rep.ind.4 <- which( AS4.final$T2D==2 
                           | AS4.final$Drugins.1==2 
                           | AS4.final$Drugdm.1==2 
                           | (AS4.final$Glu0 >= 126 & AS4.final$Glu0 < 10000) 
                           | (AS4.final$Glu120 >= 200 & AS4.final$Glu120 < 10000) ) 
length(no.self.rep.ind.4) # 128
no.self.rep.ind.4 <- setdiff(no.self.rep.ind.4,self.rep.ind.4)
length(no.self.rep.ind.4) # 126
T2Dind.4 <- union(self.rep.ind.4,no.self.rep.ind.4)
length(T2Dind.4) # 129

# Self-reported O -> Ins / Dm / Diag date 
for(i in self.rep.ind.4){
  j <- which(AS3.final$NIHID == as.character(AS4.final$NIHID[i]))
  AS3.final$time[j] <- as.integer(AS4.final$Diagdate[i]-AS4.final$baseEdate[i])
  AS3.final$status[j] <- 1
}

# Self-reported X -> time diff btw examine date 
for(i in no.self.rep.ind.4){
  j <- which(AS3.final$NIHID == as.character(AS4.final$NIHID[i]))
  AS3.final$time[j] <- as.integer(AS4.final$Edate[i]-AS4.final$baseEdate[i])
  AS3.final$status[j] <- 1
}
# normal -> not change anything
AS4.norm <- AS4.final[-T2Dind.4,]
dim(AS4.norm) # 1498 by 24
Norstage.4 <- as.character(AS4.final$NIHID[-T2Dind.4])
AS5.final <- AS5.final[which(!is.na(match(AS5.final$NIHID, Norstage.4))),]




## Stage 5

# loss of follow-up
lossind.5 <- which(is.na(match(AS4.norm$NIHID,AS5.final$NIHID)))
for(i in lossind.5){
  j <- which(AS3.final$NIHID == as.character(AS4.norm$NIHID[i]))
  AS3.final$time[j] <- as.integer(AS4.norm$Edate[i]-AS4.norm$baseEdate[i])
}


# Diabetes occurence at stage 5
self.rep.ind.5 <- which(!is.na(AS5.final$Diagdate))
length(self.rep.ind.5) # 27
no.self.rep.ind.5 <- which(AS5.final$T2D==2 | AS5.final$Drugins.1==2 | AS5.final$Drugdm.1==2 | 
                             (AS5.final$Glu0 >= 126 & AS5.final$Glu0 < 10000) |
                             (AS5.final$Glu120 >= 200 & AS5.final$Glu120 < 10000) ) 
length(no.self.rep.ind.5) # 132
no.self.rep.ind.5 <- setdiff(no.self.rep.ind.5,self.rep.ind.5)
length(no.self.rep.ind.5) # 110
T2Dind.5 <- union(self.rep.ind.5,no.self.rep.ind.5)
length(T2Dind.5) # 137

# Self-reported O -> Ins / Dm / Diag date 
for(i in self.rep.ind.5){
  j <- which(AS3.final$NIHID == as.character(AS5.final$NIHID[i]))
  AS3.final$time[j] <- as.integer(AS5.final$Diagdate[i]-AS5.final$baseEdate[i])
  AS3.final$status[j] <- 1
}

# Self-reported X -> time diff btw examine date  
for(i in no.self.rep.ind.5){
  j <- which(AS3.final$NIHID == as.character(AS5.final$NIHID[i]))
  AS3.final$time[j] <- as.integer(AS5.final$Edate[i]-AS5.final$baseEdate[i])
  AS3.final$status[j] <- 1
}

# normal -> not change anything
AS5.norm <- AS5.final[-T2Dind.5,]
dim(AS5.norm) # 2436 by 17
Norstage.5 <- as.character(AS5.final$NIHID[-T2Dind.5])
AS6.final <- AS6.final[which(!is.na(match(AS6.final$NIHID, Norstage.5))),]

## Stage 6

# loss of follow-up
lossind.6 <- which(is.na(match(AS5.norm$NIHID,AS6.final$NIHID)))

for(i in lossind.6){
  j <- which(AS3.final$NIHID == as.character(AS5.norm$NIHID[i]))
  AS3.final$time[j] <- as.integer(AS5.norm$Edate[i]-AS5.final$baseEdate[i])
}

# Diabetes occurence at stage 6
self.rep.ind.6 <- which(!is.na(AS6.final$Diagdate))
length(self.rep.ind.6) # 0
no.self.rep.ind.6 <- which(AS6.final$T2D==2 | AS6.final$Drugins.1==2 | AS6.final$Drugdm.1==2 | 
                             (AS6.final$Glu0 >= 126 & AS6.final$Glu0 < 10000) |
                             (AS6.final$Glu120 >= 200 & AS6.final$Glu120 < 10000) ) 
length(no.self.rep.ind.6) # 155
no.self.rep.ind.6 <- setdiff(no.self.rep.ind.6,self.rep.ind.6)
length(no.self.rep.ind.6) # 155
T2Dind.6 <- union(self.rep.ind.6,no.self.rep.ind.6)
length(T2Dind.6) # 155

# Self-reported O -> Ins / Dm / Diag date 
for(i in self.rep.ind.6){
  j <- which(AS3.final$NIHID == as.character(AS6.final$NIHID[i]))
  AS3.final$time[j] <- as.integer(AS6.final$Diagdate[i]-AS6.final$baseEdate[i])
  AS3.final$status[j] <- 1
}

# Self-reported X -> time diff btw examine date 
for(i in no.self.rep.ind.6){
  j <- which(AS3.final$NIHID == as.character(AS6.final$NIHID[i]))
  AS3.final$time[j] <- as.integer(AS6.final$Edate[i]-AS6.final$baseEdate[i])
  AS3.final$status[j] <- 1
}
# normal -> not change anything
AS6.norm <- AS6.final[-T2Dind.6,]
dim(AS6.norm) # 1996 by 15
Norstage.6 <- as.character(AS6.final$NIHID[-T2Dind.6])

RFSinfo <- read.csv("RFSinfo.csv") ; dim(RFSinfo)
colnames(RFSinfo)[2] <- "NIHID"
AS3.total <- inner_join(AS3.final,RFSinfo); dim(AS3.total)

# str(origin.RFS)
# str(AS3.final)
# 
# sum(!is.na(match(AS3.final$NIHID,modified.RFS$NIHID)))
# sum(!is.na(match(AS3.final$NIHID,origin.RFS$NIHID)))
# 
# AS3.o.RFS.total <- inner_join(AS3.final,origin.RFS)
# AS3.m.RFS.total <- inner_join(AS3.final,modified.RFS)
# 
# colnames(AS3.o.RFS.total)
# colnames(AS3.m.RFS.total)
# 
# data.o <- AS3.o.RFS.total
# data.m <- AS3.m.RFS.total
# 
# table(data.o$SEX)
# table(data.o$AGE)
# range(data.o$BMI.cal)
# range(data.o$income)
# table(data.o$HT)
# table(data.o$Drink)
# table(data.o$Smoke)
# table(data.o$Edu)
# table(data.o$Coffee)
# table(data.o$Redmeat1)
# table(data.o$Redmeat2)
# table(data.o$Redmeat3)
# table(data.o$Redmeat4)
# table(data.o$Redmeat5)
# 
# table(data.m$SEX)
# table(data.m$AGE)
# range(data.m$BMI.cal)
# range(data.m$income)
# table(data.m$HT)
# table(data.m$Drink)
# table(data.m$Smoke)
# table(data.m$Coffee)
# table(data.m$Redmeat1)
# table(data.m$Redmeat2)
# table(data.m$Redmeat3)
# table(data.m$Redmeat4)
# table(data.m$Redmeat5)
# 
# hist(data.m$time)
# table(data.m$time)
# hist(data.m$BMI.cal)
# hist(data.m$RFS)
# hist(data.m$Redmeat)
# hist(data.m$METs)
# 
# plot(data.m$time, data.m$BMI.cal)
# table(data.m$BMI.cal)
# plot(data.m$time, data.m$RFS)
# data.m1 <- data.m
# data.m <- data.m[-which(data.m$time == 2250),]
# data.m <- data.m[-which(data.m$time == 365),]
# plot(data.m$time, data.m$RFS)
# data.m <- data.m1
# plot(data.o$time, data.o$RFS_1st)
#data.o <- data.o[-which(data.o$time == 2250),]
#data.o <- data.o[-which(data.o$time == 365),]

##################################
########### Analysis #############
##################################

colnames(AS3.total)

# Cox with adjusting factors
cox.o <- coxph(Surv(time,status)~AGE+SEX+BMI.cal+AREA+income+Edu+Drink+Smoke+METs+Coffee+Redmeat+RFS.o, data=AS3.total); summary(cox.o)
cox.m <- coxph(Surv(time,status)~AGE+SEX+BMI.cal+AREA+income+Edu+Drink+Smoke+METs+Coffee+Redmeat+RFS.m, data=AS3.total); summary(cox.m)

# Cox without adjusting factors
cox.o.r <- coxph(Surv(time,status)~AGE+SEX+BMI.cal+AREA+RFS.o, data=AS3.total); summary(cox.o.r)
cox.m.r <- coxph(Surv(time,status)~AGE+SEX+BMI.cal+AREA+RFS.m, data=AS3.total); summary(cox.m.r)

# Logistic Regression with adjusting factors
log.o <- glm(status ~ AGE+SEX+BMI.cal+AREA+income+Edu+Drink+Smoke+METs+Coffee+Redmeat+RFS.o, data=AS3.total, family = "binomial"); summary(log.o)
log.m <- glm(status ~ AGE+SEX+BMI.cal+AREA+income+Edu+Drink+Smoke+METs+Coffee+Redmeat+RFS.m, data=AS3.total, family = "binomial"); summary(log.m)

# Logistic Regression without adjusting factors
log.o.r <- glm(status ~ AGE+SEX+BMI.cal+AREA+RFS.o, data=AS3.total, family = "binomial"); summary(log.o.r)
log.m.r <- glm(status ~ AGE+SEX+BMI.cal+AREA+RFS.m, data=AS3.total, family = "binomial"); summary(log.m.r)

plot(AS3.total$time, AS3.total$RFS.o)
plot(AS3.total$time, AS3.total$RFS.m)

# ############# Kaplan - Meier curve by RFS ########
# 
# # install.packages("survminer")
# library("survminer")
# 
# data.km <- AS3.total
# # data.km <- data.km[-which(data.km$time == 2250),]
# data.km$RFS <- ifelse(data.km$RFS.o > median(data.km$RFS.o), 1, 0)
# fit <- survfit(Surv(time, status) ~ RFS, data = data.km)
# 
# ggsurv <- ggsurvplot(
#   fit,                     # survfit object with calculated statistics.
#   data = data.km,             # data used to fit survival curves.
#   risk.table = TRUE,       # show risk table.
#   pval = TRUE,             # show p-value of log-rank test.
#   conf.int = TRUE,         # show confidence intervals for 
#   # point estimates of survival curves.
#   palette = c("#E7B800", "#2E9FDF"),
#   # survival estimates.
#   xlab = "Time in days",   # customize X axis label.
#   break.time.by = 500,     # break X axis in time intervals by 500.
#   ggtheme = theme_light(), # customize plot and risk table with a theme.
#   risk.table.y.text.col = T,# colour risk table text annotations.
#   risk.table.height = 0.25, # the height of the risk table
#   risk.table.y.text = FALSE,# show bars instead of names in text annotations
#   # in legend of risk table.
#   
#   conf.int.style = "step",  # customize style of confidence intervals
#   surv.median.line = "hv",  # add the median survival pointer.
#   legend.labs = 
#     c("Low RFS", "High RFS")    # change legend labels.
# )
# ggsurv
# 
# data.ko <- data.o
# data.ko <- data.ko[-which(data.ko$time == 2250),]
# data.ko$RFS_3rd <- ifelse(data.ko$RFS_3rd > median(data.ko$RFS_3rd), 1, 0)
# fit <- survfit(Surv(time, status) ~ RFS_3rd, data = data.ko)
# 
# ggsurv <- ggsurvplot(
#   fit,                     # survfit object with calculated statistics.
#   data = data.ko,             # data used to fit survival curves.
#   risk.table = TRUE,       # show risk table.
#   pval = TRUE,             # show p-value of log-rank test.
#   conf.int = TRUE,         # show confidence intervals for 
#   # point estimates of survival curves.
#   palette = c("#E7B800", "#2E9FDF"),
#   # survival estimates.
#   xlab = "Time in days",   # customize X axis label.
#   break.time.by = 500,     # break X axis in time intervals by 500.
#   ggtheme = theme_light(), # customize plot and risk table with a theme.
#   risk.table.y.text.col = T,# colour risk table text annotations.
#   risk.table.height = 0.25, # the height of the risk table
#   risk.table.y.text = FALSE,# show bars instead of names in text annotations
#   # in legend of risk table.
#   
#   conf.int.style = "step",  # customize style of confidence intervals
#   surv.median.line = "hv",  # add the median survival pointer.
#   legend.labs = 
#     c("Low RFS", "High RFS")    # change legend labels.
# )
# ggsurv
# 
# 
