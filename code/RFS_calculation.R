#*******************************************************************************
#**********           < RFS-SNP interaction on T2D survival >         **********
#**********           RFS calculation                                 **********
#**********                                                           **********	
#**********           Written by:				                              **********
#**********           Chanhwa Lee - chanhwa@email.unc.edu             **********
#**********                                                           **********
#**********           Version: 1.0                                    **********
#**********           Created: July 13, 2019                          **********
#**********           Modified: Oct 16, 2020                          **********
#*******************************************************************************

library(dplyr)

FFQ <- read.csv("AS3_19_FFQ.csv") ; dim(FFQ)
MM2_Diet <- read.csv("MM2_16_DIET.csv") ; dim(MM2_Diet)
RGMEAL <- data.frame(DIST_ID = MM2_Diet$DIST_ID, MEALFQ = MM2_Diet$A00_RMEALFQA)

FFQ <- FFQ[order(FFQ$DIST_ID,decreasing = F),]
RGMEAL <- RGMEAL[order(RGMEAL$DIST_ID,decreasing = F),]

FFQ <- inner_join(FFQ,RGMEAL,by="DIST_ID") ; dim(FFQ)

## 1. Original Recommended Food Score (3rd stage)

FFQ.o <- FFQ
w <- which(FFQ.o$MEALFQ > 10 | FFQ.o$AS3_F003_FQ > 10 | FFQ.o$AS3_F024_FQ > 10 | FFQ.o$AS3_F025_FQ > 10 
           | FFQ.o$AS3_F027_FQ > 10 | FFQ.o$AS3_F088_FQ > 10 | FFQ.o$AS3_F038_FQ > 10 | FFQ.o$AS3_F039_FQ > 10 
           | FFQ.o$AS3_F040_FQ > 10 | FFQ.o$AS3_F041_FQ > 10 | FFQ.o$AS3_F042_FQ > 10 | FFQ.o$AS3_F043_FQ > 10 
           | FFQ.o$AS3_F045_FQ > 10 | FFQ.o$AS3_F047_FQ > 10 | FFQ.o$AS3_F048_FQ > 10 | FFQ.o$AS3_F049_FQ > 10 
           | FFQ.o$AS3_F050_FQ > 10 | FFQ.o$AS3_F051_FQ > 10 | FFQ.o$AS3_F052_FQ > 10 | FFQ.o$AS3_F053_FQ > 10 
           | FFQ.o$AS3_F054_FQ > 10 | FFQ.o$AS3_F055_FQ > 10 | FFQ.o$AS3_F056_FQ > 10 | FFQ.o$AS3_F095_FQ > 10 
           | FFQ.o$AS3_F096_FQ > 10 | FFQ.o$AS3_F097_FQ > 10 | FFQ.o$AS3_F098_FQ > 10 | FFQ.o$AS3_F099_FQ > 10 
           | FFQ.o$AS3_F100_FQ > 10 | FFQ.o$AS3_F101_FQ > 10 | FFQ.o$AS3_F102_FQ > 10 | FFQ.o$AS3_F103_FQ > 10 
           | FFQ.o$AS3_F104_FQ > 10 | FFQ.o$AS3_F105_FQ > 10 | FFQ.o$AS3_F106_FQ > 10 | FFQ.o$AS3_F068_FQ > 10 
           | FFQ.o$AS3_F069_FQ > 10 | FFQ.o$AS3_F071_FQ > 10 | FFQ.o$AS3_F072_FQ > 10 | FFQ.o$AS3_F074_FQ > 10 
           | FFQ.o$AS3_F082_FQ > 10 | FFQ.o$AS3_F083_FQ > 10 | FFQ.o$AS3_F084_FQ > 10 | FFQ.o$AS3_F085_FQ > 10 
           | FFQ.o$AS3_F087_FQ > 10 | FFQ.o$AS3_F023_FQ > 10 | FFQ.o$AS3_F092_FQ > 10)

# Delete Missing : 183 samples are filtered.
length(w)
FFQ.o <- FFQ.o[-w,]
dim(FFQ.o)

for(i in 1:nrow(FFQ.o)){ 
  # Meals(1)
  if	(FFQ.o$MEALFQ[i]>=3)	FFQ.o$AS3_MEAL_RFS[i]=	1	else FFQ.o$AS3_MEAL_RFS[i]	=	0 # 하루 정규식사 횟수
  # Grain(1)
  if	(FFQ.o$AS3_F003_FQ[i] >=4)	FFQ.o$AS3_MIX_RFS[i]= 1	else FFQ.o$AS3_MIX_RFS[i]	=	0	# 잡곡밥 섭취빈도
  # Legumes(4)
  if	(FFQ.o$AS3_F024_FQ[i]>=4)	FFQ.o$AS3_CONG_RFS[i]	=	1		else FFQ.o$AS3_CONG_RFS[i]	=	0	# 콩/콩자반 섭취빈도
  if	(FFQ.o$AS3_F025_FQ[i]>=4)	FFQ.o$AS3_DJAN_RFS[i]	=	1		else FFQ.o$AS3_DJAN_RFS[i]	=	0	# 된장국/청국장 섭취빈도
  if	(FFQ.o$AS3_F027_FQ[i]>=4)	FFQ.o$AS3_DUBU_RFS[i]	=	1		else FFQ.o$AS3_DUBU_RFS[i]	=	0	# 두부 섭취빈도
  if	(FFQ.o$AS3_F088_FQ[i]>=4)	FFQ.o$AS3_DUYU_RFS[i]	=	1		else FFQ.o$AS3_DUYU_RFS[i]	=	0	# 두유 섭취빈도
  # Vegetable(17)
  if	(FFQ.o$AS3_F038_FQ[i]>=4)	FFQ.o$AS3_KOCB_RFS[i]	=	1		else FFQ.o$AS3_KOCB_RFS[i]	=	0	# 배추/배추국 섭취빈도
  if	(FFQ.o$AS3_F039_FQ[i]>=4)	FFQ.o$AS3_SIGM_RFS[i]	=	1		else FFQ.o$AS3_SIGM_RFS[i]	=	0	# 시금치 섭취빈도
  if	(FFQ.o$AS3_F040_FQ[i]>=4)	FFQ.o$AS3_SANC_RFS[i] =	1		else FFQ.o$AS3_SANC_RFS[i] =	0	# 상추 섭취빈도
  if	(FFQ.o$AS3_F041_FQ[i]>=4)	FFQ.o$AS3_GGAE_RFS[i] =	1		else FFQ.o$AS3_GGAE_RFS[i] =	0	# 들깻잎 섭취빈도
  if	(FFQ.o$AS3_F042_FQ[i]>=4)	FFQ.o$AS3_WRVG_RFS[i] =	1		else FFQ.o$AS3_WRVG_RFS[i] =	0	# 야채쌈/야채샐러드 섭취빈도
  if	(FFQ.o$AS3_F043_FQ[i]>=4)	FFQ.o$AS3_GRVG_RFS[i] =	1		else FFQ.o$AS3_GRVG_RFS[i] =	0	# 기타녹색채소(냉이,쑥,갓,아욱,취,무청 등) 섭취빈도도
  if	(FFQ.o$AS3_F045_FQ[i]>=4)	FFQ.o$AS3_SUKJ_RFS[i] =	1		else FFQ.o$AS3_SUKJ_RFS[i] =	0	# 콩나물/숙주나물 섭취빈도
  if	(FFQ.o$AS3_F047_FQ[i]>=4)	FFQ.o$AS3_NTAR_RFS[i]	=	1		else FFQ.o$AS3_NTAR_RFS[i]	=	0	# 느타리버섯 섭취빈도
  if	(FFQ.o$AS3_F048_FQ[i]>=4)	FFQ.o$AS3_MUSH_RFS[i] =	1		else FFQ.o$AS3_MUSH_RFS[i] =	0	# 기타버섯(양송이, 팽이, 표고) 섭취빈도
  if	(FFQ.o$AS3_F049_FQ[i]>=4)	FFQ.o$AS3_GCIP_RFS[i] =	1		else FFQ.o$AS3_GCIP_RFS[i] =	0	# 고춧잎 섭취빈도
  if	(FFQ.o$AS3_F050_FQ[i]>=4)	FFQ.o$AS3_BUCH_RFS[i] =	1		else FFQ.o$AS3_BUCH_RFS[i] =	0	# 부추/미나리 섭취빈도
  if	(FFQ.o$AS3_F051_FQ[i]>=4)	FFQ.o$AS3_OUI_RFS[i] =	1		else FFQ.o$AS3_OUI_RFS[i] =	0	# 오이 섭취빈도
  if	(FFQ.o$AS3_F052_FQ[i]>=4)	FFQ.o$AS3_CARR_RFS[i] =	1		else FFQ.o$AS3_CARR_RFS[i] =	0	# 당근/당근쥬스 섭취빈도
  if	(FFQ.o$AS3_F053_FQ[i]>=4)	FFQ.o$AS3_ONIO_RFS[i]	=	1		else FFQ.o$AS3_ONIO_RFS[i]	=	0	# 양파 섭취빈도
  if	(FFQ.o$AS3_F054_FQ[i]>=4)	FFQ.o$AS3_GOCH_RFS[i] =	1		else FFQ.o$AS3_GOCH_RFS[i] =	0	# 풋고추 섭취빈도
  if	(FFQ.o$AS3_F055_FQ[i]>=4)	FFQ.o$AS3_YOPU_RFS[i] =	1		else FFQ.o$AS3_YOPU_RFS[i] =	0	# 애호박 섭취빈도
  if	(FFQ.o$AS3_F056_FQ[i]>=4)	FFQ.o$AS3_PUMK_RFS[i] =	1		else FFQ.o$AS3_PUMK_RFS[i] =	0	# 늙은호박(호박죽)/호박즙 섭취빈도
  # Fruits(12)
  if	(FFQ.o$AS3_F095_FQ[i]>=4)	FFQ.o$AS3_STBR_RFS[i] =	1		else FFQ.o$AS3_STBR_RFS[i] =	0	# 딸기 섭취빈도
  if	(FFQ.o$AS3_F096_FQ[i]>=4)	FFQ.o$AS3_MELN_RFS[i] =	1		else FFQ.o$AS3_MELN_RFS[i] =	0	# 멜론 섭취빈도
  if	(FFQ.o$AS3_F097_FQ[i]>=4)	FFQ.o$AS3_WTML_RFS[i] =	1		else FFQ.o$AS3_WTML_RFS[i] =	0	# 수박 섭취빈도
  if	(FFQ.o$AS3_F098_FQ[i]>=4)	FFQ.o$AS3_PITC_RFS[i] =	1		else FFQ.o$AS3_PITC_RFS[i] =	0	# 복숭아/자두 섭취빈도
  if	(FFQ.o$AS3_F099_FQ[i]>=4)	FFQ.o$AS3_BANA_RFS[i] =	1		else FFQ.o$AS3_BANA_RFS[i] =	0	# 바나나 섭취빈도
  if	(FFQ.o$AS3_F100_FQ[i]>=4)	FFQ.o$AS3_PERS_RFS[i] =	1		else FFQ.o$AS3_PERS_RFS[i] =	0	# 감 섭취빈도
  if	(FFQ.o$AS3_F101_FQ[i]>=4)	FFQ.o$AS3_GYUL_RFS[i] =	1		else FFQ.o$AS3_GYUL_RFS[i] =	0	# 귤 섭취빈도
  if	(FFQ.o$AS3_F102_FQ[i]>=4)	FFQ.o$AS3_PEAR_RFS[i] =	1		else FFQ.o$AS3_PEAR_RFS[i] =	0	# 배 섭취빈도
  if	(FFQ.o$AS3_F103_FQ[i]>=4)	FFQ.o$AS3_APPL_RFS[i]	=	1		else FFQ.o$AS3_APPL_RFS[i]	=	0	# 사과/사과쥬스 섭취빈도
  if	(FFQ.o$AS3_F104_FQ[i]>=4)	FFQ.o$AS3_ORAN_RFS[i] =	1		else FFQ.o$AS3_ORAN_RFS[i] =	0	# 오렌지/오렌지쥬스 섭취빈도
  if	(FFQ.o$AS3_F105_FQ[i]>=4)	FFQ.o$AS3_GRAP_RFS[i] =	1		else FFQ.o$AS3_GRAP_RFS[i] =	0	# 포도/포도쥬스 섭취빈도
  if	(FFQ.o$AS3_F106_FQ[i]>=4)	FFQ.o$AS3_TOMA_RFS[i]	=	1		else FFQ.o$AS3_TOMA_RFS[i]	=	0	# 토마토/토마토쥬스 섭취빈도
  # Fish(5)
  if	(FFQ.o$AS3_F068_FQ[i]>=4)	FFQ.o$AS3_BLFH_RFS[i]	=	1		else FFQ.o$AS3_BLFH_RFS[i]	=	0	# 등푸른생선(고등어,꽁치,삼치) 섭취빈도
  if	(FFQ.o$AS3_F069_FQ[i]>=4)	FFQ.o$AS3_GALF_RFS[i] =	1		else FFQ.o$AS3_GALF_RFS[i] =	0	# 갈치 섭취빈도
  if	(FFQ.o$AS3_F071_FQ[i]>=4)	FFQ.o$AS3_JOGI_RFS[i] =	1		else FFQ.o$AS3_JOGI_RFS[i] =	0	# 조기 / 넙치 섭취빈도
  if	(FFQ.o$AS3_F072_FQ[i]>=4)	FFQ.o$AS3_TAEF_RFS[i] =	1		else FFQ.o$AS3_TAEF_RFS[i] =	0	# 명태, 동태, 북어 섭취빈도
  if	(FFQ.o$AS3_F074_FQ[i]>=4)	FFQ.o$AS3_ANCH_RFS[i] =	1		else FFQ.o$AS3_ANCH_RFS[i] =	0	# 멸치, 멸치볶음 섭취빈도
  # Seaweed(2)
  if	(FFQ.o$AS3_F082_FQ[i]>=4)	FFQ.o$AS3_GIM_RFS[i] =	1		else FFQ.o$AS3_GIM_RFS[i] =	0	# 김 섭취빈도
  if	(FFQ.o$AS3_F083_FQ[i]>=4)	FFQ.o$AS3_DASM_RFS[i] =	1		else FFQ.o$AS3_DASM_RFS[i] =	0	# 다시마/미역 섭취빈도
  # Diary(3)
  if	(FFQ.o$AS3_F084_FQ[i]>=4)	FFQ.o$AS3_MILK_RFS[i] =	1		else FFQ.o$AS3_MILK_RFS[i] =	0	# 우유 섭취빈도
  if	(FFQ.o$AS3_F085_FQ[i]>=4)	FFQ.o$AS3_YOGU_RFS[i]	=	1		else FFQ.o$AS3_YOGU_RFS[i]	=	0	# 요구르트, 요플레 섭취빈도
  if	(FFQ.o$AS3_F087_FQ[i]>=4)	FFQ.o$AS3_CHEZ_RFS[i] =	1		else FFQ.o$AS3_CHEZ_RFS[i] =	0	# 치즈 섭취빈도
  # Nuts(1)
  if	(FFQ.o$AS3_F023_FQ[i]>=4)	FFQ.o$AS3_NUT_RFS[i] =	1		else FFQ.o$AS3_NUT_RFS[i] =	0	# 땅콩/아몬드/잣 섭취빈도
  # Tea (1)
  if	(FFQ.o$AS3_F092_FQ[i]>=4)	FFQ.o$AS3_TEA_RFS[i] =	1		else FFQ.o$AS3_TEA_RFS[i] =	0	# 녹차 섭취빈도
}
dim(FFQ.o)

for(i in 1:nrow(FFQ.o)){
  FFQ.o$RFS.o[i] = sum(FFQ.o[i,227:273])
}

summary(FFQ.o$RFS.o)
table(FFQ.o$RFS.o)
hist(FFQ.o$RFS.o)

colnames(FFQ.o)
AS3.RFS.o <- data.frame(FFQ.o[,c(1,227:274)])
colnames(AS3.RFS.o) ; ncol(AS3.RFS.o)
write.csv(AS3.RFS.o,"AS3.RFS.o.csv")

## 2. Modified Recommended Food Score (3rd stage)

FFQ.m <- FFQ
w <- which(FFQ.m$MEALFQ > 10 | FFQ.m$AS3_F003_FQ > 10 | FFQ.m$AS3_F002_FQ > 10 | FFQ.m$AS3_F025_FQ > 10 
           | FFQ.m$AS3_F027_FQ > 10 | FFQ.m$AS3_F088_FQ > 10 | FFQ.m$AS3_F038_FQ > 10 | FFQ.m$AS3_F039_FQ > 10 
           | FFQ.m$AS3_F040_FQ > 10 | FFQ.m$AS3_F041_FQ > 10 | FFQ.m$AS3_F042_FQ > 10 | FFQ.m$AS3_F043_FQ > 10 
           | FFQ.m$AS3_F045_FQ > 10 | FFQ.m$AS3_F047_FQ > 10 | FFQ.m$AS3_F048_FQ > 10 | FFQ.m$AS3_F049_FQ > 10 
           | FFQ.m$AS3_F050_FQ > 10 | FFQ.m$AS3_F051_FQ > 10 | FFQ.m$AS3_F052_FQ > 10 | FFQ.m$AS3_F053_FQ > 10 
           | FFQ.m$AS3_F054_FQ > 10 | FFQ.m$AS3_F055_FQ > 10 | FFQ.m$AS3_F056_FQ > 10 | FFQ.m$AS3_F095_FQ > 10 
           | FFQ.m$AS3_F096_FQ > 10 | FFQ.m$AS3_F097_FQ > 10 | FFQ.m$AS3_F098_FQ > 10 | FFQ.m$AS3_F099_FQ > 10 
           | FFQ.m$AS3_F100_FQ > 10 | FFQ.m$AS3_F101_FQ > 10 | FFQ.m$AS3_F102_FQ > 10 | FFQ.m$AS3_F103_FQ > 10 
           | FFQ.m$AS3_F104_FQ > 10 | FFQ.m$AS3_F105_FQ > 10 | FFQ.m$AS3_F106_FQ > 10 | FFQ.m$AS3_F069_FQ > 10 
           | FFQ.m$AS3_F071_FQ > 10 | FFQ.m$AS3_F074_FQ > 10 | FFQ.m$AS3_F082_FQ > 10 | FFQ.m$AS3_F083_FQ > 10 
           | FFQ.m$AS3_F084_FQ > 10 | FFQ.m$AS3_F085_FQ > 10 | FFQ.m$AS3_F087_FQ > 10 | FFQ.m$AS3_F023_FQ > 10 
           | FFQ.m$AS3_F092_FQ > 10 | FFQ.m$AS3_F044_FQ > 10 | FFQ.m$AS3_F046_FQ > 10 | FFQ.m$AS3_F037_FQ > 10 
           | FFQ.m$AS3_F070_FQ > 10 | FFQ.m$AS3_F081_FQ > 10 | FFQ.m$AS3_F068_FQ > 10 | FFQ.m$AS3_F068_FQ > 10 
           | FFQ.m$AS3_F067_FQ > 10 | FFQ.m$AS3_F075_FQ > 10)



# Delete Missing : 183 samples are filtered.
length(w)
FFQ.m <- FFQ.m[-w,]
dim(FFQ.m)

for(i in 1:nrow(FFQ.m)){ 
  # Meals(1)
  if	(FFQ.m$MEALFQ[i]>=3)	FFQ.m$AS3_MEAL_RFS[i]=	1	else FFQ.m$AS3_MEAL_RFS[i]	=	0 # 하루 정규식사 횟수
  
  # Grain(2)
  if	(FFQ.m$AS3_F003_FQ[i] >=4)	FFQ.m$AS3_MIX_RFS[i]= 1	else FFQ.m$AS3_MIX_RFS[i]	=	0	# 잡곡밥 섭취빈도
  
  # if	(AS3_F003_FQ[i] >=4)	FFQ.m$AS3_MIX_RFS[i]= 1	else FFQ.m$AS3_MIX_RFS[i]	=	0	# 보리 관련 데이터 없음 -> AS3_Diet csv file not exist
  
  # Legumes(4)
  # if	(AS3_F024_FQ[i]>=4)	FFQ.m$AS3_CONG_RFS[i]	=	1		else FFQ.m$AS3_CONG_RFS[i]	=	0	# 콩/콩자반 섭취빈도
  if	(FFQ.m$AS3_F002_FQ[i]>=4)	FFQ.m$AS3_CONG_RFS[i]	=	1		else FFQ.m$AS3_CONG_RFS[i]	=	0	# 콩밥 섭취빈도
  if	(FFQ.m$AS3_F025_FQ[i]>=4)	FFQ.m$AS3_DJAN_RFS[i]	=	1		else FFQ.m$AS3_DJAN_RFS[i]	=	0	# 된장국/청국장 섭취빈도
  if	(FFQ.m$AS3_F027_FQ[i]>=4)	FFQ.m$AS3_DUBU_RFS[i]	=	1		else FFQ.m$AS3_DUBU_RFS[i]	=	0	# 두부 섭취빈도
  if	(FFQ.m$AS3_F088_FQ[i]>=4)	FFQ.m$AS3_DUYU_RFS[i]	=	1		else FFQ.m$AS3_DUYU_RFS[i]	=	0	# 두유 섭취빈도
  
  # Vegetable(20)
  if	(FFQ.m$AS3_F038_FQ[i]>=4)	FFQ.m$AS3_KOCB_RFS[i]	=	1		else FFQ.m$AS3_KOCB_RFS[i]	=	0	# 배추/배추국 섭취빈도
  if	(FFQ.m$AS3_F039_FQ[i]>=4)	FFQ.m$AS3_SIGM_RFS[i]	=	1		else FFQ.m$AS3_SIGM_RFS[i]	=	0	# 시금치 섭취빈도
  if	(FFQ.m$AS3_F040_FQ[i]>=4)	FFQ.m$AS3_SANC_RFS[i] =	1		else FFQ.m$AS3_SANC_RFS[i] =	0	# 상추 섭취빈도
  if	(FFQ.m$AS3_F041_FQ[i]>=4)	FFQ.m$AS3_GGAE_RFS[i] =	1		else FFQ.m$AS3_GGAE_RFS[i] =	0	# 들깻잎 섭취빈도
  if	(FFQ.m$AS3_F042_FQ[i]>=4)	FFQ.m$AS3_WRVG_RFS[i] =	1		else FFQ.m$AS3_WRVG_RFS[i] =	0	# 야채쌈/야채샐러드 섭취빈도
  if	(FFQ.m$AS3_F043_FQ[i]>=4)	FFQ.m$AS3_GRVG_RFS[i] =	1		else FFQ.m$AS3_GRVG_RFS[i] =	0	# 기타녹색채소(냉이,쑥,갓,아욱,취,무청 등) 섭취빈도도
  if	(FFQ.m$AS3_F045_FQ[i]>=4)	FFQ.m$AS3_SUKJ_RFS[i] =	1		else FFQ.m$AS3_SUKJ_RFS[i] =	0	# 콩나물/숙주나물 섭취빈도
  if	(FFQ.m$AS3_F047_FQ[i]>=4)	FFQ.m$AS3_NTAR_RFS[i]	=	1		else FFQ.m$AS3_NTAR_RFS[i]	=	0	# 느타리버섯 섭취빈도
  if	(FFQ.m$AS3_F048_FQ[i]>=4)	FFQ.m$AS3_MUSH_RFS[i] =	1		else FFQ.m$AS3_MUSH_RFS[i] =	0	# 기타버섯(양송이, 팽이, 표고) 섭취빈도
  if	(FFQ.m$AS3_F049_FQ[i]>=4)	FFQ.m$AS3_GCIP_RFS[i] =	1		else FFQ.m$AS3_GCIP_RFS[i] =	0	# 고춧잎 섭취빈도
  if	(FFQ.m$AS3_F050_FQ[i]>=4)	FFQ.m$AS3_BUCH_RFS[i] =	1		else FFQ.m$AS3_BUCH_RFS[i] =	0	# 부추/미나리 섭취빈도
  if	(FFQ.m$AS3_F051_FQ[i]>=4)	FFQ.m$AS3_OUI_RFS[i] =	1		else FFQ.m$AS3_OUI_RFS[i] =	0	# 오이 섭취빈도
  if	(FFQ.m$AS3_F052_FQ[i]>=4)	FFQ.m$AS3_CARR_RFS[i] =	1		else FFQ.m$AS3_CARR_RFS[i] =	0	# 당근/당근쥬스 섭취빈도
  if	(FFQ.m$AS3_F053_FQ[i]>=4)	FFQ.m$AS3_ONIO_RFS[i]	=	1		else FFQ.m$AS3_ONIO_RFS[i]	=	0	# 양파 섭취빈도
  if	(FFQ.m$AS3_F054_FQ[i]>=4)	FFQ.m$AS3_GOCH_RFS[i] =	1		else FFQ.m$AS3_GOCH_RFS[i] =	0	# 풋고추 섭취빈도
  if	(FFQ.m$AS3_F055_FQ[i]>=4)	FFQ.m$AS3_YOPU_RFS[i] =	1		else FFQ.m$AS3_YOPU_RFS[i] =	0	# 애호박 섭취빈도
  if	(FFQ.m$AS3_F056_FQ[i]>=4)	FFQ.m$AS3_PUMK_RFS[i] =	1		else FFQ.m$AS3_PUMK_RFS[i] =	0	# 늙은호박(호박죽)/호박즙 섭취빈도
  
  if	(FFQ.m$AS3_F044_FQ[i]>=4)	FFQ.m$AS3_DORA_RFS[i] =	1		else FFQ.m$AS3_DORA_RFS[i] =	0	# 도라지 섭취빈도
  if	(FFQ.m$AS3_F046_FQ[i]>=4)	FFQ.m$AS3_GOSA_RFS[i] =	1		else FFQ.m$AS3_GOSA_RFS[i] =	0	# 고사리 섭취빈도
  if	(FFQ.m$AS3_F037_FQ[i]>=4)	FFQ.m$AS3_RDSH_RFS[i] =	1		else FFQ.m$AS3_RDSH_RFS[i] =	0	# 무 섭취빈도
  
  # Fruits(12)
  if	(FFQ.m$AS3_F095_FQ[i]>=4)	FFQ.m$AS3_STBR_RFS[i] =	1		else FFQ.m$AS3_STBR_RFS[i] =	0	# 딸기 섭취빈도
  if	(FFQ.m$AS3_F096_FQ[i]>=4)	FFQ.m$AS3_MELN_RFS[i] =	1		else FFQ.m$AS3_MELN_RFS[i] =	0	# 멜론 섭취빈도
  if	(FFQ.m$AS3_F097_FQ[i]>=4)	FFQ.m$AS3_WTML_RFS[i] =	1		else FFQ.m$AS3_WTML_RFS[i] =	0	# 수박 섭취빈도
  if	(FFQ.m$AS3_F098_FQ[i]>=4)	FFQ.m$AS3_PITC_RFS[i] =	1		else FFQ.m$AS3_PITC_RFS[i] =	0	# 복숭아/자두 섭취빈도
  if	(FFQ.m$AS3_F099_FQ[i]>=4)	FFQ.m$AS3_BANA_RFS[i] =	1		else FFQ.m$AS3_BANA_RFS[i] =	0	# 바나나 섭취빈도
  if	(FFQ.m$AS3_F100_FQ[i]>=4)	FFQ.m$AS3_PERS_RFS[i] =	1		else FFQ.m$AS3_PERS_RFS[i] =	0	# 감 섭취빈도
  if	(FFQ.m$AS3_F101_FQ[i]>=4)	FFQ.m$AS3_GYUL_RFS[i] =	1		else FFQ.m$AS3_GYUL_RFS[i] =	0	# 귤 섭취빈도
  if	(FFQ.m$AS3_F102_FQ[i]>=4)	FFQ.m$AS3_PEAR_RFS[i] =	1		else FFQ.m$AS3_PEAR_RFS[i] =	0	# 배 섭취빈도
  if	(FFQ.m$AS3_F103_FQ[i]>=4)	FFQ.m$AS3_APPL_RFS[i]	=	1		else FFQ.m$AS3_APPL_RFS[i]	=	0	# 사과/사과쥬스 섭취빈도
  if	(FFQ.m$AS3_F104_FQ[i]>=4)	FFQ.m$AS3_ORAN_RFS[i] =	1		else FFQ.m$AS3_ORAN_RFS[i] =	0	# 오렌지/오렌지쥬스 섭취빈도
  if	(FFQ.m$AS3_F105_FQ[i]>=4)	FFQ.m$AS3_GRAP_RFS[i] =	1		else FFQ.m$AS3_GRAP_RFS[i] =	0	# 포도/포도쥬스 섭취빈도
  if	(FFQ.m$AS3_F106_FQ[i]>=4)	FFQ.m$AS3_TOMA_RFS[i]	=	1		else FFQ.m$AS3_TOMA_RFS[i]	=	0	# 토마토/토마토쥬스 섭취빈도
  
  # Fish(9)
  #  if	(AS3_F068_FQ[i]>=4)	FFQ.m$AS3_BLFH_RFS[i]	=	1		else FFQ.m$AS3_BLFH_RFS[i]	=	0	# 등푸른생선(고등어,꽁치,삼치) 섭취빈도
  if	(FFQ.m$AS3_F069_FQ[i]>=4)	FFQ.m$AS3_GALF_RFS[i] =	1		else FFQ.m$AS3_GALF_RFS[i] =	0	# 갈치 섭취빈도
  if	(FFQ.m$AS3_F071_FQ[i]>=4)	FFQ.m$AS3_JOGI_RFS[i] =	1		else FFQ.m$AS3_JOGI_RFS[i] =	0	# 조기 / 넙치 섭취빈도
  #  if	(AS3_F072_FQ[i]>=4)	FFQ.m$AS3_TAEF_RFS[i] =	1		else FFQ.m$AS3_TAEF_RFS[i] =	0	# 명태, 동태, 북어 섭취빈도
  if	(FFQ.m$AS3_F074_FQ[i]>=4)	FFQ.m$AS3_ANCH_RFS[i] =	1		else FFQ.m$AS3_ANCH_RFS[i] =	0	# 멸치, 멸치볶음 섭취빈도

  if	(FFQ.m$AS3_F070_FQ[i]>=4)	FFQ.m$AS3_EEL_RFS[i] =	1		else FFQ.m$AS3_EEL_RFS[i] =	0	# 장어 섭취빈도
  if	(FFQ.m$AS3_F081_FQ[i]>=4)	FFQ.m$AS3_FSCK_RFS[i] =	1		else FFQ.m$AS3_FSCK_RFS[i] =	0	# 어묵 섭취빈도
  if	(FFQ.m$AS3_F068_FQ[i]>=4)	FFQ.m$AS3_MCKR_RFS[i] =	1		else FFQ.m$AS3_MCKR_RFS[i] =	0	# 고등어 섭취빈도 >>> 따로 구분이 안되어 있고 등푸른 생선에 포함됨
  if	(FFQ.m$AS3_F068_FQ[i]>=4)	FFQ.m$AS3_MCPK_RFS[i] =	1		else FFQ.m$AS3_MCPK_RFS[i] =	0	# 꽁치 섭취빈도 >>> 따로 구분이 안되어 있고 등푸른 생선에 포함됨
  if	(FFQ.m$AS3_F067_FQ[i]>=4)	FFQ.m$AS3_SASH_RFS[i] =	1		else FFQ.m$AS3_SASH_RFS[i] =	0	# 회 섭취빈도
  if	(FFQ.m$AS3_F075_FQ[i]>=4)	FFQ.m$AS3_TUNA_RFS[i] =	1		else FFQ.m$AS3_TUNA_RFS[i] =	0	# 참치 섭취빈도
  # Seaweed(2)
  if	(FFQ.m$AS3_F082_FQ[i]>=4)	FFQ.m$AS3_GIM_RFS[i] =	1		else FFQ.m$AS3_GIM_RFS[i] =	0	# 김 섭취빈도
  if	(FFQ.m$AS3_F083_FQ[i]>=4)	FFQ.m$AS3_DASM_RFS[i] =	1		else FFQ.m$AS3_DASM_RFS[i] =	0	# 다시마/미역 섭취빈도
  # Diary(3)
  if	(FFQ.m$AS3_F084_FQ[i]>=4)	FFQ.m$AS3_MILK_RFS[i] =	1		else FFQ.m$AS3_MILK_RFS[i] =	0	# 우유 섭취빈도
  if	(FFQ.m$AS3_F085_FQ[i]>=4)	FFQ.m$AS3_YOGU_RFS[i]	=	1		else FFQ.m$AS3_YOGU_RFS[i]	=	0	# 요구르트, 요플레 섭취빈도
  if	(FFQ.m$AS3_F087_FQ[i]>=4)	FFQ.m$AS3_CHEZ_RFS[i] =	1		else FFQ.m$AS3_CHEZ_RFS[i] =	0	# 치즈 섭취빈도
  # Nuts(1)
  if	(FFQ.m$AS3_F023_FQ[i]>=4)	FFQ.m$AS3_NUT_RFS[i] =	1		else FFQ.m$AS3_NUT_RFS[i] =	0	# 땅콩/아몬드/잣 섭취빈도
  # Tea (1)
  if	(FFQ.m$AS3_F092_FQ[i]>=4)	FFQ.m$AS3_TEA_RFS[i] =	1		else FFQ.m$AS3_TEA_RFS[i] =	0	# 녹차 섭취빈도
}
dim(FFQ.m)
colnames(FFQ.m)
for(i in 1:nrow(FFQ.m)){
  FFQ.m$RFS.m[i] = sum(FFQ.m[i,227:280])
}

summary(FFQ.m$RFS.m)
table(FFQ.m$RFS.m)
hist(FFQ.m$RFS.m)

colnames(FFQ.m)
AS3.RFS.m <- data.frame(FFQ.m[,c(1,227:281)])
colnames(AS3.RFS.m) ; ncol(AS3.RFS.m)
write.csv(AS3.RFS.m,"AS3.RFS.m.csv")

####### 
ncol(AS3.RFS.m)
RFS.o <- AS3.RFS.o[,c(1,49)]; dim(RFS.o)
RFS.m <- AS3.RFS.m[,c(1,56)]; dim(RFS.m)
RFSinfo <- inner_join(RFS.o,RFS.m,by="DIST_ID") ; dim(RFSinfo)
write.csv(RFSinfo,"RFSinfo.csv")

plot(RFSinfo$RFS.o,RFSinfo$RFS.m)
AS3.RFS.o[AS3.RFS.o$DIST_ID == "NIH16K7080495",]
AS3.RFS.m[AS3.RFS.m$DIST_ID == "NIH16K7080495",]
sum(RFSinfo$RFS.o==RFSinfo$RFS.m)
