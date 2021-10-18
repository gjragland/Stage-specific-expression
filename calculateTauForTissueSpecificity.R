#R
#GJR
#4/14/2021
# calculate tissue specificity index from FlyAtlas data


#import data from local mysql database generated from:
# motif.gla.ac.uk/downloads/FlyAtlas2_21.04.18.sql
# Above linked from fly atlas website http://flyatlas.gla.ac.uk/FlyAtlas2/index.html
# The tau values are calculated following: 
#   Cridland, J.M., Majane, A.C., Sheehy, H.K. and Begun, D.J., 2020. Polymorphism and Divergence of Novel Gene Expression Patterns in Drosophila melanogaster. Genetics, 216(1), pp.79-93.
# Original development of the tau estimator:
# Yanai, I., Benjamin, H., Shmoish, M., Chalifa-Caspi, V., Shklar, M., Ophir, R., Bar-Even, A., Horn-Saban, S., Safran, M., Domany, E. and Lancet, D., 2005. Genome-wide midrange transcription profiles reveal expression level relationships in human tissue specification. Bioinformatics, 21(5), pp.650-659.

library(RMySQL)
library(tidyr)

#import info on tissues and ids
db='flyAtlas2'
mydb = dbConnect(MySQL(), user='xx', password='xx', dbname=db)
query<-'select * from Tissue;'
tissueInfo<-dbGetQuery(mydb, query)
tissueInfo<-tissueInfo[tissueInfo$Replicates != 0,]

#import FPKM data
query='select FBgn, TissueID, FPKM from GeneFPKM;'
geneFPKM<-dbGetQuery(mydb, query)
geneFPKM<-geneFPKM[geneFPKM$TissueID %in% tissueInfo$TissueID,]

#gather data to one row per flybase gene id
a<-geneFPKM %>% pivot_wider(names_from = TissueID, values_from = FPKM)
wideFPKM<-as.data.frame(a)

#categorize tissue ids to male, female, and larva
maleTissueIds<-tissueInfo$TissueID[tissueInfo$Sex=='Male']
femaleTissueIds<-tissueInfo$TissueID[tissueInfo$Sex=='Female']
larvalTissueIds<-tissueInfo$TissueID[tissueInfo$Stage=='Larval']

#set minimum FPKM value of whole bodies to '2'
# adopted from Begun dros tissue expression paper, 2020, also seems to be a feature of a Mank et al. am nat paper
# Dean, R. and Mank, J.E., 2016. Tissue specificity and sex-specific regulatory variation permit the evolution of sex-biased gene expression. The American Naturalist, 188(3), pp.E74-E84.
wideFPKM$`100`[wideFPKM$`100` <= 2]<-2
wideFPKM$`200`[wideFPKM$`200` <= 2]<-2
wideFPKM$`300`[wideFPKM$`300` <= 2]<-2

#create function to normalize expression values to whole bodies (FPKM tissue / FPKM whole body) per sex/stage
maleInd<-names(wideFPKM)[-1] %in% as.character(maleTissueIds)
femaleInd<-names(wideFPKM)[-1] %in% as.character(femaleTissueIds)
larvalInd<-names(wideFPKM)[-1] %in% as.character(larvalTissueIds)

calcTau<-function(x) {
  xi<-x/max(x)
  ( sum( 1-xi ) ) / (length(xi) - 1)
}

#longish loop, runs ~1 minute or so
tau<-vector(length=nrow(wideFPKM))
for (i in 1:nrow(wideFPKM)) {
  expressionVector<-as.numeric(wideFPKM[i,-1])
  
  #FPKM >= 2 in order to be considered expressed
  expressionVector[expressionVector < 2]<-0
  
  denom<-tail(expressionVector[larvalInd],1)
  expressionVector[larvalInd]<-expressionVector[larvalInd]/denom
  
  denom<-tail(expressionVector[maleInd],1)
  expressionVector[maleInd]<-expressionVector[maleInd]/denom
  
  denom<-tail(expressionVector[femaleInd],1)
  expressionVector[femaleInd]<-expressionVector[femaleInd]/denom
  
  #remove three entries for whole body
  expressionVector<-as.numeric(expressionVector[-c(38:40)])
  
  tau[i]<-NA
  if (any(expressionVector!=0)) {
    tau[i]<-calcTau(expressionVector)
  }
}
#note, leaves 2655 genes with 0 expression across the board, tau==NA
sum(is.na(tau))

#and 2386 genes expressed (under above filters) in only one tissue
sum(tau==1,na.rm = T)

#heavily left-skewed distribution, most genes with high tau
hist(tau)

#export flybase ids and tau estimates
out<-data.frame(FBgn=wideFPKM$FBgn,tau=tau)
write.table(out,'/media/raglandlab/ExtraDrive1/dmel/FlyAtlas2_data/tauEstimatesPerFBgn.txt',sep="\t",row.names=F,quote=F)

out<-data.frame(FBgn=wideFPKM$FBgn)
write.table(out,'/media/raglandlab/ExtraDrive1/dmel/FlyAtlas2_data/FBgn.txt',sep="\t",row.names=F,quote=F)

############################################################################
#### test whether tissue specificity is correlated with overall expression #
####### as measured as FPKM in adult males, females, and larvvae           #
############################################################################
#no evidence for simple (linear) relationship between the two, though plotting the two suggests that
# the highest expressed transcripts have roughly median or higher tau scores


plot(tau$tau[wideFPKM$`100` > 2],wideFPKM$`100`[wideFPKM$`100` > 2])
plot(tau$tau[wideFPKM$`200` > 2],wideFPKM$`200`[wideFPKM$`200` > 2])
plot(tau$tau[wideFPKM$`300` > 2],wideFPKM$`300`[wideFPKM$`300` > 2])

#remove all tau=NA entries
ind<-is.na(tau$tau)
tau<-tau[!ind,]
wideFPKM<-wideFPKM[!ind,]

#note, from the 'Tissue' mysql table, tissue ids 100, 200, and 300 are adult males, adult females, and larvae, respectively
#exclude all transcripts with FPKM < 2

cor.test(wideFPKM$`100`[wideFPKM$`100` > 2],tau$tau[wideFPKM$`100` > 2],method='pearson')
# Pearson's product-moment correlation
# 
# data:  wideFPKM$`100`[wideFPKM$`100` > 2] and tau$tau[wideFPKM$`100` > 2]
# t = 1.5458, df = 10547, p-value = 0.1222
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.00403448  0.03412330
# sample estimates:
#        cor 
# 0.01504989


cor.test(wideFPKM$`200`[wideFPKM$`200` > 2],tau$tau[wideFPKM$`200` > 2],method='pearson')
# Pearson's product-moment correlation
# 
# data:  wideFPKM$`200`[wideFPKM$`200` > 2] and tau$tau[wideFPKM$`200` > 2]
# t = 1.0952, df = 8238, p-value = 0.2734
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.009528596  0.033649436
# sample estimates:
#        cor 
# 0.01206605 

cor.test(wideFPKM$`300`[wideFPKM$`300` > 2],tau$tau[wideFPKM$`300` > 2],method='pearson')
# Pearson's product-moment correlation
# 
# data:  wideFPKM$`300`[wideFPKM$`300` > 2] and tau$tau[wideFPKM$`300` > 2]
# t = 1.6164, df = 8073, p-value = 0.106
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.003825675  0.039783528
# sample estimates:
#        cor 
# 0.01798748 


