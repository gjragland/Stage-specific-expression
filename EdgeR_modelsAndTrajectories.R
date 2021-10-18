# R
# Written by Phil Freda and Greg Ragland
# PJF 3/7/17, modified by GJR through October 2021
# Analyze RNA-seq data from the Drosophila melanogaster experiment comparing larval to adult
#       expression during and following cold exposure
# Experimental Design: 
# 2 life stages - Larvae and Adults
# Larvae - 6 lines, 4 treatments (conditions), 3 replicates per lines/treatment: 0 time point (before exposure), 30 minutes in (during exposure), 60 minutes in (pulling out of the chamber), 30 mins post exposure
# Adults - 6 lines, 4 treatments (conditions), 3 replicates per lines/treatment: 0 time point (before exposure), 30 minutes in (during exposure), 60 minutes in (pulling out of the chamber), 30 mins post exposure
# Lines are each members of one of two groups:
  # HA: High performing as adults, low performing as larvae (L380,L486, plust whatever line number is in the grand mean, i.e. no paramater (L486))
  # HL: High performing as larvae, low performing as adults (L441,832,913)
# Model: This model is a stage-specific model for analyzing differential expression

# Installing/loading software


library(extrafont)
library(edgeR)
library(stringr)
library(locfit)
library(statmod)
library(reshape)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

#Clear the workspace before re-running initial steps, some tables/data are overwritten, e.g., filtered dge objects
#rm(list = ls())

################################################################################################
####### Read in and process data, fit initial EdgeR models   ###################################
################################################################################################


# Set WD:
setwd("/media/raglandlab/ExtraDrive1/dmel/RNAseqStudy/PhilRnaSeqAnalysis/DrosDgrpRnaSeqDataCounts")

#create functions to annotate to current flybase ids/symbols/gene names
flySyms<-read.table('../FlyBase_IDsConversionAndSymbols.txt',header=T,sep="\t",row.names=NULL,stringsAsFactors=F,quote="\"",comment.char="")
idToSym <- setNames(as.list(flySyms$current_symbol), flySyms$submitted_item)
getSym<-function(id) {
  flysim<-idToSym[[as.character(id)]]
  return(flysim)
}
qw <- function(x) unlist(strsplit(x, "[[:space:]]+"))

# Create a table displaying counts for each gene in each sample using the gene-level outputs from RSEM

# create a file list containing all the files with gene results (conducted from folder containing the count data)
file_list<-list.files(pattern='*genes.results')
# check the list
file_list

# empty table to populate with sample ID's with the first parts of the filenames
# RUN THIS AGAIN WHEN REMOVING ENTRIES #########################
table<-c()

## Loops through file list and pull out ID's and pulls out expected counts for each sample ID and for each gene
# Example
# gene id     sample 1      sample 2
# gene 1         5             2
# gene 2        10            400

for (file in file_list){
  sampleID<-str_match(file,"(^[0-9]*)")[,1]
  print(sampleID)
  df<-read.table(file,header=T,sep="\t",row.names=NULL)
  df2<-df[,c(1,5)]
  colnames(df2)[colnames(df2)=="expected_count"] <- paste("droso_", sampleID,"_exp_count",sep="")
  if (length(table) ==0){table<-df2
  } else {
    table<-merge(table,df2,by="gene_id") }
}  

length(table$gene_id)
# 17559 genes

# Table of annotations (skipped for now)
# geneid.flybasemerged<-read.table("gtffile.genenameid.flybase",sep="\t",header=TRUE,row.names=NULL)
# for importing a brief description of the genes (optional)
# Merge this wih the table of counts so you have annotations with the gene ID's and counts per sample
# table2.flybase<-merge(geneid.flybasemerged,table, by="gene_id",all=TRUE)

## Removes genes that do not have a least 1 count in at least 50% of the samples
filterMinCount<- function(x) {
  pres<-x >=1
  out=F
  if ((sum(pres)/length(pres)) >= 0.5) {out=T}
  return(out)
}

filterInd.flybase<-apply(table[,(-1)],1,filterMinCount) 
table<-table[filterInd.flybase,]

nrow(table) #filtered transcript number
# 13242

## Creation of the DGElist object (A matrix of counts)

colnames(table) # Check how many columns there are for table.dge
dat.dge<-DGEList(counts=table[,-1],genes=table[,1])

setwd('../')
sampleInfo<-read.csv('AllSampleInfo.csv',stringsAsFactors = F)
sampleInfo$treatment[sampleInfo$treatment==25]<-30 #time point 30 mislabeled as '25' in the info file
sampleInfo$treatment<-as.factor(sampleInfo$treatment)
sampleInfo$stage=as.factor(sampleInfo$stage)
sampleInfo$pheno<-'empty'
lines<-c(358,380,441,486,832,913) #DGRP line numbers of samples
pheno<-c('LL','LL','HL','LL','HL','HL') #line phenotypes
for (i in 1:6) {
    sampleInfo$pheno[sampleInfo$line==lines[i]]<-pheno[i]
}

getIds<-function(x) {
    as.numeric(unlist(strsplit(x,'_'))[2])
}

sampleIds<-sapply(names(table)[-1],getIds)
tableL.dge<-DGEList(counts=table[-1][,sampleIds <= 72],genes=table[,(1)])

# Check
tableL.dge$samples
nrow(tableL.dge$samples)

tableA.dge<-DGEList(counts=table[-1][,sampleIds > 72],genes=table[,(1)])
# Check
tableA.dge$samples
nrow(tableA.dge$samples)

# Check structure of object
structure(tableL.dge)
structure(tableA.dge)

## TMM normalization
tableL.dge<- calcNormFactors(tableL.dge)
tableL.dge$samples 
tableA.dge<- calcNormFactors(tableA.dge)
tableA.dge$samples 

## Creation of MDS Plots
#pdf('Plots/MDS_Larvae_AllSamples.pdf')
plotMDS(tableL.dge)
#dev.off()
#pdf('Plots/MDS_Adults_AllSamples.pdf')
plotMDS(tableA.dge)
#dev.off()

#output library sizes for all samples
setwd('/media/raglandlab/ExtraDrive1/dmel/RNAseqStudy/PhilRnaSeqAnalysis/')
out<-rbind(tableA.dge$samples,tableL.dge$samples)
out$lib.size=round(out$lib.size)
write.table(out,'librarySizes.txt',sep="\t",row.names=T,quote=F)
median(out$lib.size) #4,808,878


######based on MDS plots and small library size, the following samples look suspect (see supplement) and are dropped:
## min libray size:
badIds<-c(2,6,13,24,56,73,80,89,96,112)
#library sizes of badIds:
tableL.dge$samples
tableL.dge$samples[ rownames(tableL.dge$samples) %in% names(sampleIds)[sampleIds %in% badIds] , ]
tableA.dge$samples[ rownames(tableA.dge$samples) %in% names(sampleIds)[sampleIds %in% badIds] , ]
mean( tableL.dge$samples[ !(rownames(tableL.dge$samples) %in% names(sampleIds)[sampleIds %in% badIds] ), 2 ] )

a<-table[-1][,!(sampleIds %in% badIds)]
table<-data.frame(gene_id=table[,1],a)
sampleIds<-sapply(names(table)[-1],getIds)
tableL.dge<-DGEList(counts=table[-1][,sampleIds <= 72],genes=table[,(1)])
tableA.dge<-DGEList(counts=table[-1][,sampleIds > 72],genes=table[,(1)])
tableAll.dge<-DGEList(counts=table[-1],genes=table[,(1)])

#re-normalize after dropping samples
tableL.dge<-calcNormFactors(tableL.dge)
tableA.dge<-calcNormFactors(tableA.dge)
tableAll.dge<-calcNormFactors(tableAll.dge)

#write.table(tableAll.dge$genes,'EdgeRFlybaseIds.txt',quote=F,row.names=F) # used to map flybase ids to gene symbols in flybase

#preliminary check for outliers (105 and 108 a ways out, but not too egregious, no clear experimental or biological reason to exclude)
#plotMDS(tableAll.dge)
#plotMDS(tableL.dge)
#plotMDS(tableA.dge)

sampleInfo<-sampleInfo[!(sampleInfo$sample_id %in% badIds),]
Info_L<-sampleInfo[sampleInfo$stage=='L',]
Info_A<-sampleInfo[sampleInfo$stage=='A',]
Info_All<-sampleInfo

# Fix column headers - change first column names to match the DGE table
Info_L$sample_id<-paste("droso_", Info_L$sample_id, "_exp_count", sep="")
Info_A$sample_id<-paste("droso_", Info_A$sample_id, "_exp_count", sep="")
Info_All$sample_id<-paste("droso_", Info_All$sample_id, "_exp_count", sep="")

# Matching the order of the DGE list
Info_L<-Info_L[order(Info_L$sample_id),]
Info_A<-Info_A[order(Info_A$sample_id),]
Info_All<-Info_All[order(Info_All$sample_id),]

# Checks to make sure - all TRUES - If we got any Falses - there is a sample out of order between the DGE list and the MDS file
Info_L$sample_id==rownames(tableL.dge$samples)
Info_A$sample_id==rownames(tableA.dge$samples)
Info_All$sample_id==rownames(tableAll.dge$samples)

############################### MODEL CREATION ###################################

############## Both Stages ###############

## Model prep: Include strings so they become factors in the model so you can see everything is in the right place
Info_All$line<-as.factor(Info_All$line) # lines need to be factors and not numbers - this is why we run with
Info_All$treatment<-as.factor(Info_All$treatment)
Info_All$pheno<-as.factor(Info_All$pheno)
designAll<-model.matrix(~treatment*line*stage,data=Info_All)

#added 11/9/20 -- list number of reps retained after removal of outliers
out<-data.frame(Info_All %>% count(interaction(line,treatment,stage)) )
#write.table(out,'nRepsAfterOutlierRemoval.txt',quote=F,row.names=F,sep="\t")

# Change rownames - add ronames to the design so we can see what the sample is called
rownames(designAll)<-colnames(tableAll.dge)
# Check - this is the original model where each line is treated separately
head(designAll)

## Dispersion (Full model) - how much of the variation is random vs. how much of it is biologically meaningful
tableAll.dge<- estimateDisp(tableAll.dge, designAll, robust=TRUE) #tagwise = TRUE for each gene separately
#tableL.dge<- estimateDisp(tableL.dge, design_test, robust=TRUE)
# check common dispersion - how much is actually random (take into account, the variability of the data)
tableAll.dge$common.dispersion
#0.05006762
sqrt(0.05006762) # take the sqrt to get the proportion (coefficient) of variation
#0.234 coefficent of variation - 23% variation among replicates


## Fitting a general linear model to data including all paramaters for stage, line, time, and all of their interactions
fitAll <- glmFit(tableAll.dge, designAll)
# Check the colnames
colnames(designAll)



################################################################################################
####### Tests for differences in expression between HA and HL phenotypes               #########
################################################################################################


colnames(designAll)
#[1] "(Intercept)"                 "treatment30"                 "treatment30a"                "treatment60"                 "line380"                     "line441"                    
#[7] "line486"                     "line832"                     "line913"                     "stageL"                      "treatment30:line380"         "treatment30a:line380"       
#[13] "treatment60:line380"         "treatment30:line441"         "treatment30a:line441"        "treatment60:line441"         "treatment30:line486"         "treatment30a:line486"       
#[19] "treatment60:line486"         "treatment30:line832"         "treatment30a:line832"        "treatment60:line832"         "treatment30:line913"         "treatment30a:line913"       
#[25] "treatment60:line913"         "treatment30:stageL"          "treatment30a:stageL"         "treatment60:stageL"          "line380:stageL"              "line441:stageL"             
#[31] "line486:stageL"              "line832:stageL"              "line913:stageL"              "treatment30:line380:stageL"  "treatment30a:line380:stageL" "treatment60:line380:stageL" 
#[37] "treatment30:line441:stageL"  "treatment30a:line441:stageL" "treatment60:line441:stageL"  "treatment30:line486:stageL"  "treatment30a:line486:stageL" "treatment60:line486:stageL" 
#[43] "treatment30:line832:stageL"  "treatment30a:line832:stageL" "treatment60:line832:stageL"  "treatment30:line913:stageL"  "treatment30a:line913:stageL" "treatment60:line913:stageL" 

source('EdgeR_LineTrajectoryContrasts.r') # to get data frame of individual trajectories for each line in each stage in each pheno, long format for ggplot
head(LineTrajectories)

#AHAt30a - AHAt0
con1  <-c(0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        (1/3),
          0,                             0,                            0,                             0,                               0,                        (1/3),
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0 )
#AHLt30a - AHLt0
con2  <-c(0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            (1/3),                         0,                               0,                        0,
          0,                             0,                            (1/3),                         0,                               0,                        (1/3),
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0 )
con<-con1-con2
#Test change from 0 to 30a in adults in HA lines to HL lines
lrt <- glmLRT(fitAll, contrast=con)
AHA30av0_AHL30av0<-as.data.frame(topTags(lrt,n=nrow(tableL.dge),sort.by=NULL))
sum(AHA30av0_AHL30av0$FDR < 0.05) #1
#So, 1 gene with a different pattern of DE during recovery in adults of HA compared to HL lines

genes<-as.character(AHA30av0_AHL30av0$genes[AHA30av0_AHL30av0$FDR < 0.05])
ind<-LineTrajectories$gene==genes[1] #FBgn0002563
lind<-LineTrajectories$stage=='Adult'
#pdf('PlotsComparingPhenos/AHA30av0_AHL30av0_g1.pdf')
ggplot(LineTrajectories[ind & lind,], aes(x = time, y = logFC)) + 
  geom_point(aes(color = pheno, shape = line) ) + geom_line(aes(color=pheno, linetype= line)) + ggtitle(paste(genes[1]))
#dev.off()
#but, plot shows that the pattern is entirely driven by one line

#LHLt30 - LHLt0
con1  <-c(0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                            (1/3),                         0,                             0,                               0,                        0,
          0,                            (1/3),                         0,                             0,                             (1/3),                      0,
          0,                             0,                            0,                             0,                               0,                        (1/3),
          0,                            (1/3),                         (1/3),                         0,                               0,                        0,
          (1/3),                         0,                            0,                             0,                               0,                        0,
          (1/3),                         0,                            0,                             (1/3),                           0,                        0 )
#LHAt30 - LHAt0
con2  <-c(0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                             (1/3),                      0,
          0,                             0,                            0,                             0,                             (1/3),                      0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                             (1/3),                      0,
          (1/3),                         0,                            0,                             (1/3),                           0,                        0,
          0,                             0,                            0,                             (1/3),                           0,                        0,
          0,                             0,                            0,                             0,                               0,                        0 )
con<-con1-con2
#Test change from 0 to 30m in larvae in HA lines to HL lines
lrt <- glmLRT(fitAll, contrast=con)
LHL30v0_LHA30v0<-as.data.frame(topTags(lrt,n=nrow(tableL.dge),sort.by=NULL))
sum(LHL30v0_LHA30v0$FDR < 0.05) #329
#So, 329 gene with a different pattern of DE during recovery in adults of HA compared to HL lines

genes<-as.character(LHL30v0_LHA30v0$genes[LHL30v0_LHA30v0$FDR < 0.05])
ind<-LineTrajectories$gene==genes[1]
lind<-LineTrajectories$stage=='Adult'
#pdf('PlotsComparingPhenos/LHL30v0_LHA30v0_g1.pdf')
ggplot(LineTrajectories[ind & lind,], aes(x = time, y = logFC)) + 
  geom_point(aes(color = pheno, shape = line) ) + geom_line(aes(color=pheno, linetype= line)) + ggtitle(paste(genes[1]))
#dev.off()
#but, once again, patterns for an arbitrary gene show that significance is driven mainly by one outlier line

################################################################################################
################################################################################################
################################################################################################


################################################################################################
####### Build reduced models that remove transcripts with complex interactions         #########
################################################################################################


#Test 3-way interaction
lrt <- glmLRT(fitAll, coef=34:48)
interactions<-as.data.frame(topTags(lrt,n=nrow(tableL.dge),sort.by=NULL))


#remove all transcripts with significant 3-way interaction
sum(interactions$FDR < 0.05) #130
ind<-interactions$FDR < 0.05
tableAll.reduced.dge<-tableAll.dge[!ind,]

#re-estimate dispersion with reduced model
designAll<-model.matrix(~stage+treatment+line+stage:treatment+
                           stage:line+treatment:line,data=Info_All)
rownames(designAll)<-colnames(tableAll.reduced.dge)
tableAll.reduced.dge<- estimateDisp(tableAll.reduced.dge, designAll, robust=TRUE) #tagwise = TRUE for each gene separately
tableAll.reduced.dge$common.dispersion

fitAll <- glmFit(tableAll.reduced.dge, designAll)
# Check the colnames
colnames(designAll) 

#Test treatment:line interaction
lrt <- glmLRT(fitAll, coef=19:33)
interactions<-as.data.frame(topTags(lrt,n=nrow(tableL.dge),sort.by=NULL))
sum(interactions$FDR < 0.05) #19
ind<-interactions$FDR < 0.05


#remove all transcripts with significant treatment:line interaction
tableAll.reduced.dge<-tableAll.reduced.dge[!ind,]

#re-estimate dispersion with reduced model
designAll<-model.matrix(~stage+treatment+line+stage:treatment+
                          stage:line,data=Info_All)
rownames(designAll)<-colnames(tableAll.reduced.dge)
tableAll.reduced.dge<- estimateDisp(tableAll.reduced.dge, designAll, robust=TRUE) #tagwise = TRUE for each gene separately
tableAll.reduced.dge$common.dispersion

fitAll <- glmFit(tableAll.reduced.dge, designAll)
# Check the colnames
colnames(designAll)

################################################################################################
################################################################################################
################################################################################################

################################################################################################
####### Build tables testing various final (reduced) model effects         #####################
################################################################################################

#Test stage:line interaction
lrt <- glmLRT(fitAll, coef=14:18)
stageByLine<-as.data.frame(topTags(lrt,n=nrow(tableL.dge),sort.by=NULL))
sum(stageByLine$FDR < 0.05) #7549

#Test stage:treatment interaction
lrt <- glmLRT(fitAll, coef=11:13)
stageByTreat<-as.data.frame(topTags(lrt,n=nrow(tableL.dge),sort.by=NULL))
sum(stageByTreat$FDR < 0.05) #880


#Test stage main effect
lrt <- glmLRT(fitAll, coef=2)
stage<-as.data.frame(topTags(lrt,n=nrow(tableL.dge),sort.by=NULL))
sum(stage$FDR < 0.05) #10524
sum(stage$FDR < 0.05 & stageByTreat$FDR > 0.05) #9780
ind<-stage$FDR < 0.05 & stageByTreat$FDR > 0.05
#write.table(stage[ind,],"Tables/StageEffectNoStageByTreat.txt",row.names=F,quote=F,sep="\t")


#Test treatment main effect
lrt <- glmLRT(fitAll, coef=3:5)
treatment<-as.data.frame(topTags(lrt,n=nrow(tableL.dge),sort.by=NULL))
sum(treatment$FDR < 0.05) #221

#stage main effects only (includes stageByLine)
sum(stage$FDR < 0.05  & stageByTreat$FDR > 0.05 ) #9780

#stage main effects or stage:line only 
sum( (stage$FDR < 0.05 | stageByLine$FDR < 0.05 )  & stageByTreat$FDR > 0.05 ) #10966

########## 4/14/2021 - check for enrichment of tissue-specific transcripts in set showing stage-specific effects

source('tissueSpecificityInFredaEtAlData.R')
# punchline1: transcripts with (stage or stageXline) effects have LOWER than expected tissue specificity
# punchline2: transcripts with (stageXtreatment) effects also have LOWER than expected tissue specificity

#treat main effects only
sum(treatment$FDR < 0.05 & stageByTreat$FDR > 0.05 ) #21

#Test stage:treatment interaction
lrt <- glmLRT(fitAll, coef=11:13)
stageByTreat<-as.data.frame(topTags(lrt,n=nrow(tableL.dge),sort.by=NULL))
sum(stageByTreat$FDR < 0.05) #880

#Create barplot of log counts of genes with stage, time, and stageXtime effects
barPlotDat<-data.frame(cat=c('stage','time','stageXtime'),logCount=log10(c(10966,21,880)))
barPlotDat$cat<-factor(barPlotDat$cat,levels=c('stage','time','stageXtime'))
#pdf('Plots/log10CountsDE_transcripts.pdf')
ggplot(data=barPlotDat, aes(x=cat, y=logCount)) +
  geom_bar(stat="identity") + theme_classic()
#dev.off()

################################################################################################
################################################################################################
################################################################################################

################################################################################################
####### Create plots of time trajectories for various sets of genes        #####################
################################################################################################

colnames(designAll)
#[1] "(Intercept)"         "stageL"              "treatment30"         "treatment30a"       
#[5] "treatment60"         "line380"             "line441"             "line486"            
#[9] "line832"             "line913"             "stageL:treatment30"  "stageL:treatment30a"
#[13] "stageL:treatment60"  "stageL:line380"      "stageL:line441"      "stageL:line486"     
#[17] "stageL:line832"      "stageL:line913"     

#build temperature trajectories per stage .... -> StageTrajectories dataframe and StageByTempMatrix matrix
source('EdgeR_StageByTimeTrajectories.r')

treatGenes<-as.character(treatment$genes[treatment$FDR < 0.05 & stageByLine$FDR > 0.05 & stageByTreat$FDR > 0.05])

ind<-StageTrajectories$gene %in% treatGenes
plotDat<-StageTrajectories[ind,]
plotDat$sym<-sapply(plotDat$gene,getSym)
plotDat<-plotDat %>%
  group_by(time,sym,gene) %>%
  summarize(meanLFC = mean(logFC, na.rm = TRUE))
plotDat<-data.frame(plotDat)

#pdf('Plots/TimeMainEffectTrajectories.pdf')
ggplot(plotDat, aes(x = time, y = meanLFC)) + scale_linetype_manual(values=rep(1:6,4)) + theme(text=element_text(family="Arial")) + scale_x_continuous(breaks=c(0,30,60,90)) +
  geom_line(aes(color=sym,linetype=sym))  + ggtitle('TimeMainEffect, mean across stages') + theme_classic() + theme(axis.text = element_text(size=20),legend.text=element_text(size=12)) + 
  labs(x='time (min)',y='Relative Expression') + geom_hline(yintercept = 0,lty=3)
#dev.off()







 
# From StageByTimeTrajectories script...
ind<-stageByTreat$FDR < 0.05
inda<-A30v0$FDR < 0.05 | A60v0$FDR < 0.05 | A90v0$FDR < 0.05
indl<-L30v0$FDR < 0.05 | L60v0$FDR < 0.05 | L90v0$FDR < 0.05
#number of genes DE across one time point in larvae and with stageXtime interaction
sum(ind & indl) #763
#number of genes DE across one time point in adults and with stageXtime interaction
sum(ind & inda) #121
out<-stageByTreat
out$sigStage<-NA
out$sigStage[ind &indl]<-'larva'
out$sigStage[ind & inda]<-'adult'
out<-out[!is.na(out$sigStage),]
#write.table(out,'Tables/sigTimeInStageAndSigStageByTime.txt',row.names=F,quote=F,sep="\t")


expMat<-expMat[ind,]

ids<-as.character(stageByTreat$genes[ind])
plotDat<-StageTrajectories[StageTrajectories$gene %in% ids,]
plotDat$Sig='none'
ids<-as.character(stageByTreat$genes[ind & inda & !indl])
plotDat$Sig[plotDat$gene %in% ids]<-'Adult'
ids<-as.character(stageByTreat$genes[ind & indl & !inda])
plotDat$Sig[plotDat$gene %in% ids]<-'Larva'
ids<-as.character(stageByTreat$genes[ind & inda & indl])
plotDat$Sig[plotDat$gene %in% ids]<-'LarvaAdult'
plotDat$direction='none'



#source('ClusterStageByTempTranscripts.r')
a<-data.frame( plotDat %>%
  group_by(stage, gene) %>%
  summarize(meanLogFC = mean(logFC, na.rm = TRUE)) )
a$direction='up'
a$direction[a$meanLogFC < 0]<-'down'


idToDirLarva <- setNames(as.list(  as.character(a$direction[a$stage=='Larva']) ),as.character(a$gene[a$stage=='Larva']))
idToDirAdult <- setNames(as.list(  as.character(a$direction[a$stage=='Adult']) ),as.character(a$gene[a$stage=='Adult']))   
#getDir<-function(id) {
#  Dir<-idToDir[[as.character(id)]]
#  return(Dir)
#}
#for (i in 1:nrow(plotDat)) {
#  if (plotDat$stage[i]=='Larva') {
#    plotDat$direction[i]<-unlist(idToDirLarva[[as.character(plotDat$gene[i])]])
#  } else {
#    plotDat$direction[i]<-unlist(idToDirAdult[[as.character(plotDat$gene[i])]])
#  } 
#}

a<-plotDat[plotDat$time==0,]
a$direction<-'down'
plotDat<-rbind(a,plotDat)

plotDatL<-plotDat
plotDatA<-plotDat
for (i in 1:nrow(plotDatL)) {
    plotDatL$direction[i]<-unlist(idToDirLarva[[as.character(plotDatL$gene[i])]])
}
for (i in 1:nrow(plotDatA)) {
  plotDatA$direction[i]<-unlist(idToDirAdult[[as.character(plotDatA$gene[i])]])
}



sePlus<-function(x) {
  if (length(x)==1) {x<-c(x,x)}
  mean(x)+2*sd(x)/sqrt(length(x))
}
seMinus<-function(x) {
  if (length(x)==1) {x<-c(x,x)}
  mean(x)-2*sd(x)/sqrt(length(x))
}

ind<-grepl('Larva',plotDat$Sig)
#pdf('Plots/TractoriesStageByTreatSigDE_time_Larvae.pdf')
ggplot(plotDatL[ind,], aes(x=time,y=logFC, color=interaction(stage,direction))) +
  stat_summary(geom="ribbon", fun.ymin="seMinus", fun.ymax="sePlus", aes(fill=interaction(stage,direction)),  alpha=0.3) +
  theme_classic()+ggtitle('TractoriesStageByTreatSigDE_time_Larvae') + coord_cartesian(ylim = c(-1, 1.5)) +
  scale_x_continuous(breaks=c(0,30,60,90)) +
  theme(axis.text = element_text(size=20)) + geom_hline(yintercept = 0,lty=3)
#dev.off()
#pdf('Plots/TractoriesStageByTreatSigDE_time_Adult.pdf')
ind<-grepl('Adult',plotDat$Sig)
ggplot(plotDatA[ind,], aes(x=time,y=logFC, color=interaction(stage,direction))) +
  stat_summary(geom="ribbon", fun.ymin="seMinus", fun.ymax="sePlus", aes(fill=interaction(stage,direction)),  alpha=0.3) +
  theme_classic()+ggtitle('TractoriesStageByTreatSigDE_time_Adult') + coord_cartesian(ylim = c(-1, 1.5)) +
  scale_x_continuous(breaks=c(0,30,60,90)) +
  theme(axis.text = element_text(size=20)) + geom_hline(yintercept = 0,lty=3)
#dev.off()
##note; very small sample for transcripts sig DE down in adults, will remove that cat from the adult plot

#Fig for supplement including 5 downreg transcripts in adults
#pdf('Plots/TractoriesStageByTreatSigDE_time_AdultFullScale.pdf')
ind<-grepl('Adult',plotDat$Sig)
ggplot(plotDatA[ind,], aes(x=time,y=logFC, color=interaction(stage,direction))) +
  stat_summary(geom="ribbon", fun.ymin="seMinus", fun.ymax="sePlus", aes(fill=interaction(stage,direction)),  alpha=0.3) +
  theme_classic()+ggtitle('TractoriesStageByTreatSigDE_time_Adult')
#dev.off()

#number transcripts in each category:
ind<-grepl('Larva',plotDat$Sig)
length(unique(plotDatL$gene[ind & plotDatL$direction=='up'])) #516 sig up larva
length(unique(plotDatL$gene[ind & plotDatL$direction=='down'])) #247 sig down larva
sigDEtimeLarvae<-unique(plotDatL$gene[ind])
#write.table(data.frame(sigDEtimeLarvae),"Tables/sigDEtimeLarveGenes.txt",sep="\t",row.names=F,quote=F)
ind<-grepl('Adult',plotDat$Sig)
length(unique(plotDatA$gene[ind & plotDatA$direction=='up'])) #116 sig up Adult
length(unique(plotDatA$gene[ind & plotDatA$direction=='down'])) #5 sig down Adult
sigDEtimeAdult<-unique(plotDatA$gene[ind])
#write.table(data.frame(sigDEtimeAdult),"Tables/sigDEtimeAdultGenes.txt",sep="\t",row.names=F,quote=F)

getSym<-function(id) {
  flysim<-idToSym[[as.character(id)]]
  return(flysim)
}

plotDat$symbol<-sapply(as.character(plotDat$gene),getSym)
syms<-c(0:4,15:18,25)
#from https://www.datanovia.com/en/blog/ggplot-colors-best-tricks-you-will-love/
cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
library(RColorBrewer)
#write.table(plotDat,'Rdatasets/plotDatStageByTime.txt',quote=F,row.names=F,sep="\t")

#hsps UP_KEYWORDS Stress Response, All FDR = 0.048, Adult FDR = 7.1E-7, Larva FDR = 1.7E-1
ppcolors<-c('#B84025','#00A29F')
ids<-qw("FBgn0263106	FBgn0001228	FBgn0001223	FBgn0001224	FBgn0001230	FBgn0013276	FBgn0013278	FBgn0013279	FBgn0000256")
#ids<-c('FBgn0001223','FBgn0001224','FBgn0001225','FBgn0001230','FBgn0013276','FBgn0013278','FBgn0051354','FBgn0013279','FBgn0263106')
#pdf('Plots/StressResponse.pdf')
subDat<-plotDat[plotDat$gene %in% ids,]
ggplot(data=subDat,aes(x=time,y=logFC,group=interaction(symbol,stage),col=symbol,lty=stage)) + theme_classic() + geom_line(lwd=2) + ggtitle('Heat Shock Response') +
  geom_hline(yintercept = 0,lty=3) + scale_linetype_manual(values=c("dashed", "solid"))+  scale_x_continuous(breaks=c(0,30,60,90)) +
  theme(legend.position=c(.1, .7)) + theme(axis.text = element_text(size=20),legend.text=element_text(size=13)) + labs(x='time (min)',y='Relative Expression') 
#dev.off()

#  	GOTERM_BP_DIRECT autophagy, All FDR = 2.1E-4, Adult FDR = ~1, Larva FDR =  1.5E-4
#ids<-c('FBgn0260945','FBgn0030960','FBgn0037363','FBgn0044452','FBgn0034366','FBgn0052672','FBgn0264490','FBgn0015277','FBgn0003943','FBgn0260935','FBgn0003997','FBgn0011706')
ids<-qw("FBgn0260945	FBgn0030960	FBgn0037363	FBgn0044452	FBgn0052672	FBgn0015279	FBgn0038129	FBgn0264490	FBgn0004387	FBgn0015277	FBgn0264357	FBgn0262516	FBgn0003943	FBgn0260935	FBgn0028582	FBgn0086656")
#pdf('Plots/autophagy.pdf')
subDat<-plotDat[plotDat$gene %in% ids,]
ggplot(data=subDat,aes(x=time,y=logFC,group=interaction(symbol,stage),col=symbol,lty=stage)) + theme_classic() + geom_line(lwd=2) + ggtitle('Autophagy') +
  geom_hline(yintercept = 0,lty=3) + scale_linetype_manual(values=c("dashed", "solid"))+  scale_x_continuous(breaks=c(0,30,60,90)) +
  theme(legend.position=c(.2, 0.3)) + theme(axis.text = element_text(size=20),legend.text=element_text(size=13)) + labs(x='time (min)',y='Relative Expression') 
#dev.off()

#GOTERM_BP_DIRECT 	response to bacterium, All FDR = 7.7E-3, Adult FDR = 1.8E-3, Larva FDR = 3.9E-2 
#ids<-c('FBgn0028484',	'FBgn0012042',	'FBgn0041581',	'FBgn0041579',	'FBgn0262081',	'FBgn0053470',	'FBgn0026760',	'FBgn0043578',	'FBgn0000276',	'FBgn0000277',	'FBgn0000278',	'FBgn0015247',	'FBgn0004240',	'FBgn0011274',	'FBgn0283461',	'FBgn0040322',	'FBgn0034329',	'FBgn0067905',	'FBgn0067903',	'FBgn0025583',	'FBgn0034328',	'FBgn0040736',	'FBgn0032006',	'FBgn0014018',	'FBgn0031975',	'FBgn0010441'
#       )
ids<-qw("FBgn0041579	FBgn0261560	FBgn0000276	FBgn0000277	FBgn0004240	FBgn0010388	FBgn0034329	FBgn0025583	FBgn0014865	FBgn0014018	FBgn0039593	FBgn0267339")
#pdf('Plots/ResponseToBacterium.pdf')
subDat<-plotDat[plotDat$gene %in% ids,]
ggplot(data=subDat,aes(x=time,y=logFC,group=interaction(symbol,stage),col=symbol,lty=stage)) + theme_classic() + geom_line(lwd=2) + ggtitle('ResponseToBacterium') +
  geom_hline(yintercept = 0,lty=3) + scale_linetype_manual(values=c("dashed", "solid"))+  scale_x_continuous(breaks=c(0,30,60,90)) +
  theme(legend.position=c(.1, .7)) + theme(axis.text = element_text(size=20),legend.text=element_text(size=13)) + labs(x='time (min)',y='Relative Expression') 
#dev.off()

#response to bacterium I
#ids<-c('FBgn0041579','FBgn0000277','FBgn0004240','FBgn0000276','FBgn0014865','FBgn0010388')
#subDat<-plotDat[plotDat$gene %in% ids,]
#ggplot(data=subDat,aes(x=time,y=logFC,group=interaction(symbol,stage),col=symbol,lty=stage)) + theme_classic() + geom_line(lwd=2) + ggtitle('') +
#  geom_hline(yintercept = 0,lty=3) + scale_linetype_manual(values=c("dashed", "solid"))+  scale_x_continuous(breaks=c(0,30,60,90)) +
#  theme(legend.position=c(.1, .7)) + theme(axis.text = element_text(size=20),legend.text=element_text(size=13)) + labs(x='time (min)',y='Relative Expression') 

#Innate Immunity I
#ids<-c('FBgn0041579',	'FBgn0000276',	'FBgn0000277',	'FBgn0004240')
#subDat<-plotDat[plotDat$gene %in% ids,]
#ggplot(data=subDat,aes(x=time,y=logFC,group=interaction(symbol,stage),col=symbol,lty=stage)) + theme_classic() + geom_line(lwd=2) + ggtitle('Innate Immunity I') +
#  geom_hline(yintercept = 0,lty=3) + scale_linetype_manual(values=c("dashed", "solid"))+ theme(text=element_text(family="Arial")) + scale_x_continuous(breaks=c(0,30,60,90)) +
#  theme(legend.position=c(.1, .7)) + theme(axis.text = element_text(size=22),legend.text=element_text(size=12)) + labs(x='time (min)',y='Relative Expression') 
 
#Innate Immunity II
#ids<-c('FBgn0028484',	'FBgn0012042',	'FBgn0041581',	'FBgn0262081',	'FBgn0053470',	'FBgn0043578',	'FBgn0000278',	'FBgn0015247',	'FBgn0011274',	'FBgn0283461',	'FBgn0040322',	'FBgn0034329',	'FBgn0067905',	'FBgn0067903',	'FBgn0025583',	'FBgn0034328',	'FBgn0040736',	'FBgn0032006',	'FBgn0014018',	'FBgn0031975',	'FBgn0010441'
#)
#subDat<-plotDat[plotDat$gene %in% ids,]
#ggplot(data=subDat,aes(x=time,y=logFC,group=interaction(symbol,stage),col=symbol,lty=stage)) + theme_classic() + geom_line(lwd=2) + ggtitle('Innate Immunity II') +
#  geom_hline(yintercept = 0,lty=3) + scale_linetype_manual(values=c("dashed", "solid"))+ theme(text=element_text(family="Arial")) + scale_x_continuous(breaks=c(0,30,60,90)) +
#  theme(legend.position=c(.1, .7)) + theme(axis.text = element_text(size=22),legend.text=element_text(size=12)) + labs(x='time (min)',y='Relative Expression') 

# 	INTERPRO Basic-leucine zipper domain , All FDR = 1.7E-2, Adult FDR = ~1, Larva FDR =  6.1E-3
ids<-c('FBgn0034534',	'FBgn0265784',	'FBgn0001291',	'FBgn0038063',	'FBgn0016694',	'FBgn0032202',	'FBgn0039209',	'FBgn0005638',	'FBgn0016076')
pdf('Plots/BasicLeucineZipper.pdf')
subDat<-plotDat[plotDat$gene %in% ids,]
ggplot(data=subDat,aes(x=time,y=logFC,group=interaction(symbol,stage),col=symbol,lty=stage)) + theme_classic() + geom_line(lwd=2) + ggtitle('Basic Leucine Zipper') +
  geom_hline(yintercept = 0,lty=3) + scale_linetype_manual(values=c("dashed", "solid"))+  scale_x_continuous(breaks=c(0,30,60,90)) +
  theme(legend.position=c(.1, .7)) + theme(axis.text = element_text(size=20),legend.text=element_text(size=13)) + labs(x='time (min)',y='Relative Expression') 
dev.off()

#KEGG_PATHWAY 	Fatty acid metabolism, , All FDR = 1.1E-4, Adult FDR = ~1, Larva FDR =  7.3E-5
ids<-qw("FBgn0033246	FBgn0263120	FBgn0035811	FBgn0029975	FBgn0036824	FBgn0027572	FBgn0038130	FBgn0031813	FBgn0039756	FBgn0086687	FBgn0032358	FBgn0034629	FBgn0027348	FBgn0261862")
pdf('Plots/FattyAcidMetabolism.pdf')
subDat<-plotDat[plotDat$gene %in% ids,]
ggplot(data=subDat,aes(x=time,y=logFC,group=interaction(symbol,stage),col=symbol,lty=stage)) + theme_classic() + geom_line(lwd=2) + ggtitle('Fatty acid metabolism') +
  geom_hline(yintercept = 0,lty=3) + scale_linetype_manual(values=c("dashed", "solid"))+  scale_x_continuous(breaks=c(0,30,60,90)) +
  theme(legend.position=c(.1, .7)) + theme(axis.text = element_text(size=20),legend.text=element_text(size=13)) + labs(x='time (min)',y='Relative Expression') 
dev.off()

######## breakpoint 10/16/2020
save.image(file='intermediateSave.Rdat')
#load it
#setwd("/media/raglandlab/ExtraDrive1/dmel/RNAseqStudy/PhilRnaSeqAnalysis")
#load('intermediateSave.Rdat')





