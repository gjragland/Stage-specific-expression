#R
#GJR, last modified 10/16/2020
#companion script to EdgeR_modelsAndTrajectories.r
#Run EdgeR contrasts to compare HL and LL phenotypes

#[1] "(Intercept)"                 "treatment30"                 "treatment30a"                "treatment60"                 "line380"                     "line441"                    
#[7] "line486"                     "line832"                     "line913"                     "stageL"                      "treatment30:line380"         "treatment30a:line380"       
#[13] "treatment60:line380"         "treatment30:line441"         "treatment30a:line441"        "treatment60:line441"         "treatment30:line486"         "treatment30a:line486"       
#[19] "treatment60:line486"         "treatment30:line832"         "treatment30a:line832"        "treatment60:line832"         "treatment30:line913"         "treatment30a:line913"       
#[25] "treatment60:line913"         "treatment30:stageL"          "treatment30a:stageL"         "treatment60:stageL"          "line380:stageL"              "line441:stageL"             
#[31] "line486:stageL"              "line832:stageL"              "line913:stageL"              "treatment30:line380:stageL"  "treatment30a:line380:stageL" "treatment60:line380:stageL" 
#[37] "treatment30:line441:stageL"  "treatment30a:line441:stageL" "treatment60:line441:stageL"  "treatment30:line486:stageL"  "treatment30a:line486:stageL" "treatment60:line486:stageL" 
#[43] "treatment30:line832:stageL"  "treatment30a:line832:stageL" "treatment60:line832:stageL"  "treatment30:line913:stageL"  "treatment30a:line913:stageL" "treatment60:line913:stageL" 



####### HA ############################

#A358t30 - A358t0
con   <-c(0,                             1,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0 )
lrt <- glmLRT(fitAll, contrast=con)
t30v0<-as.data.frame(topTags(lrt,n=nrow(tableL.dge),sort.by=NULL))

#A358t60 - A358t0
con   <-c(0,                             0,                            0,                             1,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0 )
lrt <- glmLRT(fitAll, contrast=con)
t60v0<-as.data.frame(topTags(lrt,n=nrow(tableL.dge),sort.by=NULL))

#A358t30a - A358t0
con   <-c(0,                             0,                            1,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0 )
lrt <- glmLRT(fitAll, contrast=con)
t90v0<-as.data.frame(topTags(lrt,n=nrow(tableL.dge),sort.by=NULL))

LineTrajectories<-data.frame(stage='Adult',line='L358',pheno='HA',time=30,gene=t30v0$genes,logFC=t30v0$logFC,stringsAsFactors = F)
t0<-LineTrajectories
t0$time=0
t0$logFC=0
LineTrajectories<-rbind(t0,LineTrajectories)
LineTrajectories<-rbind(LineTrajectories,data.frame( stage='Adult',line='L358',pheno='HA',time=60,gene=t60v0$genes,logFC=t60v0$logFC ) )
LineTrajectories<-rbind(LineTrajectories,data.frame( stage='Adult',line='L358',pheno='HA',time=90,gene=t90v0$genes,logFC=t90v0$logFC ) )



#A380t30 - A380t0
con   <-c(0,                             1,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               1,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0 )
lrt <- glmLRT(fitAll, contrast=con)
t30v0<-as.data.frame(topTags(lrt,n=nrow(tableL.dge),sort.by=NULL))

#A380t60 - A380t0
con   <-c(0,                             0,                            0,                             1,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          1,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0 )
lrt <- glmLRT(fitAll, contrast=con)
t60v0<-as.data.frame(topTags(lrt,n=nrow(tableL.dge),sort.by=NULL))

#A380t30a - A380t0
con   <-c(0,                             0,                            1,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        1,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0 )
lrt <- glmLRT(fitAll, contrast=con)
t90v0<-as.data.frame(topTags(lrt,n=nrow(tableL.dge),sort.by=NULL))

t0$line='L380'
LineTrajectories<-rbind(LineTrajectories,t0)
LineTrajectories<-rbind(LineTrajectories,data.frame(stage='Adult',line='L380',pheno='HA',time=30,gene=t30v0$genes,logFC=t30v0$logFC  ) )
LineTrajectories<-rbind(LineTrajectories,data.frame( stage='Adult',line='L380',pheno='HA',time=60,gene=t60v0$genes,logFC=t60v0$logFC ) )
LineTrajectories<-rbind(LineTrajectories,data.frame( stage='Adult',line='L380',pheno='HA',time=90,gene=t90v0$genes,logFC=t90v0$logFC ) )



#A486t30 - A486t0
con   <-c(0,                             1,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               1,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0 )
lrt <- glmLRT(fitAll, contrast=con)
t30v0<-as.data.frame(topTags(lrt,n=nrow(tableL.dge),sort.by=NULL))

#A486t60 - A486t0
con   <-c(0,                             0,                            0,                             1,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          1,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0 )
lrt <- glmLRT(fitAll, contrast=con)
t60v0<-as.data.frame(topTags(lrt,n=nrow(tableL.dge),sort.by=NULL))

#A486t30a - A486t0
con   <-c(0,                             0,                            1,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        1,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0 )
lrt <- glmLRT(fitAll, contrast=con)
t90v0<-as.data.frame(topTags(lrt,n=nrow(tableL.dge),sort.by=NULL))

t0$line='L486'
LineTrajectories<-rbind(LineTrajectories,t0)
LineTrajectories<-rbind(LineTrajectories,data.frame(stage='Adult',line='L486',pheno='HA',time=30,gene=t30v0$genes,logFC=t30v0$logFC  ) )
LineTrajectories<-rbind(LineTrajectories,data.frame( stage='Adult',line='L486',pheno='HA',time=60,gene=t60v0$genes,logFC=t60v0$logFC ) )
LineTrajectories<-rbind(LineTrajectories,data.frame( stage='Adult',line='L486',pheno='HA',time=90,gene=t90v0$genes,logFC=t90v0$logFC ) )


#HL ######################################

#A441t30 - A441t0
con   <-c(0,                             1,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             1,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0 )
lrt <- glmLRT(fitAll, contrast=con)
t30v0<-as.data.frame(topTags(lrt,n=nrow(tableL.dge),sort.by=NULL))

#A441t60 - A441t0
con   <-c(0,                             0,                            0,                             1,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             1,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0 )
lrt <- glmLRT(fitAll, contrast=con)
t60v0<-as.data.frame(topTags(lrt,n=nrow(tableL.dge),sort.by=NULL))

#A441t30a - A441t0
con   <-c(0,                             0,                            1,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            1,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0 )

lrt <- glmLRT(fitAll, contrast=con)
t90v0<-as.data.frame(topTags(lrt,n=nrow(tableL.dge),sort.by=NULL))

t0$line='L441'
t0$pheno='HL'
LineTrajectories<-rbind(LineTrajectories,t0)
LineTrajectories<-rbind(LineTrajectories,data.frame(stage='Adult',line='L441',pheno='HL',time=30,gene=t30v0$genes,logFC=t30v0$logFC  ) )
LineTrajectories<-rbind(LineTrajectories,data.frame( stage='Adult',line='L441',pheno='HL',time=60,gene=t60v0$genes,logFC=t60v0$logFC ) )
LineTrajectories<-rbind(LineTrajectories,data.frame( stage='Adult',line='L441',pheno='HL',time=90,gene=t90v0$genes,logFC=t90v0$logFC ) )


#A832t30 - A832t0
con   <-c(0,                             1,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             1,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0 )
lrt <- glmLRT(fitAll, contrast=con)
t30v0<-as.data.frame(topTags(lrt,n=nrow(tableL.dge),sort.by=NULL))

#A832t60 - A832t0
con   <-c(0,                             0,                            0,                             1,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             1,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0 )
lrt <- glmLRT(fitAll, contrast=con)
t60v0<-as.data.frame(topTags(lrt,n=nrow(tableL.dge),sort.by=NULL))

#A832t30a - A832t0
con   <-c(0,                             0,                            1,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            1,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0 )
lrt <- glmLRT(fitAll, contrast=con)
t90v0<-as.data.frame(topTags(lrt,n=nrow(tableL.dge),sort.by=NULL))

t0$line='L832'
LineTrajectories<-rbind(LineTrajectories,t0)
LineTrajectories<-rbind(LineTrajectories,data.frame(stage='Adult',line='L832',pheno='HL',time=30,gene=t30v0$genes,logFC=t30v0$logFC  ) )
LineTrajectories<-rbind(LineTrajectories,data.frame( stage='Adult',line='L832',pheno='HL',time=60,gene=t60v0$genes,logFC=t60v0$logFC ) )
LineTrajectories<-rbind(LineTrajectories,data.frame( stage='Adult',line='L832',pheno='HL',time=90,gene=t90v0$genes,logFC=t90v0$logFC ) )


#A913t30 - A913t0
con   <-c(0,                             1,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               1,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0 )
lrt <- glmLRT(fitAll, contrast=con)
t30v0<-as.data.frame(topTags(lrt,n=nrow(tableL.dge),sort.by=NULL))

#A913t60 - A913t0
con   <-c(0,                             0,                            0,                             1,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          1,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0 )
lrt <- glmLRT(fitAll, contrast=con)
t60v0<-as.data.frame(topTags(lrt,n=nrow(tableL.dge),sort.by=NULL))

#A913t30a - A913t0
con   <-c(0,                             0,                            1,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        1,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0 )
lrt <- glmLRT(fitAll, contrast=con)
t90v0<-as.data.frame(topTags(lrt,n=nrow(tableL.dge),sort.by=NULL))

t0$line='L913'
LineTrajectories<-rbind(LineTrajectories,t0)
LineTrajectories<-rbind(LineTrajectories,data.frame(stage='Adult',line='L913',pheno='HL',time=30,gene=t30v0$genes,logFC=t30v0$logFC  ) )
LineTrajectories<-rbind(LineTrajectories,data.frame( stage='Adult',line='L913',pheno='HL',time=60,gene=t60v0$genes,logFC=t60v0$logFC ) )
LineTrajectories<-rbind(LineTrajectories,data.frame( stage='Adult',line='L913',pheno='HL',time=90,gene=t90v0$genes,logFC=t90v0$logFC ) )




######### Larvae #####################


####### HA ############################

#A358t30 - A358t0
con   <-c(0,                             1,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             1,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0 )
lrt <- glmLRT(fitAll, contrast=con)
t30v0<-as.data.frame(topTags(lrt,n=nrow(tableL.dge),sort.by=NULL))

#A358t60 - A358t0
con   <-c(0,                             0,                            0,                             1,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             1,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0 )

lrt <- glmLRT(fitAll, contrast=con)
t60v0<-as.data.frame(topTags(lrt,n=nrow(tableL.dge),sort.by=NULL))

#A358t30a - A358t0
con   <-c(0,                             0,                            1,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            1,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0 )
lrt <- glmLRT(fitAll, contrast=con)
t90v0<-as.data.frame(topTags(lrt,n=nrow(tableL.dge),sort.by=NULL))


t0$stage='Larva'
t0$pheno='HA'
t0$line='L358'
LineTrajectories<-rbind(t0,LineTrajectories)
LineTrajectories<-rbind(LineTrajectories,data.frame(stage='Larva',line='L358',pheno='HA',time=30,gene=t30v0$genes,logFC=t30v0$logFC  ) )
LineTrajectories<-rbind(LineTrajectories,data.frame( stage='Larva',line='L358',pheno='HA',time=60,gene=t60v0$genes,logFC=t60v0$logFC ) )
LineTrajectories<-rbind(LineTrajectories,data.frame( stage='Larva',line='L358',pheno='HA',time=90,gene=t90v0$genes,logFC=t90v0$logFC ) )



#A380t30 - A380t0
con   <-c(0,                             1,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               1,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             1,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             1,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0 )
lrt <- glmLRT(fitAll, contrast=con)
t30v0<-as.data.frame(topTags(lrt,n=nrow(tableL.dge),sort.by=NULL))

#A380t60 - A380t0
con   <-c(0,                             0,                            0,                             1,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          1,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             1,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        1,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0 )
lrt <- glmLRT(fitAll, contrast=con)
t60v0<-as.data.frame(topTags(lrt,n=nrow(tableL.dge),sort.by=NULL))

#A380t30a - A380t0
con   <-c(0,                             0,                            1,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        1,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            1,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               1,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0 )
lrt <- glmLRT(fitAll, contrast=con)
t90v0<-as.data.frame(topTags(lrt,n=nrow(tableL.dge),sort.by=NULL))

t0$line='L380'
LineTrajectories<-rbind(LineTrajectories,t0)
LineTrajectories<-rbind(LineTrajectories,data.frame(stage='Larva',line='L380',pheno='HA',time=30,gene=t30v0$genes,logFC=t30v0$logFC  ) )
LineTrajectories<-rbind(LineTrajectories,data.frame( stage='Larva',line='L380',pheno='HA',time=60,gene=t60v0$genes,logFC=t60v0$logFC ) )
LineTrajectories<-rbind(LineTrajectories,data.frame( stage='Larva',line='L380',pheno='HA',time=90,gene=t90v0$genes,logFC=t90v0$logFC ) )



#A486t30 - A486t0
con   <-c(0,                             1,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               1,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             1,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             1,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0 )
lrt <- glmLRT(fitAll, contrast=con)
t30v0<-as.data.frame(topTags(lrt,n=nrow(tableL.dge),sort.by=NULL))

#A486t60 - A486t0
con   <-c(0,                             0,                            0,                             1,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          1,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             1,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        1,
          0,                             0,                            0,                             0,                               0,                        0 )
lrt <- glmLRT(fitAll, contrast=con)
t60v0<-as.data.frame(topTags(lrt,n=nrow(tableL.dge),sort.by=NULL))

#A486t30a - A486t0
con   <-c(0,                             0,                            1,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        1,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            1,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               1,                        0,
          0,                             0,                            0,                             0,                               0,                        0 )
lrt <- glmLRT(fitAll, contrast=con)
t90v0<-as.data.frame(topTags(lrt,n=nrow(tableL.dge),sort.by=NULL))

t0$line='L486'
LineTrajectories<-rbind(LineTrajectories,t0)
LineTrajectories<-rbind(LineTrajectories,data.frame(stage='Larva',line='L486',pheno='HA',time=30,gene=t30v0$genes,logFC=t30v0$logFC  ) )
LineTrajectories<-rbind(LineTrajectories,data.frame( stage='Larva',line='L486',pheno='HA',time=60,gene=t60v0$genes,logFC=t60v0$logFC ) )
LineTrajectories<-rbind(LineTrajectories,data.frame( stage='Larva',line='L486',pheno='HA',time=90,gene=t90v0$genes,logFC=t90v0$logFC ) )


#HL ######################################

#A441t30 - A441t0
con   <-c(0,                             1,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             1,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             1,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          1,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0 )
lrt <- glmLRT(fitAll, contrast=con)
t30v0<-as.data.frame(topTags(lrt,n=nrow(tableL.dge),sort.by=NULL))

#A441t60 - A441t0
con   <-c(0,                             0,                            0,                             1,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             1,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             1,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            1,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0 )
lrt <- glmLRT(fitAll, contrast=con)
t60v0<-as.data.frame(topTags(lrt,n=nrow(tableL.dge),sort.by=NULL))

#A441t30a - A441t0
con   <-c(0,                             0,                            1,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            1,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            1,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             1,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0 )
lrt <- glmLRT(fitAll, contrast=con)
t90v0<-as.data.frame(topTags(lrt,n=nrow(tableL.dge),sort.by=NULL))

t0$line='L441'
t0$pheno='HL'
LineTrajectories<-rbind(LineTrajectories,t0)
LineTrajectories<-rbind(LineTrajectories,data.frame(stage='Larva',line='L441',pheno='HL',time=30,gene=t30v0$genes,logFC=t30v0$logFC  ) )
LineTrajectories<-rbind(LineTrajectories,data.frame( stage='Larva',line='L441',pheno='HL',time=60,gene=t60v0$genes,logFC=t60v0$logFC ) )
LineTrajectories<-rbind(LineTrajectories,data.frame( stage='Larva',line='L441',pheno='HL',time=90,gene=t90v0$genes,logFC=t90v0$logFC ) )


#A832t30 - A832t0
con   <-c(0,                             1,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             1,                            0,                             0,                               0,                        0,
          0,                             1,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          1,                             0,                            0,                             0,                               0,                        0 )
lrt <- glmLRT(fitAll, contrast=con)
t30v0<-as.data.frame(topTags(lrt,n=nrow(tableL.dge),sort.by=NULL))

#A832t60 - A832t0
con   <-c(0,                             0,                            0,                             1,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             1,                               0,                        0,
          0,                             0,                            0,                             1,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            1,                             0,                               0,                        0 )
lrt <- glmLRT(fitAll, contrast=con)
t60v0<-as.data.frame(topTags(lrt,n=nrow(tableL.dge),sort.by=NULL))

#A832t30a - A832t0
con   <-c(0,                             0,                            1,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            1,                             0,                               0,                        0,
          0,                             0,                            1,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             1,                            0,                             0,                               0,                        0 )
lrt <- glmLRT(fitAll, contrast=con)
t90v0<-as.data.frame(topTags(lrt,n=nrow(tableL.dge),sort.by=NULL))

t0$line='L832'
LineTrajectories<-rbind(LineTrajectories,t0)
LineTrajectories<-rbind(LineTrajectories,data.frame(stage='Larva',line='L832',pheno='HL',time=30,gene=t30v0$genes,logFC=t30v0$logFC  ) )
LineTrajectories<-rbind(LineTrajectories,data.frame( stage='Larva',line='L832',pheno='HL',time=60,gene=t60v0$genes,logFC=t60v0$logFC ) )
LineTrajectories<-rbind(LineTrajectories,data.frame( stage='Larva',line='L832',pheno='HL',time=90,gene=t90v0$genes,logFC=t90v0$logFC ) )


#A913t30 - A913t0
con   <-c(0,                             1,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               1,                        0,
          0,                             1,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             1,                               0,                        0 )
lrt <- glmLRT(fitAll, contrast=con)
t30v0<-as.data.frame(topTags(lrt,n=nrow(tableL.dge),sort.by=NULL))

#A913t60 - A913t0
con   <-c(0,                             0,                            0,                             1,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          1,                             0,                            0,                             1,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        1 )
lrt <- glmLRT(fitAll, contrast=con)
t60v0<-as.data.frame(topTags(lrt,n=nrow(tableL.dge),sort.by=NULL))

#A913t30a - A913t0
con   <-c(0,                             0,                            1,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        1,
          0,                             0,                            1,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               0,                        0,
          0,                             0,                            0,                             0,                               1,                        0 )
lrt <- glmLRT(fitAll, contrast=con)
t90v0<-as.data.frame(topTags(lrt,n=nrow(tableL.dge),sort.by=NULL))

t0$line='L913'
LineTrajectories<-rbind(LineTrajectories,t0)
LineTrajectories<-rbind(LineTrajectories,data.frame(stage='Larva',line='L913',pheno='HL',time=30,gene=t30v0$genes,logFC=t30v0$logFC  ) )
LineTrajectories<-rbind(LineTrajectories,data.frame( stage='Larva',line='L913',pheno='HL',time=60,gene=t60v0$genes,logFC=t60v0$logFC ) )
LineTrajectories<-rbind(LineTrajectories,data.frame( stage='Larva',line='L913',pheno='HL',time=90,gene=t90v0$genes,logFC=t90v0$logFC ) )
