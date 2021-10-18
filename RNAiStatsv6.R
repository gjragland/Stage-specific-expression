setwd("~/University/Post-Doc/RNAi")

library(lme4) #package for running models
#runs a generalized linear mixed model with binomial distribution 
#(because data are binary) and with cross.date as a random effect
#(because control survival varied with cross.date)

#read in data
my.data <- read.csv("RNAiData_v6.csv", stringsAsFactors = FALSE)
str(my.data)

#make stage a factor
my.data$Stage <- factor(my.data$Stage, levels=c("Larvae", "Female", "Male"))
levels(my.data$Stage)

##subset data into the two genetic backgrounds (attP2 and attP40)
P2.data <- my.data[my.data$Insertion.site=="attP2",]
P40.data <- my.data[my.data$Insertion.site=="attP40",]

#attP2 control line subset by stage and cross.date
P2.larv.May <- P2.data[P2.data$Line=="control" & P2.data$Stage=="Larvae" & P2.data$Cross.date=="May",]
P2.larv.June <- P2.data[P2.data$Line=="control" & P2.data$Stage=="Larvae" & P2.data$Cross.date=="June",]
P2.larv.July <- P2.data[P2.data$Line=="control" & P2.data$Stage=="Larvae" & P2.data$Cross.date=="July",]
P2.larv.Feb <- P2.data[P2.data$Line=="control" & P2.data$Stage=="Larvae" & P2.data$Cross.date=="Feb",]

P2.f.May <- P2.data[P2.data$Line=="control" & P2.data$Stage=="Female" & P2.data$Cross.date=="May",]
P2.f.June <- P2.data[P2.data$Line=="control" & P2.data$Stage=="Female" & P2.data$Cross.date=="June",]
P2.f.Jan <- P2.data[P2.data$Line=="control" & P2.data$Stage=="Female" & P2.data$Cross.date=="Jan",]
P2.f.Feb <- P2.data[P2.data$Line=="control" & P2.data$Stage=="Female" & P2.data$Cross.date=="Feb",]

P2.m.May <- P2.data[P2.data$Line=="control" & P2.data$Stage=="Male" & P2.data$Cross.date=="May",]
P2.m.June <- P2.data[P2.data$Line=="control" & P2.data$Stage=="Male" & P2.data$Cross.date=="June",]
P2.m.Feb <- P2.data[P2.data$Line=="control" & P2.data$Stage=="Male" & P2.data$Cross.date=="Feb",]

#attP40 control line subset by stage and cross.date
P40.larv.May <- P40.data[P40.data$Line=="control" & P40.data$Stage=="Larvae" & P40.data$Cross.date=="May",]
P40.larv.June <- P40.data[P40.data$Line=="control" & P40.data$Stage=="Larvae" & P40.data$Cross.date=="June",]
P40.larv.July <- P40.data[P40.data$Line=="control" & P40.data$Stage=="Larvae" & P40.data$Cross.date=="July",]
P40.larv.Feb <- P40.data[P40.data$Line=="control" & P40.data$Stage=="Larvae" & P40.data$Cross.date=="Feb",]

P40.f.May <- P40.data[P40.data$Line=="control" & P40.data$Stage=="Female" & P40.data$Cross.date=="May",]
P40.f.June <- P40.data[P40.data$Line=="control" & P40.data$Stage=="Female" & P40.data$Cross.date=="June",]
P40.f.July <- P40.data[P40.data$Line=="control" & P40.data$Stage=="Female" & P40.data$Cross.date=="July",]

P40.m.May <- P40.data[P40.data$Line=="control" & P40.data$Stage=="Male" & P40.data$Cross.date=="May",]
P40.m.June <- P40.data[P40.data$Line=="control" & P40.data$Stage=="Male" & P40.data$Cross.date=="June",]
P40.m.July <- P40.data[P40.data$Line=="control" & P40.data$Stage=="Male" & P40.data$Cross.date=="July",]

#run one model for each RNAi line compared to the appropriate control
#on matched cross dates

##### attP2-CG10505 #####

P2.data[P2.data$Line=="CG10505",] #view cross.dates
#larvae June July Feb, females May June Feb, males May June Feb

#combine RNAi data with control data from appropriate cross.dates
P2.CG10505 <- rbind(P2.data[P2.data$Line=="CG10505",],
                    P2.larv.June, P2.larv.July, P2.larv.Feb,
                    P2.f.May, P2.f.June, P2.f.Feb, P2.m.May, P2.m.June, P2.m.Feb)
#make line a factor, with control as the baseline
P2.CG10505$Line <- factor(P2.CG10505$Line, levels = c("control","CG10505"))

#run model
P2.model <- glmer(cbind(alive, dead)~Line*Stage +(1|Vial), data=P2.CG10505, 
                  family=binomial(link="logit"))
summary(P2.model)

#Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
#Family: binomial  ( logit )
#Formula: cbind(alive, dead) ~ Line * Stage + (1 | Vial)
#Data: P2.CG10505
#
#AIC      BIC   logLik deviance df.resid 
#217.8    230.3   -101.9    203.8       37 
#
#Scaled residuals: 
#  Min       1Q   Median       3Q      Max 
#-2.03561 -0.38943  0.01002  0.44691  1.49870 
#
#Random effects:
#  Groups Name        Variance Std.Dev.
#Vial   (Intercept) 0.4298   0.6556  
#Number of obs: 44, groups:  Vial, 33
#
#Fixed effects:
#                        Estimate Std. Error z value Pr(>|z|)   
#(Intercept)              -0.2255     0.2244  -1.005  0.31487   
#LineCG10505               0.2223     0.3505   0.634  0.52604   
#StageFemale               0.3367     0.4197   0.802  0.42235   
#StageMale                 1.3110     0.4421   2.965  0.00302 **
#LineCG10505:StageFemale  -1.8107     0.7268  -2.491  0.01273 * 
#LineCG10505:StageMale    -2.1896     0.7335  -2.985  0.00283 **
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#
#Correlation of Fixed Effects:
#  (Intr) LnCG10505 StgFml StagMl LCG10505:SF
#LineCG10505 -0.640                                    
#StageFemale -0.535  0.342                             
#StageMale   -0.509  0.326     0.605                   
#LCG10505:SF  0.310 -0.483    -0.578 -0.353            
#LCG10505:SM  0.308 -0.479    -0.365 -0.606  0.561


##### attP2-klu #####

P2.data[P2.data$Line=="klu",] 
#larvae May July Feb, females June Jan Feb, males June Feb
P2.klu <- rbind(P2.data[P2.data$Line=="klu",],
                    P2.larv.May, P2.larv.July, P2.larv.Feb,
                    P2.f.June, P2.f.Jan, P2.f.Feb, P2.m.June, P2.m.Feb)
P2.klu$Line <- factor(P2.klu$Line, levels = c("control","klu"))

P2.model <- glmer(cbind(alive, dead)~Line*Stage +(1|Vial), data=P2.klu, 
                  family=binomial(link="logit"))
summary(P2.model)

#Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
#Family: binomial  ( logit )
#Formula: cbind(alive, dead) ~ Line * Stage + (1 | Vial)
#Data: P2.klu
#
#AIC      BIC   logLik deviance df.resid 
#235.7    248.9   -110.8    221.7       42 
#
#Scaled residuals: 
#  Min       1Q   Median       3Q      Max 
#-1.97586 -0.62812  0.02488  0.58768  1.40405 
#
#Random effects:
#  Groups Name        Variance Std.Dev.
#Vial   (Intercept) 0.3493   0.591   
#Number of obs: 49, groups:  Vial, 36
#
#Fixed effects:
#                    Estimate Std. Error z value Pr(>|z|)   
#(Intercept)         -0.23118    0.21812  -1.060   0.2892   
#Lineklu             -0.08865    0.33386  -0.266   0.7906   
#StageFemale          0.39135    0.39821   0.983   0.3257   
#StageMale            1.26846    0.44534   2.848   0.0044 **
#Lineklu:StageFemale  1.06318    0.57998   1.833   0.0368 * 
#Lineklu:StageMale    0.43839    0.63233   0.693   0.4881   
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#
#Correlation of Fixed Effects:
#   (Intr) Linekl StgFml StagMl Lnk:SF
#Lineklu     -0.653                            
#StageFemale -0.548  0.357                     
#StageMale   -0.490  0.320  0.550              
#Lnkl:StgFml  0.375 -0.576 -0.685 -0.377       
#Linkl:StgMl  0.343 -0.528 -0.386 -0.703  0.563

##### attP2-NtR #####

P2.data[P2.data$Line=="NtR",] 
#larvae July Feb, females June Feb, males June Feb
P2.NtR <- rbind(P2.data[P2.data$Line=="NtR",],
                P2.larv.July, P2.larv.Feb,
                P2.f.June, P2.f.Feb, P2.m.June, P2.m.Feb)
P2.NtR$Line <- factor(P2.NtR$Line, levels = c("control","NtR"))

P2.model <- glmer(cbind(alive, dead)~Line*Stage +(1|Vial), data=P2.NtR, 
                  family=binomial(link="logit"))
summary(P2.model)

#Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
#Family: binomial  ( logit )
#Formula: cbind(alive, dead) ~ Line * Stage + (1 | Vial)
#Data: P2.NtR
#
#AIC      BIC   logLik deviance df.resid 
#180.0    190.9    -83.0    166.0       28 
#
#Scaled residuals: 
#  Min       1Q   Median       3Q      Max 
#-1.86435 -0.37096 -0.01017  0.30159  1.35766 
#
#Random effects:
#  Groups Name        Variance Std.Dev.
#Vial   (Intercept) 0.7059   0.8402  
#Number of obs: 35, groups:  Vial, 25
#
#Fixed effects:
#                     Estimate Std. Error z value Pr(>|z|)  
#(Intercept)         -0.178500   0.323432  -0.552   0.5810  
#LineNtR             -0.001619   0.514025  -0.003   0.9975  
#StageFemale          0.553310   0.545793   1.014   0.3107  
#StageMale            1.346745   0.568095   2.371   0.0178 *
#LineNtR:StageFemale  0.273881   0.877779   0.312   0.7550  
#LineNtR:StageMale   -0.475663   0.892406  -0.533   0.5940  
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#
#Correlation of Fixed Effects:
#  (Intr) LinNtR StgFml StagMl LNR:SF
#LineNtR     -0.629                            
#StageFemale -0.593  0.372                     
#StageMale   -0.570  0.357  0.722              
#LnNtR:StgFm  0.368 -0.586 -0.620 -0.446       
#LnNtR:StgMl  0.362 -0.576 -0.459 -0.634  0.720

###### attP2-pigs #####
P2.data[P2.data$Line=="pigs",] 
#larvae Feb, females June Feb, males June Feb
P2.pigs <- rbind(P2.data[P2.data$Line=="pigs",],
                P2.larv.Feb,
                P2.f.June, P2.f.Feb, P2.m.June, P2.m.Feb)
P2.pigs$Line <- factor(P2.pigs$Line, levels = c("control", "pigs"))

P2.model <- glmer(cbind(alive, dead)~Line*Stage +(1|Vial), data=P2.pigs, 
                  family=binomial(link="logit"))
summary(P2.model) 
#error: boundary (singular) fit: see ?isSingular
#Vial has no effect

#run regular glm (no need for random vial effect)
P2.model <- glm(cbind(alive, dead)~Line*Stage, data=P2.pigs, 
                 family="binomial")
summary(P2.model)

#Call:
#  glm(formula = cbind(alive, dead) ~ Line * Stage, family = "binomial", 
#      data = P2.pigs)
#
#Deviance Residuals: 
#  Min       1Q   Median       3Q      Max  
#-2.2370  -0.5716  -0.2220   0.8495   2.3445  
#
#Coefficients:
#                     Estimate Std. Error z value Pr(>|z|)   
#(Intercept)            0.6190     0.2344   2.641  0.00827 **
#Linepigs              -0.3612     0.3272  -1.104  0.26967   
#StageFemale           -0.2826     0.3514  -0.804  0.42139   
#StageMale              0.4106     0.3813   1.077  0.28162   
#Linepigs:StageFemale  -0.2629     0.4837  -0.544  0.58673   
#Linepigs:StageMale     0.4302     0.5352   0.804  0.42151   
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#
#(Dispersion parameter for binomial family taken to be 1)
#
#Null deviance: 65.948  on 33  degrees of freedom
#Residual deviance: 45.851  on 28  degrees of freedom
#AIC: 142.92
#
#Number of Fisher Scoring iterations: 4

##### attP40-CG32533 #####

P40.data[P40.data$Line=="CG32533",] #view cross.dates
#larvae May July Feb, females May July, males May July

#combine RNAi data with control data from appropriate cross.dates
P40.CG32533 <- rbind(P40.data[P40.data$Line=="CG32533",],
                    P40.larv.May, P40.larv.July, P40.larv.Feb,
                    P40.f.May, P40.f.July, P40.m.May, P40.m.July)
#make line a factor, with control as the baseline
P40.CG32533$Line <- factor(P40.CG32533$Line, levels = c("control","CG32533"))

#run model
P40.model <- glmer(cbind(alive, dead)~Line*Stage +(1|Vial), data=P40.CG32533, 
                  family=binomial(link="logit"))
summary(P40.model)

#Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
#Family: binomial  ( logit )
#Formula: cbind(alive, dead) ~ Line * Stage + (1 | Vial)
#Data: P40.CG32533
#
#AIC      BIC   logLik deviance df.resid 
#231.9    245.0   -108.9    217.9       41 
#
#Scaled residuals: 
#  Min       1Q   Median       3Q      Max 
#-1.44306 -0.53304 -0.05996  0.41841  1.57510 
#
#Random effects:
#  Groups Name        Variance Std.Dev.
#Vial   (Intercept) 0.6065   0.7788  
#Number of obs: 48, groups:  Vial, 33
#
#Fixed effects:
#                        Estimate Std. Error z value Pr(>|z|)   
#(Intercept)              0.40905    0.27898   1.466  0.14259   
#LineCG32533             -1.13383    0.44892  -2.526  0.01155 * 
#StageFemale             -1.22740    0.46982  -2.613  0.00899 **
#StageMale                0.07103    0.46223   0.154  0.87787   
#LineCG32533:StageFemale  1.73041    0.71020   2.436  0.01483 * 
#LineCG32533:StageMale    1.14353    0.70920   1.612  0.10687   
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#
#Correlation of Fixed Effects:
#  (Intr) LnCG32533 StgFml StagMl LCG32533:SF
#LineCG32533 -0.623                                    
#StageFemale -0.595  0.372                             
#StageMale   -0.603  0.375     0.705                   
#LCG32533:SF  0.394 -0.633    -0.662 -0.467            
#LCG32533:SM  0.395 -0.633    -0.462 -0.652  0.724    

##### attP40-Clk #####

P40.data[P40.data$Line=="Clk",]
#larvae June July Feb, females May July, males May July
P40.Clk <- rbind(P40.data[P40.data$Line=="Clk",],
                     P40.larv.June, P40.larv.July, P40.larv.Feb,
                     P40.f.May, P40.f.July, P40.m.May, P40.m.July)
P40.Clk$Line <- factor(P40.Clk$Line, levels = c("control","Clk"))

P40.model <- glmer(cbind(alive, dead)~Line*Stage +(1|Vial), data=P40.Clk, 
                   family=binomial(link="logit"))
summary(P40.model)

#Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
#Family: binomial  ( logit )
#Formula: cbind(alive, dead) ~ Line * Stage + (1 | Vial)
#Data: P40.Clk
#
#AIC      BIC   logLik deviance df.resid 
#248.5    262.5   -117.2    234.5       48 
#
#Scaled residuals: 
#  Min      1Q  Median      3Q     Max 
#-1.0945 -0.5278 -0.1425  0.3880  1.4364 
#
#Random effects:
#  Groups Name        Variance Std.Dev.
#Vial   (Intercept) 1.018    1.009   
#Number of obs: 55, groups:  Vial, 39
#
#Fixed effects:
#                    Estimate Std. Error z value Pr(>|z|)    
#(Intercept)           0.5864     0.3144   1.865 0.062184 .  
#LineClk              -0.2197     0.4808  -0.457 0.647790    
#StageFemale          -1.4195     0.5429  -2.615 0.008928 ** 
#StageMale            -0.0966     0.5362  -0.180 0.857027    
#LineClk:StageFemale  -1.8460     0.8912  -2.071 0.038326 *  
#LineClk:StageMale    -3.0031     0.8744  -3.434 0.000594 ***
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#
#Correlation of Fixed Effects:
#  (Intr) LinClk StgFml StagMl LnC:SF
#LineClk     -0.650                            
#StageFemale -0.580  0.376                     
#StageMale   -0.586  0.382  0.774              
#LnClk:StgFm  0.345 -0.545 -0.601 -0.469       
#LnClk:StgMl  0.351 -0.556 -0.467 -0.610  0.680

##### attP40-Ir85a #####

P40.data[P40.data$Line=="Ir85a",]
#larvae May July, females MayJune  July, males May June July
P40.Ir85a <- rbind(P40.data[P40.data$Line=="Ir85a",],
                 P40.larv.May, P40.larv.July,
                 P40.f.May, P40.f.June, P40.f.July, P40.m.May, P40.m.June, P40.m.July)
P40.Ir85a$Line <- factor(P40.Ir85a$Line, levels = c("control","Ir85a"))

P40.model <- glmer(cbind(alive, dead)~Line*Stage +(1|Vial), data=P40.Ir85a, 
                   family=binomial(link="logit"))
summary(P40.model) 

#Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
#Family: binomial  ( logit )
#Formula: cbind(alive, dead) ~ Line * Stage + (1 | Vial)
#Data: P40.Ir85a
#
#AIC      BIC   logLik deviance df.resid 
#251.5    265.0   -118.8    237.5       44 
#
#Scaled residuals: 
#  Min      1Q  Median      3Q     Max 
#-1.7370 -0.6859 -0.2434  0.6465  2.4583 
#
#Random effects:
#  Groups Name        Variance Std.Dev.
#Vial   (Intercept) 0.4022   0.6342  
#Number of obs: 51, groups:  Vial, 33
#
#Fixed effects:
#                      Estimate Std. Error z value Pr(>|z|)  
#(Intercept)           -0.07221    0.29881  -0.242   0.8091  
#LineIr85a             -0.82062    0.41707  -1.968   0.0491 *
#StageFemale           -0.55545    0.43543  -1.276   0.2021  
#StageMale              0.53572    0.43156   1.241   0.2145  
#LineIr85a:StageFemale  0.18542    0.62271   0.298   0.7659  
#LineIr85a:StageMale    0.75462    0.60764   1.242   0.2143  
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#
#Correlation of Fixed Effects:
#  (Intr) LnIr85 StgFml StagMl LI85:SF
#LineIr85a   -0.717                             
#StageFemale -0.686  0.494                      
#StageMale   -0.692  0.495  0.711               
#LnIr85:StgF  0.480 -0.667 -0.699 -0.498        
#LnIr85:StgM  0.492 -0.686 -0.507 -0.709  0.692 

##### attP40-mthl15 #####

P40.data[P40.data$Line=="mthl15",]
#larvae July Feb, females June July, males June July
P40.mthl15 <- rbind(P40.data[P40.data$Line=="mthl15",],
                 P40.larv.July, P40.larv.Feb,
                 P40.f.June, P40.f.July, P40.m.June, P40.m.July)
P40.mthl15$Line <- factor(P40.mthl15$Line, levels = c("control","mthl15"))

P40.model <- glmer(cbind(alive, dead)~Line*Stage +(1|Vial), data=P40.mthl15, 
                   family=binomial(link="logit"))
summary(P40.model)

#Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
#Family: binomial  ( logit )
#Formula: cbind(alive, dead) ~ Line * Stage + (1 | Vial)
#Data: P40.mthl15
#
#AIC      BIC   logLik deviance df.resid 
#185.7    197.3    -85.8    171.7       32 
#
#Scaled residuals: 
#  Min      1Q  Median      3Q     Max 
#-1.7006 -0.6957 -0.1084  0.7321  1.4522 
#
#Random effects:
#  Groups Name        Variance Std.Dev.
#Vial   (Intercept) 0.3862   0.6214  
#Number of obs: 39, groups:  Vial, 26
#
#Fixed effects:
#                       Estimate Std. Error z value Pr(>|z|)   
#(Intercept)             0.75845    0.28450   2.666  0.00768 **
#Linemthl15             -1.45342    0.45813  -3.172  0.00151 **
#StageFemale            -1.33827    0.48106  -2.782  0.00540 **
#StageMale              -0.01226    0.47919  -0.026  0.97960   
#Linemthl15:StageFemale  1.00965    0.70381   1.435  0.15142   
#Linemthl15:StageMale    2.18081    0.71788   3.038  0.00238 **
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#
#Correlation of Fixed Effects:
#  (Intr) Lnmt15 StgFml StagMl L15:SF
#Linemthl15  -0.625                            
#StageFemale -0.594  0.372                     
#StageMale   -0.591  0.368  0.627              
#Lnmthl15:SF  0.407 -0.651 -0.684 -0.429       
#Lnmthl15:SM  0.399 -0.639 -0.423 -0.669  0.650

##### attP40-psh #####

P40.data[P40.data$Line=="psh",]
#larvae May June July Feb, females June July, males June July
P40.psh <- rbind(P40.data[P40.data$Line=="psh",],
                 P40.larv.May, P40.larv.June, P40.larv.July, P40.larv.Feb,
                 P40.f.June, P40.f.July, P40.m.June, P40.m.July)
P40.psh$Line <- factor(P40.psh$Line, levels = c("control","psh"))

P40.model <- glmer(cbind(alive, dead)~Line*Stage +(1|Vial), data=P40.psh, 
                   family=binomial(link="logit"))
summary(P40.model)

#Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
#Family: binomial  ( logit )
#Formula: cbind(alive, dead) ~ Line * Stage + (1 | Vial)
#Data: P40.psh
#
#AIC      BIC   logLik deviance df.resid 
#277.3    291.1   -131.6    263.3       46 
#
#Scaled residuals: 
#  Min      1Q  Median      3Q     Max 
#-1.9876 -0.5707 -0.1991  0.7728  1.6292 
#
#Random effects:
#  Groups Name        Variance Std.Dev.
#Vial   (Intercept) 0.3778   0.6147  
#Number of obs: 53, groups:  Vial, 42
#
#Fixed effects:
#                    Estimate Std. Error z value Pr(>|z|)  
#(Intercept)           0.3351     0.1953   1.716   0.0861 .
#Linepsh              -0.2784     0.2792  -0.997   0.3186  
#StageFemale          -0.9143     0.4319  -2.117   0.0342 *
#StageMale             0.4105     0.4310   0.952   0.3410  
#Linepsh:StageFemale   1.0381     0.6276   1.654   0.0981 .
#Linepsh:StageMale     0.5497     0.6484   0.848   0.3966  
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#
#Correlation of Fixed Effects:
#  (Intr) Linpsh StgFml StagMl Lnp:SF
#Linepsh     -0.699                            
#StageFemale -0.453  0.316                     
#StageMale   -0.452  0.316  0.539              
#Lnpsh:StgFm  0.312 -0.445 -0.688 -0.371       
#Lnpsh:StgMl  0.301 -0.430 -0.359 -0.665  0.531

P40.model <- glm(cbind(alive, dead)~Line*Stage, data=P40.psh, 
                   family="binomial")
summary(P40.model)
