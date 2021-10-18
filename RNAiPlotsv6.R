setwd("~/University/Post-Doc/RNAi")

my.data <- read.csv("RNAiData_v6.csv", stringsAsFactors = FALSE)
str(my.data)

P2.data <- my.data[my.data$Insertion.site=="attP2",]
P40.data <- my.data[my.data$Insertion.site=="attP40",]

P2.data$Line <- factor(P2.data$Line, levels = c("control","CG10505","klu","NtR", "pigs"))
P40.data$Line <- factor(P40.data$Line, levels = c("control","CG32533","Clk","Ir85a","mthl15", "psh"))

larv.P2.data <- P2.data[which(P2.data$Stage=="Larvae"),]
larv.P40.data <- P40.data[which(P40.data$Stage=="Larvae"),]
ad.P2.f.data <- P2.data[which(P2.data$Stage=="Female"),]
ad.P40.f.data <- P40.data[which(P40.data$Stage=="Female"),]
ad.P2.m.data <- P2.data[which(P2.data$Stage=="Male"),]
ad.P40.m.data <- P40.data[which(P40.data$Stage=="Male"),]

#points and errors for plotting relative PSurvival
P2.means <- tapply(larv.P2.data$PSurvival, larv.P2.data$Line,mean) 
P2.means <- P2.means - P2.means[1]
P2.sd <- tapply(larv.P2.data$PSurvival, larv.P2.data$Line,sd) 
P2.N <- tapply(larv.P2.data$PSurvival, larv.P2.data$Line,length)
P2.se <- P2.sd/sqrt(P2.N)

P40.means <- tapply(larv.P40.data$PSurvival, larv.P40.data$Line,mean) 
P40.means <- P40.means - P40.means[1] 
P40.sd <- tapply(larv.P40.data$PSurvival, larv.P40.data$Line,sd) 
P40.N <- tapply(larv.P40.data$PSurvival, larv.P40.data$Line,length)
P40.se <- P40.sd/sqrt(P40.N)

P2.f.means <- tapply(ad.P2.f.data$PSurvival, ad.P2.f.data$Line,mean) 
P2.f.means <- P2.f.means - P2.f.means[1]
P2.f.sd <- tapply(ad.P2.f.data$PSurvival, ad.P2.f.data$Line,sd) 
P2.f.N <- tapply(ad.P2.f.data$PSurvival, ad.P2.f.data$Line,length)
P2.f.se <- P2.f.sd/sqrt(P2.f.N)

P40.f.means <- tapply(ad.P40.f.data$PSurvival, ad.P40.f.data$Line,mean) 
P40.f.means <- P40.f.means - P40.f.means[1]
P40.f.sd <- tapply(ad.P40.f.data$PSurvival, ad.P40.f.data$Line,sd) 
P40.f.N <- tapply(ad.P40.f.data$PSurvival, ad.P40.f.data$Line,length)
P40.f.se <- P40.f.sd/sqrt(P40.f.N)


P2.m.means <- tapply(ad.P2.m.data$PSurvival, ad.P2.m.data$Line,mean) 
P2.m.means <- P2.m.means - P2.m.means[1]
P2.m.sd <- tapply(ad.P2.m.data$PSurvival, ad.P2.m.data$Line,sd) 
P2.m.N <- tapply(ad.P2.m.data$PSurvival, ad.P2.m.data$Line,length)
P2.m.se <- P2.m.sd/sqrt(P2.m.N)

P40.m.means <- tapply(ad.P40.m.data$PSurvival, ad.P40.m.data$Line,mean) 
P40.m.means <- P40.m.means - P40.m.means[1]
P40.m.sd <- tapply(ad.P40.m.data$PSurvival, ad.P40.m.data$Line,sd) 
P40.m.N <- tapply(ad.P40.m.data$PSurvival, ad.P40.m.data$Line,length)
P40.m.se <- P40.m.sd/sqrt(P40.m.N)

#single plot of relative PSurvival
par(mfrow=c(1,1), mar=c(3,3,0,0), oma=c(2,2,0,0), cex=0.75, xpd=TRUE)

plot(P2.f.means~c(1,2,3,4,5), xlim=c(0.5,9.5), ylim=c(-0.7,0.6),type="n",axes=FALSE,
     xlab="RNAi Line", ylab="Proportion survival relative to control")

axis(side=2, at=c(-0.6,-0.4,-0.2,0,0.2,0.4,0.6), cex.axis=1,las=1,
     labels = c(-0.6,-0.4,-0.2,0,0.2,0.4,0.6), pos=0.5, tck=-0.01)
lines(c(0.5,0.5),c(-0.7,0.6))
axis(side=1, at=c(1,2,3,4,5,6,7,8,9), cex.axis=1,las=1,
     labels = c("CG10505", "Clk", "klu", "Ir85a", "CG32533", "mthl15", "NtR", "pigs", "psh"), 
     pos=-0.7, tck=-0.01, font=3)
lines(c(0.5,9.5),c(-0.7,-0.7))

#zero line
lines(c(0.5,9.5), c(0,0), lty=2)

#y-axis title
mtext("Proportion Survival Relative to Control", side=2, outer=TRUE, cex=0.8, font=2)
#x-axis title
mtext("RNAi Line", side=1, outer=TRUE, cex=0.8, font=2)

#legend
legend(3,0.55, horiz=TRUE, c("Larvae",  "Adult females", "Adult males"), 
       pch=c(15,16,17), col=c("#999999", "#CC79A7", "#0072B2"))

Y.larv <- c(P2.means[2], P40.means[3], P2.means[3], P40.means[c(4,2,5)], P2.means[c(4,5)], P40.means[6])
Y.f <- c(P2.f.means[2], P40.f.means[3], P2.f.means[3], P40.f.means[c(4,2,5)], P2.f.means[c(4,5)], P40.f.means[6])
Y.m <- c(P2.m.means[2], P40.m.means[3], P2.m.means[3], P40.m.means[c(4,2,5)], P2.m.means[c(4,5)], P40.m.means[6])

se.larv <- c(P2.se[2], P40.se[3], P2.se[3], P40.se[c(4,2,5)], P2.se[c(4,5)], P40.se[6])
se.f <- c(P2.f.se[2], P40.f.se[3], P2.f.se[3], P40.f.se[c(4,2,5)], P2.f.se[c(4,5)], P40.f.se[6])
se.m <- c(P2.m.se[2], P40.m.se[3], P2.m.se[3], P40.m.se[c(4,2,5)], P2.m.se[c(4,5)], P40.m.se[6])

X <- c(1,2,3,4,5,6,7,8,9)
X.larv <- X - 0.15
X.f <- X
X.m <- X + 0.15

#error bars (standard error)
arrows(X.larv, Y.larv-se.larv, X.larv, Y.larv+se.larv,
       length=0.02, angle=90, code=3) #larvae
arrows(X.f, Y.f-se.f, X.f, Y.f+se.f,
       length=0.02, angle=90, code=3) #females
arrows(X.m, Y.m-se.m, X.m, Y.m+se.m,
       length=0.02, angle=90, code=3) #males
#means
points(X.larv, Y.larv, col="#999999",cex=1, pch=15)
points(X.f, Y.f, col="#CC79A7",cex=1, pch=16)
points(X.m, Y.m, col="#0072B2",cex=1, pch=17)

#stats
text(X.larv[c(4,5,6)], Y.larv[c(4,5,6)]+se.larv[c(4,5,6)]+0.05, "*")
text(X.f[c(1,2,3,5)], Y.f[c(1,2,3,5)]+se.f[c(1,2,3,5)]+0.05, "*")
text(X.m[c(1,2,6)], Y.m[c(1,2,6)]+se.m[c(1,2,6)]+0.05, "*")

#export as "RelativePSurvivalbyLinev6.pdf" 7 x 3.5"