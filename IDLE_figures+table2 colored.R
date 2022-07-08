# read data ===============
library(survivalROC)
library(survival)
library(pROC)
library(Hmisc)
library(timeROC)
library(e1071)
library(Matrix)
library(gplots)
rm(list=ls())

texture.roi<- read.csv("IDLE.csv")
max(texture.roi$days.to.surgery-
    texture.roi$days.to.dx)/365
# 1.558904 
# ==> all surgeries were within 1.558904 years
#     of cancer initial diagnosis
texture.roi$grade1<- 1*(texture.roi$grade==1)
texture.roi$days.pfs<- texture.roi$days.surgery.pfs+texture.roi$days.to.surgery

PPV.NPV<- function(y, score.test)
{
	# y=object from survivalROC
	# score.test = the prediction risk score used by survivalROC in y
	# sensitivity=y$TP
	# specificity=y$FP
	# S(t)=y$Survival=time dependent survival used by y
	# PPV = sensitivity*(1-S(t))/P(Test+)
	#     = y$TP*(1-y$Survival)/P(Test+)
 # NPV = specificity*S(t)/P(Test-)
 #     =(1-y$FP)*y$Survivial/P(Test-)

	# calculate P(Test+)
	x<- 1
	for (i in 2:length(y$TP))
	{	x<- c(x, mean(score.test>y$cut.values[i]))	}
	#x=P(Test+), 1-x=P(Test-)
 y$PPV<- y$TP*(1-y$Survival)/x
 y$NPV<- (1-y$FP)*y$Survival/(1-x)
 # restrict to be within 0 and 1
 y$PPV[y$PPV>1]<- 1
 y$PPV[y$PPV<0]<- 0
 y$NPV[y$NPV>1]<- 1
 y$NPV[y$NPV<0]<- 0
 x<- data.frame(cut.values=y$cut.values,
 	 sensitivity=y$TP, specificity=	1-y$FP, PPV=y$PPV, NPV=y$NPV)
 x
}



# AUC and its standard deviation ========
# bootstrap estimates of 5-year and 10-year
# time dependent AUC confidence interval

simu<- 500
yr<- c(4:6,10)


# IDLE2 =IDLE
{a.rmst<- matrix(NA, nrow=simu,ncol=length(yr))
a.rmst<- as.data.frame(a.rmst)
names(a.rmst)<- paste0("year",yr) #for rmst.ctRoi
a.staging<- a.rmst
a.grade<- a.rmst
a.rmst.roi<- a.rmst #for rmst.roi


set.seed(123)
for (i in 1:simu)
{
	ind<- sample(1:nrow(texture.roi), replace=TRUE)
	dataset<- texture.roi[ind,] #bootstrap sample
	for (j in 1:length(yr)) # year 4,5, 6
	{	
  	# rmst.ctRoi
	  y<- survivalROC(Stime=dataset$days.surgery.pfs,
		status=dataset$ind.lung.pfs,
		marker= dataset$IDLE2,
		predict.time=365*yr[j], method="KM")
    a.rmst[i,j]<- y$AUC 
    
	  # cancer staging
    y<- survivalROC(Stime=dataset$days.surgery.pfs,
		status=dataset$ind.lung.pfs,
		marker=dataset$dia.tumor,
		predict.time=365*yr[j], method="KM")
    a.staging[i,j]<- y$AUC 
    
    # tumor grade
    y<- survivalROC(Stime=dataset$days.surgery.pfs,
		status=dataset$ind.lung.pfs,
		marker=dataset$grade,
		predict.time=365*yr[j], method="KM")
    a.grade[i,j]<- y$AUC
	}
}

# IDLE2 AUC and its standard deviation
apply(a.rmst, 2, mean)
#     year4     year5     year6    year10 
# 0.7997187 0.8157354 0.8043462 0.7896338
sqrt(apply(a.rmst,2,var))
#      year4      year5      year6     year10 
# 0.03957134 0.03782348 0.03861153 0.03938712  

# confidence intervals with Bonferroni adjustment
quantile(a.rmst$year5,probs=c(0.05/2/3, 1-0.05/2/3))
# 0.8333333%  99.16667% 
#  0.7232297  0.8893349 
quantile(a.rmst$year10, probs=c(0.05/2/3, 1-0.05/2/3))
# 0.8333333%  99.16667% 
#  0.6876098  0.8695676 

# - - - - - - - - - - - - - - - -
# staging AUC and its standard deviation
sqrt(apply(a.staging,2,var))
#      year4      year5      year6     year10 
# 0.04268129 0.04245501 0.04241638 0.04143740

# confidence intervals with Bonferroni adjustment
quantile(a.staging$year5,
  probs=c(0.05/2/3, 1-0.05/2/3))
# 0.8333333%  99.16667% 
#  0.4538680  0.6503922

# - - - - - - - - - - - - - - - -
# tumor grade AUC and its standard deviation
sqrt(apply(a.grade,2,var))
#      year4      year5      year6     year10 
# 0.04385640 0.04403028 0.04353295 0.04496273 

# confidence intervals with Bonferroni adjustment
quantile(a.grade$year5,
  probs=c(0.05/2/3, 1-0.05/2/3))
# 0.8333333%  99.16667% 
#  0.4788746  0.6877301


}



# Figure 3  =======
# 5-year and 10-year AUC plot and
# Kaplan-Meier curves using texture.roi

par(mfrow=c(3,3))
{
  cut.off<- 5
# A. plot 5-year AUC in texture.roi
{dataset<- texture.roi
set.seed(123)
	y.ctRoi<- survivalROC(Stime=dataset$days.surgery.pfs,
		status=dataset$ind.lung.pfs, marker=dataset$IDLE2,
		predict.time=cut.off*365, method="KM")
set.seed(123)
	y.staging<- survivalROC(Stime=dataset$days.surgery.pfs,
		status=dataset$ind.lung.pfs, marker=dataset$dia.tumor,
		predict.time=cut.off*365, method="KM")
set.seed(123)
	y.grade<- survivalROC(Stime=dataset$days.surgery.pfs,
		status=dataset$ind.lung.pfs, marker=dataset$grade,
		predict.time=cut.off*365, method="KM")

plot(y.ctRoi$FP, y.ctRoi$TP, type="l", xlim=c(0,1), ylim=c(0,1),
	xlab="1-specificity", ylab="sensitivity", 
  col="black",lwd=2,
	cex.main=1, cex.lab=1,cex.axis=1,
  main=paste0("A. ", cut.off, 
   "-year time-dependent ROC"))
lines(y.staging$FP, y.staging$TP,lwd=2, col="cyan3")
lines(y.grade$FP, y.grade$TP,lwd=2, col="red") 
lines(c(0,1),c(0, 1), lty=2, col="grey")
#legend(0.2,0.5,
#paste("AUC=", round(y$AUC,3), sep=""))
text(0.7,0.27,paste0("IDLE: AUC=", round(y.ctRoi$AUC,3)),
  col="black",cex=0.9)
text(0.7,0.19,paste0("staging: AUC=", round(y.staging$AUC,3)),
  col="cyan3",cex=0.9)
text(0.7,0.11,paste0("grade: AUC=", round(y.grade$AUC,3)),
  col="red",cex=0.9)

print(paste0("rmst.ctRoi AUC=", round(y.ctRoi$AUC,4)))
print(paste0("staging AUC=", round(y.staging$AUC,4)))
print(paste0("grade AUC=", round(y.grade$AUC,4)))
# [1] "rmst.ctRoi AUC=0.8166"
# [1] "staging AUC=0.5614"
# [1] "grade AUC=0.5788"
}

# B-C. plot 5-year time dependent 
# PPV, NPV conditioning on sensitivity
# using texture.roi
{
dataset<- texture.roi
set.seed(123)
y.ctRoi<- survivalROC(Stime=dataset$days.surgery.pfs/365,
		status=dataset$ind.lung.pfs, marker= dataset$IDLE2,
		predict.time=cut.off, method="NNE", span = 0.07)
set.seed(123)
	y.staging<- survivalROC(Stime=dataset$days.surgery.pfs/365,
		status=dataset$ind.lung.pfs, marker=dataset$dia.tumor,
		predict.time=cut.off, method="KM")
set.seed(123)
	y.grade<- survivalROC(Stime=dataset$days.surgery.pfs/365,
		status=dataset$ind.lung.pfs, marker=dataset$grade,
	  predict.time=cut.off, method="NNE", span = 0.07)
# rmst.ctRoi
results<- PPV.NPV(y.ctRoi, dataset$IDLE2)
results<- results[-1,]
results<- results[with(results, order(sensitivity)),]
results<- results[results$sensitivity>0.6,]
results<- results[results$sensitivity<0.95,]
results.ctRoi<- results[with(results, order(sensitivity)),]

# staging
results.staging<- PPV.NPV(y.staging, dataset$dia.tumor)
round(results.staging,4)
#   cut.values sensitivity specificity    PPV    NPV
# 1       -Inf      1.0000      0.0000 0.2558    NaN
# 2          1      0.7403      0.4202 0.3050 0.8248
# 3          2      0.0859      0.8819 0.2000 0.7373
# 4          3      0.0000      1.0000    NaN 0.7442

results.grade<- PPV.NPV(y.grade, dataset$grade)
round(results.grade,4)
#   cut.values sensitivity specificity    PPV    NPV
# 1       -Inf      1.0000      0.0000 0.2558    NaN
# 2          1      0.9097      0.2717 0.3004 0.8977
# 3          2      0.3368      0.6876 0.2704 0.7510
# 4          3      0.0930      0.9655 0.4811 0.7559
# 5          4      0.0000      1.0000    NaN 0.7442

# PPV
plot(results.ctRoi$sensitivity, results.ctRoi$PPV,
    main="B. 5-year time-dependent
  positive predictive value", ylim=c(0.25,0.62), type="l", lwd=2,
  xlab="Sensitivity", ylab="", col="black")
lines(results.staging$sensitivity, results.staging$PPV,
  col="cyan3",lwd=2)
lines(results.grade$sensitivity, results.grade$PPV,
  col="red",lwd=2)
legend("topright", c("IDLE","staging","grade"),
   pch=rep(15,3), col=c("black", "cyan3","red"),
  horiz = TRUE, box.lty = 0)

# NPV
plot(results.ctRoi$sensitivity, results.ctRoi$NPV,
    main="C. 5-year time-dependent 
  negative predictive value", ylim=c(0.78,1), type="l", lwd=2,
  xlab="Sensitivity", ylab="", col="black")
lines(results.staging$sensitivity, results.staging$NPV,
  col="cyan3",lwd=2)
lines(results.grade$sensitivity, results.grade$NPV,
  col="red",lwd=2)
legend("topright", c("IDLE","staging","grade"),
   pch=rep(15,3), col=c("black", "cyan3","red"),
  horiz = TRUE, box.lty = 0)
}

cut.off<- 10
# 10-year prediction
# D. plot 10-year AUC in texture.roi
{dataset<- texture.roi
set.seed(123)
	y.ctRoi<- survivalROC(Stime=dataset$days.surgery.pfs,
		status=dataset$ind.lung.pfs, marker=dataset$IDLE2,
		predict.time=cut.off*365, method="KM")
set.seed(123)
	y.staging<- survivalROC(Stime=dataset$days.surgery.pfs,
		status=dataset$ind.lung.pfs, marker=dataset$dia.tumor,
		predict.time=cut.off*365, method="KM")
set.seed(123)
	y.grade<- survivalROC(Stime=dataset$days.surgery.pfs,
		status=dataset$ind.lung.pfs, marker=dataset$grade,
		predict.time=cut.off*365, method="KM")

plot(y.ctRoi$FP, y.ctRoi$TP, type="l", xlim=c(0,1), ylim=c(0,1),
	xlab="1-specificity", ylab="sensitivity",lwd=2,
	cex.main=1, cex.lab=1,cex.axis=1,
  main=paste0("D. ", cut.off, 
   "-year time-dependent ROC"))
lines(y.staging$FP, y.staging$TP,lwd=2, col="cyan3")
lines(y.grade$FP, y.grade$TP,lwd=2, col="red") 
lines(c(0,1),c(0, 1), lty=2, col="grey")
#legend(0.2,0.5,
#paste("AUC=", round(y$AUC,3), sep=""))
text(0.7,0.27,paste0("IDLE: AUC=", round(y.ctRoi$AUC,3)),
  col="black",cex=0.9)
text(0.7,0.19,paste0("staging: AUC=", round(y.staging$AUC,3)),
  col="cyan3",cex=0.9)
text(0.7,0.11,paste0("grade: AUC=", round(y.grade$AUC,3)),
  col="red",cex=0.9)
}

# E-F. plot 10-year time-dependent 
# PPV, NPV conditioning on sensitivity
# using texture.roi
{
dataset<- texture.roi
set.seed(123)
y.ctRoi<- survivalROC(Stime=dataset$days.surgery.pfs/365,
		status=dataset$ind.lung.pfs, marker= dataset$IDLE2,
		predict.time=cut.off, method="NNE", span = 0.07)
set.seed(123)
	y.staging<- survivalROC(Stime=dataset$days.surgery.pfs/365,
		status=dataset$ind.lung.pfs, marker=dataset$dia.tumor,
		predict.time=cut.off, method="KM")
set.seed(123)
	y.grade<- survivalROC(Stime=dataset$days.surgery.pfs/365,
		status=dataset$ind.lung.pfs, marker=dataset$grade,
	  predict.time=cut.off, method="NNE", span = 0.05)

# rmst.ctRoi
results<- PPV.NPV(y.ctRoi, dataset$IDLE2)
results<- results[-1,]
results<- results[with(results, order(sensitivity)),]
results<- results[results$sensitivity>0.6,]
results<- results[results$sensitivity<0.95,]
results.ctRoi<- results[with(results, order(sensitivity)),]

# staging
results.staging<- PPV.NPV(y.staging, dataset$dia.tumor)
round(results.staging,4)
#   cut.values sensitivity specificity    PPV    NPV
# 1       -Inf      1.0000      0.0000 0.3150    NaN
# 2          1      0.6450      0.3902 0.3273 0.7051
# 3          2      0.0897      0.8808 0.2571 0.6778
# 4          3      0.0000      1.0000    NaN 0.6850

results.grade<- PPV.NPV(y.grade, dataset$grade)
round(results.grade,4)
#   cut.values sensitivity specificity    PPV    NPV
# 1       -Inf      1.0000      0.0000 0.3150    NaN
# 2          1      0.8859      0.2764 0.3602 0.8405
# 3          2      0.3415      0.6918 0.3375 0.6955
# 4          3      0.0756      0.9626 0.4815 0.6936
# 5          4      0.0000      1.0000    NaN 0.6850

# PPV
plot(results.ctRoi$sensitivity, results.ctRoi$PPV,
    main="E. 10-year time-dependent positive
  predictive value (sample 2)", ylim=c(0.30,0.7), type="l", lwd=2,
  xlab="Sensitivity", ylab="", col="black")
lines(results.staging$sensitivity, results.staging$PPV,
  col="cyan3",lwd=2)
lines(results.grade$sensitivity, results.grade$PPV,
  col="red",lwd=2)
legend("topright", c("IDLE","staging","grade"),
   pch=rep(15,3), col=c("black", "cyan3","red"),
  horiz = TRUE, box.lty = 0)

# NPV
plot(results.ctRoi$sensitivity, results.ctRoi$NPV,
    main="F. 10-year time-dependent 
  negative predictive value", ylim=c(0.65,1), type="l", lwd=2,
  xlab="Sensitivity", ylab="", col="black")
lines(results.staging$sensitivity, results.staging$NPV,
  col="cyan3",lwd=2)
lines(results.grade$sensitivity, results.grade$NPV,
  col="red",lwd=2)
legend("topright", c("IDLE","staging","grade"),
   pch=rep(15,3), col=c("black", "cyan3","red"),
  horiz = TRUE, box.lty = 0)
}


# G-L. plot Kaplan-Meier curves and 
# hazards ratios 
dataset<- texture.roi
library(survminer) 
# G. rmst.ctRoi
{
x<- coxph(Surv(days.surgery.pfs/365, ind.lung.pfs) ~
		    (IDLE2>=0.271),
	      data = texture.roi, ties = "breslow")
summary(x)
# n= 182, number of events= 54 
#                          coef exp(coef) se(coef)     z Pr(>|z|)    
# IDLE2 >= 0.2714301TRUE 1.7305    5.6435   0.2779 6.228 4.73e-10 ***
#                        exp(coef) exp(-coef) lower .95 upper .95
# IDLE2 >= 0.2714301TRUE     5.643     0.1772     3.274     9.729

x<- survfit(Surv(days.surgery.pfs/365, ind.lung.pfs) ~ 
    (IDLE2>=0.271),
    data = dataset, na.action=na.exclude,
    conf.int = 0.95, conf.type="log-log")
x
#                        n events median 0.95LCL 0.95UCL
# IDLE2 >= 0.271=FALSE 134     24     NA      NA      NA
# IDLE2 >= 0.271=TRUE   48     30   3.25    2.24    6.32
n1<- as.character(c(134,123,108,102,92,78,30))
n2<- as.character(c(48,32,19,15,8,6,3))
plot(x,conf.int=F, lty=rep(1,2),
 log=FALSE, col=c("black", "red"),
	main="G. IDLE", mark.time=TRUE,
	xlab= "Years since initial surgery",
	ylab="Progression free survival",
	cex.lab=1,cex.main=1,cex.axis=1, lwd=3)
legend(6,0.78,c("IDLE<0.271",  "IDLE>0.271"),
	lty=rep(1,2),col=c("black", "red"),#lwd=3,
   pch=rep(15,4),	cex=0.8,box.lty=0)
text(9.5,0.5, "HR=5.643 \n(p<0.0001)")
text(1.5, 0.2, "Number at risk:", cex=0.8)
for (i in 1:length(n1))
{
  text(2*(i-1), 0.10, n1[i], col="black", cex=0.8)
  text(2*(i-1), 0.03, n2[i], col="red", cex=0.8)
}
}

# H. TNM staging
{
#cancer staging KM curve 
dataset<- texture.roi
x<- coxph(Surv(days.surgery.pfs/365, ind.lung.pfs) ~
		    staging,
	      data = dataset, ties = "breslow")
summary(x)
# n= 182, number of events= 54 
#                coef exp(coef) se(coef)      z Pr(>|z|)
# stagingT1b  0.27733   1.31960  0.29389  0.944    0.3454
# stagingT1c -0.08912   0.91473  0.50265 -0.177    0.8593
# 
#            exp(coef) exp(-coef) lower .95 upper .95
# stagingT1b    1.3196     0.7578    0.7418     2.347
# stagingT1c    0.9147     1.0932    0.3415     2.450

x<- survfit(Surv(days.surgery.pfs/365, ind.lung.pfs) ~ 
    staging,
    data = dataset, na.action=na.exclude,
    conf.int = 0.95, conf.type="log-log")
x
#              n events median 0.95LCL 0.95UCL
# staging=T1a 69     19     NA      NA      NA
# staging=T1b 93     30     NA      NA      NA
# staging=T1c 20      5     NA    6.98      NA

n1<- c(69,59,56,52,46,38,11)
n2<- c(93,79,55,50,42,35,16)
n3<- c(20,17,16,15,12,11,6)
plot(x,conf.int=F, lty=rep(1,2),
 log=FALSE, col=c("green", "black", "red"),
	main="H. TNM staging", mark.time=TRUE,
	xlab= "Years since initial surgery",
	ylab="Progression free survival",
	cex.lab=1,cex.main=1,cex.axis=1, lwd=3)
legend("topright",c("T1a", "T1b", "T1c"),
  horiz=T, lty=rep(1,3),
  col=c("green","black", "red"),#lwd=3,
   pch=rep(15,4),	cex=1,box.lty=0)
text(6,0.5,"T1b vs T1a: HR=1.319 (p=0.3454)")
text(6,0.4,"T1c vs T1a: HR=0.914 (p=0.8593)")
text(1.5, 0.27, "Number at risk:", cex=0.8)
for (i in 1:length(n1))
{
  text(2*(i-1), 0.17, n1[i], col="green", cex=0.8)
  text(2*(i-1), 0.10, n2[i], col="black", cex=0.8)
  text(2*(i-1), 0.03, n3[i], col="red", cex=0.8)
}


}

# I. Tumor grade
{
texture.roi$grade34<- 1*(texture.roi$grade>=3)
dataset<- texture.roi
x<- coxph(Surv(days.surgery.pfs/365, ind.lung.pfs) ~
		    grade34, data = dataset, subset = grade>0,
        ties = "breslow")
x<- summary(x)
x$coefficients
# n= 179, number of events= 52 
#           coef exp(coef) se(coef)     z Pr(>|z|)
# grade34 0.1823    1.2000   0.2916 0.625    0.531899
#         exp(coef) exp(-coef) lower .95 upper .95
# grade34       1.2     0.8334    0.6775     2.125

x<- survfit(Surv(days.surgery.pfs/365, ind.lung.pfs) ~ 
    grade34, data = dataset, subset = grade>0,
  na.action=na.exclude,
    conf.int = 0.95, conf.type="log-log")
x
#             n events median 0.95LCL 0.95UCL
# grade34=0 121     34     NA      NA      NA
# grade34=1  58     18     NA      NA      NA

n1<- c(121,105,88,82,72,60,18)# with grade=0 excluded
n2<- c(58,47,37,34,27,23,14)

plot(x,conf.int=F, lty=rep(1,2),
 log=FALSE, col=c( "black", "red"),
	main="I. Tumor grade", mark.time=TRUE,
	xlab= "Years since initial surgery",
	ylab="Progression free survival",
	cex.lab=1,cex.main=1,cex.axis=1, lwd=3)
legend("topright",c("grade=1,2", 
  "grade=3,4"), horiz = T,
	lty=rep(1,2),col=c("black", "red"),#lwd=3,
   pch=rep(15,4),	cex=1,box.lty=0)
text(6, 0.5, "HR=1.200 (p=0.5319)")
text(1.5, 0.20, "Number at risk:", cex=0.8)
for (i in 1:length(n1))
{
  text(2*(i-1), 0.10, n1[i], col="black", cex=0.8)
  text(2*(i-1), 0.03, n2[i], col="red", cex=0.8)
}

}

par(mfrow=c(1,1))
}



# Figure 4 heatmap ======
data4B<- read.csv("Figure 4B data.csv")

# Figure 4B
data1<- data4B[, c(-1,-2)]
data1<- as.matrix(data1)
data1<- t(data1)
col.pfs<- rep("red", nrow(data4B))
col.pfs[data4B$ind.lung.pfs==0]<- "green"

library(gplots)
png(filename= 'CT+HEpredictors.png',
	    width=1200, height=900)
heatmap.2(data1,
  scale = "row",
  ColSideColors=col.pfs,
	col = bluered,
  trace = "none",
	density.info = "none",
	cexRow=1,cexCol=1.5,
	margins=c(5,10),
	keysize=0.5  )
 legend("topright",  legend=c("no progress","progressed"), 
         fill=c("green","red"), cex=1, box.lty=0)
graphics.off()


# Table 2 ========

x<- coxph(Surv(days.surgery.pfs/365, ind.lung.pfs) ~ 
   IDLE2h+staging+grade34+age.dx+chemo+radiotherapy,
	data = texture.roi, na.action = na.exclude, method = "efron", robust = F)
x<- summary(x)
round(x$coefficients,4)
#                 coef exp(coef) se(coef)       z Pr(>|z|)
# IDLE2h        1.7353    5.6708   0.2975  5.8321   0.0000
# stagingT1b   -0.1433    0.8665   0.3157 -0.4539   0.6499
# stagingT1c   -0.2603    0.7708   0.5506 -0.4727   0.6364
# grade34      -0.0184    0.9818   0.3013 -0.0610   0.9513
# age.dx        0.0313    1.0318   0.0295  1.0609   0.2888
# chemo        -0.4005    0.6700   0.4440 -0.9020   0.3671
# radiotherapy  0.3335    1.3959   0.6296  0.5297   0.5963
round(x$conf.int,4)
x1<- cbind(x$conf.int[,-2], x$coefficients[, 5])
round(x1, 4)
#              exp(coef) lower .95 upper .95       
# IDLE2h          5.6708    3.1650   10.1605 0.0000
# stagingT1b      0.8665    0.4667    1.6087 0.6499
# stagingT1c      0.7708    0.2620    2.2680 0.6364
# grade34         0.9818    0.5440    1.7720 0.9513
# age.dx          1.0318    0.9738    1.0934 0.2888
# chemo           0.6700    0.2806    1.5996 0.3671
# radiotherapy    1.3959    0.4064    4.7945 0.5963


