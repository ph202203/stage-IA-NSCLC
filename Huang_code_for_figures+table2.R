
library(survivalROC)
library(survival)
library(pROC)
library(Hmisc)
library(timeROC)
library(e1071)
library(Matrix)
library(gplots)
rm(list=ls())

roi.stage1<- read.csv("IDLE1.csv")
texture.roi<- read.csv("IDLE2.csv")
max(texture.roi$days.to.surgery-
    texture.roi$days.to.dx)/365
# 1.558904 
# ==> all surgeries were within 1.558904 years
#     of cancer initial diagnosis
roi.stage1$grade1<- 1*(roi.stage1$grade==1)
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


# IDLE1
{a.rmst<- matrix(NA, nrow=simu,ncol=length(yr))
a.rmst<- as.data.frame(a.rmst) 
names(a.rmst)<- paste0("year",yr)
a.staging<- a.rmst
a.grade<- a.rmst

set.seed(123)
for (i in 1:simu)
{
	ind<- sample(1:nrow(roi.stage1), replace=TRUE)
	dataset<- roi.stage1[ind,] #bootstrap sample
	for (j in 1:length(yr)) # year 4,5, 6
	{	
  	# rmst.roi
	  y<- survivalROC(Stime=dataset$days.surgery.pfs,
		status=dataset$ind.lung.pfs,
		marker= dataset$IDLE1,
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

# year 4,5,6, 10 rmst.roi AUC standard deviation
sqrt(apply(a.rmst,2,var))
#      year4      year5      year6     year10 
# 0.03725680 0.03567214 0.03464121 0.03588640 
# - - - - - - - - - - - - - - - -
# year 4,5,6,10 staging AUC and its standard deviation
apply(a.staging, 2, mean)
#     year4     year5     year6    year10 
# 0.5566203 0.5452142 0.5423832 0.5156950
sqrt(apply(a.staging,2,var))
#      year4      year5      year6     year10 
# 0.03629808 0.03549619 0.03513032 0.03597050 

# - - - - - - - - - - - - - - - -
# year 4,5,6,10 tumor grade AUC and its standard deviation
apply(a.grade, 2, mean)
#     year4     year5     year6    year10 
# 0.5453314 0.5314251 0.5423595 0.5441878 
sqrt(apply(a.grade,2,var))
#      year4      year5      year6     year10 
# 0.03687105 0.03592194 0.03604158 0.03623576 
# - - - - - - - - - - - - - - - -
}

# IDLE2
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
    
    # rmst.roi
	  y<- survivalROC(Stime=dataset$days.surgery.pfs,
		status=dataset$ind.lung.pfs,
		marker= dataset$IDLE1,
		predict.time=365*yr[j], method="KM")
    a.rmst.roi[i,j]<- y$AUC 
	  
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
quantile(a.rmst$year5,
  probs=c(0.05/2/3, 1-0.05/2/3))
# 0.8333333%  99.16667% 
#  0.7232297  0.8893349 
quantile(a.rmst$year10,
  probs=c(0.05/2/3, 1-0.05/2/3))
# 0.8333333%  99.16667% 
#  0.6876098  0.8695676 

# - - - - - - - - - - - - - - - -
# IDLE1 AUC and its standard deviation
sqrt(apply(a.rmst.roi,2,var))
#      year4      year5      year6     year10 
# 0.04999770 0.04965809 0.04880446 0.04712344 

# confidence intervals with Bonferroni adjustment
quantile(a.rmst.roi$year5,
  probs=c(0.05/2/3, 1-0.05/2/3))
# 0.8333333%  99.16667% 
#  0.5478984  0.7656477 
quantile(a.rmst.roi$year10,
  probs=c(0.05/2/3, 1-0.05/2/3))
# 0.8333333%  99.16667% 
#  0.5241175  0.7475834 

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

par(mfrow=c(4,3))
{
  cut.off<- 5
# A. plot 5-year AUC in texture.roi
{dataset<- texture.roi
set.seed(123)
	y.ctRoi<- survivalROC(Stime=dataset$days.surgery.pfs,
		status=dataset$ind.lung.pfs, marker=dataset$IDLE2,
		predict.time=cut.off*365, method="KM")
set.seed(123)
	y.roi<- survivalROC(Stime=dataset$days.surgery.pfs,
		status=dataset$ind.lung.pfs, marker=dataset$IDLE1,
		predict.time=cut.off*365, method="KM")
set.seed(123)
	y.staging<- survivalROC(Stime=dataset$days.surgery.pfs,
		status=dataset$ind.lung.pfs, marker=dataset$dia.tumor,
		predict.time=cut.off*365, method="KM")
set.seed(123)
	y.grade<- survivalROC(Stime=dataset$days.surgery.pfs,
		status=dataset$ind.lung.pfs, marker=dataset$grade,
		predict.time=cut.off*365, method="KM")

plot(y.roi$FP, y.roi$TP, type="l", xlim=c(0,1), ylim=c(0,1),
	xlab="1-specificity", ylab="sensitivity", 
  col="darkorange4",lwd=2,
	cex.main=1, cex.lab=1,cex.axis=1,
  main=paste0("A. ", cut.off, 
   "-year time-dependent ROC \n(sample 2)"))
lines(y.ctRoi$FP, y.ctRoi$TP,lwd=2, col="black")
lines(y.staging$FP, y.staging$TP,lwd=2, col="cyan3")
lines(y.grade$FP, y.grade$TP,lwd=2, col="red") 
lines(c(0,1),c(0, 1), lty=2, col="grey")
#legend(0.2,0.5,
#paste("AUC=", round(y$AUC,3), sep=""))
text(0.7,0.35,paste0("IDLE1: AUC=", round(y.roi$AUC,3)),
  col="darkorange4",cex=0.9)
text(0.7,0.27,paste0("IDLE2: AUC=", round(y.ctRoi$AUC,3)),
  col="black",cex=0.9)
text(0.7,0.19,paste0("staging: AUC=", round(y.staging$AUC,3)),
  col="cyan3",cex=0.9)
text(0.7,0.11,paste0("grade: AUC=", round(y.grade$AUC,3)),
  col="red",cex=0.9)

}

# B-C. plot 5-year time dependent 
# PPV, NPV conditioning on sensitivity
# using texture.roi
{
dataset<- texture.roi
set.seed(123)
y.roi<- survivalROC(Stime=dataset$days.surgery.pfs/365,
		status=dataset$ind.lung.pfs, marker= dataset$IDLE1,
		predict.time=cut.off, method="NNE", span = 0.09)
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
# IDLE1
results<- PPV.NPV(y.roi, dataset$IDLE1)
results<- results[-1,]
results<- results[with(results, order(sensitivity)),]
results<- results[results$sensitivity>0.6,]
results<- results[results$sensitivity<0.95,]
results.roi<- results[with(results, order(sensitivity)),]

# IDLE2
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
    main="B. 5-year time-dependent positive
  predictive value (sample 2)", ylim=c(0.25,0.62), type="l", lwd=2,
  xlab="Sensitivity", ylab="", col="black")
lines(results.roi$sensitivity, results.roi$PPV,
  col="darkorange4",lwd=2)
lines(results.staging$sensitivity, results.staging$PPV,
  col="cyan3",lwd=2)
lines(results.grade$sensitivity, results.grade$PPV,
  col="red",lwd=2)
legend("topright", c("IDLE1","IDLE2","staging","grade"),
   pch=rep(15,3), col=c("darkorange4","black", "cyan3","red"),
  horiz = TRUE, box.lty = 0)

# NPV
plot(results.ctRoi$sensitivity, results.ctRoi$NPV,
    main="C. 5-year time-dependent negative
  predictive value (sample 2)", ylim=c(0.78,1), type="l", lwd=2,
  xlab="Sensitivity", ylab="", col="black")
lines(results.roi$sensitivity, results.roi$NPV,
  col="darkorange4",lwd=2)
lines(results.staging$sensitivity, results.staging$NPV,
  col="cyan3",lwd=2)
lines(results.grade$sensitivity, results.grade$NPV,
  col="red",lwd=2)
legend("topright", c("IDLE1","IDLE2","staging","grade"),
   pch=rep(15,3), col=c("darkorange4","black", "cyan3","red"),
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
	y.roi<- survivalROC(Stime=dataset$days.surgery.pfs,
		status=dataset$ind.lung.pfs, marker=dataset$IDLE1,
		predict.time=cut.off*365, method="KM")
set.seed(123)
	y.staging<- survivalROC(Stime=dataset$days.surgery.pfs,
		status=dataset$ind.lung.pfs, marker=dataset$dia.tumor,
		predict.time=cut.off*365, method="KM")
set.seed(123)
	y.grade<- survivalROC(Stime=dataset$days.surgery.pfs,
		status=dataset$ind.lung.pfs, marker=dataset$grade,
		predict.time=cut.off*365, method="KM")

plot(y.roi$FP, y.roi$TP, type="l", xlim=c(0,1), ylim=c(0,1),
	xlab="1-specificity", ylab="sensitivity", 
  col="darkorange4",lwd=2,
	cex.main=1, cex.lab=1,cex.axis=1,
  main=paste0("D. ", cut.off, 
   "-year time-dependent ROC \n(sample 2)"))
lines(y.ctRoi$FP, y.ctRoi$TP,lwd=2, col="black")
lines(y.staging$FP, y.staging$TP,lwd=2, col="cyan3")
lines(y.grade$FP, y.grade$TP,lwd=2, col="red") 
lines(c(0,1),c(0, 1), lty=2, col="grey")
text(0.7,0.35,paste0("IDLE1: AUC=", round(y.roi$AUC,3)),
  col="darkorange4",cex=0.9)
text(0.7,0.27,paste0("IDLE2: AUC=", round(y.ctRoi$AUC,3)),
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
y.roi<- survivalROC(Stime=dataset$days.surgery.pfs/365,
		status=dataset$ind.lung.pfs, marker= dataset$IDLE1,
		predict.time=cut.off, method="NNE", span = 0.09)
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
# IDLE1
results<- PPV.NPV(y.roi, dataset$IDLE1)
results<- results[-1,]
results<- results[with(results, order(sensitivity)),]
results<- results[results$sensitivity>0.6,]
results<- results[results$sensitivity<0.95,]
results.roi<- results[with(results, order(sensitivity)),]

# IDLE2
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
lines(results.roi$sensitivity, results.roi$PPV,
  col="darkorange4",lwd=2)
lines(results.staging$sensitivity, results.staging$PPV,
  col="cyan3",lwd=2)
lines(results.grade$sensitivity, results.grade$PPV,
  col="red",lwd=2)
legend("topright", c("IDLE1","IDLE2","staging","grade"),
   pch=rep(15,3), col=c("darkorange4","black", "cyan3","red"),
  horiz = TRUE, box.lty = 0)

# NPV
plot(results.ctRoi$sensitivity, results.ctRoi$NPV,
    main="F. 10-year time-dependent negative
  predictive value (sample 2)", ylim=c(0.65,1), type="l", lwd=2,
  xlab="Sensitivity", ylab="", col="black")
lines(results.roi$sensitivity, results.roi$NPV,
  col="darkorange4",lwd=2)
lines(results.staging$sensitivity, results.staging$NPV,
  col="cyan3",lwd=2)
lines(results.grade$sensitivity, results.grade$NPV,
  col="red",lwd=2)
legend("topright", c("IDLE1","IDLE2","staging","grade"),
   pch=rep(15,3), col=c("darkorange4","black", "cyan3","red"),
  horiz = TRUE, box.lty = 0)

}


# G-L. plot Kaplan-Meier curves and 
# hazards ratios 
dataset<- texture.roi
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
# use ggsurvplot to find number at risk n1 and n2
n1<- as.character(c(134,123,108,102,92,78,30))
n2<- as.character(c(48,32,19,15,8,6,3))
plot(x,conf.int=F, lty=rep(1,2),
 log=FALSE, col=c("black", "red"),
	main="G. IDLE2 (sample 2)", mark.time=TRUE,
	xlab= "Years since initial surgery",
	ylab="Progression free survival",
	cex.lab=1,cex.main=1,cex.axis=1, lwd=3)
legend(6,0.78,c("IDLE2<0.271",  "IDLE2>0.271"),
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
	main="H. TNM staging (sample 2)", mark.time=TRUE,
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

n1<- c(121,105,88,82,72,60,18)
n2<- c(58,47,37,34,27,23,14)

plot(x,conf.int=F, lty=rep(1,2),
 log=FALSE, col=c( "black", "red"),
	main="I. Tumor grade (sample 2)", mark.time=TRUE,
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

dataset<- roi.stage1

# J-L. plot Kaplan-Meier curves and 
# hazards ratios using roi.stage1
# J. rmst.roi 
{x<- coxph(Surv(days.surgery.pfs/365, ind.lung.pfs) ~
		    (IDLE1>=0.328),
	      data = dataset, ties = "breslow")
summary(x)
# n= 288, number of events= 95 
#                      coef exp(coef) se(coef)     z Pr(>|z|)    
# IDLE1 >= 0.328TRUE 1.2039    3.3329   0.2443 4.927 8.36e-07 ***
#                    exp(coef) exp(-coef) lower .95 upper .95
# IDLE1 >= 0.328TRUE     3.333        0.3     2.065      5.38
x<- survfit(Surv(days.surgery.pfs/365, ind.lung.pfs) ~ 
    (IDLE1>=0.328),
    data = dataset, na.action=na.exclude,
    conf.int = 0.95, conf.type="log-log")
x
#                        n events median 0.95LCL 0.95UCL
# IDLE1 >= 0.328=FALSE 254     73     NA      NA      NA
# IDLE1 >= 0.328=TRUE   34     22   3.06    1.92    10.3
# use ggsurvplot to find number at risk n1 and n2
n1<- as.character(c(254,217,188,174,154,112,41))
n2<- as.character(c(34,22,13,11,10,8,2))
plot(x,conf.int=F, lty=rep(1,2),
 log=FALSE, col=c("black", "red"),
	main="J. IDLE1 (sample 1)", mark.time=TRUE,
	xlab= "Years since initial surgery",
	ylab="Progression free survival",
	cex.lab=1,cex.main=1,cex.axis=1, lwd=3)
legend("topright",c("IDLE1<0.328",  "IDLE1>0.328"),
	lty=rep(1,2),col=c("black", "red"),#lwd=3,
  horiz = T, pch=rep(15,4),	cex=0.8,box.lty=0)
text(8,0.6, "HR=3.333 (p<0.0001)")
text(2, 0.2, "Number at risk:", cex=0.8)
for (i in 1:length(n1))
{
  text(2*(i-1), 0.10, n1[i], col="black", cex=0.8)
  text(2*(i-1), 0.03, n2[i], col="red", cex=0.8)
}
}

# K. cancer staging in sample 1 (roi.stage1)
{
#cancer staging KM curve 
dataset<- roi.stage1
x<- coxph(Surv(days.surgery.pfs/365, ind.lung.pfs) ~
		    staging,
	      data = dataset, ties = "breslow")
summary(x)
#  n= 288, number of events= 95 
#              coef exp(coef) se(coef)     z Pr(>|z|)
# stagingT1b 0.0672    1.0695   0.2306 0.291    0.7707
# stagingT1c 0.1720    1.1877   0.3079 0.559    0.5763
# 
#            exp(coef) exp(-coef) lower .95 upper .95
# stagingT1b     1.070      0.935    0.6806     1.681
# stagingT1c     1.188      0.842    0.6496     2.172

x<- survfit(Surv(days.surgery.pfs/365, ind.lung.pfs) ~ 
    staging,
    data = dataset, na.action=na.exclude,
    conf.int = 0.95, conf.type="log-log")
x
#               n events median 0.95LCL 0.95UCL
# staging=T1a  95     31     NA      NA      NA
# staging=T1b 149     48     NA      NA      NA
# staging=T1c  44     16     NA    9.36      NA

n1<- c(95,80,74,69,61,45,12)
n2<- c(149,124,97,87,78,56,23)
n3<- c(44,35,30,29,25,19,8)
plot(x,conf.int=F, lty=rep(1,2),
 log=FALSE, col=c("green", "black", "red"),
	main="K. TNM staging (sample 1)", mark.time=TRUE,
	xlab= "Years since initial surgery",
	ylab="Progression free survival",
	cex.lab=1,cex.main=1,cex.axis=1, lwd=3)
legend("topright",c("T1a", "T1b", "T1c"),
  horiz=T, lty=rep(1,3),
  col=c("green","black", "red"),#lwd=3,
   pch=rep(15,4),	cex=1,box.lty=0)
text(6,0.5,"T1b vs T1a: HR=1.070 (p=0.7710)")
text(6,0.4,"T1c vs T1a: HR=1.188 (p=0.5763)")
text(1.5, 0.27, "Number at risk:", cex=0.8)
for (i in 1:length(n1))
{
  text(2*(i-1), 0.17, n1[i], col="green", cex=0.8)
  text(2*(i-1), 0.10, n2[i], col="black", cex=0.8)
  text(2*(i-1), 0.03, n3[i], col="red", cex=0.8)
}


}

# L. Tumor grade
{
roi.stage1$grade34<- 1*(roi.stage1$grade>=3)
dataset<- roi.stage1
x<- coxph(Surv(days.surgery.pfs/365, ind.lung.pfs) ~
		    grade34, data = dataset,subset = grade>0,
     ties = "breslow")
summary(x)
# n= 281, number of events= 91 
#           coef exp(coef) se(coef)     z Pr(>|z|)
# grade34 0.0988    1.1038   0.2155 0.458    0.6466823
#         exp(coef) exp(-coef) lower .95 upper .95
# grade34     1.104     0.9059    0.7235     1.684

x<- survfit(Surv(days.surgery.pfs/365, ind.lung.pfs) ~ 
    grade34,  data = dataset, subset=grade>0,
    na.action=na.exclude,
    conf.int = 0.95, conf.type="log-log")
x
#             n events median 0.95LCL 0.95UCL
# grade34=0 177     56     NA      NA      NA
# grade34=1 104     35     NA    11.2      NA

n1<- c(177,148,125,117,105,76,24)
n2<- c(104,85,71,64,56,42,18)

plot(x,conf.int=F, lty=rep(1,2),
 log=FALSE, col=c( "black", "red"),
	main="L. Tumor grade (sample 1)", mark.time=TRUE,
	xlab= "Years since initial surgery",
	ylab="Progression free survival",
	cex.lab=1,cex.main=1,cex.axis=1, lwd=3)
legend("topright",c("grade=1,2", 
  "grade=3,4"), horiz = T,
	lty=rep(1,2),col=c("black", "red"),#lwd=3,
   pch=rep(15,4),	cex=0.8,box.lty=0)
text(6, 0.5, "HR=1.104 (p=0.6467)")
text(3, 0.20, "Number at risk:", cex=0.8)
for (i in 1:length(n1))
{
  text(2*(i-1), 0.10, n1[i], col="black", cex=0.8)
  text(2*(i-1), 0.03, n2[i], col="red", cex=0.8)
}

}

par(mfrow=c(1,1))
}



# Figure 4 heatmap ======
data4A<- read.csv("Figure 4A data.csv")
data4B<- read.csv("Figure 4B data.csv")

# Figure 4A
data1<- data4A[, c(-1,-2)]
data1<- as.matrix(data1)
data1<- t(data1)
col.pfs<- rep("red", nrow(data4A))
col.pfs[data4A$ind.lung.pfs==0]<- "green"

png(filename='HEpredictors.png',
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
         fill=c("green","red"), cex=0.8, box.lty=0)
graphics.off()




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


# Figure S1 ==========
# plot Kaplan-Meier curves of top variables
{
par(mfrow=c(3,3))
dataset<- roi.stage1
# grade=1  
{x<- survfit(Surv(days.surgery.pfs/365, ind.lung.pfs) ~ 
    grade1,
    data = dataset, na.action=na.exclude,
    conf.int = 0.95, conf.type="log-log")
x
#            n events median 0.95LCL 0.95UCL
# grade1=0 239     88     NA      NA      NA
# grade1=1  49      7     NA      NA      NA
# use ggsurvplot in survminer package to 
# find number at risk n1 and n2
n1<- as.character(c(49,44,42,41,40,32,7))
n2<- as.character(c(239,195,159,144,124,88,36))
plot(x,conf.int=F, lty=rep(1,2),
 log=FALSE, col=c("red","black"),
	main="Regional highest grade (sample 1)", mark.time=TRUE,
	xlab= "Years since initial surgery",
	ylab="Progression free survival",
	cex.lab=1,cex.main=1,cex.axis=1, lwd=3)
legend(1,0.6,c("Grade=1",  "Grade>1"),horiz = T,
	lty=rep(1,2),col=c("black", "red"),#lwd=3,
   pch=rep(15,4),	box.lty=0, cex=0.8)
text(6,0.38, "HR=3.135 (p=0.0036)")#1/0.319=3.135
text(2, 0.2, "Number at risk:", cex=0.8)
for (i in 1:length(n1))
{
  text(2*(i-1), 0.10, n1[i], col="black", cex=0.8)
  text(2*(i-1), 0.03, n2[i], col="red", cex=0.8)
}

}

# not AIS
{dataset$nonAIS<- 1-dataset$AIS
  x<- survfit(Surv(days.surgery.pfs/365, ind.lung.pfs) ~ 
    nonAIS,
    data = dataset, na.action=na.exclude,
    conf.int = 0.95, conf.type="log-log")
x
#            n events median 0.95LCL 0.95UCL
# nonAIS=0  93     19     NA      NA      NA
# nonAIS=1 195     76     NA    10.3      NA
# use ggsurvplot in survminer package to 
# find number at risk n1 and n2
n1<- as.character(c(93,81,75,71,66,50,16))
n2<- as.character(c(195, 158,126,114,98,70,27))
plot(x,conf.int=F, lty=rep(1,2),
 log=FALSE, col=c("black", "red"),
	main="AIS (sample 1)", mark.time=TRUE,
	xlab= "Years since initial surgery",
	ylab="Progression free survival",
	cex.lab=1,cex.main=1,cex.axis=1, lwd=3)
legend("topright",c("AIS",  "not AIS"),horiz = T,
	lty=rep(1,2),col=c("black", "red"),#lwd=3,
   pch=rep(15,4),	box.lty=0, cex=0.8)
text(5,0.4, "HR=2.208 (p=0.0020)")#1/0.453=2.208
text(2, 0.2, "Number at risk:", cex=0.8)
for (i in 1:length(n1))
{
  text(2*(i-1), 0.10, n1[i], col="black", cex=0.8)
  text(2*(i-1), 0.03, n2[i], col="red", cex=0.8)
}

}

# tumorD
{x<- survfit(Surv(days.surgery.pfs/365, ind.lung.pfs) ~ 
    (tumorD>13.25),
    data = dataset, na.action=na.exclude,
    conf.int = 0.95, conf.type="log-log")
x
#                        n events median 0.95LCL 0.95UCL
# tumorD > 13.25=FALSE 148     39     NA      NA      NA
# tumorD > 13.25=TRUE  140     56     NA    10.2      NA
# use ggsurvplot in survminer package to 
# find number at risk n1 and n2
n1<- as.character(c(148,128,116,108,96,70,21))
n2<- as.character(c(140,111,85,77,68,50,22))
plot(x,conf.int=F, lty=rep(1,2),
 log=FALSE, col=c("black", "red"),
	main="Largest invasive tumor size (sample 1)", mark.time=TRUE,
	xlab= "Years since initial surgery",
	ylab="Progression free survival",
	cex.lab=1,cex.main=1,cex.axis=1, lwd=3)
legend("topright",c("<13.2mm",  ">=13.2mm"),horiz = T,
	lty=rep(1,2),col=c("black", "red"),#lwd=3,
   pch=rep(15,4),	cex=0.8,box.lty=0)
text(5,0.4, "HR=1.738 (p=0.0081)")
text(2, 0.2, "Number at risk:", cex=0.8)
for (i in 1:length(n1))
{
  text(2*(i-1), 0.10, n1[i], col="black", cex=0.8)
  text(2*(i-1), 0.03, n2[i], col="red", cex=0.8)
}
}

# age.surgery 
{x<- survfit(Surv(days.surgery.pfs/365, ind.lung.pfs) ~ 
    (age.surgery>65.5192),
    data = dataset, na.action=na.exclude,
    conf.int = 0.95, conf.type="log-log")
x
#                               n events median 0.95LCL 0.95UCL
# age.surgery > 65.5192=FALSE 148     41     NA      NA      NA
# age.surgery > 65.5192=TRUE  140     54     NA    10.3      NA
# use ggsurvplot in survminer package to 
# find number at risk n1 and n2
n1<- as.character(c(148,128,113,109,95,74,27))
n2<- as.character(c(140,111,88,76,69,46,16))
plot(x,conf.int=F, lty=rep(1,2),
 log=FALSE, col=c("black", "red"),
	main="Age at surgery (sample 1)", mark.time=TRUE,
	xlab= "Years since initial surgery",
	ylab="Progression free survival",
	cex.lab=1,cex.main=1,cex.axis=1, lwd=3)
legend("topright",c("<=65",  ">65"),horiz = T,
	lty=rep(1,2),col=c("black", "red"),#lwd=3,
   pch=rep(15,4),	cex=0.8,box.lty=0)
text(5,0.4, "HR=1.611 (p=0.0216)")
text(2, 0.2, "Number at risk:", cex=0.8)
for (i in 1:length(n1))
{
  text(2*(i-1), 0.10, n1[i], col="black", cex=0.8)
  text(2*(i-1), 0.03, n2[i], col="red", cex=0.8)
}
}

# Pneumonia
{x<- survfit(Surv(days.surgery.pfs/365, ind.lung.pfs) ~ 
    Pneumonia,
    data = dataset, na.action=na.exclude,
    conf.int = 0.95, conf.type="log-log")
x
#               n events median 0.95LCL 0.95UCL
# Pneumonia=0 256     77     NA      NA      NA
# Pneumonia=1  32     18   9.36    2.44      NA
# use ggsurvplot in survminer package to 
# find number at risk n1 and n2
n1<- as.character(c(256,216,185,170,152,113,41))
n2<- as.character(c(32,23,16,15,12,7,2))
plot(x,conf.int=F, lty=rep(1,2),
 log=FALSE, col=c("black", "red"),
	main="Pneumonia (sample 1)", mark.time=TRUE,
	xlab= "Years since initial surgery",
	ylab="Progression free survival",
	cex.lab=1,cex.main=1,cex.axis=1, lwd=3)
legend("topright",c("Not present", "present"),
  horiz = T, lty=rep(1,2),col=c("black", "red"),#lwd=3,
   pch=rep(15,4),	cex=0.8,box.lty=0)
text(5,0.4, "HR=2.397 (p=0.0009)")
text(2, 0.2, "Number at risk:", cex=0.8)
for (i in 1:length(n1))
{
  text(2*(i-1), 0.10, n1[i], col="black", cex=0.8)
  text(2*(i-1), 0.03, n2[i], col="red", cex=0.8)
}
}

# cigsmok
{x<- survfit(Surv(days.surgery.pfs/365, ind.lung.pfs) ~ 
    cigsmok,
    data = dataset, na.action=na.exclude,
    conf.int = 0.95, conf.type="log-log")
x
#             n events median 0.95LCL 0.95UCL
# cigsmok=0 128     34     NA      NA      NA
# cigsmok=1 160     61     NA    11.2      NA
# use ggsurvplot in survminer package to 
# find number at risk n1 and n2
n1<- as.character(c(128,109,94,90,82,56,18))
n2<- as.character(c(160,130,107,95,82,64,25))
plot(x,conf.int=F, lty=rep(1,2),
 log=FALSE, col=c("black", "red"),
	main="Smoking status (sample 1)", mark.time=TRUE,
	xlab= "Years since initial surgery",
	ylab="Progression free survival",
	cex.lab=1,cex.main=1,cex.axis=1, lwd=3)
legend("topright",c("Former",  "Current"),horiz = T,
	lty=rep(1,2),col=c("black", "red"),#lwd=3,
   pch=rep(15,4),	cex=0.8,box.lty=0)
text(5,0.4, "HR=1.565 (p=0.0365)")
text(2, 0.2, "Number at risk:", cex=0.8)
for (i in 1:length(n1))
{
  text(2*(i-1), 0.10, n1[i], col="black", cex=0.8)
  text(2*(i-1), 0.03, n2[i], col="red", cex=0.8)
}
}

# pre.sq.is
{x<- survfit(Surv(days.surgery.pfs/365, ind.lung.pfs) ~ 
    pre.sq.is,
    data = dataset, na.action=na.exclude,
    conf.int = 0.95, conf.type="log-log")
x
#               n events median 0.95LCL 0.95UCL
# pre.sq.is=0 267     86     NA      NA      NA
# pre.sq.is=1  21      9     NA    1.52      NA
# use ggsurvplot in survminer package to 
# find number at risk n1 and n2
n1<- as.character(c(267,227,191,176,157,115,43))
n2<- as.character(c(21,12,10,9,7,5,0))
plot(x,conf.int=F, lty=rep(1,2),
 log=FALSE, col=c("black", "red"),
	main="Squamous in situ in 
premalignant region (sample 1)", mark.time=TRUE,
	xlab= "Years since initial surgery",
	ylab="Progression free survival",
	cex.lab=1,cex.main=1,cex.axis=1, lwd=3)
legend("topright",c("Not present",  "present"),horiz = T,
	lty=rep(1,2),col=c("black", "red"),#lwd=3,
   pch=rep(15,4),	cex=1,box.lty=0)
text(5,0.4, "HR=1.835 (p=0.0836)")
text(2, 0.2, "Number at risk:", cex=0.8)
for (i in 1:length(n1))
{
  text(2*(i-1), 0.10, n1[i], col="black", cex=0.8)
  text(2*(i-1), 0.03, n2[i], col="red", cex=0.8)
}
}


# pkyr
{x<- survfit(Surv(days.surgery.pfs/365, ind.lung.pfs) ~ 
    (pkyr>139),
    data = dataset, na.action=na.exclude,
    conf.int = 0.95, conf.type="log-log")
x
#                    n events median 0.95LCL 0.95UCL
# pkyr > 139=FALSE 279     89     NA      NA      NA
# pkyr > 139=TRUE    9      6   4.07    1.05      NA
# use ggsurvplot in survminer package to 
# find number at risk n1 and n2
n1<- as.character(c(279,232,196,183,162,118,43))
n2<- as.character(c(9,7,5,2,2,2,0))
plot(x,conf.int=F, lty=rep(1,2),
 log=FALSE, col=c("black", "red"),
	main="Smoking pack-years (sample 1)", mark.time=TRUE,
	xlab= "Years since initial surgery",
	ylab="Progression free survival",
	cex.lab=1,cex.main=1,cex.axis=1, lwd=3)
legend("topright",c("<=139",  ">139"),horiz = T,
	lty=rep(1,2),col=c("black", "red"),#lwd=3,
   pch=rep(15,4),	cex=0.8,box.lty=0)
text(9,0.5, "HR=2.627 (p=0.0224)")
text(2, 0.2, "Number at risk:", cex=0.8)
for (i in 1:length(n1))
{
  text(2*(i-1), 0.10, n1[i], col="black", cex=0.8)
  text(2*(i-1), 0.03, n2[i], col="red", cex=0.8)
}

}
par(mfrow=c(1,1))
}


# Figure S2 ========

{
par(mfrow=c(3,3))
dataset<- texture.roi
# days.to.surgery 
{x<- survfit(Surv(days.surgery.pfs/365, ind.lung.pfs) ~ 
    (days.to.surgery>88),
    data = dataset, na.action=na.exclude,
    conf.int = 0.95, conf.type="log-log")
x
#                              n events median 0.95LCL 0.95UCL
# days.to.surgery > 88=FALSE  81     14     NA      NA      NA
# days.to.surgery > 88=TRUE  101     40     NA    6.75      NA
# use ggsurvplot in survminer package to 
# find number at risk n1 and n2
n1<- as.character(c(81,72,67,64,57,53,22))
n2<- as.character(c(101,83,60,53,43,31,11))
plot(x,conf.int=F, lty=rep(1,2),
 log=FALSE, col=c("black", "red"),
	main="Days since LDCT screening date to 
  the surgery date (sample 2)", mark.time=TRUE,
	xlab= "Years since initial surgery",
	ylab="Progression free survival",
	cex.lab=1,cex.main=1,cex.axis=1, lwd=3)
legend("topright",c("<88 days ",  ">88 days"),horiz = T,
	lty=rep(1,2),col=c("black", "red"),#lwd=3,
   pch=rep(15,4),	cex=0.8,box.lty=0)
text(4,0.4, "HR=2.736 (p=0.0012)")
text(2, 0.2, "Number at risk:", cex=0.8)
for (i in 1:length(n1))
{
  text(2*(i-1), 0.10, n1[i], col="black", cex=0.8)
  text(2*(i-1), 0.03, n2[i], col="red", cex=0.8)
}

x<- coxph(Surv(days.surgery.pfs, ind.lung.pfs) ~ 
    (days.to.surgery>88)+staging+grade+ AIS+
    age.surgery+loc1.sd.v, 
     data = texture.roi, na.action = na.exclude, 
     method = "efron", robust = F)
summary(x)
# n= 182, number of events= 54 
#                              coef exp(coef) se(coef)      z Pr(>|z|)   
# days.to.surgery > 88TRUE  0.88069   2.41257  0.32104  2.743  0.00608 **
# stagingT1b               -0.20286   0.81639  0.32174 -0.631  0.52835   
# stagingT1c               -0.61150   0.54254  0.53101 -1.152  0.24949   
# grade                     0.08319   1.08675  0.18575  0.448  0.65425   
# AIS                      -0.87778   0.41570  0.35479 -2.474  0.01336 * 
# age.surgery               0.05215   1.05354  0.02911  1.792  0.07317 . 
# loc1.sd.v                 1.98149   7.25357  0.77918  2.543  0.01099 * 
}

# grade1
{x<- survfit(Surv(days.surgery.pfs/365, ind.lung.pfs) ~ 
    grade1,
    data = dataset, na.action=na.exclude,
    conf.int = 0.95, conf.type="log-log")
x
#            n events median 0.95LCL 0.95UCL
# grade1=0 144     50     NA      NA      NA
# grade1=1  38      4     NA      NA      NA
# use ggsurvplot in survminer package to 
# find number at risk n1 and n2
n2<- as.character(c(144,120,93,84,68,56,26))
n1<- as.character(c(38,35,34,33,32,28,7))
plot(x,conf.int=F, lty=rep(1,2),
 log=FALSE, col=c("red","black"),
	main="Regional highest grade (sample 2)", mark.time=TRUE,
	xlab= "Years since initial surgery",
	ylab="Progression free survival",
	cex.lab=1,cex.main=1,cex.axis=1, lwd=3)
legend(1,0.6,c("Grade=1",  "Grade>1"),horiz = T,
	lty=rep(1,2),col=c("black", "red"),#lwd=3,
   pch=rep(15,4),	box.lty=0, cex=0.8)
text(6,0.38, "HR=4.036 (p=0.0073)")#1/0.2478=4.036
text(2, 0.2, "Number at risk:", cex=0.8)
for (i in 1:length(n1))
{
  text(2*(i-1), 0.10, n1[i], col="black", cex=0.8)
  text(2*(i-1), 0.03, n2[i], col="red", cex=0.8)
}

}


# loc1.sd.v
{x<- survfit(Surv(days.surgery.pfs/365, ind.lung.pfs) ~ 
    (loc1.sd.v>0.9734),
    data = dataset, na.action=na.exclude,
    conf.int = 0.95, conf.type="log-log")
x
#                            n events median 0.95LCL 0.95UCL
# loc1.sd.v > 0.9734=FALSE 130     30     NA      NA      NA
# loc1.sd.v > 0.9734=TRUE   52     24     NA    3.57      NA
# use ggsurvplot in survminer package to 
# find number at risk n1 and n2
n1<- as.character(c(130,113,95,89,75,63,26))
n2<- as.character(c(52,42,32,28,25,21,7))
plot(x,conf.int=F, lty=rep(1,2),
 log=FALSE, col=c("black", "red"),
	main="LocSd100V (sample 2)", mark.time=TRUE,
	xlab= "Years since initial surgery",
	ylab="Progression free survival",
	cex.lab=1,cex.main=1,cex.axis=1, lwd=3)
legend("topright",c("<=0.9734 ",  ">0.9734"),horiz = T,
	lty=rep(1,2),col=c("black", "red"),#lwd=3,
   pch=rep(15,4),	cex=0.8,box.lty=0)
text(4,0.4, "HR=2.226 (p=0.0035)")
text(2, 0.2, "Number at risk:", cex=0.8)
for (i in 1:length(n1))
{
  text(2*(i-1), 0.10, n1[i], col="black", cex=0.8)
  text(2*(i-1), 0.03, n2[i], col="red", cex=0.8)
}
}

# rms1
{x<- survfit(Surv(days.surgery.pfs/365, ind.lung.pfs) ~ 
    (rms1>1.0762),
    data = dataset, na.action=na.exclude,
    conf.int = 0.95, conf.type="log-log")
x
#                       n events median 0.95LCL 0.95UCL
# rms1 > 1.0762=FALSE 173     48     NA      NA      NA
# rms1 > 1.0762=TRUE    9      6   2.44   0.249      NA
# use ggsurvplot in survminer package to 
# find number at risk n1 and n2
n1<- as.character(c(173,150,125,115,98,82,31))
n2<- as.character(c(9,5,2,2,2,2,2))
plot(x,conf.int=F, lty=rep(1,2),
 log=FALSE, col=c("black", "red"),
	main="Root mean square (sample 2)", mark.time=TRUE,
	xlab= "Years since initial surgery",
	ylab="Progression free survival",
	cex.lab=1,cex.main=1,cex.axis=1, lwd=3)
legend("topright",c("Low ",  "High"),horiz = T,
	lty=rep(1,2),col=c("black", "red"),#lwd=3,
   pch=rep(15,4),	cex=0.8,box.lty=0)
text(8,0.5, "HR=4.048 (p=0.0013)")
text(2, 0.2, "Number at risk:", cex=0.8)
for (i in 1:length(n1))
{
  text(2*(i-1), 0.10, n1[i], col="black", cex=0.8)
  text(2*(i-1), 0.03, n2[i], col="red", cex=0.8)
}
}


# age.surgery
{x<- survfit(Surv(days.surgery.pfs/365, ind.lung.pfs) ~ 
    (age.surgery>65.4699),
    data = dataset, na.action=na.exclude,
    conf.int = 0.95, conf.type="log-log")
x
#                              n events median 0.95LCL 0.95UCL
# age.surgery > 65.4699=FALSE 99     22     NA      NA      NA
# age.surgery > 65.4699=TRUE  83     32     NA    6.28      NA
# use ggsurvplot in survminer package to 
# find number at risk n1 and n2
n1<- as.character(c(99,90,79,76,64,55,20))
n2<- as.character(c(83,65,48,41,36,29,13))
plot(x,conf.int=F, lty=rep(1,2),
 log=FALSE, col=c("black", "red"),
	main="Age at surgery (sample 2)", mark.time=TRUE,
	xlab= "Years since initial surgery",
	ylab="Progression free survival",
	cex.lab=1,cex.main=1,cex.axis=1, lwd=3)
legend("topright",c("<=65",  ">65"),horiz = T,
	lty=rep(1,2),col=c("black", "red"),#lwd=3,
   pch=rep(15,4),	cex=0.8,box.lty=0)
text(6,0.4, "HR=2.106 (p=0.0073)")
text(2, 0.2, "Number at risk:", cex=0.8)
for (i in 1:length(n1))
{
  text(2*(i-1), 0.10, n1[i], col="black", cex=0.8)
  text(2*(i-1), 0.03, n2[i], col="red", cex=0.8)
}
}


# glrlmaskewness 
{x<- survfit(Surv(days.surgery.pfs/365, ind.lung.pfs) ~ 
    (glrlmaskewness >1.6533),
    data = dataset, na.action=na.exclude,
    conf.int = 0.95, conf.type="log-log")
x
#                                 n events median 0.95LCL 0.95UCL
# glrlmaskewness > 1.6533=FALSE  42      5     NA      NA      NA
# glrlmaskewness > 1.6533=TRUE  140     49     NA      NA      NA
# use ggsurvplot in survminer package to 
# find number at risk n1 and n2
n1<- as.character(c(42,39,35,35,32,25,10))
n2<- as.character(c(140,116,92,82,68,59,23))
plot(x,conf.int=F, lty=rep(1,2),
 log=FALSE, col=c("black", "red"),
	main="GLRLM skewness (sample 2)", mark.time=TRUE,
	xlab= "Years since initial surgery",
	ylab="Progression free survival",
	cex.lab=1,cex.main=1,cex.axis=1, lwd=3)
legend("topright",c("<=1.6533 ",  ">1.6533"),horiz = T,
	lty=rep(1,2),col=c("black", "red"),#lwd=3,
   pch=rep(15,4),	cex=0.8,box.lty=0)
text(4,0.5, "HR=3.490 (p=0.0078)")
text(2, 0.2, "Number at risk:", cex=0.8)
for (i in 1:length(n1))
{
  text(2*(i-1), 0.10, n1[i], col="black", cex=0.8)
  text(2*(i-1), 0.03, n2[i], col="red", cex=0.8)
}
}

# glrlmakurtosis
{x<- survfit(Surv(days.surgery.pfs/365, ind.lung.pfs) ~ 
    (glrlmakurtosis >1.9048),
    data = dataset, na.action=na.exclude,
    conf.int = 0.95, conf.type="log-log")
x
#                                 n events median 0.95LCL 0.95UCL
# glrlmakurtosis > 1.9048=FALSE  42      6     NA      NA      NA
# glrlmakurtosis > 1.9048=TRUE  140     48     NA      NA      NA
# use ggsurvplot in survminer package to 
# find number at risk n1 and n2
n1<- as.character(c(42,39,33,32,30,24,12))
n2<- as.character(c(140,116,94,85,70,60,21))
plot(x,conf.int=F, lty=rep(1,2),
 log=FALSE, col=c("black", "red"),
	main="GLRLM kurtosis (sample 2)", mark.time=TRUE,
	xlab= "Years since initial surgery",
	ylab="Progression free survival",
	cex.lab=1,cex.main=1,cex.axis=1, lwd=3)
legend("topright",c("<=1.9048 ",  ">1.9048"),horiz = T,
	lty=rep(1,2),col=c("black", "red"),#lwd=3,
   pch=rep(15,4),	cex=0.8,box.lty=0)
text(4,0.5, "HR=2.737 0.0201)")
text(2, 0.2, "Number at risk:", cex=0.8)
for (i in 1:length(n1))
{
  text(2*(i-1), 0.10, n1[i], col="black", cex=0.8)
  text(2*(i-1), 0.03, n2[i], col="red", cex=0.8)
}
}


# tumorD
{x<- survfit(Surv(days.surgery.pfs/365, ind.lung.pfs) ~ 
    (tumorD>12.3),
    data = dataset, na.action=na.exclude,
    conf.int = 0.95, conf.type="log-log")
x
#                      n events median 0.95LCL 0.95UCL
# tumorD > 12.3=FALSE 97     22     NA      NA      NA
# tumorD > 12.3=TRUE  85     32     NA    6.98      NA
# use ggsurvplot in survminer package to 
# find number at risk n1 and n2
n1<- as.character(c(97,84,78,73,64,54,18))
n2<- as.character(c(85,71,49,44,36,30,15))
plot(x,conf.int=F, lty=rep(1,2),
 log=FALSE, col=c("black", "red"),
	main="Largest invasive tumor size (sample 2)", mark.time=TRUE,
	xlab= "Years since initial surgery",
	ylab="Progression free survival",
	cex.lab=1,cex.main=1,cex.axis=1, lwd=3)
legend("topright",c("<12.3mm",  ">=12.3mm"),horiz = T,
	lty=rep(1,2),col=c("black", "red"),#lwd=3,
   pch=rep(15,4),	cex=0.8,box.lty=0)
text(6,0.4, "HR=1.913 (p=0.0196)")
text(2, 0.2, "Number at risk:", cex=0.8)
for (i in 1:length(n1))
{
  text(2*(i-1), 0.10, n1[i], col="black", cex=0.8)
  text(2*(i-1), 0.03, n2[i], col="red", cex=0.8)
}
}


# AIS
{x<- survfit(Surv(days.surgery.pfs/365, ind.lung.pfs) ~ 
    AIS,
    data = dataset, na.action=na.exclude,
    conf.int = 0.95, conf.type="log-log")
x
#         n events median 0.95LCL 0.95UCL
# AIS=0 113     41     NA      NA      NA
# AIS=1  69     13     NA      NA      NA
# use ggsurvplot in survminer package to 
# find number at risk n1 and n2
n2<- as.character(c(113,94,71,64,51,42,21))
n1<- as.character(c(69,61,56,53,49,42,12))
plot(x,conf.int=F, lty=rep(1,2),
 log=FALSE, col=c("red","black"),
	main="AIS (sample 2)", mark.time=TRUE,
	xlab= "Years since initial surgery",
	ylab="Progression free survival",
	cex.lab=1,cex.main=1,cex.axis=1, lwd=3)
legend("topright",c("AIS", "Not AIS"),horiz = T,
	lty=rep(1,2),col=c("black", "red"),#lwd=3,
   pch=rep(15,4),	box.lty=0, cex=0.8)
text(5,0.4, "HR=2.237 (p=0.0114)")#1/0.447=2.237
text(2, 0.2, "Number at risk:", cex=0.8)
for (i in 1:length(n1))
{
  text(2*(i-1), 0.10, n1[i], col="black", cex=0.8)
  text(2*(i-1), 0.03, n2[i], col="red", cex=0.8)
}

}

par(mfrow=c(1,1))
}



# Figure S3 ========

# rmst.ct prediction 
{par(mfrow=c(1,3))
# 5- and 10-year AUC of rmst.ct in texture.roi
{dataset<- texture.roi
set.seed(123)
 y5.ct<- survivalROC(Stime=dataset$days.pfs,
		status=dataset$ind.lung.pfs, marker=-dataset$rmst.ct,
		predict.time=5*365, method="KM")
 y10.ct<- survivalROC(Stime=dataset$days.pfs,
		status=dataset$ind.lung.pfs, marker=-dataset$rmst.ct,
		predict.time=10*365, method="KM")

plot(y5.ct$FP, y5.ct$TP, type="l", xlim=c(0,1), ylim=c(0,1),
	xlab="1-specificity", ylab="sensitivity", 
  col="black",lwd=2,
	cex.main=1, cex.lab=1,cex.axis=1,
  main="Time-dependent ROC from\nCT features (sample 2)")
lines(y10.ct$FP, y10.ct$TP,lwd=2, col="red") 
lines(c(0,1),c(0, 1), lty=2, col="grey")
#legend(0.2,0.5,
#paste("AUC=", round(y$AUC,3), sep=""))
text(0.7,0.25,paste0("5-year AUC=", round(y5.ct$AUC,3)),
  col="black",cex=0.9)
text(0.7,0.15,paste0("10-year AUC=", round(y10.ct$AUC,3)),
  col="red",cex=0.9)

print(paste0("5-year rmst.ct AUC=", round(y5.ct$AUC,4)))
print(paste0("10-year rmst.ct AUC=", round(y10.ct$AUC,4)))
# "5-year rmst.ct AUC=0.6733"
# "10-year rmst.ct AUC=0.632"
}

# 5- and 10-year time dependent 
# PPV, NPV conditioning on sensitivity
# of rmst.ct using texture.roi
{
dataset<- texture.roi
set.seed(123)
y5.ct<- survivalROC(Stime=dataset$days.pfs/365,
		status=dataset$ind.lung.pfs, marker= -dataset$rmst.ct,
		predict.time=5, method="NNE", span = 0.07)
set.seed(123)
y10.ct<- survivalROC(Stime=dataset$days.pfs/365,
		status=dataset$ind.lung.pfs, marker= -dataset$rmst.ct,
		predict.time=10, method="NNE", span = 0.07)
# 5-year
results<- PPV.NPV(y5.ct, -dataset$rmst.ct)
results<- results[-1,]
results<- results[with(results, order(sensitivity)),]
results<- results[results$sensitivity>0.6,]
results<- results[results$sensitivity<0.95,]
results.5<- results[with(results, order(sensitivity)),]

# 10-year
results<- PPV.NPV(y10.ct, -dataset$rmst.ct)
results<- results[-1,]
results<- results[with(results, order(sensitivity)),]
results<- results[results$sensitivity>0.6,]
results<- results[results$sensitivity<0.95,]
results.10<- results[with(results, order(sensitivity)),]

# PPV
plot(results.5$sensitivity, results.5$PPV,
    main="Time-dependent positive predictive
   value from CT features (sample 2)", ylim=c(0.25,0.45), type="l", lwd=2,
  xlab="Sensitivity", ylab="", col="black")
lines(results.10$sensitivity, results.10$PPV,
  col="red",lwd=2)
legend("topright", c("year 5","year 10"),
   pch=rep(15,3), col=c("black", "red"),
  horiz = TRUE, box.lty = 0)

# NPV
plot(results.5$sensitivity, results.5$NPV,
    main="Time-dependent negative predictive
   value from CT features (sample 2)", ylim=c(0.72,1), type="l", lwd=2,
  xlab="Sensitivity", ylab="", col="black")
lines(results.10$sensitivity, results.10$NPV,
  col="red",lwd=2)
legend("topright", c("year 5","year 10"),
   pch=rep(15,3), col=c("black", "red"),
  horiz = TRUE, box.lty = 0)


}
par(mfrow=c(1,1))
}


# Table 2 ========

# 4 Cox PH regression models that adjust 
# for age and treatment
{
roi.stage1$grade34<- 1*(roi.stage1$grade>=3)
texture.roi$grade34<- 1*(texture.roi$grade>=3)

# use binary rmst.roi and rmst.roiCT high subgroups. 
# The threshold values are the same as used in 
# Kaplan-Meier curves:
roi.stage1$IDLE1h<- 1*(roi.stage1$IDLE1>=0.328)
texture.roi$IDLE1h<- 1*(texture.roi$IDLE1>=0.328)
texture.roi$IDLE2h<- 1*(texture.roi$IDLE2>=0.271)


x<- coxph(Surv(days.surgery.pfs/365, ind.lung.pfs) ~ 
   IDLE1h+staging+grade34+age.dx+chemo+radiotherapy,
	data = roi.stage1, na.action = na.exclude, method = "efron", robust = F)
x<- summary(x)$coefficients
round(x,4)
#                 coef exp(coef) se(coef)       z Pr(>|z|)
# IDLE1h        1.2541    3.5046   0.2631  4.7672   0.0000
# stagingT1b   -0.0039    0.9961   0.2396 -0.0162   0.9870
# stagingT1c   -0.1933    0.8242   0.3318 -0.5826   0.5602
# grade34       0.0144    1.0145   0.2206  0.0651   0.9481
# age.dx        0.0145    1.0146   0.0200  0.7263   0.4676
# chemo        -0.1347    0.8739   0.2935 -0.4590   0.6462
# radiotherapy  0.6956    2.0049   0.4232  1.6437   0.1002


# model 2
x<- coxph(Surv(days.surgery.pfs/365, ind.lung.pfs) ~ 
   IDLE1h+staging+grade34+age.dx+chemo+radiotherapy,
	data = texture.roi, na.action = na.exclude, method = "efron", robust = F)
x<- summary(x)$coefficients
round(x,4)
#                 coef exp(coef) se(coef)       z Pr(>|z|)
# IDLE1h        1.3357    3.8028   0.3729  3.5824   0.0003
# stagingT1b    0.0004    1.0004   0.3248  0.0013   0.9990
# stagingT1c   -0.3470    0.7068   0.5334 -0.6506   0.5153
# grade34       0.0936    1.0982   0.3073  0.3048   0.7605
# age.dx        0.0446    1.0456   0.0286  1.5587   0.1191
# chemo        -0.4427    0.6423   0.4335 -1.0212   0.3072
# radiotherapy  0.9376    2.5538   0.6302  1.4879   0.1368



# model 3
x<- coxph(Surv(days.surgery.pfs/365, ind.lung.pfs) ~ 
   IDLE2h+staging+grade34+age.dx+chemo+radiotherapy,
	data = texture.roi, na.action = na.exclude, method = "efron", robust = F)
x<- summary(x)$coefficients
round(x,4)
#                 coef exp(coef) se(coef)       z Pr(>|z|)
# IDLE2h        1.7353    5.6708   0.2975  5.8321   0.0000
# stagingT1b   -0.1433    0.8665   0.3157 -0.4539   0.6499
# stagingT1c   -0.2603    0.7708   0.5506 -0.4727   0.6364
# grade34      -0.0184    0.9818   0.3013 -0.0610   0.9513
# age.dx        0.0313    1.0318   0.0295  1.0609   0.2888
# chemo        -0.4005    0.6700   0.4440 -0.9020   0.3671
# radiotherapy  0.3335    1.3959   0.6296  0.5297   0.5963



x<- coxph(Surv(days.surgery.pfs/365, ind.lung.pfs) ~ 
   IDLE1h+IDLE2h+staging+grade34+age.dx+chemo+radiotherapy,
	data = texture.roi, na.action = na.exclude, method = "efron", robust = F)
x<- summary(x)$coefficients
round(x,4)
#                 coef exp(coef) se(coef)       z Pr(>|z|)
# IDLE1h        0.7746    2.1697   0.3949  1.9617   0.0498
# IDLE2h        1.6251    5.0792   0.3115  5.2173   0.0000
# stagingT1b   -0.2653    0.7670   0.3294 -0.8052   0.4207
# stagingT1c   -0.4111    0.6629   0.5543 -0.7417   0.4582
# grade34      -0.0756    0.9272   0.3041 -0.2487   0.8036
# age.dx        0.0155    1.0156   0.0312  0.4975   0.6188
# chemo        -0.4629    0.6295   0.4427 -1.0455   0.2958
# radiotherapy  0.3634    1.4382   0.6352  0.5721   0.5673

}


