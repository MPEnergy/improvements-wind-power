rm(list = ls())
library(R.matlab)
library(chron)
library(glmnet)
library(forecast)
library("astsa")
library(lmtest)
require(rugarch)
library("fGarch")
require(tSeries)
library("ggpubr")

data <- readMat('Hourly_residuals_sigmaGARCH.mat')



d<-data$dat
newd<-strftime(as.POSIXct( ( d - 719529) * 86400,origin = "1970-01-01",tz = "UTC"), format = '%Y-%m-%d %H:%M', 
               tz = 'UTC', usetz = FALSE)
library(lubridate)
timeDate <- as.POSIXct(newd,tz="UTC")
v <- data$val

xh<-format(strptime(timeDate,format="%Y-%m-%d %H:%M:%S"),'%H')
xm<-format(strptime(timeDate,format="%Y-%m-%d %H:%M:%S"),'%M')

xh<-as.numeric(xh)
xm<-as.numeric(xm)

idx = xm==59
timeDate[idx]=timeDate[idx]+minutes(1)


st<-20
v1 <- v[,st]
#acf2(v1^2,main=paste('Squared Residuals in Point', st))
#dev.copy(jpeg,paste("C:/Users/u6032456/Documents/Thesis/PhDMeeting/article3/Step4_ResidualsBSB/P_", st, "_var2_sigmaBSB.jpg"))
#dev.off()


modm11 = ugarchspec(variance.model=list(model="gjrGARCH", garchOrder=c(1,1),submodel=NULL), 
                    mean.model=list(armaOrder=c(1,0), archm=F, archpow=2, 
                                    include.mean=F),distribution.model = "std")

#spec1 <- ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = TRUE),distribution.model = "nig")

gfit <- ugarchfit(spec = modm11, data = coredata(v1), solver = "hybrid")


#gfit <- garchFit(formula = ~garch(1,1), data=coredata(v1), trace=FALSE, cond.dist="QMLE")
#summary(gfit)
#xi=v1/gfit@h.t

r = residuals(gfit, standardize=TRUE)

acf2(coredata(r))
dev.copy(jpeg,paste("C:/article3/Step5_ResidualsGARCH/P_", st, "_var_GARCH1_1.jpg"))
dev.off()



# Make qq plots of the residuals

isnormal<-matrix(nrow=85,ncol=1)
for (st in 1:85)
{
  v1 <- v[,st]
  gfit <- garchFit(formula = ~arma(1,0)+garch(1,1), data=coredata(v1), trace=FALSE, cond.dist="QMLE")
  r = residuals(gfit, standardize=TRUE)
  t=ks.test(r,"pnorm", mean=mean(r), sd=sd(r))
  isnormal[st]<-t$p.value>=0.05
}




params <- matrix(nrow=85, ncol=7)
pval<- matrix(nrow=85, ncol=7)
R <- matrix(nrow=38664, ncol=85)
for(st in 1:85){
  v1 <- v[,st]
  #gfit <- garchFit(formula = ~arma(1,0)+garch(1,1), data=v1, trace=FALSE)
  gfit <- ugarchfit(spec = modm11, data = coredata(v1), solver = "hybrid")
  r = residuals(gfit, standardize=TRUE)
  R[,st]<-r
  #pval[st,] <- gfit@fit$coef[,4]
  params[st,]<- gfit@fit$coef
  #acf2(coredata(r),main=paste('Residuals in Point', st))
  #dev.copy(jpeg,paste("C:/article3/Step5_ResidualsGARCH/P_", st, "_var_AR1_GFJGARCH1_1.jpg"))
  #dev.off()
  #acf2(coredata(r)^2,main=paste('Squared Residuals in Point', st))
  #dev.copy(jpeg,paste("C:/article3/Step5_ResidualsGARCH/P_", st, "_var2_AR1_GJRGARCH1_1.jpg"))
  #dev.off()
  #ggqqplot(r,main=paste("Point", st))
  #dev.copy(jpeg,paste("C:/article3/Step5_ResidualsGARCH/P_", st, "_QQplot.jpg"))
  #dev.off()
  
  #ggdensity(r, fill = "lightgray")
  #dev.copy(jpeg,paste("C:/article3/Step5_ResidualsGARCH/P_", st, "_Density.jpg"))
  #dev.off()
}

library("xlsx")
write.csv(xi, file = "C:/article3/Residuals_AR(1)_GJRGARCH(1,1).csv")




# Make QQ plots for t-student distribution of residuals
library(fitdistrplus)
library(metRology)
library(ggplot2)
library(qqplotr)

for (st in 3:85){
  X<-R[,st]
  #par.est = fitdist(as.numeric(X), distr = "t.scaled", start = list(mean = mean(X), sd = sd(X), df = 5.780090))
  
  my_title <- paste("Q-Q plot of residuals in Point",st)
  #ggplot(mapping = aes(sample = X)) +
  #  stat_qq_band(distribution = "t.scaled", dparams = par.est$estimate) +
  #  stat_qq_point(distribution = "t.scaled", dparams = par.est$estimate, size = 2) +  stat_qq_line(distribution = "t.scaled", dparams = par.est$estimate) +
  #  xlab("Theoretical Quantiles") + ylab("Sample Quantiles") +
  #  ggtitle(my_title) +
  #  theme(plot.title = element_text(hjust = 0.5))
  
  
  dsstd(X, mean = mean(X), sd = sd(X), nu = params[st,7], xi = params[st,6])
  psstd(X, mean = mean(X), sd = sd(X), nu = params[st,7], xi = params[st,6])
  qsstd(X, mean = mean(X), sd = sd(X), nu = params[st,7], xi = params[st,6])
  rsstd(X, mean = mean(X), sd = sd(X), nu = params[st,7], xi = params[st,6])
  
  ggplot(mapping = aes(sample = X)) +
    stat_qq_band(distribution = "sstd", dparams = list(mean = mean(X), sd = sd(X), nu = gfit@fit$coef[7], xi = gfit@fit$coef[6])) +
    stat_qq_point(distribution = "sstd", dparams = list(mean = mean(X), sd = sd(X), nu = gfit@fit$coef[7], xi = gfit@fit$coef[6]), size = 2) +  stat_qq_line(distribution = "sstd", dparams = list(mean = mean(X), sd = sd(X), nu = gfit@fit$coef[7], xi = gfit@fit$coef[6])) +
    xlab("Theoretical Quantiles") + ylab("Sample Quantiles") +
    ggtitle(my_title) +
    theme(plot.title = element_text(hjust = 0.5))
  
  
  dev.copy(jpeg,paste("C:/article3/Step5_ResidualsGARCH/P_", st, "QQplot.jpg"))
  dev.off()
  
}



# Make some fitted histograms
st<-85
X<-R[,st]
hist(X,100,prob=TRUE,main=paste('Histogram of Residuals in Point', st),fill = "lightgray")
box()
A <- seq(min(X), max(X), length = 100)
D<-dsstd(A, mean = mean(X), sd = sd(X), nu = params[st,7], xi = params[st,6])
lines(A,D,col="blue", lwd=2)
dev.copy(jpeg,paste("C:/article3/Step5_ResidualsGARCH/P_", st, "Histogram.jpg"))
dev.off()









Box.test(r, lag = 1, type = "Ljung-Box")





tsdisplay(r, lag.max = 200)

tsdisplay((r)**2, lag.max = 200)



