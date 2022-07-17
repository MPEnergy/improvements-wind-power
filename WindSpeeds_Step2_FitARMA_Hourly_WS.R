rm(list = ls())
library(R.matlab)
library(chron)
library(glmnet)
library(forecast)
library("astsa")
library(lmtest)



data <- readMat('Hourly_deseasonalized_test.mat')


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


# See the autocorrelation
for (st in 1:85){
v1 <- v[,st]
require(xts)
aux<-data.frame(dates=c(timeDate),speeds=v1)
ws<-xts(aux$speeds, order.by=as.POSIXct(aux$dates))
WS<-coredata(ws)
acf2(v1, main=paste("Point ",st))
dev.copy(jpeg,paste("C:/Users/u6032456/Documents/Thesis/PhDMeeting/article2/ARMA/NewSeasonality/P_", st, "_acf_pacf_deseasonalized.jpg"))
dev.off()
}



#Test some orders for ARIMA
st<-20
v1 <- v[,st]
aux<-data.frame(dates=c(timeDate),speeds=v1)
ws<-xts(aux$speeds, order.by=as.POSIXct(aux$dates))
WS<-coredata(ws)
obj <- auto.arima(ws)
p<-obj$arma[1] 
q<-obj$arma[2] #1
d<-obj$arma[6]  #2
fitarima <- arima(ws, order=c(p,d,q))
resarima <- residuals(fitarima)
acf2(resarima^2)


# Fit AR(2)
st<-1
v1 <- v[,st]
aux<-data.frame(dates=c(timeDate),speeds=v1)
ws<-xts(aux$speeds, order.by=as.POSIXct(aux$dates))
WS<-coredata(ws)
fit <- arima(v1, order=c(2,0,24),optim.control=list(maxit = 1000))
res<- residuals(fit)
acf2(res)
dev.copy(jpeg,paste("C:/Users/u6032456/Documents/Thesis/PhDMeeting/article2/Removing_Var/P_", st, "_acf_pacf_novar2_ARMA124.jpg"))
dev.off()

# Fit several models at a time and compare
aicvector<-matrix(nrow=85, ncol=2)
coefs<-matrix(nrow=85, ncol=28)
residuos<-matrix(nrow=38664, ncol=85)
for (st in 60:85){
  v1<-v[,st]
  tryCatch({
    fit24<-arima(v1, order=c(24,0,0),optim.control=list(maxit = 1000))
    res24<-residuals(fit24)
    acf2(res24,main=paste('Residuals in Point', st))
    dev.copy(jpeg,paste("C:/Users/u6032456/Documents/Thesis/PhDMeeting/article2/ARMA/NewSeasonality/AR(24)/P_", st, "_acf_pacf_var_AR24.jpg"))
    dev.off()
    acf2(res24^2,main=paste('Residuals in Point', st))
    dev.copy(jpeg,paste("C:/Users/u6032456/Documents/Thesis/PhDMeeting/article2/ARMA/NewSeasonality/AR(24)/P_", st, "_acf_pacf_var2_AR24.jpg"))
    dev.off()
    aicvector[st,1]<-fit24$aic
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  tryCatch({
    fit224<-arima(v1, order=c(2,0,24), optim.control=list(maxit = 1000))
    res224<-residuals(fit224)
    acf2(res224,main=paste('Residuals in Point', st))
    dev.copy(jpeg,paste("C:/Users/u6032456/Documents/Thesis/PhDMeeting/article2/ARMA/NewSeasonality/ARMA(2,24)/P_", st, "_acf_pacf_var_ARMA224.jpg"))
    dev.off()
    acf2(res224^2,main=paste('Residuals in Point', st))
    dev.copy(jpeg,paste("C:/Users/u6032456/Documents/Thesis/PhDMeeting/article2/ARMA/NewSeasonality/ARMA(2,24)/P_", st, "_acf_pacf_var2_ARMA224.jpg"))
    dev.off()
    aicvector[st,2]<-fit224$aic
    residuos[,st]=res224
    p = coeftest(fit224)
    for (i in 1:27){coefs[st,i]<-fit224$coef[i]}
    coefs[st,28]<-p[27,4]
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


coefs<-matrix(nrow = 85, ncol = 25)
residuos<-matrix(nrow=38664, ncol=85)
for (st in 1:85){
  v1 <- v[,st]
  fit <- arima(v1, order=c(1,0,24))
  res<- residuals(fit)
  residuos[,st]=res
  for (i in 1:25){coefs[st,i]<-fit$coef[i]}
  
 # acf2(res,  main=paste('Residuals in Point', st))
 # dev.copy(jpeg,paste("C:/Users/u6032456/Documents/Thesis/PhDMeeting/article2/ARMA/ARMA(1,24)/P_", st, "_acf_pacf_ARMA124.jpg"))
 # dev.off()
  
}

for (st in 1:85){
  r<-residuos[,st]
  acf2(r^2, main=paste('Squared Residuals in Point', st))
  dev.copy(jpeg,paste("C:/Users/u6032456/Documents/Thesis/PhDMeeting/article2/ARMA/ARMA(1,24)/P_", st, "_acf_pacf_res2.jpg"))
  dev.off()
}


trend = ma(log(ws), order = 24, centre = T)
plot(as.ts(log(ws)))
lines(trend)
plot(as.ts(trend))

detrend_ws = ws - trend
plot(as.ts(detrend_ws))









# Point 1
st <- 1
v1 <- v[,st]
require(xts)
aux<-data.frame(dates=c(timeDate),speeds=v1)
ws<-xts(aux$speeds, order.by=as.POSIXct(aux$dates))
WS<-coredata(ws)
acf2(WS, main=paste("Point ",st))

# Save
dev.copy(pdf,'C:/Users/u6032456/Desktop/Thesis/PhDMeeting/article1/GridData_TSAnalysis/P85_acf_pacf.pdf')
dev.off()

# Stationarity tests
library(tseries)
ADF<-adf.test(ws)
pADF<-ADF$p.value
KPSS<-kpss.test(ws, null="Trend")
pKPSS<-KPSS$p.value

# Fit ARIMA
obj <- auto.arima(ws)
p<-obj$arma[1] 
q<-obj$arma[2] #1
d<-obj$arma[6]  #2
fitarima <- arima(ws, order=c(p,d,q))
resarima <- residuals(fitarima)

# Fit AR(2)
fitar <- arima(WS, order=c(2,0,0))
coeftest(fitar)
res <- residuals(fitar)



 # Ploting the fit
ts.plot(WS)
AR_fit <- WS - res
points(AR_fit, type = "l", col = 2, lty = 2)
acf2(coredata(res), main=paste("Point ",st, "residuals"))
dev.copy(pdf,'C:/Users/u6032456/Desktop/Thesis/PhDMeeting/article1/GridData_TSAnalysis/P85_acf_pacf_res_ar7.pdf')
dev.off()

epsilon <- res
hist(epsilon, main=paste("Histogram of residuals in Point ",st))
dev.copy(pdf,'C:/Users/u6032456/Desktop/Thesis/PhDMeeting/article1/GridData_TSAnalysis/P85_hist_res.pdf')
dev.off()

## Check stationatiry
library(aTSA)
st = 47
v1 <- v[,st]
require(xts)
aux<-data.frame(dates=c(timeDate),speeds=v1)
ws<-xts(aux$speeds, order.by=as.POSIXct(aux$dates))

tseries::adf.test(ws)
tseries::kpss.test(ws)
forecast::ndiffs(ws,test="kpss")






isST<- vector()
ps<- matrix(nrow = 85, ncol = 2)
for (st in 1:85)
{
v1 <- v[,st]
require(xts)
aux<-data.frame(dates=c(timeDate),speeds=v1)
ws<-xts(aux$speeds, order.by=as.POSIXct(aux$dates))
ADF<-adf.test(ws)
pADF<-ADF$p.value
KPSS<-kpss.test(diff(ws,1), null = c("Level", "Trend"), lshort = TRUE)
pKPSS<-KPSS$p.value

if (pKPSS>=0.05)
{isST[st]=1
}
else
{isST[st]=0
}
 ps[st,1] = pADF
 ps[st,2]=pKPSS
}



## Store the residuals after AR(2)/AR(1)
library("xlsx")

st<-1
aux<-data.frame(dates=c(timeDate),speeds=v1)
ws<-xts(aux$speeds, order.by=as.POSIXct(aux$dates))
WS<-coredata(ws)
fitar <- arima(WS, order=c(2,0,0))
Epsilon <- residuals(fitar)

for (st in 2:85)
{
v1 <- v[,st]
#require(xts)
aux<-data.frame(dates=c(timeDate),speeds=v1)
ws<-xts(aux$speeds, order.by=as.POSIXct(aux$dates))
WS<-coredata(ws)

if (st ==14 || st==34 || st==35 || st==36)
{fitar <- arima(WS, order=c(1,0,0))
res <- residuals(fitar)}
else
{fitar <- arima(WS, order=c(2,0,0))
res <- residuals(fitar)}


Epsilon<-cbind(Epsilon,res)
}

write.xlsx(Epsilon, file = "C:/matlab/code/PMT_NP_Wind_TrabajoInvestigacion/Residuals_ar2.xlsx",
           sheetName = "Epsilon", append = TRUE)

# Test for normality
isNorm<-vector()
require(xts)
library("ggpubr")
for (st in 1:85)
{S=shapiro.test(Epsilon[,st])
if(S$p>0.05)
{isNorm[st]=1}
else
{isNorm[st]=0}
}
isNorm=which(isNorm==1)

ggqqplot(Epsilon[,30],main="Point 30")
dev.copy(pdf,'C:/Users/u6032456/Desktop/Thesis/PhDMeeting/article1/GridData_TSAnalysis/Residuals_p30_Q-Qplot.pdf')
dev.off()


# Study further the residuals
st<-1
epsilon<-Epsilon[,st]
aux<-data.frame(dates=c(timeDate),speeds=epsilon)
e<-xts(aux$speeds, order.by=as.POSIXct(aux$dates))
e<-coredata(e)
fitar <- arima(e, order=c(4,0,4))
res <- residuals(fitar)
ggqqplot(e,main="Point 1")

# Test normality for final residuals
rm(list = ls())
library(R.matlab)

residuals <- readMat('C:/matlab/code/PMT_NP_Wind_TrabajoInvestigacion/Residuals_final.mat')

d<-residuals$dates
newd<-strftime(as.POSIXct( ( d - 719529) * 86400,origin = "1970-01-01",tz = "UTC"), format = '%Y-%m-%d', 
               tz = 'UTC', usetz = FALSE)
library(lubridate)
timeDate <- as.POSIXct(newd,tz="UTC")
e <- residuals$e

st<-21
epsilon<-e[,st]
aux<-data.frame(dates=c(timeDate),speeds=epsilon)
r<-xts(aux$speeds, order.by=as.POSIXct(aux$dates))
r<-coredata(r)
ggqqplot(r,main="Point 1")


isNorm<-vector()
for (st in 1:85)
{epsilon<-e[,st]
aux<-data.frame(dates=c(timeDate),speeds=epsilon)
r<-xts(aux$speeds, order.by=as.POSIXct(aux$dates))
r<-coredata(r)
S=shapiro.test(r)
if(S$p>0.05)
{isNorm[st]=1}
else
{isNorm[st]=0}
}
isNorm=which(isNorm==1)









