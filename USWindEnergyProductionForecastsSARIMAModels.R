#LIBRERIAS Y PAQUETES ----------------------------------------------------------
rm(list=ls(all=TRUE))
library(forecast)
library(car)
library(carData)
library(TSA)
library(fANCOVA)
library("ggplot2")
library("showtext")
library(pdR)
library(timsac);library(lmtest) 
library(uroot)
windowsFonts(A = windowsFont("Times New Roman"))
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-regexponencial.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funciones-BP.LB.test-pruebaDW1.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funciones-Criterios.Informacion-Calidad.Intervalos.R")

#Lectura de datos--------------------------------------------------------------

Data=read.table(file.choose(),header=T,sep=",")
Data = Data$wind_power
yt=ts(Data,freq=12,start=c(2001,1))
lnyt=log(yt)

#Grafica de la serie --------------------------------------------------------------

win.graph() 
plot(yt,xlab = "Time",ylab='Yt')
win.graph() 
plot(lnyt,xlab = "Time",ylab='Log(Yt)')

#Identificacion usando solo los primeros 254 datos--------------------------------------------------------------

m=12
n=length(yt)-m
t=1:n
yt2=ts(yt[t],frequency=12,start=c(2001,1))
lnyt2=log(yt2)
dif1dif4lnyt2=diff(diff(lnyt2,12),1) #diferencia regular y estacional de orden 1 para logYt
tnuevo=(n+1):length(yt)
ytf=ts(yt[tnuevo],freq=12,start=c(2022,3))

#Graficos de las diferencias--------------------------------------------------------------
win.graph()
par(mfrow = c(2, 2))
plot(diff(lnyt2),ylab=expression(paste(nabla,sep="",log(Y[t]))))
abline(h=mean(diff(lnyt2)))
plot(diff(lnyt2,12),ylab=expression(paste(nabla[12],sep="",log(Y[t]))))
abline(h=mean(diff(lnyt2,12)))
plot(diff(diff(lnyt2,12),1),ylab=expression(paste(nabla,sep="",nabla[12],sep="",log(Y[t]))))
abline(h=mean(diff(diff(lnyt2,12),1)))
plot(diff(diff(lnyt2,1),12),ylab=expression(paste(nabla[12],sep="",nabla,sep="",log(Y[t]))))
abline(h=mean(diff(diff(lnyt2,1),12)))
par(mfrow = c(1, 1))
#ACF y PACF de la diferencia mixta de logYt, pero con identificacion de los k multiplos de s=12 --------------------------------------------------------------
win.graph()
par(mfrow = c(1, 2))
acf(as.numeric(dif1dif4lnyt2),ci.type="ma",lag.max=36,
    main=expression(paste("ACF ",sep=" ",nabla,sep="",nabla[12],sep="",log(Y[t]))),lwd=2)
abline(v=seq(12,36,by=12),col=2,lty=2)
pacf(as.numeric(dif1dif4lnyt2),lag.max=36,main=expression(paste("PACF ",sep=" ",nabla,sep="",nabla[12],sep="",log(Y[t]))),lwd=2)
abline(v=seq(12,36,by=12),col=2,lty=2)


#Test HEGY para la serie del logaritmo s ́olo con los primeros n=254 datos--------------------------------------------------------------
HEGY.test(wts=lnyt2,itsd=c(0,0,c(0)),selectlags=list(mode="aic", Pmax=12))$stats

#Identificacion de un ARMA estacionario sobre diferencia mixta de logYt armasubsets() --------------------------------------------------------------

win.graph(width=8,height=5)
plot(armasubsets(y=dif1dif4lnyt2,nar=12,nma=12,y.name='Y',ar.method='ml'))
win.graph(width=15,height=15)
plot(armasubsets(y=dif1dif4lnyt2,nar=24,nma=24,y.name='Y',ar.method='ols'))

#Auto.arima() --------------------------------------------------------------

auto.arima(lnyt2,ic="aic",seasonal.test="ocsb")
auto.arima(lnyt2,ic="aic",seasonal.test="ch")
auto.arima(lnyt2,ic="aic",seasonal.test="seas")
auto.arima(lnyt2,ic="bic",seasonal.test="ocsb")
auto.arima(lnyt2,ic="bic",seasonal.test="ch")
auto.arima(lnyt2,ic="bic",seasonal.test="seas")

#-----MODELO 1:ARIMA(2,1,1)(0,1,1)[12]------------------------------------------------------
modelo1=Arima(lnyt2,order=c(2,1,1),seasonal=list(order=c(0,1,1)),method="ML")
k1=length(coef(modelo1)[coef(modelo1)!=0]);k1 #n ́umero de par ́ametros del modelo1
coeftest(modelo1)
ythat1=exp(modelo1$fitted)*exp(modelo1$sigma2/2)
seudores1=yt2-ythat1

#-----MODELO 2 ARIMA(4,1,1)(2,0,0)[12] con deriva ------------------------------------------------------
modelo2=Arima(lnyt2,order=c(4,1,1),seasonal=list(order=c(2,0,0)),include.drift=T,method="ML")
k2=length(modelo2$coef[modelo2$coef!=0]);k2 
coeftest(modelo2)
ythat2=exp(modelo2$fitted)*exp(modelo2$sigma2/2)
seudores2=yt2-ythat2

#-----MODELO 3 armasubsets() ARIMA (1,1,7)(1,1,2)[12] ------------------------------------------------------
modelo3=Arima(lnyt2,order=c(1,1,7),seasonal=list(order=c(1,1,2)),fixed=c(NA,NA,NA,rep(0,4),NA,NA,0,NA),method="ML")
k3=length(modelo3$coef[modelo3$coef!=0]);k3
coeftest(modelo3)
ythat3=exp(modelo3$fitted)*exp(modelo3$sigma2/2)
seudores3=yt2-ythat3

#-----MODELO 4 armasubsets() ARIMA (5,1,5)(1,1,1)[12] ------------------------------------------------------
modelo4=Arima(lnyt2,order=c(1,1,2),seasonal=list(order=c(1,1,1)),method="ML")
k4=length(modelo4$coef[modelo4$coef!=0]);k4
coeftest(modelo4)
ythat4=exp(modelo4$fitted)*exp(modelo4$sigma2/2)
seudores4=yt2-ythat4

#------------------------------------------------------------------
#INDICADORES-PRUEBAS

#AIC-BIC ------------------------------------------------------------------
Criterios1=exp.crit.inf.resid(residuales=seudores1,n.par=k1); Criterios1
Criterios2=exp.crit.inf.resid(residuales=seudores2,n.par=k2); Criterios2 
Criterios3=exp.crit.inf.resid(residuales=seudores3,n.par=k3); Criterios3 
Criterios4=exp.crit.inf.resid(residuales=seudores4,n.par=k4); Criterios4 

#Ljun-Box------------------------------------------------------------------
BP.LB.test(residuals(modelo1),maxlag=36,type="Ljung")
BP.LB.test(residuals(modelo2),maxlag=36,type="Ljung")
BP.LB.test(residuals(modelo3),maxlag=36,type="Ljung")
BP.LB.test(residuals(modelo4),maxlag=36,type="Ljung")

#Normalidad------------------------------------------------------------------
shapiro.test(residuals(modelo1))
shapiro.test(residuals(modelo2))
shapiro.test(residuals(modelo3))
shapiro.test(residuals(modelo4))

#Pronostico------------------------------------------------------------------
predmod1=ts(exp(as.data.frame(forecast(modelo1,h=m,level=95)))*exp(modelo1$sigma2/2),freq=12,start=c(2022,3))
predmod1 #pron ́osticos puntuales e I.P
ytpron1=predmod1[,1] #pron ́osticos puntuales

predmod2=ts(exp(as.data.frame(forecast(modelo2,h=m,level=95)))*exp(modelo2$sigma2/2),freq=12,start=c(2022,3))
predmod2 #pron ́osticos puntuales e I.P
ytpron2=predmod2[,1] #pron ́osticos puntuales

predmod3=ts(exp(as.data.frame(forecast(modelo3,h=m,level=95)))*exp(modelo3$sigma2/2),freq=12,start=c(2022,3))
predmod3 #pron ́osticos puntuales e I.P
ytpron3=predmod3[,1] #pron ́osticos puntuales

predmod4=ts(exp(as.data.frame(forecast(modelo4,h=m,level=95)))*exp(modelo4$sigma2/2),freq=12,start=c(2022,3))
predmod4 #pron ́osticos puntuales e I.P
ytpron4=predmod4[,1] #pron ́osticos puntuales


#Medidas de pronostico ------------------------------------------------------------------
accuracy(ytpron1,ytf)
AmpCobmodelo1=amplitud.cobertura(real=ytf,LIP=predmod1[,2],LSP=predmod1[,3]) 
accuracy(ytpron2,ytf)
AmpCobmodelo2=amplitud.cobertura(real=ytf,LIP=predmod2[,2],LSP=predmod2[,3]) 
accuracy(ytpron3,ytf)
AmpCobmodelo3=amplitud.cobertura(real=ytf,LIP=predmod3[,2],LSP=predmod3[,3])
accuracy(ytpron4,ytf)
AmpCobmodelo4=amplitud.cobertura(real=ytf,LIP=predmod4[,2],LSP=predmod4[,3]) 

#------------------------------------------------------------------
#GRAFICAS
#Ajustes
win.graph()
par(mfrow = c(2, 2))
plot(yt, ylab = "Wind Energy", main = "Model 1")
lines(ythat1,col=2)
legend("topleft",legend=c("Original","Fit model 1"),col=c(1,2,4),lty=1,lwd=2)
plot(yt, ylab = "Wind Energy", main = "Model 2")
lines(ythat2,col=2)
legend("topleft",legend=c("Original","Fit model 2"),col=c(1,2,4),lty=1,lwd=2)
plot(yt, ylab = "Wind Energy", main = "Model 3")
lines(ythat3,col=2)
legend("topleft",legend=c("Original","Fit model 3"),col=c(1,2,4),lty=1,lwd=2)
plot(yt, ylab = "Wind Energy", main = "Model 4")
lines(ythat4,col=2)
legend("topleft",legend=c("Original","Fit model 4"),col=c(1,2,4),lty=1,lwd=2)

#Analisis de residuales
#Graficas de residuales vs tiempo de los modelos 1 al 4------------------------------------------------------------------
win.graph()
par(mfrow = c(2, 2))
plot(residuals(modelo1),ylab="Residuales",main="Model 1")
abline(h=c(-2*sqrt(modelo1$sigma2),0,2*sqrt(modelo1$sigma2)),lty=2,col=2)
plot(residuals(modelo2),ylab="Residuales",main="Model 2")
abline(h=c(-2*sqrt(modelo2$sigma2),0,2*sqrt(modelo2$sigma2)),lty=2,col=2)
plot(residuals(modelo3),ylab="Residuales",main="Model 3")
abline(h=c(-2*sqrt(modelo3$sigma2),0,2*sqrt(modelo3$sigma2)),lty=2,col=2)
plot(residuals(modelo4),ylab="Residuales",main="Model 4")
abline(h=c(-2*sqrt(modelo4$sigma2),0,2*sqrt(modelo4$sigma2)),lty=2,col=2)

#Graficas de residuales vs ajustados de los modelos 1 al 4------------------------------------------------------------------
win.graph()
par(mfrow = c(2, 2))
plot(modelo1$fitted,residuals(modelo1),ylab="Residuales",main="Model 1")
abline(h=c(-2*sqrt(modelo1$sigma2),0,2*sqrt(modelo1$sigma2)),lty=2,col=2)
plot(modelo2$fitted,residuals(modelo2),ylab="Residuales",main="Model 2")
abline(h=c(-2*sqrt(modelo2$sigma2),0,2*sqrt(modelo2$sigma2)),lty=2,col=2)
plot(modelo3$fitted,residuals(modelo3),ylab="Residuales",main="Model 3")
abline(h=c(-2*sqrt(modelo3$sigma2),0,2*sqrt(modelo3$sigma2)),lty=2,col=2)
plot(modelo4$fitted,residuals(modelo4),ylab="Residuales",main="Model 4")
abline(h=c(-2*sqrt(modelo4$sigma2),0,2*sqrt(modelo4$sigma2)),lty=2,col=2)

#Gráficas de ACF y PACf residuales de los modelos 1 a #Graficas de residuales vs tiempo de los modelos 1 al 4------------------------------------------------------------------4
win.graph()
par(mfrow = c(2, 2))
acf(as.numeric(residuals(modelo1)),ci.type="ma",lag.max=36,main="ACF residuals Model 1")
acf(as.numeric(residuals(modelo2)),ci.type="ma",lag.max=36,main="ACF residuals Model 2")
acf(as.numeric(residuals(modelo3)),ci.type="ma",lag.max=36,main="ACF residuals Model 3")
acf(as.numeric(residuals(modelo4)),ci.type="ma",lag.max=36,main="ACF residuals Model 4")

win.graph()
par(mfrow = c(2, 2))
pacf(as.numeric(residuals(modelo1)),lag.max=36,main="PACF residuals Model 1")
pacf(as.numeric(residuals(modelo2)),lag.max=36,main="PACF residuals Model 2")
pacf(as.numeric(residuals(modelo3)),lag.max=36,main="PACF residuals Model 3")
pacf(as.numeric(residuals(modelo4)),lag.max=36,main="PACF residuals Model 4")

#Tabla resumen criterios de información ------------------------------------------------------------------
tablacriterios=rbind(Criterios1,Criterios2,Criterios3,Criterios4)
rownames(tablacriterios)=paste0("Model ",1:4)
tablacriterios

#Tabla resumen de medidas de precision pronosticos------------------------------------------------------------------
comp=rbind(accuracy(ytpron1,ytf),accuracy(ytpron2,ytf),accuracy(ytpron3,ytf),accuracy(ytpron4,ytf))[,c(2,3,5)]
otros=rbind(AmpCobmodelo1,AmpCobmodelo2,AmpCobmodelo3,AmpCobmodelo4)
tablaprecis=cbind(comp,otros)
rownames(tablaprecis)=paste0("Model ",1:4)
tablaprecis

#Comparación Gráfica de los pronósticos de los 4 modelos------------------------------------------------------------------
win.graph()
plot(ytf,type="b",ylab="Wind Energy",pch=19,col=1,lwd=2,xaxt="n")
#plot(ytnuevo,type="b",ylab="Wind Energy",col=1,pch=19,ylim=c(min(ytf,ytpron1,ytpron2,ytpron3,ytpron4),max(ytnuevo,ytpron1,ytpron2,ytpron3,ytpron4)+25),lwd=2,xaxt="n", main="Comparative graph of the forecasts of the 4 models")
axis(1,at=time(ytf),labels=c("22-III","22-IV","22-V","22-VI","22-VII","22-VIII","22-IX","22-X","22-XI","22-XII","23-I","23-II"),cex.axis=0.7)
lines(ytpron1,col=2,pch=1,type="b",lwd=2)
lines(ytpron2,col=3,pch=2,type="b",lwd=2)
lines(ytpron3,col=4,pch=3,type="b",lwd=2)
lines(ytpron4,col=5,pch=4,type="b",lwd=2)
legend("bottomleft",legend=c("Real","Model 1","Model 2","Model 3","Model 4"),pch=c(19,1:4),col=c(1:5),lwd=2)
#tabla resumen de pronosticos model 4
tablar.pronosticos=cbind(ytf,ytpron3)
tablar.pronosticos

#--------------------------------------------------------------
#-------------------------------------------------------------
#Prediction model for 2 next years 
n_S=length(yt) #Using all data available

#Fitting
t_S=1:n_S
yt_S=ts(yt[t_S],freq=m,start=c(2001,1))
mes_S=seasonaldummy(yt_S) 

#PARA LOS PRONOSTICOS
tnuevo_S=(n_S+1):(length(yt)+24)
ytnuevo_S=ts(yt[tnuevo_S],freq=m,start=c(2023,2)) 
mesnuevo=seasonaldummy(yt_S,h=m)

#Initialize model
#h=24 next 24 months
mod_S=Arima(log(yt_S),order=c(1,1,7),seasonal=list(order=c(1,1,2)),fixed=c(NA,NA,NA,rep(0,4),NA,NA,0,NA),method="ML")
kS=length(mod_S$coef[mod_S$coef!=0]);kS
coeftest(mod_S)
ythatS=exp(mod_S$fitted)*exp(mod_S$sigma2/2)
seudores3=yt2-ythatS

#Pronosticos
predmodS=ts(exp(as.data.frame(forecast(mod_S,h=24,level=95)))*exp(mod_S$sigma2/2),freq=12,start=c(2023,3))
predmodS #pron ́osticos puntuales e I.P
ytpronS=predmodS[,1] #pron ́osticos puntuales


#Graph predictions
months_labels <- c("23-Mar", "23-Apr", "23-May", "23-Jun", "23-Jul", "23-Aug", "23-Sep", "23-Oct", "23-Nov", "23-Dec", "24-Jan", "24-Feb", "24-Mar", "24-Apr", "24-May", "24-Jun", "24-Jul", "24-Aug", "24-Sep", "24-Oct", "24-Nov", "24-Dec", "25-Jan", "25-Feb")
win.graph()
plot(ytpronS,
     xaxt = "n", 
     xlab = "Time",
     ylab = "Wind Energy (MW)",
     main = "Forecasts Wind Energy Production in US 2023 - 2025",
     type = "l", col = "blue", lwd = 2)
axis(1, at = time(ytpronS), labels = months_labels,las = 2)
grid()