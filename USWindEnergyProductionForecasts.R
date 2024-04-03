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
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/
Funciones-Criterios.Informacion-Calidad.Intervalos.R")

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
plot(diff(lnyt2),ylab=expression(paste(nabla,sep="",log(Y[t]))))
abline(h=mean(diff(lnyt2)))
win.graph()
plot(diff(lnyt2,12),ylab=expression(paste(nabla[12],sep="",log(Y[t]))))
abline(h=mean(diff(lnyt2,12)))
win.graph()
plot(diff(diff(lnyt2,12),1),ylab=expression(paste(nabla,sep="",nabla[12],sep="",log(Y[t]))))
abline(h=mean(diff(diff(lnyt2,12),1)))
win.graph()
plot(diff(diff(lnyt2,1),12),ylab=expression(paste(nabla[12],sep="",nabla,sep="",log(Y[t]))))
abline(h=mean(diff(diff(lnyt2,1),12)))

#ACF y PACF de la diferencia mixta de logYt, pero con identificacion de los k multiplos de s=12 --------------------------------------------------------------

win.graph(width=10,height=5)
acf(as.numeric(dif1dif4lnyt2),ci.type="ma",lag.max=36,
    main=expression(paste("ACF de",sep=" ",nabla,sep="",nabla[12],sep="",log(Y[t]))),lwd=3)
abline(v=seq(12,36,by=12),col=2,lty=2)

win.graph(width=10,height=5)
pacf(as.numeric(dif1dif4lnyt2),lag.max=36,main=expression(paste("PACF de",sep=" ",nabla,sep="",nabla[12],sep="",log(Y[t]))),lwd=3)
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

#-----MODELO 1:ARIMA(1,1,1)(0,1,1)[12]------------------------------------------------------
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
#Modelo 1
win.graph()
plot(yt,ylab=expression(Y[t]))
lines(ythat1,col=2)
legend("topleft",legend=c("Original","Ajuste modelo 1"),col=c(1,2,4),lty=1,lwd=2)
#Modelo 2
win.graph()
plot(yt,ylab=expression(Y[t]))
lines(ythat2,col=2)
legend("topleft",legend=c("Original","Ajuste modelo 2"),col=c(1,2,4),lty=1,lwd=2)
#Modelo 3
win.graph()
plot(yt,ylab=expression(Y[t]))
lines(ythat3,col=2)
legend("topleft",legend=c("Original","Model 3 adjustment"),col=c(1,2,4),lty=1,lwd=2)
#Modelo 4
win.graph()
plot(yt,ylab=expression(Y[t]))
lines(ythat4,col=2)
legend("topleft",legend=c("Original","Ajuste modelo 4"),col=c(1,2,4),lty=1,lwd=2)

#Analisis de residuales
#Graficas de residuales vs tiempo de los modelos 1 al 4------------------------------------------------------------------
win.graph()
plot(residuals(modelo1))
abline(h=c(-2*sqrt(modelo1$sigma2),0,2*sqrt(modelo1$sigma2)),lty=2,col=2)
win.graph()
plot(residuals(modelo2))
abline(h=c(-2*sqrt(modelo2$sigma2),0,2*sqrt(modelo2$sigma2)),lty=2,col=2)
win.graph()
plot(residuals(modelo3))
abline(h=c(-2*sqrt(modelo3$sigma2),0,2*sqrt(modelo3$sigma2)),lty=2,col=2)
win.graph()
plot(residuals(modelo4))
abline(h=c(-2*sqrt(modelo4$sigma2),0,2*sqrt(modelo4$sigma2)),lty=2,col=2)

#Graficas de residuales vs ajustados de los modelos 1 al 4------------------------------------------------------------------
win.graph()
plot(modelo1$fitted,residuals(modelo1))
abline(h=c(-2*sqrt(modelo1$sigma2),0,2*sqrt(modelo1$sigma2)),lty=2,col=2)
win.graph()
plot(modelo2$fitted,residuals(modelo2))
abline(h=c(-2*sqrt(modelo2$sigma2),0,2*sqrt(modelo2$sigma2)),lty=2,col=2)
win.graph()
plot(modelo3$fitted,residuals(modelo3))
abline(h=c(-2*sqrt(modelo3$sigma2),0,2*sqrt(modelo3$sigma2)),lty=2,col=2)
win.graph()
plot(modelo4$fitted,residuals(modelo4))
abline(h=c(-2*sqrt(modelo4$sigma2),0,2*sqrt(modelo4$sigma2)),lty=2,col=2)

#Gráficas de ACF y PACf residuales de los modelos 1 a #Graficas de residuales vs tiempo de los modelos 1 al 4------------------------------------------------------------------4
win.graph()
acf(as.numeric(residuals(modelo1)),ci.type="ma",lag.max=36,main="ACF residuos Modelo 1")
win.graph()
acf(as.numeric(residuals(modelo2)),ci.type="ma",lag.max=36,main="ACF residuos Modelo 2")
win.graph()
acf(as.numeric(residuals(modelo3)),ci.type="ma",lag.max=36,main="ACF residuos Modelo 3")
win.graph()
acf(as.numeric(residuals(modelo4)),ci.type="ma",lag.max=36,main="ACF residuos Modelo 4")

win.graph()
pacf(as.numeric(residuals(modelo1)),lag.max=36,main="PACF residuos Modelo 1")
win.graph()
pacf(as.numeric(residuals(modelo2)),lag.max=36,main="PACF residuos Modelo 2")
win.graph()
pacf(as.numeric(residuals(modelo3)),lag.max=36,main="PACF residuos Modelo 3")
win.graph()
pacf(as.numeric(residuals(modelo4)),lag.max=36,main="PACF residuos Modelo 4")

#Tabla resumen criterios de información ------------------------------------------------------------------
tablacriterios=rbind(Criterios1,Criterios2,Criterios3,Criterios4)
rownames(tablacriterios)=paste0("Modelo",1:4)
tablacriterios

#Tabla resumen de medidas de precision pronosticos------------------------------------------------------------------
comp=rbind(accuracy(ytpron1,ytf),accuracy(ytpron2,ytf),accuracy(ytpron3,ytf),accuracy(ytpron4,ytf))[,c(2,3,5)]
otros=rbind(AmpCobmodelo1,AmpCobmodelo2,AmpCobmodelo3,AmpCobmodelo4)
tablaprecis=cbind(comp,otros)
rownames(tablaprecis)=paste0("Modelo",1:4)
tablaprecis

#Comparación Gráfica de los pronósticos de los 4 modelos------------------------------------------------------------------
win.graph()
plot(ytf,type="b",pch=19,col=1,lwd=2,xaxt="n")
axis(1,at=time(ytf),labels=c("Mar-22","Apr-22","May-22","Jun-22","Jul-22","Aug-22","Sep-22","Oct-22","Nov-22","Dec-22","Jan-23","Feb-23"))
lines(ytpron1,pch=1,col="brown",type="b",lwd=2)
lines(ytpron2,pch=2,type="b",col=2,lwd=2)
lines(ytpron3,pch=3,type="b",col=3,lwd=2)
lines(ytpron4,pch=4,type="b",col=4,lwd=2)
legend("bottomleft",legend=c("Real",paste0("Model",1:4)),pch=c(19,1:4),col=c(1,"brown",2:4),lwd=2)