#LIBRERIAS Y PAQUETES ----------------------------------------------------------
rm(list=ls(all=TRUE))
library(TSA)
library(forecast)
library(fANCOVA)
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-regexponencial.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funciones-Criterios.Informacion-Calidad.Intervalos.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-Descomp.Loess.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-SuavizamientoEstacional.R")

#Lectura de datos--------------------------------------------------------------

Data=read.table(file.choose(),header=T,sep=",")
Data = Data$wind_power
yt=ts(Data,freq=12,start=c(2001,1))
lnyt=log(yt)

#Grafica de la serie --------------------------------------------------------------

win.graph()
par(mfrow = c(1, 2))
plot(yt, xlab = "Time", ylab = "Wind Energy (MW)", type = "l", col = "blue", lwd = 2)
title(main = "Wind Energy Production 2001-2023", cex.main = 1.2)
grid()
plot(lnyt, xlab = "Time", ylab = "Log(Wind Energy)", type = "l", col = "blue", lwd = 2)
title(main = "Wind Energy Production 2001-2023 in Log Scale", cex.main = 1.2)
grid()
par(cex.lab = 1.2, cex.axis = 1.2)

#--ANALISIS DESCRIPTIVO --------------------------------------------------------

#Realizando descompose
win.graph()
descom.multiplicativo= decompose(lnyt,type = "multiplicative")
plot(descom.multiplicativo)

#Graficando la tendencia 
win.graph()
Tt.log=decompose(lnyt)$trend
plot(Tt.log,ylim=c(min(lnyt),max(lnyt)),main="Trend of Time Serie")

win.graph()
boxplot(lnyt ~ cycle(yt), 
        main = "Monthly Box-Plots of Wind Energy", 
        xlab = "Month", 
        ylab = "Log Wind Energy (MW)", 
        col = "lightblue",  
        border = "blue",    
        ylim = c(min(lnyt), max(lnyt)),  
        cex.axis = 1.2,     
        cex.lab = 1.2,
        names = month.abb)      

#DEFINIENDO VARIABLES Y CREACION DE LA MATRIZ NO MOVER -------------------------
m=12
n=length(yt)-m

#PARA EL AJUSTE

t=1:n
t2=t^2
yt_a=ts(yt[t],freq=m,start=c(2001,1))
mes=seasonaldummy(yt_a) 

#Matriz de diseño en el ajuste
X1=data.frame(t,t2,mes)
head(X1) 

#PARA LOS PRONOSTICOS

tnuevo=(n+1):length(yt) 
t2nuevo=tnuevo^2
ytnuevo=ts(yt[tnuevo],freq=m,start=c(2022,3)) 
mesnuevo=seasonaldummy(yt_a,h=m)

#Matriz de diseño en los pronosticos
X1nuevo=data.frame(t=tnuevo,t2=t2nuevo,mesnuevo)
head(X1nuevo) 


#-----MODELO 1:COMPLETAMENTE MULTIPLICATIVA ------------------------------------------------------

mod1=lm(log(yt_a)~.,data=X1)
summary(mod1)

#Calculo valores ajustados del modelo 1
ythatmod1=ts(exp(fitted(mod1))*exp(summary(mod1)$sigma^2/2),freq=m,start=start(yt_a))

#Calculo de los criterios AIC y BIC en modelo 1
nparmod1=length(coef(mod1)[coef(mod1)!=0]);nparmod1 
res.orig.mod1=yt_a-ythatmod1
Criterios1= exp.crit.inf.resid(residuales=res.orig.mod1,n.par=nparmod1);Criterios1

#Pronosticos del modelo 1 en la escala original
pronos1=exp(predict(mod1,newdata=X1nuevo,interval="prediction",level=0.95))*exp(summary(mod1)$sigma^2/2)
pronos1=ts(pronos1,freq=m,start=start(ytnuevo))
pronos1
ytpron1=pronos1[,1]

#precision pronosticos puntuales modelo 1
accuracy(ytpron1,ytnuevo)

#precision pronosticos por I.P modelo 1
amplcobmod1=amplitud.cobertura(real=ytnuevo,LIP=pronos1[,2],LSP=pronos1[,3]);amplcobmod1

#-----MODELO 2:PARCIALMENTE MULTIPLIVATIVA ------------------------------------------------------

parametros.mod2=c(paste0("beta",0:2),paste0("delta",1:11)) 
parametros.mod2

mod2=regexponencial(respuesta=yt_a,data=X1,names.param=parametros.mod2)
summary(mod2)

#Calculo valores ajustados del modelo 2
ythatmod2=ts(fitted(mod2),freq=m,start=start(yt))

#Calculo de los criterios AIC y BIC en modelo 2
nparmod2=length(coef(mod2)[coef(mod2)!=0]);nparmod2          
Criterios2=exp.crit.inf.resid(residuales=residuals(mod2),n.par=nparmod2); Criterios2

#Pronosticos del modelo 2 en la escala original, solo son de tipo puntual por ser modelo no lineal
pronos2=predict(mod2,newdata=X1nuevo,interval="prediction",level=0.95)
ytpron2=ts(pronos2,freq=m,start=start(ytnuevo))
ytpron2 

#precision pronosticos puntuales modelo 2
accuracy(ytpron2,ytnuevo) 

#-----MODELO 3:SEHW ------------------------------------------------------------

mod3=SuavizamientoEstacional(yt_a,seasonal="multiplicative",h=m)
str(mod3) 

#Calculo de AIC y BIC 
s=m
p3=(s-1)+2
Criterios3=exp.crit.inf.resid(residuales=residuals(mod3),n.par=p3);Criterios3
MSE3=mod3$MSE 

#Pronosticos del modelo 3 
pronos3=mod3$forecast
pronos3
ytpron3=pronos3[,1] #solo los pronosticos puntuales del suavizamiento

#precision pronosticos puntuales modelo 3
accuracy(ytpron3,ytnuevo) 

#Precision pronosticos por I.P Modelo 3
amplcobmod3=amplitud.cobertura(real=ytnuevo,LIP=pronos3[,2],LSP=pronos3[,3]);amplcobmod3

#-----MODELO 4:DLL(AICC) -------------------------------------------------------

mod4=Descomp.Loess(serie.ajuste=yt_a,h=m,tipo.descomp="multiplicative",grado=1,criterio="aicc")
str(mod4) 

#Calculo AIC y BIC 
Criterios4=exp.crit.inf.resid(residuales=residuals(mod4),n.par=mod4$p);Criterios4
MSE4=mod4$MSE 

#Pronosticos del modelo 4 
pronos4=mod4$tablapron
pronos4

#Precision pronosticos puntuales
ytpron4=mod4$ytpron
accuracy(ytpron4,ytnuevo)

#Graficando St estimada por el filtro de descomposicion
win.graph()
plot(mod4$St,ylab=expression(hat(S)[t]), main="Gráfica de la estimación de la estacionalidad") 

#--TABLAS DE LAS MEDIDAS DE AJUSTE Y PRONOSTICOS DE LOS 4 MODELOS --------------

#Tabulando medidas de ajuste
tabla1.criterios=rbind(Criterios1,Criterios2,Criterios3,Criterios4)
rownames(tabla1.criterios)=c("Modelo 1","Modelo 2","Modelo 3","Modelo 4")
tabla1.criterios

#tabla resumen de pronosticos 
tabla.pronosticos=cbind(pronos1,ytpron2,pronos3,ytpron4)
tabla.pronosticos

#Tabulando medidas de pronosticos
precision.puntuales=rbind(accuracy(ytpron1,ytnuevo), accuracy(ytpron2,ytnuevo), accuracy(ytpron3,ytnuevo), accuracy(ytpron4,ytnuevo))[,c(2,3,5)]
precision.intervalos=rbind(amplcobmod1,c(NA,NA),amplcobmod3,c(NA,NA))
tabla2.precision=cbind(precision.puntuales,precision.intervalos)
rownames(tabla2.precision)=c("Modelo 1","Modelo 2","Modelo 3","Modelo 4")
tabla2.precision

#--GRAFICAS VALORES AJUSTADOS PARA LOS 4 MODELOS -------------------------------

win.graph()
par(mfrow = c(2, 2))
plot(yt, ylab = "Wind Energy", main = "Model 1")
lines(ythatmod1, col = 2, lwd = 2)
legend("topleft", legend = c("Original", "Model 1"), lty = 1, col = c(1, 2))
plot(yt, ylab = "Wind Energy", main = "Model 2")
lines(ythatmod2, col = 2, lwd = 2)
legend("topleft", legend = c("Original", "Model 2"), lty = 1, col = c(1, 2))
plot(yt, ylab = "Wind Energy", main = "Model 3")
lines(fitted(mod3), col = 2, lwd = 2)
legend("topleft", legend = c("Original", "Model 3"), lty = 1, col = c(1, 2))
plot(yt, ylab = "Wind Energy", main = "Model 4")
lines(fitted(mod4), col = 2, lwd = 2)
legend("topleft", legend = c("Original", "Model 4"), lty = 1, col = c(1, 2))
par(mfrow = c(1, 1))

#--GRAFICAS RESIDUOS DE AJUSTE PARA LOS 4 MODELOS ------------------------------
#Residuales vs Time
win.graph()
par(mfrow = c(2, 2))
plot.ts(residuals(mod1),ylab="Residuales",main="Model 1")
abline(h=c(-2*summary(mod1)$sigma,0,2*summary(mod1)$sigma),lty=1,col=2)
plot.ts(residuals(mod2),ylab="Residuales",main="Model 2")
abline(h=c(-2*summary(mod2)$sigma,0,2*summary(mod2)$sigma),lty=1,col=2)
plot(residuals(mod3),ylim=c(min(-2*sqrt(MSE3),residuals(mod3)),max(2*sqrt(MSE3),residuals(mod3))),ylab="Residuales",main="Model 3")
abline(h=c(-2*sqrt(MSE3),0,2*sqrt(MSE3)),lty=1,col=2)
plot(residuals(mod4),ylim=c(min(-2*sqrt(MSE4),residuals(mod4)),max(2*sqrt(MSE4),residuals(mod4))),ylab="Residuales",main="Model 4")
abline(h=c(-2*sqrt(MSE4),0,2*sqrt(MSE4)),lty=1,col=2)
par(mfrow = c(1, 1))

#Residuales vs fitted values
win.graph()
par(mfrow = c(2, 2))
plot(fitted(mod1),residuals(mod1),ylab="Residuals vs fitted values",main="Model 1")
abline(h=c(-2*summary(mod1)$sigma,0,2*summary(mod1)$sigma),lty=1,col=2)
plot(fitted(mod2),residuals(mod2),ylab="Residuals vs fitted values",main="Model 2")
abline(h=c(-2*summary(mod2)$sigma,0,2*summary(mod2)$sigma),lty=1,col=2)
plot(as.numeric(fitted(mod3)),residuals(mod3),ylim=c(min(-2*sqrt(MSE3),residuals(mod3)),max(2*sqrt(MSE3),residuals(mod3))),ylab="Residuals vs fitted values",main="Model 3")
abline(h=c(-2*sqrt(MSE3),0,2*sqrt(MSE3)),lty=1,col=2)
plot(as.numeric(fitted(mod4)),residuals(mod4),ylim=c(min(-2*sqrt(MSE4),residuals(mod4)),max(2*sqrt(MSE4),residuals(mod4))),ylab="Residuals vs fitted values",main="Model 4")
abline(h=c(-2*sqrt(MSE4),0,2*sqrt(MSE4)),lty=1,col=2)

#--GRAFICA COMPARATIVA DE LOS PRONOSTICOS DE LOS 4 MODELOS ---------------------
win.graph()
plot(ytnuevo,type="b",ylab="Wind Energy",col=1,pch=19,ylim=c(min(ytnuevo,ytpron1,ytpron2,ytpron3,ytpron4),max(ytnuevo,ytpron1,ytpron2,ytpron3,ytpron4)+25),lwd=2,xaxt="n", main="Comparative graph of the forecasts of the 4 models")
axis(1,at=time(ytnuevo),labels=c("22-III","22-IV","22-V","22-VI","22-VII","22-VIII","22-IX","22-X","22-XI","22-XII","23-I","23-II"),cex.axis=0.7)
lines(ytpron1,col=2,pch=1,type="b",lwd=2)
lines(ytpron2,col=3,pch=2,type="b",lwd=2)
lines(ytpron3,col=4,pch=3,type="b",lwd=2)
lines(ytpron4,col=5,pch=4,type="b",lwd=2)
legend("bottomleft",legend=c("Real","Model 1","Model 2","Model 3","Model 4"),pch=c(19,1:4),col=c(1:5),lwd=2)

#tabla resumen de pronosticos model 4
tablar.pronosticos=cbind(ytnuevo,ytpron4)
tablar.pronosticos


#--------------------------------------------------------------
#-------------------------------------------------------------
#Prediction model for 2 next years 
n_DLL=length(yt) #Using all data available

#Fitting
t_DLL=1:n_DLL
yt_DLL=ts(yt[t_DLL],freq=m,start=c(2001,1))
mes_DLL=seasonaldummy(yt_DLL) 

#PARA LOS PRONOSTICOS
tnuevo_DLL=(n_DLL+1):(length(yt)+24)
ytnuevo=ts(yt[tnuevo_DLL],freq=m,start=c(2023,2)) 
mesnuevo=seasonaldummy(yt_DLL,h=m)

#Initialize model
#h=24 next 24 months
mod_DLL=Descomp.Loess(serie.ajuste=yt_DLL,h=24,tipo.descomp="multiplicative",grado=1,criterio="aicc")
str(mod_DLL) 

#Pronosticos del modelo DLL
predictions_next_24_months=mod_DLL$tablapron
predictions <- predictions_next_24_months[, "Pron_serie"]
predictions
predictions_ts <- ts(predictions, start = c(2023, 3), frequency = 12)


#Graph predictions
months_labels <- c("23-Mar", "23-Apr", "23-May", "23-Jun", "23-Jul", "23-Aug", "23-Sep", "23-Oct", "23-Nov", "23-Dec", "24-Jan", "24-Feb", "24-Mar", "24-Apr", "24-May", "24-Jun", "24-Jul", "24-Aug", "24-Sep", "24-Oct", "24-Nov", "24-Dec", "25-Jan", "25-Feb")
win.graph()
plot(predictions_ts,
     xaxt = "n", 
     xlab = "Time",
     ylab = "Wind Energy (MW)",
     main = "Forecasts Wind Energy Production in US 2023 - 2025",
     type = "l", col = "blue", lwd = 2)
axis(1, at = time(predictions_ts), labels = months_labels,las = 2)
grid()