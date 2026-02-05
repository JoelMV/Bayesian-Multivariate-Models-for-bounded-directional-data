library(armspp)
library(TeachingDemos)
library(tictoc)
library(mvtnorm)
library(coda)
library(sns)

GammaP<-function(Th,A1,A2,B1){
  (B1^A1)*(cos(Th)^(A1-1))*(sin(Th)^(A2-1))/(beta(A1,A2)*(B1*cos(Th)+sin(Th))^(A1+A2))     
}
GammaP2<-function(Th,A1,A2,B2){
  (B2^A2)*(cos(Th)^(A1-1))*(sin(Th)^(A2-1))/(beta(A1,A2)*(cos(Th)+B2*sin(Th))^(A1+A2))     
}
qgammap<-function(U,A1,A2,B1){
  X1<-rgamma(100000,A1,B1)
  X2<-rgamma(100000,A2,1)
  T1<-atan2(X2,X1)
  q<-sapply(U, function(x) quantile(T1,x,names=FALSE ) )
  return(q)
}
qgammap2<-function(U,A1,A2,B2){
  X1<-rgamma(100000,A1,1)
  X2<-rgamma(100000,A2,B2)
  T1<-atan2(X2,X1)
  q<-sapply(U, function(x) quantile(T1,x,names=FALSE ) )
  return(q)
}
pgammap2<-function(U,A1,A2,B2){
  p<-sapply(U, function(u) integrate(function(s) GammaP2(s,A1,A2,B2),lower = 0,upper = u)$value)
  p<-replace(p,which(p>0.9998),0.9998)
  p<-replace(p,which(p<0.0002),0.0002)
  return(p)
}
logdctgammay.x<-function(x,y,A1,A2,B1,A3,A4,B2,V,R){
  dmvt(c(qt(pgammap2(x,A1,A2,B1),V),qt(pgammap2(y,A3,A4,B2),V)),sigma = matrix(c(1,R,R,1),ncol = 2),df=V)+log(GammaP2(y,A3,A4,B2))-sum(log(dt(qt(pgammap2(x,A1,A2,B1),V),V)),log(dt(qt(pgammap2(y,A3,A4,B2),V),V)))
}
logdctgammax.y<-function(x,y,A1,A2,B1,A3,A4,B2,V,R){
  dmvt(c(qt(pgammap2(x,A1,A2,B1),V),qt(pgammap2(y,A3,A4,B2),V)),sigma = matrix(c(1,R,R,1),ncol = 2),df=V)+log(GammaP2(x,A1,A2,B1))-sum(log(dt(qt(pgammap2(x,A1,A2,B1),V),V)),log(dt(qt(pgammap2(y,A3,A4,B2),V),V)))
}
rcopT<-function(n,rho,v){
  V<- rmvt(n,sigma = matrix(c(1,rho,rho,1),ncol=2),df=v)
  U1<- pt(V[,1],df=v); U2<-pt(V[,2],df=v)
  return(matrix(data=c(U1,U2),ncol = 2))
}
logmarg_A1<-function(A,T2,ni){
  (A-1)*T2+dgamma(A,a,b,log = TRUE)-ni*lgamma(A)
}
logmarg_A2<-function(A,T3,T4,ni){
  lgamma(ni*A+a) - ni*lgamma(A) + (A-1)*T4 - (ni*A+a)*log(b+T3)+ dgamma(A,a,b,log = TRUE)  
}
logver_v<-function(V){
  sum(dmvt(matrix(data = c(qt(F1,V),qt(F2,V)),ncol = 2),sigma = Sigma,df=V))+
    -sum(log(dt(qt(F1,V),V))+log(dt(qt(F2,V),V)))+dgamma(V-2,a,b,log = TRUE)
}
logver_rho<-function(p){
  sum(dmvt(matrix(data = c(qt(F1,v_0),qt(F2,v_0)),ncol = 2),sigma = matrix(c(1,p,p,1),ncol = 2),df=v_0))
}

gl<- 4; rho<- 0.7; alpha1<-2; alpha2<-2; beta1<-1; alpha3<-0.5; alpha4<-0.5; beta2<-1

set.seed(7) # seed 7
n=20;  U<-rcopT(n,rho,gl); U1<-U[,1]; U2<-U[,2]; 
theta1<-qgammap2(U1,alpha1,alpha2,beta1); theta2<-qgammap2(U2,alpha3,alpha4,beta2)
a<-1; b<-0.05; cal<-100000; lag<-170; muestra<-1000 ; N=cal+lag*muestra
pmt<-matrix(data = NA,ncol = 8,nrow = muestra); Int20<-matrix(data = NA,ncol = 2,nrow = 8)
alpha1_0<-4; alpha2_0<-4; beta1_0<-4; alpha3_0<-4; alpha4_0<-4; beta2_0<-4; v_0<-4; rho_0<- 0; j<-1
for (i in 1:N) {
  r1<-rgamma(n,sum(alpha1_0+alpha2_0),cos(theta1)+beta1_0*sin(theta1) )
  r2<-rgamma(n,sum(alpha3_0+alpha4_0),cos(theta2)+beta2_0*sin(theta2) )
  
  T2<-sum(log(r1*cos(theta1))); S2<-sum(log(r2*cos(theta2)))
  T3<-sum(r1*sin(theta1))     ; S3<-sum(r2*sin(theta2))
  T4<-sum(log(r1*sin(theta1))); S4<-sum(log(r2*sin(theta2)))
  
  alpha1_0<-arms(1,function(x) logmarg_A1(x,T2,n),0.000001,200,metropolis = TRUE)
  alpha2_0<-arms(1,function(x) logmarg_A2(x,T3,T4,n),0.000001,200,metropolis = TRUE)
  beta1_0 <-rgamma(1,n*alpha2_0+a,T3+b)
  alpha3_0<-arms(1,function(x) logmarg_A1(x,S2,n),0.000001,200,metropolis = TRUE)
  alpha4_0<-arms(1,function(x) logmarg_A2(x,S3,S4,n),0.000001,200,metropolis = TRUE)
  beta2_0 <-rgamma(1,n*alpha4_0+a,S3+b)

  if (i>cal && (cal-i)%%lag==0) {
    pmt[j,1]<-alpha1_0
    pmt[j,2]<-alpha2_0
    pmt[j,3]<-beta1_0
    pmt[j,4]<-alpha3_0
    pmt[j,5]<-alpha4_0
    pmt[j,6]<-beta2_0
    j<-j+1
  }
}
TM2<-3; CAL2<-0; LAG2<-1; N2<-CAL2+TM2*LAG2
for (k in 1:muestra) {
  alpha1_0<-pmt[k,1]
  alpha2_0<-pmt[k,2]
  beta1_0 <-pmt[k,3]
  alpha3_0<-pmt[k,4]
  alpha4_0<-pmt[k,5]
  beta2_0 <-pmt[k,6]
  F1<-pgammap2(theta1,alpha1_0,alpha2_0,beta1_0)
  F2<-pgammap2(theta2,alpha3_0,alpha4_0,beta2_0)
  
  for (i2 in 1:N2) {
    Sigma<-matrix(c(1,rho_0,rho_0,1),ncol = 2)
    v_0  <-arms(1,logver_v,2.0001,50,metropolis = TRUE)
    rho_0<-arms(1,logver_rho,-1,1,metropolis = TRUE)
  }
  
  pmt[k,7]<-v_0
  pmt[k,8]<-rho_0
  
}
for (i in 1:8) {
  Int20[i,]<-c(round(emp.hpd(pmt[,i],0.95),digits = 2))
}
#Credibility Intervals for n=20
#Int20

n=50;  U<-rcopT(n,rho,gl); U1<-U[,1]; U2<-U[,2]; 
theta1<-qgammap2(U1,alpha1,alpha2,beta1); theta2<-qgammap2(U2,alpha3,alpha4,beta2)
a<-1; b<-0.05; cal<-100000; lag<-170; muestra<-1000 ; N=cal+lag*muestra
pmt<-matrix(data = NA,ncol = 8,nrow = muestra); Int50<-matrix(data = NA,ncol = 2,nrow = 8)
alpha1_0<-4; alpha2_0<-4; beta1_0<-4; alpha3_0<-4; alpha4_0<-4; beta2_0<-4; v_0<-4; rho_0<- 0; j<-1
for (i in 1:N) {
  r1<-rgamma(n,sum(alpha1_0+alpha2_0),cos(theta1)+beta1_0*sin(theta1) )
  r2<-rgamma(n,sum(alpha3_0+alpha4_0),cos(theta2)+beta2_0*sin(theta2) )
  
  T2<-sum(log(r1*cos(theta1))); S2<-sum(log(r2*cos(theta2)))
  T3<-sum(r1*sin(theta1))     ; S3<-sum(r2*sin(theta2))
  T4<-sum(log(r1*sin(theta1))); S4<-sum(log(r2*sin(theta2)))
  
  alpha1_0<-arms(1,function(x) logmarg_A1(x,T2,n),0.000001,200,metropolis = TRUE)
  alpha2_0<-arms(1,function(x) logmarg_A2(x,T3,T4,n),0.000001,200,metropolis = TRUE)
  beta1_0 <-rgamma(1,n*alpha2_0+a,T3+b)
  alpha3_0<-arms(1,function(x) logmarg_A1(x,S2,n),0.000001,200,metropolis = TRUE)
  alpha4_0<-arms(1,function(x) logmarg_A2(x,S3,S4,n),0.000001,200,metropolis = TRUE)
  beta2_0 <-rgamma(1,n*alpha4_0+a,S3+b)
  
  if (i>cal && (cal-i)%%lag==0) {
    pmt[j,1]<-alpha1_0
    pmt[j,2]<-alpha2_0
    pmt[j,3]<-beta1_0
    pmt[j,4]<-alpha3_0
    pmt[j,5]<-alpha4_0
    pmt[j,6]<-beta2_0
    j<-j+1
  }
}
TM2<-3; CAL2<-0; LAG2<-1; N2<-CAL2+TM2*LAG2
for (k in 1:muestra) {
  alpha1_0<-pmt[k,1]
  alpha2_0<-pmt[k,2]
  beta1_0 <-pmt[k,3]
  alpha3_0<-pmt[k,4]
  alpha4_0<-pmt[k,5]
  beta2_0 <-pmt[k,6]
  F1<-pgammap2(theta1,alpha1_0,alpha2_0,beta1_0)
  F2<-pgammap2(theta2,alpha3_0,alpha4_0,beta2_0)
  
  for (i2 in 1:N2) {
    Sigma<-matrix(c(1,rho_0,rho_0,1),ncol = 2)
    v_0  <-arms(1,logver_v,2.0001,50,metropolis = TRUE)
    rho_0<-arms(1,logver_rho,-1,1,metropolis = TRUE)
  }
  
  pmt[k,7]<-v_0
  pmt[k,8]<-rho_0
  
}
for (i in 1:8) {
  Int50[i,]<-c(round(emp.hpd(pmt[,i],0.95),digits = 2))
}
#Credibility Intervals for n=50
#Int50

n=100;  U<-rcopT(n,rho,gl); U1<-U[,1]; U2<-U[,2]; 
theta1<-qgammap2(U1,alpha1,alpha2,beta1); theta2<-qgammap2(U2,alpha3,alpha4,beta2)
a<-1; b<-0.05; cal<-100000; lag<-170; muestra<-1000 ; N=cal+lag*muestra
pmt<-matrix(data = NA,ncol = 8,nrow = muestra); Int100<-matrix(data = NA,ncol = 2,nrow = 8)
alpha1_0<-4; alpha2_0<-4; beta1_0<-4; alpha3_0<-4; alpha4_0<-4; beta2_0<-4; v_0<-4; rho_0<- 0; j<-1
for (i in 1:N) {
  r1<-rgamma(n,sum(alpha1_0+alpha2_0),cos(theta1)+beta1_0*sin(theta1) )
  r2<-rgamma(n,sum(alpha3_0+alpha4_0),cos(theta2)+beta2_0*sin(theta2) )
  
  T2<-sum(log(r1*cos(theta1))); S2<-sum(log(r2*cos(theta2)))
  T3<-sum(r1*sin(theta1))     ; S3<-sum(r2*sin(theta2))
  T4<-sum(log(r1*sin(theta1))); S4<-sum(log(r2*sin(theta2)))
  
  alpha1_0<-arms(1,function(x) logmarg_A1(x,T2,n),0.000001,200,metropolis = TRUE)
  alpha2_0<-arms(1,function(x) logmarg_A2(x,T3,T4,n),0.000001,200,metropolis = TRUE)
  beta1_0 <-rgamma(1,n*alpha2_0+a,T3+b)
  alpha3_0<-arms(1,function(x) logmarg_A1(x,S2,n),0.000001,200,metropolis = TRUE)
  alpha4_0<-arms(1,function(x) logmarg_A2(x,S3,S4,n),0.000001,200,metropolis = TRUE)
  beta2_0 <-rgamma(1,n*alpha4_0+a,S3+b)
  
  if (i>cal && (cal-i)%%lag==0) {
    pmt[j,1]<-alpha1_0
    pmt[j,2]<-alpha2_0
    pmt[j,3]<-beta1_0
    pmt[j,4]<-alpha3_0
    pmt[j,5]<-alpha4_0
    pmt[j,6]<-beta2_0
    j<-j+1
  }
}
TM2<-3; CAL2<-0; LAG2<-1; N2<-CAL2+TM2*LAG2
for (k in 1:muestra) {
  alpha1_0<-pmt[k,1]
  alpha2_0<-pmt[k,2]
  beta1_0 <-pmt[k,3]
  alpha3_0<-pmt[k,4]
  alpha4_0<-pmt[k,5]
  beta2_0 <-pmt[k,6]
  F1<-pgammap2(theta1,alpha1_0,alpha2_0,beta1_0)
  F2<-pgammap2(theta2,alpha3_0,alpha4_0,beta2_0)
  
  for (i2 in 1:N2) {
    Sigma<-matrix(c(1,rho_0,rho_0,1),ncol = 2)
    v_0  <-arms(1,logver_v,2.0001,50,metropolis = TRUE)
    rho_0<-arms(1,logver_rho,-1,1,metropolis = TRUE)
  }
  
  pmt[k,7]<-v_0
  pmt[k,8]<-rho_0
  
}
for (i in 1:8) {
  Int100[i,]<-c(round(emp.hpd(pmt[,i],0.95),digits = 2))
}
#Credibility Intervals for n=100
#Int100

n=500;  U<-rcopT(n,rho,gl); U1<-U[,1]; U2<-U[,2]; 
theta1<-qgammap2(U1,alpha1,alpha2,beta1); theta2<-qgammap2(U2,alpha3,alpha4,beta2)
a<-1; b<-0.05; cal<-100000; lag<-170; muestra<-1000 ; N=cal+lag*muestra
pmt<-matrix(data = NA,ncol = 8,nrow = muestra); Int500<-matrix(data = NA,ncol = 2,nrow = 8)
alpha1_0<-4; alpha2_0<-4; beta1_0<-4; alpha3_0<-4; alpha4_0<-4; beta2_0<-4; v_0<-4; rho_0<- 0; j<-1
for (i in 1:N) {
  r1<-rgamma(n,sum(alpha1_0+alpha2_0),cos(theta1)+beta1_0*sin(theta1) )
  r2<-rgamma(n,sum(alpha3_0+alpha4_0),cos(theta2)+beta2_0*sin(theta2) )
  
  T2<-sum(log(r1*cos(theta1))); S2<-sum(log(r2*cos(theta2)))
  T3<-sum(r1*sin(theta1))     ; S3<-sum(r2*sin(theta2))
  T4<-sum(log(r1*sin(theta1))); S4<-sum(log(r2*sin(theta2)))
  
  alpha1_0<-arms(1,function(x) logmarg_A1(x,T2,n),0.000001,200,metropolis = TRUE)
  alpha2_0<-arms(1,function(x) logmarg_A2(x,T3,T4,n),0.000001,200,metropolis = TRUE)
  beta1_0 <-rgamma(1,n*alpha2_0+a,T3+b)
  alpha3_0<-arms(1,function(x) logmarg_A1(x,S2,n),0.000001,200,metropolis = TRUE)
  alpha4_0<-arms(1,function(x) logmarg_A2(x,S3,S4,n),0.000001,200,metropolis = TRUE)
  beta2_0 <-rgamma(1,n*alpha4_0+a,S3+b)
  
  if (i>cal && (cal-i)%%lag==0) {
    pmt[j,1]<-alpha1_0
    pmt[j,2]<-alpha2_0
    pmt[j,3]<-beta1_0
    pmt[j,4]<-alpha3_0
    pmt[j,5]<-alpha4_0
    pmt[j,6]<-beta2_0
    j<-j+1
  }
}
TM2<-3; CAL2<-0; LAG2<-1; N2<-CAL2+TM2*LAG2
for (k in 1:muestra) {
  alpha1_0<-pmt[k,1]
  alpha2_0<-pmt[k,2]
  beta1_0 <-pmt[k,3]
  alpha3_0<-pmt[k,4]
  alpha4_0<-pmt[k,5]
  beta2_0 <-pmt[k,6]
  F1<-pgammap2(theta1,alpha1_0,alpha2_0,beta1_0)
  F2<-pgammap2(theta2,alpha3_0,alpha4_0,beta2_0)
  
  for (i2 in 1:N2) {
    Sigma<-matrix(c(1,rho_0,rho_0,1),ncol = 2)
    v_0  <-arms(1,logver_v,2.0001,50,metropolis = TRUE)
    rho_0<-arms(1,logver_rho,-1,1,metropolis = TRUE)
  }
  
  pmt[k,7]<-v_0
  pmt[k,8]<-rho_0
  
}
for (i in 1:8) {
  Int500[i,]<-c(round(emp.hpd(pmt[,i],0.95),digits = 2))
}
#Credibility Intervals for n=50
#Int500

hist(pmt[,1],xlab = expression(alpha[11]),col = 'red',main = '',border = 'white',freq = F,breaks = 8)
hist(pmt[,2],xlab = expression(alpha[12]),col = 'red',main = '',border = 'white',freq = F,breaks = 8)
hist(pmt[,3],xlab = expression(beta[1]),col = 'red',main = '',border = 'white',freq = F,breaks = 8)
hist(pmt[,4],xlab = expression(alpha[21]),col = 'red',main = '',border = 'white',freq = F,breaks = 9) 
hist(pmt[,5],xlab = expression(alpha[22]),col = 'red',main = '',border = 'white',freq = F,breaks = 8)
hist(pmt[,6],xlab = expression(beta[2]),col = 'red',main = '',border = 'white',freq = F,breaks = 8)
hist(pmt[,7],xlab = expression(nu),col = 'red',main = '',border = 'white',freq = F,breaks = 8)
hist(pmt[,8],xlab = expression(rho[12]),col = 'red',main = '',border = 'white',freq = F,breaks = 8)


ConfidenceIntervals<-cbind(Int20,Int50,Int100,Int500)

ConfidenceIntervals
geweke.diag(pmt)
# 1.27170 -0.09260 -0.60030  0.09476 -0.24904 -0.03074  0.61097 -0.04343 
ess(pmt)
# 1251.683 1000.000 1000.000 1000.000 1000.000  866.345 1000.000 1128.055
