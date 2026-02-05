# LIBRARIES and FUNCTIONS ####
library(mvtnorm)
library(LaplacesDemon)
library(TeachingDemos)
library(armspp)
library(TeachingDemos)
library(tictoc)
library(mvtnorm)
library(rgl)
library(mvnormalTest)
library(psych)

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
rcopG<-function(N,SIGMA) {
  V<-rmvnorm(N,sigma = SIGMA)
  U<-pnorm(V)
  return(U)
}
logmarg_A11<-function(A){
  (A-1)*T2+dgamma(A,a,b,log = TRUE)-n*lgamma(A)
}
logmarg_A21<-function(A){
  (A-1)*S2+dgamma(A,a,b,log = TRUE)-n*lgamma(A)
}
logmarg_A31<-function(A){
  (A-1)*R2+dgamma(A,a,b,log = TRUE)-n*lgamma(A)
}
logmarg_A12<-function(A){
  lgamma(n*A+a) - n*lgamma(A) + (A-1)*T4 - (n*A+a)*log(b+T3)+ dgamma(A,a,b,log = TRUE)  
}
logmarg_A22<-function(A){
  lgamma(n*A+a) - n*lgamma(A) + (A-1)*S4 - (n*A+a)*log(b+S3)+ dgamma(A,a,b,log = TRUE)  
}
logmarg_A32<-function(A){
  lgamma(n*A+a) - n*lgamma(A) + (A-1)*R4 - (n*A+a)*log(b+R3)+ dgamma(A,a,b,log = TRUE)  
}

# - ####

set.seed(11) #seed: 11

rho12<-  0.75; rho13<-  -0.75; rho23<-  -0.75; S<-matrix(c(1,rho12,rho13,rho12,1,rho23,rho13,rho23,1),ncol = 3)
n=500 ; U<-rcopG(n,S); U1<-U[,1]; U2<-U[,2]; U3<-U[,3]
alpha11<-2; alpha12<-2; beta1<-2; alpha21<-0.5; alpha22<-3; beta2<-1; alpha31<-3; alpha32<-5; beta3<-3
theta1<-qgammap2(U1,alpha11,alpha12,beta1); theta2<-qgammap2(U2,alpha21,alpha22,beta2); theta3<-qgammap2(U3,alpha31,alpha32,beta3)
plot3d(theta1,theta2,theta3,type = 's',col='red',size = 1,xlab = expression(theta[1]),ylab = expression(theta[2]),zlab = expression(theta[3]))

DatMat<-matrix(data = c(cos(theta1),sin(theta2)*cos(theta1),sin(theta1)*sin(theta2)*cos(theta3),sin(theta1)*sin(theta2)*sin(theta3)),ncol = 4)
mvnormalTest::mardia(DatMat)$mv.test
psych::mardia(DatMat)
muestra=1000

cal<-50000; lag<-50; N=cal+lag*muestra
alpha11_0<-4; alpha12_0<-4; beta1_0<-4; alpha21_0<-4; alpha22_0<-4; beta2_0<-4; alpha31_0<-4; alpha32_0<-4; beta3_0<-4
a<-1; b<-0.5; 
MuestraFinal<-matrix(data = NA,ncol = muestra,nrow = 12)
Intervalos.5<-matrix(data = NA,ncol = 2,nrow = 12)
j<-1
for (i in 1:N) {
  r1<-rgamma(n,sum(alpha11_0+alpha12_0),cos(theta1)+beta1_0*sin(theta1) )
  r2<-rgamma(n,sum(alpha21_0+alpha22_0),cos(theta2)+beta2_0*sin(theta2) )
  r3<-rgamma(n,sum(alpha31_0+alpha32_0),cos(theta3)+beta3_0*sin(theta3) )
  
  T2<-sum(log(r1*cos(theta1))); S2<-sum(log(r2*cos(theta2))); R2<-sum(log(r3*cos(theta3))); 
  T3<-sum(r1*sin(theta1))     ; S3<-sum(r2*sin(theta2));      R3<-sum(r3*sin(theta3))
  T4<-sum(log(r1*sin(theta1))); S4<-sum(log(r2*sin(theta2))); R4<-sum(log(r3*sin(theta3)));
  
  alpha11_0<-arms(1,logmarg_A11,0.000001,200,metropolis = TRUE)
  alpha12_0<-arms(1,logmarg_A12,0.000001,200,metropolis = TRUE)
  beta1_0 <-rgamma(1,n*alpha12_0+a,T3+b)
  alpha21_0<-arms(1,logmarg_A21,0.000001,200,metropolis = TRUE)
  alpha22_0<-arms(1,logmarg_A22,0.000001,200,metropolis = TRUE)
  beta2_0 <-rgamma(1,n*alpha22_0+a,S3+b)
  alpha31_0<-arms(1,logmarg_A31,0.000001,200,metropolis = TRUE)
  alpha32_0<-arms(1,logmarg_A32,0.000001,200,metropolis = TRUE)
  beta3_0 <-rgamma(1,n*alpha32_0+a,R3+b)
  if (i>cal && (cal-i)%%lag==0) {
    MuestraFinal[1,j]<-alpha11_0
    MuestraFinal[2,j]<-alpha12_0
    MuestraFinal[3,j]<-beta1_0
    MuestraFinal[4,j]<-alpha21_0
    MuestraFinal[5,j]<-alpha22_0
    MuestraFinal[6,j]<-beta2_0
    MuestraFinal[7,j]<-alpha31_0
    MuestraFinal[8,j]<-alpha32_0
    MuestraFinal[9,j]<-beta3_0
    j<-j+1
  }
}
nu0<-4
Sigma0<-diag(0.0001,nrow = 3)
for (k in 1:muestra) {
  alpha11_0<-MuestraFinal[1,k]
  alpha12_0<-MuestraFinal[2,k]
  beta1_0  <-MuestraFinal[3,k]
  alpha21_0<-MuestraFinal[4,k]
  alpha22_0<-MuestraFinal[5,k]
  beta2_0  <-MuestraFinal[6,k]
  alpha31_0<-MuestraFinal[7,k]
  alpha32_0<-MuestraFinal[8,k]
  beta3_0  <-MuestraFinal[9,k]
  
  F1<-pgammap2(theta1,alpha11_0,alpha12_0,beta1_0)
  F2<-pgammap2(theta2,alpha21_0,alpha22_0,beta2_0)
  F3<-pgammap2(theta3,alpha31_0,alpha32_0,beta3_0)
  
  Y1<-qnorm(F1); Y2<-qnorm(F2); Y3<-qnorm(F3)
  
  S.mu<-Sigma0
  for (j in 1:n) {
    S.mu<-S.mu+c(Y1[j],Y2[j],Y3[j])%*%t(c(Y1[j],Y2[j],Y3[j]))
  }
  
  R.posterior<-cov2cor(rinvwishart(nu0+n,S.mu))
  
  MuestraFinal[10,k]<-R.posterior[1,2]
  MuestraFinal[11,k]<-R.posterior[1,3]
  MuestraFinal[12,k]<-R.posterior[2,3]
}
for (i in 1:12) {
  Intervalos.5[i,]<-round(emp.hpd(MuestraFinal[i,],0.95),digits = 2)
}

a<-1; b<-0.05;
MuestraFinal<-matrix(data = NA,ncol = muestra,nrow = 12)
Intervalos.05<-matrix(data = NA,ncol = 2,nrow = 12)
j<-1
for (i in 1:N) {
  r1<-rgamma(n,sum(alpha11_0+alpha12_0),cos(theta1)+beta1_0*sin(theta1) )
  r2<-rgamma(n,sum(alpha21_0+alpha22_0),cos(theta2)+beta2_0*sin(theta2) )
  r3<-rgamma(n,sum(alpha31_0+alpha32_0),cos(theta3)+beta3_0*sin(theta3) )
  
  T2<-sum(log(r1*cos(theta1))); S2<-sum(log(r2*cos(theta2))); R2<-sum(log(r3*cos(theta3))); 
  T3<-sum(r1*sin(theta1))     ; S3<-sum(r2*sin(theta2));      R3<-sum(r3*sin(theta3))
  T4<-sum(log(r1*sin(theta1))); S4<-sum(log(r2*sin(theta2))); R4<-sum(log(r3*sin(theta3)));
  
  alpha11_0<-arms(1,logmarg_A11,0.000001,200,metropolis = TRUE)
  alpha12_0<-arms(1,logmarg_A12,0.000001,200,metropolis = TRUE)
  beta1_0 <-rgamma(1,n*alpha12_0+a,T3+b)
  alpha21_0<-arms(1,logmarg_A21,0.000001,200,metropolis = TRUE)
  alpha22_0<-arms(1,logmarg_A22,0.000001,200,metropolis = TRUE)
  beta2_0 <-rgamma(1,n*alpha22_0+a,S3+b)
  alpha31_0<-arms(1,logmarg_A31,0.000001,200,metropolis = TRUE)
  alpha32_0<-arms(1,logmarg_A32,0.000001,200,metropolis = TRUE)
  beta3_0 <-rgamma(1,n*alpha32_0+a,R3+b)
  if (i>cal && (cal-i)%%lag==0) {
    MuestraFinal[1,j]<-alpha11_0
    MuestraFinal[2,j]<-alpha12_0
    MuestraFinal[3,j]<-beta1_0
    MuestraFinal[4,j]<-alpha21_0
    MuestraFinal[5,j]<-alpha22_0
    MuestraFinal[6,j]<-beta2_0
    MuestraFinal[7,j]<-alpha31_0
    MuestraFinal[8,j]<-alpha32_0
    MuestraFinal[9,j]<-beta3_0
    j<-j+1
  }
}
nu0<-4
Sigma0<-diag(0.0001,nrow = 3)
for (k in 1:muestra) {
  alpha11_0<-MuestraFinal[1,k]
  alpha12_0<-MuestraFinal[2,k]
  beta1_0  <-MuestraFinal[3,k]
  alpha21_0<-MuestraFinal[4,k]
  alpha22_0<-MuestraFinal[5,k]
  beta2_0  <-MuestraFinal[6,k]
  alpha31_0<-MuestraFinal[7,k]
  alpha32_0<-MuestraFinal[8,k]
  beta3_0  <-MuestraFinal[9,k]
  
  F1<-pgammap2(theta1,alpha11_0,alpha12_0,beta1_0)
  F2<-pgammap2(theta2,alpha21_0,alpha22_0,beta2_0)
  F3<-pgammap2(theta3,alpha31_0,alpha32_0,beta3_0)
  
  Y1<-qnorm(F1); Y2<-qnorm(F2); Y3<-qnorm(F3)
  
  S.mu<-Sigma0
  for (j in 1:n) {
    S.mu<-S.mu+c(Y1[j],Y2[j],Y3[j])%*%t(c(Y1[j],Y2[j],Y3[j]))
  }
  
  R.posterior<-cov2cor(rinvwishart(nu0+n,S.mu))
  
  MuestraFinal[10,k]<-R.posterior[1,2]
  MuestraFinal[11,k]<-R.posterior[1,3]
  MuestraFinal[12,k]<-R.posterior[2,3]
}
for (i in 1:12) {
  Intervalos.05[i,]<-round(emp.hpd(MuestraFinal[i,],0.95),digits = 2)
}

hist(MuestraFinal[1,],xlab = expression(alpha[11]),col = 'red',main = '',border='white',freq=F,breaks = 10)
hist(MuestraFinal[2,],xlab = expression(alpha[12]),col = 'red',main = '',border='white',freq=F,breaks = 10)
hist(MuestraFinal[3,],xlab = expression(beta[1]),col = 'red',main = '',border='white',freq=F,breaks = 10)
hist(MuestraFinal[4,],xlab = expression(alpha[21]),col = 'red',main = '',border='white',freq=F,breaks = 10)
hist(MuestraFinal[5,],xlab = expression(alpha[22]),col = 'red',main = '',border='white',freq=F,breaks = 10)
hist(MuestraFinal[6,],xlab = expression(beta[2]),col = 'red',main = '',border='white',freq=F,breaks = 10)
hist(MuestraFinal[7,],xlab = expression(alpha[31]),col = 'red',main = '',border='white',freq=F,breaks = 10)
hist(MuestraFinal[8,],xlab = expression(alpha[32]),col = 'red',main = '',border='white',freq=F,breaks = 10)
hist(MuestraFinal[9,],xlab = expression(beta[3]),col = 'red',main = '',border='white',freq=F,breaks = 10)
hist(MuestraFinal[10,],xlab = expression(rho[12]),col = 'red',main = '',border='white',freq=F,breaks = 10)
hist(MuestraFinal[11,],xlab = expression(rho[13]),col = 'red',main = '',border='white',freq=F,breaks = 10)
hist(MuestraFinal[12,],xlab = expression(rho[23]),col = 'red',main = '',border='white',freq=F,breaks = 10)

#a<-1; b<-0.01;
#MuestraFinal<-matrix(data = NA,ncol = muestra,nrow = 12)
#Intervalos.01<-matrix(data = NA,ncol = 2,nrow = 12)
#j<-1
#for (i in 1:N) {
  r1<-rgamma(n,sum(alpha11_0+alpha12_0),cos(theta1)+beta1_0*sin(theta1) )
  r2<-rgamma(n,sum(alpha21_0+alpha22_0),cos(theta2)+beta2_0*sin(theta2) )
  r3<-rgamma(n,sum(alpha31_0+alpha32_0),cos(theta3)+beta3_0*sin(theta3) )
  
  T2<-sum(log(r1*cos(theta1))); S2<-sum(log(r2*cos(theta2))); R2<-sum(log(r3*cos(theta3))); 
  T3<-sum(r1*sin(theta1))     ; S3<-sum(r2*sin(theta2));      R3<-sum(r3*sin(theta3))
  T4<-sum(log(r1*sin(theta1))); S4<-sum(log(r2*sin(theta2))); R4<-sum(log(r3*sin(theta3)));
  
  alpha11_0<-arms(1,logmarg_A11,0.000001,200,metropolis = TRUE)
  alpha12_0<-arms(1,logmarg_A12,0.000001,200,metropolis = TRUE)
  beta1_0 <-rgamma(1,n*alpha12_0+a,T3+b)
  alpha21_0<-arms(1,logmarg_A21,0.000001,200,metropolis = TRUE)
  alpha22_0<-arms(1,logmarg_A22,0.000001,200,metropolis = TRUE)
  beta2_0 <-rgamma(1,n*alpha22_0+a,S3+b)
  alpha31_0<-arms(1,logmarg_A31,0.000001,200,metropolis = TRUE)
  alpha32_0<-arms(1,logmarg_A32,0.000001,200,metropolis = TRUE)
  beta3_0 <-rgamma(1,n*alpha32_0+a,R3+b)
  if (i>cal && (cal-i)%%lag==0) {
    MuestraFinal[1,j]<-alpha11_0
    MuestraFinal[2,j]<-alpha12_0
    MuestraFinal[3,j]<-beta1_0
    MuestraFinal[4,j]<-alpha21_0
    MuestraFinal[5,j]<-alpha22_0
    MuestraFinal[6,j]<-beta2_0
    MuestraFinal[7,j]<-alpha31_0
    MuestraFinal[8,j]<-alpha32_0
    MuestraFinal[9,j]<-beta3_0
    j<-j+1
  }
}
#nu0<-4
#Sigma0<-diag(0.0001,nrow = 3)
#for (k in 1:muestra) {
  alpha11_0<-MuestraFinal[1,k]
  alpha12_0<-MuestraFinal[2,k]
  beta1_0  <-MuestraFinal[3,k]
  alpha21_0<-MuestraFinal[4,k]
  alpha22_0<-MuestraFinal[5,k]
  beta2_0  <-MuestraFinal[6,k]
  alpha31_0<-MuestraFinal[7,k]
  alpha32_0<-MuestraFinal[8,k]
  beta3_0  <-MuestraFinal[9,k]
  
  F1<-pgammap2(theta1,alpha11_0,alpha12_0,beta1_0)
  F2<-pgammap2(theta2,alpha21_0,alpha22_0,beta2_0)
  F3<-pgammap2(theta3,alpha31_0,alpha32_0,beta3_0)
  
  Y1<-qnorm(F1); Y2<-qnorm(F2); Y3<-qnorm(F3)
  
  S.mu<-Sigma0
  for (j in 1:n) {
    S.mu<-S.mu+c(Y1[j],Y2[j],Y3[j])%*%t(c(Y1[j],Y2[j],Y3[j]))
  }
  
  R.posterior<-cov2cor(rinvwishart(nu0+n,S.mu))
  
  MuestraFinal[10,k]<-R.posterior[1,2]
  MuestraFinal[11,k]<-R.posterior[1,3]
  MuestraFinal[12,k]<-R.posterior[2,3]
}
#for (i in 1:12) {
  Intervalos.01[i,]<-round(emp.hpd(MuestraFinal[i,],0.95),digits = 2)
}

#Intervalos<-cbind(c(alpha11,alpha12,beta1,alpha21,alpha22,beta2,alpha31,alpha32,beta3,rho12,rho13,rho23),Intervalos.5,Intervalos.05,Intervalos.01)
