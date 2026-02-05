## ####

library(armspp)
library(TeachingDemos)
library(tictoc)
library(mvtnorm)
library(ks)

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
rcopG<-function(n,S) {
  V<-rmvnorm(n,sigma = S)
  U<-pnorm(V)
  return(U)
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
logdctgammay.x<-function(x,y,A1,A2,B1,A3,A4,B2,V,R){
  dmvt(c(qt(pgammap2(x,A1,A2,B1),V),qt(pgammap2(y,A3,A4,B2),V)),sigma = matrix(c(1,R,R,1),ncol = 2),df=V)+log(GammaP2(y,A3,A4,B2))-sum(log(dt(qt(pgammap2(x,A1,A2,B1),V),V)),log(dt(qt(pgammap2(y,A3,A4,B2),V),V)))
}
logdctgammax.y<-function(x,y,A1,A2,B1,A3,A4,B2,V,R){
  dmvt(c(qt(pgammap2(x,A1,A2,B1),V),qt(pgammap2(y,A3,A4,B2),V)),sigma = matrix(c(1,R,R,1),ncol = 2),df=V)+log(GammaP2(x,A1,A2,B1))-sum(log(dt(qt(pgammap2(x,A1,A2,B1),V),V)),log(dt(qt(pgammap2(y,A3,A4,B2),V),V)))
}
logdcggammay.x<-function(x,y,A1,A2,B1,A3,A4,B2,R){
  dmvnorm(c(qnorm(pgammap2(x,A1,A2,B1)),qnorm(pgammap2(y,A3,A4,B2))),sigma = matrix(c(1,R,R,1),ncol = 2),log = TRUE)+log(GammaP2(y,A3,A4,B2))-sum(log(dnorm(qnorm(pgammap2(x,A1,A2,B1)))),log(dnorm(qnorm(pgammap2(y,A3,A4,B2)))))
}
dcggp<-function(x,y,A1,A2,B1,A3,A4,B2,R){
  dmvnorm(c(qnorm(pgammap2(x,A1,A2,B1)),qnorm(pgammap2(y,A3,A4,B2))),sigma = matrix(c(1,R,R,1),ncol = 2))*GammaP2(x,A1,A2,B1)*GammaP2(y,A3,A4,B2)/(dnorm(qnorm(pgammap2(x,A1,A2,B1)))*dnorm(qnorm(pgammap2(y,A3,A4,B2))))
}
dctgp<-function(x,y,A1,A2,B1,A3,A4,B2,V,R){
  dmvt(c(qt(pgammap2(x,A1,A2,B1),V),qt(pgammap2(y,A3,A4,B2),V)),sigma = matrix(c(1,R,R,1),ncol = 2),df=V,log = FALSE)*GammaP2(x,A1,A2,B1)*GammaP2(y,A3,A4,B2)/(dt(qt(pgammap2(x,A1,A2,B1),V),V)*dt(qt(pgammap2(y,A3,A4,B2),V),V))
}
## Principal Program ####
#SubjectData <- read.csv("~/Downloads/SubjectData.csv")
SubjectData <- read.csv("~/Library/CloudStorage/Dropbox/Escuela/Paper GammaP Multivariado/SubjectData.csv")

theta1<-(SubjectData$LFA[-c(1,18)]*pi/180)[1:73]
theta2<-(SubjectData$RFA[-c(1,18)]*pi/180)[1:73]

n<-length(theta1)
plot(theta1,theta2,pch=16,col=rgb(1.0,0.0,0.0,0.9),xlab = 'LFA',ylab = 'RFA')

ejex<-seq(0.001,pi/2-0.001,length=1000)

set.seed(1) # 1
# INFERENCIA CON COPULA T ####

a<-1; b<-0.05; cal<-250000; lag<-50; muestra<-1000; N=cal+lag*muestra

pmt<-matrix(data = NA,ncol = 8,nrow = muestra)
IntCt<-matrix(data = NA,ncol = 2,nrow = 8)

alpha1_0<-4; alpha2_0<-4; beta1_0<-4; alpha3_0<-4; alpha4_0<-4; beta2_0<-4
v_0<-4; rho_0<-0; j<-1

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
  
  #Sigma<-matrix(c(1,rho_0,rho_0,1),ncol = 2)
  #v_0<-arms(1,logver_v,2.0001,50,metropolis = TRUE)
  #rho_0<-arms(1,logver_rho,-1,1,metropolis = TRUE)
  
  for (i2 in 1:N2) {
    Sigma<-matrix(c(1,rho_0,rho_0,1),ncol = 2)
    v_0  <-arms(1,logver_v,2.0001,50,metropolis = TRUE)
    rho_0<-arms(1,logver_rho,-1,1,metropolis = TRUE)
    #if(i2>CAL2 && (i2-CAL2)%%LAG2==0){
    #  rmt[j,1]<-v_0
    #  rmt[j,2]<-rho_0
    #  j<-j+1
    #}
  }
  
  pmt[k,7]<-v_0
  pmt[k,8]<-rho_0
  
}
for (i in 1:8) {
  IntCt[i,]<-c(emp.hpd(pmt[,i],0.95))
}  

f_q1<-rep(0,length=length(ejex));f_q2<-rep(0,length=length(ejex));f_m<-rep(0,length=length(ejex))
g_q1<-rep(0,length=length(ejex));g_q2<-rep(0,length=length(ejex));g_m<-rep(0,length=length(ejex))
q1<-rep(NA,muestra); q2<-rep(NA,muestra)

for (i in 1:length(ejex)) {
  for (l in 1:muestra) {
    q1[l]<-GammaP2(ejex[i],pmt[l,1] , pmt[l,2],pmt[l,3])
    q2[l]<-GammaP2(ejex[i],pmt[l,4] , pmt[l,5],pmt[l,6])
  }
  f_q1[i]<-quantile( q1,0.025 )
  f_q2[i]<-quantile( q1,0.975 )
  f_m[i]<-quantile( q1,0.5 )
  g_q1[i]<-quantile( q2,0.025 )
  g_q2[i]<-quantile( q2,0.975 )
  g_m[i]<-quantile( q2,0.5 )
}

hist(theta1,freq = FALSE,border = 'white',col = 'gray70',
     xlim = c(0,pi/2),main='',xlab = expression(theta[1]))
lines(ejex,f_q1,col='blue',lwd=2,lty=2)
lines(ejex,f_q2,col='blue',lwd=2,lty=2)
lines(ejex,f_m,col='red',lwd=2)

hist(theta2,freq = FALSE,border = 'white',col = 'gray70',
     xlim = c(0,pi/2),breaks = 7,main='',xlab = expression(theta[2]))
lines(ejex,g_q1,col='blue',lwd=2,lty=2)
lines(ejex,g_q2,col='blue',lwd=2,lty=2)
lines(ejex,g_m,col='red',lwd=2)

n.pred<-500
ind.pred<-sample(1000,n.pred)
theta1.predct<-rep(NA,n.pred)
theta2.predct<-rep(NA,n.pred)

for (i in 1:n.pred) {
  u<-rcopT(1,pmt[ind.pred[i],8],pmt[ind.pred[i],7])
  theta1.predct[i]<-qgammap2(u[1],pmt[ind.pred[i],1],pmt[ind.pred[i],2],pmt[ind.pred[i],3]) 
  theta2.predct[i]<-qgammap2(u[2],pmt[ind.pred[i],4],pmt[ind.pred[i],5],pmt[ind.pred[i],6])
}

df <- data.frame(x = theta1, y = theta2)
df.predct <- data.frame(x = theta1.predct, y = theta2.predct)

cpo.t<-rep(NA,n)
for (i in 1:n) {
  cpo.t[i]<-0
  for (j in 1:muestra) {
    cpo.t[i]<-cpo.t[i]+1/dctgp(theta1[i],theta2[i],pmt[j,1],pmt[j,2],pmt[j,3],pmt[j,4],pmt[j,5],pmt[j,6],pmt[j,7],pmt[j,8])
  }
  cpo.t[i]<-n/cpo.t[i]
}
lpml.t<-sum(log(cpo.t))

rho.t<-pmt[,8]
tau.t<-2*asin(rho.t)/pi

emp.hpd(tau.t)

Mpred_t<-matrix(data = NA,nrow = n,ncol = muestra)

lppd_ct<-0
pwaic2_ct<-0
for (i in 1:n) {
  for (j in 1:muestra) {
    Mpred_t[i,j]<-dctgp(theta1[i],theta2[i],pmt[j,1],pmt[j,2],pmt[j,3],pmt[j,4],pmt[j,5],pmt[j,6],pmt[j,7],pmt[j,8])
  }
  lppd_ct<-lppd_ct+log(mean(Mpred_t[i,]))
  pwaic2_ct<-pwaic2_ct+var(log(Mpred_t[i,]))
}

lppd_ct
pwaic2_ct

waic_ct<- -2*(lppd_ct-pwaic2_ct)

# INFERENCIA CON COPULA GAUSSIANA ####

a<-1; b<-0.05
cal<-100000; lag<-50; muestra<-1000; N=cal+lag*muestra

pmt<-matrix(data = NA,ncol = 7,nrow = muestra)
IntCg<-matrix(data = NA,ncol = 2,nrow = 8)

alpha1_0<-4; alpha2_0<-4; beta1_0<-4; alpha3_0<-4; alpha4_0<-4; beta2_0<-4
rho_0<- 0
j<-1

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
for (k in 1:muestra) {
  alpha1_0<-pmt[k,1]
  alpha2_0<-pmt[k,2]
  beta1_0 <-pmt[k,3]
  alpha3_0<-pmt[k,4]
  alpha4_0<-pmt[k,5]
  beta2_0 <-pmt[k,6]
  F1<-pgammap2(theta1,alpha1_0,alpha2_0,beta1_0)
  F2<-pgammap2(theta2,alpha3_0,alpha4_0,beta2_0)
  
  logver_rho<-function(p){
    sum(dmvnorm( matrix(data = c(qnorm(F1),qnorm(F2)),ncol = 2),mean = c(0,0), sigma = matrix(c(1,p,p,1),ncol = 2),log=TRUE))
  }
  
  rho.0<-arms(1,logver_rho,-1,1,metropolis = TRUE)
  
  pmt[k,7]<-rho.0
  
}

for (i in 1:7) {
  IntCg[i,]<-c(emp.hpd(pmt[,i],0.95))
}  

## Prediciivas marginales

f_q1<-rep(0,length=length(ejex));f_q2<-rep(0,length=length(ejex));f_m<-rep(0,length=length(ejex))
g_q1<-rep(0,length=length(ejex));g_q2<-rep(0,length=length(ejex));g_m<-rep(0,length=length(ejex))
q1<-rep(NA,muestra); q2<-rep(NA,muestra)

for (i in 1:length(ejex)) {
  for (l in 1:muestra) {
    q1[l]<-GammaP2(ejex[i],pmt[l,1] , pmt[l,2],pmt[l,3])
    q2[l]<-GammaP2(ejex[i],pmt[l,4] , pmt[l,5],pmt[l,6])
  }
  f_q1[i]<-quantile( q1,0.025 )
  f_q2[i]<-quantile( q1,0.975 )
  f_m[i]<-quantile( q1,0.5 )
  g_q1[i]<-quantile( q2,0.025 )
  g_q2[i]<-quantile( q2,0.975 )
  g_m[i]<-quantile( q2,0.5 )
}

hist(theta1,freq = FALSE,border = 'white',col = 'gray70',
     xlim = c(0,pi/2),main='',xlab = expression(theta[1]))
lines(ejex,f_q1,col='blue',lwd=2,lty=2)
lines(ejex,f_q2,col='blue',lwd=2,lty=2)
lines(ejex,f_m,col='red',lwd=2)

hist(theta2,freq = FALSE,border = 'white',col = 'gray70',
     xlim = c(0,pi/2),breaks = 7,main='',xlab = expression(theta[2]))
lines(ejex,g_q1,col='blue',lwd=2,lty=2)
lines(ejex,g_q2,col='blue',lwd=2,lty=2)
lines(ejex,g_m,col='red',lwd=2)


#f_q1<-rep(0,length=length(ejex)); f_q2<-rep(0,length=length(ejex))
#g_q1<-rep(0,length=length(ejex)); g_q2<-rep(0,length=length(ejex))
#q1<-rep(NA,muestra); q2<-rep(NA,muestra)
#for (i in 1:length(ejex)) {
#  for (l in 1:muestra) {
#    q1[l]<-GammaP2(ejex[i],pmt[l,1] , pmt[l,2],pmt[l,3])
#    q2[l]<-GammaP2(ejex[i],pmt[l,4] , pmt[l,5],pmt[l,6])
#  }
#  f_q1[i]<-quantile( q1,0.025 )
#  f_q2[i]<-quantile( q1,0.975 )
#  g_q1[i]<-quantile( q2,0.025 )
#  g_q2[i]<-quantile( q2,0.975 )
#}

#hist(theta1,freq = FALSE,border = 'white',col = 'gray80',xlim = c(0,pi/2))#,ylim = c(0,6))
#lines(ejex,f_q1,col='blue',lwd=3)
#lines(ejex,f_q2,col='blue',lwd=3)

#hist(theta2,freq = FALSE,border = 'white',col = 'gray80',xlim = c(0,pi/2))#,ylim = c(0,8))
#lines(ejex,g_q1,col='blue',lwd=3)
#lines(ejex,g_q2,col='blue',lwd=3)

n.pred<-500
ind.pred<-sample(1000,n.pred)

theta1.predcg<-rep(NA,n.pred); theta2.predcg<-rep(NA,n.pred)

for (i in 1:n.pred) {
  u<-rcopG(1,matrix(c(1,pmt[ind.pred[i],7],pmt[ind.pred[i],7],1),ncol = 2))
  theta1.predcg[i]<-qgammap2(u[1],pmt[ind.pred[i],1],pmt[ind.pred[i],2],pmt[ind.pred[i],3]) 
  theta2.predcg[i]<-qgammap2(u[2],pmt[ind.pred[i],4],pmt[ind.pred[i],5],pmt[ind.pred[i],6])
}

df.predcg <- data.frame(x = theta1.predcg, y = theta2.predcg)

cpo.g<-rep(NA,n)
for (i in 1:n) {
  cpo.g[i]<-0
  for (j in 1:muestra) {
    cpo.g[i]<-cpo.g[i]+1/dcggp(theta1[i],theta2[i],pmt[j,1],pmt[j,2],pmt[j,3],pmt[j,4],pmt[j,5],pmt[j,6],pmt[j,7])
  }
  cpo.g[i]<-n/cpo.g[i]
}
lpml.g<-sum(log(cpo.g))

rho.g<-pmt[,7]
tau.g<-2*asin(rho.g)/pi

hist(tau.g)
emp.hpd(tau.g)

Mpred_g<-matrix(data = NA,nrow = n,ncol = muestra)

lppd_cg<-0
pwaic2_cg<-0
for (i in 1:n) {
  for (j in 1:muestra) {
    Mpred_g[i,j]<-dcggp(theta1[i],theta2[i],pmt[j,1],pmt[j,2],pmt[j,3],pmt[j,4],pmt[j,5],pmt[j,6],pmt[j,7])
  }
  lppd_cg<-lppd_cg+log(mean(Mpred_g[i,]))
  pwaic2_cg<-pwaic2_cg+var(log(Mpred_g[i,]))
}

lppd_cg
pwaic2_cg

waic_cg<- -2*(lppd_cg-pwaic2_cg)

# COMPARACIONES ####

mean(tau.g)
mean(tau.t)

emp.hpd(tau.t)
emp.hpd(tau.g)

lpml.t #-59.633
lpml.g #-60.3493

waic_ct #-262.8532
waic_cg #-261.5193

# INTERVALOS

round(cbind(IntCg,IntCt),digits = 2)

# PREDICTIVAS KDE ####

H1.ct <- Hpi(x=df.predct) # H1.ct<-matrix(c(0.01,0.001,0.001,0.02),nrow = 2)
H2.ct <- Hpi.diag(x=df.predct)

fhat1.ct <- kde(x=df.predct, H=H1.ct,compute.cont=T)
fhat2.ct <- kde(x=df.predct, H=H2.ct)

plot(fhat1.ct,xlab = 'LFA',ylab = 'RFA',main='CopT, H no diagonal',xlim = c(0.1,0.9),ylim = c(0.15,0.85),lwd=3,display="filled.contour",cont=c(10,50,75,95),drawlabels=F,col=c('white','gray90','gray70','gray45','gray25')) #points(df.predct)
points(theta1,theta2,pch=16,col='red')

plot(fhat2.ct,xlab = 'LFA',ylab = 'RFA',main='CopT, H diagonal',xlim = c(0.15,0.8),ylim = c(0.15,0.8),lwd=3,display="filled.contour",cont=c(10,50,75,95),drawlabels=F,col=c('white','gray90','gray70','gray45','gray25')) #points(df.predct)
points(theta1,theta2,pch=16,col='red')

H1.cg <- Hpi(x=df.predcg)
H2.cg <- Hpi.diag(x=df.predcg)

fhat1.cg <- kde(x=df.predcg, H=H1.cg)
fhat2.cg <- kde(x=df.predcg, H=H2.cg)

plot(fhat1.cg,xlab = 'LFA',ylab = 'RFA',main='CopG, H no diagonal',xlim = c(0.1,0.9),ylim = c(0.15,0.85),lwd=3,display="filled.contour",cont=c(10,50,75,95),drawlabels=F,col=c('white','gray90','gray70','gray45','gray25')) #points(df.predcg)
points(theta1,theta2,pch=16,col='red')

plot(fhat2.cg,xlab = 'LFA',ylab = 'RFA',main='CopG, H diagonal',xlim = c(0.15,0.8),ylim = c(0.15,0.8),lwd=3,display="filled.contour",cont=c(10,50,75,95),drawlabels=F,col=c('white','gray90','gray70','gray45','gray25')) #points(df.predcg)
points(theta1,theta2,pch=16,col='red')

