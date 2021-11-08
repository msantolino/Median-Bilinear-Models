################
## packages
################

library(gnm)
library(quantreg)
library(fastDummies)

######################################
### Calibration of bilinear models
######################################

################################
# Mortality data Spain 1908-2016
################################

exposure<-read.table("Exposures.txt",header=T,quote="")
colnames(exposure)<-c('Year','Age','ExpFem','ExpMal','ExpTot')
mortality<-read.table("Mx.txt",header=T,quote="")
colnames(mortality)<-c('Year','Age','MortFem','MortMal','MortTot')

dataspain<-data.frame(exposure,mortality[,3:5])

# Males 0-100 years
bdspainM<-subset(dataspain,select=c(Year,Age,ExpMal,MortMal))
bdspain<-bdspainM[complete.cases(bdspainM),]
bdspainRH<-data.frame(bdspain$MortMal,bdspain$Age,bdspain$Year)
colnames(bdspainRH)<-c("rate","age","year")
bdspainR<-subset(bdspainRH, age<101)


#binary variables
bdspainR$age<-as.factor(bdspainR$age)
bdspainR$year<-as.factor(bdspainR$year)
bddich<-dummy_cols(bdspainR)

varD<-cbind(log(bddich$rate))#dependent variable
varXage<-as.matrix(bddich[,4:104])#explanatory variables
nyear<-105+length(table(bdspainR$year))-1
varXyear<-as.matrix(bddich[,105:nyear])

##########################
# Mean SVD model
##########################
nag<-101
nyear<-109

reg_insv <- lm(log(rate) ~-1+ factor(age), data = bdspainR)

yest<-reg_insv$fitted.values
yob<-log(bdspainR$rate)
error<-yob-yest

A<-matrix(c(error), nrow=nag, ncol=nyear)
B<-svd(A)

axM<-reg_insv$coef[1:nag]
bxM<-B$u[,1]
gtM<-(B$v[,1])*B$d[1]

axMsv<-axM+mean(gtM)
bxMsv<-bxM/sum(bxM) #restriction: sum=1
gtMsv<-(gtM-mean(gtM))*sum(bxM) #restriction: sum=0


##########################
# Mean MV model
##########################

reg_init <- gnm(log(rate) ~ age + Mult(age, year), family = "gaussian", data = bdspainR)
betaLC<-coef(reg_init)

axM<-betaLC[1:101]
axM[2:101]<-betaLC[2:101]+axM[1]
bxM<-betaLC[102:202]
gtM<-betaLC[203:length(betaLC)]

axMmv<-axM+mean(gtM)*bxM
bxMmv<-bxM/sum(bxM) #restriction: sum=1
gtMmv<-(gtM-mean(gtM))*sum(bxM) #restriction: sum=0


##########################
# Med A-AS model
##########################


t1<-coefLC<-c(axMmv,bxMmv,gtMmv) #initial values
QLC50<-nlrqB(LCmodel,t1,theta=0.5)

axMEa<-QLC50$coef[1:nag]+mean(QLC50$coef[(nag+nag+1):length(t1)])*QLC50$coef[(1+nag):(2*nag)]
bxMEa<-QLC50$coef[(1+nag):(2*nag)]/sum(QLC50$coef[(1+nag):(2*nag)]) #restriction: sum=1
gtMEa<-(QLC50$coef[(nag+nag+1):length(QLC50$coef)]-mean(QLC50$coef[(nag+nag+1):length(QLC50$coef)]))*sum(QLC50$coef[(1+nag):(2*nag)]) #restriction: sum=0


##########################
#  Median B-MV model
##########################

# estimate a_x
tau<-0.5
coefA<-EM.qrg(varD,varXage,tau)$theta[1:nag]
a_est<-cbind(rep(coefA,nyear))
yadj<-varD-a_est

# estimate gamma_y
coefY<-EM.qrg(yadj,varXyear,tau)$theta[1:nyear]
A<-matrix(rep(coefY,each=nag*nag),ncol=nag,byrow=T)
XAd<-varXage*A  

# estimate b_x
coefAY<-EM.qrg(yadj,XAd,tau)$theta[1:nag]

# fitted y
a_est<-cbind(rep(coefA,nyear))
b_est<-cbind(rep(coefAY,nyear))
g_est<-cbind(rep(coefY,each=nag))
ypred<-a_est+b_est*g_est

# objective function
error<-varD-ypred
objec<-ifelse(error>=0, error*tau,error*(tau-1))
minobj<-sum(objec)
epsilon<-minobj

# loop
max_it<-100
iter <- 0;
while (epsilon>0.001 & iter < max_it){
      print(iter) 
	iter <- iter + 1

# estimate a_x
varXageL<-varXage*matrix(rep(coefA,nyear*nag),ncol=nag,byrow=F)
b_est<-cbind(rep(coefAY,nyear))
g_est<-cbind(rep(coefY,each=nag))
yL<-varD-b_est*g_est
coefAL<-EM.qrg(yL,varXageL,tau,error=0.1)$theta[1:nag]
a_estL<-cbind(rep(coefA*coefAL,nyear))

yadjL<-varD-a_estL

# estimate gamma_y
varXyearL<-varXyear*matrix(rep(coefY,each=nyear*nag),ncol=nyear,byrow=T)
varXyearL2<-varXyear*matrix(rep(coefAY,nyear*nyear),ncol=nyear,byrow=F)
coefYL<-EM.qrg(yadjL,varXyearL*varXyearL2,tau)$theta[1:nyear]
g_estL<-cbind(rep(coefY*coefYL,each=nag))

# estimate beta_x
XAdL<-varXage*matrix(rep(coefAY,nyear*nag),ncol=nag,byrow=F)
AL<-matrix(rep(coefYL*coefY,each=nag*nag),ncol=nag,byrow=T)
coefAYL<-EM.qrg(yadjL,XAdL*AL,tau)$theta[1:nag]
b_estL<-cbind(rep(coefAYL*coefAY,nyear))

ypredL<-a_estL+b_estL*g_estL

errorL<-varD-ypredL
objecL<-ifelse(errorL>=0, errorL*tau,errorL*(tau-1))
minobjL<-sum(objecL)
epsilon<-minobj-minobjL
print(epsilon)
parameters<-c(coefA,coefAY,coefY)

#step
coefA<-coefA*coefAL
coefY<-coefYL*coefY
coefAY<-coefAY*coefAYL
minobj<-minobjL
}

axMEmv<-coefA+mean(coefY)*coefAY
bxMEmv<-coefAY/sum(coefAY) #restriction: sum=1
gtMEmv<-(coefY-mean(coefY))*sum(coefAY) #restriction: sum=0



##########################
#  Median B-MP model
##########################

# estimate a_x
coefA<-rq.fit.fnb(varXage,varD,tau)$coeff
n<-nrow(varXage)/length(coefA)
a_est<-cbind(rep(coefA,n))
yadj<-varD-a_est

# estimate gamma_y
coefY<-rq.fit.fnb(varXyear,yadj,tau)$coeff
n2<-nrow(varXyear)/length(coefY)
A<-matrix(rep(coefY,each=n2*n2),ncol=n2,byrow=T)
XAd<-varXage*A  

# estimate b_x
coefAY<-rq.fit.fnb(XAd,yadj,tau)$coeff

# fitted y
a_est<-cbind(rep(coefA,n))
b_est<-cbind(rep(coefAY,n))
g_est<-cbind(rep(coefY,each=n2))
ypred<-a_est+b_est*g_est

# objective function
error<-varD-ypred
objec<-ifelse(error>=0, error*tau,error*(tau-1))
minobj<-sum(objec)
epsilon<-minobj

# loop
max_it<-100
iter <- 0;
while (epsilon>0.01 & iter < max_it){
      print(iter) 
	iter <- iter + 1

# estimate a_x
varXageL<-varXage*matrix(rep(coefA,n*n2),ncol=n2,byrow=F)
b_est<-cbind(rep(coefAY,n))
g_est<-cbind(rep(coefY,each=n2))
yL<-varD-b_est*g_est
coefAL<-rq.fit.fnb(varXageL,yL,tau)$coeff
a_estL<-cbind(rep(coefA*coefAL,n))

yadjL<-varD-a_estL

# estimate gamma_y
varXyearL<-varXyear*matrix(rep(coefY,each=n*n2),ncol=n,byrow=T)
varXyearL2<-varXyear*matrix(rep(coefAY,n*n),ncol=n,byrow=F)
coefYL<-rq.fit.fnb(varXyearL*varXyearL2,yadjL,tau)$coeff
g_estL<-cbind(rep(coefY*coefYL,each=n2))

# estimate beta_x
XAdL<-varXage*matrix(rep(coefAY,n*n2),ncol=n2,byrow=F)
AL<-matrix(rep(coefYL*coefY,each=n2*n2),ncol=n2,byrow=T)
coefAYL<-rq.fit.fnb(XAdL*AL,yadjL,tau)$coeff
b_estL<-cbind(rep(coefAYL*coefAY,n))

ypredL<-a_estL+b_estL*g_estL

errorL<-varD-ypredL
objecL<-ifelse(errorL>=0, errorL*tau,errorL*(tau-1))
minobjL<-sum(objecL)
epsilon<-minobj-minobjL
print(epsilon)
parameters<-c(coefA,coefAY,coefY)

#step
coefA<-coefA*coefAL
coefY<-coefYL*coefY
coefAY<-coefAY*coefAYL
minobj<-minobjL
}


axMEb<-coefA+mean(coefY)*coefAY
bxMEb<-coefAY/sum(coefAY) #restriction: sum=1
gtMEb<-(coefY-mean(coefY))*sum(coefAY) #restriction: sum=0







