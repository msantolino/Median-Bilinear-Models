################################
########### Auxiliary functions
################################

#####################################
# Lee-Carter model specification
#####################################
LCmodel<-function(t,nag=101,ntiem=109)
{
ypred<-cbind(rep(t[1:nag],ntiem))+cbind(rep(t[(1+nag):(2*nag)],ntiem))*cbind(rep(t[(nag+nag+1):length(t)],each=nag))
A<-varD-ypred
#gradiente
gradA<- varXage
gradY<- varXyear*matrix(rep(t[(1+nag):(2*nag)],ntiem*ntiem),ncol=ntiem,byrow=F)
gradAY<- varXage*matrix(rep(t[(2*nag+1):length(t)],each=nag*nag),ncol=nag,byrow=T)
B<-cbind(gradA,gradAY,gradY)
C=list(f=A,J=B)
return(C)
}


######################################################################################
# Functions nlrqB, rho.rq, mek.rq and model.step.rq are derived from
# Koenker, R. (2020). Non linear quantile regression.
# http://www.econ.uiuc.edu/~roger/research/nlrq/nlrq.html, Accessed 16 February, 2021.
# Functions are adapted to deal with non-full-rank design matrices.
######################################################################################

nlrqB<-
function(model, t, k = 3, theta = .5, big=1e+20, nit.max = 100,
eps = 1e-07, beta = 0.97)
{
#function to compute nonlinear rq estimate
# t is the initial value of the unknown parameter
# model is a user-provided function which returns components
# f=(f_i (x_i , t)
# J=(grad f_i )
# theta is the desired quantile
# k is the number of Meketon steps per iteration
# eps and eta are tolerance parameters
#
#function returns
# coef is the value of the parameter at the solution
# obj is the value of the objective function at the solution
# nit is the number of "Meketon steps" taken
m <- model(t)
n <- length(m$f)
w <- rep(0, n)
snew <- sum(rho.rq(m$f,theta))
sold <- big
nit <- 0
while(sold - snew > eps & nit < nit.max) {
z <- mek.rq(m$J, m$f, k, w, theta=theta, int = F,
eps = eps, beta = beta)
step <- z$coef
t0 <- t
l <- nlminb(start = 1, objective = model.step.rq, lower = 0,
upper = 1, model = model, t0 = t, theta=theta,
step = step)$par
t <- t0 + l * step
m <- model(t)
sold <- snew
snew <- sum(rho.rq(m$f,theta))
w <- lsfit(m$J, z$w, int = F)$resid
w1 <- max(pmax(w,0))
if(w1>theta)
w <- w*theta/(w1 + eps)
w0 <- max(pmax(-w,0))
if(w0>1-theta)
w <- w*(1-theta)/(w0 + eps)
#print(c(t, l, sum(rho.rq(m$f,theta))))
nit <- nit+1
}
return(list(coef=t,obj=snew,nit=nit))
}

rho.rq<-
function(u,theta){theta*pmax(u,0)+(theta - 1)*pmin(u,0)}

mek.rq<-
function(x, y, kmax = 1000, w, theta=.5, int = T, big=1e+20,
eps = 1e-06, beta = 0.97)
{
if(int == T)
x <- cbind(1, x)
yw <- big
k <- 1
while(k <= kmax & yw - crossprod(y, w) > eps) {
d <- pmin(theta - w, 1 - theta + w)
z <- lsfit(x, y, d^2, int = F)
yw <- sum(rho.rq(z$resid,theta))
k <- k + 1
s <- z$resid * d^2
alpha <- max(eps, pmax(s/(theta - w), -s/(1 - theta + w)))
w <- w + (beta/alpha) * s
}
coef <- z$coef
return(list(coef=coef,w=w))
}

model.step.rq<-
function(lambda, t0, step, model, theta)
{
sum(rho.rq(model(t0 + lambda * step)$f, theta))
}



######################################################################################
# Function EM.qrg is derived from the R package ALDqr.
# ALDqr: Quantile Regression Using Asymmetric Laplace Distribution
# Luis Benites Sanchez, Christian E. Galarza, Victor H. Lachos
# The function EM.qr of this package is modified to deal with sparse matrices.
######################################################################################

EM.qrg<-function (y, x = NULL, tau = NULL, error = 0.01, iter = 2000, 
    envelope = FALSE) 
{
    logVQR <- function(y, x, tau, theta) {
        p <- ncol(x)
        n <- nrow(x)
        beta <- theta[1:p]
        sigma <- theta[p + 1]
        mu <- x %*% beta
        muc <- (y - mu)/sigma
        Ind <- (muc < 0) + 0
        logver <- sum(-log(sigma) + log(tau * (1 - tau)) - muc * 
            (tau - Ind))
        return(logver)
    }
    p <- ncol(x)
    n <- nrow(x)
    reg <- lm(y ~ x[, 2:p])
    taup2 <- (2/(tau * (1 - tau)))
    thep <- (1 - 2 * tau)/(tau * (1 - tau))
    beta <- as.vector(coefficients(reg), mode = "numeric")
    sigma <- sqrt(sum((y - x %*% beta)^2)/(n - p))
    lk = lk1 = lk2 <- logVQR(y, x, tau, c(beta, sigma))
    teta_velho <- matrix(c(beta, sigma), ncol = 1)
    cont <- 0
    criterio <- 1
    while (criterio > error) {
        print(criterio)
        cont <- (cont + 1)
        muc <- (y - x %*% beta)
        delta2 <- (y - x %*% beta)^2/(taup2 * sigma)
        gamma2 <- (2 + thep^2/taup2)/sigma
        vchpN <- besselK(sqrt(delta2 * gamma2), 0.5 - 1)/(besselK(sqrt(delta2 * 
            gamma2), 0.5)) * (sqrt(delta2/gamma2))^(-1)
        vchp1 <- besselK(sqrt(delta2 * gamma2), 0.5 + 1)/(besselK(sqrt(delta2 * 
            gamma2), 0.5)) * (sqrt(delta2/gamma2))
        sigma <- sum(vchpN * muc^2 - 2 * muc * thep + vchp1 * 
            (thep^2 + 2 * taup2), na.rm=TRUE)/(3 * n * taup2)
        ym<-y-(thep/vchpN)
        beta<-cbind(slm(ym~-1+x, weights=c(vchpN))$coef);
	  teta_novo <- matrix(c(beta, sigma), ncol = 1)
        criterio <- sqrt(sum((teta_velho - teta_novo)^2))
        lk3 <- logVQR(y, x, tau, c(beta, sigma))
        if (cont < 2) 
            criterio <- abs(lk2 - lk3)/abs(lk3)
        else {
            tmp <- (lk3 - lk2)/(lk2 - lk1)
            tmp2 <- lk2 + (lk3 - lk2)/(1 - tmp)
            criterio <- abs(tmp2 - lk3)
        }
        lk2 <- lk3
        if (cont == iter) {
            break
        }
        teta_velho <- teta_novo
    }
    Weights <- vchpN * vchp1
    logver <- logVQR(y, x, tau, teta_novo)
    return(list(theta = teta_novo,logver = logver, 
        iter = cont, Weights = Weights, di = abs(muc)/sigma))
}

