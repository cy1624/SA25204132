## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
set.seed(1)
require(magrittr)
require(Rcpp)
require(microbenchmark)
require(stats)
require(rbenchmark)
require(RcppArmadillo)
require(snowfall)
require(boot)
require(lpSolve)
require(DAAG)
require(bootstrap)
require(gtools)

## -----------------------------------------------------------------------------
set.seed(0)
rRayleigh <- function(n, sigma){
  y <- runif(n)
  x <- sqrt(-2*sigma^2*log(1-y))
  x
}

## ----echo=FALSE---------------------------------------------------------------
sigmas <- c(.5,1,2,3)

for (sigma in sigmas) {
  x <- rRayleigh(100, sigma = sigma)
  hist(x, main = paste0("sigma=", sigma))
  abline(v=sigma, col="red")
}

## -----------------------------------------------------------------------------
n <- 1000
y <- runif(n)
p <- c(.1,.2,.2,.2,.3)
cp <- cumsum(p)
x <- findInterval(y, cp)

## -----------------------------------------------------------------------------
phat <- as.vector(table(x)/n)
data.frame(phat,p)


## -----------------------------------------------------------------------------
rbeta2 <- function(n, a, b){
  if(a>1 && b>1){
    c <- (a-1)^(a-1)*(b-1)^(b-1)/beta(a,b)/(a+b-2)^(a+b-2)
  } else{
    c <- 1
  }
  y <- runif(n*(c+0.1))
  r <- runif(length(y))
  x <- y[r<dbeta(y,a,b)/c]
  if(length(x)>=n){
    return(x[1:n])
  } else{
    return(c(x, rbeta2(n-length(x), a, b)))
  }
}

## -----------------------------------------------------------------------------
x <- rbeta2(1000,3,2)

## -----------------------------------------------------------------------------
x_seq <- seq(0,1,0.02)
y_seq <- sapply(x_seq, function(x) dbeta(x,3,2))


hist(x, freq = F, ylim = c(0,1.8),
     main = "Beta(3,2)直方图和理论密度函数（红色）")
lines(x_seq,y_seq,col="red")

## -----------------------------------------------------------------------------
rmixnorm <- function(n,p){
  x1 <- rnorm(n,0,1)
  x2 <- rnorm(n,3,1)
  r <- sample(0:1,n, replace = T, prob = c(1-p,p))
  return(r*x1+(1-r)*x2)
}

## -----------------------------------------------------------------------------
x_seq <- seq(-3,6,0.02)
y_seq <- sapply(x_seq, 
                function(x) (dnorm(x,0,1)+dnorm(x,3,1))/2)

n <- 1000
x <- rmixnorm(n,p=0.5)
hist(x, freq = F,
     main = "p=0.5时的直方图和理论密度函数（红色）")
lines(x_seq,y_seq,col="red")

## ----echo=FALSE---------------------------------------------------------------

for(p in seq(0.1,0.9,0.1)){
  x_seq <- seq(-3,6,0.02)
  y_seq <- sapply(x_seq, 
                function(x) p*dnorm(x,0,1)+(1-p)*dnorm(x,3,1))
  x <- rmixnorm(n,p=p)
  hist(x, freq = F, ylim = c(0,max(y_seq)+0.02),
      main = paste("p=",p))
  lines(x_seq,y_seq,col="red")
}

## -----------------------------------------------------------------------------
n <- 1000
Gam <- rgamma(n, 4, 2)
y <- rexp(n, Gam)

## -----------------------------------------------------------------------------
set.seed(0)
pbeta1 <- function(n,a,b,x){
  r <- rbeta(n,a,b)
  phat <- sapply(x, function(k) sum(r<k)/n)
  return(phat)
}

## -----------------------------------------------------------------------------
p1 <- pbeta1(10000,3,3,(1:9)/10) 
p2 <- pbeta((1:9)/10,3,3) 
p1 %>% round(.,4)# 蒙特卡洛的结果
p2 %>% round(.,4)# pbeta的结果
round(abs(p1-p2)/p2,4) # 蒙特卡洛估计的相对偏差

## -----------------------------------------------------------------------------
n <- 1000

hat1 <- replicate(100, expr = {
  x <- runif(n)
  y <- sqrt(1-2*log(x))
  mean(y/sqrt(2*pi*exp(1)))
})

hat2 <- replicate(100, expr = {
    x <- rnorm(n*10)
    x <- x[x>=1][1:n]
    mean(x^2*(1-pnorm(1)))
})

c(mean(hat1),mean(hat2))
c(sd(hat1),sd(hat2))

## -----------------------------------------------------------------------------
x <- seq(1,10,0.02)
y1 <- sapply(x, function(x){
  x^2/sqrt(2*pi)*exp(-x^2/2)
})
y2 <- sapply(x, function(x){
  x*exp(-x^2/2+1/2)
})
y3 <- sapply(x, function(x){
  1/sqrt(2*pi)*exp(-x^2/2)/(1-pnorm(1))
})

plot(x,y1,ylim = c(0,1.6), type = "l",ylab="")
lines(x,y2,col="red")
lines(x,y3,col="lightgreen")
legend("topright", legend = c("g(x)","f1(x)","f2(x)"),
       col = c("black","red","lightgreen"),lty = rep(1,3))

## -----------------------------------------------------------------------------
x <- c(1,2,4,6,8)*1e4
t <- benchmark(replications = 1000,
          sort(sample(1:x[1],x[1])),
          sort(sample(1:x[2],x[2])),
          sort(sample(1:x[3],x[3])),
          sort(sample(1:x[4],x[4])),
          sort(sample(1:x[5],x[5])),
          columns = c("user.self"))
tn <- t$user.self

## -----------------------------------------------------------------------------
an <- x*log(x)
fit <- lm(tn~an)
summary(fit)
plot(an,tn)
lines(an,fit$fitted.values,col="red")

## ----warning=F----------------------------------------------------------------
set.seed(0)
m <- 10000
k <- 5
g <- function(x) {
exp(-x- log(1+x^2)) * (x > 0) * (x < 1)
}
F_inv <- function(x){ # 逆变换法
  -log(1-(1-exp(-1))*x)
}
rf <- function(n,a0,a1){
    u <- runif(n)
    x <- -log(exp(-a0)-(exp(-a0)-exp(-a1))*u)
    x
}
fg <- function(x){ # g/f
  g(x)*exp(x)*(1-exp(-1))
}

j <- seq(0,1,length.out=k+1) %>% F_inv
test <- function(){mean(sapply(1:k,FUN = function(i){
    x <- rf(m,j[i],j[i+1])
    mean(fg(x))
}))}



sfInit(parallel = T,cpus = 8)
sfExport("j","k","m","g","fg","F_inv","rf","test")
result <- sfSapply(1:1000,fun = function(k) test())
sfStop()

mean(result)
sd(result)

## -----------------------------------------------------------------------------
n <- 1:5 *10
m <- 1000
mu0 <- 500
mu <- c(seq(450, 650, 10))
sigma <- 100

powers <- sapply(n, FUN = function(n){
  sapply(mu, FUN = function(mu){
    pvalues <- replicate(m, expr = {
    x <- rnorm(n, mean = mu, sd = sigma)
    ttest <- t.test(x,
    alternative = "greater", mu = mu0)
    ttest$p.value } )
    mean(pvalues <= .05)
  })
})

## -----------------------------------------------------------------------------
plot(mu,powers[,1], type="l",col=rainbow(5)[1],ylab = "power")
abline(v=500,h=0.05)
for(i in 2:5){
  lines(mu,powers[,i],col=rainbow(5)[i])
}
legend("bottomright",
       legend = c("n=10","n=20","n=30","n=40","n=50"),
       col = rainbow(5),lty=rep(1,5))

## -----------------------------------------------------------------------------
cover1 <- function(x,mu0,alpha=0.05){# 判断均值是否落在置信区间内
  n <- length(x)
  return((mu0>mean(x)-qt(1-alpha/2,n-1)*sd(x)/sqrt(n))*
           (mu0<mean(x)+qt(1-alpha/2,n-1)*sd(x)/sqrt(n)))
}
cover2 <- function(x,sigma2,alpha=0.05){# 判断方差是否落在置信区间内
  n <- length(x)
  S <- (n-1) * var(x) / qchisq(alpha, df = n-1)
  return(sigma2<S)
}

test <- replicate(50000,expr = {
  x <- rchisq(20,2)
  c(cover1(x,2),cover2(x,4))
})

rowMeans(test)


## -----------------------------------------------------------------------------
n <- 20
m <- 10000
fs <- c(function(.) rchisq(.,df=1),
        function(.) runif(.,min=0,max=2),
        function(.) rexp(.,rate=1))
mu0 <- 1
result <- sapply(fs, function(f){
    replicate(m, expr = {
      x <- f(n)
      cover1(x,mu0)
    }) %>% mean
})

result


## ----echo=FALSE---------------------------------------------------------------
x <- seq(0,5,.05)
y1 <- dnorm(x,1)
y2 <- dchisq(x,1)
y3 <- dunif(x,0,2)
y4 <- dexp(x,1)
plot(x,y1,type="l",ylim=c(0,2),ylab = "density")
lines(x,y2,col="red")
lines(x,y3,col="green")
lines(x,y4,col="blue")
legend("topright",legend = c("N(0,1)",
  expression({chi^2}(1)),"U(0,2)","Exp(1)"),
  col = c("black","red","green","blue"),
  lty = rep(1,5))


## -----------------------------------------------------------------------------
set.seed(0)
test <- function(r, r_true){
  t <- table(r, r_true)
  return(list(FWER=t[1,1]>0,
              FDR=t[1,1]/sum(t[1,]),
              TPR=t[1,2]/sum(t[,2])))
}

## -----------------------------------------------------------------------------
n0 <- 950
n1 <- 50
alpha <- 0.1
m <- 10000
H <- c(rep(0,n0),rep(1,n1))

result <- replicate(m, expr = {
  p <- c(runif(n0), rbeta(n1,0.1,1))
  p_Bonf <- p.adjust(p, method = "bonf")
  r_Bonf <- p_Bonf > alpha
  p_BH <- p.adjust(p, method = "BH")
  r_BH <- p_BH > alpha
  c(unlist(test(r_Bonf,H)), unlist(test(r_BH,H)))
})

matrix(rowMeans(result), nrow = 3,
       dimnames = list(c('FWER', 'FDR', 'TPR'),
       c('Bonferroni correction', 'B-H correction')))


## -----------------------------------------------------------------------------
data <- boot::aircondit %>% unlist
lambda_hat <- 1/mean(data)
lambda_hat

## -----------------------------------------------------------------------------
lambda_boot <- replicate(10000, expr = {
  1/mean(sample(data, replace = T))
})
mean(lambda_boot)-lambda_hat # 偏差
sd(lambda_boot) # 标准差


## -----------------------------------------------------------------------------
data <- bootstrap::scor
n <- dim(data)[1]
Cov <- cov(data)
Eig <- eigen(Cov)$values
theta <- Eig[1]/sum(Eig)

## -----------------------------------------------------------------------------
theta_boot <- replicate(10000, expr = {
  i <- sample(1:n, replace = T)
  Cov <- cov(data[i,])
  Eig <- eigen(Cov)$values
  Eig[1]/sum(Eig)
})

mean(theta_boot)-theta # 偏差
sd(theta_boot) # 标准差

## -----------------------------------------------------------------------------
set.seed(0)
estimate <- function(data){
  n <- dim(data)[1]
  Cov <- cov(data)
  Eig <- eigen(Cov)$values
  Eig[1]/sum(Eig)
}

## -----------------------------------------------------------------------------
data <- bootstrap::scor
n <- dim(data)[1]


theta <- estimate(data)
theta_jack <- sapply(1:n, function(i) estimate(data[-i,]))

## -----------------------------------------------------------------------------
data <- DAAG::ironslag
n <- dim(data)[1]
index <- gtools::permutations(n,2)

res <- apply(index, MARGIN = 1, FUN = function(i){
    fit1 <- lm(magnetic~chemical, data[-i,])
    fit2 <- lm(magnetic~chemical+I(chemical^2), data[-i,])
    fit3 <- lm(log(magnetic)~chemical, data[-i,])
    fit4 <- lm(log(magnetic)~log(chemical), data[-i,])
    return(c(predict(fit1, newdata=data[i,]) - data[i,'magnetic'],
           predict(fit2, newdata=data[i,]) - data[i,'magnetic'],
           exp(predict(fit3, newdata=data[i,])) - data[i,'magnetic'],
           exp(predict(fit4, newdata=data[i,])) - data[i,'magnetic']))
  })

## -----------------------------------------------------------------------------
in_interval <- function(theta, interval){
  if(theta<interval[1]){
    return(-1)
  } else if(theta>interval[2]){
    return(1)
  } else return(0)
}

## -----------------------------------------------------------------------------
# 参数设置
R <- 200
m <- 1000
mu <- 0

## -----------------------------------------------------------------------------
n <- 20
result1 <- replicate(m, expr = {
  x <- rnorm(n,mu,1)
  x_boot <- boot(x, statistic = function(x,i) mean(x[i]), R=R)
  interval <- boot.ci(x_boot, type = c("norm",'basic','perc'))
  return(c(in_interval(mu,interval$normal[2:3]),
           in_interval(mu,interval$basic[4:5]),
           in_interval(mu,interval$percent[4:5])))
})
n <- 50
result2 <- replicate(m, expr = {
  x <- rnorm(n,mu,1)
  x_boot <- boot(x, statistic = function(x,i) mean(x[i]), R=R)
  interval <- boot.ci(x_boot, type = c("norm",'basic','perc'))
  return(c(in_interval(mu,interval$normal[2:3]),
           in_interval(mu,interval$basic[4:5]),
           in_interval(mu,interval$percent[4:5])))
})
n <- 100
result3 <- replicate(m, expr = {
  x <- rnorm(n,mu,1)
  x_boot <- boot(x, statistic = function(x,i) mean(x[i]), R=R)
  interval <- boot.ci(x_boot, type = c("norm",'basic','perc'))
  return(c(in_interval(mu,interval$normal[2:3]),
           in_interval(mu,interval$basic[4:5]),
           in_interval(mu,interval$percent[4:5])))
})
n <- 200
result4 <- replicate(m, expr = {
  x <- rnorm(n,mu,1)
  x_boot <- boot(x, statistic = function(x,i) mean(x[i]), R=R)
  interval <- boot.ci(x_boot, type = c("norm",'basic','perc'))
  return(c(in_interval(mu,interval$normal[2:3]),
           in_interval(mu,interval$basic[4:5]),
           in_interval(mu,interval$percent[4:5])))
})

# n=20的结果（-1代表落在区间左侧，0代表落在区间内，1代表落在区间右侧）
(apply(result1,1, table)/m) %>% `colnames<-`(.,c("norm",'basic','perc'))
# n=50的结果
(apply(result2,1, table)/m) %>% `colnames<-`(.,c("norm",'basic','perc'))
# n=100的结果
(apply(result3,1, table)/m) %>% `colnames<-`(.,c("norm",'basic','perc'))
# n=200的结果
(apply(result4,1, table)/m) %>% `colnames<-`(.,c("norm",'basic','perc'))

## -----------------------------------------------------------------------------
# 参数设置
lambda <- 2
R <- 1000
m <- 1000

## -----------------------------------------------------------------------------
n <- 5
bias_n <- lambda/(n-1)
sd_n <- lambda*n/(n-1)/sqrt(n-2)
result1 <- replicate(m, expr = {
  x <- rexp(n, rate = lambda)
  x_boot <- boot(x,statistic = function(x,i) 1/mean(x[i]),R=R,parallel = 'multicore')$t
  return(c(mean(x_boot)-1/mean(x), sd(x_boot)))
})
# n=5
data.frame(true=c(bias_n,sd_n),bootstrap=rowMeans(result1),
           row.names = c('bias','sd'))

## -----------------------------------------------------------------------------
n <- 10
bias_n <- lambda/(n-1)
sd_n <- lambda*n/(n-1)/sqrt(n-2)
result2 <- replicate(m, expr = {
  x <- rexp(n, rate = lambda)
  x_boot <- boot(x,statistic = function(x,i) 1/mean(x[i]),R=R,parallel = 'multicore')$t
  return(c(mean(x_boot)-1/mean(x), sd(x_boot)))
})
# n=10
data.frame(true=c(bias_n,sd_n),bootstrap=rowMeans(result2),
           row.names = c('bias','sd'))

## -----------------------------------------------------------------------------
n <- 20
bias_n <- lambda/(n-1)
sd_n <- lambda*n/(n-1)/sqrt(n-2)
result3 <- replicate(m, expr = {
  x <- rexp(n, rate = lambda)
  x_boot <- boot(x,statistic = function(x,i) 1/mean(x[i]),R=R,parallel = 'multicore')$t
  return(c(mean(x_boot)-1/mean(x), sd(x_boot)))
})
# n=20
data.frame(true=c(bias_n,sd_n),bootstrap=rowMeans(result3),
           row.names = c('bias','sd'))

## -----------------------------------------------------------------------------
set.seed(0)
CM <- function(x,y){
  x <- sort(x)
  y <- sort(y)
  n <- length(x)
  m <- length(y)
  r <- rank(c(x,y))
  rx <- r[1:n]
  ry <- r[-(1:n)]
  U <- n*sum((rx-1:n)^2)+m*sum((ry-1:m)^2)
  W <- U/n/m/(n+m)-(4*m*n-1)/6/(n+m)
  return(W)
}

CM.boot <- function(z,ix,nx){
  z <- z[ix]
  CM(z[1:nx],z[-(1:nx)])
}

## -----------------------------------------------------------------------------
attach(chickwts)
x <- sort(weight[feed == "soybean"])
y <- sort(weight[feed == "linseed"])
detach(chickwts)

boot.obj <- boot(data = c(x,y), statistic = CM.boot, R=999,
                 sim = "permutation", nx=length(x))
ts <- c(boot.obj$t0,boot.obj$t)
p.value <- mean(ts>=ts[1])
hist(ts)
abline(v=ts[1], col='red')

## -----------------------------------------------------------------------------
K <- function(x, y) {
 X <- x- mean(x)
 Y <- y- mean(y)
 outx <- sum(X > max(Y)) + sum(X < min(Y))
 outy <- sum(Y > max(X)) + sum(Y < min(X))
 return(as.integer(max(c(outx, outy))))
 }

K.boot <- function(z, ix, nx){
  z <- z[ix]
  return(K(z[1:nx],z[-(1:nx)]))
}

K.test <- function(x,y,...){
  boot.obj <- boot(data = c(x,y), statistic = K.boot,...,
                 sim = "permutation", nx=length(x))
  ts <- c(boot.obj$t0,boot.obj$t)
  mean(ts>=ts[1])
} 

## ----warning=FALSE------------------------------------------------------------
n1 <- 20
n2 <- 30
mu1 <- mu2 <- 0
sigma1 <- sigma2 <- 1
m <- 10000

sfInit(parallel = T,cpus = 8)
sfExportAll()
sfLibrary(boot)
ps1 <- sfSapply(1:m, fun = function(k){
  x <- rnorm(n1, mu1, sigma1)
  y <- rnorm(n2, mu2, sigma2)
  x <- x- mean(x) #centered by sample mean
  y <- y- mean(y)
  K.test(x, y, R=999)
})
sfStop()

n1 <- 20
n2 <- 50
mu1 <- mu2 <- 0
sigma1 <- sigma2 <- 1
m <- 10000

sfInit(parallel = T,cpus = 8)
sfExportAll()
sfLibrary(boot)
ps2 <- sfSapply(1:m, fun = function(k){
  x <- rnorm(n1, mu1, sigma1)
  y <- rnorm(n2, mu2, sigma2)
  x <- x- mean(x) #centered by sample mean
  y <- y- mean(y)
  K.test(x, y, R=999)
})
sfStop()

## -----------------------------------------------------------------------------
set.seed(0)
RandomWalkMCMC <- function(N, sigma2, X1=NULL){
  x <- numeric(N)
  r <- numeric(N)
  sigma <- sqrt(sigma2)
  if(is.null(X1)){
    x[1] <- rnorm(1,0,sigma)
  } else {
    x[1] <- X1
  }
  
  r[1] <- T
  for(i in 2:N){
    y <- rnorm(1, x[i-1], sigma)
    u <- runif(1)
    #a <- exp(abs(x[i-1])-abs(y))
    a <- dnorm(y,0,1)/dnorm(x[i-1],0,1)
    r[i] <- u<a
    x[i] <- r[i]*y + (1-r[i])*x[i-1]
  }
  return(list(x=x, acceptance=r))
}

## -----------------------------------------------------------------------------
N <- 5000
rst1 <- RandomWalkMCMC(N=N, sigma2 = 0.1)
rst2 <- RandomWalkMCMC(N=N, sigma2 = 0.5)
rst3 <- RandomWalkMCMC(N=N, sigma2 = 2)
rst4 <- RandomWalkMCMC(N=N, sigma2 = 16)


plot(rst1$x, main = 'sigma^2 = 0.1',type = 'l',ylab = 'x')
plot(rst2$x, main = 'sigma^2 = 0.5',type = 'l',ylab = 'x')
plot(rst3$x, main = 'sigma^2 = 2',type = 'l',ylab = 'x')
plot(rst4$x, main = 'sigma^2 = 16',type = 'l',ylab = 'x')

# 接受率
c(mean(rst1$acceptance),mean(rst2$acceptance),
  mean(rst3$acceptance),mean(rst4$acceptance))

## -----------------------------------------------------------------------------
N <- 3000
a <- 1
b <- 2
n <- 15
x <- numeric(N)
y <- numeric(N)
x[1] <- 1

for(i in 1:N){
  y[i] <- rbeta(1,x[i]+a,n-x[i]+b)
  x[i+1] <- rbinom(1,n,y[i])
}
x <- x[1:N]

## -----------------------------------------------------------------------------
Gelman.Rubin <- function(psi) {
  # psi[i,j] is the statistic psi(X[i,1:j])
  # for chain in i-th row of X
  psi <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)
  psi.means <- rowMeans(psi)
  B <- n * var(psi.means)
  #row means
  #between variance est.
  psi.w <- apply(psi, 1, "var") #within variances
  W <- mean(psi.w)
  #within est.
  v.hat <- W*(n-1)/n + (B/n)
  r.hat <- v.hat / W
  return(r.hat)
 }

## ----echo=FALSE---------------------------------------------------------------
set.seed(123)
N <- 15000

sigma2 <- .1
X1 <- matrix(c(RandomWalkMCMC(N=N,sigma2 = sigma2,X1=-10)$x,
              RandomWalkMCMC(N=N,sigma2 = sigma2,X1=-5)$x,
              RandomWalkMCMC(N=N,sigma2 = sigma2,X1=5)$x,
              RandomWalkMCMC(N=N,sigma2 = sigma2,X1=10)$x), nrow=4, byrow = T)
N <- length(X1[1,])
psi1 <- t(apply(X1, 1, cumsum))
for (i in 1:nrow(psi1)){
   psi1[i,] <- psi1[i,] / (1:ncol(psi1))
}

R1 <- sapply(2:N, FUN = function(k) Gelman.Rubin(psi1[,1:k]))

sigma2 <- .5
X2 <- matrix(c(RandomWalkMCMC(N=N,sigma2 = sigma2,X1=-10)$x,
              RandomWalkMCMC(N=N,sigma2 = sigma2,X1=-5)$x,
              RandomWalkMCMC(N=N,sigma2 = sigma2,X1=5)$x,
              RandomWalkMCMC(N=N,sigma2 = sigma2,X1=10)$x), nrow=4, byrow = T)
N <- length(X2[1,])
psi2 <- t(apply(X2, 1, cumsum))
for (i in 1:nrow(psi2)){
   psi2[i,] <- psi2[i,] / (1:ncol(psi2))
}

R2 <- sapply(2:N, FUN = function(k) Gelman.Rubin(psi2[,1:k]))



sigma2 <- 2
X3 <- matrix(c(RandomWalkMCMC(N=N,sigma2 = sigma2,X1=-10)$x,
              RandomWalkMCMC(N=N,sigma2 = sigma2,X1=-5)$x,
              RandomWalkMCMC(N=N,sigma2 = sigma2,X1=5)$x,
              RandomWalkMCMC(N=N,sigma2 = sigma2,X1=10)$x), nrow=4, byrow = T)
N <- length(X1[1,])
psi3 <- t(apply(X3, 1, cumsum))
for (i in 1:nrow(psi3)){
   psi3[i,] <- psi3[i,] / (1:ncol(psi3))
}

R3 <- sapply(2:N, FUN = function(k) Gelman.Rubin(psi3[,1:k]))



sigma2 <- 16
X4 <- matrix(c(RandomWalkMCMC(N=N,sigma2 = sigma2,X1=-10)$x,
              RandomWalkMCMC(N=N,sigma2 = sigma2,X1=-5)$x,
              RandomWalkMCMC(N=N,sigma2 = sigma2,X1=5)$x,
              RandomWalkMCMC(N=N,sigma2 = sigma2,X1=10)$x), nrow=4, byrow = T)
N <- length(X4[1,])
psi4 <- t(apply(X4, 1, cumsum))
for (i in 1:nrow(psi4)){
   psi4[i,] <- psi4[i,] / (1:ncol(psi4))
}

R4 <- sapply(2:N, FUN = function(k) Gelman.Rubin(psi4[,1:k]))

## -----------------------------------------------------------------------------
b <- 200 # burn-in

plot(R1[-(1:b)],type='l',ylab='R',
     main = 'sigma^2=0.1')
abline(h=1.2, col='red')

plot(R2[-(1:b)],type='l',ylab='R',
     main = 'sigma^2=0.5')
abline(h=1.2, col='red')

plot(R3[-(1:b)],type='l',ylab='R',
     main = 'sigma^2=2')
abline(h=1.2, col='red')

plot(R4[-(1:b)],type='l',ylab='R',
     main = 'sigma^2=16')
abline(h=1.2, col='red')


## -----------------------------------------------------------------------------

plot(psi1[1,-(1:b)],type = "l", ylab=bquote(psi),
     ylim=c(min(psi1[,-(1:b)]),max(psi1[,-(1:b)])),
     main = 'sigma^2=0.1')
lines(psi1[2,-(1:b)],col=2)
lines(psi1[3,-(1:b)],col=3)
lines(psi1[4,-(1:b)],col=4)
abline(v=min(which(R1[-(1:b)]<1.2)),lwd=1.5)

plot(psi2[1,-(1:b)],type = "l", ylab=bquote(psi),
     ylim=c(min(psi2[,-(1:b)]),max(psi2[,-(1:b)])),
     main = 'sigma^2=0.5')
lines(psi2[2,-(1:b)],col=2)
lines(psi2[3,-(1:b)],col=3)
lines(psi2[4,-(1:b)],col=4)
abline(v=min(which(R2[-(1:b)]<1.2)),lwd=1.5)

plot(psi3[1,-(1:b)],type = "l", ylab=bquote(psi),
     ylim=c(min(psi3[,-(1:b)]),max(psi3[,-(1:b)])),
     main = 'sigma^2=2')
lines(psi3[2,-(1:b)],col=2)
lines(psi3[3,-(1:b)],col=3)
lines(psi3[4,-(1:b)],col=4)
abline(v=min(which(R3[-(1:b)]<1.2)),lwd=1.5)

plot(psi4[1,-(1:b)],type = "l", ylab=bquote(psi),
     ylim=c(min(psi4[,-(1:b)]),max(psi4[,-(1:b)])),
     main = 'sigma^2=16')
lines(psi4[2,-(1:b)],col=2)
lines(psi4[3,-(1:b)],col=3)
lines(psi4[4,-(1:b)],col=4)
abline(v=min(which(R4[-(1:b)]<1.2)),lwd=1.5)


## -----------------------------------------------------------------------------
Gibbs <- function(N,a,b,n,X1){
  x <- numeric(N)
  y <- numeric(N)
  x[1] <- X1
  
  for(i in 1:N){
    y[i] <- rbeta(1,x[i]+a,n-x[i]+b)
    x[i+1] <- rbinom(1,n,y[i])
  }
  x <- x[1:N]
  return(data.frame(x=x,y=y))
}

## ----echo=FALSE---------------------------------------------------------------
set.seed(123)
N <- 3000
a <- 1
b <- 2
n <- 15


data1 <- Gibbs(N,a,b,n,X1=0)
data2 <- Gibbs(N,a,b,n,X1=8)
data3 <- Gibbs(N,a,b,n,X1=15)

X <- matrix(c(data1$x,data2$x,data3$x), nrow=3, byrow = T)
psiX <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psiX)){
   psiX[i,] <- psiX[i,] / (1:ncol(psiX))
}

RX <- sapply(2:N, FUN = function(k) Gelman.Rubin(psiX[,1:k]))

Y <- matrix(c(data1$y,data2$y,data3$y), nrow=3, byrow = T)
psiY <- t(apply(Y, 1, cumsum))
for (i in 1:nrow(psiY)){
   psiY[i,] <- psiY[i,] / (1:ncol(psiY))
}

RY <- sapply(2:N, FUN = function(k) Gelman.Rubin(psiY[,1:k]))

## -----------------------------------------------------------------------------
b <- 500

plot(RX[-(1:b)],type = 'l',ylab = 'RX')
abline(h=1.2,col='red')

plot(RY[-(1:b)],type = 'l',ylab='RY')
abline(h=1.2,col='red')

plot(psiX[1,-(1:b)],type = "l", ylab='psiX',
     ylim=c(min(psiX[,-(1:b)]),max(psiX[,-(1:b)])))
lines(psiX[2,-(1:b)],col=2)
lines(psiX[3,-(1:b)],col=3)
abline(v=min(which(RX[-(1:b)]<1.2)),lwd=1.5)

plot(psiY[1,-(1:b)],type = "l", ylab='psiY',
     ylim=c(min(psiY[,-(1:b)]),max(psiY[,-(1:b)])))
lines(psiY[2,-(1:b)],col=2)
lines(psiY[3,-(1:b)],col=3)
abline(v=min(which(RY[-(1:b)]<1.2)),lwd=1.5)


## -----------------------------------------------------------------------------
estimate <- function(N,b1,b2,b3,f0){
    x1 <- rpois(N,1)
    x2 <- rexp(N,1)
    x3 <- rbinom(N,1,0.5)
    uniroot(f = function(alpha){
        r <- rbinom(N,1,prob = 1/(1+exp(-(alpha+b1*x1+b2*x2+b3*x3))))
        mean(r)-f0
    }, interval = c(-30,0))$root
    
}

## -----------------------------------------------------------------------------
N <- 1e6
b1 <- b2 <- 1
b3 <- -1

alpha <- c(estimate(N,b1,b2,b3,f0 = 0.1),
estimate(N,b1,b2,b3,f0 = 0.01),
estimate(N,b1,b2,b3,f0 = 0.001),
estimate(N,b1,b2,b3,f0 = 0.0001))

alpha

## -----------------------------------------------------------------------------
plot(-log(c(.1,.01,.001,.0001)),alpha,
     xlab = '-logf0',ylab = bquote(alpha))

## -----------------------------------------------------------------------------
set.seed(0)
obj <- c(4,2,9)
c.mat <- rbind(matrix(c(2,1,1,-1,1,3),nrow = 2),
               diag(1,3,3))
c.dir <- c(rep('<=',2), rep('>=',3))
c.rhs <- c(2,3,0,0,0)

result <- lp('min',obj, c.mat,c.dir,c.rhs)
# 最小值
result$objval
# x,y,z的取值
result$solution

## -----------------------------------------------------------------------------
u <- c(11,8,27,13,16,0,23,10,24,2)
v <- u+1
n <- length(u)
likelihood_o <- function(lambda){
  sum(log(exp(-lambda*u)-exp(-lambda*v)))
}
dlo <- function(lambda){
  sum( (v*exp(-lambda*v)-u*exp(-lambda*u))/(exp(-lambda*u)-exp(-lambda*v)) )
}
ddlo <- function(lambda){
  -sum((u-v)^2*exp(-(u+v)*lambda)/(exp(-u*lambda)-exp(-v*lambda))^2)
}

## -----------------------------------------------------------------------------
lambda_old  <- 0
lambda_new <- 10
iter <- 0
max_iter <- 200
ts <- NULL
lambda_NT <- numeric(max_iter)
while (abs((lambda_old-lambda_new))>1e-8 && iter<max_iter) {
    lambda_old <- lambda_new
    lambda_NT[iter+1] <- lambda_new
    t <- 1
    delta <- dlo(lambda_new)/ddlo(lambda_new)
    while(lambda_new-t*dlo(lambda_new)/ddlo(lambda_new)<0) t <- t/2
    lambda_new <- lambda_new-t*dlo(lambda_new)/ddlo(lambda_new)
    iter <- iter+1
    ts <- c(ts,t)
}
lambda_NT <- lambda_NT[1:iter]
lambda_new

## -----------------------------------------------------------------------------
lambda_old  <- 0
lambda_new <- 10
iter <- 0
max_iter <- 200
lambda_EM <- numeric(max_iter)
while (abs((lambda_old-lambda_new))>1e-8 && iter<max_iter) {
  lambda_old <- lambda_new
  lambda_EM[iter+1] <- lambda_new
  lambda_new <- n/(n/lambda_new-dlo(lambda_new))
  iter <- iter+1
}
lambda_EM <- lambda_EM[1:iter]
lambda_new

## ----echo=FALSE---------------------------------------------------------------
plot(log10(abs(lambda_NT-rev(lambda_NT)[1])),
     ylab = 'log(x-x*)',main = 'Newton-Raphson')

plot(log10(abs(lambda_EM-rev(lambda_EM)[1])),
     ylab = 'log(x-x*)',main = 'EM')

## -----------------------------------------------------------------------------
set.seed(0)
lapply2 <- function(x, f, ...){
  m <- length(x)
  n <- x[[1]] %>% unlist %>% f %>% length
  unlist(Map(function(z) vapply(z, f, rep(0,n), ...)
             , x)) %>% structure('dim'=c(n,m))
}

## -----------------------------------------------------------------------------
chisq.test2 <- function(x,y){
  A <- table(x,y) %>% as.matrix
  n <- sum(A)
  nc <- ncol(A)
  nr <- nrow(A)
  Ac <- matrix(rep(colSums(A),each=nr),nr,nc)
  Ar <- matrix(rep(rowSums(A),nc),nr,nc)
  n*(sum(A*A/Ac/Ar)-1)
}

## -----------------------------------------------------------------------------
table2 <- function(x,y){
  x0 <- unique(x) %>% sort
  y0 <- unique(y) %>% sort
  T <- sapply(y0, function(Y){
    tmp <- x[y==Y]
    vapply(x0, function(X){
        sum(tmp==X)
    },FUN.VALUE = 0L)
})
  class(T) <- 'table'
  T
}

## -----------------------------------------------------------------------------
chisq.test3 <- function(x,y){
  A <- table2(x,y) %>% as.matrix
  n <- sum(A)
  nc <- ncol(A)
  nr <- nrow(A)
  Ac <- matrix(rep(colSums(A),each=nr),nr,nc)
  Ar <- matrix(rep(rowSums(A),nc),nr,nc)
  n*(sum(A*A/Ac/Ar)-1)
}

## -----------------------------------------------------------------------------
set.seed(0)
cppFunction('
NumericVector MCMCC(int N) {
  auto f = [](double x) -> double {
        return pow(2+x, 125) * pow(1-x, 38) * pow(x, 34) * (x > 0) * (x < 1);
    };
	NumericVector X(N);
	X[0] = 0.6;
	double u = 0,y = 0;
	for(int i = 1; i < N; i++) {
    y = rnorm(1,X[i-1],pow(0.5,2))[0];
		u = runif(1)[0];
		X[i]=((u<f(y)/f(X[i-1]))?y:X[i-1]);
  	}
	return(X);
}')

## -----------------------------------------------------------------------------
f <- function(x){
  (2+x)^125*(1-x)^38*x^34*(x<1)*(x>0)/2.357695e+28
}

MCMCR <- function(N){
  X <- numeric(N)
  X[1] <- 0.6
  u <- runif(N)
  for(i in 2:N){
    y <- rnorm(1,X[i-1],0.5^2)
    if(u[i]<f(y)/f(X[i-1])){
      X[i] <- y
    } else X[i] <- X[i-1]
  }
  X
}

## -----------------------------------------------------------------------------
N <- 20000
b <- 7000

xC <- MCMCC(N)
hist(xC[-(1:b)], freq = F,
     main = 'simulation with Rcpp')
lines(seq(0,1,0.01),f(seq(0,1,0.01)),col='red')

xR <- MCMCR(N)
hist(xR[-(1:b)], freq = F,
     main = 'simulaton with R')
lines(seq(0,1,0.01),f(seq(0,1,0.01)),col='red')

## -----------------------------------------------------------------------------
qqplot(xC[-(1:b)],xR[-(1:b)])

## -----------------------------------------------------------------------------
microbenchmark(MCMCC(N),MCMCR(N))

