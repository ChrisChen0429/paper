## Solutions checked on 2.14.0 on Ubuntu. I've used
## `##' before comments, and `#' before expected R output
## Simon Wood Nov 2012

## Q1. Note that `machine precision' is slightly more slippery than
##     is suggested here. The question focuses on the precision implied
##     by the storage of floating point numbers in main memory. The part 
##     of the CPU which carries out floating point operations often 
##     works at a higher precision, carrying more significant figures 
##     in intermendiate quantities. Compilers sometimes use these extra 
##     precision floating point storage locations as temporary storage 
##     to speed up code, which can mean that a program has some of its
##     floating point numbers stored at higher precision than others. 
##     The upshot is that it is usually a bad idea to test floating point
##     numbers for exact equality, even if they should be exactly equal. 
##     Do read Section 1.1 of Press et al. (2007) Numerical Recipes, CUP.
## a)
eps <- 1
x <- 1 
while (x+eps != x) eps <- eps/2
eps/x
# 1.110223e-16

## b)
eps+x==x
# [1] TRUE
eps1 <- eps*1.01
eps1+x==x
# [1] FALSE
## ... confirming the statement

## c) 
.Machine$double.eps

## d) 
eps <- 1
x <- .125
while (x+eps != x) eps <- eps/2
eps/x
#[1] 1.110223e-16

## e) 
eps <- 1
x <- pi
while (x+eps != x) eps <- eps/2
eps/x
#[1] 7.067899e-17

## f)
eps <- 1
x <- 1 
while (x+eps != x) eps <- eps/2
log(1/(2*eps))/log(10) 
#[1] 15.65356

## Q2.
## a)
X <- matrix(runif(100000),1000,100)
z <- rep(0,1000)
system.time(
for (i in 1:1000) {
  for (j in 1:100) z[i] <- z[i] + X[i,j]
})
system.time(z1 <- rowSums(X))
system.time(z2 <- apply(X,1,sum))
range(z-z1)
range(z-z2)

## b)
n <- 100000
z <- rnorm(n)
proc.time()
zneg <- 0;j <- 1
for (i in 1:n) {
  if (z[i]<0) {
    zneg[j] <- z[i]
    j <- j + 1
  }
}
proc.time()
zneg1 <- z[z<0] ## fast alternative
proc.time()
range(zneg1-zneg)

## Q3. 

## basic manipulations....
set.seed(1)
n <- 1000
A <- matrix(runif(n*n),n,n)
x <- runif(n)

xtAx <- t(x)%*%A%*%x

trA <- sum(diag(A)) 

AtXA <- t(A)%*%diag(x)%*%A

AtXA.2 <- t(A)%*%(x*A) ## alternative to above, using `recycling' 

## Q4 
## a)
set.seed(0)
n <- 1000
A <- matrix(runif(n*n),n,n)
x.true <- runif(n)
y <- A%*%x.true

## b)
proc.time()
Ai <- solve(A)
x1 <- Ai%*%y
proc.time() ## took about 4 seconds
mean(abs(x1-x.true))
#[1] 2.869841e-11

## c)
system.time(x2 <- solve(A,y)) ## .8 seconds
mean(abs(x2-x.true))
#[1] 1.356142e-12

## d) It seems that method c) is both faster and more accurate (which is what 
##    theory would confirm).

## Q5. Many solutions are possible. Here is one that does the job.

obs.cdf <- function(x,plot.cdf=TRUE) {
  sx <- sort.int(x,index.return=TRUE) ## order the observations
  n <- length(x)
  sx$F <- 1:n/n  ## \hat F for the ordered data
  if (plot.cdf) {
    dx <- (max(x)-min(x))/20 ## how far to plot beyond data range
    xp <- c(sx$x[1]-dx,sx$x,sx$x[n]+dx) ## the x values for plotting
    xp <- rep(xp,rep(3,n+2))[-1] ## repeat the x values for step function plot
    Fp <- c(0,sx$F,1)  ## \hat F values for plotting
    Fp <- rep(Fp,rep(3,n+2))[-3*(n+2)] ## repeat in same way as x
    plot(xp,Fp,type="l",xlab="x",ylab="F(x)") ## plot!
  }
  sx$F[sx$ix] <- sx$F  ## \hat F in original data order
  sx$F ## return \hat F
}
obs.cdf(rnorm(30))

## Q6... NA

## Q7. 
## a)
rb <- function(x,z) {
  100*(z-x^2)^2 + (1-x)^2 
}

## b)
n <- 100
x <- seq(-1.5,1.5,length=n)
z <- seq(-.5,1.5,length=n)
f <- outer(x,z,rb)
contour(x,z,matrix(f,n,n))

## c)
contour(x,z,matrix(log10(f),n,n),levels=(1:10/2))

## d) 
rb.grad <- function(x,z) {
## gradient of rosenbrocks function at a single point
  g <- rep(NA,2)
  g[2] <- 200*(z-x^2)
  g[1] <- 400*(x^3-z*x) + 2*(x-1)
  g
}

## e)
x0 <- .5; z0 <- 1; eps <- 1e-7
f0 <- rb(x0,z0)
g <- g0 <- rb.grad(x0,z0) ## exact gradiant
## put FD gradiant in g0
g0[1] <- (rb(x0+eps,z0)-f0)/eps
g0[2] <- (rb(x0,z0+eps)-f0)/eps
g;g0 ## compare

## f) 
rb.hess <- function(x,z) {
## Hessian of Rosenbrock's function, at a single point
  H <- matrix(NA,2,2)
  H[2,2] <- 200
  H[1,1] <- 1200*x^2 - 400*z + 2
  H[2,1] <- H[1,2] <- -400*x
  H
}

## g)
H <- H0 <- rb.hess(x0,z0)
H0[,1] <- (rb.grad(x0+eps,z0)-g)/eps
H0[,2] <- (rb.grad(x0,z0+eps)-g)/eps
H;H0

## h)

ad.quad <- function(x,z,col=2,trans=log10,lev=1:10/2) {
## function to add a quadratic approximation to Rosenbrocks
## to a plot of the function, based on a Taylor expansion
## at x,z
  f0 <- rb(x,z)
  g <- rb.grad(x,z)
  H <- rb.hess(x,z)
  qap <- function(x,z,f0,g,H) {
    X <- matrix(c(x,z),length(x),2)
    f0 + X%*%g + rowSums((X%*%H)*X)/2
  }
  n <- 100
  xx <- seq(-1.5,1.5,length=n)
  zz <- seq(-.5,1.5,length=n)
  q <- trans(outer(xx-x,zz-z,qap,f0=f0,g=g,H=H)) 
  contour(xx,zz,matrix(q,n,n),col=col,add=TRUE,levels=lev)
}

## i)

contour(x,z,matrix(log10(f),n,n),levels=(1:10/2),lwd=2)
ad.quad(-1,0.5);points(-1,.05,col=2,pch=20,cex=2)
ad.quad(0,0,col=3);points(0,0,col=3,pch=20,cex=2)
ad.quad(1,1,col=4);points(1,1,col=4,pch=20,cex=2)

## j)
contour(x,z,matrix(log10(f),n,n),levels=(1:10/2),lwd=2)
ad.quad(0.5,0.5);points(.5,.5,col=2,pch=20,cex=2)

## Q8.
## a) 
?optim

## b) versions of rb and rb.grad for `optim' use

rbo <- function(x) {
  100*(x[2]-x[1]^2)^2 + (1-x[1])^2 
}

rbg <- function(x) rb.grad(x[1],x[2])

## c)
optim(c(-.5,1),rbo) ## reports successful convergence

## d)
optim(c(-.5,1),rbo,method="BFGS")
## ... modest improvement

## e)
optim(c(-.5,1),rbo,rbg,method="BFGS")
## ... takes fewer evaluations to converge to much higher accuracy

## f) 
optim(c(-.5,1),rbo,rbg,method="CG",control=list(maxit=11000))
## ... quite slow convergence in this case.

