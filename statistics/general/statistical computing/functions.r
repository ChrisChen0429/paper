## Functions for sparse Markov linear models:
##   linear.basis      (evaluate piecewise linear basis functions)
##   derivation.matrix (approximate squared derivative penalties)
##   qinv              (obtain elements of the inverse of a sparse matrix)
##
## Finn Lindgren 24/12/2014
## For use in APTS 2014/15 assessment
## See the comments preceding each function for explanations

## linear.basis
##
## Evaluate piecewise linear basis functions
##
##   A <-  linear.basis(x, knots)
## returns a length(x)-by-length(knots) matrix with length(knots)
## piecewise linear basis functions evaluated at the locations in the
## vector 'x'.  'knots' must be given as an increasing sequence.
## The first and last basis functions are defined to be constant
## outside of range(knots)
##
## Runtime ~ O(nx + nk) for optimal implementation

linear.basis <- function(x, knots) {
  require(Matrix)
  n <- length(x)
  nk <- length(knots)
  if (nk < 2) {
    stop("'knots' must have length >= 2")
  }
  A <- Matrix(0.0, n, nk)
  ## Left half-interval
  ii <- which(x <= knots[1])
  if (length(ii) > 0) {
    A[ii, 1] <- 1.0
  }
  ## Interior intervals
  dk <- diff(knots)
  for (k in seq_len(nk - 1)) {
    ii <- which((x >= knots[k]) & (x < knots[k + 1]))
    if (length(ii) > 0) {
      A[ii, k] <- (knots[k + 1] - x[ii]) / dk[k]
      A[ii, k + 1] <- (x[ii] - knots[k]) / dk[k]
    }
  }
  ## Right half-interval
  ii <- which(x >= knots[nk])
  if (length(ii) > 0) {
    A[ii, nk] <- 1.0
  }
  A
}

## derivation.matrix
##
## Construct a derivation operator matrix, suitable for integrated
## squared derivative penalty matrix construction.
##
##   der <- derivation.matrix(knots, order)
## returns a  list(D, dx)  such that
##   sum((der$D %*% f)^2 * der$dx)
## is an approximation of
##   \int (d^o f / dx^o)^2 dx, where o=order,
## and the vector f is evaluated at the locations given by 'knots'.
## knots must be given as an increasing sequence.
##
## The corresponding quadratic form matrix is given by
##   DD <- t(der$D) %*% Diagonal(length(der$dx), der$dx) %*% der$D
## so that
##   sum((der$D %*% f)^2 * der$dx) = t(f) %*% DD %*% f
##
## Runtime ~ O(nk * order) for optimal implementation

derivation.matrix <- function(knots, order=2) {
  nk <- length(knots)
  if (!((1 <= order) && (order < nk))) {
    stop(paste("order = ", order, ", but must be in [1, ", nk-1, "]", sep=""))
  }
  dx <- diff(knots)
  D <- (1 / dx) * diff(Diagonal(nk))
  for (k in seq_len(order-1)) {
    ii <- seq_len(length(dx) - 1)
    dx <- (dx[ii] + dx[ii + 1]) / 2
    D <- (1 / dx) * diff(D)
  }
  ## dx is here constructed heuristically to give the correct total
  ## penalty for simple polynomials.
  list(D=D, dx=dx * diff(range(knots)) / sum(dx))
}

## qinv
##
## Calculate variances and neighbour covariances
## Uses recursion formulas based on the Cholesky decomposition,
## see Takahashi (1973), Rue, Martino (2007), Rue, Martino, Chopin (2009)
##
## This is a simple and slow R implementation of INLA::inla.qinv(),
## see http://r-inla.org/. For a 1D Markov model, the runtime
## of qinv(Q), ncol(Q)==1,000 is approximately the same as the runtime
## of inla.qinv(Q), ncol(Q)==1,000,000
##
## The call
##   S <- qinv(Q)
## gives the elements needed to calculate
##   diag(Q^-1) = diag(S),
## and
##   dlog(det(Q))/dt = tr(Q^-1 dQ/dt) = tr(S dQ/dt),
## where dQ/dt is a derivative of Q with respect to a parameter t
## (assuming that Q_ij == 0 implies (dQ/dt)_ij == 0)
##
## The call
##   S <- qinv(Q, sparse=FALSE)
## gives the full inverse of Q.  Use only to check that
##   S == solve(Q)
## up to small numerical deviations.
##
## Runtime ~ O(nnz(L)) for optimal implementation, where Q=LL'

qinv <- function(Q, sparse=TRUE, perm=TRUE) {
  require("Matrix")
  n <- nrow(Q)
  ## P Q P' = L L'
  if (!is(Q, "Matrix")) {
    ch <- expand(Matrix::Cholesky(as(Q, "Matrix"), perm=perm))
  } else {
    ch <- expand(Matrix::Cholesky(Q, perm=perm))
  }
  ## S = P Q^-1 P'
  ## S_ij = I(i=j) / L_ii^2 - sum_{k=i+1}^n L_ki S_kj / L_ii,
  ##   j \geq i, i=n,...,1
  S <- Diagonal(n, 1 / diag(ch$L)^2)
  for (i in (n-1):1) {
    kk <- ((i+1):n)[ch$L[(i+1):n, i] != 0.0]
    Lki <- ch$L[kk, i]
    jj <- n:i
    if (sparse) { ## Only compute elements needed for the variances
      jj <- jj[ch$L[jj, i] != 0.0]
    }
    for (j in jj) {
      S[j, i] <- S[i, j] <- S[i, j] - sum(Lki * S[kk, j]) / ch$L[i, i]
    }
  }
  ## Q^-1 = P' S P
  t(ch$P) %*% S %*% ch$P
}
