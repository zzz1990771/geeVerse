#' Cross-Validation for Generalized Estimating Equations (GEE)
#'
#' This function performs k-fold cross-validation for model selection in the context
#' of Generalized Estimating Equations (GEE). It is designed to evaluate the performance
#' of different models specified by a range of lambda values, choosing the one that
#' minimizes the cross-validation criterion.
#'
#' @param formula an object of class \code{"formula"} (or one that can be coerced to that class):
#' a symbolic description of the model to be fitted.
#' @param id a vector which identifies the cluster/group for each observation.
#' @param data an optional data frame containing the variables in the model.
#' @param family a description of the error distribution and link function to be used in the model.
#' @param scale.fix logical; if \code{TRUE}, the scale parameter is fixed to \code{scale.value}.
#' @param scale.value the value of the scale parameter when \code{scale.fix} is \code{TRUE}.
#' @param fold the number of folds to be used in the cross-validation.
#' @param lambda.vec a vector of lambda values for which the cross-validation error will be calculated.
#' @param pindex an optional numeric vector specifying a parameter index.
#' @param eps the threshold for convergence criteria.
#' @param maxiter the maximum number of iterations for the convergence of the algorithm.
#' @param tol the tolerance level for the convergence of the algorithm.
#'
#' @return An object of class \code{"CVfit"}, which is a list containing:
#' \describe{
#'   \item{\code{fold}}{The number of folds used in the cross-validation.}
#'   \item{\code{lam.vect}}{The vector of lambda values tested.}
#'   \item{\code{cv.vect}}{The cross-validation error for each lambda.}
#'   \item{\code{lam.opt}}{The lambda value that resulted in the minimum cross-validation error.}
#'   \item{\code{cv.min}}{The minimum cross-validation error.}
#'   \item{\code{call}}{The matched call.}
#' }
#'
#' @importFrom stats gaussian
#' @importFrom stats model.extract
#' @importFrom stats model.matrix
#'
#' @export
CVfit <-
function(formula, id, data, family, scale.fix, scale.value, fold,
lambda.vec, pindex, eps, maxiter, tol) {

call_cv <- match.call()
mf <- match.call(expand.dots=FALSE)

mf$family <- mf$link <- mf$varfun <-
mf$scale.fix <- mf$scale.value <-
mf$fold <- mf$lambda.vec <-mf$pindex <-
mf$eps <- mf$maxiter <- mf$tol <- NULL

if(is.null(mf$id)) mf$id <- as.name("id")

mf[[1]] <- as.name("model.frame")
mf <- eval(mf, parent.frame())
Terms <- attr(mf, "terms")
y <- model.extract(mf, "response")
X <- model.matrix(Terms, mf)
id <- model.extract(mf, id)

if(missing(family)) family=gaussian(link="identity")

#if (is.character(family)) family <- get(family)
#if (is.function(family))  family <- family()

if(missing(pindex)) pindex <- NULL
if(missing(scale.fix))   scale.fix <- TRUE
if(missing(scale.value)) scale.value <- 1
if(missing(eps)) eps <- 10^-6
if(missing(maxiter)) maxiter <- 30
if(missing(tol)) tol <- 10^-3

N <- length(unique(id))
K <- dim(X)[2]-1
nx <- dim(X)[2]
nt <- dim(X)[1]/N
nt <- rep(nt,N)

#pay attention to this part
lam.min <- -1
cv.min <- Inf
cv.vect <- NULL

for (j in 1:length(lambda.vec))   {
#get one of the lambda's
lam.temp <- lambda.vec[j]
#initial value for cv.values is 0.
cv.value <- 0

for (k in 1:fold) {

#select the index that will be omitted.
index.cv <- ((k-1)*nt[1]*(N/fold)+1):(k*nt[1]*(N/fold))

#training part#
y.train <- y[-index.cv]
if (colnames(X)[1]== "(Intercept)") x.train <- X[-index.cv,-1] else x.train <- X[-index.cv,]
id.train <- id[-index.cv]

#compute beta.train
data.train=data.frame("id"=id.train,"y"=y.train, x.train)
mm<-match.call(expand.dots = FALSE)
mm$formula<-formula
mm$id<-id.train
mm$data<-data.train
mm$na.action<-"na.omit"
mm$family<-family
mm$corstr<-"independence"
mm$Mv<-NULL
mm$beta_int<-rep(0,nx)  ##NULL##
mm$R<-NULL
mm$scale.fix<-scale.fix
mm$scale.value<-scale.value
mm$lambda<-lam.temp
mm$pindex<-pindex
mm$eps<-eps
mm$maxiter<-maxiter
mm$tol<-tol
mm$silent<-TRUE
mm$lambda.vec<-NULL
mm$fold<-NULL
#mm$link <- NULL
#mm$varfun<NULL

mm[[1]]<-as.name("PGEE")

beta.train <- eval(mm, parent.frame())$coefficients

#testing part##
y.cv<-y[index.cv]
x.cv<-X[index.cv,]
id.cv<-id[index.cv]
yy=y.cv
eta=x.cv%*%beta.train
mu=family$linkinv(eta)
##family$dev.resids gives the square of the residuals
cv.value<-cv.value+sum((family$dev.resids(yy,mu,wt=1)))
} #k

cv.vect<-c(cv.vect, cv.value)

if(cv.value<cv.min) {
lam.min<-lam.temp
cv.min<-cv.value
}

} #j

out<-list()
attr(out, "class") <- c("CVfit")
out$fold=fold
out$lam.vect=lambda.vec
out$cv.vect=cv.vect
out$lam.opt=lam.min
out$cv.min=cv.min
out$call <- call_cv
out
}
