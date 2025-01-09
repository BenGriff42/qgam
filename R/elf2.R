elf2 <- function(theta = NULL, link = "identity", qu, co){

  # Some checks
  if (!is.na(qu) && (findInterval(qu, c(0, 1) ) != 1)) stop("qu should be in (0, 1)")

  linktemp <- substitute(link)
  if (!is.character(linktemp)) linktemp <- deparse(linktemp)
  if (linktemp %in% c("log", "identity", "sqrt")) stats <- make.link(linktemp)
  else if (is.character(link)) {
    stats <- make.link(link)
    linktemp <- link
  } else {
    if (inherits(link, "link-glm")) {
      stats <- link
      if (!is.null(stats$name))
        linktemp <- stats$name
    }
    else stop(linktemp, " link not available for elf family; available links are \"identity\", \"log\" and \"sqrt\"")
  }
  ## Theta <-  NULL;
  n.theta <- 1
  if ( !is.null(theta) ) {

    iniTheta <- theta ## fixed log theta supplied

    n.theta <- 0 ## signal that there are no theta parameters to estimate

  } else iniTheta <- 0 ## inital log theta value

  env <- new.env(parent = environment(elf)) #.GlobalEnv) ##########!!!!!!!!!!!!!!!!~########################

  assign(".Theta", iniTheta, envir = env)
  getTheta <- function(trans=FALSE) if (trans) exp(get(".Theta")) else get(".Theta")
  putTheta <- function(theta) assign(".Theta", theta, envir=environment(sys.function()))

  assign(".qu", qu, envir = env)
  getQu <- function( ) get(".qu")
  putQu <- function(qu) assign(".qu", qu, envir=environment(sys.function()))

  assign(".co", co, envir = env)
  getCo <- function( ) get(".co")
  putCo <- function(co) assign(".co", co, envir=environment(sys.function()))

  # variance <- function(mu) exp(get(".Theta"))  ##### XXX ##### Necessary?

  validmu <- function(mu) all( is.finite(mu) )

  dev.resids <- function(y, mu, wt, theta=NULL) {        ##### XXX #####
    if( is.null(theta) ) theta <- get(".Theta")
    tau <- get(".qu")
    co <- get(".co")

    mu <- drop(mu)

    sig <- exp(theta)
    lam <- mean(co / sig)
    sig <- co / lam

    term <-
      (1 - tau)*lam*log1p(-tau) +
      lam*tau*log(tau) -
      (1 - tau)*(y - mu) +
      lam*log1pexp((y - mu) / lam)

    2 * wt * term / sig
  }

  Dd <- function(y, mu, theta, wt, level=0) {

    tau <- get(".qu")
    co <- get(".co")
    mu <- drop(mu)

    ## derivatives of the deviance...
    sig <- exp(theta)
    lam <- mean(co / sig)
    sig <- co / lam

    z <- (y - mu) / lam

    der <- sigmoid(z, deriv = TRUE)

    r <- list()
    ## get the quantities needed for IRLS.
    ## Dmu2eta2 is deriv of D w.r.t mu twice and eta twice,
    ## Dmu is deriv w.r.t. mu once, etc...
    r$Dmu <- -2 * wt * ( (der$D0 - 1 + tau) / sig )
    r$Dmu2 <- 2 * wt * ( der$D1 / (sig * lam) )
    # r$EDmu2 <- 2 * wt * ((1-tau)*tau / (lam + 1)) / sig^2 ## exact (or estimated) expected weight #### XXX ####
    r$EDmu2 <- r$Dmu2 # It make more sense using the observed information everywhere
    if (level > 0) { ## quantities needed for first derivatives

      r$Dth <- -dev.resids(y, mu = mu, wt = wt, theta = theta)
      r$Dmuth <- -r$Dmu
      r$Dmu3 <- -(2 * wt * der$D2) / (sig * lam^2)
      r$Dmu2th <- -r$Dmu2
    }
    if (level > 1) { ## whole damn lot
      r$Dmu4 <- (2 * wt * der$D3) / (sig * lam^3)
      r$Dth2 <- -r$Dth
      r$Dmuth2 <- r$Dmu
      r$Dmu2th2 <- r$Dmu2
      r$Dmu3th <- -r$Dmu3
    }
    r
  }

  aic <- function(y, mu, theta=NULL, wt, dev) {
    if (is.null(theta)) theta <- get(".Theta")
    sig <- exp(theta)
    tau <- get(".qu")
    co <- get(".co")
    lam <- mean(co / sig)
    sig <- co / lam

    mu <- drop(mu)

    term <-
      -(1 - tau) * (y - mu) / sig +
      lam * log1pexp( (y - mu) / lam ) / sig +
      log(lam * beta(lam * (1 - tau) / sig, tau * lam / sig))

    2 * sum(term * wt)
  }

  ls <- function(y, w, theta, scale) {
    tau <- get(".qu")
    co <- get(".co")
    sig <- exp(theta)
    lam <- mean(co / sig)
    sig <- co / lam

    ## the log saturated likelihood function.
    ls <- sum(w * (
      (1 - tau) * lam * log1p(-tau) / sig +
        lam * tau * log(tau) / sig -
        log(lam * beta(lam * (1 - tau) / sig, lam * tau / sig))
    ))

    lsth <-
      (
        lam * (1 - tau) * digamma(lam * (1 - tau) / sig) +
          lam * tau * digamma(lam * tau / sig) -
          lam * digamma(lam / sig)
      ) / sig

    lsth2 <-
      -lsth -
      (
        lam^2 * (1 - tau)^2 * trigamma(lam * (1 - tau) / sig) +
          lam^2 * tau*2 * trigamma(lam * tau / sig) -
          lam^2 * trigamma(lam / sig)
      ) / sig

    list(ls = ls, ## saturated log likelihood
         lsth1 = lsth, ## first deriv vector w.r.t theta - last element relates to scale, if free
         lsth2 = lsth2) ## Hessian w.r.t. theta, last row/col relates to scale, if free
  }

  initialize <- expression({

    mustart <- quantile(y, family$getQu()) + y * 0 # this ---> y <--- is very bad idea

  })

  #postproc <- expression({  ####### XXX ??? #######
  #  object$family$family <-
  #    paste("elf(",round(object$family$getTheta(TRUE),3),")",sep="")
  #})

  #   rd <- function(mu,wt,scale) {  ####### XXX TODO #######
  #     Theta <- exp(get(".Theta"))
  #     rnbinom(mu,size=Theta,mu=mu)
  #   }
  #
  #   qf <- function(p,mu,wt,scale) {  ####### XXX TODO #######
  #     Theta <- exp(get(".Theta"))
  #     qnbinom(p,size=Theta,mu=mu)
  #   }

  get.null.coef <- function(G,start=NULL,etastart=NULL,mustart=NULL,...) {
    ## Get an estimate of the coefs corresponding to maximum reasonable deviance...
    y <- G$y
    weights <- G$w
    nobs <- G$n ## ignore codetools warning!!
    ##start <- etastart <- mustart <- NULL
    family <- G$family
    eval(family$initialize) ## have to do this to ensure y numeric
    y <- as.numeric(y)
    mum <- quantile(y, get(".qu")) + 0*y
    etam <- family$linkfun(mum)
    null.coef <- qr.coef(qr(G$X), etam)
    null.coef[is.na(null.coef)] <- 0;
    ## get a suitable function scale for optimization routines
    null.scale <- sum(family$dev.resids(y, mum, weights))/nrow(G$X)
    list(null.coef=null.coef,null.scale=null.scale)
  }


  #  environment(rd)<- environment(qf) <- environment(variance) <-
  environment(dev.resids) <- environment(ls) <- environment(aic) <- environment(Dd) <-
    environment(getTheta) <-
    environment(putTheta) <- environment(putCo) <- environment(getCo) <-
    environment(putQu) <- environment(getQu) <- environment(get.null.coef) <- env

  structure(list(family = "elf", link = linktemp, linkfun = stats$linkfun,
                 linkinv = stats$linkinv, dev.resids = dev.resids,Dd=Dd,
                 #variance=variance,
                 aic = aic, mu.eta = stats$mu.eta, initialize = initialize,
                 #postproc=postproc,
                 ls=ls,
                 validmu = validmu, valideta = stats$valideta, n.theta=n.theta,
                 ini.theta = iniTheta, putTheta=putTheta,getTheta=getTheta,
                 putQu=putQu, getQu=getQu,
                 putCo=putCo,getCo=getCo, get.null.coef=get.null.coef,
                 use.wz=TRUE
                 #, rd=rd,qf=qf
  ),
  class = c("extended.family","family"))
} ## elf
