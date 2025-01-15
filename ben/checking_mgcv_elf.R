load("ben/fitting_stuff.RData")
load("ben/mObj.RData")

# old elf version
mObj.e <- mObj
mObj.e$family <- qgam.old::elf(theta = mObj$family$getTheta(), co = mObj$family$getCo(), qu = mObj$family$getQu())

# New call_lists
call_list1 <-  call_list
call_list1$G <- quote(mObj1)
call_list2 <-  call_list
call_list2$G <- quote(mObj2)

theta0 <- mObj$family$getTheta()
co0 <- mObj$family$getCo()

################ define lam1 and sig1 ##################

co <- co0
theta <- theta0

mObj1 <- mObj
mObj1$family$putTheta(theta)
mObj1$family$putCo(co)

debugonce(mgcv:::gam.fit4)
m1 <- do.call(gam, call_list1)

mObj2 <- mObj.e
mObj2$family$putTheta(theta)
mObj2$family$putCo(co)

debugonce(mgcv:::gam.fit4)
m2 <- do.call(gam, call_list2)

m1$iter
m2$iter

plot(m1$fitted.values, m2$fitted.values)
abline(0,1,col=2)
################ define lam2 and sig2 ##################


lam2 <- 0.85
sig2 <- exp(1.097)

lam2 <- 1.2
sig2 <- exp(1)

lam1 <- lam2 * sig2
sig1 <- sig2

mObj1 <- mObj
mObj1$family$putTheta(log(sig1))
mObj1$family$putCo(lam1)

m1 <- do.call(gam, call_list1)

mObj2 <- mObj.e
mObj2$family$putTheta(log(sig2))
mObj2$family$putCo(lam2 * sig2)

m2 <- do.call(gam, call_list2)


################# check plot ######################

l <- function(y, sig, lam, tau) {
  ((tau - 1) * y + lam * log(1 + exp(y / lam))) / sig
}
ol <- function(y, sig, lam, tau) {
  (tau - 1) * y / sig + lam * log(1 + exp(y / (lam * sig)))
}

y <- seq(-10, 10, by = 0.01)

tau <- mObj$family$getQu()
plot(y, pinLoss(y / exp(theta), 0, tau, add = F), type = "l")
lines(y, l(y, exp(theta), co, tau), col = 2)
lines(y, ol(y, exp(theta), co / exp(theta), tau), col = 3, lty = 2)

############## dDeta #########################

# gam.fit4 uses dDeta with deriv 0. I think these should be the same

