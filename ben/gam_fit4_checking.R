load("ben/gam_fit4_output.RData")
fam1 <- aaa$family

load("ben/gam_fit4_output2.RData")
fam2 <- aaa$family

y <- aaa$y
mu <- aaa$mu
weights <- aaa$weights
theta <- aaa$theta


dd1 <- dDeta(y,mu,weights,theta,fam1,0)
sum(!is.finite(dd1$Deta.Deta2))

r <- fam1$Dd(y,mu,theta,weights,0)
sum(!is.finite(r$Dmu/r$Dmu2))
sum(!is.finite(r$Dmu))
sum(!is.finite(1/r$Dmu2))
sum(r$Dmu2 == 0)


lam <- fam1$getCo()
sig <- exp(fam1$getTheta())
der <- sigmoid((y - mu) / lam, deriv = TRUE)
sum(( der$D1 ) == 0)

dl <- dlogis(y-mu, 0, lam)
sum( dl == 0)



dd2 <- dDeta(y,mu,weights,theta,fam2,0)
sum(!is.finite(dd2$Deta.Deta2))
