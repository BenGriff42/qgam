# install on master
# load_all on new_elf

# use qgam.old::qgam for original elf
# use qgam::qgam for new elf

# set.seed(4)
n <- 1e4
qu <- 0.1

dat <- mgcv::gamSim(1,n=n,dist="normal",scale=10)
form <- y ~ s(x0) + s(x1) + s(x2) + s(x3)

fit.elf <- qgam.old::qgam(form = form, data = dat, qu = qu, discrete = TRUE)
fit.nelf <- qgam::qgam(form = form, data = dat, qu = qu, discrete = TRUE)

sum( (qnorm(qu, dat$f, 10) - fit.elf$fitted.values)^2 )
sum( (qnorm(qu, dat$f, 10) - fit.nelf$fitted.values)^2 )

plot(fit.elf$fitted.values, fit.nelf$fitted.values)
plot(qnorm(qu, dat$f, 10), fit.nelf$fitted.values)
abline(0,1,col=2)

c(fit.elf$aic, fit.nelf$aic)
c(fit.elf$family$getTheta(), fit.nelf$family$getTheta())
c(fit.elf$family$getCo(), fit.nelf$family$getCo())


dat$y <- dat$y * ((dat$x2) * 3 / 5 + 0.7)
form <- list(y ~ s(x0) + s(x1) + s(x2) + s(x3), ~ s(x2))

fit.elf <- qgam.old::qgam(form = form, data = dat, qu = qu, discrete = TRUE)
fit.nelf <- qgam::qgam(form = form, data = dat, qu = qu, discrete = TRUE)

pinLoss(dat$y, fit.elf$fitted.values, qu = qu)
pinLoss(dat$y, fit.nelf$fitted.values, qu = qu)

plot(fit.elf$fitted.values, fit.nelf$fitted.values)
abline(0,1,col=2)

fit.elf$aic
fit.nelf$aic
