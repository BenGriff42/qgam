n <- 1e4
qu <- 0.5

dat <- mgcv::gamSim(1, n = n, dist = "normal", scale = 10)
form <- y ~ s(x0) + s(x1) + s(x2)

co <- 0.85
sig <- exp(1.097) #* 0.1

# fit1 <- gam(form, family = qgam::elf(theta = log(sig), co = co, qu = qu), data = dat,
#             control = list(trace = TRUE),
#             sp = c(1.052212, 0.2049904, 0.002499400))
# fit2 <- gam(form, family = qgam.old::elf(theta = log(sig), co = co, qu = qu), data = dat,
#             control = list(trace = TRUE),
#             sp = c(1.052212, 0.2049904, 0.002499400))
# 
# plot(fit1$fitted.values, fit2$fitted.values)
# abline(0, 1, col = 2)


fit1 <- gam(form, family = qgam::elf( co = co, qu = qu), data = dat,
            control = list(trace = TRUE))

fit2 <- gam(form, family = qgam.old::elf(co = co, qu = qu), data = dat,
            control = list(trace = TRUE))

plot(fit1$fitted.values, fit2$fitted.values)
abline(0, 1, col = 2)
