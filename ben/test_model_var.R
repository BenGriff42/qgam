library(MASS)
data(mcycle)

fit1 <- qgam::qgam(accel~s(times, k=20, bs="ad"), 
                   data = mcycle, 
                   qu = 0.8)
fit2 <- qgam.old::qgam(accel~s(times, k=20, bs="ad"), 
                       data = mcycle, 
                       qu = 0.8)
xSeq <- data.frame(cbind("accel" = rep(0, 1e3), "times" = seq(2, 58, length.out = 1e3)))
pred1 <- predict(fit1, newdata = xSeq, se=TRUE)
pred2 <- predict(fit2, newdata = xSeq, se=TRUE)
plot(mcycle$times, mcycle$accel, xlab = "Times", ylab = "Acceleration", ylim = c(-150, 80))
lines(xSeq$times, pred1$fit, lwd = 1)
lines(xSeq$times, pred1$fit + 2*pred1$se.fit, lwd = 1, lty = 2)
lines(xSeq$times, pred1$fit - 2*pred1$se.fit, lwd = 1, lty = 2)   
lines(xSeq$times, pred2$fit, lwd = 1, col = 2)
lines(xSeq$times, pred2$fit + 2*pred2$se.fit, lwd = 1, lty = 2, col = 2)
lines(xSeq$times, pred2$fit - 2*pred2$se.fit, lwd = 1, lty = 2, col = 2)   

fit3 <- qgam::qgam(list(accel ~ s(times, k=20, bs="ad"), ~ s(times)),
                   data = mcycle, 
                   qu = 0.8)

fit4 <- qgam.old::qgam(list(accel ~ s(times, k=20, bs="ad"), ~ s(times)),
                       data = mcycle, 
                       qu = 0.8)
pred3 <- predict(fit3, newdata = xSeq, se=TRUE)
pred4 <- predict(fit4, newdata = xSeq, se=TRUE)
plot(mcycle$times, mcycle$accel, xlab = "Times", ylab = "Acceleration", ylim = c(-150, 80))
lines(xSeq$times, pred3$fit, lwd = 1)
lines(xSeq$times, pred3$fit + 2*pred3$se.fit, lwd = 1, lty = 2)
lines(xSeq$times, pred3$fit - 2*pred3$se.fit, lwd = 1, lty = 2)   
lines(xSeq$times, pred4$fit, lwd = 1, col = 2)
lines(xSeq$times, pred4$fit + 2*pred4$se.fit, lwd = 1, lty = 2, col = 2)
lines(xSeq$times, pred4$fit - 2*pred4$se.fit, lwd = 1, lty = 2, col = 2)   





set.seed(651)
n <- 2000
x <- seq(-4, 3, length.out = n)
X <- cbind(1, x, x^2)
beta <- c(0, 1, 1)
sigma =  1.2 + sin(2*x)
f <- drop(X %*% beta)
dat <- f + rnorm(n, 0, sigma)
dataf <- data.frame(cbind(dat, x))
names(dataf) <- c("y", "x")

fit1 <- qgam::qgam(list(y~s(x, k = 30, bs = "cr"), ~ s(x, k = 30, bs = "cr")), 
            data = dataf, qu = 0.95)
fit2 <- qgam.old::qgam(list(y~s(x, k = 30, bs = "cr"), ~ s(x, k = 30, bs = "cr")), 
                   data = dataf, qu = 0.95)

plot(x, dat, col = "grey", ylab = "y")
tmp1 <- predict(fit1, se = TRUE)
tmp2 <- predict(fit2, se = TRUE)
lines(x, tmp1$fit)
lines(x, tmp1$fit + 2 * tmp1$se.fit, lty = 2)
lines(x, tmp1$fit - 2 * tmp1$se.fit, lty = 2)
lines(x, tmp2$fit, col = 2)
lines(x, tmp2$fit + 2 * tmp2$se.fit, col = 2, lty = 2)
lines(x, tmp2$fit - 2 * tmp2$se.fit, col = 2, lty = 2)

