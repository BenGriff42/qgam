co <- 0.8
sig <- exp(1.2)
qu <- 0.1

# co <- co0
# sig <- exp(theta0)

fam1 <- elf(theta = log(sig), co = co, qu = qu)
fam2 <- qgam.old::elf(theta = log(sig), co = co, qu = qu)

y <- seq(-10, 10, by = 0.01)
mu <- y*0

# dev.resids
plot(y,fam1$dev.resids(y, mu, 1),type = "l")
lines(y,fam2$dev.resids(y, mu, 1),col = 2, lty = 2)

# Dd
Dd1 <- fam1$Dd(y, mu, theta = log(sig), wt = 1, level = 2)
Dd2 <- fam2$Dd(y, mu, theta = log(sig), wt = 1, level = 2)

plot(y,Dd1$Dmu,type = "l")
lines(y,Dd2$Dmu,col = 2, lty = 2)

plot(y,Dd1$Dmu2,type = "l")
lines(y,Dd2$Dmu2,col = 2, lty = 2)

plot(y,Dd1$Dmu3,type = "l")
lines(y,Dd2$Dmu3,col = 2, lty = 2)

plot(y,Dd1$Dmu4,type = "l")
lines(y,Dd2$Dmu4,col = 2, lty = 2)

plot(y,Dd1$Dth,type = "l")
lines(y,Dd2$Dth,col = 2, lty = 2)

plot(y,Dd1$Dth2,type = "l")
lines(y,Dd2$Dth2,col = 2, lty = 2)

plot(y,Dd1$Dmuth,type = "l")
lines(y,Dd2$Dmuth,col = 2, lty = 2)

plot(y,Dd1$Dmu2th,type = "l")
lines(y,Dd2$Dmu2th,col = 2, lty = 2)

plot(y,Dd1$Dmuth2,type = "l")
lines(y,Dd2$Dmuth2,col = 2, lty = 2)

plot(y,Dd1$Dmu2th2,type = "l")
lines(y,Dd2$Dmu2th2,col = 2, lty = 2)

plot(y,Dd1$Dmu3th,type = "l")
lines(y,Dd2$Dmu3th,col = 2, lty = 2)

#############################

dd1 <- dDeta(y, y*0, wt = 1, theta = log(sig),fam = fam1, deriv = 0)
dd2 <- dDeta(y, y*0, wt = 1, theta = log(sig),fam = fam2, deriv = 0)

str(dd1)
str(dd2)

sum(abs(dd1$Deta - dd2$Deta))
sum(abs(dd1$Deta2 - dd2$Deta2))
sum(abs(dd1$EDeta2 - dd2$EDeta2))
sum(abs(dd1$Deta.Deta2 - dd2$Deta.Deta2))
sum(abs(dd1$Deta.EDeta2 - dd2$Deta.EDeta2))

#################################
n <- 1e4
dat <- mgcv::gamSim(1, n = n, dist = "normal", scale = 10)
y <- dat$y
f <- dat$f


fam1$ls(y, f, theta = log(sig), scale = 1)
fam2$ls(y,  f, theta = log(sig), scale = 1)

thetas <- seq(-2, 2, length = 100)
ls1 <- ls2 <- rep(NA, 100)
for(i in 1:100){
  ls1[i] <- fam1$ls(y, f, theta = thetas[i], scale = 1)$ls
  ls2[i] <- fam2$ls(y, f, theta = thetas[i], scale = 1)$ls
}

plot(thetas, ls1, type = "l")
lines(thetas, ls2, col = 2, lty = 2)


thetas <- seq(-2, 2, length = 100)
ls1 <- lsth <- rep(NA, 100)
for(i in 1:100){
  ls1[i] <- fam1$ls(y, f, theta = thetas[i], scale = 1)$ls
  lsth[i] <- fam1$ls(y, f, theta = thetas[i], scale = 1)$lsth1
}

plot(thetas, lsth, type = "l")
lines((thetas[1:99] + thetas[2:100]) / 2,(ls1[2:100] - ls1[1:99]) / (thetas[2] - thetas[1]), col = 2, lty = 2)


thetas <- seq(-2, 2, length = 100)
lsth <- lsth2 <- rep(NA, 100)
for(i in 1:100){
  lsth[i] <- fam1$ls(y, f, theta = thetas[i], scale = 1)$lsth1
  lsth2[i] <- fam1$ls(y, f, theta = thetas[i], scale = 1)$lsth2
}

plot(thetas, lsth2, type = "l")
lines((thetas[1:99] + thetas[2:100]) / 2,(lsth[2:100] - lsth[1:99]) / (thetas[2] - thetas[1]), col = 2, lty = 2)

# elf second gradient is positive, clearly should be negative