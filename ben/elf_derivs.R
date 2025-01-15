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

