bp <- 1.86191; #for patient
bn <- 1.80954; #for healthy

tt <- seq(0,0.5,0.01);

msdp <- tt^(1-1/bp);
msdn <- tt^(1-1/bn);

plot(x=tt, y=msdp, type="l", xlab="t", ylab="MSD", ylim=c(0,0.8), col = 1, lty = 1);
par(new=TRUE);
plot(x=tt, y=msdn, type="l", xlab="t", ylab="MSD", ylim=c(0,0.8), col = 1, lty = 2);
legend("topleft", legend=c("patiehnt","healthy"),col=c(1,1), lty = c(1,2));

aa <- seq(0.5,1,0.01);
dfp <- 2*(aa/(1-1/bp)-1);
dfn <- 2*(aa/(1-1/bn)-1);

plot(x=aa, y=dfp, type="l", xlab="alpha", ylab="df", ylim=c(0,2.5), col = 1, lty = 1);
par(new=TRUE);
plot(x=aa, y=dfn, type="l", xlab="alpha", ylab="df", ylim=c(0,2.5), col = 1, lty = 2);
legend("topleft", legend=c("patient","healthy"),col=c(1,1), lty = c(1,2));

b30 <- 1.286723468;
b15 <- 1.324001791;
b10 <- 1.492210317;



msd30 <- tt^(1-1/b30);
msd15 <- tt^(1-1/b15);
msd10 <- tt^(1-1/b10);

plot(x=tt, y=msd30, type="l", xlab="t", ylab="MSD", ylim=c(0,0.9), col = 1, lty = 1);
par(new=TRUE);
plot(x=tt, y=msd15, type="l", xlab="t", ylab="MSD", ylim=c(0,0.9), col = 1, lty = 2);
par(new=TRUE);
plot(x=tt, y=msd10, type="l", xlab="t", ylab="MSD", ylim=c(0,0.9), col = 1, lty = 3);
legend("topleft", legend=c("30","15,15","10,10,10"),col=c(1,1,1), lty = c(1,2,3));

df30 <- 2*(aa/(1-1/b30)-1);
df15 <- 2*(aa/(1-1/b15)-1);
df10 <- 2*(aa/(1-1/b10)-1);

plot(x=aa, y=df30, type="l", xlab="alpha", ylab="df", ylim=c(0,7), col = 1, lty = 1);
par(new=TRUE);
plot(x=aa, y=df15, type="l", xlab="alpha", ylab="df", ylim=c(0,7), col = 1, lty = 2);
par(new=TRUE);
plot(x=aa, y=df10, type="l", xlab="alpha", ylab="df", ylim=c(0,7), col = 1, lty = 3);
legend("topleft", legend=c("30","15,15","10,10,10"),col=c(1,1,1), lty = c(1,2,3));
