zx=1

if (zx==1){
  fname <- "BRD3179_50";
}else if (zx==2){
  fname <- "BRD3179N_50";
}else if (zx==3){
  fname <- "BRD3187_50";
}else if (zx==4){
  fname <- "BRD3187N_50";
}

param_l <- read.csv(paste(paste("inf_param_df_", fname, sep=""),"_l.csv", sep =""), header = TRUE);
theta <- param_l[1,2];
beta <- param_l[2,2];
nsc <- param_l[3,2];
kbt <- param_l[4,2];
lsc <- param_l[5:length(param_l[,2]),2];

dat <- read.csv(paste(paste("dis_data_", fname, sep=""),"_l.csv", sep =""), header = TRUE);
x<- as.numeric(10^dat[1:length(dat[,1]),2]);

ss <- x[2]-x[1];
ssr <- 21.8/0.3;

nx <- ss/ssr;
ssq <- nx*ss;
kk <- kbt/ss^2;

kf <- 0.01*kbt;
r0 <- 1.6*ss;

epsilon1=0.01;
epsilon2=0.01;

dx <- ssr;
x <- seq(ssr,ss,dx);
dy <- ss;
y <- seq(ss,ssq,dy);

phi <- rep(kk/2*(ss-ssr)^2,length(x));
#phi <- -kf*r0^2/2*log(1-x^2/r0^2);
for (i in 1:length(x)){
  if (x[i]<2^(1/6)*ssr){
    phi[i] <- phi[i]+epsilon1*4*kbt*(ssr^12/x[i]^12-ssr^6/x[i]^6)
  }
}

#phiy <- rep(kk/2*ss^2,length(x));
phiy <- epsilon2*(4.*kbt*(ssr^12/y^12-ssr^6/y^6));

Zx <- sum(exp(-phi/kbt)*dx);
Zy <- sum(exp(-phiy/kbt)*dy);
  
rhox <- exp(-phi/kbt)*dx/Zx;
rhoy <- exp(-phiy/kbt)*dy/Zy;

fex <- sum(phi*rhox*dx)+kbt*sum(rhox*log(rhox)*dx);
fey <- sum(phiy*rhoy*dy)+kbt*sum(rhoy*log(rhoy)*dy);

entx <- sum(rhox*log(rhox)*dx/kbt);
enty <- sum(rhoy*log(rhoy)*dy/kbt);

xx <- rep(0,length(rhox));
cumx <- rep(0,length(rhox));
cumy <- rep(0,length(rhoy));
for (i in 1:length(rhox)){
  xx[i]=i*dx/(nx*dx);
  cumx[i] <- sum(rhox[1:i]);
  cumy[i] <- sum(rhoy[1:i]);
}
plot(x=xx, y=cumx, xlab="scaled distance", ylab="cummurative probability density", xlim=c(0,1), ylim=c(0,1), col=1, type="l", lty = 1);
par(new=T);
plot(x=xx, y=cumy, xlab="scaled distance", ylab="cummurative probability density", xlim=c(0,1), ylim=c(0,1), col=2, type="l", lty = 2);
legend("topleft", legend=c("No extrusion","Extrusion"),col=c(1,2), lty = c(1,2));
plot(x=xx, y=cumx-cumy, xlab="scaled distance", ylab="difference of cummurative probability density", xlim=c(0,1), col=1, type="l", lty = 1);
