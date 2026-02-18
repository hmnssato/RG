k=1
ndat=51

if (k==1){
  dat <- read.csv("mean_contact_frequency.BRD3179_40000.csv",header = TRUE);
  rd1 <- fread("norm_contact.BRD3179_40000.csv")
  nn <- "BRD3179_40000"
}else if (k==2){
  dat <- read.csv("mean_contact_frequency.BRD3179N_40000.csv",header = TRUE);
  rd1 <- fread("norm_contact.BRD3179N_40000.csv")
  nn <- "BRD3179N_40000"
}else if (k==3){
  dat <- read.csv("mean_contact_frequency.BRD3187_40000.csv",header = TRUE);
  rd1 <- fread("norm_contact.BRD3187_40000.csv")
  nn <- "BRD3187_40000"
}else if (k==4){
  dat <- read.csv("mean_contact_frequency.BRD3187N_40000.csv",header = TRUE);
  rd1 <- fread("norm_contact.BRD3187N_40000.csv")
  nn <- "BRD3187N_40000"
}
dat <- log10(dat);
dat<-dat[order(dat[,1]),];
plot(dat);
y<-10^dat;
x<-seq(y[5,1],10^6.7,(10^6.7-y[5,1])/(ndat-1));

bb_number <- 1000
rd1 <- data.frame(rd1)
ll <- length(rd1[,1]) 
res_dat <- log10(x)
for (i in 1:bb_number){
  cc <- sample(c(1:ll),6200, replace=TRUE)
  xx <- rd1[cc,3]
  yy <- rd1[cc,4]
  zz <- cbind(xx,yy)
  
  zz2 <- zz[order(zz[,1]),]
  #semi log scale: only CP is log scale
  sp <- smooth.spline(zz2[,1],log10(zz2[,2]),spar=0)
  resp <- predict(sp, x)
  
  res_dat <- cbind(res_dat, resp$y)
}

res_dat <- data.frame(res_dat)
write.csv(res_dat, paste(paste(paste("boot_dat2","_", sep=""),nn, sep=""),".csv",sep=""))

rr_mean=rep(0, length(for_rr_mean[,1]))
for (i in 1:length(for_rr_mean[,1])){
  rr_mean[i]=mean(for_rr_mean[i,])
}
res_dat_mod <- cbind(log10(x), log10(rr_mean))
write.csv(res_dat_mod, paste(paste(paste("mean_boot_dat","_", sep=""),nn, sep=""),".csv",sep=""))
