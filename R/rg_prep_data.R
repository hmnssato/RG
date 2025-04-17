library(GenKern);

#generate global data
f1<-function(k){
  if (k==1){
    dat <- read.csv("mean_contact_frequency.BRD3179_40000.csv",header = TRUE);
  }else if (k==2){
    dat <- read.csv("mean_contact_frequency.BRD3179N_40000.csv",header = TRUE);
  }else if (k==3){
    dat <- read.csv("mean_contact_frequency.BRD3187_40000.csv",header = TRUE);
  }else if (k==4){
    dat <- read.csv("mean_contact_frequency.BRD3187N_40000.csv",header = TRUE);
  }
  dat<-dat[order(dat[,1]),];
  plot(dat);
  datt<-log10(dat);
  plot(datt);
  ll<-nearest(datt[,1], 5.2);
  #ll<-nearest(datt[,1], 5.6);
  kk<-nearest(datt[,1], 6.7);
  x<-datt[ll:kk,1];
  y<-datt[ll:kk,2];
  llx<-stats::smooth.spline(x,y,spar=0);
  xx<-log10(seq(dat[ll,1],dat[kk,1],(dat[kk,1]-dat[ll,1])/49));
  yy<-stats::predict(llx,xx)$y;
  datl<-cbind(xx,yy);
  plot(datl, xlab="log of distance (bp)", ylab="log of CP");
  if (k==1){
    write.csv(datl,"dis_data_BRD3179_50_g.csv");
  }else if (k==2){
    write.csv(datl,"dis_data_BRD3179N_50_g.csv");
  }else if (k==3){
    write.csv(datl,"dis_data_BRD3187_50_g.csv");
  }else if (k==4){
    write.csv(datl,"dis_data_BRD3187N_50_g.csv");
  }
}

f1(4);

#generate semi-local data (3, 5)
f2<-function(k){
  if (k==1){
    dat <- read.csv("mean_contact_frequency.BRD3179_40000.csv",header = TRUE);
  }else if (k==2){
    dat <- read.csv("mean_contact_frequency.BRD3179N_40000.csv",header = TRUE);
  }else if (k==3){
    dat <- read.csv("mean_contact_frequency.BRD3187_40000.csv",header = TRUE);
  }else if (k==4){
    dat <- read.csv("mean_contact_frequency.BRD3187N_40000.csv",header = TRUE);
  }
  dat<-dat[order(dat[,1]),];
  datt<-log10(dat);
  ll<-nearest(datt[,1], 5.2);
  kk<-nearest(datt[,1], 6.7);
  x<-datt[ll:kk,1];
  y<-datt[ll:kk,2];
  llx<-stats::smooth.spline(x,y,spar=0);
  ss <- 5000;
  mm <- nearest(datt[,1],log10(500000));
  xx <- log10(10^datt[mm,1]+c(-2:2)*ss);
  yy<-stats::predict(llx,xx)$y;
  datl<-cbind(xx,yy);
  xxx <- log10(10^datt[mm,1]+c(-1:1)*2*ss);
  yyy<-stats::predict(llx,xxx)$y;
  datlg<-cbind(xxx,yyy);
  plot(datl, xlab="distance log(bp)", ylab="log(CP)");
  plot(datlg, xlab="distance log(bp)", ylab="log(CP)");
  if (k==1){
    write.csv(datl,"dis_data_BRD3179_5_sl.csv");
    write.csv(datlg,"dis_data_BRD3179_3_sl.csv");
  }else if (k==2){
    write.csv(datl,"dis_data_BRD3179N_5_sl.csv");
    write.csv(datlg,"dis_data_BRD3179N_3_sl.csv");
  }else if (k==3){
    write.csv(datl,"dis_data_BRD3187_5_sl.csv");
    write.csv(datlg,"dis_data_BRD3187_3_sl.csv");
  }else if (k==4){
    write.csv(datl,"dis_data_BRD3187N_5_sl.csv");
    write.csv(datlg,"dis_data_BRD3187N_3_sl.csv");
  }
}

f2(1);

#generate semi-local data (30)
f3<-function(k){
  if (k==1){
    dat <- read.csv("mean_contact_frequency.BRD3179_40000.csv",header = TRUE);
  }else if (k==2){
    dat <- read.csv("mean_contact_frequency.BRD3179N_40000.csv",header = TRUE);
  }else if (k==3){
    dat <- read.csv("mean_contact_frequency.BRD3187_40000.csv",header = TRUE);
  }else if (k==4){
    dat <- read.csv("mean_contact_frequency.BRD3187N_40000.csv",header = TRUE);
  }
  dat<-dat[order(dat[,1]),];
  datt<-log10(dat);
  ll<-nearest(datt[,1], 5.2);
  kk<-nearest(datt[,1], 6.7);
  x<-datt[ll:kk,1];
  y<-datt[ll:kk,2];
  llx<-stats::smooth.spline(x,y,spar=0);
  ss <- 5000;
  mm <- nearest(datt[,1],log10(500000));
  xx <- log10(10^datt[mm,1]+c(-15:14)*ss);
  yy<-stats::predict(llx,xx)$y;
  datl<-cbind(xx,yy);
  plot(datl, xlab="distance log(bp)", ylab="log(CP)");
  if (k==1){
    write.csv(datl,"dis_data_BRD3179_30_sl.csv");
  }else if (k==2){
    write.csv(datl,"dis_data_BRD3179N_30_sl.csv");
  }else if (k==3){
    write.csv(datl,"dis_data_BRD3187_30_sl.csv");
  }else if (k==4){
    write.csv(datl,"dis_data_BRD3187N_30_sl.csv");
  }
}

f3(1);
