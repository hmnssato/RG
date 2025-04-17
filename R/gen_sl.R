library(reticulate);
use_python("set your python folder")

#generate semi-local cp (3,5)

zx=1;

if (zx==1){
  fname <- "BRD3179_5";
  ffname <-  "BRD3179_3";
}else if (zx==2){
  fname <- "BRD3179N_5";
  ffname <-  "BRD3179N_3";
}else if (zx==3){
  fname <- "BRD3187_5";
  ffname <-  "BRD3187_3";
}else if (zx==4){
  fname <- "BRD3187N_5";
  ffname <-  "BRD3187N_3";
}

param_sl <- read.csv(paste(paste("inf_param_df_", fname, sep=""),"_sl.csv", sep =""), header = TRUE);
theta_sl <- param_sl[1,2];
beta_sl <- param_sl[2,2];
nsc_sl <- param_sl[3,2];
kbt_sl <- param_sl[4,2];
lsc_sl <- param_sl[5:length(param_sl[,2]),2];

dat_sl <- read.csv(paste(paste("dis_data_", fname, sep=""),"_sl.csv", sep =""), header = TRUE);
x<- as.numeric(10^dat_sl[,2]);

dat_sl1 <- t(read.csv(paste(paste("yre_df_", fname, sep=""),"_sl.csv", sep =""), header = TRUE));
yre <- as.numeric(dat_sl1[2:(length(dat_sl1[,2])),2]);
mi_yre=min(yre);

Nx <- length(x);
ss <- abs(x[2]-x[1]);
kk <- kbt_sl/ss^2;

recr_sl<-matrix(rep(0,Nx*Nx),Nx,Nx);
recr_sl1<-matrix(rep(0,Nx*Nx),Nx,Nx);

rate<-1000;
rrr <- rep(0,rate);

lls = stats::smooth.spline(x-x[1],lsc_sl*(x[2]-x[1]),spar = 0);

for (k in 1:Nx){
  for (l in 1:Nx){
    if(k!=l){
      lx=abs(x[k]-x[l]);
      scc_l = stats::predict(lls,lx)$y;
      lxx = scc_l;
      r1 = lx - lxx;
      r2 = lx + lxx;
      del_r = (r2 - r1) / (rate-1);  
      for (m in 1:rate){
        rrr[m] = r1+del_r*m;
      }
      numpy <- import("numpy");
      kkk <- numpy$fft$fftfreq(length(rrr), del_r);
      kx <- kkk[11:(rate-10)]
      ggg = 1 / (4 * kk * ((kx / lxx^nsc_sl) ^ 2 / 2 - (kx / lxx^nsc_sl) ^ 4 / 24)) * (sin(kx / 2) / sin(kx / lxx^nsc_sl / 2)) ^ beta_sl;
      gy = numpy$fft$ifft(ggg);
      gyc= 1 / lxx^((2 * theta_sl + 1)*nsc_sl) * Re(gy);
      recr_sl[k,l] = gyc[floor(rate/2)-10];
    }else{
      recr_sl[k,l] = 10^(mi_yre-5);
    }
  }
}

recp_sl_max=10^(mi_yre-5);
for (k in 1:Nx){
  for (l in 1:Nx){
    if (is.nan(recr_sl[k,l])==FALSE && recr_sl[k,l]>recp_sl_max){
      recp_sl_max <- recr_sl[k,l];
    }
  }
}

recp_sl_max=max(c(10^(-10),recp_sl_max))

diag(recr_sl)=rep(recp_sl_max,length(recr_sl[1,]));

for (k in 1:Nx){
  for (l in 1:Nx){
    recr_sl1[k,Nx-l+1]=recr_sl[k,l]
  }
}

recr_sl_log <- log10(recr_sl1);

write.csv(recr_sl1,paste(paste("cp_",fname,sep=""),"_sl.csv",sep=""));
write.csv(recr_sl_log,paste(paste("cp_",fname,sep=""),"_sl_log.csv",sep=""));

write.csv(cbind(c("theta", theta_sl),c("beta", beta_sl),c("n", nsc_sl),c("kbt", kbt_sl)), paste(paste("param_",fname,sep=""),"_sl.csv",sep = ""));
write.csv(lsc_sl,paste(paste("lsc_",fname,sep=""),"_sl.csv",sep = ""));

param_sl <- read.csv(paste(paste("inf_param_df_", ffname, sep=""),"_sl.csv", sep =""), header = TRUE);
theta_sl <- param_sl[1,2];
beta_sl <- param_sl[2,2];
nsc_sl <- param_sl[3,2];
kbt_sl <- param_sl[4,2];
lsc_sl <- param_sl[5:length(param_sl[,2]),2];

dat_sl <- read.csv(paste(paste("dis_data_", ffname, sep=""),"_sl.csv", sep =""), header = TRUE);
x<- as.numeric(10^dat_sl[,2]);

dat_sl1 <- t(read.csv(paste(paste("yre_df_", ffname, sep=""),"_sl.csv", sep =""), header = TRUE));
yre <- as.numeric(dat_sl1[2:(length(dat_sl1[,2])),2]);
mi_yre=min(yre);

Nx <- length(x);
ss <- abs(x[2]-x[1]);
kk <- kbt_sl/ss^2;

recr_sl<-matrix(rep(0,Nx*Nx),Nx,Nx);
recr_sl1<-matrix(rep(0,Nx*Nx),Nx,Nx);

rate<-1000;
rrr <- rep(0,rate);

library(Hmisc);

for (k in 1:Nx){
  for (l in 1:Nx){
    if(k!=l){
      lx=abs(x[k]-x[l]);
      scc_l = approxExtrap(x-x[1],lsc_sl*(x[2]-x[1]),lx)$y;
      lxx = scc_l;
      r1 = lx - lxx;
      r2 = lx + lxx;
      del_r = (r2 - r1) / (rate-1);  
      for (m in 1:rate){
        rrr[m] = r1+del_r*m;
      }
      numpy <- import("numpy");
      kkk <- numpy$fft$fftfreq(length(rrr), del_r);
      kx <- kkk[11:(rate-10)]
      ggg = 1 / (4 * kk * ((kx / lxx^nsc_sl) ^ 2 / 2 - (kx / lxx^nsc_sl) ^ 4 / 24)) * (sin(kx / 2) / sin(kx / lxx^nsc_sl / 2)) ^ beta_sl;
      gy = numpy$fft$ifft(ggg);
      gyc= 1 / lxx^((2 * theta_sl + 1)*nsc_sl) * Re(gy);
      recr_sl[k,l] = gyc[floor(rate/2)-10];
    }else{
      recr_sl[k,l] = 10^(mi_yre-5);
    }
  }
}

recp_sl_max=10^(mi_yre-5);
for (k in 1:Nx){
  for (l in 1:Nx){
    if (is.nan(recr_sl[k,l])==FALSE && recr_sl[k,l]>recp_sl_max){
      recp_sl_max <- recr_sl[k,l];
    }
  }
}

recp_sl_max=max(c(10^(-10),recp_sl_max))

diag(recr_sl)=rep(recp_sl_max,length(recr_sl[1,]));

for (k in 1:Nx){
  for (l in 1:Nx){
    recr_sl1[k,Nx-l+1]=recr_sl[k,l]
  }
}

recr_sl_log <- log10(recr_sl1);

write.csv(recr_sl1,paste(paste("cp_",ffname,sep=""),"_sl.csv",sep=""));
write.csv(recr_sl_log,paste(paste("cp_",ffname,sep=""),"_sl_log.csv",sep=""));

write.csv(cbind(c("theta", theta_sl),c("beta", beta_sl),c("n", nsc_sl),c("kbt", kbt_sl)), paste(paste("param_",ffname,sep=""),"_sl.csv",sep = ""));
write.csv(lsc_sl,paste(paste("lsc_",ffname,sep=""),"_sl.csv",sep = ""));
