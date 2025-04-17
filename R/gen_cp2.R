library(reticulate);
use_python("set your python folder")

#generate semi-local cp (15, 15)

zx=1;

if (zx==1){
  fname <-  "BRD3179_30";
}else if (zx==2){
  fname <-  "BRD3179N_30";
}else if (zx==3){
  fname <-  "BRD3187_30";
}else if (zx==4){
  fname <-  "BRD3187N_30";
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
yre <- as.numeric(dat_sl1[2:length(dat_sl1[,2]),2]);
min_yre=min(yre);


ss <- x[2]-x[1];
kk <- kbt_sl/ss^2;
Nx <- length(x);

recr_sl <- matrix(rep(0,Nx*Nx),Nx,Nx);
recr_sl1 <- matrix(rep(0,Nx*Nx),Nx,Nx);
recr_sl1_log <- matrix(rep(0,Nx*Nx),Nx,Nx);

lls = stats::smooth.spline(x-x[1],lsc_sl*(x[2]-x[1]),spar = 0);

rate <- 1000;
rrr <- rep(0,rate);

for (k in 1:floor(Nx*0.5)){
  for (l in 1:floor(Nx*0.5)){
    if(k!=l){
      lx=abs(x[k]-x[l]);
      scc_sl = stats::predict(lls,lx)$y;
      lxx = scc_sl;
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
      recr_sl[k,l] = 10^(min_yre-5);
    }
  }
}

for (k in (floor(Nx*0.5)+1):Nx){
  for (l in (floor(Nx*0.5)+1):Nx){
    if(k!=l){
      lx=abs(x[k]-x[l]);
      scc_sl = stats::predict(lls,lx)$y;
      lxx = scc_sl;
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
      recr_sl[k,l] = 10^(min_yre-5);
    }
  }
}

recp_sl_max=10^(min_yre-5);
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

for (k in 1:Nx){
  for (l in 1:Nx){
    recr_sl1[k,Nx-l+1]=recr_sl[k,l]
    if (recr_sl1[k,Nx-l+1]>0){
      recr_sl1_log[k,Nx-l+1] <-  log10(recr_sl1[k,Nx-l+1]);
    }else
      recr_sl1_log[k,Nx-l+1] <- -10;
  }
}

library("tagcloud")
library(reshape2)
library(dplyr)

recp_melt = recr_sl1 %>% melt() 
color_recp <- smoothPalette(recp_melt$value, pal="Reds");
png("cp_gen_xr2_sl_l.png", width=500, height=500);
plot(recp_melt$Var1,recp_melt$Var2,col=color_recp,pch=16, ann=F, asp = 1,mar = c(0,0,0,0))
dev.off()

recp_melt = recr_sl1_log %>% melt() 
color_recp <- smoothPalette(recp_melt$value, pal="Reds");
png("cp_gen_xr2_sl_log.png", width=500, height=500);
plot(recp_melt$Var1,recp_melt$Var2,col=color_recp,pch=16, ann=F, asp = 1,mar = c(0,0,0,0))
dev.off()

xr<-matrix(rep(0,((Nx-1)*Nx+Nx-1-Nx+1)*2),(Nx-1)*Nx+Nx-1-Nx+1,2);
i=1;
for (k in 1:Nx){
  for (l in 1:Nx){
    if (k!=l){
      xr[i,1]=abs(x[k]-x[l]);
      xr[i,2]=recr_sl[k,l];
      i=i+1;
    }
  }
}

colnames(xr) <- c("dis","cp");
xr <- as.data.frame(xr);
yr<-xr[!duplicated(xr[c("dis","cp")]),];
zr<-yr[2,];
for(i in 3:length(yr[,1])){
    if (yr[i,2]!=0){
      zr<-rbind(zr,yr[i,]);
    }
}

zr<-zr[order(zr[,1]),];

lzr<-log10(zr);
sx<-(zr[length(zr[,1]),1]-zr[1,1])/(30-1);
xxx<-seq(zr[1,1],zr[length(zr[,1]),1],sx);
xxl<-log10(xxx);
llx = stats::smooth.spline(lzr[,1],lzr[,2],spar = 0);
yyl<-stats::predict(llx,xxl)$y;
xr2_log <- cbind(xxl,yyl);
write.csv(xr2_log,"dis_data_xr2_sl_log.csv");
