# one point loop extrusion simulation (moving)
# "image magick" have to be installed

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

dat_1 <- t(read.csv(paste(paste("yre_df_", fname, sep=""),"_l.csv", sep =""), header = TRUE));
yre <- as.numeric(dat_1[2:length(dat_1[,2]),2]);

ssr <- 21.8/0.3;
ss <- 5000;
kk <- kbt/ss^2;
Nx <- length(x);

ssq <- 40*ss;

Iter <- 10000;
imax <- 100;
alpha <- 1-1/beta;
aloc=14;
bloc=16;
bound <- round(Nx*0.4,0);
upb <- round(Nx*0.9,0);
lowb <- round(Nx*0.1,0);

gm <- 2;
D <- (2*kbt*gm)^0.5/gm;
t_lj <- (ss^2/kbt)^alpha;

epsilon1 <- 0.01;
epsilon2 <- 0.01;

tmax <- dt*imax;
print(tmax);

diff_alpha <- dt^(alpha/2);
print(diff_alpha);

aax <- rep(0,imax);
bbx <- rep(0,imax);

print(Sys.time());
st1 <- Sys.time();

library(Rcpp);
library (inline);
library (RcppArmadillo);
sourceCpp('sim_arma_1.cpp');
recpx <- rcpp_sim(imax, kk, aloc, bloc, bound, Iter, epsilon1, epsilon2, ssr, ss, ssq, alpha, beta, D, gm, dt, kbt, upb, lowb, x);

imx <- array(rep(0,imax*Nx),dim=c(imax,Nx));
imy <- array(rep(0,imax*Nx),dim=c(imax,Nx));

for (i in 0:(imax-1)){
  for (n in 0:(Nx-1)){
    imx[i+1,n+1] <- recpx[i*Nx+n+1];
    imy[i+1,n+1] <- recpx[(imax-1)*Nx+Nx+i*Nx+n+1];
  }
}

for (i in 0:(imax-1)){
  aax[i+1] <- recpx[(imax-1)*Nx+Nx+(imax-1)*Nx+Nx+i+1]+1;
  bbx[i+1] <- recpx[(imax-1)*Nx+Nx+(imax-1)*Nx+Nx+imax+i+1]+1;
}

rate=1000;

kkk <- rep(0,rate);
kx <- rep(0,rate-20);
rx <- rep(0,rate-20);
rrr <- rep(0,rate);
recr <- array(rep(0,imax*Nx*Nx),dim=c(imax,Nx,Nx));
recp <- array(rep(0,imax*Nx*Nx),dim=c(imax,Nx,Nx));

min_yre=min(yre);

lls = stats::smooth.spline(x-x[1],lsc*(x[2]-x[1]),spar = 0);

library(reticulate);
use_python("C:/Users/Masamichi Sato/PycharmProjects/RG/venv/Scripts/")

for (i in 1:imax){
  for (k in 1:Nx){
    for (l in 1:Nx){
      if(k!=l){
        if (k==aax[i] && l==bbx[i]){
          lx=abs(imy[i,k]-imy[i,l]);
        }else if (k==bbx[i] && l==aax[i]){
          lx=abs(imy[i,k]-imy[i,l]);
        }else{
          lx=abs(imx[i,k]-imx[i,l]);
        }
        scc = stats::predict(lls,lx)$y;#mean(lsc);
        lxx = scc;
        r1 = lx - lxx;
        r2 = lx + lxx;
        del_r = (r2 - r1) / (rate-1);  
        numpy <- import("numpy");
        kkk <- numpy$fft$fftfreq(length(rrr), del_r);
        kx <- kkk[11:(rate-10)]
        ggg = 1 / (4 * kk * ((kx / lxx^nsc) ^ 2 / 2 - (kx / lxx^nsc) ^ 4 / 24)) * (sin(kx / 2) / sin(kx / lxx^nsc / 2)) ^ beta;
        gy = numpy$fft$ifft(ggg);
        gyc= 1 / lxx^((2 * theta + 1)*nsc) * Re(gy);
        recr[i,k,l] = gyc[floor(rate/2)-10];
      }else{
        recr[i,k,l] = 10^(min_yre-5);
      }
    }
  }
}


recp_max = max(recr);
recp_max=max(c(10^(-10),recp_max))

for (i in 1:imax){
  diag(recr[i,,])=rep(recp_max,length(recr[1,1,]));
}

for (i in 1:imax){
  for (k in 1:Nx){
    for (l in 1:Nx){
      recp[i,k,Nx-l+1]=recr[i,k,l]
    }
  }
}

recb <- log10(recp);

print(Sys.time());
st2 <- Sys.time();

print(st2-st1);

library("tagcloud")
library(reshape2)
library(dplyr)
library("spatstat")
library("plotrix")

for (i in 1:imax){
  recp_melt = recp[i,,] %>% melt() 
  color_recp <- smoothPalette(recp_melt$value, pal="Reds");
  png(sprintf("img/recpax%03d.png",i), width=500, height=500);
  plot(recp_melt$Var1,recp_melt$Var2,col=color_recp,pch=16, ann=F, asp = 1,mar = c(0,0,0,0))
  x <- recp_melt$value;
  n.colors <- 5;
  # interpolate colours:
  cc <- smoothPalette(sort(x, decreasing = F), pal="Reds");
  palette <- colorRampPalette(cc, space = "rgb")(n.colors)
  color.map <- colourmap( palette , range=range(x) )
  color.range <- color.map( seq(min(x), max(x), length.out = n.colors) )
  #the labels of the legend
  col.labels <- round(seq(min(x),max(x),length=5), 4);
  color.legend( xl =50 , yb = 40, xr = 52, yt = 50 , # the coordinates
                legend = col.labels , gradient="y", 
                rect.col=color.range, align="rb")
  dev.off()
}
command <- sprintf("magick convert -delay 10 -loop 1 img/recpax*.png movierecpax.gif");
system(command);


for (i in 1:imax){
  recb_melt = recb[i,,] %>% melt() 
  color_recp <- smoothPalette(recb_melt$value, pal="Reds");
  png(sprintf("img/recpbx%03d.png",i), width=500, height=500);
  plot(recb_melt$Var1,recp_melt$Var2,col=color_recp,pch=16, ann=F, asp = 1,mar = c(0,0,0,0))
  x <- recb_melt$value;
  n.colors <- 5;
  # interpolate colours:
  cc <- smoothPalette(sort(x, decreasing = F), pal="Reds");
  palette <- colorRampPalette(cc, space = "rgb")(n.colors)
  color.map <- colourmap( palette , range=range(x) )
  color.range <- color.map( seq(min(x), max(x), length.out = n.colors) )
  #the labels of the legend
  col.labels <- round(seq(min(x),max(x),length=5), 4);
  color.legend( xl =50 , yb = 40, xr = 52, yt = 50 , # the coordinates
                legend = col.labels , gradient="y", 
                rect.col=color.range, align="rb")
  dev.off()
}
command <- sprintf("magick convert -delay 10 -loop 1 img/recpbx*.png movierecpbx.gif");
system(command);

