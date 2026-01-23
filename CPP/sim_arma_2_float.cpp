//#include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]

// [[Rcpp::export]]
NumericVector rcpp_sim(int imax, double kk, int aloc1, int bloc1, int aloc2, int bloc2,
                       double ploc1, double ploc2, double ploc3, double ploc4,
                       double epsilon1, double epsilon2, double ssr, double ss, double ssq, double alpha,double beta,
                       double D, double gm, double dt, double kbt,
                       int bound1, int bound2,
                       NumericVector x){
  int Nx=x.size();
  arma::mat xx=arma::zeros(imax,Nx);
  arma::mat yy1=arma::zeros(imax,Nx);
  arma::mat yy2=arma::zeros(imax,Nx);
  arma::mat imx=arma::zeros(imax,Nx);
  arma::mat imy1=arma::zeros(imax,Nx);
  arma::mat imy2=arma::zeros(imax,Nx);
  arma::vec ax1=arma::zeros(imax); //initial location of loop extrusion occurring, i
  arma::vec bx1=arma::zeros(imax); //terminal location of loop extrusion occurring, j
  arma::vec ax2=arma::zeros(imax); //initial location of loop extrusion occurring, i
  arma::vec bx2=arma::zeros(imax); //terminal location of loop extrusion occurring, j
  arma::vec dphi=arma::zeros(Nx);
  arma::vec dphiy1=arma::zeros(Nx);
  arma::vec dphiy2=arma::zeros(Nx);
  arma::vec rr=arma::zeros(imax*Nx*Nx);
  arma::vec zz=arma::zeros((imax-1)*Nx+Nx+(imax-1)*Nx+Nx+(imax-1)*Nx+Nx+imax+imax+imax+imax+1);

  double xi;
  double ya1;
  double ya2;
  double nr1;
  double nr2;
  
  int az1;
  int bz1;
  int az2;
  int bz2;
  
  int upb1;
  int lowb1;
  int upb2;
  int lowb2;
  
  double pp;
  
  upb1 = bloc1+bound1-1;
  lowb1 = aloc1-bound1-1;
  upb2 = bloc2+bound2-1;
  lowb2 = aloc2-bound2-1;

  az1 = aloc1-1;
  bz1 = bloc1-1;
  az2 = aloc2-1;
  bz2 = bloc2-1;
  
  for(int n=0; n<Nx; n++){
    xx(0,n) = x[n];
    yy1(0,n) = x[n];
    yy2(0,n) = x[n];
    imx(0,n)=imx(0,n)+xx(0,n);
    imy1(0,n)=imy1(0,n)+yy1(0,n);
    imy2(0,n)=imy2(0,n)+yy2(0,n);
  }
  
  ax1[0] = az1;
  bx1[0] = bz1;
  ax2[0] = az2;
  bx2[0] = bz2;
  for (int i=1; i<imax; i++){
    for (int l=0; l<Nx; l++){
      if (l != Nx-1){
        xi = abs(xx(i-1,l+1)-xx(i-1,l));
        dphi[l] = kk*abs(xi-ssr);
        if (xi < pow(2.,(1/6)) * ssr){
          dphi[l]=dphi[l]+epsilon1*(4.*kbt*(-12*pow(ssr,12)/pow(xi,13)+6*pow(ssr,6)/pow(xi,7))+30*kbt*xi/(1-pow(xi/1.6/ssr,2));
        }
      }
    }
    if (bx1[i-1]!=ax1[i-1]){
      ya1 = abs(yy1(i-1,bx1[i-1])-yy1(i-1,ax1[i-1]));
      dphiy1[ax1[i-1]]=kk*abs(ya1-abs(ax1[i-1]-bx1[i-1])*ss);
      dphiy1[bx1[i-1]]=kk*abs(ya1-abs(ax1[i-1]-bx1[i-1])*ss);
      pp = arma::randu();
      if (pp >= 0 && pp < ploc1 && ya1 < ssq){
        dphiy1[ax1[i-1]]=dphiy1[ax1[i-1]]+epsilon2*(4*kbt*(-12*pow(ssr,12)/pow(ya1,13)+6*pow(ssr,6)/pow(ya1,7)));
        dphiy1[bx1[i-1]]=dphiy1[bx1[i-1]]+epsilon2*(4*kbt*(-12*pow(ssr,12)/pow(ya1,13)+6*pow(ssr,6)/pow(ya1,7)));
        ax1[i] = ax1[i-1]+1;
        bx1[i] = bx1[i-1];
        if (ax1[i]>upb1) {ax1[i] = upb1;}
      }else if (pp >= ploc1 && pp < ploc2 && ya1 < ssq){
        dphiy1[ax1[i-1]]=dphiy1[ax1[i-1]]+epsilon2*(4*kbt*(-12*pow(ssr,12)/pow(ya1,13)+6*pow(ssr,6)/pow(ya1,7)));
        dphiy1[bx1[i-1]]=dphiy1[bx1[i-1]]+epsilon2*(4*kbt*(-12*pow(ssr,12)/pow(ya1,13)+6*pow(ssr,6)/pow(ya1,7)));
        ax1[i] = ax1[i-1]-1;
        bx1[i] = bx1[i-1];
        if (ax1[i]<lowb1) {ax1[i] = lowb1;}
      }else if (pp >= ploc2 && pp < ploc3 && ya1 < ssq){
        dphiy1[ax1[i-1]]=dphiy1[ax1[i-1]]+epsilon2*(4*kbt*(-12*pow(ssr,12)/pow(ya1,13)+6*pow(ssr,6)/pow(ya1,7)));
        dphiy1[bx1[i-1]]=dphiy1[bx1[i-1]]+epsilon2*(4*kbt*(-12*pow(ssr,12)/pow(ya1,13)+6*pow(ssr,6)/pow(ya1,7)));
        ax1[i] = ax1[i-1];
        bx1[i] = bx1[i-1]+1;
        if (bx1[i]>upb1) {bx1[i] = upb1;}
      }else if (pp >= ploc3 && pp < ploc4 && ya1 < ssq){
        dphiy1[ax1[i-1]]=dphiy1[ax1[i-1]]+epsilon2*(4*kbt*(-12*pow(ssr,12)/pow(ya1,13)+6*pow(ssr,6)/pow(ya1,7)));
        dphiy1[bx1[i-1]]=dphiy1[bx1[i-1]]+epsilon2*(4*kbt*(-12*pow(ssr,12)/pow(ya1,13)+6*pow(ssr,6)/pow(ya1,7)));
        ax1[i] = ax1[i-1];
        bx1[i] = bx1[i-1]-1;
        if (bx1[i]<lowb1) {bx1[i] = lowb1;}
      }else if (pp >= ploc4 && pp < 1 && ya1 < ssq){
        dphiy1[ax1[i-1]]=dphiy1[ax1[i-1]]+epsilon2*(4.*kbt*(-12*pow(ssr,12)/pow(ya1,13)+6*pow(ssr,6)/pow(ya1,7)));
        dphiy1[bx1[i-1]]=dphiy1[bx1[i-1]]+epsilon2*(4.*kbt*(-12*pow(ssr,12)/pow(ya1,13)+6*pow(ssr,6)/pow(ya1,7)));
        ax1[i] = ax1[i-1];
        bx1[i] = bx1[i-1];
      }else{
        dphiy1[ax1[i-1]]=dphi[ax1[i-1]];
        dphiy1[bx1[i-1]]=dphi[bx1[i-1]];
        ax1[i] = ax1[i-1];
        bx1[i] = bx1[i-1];
      }
    }else{
      dphiy1[ax1[i-1]]=dphiy1[ax1[i-2]];
      dphiy1[bx1[i-1]]=dphiy1[bx1[i-2]];
      ax1[i] = ax1[i-2];
      bx1[i] = bx1[i-2];
    }      
    if (ax1[i]==bx1[i]){
      ax1[i]=ax1[i-1];
      bx1[i]=bx1[i-1];
    }
    if (bx2[i-1]!=ax2[i-1]){
      ya2 = abs(yy2(i-1,bx2[i-1])-yy2(i-1,ax2[i-1]));
      pp = arma::randu();
      dphiy2[ax2[i-1]]=kk*abs(ya2-abs(ax2[i-1]-bx2[i-1])*ss);
      dphiy2[bx2[i-1]]=kk*abs(ya2-abs(ax2[i-1]-bx2[i-1])*ss);
      if (pp >= 0 && pp < ploc1 && ya2 < ssq){
        dphiy2[ax2[i-1]]=dphiy2[ax2[i-1]]+epsilon2*(4.*kbt*(-12*pow(ssr,12)/pow(ya2,13)+6*pow(ssr,6)/pow(ya2,7)));
        dphiy2[bx2[i-1]]=dphiy2[bx2[i-1]]+epsilon2*(4.*kbt*(-12*pow(ssr,12)/pow(ya2,13)+6*pow(ssr,6)/pow(ya2,7)));
        ax2[i] = ax2[i-1]+1;
        bx2[i] = bx2[i-1];
        if (ax2[i]>upb2) {ax2[i] = upb2;}
      }else if (pp >= ploc1 && pp < ploc2 && ya2 < ssq){
        dphiy2[ax2[i-1]]=dphiy2[ax2[i-1]]+epsilon2*(4.*kbt*(-12*pow(ssr,12)/pow(ya2,13)+6*pow(ssr,6)/pow(ya2,7)));
        dphiy2[bx2[i-1]]=dphiy2[bx2[i-1]]+epsilon2*(4.*kbt*(-12*pow(ssr,12)/pow(ya2,13)+6*pow(ssr,6)/pow(ya2,7)));
        ax2[i] = ax2[i-1]-1;
        bx2[i] = bx2[i-1];
        if (ax2[i]<lowb2) {ax2[i] = lowb2;}
      }else if (pp >= ploc2 && pp < ploc3 && ya2 < ssq){
        dphiy2[ax2[i-1]]=dphiy2[ax2[i-1]]+epsilon2*(4.*kbt*(-12*pow(ssr,12)/pow(ya2,13)+6*pow(ssr,6)/pow(ya2,7)));
        dphiy2[bx2[i-1]]=dphiy2[bx2[i-1]]+epsilon2*(4.*kbt*(-12*pow(ssr,12)/pow(ya2,13)+6*pow(ssr,6)/pow(ya2,7)));
        ax2[i] = ax2[i-1];
        bx2[i] = bx2[i-1]+1;
        if (bx2[i]>upb2) {bx2[i] = upb2;}
      }else if (pp >= ploc3 && pp < ploc4 && ya2 < ssq){
        dphiy2[ax2[i-1]]=dphiy2[ax2[i-1]]+epsilon2*(4.*kbt*(-12*pow(ssr,12)/pow(ya2,13)+6*pow(ssr,6)/pow(ya2,7)));
        dphiy2[bx2[i-1]]=dphiy2[bx2[i-1]]+epsilon2*(4.*kbt*(-12*pow(ssr,12)/pow(ya2,13)+6*pow(ssr,6)/pow(ya2,7)));
        ax2[i] = ax2[i-1];
        bx2[i] = bx2[i-1]-1;
        if (bx2[i]<lowb2) {bx2[i] = lowb2;}
      }else if (pp >= ploc4 && pp < 1 && ya2 < ssq){
        dphiy2[ax2[i-1]]=dphiy2[ax2[i-1]]+epsilon2*(4.*kbt*(-12*pow(ssr,12)/pow(ya2,13)+6*pow(ssr,6)/pow(ya2,7)));
        dphiy2[bx2[i-1]]=dphiy2[bx2[i-1]]+epsilon2*(4.*kbt*(-12*pow(ssr,12)/pow(ya2,13)+6*pow(ssr,6)/pow(ya2,7)));
        ax2[i] = ax2[i-1];
        bx2[i] = bx2[i-1];
      }else{
        dphiy2[ax2[i-1]]=dphi[ax2[i-1]];
        dphiy2[bx2[i-1]]=dphi[bx2[i-1]];
        ax2[i] = ax2[i-1];
        bx2[i] = bx2[i-1];
      }
    }else{
      dphiy2[ax2[i-1]]=dphiy2[ax2[i-2]];
      dphiy2[bx2[i-1]]=dphiy2[bx2[i-2]];
      ax2[i] = ax2[i-2];
      bx2[i] = bx2[i-2];
    }      
    if (ax2[i]==bx2[i]){
      ax2[i]=ax2[i-1];
      bx2[i]=bx2[i-1];
    }
    if (ax1[i]>Nx-1){
      ax1[i] = Nx-1;
    }
    if (bx1[i]>Nx-1){
      bx1[i] = Nx-1;
    }
    if (ax2[i]>Nx-1){
      ax2[i] = Nx-1;
    }
    if (bx2[i]>Nx-1){
      bx2[i] = Nx-1;
    }
    if (ax1[i]<1){
      ax1[i] = 1;
    }
    if (bx1[i]<1){
      bx1[i] = 1;
    }
    if (ax2[i]<1){
      ax2[i] = 1;
    }
    if (bx2[i]<1){
      bx2[i] = 1;
    }
    
    for (int n=0; n<Nx; n++){
      nr1 = arma::randn();
      nr2 = arma::randn();
      xx(i,n) = xx(i-1,n)-1/gm*dphi[n]*dt+D*pow(dt,0.5)*nr1;  
      if (n != ax1[i] && n != bx1[i]){
        yy1(i,n) = xx(i,n);
      }
      if (n != ax2[i] && n != bx2[i]){
        yy2(i,n) = xx(i,n);
      }
      if (ax1[i-1]!=bx1[i-1]){
        if (n==ax1[i]){
          yy1(i,ax1[i]) = yy1(i-1,ax1[i-1])-1/gm*dphiy1[ax1[i-1]]*dt+D*pow(dt,0.5)*nr1;  
        }
        if (n==bx1[i]){
          yy1(i,bx1[i]) = yy1(i-1,bx1[i-1])-1/gm*dphiy1[bx1[i-1]]*dt+D*pow(dt,0.5)*nr1;  
        }
      }else{
        yy1(i,ax1[i]) = yy1(i-2,ax1[i-2]);
        yy1(i,bx1[i]) = yy1(i-2,bx1[i-2]);
      }
      if (ax2[i-1]!=bx2[i-1]){
        if (n==ax2[i]){
          yy2(i,ax2[i]) = yy2(i-1,ax2[i-1])-1/gm*dphiy2[ax2[i-1]]*dt+D*pow(dt,0.5)*nr2;  
        }
        if (n==bx2[i]){
          yy2(i,bx2[i]) = yy2(i-1,bx2[i-1])-1/gm*dphiy2[bx2[i-1]]*dt+D*pow(dt,0.5)*nr2;  
        }
      }else{
        yy2(i,ax2[i]) = yy2(i-2,ax2[i-2]);
        yy2(i,bx2[i]) = yy2(i-2,bx2[i-2]);
      }
    }
    for (int n=0; n<Nx; n++){
      imx(i,n)=imx(i,n)+xx(i,n);
      imy1(i,n)=imy1(i,n)+yy1(i,n);
      imy2(i,n)=imy2(i,n)+yy2(i,n);
    }
  }

  for(int i=0; i<imax; i++){
    for (int n=0; n<Nx; n++){
      zz[i*Nx+n]=imx(i,n);
      zz[(imax-1)*Nx+Nx+i*Nx+n]=imy1(i,n);
      zz[(imax-1)*Nx+Nx+(imax-1)*Nx+Nx+i*Nx+n]=imy2(i,n);
    }
  }
      
  for(int i=0; i<imax; i++){
    zz[(imax-1)*Nx+Nx+(imax-1)*Nx+Nx+(imax-1)*Nx+Nx+i]=ax1[i];
    zz[(imax-1)*Nx+Nx+(imax-1)*Nx+Nx+(imax-1)*Nx+Nx+imax+i]=bx1[i];
    zz[(imax-1)*Nx+Nx+(imax-1)*Nx+Nx+(imax-1)*Nx+Nx+imax+imax+i]=ax2[i];
    zz[(imax-1)*Nx+Nx+(imax-1)*Nx+Nx+(imax-1)*Nx+Nx+imax+imax+imax+i]=bx2[i];
  }

  Rcpp::NumericVector zzz = NumericVector(zz.begin(), zz.end());

  return zzz;
}

