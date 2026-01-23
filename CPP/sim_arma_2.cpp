//#include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]

// [[Rcpp::export]]
NumericVector rcpp_sim(int imax, double kk, int aloc1, int bloc1, int aloc2, int bloc2,
                       double epsilon1, double epsilon2, int Iter, double ssr, double ss, double ssq, int bound, double alpha, double beta,
                       double D, double gm, double dt, double kbt,
                       int upb, int lowb,
                       NumericVector x){
  int Nx=x.size();
  arma::mat xx=arma::zeros(imax,Nx);
  arma::mat yy1=arma::zeros(imax,Nx);
  arma::mat yy2=arma::zeros(imax,Nx);
  arma::mat xxx=arma::zeros(imax,Nx);
  arma::mat yyy1=arma::zeros(imax,Nx);
  arma::mat yyy2=arma::zeros(imax,Nx);
  arma::mat imx=arma::zeros(imax,Nx);
  arma::mat imy1=arma::zeros(imax,Nx);
  arma::mat imy2=arma::zeros(imax,Nx);
  arma::vec ax1=arma::zeros(imax); //initial location of loop extrusion occurring, i
  arma::vec bx1=arma::zeros(imax); //terminal location of loop extrusion occurring, j
  arma::vec axx1=arma::zeros(imax);
  arma::vec bxx1=arma::zeros(imax);
  arma::vec ax2=arma::zeros(imax); //initial location of loop extrusion occurring, i
  arma::vec bx2=arma::zeros(imax); //terminal location of loop extrusion occurring, j
  arma::vec axx2=arma::zeros(imax);
  arma::vec bxx2=arma::zeros(imax);
  arma::vec dphi=arma::zeros(Nx);
  arma::vec dphiy1=arma::zeros(Nx);
  arma::vec dphiy2=arma::zeros(Nx);
  arma::vec rr=arma::zeros(imax*Nx*Nx);
  arma::vec zz=arma::zeros((imax-1)*Nx+Nx+(imax-1)*Nx+Nx+(imax-1)*Nx+Nx+imax+imax+imax+imax+1);
  arma::vec aax1=arma::zeros(imax);
  arma::vec bbx1=arma::zeros(imax);
  arma::vec aax2=arma::zeros(imax);
  arma::vec bbx2=arma::zeros(imax);
  
  double xi;
  double ya1;
  double ya2;
  double nr;
  
  int az1;
  int bz1;
  int az2;
  int bz2;
  
  upb = upb-1;
  lowb = lowb-1;
  
  az1 = aloc1-1;
  bz1 = bloc1-1;
  az2 = aloc2-1;
  bz2 = bloc2-1;
  
  for (int j=0; j<Iter; j++){
    for(int n=0; n<Nx; n++){
      xx(0,n) = x[n];
      yy1(0,n) = x[n];
      yy2(0,n) = x[n];
      xxx(0,n)=xxx(0,n)+xx(0,n);
      yyy1(0,n)=yyy1(0,n)+yy1(0,n);
      yyy2(0,n)=yyy2(0,n)+yy2(0,n);
    }
    
    ax1[0] = az1;
    bx1[0] = bz1;
    ax2[0] = az2;
    bx2[0] = bz2;
    axx1[0]=axx1[0]+ax1[0];
    bxx1[0]=bxx1[0]+bx1[0];
    axx2[0]=axx2[0]+ax2[0];
    bxx2[0]=bxx2[0]+bx2[0];
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
      ya1 = abs(yy1(i-1,bx1[i-1])-yy1(i-1,ax1[i-1]));
      dphiy1[ax1[i-1]]=kk*abs(ya1-abs(ax1[i-1]-bx1[i-1])*ss);
      dphiy1[bx1[i-1]]=kk*abs(ya1-abs(ax1[i-1]-bx1[i-1])*ss);
      if (ya1 < ssq){
        dphiy1[ax1[i-1]]=dphiy1[ax1[i-1]]+epsilon2*(4*kbt*(-12*pow(ssr,12)/pow(ya1,13)+6*pow(ssr,6)/pow(ya1,7)));
        dphiy1[bx1[i-1]]=dphiy1[bx1[i-1]]+epsilon2*(4*kbt*(-12*pow(ssr,12)/pow(ya1,13)+6*pow(ssr,6)/pow(ya1,7)));
        if (ax1[i-1] > az1 - bound && bx1[i-1] == bz1){
          ax1[i] = ax1[i-1]-1;
          bx1[i] = bx1[i-1];
        }else if (ax1[i-1] == az1 - bound && bx1[i-1] > bz1 - bound){
            bx1[i]=bx1[i-1]-1;
            ax1[i]=ax1[i-1];
        }
      }else{
        ax1[i] = ax1[i-1];
        bx1[i] = bx1[i-1];
      }
      if (ax1[i]==bx1[i]){
        ax1[i]=ax1[i-1];
        bx1[i]=bx1[i-1];
      }
      if (ax1[i]<lowb) {ax1[i] = lowb;}
      if (ax1[i]>upb) {ax1[i] = upb;}
      if (bx1[i]<lowb) {bx1[i] = lowb;}
      if (bx1[i]>upb) {bx1[i] = upb;}
      
      ya2 = abs(yy2(i-1,bx2[i-1])-yy2(i-1,ax2[i-1]));
      dphiy2[ax2[i-1]]=kk*abs(ya2-abs(ax2[i-1]-bx2[i-1])*ss);
      dphiy2[bx2[i-1]]=kk*abs(ya2-abs(ax2[i-1]-bx2[i-1])*ss);
      if (ya2 < ssq){
        dphiy2[ax2[i-1]]=dphiy2[ax2[i-1]]+epsilon2*(4*kbt*(-12*pow(ssr,12)/pow(ya2,13)+6*pow(ssr,6)/pow(ya2,7)));
        dphiy2[bx2[i-1]]=dphiy2[bx2[i-1]]+epsilon2*(4*kbt*(-12*pow(ssr,12)/pow(ya2,13)+6*pow(ssr,6)/pow(ya2,7)));
        if (ax2[i-1] == az2 && bx2[i-1] < bz2 + bound){
          ax2[i] = ax2[i-1];
          bx2[i] = bx2[i-1]+1;
        }else if (ax2[i-1] < az2 + bound && bx2[i-1] == bz2 + bound){
          ax2[i] = ax2[i-1]+1;
          bx2[i] = bx2[i-1];
        }
      }else{
        ax2[i] = ax2[i-1];
        bx2[i] = bx2[i-1];
      }
      if (ax2[i]==bx2[i]){
        ax2[i]=ax2[i-1];
        bx2[i]=bx2[i-1];
      }
      if (ax2[i]<lowb) {ax2[i] = lowb;}
      if (ax2[i]>upb) {ax2[i] = upb;}
      if (bx2[i]<lowb) {bx2[i] = lowb;}
      if (bx2[i]>upb) {bx2[i] = upb;}      
      
      for (int n=0; n<Nx; n++){
        nr = arma::randn();
        xx(i,n) = xx(i-1,n)-1/gm*dphi[n]*dt+D*pow(dt,0.5)*nr;  
        if (n != ax1[i] && n != bx1[i]){
          yy1(i,n) = xx(i,n);
        }
        if (n != ax2[i] && n != bx2[i]){
          yy2(i,n) = xx(i,n);
        }
        if (n==ax1[i]){
          yy1(i,ax1[i]) = yy1(i-1,ax1[i-1])-1/gm*dphiy1[ax1[i-1]]*dt+D*pow(dt,0.5)*nr;  
        }
        if (n==bx1[i]){
          yy1(i,bx1[i]) = yy1(i-1,bx1[i-1])-1/gm*dphiy1[bx1[i-1]]*dt+D*pow(dt,0.5)*nr;  
        }
        if (n==ax2[i]){
          yy2(i,ax2[i]) = yy2(i-1,ax2[i-1])-1/gm*dphiy2[ax2[i-1]]*dt+D*pow(dt,0.5)*nr;  
        }
        if (n==bx2[i]){
          yy2(i,bx2[i]) = yy2(i-1,bx2[i-1])-1/gm*dphiy2[bx2[i-1]]*dt+D*pow(dt,0.5)*nr;  
        }
      }
      for (int n=0; n<Nx; n++){
        xxx(i,n)=xxx(i,n)+xx(i,n);
        yyy1(i,n)=yyy1(i,n)+yy1(i,n);
        yyy2(i,n)=yyy2(i,n)+yy2(i,n);
      }
      axx1[i]=axx1[i]+ax1[i];
      bxx1[i]=bxx1[i]+bx1[i];
      axx2[i]=axx2[i]+ax2[i];
      bxx2[i]=bxx2[i]+bx2[i];
    }
  }

  for(int i=0; i<imax; i++){
    for (int n=0; n<Nx; n++){
      imx(i,n)=xxx(i,n)/Iter;
      imy1(i,n)=yyy1(i,n)/Iter;
      imy2(i,n)=yyy2(i,n)/Iter;
    }
    aax1[i]=round(axx1[i]/Iter);
    bbx1[i]=round(bxx1[i]/Iter);
    if (aax1[i]>Nx-1){
      aax1[i] = Nx-1;
    }
    if (bbx1[i]>Nx-1){
      bbx1[i] = Nx-1;
    }
    aax2[i]=round(axx2[i]/Iter);
    bbx2[i]=round(bxx2[i]/Iter);
    if (aax2[i]>Nx-1){
      aax2[i] = Nx-1;
    }
    if (bbx2[i]>Nx-1){
      bbx2[i] = Nx-1;
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
    zz[(imax-1)*Nx+Nx+(imax-1)*Nx+Nx+(imax-1)*Nx+Nx+i]=aax1[i];
    zz[(imax-1)*Nx+Nx+(imax-1)*Nx+Nx+(imax-1)*Nx+Nx+imax+i]=bbx1[i];
    zz[(imax-1)*Nx+Nx+(imax-1)*Nx+Nx+(imax-1)*Nx+Nx+imax+imax+i]=aax2[i];
    zz[(imax-1)*Nx+Nx+(imax-1)*Nx+Nx+(imax-1)*Nx+Nx+imax+imax+imax+i]=bbx2[i];
  }

  Rcpp::NumericVector zzz = NumericVector(zz.begin(), zz.end());

  return zzz;
}

