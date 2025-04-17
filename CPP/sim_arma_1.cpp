//#include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]

// [[Rcpp::export]]
NumericVector rcpp_sim(int imax, double kk, int aloc, int bloc, int bound,
                       int Iter, double epsilon1, double epsilon2, double ssr, double ss, double ssq, double alpha, double beta,
                       double D, double gm, double dt, double kbt,
                       int upb, int lowb,
                       NumericVector x){
  int Nx=x.size();
  arma::mat xx=arma::zeros(imax,Nx);
  arma::mat yy=arma::zeros(imax,Nx);
  arma::mat xxx=arma::zeros(imax,Nx);
  arma::mat yyy=arma::zeros(imax,Nx);
  arma::mat imx=arma::zeros(imax,Nx);
  arma::mat imy=arma::zeros(imax,Nx);
  arma::vec ax=arma::zeros(imax); //initial location of loop extrusion occurring, i
  arma::vec bx=arma::zeros(imax); //terminal location of loop extrusion occurring, j
  arma::vec axx=arma::zeros(imax);
  arma::vec bxx=arma::zeros(imax);
  arma::cube recp=arma::zeros(imax,Nx,Nx);
  arma::vec dphi=arma::zeros(Nx);
  arma::vec dphiy=arma::zeros(Nx);
  arma::vec rr=arma::zeros(imax*Nx*Nx);
  arma::vec zz=arma::zeros((imax-1)*Nx+Nx+(imax-1)*Nx+Nx+imax+imax+1);
  arma::vec aax=arma::zeros(imax);
  arma::vec bbx=arma::zeros(imax);
  
  double xi;
  double ya;
  double nr;
  
  int az;
  int bz;
  
  upb = upb-1;
  lowb = lowb-1;
  
  az = aloc-1;
  bz = bloc-1;
  
  for (int j=0; j<Iter; j++){
    for(int n=0; n<Nx; n++){
      xx(0,n) = x[n];
      yy(0,n) = x[n];
      xxx(0,n)=xxx(0,n)+xx(0,n);
      yyy(0,n)=yyy(0,n)+yy(0,n);
    }
    
    ax[0] = az;
    bx[0] = bz;
    axx[0]=axx[0]+ax[0];
    bxx[0]=bxx[0]+bx[0];
    for (int i=1; i<imax; i++){
      for (int l=0; l<Nx; l++){
        if (l != Nx-1){
          xi = abs(xx(i-1,l+1)-xx(i-1,l));
          dphi[l]=kk*abs(xi-ssr);
          if (xi < pow(2,1/6) * ssr){
            dphi[l]=dphi[l]+epsilon1*(4*kbt*(-12*pow(ssr,12)/pow(xi,13)+6*pow(ssr,6)/pow(xi,7)));
          }
        }
      }
      ya = abs(yy(i-1,bx[i-1])-yy(i-1,ax[i-1]));
      dphiy[ax[i-1]]=kk*abs(ya-abs(ax[i-1]-bx[i-1])*ss);
      dphiy[bx[i-1]]=kk*abs(ya-abs(ax[i-1]-bx[i-1])*ss);
      if (ya < ssq){
        dphiy[ax[i-1]]=dphiy[ax[i-1]]+epsilon2*(4*kbt*(-12*pow(ssr,12)/pow(ya,13)+6*pow(ssr,6)/pow(ya,7)));
        dphiy[bx[i-1]]=dphiy[bx[i-1]]+epsilon2*(4*kbt*(-12*pow(ssr,12)/pow(ya,13)+6*pow(ssr,6)/pow(ya,7)));
        if (ax[i-1] == az && bx[i-1] < bz + bound + 1){
          ax[i] = ax[i-1];
          bx[i] = bx[i-1]+1;
        }
        if (ax[i-1] < az + bound && bx[i-1] == bz + bound + 1){
          ax[i] = ax[i-1]+1;
          bx[i] = bx[i-1];
        }
      }else{
        ax[i] = ax[i-1];
        bx[i] = bx[i-1];
      }
      if (ax[i]==bx[i]){
        ax[i]=ax[i-1];
        bx[i]=bx[i-1];
      }
      if (ax[i]<lowb) {ax[i] = lowb;}
      if (ax[i]>upb) {ax[i] = upb;}
      if (bx[i]<lowb) {bx[i] = lowb;}
      if (bx[i]>upb) {bx[i] = upb;}
      for (int n=0; n<Nx; n++){
        nr = arma::randn();
        xx(i,n) = xx(i-1,n)-1/gm*dphi[n]*dt+D*pow(dt,0.5)*nr;  
        if (n != ax[i] && n != bx[i]){
          yy(i,n) = xx(i,n);
        }
        if (n==ax[i]){
          yy(i,ax[i]) = yy(i-1,ax[i-1])-1/gm*dphiy[ax[i-1]]*dt+D*pow(dt,0.5)*nr;  
        }
        if (n==bx[i]){
          yy(i,bx[i]) = yy(i-1,bx[i-1])-1/gm*dphiy[bx[i-1]]*dt+D*pow(dt,0.5)*nr;  
        }
      }
      for (int n=0; n<Nx; n++){
        xxx(i,n)=xxx(i,n)+xx(i,n);
        yyy(i,n)=yyy(i,n)+yy(i,n);
      }
      axx[i]=axx[i]+ax[i];
      bxx[i]=bxx[i]+bx[i];
    }
  }

  for(int i=0; i<imax; i++){
    for (int n=0; n<Nx; n++){
      imx(i,n)=xxx(i,n)/Iter;
      imy(i,n)=yyy(i,n)/Iter;
    }
    aax[i]=round(axx[i]/Iter);
    bbx[i]=round(bxx[i]/Iter);
    if (aax[i]>Nx-1){
      aax[i] = Nx-1;
    }
    if (bbx[i]>Nx-1){
      bbx[i] = Nx-1;
    }
  }
  
  for(int i=0; i<imax; i++){
    for (int n=0; n<Nx; n++){
      zz[i*Nx+n]=imx(i,n);
      zz[(imax-1)*Nx+Nx+i*Nx+n]=imy(i,n);
    }
  }
      
  for(int i=0; i<imax; i++){
    zz[(imax-1)*Nx+Nx+(imax-1)*Nx+Nx+i]=aax[i];
  }
      
  for(int i=0; i<imax; i++){
    zz[(imax-1)*Nx+Nx+(imax-1)*Nx+Nx+imax+i]=bbx[i];
  }
      
  Rcpp::NumericVector zzz = NumericVector(zz.begin(), zz.end());

  return zzz;
}
