//#include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]

// [[Rcpp::export]]
NumericVector rcpp_sim(int imax, double kk,
                       double ploc1, double ploc2, double ploc3, double ploc4,
                       int aloc, int bloc,
                       double epsilon1, double epsilon2,int Iter,double ssr, double ss, double ssq,double alpha,double beta,
                       double D, double gm, double dt, double kbt,
                       int upb, int lowb,
                        NumericVector x){

  int Nx=x.size();
  arma::mat xx=arma::zeros(imax,Nx);
  arma::mat yy=arma::zeros(imax,Nx);
  arma::vec ax=arma::zeros(imax); //initial location of loop extrusion occurring, i
  arma::vec bx=arma::zeros(imax); //terminal location of loop extrusion occurring, j
  arma::cube recp=arma::zeros(imax,Nx,Nx);
  arma::vec dphi=arma::zeros(Nx);
  arma::vec dphiy=arma::zeros(Nx);
  arma::vec rr=arma::zeros(imax*Nx*Nx);
  arma::vec zz=arma::zeros((imax-1)*Nx+Nx+(imax-1)*Nx+Nx+imax+imax+1);

  double xi;
  double ya;
  double pp;
  double nr;
  
  int az;
  int bz;
  
  upb = upb-1;
  lowb = lowb-1;
  
  az = aloc-1;
  bz = bloc-1;
  
  for(int n=0; n<Nx; n++){
    xx(0,n) = x[n];
    yy(0,n) = x[n];
  }
  
  ax[0] = az;
  bx[0] = bz;
  for (int i=1; i<imax; i++){
    for (int l=0; l<Nx; l++){
      if (l != Nx-1){
        xi = abs(xx(i-1,l+1)-xx(i-1,l));
        dphi[l] = kk*abs(xi-ssr);
        if (xi < pow(2.,(1/6)) * ssr){
          dphi[l]=dphi[l]+epsilon1*(4.*kbt*(-12*pow(ssr,12)/pow(xi,13)+6*pow(ssr,6)/pow(xi,7)));
        }
      }
    }
    ya = abs(yy(i-1,bx[i-1])-yy(i-1,ax[i-1]));
    dphiy[ax[i-1]]=kk*abs(ya-abs(ax[i-1]-bx[i-1])*ss);
    dphiy[bx[i-1]]=kk*abs(ya-abs(ax[i-1]-bx[i-1])*ss);
    pp = arma::randu();
    if (pp >= 0 && pp < ploc1 && ya < ssq){
      dphiy[ax[i-1]]=dphiy[ax[i-1]]+epsilon2*(4.*kbt*(-12*pow(ssr,12)/pow(ya,13)+6*pow(ssr,6)/pow(ya,7)));
      dphiy[bx[i-1]]=dphiy[bx[i-1]]+epsilon2*(4.*kbt*(-12*pow(ssr,12)/pow(ya,13)+6*pow(ssr,6)/pow(ya,7)));
      ax[i] = ax[i-1]+1;
      bx[i] = bx[i-1];
      if (ax[i]>upb) {ax[i] = upb;}
    }else if (pp >= ploc1 && pp < ploc2 && ya < ssq){
      dphiy[ax[i-1]]=dphiy[ax[i-1]]+epsilon2*(4.*kbt*(-12*pow(ssr,12)/pow(ya,13)+6*pow(ssr,6)/pow(ya,7)));
      dphiy[bx[i-1]]=dphiy[bx[i-1]]+epsilon2*(4.*kbt*(-12*pow(ssr,12)/pow(ya,13)+6*pow(ssr,6)/pow(ya,7)));
      ax[i] = ax[i-1]-1;
      bx[i] = bx[i-1];
      if (ax[i]<lowb) {ax[i] = lowb;}
    }else if (pp >= ploc2 && pp < ploc3 && ya < ssq){
      dphiy[ax[i-1]]=dphiy[ax[i-1]]+epsilon2*(4.*kbt*(-12*pow(ssr,12)/pow(ya,13)+6*pow(ssr,6)/pow(ya,7)));
      dphiy[bx[i-1]]=dphiy[bx[i-1]]+epsilon2*(4.*kbt*(-12*pow(ssr,12)/pow(ya,13)+6*pow(ssr,6)/pow(ya,7)));
      ax[i] = ax[i-1];
      bx[i] = bx[i-1]+1;
      if (bx[i]>upb) {bx[i] = upb;}
    }else if (pp >= ploc3 && pp < ploc4 && ya < ssq){
      dphiy[ax[i-1]]=dphiy[ax[i-1]]+epsilon2*(4.*kbt*(-12*pow(ssr,12)/pow(ya,13)+6*pow(ssr,6)/pow(ya,7)));
      dphiy[bx[i-1]]=dphiy[bx[i-1]]+epsilon2*(4.*kbt*(-12*pow(ssr,12)/pow(ya,13)+6*pow(ssr,6)/pow(ya,7)));
      ax[i] = ax[i-1];
      bx[i] = bx[i-1]-1;
      if (bx[i]<lowb) {bx[i] = lowb;}
    }else if (pp >= ploc4 && pp < 1 && ya < ssq){
      dphiy[ax[i-1]]=dphiy[ax[i-1]]+epsilon2*(4.*kbt*(-12*pow(ssr,12)/pow(ya,13)+6*pow(ssr,6)/pow(ya,7)));
      dphiy[bx[i-1]]=dphiy[bx[i-1]]+epsilon2*(4.*kbt*(-12*pow(ssr,12)/pow(ya,13)+6*pow(ssr,6)/pow(ya,7)));
      ax[i] = ax[i-1];
      bx[i] = bx[i-1];
    }else{
      dphiy[ax[i-1]]=dphi[ax[i-1]];
      dphiy[bx[i-1]]=dphi[bx[i-1]];
      ax[i] = ax[i-1];
      bx[i] = bx[i-1];
    }
    if (ax[i]==bx[i]){
      ax[i]=ax[i-1];
      bx[i]=bx[i-1];
    }
    for (int n=0; n<Nx; n++){
      nr = arma::randn();
      xx(i,n) = xx(i-1,n)-1/gm*dphi[n]*dt+D*pow(dt,0.5)*nr;  
        yy(i,n) = xx(i,n);
      if (ax[i-1]!=bx[i-1]){
        if (n==ax[i]){
          yy(i,ax[i]) = yy(i-1,ax[i-1])-1/gm*dphiy[ax[i-1]]*dt+D*pow(dt,0.5)*nr;  
        }
        if (n==bx[i]){
          yy(i,bx[i]) = yy(i-1,bx[i-1])-1/gm*dphiy[bx[i-1]]*dt+D*pow(dt,0.5)*nr;  
        }
      }else{
        yy(i,ax[i]) = yy(i-2,ax[i-2]);
        yy(i,bx[i]) = yy(i-2,bx[i-2]);
      }
    }
  }

  for(int i=0; i<imax; i++){
    for (int n=0; n<Nx; n++){
      zz[i*Nx+n]=xx(i,n);
      zz[(imax-1)*Nx+Nx+i*Nx+n]=yy(i,n);
    }
  }
      
  for(int i=0; i<imax; i++){
    zz[(imax-1)*Nx+Nx+(imax-1)*Nx+Nx+i]=ax[i];
  }
      
  for(int i=0; i<imax; i++){
    zz[(imax-1)*Nx+Nx+(imax-1)*Nx+Nx+imax+i]=bx[i];
  }
      
  Rcpp::NumericVector zzz = NumericVector(zz.begin(), zz.end());

  return zzz;
}
