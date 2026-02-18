#read recp from rgsim_ver1.R

library(data.table)

rd1<-fread("norm_contact.BRD3179_40000.csv")
rr<-rd1$norm_contact
ss<-rep(0,length(rr))
k=1
m<-max(rr)
for (i in 1:length(rr)){
  if (rr[i]<m){
    ss[k]=rr[i]
    k=k+1
  } 
}

rd2<-read.csv("mean_contact_frequency.BRD3179_40000.csv")

x1<-log10(recp)
x2<-log10(rd1$norm_contact)
x3<-log10(rd2$mean_norm_contact)

hist(x1, xlab="log CP", main="Histogram of simulated CP", cex.main=2, cex.lab=1.5,cex.axis=2)
hist(x2, xlab="log CP", main="Histogram of original CP data", cex.main=2, cex.lab=1.5,cex.axis=2)
hist(x3, xlab="log CP", main="Histogram of original CP data", cex.main=2, cex.lab=1.5,cex.axis=2)