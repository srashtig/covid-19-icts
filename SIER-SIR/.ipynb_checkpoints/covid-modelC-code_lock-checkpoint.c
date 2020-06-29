#include<stdio.h>
#include<math.h>

main()
{
 long N,n,i;
 double t,t0,T,tl,td,tr,s,e,iu,id,rd,ru,s0,e0,iu0,id0,rd0,ru0;
 double ss,ee,iiu,iid,rrd,rru,ss1,ee1,iiu1,iid1,rrd1,rru1,ss2,ee2,iiu2,iid2,rrd2,rru2;
 double bet,sig,gamd,gamr,p,lam,lam1,lam2,betf,tau,tint,bet0;
 double h,ks1,ks2,ks3,ks4,ke1,ke2,ke3,ke4,kid1,kid2,kid3,kid4,kiu1,kiu2,kiu3,kiu4;
 double krd1,krd2,krd3,krd4,kru1,kru2,kru3,kru4;

 FILE *f1;
 f1=fopen("seir-p_newIC_largetime_p0.0033_lock.dat","w");
//  f1=fopen("confirmedcases_model_covid19_UK_p0.0018.dat","w");

 N=60000000; //approx population of UK
 T=1000; //no. of days
 h=0.1; //increment in time

 sig=1.0/5;
// gamd=1.0/7;
 gamr=1.0/15;

 tl=1.0/sig;
 //td=1.0/gamd;
 tr=1.0/gamr;

 p=0.0;
 bet0=0.04369/(0.07257-0.015238*p);
 p=0.0033;
 betf=0.2;

 tau=4.0;
 tint=40.0;

 lam=0.22;
 lam1=-0.4868;
 lam2=-0.14273;

 //initialisation
 t0=20;

 rrd=1.628;
 iiu=lam/(p*gamr)*rrd;;
 rru=(1-p)*gamr/lam*iiu;
 ee=bet0/(sig+lam)*iiu;
 ss=-bet0/lam*iiu;
 id0=0.0;

 e0=ee*exp(lam*t0);//+ee1*exp(lam1*t0)+ee2*exp(lam2*t0);
 iu0=iiu*exp(lam*t0);//+iiu1*exp(lam1*t0)+iiu2*exp(lam2*t0);
 rd0=rrd*exp(lam*t0);//+rrd1*exp(lam1*t0)+rrd2*exp(lam2*t0);
 ru0=rru*exp(lam*t0);//+rru1*exp(lam1*t0)+rru2*exp(lam2*t0);
 s0=ss*exp(lam*t0);//+ss1*exp(lam1*t0)+ss2*exp(lam2*t0);
 s0=N+s0;

 printf("%lf %lf %lf %lf %lf %lf %lf\n",bet,s0,e0,iu0,rd0,ru0,ss);

// R-K integration of the dynamics of the modified-SEIR model

 n=(T-t0)/h;

 for(i=1;i<n+1;i++)
 {
  if((i*h+t0)<tint) {bet=bet0;}
  else {bet = betf+(bet0-betf)*exp(-((i*h+t0)-tint)/tau);}

  ks1=-bet/N*(iu0)*s0;        //R-K step 1
  ke1=bet/N*(iu0)*s0 - sig*e0;
  kiu1=sig*e0 - gamr*iu0;
  krd1=p*gamr*iu0;
  kru1=(1.0-p)*gamr*iu0;

  ks2=-bet/N*((iu0+h*kiu1/2))*(s0+h*ks1/2); //R-K step 2
  ke2=bet/N*((iu0+h*kiu1/2))*(s0+h*ks1/2) - sig*(e0+h*ke1/2);
  kiu2=sig*(e0+h*ke1/2) - gamr*(iu0+h*kiu1/2);
  krd2=p*gamr*(iu0+h*kiu1/2);
  kru2=(1-p)*gamr*(iu0+h*kiu1/2);

  ks3=-bet/N*((iu0+h*kiu2/2))*(s0+h*ks2/2); //R-K step 3
  ke3=bet/N*((iu0+h*kiu2/2))*(s0+h*ks2/2) - sig*(e0+h*ke2/2);
  kiu3=sig*(e0+h*ke2/2) - gamr*(iu0+h*kiu2/2);
  krd3=p*gamr*(iu0+h*kiu2/2);
  kru3=(1.0-p)*gamr*(iu0+h*kiu2/2);

  ks4=-bet/N*((iu0+h*kiu3))*(s0+h*ks3); //R-K step 4
  ke4=bet/N*((iu0+h*kiu3))*(s0+h*ks3) - sig*(e0+h*ke3);
  kiu4=sig*(e0+h*ke3) - gamr*(iu0+h*kiu3);
  krd4=p*gamr*(iu0+h*kiu3);
  kru4=(1.0-p)*gamr*(iu0+h*kiu3);

  s=s0+h/6*(ks1+2*ks2+2*ks3+ks4);
  e=e0+h/6*(ke1+2*ke2+2*ke3+ke4);
  iu=iu0+h/6*(kiu1+2*kiu2+2*kiu3+kiu4);
  rd=rd0+h/6*(krd1+2*krd2+2*krd3+krd4);
  ru=ru0+h/6*(kru1+2*kru2+2*kru3+kru4);

  fprintf(f1,"%lf %lf %lf %lf %lf %lf\n",t0+i*h,s,e,iu,rd,ru);

  s0=s;
  e0=e;
  iu0=iu;
  rd0=rd;
  ru0=ru;
 }

 fclose(f1);
}
