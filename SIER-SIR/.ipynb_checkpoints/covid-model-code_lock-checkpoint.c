#include<stdio.h>
#include<math.h>

main()
{
 long N,n,i;
 double t,t0,T,tl,td,tr,s,e,iu,id,rd,ru,s0,e0,iu0,id0,rd0,ru0,tau,tint;
 double ss,ee,iiu,iid,rrd,rru,ss1,ee1,iiu1,iid1,rrd1,rru1,ss2,ee2,iiu2,iid2,rrd2,rru2,ssalt;
 double bet,sig,gamd,gamr,p,lam,lam1,lam2,betf;
 double h,ks1,ks2,ks3,ks4,ke1,ke2,ke3,ke4,kid1,kid2,kid3,kid4,kiu1,kiu2,kiu3,kiu4;
 double krd1,krd2,krd3,krd4,kru1,kru2,kru3,kru4;

 FILE *f1;
// f1=fopen("confirmedcases_model_covid19_spain_lin_correct_prot1.dat","w");
  f1=fopen("confirmedcases_model_covid19_India_p0.001_lock_betf0.1.dat","w");

 N=130*10000000; //approx population of India
 T=150; //no. of days
 h=0.1; //increment in time

 sig=1.0/5;
 gamd=1.0/7;
 gamr=1.0/15;

 tl=1.0/sig;
 td=1.0/gamd;
 tr=1.0/gamr;

 tau=3;
 tint=25;

 p=0.001;
 bet=0.02776/(0.0628-0.015238*p);
 betf=0.1;

/* p=0.1;
 bet=0.5816/(0.1686-0.0152*p);*/

 lam=0.1713;
 lam1=-0.4868;
 lam2=-0.14273;

 //initialisation, set 1
 /*t0=2;
 rd0=2.0;
 id0=rd0*td/t0;
 iu0=rd0*(1-p)/p*((t0+td)/(t0+tr))*tr/t0;
 ru0=rd0*(1-p)/p*((t0+td)/(t0+tr));
 e0=rd0/p*(t0+td)*tl/(t0*t0);
 s0=N-(e0+iu0+id0+rd0+ru0);*/

 //initialisation, set 2
 t0=15.0;

 rrd=8.6;
 iid=lam*td*rrd;
 iiu=(1-p)/p*(gamd+lam)/(gamr+lam)*iid;
 rru=gamr/lam*iiu;
 ee=(gamd+lam)/(p*sig)*iid;
 ss=-(ee+iiu+iid+rrd+rru);
 ssalt=-(sig+lam)/lam*ee;

 printf("%lf %lf %lf %lf %lf %lf %lf\n",rrd,iid,iiu,rru,ee,ss,ssalt);

 rrd1=0;
 iid1=lam1*td*rrd1;
 iiu1=(1-p)/p*(gamd+lam1)/(gamr+lam1)*iid1;
 rru1=gamr/lam1*iiu1;
 ee1=(gamd+lam1)/(p*sig)*iid1;
 ss1=-(ee1+iiu1+iid1+rrd1+rru1);

 rrd2=0;
 iid2=lam2*td*rrd2;
 iiu2=(1-p)/p*(gamd+lam2)/(gamr+lam2)*iid2;
 rru2=gamr/lam2*iiu2;
 ee2=(gamd+lam2)/(p*sig)*iid2;
 ss2=-(ee2+iiu2+iid2+rrd2+rru2);

 e0=ee*exp(lam*t0)+ee1*exp(lam1*t0)+ee2*exp(lam2*t0);
 iu0=iiu*exp(lam*t0)+iiu1*exp(lam1*t0)+iiu2*exp(lam2*t0);
 id0=iid*exp(lam*t0)+iid1*exp(lam1*t0)+iid2*exp(lam2*t0);
 rd0=rrd*exp(lam*t0)+rrd1*exp(lam1*t0)+rrd2*exp(lam2*t0);
 ru0=rru*exp(lam*t0)+rru1*exp(lam1*t0)+rru2*exp(lam2*t0);
 s0=ss*exp(lam*t0)+ss1*exp(lam1*t0)+ss2*exp(lam2*t0);
 s0=N-s0;

 printf("%ld	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf\n",N,bet,s0,e0,iu0,id0,rd0,ru0,ss);

// R-K integration of the dynamics of the modified-SEIR model

 n=(T-t0)/h;

 for(i=1;i<n+1;i++)
 {
  if(i*h<tint) {bet=bet;}
  else {bet = betf+(bet-betf)*exp(-(i*h-tint)/tau);}

  ks1=-bet/N*(iu0+id0)*s0;        //R-K step 1
  ke1=bet/N*(iu0+id0)*s0 - sig*e0;
  kiu1=(1.0-p)*sig*e0 - gamr*iu0;
  kid1=p*sig*e0 - gamd*id0;
  krd1=gamd*id0;
  kru1=gamr*iu0;

  ks2=-bet/N*((iu0+h*kiu1/2)+(id0+h*kid1/2))*(s0+h*ks1/2); //R-K step 2
  ke2=bet/N*((iu0+h*kiu1/2)+(id0+h*kid1/2))*(s0+h*ks1/2) - sig*(e0+h*ke1/2);
  kiu2=(1.0-p)*sig*(e0+h*ke1/2) - gamr*(iu0+h*kiu1/2);
  kid2=p*sig*(e0+h*ke1/2) - gamd*(id0+h*kid1/2);
  krd2=gamd*(id0+h*kid1/2);
  kru2=gamr*(iu0+h*kiu1/2);

  ks3=-bet/N*((iu0+h*kiu2/2)+(id0+h*kid2/2))*(s0+h*ks2/2); //R-K step 3
  ke3=bet/N*((iu0+h*kiu2/2)+(id0+h*kid2/2))*(s0+h*ks2/2) - sig*(e0+h*ke2/2);
  kiu3=(1.0-p)*sig*(e0+h*ke2/2) - gamr*(iu0+h*kiu2/2);
  kid3=p*sig*(e0+h*ke2/2) - gamd*(id0+h*kid2/2);
  krd3=gamd*(id0+h*kid2/2);
  kru3=gamr*(iu0+h*kiu2/2);

  ks4=-bet/N*((iu0+h*kiu3)+(id0+h*kid3))*(s0+h*ks3); //R-K step 4
  ke4=bet/N*((iu0+h*kiu3)+(id0+h*kid3))*(s0+h*ks3) - sig*(e0+h*ke3);
  kiu4=(1.0-p)*sig*(e0+h*ke3) - gamr*(iu0+h*kiu3);
  kid4=p*sig*(e0+h*ke3) - gamd*(id0+h*kid3);
  krd4=gamd*(id0+h*kid3);
  kru4=gamr*(iu0+h*kiu3);

  s=s0+h/6*(ks1+2*ks2+2*ks3+ks4);
  e=e0+h/6*(ke1+2*ke2+2*ke3+ke4);
  iu=iu0+h/6*(kiu1+2*kiu2+2*kiu3+kiu4);
  id=id0+h/6*(kid1+2*kid2+2*kid3+kid4);
  rd=rd0+h/6*(krd1+2*krd2+2*krd3+krd4);
  ru=ru0+h/6*(kru1+2*kru2+2*kru3+kru4);

  fprintf(f1,"%lf	%lf	%lf	%lf	%lf	%lf	%lf\n",t0+i*h,s,e,iu,id,rd,ru);

  s0=s;
  e0=e;
  iu0=iu;
  id0=id;
  rd0=rd;
  ru0=ru;
 }

 fclose(f1);
}
