#include <stdio.h>
#include <math.h>
#ifndef PI                      /* as in stroustrup */
#define PI  M_PI
#endif

#define FULLOUTPUT

#include "ch1.h"
#include "init.c"
#define NZ 1
#define NB 2

#define ARRAY(x) float x[NZ][NYP][NX]
#define VECTOR(x) float x[NZ][NYP]
#define ARRAY0(x) float x[NYP][NX]
#define VECTOR0(x) float x[NYP]

#define LOOP(N1,N2,N3) for(k=0;k<N3;k++) for(j=0;j<N2;j++) for(i=0;i<N1;i++)
#define LOOPXY(N1,N2) for(j=0;j<N2;j++) for(i=0;i<N1;i++)
#define LOOPYZ(N2,N3) for(k=0;k<N3;k++) for(j=0;j<N2;j++)
#define LOOPZ(N3) for(k=0;k<N3;k++)
#define LOOPY(N2) for(j=0;j<N2;j++)
#define LOOPX(N1) for(i=0;i<N1;i++)

int qpsi();
int ubqby();
int rhsq();
int rhsu();
int fricq();
int fricu();

ARRAY(psi);
//ARRAY(q);
ARRAY(qreal);
ARRAY(work);
ARRAY(dqdt);
ARRAY(dq1dt);
ARRAY(dq2dt);
ARRAY0(wv);
ARRAY0(trnc);
ARRAY0(filter);
ARRAY(u);
ARRAY(v);

float dx;
float fr1;
float fr2;
//float beta;
//float r,r1,kappa;
//float width,f0;
float k0;

#define writef2(yvar,nnx,nny,nnz) fwrite(yvar,sizeof(float),nnx*nny*nnz,stdout)

main(argc,argv)
int argc;
char *argv[];
{
  int i,j,k;
  float t=0;
  float dt,dt2rk,dt0,dt1,dt2;
  int icase=0;
  int kmax2=(NX/2-1)*(NX/2-1);
  float cphi=0.65*NX/2;
  float qmax;
  float u1max,u2max;
  int tcount=0;
  float *dqpntr,*dupntr;
  FILE *dbg;
  
  init();
  fr1=1.0/rd/rd/(1.0+del);fr2=del*fr1;
  dt=1.0/dtinv;
  k0=PI/width;

  fprintf(stderr,"fr1=%f, fr2=%f, beta=%f, dt=%g, tmax=%g\n",fr1,fr2,beta,dt,tmax);
  /* scanf("%f %f %f",&r1,&u1max,&u2max); */

  dt2rk=dt/2;
  dt0=23.0*dt/12.0;
  dt1=-4.0*dt/3.0;
  dt2=5.0*dt/12.0;

  LOOPXY(NX,NYP) wv[j][i]=j*j+(i>>1)*(i>>1);
  LOOPY(NYP)wv[j][1]=j*j+NX*NX/4;
  LOOPXY(NX,NYP) trnc[j][i]=(wv[j][i]<kmax2)?1.0/(NX*NY/4):0.0;
  LOOPXY(NX,NYP) filter[j][i]=(wv[j][i]<cphi*cphi)?1.0:
    exp(-18.0*pow(2*PI/NX*(sqrt((double)wv[j][i])-cphi),(double)7.0));
  LOOPY(NYP) filter[j][1]=0;
  LOOPXY(NX,NYP) wv[j][i] *= k0*k0;

  sinetrans(q,NX,NYP,NZ);

  while(t<=tmax){
    //    fprintf(stderr,"calc t=%e\n",t);
    /* calculate psi trans */
    qpsi();
    
    /* calculate q real */
    LOOP(NX,NYP,NZ) qreal[k][j][i]=q[k][j][i]*trnc[j][i];
    invsinetrans(qreal,NX,NYP,NZ);
    
    qmax=0; LOOP(NX,NY,NZ) if(qmax<qreal[k][j][i])qmax=qreal[k][j][i];
    //    fprintf(stderr,"max q %g\n",qmax);
    if(qmax>1000)exit(1);

    switch(icase){
    case 0:
      dqpntr=&dq2dt[0][0][0];
      break;
    case -1:
      dqpntr=&dq1dt[0][0][0];
      break;
    default:
      dqpntr=&dqdt[0][0][0];
      break;
    };
    
    /* calculate rhs of q and ubar (advection) */
    rhsq(dqpntr);

    /* output step */
    if((tcount%tpl)==0 && icase<1){
    /*
    if(1){
    */
      fprintf(stderr,"writing t = %f\n",t);
      writef2(&t,1,1,1);
      LOOP(NX,NYP,NZ) work[k][j][i]=psi[k][j][i]*trnc[j][i];
      invsinetrans(work,NX,NYP,NZ);
      writef2(&work[0][0][0],NX,NYP,NZ);
      writef2(&qreal[0][0][0],NX,NYP,NZ);
      writef2(&u[0][0][0],NX,NYP,NZ);
      writef2(&v[0][0][0],NX,NYP,NZ);
      /*
      writef2(dupntr,1,NYP,NZ);
      */
      fflush(stdout);
    };
    
    switch(icase){
    case 0:
      LOOP(NX,NYP,NZ) q[k][j][i] += dt*dq2dt[k][j][i];
      t += dt;tcount++;
      icase=1;
      break;
    case 1:
      LOOP(NX,NYP,NZ) q[k][j][i] += dt2rk*(dqdt[k][j][i]-dq2dt[k][j][i]);
      icase=-1;
      break;
    case -1:
      LOOP(NX,NYP,NZ) q[k][j][i] += dt*dq1dt[k][j][i];
      t += dt;tcount++;
      icase=3;
      break;
    case 3:
      LOOP(NX,NYP,NZ) q[k][j][i] += dt2rk*(dqdt[k][j][i]-dq1dt[k][j][i]);
      icase= -2;
      break;
    case -2:
      LOOP(NX,NYP,NZ) {
	q[k][j][i] += dt0*dqdt[k][j][i]+dt1*dq1dt[k][j][i]+dt2*dq2dt[k][j][i];
	dq2dt[k][j][i]=dq1dt[k][j][i];
	dq1dt[k][j][i]=dqdt[k][j][i];
      };

      t += dt;tcount++;
      break;
    };
    LOOP(NX,NYP,NZ) q[k][j][i] *= filter[j][i];
    //    fprintf(stderr,"here %g %g %g\n",t,tmax,q[0][32][0]);
  };

/*
dbg=fopen("dbg","w");
fwrite(&qreal[0][0][0],sizeof(float),NX*NYP*NZ,dbg);
fwrite(&filter[0][0],sizeof(float),NX*NYP,dbg);
fwrite(&q[0][0],sizeof(float),NX*NYP*NZ,dbg);
LOOP(NZ) qreal[k][j][i]=q[k][j][i]*trnc[j][i];
invsinetrans(qreal,NX,NYP,NZ);
fwrite(&qreal[0][0][0],sizeof(float),NX*NYP*NZ,dbg);
LOOP(NZ) qreal[k][j][i]=q[k][j][i]*trnc[j][i]*filter[j][i];
invsinetrans(qreal,NX,NYP,NZ);
fwrite(&qreal[0][0][0],sizeof(float),NX*NYP*NZ,dbg);
fclose(dbg);
*/

exit(0);
}

int qpsi()
{
int i,j;
float det;
for(j=1;j<NYP;j++){
  for(i=1;i<NX;i++){
    det=(wv[j][i]+fr1)*(wv[j][i]+fr2)-fr1*fr2;
    psi[0][j][i]= -((wv[j][i]+fr2)*q[0][j][i]+fr1*q[1][j][i])/det;
    psi[1][j][i]= -((wv[j][i]+fr1)*q[1][j][i]+fr2*q[0][j][i])/det;
  };
};
for(i=1;i<NX;i++){
  psi[0][0][i]=0.0;
  psi[1][0][i]=0.0;
};
}

int rhsq(f)
ARRAY(f); /* transform of rhs (output) */
{
  int i,j,k;
  float temp;
  /* calculate u */
  LOOP(NX,NYP,NZ) u[k][j][i] = -psi[k][j][i]*trnc[j][i];
  ddy(u,k0,NX,NYP,NZ);
  invcosinetrans(u,NX,NYP,NZ);
  
  /* calculate v */
  LOOP(NX,NYP,NZ) v[k][j][i]=psi[k][j][i]*trnc[j][i];
  ddx(v,k0,NX,NYP,NZ);
  invsinetrans(v,NX,NYP,NZ);
  
  LOOP(NX,NYP,NZ) f[k][j][i] = -(u[k][j][i]+ubar[k][j])*qreal[k][j][i];
  sinetrans(f,NX,NYP,NZ);
  ddx(f,k0,NX,NYP,NZ);
  
  LOOP(NX,NYP,NZ) work[k][j][i] = -v[k][j][i]*qreal[k][j][i];
  /* remove mean eddy flux */
  LOOPYZ(NYP,NZ){
    temp=0.0;
    LOOPX(NX) temp += work[k][j][i];
    temp=-temp/NX;
    vpqpbar[k][j]= temp;
    LOOPX(NX) work[k][j][i] += temp;
  };
  cosinetrans(work,NX,NYP,NZ);
  ddy(work,k0,NX,NYP,NZ); /*gives negative of result*/
  LOOP(NX,NYP,NZ) f[k][j][i] -= work[k][j][i];
  
  /* mean PV advection */
  LOOP(NX,NYP,NZ) work[k][j][i] = -v[k][j][i]*qbary[k][j];
  sinetrans(work,NX,NYP,NZ);
  LOOP(NX,NYP,NZ) f[k][j][i] += work[k][j][i];
  
  /* friction */
  LOOPXY(NX,NYP) {
    f[0][j][i] -= r1*fr1*(q[0][j][i]-q[1][j][i]
		+(fr1+fr2)*(psi[0][j][i]-psi[1][j][i]))
                +kappa*wv[j][i]*q[0][j][i];
    f[1][j][i] -= r1*fr2*(q[1][j][i]-q[0][j][i]
                +(fr1+fr2)*(psi[1][j][i]-psi[0][j][i]))
                +r*fr2*(q[1][j][i]+fr2*(psi[1][j][i]-psi[0][j][i]))
                +kappa*wv[j][i]*q[1][j][i];
  };
}

int ubqby()
{
  int j,k;
  /* filter ubar */
  LOOPZ(NZ) sinft(&ubar[k][0],NY,0);
  LOOPYZ(NYP,NZ) qbary[k][j]=ubar[k][j];
  LOOPYZ(NYP,NZ) ubar[k][j] *= filter[j][0]*trnc[j][0]*NX/2;
  LOOPZ(NZ) sinft(&ubar[k][0],NY,0);

  LOOPYZ(NYP,NZ) qbary[k][j] *= j*j*k0*k0*trnc[j][0]*NX/2;
  LOOPZ(NZ) sinft(&qbary[k][0],NY,0);
  LOOPY(NYP) {
    qbary[0][j]=beta+qbary[0][j]+fr1*(ubar[0][j]-ubar[1][j]);
    qbary[1][j]=beta+qbary[1][j]+fr2*(ubar[1][j]-ubar[0][j]);
  };
}

int rhsu(f)
VECTOR(f);
{
  int j,k;

  /* calculate meridional circ */
  LOOPY(NYP) {
    ru[0][j] = forcing[0][j]-fr1*r1*(ubar[0][j]-ubar[1][j]);
    ru[1][j] = -fr2*r1*(ubar[1][j]-ubar[0][j])-fr2*r*ubar[1][j];
  };
  LOOPYZ(NYP,NZ) ru[k][j] += vpqpbar[k][j]-kappa*(qbary[k][j]-beta);
  
  LOOPY(NYP) phi[j]=ru[0][j]-ru[1][j];
  sinft(phi, NY, 0);
  for(j=1;j<NYP;j++) phi[j]= -phi[j]/(j*j*k0*k0+fr1+fr2)*trnc[j][0]*NX/2;
  sinft(phi,NY,0);
  /* calculate dubardt */
  LOOPY(NYP) {
    f[0][j]=ru[0][j]+fr1*phi[j];
    f[1][j]=ru[1][j]-fr2*phi[j];
  };
}


sinetrans(f,nx,nyp,nz)
     ARRAY(f);
     int nx,nyp,nz;
{
  int i,j,k;
  for(k=0;k<nz;k++){
    for(i=0;i<nx;i++)sinft(&f[k][0][i],nyp-1,SHIFT);
    for(j=0;j<nyp;j++)realft(&f[k][j][0],nx,1);
  };
}

cosinetrans(f,nx,nyp,nz)
     ARRAY(f);
     int nx,nyp,nz;
{
  int i,j,k;
  for(k=0;k<nz;k++){
    for(i=0;i<nx;i++)cosft(&f[k][0][i],nyp-1,SHIFT);
    for(j=0;j<nyp;j++)realft(&f[k][j][0],nx,1);
  };
}

invsinetrans(f,nx,nyp,nz)
     ARRAY(f);
     int nx,nyp,nz;
{
  int i,j,k;
  for(k=0;k<nz;k++){
    for(j=0;j<nyp;j++)realft(&f[k][j][0],nx,-1);
    for(i=0;i<nx;i++)sinft(&f[k][0][i],nyp-1,SHIFT);
  };
}

invcosinetrans(f,nx,nyp,nz)
     ARRAY(f);
     int nx,nyp,nz;
{
  int i,j,k;
  for(k=0;k<nz;k++){
    for(j=0;j<nyp;j++)realft(&f[k][j][0],nx,-1);
    for(i=0;i<nx;i++)cosft(&f[k][0][i],nyp-1,SHIFT);
  };
}

ddx(f,k0,nx,nyp,nz)
     ARRAY(f);
     float k0;
     int nx,nyp,nz;
{
  int i,j,k,i1;
  float temp;
  for(k=0;k<nz;k++)
    for(j=0;j<nyp;j++){
      f[k][j][0]=0.0;
      f[k][j][1]=0.0;
      for(i=1;i<nx>>1;i++){
	i1=i<<1;
	temp=-i*k0*f[k][j][i1];
	f[k][j][i1]=i*k0*f[k][j][i1+1];
	f[k][j][i1+1]=temp;
      };
    };
}

ddy(f,k0,nx,nyp,nz)
     ARRAY(f);
     float k0;
     int nx,nyp,nz;
{
  int i,j,k;
  for(k=0;k<nz;k++)
    for(j=0;j<nyp;j++){
      for(i=1;i<nx;i++){
	f[k][j][i]=j*k0*f[k][j][i];
      };
    };
}
