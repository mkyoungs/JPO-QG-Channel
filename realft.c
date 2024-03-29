#include <math.h>

void realft(data, n, isign)
float data[];
int n;
int isign;
{
/*  void four1(float data[], int nn, int isign); */
  unsigned long i,i1,i2,i3,i4,np1;
  float c1=0.5,c2,h1r,h1i,h2r,h2i;
  double wr,wi,wpr,wpi,wtemp,theta;
  
  theta=3.141592653589793/(double) (n>>1);
  if (isign == 1)
    {
      c2 = -0.5;
      four1(data,n>>1,1);
    }
  else
    {
      c2=0.5;
      theta = -theta;
    }
  wtemp=sin(0.5*theta);
  wpr = -2.0*wtemp*wtemp;
  wpi=sin(theta);
  wr=1.0+wpr;
  wi=wpi;
  np1=n+1;
  for (i=2;i<=(n>>2);i++)
    {
      i4=1+(i3=np1-(i2=1+(i1=i+i-2)));
      h1r=c1*(data[i1]+data[i3]);
      h1i=c1*(data[i2]-data[i4]);
      h2r = -c2*(data[i2]+data[i4]);
      h2i=c2*(data[i1]-data[i3]);
      data[i1]=h1r+wr*h2r-wi*h2i;
      data[i2]=h1i+wr*h2i+wi*h2r;
      data[i3]=h1r-wr*h2r+wi*h2i;
      data[i4] = -h1i+wr*h2i+wi*h2r;
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
    }
  if (isign == 1)
    {
      data[0] = (h1r=data[0])+data[1];
      data[1] = h1r-data[1];
    }
  else
    {
      data[0]=c1*((h1r=data[0])+data[1]);
      data[1]=c1*(h1r-data[1]);
      four1(data,n>>1,-1);
    }
}
