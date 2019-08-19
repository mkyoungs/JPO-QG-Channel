#include <math.h>
#include <stdio.h>

void sinft(y, n, ishift)
float y[];
int n;
int ishift;
{
/*  void realftshift(float data[], int n, int isign, int ishift); */
  int j;
  float sum,y1,y2;
  double theta,wi=0.0,wr=1.0,wpi,wpr,wtemp;
 
  
  theta=3.14159265358979/(double) n;
  wtemp=sin(0.5*theta);
  wpr = -2.0*wtemp*wtemp;
  wpi=sin(theta);
  y[0]=0.0;
  for (j=1;j<=(n>>1);j++)
    { 
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
      y1=wi*(y[j<<ishift]+y[(n-j)<<ishift]);
      y2=0.5*(y[j<<ishift]-y[(n-j)<<ishift]);
      y[j<<ishift]=y1+y2;
      y[(n-j)<<ishift]=y1-y2;

    }
  realftshift(y,n,1,ishift);
  y[0]*=0.5;
  sum=y[1<<ishift]=0.0;
  for (j=0;j<n-1;j+=2)
    {
      sum += y[j<<ishift];
      y[j<<ishift]=y[(j+1)<<ishift];
      y[(j+1)<<ishift]=sum;
      //fprintf(stderr,"%g",y[j]);
    }
}
