#include <math.h>
#ifndef PI                      /* as in stroustrup */
#define PI  M_PI
#endif

void cosft(y, n, ishift)
float y[];
int n;
int ishift;
{
/*  void realftshift(float data[], int n, int isign, int ishift); */
  int j;
  float sum,y1,y2;
  double theta,wi=0.0,wpi,wpr,wr=1.0,wtemp;
  
  theta=PI/n;
  wtemp=sin(0.5*theta);
  wpr = -2.0*wtemp*wtemp;
  wpi=sin(theta);
  sum=0.5*(y[0]-y[n<<ishift]);
  y[0]=0.5*(y[0]+y[n<<ishift]);
  for (j=1;j<(n>>1);j++)
    {
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
      y1=0.5*(y[j<<ishift]+y[(n-j)<<ishift]);
      y2=(y[j<<ishift]-y[(n-j)<<ishift]);
      y[j<<ishift]=y1-wi*y2;
      y[(n-j)<<ishift]=y1+wi*y2;
      sum += wr*y2;
    }
  realftshift(y,n,1,ishift);
  y[n<<ishift]=y[1<<ishift];
  y[1<<ishift]=sum;
  for (j=3;j<n;j+=2)
    {
      sum += y[j<<ishift];
      y[j<<ishift]=sum;
    }
}
