#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define iter 1000000
#define desv 0.10
#define sigma 0.01
#define N 6
#define v 5.0
#define scale 10.0
#define pi 3.1415

float modelo(float x_ahora, float y_ahora, float x_est, float y_est);
float rand_norm();
float likelihood(float *y_calc, float *y_teo);
float alpha(float *y_ahora, float *y_antes, float *y_teo);

int main(){

  FILE *fp;
  fp=fopen("M.txt", "w");
  fprintf(fp, "%i \n %f \n", iter, scale);

  float time_obs[N]={3.12, 2.98, 2.84, 3.26, 3.12, 2.98};
  float posx[N]={3, 4, 5, 3, 4, 5};
  float posy[N]={15, 15, 15, 16, 16, 16};

  int i;
  int j;
  for(i=0;i<N;i++){
    time_obs[i]=time_obs[i]/scale;
    posx[i]=posx[i]/scale;
    posy[i]=posy[i]/scale;
  }

  float *xs;
  float *ys;
  float *like;
  xs=malloc(iter*sizeof(float));
  ys=malloc(iter*sizeof(float));

  xs[0]=10.0/scale;
  ys[0]=8.0/scale;

  float *time_now;
  float *time_bef;
  time_now=malloc(N*sizeof(float));
  time_bef=malloc(N*sizeof(float));

  for(i=0;i<N;i++){
    time_bef[i]=modelo(xs[0], ys[0], posx[i], posy[i]);
  }

  like=malloc(iter*sizeof(float));
  like[0]=likelihood(time_bef, time_obs);

  float x_temp;
  float y_temp;

  for(j=1;j<iter;j++){
    x_temp=xs[j-1]+rand_norm();
    y_temp=ys[j-1]+rand_norm();

    for(i=0;i<N;i++){
      time_now[i]=modelo(x_temp, y_temp, posx[i], posy[i]);
    }

    float r=alpha(time_now, time_bef, time_obs);

    if(r>=1.0){
      xs[j]=x_temp;
      ys[j]=y_temp;
      like[j]=likelihood(time_now, time_obs);
    }
    else{
      xs[j]=xs[j-1];
      ys[j]=ys[j-1];
      like[j]=like[j-1];
    }

    for(i=0;i<N;i++){
      time_bef[i]=time_now[i];
    }
  }

  for(i=0;i<iter;i++){
    printf("%f \n", xs[i]);
  }

  for(i=0;i<iter;i++){
    printf("%f \n", ys[i]);
  }

  for(i=0;i<iter;i++){
    printf("%f \n", like[i]);
  }

  return 0;
}

float modelo(float x_ahora, float y_ahora, float x_est, float y_est){
  float r=pow(pow(x_ahora-x_est,2)+pow(y_ahora-y_est,2),0.5);
  return r/v ;
}

float rand_norm(){

  float phi;
  float gamma;
  float r;
  float doris;
  phi = drand48()*2*pi;
  gamma=-log(drand48());
  r=sigma*sqrt(2*gamma);
  doris=r*cos(phi);

  return doris;
}

float likelihood(float *y_calc, float *y_teo){
  int i;
  float chi2=0;
  float hood;
  for(i=0;i<N;i++){
    chi2=chi2+pow((y_calc[i]-y_teo[i])*scale/desv,2);
  }

  hood=exp(-chi2/2.0);
  return hood;
}

float alpha(float *y_ahora, float *y_antes, float *y_teo){
  int i;
  float chi_now=0;
  float chi_bef=0;
  float alfa;

  for(i=0;i<N;i++){
    chi_now=chi_now+pow((y_ahora[i]-y_teo[i])*scale/desv,2);
    chi_bef=chi_bef+pow((y_antes[i]-y_teo[i])*scale/desv,2);
  }

  alfa=exp((-chi_now+chi_bef)/200.0);
  return alfa;
}
