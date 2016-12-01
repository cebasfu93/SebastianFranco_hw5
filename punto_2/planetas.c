#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define iter 1000000
#define G 6.67408e-11
#define N 8
#define au 1.496e11
#define yr 3.154e7
#define pi 3.1415
#define sigma 0.01

float modelo(float x, float m, float b);
float rand_norm();
float likelihood(float *y_calc, float *y_teo);
float alpha(float *y_ahora, float *y_antes, float *y_teo);

int main(){

  FILE *fp;
  fp=fopen("M.txt", "w");
  fprintf(fp, "%i \n", iter);

  int i;
  int j;

  float posx[N]={0.324190175, -0.701534590, -0.982564148, 1.104185888, 3.266443877, -9.218802228, 19.930781147, 24.323085642};
  float posy[N]={0.090955208, -0.168809218, -0.191145980, -0.826097003, -3.888055863, 1.788299816, -2.555241579, -17.606227355};
  float posz[N]={-0.022020510, 0.037947785, -0.000014724, 3.260215854, -0.057015321, 0.335737817, -0.267710968, -0.197974999};
  float velx[N]={-4.627851589, 1.725066954, 1.126784520, 4.524583075, 2.076140727, -0.496457364, 0.172224285, 0.664855006};
  float vely[N]={10.390063716, -7.205747212, -6.187988860, 4.524583075, 1.904040630, -2.005021061, 1.357933443, 0.935497107};
  float velz[N]={1.273504997, -0.198268558, 0.000330572, 0.014760239, -0.054374153, 0.054667082, 0.002836325, -0.034716967};

  float *pos;
  float *vel;
  pos=malloc(N*sizeof(float));
  vel=malloc(N*sizeof(float));

  for(i=0; i<N; i++){
    posx[i]=posx[i]*au;
    posy[i]=posy[i]*au;
    posz[i]=posz[i]*au;
    velx[i]=velx[i]*au/yr;
    vely[i]=vely[i]*au/yr;
    velz[i]=velz[i]*au/yr;
  }

  for(i=0;i<N;i++){
    pos[i]=pow(pow(posx[i],2)+pow(posy[i],2)+pow(posz[i],2),0.5);
    vel[i]=pow(pow(velx[i],2)+pow(vely[i],2)+pow(velz[i],2),0.5);
  }

  float *y_obs;
  float *x_obs;
  y_obs=malloc(N*sizeof(float));
  x_obs=malloc(N*sizeof(float));
  for(i=0;i<N;i++){
    y_obs[i]=log10(pow(vel[i],2))/10.0;
    x_obs[i]=log10(pos[i])/10.0;
  }

  srand48(1);

  float *pends;
  float *cortes;
  float *like;
  pends=malloc(iter*sizeof(float));
  cortes=malloc(iter*sizeof(float));
  like=malloc(iter*sizeof(float));

  pends[0]=0;
  cortes[0]=1;


  float like_now;
  float like_bef;
  float *y_now;
  float *y_bef;
  y_now=malloc(N*sizeof(float));
  y_bef=malloc(N*sizeof(float));

  for(i=0;i<N;i++){
    y_bef[i]=modelo(x_obs[i], pends[0], cortes[0]);
  }

  like[0]=likelihood(y_bef, y_obs);

  float pend_temp;
  float cort_temp;
  for(i=1; i<iter; i++){

    pend_temp = pends[i-1]+rand_norm();
    cort_temp = cortes[i-1]+rand_norm();

    for(j=0;j<N;j++){
      y_now[j]=modelo(x_obs[j], pend_temp, cort_temp);
    }

    float r=alpha(y_now, y_bef, y_obs);
    /*printf("%f %f %f\n", r, pends[i-1], cortes[i-1]);*/
    if(r>=1.0){
      pends[i]=pend_temp;
      cortes[i]=cort_temp;
      like[i]=likelihood(y_now, y_obs);
    }

    else{
      pends[i]=pends[i-1];
      cortes[i]=cortes[i-1];
      like[i]=like[i-1];
    }

    for(j=0;j<N;j++){
      y_bef[j]=y_now[j];
    }
  }

  for(i=0;i<iter;i++){
    printf("%f \n", pends[i]);
  }
  for(i=0;i<iter;i++){
    printf("%f \n", cortes[i]);
  }
  for(i=0;i<iter;i++){
    printf("%f \n", like[i]);
  }

  return 0;
}

float modelo(float x, float m, float b){
  return x*m+b;
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
    chi2=chi2+pow(y_calc[i]-y_teo[i],2);
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
    chi_now=chi_now+pow(y_ahora[i]-y_teo[i],2);
    chi_bef=chi_bef+pow(y_antes[i]-y_teo[i],2);
  }

  alfa=exp((-chi_now+chi_bef)/2000.0);
  return alfa;
}
