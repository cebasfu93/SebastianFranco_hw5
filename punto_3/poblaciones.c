#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define pi 3.1415
#define iter 1000000
#define N 96
#define sigma 0.01

float rand_norm();
float * modelo(float a, float b, float d, float g, float x_ini, float y_ini, float * time_teo);
float likelihood(float * y_calc, float * y_teo);
float alpha(float * y_ahora, float * y_antes, float * y_teo);

int main(int argc, char const *argv[]){
  int i;
  int j;

  FILE *in;
  float n1,n2,n3;
  char bas1[200];
  int res;
  float *time_obs;
  float *x_obs; //x=prey
  float *y_obs; //y=predator
  time_obs=malloc(N*sizeof(float));
  x_obs=malloc(2*N*sizeof(float));

  in = fopen("lotka_volterra_obs.dat", "r");

  for(i=0;i<5;i++){
    res = fscanf(in, "%s \n", &bas1);
  }

  for(i=0;i<N;i++){
    res = fscanf(in, "%f %f %f\n", &n1, &n2, &n3);
    time_obs[i]=n1;
    x_obs[i]=n2;
    x_obs[N+i]=n3;
  }
  fclose(in);

  srand48(1);
  float *alfas, *betas, *deltas, *gammas, *like;
  alfas=malloc(iter*sizeof(float));
  betas=malloc(iter*sizeof(float));
  deltas=malloc(iter*sizeof(float));
  gammas=malloc(iter*sizeof(float));
  like=malloc(iter*sizeof(float));

//Despues de correr muchas veces el makefile escogi estos valores iniciales para que las graficas se vieran descentes
//Si se cambian, igual converge a los mismos valores
  alfas[0]=28;
  betas[0]=7;
  deltas[0]=2;
  gammas[0]=6;

  float like_now, like_bef;
  float *data_now, *data_bef;
  data_now=malloc(2*N*sizeof(float));
  data_bef=malloc(2*N*sizeof(float));

  data_bef=modelo(alfas[0], betas[0], deltas[0], gammas[0], x_obs[0], x_obs[N], time_obs);

  like[0]=likelihood(data_bef, x_obs);

  float a_temp, b_temp, g_temp, d_temp;

  for(j=1;j<iter;j++){

    a_temp=alfas[j-1]+rand_norm();
    b_temp=betas[j-1]+rand_norm();
    g_temp=gammas[j-1]+rand_norm();
    d_temp=deltas[j-1]+rand_norm();

    data_now=modelo(a_temp, b_temp, d_temp, g_temp, x_obs[0], x_obs[N], time_obs);

    float r=alpha(data_now, data_bef, x_obs);

    if(r>=1.0){
      alfas[j]=a_temp;
      betas[j]=b_temp;
      gammas[j]=g_temp;
      deltas[j]=d_temp;
      like[j]=likelihood(data_now, x_obs);
    }

    else{
      alfas[j]=alfas[j-1];
      betas[j]=betas[j-1];
      gammas[j]=gammas[j-1];
      deltas[j]=deltas[j-1];
      like[j]=like[j-1];
    }
    for(i=0;i<2*N;i++){
      data_bef[i]=data_now[i];
    }
  }

  for(i=0;i<iter;i++){
    printf("%f %f %f %f %f \n", alfas[i], betas[i], gammas[i], deltas[i], like[i]);
  }

  return 0;
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

float * modelo(float a, float b, float d, float g, float x_ini, float y_ini, float *time_teo){
  int i;
  float *x_calc;
  x_calc=malloc(2*N*sizeof(float));
  x_calc[0]=x_ini;
  x_calc[N]=y_ini;
  for(i=1;i<N;i++){
    float delt=time_teo[i]-time_teo[i-1];
    x_calc[i]=x_calc[i-1]+delt*x_calc[i-1]*(a-b*x_calc[N+i-1]);
    x_calc[N+i]=x_calc[N+i-1]-delt*x_calc[N+i-1]*(g-d*x_calc[i-1]);
  }
  return x_calc;

}

float alpha(float *y_ahora, float *y_antes, float *y_teo){
  int i;
  float chi_now=0;
  float chi_bef=0;
  float alfa;

  for(i=0;i<2*N;i++){
    chi_now=chi_now+pow(y_ahora[i]-y_teo[i],2);
    chi_bef=chi_bef+pow(y_antes[i]-y_teo[i],2);
  }

  alfa=exp((-chi_now+chi_bef)/200.0);
  return alfa;
}

float likelihood(float *y_calc, float *y_teo){
  int i;
  float chi2=0;
  float hood;
  for(i=0;i<2*N;i++){
    chi2=chi2+pow(y_calc[i]-y_teo[i],2);
  }

  hood=exp(-chi2/2.0);
  return hood;
}
