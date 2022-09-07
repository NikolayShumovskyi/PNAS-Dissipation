#include "rng.h"
#include <math.h>
#define mult 314159269   /* multiplier */
#define add 907633385     /* additive constant for rng */
#define big 4294967296.0
static unsigned int rn; /* seed */
static double fact; /* normalization constant */
static double pi;

void rng_init(unsigned int value)
{
fact=1/big;
pi=8*atan((double)1)/big;
rn=value;
return;
}

unsigned int rng_get_seed(void)
{
return rn;
}

unsigned int rng_int(void)
{
  rn=rn*mult+add;
  return rn;
}

double rng(void)
{
  rn=rn*mult+add;
  return fact*rn;
}
double rng_gauss(double d)
{
  double r, phi;  
  rn=rn*mult+add;
  r=sqrt(-2*log(rn*fact))*d;
  rn=rn*mult+add;
  phi=rn*pi;
  return r*cos(phi);
}
