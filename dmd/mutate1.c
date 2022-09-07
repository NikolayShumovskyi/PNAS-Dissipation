#include <stdio.h> /* standart input,output */
#include <ctype.h> /* for charachter recognition */
#include <stdlib.h> /* for conversion from char to dec */
#include <strings.h>
#include <time.h>
#include <math.h>
#include "bcp.h"
#include "mutate.h"
#include "rng.h"

enum mutate_keys 
{MIN_KEY=0,GROUP_KEY=MIN_KEY,IN_KEY,OUT_KEY,X_KEY,Y_KEY,Z_KEY,IN_T_NEW_KEY,
OUT_T_NEW_KEY,END_KEY,MAX_KEY};

static int nmut=0;
static int nT=0;
static int ngroup=0;
static int n_atom;
static double next_time=0;
static double init_time=0;
static int good=0;
static char ** keywords;
static double L[3];
static double L2[3];
static double X0[3]={0,0,0};
static double A[3]={0,0,0};
static double RHS=1;
static int FROM[1000];
static int FROM_T[1000];
static int TO[1000];
static double T[1000];
static int IO[1000];
static int IO_T[1000];
static int *c=0;
static double *vx;
static double *vy;
static double *vz;
char keyword[60];

int mutate_init_param(void)
{
  int i;
  for(i=0;i<3;i++)
    {
      A[i]=0;
      X0[i]=0;
    }
  nmut=0;
  nT=0;
  RHS=1;
}

int init_mutate_keywords(void)
{ 
  int i;   
  keywords=(char **)malloc(MAX_KEY*sizeof(unsigned char *));
  if(!keywords) return 0;
  keywords[GROUP_KEY]="SHAPE";
  keywords[X_KEY]="X";
  keywords[Y_KEY]="Y";
  keywords[Z_KEY]="Z";
  keywords[IN_KEY]="IN";
  keywords[OUT_KEY]="OUT";
  keywords[IN_T_NEW_KEY]="IN_T_NEW";
  keywords[OUT_T_NEW_KEY]="OUT_T_NEW";
  keywords[END_KEY]="END";
   mutate_init_param();
  return 1;
}
extern int mutate(void)
{
  if(good)
    {  
      atom * a=(atom *)get_atom();
      int i;
      int k=0;
      update_atoms();
      for(i=0;i<n_atom;i++)
	if(a[i].c!=c[i])
	  {
	    k++;
	    a[i].c=c[i];
	  }
      for(i=0;i<n_atom;i++)
	{
	  if(a[i].v.x==a[i].u.x)
	    {
	      k++;
	      a[i].v.x=a[i].u.x;
	    }
	  if(a[i].v.y==a[i].u.y)
	    {
	      k++;
	      a[i].v.y=a[i].u.y;
	    }
	  if(a[i].v.z==a[i].u.z)
	    {
	      k++;
	      a[i].v.z=a[i].u.z;
	    }
	}
      free(c);
      c=0;
      good=0;
      return k;
    }
  return 0;
}

int count(void)
{
  iatom * a=(iatom *)get_atom();
  int n_atom=get_atom_number();
  int n_dim=get_dimension();
  dimensions * bound=get_bounds();
  int i,j,k;
  int l=0;
  double corr=get_corr();
  for(i=0;i<n_dim;i++)
    {
      L[i]=bound[i].length;
      L2[i]=L[i]*0.5;
    }
  moveatoms();
  for(i=0;i<n_atom;i++)
    for(j=0;j<nmut;j++)
      if(FROM[j]==a[i].c)
	{
	  double s=0;
	  for(k=0;k<n_dim;k++)
	    {
	      double dx=a[i].q[k]-X0[k];
	      if(dx>L2[k])dx-=L[k];
	      if(dx<-L2[k])dx+=L[k];
	      s+=dx*dx*A[k];
	    }
	  if(s*IO[j]<RHS*IO[j])
	    {
	      c[i]=TO[j];
	      l++;
	    }
	}

  for(i=0;i<n_atom;i++)
    {
      a[i].u.x=a[i].v.x;
      a[i].u.y=a[i].v.y;
      a[i].u.z=a[i].v.z;
      for(j=0;j<nT;j++)
	if((FROM_T[j]==a[i].c)||!FROM_T[j])
	  {
	    double s=0;
	    for(k=0;k<n_dim;k++)
	      {
		double dx=a[i].q[k]-X0[k];
		if(dx>L2[k])dx-=L[k];
		if(dx<-L2[k])dx+=L[k];
		s+=dx*dx*A[k];
	      }
	    if(s*IO_T[j]<RHS*IO[j])
	      {
		a[i].u.x=rng_gauss(T[j]/a.m[i])*corr;
		a[i].u.y=rng_gauss(T[j]/a.m[i])*corr;
		a[i].u.z=rng_gauss(T[j]/a.m[i])*corr;
		l++;
	      }
	  }
    }
  return l;
}


extern int mutate_init(char * mut_name)
{
  FILE * infile;
  FILE * save=NULL;
  char value1[60];
  char value2[60];
  double dummy1,dummy2;   
  atom * a;
  int i,nt;
  n_atom=get_atom_number();
  good=0;
  init_param();
  ngroup=0;
  infile=fopen(mut_name,"r");
  if(!infile) return good;

  
  if(!init_mutate_keywords())
    return good;
  if(!c)
    {
      c=(int *)malloc(n_atom*sizeof(int));
      a=get_atom();
      for(i=0;i<n_atom;i++)
	c[i]=a[i].c;
    }
 	      
  while(!feof(infile))
    {
      enum mutate_keys key;
      fscanf(infile,"%s%s%s",keyword,value1,value2);
      for(key=MIN_KEY;key<MAX_KEY;key++)
	if(!strcmp(keyword,keywords[key])){
	  switch(key){
	  case GROUP_KEY:{
	    if(nmut+nT){
	      if(count())ngroup++;
	      mutate_init_param(); 
	    }
	    break;
	  }	      
	  case X_KEY:{
	    dummy1=atof(value1);
	    dummy2=atof(value2);
	    X0[0]=dummy1;
	    if(!dummy2)
	      A[0]=0;
	    else
	      A[0]=1/(dummy2*dummy2);
	    break;
	  }
	  case Y_KEY:{
	    dummy1=atof(value1);
	    dummy2=atof(value2);
	    X0[1]=dummy1;
	    if(!dummy2)
	      A[1]=0;
	    else
	      A[1]=1/(dummy2*dummy2);
	    break;
	  }
	  case Z_KEY:{
	    dummy1=atof(value1);
	    dummy2=atof(value2);
	    X0[2]=dummy1;
	    if(!dummy2)
	      A[2]=0;
	    else
	      A[2]=1/(dummy2*dummy2);
	    break;
	  }
	  case IN_KEY:{
	    dummy1=atof(value1);
	    dummy2=atof(value2);
	    FROM[nmut]=(int)dummy1;
	    TO[nmut]=(int)dummy2;
	    IO[nmut]=1;
	    nmut++;
	    break;
	    }
	  case OUT_KEY:{
	    dummy1=atof(value1);
	    dummy2=atof(value2);
	    FROM[nmut]=(int)dummy1;
	    TO[nmut]=(int)dummy2;
	    IO[nmut]=-1;
	    nmut++;
	    break;
	  }
	  case IN_T_NEW_KEY:{
	    dummy1=atof(value1);
	    dummy2=atof(value2);
	    FROM_T[nT]=(int)dummy1;
	    T[nT]=dummy2;
	    IO_T[nT]=1;
	    nT++;
	    break;
	  }
	  case OUT_T_NEW_KEY:{
	    dummy1=atof(value1);
	    dummy2=atof(value2);
	    FROM_T[nT]=(int)dummy1;
	    T[nT]=(int)dummy2;
	    IO_T[nT]=-1;
	    nT++;
	    break;
	  }

	  case END_KEY:{
	    if(nmut+nT){
	      if(count())ngroup++;
	      mutate_init_param();
	    }
	    goto finish;
	  }
	  }
	}
    }
 finish:
  fclose(infile);
  if(ngroup)
    {
      ngroup=0;
      good=1;
      return good;
    }
  free(c);
  c=0;
  return good; 
    }

