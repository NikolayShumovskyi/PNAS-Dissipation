#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include "nucleation.h"
#include "bcp.h"
#include "search.h"
static int at1=1;
static double qmin;
static double qmax;
static int maxN=12;
static int Neib=12;
static int n_atom;
static int ncl;
static int maxcl;
static int maxcs;
static int good;
static double dt_dn;
static FILE *nuc =0;
static FILE *maxclust =0;
static FILE *density=0;
static double nuctime;
static double **coef;
static double norm[7]; 
static double pi;
static int mp[7];
static double sm;
static moved_atom * a;
static dimensions * bound;
static int nAtoms;
static double dim;
static double range;
static double start_time;
static double finish_time;
static double rbin;
static double xc,yc,zc;
void spherh6_init(void)
{
  mp[0]=4;
  mp[1]=3;
  mp[2]=3;
  mp[3]=2;
  mp[4]=2;
  mp[5]=1;
  mp[6]=1;
  pi=4*atan((double)1);
  norm[0]=sqrt(13/pi)/32;
  norm[1]=-sqrt(273*0.5/pi)/16;
  norm[2]=sqrt(1365/pi)/64;
  norm[3]=-sqrt(1365/pi)/32;
  norm[4]=3*sqrt(91*0.5/pi)/32;
  norm[5]=-3*sqrt(1001/pi)/32;
  norm[6]=sqrt(3003/pi)/64;

  coef=(double **)malloc(7*sizeof(double*));
  coef[0]=(double *)malloc((4+3+3+2+2+1+1)*sizeof(double));
  coef[1]=coef[0]+4;
  coef[2]=coef[1]+3;
  coef[3]=coef[2]+3;
  coef[4]=coef[3]+2;
  coef[5]=coef[4]+2;
  coef[6]=coef[5]+1;
  coef[0][0]=231;
  coef[0][1]=-315;
  coef[0][2]=105;
  coef[0][3]=-5;
  coef[1][0]=33;
  coef[1][1]=-30;
  coef[1][2]=5;
  coef[2][0]=33;
  coef[2][1]=-18;
  coef[2][2]=1;
  coef[3][0]=11;
  coef[3][1]=-3;
  coef[4][0]=11;
  coef[4][1]=-1;
  coef[5][0]=1;
  coef[6][0]=1;
}
double spherh6(int m,double theta)
{
  double s=sin(theta);
  double c=cos(theta);
  int i,n=abs(m);
  double a=norm[n]*pow(s,(double)n);
  double p=0;
  if(n&1)
    {    
      a*=c;
      if(m<0)a=-a;
    }
  c*=c;
  for(i=0;i<mp[n];i++)
    p=p*c+coef[n][i];
  sm=p*a;
  //  printf("%d %lf %lf %lf %lf\n",m,theta,sm,p,c);
  return sm;
}
double spherh6_re(int m,double theta,double phi)
{
  spherh6(m,theta);
  return  sm*cos(phi*m);
}
double spherh6_im(int m,double theta,double phi)
{
  return sm*sin(phi*m);
}
double get_q6(int n, double *x, double *y, double *z)
{
  double s,a,b,r,rho,theta,phi;
  int i,m;
  s=0;
  for(m=-6;m<=6;m++)
    {
      a=0;
      b=0;
      for(i=0;i<n;i++)
	{
	  rho=x[i]*x[i]+y[i]*y[i];
	  r=sqrt(rho+z[i]*z[i]);
	  rho=sqrt(rho);
	  theta=acos(z[i]/r);
	  if((!y[i])&&(x[i]<=0))
	    phi=pi;
	  else
	    phi=2*atan(y[i]/(x[i]+rho));
	  a+=spherh6_re(m,theta,phi);
	  b+=spherh6_im(m,theta,phi);
	}
      a/=n;
      b/=n;
      s+=a*a+b*b;
    }
  return sqrt(s*pi*4/13);
}

int compare1(const void *a, const void *b)
{
  return (*(int*)b -*(int*)a);
}

int open_q6(char * fname)
{
    FILE * ff;
    char dummy[100];    
    if(good)
      {
	if(nuc)fclose(nuc);
        nuc=0;
	if(maxclust)fclose(maxclust);
        maxclust=0;
	if(density)fclose(density);
        density=0;
      }
    good=0;
    ff=fopen(fname,"r");
    if(!ff)return good;
    fscanf(ff,"%s%lf",dummy,&dt_dn);
    fscanf(ff,"%s%d",dummy,&at1);
    fscanf(ff,"%s%lf",dummy, &qmin);
    fscanf(ff,"%s%lf",dummy, &qmax);
    fscanf(ff,"%s%d",dummy,&maxN);
    fscanf(ff,"%s%d",dummy,&maxcs);
    fscanf(ff,"%s%lf",dummy,&range);  
    fscanf(ff,"%s%d",dummy,&Neib);  
    fscanf(ff,"%s%s",dummy,dummy);
    nuc=fopen(dummy,"w");
    fscanf(ff,"%s%lf",dummy,&start_time);
    fscanf(ff,"%s%lf",dummy,&finish_time);
    fscanf(ff,"%s%s",dummy,dummy);
    //printf("%s\n",dummy);
    maxclust=fopen(dummy,"w");
    fscanf(ff,"%s%lf",dummy,&rbin);
    if(feof(ff))rbin=0;
    if(rbin)
      {
 	fscanf(ff,"%s%s",dummy,dummy);
	density=fopen(dummy,"w");
      }
    n_atom=get_atom_number();
    a=(moved_atom *)get_atom();
    bound =get_bounds();
    xc=-bound[0].length;
    yc=-bound[1].length;
    zc=-bound[2].length;
    nAtoms=get_atom_number();
    dim=get_dimension();
    if(bound[2].period==1)
      dim=2;
    else 
      dim=3;
    if((at1>0)&&(qmin<qmax)&&(maxN>0)&&nuc&&(dt_dn>0)&&(maxcs>0)&&(dim==3))good=1;
    fclose(ff);
    if(good)
      {
	nuctime=get_time();
	spherh6_init();
      }
    printf("q6:%d\n",good);
    return good;
}

int compare_i(const void *a, const void *b)
{
  if(**(double**)b < **(double**)a)return 1;
     else if(**(double**)b == **(double**)a)return 0;
	     else return -1;
}


int find_neib(int i,int at1,int maxN,int * list,double * q6)
{
  int neib;
  q6[i]=0;
  if ((a[i].c-at1)*at1)
    return -1;
  else
    {
      int j,k,l;
      double lx=bound[0].length;
      double ly=bound[1].length;
      double lz=bound[2].length;
      double lx2=lx/2;
      double ly2=ly/2;
      double lz2=lz/2;
      double **d;
      double *d0;
      moveatoms();
      d=(double **)malloc(n_atom*sizeof(double*));
      d0=(double *)malloc(n_atom*2*sizeof(double));
      d[0]=d0;
      for(j=1;j<n_atom;j++)
	d[j]=d[j-1]+2;
      k=0;
      //	printf("q6AA\n");
      for(j=0;j<n_atom;j++)
	if((i!=j)&&(!((a[j].c-at1)*at1)))
	  {
	    double rx,ry,rz;
	    int dd;
	    rx=a[i].r.x-a[j].r.x;
	    ry=a[i].r.y-a[j].r.y;
	    rz=a[i].r.z-a[j].r.z;
	    if(rx<-lx2)rx+=lx;
	    if(rx>lx2)rx-=lx;
	    if(ry<-ly2)ry+=ly;
	    if(ry>ly2)ry-=ly;
	    if(rz<-lz2)rz+=lz;
	    if(rz>lz2)rz-=lz;
	    d[k][0]=rx*rx+ry*ry+rz*rz;
	    d[k][1]=j;
	    k++;
	  }
      //	printf("q6AB\n");
      qsort(d,k,sizeof(double*),compare_i);
      for(j=0;j<maxN;j++)
	{
	  list[j]=(int)d[j][1];
	  //  printf("%lf %d\n",sqrt(d[j][0]),list[j]);
	}
      free(d0);
      free(d);
      if(k<maxN) 
	return -1;
      else
	{
	  double *x;
          double *y;
          double *z;
	  x=(double *)malloc(maxN*sizeof(double));
	  y=(double *)malloc(maxN*sizeof(double));
	  z=(double *)malloc(maxN*sizeof(double));
	  //printf("D\n");
	  neib=0;
	  for(l=0;l<maxN;l++)
	    {
	      double rx,ry,rz;
	      int dd;
	      j=list[l];
	      rx=a[i].r.x-a[j].r.x;
	      ry=a[i].r.y-a[j].r.y;
	      rz=a[i].r.z-a[j].r.z;
	      if(rx<-lx2)rx+=lx;
	      if(rx>lx2)rx-=lx;
	      if(ry<-ly2)ry+=ly;
	      if(ry>ly2)ry-=ly;
	      if(rz<-lz2)rz+=lz;
	      if(rz>lz2)rz-=lz;
	      x[l]=rx;
	      y[l]=ry;
	      z[l]=rz;
	      if( rx*rx+ry*ry+rz*rz<=range*range)neib++;
	      // printf("%lf %lf %lf\n",x[l],y[l],z[l]);
	    }
	  q6[i]=get_q6(maxN,x,y,z);
	  // printf("q6=%d %lf\n",neib,q6[i]);
	  // if((neib==12)&&(q6[i]<0.44)&&(q6[i]>0.43))
	  // for(l=0;l<maxN;l++)
	  // printf("%2d %lf %lf %lf\n",l,x[l],y[l],z[l]);
	  free(x);
	  free(y);
	  free(z);
	  if(neib>=Neib)
	    return maxN;
	  else
	    return -1;
	}
      free(list);
    }
}


void save_density(void)
{
  int j,k,l,nbin,ntot;
  double lx=bound[0].length;
  double ly=bound[1].length;
  double lz=bound[2].length;
  double lx2=lx/2;
  double ly2=ly/2;
  double lz2=lz/2;
  double range=lx2;
  int * neib;  
  if(ly2<range)range=ly2;
  if(lz2<range)range=lz2;
  nbin=(int)(range/rbin);
  moveatoms();
  neib=(int *)malloc(nbin*sizeof(int));
  for(j=0;j<nbin;j++)
    neib[j]=0;;
  for(j=0;j<n_atom;j++)
    if(!((a[j].c-at1)*at1))
      {
	double rx,ry,rz;
	int dd;
	rx=xc-a[j].r.x;
	ry=yc-a[j].r.y;
	rz=zc-a[j].r.z;
	if(rx<-lx2)rx+=lx;
	if(rx>lx2)rx-=lx;
	if(ry<-ly2)ry+=ly;
	if(ry>ly2)ry-=ly;
	if(rz<-lz2)rz+=lz;
	if(rz>lz2)rz-=lz;
	k=(int)(sqrt(rx*rx+ry*ry+rz*rz)/rbin);
	if(k<nbin)
	  neib[k]++;
      }
  ntot=0;
for(j=0;j<nbin;j++)
  {
    double vol=pow((j+1)*rbin,3.0);
    double dvol=vol-pow(j*rbin,3.0);
    vol*=4*pi/3; 
    dvol*=4*pi/3;
    ntot+=neib[j];
    fprintf(density,"%13.5le %12.5le %12.5le %lf\n",nuctime,ntot/vol,neib[j]/dvol,(j+1)*rbin);
  }
  fprintf(density,"& %lf %lf %lf %lf\n",nuctime,xc,yc,zc);
  fflush(density);
  free(neib);
}



void compute_q6(void)
{
  int *clsize;
  int found,max_cls,max_cln;
  int *found_atoms;
  int *found_clusters; 
  int *burning; 
  int **matrix;
  double * q6;
  int i,j,k,l,begin,nburn,nburned,ntot;
  if(good)
    if(get_time()-nuctime>dt_dn)
      {
	nuctime=get_time();
        bound=get_bounds();
	found_clusters=(int *)malloc((n_atom+1)*sizeof(int));
	found_atoms=(int *)malloc((n_atom+1)*sizeof(int));
	burning=(int *)malloc(n_atom*sizeof(int));
	clsize=(int *)malloc(n_atom*sizeof(int));
	//clsize_t=(int *)malloc(n_atom*sizeof(int));
	q6=(double *)malloc(n_atom*sizeof(double));
	matrix=(int **)malloc((n_atom+1)*sizeof(int*));
	matrix[0]=(int *)malloc((maxN+1)*(n_atom+1)*sizeof(int*));
	for(i=1;i<n_atom+1;i++)	
	  matrix[i]=matrix[i-1]+maxN+1;
	
	for(i=0;i<n_atom;i++)
	  {	
	    clsize[i]=0;
	    //clsize_t[i]=0;
	  }
	//printf("q6A\n");
	for(i=0;i<n_atom;i++)
	  matrix[i][0]=find_neib(i,at1,maxN,matrix[i]+1,q6);
	//printf("q6B\n");
	for(i=0;i<n_atom;i++)
	  if(matrix[i][0]>-1)
	    {
	      if((q6[i]<qmin)||(q6[i]>qmax))
		matrix[i][0]=-1;
	      else
		{
		  k=0;
		  for(j=1;j<=matrix[i][0];j++)
		    {
		      l=matrix[i][j];
		      if((q6[l]>=qmin)&&(q6[l]<=qmax))
			{
			  k++;
			  matrix[i][k]=l; 
			}
		    }
		  matrix[i][0]=k;
		}
	    }
	free(q6);
	//printf("q6C\n");
	begin=0;
	ncl=0;
	found=0;
        max_cls=0;
        max_cln=0;
	matrix[n_atom][0]=-1;
	while(begin<n_atom)
	  { 
	    for(i=begin;matrix[i][0]<0;i++)
	      if(i==n_atom)goto done;
	    begin=i;
	    //printf("q6D %d %d %lf %lf %lf\n",i,matrix[i][0],q6[i],qmin,qmax);
	    burning[0]=begin;
	    matrix[i][0]=~matrix[i][0];
	    nburn=1;
	    found_clusters[ncl]=found;
	    nburned=0;
	    ncl++;
	    while (nburn)
	      {
		nburn--;
		j=burning[nburn];
		found_atoms[found]=j;
		found++;
		nburned++;
		for(k=1;k<=~matrix[j][0];k++)
		  {
		    i=matrix[j][k];
		    if(matrix[i][0]>=0)
		      {
			matrix[i][0]=~matrix[i][0];
			burning[nburn]=i;
			nburn++;
			//ntot++;
		      }
		  }
	      } 
	    if(nburned>max_cls)
	      {
		max_cls=nburned;
		max_cln=ncl-1;
	      }
	    clsize[ncl]=nburned;
	    //clsize_t[ncl]=ntot;
	  }
      done:
	found_clusters[ncl]=found;  
	qsort(clsize,ncl+1,sizeof(int),compare1);
	//qsort(clsize_t,ncl+1,sizeof(int),compare1);
	maxcl=clsize[0];
	fprintf(nuc,"%lf %d %d",nuctime,maxcl,ncl);
	for(i=1;i<ncl;i++)
	  fprintf(nuc," %d",clsize[i]);
	fprintf(nuc,"\n");
	fflush(nuc);
	//printf("q6E %d %d %d\n",max_cln,found_clusters[max_cln],found_clusters[max_cln+1] );
	if((nuctime>=start_time)&&(nuctime<=finish_time))
	  {
	    if(ncl)
	    for(i=found_clusters[max_cln];i<found_clusters[max_cln+1];i++)
	      fprintf(maxclust,"%lf %lf %lf\n",a[found_atoms[i]].r.x,a[found_atoms[i]].r.y,a[found_atoms[i]].r.z);
	    fprintf(maxclust,"& %lf\n",nuctime);
	    fflush(maxclust);
	  }
	if((nuctime>=start_time)&&(nuctime<=finish_time)&&rbin)
	  if(ncl)
	    {
	      double lx=bound[0].length;
	      double ly=bound[1].length;
	      double lz=bound[2].length;
	      double lx2=lx/2;
	      double ly2=ly/2;
	      double lz2=lz/2;
	      double sx=0;
	      double sy=0;
	      double sz=0;
	      double x0=0;
	      double y0=0;
	      double z0=0;
	      k=found_clusters[max_cln];
	      j=found_atoms[k];
	      xc=a[j].r.x;
	      yc=a[j].r.y;
	      zc=a[j].r.z;
		for(i=k;i<found_clusters[max_cln+1];i++)
		  {
		    j=found_atoms[i];
		    x0=a[j].r.x;
		    y0=a[j].r.y;
		    z0=a[j].r.z;
		    if(x0-xc>lx2)x0-=lx;
		    if(x0-xc<-lx2)x0+=lx;
		    if(y0-yc>ly2)y0-=ly;
		    if(y0-yc<-ly2)y0+=ly;
		    if(z0-zc>lz2)z0-=lz;
		    if(z0-zc<-lz2)z0+=lz;
		    sx+=x0;
		    sy+=y0;
		    sz+=z0;
		  }
		xc=sx/maxcl;
		yc=sy/maxcl;
		zc=sz/maxcl;
	    }
	if(xc!=-bound[0].length)
	  save_density();
	if(maxcl>=maxcs)set_max_time(0.0);
	free(burning);
	free(clsize);
	free(matrix[0]);
	free(matrix);
        free(found_atoms);
        free(found_clusters);
      }
}

