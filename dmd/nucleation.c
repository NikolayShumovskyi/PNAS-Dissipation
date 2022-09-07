#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include "nucleation.h"
#include "bcp.h"
#include "search.h"
static int ct=1;
static int minN=10;
static int maxN=12;
static int n_atom;
static int ncl;
static int maxcl;
static int maxcs;
static int good;
static double dt_dn;
static FILE **nuc;
static FILE **maxclust;
static double nuctime;
static double start_time;
static double finish_time;
static moved_atom * a;
int compare(const void *a, const void *b)
{
  return (*(int*)b -*(int*)a);
}

int open_nucleation(char * fname)
{
    FILE * ff;
    int i;
    char fname1[100];
    char dummy[100];
    char dummy1[100];
    if(good)
      {
	for(i=minN;i<=maxN;i++)
	if(nuc[i])fclose(nuc[i]);
        free(nuc);
	for(i=minN;i<=maxN;i++)
	if(maxclust[i])fclose(maxclust[i]);
        free(maxclust);
      }
    good=0;
    ff=fopen(fname,"r");
    if(!ff)return good;
    fscanf(ff,"%s%lf",dummy,&dt_dn);
    fscanf(ff,"%s%d",dummy,&ct);
    fscanf(ff,"%s%d",dummy, &minN);
    fscanf(ff,"%s%d",dummy,&maxN);
    fscanf(ff,"%s%d",dummy,&maxcs);
    fscanf(ff,"%s%s",dummy,dummy1);
    fscanf(ff,"%s%lf",dummy,&start_time);
    fscanf(ff,"%s%lf",dummy,&finish_time);
    fscanf(ff,"%s%s",dummy,dummy);
    fclose(ff);
    if((ct>0)&&(minN>0)&&(maxN>=minN)&&(dt_dn>0)&&(maxcs>0))
      {
	a=(moved_atom *)get_atom();
	n_atom=get_atom_number();
	nuc=(FILE**)malloc((maxN+1)*sizeof(FILE*));
	maxclust=(FILE**)malloc((maxN+1)*sizeof(FILE*));
	if(nuc)
	  for(i=minN;i<=maxN;i++)
	    {
	      sprintf(fname1,"%s-%02d",dummy1,i);
	      nuc[i]=fopen(fname1,"w");
	    }
	if(maxclust)
	for(i=minN;i<=maxN;i++)
	  {
	    sprintf(fname1,"%s-%02d",dummy,i);
	    maxclust[i]=fopen(fname1,"w");
	  }
	nuctime=get_time();
	good=1;
      }
    return good;
}

void compute_nucleation(void)
{
  int *clsize;
  int *clsize_t;
  int *burning;
  int found,max_cls,max_cln,N;
  int *found_atoms;
  int *found_clusters; 
  int **matrix;
  int i,j,k,begin,nburn,nburned,ntot;
  if(good)
    if(get_time()-nuctime>dt_dn)
      {
	nuctime=get_time();
	burning=(int *)malloc(n_atom*sizeof(int));
	clsize=(int *)malloc(n_atom*sizeof(int));
	found_clusters=(int *)malloc((n_atom+1)*sizeof(int));
	found_atoms=(int *)malloc((n_atom+1)*sizeof(int));
	clsize_t=(int *)malloc(n_atom*sizeof(int));
	matrix=(int **)malloc((n_atom+1)*sizeof(int*));
	matrix[0]=(int *)malloc((maxN+1)*(n_atom+1)*sizeof(int*));
	for(i=1;i<n_atom+1;i++)	
	  matrix[i]=matrix[i-1]+maxN+1;
	for(i=0;i<n_atom;i++)
	  matrix[i][0]=find_shell(i,ct,matrix[i]+1);
	moveatoms();
	for(N=minN;N<=maxN;N++)
	  {
	    for(i=0;i<n_atom;i++)
	      {	
		clsize[i]=0;
		clsize_t[i]=0;
	      }    
	    begin=0;
	    ncl=0;
	    found=0;
	    max_cls=0;
	    max_cln=0;
	    matrix[n_atom][0]=-1;
	    while(begin<n_atom)
	      { for(i=begin;matrix[i][0]<N;i++)
		if(i==n_atom)goto done;
	      begin=i;
	      burning[0]=begin;
	      matrix[i][0]=-matrix[i][0];
	      nburn=1;
	      ntot=1;
	      nburned=0;
	      found_clusters[ncl]=found;
	      ncl++;
	      while (nburn)
		{
		  nburn--;
		  j=burning[nburn];
		  found_atoms[found]=j;
		  found++;
		  nburned++;
		  for(k=1;k<=-matrix[j][0];k++)
		    {
		      i=matrix[j][k];
		      if(matrix[i][0]>=N)
			{
			  matrix[i][0]=-matrix[i][0];
			  burning[nburn]=i;
			  nburn++;
			  ntot++;
			}
		      else if(matrix[i][0]>0)
			{
			  matrix[i][0]=-matrix[i][0];
			  ntot++;
			}
		    }
		} 
	      if(nburned>max_cls)
		{
		  max_cls=nburned;
		  max_cln=ncl-1;
		}
	      clsize[ncl]=nburned;
	      clsize_t[ncl]=ntot;
	      }
	  done:
	    for(i=0;i<=n_atom;i++)
	      if(matrix[i][0]<0)
		matrix[i][0]=-matrix[i][0];
	    found_clusters[ncl]=found;
	    qsort(clsize,ncl+1,sizeof(int),compare);
	    qsort(clsize_t,ncl+1,sizeof(int),compare);
	    maxcl=clsize[0];
	    if(maxcl!=max_cls)printf("error_cls\n");
	    fprintf(nuc[N],"%lf %d %d %d",nuctime,maxcl,clsize_t[0],ncl);
	    for(i=1;i<=ncl;i++)
	      fprintf(nuc[N]," %d %d",clsize[i],clsize_t[i]);
	    fprintf(nuc[N],"\n");
	    fflush(nuc[N]);
	    if((nuctime>=start_time)&&(nuctime<=finish_time))
	      {
		if(ncl)
		for(i=found_clusters[max_cln];i<found_clusters[max_cln+1];i++)
		  fprintf(maxclust[N],"%lf %lf %lf\n",a[found_atoms[i]].r.x,a[found_atoms[i]].r.y,a[found_atoms[i]].r.z);
		fprintf(maxclust[N],"& %lf\n",nuctime);
		fflush(maxclust[N]);
	      }
	  }
	if(maxcl>=maxcs)set_max_time(0.0);
	free(burning);
	free(clsize);
	free(clsize_t);
	free(matrix[0]);
	free(matrix);
	free(found_atoms);
	free(found_clusters);
      }  
}

int get_nucleus(void)
{return maxcl;}
