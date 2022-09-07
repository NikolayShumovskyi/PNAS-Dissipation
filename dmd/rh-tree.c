#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include "rng.h"
typedef struct atom_list{
  int t;
  int c;
  struct atom_list * left;
  struct atom_list * right;
  struct atom_list * parent;
} alist;  
int level;

alist * treeadd(alist * p, alist * npt)
{  
  if(p==NULL)
      p=npt;
  else{
    npt->parent=p;
    if (npt->t<p->t)
      p->left=treeadd(p->left,npt);
    else if (npt->t>p->t)
      p->right=treeadd(p->right,npt);
    else if(rng()<0.5)
      p->left=treeadd(p->left,npt);
    else
      p->right=treeadd(p->right,npt);
  }
  return p;
}

alist * treeleft(alist *p)
{
  if(p->left!=NULL)p=treeleft(p->left);
  return p;
}

alist * treeright(alist *p)
{
  level++;
  if(p->right!=NULL)p=treeright(p->right);
  return p;
}


alist * treekill(alist * root,alist * npt)
{
  int branch=0;
  alist * parent=npt->parent;
  alist * left=npt->left;
  alist * right=npt->right;
  alist *root1=root;
  if(parent!=NULL)
    {
      if(parent->left==npt)
	branch=-1;
      else 
	branch=1;
    }

  if((left==NULL)&&(right==NULL))
    {
      if(branch>0)parent->right=NULL;
      else if(branch<0)parent->left=NULL;
      else root1=NULL;
    }
  else if (left==NULL) 
    {
      if(branch>0)parent->right=right;
      else if(branch<0)parent->left=right;
      else root1=right;
      right->parent=parent;
    }
  else if (right==NULL)
    {
      if(branch>0)parent->right=left;
      else if(branch<0)parent->left=left;
      else root1=left;
      left->parent=parent;
    }
  else if(rng()<0.5)
    {
      alist * right1=treeright(left);
      if(right1==left)
	{
	  if(branch>0)parent->right=left;
	  else if(branch<0)parent->left=left;
	  else root1=left;
	  left->parent=parent;
	  left->right=right;
	  right->parent=left;
	}
      else
	{
	  alist * parent1=right1->parent;
          alist * left1=right1->left;
	  parent1->right=left1;
	  if(left1!=NULL)left1->parent=parent1; 
	  if(branch>0)parent->right=right1;
	  else if(branch<0)parent->left=right1;
	  else root1=right1;
	  right1->parent=parent;
	  right1->right=right;
          right1->left=left;
          left->parent=right1;
          right->parent=right1; 
	}
    }
  else
   {
      alist * left1=treeleft(right);
      if(left1==right)
	{
	  if(branch>0)parent->right=right;
	  else if(branch<0)parent->left=right;
	  else root1=right;
	  right->parent=parent;
	  right->left=left;
          left->parent=right;
	}
      else
	{
	  alist * parent1=left1->parent;
          alist * right1=left1->right;
	  parent1->left=right1;
	  if(right1!=NULL)right1->parent=parent1; 
	  if(branch>0)parent->right=left1;
	  else if(branch<0)parent->left=left1;
	  else root1=left1;
	  left1->parent=parent;
	  left1->right=right;
          left1->left=left;
          left->parent=left1;
          right->parent=left1; 
	}
    }
  return root1;
}

void  treeprint(alist * p)
{ 
  if(p!=NULL)
    {
      treeprint(p->left);
      {
	int cright;
        int cleft;
        int cparent;
	if(p->left!=NULL)
	  cleft=p->left->c;
	else
	  cleft=-1;
	
	if(p->right!=NULL)
	  cright=p->right->c;
	else
	  cright=-1;
	
	if(p->parent!=NULL)
	  cparent=p->parent->c;
	else
	  cparent=-1;
	
	printf("%d %d %d %d %d\n",p->t,p->c,cleft,cright,cparent);
	treeprint(p->right);
      }
    }
}



main()
{
  alist * tree;
  alist *tp;
  alist *root;
  unsigned int rn,rn1;
  double avnrich=0;
  int i,j,k,l,m,n,nup,nrun,irun,nbin,old,dj,jcur,jold,n2,dh;
  double fact,*t,bin;
  int * hh, * rich, * index, * nh ;
  int nrich,step;
  int hth;
  double pth; 
  double *disth;
  int *distj;
  int maxh,r,shift,nbinp;
  double p;
  int mold;
  int nup1;
  double avD=0;
  int nD=0;
  int ip;
  int * dist, *dist1,*dist2;
  FILE *ff;
  char fname[80],fnamet[80];
  printf("out name d?\n");
  scanf("%s",fname);
  printf("out name t?\n");
  scanf("%s",fnamet);
  ff=fopen(fnamet,"w");
  printf("what is n?\n");
  scanf ("%d",&n);
  n2=n/2;
  printf("what is nbin in powers of two?\n");
  scanf ("%d",&nbin);
  shift=32-nbin;
  dh=1<<shift;
  nbin=1<<nbin;
  printf("what is bin for power law distributions?\n");
  scanf ("%lf",&bin);
  bin=1/log(bin);
  nbinp=log(DBL_MAX)*bin;
  printf("what is nrun?\n");
  scanf ("%d",&nrun);
  printf("what is p min?\n");
  scanf ("%lf",&p);
  printf("what is p threshold?\n");
  scanf ("%lf",&pth);
  printf("what is step?\n");
  scanf ("%d",&step);
  ip=p*pow(2.0,31);
  hth=pth*pow(2.0,31);
  printf("%d %d\n",ip,hth);
  printf("what is rn?\n");
  scanf ("%u",&rn);
  rng_init(rn);
  dist=(int *)malloc(nbin*sizeof(int));
  distj=(int *)malloc(n*sizeof(int));
  dist2=(int *)malloc(n*sizeof(int));
  dist1=(int *)malloc(nbin*sizeof(int));
  disth=(double *)malloc(nbin*2*sizeof(double));
  hh=(int *)malloc(n*sizeof(double));
  tree=(alist *)malloc(n*sizeof(alist));
  nh=(int *)malloc(n*sizeof(double));
  rich=(int *)malloc(n*sizeof(double));
  index=(int *)malloc(n*sizeof(double));
  nrich=0;
  /*t=(double *)malloc(nrun*sizeof(double));*/
  root=0;
  for(i=0;i<n;i++)
    {
      hh[i]=0;
      nh[i]=0;
      index[i]=-1;
      distj[i]=0;
      dist2[i]=0;
      tree[i].c=i;
      tree[i].t=0;
      tree[i].left=NULL;
      tree[i].right=NULL;
      tree[i].parent=NULL;
      root=treeadd(root,tree+i);
  }

 

  for(i=0;i<nbin;i++)
    disth[i]=0;
  for(i=0;i<nbin;i++)
    dist[i]=0;
  for(i=0;i<nbinp;i++)
    dist1[i]=0;
  nup=1;

  jold=0;
  jcur=0;
  for(irun=0;1;irun++)
    {
      level=0;
      tp=treeright(root);
      printf("%d %d\n",level,nrich);
      //     printf("& %d %d %d\n", tp->t,tp->c,root->c);
      //treeprint(root);
      maxh=tp->t;
      j=tp->c;     

      if(!nrich)
	{	   
	  if(!nup)
	    {
	      maxh-=dh;
	      k=0;
	      for(i=0;i<n;i++)
	       if(hh[i]>maxh)k++;
             dist2[k]++;
	     /* fprintf(ff,"%d\n",k);*/ 
	     maxh+=dh;
	   }
	}
      else
	{
	  if(!nup)
	    {
	      maxh-=dh;
	      k=0;
	      if(maxh<hth)
		{
		  for(i=0;i<n;i++)
		   if(hh[i]>maxh)k++;
		}
	     else
	       for(i=0;i<nrich;i++)
		 if(hh[rich[i]]>maxh)k++;
	      dist2[k]++;
	      /* fprintf(ff,"%d\n",k);*/ 
	      maxh+=dh;
	   }
	}
      
      dj=j-jold;
      if(dj>n2)dj-=n;
      if(dj<-n2)dj+=n;
      jcur+=dj;
      jold=j;
      distj[abs(dj)]++;
      r=rng_int()>>2;
      hh[j]-=r;
      hh[j]-=r;
      root=treekill(root,tree+j);
      // printf("kill j\n");
      // treeprint(root);
      tree[j].t=hh[j];
      tree[j].parent=NULL;
      tree[j].right=NULL;
      tree[j].left=NULL;
      root=treeadd(root,tree+j);
      //printf("add j\n");
      //treeprint(root);



      if((index[j]>=0)&&(hh[j]<hth))
	{
	  int i1=index[j];
	  index[j]=-1;
	  nrich--;
	  if(i1<nrich)
	    {
	      int i2=rich[nrich];
	      rich[i1]=i2;
	      index[i2]=i1;
	    }
	}

      i=j-1;
      if(i<0)i=n-1;
      if(!nh[i])
	{
	  nh[i]=1;
	  if(nup)
	    nup++;
	  else
	    nup1++; 
	}

      hh[i]+=r; 
      root=treekill(root,tree+i);
      // printf("kill j-1\n");
      //treeprint(root);
      tree[i].t=hh[i];
      tree[i].parent=NULL;
      tree[i].right=NULL;
      tree[i].left=NULL;
      root=treeadd(root,tree+i);
      //printf("add j-1\n");
      //treeprint(root);

      /*      printf("%d %d\n",hh[i],hth); */
      if((hh[i]>=hth)&&(index[i]<0))
	    {
	      rich[nrich]=i;
	      index[i]=nrich;
	      nrich++;
	    }

      i=j+1;
      if(i==n)i=0;
      if(!nh[i])
	{
	  nh[i]=1;
	  if(nup)
	    nup++;
	  else
	    nup1++; 
	}

      hh[i]+=r;

      root=treekill(root,tree+i);
      // printf("kill j+1\n");
      //treeprint(root);
      tree[i].t=hh[i];
      tree[i].parent=NULL;
      tree[i].right=NULL;
      tree[i].left=NULL;
      root=treeadd(root,tree+i);


      if((hh[i]>=hth)&&(index[i]<0))
	{
	  rich[nrich]=i;
	  index[i]=nrich;
	  nrich++;
	}

      if(nup==n)
	{
	  nup=0;
	  //printf("%d %d\n%",irun,nrich);
	  old=irun;
	  m=0;
	  mold=0;
	  nup1=0;
	  for(i=0;i<n;i++)
	    nh[i]=0;
	  /*t[0]=0;*/
	}
      if(!nup)
	{
	  if(nup1==n)
	    {
	      nup1=0;
	      avD+=m-mold;
	      mold=m;
	      //printf("%d %d\n%",irun,nrich);
	      nD++;
	      for(i=0;i<n;i++)
		nh[i]=0;
	    }
	  avnrich+=nrich;
	  if(!(m%step))
	    {
	    for(i=0;i<n;i++)
	      {
		k=((unsigned int)hh[i])>>shift;
		disth[k]++;
	      }
	
	    }
	  k=maxh>>shift;
	  dist[k]++;
	  /* fprintf(ff,"%d\n",dj);*/ 
	  if(maxh<ip)
	      {
		k=irun-old;
		k=log((double)k)*bin;
                dist1[k]++;
		old=irun;
	      }
	  m++;
	  if(m>nrun)break;
	  /* t[m]=t[m-1]+exp(p/maxh);*/
	  /* fprintf(ff,"%d %lf %d\n",irun,maxh,j);*/
	  /* fprintf(ff,"%le %d\n",t,irun);*/
	  
	}
      }
  avnrich/=nrun;
  avD/=nD;
  fclose(ff);
  ff=fopen(fname,"w");
  for(i=0;i<nbin;i++)
    if(disth[i])
      fprintf(ff,"%lf %lf\n",(i<nbin*0.5)?(2.0*i)/nbin:(i-nbin)/(nbin*0.5),disth[i]);
    fprintf(ff,"&\n");  
    for(i=0;i<nbin;i++)
      if(dist[i])
	fprintf(ff,"%lf %d\n",(2.0*i/nbin),dist[i]);
    fprintf(ff,"&\n");
    for(i=0;i<nbin;i++)
      if(dist1[i])
      fprintf(ff,"%d %d\n",i,dist1[i]);
    fprintf(ff,"&\n");
    for(i=0;i<n;i++)
    if(distj[i])
      fprintf(ff,"%d %d\n",i,distj[i]);
    fprintf(ff,"&d\n"); 
    for(i=0;i<n;i++)
      if(dist2[i])
      fprintf(ff,"%d %d\n",i,dist2[i]);
    fprintf(ff,"&d\n"); 
    for(i=0;i<n;i++)
    fprintf(ff,"%d\n",hh[i]);
    printf("%lf %lf\n",avnrich,avD);
    fclose(ff);
}







