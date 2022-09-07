#include <math.h>
#include <stdio.h>
#include <float.h>
#include "bcp.h"
#include "stdlib.h"
#include "search.h"
#include "controls.h"


static tsearch search;
static atom * a;
static well_type ** ecoll;
static FILE * fp;

void initsearch(void)
{ 
  alist * inpt;
  int i;
  int n=search.n;
  alist * tc=search.storage;
  alist ** ptc=search.begin;
  /*  fp=fopen("search_db","w");*/
  for (i=0;i<n;i++)
    {
    search.collp[i]=-1;
    search.collq[i]=-1;
    search.atom_change[i]=0;
    search.atom_storage[i].n=i;
    search.atom_storage[i].next=NULL;
    }
  search.np=0;
  search.nq=0;
  search.ncha=0;

  for(i=0;i<search.maxfree;i++)ptc[i]=&tc[i];
  ptc[search.maxfree]=NULL;
  search.free=(search.begin)+2*search.n;
  for(i=0;i<2*search.n;i++)
    {
      inpt=search.begin[i];
      inpt->t=DBL2;
      inpt->p=-1;
      inpt->q=-1; 
      inpt->c=-1;
      inpt->pprev=NULL;
      inpt->qprev=NULL;
      inpt->pnext=NULL;
      inpt->qnext=NULL;
    }

  for(i=0;i<n;i++)
    {
      search.atom_collisionp[i]=search.begin[i];
      search.atom_collisionq[i]=search.begin[i+search.n];
    }


  for(i=0;i<=search.z;i++)
    {
    search.cell_atoms[i]=NULL;
    search.cell_champion[i]=NULL;
    } 
}


int realloc_search (void)
{
  int i,j,k,level;
  int maxadd;
  dimensions * bound=get_bounds();
  if ((search.x==bound[0].period)&&(search.y==search.x*bound[1].period)
      &&(search.z==search.y*bound[2].period))return 1;
  search.x=bound[0].period;
  search.y=search.x*bound[1].period;
  if(search.z==search.y*bound[2].period)return 1;
  search.z=search.y*bound[2].period;
  maxadd=search.z;
   k=0;
  j=maxadd;
  level=0;
  do
    { j=(j>>1)+(j&1);
      k+=j;
      level++;
    }while(j>1);
  search.final=level;

if(maxadd>search.maxalloc)
  {
    free(search.cell_atoms);
    free(search.cell_champion);
    if(!(search.cell_atoms=(clist **)malloc((maxadd+1)*sizeof(clist *))))return 0;
    if(!(search.cell_champion=(alist **)malloc((maxadd+1)*sizeof(alist*))))return 0;
    free(search.olymp[0]);
    free(search.olymp);
    search.olymp=(size_al **)malloc((level+1)*sizeof(size_al *));
    if(!(search.olymp))return 0;
    search.olymp[0]=(size_al *)malloc((k+1)*sizeof(size_al));
    if(!(search.olymp[0]))return 0;
    j=maxadd;
    level=0;
    
    while(j>1)
      {
	j=(j>>1)+(j&1);
	search.olymp[level+1]=search.olymp[level]+j;
	level++;
      }
    
    search.change=(size_al *)malloc(46*sizeof(size_al));
    if(!(search.change))return 0;
    search.maxalloc=maxadd; 
  }
  return 1;
}

int allocsearch(int n)
{
 int i,j,k,level;
 int maxadd;
 dimensions * bound=get_bounds();
 a=get_atom();
 ecoll=get_ecoll();
 search.x=bound[0].period;
 search.y=search.x*bound[1].period;
 search.z=search.y*bound[2].period;
 maxadd=search.z;
 search.maxalloc=search.z;
 search.n=n;
 search.nhalf =n>>1;     
 if(maxadd>MAXADD)return 0;
 
 if(!(search.collp=(int *)malloc(n*sizeof(int))))return 0;
 if(!(search.collq=(int *)malloc(n*sizeof(int))))return 0;
 
 if(!(search.atomp=(int *)malloc(n*sizeof(int))))return 0;    
 if(!(search.atomq=(int *)malloc(n*sizeof(int))))return 0;

 if(!(search.atom_change=(int *)malloc(n*sizeof(int))))return 0;    
 if(!(search.changed_atoms=(int *)malloc(n*sizeof(int))))return 0;
 
 if(!(search.cell_atoms=(clist **)malloc((maxadd+1)*sizeof(clist *))))return 0;
 if(!(search.cell_champion=(alist **)malloc((maxadd+1)*sizeof(alist*))))return 0;
 if(!(search.atom_storage=(clist *)malloc(n*sizeof(clist))))return 0;

 search.maxfree=NFREE;
 if(!(search.storage =(alist *)malloc(search.maxfree*sizeof(alist))))return 0;
 if(!(search.begin=(alist **)malloc((search.maxfree+1)*sizeof(alist *))))return 0;
 if(!(search.atom_collisionp=(alist **)malloc((n)*sizeof(alist *))))return 0;
 if(!(search.atom_collisionq=(alist **)malloc((n)*sizeof(alist *))))return 0;
 
 k=0;
 j=maxadd;
 level=0;
 
 do
 { j=(j>>1)+(j&1);
   k+=j;
   level++;
  }while(j>1);
  
  search.final=level;
  search.olymp=(size_al **)malloc((level+1)*sizeof(size_al *));
  if(!(search.olymp))return 0;
  search.olymp[0]=(size_al *)malloc((k+1)*sizeof(size_al));
    if(!(search.olymp[0]))return 0;
  j=maxadd;
  level=0;
  
  while(j>1)
  {
   j=(j>>1)+(j&1);
   search.olymp[level+1]=search.olymp[level]+j;
   level++;
   }
  
   search.change=(size_al *)malloc(46*sizeof(size_al));
    if(!(search.change))return 0;
 return 1;
}
/* used in find_atoms
 when atom moves from one cell to another */
int find_collisions(int p1)
{
  int np=0;
  alist * pt=search.atom_collisionp[p1];
  *(--search.free)=pt;
  pt=pt->pnext;
  pt->pprev=NULL;
  search.atom_collisionp[p1]=pt;
  if(!search.atom_change[p1])
    {
      search.atom_change[p1]=1;

      search.changed_atoms[search.ncha++]=p1;
    }
  while(pt->p==p1)
    {
      int q=pt->q;
      if(search.collp[q]==-1)search.atomp[np++]=q;
      search.collp[q]=-2;
      pt=pt->pnext;
    }
  pt=search.atom_collisionq[p1];
  while(pt->q==p1)
    {
      int q=pt->p;
      if(search.collp[q]==-1)search.atomp[np++]=q;
      search.collp[q]=-2;
      pt=pt->qnext;
    }
  return np;
}

/* used in squeeze1 */ 
int kill_collisions(int p1, int * collp, int * atomp)
{
  int np=0;
  alist ** free=search.free;
  alist * pt=search.atom_collisionp[p1];
  int ncha=search.ncha;
  int *changed_atoms=search.changed_atoms;
  int *atom_change=search.atom_change;
  if(!atom_change[p1])
    {
      atom_change[p1]=1;
      changed_atoms[ncha++]=p1;
    } 
  while(pt->p==p1)
    {
      alist * qprev=pt->qprev;
      int q=pt->q;
      if(q<search.n)
	{
	  if(collp[q]==-1)
	    {
	      collp[q]=pt->c;
	      atomp[np++]=q;
	    }
	  if(qprev==NULL)
	    {
	      search.atom_collisionq[q]=pt->qnext;
	      search.atom_collisionq[q]->qprev=NULL;
	    }
	  else
	    { 
	      qprev->qnext=pt->qnext;
	      pt->qnext->qprev=qprev;
	    }
	}
      *(--free)=pt;
      pt=pt->pnext;
      pt->pprev=NULL;
    }
  search.atom_collisionp[p1]=pt;
  pt=search.atom_collisionq[p1];
  while(pt->q==p1)
    {
      alist * pprev=pt->pprev;
      int p=pt->p;
      if(collp[p]==-1)
	{
	  collp[p]=pt->c;
          atomp[np++]=p;
	}
	  if(pprev==NULL)
	    {
	      search.atom_collisionp[p]=pt->pnext;
	      search.atom_collisionp[p]->pprev=NULL;
	      if(!atom_change[p])
		{
		  atom_change[p]=1;
		  changed_atoms[ncha++]=p;
		} 
	    }
	  else
	    { 
	      pprev->pnext=pt->pnext;
	      pt->pnext->pprev=pprev;
	    }
      *(--free)=pt;
      pt=pt->qnext;
      pt->qprev=NULL;
    }
  search.atom_collisionq[p1]=pt;
  search.ncha=ncha;
  search.free=free;
  return np;
}

/* used in new_loc in bcp.c */
void find_atoms(int p1, size_al address1)
{ 
  clist * inpt;
  int np;
  int nch=search.nch;
  int i,j,k,i1,j1,k1,i2,j2,k2,address=a[p1].add;
  int addressz,addressy;
  clist * old_rec=search.cell_atoms[address];
  if(old_rec->n==p1)
    search.cell_atoms[address]=search.cell_atoms[address]->next;
  else
    {
      clist *  prev_rec=old_rec;
      old_rec=prev_rec->next;
      while(old_rec->n!=p1)
	{
	  prev_rec=old_rec;

	  old_rec=prev_rec->next;
	}
      prev_rec->next=old_rec->next;
    }
  old_rec->next=search.cell_atoms[address1];
  search.cell_atoms[address1]=old_rec;
  a[p1].add=address1;
  search.cell_champion[address]=NULL;
  search.change[nch]=address;
  nch++;
  search.cell_champion[address1]=NULL;
  search.change[nch]=address1;
  nch++;
  search.nch=nch;
  np=find_collisions(p1);  

  i1=a[p1].i.x.i-1;
  j1=(a[p1].i.y.i-1)*search.x;
  i2=i1+2;
  j2=j1+(search.x<<1);
  if(search.z==search.y)
  {k1=a[p1].i.z.i*search.y;k2=k1;}
  else
  {k1=(a[p1].i.z.i-1)*search.y;k2=k1+(search.y<<1);}
  
  for(k=k1;k<=k2;k+=search.y)
    { 
      addressz=k;
      if(addressz<0)addressz+=search.z;
      if(addressz==search.z)addressz=0;
      for(j=j1;j<=j2;j+=search.x)
	{ 
	  addressy=j; 
	  if(addressy<0)addressy+=search.y;
	  if(addressy==search.y)addressy=0;
	  addressy+=addressz;
	  for(i=i1;i<=i2;i++)
	    { 
	      address=i; 
	      if(address<0)address+=search.x;
	      if(address==search.x)address=0;
	      address+=addressy;
	      for(inpt=search.cell_atoms[address];inpt;inpt=inpt->next)
	      {
		int p2=inpt->n;
		if((search.collp[p2]==-1)&&(p1!=p2))
		  { 
		    search.atomp[np++]=p2;
                    search.collp[p2]=ecoll[a[p1].c][a[p2].c];
		  }

	      }
	    } 

	}
    }  

 search.nq=0;
 search.np=np;
}

/* find atoms in ajacent cells to a given atom; 
   returns number of such atoms */
int find_neighbors(int p, int * neib)
{ 
 clist * inpt;
 int nn=0;
 int i,j,k,i1,j1,k1,i2,j2,k2,address;
 int addressz,addressy;
  
  i1=a[p].i.x.i-1;
  j1=(a[p].i.y.i-1)*search.x;
  i2=i1+2;
  j2=j1+(search.x<<1);
  if(search.z==search.y)
  {k1=a[p].i.z.i*search.y;k2=k1;}
  else
  {k1=(a[p].i.z.i-1)*search.y;k2=k1+(search.y<<1);}
  
  for(k=k1;k<=k2;k+=search.y)
    { 
      addressz=k;
      if(addressz<0)addressz+=search.z;
      if(addressz==search.z)addressz=0;
      for(j=j1;j<=j2;j+=search.x)
	{ 
	  addressy=j; 
	  if(addressy<0)addressy+=search.y;
	  if(addressy==search.y)addressy=0;
	  addressy+=addressz;
	  for(i=i1;i<=i2;i++)
	    { 
	      address=i; 
	      if(address<0)address+=search.x;
	      if(address==search.x)address=0;
	      address+=addressy;
	      for(inpt=search.cell_atoms[address];inpt;inpt=inpt->next)
		    neib[nn++]=inpt->n;
	    } 
	}  
    }
 return nn;
}  
/* list atoms in a certain cell, returns number of atoms in a cell */
int list_atoms(size_al address,int * atomx)
{ 
	int i=0;
	clist *inpt1;
	for( inpt1=search.cell_atoms[address];inpt1;inpt1=inpt1->next)
 	{ atomx[i++]=inpt1->n;}
 return i;	
}

int bond(int i, int j)
{   
    alist *inpt1;
	for( inpt1=search.atom_collisionp[i];inpt1->p==i;inpt1=inpt1->pnext)
	if(j==inpt1->q) return inpt1->c;
	for( inpt1=search.atom_collisionq[i];inpt1->q==i;inpt1=inpt1->qnext)
	if(j==inpt1->p) return inpt1->c;
	return -ecoll[a[i].c][a[j].c];
}

void add_atom2cell(int i)
{
  size_al address=a[i].add;
  clist * pt=search.cell_atoms[address];
  search.cell_atoms[address]=&search.atom_storage[i];
  search.cell_atoms[address]->next=pt;
}

alist * local_champion(size_al address)
{
  double t;
  int i;
  clist *  inpt=search.cell_atoms[address];
  if(!inpt)return *search.begin;
  i=inpt->n;
  t=search.atom_collisionp[i]->t;
  inpt=inpt->next;
  while(inpt)
    {
      int j=inpt->n;
      double t1=search.atom_collisionp[j]->t;
      if(t1<t)
	{
	  t=t1;
	  i=j;
	}
      inpt=inpt->next;
    }
  return search.atom_collisionp[i];
}

int init_tables(void)
{
int i0,j0,k0,address,ix,iy,iz,level;
  int i,j,k,i1,j1,k1,i2,j2,k2;
  int maxadd=search.z;
  size_al address1;
  size_al address2;
  int n=search.n-1;
  int n1=n+1;
  double t;
  int p,q,ct;
   
  int addressz,addressy; 
  
  initsearch();
  for(i=0;i<=n;i++)
    {
      add_atom2cell(i);
      q=twall(i,&t);
      bubble(t,i,q,-1);
    }


  a[n1].c=0;
  for (k0=0;k0<maxadd;k0++)
  if(search.np=list_atoms(k0,search.atomp))
  {
   i0=search.atomp[0];
   i1=a[i0].i.x.i-1;
   i2=i1+2;
   j1=(a[i0].i.y.i-1)*search.x;
   j2=j1+(search.x<<1);
  if(search.z==search.y)
  {k1=a[i0].i.z.i*search.y;k2=k1;}
  else
  {k1=(a[i0].i.z.i-1)*search.y;k2=k1+(search.y<<1);}
   for(k=k1;k<=k2;k+=search.y)
   { 
     addressz=k;
     if(addressz<0)addressz+=search.z;
     if(addressz==search.z)addressz=0;
     for(j=j1;j<=j2;j+=search.x)
	 { 
	  addressy=j; 
	  if(addressy<0)addressy+=search.y;
	  if(addressy==search.y)addressy=0;
	  addressy+=addressz;
	  for(i=i1;i<=i2;i++)
	  { 
	   int iq;
	   int ip;
	   address=i; 
	   if(address<0)address+=search.x;
	   if(address==search.x)address=0;
	   address+=addressy;
 	   search.nq=list_atoms(address,search.atomq);
 	   for (iq=0;iq<search.nq;iq++)
 	   for (ip=0;ip<search.np;ip++)
 	   {
  	    i0=search.atomp[ip]; 
  	    j0=search.atomq[iq];
	    if(i0<j0)
	    {

	     int ct=collision_type(i0,j0);
	     if(ct<0) return (int) ct;
	     add_potential(ct);
	     if(tball(i0,j0,ct,&t)) 
	       {
		 p=i0;
		 q=j0;
	       }
	     else
	       {
		 p=-1;
		 q=-1;
		 ct=-1;
	       }  
	     bubble(t,p,q,ct);
	    }
	   }
	  }
	 }
   }
  }

  for (k0=0;k0<maxadd;k0++)
      search.cell_champion[k0]=local_champion(k0);

  k=maxadd;
  k=(k>>1)+(k&1);
  for(i=0;i<k;i++)
    {
      address1=i<<1;
      address2=address1+1;
      search.olymp[0][i]=address1;
      if((address2<maxadd)&&(search.cell_champion[address1]->t
			     >search.cell_champion[address2]->t))
	search.olymp[0][i]=address2;
    }

  for(level=1;level<search.final;level++)
    { j=(k>>1)+(k&1);
    for(i=0;i<j;i++)
      {
	size_al i1=i<<1;
	size_al i2=i1+1;
	search.olymp[level][i]=search.olymp[level-1][i1];
	if(i2<k)
	  {
	    address1=search.olymp[level-1][i1];
	    address2=search.olymp[level-1][i2];
   	if(search.cell_champion[address1]->t
	   >search.cell_champion[address2]->t) 
	  search.olymp[level][i]=address2;
	  }
      }
    k=j;
    }  
  return 1;  
}

int bubble (double t, int p1, int q1, int c)
{
  int p,q;
  alist *npt,*pt,*inpt;
  if(q1<p1)
    {
      p=q1;
      q=p1;
    }
  else
    {
      p=p1;
      q=q1;
    }
  if((q<search.n)&&(q-p>search.nhalf))
    {
      int r=q;
      q=p;
      p=r;
    }

  if (p>=0)
    {
      npt=*(search.free);
      if(!npt){writetext(get_text_name());exit(0);}
      search.free++;
      npt->t=t;
      npt->p=p;
      npt->q=q;
      npt->c=c;
      if(q<search.n)
	{
	  inpt=search.atom_collisionq[q];	
	  npt->qnext=inpt;
	  inpt->qprev=npt;          
	  npt->qprev=NULL;
	  search.atom_collisionq[q]=npt;
	}
      inpt=search.atom_collisionp[p];
      if(inpt->t >= t)
	{
	  npt->pnext=inpt;
	  inpt->pprev=npt;          
	  npt->pprev=NULL;
	  search.atom_collisionp[p]=npt;
	  return p;
	}
      else 
	{
	  pt=inpt->pnext;
	  while ( pt->t<t)
	    pt=pt->pnext;
	  
         pt->pprev->pnext=npt;
         npt->pprev=pt->pprev;
         npt->pnext=pt;
         pt->pprev=npt;
	}
    }
    return -1;
}


void olymp_sort(void)
{
  size_al *change=search.change;
  size_al * curr, *next, *prev;
  size_al i,j,k,address1,address2;
  int ch,nch1,nch=search.nch;
  int level;
  nch1=0;
  curr=search.olymp[0];
  next=search.olymp[1];
  for(ch=0;ch<nch;ch++)
   {
    i=change[ch];
    search.cell_champion[i]=local_champion(i);
    address1=i>>1;     
    if(curr[address1]<MAXADD)
	{
	  curr[address1]=MAXADD;
	  change[nch1]=address1;
	  nch1++;
	}
   }
  nch=nch1;
  nch1=0;

  for(ch=0;ch<nch;ch++)
   {
    i=change[ch];
    address1=i<<1;
    address2=address1+1;
    curr[i]=address1;
    if((address2<search.z)&&(search.cell_champion[address1]->t>search.cell_champion[address2]->t))
      curr[i]=address2;
    j=i>>1;
    if(next[j]<MAXADD)
    {
      next[j]=MAXADD;
      change[nch1]=j;
      nch1++;
    } 
   }
 nch=nch1;
 k=search.z; 
 for(level=2;level<search.final;level++) 
  { 
   k=(k>>1)+(k&1);
   nch1=0;
   prev=curr;
   curr=next;
   next=search.olymp[level];
   for(ch=0;ch<nch;ch++)
    { 
     i=change[ch];
     address1=i<<1;
     address2=address1+1;
     address1=prev[address1];
     curr[i]=address1;
     if(address2<k)
      { 
      	address2=prev[address2];
        if(search.cell_champion[address1]->t>search.cell_champion[address2]->t)
  	      curr[i]=address2;
  	   }
     j=i>>1;
     if(next[j]<MAXADD)
      {
        next[j]=MAXADD;
        change[nch1]=j;
        nch1++;
      } 
       
    }  
   nch=nch1;
 } 
 next[0]=(search.cell_champion[curr[0]]->t>search.cell_champion[curr[1]]->t) ? curr[1]:curr[0];
 
 return;  
}

int get_free(void)
{
return (int)(search.free-search.begin);
}
int get_maxfree(void)
{
return search.maxfree;
}
void set_maxfree(int a)
{
search.maxfree=a;
}


void update_table(int p1, int q1, int ct1)
{
  int i,p2,ct;
  int add;
  int ncha=search.ncha;
  double t;
  int * changed=search.changed_atoms;
  int * change=search.atom_change;
  int p=p1;
  int q=twall(p1,&t);
  if((p=bubble(t,p,q,-1))>=0)
    { 
      if(!change[p])
	{ 
	  change[p]=1;
	  changed[ncha++]=p;
	}
    }
  search.collp[p1]=-1;    
  if (q1<search.n)
    {
      search.collp[q1]=ct1; /*to make sure that we compute 
			      the collision of atoms p1 and q1, when we will look through atomq 
			      which is different from the old type ct1 in the list collq */
      search.collq[p1]=-1;
      search.collp[p1]=-1;
      search.collq[q1]=-1;
      p=q1;
      ct=-1;
      q=twall(q1,&t);
      if((p=bubble(t,p,q,-1))>=0)
	{ 
	  if(!change[p])
	    { 
	      change[p]=1;
	      changed[ncha++]=p;
	    }
	}   
    }
  for(i=0;i<search.np;i++)
    {
      p2=search.atomp[i];
      ct=search.collp[p2];
      search.collp[p2]=-1; 
      if(ct>=0) /* the collision is not recalculated if ct is -2
		   which means that the atoms will collide as before */
	if(tball(p2,p1,ct,&t))
	  {  
	    if((p=bubble(t,p1,p2,ct))>=0)
              { 
		if(!change[p])
		  { 
		    change[p]=1;
		    changed[ncha++]=p;
		  }   
	      }
	  }
    }
  if (q1<search.n)
    for(i=0;i<search.nq;i++)
      {
	p2=search.atomq[i];
	ct=search.collq[p2];
	search.collq[p2]=-1; 
	if(ct>=0)
	  if(tball(p2,q1,ct,&t))
	    {
	      if((p=bubble(t,q1,p2,ct))>=0)
		{ 
		  if(!change[p])
		    { 
		      change[p]=1;
		      changed[ncha++]=p;
		    }  
		}   
	    }
      }
  for(i=0;i<ncha;i++)
    {
      int j=changed[i];
      if(change[j])
	{
	  change[j]=0;
	  if(search.cell_champion[a[j].add])
	    {
	      search.change[search.nch++]=a[j].add;
	      search.cell_champion[a[j].add]=NULL;
	    }
	}
    }
  search.ncha=0;
/*  if((p1==491))
    print_collisions(p1,q1);
  if((q1==491))
    print_collisions(q1,p1);
  if((p1==990))
    print_collisions(p1,q1);
  if((q1==990))
    print_collisions(q1,p1);
*/

  olymp_sort();
}



int squeeze1(int p1, int * collp, int * atomp)
{ 
 clist *inpt;
 int i,j,k,i1,j1,k1,i2,j2,k2,address;
 int addressz,addressy;
 int np=kill_collisions(p1,collp,atomp);
 i1=a[p1].i.x.i-1;
 j1=(a[p1].i.y.i-1)*search.x;
 i2=i1+2;
 j2=j1+(search.x<<1);
 
 if(search.z==search.y){
   k1=a[p1].i.z.i*search.y;
   k2=k1;
 }
 else{
   k1=(a[p1].i.z.i-1)*search.y;
   k2=k1+(search.y<<1);
 }
 
 for(k=k1;k<=k2;k+=search.y){ 
   addressz=k;
   if(addressz<0)
     addressz+=search.z;
   if(addressz==search.z)
     addressz=0;
   
   for(j=j1;j<=j2;j+=search.x){ 
     addressy=j; 
     if(addressy<0)
       addressy+=search.y;
     if(addressy==search.y)
       addressy=0;
     addressy+=addressz;
     
     for(i=i1;i<=i2;i++){ 
       address=i; 
       if(address<0)
	 address+=search.x;
       if(address==search.x)
	 address=0;
       address+=addressy;
       
       inpt=search.cell_atoms[address];
       while(inpt)
	 { 
	   int p=inpt->n;
	   if(collp[p]==-1)
	     {
	       atomp[np++]=p;
	       collp[p]=ecoll[a[p].c][a[p1].c];
	     }  
	   inpt=inpt->next;
	 }	   	
     }
   }
 }  
 return np;
}

int squeeze_table(int * p1, int * q1, double * timea)
{ 
  alist *inpt=search.cell_champion[search.olymp[search.final-1][0]];
  search.nch=0;
  search.ncha=0;
  *timea=inpt->t;
  *p1=inpt->p;
  *q1=inpt->q;
/*  if((*p1==491))
    print_collisions(*p1,*q1);
  if((*q1==491))
    print_collisions(*q1,*p1);
  if((*p1==990))
    print_collisions(*p1,*q1);
  if((*q1==990))
    print_collisions(*q1,*p1);
*/
  if((*q1)>=search.n)
  {
   return (int)(inpt->c);
  }
 search.np=squeeze1(*p1,search.collp, search.atomp);
 search.nq=squeeze1(*q1,search.collq, search.atomq);
  return (int)(inpt->c);
}

  
int * get_collp(void){return search.collp;}
int * get_collq(void){return search.collq;}
int * get_atomp(void){return search.atomp;}
int * get_atomq(void){return search.atomq;}
int get_np(void){return search.np;}
int get_nq(void){return search.nq;}

int pairs(int (*do_something)(int,int,int))
{ 
  int n_bonds=0;
  int i;
  alist *inpt1;
  /*  for (i=0;i<search.n;i++)
    for( inpt1=search.atom_collisionp[i];inpt1->p==i;inpt1=inpt1->pnext)
      if(inpt1->q<search.n)
      n_bonds+=do_something(inpt1->p,inpt1->q,inpt1->c);*/
  for (inpt1=search.begin[1];inpt1<*search.free;inpt1++)
    if(inpt1->q<search.n)
      n_bonds+=do_something(inpt1->p,inpt1->q,inpt1->c);
  return n_bonds;
}  

void print_collisions(int p1,int q1)
{
  int np=0;
  alist * pt=search.atom_collisionp[p1];
      fprintf(fp,"collision %lf %d %d %d\n",get_time(),get_ll(),p1,q1);
  while(pt->p==p1)
    {
      fprintf(fp,"%10.6lf %3d %3d %2d\n",pt->t,pt->p,pt->q,pt->c);
      pt=pt->pnext;
    }
  pt=search.atom_collisionq[p1];
  fprintf(fp,"--------\n");
  while(pt->q==p1)
    {
      fprintf(fp,"%10.6lf %3d %3d %2d\n",pt->t,pt->p,pt->q,pt->c);
      pt=pt->qnext;
    }
fflush(fp);
}

int find_shell(int p1, int ct, int * neib)
{ 
  
  int nn=0;
  alist * pt=search.atom_collisionp[p1];
  while(pt)
    {
      if(pt->c==ct)
	{
	  neib[nn]=pt->q;
	  nn++;
	}
      pt=pt->pnext;
    }
  pt=search.atom_collisionq[p1];
  while(pt)
    {
      if(pt->c==ct)
	{
	  neib[nn]=pt->p;
	  nn++;
	}
      pt=pt->qnext;
    }
  return nn;
}
