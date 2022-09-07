//Metropolis method for Ising model on square lattice
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define mult (314159261)
#define mult1 (314159269)
#define add (907633385)
#define big 4294967296.0
//saves bmp picture 
void savebmp(char * fname, short * s,int step){ 
  unsigned char * color;
  int bmp[269];
  FILE *ff;
  int i,lx,ly,ns;
  lx=step;
  ly=step;
  ns=step*step;
  ff=fopen(fname,"wb");
  color=(unsigned char*) malloc(2*sizeof(unsigned char));
  color[0]='B';
  color[1]='M';
  fwrite(&color[0],1,2,ff);
  for(i=0;i<269;i++)
    bmp[i]=0;
  bmp[0]=lx*ly*3+1078;
  bmp[2]=1078;
  bmp[3]=40;
  bmp[4]=lx;
  bmp[5]=ly;
  bmp[6]=524289;
  bmp[8]=lx*ly*3;
  bmp[9]=2834;
  bmp[10]=2834;
  bmp[11]=256;
  bmp[12]=256;
    for(i=0;i<256;i++)
      bmp[13+i]=i*65793;
  color[0]=255;
  color[1]=0;
  fwrite(&bmp[0],4,269,ff);
  for(i=0;i<ns;i++)
    fwrite(&color[(s[i]+1)>>1],1,1,ff);
  fclose(ff);
  free(color);
}

void savebit(char * fname, short * s,int step){ 
  unsigned char * color;
  FILE *ff;
  int i,j,k,l,n,m;
  ff=fopen(fname,"wb");
  n=step>>3;
  color=(unsigned char *)malloc(n*sizeof(unsigned char));
  k=0;
  l=0;
  
  for(i=0;i<n;i++)
    {
      color[i]=0;
      for(j=0;j<8;j++)
	{
	  if(s[k]>0)
	    {
	      color[i]|=1<<j;
	      l++;
	    }
	  k++;
      }
    }
  // printf("sum %d\n",l);
  //for(i=0;i<n;i++)
  //  printf("%d %d\n",i, (int)color[i]);

  m=fwrite(color,1,n,ff);
  fclose(ff);
  free(color);
}


void main()
{
  unsigned int rn,rn1; //random number generator seeds
  int i,j,k,l,m,l1,j1; //dummy interger indexes
  int j2,flip;
  int nrun,irun; //number of runs,run id.
  int n; //how many Metropolis steps are taken between averages
  double fn; //fn=1/n
  int step; //lattice size L
  int step1,step2,step3;//number of points in a lattice, LxL
  int mag,mag1; //magnetization
  int en,en1; //energy
  double U,Ch;
  double time,dt;
  int de=0; //change in energy
  double avmag,aven,amag; //average absolute value of magnetization, av. energy
  //average magnetization
  int inde; //index of Metropolis probabilities exp(\Delta E/kT)
  double p; //Metropolis p; 
  double H,H0,dH; //current magnetic field, initial H, step in H
  double T,T0,dT;//current temperature, initial T, step in T
  int iH,nH; //index of a mag. field, number of different M.f. +1
  int iT,nT; //index  of a temperature, number of different temperatures+1
  short * s; //array of sipns
  int *target;
  int *pairlist;
  int *index;
  int npairs;
  double *dist_m,bin,bin1;
  int nbin,nbin2;
  FILE *ff,*ff1; //output file
  char fname[80],fnamebmp[80],fname1[80]; //name of files
  int dz[6]; //change of addresses for nearest neighbors 
  double ex[8]; //array of Metropolis probabilities
  double kaw[20];
  double pr,E,fr;
  double fact,fact1,fact6;//factors for getting r.n. in a given range 
  double tsm,bst,ts,tc;
  int lc,nc;
  printf("out name?\n"); //input
  scanf("%s",fname);
  ff=fopen(fname,"w");  
  printf("what is L?\n");
  scanf ("%d",&step);
  printf("what is nrun?\n");
  scanf ("%d",&nrun);
  printf("what is average length?\n");
  scanf ("%d",&n);
  /* n tells how many MC steps is taken between two outputs */
  fn=((double)1)/n;
  printf("what is T0?\n");
  scanf ("%lf",&T0);
  printf("what is dT\n");
  scanf ("%lf",&dT);
  printf("what is nT?\n");
  scanf ("%d",&nT);
  printf("what is H0?\n");
  scanf ("%lf",&H0);
  printf("what is dH\n");
  scanf ("%lf",&dH);
  printf("what is nH?\n");
  scanf ("%d",&nH);
  printf("what is n bins?\n");
  scanf ("%d",&nbin);
  printf("what is probability of reaction?\n");
  scanf ("%lf",&pr); 
  printf("what is energy of reaction?\n");
  scanf ("%lf",&E);
  printf("initial fraction of spins up?\n");
  scanf ("%lf",&fr);
  printf("minimal ubdate time?\n");
  scanf ("%lf",&tsm);
  printf("logarithmic base of saving time?\n");
  scanf ("%lf",&bst);
  printf("cycle length for saving time?\n");
  scanf ("%d",&lc);
  tc=tsm*pow(bst,(double)lc);
    printf("what is rn?\n");
  scanf ("%d",&rn); //end of inputs

  /* step1 is the total number of spins */ 
  step1=step*step;
  step2=step1*step;
  step3=step2*3;
  /*dz is the array of address differences of the nearest neighbors*/
  dz[0]=1;
  dz[1]=step;
  dz[2]=step1;
  dz[3]=-1;
  dz[4]=-step;
  dz[5]=-step1;
  /* array of spins */
  s=(short *)malloc(step2*sizeof(short));
  target=(int *)malloc(step2*sizeof(int));
  index=(int *)malloc(step3*sizeof(int));
  pairlist=(int *)malloc(step3*sizeof(int));
  nbin2=2*nbin+1;
  bin =step2/(double)nbin;
  bin1=1/bin;
  dist_m=(double *)malloc(nbin2*sizeof(double));
  /* rundom number generator coefficients */
  fact= ((double)step2)/big; //factor for selecting random spin;
  fact6= ((double)6)/big; //factor for selecting random spin;
  fact1= ((double)1)/big; //factor for selecting selecting p; 


  for(iH=0;iH<=nH;iH++)//loop in magnetic fields
    {
      H=H0+dH*iH;
      for(iT=0;iT<=nT;iT++) //loop in temperatures
	{
	  T=T0+dT*iT;
	  
	  /* probabilities to increase the potential energy */
	  /* we assume that -2<H<2, thus if energy change due to spin interactions
	     is negative, the gross change in energy is also negative */
	  //odd inceces correspond to switching from negative to positive orientation 
	  for(i=0;i<nbin2;i++)
	    dist_m[i]=0;
	  ex[0]=exp(-2*H/T+E/T);
	  ex[2]=exp(-4/T-2*H/T+E/T);
	  ex[4]=exp(-8/T-2*H/T+E/T);
	  ex[6]=exp(-12/T-2*H/T+E/T);
	  ex[1]=exp(2*H/T+E/T);
	  ex[3]=exp(-4/T+2*H/T+E/T);
	  ex[5]=exp(-8/T+2*H/T+E/T);
	  ex[7]=exp(-12/T+2*H/T+E/T);
	  for(i=0;i<=10;i++)
	    kaw[i]=exp(-(2*i)/T);
	  /* spin initialization */
	  
	  for(i=0;i<step2;i++)
	    {
	      s[i]=-1;
	      target[i]=i;
	    }
	  m=step2*fr;
	  l=step2;
	  for(i=0;i<m;i++)
	    {
	      rn=rn*mult+add;
	      k=l*(fact1*rn);
	      j=target[k];
	      l--;
	      target[k]=target[l];
	      target[l]=j;
	      s[j]=1;
	    }
	  
	  
	  
	  /* en -energy, mag- magnetization initialization */
	  en=0;
	  mag=0;
	  for(i=0;i<step2;i++)
	    {
	      int eni=0; 
	      for(k=0;k<6;k++)
		{
		  j=i+dz[k];
		  if(j>=step2)j-=step2;
		  if(j<0)j+=step2;
		  eni+=s[j];
		}
	      if(s[i]<0)
		en-=eni;
	      else 
		en+=eni;
	      mag+=s[i];
	    }
	  en>>=1;
	  npairs=0;
	  for(i=0;i<step2;i++)
	    for(k=0;k<3;k++)
	      {
		j=i+dz[k];
		if(j>=step2)j-=step2;
		if(s[i]!=s[j])
		  {
		    l=3*i+k;
		    pairlist[npairs]=l;
		    index[l]=npairs;
		    npairs++;
		  }
	      }
	      time=0;
	      ts=tsm;
	      nc=0;
	      for(irun=0;irun<nrun;irun++)
	    {
	      aven=0;
	      avmag=0;
	      amag=0;
	      Ch=0;

	      for(m=0;m<n;m++)
		{
		  dt=step3/(double)npairs;
		  rn=rn*mult+add;
		  if(rn*fact1<pr*dt)
		    {  
		      de=0;
		      flip=0;
		      rn=rn*mult+add;
		      i=rn*fact;
		      /* select at random i-th spin to flip */ 
		      for(k=0;k<6;k++)
			{
			  j=i+dz[k];
			  if(j>=step2)j-=step2;//taking into account
			  if(j<0)j+=step2;//periodic boundaries
			  de+=s[j];
			}
		      /* de is the magnetization of neiboring spins */
		      inde=de;
		      if(s[i]<0){de=-de;inde=-inde+1;}
		      /* de is the energy of spin interaction */
		      /*inde is the index of the Boltzmann factor */
		      //printf("%d\n",inde);
		      if(inde<0)
			flip=1; //accept spin flip
		      else
			{
			  rn=rn*mult+add; //apply Metropolis criyterion;
			  p=rn*fact1;
			  if(p<ex[inde])
			    flip=1; //accept spin flip
			}
		      if(flip)
			{
			  for(k=0;k<3;k++)
			    {
			      l=i*3+k;
			      j=i+dz[k];
			      if(j>=step2)j-=step2;//taking into account
			      if(s[i]!=s[j])
				{
				  npairs--;
				  l1=pairlist[npairs];
				  j1=index[l];
				  pairlist[j1]=l1;
				  index[l1]=j1;
				}
			      else
				{
				  pairlist[npairs]=l;
				  index[l]=npairs;
				  npairs++;
				}
			    }
			  for(k=0;k<3;k++)
			    {
			      j=i-dz[k];
			      if(j<0)j+=step2;//taking into account
			      l=j*3+k;
			      if(s[i]!=s[j])
				{
				  npairs--;
				  l1=pairlist[npairs];
				  j1=index[l];
				  pairlist[j1]=l1;
				  index[l1]=j1;
				}
			      else
				{
				  pairlist[npairs]=l;
				  index[l]=npairs;
				  npairs++;
				}
			    }
			  en-=de+de; 
			  s[i]=-s[i];
			  mag+=s[i]+s[i];
			}
		    }
		  else //kawasaki;
		    {
		      flip=0;
		      rn=rn*mult+add;
		      l=(rn*fact1)*npairs;
		      l=pairlist[l];
		      i=l/3;
		      k=l%3;
		      j=i+dz[k];
		      if(j>=step2)j-=step2;//taking into account		
		      de=0;
		      for(k=0;k<6;k++)
			{
			  l=i+dz[k];
			  if(l>=step2)l-=step2;//taking into account
			  if(l<0)l+=step2;//periodic boundaries
			  if(l!=j)
			    de+=s[l];
			  l=j+dz[k];
			  if(l>=step2)l-=step2;//taking into account
			  if(l<0)l+=step2;//periodic boundaries
			  if(l!=i)
			    de-=s[l];
			}
		      if(s[i]<0)
			de=-de;
		      if(de<0)
			flip=1;
		      else
			{
			  rn=rn*mult+add; //apply Metropolis criyterion;
			  p=rn*fact1;
			  if(p<kaw[de]) 
			    flip=1;            //accept spin flip
			}
		      if(flip)
			{
			  for(k=0;k<3;k++)
			    {
			      l=i*3+k;
			      j2=i+dz[k];
			      if(j2>=step2)j2-=step2;//taking into account
			      if(j2==j)continue;			      
			      if(s[i]!=s[j2])
				{
				  npairs--;
				  l1=pairlist[npairs];
				  j1=index[l];
				  pairlist[j1]=l1;
				  index[l1]=j1;
				}
			      else
				{
				  pairlist[npairs]=l;
				  index[l]=npairs;
				  npairs++;
				}
			    }
			  for(k=0;k<3;k++)
			    {
			      j2=i-dz[k];
			      if(j2<0)j2+=step2;//taking into account
			      l=j2*3+k;
			      if(s[i]!=s[j2])
				{
				  npairs--;
				  l1=pairlist[npairs];
				  j1=index[l];
				  pairlist[j1]=l1;
				  index[l1]=j1;
				}
			      else
				{
				  pairlist[npairs]=l;
				  index[l]=npairs;
				  npairs++;
				}
			    }

			  for(k=0;k<3;k++)
			    {
			      l=j*3+k;
			      j2=j+dz[k];
			      if(j2>=step2)j2-=step2;//taking into account
			      
			      if(s[j]!=s[j2])
				{
				  npairs--;
				  l1=pairlist[npairs];
				  j1=index[l];
				  pairlist[j1]=l1;
				  index[l1]=j1;
				}
			      else
				{
				  pairlist[npairs]=l;
				  index[l]=npairs;
				  npairs++;
				}
			    }
			  for(k=0;k<3;k++)
			    {
			      j2=j-dz[k];
			      if(j2<0)j2+=step2;//taking into account
			      l=j2*3+k;
			      if(j2==i)continue;			      
			      if(s[j]!=s[j2])
				{
				  npairs--;
				  l1=pairlist[npairs];
				  j1=index[l];
				  pairlist[j1]=l1;
				  index[l1]=j1;
				}
			      else
				{
				  pairlist[npairs]=l;
				  index[l]=npairs;
				  npairs++;
				}
			    }
			  s[i]=-s[i];
			  s[j]=-s[j];
			  en-=de+de;			  
			}
		    }
		  U=-en-mag*H;
		  k=mag*bin1+nbin;
		  dist_m[k]++;
		  Ch+=U*U;
		  aven+=U;
		  amag+=mag;
		  avmag+=abs(mag);
		  time+=dt;
		  if(time >=ts)
		    {
		      printf("%d %lf\n",nc,ts);
		      sprintf(fname1,"%s-%03d-%03d-%03d.mag",fname,iT,iH,nc);
		      ff1=fopen(fname1,"w");
		      for(i=0;i<nbin2;i++)
			if(dist_m[i])
			  {
			    fprintf(ff1,"%lf %lf\n",(i-nbin+0.5)/nbin,dist_m[i]);
			    dist_m[i]=0;
			  }
		      fclose(ff1);
		      sprintf(fnamebmp,"%s-%03d-%03d-%03d.bmp",fname,iT,iH,nc);
		      savebmp(fnamebmp,s,step);
		      sprintf(fnamebmp,"%s-%03d-%03d-%03d.bit",fname,iT,iH,nc);
		      savebit(fnamebmp,s,step2);
		      if(nc<=lc)
			ts*=bst;
		      else
			ts+=tc;
		      nc++;
		    }
		}
	      aven*=fn;
	      Ch*=fn;
	      Ch-=aven*aven;
	      fprintf(ff,"%lf %lf %lf %lf %lf %lf %lf\n",time/n,aven/step2,Ch/step2,avmag*fn/step2,amag*fn/step2,T,H);
	      fflush(ff);
	      //	printf("%d %lf %lf %d %d %lf %lf\n",irun,aven,avmag*fn,mag,en,T,H);
	    }
	  // fprintf(ff,"&\n");
	}
    }
  /* 
     simple spin output
     for(k=0;k<step;k++)
     {
     for(i=0;i<step;i++)
     {
     printf("%1d",(s[j]+1)>>1);
     j++;
     }
     printf("\n");
     } */
      fclose(ff);
}



