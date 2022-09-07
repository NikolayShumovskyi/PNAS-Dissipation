#include <stdio.h> /* standart input,output */
#include <ctype.h> /* for charachter recognition */
#include <stdlib.h> /* for conversion from char to dec */
#include <strings.h>
#include <time.h>
#include <math.h>
#include "bcp.h"
#include "controls.h"
#include "schedule.h"
#include "mutate.h"
#include "nucleation.h"
#include "q6.h"
enum schedule_keys 
{MIN_KEY=0,TIME_KEY=MIN_KEY,LX_KEY,LY_KEY,LZ_KEY,
N_COLL_T_KEY,N_COLL_V_KEY,HEAT_X_C_KEY,T_0_KEY,T_NEW_KEY,
 DV_DP_KEY,DX_DP_KEY,DY_DP_KEY,DZ_DP_KEY,P_0_KEY,P_X_KEY,P_Y_KEY,P_Z_KEY,MUTATE_KEY,NUC_KEY,Q6_KEY,ENERGY_KEY,END_KEY,MAX_KEY};

static int save_pos=0;
static double next_time=0;
static double init_time=0;
static int good=0;
static char ** keywords;
static char sched_name[100];
static double maxtime;
static double param[4];
static int n_mes=0;
static double L[3];
static int n_coll_T;
static int n_coll_V=10;
static double heat_x_c=0;
static double T0;
static double Tnew;
static double dv_dp=0;
static double dx_dp=0;
static double dy_dp=0;
static double dz_dp=0;
static double P0=0;
static double PX=0;
static double PY=0;
static double PZ=0;
  char keyword[60];
extern int schedule_ok(void)
{
  if(good)
  {
    if(get_time()>=next_time)
      return 1;
  }
  return 0;
}

void init_schedule_param(void)
{
  param[0]=0;
  param[1]=0;
  param[2]=0;
  param[3]=0;
  n_mes=0;
  return;
}



int init_schedule_keywords(void)
{ 
  int i;   
  keywords=(char **)malloc(MAX_KEY*sizeof(unsigned char *));
  if(!keywords) return 0;
  keywords[TIME_KEY]="TIME";
  keywords[LX_KEY]="LX";
  keywords[LY_KEY]="LY";
  keywords[LZ_KEY]="LZ";
  keywords[N_COLL_T_KEY]="N_COLL_T";
  keywords[N_COLL_V_KEY]="N_COLL_V";
  keywords[HEAT_X_C_KEY]="HEAT_X_C";
  keywords[T_0_KEY]="T_0";
  keywords[T_NEW_KEY]="T_NEW";
  keywords[P_0_KEY]="P_0";
  keywords[P_X_KEY]="P_X";
  keywords[P_Y_KEY]="P_Y";
  keywords[P_Z_KEY]="P_Z";
  keywords[DV_DP_KEY]="DV_DP";
  keywords[DX_DP_KEY]="DX_DP";
  keywords[DY_DP_KEY]="DY_DP";
  keywords[DZ_DP_KEY]="DZ_DP";
  keywords[MUTATE_KEY]="MUTATE";
  keywords[NUC_KEY]="NUC_NAME";
  keywords[Q6_KEY]="Q6_NAME";
  keywords[ENERGY_KEY]="ENERGY";
  keywords[END_KEY]="END";
  return 1;
}

extern int schedule_init(void)
{
  FILE * infile;
  FILE * save=NULL;
  char value[60];
  double dummy;   
  int if_energy=0;

  int i,nt;
  good=0;
  nt=0;
  printf("we do not use schedule\n");
  if(yes())return good;
  n_coll_T= get_atom_number();
  for(i=0;i<3;i++)
    L[i]=0;
    Tnew=get_temp();
    T0=get_temp();

  printf("What is schedule file name?\n");
  scanf("%s",sched_name);
  infile=fopen(sched_name,"r");
  if(!infile)return good;
  
  if(!init_schedule_keywords())return good;
  init_schedule_param();
  while(!feof(infile))
    {
      enum schedule_keys key;
      fscanf(infile,"%s%s",keyword,value);
      for(key=MIN_KEY;key<MAX_KEY;key++)
	if(!strcmp(keyword,keywords[key]))
	  {
	    switch(key){
	    case TIME_KEY:{
	      dummy=atof(value);
              if(dummy<0)return good; 
	      if(!nt)
		init_time=dummy;
	      else 
		next_time=dummy;
	      if(nt==1)goto finish;
	      nt++;
	      break;
	    }
	    case LX_KEY:{nt=1;dummy=atof(value);if(dummy>0)L[0]=dummy;break;}
	    case LY_KEY:{nt=1;dummy=atof(value);if(dummy>0)L[1]=dummy;break;}
	    case LZ_KEY:{nt=1;dummy=atof(value);if(dummy>0)L[2]=dummy;break;}
            case N_COLL_T_KEY:
	      {nt=1;dummy=atof(value);if(dummy>1)n_coll_T=dummy;break;}
            case N_COLL_V_KEY:
	      {nt=1;dummy=atof(value);if(dummy>1)n_coll_V=dummy;break;}
	    case HEAT_X_C_KEY:
	      {nt=1;dummy=atof(value);if(dummy>=0)heat_x_c=dummy;break;}
	    case T_0_KEY:{nt=1;dummy=atof(value);if(dummy>0)T0=dummy;break;}
	    case P_0_KEY:{nt=1;dummy=atof(value);P0=dummy;break;}
	    case P_X_KEY:{nt=1;dummy=atof(value);PX=dummy;break;}
	    case P_Y_KEY:{nt=1;dummy=atof(value);PY=dummy;break;}
	    case P_Z_KEY:{nt=1;dummy=atof(value);PZ=dummy;break;}
	    case T_NEW_KEY:
	      {nt=1;dummy=atof(value);if(dummy>0)Tnew=dummy;break;}
	    case ENERGY_KEY:
	      {nt=1;dummy=atof(value);if_energy=set_energy(dummy);break;}
	    case DV_DP_KEY:
	      {nt=1;dummy=atof(value);if(dummy>=0)dv_dp=dummy;break;}
	    case DX_DP_KEY:
	      {nt=1;dummy=atof(value);if(dummy>=0)dx_dp=dummy;break;}
	    case DY_DP_KEY:
	      {nt=1;dummy=atof(value);if(dummy>=0)dy_dp=dummy;break;}
	    case DZ_DP_KEY:
	      {nt=1;dummy=atof(value);if(dummy>=0)dz_dp=dummy;break;}
	    case MUTATE_KEY:
	      {nt=1;mutate_init(value);break;}
	    case NUC_KEY:
	      {nt=1;open_nucleation(value);break;}
	    case Q6_KEY:
	      {nt=1;open_q6(value);break;}
            case END_KEY: goto finish;   
	    }
	    break;
	  }
    }
 finish:
  if (next_time<0)return good;
  maxtime=next_time;
  while(!feof(infile))
      {
	fscanf(infile,"%s %s",keyword, value);
	if(feof(infile))break;
	if(!strcmp(keyword,keywords[END_KEY]))break;
	if(!strcmp(keyword,keywords[TIME_KEY]))
	  {
	    dummy=atof(value);
	    if(dummy>maxtime)
	      maxtime=dummy;
            else
	      {
		printf("Shcedule is not correct: TIME=%s\n",value);
		return good;
	      }
	  }   
      }
  fclose(infile);
  good=1;
  
  set_time(init_time);
  set_max_time(maxtime);
  set_delta_ll(n_coll_T);
  if(!if_energy)set_temp(Tnew);
  set_coeff(heat_x_c);
  set_temp_limit(T0);
  set_max_gap_mes(n_coll_V);
  set_L_limit(L);
  set_press_limit(P0,PX,PY,PZ);
  set_press_coeff(dv_dp,dx_dp,dy_dp,dz_dp); 
  return good; 
}

extern int read_schedule(void)
{
  char value[60];
  int nt=0;
  double new_next_time,newmaxtime,dummy;
  int found=0;
  FILE * infile=fopen(sched_name,"r");
  if(infile)
    {
  newmaxtime=next_time;
  while(!feof(infile))
      {
	fscanf(infile,"%s %s",keyword, value);
	if(!strcmp(keyword,keywords[END_KEY]))break;
	if(!strcmp(keyword,keywords[TIME_KEY]))
	  {
	    dummy=atof(value);
	    if(dummy>=next_time)
	      {found=1;break;}
	  }   
      }
  if(!found)
    {
      good=0; 
      fclose(infile);
      maxtime=0; 
      set_max_time(maxtime);
    }
  if(dummy>next_time)
    {
      fclose(infile);
      if(dummy>get_max_time())
	{maxtime=dummy;
	 set_max_time(maxtime);}
      return dummy;
    }
  while(!feof(infile))
    {
      enum schedule_keys key;
      fscanf(infile,"%s%s",keyword,value);
      for(key=MIN_KEY;key<MAX_KEY;key++)
	if(!strcmp(keyword,keywords[key]))
	  {
	    switch(key){
	    case TIME_KEY:{
	      dummy=atof(value);
	      new_next_time=dummy;
	      goto finish;
	    }
	    case LX_KEY:{dummy=atof(value);if(dummy>0){nt=1;L[0]=dummy;}break;}
	    case LY_KEY:{dummy=atof(value);if(dummy>0){nt=1;L[1]=dummy;}break;}
	    case LZ_KEY:{dummy=atof(value);if(dummy>0){nt=1;L[2]=dummy;}break;}
            case N_COLL_T_KEY:
	      {
		dummy=atof(value);
		if(dummy>1)
		 {
		   n_coll_T=dummy;
                   set_delta_ll(n_coll_T);
		 }   
		break;
	      }
            case N_COLL_V_KEY:
	      {
		dummy=atof(value);
		if(dummy>1)
		 {
		   n_coll_V=dummy;
                   set_max_gap_mes(n_coll_V);
		 }   
		break;
	      }
	    case HEAT_X_C_KEY:
	      {
		dummy=atof(value);
		if(dummy>=0)
		  {
		    heat_x_c=dummy;
		    set_coeff(heat_x_c);
		  }
		    break;
	      }
	    case T_0_KEY:
	      {
		dummy=atof(value);
		if(dummy>0)
		  {
		    T0=dummy;
		    set_temp_limit(T0);
		  }
		break;
	      }
	    case P_0_KEY:
	      {
		dummy=atof(value);P0=dummy;
		set_press_limit(P0,PX,PY,PZ);
		break;
	      }
	    case P_X_KEY:
	      {
		dummy=atof(value);PX=dummy;
		set_press_limit(P0,PX,PY,PZ);
		break;
	      }
	    case P_Y_KEY:
	      {
		dummy=atof(value);PY=dummy;
		set_press_limit(P0,PX,PY,PZ);
		break;
	      }
	    case P_Z_KEY:
	      {
		dummy=atof(value);PZ=dummy;
		set_press_limit(P0,PX,PY,PZ);
		break;
	      }
	    case T_NEW_KEY:
	      {
		dummy=atof(value);
		if(dummy>0)
		  {
		    Tnew=dummy;
		    set_temp(Tnew);
		  }
		break;
	      }
	    case ENERGY_KEY:
	      {
		dummy=atof(value);
		set_energy(dummy);		  
		break;
	      }
	    case DV_DP_KEY:
	      {
		dummy=atof(value);
		if(dummy>=0)
		  {
		    dv_dp=dummy;
		    set_press_coeff(dv_dp,dx_dp,dy_dp,dz_dp);
		  }
		break;
	      }
	    case DX_DP_KEY:
	      {
		dummy=atof(value);
		if(dummy>=0)
		  {
		    dx_dp=dummy;
		    set_press_coeff(dv_dp,dx_dp,dy_dp,dz_dp);
		  }
		break;
	      }
	    case DY_DP_KEY:
	      {
		dummy=atof(value);
		if(dummy>=0)
		  {
		    dy_dp=dummy;
		    set_press_coeff(dv_dp,dx_dp,dy_dp,dz_dp);
		  }
		break;
	      }
	    case DZ_DP_KEY:
	      {
		dummy=atof(value);
		if(dummy>=0)
		  {
		    dz_dp=dummy;
		    set_press_coeff(dv_dp,dx_dp,dy_dp,dz_dp);
		  }
		break;
	      }
	    case MUTATE_KEY:
	      {nt=1;mutate_init(value);break;}
	    case NUC_KEY:
	      {nt=1;open_nucleation(value);break;}
	    case Q6_KEY:
	      {nt=1;open_q6(value);break;}
            case END_KEY: goto finish;   
	    }
	    break;
	  }
    }
 finish:
  if(new_next_time<=next_time)
    {
      good=0; 
      fclose(infile);
      maxtime=0; 
      set_max_time(maxtime);
      next_time=0;
      return good;
    }
  if(new_next_time>get_max_time())
    {
      maxtime=new_next_time;
      set_max_time(maxtime);
    }
     next_time=new_next_time;
  fclose(infile);
  if(nt)  
    set_L_limit(L);
    }
  return good;
}














