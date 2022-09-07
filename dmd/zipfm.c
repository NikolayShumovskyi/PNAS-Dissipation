#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define alphabet "abcdefghijklmnopqrstuvwxyz"
#define Alphabet "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
#define MAX_SIZE 0x10000

struct node{
  unsigned char *word;
  long int  count;
  long int  order;
  long int  rank;
  struct node *left;
  struct node *right;
  struct node *more;
  struct node *less;
  struct node *equal;
};

struct node ** order, ** rank;
short alv [256];
FILE *infile, *outfile;
unsigned char infilename[30], outfilename[30];
unsigned char *t,*PEOF;
long int  word_number,word_order,word_rank;

void fillalphabet (unsigned char * c, unsigned char * C);
long int  getword( unsigned char ** p);
long int  wordcmp(unsigned char * p1, unsigned char * p2);
struct node * treeadd( struct node * p, unsigned char * w, long n);
struct node * treevine( struct node * p, struct node * q);
struct node * treevine1( struct node * p, struct node * q);
void  treeprint(struct node * p);
void wordprint( struct node  * p);
long wordlength( struct node  * p);
void wordprintf( struct node  * p,FILE *ff);
void rankprintf( struct node  * p,FILE *ff);
void  treesort(struct node * p);
void  treerank(struct node * p);
unsigned char * wordalloc(unsigned char * w, long n);

void fillalphabet (unsigned char * c, unsigned char * C)
{
  long int  i;
  unsigned char *d;
  for(i=0;i<256;i++)alv[i]=0;
  i=0;
  for(d=c;*d;d++)alv[*d]=++i; 
  i=0;
  for(d=C;*d;d++)alv[*d]=++i;
}

long int  getword( unsigned char ** p)
{ unsigned char *d;
  long int  i;
  d=*p;
  /*if (d>PEOF-15)
  {printf("yyy");}*/
  while(d<PEOF)
    {
      if(alv[*d])break;
      d++;
    }
  *p=d;
  if(d==PEOF)return 0;
  i=0;
  while(alv[*d])
    {i++;if(d==PEOF)return -i;d++;}
  if(*d=='-')
    {
      d++; 
      if(*d=='\n')
	{ d++;
          while(alv[*d]){if(d==PEOF)return 0;d++;}
          *p=d;
	  i=getword(p);
        }
    }
  return i;
}

long int  wordcmp(unsigned char * p1, unsigned char * p2)
{ 

  long int  i;
  i=0;

  while((alv[p1[i]]==alv[p2[i]])&(alv[p1[i]]>0))i++;
  if (alv[p1[i]]>alv[p2[i]])return 1;
  if (alv[p1[i]]<alv[p2[i]])return -1;
  return 0; 

}

unsigned char * wordalloc(unsigned char * w, long n)
	{ 
			long i;	
			unsigned char * c = (unsigned char *)malloc((n+1)*sizeof(unsigned char));
 			if (c==NULL){printf("no memomory %ld \n",word_number);return NULL;}
  	 c=memcpy(c,w,n);
  	 c[n]=0;
  	 return c;
  }

struct node * treeadd( struct node * p, unsigned char * w, long n)
{
  long int  cond;
  if(p==NULL)
    {
      p=(struct node * )malloc(sizeof(struct node));
     	if (p==NULL){printf("no memory %ld \n",word_number);return NULL;}
      else
      {
     p->word=wordalloc(w,n);
     if(p->word!=NULL)
	{  
	  p->count=1;
	  p->left=NULL;
	  p->right=NULL;
          p->more=NULL;
          p->less=NULL;
          p->equal=NULL;
          word_number++;
	}
    }
    }    
  else if ((cond=wordcmp(w,p->word))==0)p->count++;
  else if (cond<0)p->left=treeadd(p->left,w,n);
  else p->right=treeadd(p->right,w,n);
  return p;
  
}
struct node * treevine( struct node * p, struct node * q)
{ int cond;

  if (p==NULL)p=q;
  else if ((cond=q->count - p->count)>0) p->more=treevine(p->more,q);
  else if(cond<0)
    p->less=treevine(p->less,q);
    else {q->equal=p->equal; p->equal=q;}
  return p;
}


struct node * treevine1( struct node * p, struct node * q)
{ int cond;

  if (p==NULL)p=q;
     
  else {wordprint(p); if ((cond=q->count - p->count)>0) p->more=treevine1(p->more,q);
  else if(cond<0)
    p->less=treevine1(p->less,q);
    else { q->equal=p->equal; p->equal=q; }}
   
  return p;
}


void  treeprint(struct node * p)
{ 
  unsigned char *c;
  if(p!=NULL)
    {
      treeprint(p->left);
      printf("%6ld ",p->count);
      c=p->word;
      while(alv[*c]){printf("%c",*c);c++;}
      printf("\n");
      treeprint(p->right);
    }
}
void wordprint( struct node  * p)
{ 
  unsigned char *c;
  unsigned char d;
  c=p->word;
  if (p==NULL){printf("No word"); return;}
  printf("%6ld %6ld %6ld ",p->rank,p->count,p->order);
 /* printf("%lx %lx %lx ",p->less,p->more,p->word);*/
  while(alv[*c]){printf("%c",*c);c++;}
  printf("\n");
 /* scanf("%c",&d);*/
}

long wordlength( struct node  * p)
{ 
  unsigned char *c;
  long j;
  c=p->word;
  j=1;
  while(alv[*c]){j++;c++;}
return j;
}
void wordprintf( struct node  * p,FILE *ff)
{ 
  unsigned char *c;
  const unsigned char d=9;
  fprintf(ff,"%6ld%c%6ld%c%6ld%c",p->rank,d,p->count,d,p->order,d);
  c=p->word;
  while(alv[*c]){fprintf(ff,"%c",*c);c++;}
  fprintf(ff,"\n");
}
void rankprintf( struct node  * p,FILE *ff)
{ 
  const unsigned char c=9;
  fprintf(ff,"%6ld%c%6ld\n",p->rank+1,c,p->count);
}
void  treesort(struct node * p)
{ 
  unsigned char *c;
  if(p!=NULL)
    {
      treesort(p->left);
      order[word_order]=p;
      p->order=word_order;
      word_order++;
      treesort(p->right);
    }
}

void  treerank(struct node * p)
{ 
  unsigned char *c;
  if(p!=NULL)
    { struct node *q; 
      treerank(p->more);
      for(q=p;q;q=q->equal)
      {
      rank[word_rank]=q;
      q->rank=word_rank;
      word_rank++;
      }
      treerank(p->less);
    }
}


main()
{
  long int i,n,lt,cur_size,j;
  unsigned char * p, * pold;
  struct node * root, * root1;
  root= NULL;
  fillalphabet(alphabet,Alphabet);
  t= (unsigned char * )malloc(MAX_SIZE * sizeof(unsigned char));
  printf("%d\n",(int)(MAX_SIZE));
  printf("What is in file name\n");
  scanf("%s",&infilename);
 
  infile=fopen(&infilename,"r");

  n=0;
  i=0; 
  word_number=0;
  
  do
  	{
  	cur_size=MAX_SIZE-n;
  	p=t+n;
  	lt=fread(p,sizeof(unsigned char),cur_size,infile);
    PEOF=p+lt-1;
  	p=t;
  	while((n=getword(&p))>0)
    	{
      	root=treeadd(root,p,n);
      	if(root==NULL){fclose(infile);return;}
      	p+=n;
   		}
  	if (n<0)
    	{
    		n=-n;
    		p=PEOF-n+1;
     		t=memcpy(t,p,n);
    	}
  	i+=cur_size;
  	printf("char. read=%ld, words found=%ld\n",i,word_number);     
   	}
  	while(lt==cur_size); 
  fclose(infile); 
  printf("Done! Words found=%ld\n",word_number);      
  order=(struct node **)malloc(word_number*sizeof(struct node *));
  if (order==NULL){printf("no memory\n");return;}
  
  rank=(struct node **)malloc(word_number*sizeof(struct node *));
  if (rank==NULL){printf("no memory\n");return;}
  printf("Sorting alphabetically...\n");
  word_order=0;
  treesort(root);
  printf("Alphabeticall sort done\n");
  root1=NULL;
  j=0;
  printf("Now sorting according to word frequency...\n");
  for(i=0;i<word_number;i++)
    {
      if(order[i]==NULL){printf("error\n");return;}
      root1=treevine(root1,order[i]);
      j+=wordlength(order[i])*order[i]->count;
    }
  printf("Total number of letters=%ld \n",j);
  word_rank=0;
  treerank(root1);
  printf("Sorting done\n\n");
  printf("What is out file name for Zipf analysis?\n");
  scanf("%s",&outfilename);
  outfile=fopen(&outfilename,"w");
   { 
      long int old_rank=0;
      for(i=0;i<word_number;i++)
      	if (rank[i]->count != old_rank)
      		{
     			 rankprintf(rank[i],outfile);
      		 old_rank=rank[i]->count;
      		}
   }
  fclose(outfile);
  { 
    unsigned char c;
    c=' ';
    while((c!='y')&& (c!='n'))
    {
    printf("List alphabetically ? y/n\n");
    scanf("%c",&c);
    }
    if (c=='y')
    {
      printf("What is out file name\n");
      scanf("%s",&outfilename);
      outfile=fopen(&outfilename,"w");
      fprintf(outfile,"Rank      Count     Order\n"); 
      for(i=0;i<word_number;i++)
      wordprintf(order[i],outfile); 
      fclose(outfile);
    }
    c=' ';
    while((c!='y')&& (c!='n'))
    {
    printf("List according to frequency? y/n\n");
    scanf("%c",&c);
    }
    if (c=='y')
    {
      printf("What is out file name\n");
      scanf("%s",&outfilename);
      outfile=fopen(&outfilename,"w");
      fprintf(outfile,"Rank      Count     Order\n"); 
      for(i=0;i<word_number;i++)
      wordprintf(rank[i],outfile);
      fclose(outfile); 
    }
   } 
}










