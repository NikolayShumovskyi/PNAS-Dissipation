/*Rules : 
1. Every cell has a list of atoms (clist). 

2. Every atom has a neighbor list of atoms. These are atoms that
were located at time t0 within the spere of radius maxrb +2eps.
the eps is passed from make_key_system into search as argument of init_tables.
The coordinates of each atom is saved in a list of atoms  
search.atoms. It also has the time of exit of the circle of radius
eps. the collision type with this circle has type n-1+4. 
Once this collision happens, the neiborlist is updated. when the
neibor list is updataed all the atoms from the list with external collisions
with this atoms will be removed from the collision tables, because they are added to
the collision tables only if their collsion times are smaller than the time of the 
update of the list. 
when the atom is moving from one cell to another no event updates are taken place except
that the the event ofthis collision is removed from the queue.


 
Every atom is represented in the collision list of a cell . 


2. Each atom has a collision list (alist) which has 8 elements: 
first atom p, next atom q>p, collision type c, time t, nextp, nextq,
prevp, prevq. nextp points to the next element involved with p;
the same for q and the nextq points to the next element involved with q. 
Each atom list starts with the wall collision (q>=n_atom). Scanning
of the list for say atom m is going in the following way :
if m==p we use nextp or prev p, otherwise we use nextq or prevq 

3. If the collision is of internal type, it is always listed in the
collision list.

4. If the collision is of external type it is listed if and only if
its time is smaller than both of the times of leaving the cell for each
atom.

6. matrix of collisions shows the types of external colisions for atom
pairs which are not represented in collision table. External collisions
can be of finite energy. Historically, for private use, if their energy 
is positive , it means attraction (gain of kinetic energy) For public
use we use conventional negative potential energies for attraction.
The numbers of external collisions are numbered from 0 to nen1.
From 0 to nen-1 - collisions with finite energies;
from nen to nen0-1 collisions with infinite energy: ellastic repulsions)
from nen0 to nen1-1 collisions of bonded pairs. 
It is a triangular array of int ** excoll

7. If the collision is of the wall type (q>=n_atom), 
then find_atoms from new_loc 
is called. This function is removing atom from one cell and adding
it to another. From this function find_collisions is called
in which we go 
through the collision list of the atom, and do not clean up the 
lists of the neighboring atoms which must be in
the intersection of the neighboring cells of the old and the new
cell, constructing the array collp[], which originally filled with -1,
indicating that an atom is not present in the cell, and fill this
array with -2 which means that the atoms will collide as before.
Then in function find_atoms we go though the cell list of 27 
neighboring cells of the new cell and add  collp[] if it is equal to -1.  
Such element of collp collq is filled with excoll[][] of atoms type.

8. In squeeze1, in function kill_collisions 
we scan the collision lists of atom lists of p and q
and completly clean them up cleaning also the records in the lists of 
the atoms involved in these collisions and construct an arrays 
collp[] and collq[] which
originally filled with -1 indicating that an atom is not present in
the cells. While scanning we delete the collisions in which given
atoms (p or q) take part and put the values of collisions in the collp
or collq memorizing the values of atoms into arrays atomp and atomq
and amounts nq and np of colliding atoms. After this in function 
add_extern_atom we go through
the cell lists of 27 neighboring cells and add those atoms
into collp and collq if the corresponding element in collp and collq is -1,
and such element of collp collq is filled with excoll[][] of atoms
type. These values are ovewritten by a type of particle collision if
it is found later.  
*/

/* cell list of atoms present in the cell togther with the data on
   wall collisions */
typedef struct cell_list
{ 
  int n;
  struct cell_list * next;
} clist;
typedef struct neib_list
{ 
  int n;
  int c;
  struct neib_list * next;
} nlist;

/* atom list of collisions */
typedef struct atom_list{
  double t;
  int p;
  int q;
  struct atom_list * pnext;
  struct atom_list * qnext;
  struct atom_list * pprev;
  struct atom_list * qprev;
  struct atom_list * left;
  struct atom_list * right;
  struct atom_list * parent;

} alist;  

typedef struct fixed_atom{
  double x;
  double y;
  double z;
  double t;
}fatom;


/**********************************************************************
 * The structure for effective search of minimal collision   
 * each cell has its  cellist which is the cell list
 * of atoms (clist). Each atom has two lists of collisions (alist).
 * atom_collisionp and atom_cllisionq; In the atom_collisionp
 * the collisions in which the second partner has number larger
 * than this atom. This list is sorted. Other list atom_cllisionq
 * has collisions in which the second partner has smaller number
 * than this atom. This list will be kept unsorted.
 * The nodes of alist are kept in "storage". To indicate which of 
 * nodes are free index array "begin" and a current free address  
 * "free" is used. In update_table, we make bubble sort for each
 * atom involved in atomq and atomp. If the topmost record, corresponding
 * to the next collision of this atom is changed we determine
 * if this atom is the smallest partner construct the array of
 * atoms changed_atom, analgous to atomp, which will keep the the atom numbers
 * which soonest collision has been changed. Next we go through this
 * list and find cells wich contain atom with p-participant of these 
 * collisions. Each cell has a pointer to the soonest collision
 * record cell_champion.
 * We define the affected cells in function
 * cell_change. 
 * The champions of each cell is determined by the soonest collision
 * in the collision list of cell_atom[i] in function local_champion
 * (i<maxadd,    
 * the number of cells in the system) form an olympic table "olymp" 
 * containing the cell numbers who win in the levels from zero to   
 * final. If the cell is empty, it has cell_atoms[i]=NULL. Local          
 * championships  are held by bubble sort.                          
 **********************************************************************/  


typedef struct
{
  clist * cell_storage;
  clist ** cell_atoms;
  nlist * neib_storage;
  nlist ** neib;
  nlist **neib_begin;
  nlist **neib_free;
  alist **begin;
  alist **free;
  alist **atom_collisionp;
  alist **atom_collisionq;
  alist * storage;
  alist * root;
  alist * left;
  fatom * fixed_atoms;

  int maxfree;
  int x;
  int y;
  int z;
  int * atomp;
  int * atomq;
  int * collp;
  int * collq;
  int np;
  int nq;
  int n;
  int nhalf;
  int maxalloc;
  double eps2;
  double rneib;
  //double b_d;
  //double b_n;
  //int t_l;
  //double l_d;
  //double l_n;
} tsearch;

extern void change_cell(int p1, size_al address1);
extern void change_neib(int p1); 
extern int allocsearch(int n);
extern void initsearch(void);
extern int init_tables(void);
extern void set_maxfree(int a);
extern int get_maxfree(void);
extern double get_eps2(void);
extern int get_free(void);
extern int squeeze1(int p1, int * collp, int * atomp);
extern void update_table(int p1,int q1, int ct1);
extern int squeeze_table(int * p1, int * q1, double * timea);
extern int check_atoms(int add);
extern int bond(int i, int j);
extern int find_neighbors(int p, int * neib);
extern fatom * get_fixed_atom(int i);
  
extern int * get_collp(void);
extern int * get_collq(void);
extern int * get_atomp(void);
extern int * get_atomq(void);
extern int get_np(int i);
extern int get_nq(int i);
extern int pairs(int (*)(int,int,int));

extern void add_atom2cell(int i);
int find_collisions(int p1);
void print_collisions(int p1,int q1);
extern int realloc_search (void);
extern int find_shell(int p1, int ct, int * neib);
void set_rneib(double eps1, double maxrb1);
double get_eps2(void);
int kill_neib_list(int i);
int kill_neib(int i, int j);
int kill_collisions(int p1);
int get_coll_type(int i,int j);
int set_coll_type(int i,int j,int ct);
extern void clean_neib(void);
extern void update_neib(int i1,int i2);
int bubble (double t, int p1, int q1);
