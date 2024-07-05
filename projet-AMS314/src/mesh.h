/* 


  Mesh Structures that will be used in the project 
  
  
  Provided functions : 
       
  

*/
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <libmesh6.h>
#include <lplib3.h>
#include <math.h>
#include <time.h>



typedef double double3[3];
typedef double double2[2];
typedef int    int5[5];
typedef int    int4[4]; 
typedef int    int3[3];
typedef int    int1[1];
typedef unsigned long long int u64;


/* 
  
  A simple mesh structure to hold the geomety/connectivity 
  Note that 2D mesh are stored as 3D Mesh, z-coordinate will be set to zero

  Provided function msh_init, msh_read, msh_write, msh_check
  The mesh files are based on the meshb format (https://github.com/LoicMarechal/libMeshb), can be viewed with ViZiR4, see README to install it. 

  The following functions has to be implemented  : 
    msh_bondingbox  : compute the bounding box of mesh
    msh_reorder     : reorder an input mesh according to Z-ordering 
    msh_reorder     : simple smoothing algorithm of volume points   
    msh_neighbors   : build the set of Tris surrounding elements, efficient implementation with Hash Tab
    msh_neighborsQ2 : a quadratic setting of the neigbhoring structures 

 

*/

typedef struct mesh_vertex
{  
  double Crd[2];
  
  long long int icrit; /* sorting creteria, to be used with qsort  */
  int idxNew;          /* new if after sort */
  int idxOld;          /* initial id  */ 
  
} Vertex; 


typedef struct mesh_edge
{
  int Ver[2]; /* index of the two vertices  */
  int Voi[2]; /* neigbhoring structures */
  int Ref;
  
  long long int icrit; /* sorting creteria, to be used with qsort  */
  
} Edge;


typedef struct mesh_triangle
{
  int Ver[3]; /* index of the three vertices  */
  int Voi[3]; /* neigbhoring structures */
  int Ref;
  
  long long int icrit; /* sorting creteria, to be used with qsort  */
  
} Triangle;


typedef struct t_mesh
{
  int Dim;
  int NbrVer, NbrTri, NbrEfr, NbrEdg;
  
  Vertex      *Ver;   /* list of vertices */
  Triangle    *Tri;   /* list of triangles */
  Edge        *Efr;   /* list of boundary edges */
  Edge        *Edg;   /* list of all the mesh egdes */
  
  /* bounding box */
  double bb[4];
  
} Mesh; 



/* Provided functions */
Mesh * msh_init();
Mesh * msh_read(char *file, int readEfr);
int    msh_write(Mesh *msh, char *file); 

double * sol_read(char *file, int mshDim, int mshNbrSol);



/* functions to be implemented  */
double aire(Mesh* msh, int iTri);
double    Q1(Mesh* msh, int iTri);
double    Q2(Mesh* msh, int iTri);
int    msh_boundingbox(Mesh *msh);         /* compute the bouding box of the mesh                         */
int    msh_reorder(Mesh *msh);             /* perform a mesh using morton curve/octree-based              */
int    msh_smooth(Mesh *msh, int nbrStep); /* a simple mesh smoohting algorithm                           */
int    msh_neighbors(Mesh *msh);           /* build TriVois with a hash table                             */
int    msh_neighborsQ2(Mesh *msh);         /* build TriVois with the naive quadratic approach             */
int    msh_neighbors_efr(Mesh *msh);  

int    nb_boundary_edge(Mesh *msh);
/* a provided simple hash table data structure */

typedef struct mesh_hash_table
{
  int  SizHead;   /* Maxmimum entries, the key is in [0,SizHead-1]*/
  int  NbrObj;    /* Number of object in the hash tables */
  int  NbrMaxObj; /* Maximum of object that can be store in the hash tab */ 
  int  *Head ;    /* Head[key%(SizHead)] = link to the first object having this key  in the LstObj list */
  int5 *LstObj;   /* List of objects in the Hash Tab */
  int NbrBoundaryEdge;
  int NbrEdge;
  
  /* LstObj[id][0:1] = ip1-ip2,     the 2 points defining the edge  */
  /* LstObj[id][2:3] = iTri1,iTri2, the two neighboring triangles having ip1-ip2 as points */
  /* LstObj[id][4]   = idnxt,       the link to the next element in collision, if = 0 last element of the list */
  
} HashTable;

/* Implementing the following function should be necessary */
HashTable * hash_init(int SizHead, int NbrMaxObj);          /* alloc and set htable ==> allocate Head, LstObj */

int hash_find(HashTable *hsh, int ip1, int ip2);            /* return the id found (in LstObj ), if 0 the object is not in the list */
int hash_add (HashTable *hsh, int ip1, int ip2, int iTri);  /* ==> add this entry in the hash tab */

HashTable* hashTable_Voi(Mesh *msh);

typedef struct mesh_doubly_linked_list DL; 
struct mesh_doubly_linked_list
{
  /* data */
  DL * previous_node;
  DL * next_node;
  int value;
} ;

DL * dl_init(int value);
int dl_append(DL* dl,int value);

typedef struct mesh_stack
{
  DL * top;
  int NbrObj;
} Stack;

Stack* stack_init();
int stack_push(Stack* sk, int value);
int stack_pop(Stack* sk);

typedef struct edge_hash_table
{
  /* data */
  int SizeHead;
  int NbrObj;
  int NbrMaxObj;
  int *Head;
  int3 *LstObj;
} EdgeHashTable;

EdgeHashTable* edge_hash_init(int SizeHead, int NbrMaxObj);
int edge_hash_find(EdgeHashTable *hsh, int ip1, int ip2);
int edge_hash_add(EdgeHashTable *hsh, int ip1, int ip2);

EdgeHashTable* hashTable_Efr(Mesh *msh);
/* Fonction used for adaptation */


typedef struct t_vector
{
  DL* head;
  DL* tail;
  int length;
} Vector;

Vector * vector_init();
int vector_append(Vector* v, int val);
DL* vector_elementAt(Vector* v, int pos);
int vector_inser(Vector* v, int pos, int val);
int vector_remove(Vector* v, int pos);
int vector_valueAt(Vector* v, int pos);
/* Writes a 2d metric tensor field in a file (extension should be .sol)*/
int msh_write2dmetric(char *file, int nmetric, double3 *metric);  
int msh_write2dfield_Triangles(char *file, int nfield, double *field);
int msh_write2dfield_Vertices(char *file, int nfield, double *field);


double * colorisation(Mesh *msh);

double2* points_unif(double xmin,double ymin,double xmax,double ymax,int Nx, int Ny);
//int localiser(Mesh *msh,int K, double2 P);

typedef struct t_circle
{
  double r;
  double* centre;
} Circle;

Circle* circle_init();

typedef struct t_set
{
  int size;
  DL* elem;
}Set;

Set* set_init();
int set_add(Set* s, int value);
int set_pop(Set* s);
int set_remove(Set* s, int value);
int set_contains(Set* s, int value);

typedef struct t__dinamic_mesh
{
  int Dim;
  int NbrVerMax, NbrTriMax;
  Mesh* msh;
}DMesh;

DMesh* dMesh_init(Mesh* msh, int lambda);
int dMesh_resize(DMesh* dMsh, int NbrVerMax, int NbrTriMax);


void test(Mesh* msh);

int Triangulation(double2* Ps, int NbrP, DMesh* dmsh);

Mesh* toMesh(DMesh* dmsh);

int* get_corner_index(Mesh* msh);

typedef struct t_photo
{
  Mesh* msh;
  int size;
  int* corners;
  double* data;
}Photo;


Photo* photo_init();

Photo * photo_read(char* file_sol, char *file_msh, int readEfr);

int save_photo(Photo* pht, char* file_name);

//double2* photo_gradient(Photo* pht); 

typedef struct t_compressed_images
{
  int size;
  double xmin,xmax,ymin,ymax;
  double2* pix;
  double* data;
}cPhoto;

cPhoto* cPhoto_init();

cPhoto* compress_1(Photo* pht, double level);

Photo* decode_cPhoto(cPhoto* cpt);