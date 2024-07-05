#include <mesh.h>
#define EMPTYOBJECT 0
#define PI 3.1415926535897932384626433832795

int tri2edg[3][2] = {{1,2},{2,0},{0,1}};


Mesh * msh_init()
{
  Mesh *msh = malloc(sizeof(Mesh));
  if ( ! msh ) return NULL;
  
  msh->Dim    = 0;
  msh->NbrVer = 0;
  msh->NbrTri = 0;
  msh->NbrEfr = 0;
  msh->NbrEdg = 0;
  
  msh->Ver = NULL;
  msh->Tri = NULL;  
  msh->Efr = NULL;  
  msh->Edg = NULL;  
  
  msh->bb[0] = 0.0; /* xmin  */
  msh->bb[1] = 0.0; /* xmax  */
  msh->bb[2] = 0.0; /* ymin  */
  msh->bb[3] = 0.0; /* ymax  */
  
  return msh;
}  
  


Mesh * msh_read(char *file, int readEfr)
{
  char   InpFil[1024];
  float  bufFlt[2];
  double bufDbl[2];
  int    i, bufTri[4], bufEfr[3];
  int    FilVer, ref; 
  
  int fmsh = 0;
  
  if ( ! file ) return NULL;
  
  Mesh * msh = msh_init();
    
  //--- set file name 
  strcpy(InpFil,file);
  if ( strstr(InpFil,".mesh") ) {
    if ( !(fmsh = GmfOpenMesh(InpFil, GmfRead, &FilVer, &msh->Dim)) ) {
      return NULL;
    }    
  }
  else {
    strcat(InpFil,".meshb");
    if ( !(fmsh = GmfOpenMesh(InpFil, GmfRead, &FilVer, &msh->Dim)) ) {
      strcpy(InpFil,file);
      strcat(InpFil,".mesh");
      if ( !(fmsh = GmfOpenMesh(InpFil, GmfRead, &FilVer, &msh->Dim)) ) {
        return NULL;
      }    
    } 
  }
  
  printf(" File %s opened Dimension %d Version %d \n", InpFil, msh->Dim, FilVer);
  
  msh->NbrVer = GmfStatKwd(fmsh, GmfVertices);
  msh->NbrTri = GmfStatKwd(fmsh, GmfTriangles);
  
  //--- allocate arrays 
  msh->Ver = calloc( (msh->NbrVer+1), sizeof(Vertex)   );
  msh->Tri = calloc( (msh->NbrTri+1), sizeof(Triangle) );  
  
  
  //--- read vertices   
  GmfGotoKwd(fmsh, GmfVertices);
  if ( msh->Dim == 2 ) {
    if ( FilVer == GmfFloat ) {		// read 32 bits float
      for (i=1; i<=msh->NbrVer; ++i) {
        GmfGetLin(fmsh, GmfVertices, &bufFlt[0], &bufFlt[1], &ref);
        msh->Ver[i].Crd[0] = (double)bufFlt[0];
        msh->Ver[i].Crd[1] = (double)bufFlt[1];
      }
    }
    else  {	// read 64 bits float
      for (i=1; i<=msh->NbrVer; ++i) {
        GmfGetLin(fmsh, GmfVertices, &bufDbl[0], &bufDbl[1], &ref);
        msh->Ver[i].Crd[0] = bufDbl[0];
        msh->Ver[i].Crd[1] = bufDbl[1];
      }  
    }
  }
  else {
    fprintf(stderr,"  ## ERROR: 3D is not implemented\n");
    exit(1);
  }
  
  
  //--- read triangles   
  GmfGotoKwd(fmsh, GmfTriangles);
  for (i=1; i<=msh->NbrTri; ++i) {
    GmfGetLin(fmsh, GmfTriangles, &bufTri[0], &bufTri[1], &bufTri[2], &bufTri[3]);
    msh->Tri[i].Ver[0] = bufTri[0];
    msh->Tri[i].Ver[1] = bufTri[1];
    msh->Tri[i].Ver[2] = bufTri[2];
    msh->Tri[i].Ref    = bufTri[3];
  }
  
  
  //--- read boundary edges
  if ( readEfr == 1 ) {
    msh->NbrEfr = GmfStatKwd(fmsh, GmfEdges);
    msh->Efr    = calloc( (msh->NbrEfr+1), sizeof(Edge) );  

    GmfGotoKwd(fmsh, GmfEdges);
    for (i=1; i<=msh->NbrEfr; ++i) {
      GmfGetLin(fmsh, GmfEdges, &bufEfr[0], &bufEfr[1], &bufEfr[2]);
      msh->Efr[i].Ver[0] = bufEfr[0];
      msh->Efr[i].Ver[1] = bufEfr[1];
      msh->Efr[i].Ref    = bufEfr[2];
    }
  }
  
  
  
  GmfCloseMesh(fmsh);
  
  return msh;
  
}



double * sol_read(char *file, int mshDim, int mshNbrSol)
{
  char   InpFil[1024];
  int    FilVer, SolTyp, NbrTyp, SolSiz, TypTab[ GmfMaxTyp ]; 
  float  bufFlt;
  double bufDbl;
  int    i, dim, nbrSol;
  
  int fsol = 0;
  
  if ( ! file ) return NULL;
  
  double * sol = NULL;
    
    
  //--- set file name 
  strcpy(InpFil, file);
  if ( strstr(InpFil,".sol") ) {
    if ( !(fsol = GmfOpenMesh(InpFil,GmfRead,&FilVer,&dim)) ) {
      return NULL;
    }    
  }
  else {
    strcat(InpFil,".solb");
    if ( !(fsol = GmfOpenMesh(InpFil,GmfRead,&FilVer,&dim)) ) {
      strcpy(InpFil,file);
      strcat(InpFil,".sol");
      if ( !(fsol = GmfOpenMesh(InpFil,GmfRead,&FilVer,&dim)) ) {
        return NULL;
      }    
    } 
  }
  
  printf(" File %s opened Dimension %d Version %d \n", InpFil, dim, FilVer);
  
  SolTyp = GmfSolAtVertices;		// read only sol at vertices
  nbrSol = GmfStatKwd(fsol, SolTyp, &NbrTyp, &SolSiz, TypTab);
	
	
  if ( nbrSol == 0 ) {
    printf("  ## WARNING: No SolAtVertices in the solution file !\n");
    return NULL;
  }
  if ( dim != mshDim ) {
    printf("  ## WARNING: WRONG DIMENSION NUMBER. IGNORED\n");
    return NULL;
  }
  if ( nbrSol != mshNbrSol ) {
    printf("  ## WARNING: WRONG SOLUTION NUMBER. IGNORED\n");
    return NULL;
  }
  if (  NbrTyp != 1 ) {
    printf("  ## WARNING: WRONG FIELD NUMBER. IGNORED\n");
    return NULL;
  }
  if ( TypTab[0] != GmfSca ) {
    printf("  ## WARNING: WRONG FIELD TYPE. IGNORED\n");
    return NULL;
  }
	
	sol = (double *)calloc(nbrSol+1, sizeof(double));
	
	
  GmfGotoKwd(fsol, SolTyp);

  for (i=1; i<=nbrSol; ++i) {
		if ( FilVer == GmfFloat ) {
	    GmfGetLin(fsol, SolTyp, &bufFlt);
	    sol[i] = (double)bufFlt;
		}
		else {
	    GmfGetLin(fsol, SolTyp, &bufDbl);
	    sol[i] = bufDbl;
		}
  }
	
  if ( !GmfCloseMesh(fsol) ) {
    fprintf(stderr, "  ## ERROR: Cannot close solution file %s ! \n", InpFil);
    //myexit(1);
  }
	
  return sol;	
}




int compar_vertex(const void *a, const void *b)
{
  Vertex *va = (Vertex *) a;
  Vertex *vb = (Vertex *) b;
  return ( vb->icrit - va->icrit );
}

int compar_triangle(const void *a, const void *b)
{
  Triangle *va = (Triangle *) a;
  Triangle *vb = (Triangle *) b;
  return ( vb->icrit - va->icrit );
}

double aire(Mesh* msh, int iTri){
  double x[3];
  double y[3];
  Triangle* K = (Triangle*) &msh->Tri[iTri];
  for(int i=0; i<3; i++){
    int I = K->Ver[i];
    Vertex* Si = (Vertex*) &msh->Ver[I];
    x[i] = Si->Crd[0];
    y[i] = Si->Crd[1];
  }
  
  return 0.5 * ( (x[1]-x[0])*(y[2]-y[0]) - (x[2]-x[0])*(y[1]-y[0]) );
}

double    Q1(Mesh* msh, int iTri){
  const double alpha = sqrt(3.0)/12.0;
  double x[3];
  double y[3];

  Triangle* K = (Triangle*)&msh->Tri[iTri];
  for(int i=0; i<3; i++){
    int I = K->Ver[i];
    Vertex* Si = (Vertex*)&msh->Ver[I];
    x[i] = Si->Crd[0];
    y[i] = Si->Crd[1];
    //printf("xi = %f yi = %f\n",x[i],y[i]);
  }
  double aire_k = (0.5 * ( (x[1]-x[0])*(y[2]-y[0]) - (x[2]-x[0])*(y[1]-y[0]) ));
  double sum_l = 0;
  for(int i=0; i<3; i++){
    int j = (i+1) % 3;
    sum_l += pow((x[i]-x[j]),2) + pow((y[i]-y[j]),2);
  }
  //printf("sum_l = %f aire = %f\n",sum_l,aire_k);
  return alpha*sum_l/aire_k;
}
double    Q2(Mesh* msh, int iTri){
  const double alpha = sqrt(3.0)/6.0;
  double x[3];
  double y[3];

  Triangle* K = (Triangle*)&msh->Tri[iTri];
  for(int i=0; i<3; i++){
    int I = K->Ver[i];
    Vertex* Si = (Vertex*)&msh->Ver[I];
    x[i] = Si->Crd[0];
    y[i] = Si->Crd[1];
  }
  double aire_k = (0.5 * ( (x[1]-x[0])*(y[2]-y[0]) - (x[2]-x[0])*(y[1]-y[0]) ));
  double max_l = 0;
  double sum_l = 0;
  for(int i=0; i<3; i++){
    int j = (i+1) % 3;
    double dl = sqrt(pow((x[i]-x[j]),2.0) + pow((y[i]-y[j]),2.0));
    sum_l += dl;
    if(dl>max_l) max_l = dl;
  }
  double rho = 2.0 * aire_k / sum_l;
  return alpha*max_l/rho;
}

int msh_reorder(Mesh *msh)
{
  int iVer;  // iTri, 
  
  if ( ! msh            ) return 0;
  if ( msh->NbrVer <= 0 ) return 0;
  
  /* compute bonding box */
  for (iVer=1; iVer<=msh->NbrVer; iVer++) {
    /* todo msh->bb : used to compute the Z-curve index */
    
  }
  
  for (iVer=1; iVer<=msh->NbrVer; iVer++) {
    msh->Ver[iVer].icrit  = rand();   /* change the randon  by Z  order */
    msh->Ver[iVer].idxNew = iVer;
    msh->Ver[iVer].idxOld = iVer;
  }
  //qsort(&msh->Ver[1], msh->NbrVer, sizeof(Vertex), compar_vertex); 

  
  /* update idxNew for vertices */
   
  
  /* re-assign triangles ids */
  
  
  /* sort triangles */
  //for (iTri=1; iTri<=msh->NbrTri; iTri++) {
  //  msh->Tri[iTri].icrit = rand();   /* change the randon  by an improved order */  
  //}
  //qsort(&msh->Tri[1], msh->NbrTri, sizeof(Triangle), compar_triangle);
  
  
  return 1;
}




int msh_write(Mesh *msh, char *file)
{
  int iVer, iTri;
  int FilVer = 2;
  
  if ( ! msh  ) return 0;
  if ( ! file ) return 0;
  
  int fmsh = GmfOpenMesh(file, GmfWrite, FilVer, msh->Dim);
  if ( fmsh <=  0 ) {
    printf("  ## ERROR: CANNOT CREATE FILE \n");
    return 0;
  }
  
  GmfSetKwd(fmsh, GmfVertices, msh->NbrVer);
  for (iVer=1; iVer<=msh->NbrVer; iVer++) 
    GmfSetLin(fmsh, GmfVertices, msh->Ver[iVer].Crd[0], msh->Ver[iVer].Crd[1], 0); 
  
  GmfSetKwd(fmsh, GmfTriangles, msh->NbrTri);
  for (iTri=1; iTri<=msh->NbrTri; iTri++)  
    GmfSetLin(fmsh, GmfTriangles, msh->Tri[iTri].Ver[0], msh->Tri[iTri].Ver[1], msh->Tri[iTri].Ver[2], msh->Tri[iTri].Ref);  
     
  GmfCloseMesh(fmsh);
  
  return 1;
}





int msh_neighborsQ2(Mesh *msh)
{
  int iTri, iEdg, jTri, jEdg, ip1, ip2, jp1, jp2;
  
  if ( ! msh ) return 0;
  
  for (iTri=1; iTri<=msh->NbrTri; iTri++) {
    //printf("voi = %d %d %d\n",msh->Tri[iTri].Voi[0],msh->Tri[iTri].Voi[1],msh->Tri[iTri].Voi[2]);
    for (iEdg=0; iEdg<3; iEdg++) {
      ip1 = msh->Tri[iTri].Ver[tri2edg[iEdg][0]];
      ip2 = msh->Tri[iTri].Ver[tri2edg[iEdg][1]];
			
      // find the Tri different from iTri that has ip1, ip2 as vertices 
      for (jTri=1; jTri<=msh->NbrTri; jTri++) {
        if ( iTri == jTri ) 
          continue;
        
        for (jEdg=0; jEdg<3; jEdg++) {
          jp1 = msh->Tri[jTri].Ver[tri2edg[jEdg][0]];
          jp2 = msh->Tri[jTri].Ver[tri2edg[jEdg][1]];
          
          // compare the 4 points 
          if(((ip1 == jp1)&&(ip2 == jp2))||((ip1 == jp2)&&(ip2 == jp1))){
            msh->Tri[iTri].Voi[iEdg] = jTri;
            msh->Tri[jTri].Voi[jEdg] = iTri;
            //printf("Tri%d <-> Tri%d with Edge(%d,%d)\n",iTri,jTri,ip1,ip2);
            break;
          }
        }
      }
      
    }
  }
  
  return 1;
}

HashTable* hashTable_Voi(Mesh *msh){
  int iTri, iEdg, ip1, ip2;
  
  if ( ! msh ) return 0;
  
  /* initialize HashTable and set the hash table */
  const int SizeHead = 2*msh->NbrVer;
  const int NbrMaxObj = 3*msh->NbrTri;
  HashTable * hsh = hash_init(SizeHead,NbrMaxObj);

  for(iTri=1; iTri<=msh->NbrTri; iTri++){
    
    for (iEdg=0; iEdg<3; iEdg++){
      ip1 = msh->Tri[iTri].Ver[tri2edg[iEdg][0]];
      ip2 = msh->Tri[iTri].Ver[tri2edg[iEdg][1]];
      hash_add(hsh,ip1,ip2,iTri);
    }
  }

  return hsh;
}



int msh_neighbors(Mesh *msh)
{
  int iTri, iEdg, ip1, ip2;
  
  if ( ! msh ) return 0;
  
  /* initialize HashTable and set the hash table */
  
  HashTable * hsh = hashTable_Voi(msh);
  EdgeHashTable* ehsh = hashTable_Efr(msh);

  for(iTri=1; iTri<=msh->NbrTri; iTri++){
    
    for (iEdg=0; iEdg<3; iEdg++){
      if(msh->Tri[iTri].Voi[iEdg] > 0) continue;
      ip1 = msh->Tri[iTri].Ver[tri2edg[iEdg][0]];
      ip2 = msh->Tri[iTri].Ver[tri2edg[iEdg][1]];
      int index_key = hash_find(hsh, ip1, ip2);
      int* tmp = (int*)&hsh->LstObj[index_key];
      int t1 = tmp[2];
      int t2 = tmp[3];

      if(edge_hash_find(ehsh, ip1,ip2)==0){ msh->Tri[iTri].Voi[iEdg] = t1 + t2 - iTri; 
      //printf("T%d <--> T%d -- E%d\n",iTri,t1 + t2 - iTri,iEdg);
      }

      //msh->Tri[t2].Voi[iEdg] = t1;
    }
  }

  return 1;

}



static inline int key_function(int ip1, int ip2){
  return ip1+ip2; //ip1<ip2?ip1:ip2; //On peut esayer d'autre façon.
}


HashTable * hash_init(int SizHead, int NbrMaxObj)
{
	HashTable *hsh = NULL;
	
	// to be implemented

	// allocate hash table
	
	// initialize hash table
	
	// allocate Head, LstObj

  hsh = (HashTable *) malloc(sizeof(HashTable));
  if ( ! hsh ) return NULL;

  hsh->SizHead = SizHead;
  hsh->NbrMaxObj = NbrMaxObj;
  hsh->NbrObj = 0;
  hsh->Head = (int*) malloc(sizeof(int) * SizHead);
  memset(hsh->Head,EMPTYOBJECT,sizeof(int) * SizHead);
  hsh->LstObj = (int5*) calloc(NbrMaxObj + 1, sizeof(int5));
  hsh->NbrBoundaryEdge = 0;
  hsh->NbrEdge = 0;
	
	
  return hsh;
}


int hash_find(HashTable *hsh, int ip1, int ip2)
{
  
	// to be implemented
	
	// return the id found (in LstObj ), if 0 the object is not in the list
  int key = key_function(ip1,ip2) % hsh->SizHead;

  int head = hsh->Head[key];

  if(head == EMPTYOBJECT) return 0;

  int jp1 = ip1<ip2?ip1:ip2; //jp1 = min(ip1,ip2)
  int jp2 = ip1 + ip2 - jp1;//jp2 = max(ip1,ip2)

  while(head != EMPTYOBJECT){
    int* tmp = (int*)&hsh->LstObj[head];
    if((tmp[0] == jp1)&&(tmp[1] == jp2)) {return head;}
    head = tmp[4];
  }

  return 0;

}


int hash_add(HashTable *hsh, int ip1, int ip2, int iTri)
{

  // to be implemented
	
  // ===> add this entry in the hash tab 

  int key = key_function(ip1,ip2) % hsh->SizHead;
  int jp1 = ip1<ip2?ip1:ip2; //jp1 = min(ip1,ip2)
  int jp2 = ip1 + ip2 - jp1;//jp2 = max(ip1,ip2)
  int head = hsh->Head[key];

  //if(iTri == 57) printf("Tri57\n");
  //if(iTri == 57) printf("jp1 = %d jp2 = %d\n",jp1,jp2);
  //if(iTri == 57) printf("head = %d\n",head);

  if(head  == EMPTYOBJECT){ //aucun objet pour cette cle
    //if(iTri == 57) printf("head  == EMPTYOBJECT\n");
    int* tmp = (int*)&hsh->LstObj[hsh->NbrObj + 1];
    tmp[0] = jp1;
    tmp[1] = jp2;
    tmp[2] = iTri;
    tmp[3] = EMPTYOBJECT;
    tmp[4] = EMPTYOBJECT;
    hsh->Head[key] = hsh->NbrObj + 1;
    //if(iTri == 57) printf("head = hsh->NbrObj + 1");
    hsh->NbrObj++;
    hsh->NbrBoundaryEdge++;
    hsh->NbrEdge++;
    return 1;
  }

  
  int headNext = head;
  //if(iTri == 57) printf("begin while\n");
  while (headNext != EMPTYOBJECT)
  {
    int* tmp = (int*)&hsh->LstObj[headNext];
    if((tmp[1] == jp2)&&(tmp[0] == jp1)){
      //if(iTri == 57) printf("(tmp[1] == jp2)&&(tmp[0] == jp1)\n");
      //if(iTri == 57) printf("tmp[2] = %d tmp[3] = %d\n",tmp[2],tmp[3]);
      tmp[3] = iTri;
      hsh->NbrBoundaryEdge--;
      return 1;
    }
    head = headNext;
    headNext = hsh->LstObj[head][4];
  }
  //if(iTri == 57) printf("end while\n");

  int* tmp = (int*)&hsh->LstObj[head];
  if((tmp[1] == jp2)&&(tmp[0] == jp1)){
    //if(iTri == 57) printf("(tmp[1] == jp2)&&(tmp[0] == jp1)\n");
    //if(iTri == 57) printf("tmp[2] = %d tmp[3] = %d\n",tmp[2],tmp[3]);
    tmp[3] = iTri;
    hsh->NbrBoundaryEdge--;
    return 1;
  }

  if(hsh->NbrObj == hsh->NbrMaxObj + 1) {
    printf("NbrObj > NbrMaxObj\n");
    return 0;
  }


  //if(iTri == 57) printf("tmp2\n");
  int* tmp2 = (int*)&hsh->LstObj[hsh->NbrObj + 1];
  tmp2[0] = jp1;
  tmp2[1] = jp2;
  tmp2[2] = iTri;
  tmp2[3] = EMPTYOBJECT;
  tmp2[4] = EMPTYOBJECT;
  
  tmp[4] = hsh->NbrObj + 1;
  //if(iTri == 57) printf("tmp[4] = %d\n",hsh->NbrObj + 1);
  hsh->NbrBoundaryEdge++;
  hsh->NbrEdge++;
  hsh->NbrObj++;
  
	return 1;
}


DL* dl_init(int value){
  DL* res = (DL*)malloc(sizeof(DL));
  res->previous_node = NULL;
  res->next_node = NULL;
  res->value = value;
  return res;
}

int dl_append(DL* dl,int value){
  if(dl==NULL) return 0;
  DL* newDl = dl_init(value);
  if(newDl==NULL) return 0;
  DL* next = dl;
  while (next->next_node != NULL)
  {
    ;
  }
  next->next_node = newDl;
  newDl->previous_node = next;
  return 1;  

}

Stack* stack_init(){
  Stack* res = (Stack*)malloc(sizeof(Stack));
  res->NbrObj = 0;
  res->top = NULL;
  return res;
}

int stack_push(Stack* sk, int value){
  if(sk==NULL) return 0;
  if(sk->top==NULL){
    DL* dl = dl_init(value);
    if(dl==NULL) return 0;
    sk->top = dl;
    sk->NbrObj ++;
    return 1;
    }

  int res = dl_append(sk->top,value);
  if(res ==0) return 0;
  sk->top = sk->top->next_node;
  sk->NbrObj++;
  return 1;
}
int stack_pop(Stack* sk){
  if(sk==NULL){
    printf("NULL Stack\n");
    return 0;
  };
  if(sk->NbrObj==0){
    printf("Empty Stack\n");
    return 0;
  }
  DL* old_top = sk->top;
  int res = old_top->value;
  sk->top = old_top->previous_node;
  if(sk->top != NULL){
    sk->top->next_node=NULL;
  }
  free(old_top);
  sk->NbrObj--;
  return res;
}

EdgeHashTable* edge_hash_init(int SizeHead, int NbrMaxObj){
  EdgeHashTable *hsh = NULL;
	
  hsh = (EdgeHashTable *) malloc(sizeof(EdgeHashTable));
  
  if ( ! hsh ) return NULL;
  
  hsh->SizeHead = SizeHead;
  hsh->NbrMaxObj = NbrMaxObj;
  hsh->NbrObj = 0;
  hsh->Head = (int*) malloc(sizeof(int) * SizeHead);
  memset(hsh->Head,EMPTYOBJECT,sizeof(int) * SizeHead);
  hsh->LstObj = (int3*) calloc(NbrMaxObj + 1, sizeof(int3));
	
  return hsh;
}
int edge_hash_find(EdgeHashTable *hsh, int ip1, int ip2){
  int key = key_function(ip1,ip2) % hsh->SizeHead;

  int head = hsh->Head[key];

  if(head == EMPTYOBJECT) return 0;

  int jp1 = ip1<ip2?ip1:ip2; //jp1 = min(ip1,ip2)
  int jp2 = ip1 + ip2 - jp1;//jp2 = max(ip1,ip2)

  while(head != EMPTYOBJECT){
    int* tmp = (int*)&hsh->LstObj[head];
    if((tmp[0] == jp1)&&(tmp[1] == jp2)) {return head;}
    head = tmp[2];
  }

  return 0;

}

int edge_hash_add(EdgeHashTable *hsh, int ip1, int ip2){
  int key = key_function(ip1,ip2) % hsh->SizeHead;
  int jp1 = ip1<ip2?ip1:ip2; //jp1 = min(ip1,ip2)
  int jp2 = ip1 + ip2 - jp1;//jp2 = max(ip1,ip2)
  int head = hsh->Head[key];


  if(head  == EMPTYOBJECT){ //aucun objet pour cette cle
;
    int* tmp = (int*)&hsh->LstObj[hsh->NbrObj + 1];
    tmp[0] = jp1;
    tmp[1] = jp2;
    tmp[2] = EMPTYOBJECT;
  
    hsh->Head[key] = hsh->NbrObj + 1;

    hsh->NbrObj++;
    return 1;
  }

  
  int headNext = head;

  while (headNext != EMPTYOBJECT)
  {
    int* tmp = (int*)&hsh->LstObj[headNext];
    if((tmp[1] == jp2)&&(tmp[0] == jp1)){
      return 1;
    }
    head = headNext;
    headNext = hsh->LstObj[head][2];
  }


  int* tmp = (int*)&hsh->LstObj[head];
  if((tmp[1] == jp2)&&(tmp[0] == jp1)){
    return 1;
  }

  if(hsh->NbrObj == hsh->NbrMaxObj + 1) {
    printf("NbrObj > NbrMaxObj\n");
    return 0;
  }


  //if(iTri == 57) printf("tmp2\n");
  int* tmp2 = (int*)&hsh->LstObj[hsh->NbrObj + 1];
  tmp2[0] = jp1;
  tmp2[1] = jp2;
  tmp2[2] = EMPTYOBJECT;
  
  tmp[2] = hsh->NbrObj + 1;
  //if(iTri == 57) printf("tmp[4] = %d\n",hsh->NbrObj + 1);

  hsh->NbrObj++;
  
	return 1;
}

EdgeHashTable* hashTable_Efr(Mesh *msh){

  int iEdg, ip1, ip2;
  
  if ( ! msh ) return 0;

  /* initialize HashTable and set the hash table */
  const int SizeHead = 2*msh->NbrVer;
  const int NbrMaxObj = msh->NbrEfr;

  EdgeHashTable * hsh = edge_hash_init(SizeHead,NbrMaxObj);
 
  for(iEdg = 0; iEdg < msh->NbrEfr; iEdg++){
    Edge* tmp = &msh->Efr[iEdg];
    ip1 = tmp->Ver[0];
    ip2 = tmp->Ver[1];
    edge_hash_add(hsh,ip1,ip2);
  }

  return hsh;
}

Vector* vector_init(){
  Vector* res = (Vector*) malloc(sizeof(Vector));
  res->length = 0;
  return res;
}

int vector_append(Vector* v, int val){
  if(! v) return 0;
  if(v->length==0){
    DL* dl = dl_init(val);
    if(! dl) return 0;
    v->head = dl;
    v->tail = dl;
    v->length ++;
    return 1;
  }
  if(dl_append(v->tail, val)){
    v->tail = v->tail->next_node;
    v->length ++;
    return 1;
  }
  
  return 0;
}
DL* vector_elementAt(Vector* v, int pos){
  if( ! v) return NULL;
  if(pos < 0){
    printf("pos < 0\n");
    return NULL;
  }
  if(pos >= v->length){
    printf("pos >= length\n");
    return NULL;
  }
  if(pos < v->length / 2){
    DL* res = v->head;
    for(int i = 0; i < pos; i++){
      res = res->next_node;
    }
    return res;
  }
  DL * res = v->tail;
  for(int i = 0; i < v->length - 1 - pos; i++){
    res = res->previous_node;
  }
  return res;
}
int vector_inser(Vector* v, int pos, int val){
  if(! v) return 0;
  if(pos < 0){
    printf("pos < 0\n");
    return 0;
  }
  if(pos > v->length){
    printf("pos > length\n");
    return 0;
  }
  if(pos == 0){
    DL * dl = dl_init(val);
    if( ! dl) return 0;
    dl->next_node = v->head;
    v->head->previous_node = dl;
    v->head = dl;
    v->length ++;
    return 1;
  }
  if(pos == v->length){
    return vector_append(v,val);
  }
  DL* tmp = vector_elementAt(v,pos - 1);
  DL* dl = dl_init(val);
  if( ! dl) return 0;
  dl->next_node = tmp->next_node;
  tmp->next_node->previous_node = dl;
  dl->previous_node = tmp;
  tmp->next_node = dl;
  v->length ++;
  return 1;
}

int vector_remove(Vector* v, int pos){
  if(! v) return 0;
  if(pos < 0){
    printf("pos < 0\n");
    return 0;
  }
  if(pos >= v->length){
    printf("pos >= length\n");
    return 0;
  }
  if(pos == 0){
    DL* tmp =v->head->next_node;
    tmp->previous_node = NULL;
    free(v->head);
    v->head = tmp;
    v->length --;
    return 1;
  }
  if(pos == v->length - 1){
    DL* tmp = v->tail->previous_node;
    tmp->next_node = NULL;
    free(v->tail);
    v->tail = tmp;
    v->length --;
    return 1;
  }
  DL* tmp = vector_elementAt(v,pos);
  tmp->previous_node->next_node = tmp->next_node;
  tmp->next_node->previous_node = tmp->previous_node;
  free(tmp);
  v->length --;
  return 1;
}

int vector_valueAt(Vector* v, int pos){
  if(! v){
    printf("V = NULL\n");
    return 0;
  }
  DL* tmp = vector_elementAt(v, pos);
  if(! tmp){
    printf("V[pos] = NULL\n");
    return 0;
  }
  return tmp->value;
}

double * colorisation(Mesh *msh){
  const int nbrTri = msh->NbrTri + 1;
  double* res = (double*) calloc(nbrTri, sizeof(double));
  double color = 1.0;

  HashTable* hsh = hashTable_Voi(msh);
  Stack * sk = stack_init();
  for(int index_efr=0; index_efr < msh->NbrEfr; index_efr++){
    Edge* tmp = &msh->Efr[index_efr];
    int ip1 = tmp->Ver[0];
    int ip2 = tmp->Ver[1];
    int index_key = hash_find(hsh,ip1,ip2);
    int iTris[2] = {hsh->LstObj[index_key][2],hsh->LstObj[index_key][3]};
    for(int index_tri = 0; index_tri < 2; index_tri++){
      int iTri = iTris[index_tri];
      if(iTri==0) continue;
      if(res[iTri]>0) continue;
      //printf("iTri = %d\n",iTri);
      stack_push(sk,iTri);
      
      while (sk->NbrObj>0)
      {
        int top = stack_pop(sk);
        //printf("top = %d\n",top);
        res[top] = color;
        //printf("T%d = %d\n",top,color);
        Triangle* T = &msh->Tri[top];
        for(int index_voi = 0; index_voi < 3; index_voi++){
          int voiT = T->Voi[index_voi];
          if(voiT==0) continue;
          //printf("voi = %d avec c = %d\n",voiT,res[voiT]);
          if(res[voiT] > 0) continue;
          
          stack_push(sk,voiT);
        }
      }
      color += 1.0;
      
    }
    
  }

  return res;

  

}


int    nb_boundary_edge(Mesh *msh){
  //printf("###############################################\n");
  HashTable * hsh = hashTable_Voi(msh);
  //printf("###############################################\n");
  printf("nombre total d’arêtes = %d\n",hsh->NbrEdge);
  int s = 0;
  int nbObj = hsh->NbrObj + 1;
  const int nbTri = msh->NbrTri;
  int info[nbTri+1];
  memset(info,0,sizeof(int) * (nbTri+1));
  for(int i=1; i<nbObj; i++){
    int* tmp = (int*) &hsh->LstObj[i];
    if((tmp[3] == EMPTYOBJECT)&&(tmp[2] != EMPTYOBJECT)) s++;
    info[tmp[2]]++;
    info[tmp[3]]++;
    //if((tmp[3] == 57)||(tmp[2] == 57)) printf("S1 = %d S2 = %d T1 = %d, T2 = %d\n",tmp[0],tmp[1],tmp[2],tmp[3]);
  }
  for(int i=0;i<nbTri+1;i++){
    //printf("nb of Tri%d = %d\n",i,info[i]);
  }
  /*int head = 62;
  
  while(head != EMPTYOBJECT){
    printf("head = %d\n",head);
    int * tmp = (int*)&hsh->LstObj[head];
    printf("L62 : %d %d %d %d %d\n",tmp[0],tmp[1],tmp[2],tmp[3],tmp[4]);
    head = tmp[4];
  }*/
 
  return hsh->NbrBoundaryEdge;


}



int msh_write2dfield_Vertices(char *file, int nfield, double *field) 
{
  int iVer;
  
  int fmsh = GmfOpenMesh(file, GmfWrite, GmfDouble, 2);
  if ( fmsh <=  0 ) {
    printf("  ## ERROR: CANNOT CREATE FILE \n");
    return 0;
  }
  
  int sizfld[1];
  sizfld[0] = GmfSca;
  
  GmfSetKwd(fmsh, GmfSolAtVertices, nfield, 1, sizfld);
  
  for (iVer=1; iVer<=nfield; iVer++) 
    GmfSetLin(fmsh, GmfSolAtVertices, &field[iVer]); 
  
  GmfCloseMesh(fmsh);
  
  return 1;
}



int msh_write2dfield_Triangles(char *file, int nfield, double *field) 
{
  int iTri;
  
  int fmsh = GmfOpenMesh(file, GmfWrite, GmfDouble, 2);
  if ( fmsh <=  0 ) {
    printf("  ## ERROR: CANNOT CREATE FILE \n");
    return 0;
  }
  
  int sizfld[1];
  sizfld[0] = GmfSca;
  
  GmfSetKwd(fmsh, GmfSolAtTriangles, nfield, 1, sizfld);
  
  for (iTri=1; iTri<=nfield; iTri++) 
    GmfSetLin(fmsh, GmfSolAtTriangles, &field[iTri]); 
  
  GmfCloseMesh(fmsh);
  
  return 1;
}



int msh_write2dmetric(char *file, int nmetric, double3 *metric) 
{  
  int iVer;
  
  int fmsh = GmfOpenMesh(file, GmfWrite, GmfDouble, 2);
  if ( fmsh <=  0 ) {
    printf("  ## ERROR: CANNOT CREATE FILE \n");
    return 0;
  }
  
  int sizfld[1];
  sizfld[0] = GmfSymMat;
  
  GmfSetKwd(fmsh, GmfSolAtVertices, nmetric, 1, sizfld);
  
  for (iVer=1; iVer<=nmetric; iVer++) 
    GmfSetLin(fmsh, GmfSolAtVertices, &metric[iVer][0], &metric[iVer][1], &metric[iVer][2]); 
  
  GmfCloseMesh(fmsh);
  
  return 1;
}



double2* points_unif(double xmin,double ymin,double xmax,double ymax,int Nx, int Ny){
  if( (xmax<=xmin) || (ymax<=ymin) || (Nx <= 0) || (Ny <= 0)){
    printf("Error\n");
    return NULL;
  }

  double dx = (xmax - xmin)/ (double)(Nx - 1);
  double dy = (ymax - ymin)/ (double)(Ny - 1);

  double2* res = (double2 *)malloc(sizeof(double2) * Nx * Ny);
  if(res==NULL) printf("malloc error\n");
  
  int index_res = 0;
  for(int index_x = 0; index_x < Nx; index_x++){
    for(int index_y = 0; index_y < Ny; index_y++){
      res[index_res][0] = xmin + ((double)index_x)*dx;
      res[index_res][1] = ymin + ((double)index_y)*dy;
      index_res ++;
    }
  }
  
  return res;
}

static inline double aire_signe(double2* s){
  //s[i][0] --> xi; s[j][1] --> yj;
  //(x[1]-x[0])*(y[2]-y[0]) - (x[2]-x[0])*(y[1]-y[0])
  return 0.5 * ( (s[1][0]-s[0][0])*(s[2][1]-s[0][1]) - (s[2][0]-s[0][0])*(s[1][1]-s[0][1]) );
}

static inline double2* get_sommets_T(Triangle* T, Mesh* msh){
  double2* res = (double2*) malloc(sizeof(double2)*3);
  int* ver = (int*) T->Ver;
  //printf("sommets = \n");
  for(int q=0; q<3; q++){
    //printf("%d, ",ver[q]);
    double* ver_q = (double*) msh->Ver[ver[q]].Crd;
    res[q][0] = ver_q[0];
    res[q][1] = ver_q[1];
  }
  //printf("\n");
  return res;
}

static inline int* signe_P(Triangle* T, Mesh * msh, double2 P){
  double aire_T = aire_signe(get_sommets_T(T,msh));
  if(aire_T==0){
    printf("aire nulle\n");
    return NULL;
  }
  //printf("aire = %f\n",aire_T);
  int* res = (int*) malloc(sizeof(int) * 3);
  for(int q=0; q<3; q++){
    double2* s_T = get_sommets_T(T,msh);
    s_T[q][0] = P[0];
    s_T[q][1] = P[1];
    double aire_q = aire_signe(s_T);
    //printf("x = %f %f %f\ny = %f %f %f\n",s_T[0][0],s_T[1][0],s_T[2][0],s_T[0][1],s_T[1][1],s_T[2][1]);
    //printf("a = %f, b = %f, a*b>=0 = %d, sgn = %d\n",aire_T,aire_q,aire_q * aire_T >= 0,2*(aire_q * aire_T >= 0) - 1);
    res[q] = 2*(aire_q * aire_T >= 0) - 1; // res = -1 si aire_q et aire_T n'ont pas le meme signe. res = 1 sinon.
  }
  return res;
}

static inline int eva_sgn(int* sgn_P){
  int* sigma = malloc(sizeof(int)*3);
  srand((unsigned)time(NULL));
  sigma[0] = rand()%3; //x = 0, 1 ,2
  const int s = sigma[0] + 1;
  sigma[1] = (s + rand()%2) %3;
  sigma[2] = 3 - sigma[1] - sigma[0];
  //printf("sig = %d %d %d\n",sigma[0],sigma[1],sigma[2]);
  for(int q=0; q<3;q++){
    int sq = sigma[q]; //sq = 0 , 1, 2 
    //printf("sgn[%d] = %d\n",sq,sgn_P[sq]);
    if(sgn_P[sq]<0) return sq;
  }
  return -1; //P est dans interieur de T
}

static inline int localiser(Mesh *msh, int K, double2 P){
  if((msh==NULL) || (P == NULL)){
    //printf("Error\n");
    return 0;
  }
  if((K<=0)||(K>=msh->NbrTri)){
    printf("Error\n");
    return 0;
  }
  Triangle* T = (Triangle*) &msh->Tri[K];
  //printf("T%d avec S%d, %d, %d\n",K,T->Ver[0],T->Ver[1],T->Ver[2]);

  int* sgn_P = signe_P(T,msh,P);
  //printf("@sgn = %d",sgn_P);
  //printf("T%d\n",K);
  //printf("sgn = %d %d %d\n",sgn_P[0],sgn_P[1],sgn_P[2]);

  int direction = eva_sgn(sgn_P);

  //printf("direction = %d\n",direction);

  //printf("voi = %d %d %d\n",T->Voi[0],T->Voi[1],T->Voi[2]);

  while (direction >= 0)
  {
    K = T->Voi[direction];
    //printf("T%d\n",K);
    T = (Triangle*) &msh->Tri[K];
    //printf("T%d avec S%d, %d, %d\n",K,T->Ver[0],T->Ver[1],T->Ver[2]);
    sgn_P = signe_P(T,msh,P);
    direction = eva_sgn(sgn_P);
    

    //printf("sgn = %d %d %d\n",sgn_P[0],sgn_P[1],sgn_P[2]);

    //printf("direction = %d\n",direction);

    //system("read -p 'Press Enter to continue...' var");
  }

  return K;
  
}

Circle* circle_init(){
  Circle* res = (Circle*) malloc(sizeof(Circle));
  res->r = 0;
  res->centre = (double*)malloc(2 * sizeof(double));
  return res;
}

Set* set_init(){
  Set* res = (Set*)malloc(sizeof(Set));
  res->size = 0;
  return res;
}
int set_add(Set* s, int value){
  if(! s) return 0;
  if(s->size==0){
    DL* dl = dl_init(value);
    if(! dl) return 0;
    s->elem = dl;
    s->size ++;
    return 1;
  }
  DL* elem = s->elem;
  if(! elem) return 0;
  while(elem->next_node != NULL){
    if(elem->value == value) return 1;
    elem = elem->next_node;
  }
  if(elem->value == value) return 1;
  DL * new_dl = dl_init(value);
  if(! new_dl) return 0;
  elem->next_node = new_dl;
  new_dl->previous_node = elem;
  s->size ++;
  return 1;
}

int set_pop(Set* s){
  if(! s){
    printf("set = NULL\n");
    return 0;
  }
  if(s->size == 0){
    printf("size = 0\n");
    return 0;
  }
  if(! s->elem){
    printf("elem = NULL\n");
    return 0;
  }
  int res = s->elem->value;
  DL* tmp = s->elem->next_node;
  free(s->elem);
  s->size --;
  s->elem = tmp;
  if(s->elem != NULL) s->elem->previous_node = NULL;
  return res;
}

int set_remove(Set* s, int value){
  if(! s) return 0;
  DL* elem = s->elem;
  while (elem != 0)
  {
    if(elem->value == value){
      DL* tmp1 = elem->previous_node;
      DL* tmp2 = elem->next_node;
      if(tmp1 != 0) tmp1->next_node = tmp2;
      if(tmp2 != 0) tmp2->previous_node = tmp1;
      if(s->elem == elem) s->elem = tmp2;
      free(elem);
      s->size --;
      return 1;
    }
    elem = elem->next_node;
  }

  //printf("Il n'y existe pas cet element\n");
  return 0;
}

int set_contains(Set* s, int value){
  if(! s) return 0;
  if(s->size ==0) return 0;
  for(DL* it = s->elem; it != NULL; it = it->next_node){
    if(it->value == value) return 1;
  }
  return 0;
}



static inline Circle* cercle_circonscrit(double2* s){
  if(s==NULL){
    printf("Null pointer\n");
    return NULL;
  }
  //printf("(%f,%f), (%f,%f), (%f,%f)\n",s[0][0],s[0][1],s[1][0],s[1][1],s[2][0],s[2][1]);
  double aire_k = aire_signe(s);
  if(aire_k ==0){
    printf("aire = 0\n");
    return NULL;
  }
  Circle * res = circle_init();
  double r = 1;
  for (int n = 0; n < 3; n++){
    r *=  (pow((s[(n + 1) % 3][0] - s[n][0]),2) + pow((s[(n + 1) % 3][1] - s[n][1]),2));
  }
  r = sqrt(r) / (4.0 * aire_k);
  res->r = fabs(r);
  double A = pow(s[0][0],2) + pow(s[0][1],2);
  double B = pow(s[1][0],2) + pow(s[1][1],2);
  double C = pow(s[2][0],2) + pow(s[2][1],2);

  res->centre[0] = ((A - B) * (s[0][1] - s[2][1]) - (s[0][1] - s[1][1]) * (A - C)) / (4.0 * aire_k);
  res->centre[1] = ((s[0][0] - s[1][0]) * (A - C)  - (A - B) * (s[0][0] - s[2][0])) / (4.0 * aire_k);
  return res;
}

static inline Vector* get_cavity(double2 P, int K_P, Mesh* msh){
  Stack* sk = stack_init();
  const int nbrTri = msh->NbrTri;
  int* visited = (int*) calloc(nbrTri+1, sizeof(int));
  stack_push(sk,K_P);
  Vector* res = vector_init();
  while (sk->NbrObj>0)
  {
    int K = stack_pop(sk);
    
    //printf("getCav T = %d\n",K);
    if(visited[K] == 1) continue;
    visited[K] = 1;
    Triangle* T = (Triangle*)&msh->Tri[K];
    Circle* circ = cercle_circonscrit(get_sommets_T(T,msh));
    //printf("O = (%f,%f), r = %f\n",circ->centre[0],circ->centre[1],circ->r);
    //printf("d = %f\n",sqrt(pow(P[0] -circ->centre[0],2) + pow(P[1] -circ->centre[1],2)));
    if((pow(P[0] -circ->centre[0],2) + pow(P[1] -circ->centre[1],2)) > pow(circ->r,2)) continue;
    vector_append(res, K);
    for(int q=0; q<3; q++){
      int K_voi = T->Voi[q];
      //printf("K_voi = %d, visted[K_voi] = %d\n",K_voi,visited[K_voi]);
      if((K_voi >0)&&(visited[K_voi] == 0)){ 
        stack_push(sk,K_voi);
        //printf("K_voi = %d, visted[K_voi] = %d\n",K_voi,visited[K_voi]);
      }
    }
  }
  free(visited);
  free(sk);
  return res;
}

void test(Mesh* msh){
  double2 P;
  P[0] =15.0/90.0;
  P[1] =15.0/90.0;
  printf("P = (%f,%f)\n",P[0],P[1]);
  int K = localiser(msh,1,P);
  printf("K_P = %d\n",K);
  Vector* cav = get_cavity(P,K,msh);
  printf("La cavité de P est :\n");
  for(int i = 0; i < cav->length; i++){
    printf("%d, ",vector_valueAt(cav,i));
  }
  printf("\n");
  
}
DMesh* dMesh_init(Mesh* msh, int lambda){
  if(lambda<1){{
    printf("Lambda must be >= 1\n");
    return NULL;
  }}
  if( ! msh ) return NULL;
  DMesh * res = (DMesh*)malloc(sizeof(DMesh));
  if(! res) return NULL;
  res->Dim = msh->Dim;
  res->NbrVerMax = (msh->NbrVer + 1) * lambda;
  res->NbrTriMax = (msh->NbrTri + 1) * lambda;
  msh->Ver = (Vertex *)realloc(msh->Ver,sizeof(Vertex)*res->NbrVerMax);
  msh->Tri = (Triangle*)realloc(msh->Tri,sizeof(Triangle)*res->NbrTriMax);
  res->msh = msh;
  return res;
}

int dMesh_resize(DMesh* dMsh, int NbrVerMax, int NbrTriMax){
  if(! dMsh) return 0;
  if((NbrVerMax<=0)||(NbrTriMax<=0)){
    printf("NbrVerMax and NbrTriMax must be >= 1\n");
    return 0;
  }

  //printf("newNbrVerMax= % d,  newNbrTriMax= %d\n",NbrVerMax,NbrTriMax);
  //printf("Resize NbrVerMax %d -> %d\n",dMsh->NbrVerMax,NbrVerMax);
  //printf("@Vertex = %p\n",dMsh->msh->Ver);
  dMsh->msh->Ver = (Vertex*)realloc(dMsh->msh->Ver,sizeof(Vertex) * (NbrVerMax));
  //printf("@Vertex = %p\n",dMsh->msh->Ver);
  if(! dMsh->msh->Ver){
    printf("Ver = NULL\n");
    return 0;
  }
  dMsh->NbrVerMax = NbrVerMax;

  //printf("Resize NbrTriMax %d -> %d\n",dMsh->NbrTriMax,NbrTriMax);
  dMsh->msh->Tri = (Triangle*)realloc(dMsh->msh->Tri,sizeof(Triangle) * (NbrTriMax));
  if(! dMsh->msh->Tri){
    printf("Tri = NULL\n");
    return 0;
  }
  dMsh->NbrTriMax = NbrTriMax;
  return 1;
}
static inline HashTable* hsh_cav(Vector* cav, Mesh* msh){
  if ( ! msh ) return NULL;
  if ( ! cav ) return NULL; 
  int iTri, iEdg, ip1, ip2;
  
  /* initialize HashTable and set the hash table */
  const int nbrTri = cav->length;
  const int SizeHead = 2*3*nbrTri;
  const int NbrMaxObj = 3*nbrTri;
  HashTable * hsh = hash_init(SizeHead,NbrMaxObj);

  for(int i = 0; i < cav->length; i++)
  {
    iTri = vector_valueAt(cav, i);
    for (iEdg=0; iEdg<3; iEdg++){
      ip1 = msh->Tri[iTri].Ver[tri2edg[iEdg][0]];
      ip2 = msh->Tri[iTri].Ver[tri2edg[iEdg][1]];
      hash_add(hsh,ip1,ip2,iTri);
    }
  }

  return hsh;
}


static inline int reinit_voi(Vector* cav, Mesh* msh){
  if(! cav){
    printf("cav ==Null\n");
    return 0;
  }
  if(! msh){
    printf("msh = Null\n");
    return 0;
  }
  Vector* tmp = vector_init();
  if(! tmp) return 0;
  for(DL* it = cav->head; it != NULL; it = it->next_node){
    if(it->value == 0) continue;
    vector_append(tmp,it->value);
    for(int i = 0; i <3; i++){
      int voi = msh->Tri[it->value].Voi[i];
      if(voi == 0) continue;
      vector_append(tmp,voi);
    }
  }
  for(DL* it = tmp->head; it != NULL; it= it->next_node){
    for(int i =0; i<3; i++){
      msh->Tri[it->value].Voi[i] = 0;
    }
  }
  return 1;
}

static inline int* order_Tri(int* ips, Mesh* msh){
  if(! msh){
    printf("msh = Null\n");
    return NULL;
  }
  if(! ips){
    printf("ips = Null\n");
    return NULL;
  }
  
  int* res = (int*)malloc(sizeof(int)*3);
  if(! res) return NULL;
  double2* P = (double2*)malloc(sizeof(double2)* 3);
  for(int i = 0; i < 3 ; i++){
    res[i] = i;
    P[i][0] = msh->Ver[ips[i]].Crd[0];
    P[i][0] = msh->Ver[ips[i]].Crd[1];
  }

  if(aire_signe(P)< 0){
    res[0] = 1;
    res[1] = 0;
  }
  return res;
}

static inline int Bowyer_Watson(double2 P, DMesh* dMsh, Set* vertex_deleted){
  if(! P) return 0;
  if(! dMsh) return 0;
  if(! dMsh->msh) return 0;
  
  //printf("P = (%f,%f)\n",P[0],P[1]);

  //Localiser le point P.
  int K_P = localiser(dMsh->msh, 1, P);

  //printf("loc P = %d\n",K_P);
  //printf("S = %d %d %d\n",dMsh->msh->Tri[K_P].Ver[0],dMsh->msh->Tri[K_P].Ver[1],dMsh->msh->Tri[K_P].Ver[2]);

  //Trouver la cavite de P
  Vector* cav = get_cavity(P,K_P,dMsh->msh);
  if(!cav){
    printf("cav = NULL\n");
    return 0;
  }
  HashTable* hcav = hsh_cav(cav, dMsh->msh);
  if(!hcav){
    printf("hcav = NULL\n");
    return 0;
  }

  /*printf("cav = ");
  for(int q=0; q<cav->length;q++){
    printf("%d, ",vector_valueAt(cav,q));
  }
  printf("\n");*/

  if(reinit_voi(cav,dMsh->msh) == 0){
    printf("error reinit_voi \n");
      return 0;
  }

  Vector* edges_bord_cav = vector_init();

  for (DL* it = cav->head; it != NULL; it = it->next_node)
  {
    int K = it->value;
    /*if(! reinit_voi(K,dMsh->msh)){ 
      printf("error reinit_voi for T%d\n",K);
      return 0;
    }*/
    for(int q = 0; q < 3; q++){
      int ip1 = dMsh->msh->Tri[K].Ver[q];
      int ip2 = dMsh->msh->Tri[K].Ver[(q+1)%3];
      int tmp = hash_find(hcav,ip1,ip2);
      if(! tmp){
        printf("hsh error\n");
        return 0;
      }
      int* obj = hcav->LstObj[tmp];
      if((obj[2]==0||obj[3]==0)){
        vector_append(edges_bord_cav, ip1);
        vector_append(edges_bord_cav, ip2);
      }else{
        set_add(vertex_deleted, ip1);
        set_add(vertex_deleted, ip2);
      }
    }
    
  }

  if(edges_bord_cav->length < 2 * cav->length){
    printf("nbr_edges_bord_cav < nbrTri_cav\n");
    return 0;
  }


  for(DL* it = edges_bord_cav->head; it != NULL; it = it->next_node){
    int ip = it->value;
    set_remove(vertex_deleted, ip);
  }

  int numero_P = 0;
  

  if(vertex_deleted->size>0){
    numero_P = set_pop(vertex_deleted);
    dMsh->msh->NbrVer--;
  }else{
    numero_P = dMsh->msh->NbrVer + 1;
    if(dMsh->NbrVerMax <= numero_P){
      dMesh_resize(dMsh,2*dMsh->NbrVerMax,dMsh->NbrTriMax);
    }
  }

  if(! numero_P){
    printf("numero_P = 0\n");
    return 0;
  }

  //printf("numero P = %d\n",numero_P);

  dMsh->msh->Ver[numero_P].Crd[0] = P[0];
  dMsh->msh->Ver[numero_P].Crd[1] = P[1];
  //printf("P%d = (%f,%f)\n",numero_P,dMsh->msh->Ver[numero_P].Crd[0],dMsh->msh->Ver[numero_P].Crd[1]);
  dMsh->msh->NbrVer++;


  DL* tri_cav = cav->head;
  int numero_T = 0;

  for(DL* it = edges_bord_cav->head; it != NULL; it= it->next_node){
    int ip1 = it->value;
    it = it->next_node;
    if(! it){
      printf("iterator edges_bord_cav error\n");
      return 0;
    }
    int ip2 = it->value;
    //printf("ip1 = %d, ip2 = %d\n",ip1,ip2);
    if(! tri_cav){
      numero_T = dMsh->msh->NbrTri + 1;
      if(dMsh->NbrTriMax<=numero_T){
        dMesh_resize(dMsh,dMsh->NbrVerMax,2*dMsh->NbrTriMax);
      }
      dMsh->msh->NbrTri ++;
    }else{
      numero_T = tri_cav->value;
      tri_cav = tri_cav->next_node;
    }

     if(! numero_T){
     printf("numero_T = 0\n");
     return 0;
     }

     //printf("numero_T = %d\n",numero_T);

     int* ips = (int*)malloc(sizeof(int)*3);
     ips[0] = ip1;
     ips[1] = ip2;
     ips[2] = numero_P;
     int* order = order_Tri(ips, dMsh->msh);
     if(! order){
      printf("order = Null\n");
      return 0;
     }
     for(int j = 0; j < 3; j++){
      dMsh->msh->Tri[numero_T].Ver[j] = ips[order[j]];
      dMsh->msh->Tri[numero_T].Voi[j] = 0;
     }
     
     
  }
  return msh_neighbors(dMsh->msh);
}

int Triangulation(double2* Ps, int NbrP, DMesh* dmsh){
  if(! Ps){
    printf("Ps = Null\n");
    return 0;
  }
  if(NbrP<=0) return 0;
  if(! dmsh){
    printf("dMsh = Null\n");
    return 0;
  }

  Set* vertex_deleted = set_init();

  for(int i = 0; i < NbrP; i++){
    //printf("P_%d  = (%f,%f)\n",i,Ps[i][0],Ps[i][1]);
    if(! Bowyer_Watson(Ps[i], dmsh, vertex_deleted)){
      printf("error Bowyer_Watson for P_%d  = (%f,%f)\n",i,Ps[i][0],Ps[i][1]);
      return 0;
    }
    
  }

  if(vertex_deleted->size>0){
    printf("vertex_deleted > 0");
    return 0;
  }
  

  return 1;
}

Mesh* toMesh(DMesh* dmsh){
  Set* s = set_init();
  set_add(s,1);
  set_add(s,2);
  set_add(s,3);
  set_add(s,4);

  Mesh* msh = msh_init();
  if(! msh) return NULL;
  msh->Dim = dmsh->Dim;
  msh->NbrVer = dmsh->msh->NbrVer - 4;
  msh->NbrTri = 0;
  msh->Ver = calloc( (msh->NbrVer+1), sizeof(Vertex));
  msh->Tri = calloc( (dmsh->msh->NbrTri+1), sizeof(Triangle) ); 

  
  for(int i = 1; i <= msh->NbrVer; i++){
    for(int j = 0; j < 2 ; j++){
      msh->Ver[i].Crd[j] = dmsh->msh->Ver[i+4].Crd[j];
    }
  }

  int index_tri = 0;
  for(int i = 1; i <= dmsh->msh->NbrTri; i++){
    Triangle* tmp = (Triangle*)&dmsh->msh->Tri[i];
    int tri_bord = 0;
    for(int j = 0; j < 3; j++){
      tri_bord += set_contains(s,tmp->Ver[j]);
    }
    if(tri_bord == 0){
      index_tri ++;
      for(int j = 0; j < 3; j++){
        msh->Tri[index_tri].Ver[j] = tmp->Ver[j] - 4;
      }
    }
  }

  free(s);

  msh->NbrTri = index_tri;

  return msh;
}

static inline Mesh* simple_mesh_init(double xmin, double ymin, double xmax, double ymax){
  Mesh* msh = msh_init();
  if(! msh) return NULL;
  const double dX = 2.5*(xmax - xmin );
  const double dY = 2.5*(ymax - ymin );

  msh->Dim = 2;
  msh->NbrVer = 4;
  msh->NbrTri = 2;
  msh->Ver = calloc( (msh->NbrVer+1), sizeof(Vertex)   );
  msh->Tri = calloc( (msh->NbrTri+1), sizeof(Triangle) ); 
  msh->Ver[1].Crd[0] = xmin - dX;
  msh->Ver[1].Crd[1] = ymin - dY;

  msh->Ver[2].Crd[0] = xmax + dX;
  msh->Ver[2].Crd[1] = ymin - dY;

  msh->Ver[3].Crd[0] = xmax + dX;
  msh->Ver[3].Crd[1] = ymax + dY;

  msh->Ver[4].Crd[0] = xmin - dX;
  msh->Ver[4].Crd[1] = ymax + dY;

  msh->Tri[1].Ver[0] = 1;
  msh->Tri[1].Ver[1] = 2;
  msh->Tri[1].Ver[2] = 4;

  msh->Tri[2].Ver[0] = 2;
  msh->Tri[2].Ver[1] = 3;
  msh->Tri[2].Ver[2] = 4;
  
  msh_neighbors(msh);

  return msh;
}

int* get_corner_index(Mesh* msh){
  if(! msh){
    printf("msh null\n");
    return NULL;
  }
  int* res = (int*)calloc(4,sizeof(int));
  if(! res) return NULL;
  const float L_max = 100000;
  float lmax = -L_max;
  float lmin = L_max;
  float hmax = -L_max;
  float hmin = L_max;

  for(int i = 1; i <= msh->NbrVer; i++){
    Vertex* Si = (Vertex*)&msh->Ver[i];
    float l = Si->Crd[0];
    float h = Si->Crd[1];
    //printf("l h = %f %f\n",l,h);

    if(l>=lmax){
      lmax = l;
      if(h >= hmax){
        hmax = h;
        res[2] = i;
        
      }else if(h <= hmin){
        hmin = h;
        res[1] = i;
        
      }
      
    }
    if(l<=lmin){
      lmin = l;
      //printf("h = %f, hmax = %f\n",h,hmax);
      if(h >= hmax){
        hmax = h;
        res[3] = i;
        
      }else if(h <= hmin){
        hmin = h;
        res[0] = i;
        
      }
    }
  }

  //printf("get_corner_index:\nindex = %d %d %d %d\n",res[0],res[1],res[2],res[3]);

  return res;
}

Photo* photo_init(){
  Photo* res = (Photo*)malloc(sizeof(Photo));
  res->size = 0;
  return res;
};

Photo* photo_read(char* file_sol, char* file_msh, int readEfr){
  Photo* photo = photo_init();
  if(! photo) return NULL;
  Mesh* msh = msh_read(file_msh,readEfr);
  //printf("NbrTri = %d\n",msh->NbrTri);
  //printf("@ T_9034 = %p\n",&msh->Tri[9034]);
  //printf("S3 T_9034 = %d\n",msh->Tri[9034].Ver[2]);
  if(! msh){
    printf("msh read error\n");
    free(photo);
    return NULL;
  }
  int * corners = get_corner_index(msh);
  if(! corners){
    printf("corners null\n");
    free(photo);
    free(msh);
    return NULL;
  }

  photo->corners = corners;

  double* sol = sol_read(file_sol,msh->Dim,msh->NbrVer);
  if(! sol){
    printf("sol read error\n");
    free(photo);
    free(msh);
    return NULL;
  }
  photo->msh = msh;
  photo->size = msh->NbrVer;
  photo->data = sol;
  return photo;
}

static inline double* gradient_tri(double* vals, double2* S){
  double* res = (double*)malloc(sizeof(double)*2);
  if(! res) return NULL;
  double Ak = aire_signe(S);
  if(Ak == 0){
    printf("aire_K = 0\n");
    free(res);
    return NULL;
  }
  for(int i = 0; i < 3; i++){
    double dx = (S[(i+1)%3][1] - S[(i+2)%3][1])*0.5/Ak;
    double dy = (S[2 - (2*i)%3][0] - S[(1+i)%3][0])*0.5/Ak;
    res[0] += vals[i]*dx;
    res[1] += vals[i]*dy;
  }
  return res;
}

int save_photo(Photo* pht, char* file_name){
  if(!pht){
    printf("photo null\n");
    return 0;
  }
  if(! file_name){
    printf("file name error\n");
    return 0;
  }
  if(pht->size == 0){
    printf("size = 0\n");
    return 0;
  }
  if(! pht->msh){
    printf("msh null\n");
    return 0;
  }
  if(! pht->data){
    printf("data null\n");
    return 0;
  }

  char* name_msh = (char*)malloc(sizeof(char)*(strlen(".meshb") + strlen(file_name)));
  if(! name_msh) return 0;
  sprintf(name_msh,"%s.meshb",file_name);
  char* name_sol = (char*)malloc(sizeof(char)*(strlen(".solb") + strlen(file_name)));
  if(! name_sol) return 0;
  sprintf(name_sol,"%s.solb",file_name);


  int res0 = msh_write(pht->msh,name_msh);
  int res1 = msh_write2dfield_Vertices(name_sol, pht->msh->NbrVer, pht->data);

  return res0 * res1;
  
}

static inline double2* photo_gradient(Photo* pht){
  if(! pht) {
    printf("photo Null\n");
    return NULL;
  }
  if(pht->size == 0){
    return NULL;
  }
  if(! pht->msh){
    printf("msh Null\n");
    return NULL;
  }
  if(! pht->data){
    printf("data = Null\n");
    return NULL;
  }
  double2* grad_si = (double2*)calloc((pht->msh->NbrVer +1),sizeof(double2));
  if(! grad_si) return NULL;
  //printf("1\n");
  double2* grad_tri = (double2*)malloc(sizeof(double2)*(pht->msh->NbrTri + 1));
  //printf("2\n");
  if(!grad_tri) return NULL;
  double* aire_tri = (double*)malloc(sizeof(double)*(pht->msh->NbrTri + 1));
  //printf("3\n");
  if(! aire_tri) {
    free(grad_tri);
    return NULL;
  }
  double* aire_si = (double*)calloc((pht->msh->NbrVer + 1),sizeof(double));
  if(! aire_si){
    free(grad_tri);
    return NULL;
  }
  for(int i=1; i <= pht->msh->NbrTri; i ++){
    Triangle* T = (Triangle*)&pht->msh->Tri[i];
    double* vals = (double*)malloc(sizeof(double)*3);
    //printf("Tri = %d, NbrTri = %d\n",i,pht->msh->NbrTri);
    if(! vals){
      printf("vals malloc error\n");
      free(grad_tri);
      return NULL;
    }
    for(int j = 0; j < 3; j++){
      vals[j] = pht->data[T->Ver[j]];
    }
    double* grad_i = gradient_tri(vals,get_sommets_T(T,pht->msh));
    if(! grad_i){
      printf("grad_i error\n");
      free(vals);
      free(grad_tri);
      return NULL;
    }
    grad_tri[i][0] = grad_i[0];
    grad_tri[i][1] = grad_i[1];
    double aire_T = aire_signe(get_sommets_T(T,pht->msh));
    aire_tri[i] = aire_T;
    for(int j = 0; j < 3; j++){
      aire_si[T->Ver[j]] += fabs(aire_T);
    }
  }
  //printf("5\n");

  for(int i = 1; i <= pht->msh->NbrTri; i++){
    Triangle* T = (Triangle*)&pht->msh->Tri[i];
    //printf("6\n");
    for(int j = 0; j < 3; j++){
      int si = T->Ver[j];
      for(int z = 0; z < 2; z++){
        grad_si[si][z] += fabs(aire_tri[i]) * grad_tri[i][z] / aire_si[si];
      }
    }
  }
  //printf("7\n");

  return grad_si;
}

static inline double norm2(double x, double y){
  return sqrt(x*x + y*y);
}

cPhoto* cPhoto_init(){
  cPhoto* res = (cPhoto*)malloc(sizeof(cPhoto));
  if(! res) return NULL;
  res->size = 0;
  res->xmin = 0;
  res->xmax = 0;
  res->ymax = 0;
  res->ymax = 0;
  return res;
}

cPhoto* compress_1(Photo* pht, double level){
  if((level<=0)||(level>1)){
    printf("level error\n");
    return NULL;
  }
  if(! pht){
    printf("photo null\n");
    return NULL;
  }
  if(pht->size == 0){
    printf("photo size 0\n");
    return NULL;
  }
  cPhoto* res = cPhoto_init();
  if(! res){
    printf("cPhoto_init = NULL\n");
    return NULL;
  }

  res->xmin = pht->msh->Ver[pht->corners[0]].Crd[0];
  res->ymin = pht->msh->Ver[pht->corners[0]].Crd[1];
  res->xmax = pht->msh->Ver[pht->corners[2]].Crd[0];
  res->ymax = pht->msh->Ver[pht->corners[2]].Crd[1];

  //printf("xmin xmax : %f %f; ymin ymax : %f %f\n",res->xmin,res->xmax,res->ymin,res->ymax);

  int* mark_si = (int*)calloc(pht->size + 1, sizeof(int));
  if(! mark_si) return NULL;
  double2* grad_si = photo_gradient(pht);
  if(! grad_si){
    printf("grad_si = NULL\n");
    free(mark_si);
    free(res);
    return NULL;
  }
  double G_max = -1;
  double G_min = -1;
  double* Gs = (double*)malloc(sizeof(double)*(pht->size+1));
  if(!Gs) return NULL;
  for(int i = 1; i <= pht->size; i++){
    double* grad_i = grad_si[i];
    double G_i = norm2(grad_i[0],grad_i[1]);
    if(G_max<G_i) G_max = G_i;
    if(G_min<0) G_min = G_i;
    if(G_min>G_i) G_min = G_i;
    Gs[i] = G_i;
  }

  int Nbr_p = 0;
  //srand((unsigned)time(NULL));
  //srand(7);
  //const double db = (res->xmax - res->xmin) / (float)pht->size;
  //const double p = 0.25;
  for(int i = 1; i <= pht->size; i++){
    
    //double p = level * 2.0 * atan(norm2(grad_i[0],grad_i[1])) / PI;
    
    //double omega = (double)rand()/RAND_MAX;
    //int mark_i = (omega < p);
    double omega = 0.5;
    if(G_max > G_min) omega = (Gs[i] - G_min) / (G_max - G_min);
    //printf("omega = %f\n",omega);
    int mark_i = (omega > level);
    Nbr_p += mark_i;
    mark_si[i] = mark_i;
  }

  for(int q = 0; q<4; q++){
    int cq = pht->corners[q];
    if(mark_si[cq] == 0){
      mark_si[cq] = 1;
      Nbr_p ++;
    }
  }

  double2* pix = (double2*)malloc(sizeof(double2)*(Nbr_p));
  if(! pix) return NULL;
  double* data = (double*)calloc((Nbr_p+1),sizeof(double));
  if(! data) return NULL;

  res->size = Nbr_p;

  int index_p = 0;
  /*for(int q = 0; q<4; q++){
    int cq = pht->corners[q];
    printf("q = %d, cq = %f %f\n",cq, pht->msh->Ver[cq].Crd[0],pht->msh->Ver[cq].Crd[1]);
  }*/

  for(int i = 1; i <= pht->size; i++){
    /*int pq = 0;
    for(int q = 0; q<4; q++){
    int cq = pht->corners[q];
    if(i == cq){
      //printf("cq = %f %f\n",pht->msh->Ver[cq].Crd[0],pht->msh->Ver[cq].Crd[1]);
      pq = 1;
      break;
    }
  }*/
    if(mark_si[i]){
      pix[index_p][0] = pht->msh->Ver[i].Crd[0];
      pix[index_p][1] = pht->msh->Ver[i].Crd[1];
      data[index_p+1] = pht->data[i];
      //if(pq == 1) printf("i = %d, caq = %f %f\n",i,pht->msh->Ver[i].Crd[0],pht->msh->Ver[i].Crd[1]);
      index_p ++;
    }
  }

  res->data = data;
  res->pix = pix;

  return res;
}

Photo* decode_cPhoto(cPhoto* cpt){
  if(! cpt){
    printf("cpt null\n");
    return NULL;
  }
  if(cpt->size == 0){
    printf("size = 0\n");
    return NULL;
  }
  if(! cpt->data){
    printf("data null\n");
    return NULL;
  }
  if(! cpt->pix){
    printf("pix null\n");
    return NULL;
  }

  Photo* res = photo_init();
  if(!res) return NULL;

  Mesh* msh = simple_mesh_init(cpt->xmin,cpt->ymin,cpt->xmax,cpt->ymax);
  if(! msh){
    printf("msh null\n");
    return NULL;
  }

  DMesh* dmsh = dMesh_init(msh,1);

  if(! dmsh){
    printf("dmsh null\n");
    return NULL;
  }

  Triangulation(cpt->pix,cpt->size,dmsh);

  //msh_write(dmsh->msh,"res_test0.meshb");

  Mesh* msh_deco = toMesh(dmsh);
  if(! msh_deco){
    printf("msh_deco error\n");
    return NULL;
  }

  res->data = cpt->data;
  res->msh = msh_deco;

  res->size = cpt->size;
   
  return res;
}