#include <mesh.h>
#include <stdio.h>
#include <time.h>

int TP1(int argc, char *argv[]){
  double to,ti;
  
  if ( argc < 2 ) {
    printf(" usage : mesh file \n");
    return 0;
  }
  
  /* read a mesh */
  to =  GetWallClock();
  int readEfr = 1;
  Mesh * msh = msh_read(argv[1], readEfr);
  printf("NbrEfr = %d\n",msh->NbrEfr);
  ti =  GetWallClock();
  
  if ( ! msh ) return 0;
  
  printf("  Vertices   %10d \n", msh->NbrVer);
  printf("  Triangles  %10d \n", msh->NbrTri);
  printf("  time to read the mesh %lg (s) \n",ti-to);
  
  /* re-order a mesh */ 
  to =  GetWallClock();
  msh_reorder(msh);
  ti =  GetWallClock();  
  printf("  time to re-order the mesh  %lg (s) \n",ti-to);
   
  /* create neigbhors Q2 version */
  to =  GetWallClock();
  //msh_neighborsQ2(msh);  %Pour mesurer correctement le temps d'excutation de msh_neighbors avec Hsh, il faut commenter cette ligne
  ti =  GetWallClock();
  printf("  time q2 neigh.        %lg (s) \n",ti-to);
  
  /* create neigbhors with hash table */
  to =  GetWallClock();
  msh_neighbors(msh);
  ti =  GetWallClock();
  printf("  time hash tab neigh.  %lg (s) \n",ti-to);
  
  /* write reordered mesh */
  to =  GetWallClock();
  msh_write(msh,"output.meshb");
  ti =  GetWallClock();
  
  int      iTri, iVer;
  double  *Qua = (double  *)malloc(sizeof(double ) * (msh->NbrTri+1));
  double  *Qua2 = (double *)malloc(sizeof(double ) * (msh->NbrTri+1));
  double3 *Met = (double3 *)malloc(sizeof(double3) * (msh->NbrVer+1));
  

  printf("Q1 = %f, Q2 = %f\n",Q1(msh, 1),Q2(msh, 1));
  
  for (iTri=1; iTri<=msh->NbrTri; iTri++) {
    Qua[iTri]  = Q1(msh, iTri);

    Qua2[iTri] = Q2(msh, iTri);
  } 
  
  for (iVer=1; iVer<=msh->NbrVer; iVer++) {
  	 Met[iVer][0] = 1.;
  	 Met[iVer][1] = 0.;
  	 Met[iVer][2] = 1.;
  } 

  double * res = colorisation(msh); //colorisation d'un maillage
  
  msh_write2dfield_Triangles("quality.solb", msh->NbrTri, Qua);
  msh_write2dfield_Triangles("quality2.solb", msh->NbrTri, Qua2);
  msh_write2dfield_Triangles("color.solb", msh->NbrTri, res);
 
  msh_write2dmetric("metric.solb", msh->NbrVer, Met);

  return 0;
}

void TP2_1(){
  Mesh* msh = msh_init();
  const double xmin = 0;
  const double ymin = 0;
  const double xmax = 5;
  const double ymax = 5;

  msh->Dim = 2;
  msh->NbrVer = 4;
  msh->NbrTri = 2;
  msh->Ver = calloc( (msh->NbrVer+1), sizeof(Vertex)   );
  msh->Tri = calloc( (msh->NbrTri+1), sizeof(Triangle) ); 
  msh->Ver[1].Crd[0] = xmin - 1;
  msh->Ver[1].Crd[1] = ymin - 1;

  msh->Ver[2].Crd[0] = xmax + 1;
  msh->Ver[2].Crd[1] = ymin - 1;

  msh->Ver[3].Crd[0] = xmax + 1;
  msh->Ver[3].Crd[1] = ymax + 1;

  msh->Ver[4].Crd[0] = xmin - 1;
  msh->Ver[4].Crd[1] = ymax + 1;

  msh->Tri[1].Ver[0] = 1;
  msh->Tri[1].Ver[1] = 2;
  msh->Tri[1].Ver[2] = 4;

  msh->Tri[2].Ver[0] = 2;
  msh->Tri[2].Ver[1] = 3;
  msh->Tri[2].Ver[2] = 4;


  const int Nx = 128;
  const int Ny = 64;
  
  double2* Ps = points_unif(xmin,ymin,xmax,ymax,Nx,Ny);

  DMesh* dmsh = dMesh_init(msh,1);
  //dMesh_resize(dmsh, Nx * Ny + 5,Nx * Ny * 3);
  printf("%d\n",Triangulation(Ps,Nx * Ny,dmsh));

  Mesh* res = toMesh(dmsh);

  msh_write(dmsh->msh,"res0.meshb");
  msh_write(res,"res.meshb");
}

void TP2_2(){

  //Photo* photo = photo_read("data/joconde.lowres.sol","data/joconde.lowres.mesh",1);
  Photo* photo = photo_read("data/joconde.sol","data/joconde.mesh",1);
  printf("size_init = %d\n",photo->size);
  double level = 0.12;
  //level = 0.2;
  //level = 0.175;
  //level = 0.15;
  level = 0.125;
  //level = 0.1;
  printf("level = %f\n",level);
  cPhoto* cpt = compress_1(photo,level);
  printf("size_compressed = %d ~= %f %%\n",cpt->size,(double)cpt->size * 100.0 / (double)photo->size);
  Photo* photo_deco = decode_cPhoto(cpt);
 
  char* s = "res_test";
  printf("%d\n",save_photo(photo_deco,s));
}

int main(int argc, char *argv[])
{ 
  //TP1(argc,argv);
  //TP2_1();
  TP2_2();
   
  return 0;
}

