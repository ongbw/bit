/*---------------------------------------------------------------------------*/
// treecode_mp.cpp
// checking output of treecode-mp to treecode
/*---------------------------------------------------------------------------*/

#include "tree3d.h" // tree code structures
#include "problem.h"
#include "constants.h"
#include <omp.h>
#include <time.h>
#include <limits>


int main()  {

  int npart = 16e3; // number of particles

  PARTICLE* part;  // initialize Particle data type
  part = new PARTICLE[npart]; 
 
  //Tree Related Data Type
  int numtree=0;    // Number of trees
  int maxmemb=100;  // Maximum number of particles per bottom level cluster
  int maxtrees;     // Maximum number of clusters
  maxtrees = (int) 3*npart/maxmemb + 10;
  TREE* tree;       // initialize tree data type
  tree = new TREE[maxtrees];

  // Domain setup
  double dlength = 1.0;     // Set the domain length
  double dheight = 1.0;     // Set the domain height
  double ddepth = 1.0;      // Set the domain width
  double density = 1.0;    // plasma density

  // Initialize all tree pointers to NULL
  for (int i = 0; i < maxtrees; i++)  {
    tree[i].members = NULL;
    tree[i].moments = NULL;
  }

 // Initialize the particles in the domain
  boxparticles(npart, dlength, dheight, ddepth, part, density);

  string filename;
  
  // variables for computing force
  double xfor;
  double yfor;
  double zfor;

  //store the treecode solution to the variables fx, fy, fz
  double* fx = new double[npart]; // store treecode solution
  double* fy = new double[npart]; // store treecode solution
  double* fz = new double[npart]; // store treecode solution
  

  int nterms=5;
  // create tree
  numtree = treemake(npart,nterms,maxmemb,part,tree);
  treemoments(numtree,nterms,tree,part);


  double clocktime;

  clocktime = omp_get_wtime();
  // directsum
  for (int iii=0; iii < npart; iii++) {
    xfor = 0.0; yfor = 0.0; zfor = 0.0;
    treeforce(0,iii,nterms,part[iii].x,part[iii].y,part[iii].z,
	      ACC,&xfor, &yfor, &zfor,tree,part);
    fx[iii]=xfor*charge_constant*part[iii].tot_charge;
    fy[iii]=yfor*charge_constant*part[iii].tot_charge;
    fz[iii]=zfor*charge_constant*part[iii].tot_charge;
  }
  cout << "single core = " << omp_get_wtime()-clocktime << "s" << endl; 


  // treecode_mp
  double* fx2 = new double[npart]; // store treecode-mp solution
  double* fy2 = new double[npart]; // store treecode-mpsolution
  double* fz2 = new double[npart]; // store treecode-mp solution
  
#pragma omp parallel 
  {
    clocktime = omp_get_wtime();
#pragma omp for private(xfor,yfor,zfor)
    for (int iii=0; iii<npart; iii++)  {
      xfor = 0.0; yfor = 0.0; zfor = 0.0;
      treeforce(0,iii,nterms,part[iii].x,part[iii].y,part[iii].z,
		ACC,&xfor, &yfor, &zfor,tree,part);
      fx2[iii]=xfor*charge_constant*part[iii].tot_charge;
      fy2[iii]=yfor*charge_constant*part[iii].tot_charge;
      fz2[iii]=zfor*charge_constant*part[iii].tot_charge;
    }
  }

 cout << "multi-core timing (" << omp_get_max_threads() <<" thread(s)) = " 
       << omp_get_wtime()-clocktime << "s" << endl; 

  // check and output error
  double errx;
  double erry;
  double errz;

  errx = 0.0;  erry = 0.0; errz = 0.0;
  for (int i = 0; i<npart; i++) {
    if (abs(fx2[i]-fx[i])>errx) 
      errx = abs(fx2[i]-fx[i]);
    if (abs(fy2[i]-fy[i])>erry) 
      erry = abs(fy2[i]-fy[i]);
    if (abs(fz2[i]-fz[i])>errz) 
      errz = abs(fz2[i]-fz[i]);
  }

  double err_tot = abs(errx) + abs(erry) + abs(errz);

  if (abs(err_tot) < 1.0e-15){
    cout << "Test PASSED" << endl;
    return 0;
  }  else {
    cout << "Test FAILED" << endl;
    cout << "ERROR: treecode and treecode-mp not within machine precision";
    return 1;
  }

  if (part != NULL)
    delete [] part;

  delete [] fx;
  delete [] fy;
  delete [] fz;
  delete [] fx2;
  delete [] fy2;
  delete [] fz2;

} // end main


