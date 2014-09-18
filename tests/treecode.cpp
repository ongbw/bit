/*---------------------------------------------------------------------------*/
// treecode.cpp
// bwo, May 2010
// convergence of treecode-mp to directsum-mp
/*---------------------------------------------------------------------------*/

#include "tree3d.h" // tree code structures
#include "problem.h"
#include "constants.h"
#include <omp.h>
#include <time.h>

int main()  {

  int npart = 16000; // number of particles

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


  double* fx = new double[npart]; // store directsum solution
  double* fy = new double[npart]; // store directsum solution
  double* fz = new double[npart]; // store directsum solution
  
  
  // directsum_mp
#pragma omp parallel 
  {
#pragma omp for private(xfor,yfor,zfor)
    for (int iii=0; iii<npart; iii++)  {
      xfor = 0.0; yfor = 0.0; zfor = 0.0;
      force_particle_ds(iii,npart,part,&xfor, &yfor, &zfor);
      fx[iii]=xfor*charge_constant*part[iii].tot_charge;
      fy[iii]=yfor*charge_constant*part[iii].tot_charge;
      fz[iii]=zfor*charge_constant*part[iii].tot_charge;
    }
  }
  

  // treecode_mp
  double* fx2 = new double[npart]; // store treecode solution
  double* fy2 = new double[npart]; // store treecode solution
  double* fz2 = new double[npart]; // store treecode solution

  double errx;
  double erry;
  double errz;

  double errstore;
  double errnew;

  // check convergence in the number of terms
  int failflag=0;
  for (int nterms = 1; nterms<5; nterms++) {
    // create tree
    tree = new TREE[maxtrees];
    numtree = treemake(npart,nterms,maxmemb,part,tree);
    treemoments(numtree,nterms,tree,part);

    // calculate new force
#pragma omp parallel 
    {
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

    // check and output error
    errx = 0.0;  erry = 0.0; errz = 0.0;
    for (int i = 0; i<npart; i++) {
      if (abs(fx2[i]-fx[i])>errx) 
	errx = abs(fx2[i]-fx[i]);
      if (abs(fy2[i]-fy[i])>erry) 
	erry = abs(fy2[i]-fy[i]);
      if (abs(fz2[i]-fz[i])>errz) 
	errz = abs(fz2[i]-fz[i]);
    }
    cout << "nterms = " << nterms 
	 << ", errx = " << errx
	 << ", erry = " << erry
	 << ", errz = " << errz
	 << endl;

    errnew = abs(errx) + abs(erry) + abs(errz);
    if (nterms == 1)
      errstore=errnew;
    else
      if (errnew > errstore)
	failflag=1;
      else
	errstore = errnew;


    // normally, one doesn't create multiple trees. this was just
    // easier than reinitializing structure components which were
    // changing with nterms
    if (tree != NULL)
      delete [] tree;

  } // end nterms


  if (part != NULL)
    delete [] part;
  
  
  delete [] fx;
  delete [] fy;
  delete [] fz;
  delete [] fx2;
  delete [] fy2;
  delete [] fz2;

  if (failflag == 1){
    cout << "ERROR: treecode not converging to directsum";
    cout << "Test FAILED" << endl;
  }  else {
    cout << "Test PASSED" << endl;
  }

  return failflag;

} // end main

