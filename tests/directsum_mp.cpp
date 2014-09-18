/*---------------------------------------------------------------------------*/
// directsum_mp.cpp
// checking output of directsum-mp to directsum
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
 
  // Domain setup
  double dlength = 1.0;     // Set the domain length
  double dheight = 1.0;     // Set the domain height
  double ddepth = 1.0;      // Set the domain width
  double density = 1.0;    // plasma density

 // Initialize the particles in the domain
  boxparticles(npart, dlength, dheight, ddepth, part, density);

  string filename;
  
  //local variables for computing force
  double xfor;
  double yfor;
  double zfor;

  double clocktime;


  //store the directsum solution to the variables fx, fy, fz
  double* fx = new double[npart]; // store directsum solution
  double* fy = new double[npart]; // store directsum solution
  double* fz = new double[npart]; // store directsum solution
  

  clocktime = omp_get_wtime();
  // directsum
  for (int iii=0; iii < npart; iii++) {
    xfor = 0.0; yfor = 0.0; zfor = 0.0;
    force_particle_ds(iii,npart,part,&xfor, &yfor, &zfor);
    fx[iii]=xfor*charge_constant*part[iii].tot_charge;
    fy[iii]=yfor*charge_constant*part[iii].tot_charge;
    fz[iii]=zfor*charge_constant*part[iii].tot_charge;
  }
  cout << "single-core timing = " << omp_get_wtime()-clocktime << "s" << endl; 


  // directsum_mp
  // treecode_mp
  double* fx2 = new double[npart]; // store directsum-mp solution
  double* fy2 = new double[npart]; // store directsum-mpsolution
  double* fz2 = new double[npart]; // store directsum-mp solution
  
#pragma omp parallel 
  {
    clocktime = omp_get_wtime();
#pragma omp for private(xfor,yfor,zfor)
    for (int iii=0; iii<npart; iii++)  {
      xfor = 0.0; yfor = 0.0; zfor = 0.0;
      force_particle_ds(iii,npart,part,&xfor, &yfor, &zfor);
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
    cout << "ERROR: directsum and directsum-mp not within machine precision";
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


