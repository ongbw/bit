/*---------------------------------------------------------------------------*/
// kernel_convergence.cpp
// checking convergence of regularized kernel to non-regularized kernel
/*---------------------------------------------------------------------------*/

#include "tree3d.h" // tree code structures
#include "problem.h"
#include "constants.h"
#include <time.h>
#include <limits>

int main()  {


  double dx = 0.1;
  double dy = 0.2;
  double dz = 0.3;
  
  //store the unregularized solution to the variables fx, fy, fz
  double fx, fy, fz;
  double rs = dx*dx + dy*dy + dz*dz;
  printf("computing gradient of unregularized kernel\n");

  eval_grad_kernel_noreg(dx,dy,dz,rs,&fx,&fy,&fz);
  
  // computing regularized solution
  double fx2, fy2, fz2;
  rs = rs + DEL;
  

  //for (int n = 0; n < 20; n++) {
  printf("testing kernel order %d\n",KERNEL_ORDER);
  eval_grad_kernel_reg(dx,dy,dz,rs,&fx2,&fy2,&fz2);

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

  // for KERNEL_ORDER = 0, DEL = 0.01,
  // errx = 1.876829391553e-01, erry = 3.753658783106e-01,  errz = 5.630488174659e-01 
  
  double err_tot = abs(errx) + abs(erry) + abs(errz);

  if (abs(err_tot) < 1.2){
    cout << "Test PASSED" << endl;
    return 0;
  }  else {
    cout << "Test FAILED" << endl;
    cout << "ERROR: kernel regularzation not within tolerance";
    return 1;
  }  
} // end main


