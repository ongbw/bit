/*---------------------------------------------------------------------------*/
//  constants.h                                                              //
/*---------------------------------------------------------------------------*/
 
#ifndef constants_h
#define constants_h

const int MAXPART =1000000;       // Maximum number of particles
const int NUMCHILD=8;           // Number of children per tree
const int MAXTERM=12;           // Maximum number of Taylor expansion terms
const int MAXCHILD=100000;       // Maximum number of trees
const double pi=3.14159265359;
const double inv2pi=0.5/pi;


// Physical constants
const double perm=8.85418781762e-12;
const double elec=1.60217733e-19;
const double invmass=1.0/9.1093897e-31;
const double mu=4e-7*pi;
const double c=3.0e8;
const double charge_constant=inv2pi/perm;


/*
// normalized constants
const double perm=1.0;
const double elec=1.0;
const double invmass=1.0;  
const double charge_constant = 1.0;
const double mu=1.0;
const double c = 0.1;
*/
                            
const double eps=1e-14;  // small #check distance from particles.
const double Bfield=10.0;

const double endcap = -50.0;
// Seed the ranf() random number generator
static unsigned int seed = 0x11111111;  

#endif //!defined(constants.h)
