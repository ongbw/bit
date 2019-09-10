/*---------------------------------------------------------------------------*/
//  tree3d_v2.h                                                            //
//  original 2D Code by Jerry Emhoff                                         //
//  Modified by Andrew Christlieb (6/4/07)                                   //
//  modified by bwo (16/07/08)                                               //
//  modified by bwo (02/16/10)                                               //
//    -- using pointers for moments
//  modified by ongbw (09/10/19) - two parameter family of regularization //
/*---------------------------------------------------------------------------*/
//
// Global Variables
// bdy_flg: flag for free space problem or bounded problem
//             1-free space (to be coded first)

#ifndef tree3D_h
#define tree3D_h 1

#define bound_flg 0
#define panelreg 1e-1

#include <fstream>             // File stream library (std. C)
#include <iostream>            // File stream library (C++)
#include <iomanip>
#include <cmath>              //  Math library
#include <ctime>
#include <cstring>
#include <cstdlib>
#include "constants.h"
#include <sstream>
using namespace std;

// Structure for particle quantities
struct PARTICLE{
  double charge;  // Charge on the species (i.e. +1, -1, +2, etc.)
  double x;       // x-position of the particle
  double y;       // y-position
  double z;       // z-position
  double u;       // Velocity in the x-direction
  double v;       // Velocity in the y-direction
  double w;       // Velocity in the z-direction
  double xforce;  // Force in the x-direction
  double yforce;  // Force in the y-direction
  double zforce;  // Force in the z-direction
  // Particle weight - i.e. the amount of real particles represented
  double weight;
  double tot_charge;  // Charge x Weight
};


// Structure for tree properties
struct TREE {
  int nummemb;      // Number of member particles
  double x0;        // Left x-limit
  double x1;        // Right x-limit
  double y0;        // Left y-limit
  double y1;        // Right y-limit
  double z0;        // Bottom z-limit
  double z1;        // Top z-limit
  double xmid;      // x midpoint
  double ymid;      // y midpoint
  double zmid;      // z midpoint
  double aspectxy;  // xy aspect ratio
  double aspectxz;  // xz aspect ratio
  double aspectyz;  // yz aspect ratio
  double sqradius;  // square of the radius from center of mass to origin
  int children[NUMCHILD]; // Child trees of the current tree
  int *members;     // Index of member particles
  double ***moments;  // Moments of the cluster

  TREE ()
  : nummemb(0),
    members(0x0),
    moments(0x0) { }
};

// Structure for panel
struct PANEL{
  double x0;      // Starting x-point
  double x1;      // Ending x-point
  double y0;      // Starting y-point
  double y1;      // Ending y-point
  double z0;      // Starting z-point
  double z1;      // Ending z-point
  double str;     // Panel strength   : for boundary conditions, sigma
  double xnrm;    // X-normal direction
  double ynrm;    // Y-normal direction
  double znrm;    // Z-normal direction
  double refpot;  // Reference potential or slope ??? what is this?
  double midx;    // x midpoint
  double midy;    // y midpoint
  double midz;    // z midpoint
  int type;       // Type of panel - Dirichlet or Neumann
  int orientation; // 1: xy-panels, 2: xz panels, : yzpanels
  double area;    // area of panel
};

struct BOUNDARY{
  int type;
  double length;
  double width;
  double x0;
  double y0;
  double z0;
  double x_normal;
  double y_normal;
  double z_normal;
  double refpot;
  int direction;
};

////////////////////////////
// Treecode functions

// Creates the tree
int treemake(int num, int nterms, int maxmemb, PARTICLE* part, TREE* tree);
// Deletes the tree particle arrays
void treeempty(int numtree, TREE* tree);
// Initialize the tree
void treeinit(TREE* tree1, double x0, double x1, double y0,
	      double y1, double z0, double z1);
// Divide a cluster into two pieces
int treedivide(int dir, TREE* split, TREE* c1, TREE* c2, PARTICLE* part);
// Shrink a cluster around its members
void treeshrink(TREE* tree1, PARTICLE* part);
// Calculates the moments of all tree clusters
void treemoments(int numtree, int nterms, TREE* tree, PARTICLE* part);


/*****************************************************************************/

// create random double between 0 and 1
double RandL(void);


////////////////////////////
// particle initializations

// initialize mesh
void init_mesh(int Nx, int Ny, int Nz,
	       double dlength, double dheight, double ddepth,
	       PARTICLE* part);

// Initializes particles in the box domain
void boxparticles(int npart, double dlength, double dheight,
                  double ddepth, PARTICLE* part, double density);
int part_vel(PARTICLE* part,int hot_or_cold,int s, int cnt);

//find bounding box
int bounding_box(int npart, PARTICLE* part,double &Mx, double &mx, double &My,
		 double &my, double &Mz, double &mz);


////////////////////////////
// boundary initialization
int setboundaries(int xpanels, int ypanels, int zpanels, double dlength,
                  double dheight, double ddepth, PANEL* pan);



/////////////////////////////////////
// Force / electric field / potential calculation

// evaluate gradient of unregularized kernel
void eval_grad_kernel_noreg(double dx, double dy, double dz, double rs,
			    double *xfor, double *yfor, double *zfor);

// evaluate gradient of regularized kernel
void eval_grad_kernel_reg(double dx, double dy, double dz, double rs,
			  double *xfor, double *yfor, double *zfor);

// Computes the force at a point (x,y,z) due to particles using direct sum, no regularization
void force_point_ds_noreg(double x, double y, double z, int npart, PARTICLE* part,
			  double *xfor, double *yfor, double *zfor);

// Computes the force on particle part[ind] due to particles using direct sum, no regularization
void force_particle_ds_noreg(int ind, int npart, PARTICLE* part,
			     double *xfor, double *yfor, double *zfor);

// Computes the force at a point (x,y,z) due to particles using direct sum, regularized kernel
void force_point_ds_reg(double x, double y, double z, int npart, PARTICLE* part,
			double *xfor, double *yfor, double *zfor);

// Computes the force on a particle, part[ind], due to particles using direct sum, regularized kernel
void force_particle_ds_reg(int ind, int npart, PARTICLE* part,
			   double *xfor, double *yfor, double *zfor);

// calculates the electrostatic force on all particles, storing force in
// part[j].xforce, part[j].yforce, part[j].zforce using direct sum
void freespace_particles_ds_noreg(int npart, PARTICLE* part);

// calculates the electrostatic force on all particles, storing force in
// part[j].xforce, part[j].yforce, part[j].zforce using direct sum
void freespace_particles_ds_reg(int npart, PARTICLE* part);


// Compute the force at a point due to a cluster of particles
void treetaylorforce(int nterms, double dx, double dy, double dz, double rs,
                     double *xfor, double *yfor, double* zfor,
		     double*** moment);

// Main algorithm for computing the force at a point due to all particles
// i.e. decision on when to use direct sum or treecode
double treeforce(int ind, int p, int nterms, double x, double y, double z,
		 double acc, double *xfor, double *yfor, double *zfor,
		 TREE* tree, PARTICLE* part);





// Computes shielded force on a particle due to other particles using direct sum
void shielded_force_particle_ds(int ind, int npart, PARTICLE* part,
				double *xfor, double*yfor, double*zfor,
				double alpha);

// Computes shielded pot on a particle due to other particles using direct sum
void shielded_pot_particle_ds(int ind, int npart, PARTICLE* part,
			      double *pot, double alpha);


// Computes the force at a point due to all particles using direct sum
void shielded_force_point_ds(double x, double y, double z,
			     int npart, PARTICLE* part,
			     double* xfor, double* yfor, double* zfor,
			     double alpha);



// Computes the potential at a point due to all particles using direct sum
double potreal_point(double x, double y, double z, int npart,
		   PARTICLE* part);




// Compute the shielded force at a point due to a cluster of particles
void shielded_treetaylorforce(int nterms, double dx, double dy, double dz,
			      double rs, double *xfor, double *yfor,
			      double* zfor, double*** moment, double kappa);
void shielded_treetaylorpot(int nterms, double dx, double dy, double dz,
			    double rs, double *pot,
			    double*** moment, double kappa);

/*
void streetaylorforce(int nterms, double dx, double dy, double dz,
		      double rs, double *xfor, double *yfor,
		      double* zfor, double*** moment, double kappa);
*/



// Main algorithm for computing the shielded force at a point due to
// all particles, i.e. decision on when to use direct sum or treecode
double shielded_treeforce(int ind, int p, int nterms,
			  double x, double y, double z, double alpha,
			  double acc, double *xfor, double *yfor,
			  double *zfor, TREE* tree, PARTICLE* part);

double shielded_treepot(int ind, int p, int nterms,
			double x, double y, double z, double alpha,
			double acc, double *pot, TREE* tree, PARTICLE* part,
			int Nx, int Ny, int Nz,
			double dlength, double dheight, double ddepth);



// Compute the potential at a point due to a cluster of particles
void treetaylorpotential(int nterms, double dx, double dy, double dz, double rs,
			 double* pot, double*** moment);

// Main algorithm for computing the  force at a point due to all particles
// i.e. decision on when to use direct sum or treecode
double treepotential(int ind, int p, int nterms, double x, double y, double z,
		     double acc, double* pot,
		     TREE* tree, PARTICLE* part);




// calculates the shielded electrostatic force on each particle,
// and assigns it to  part[j].xforce, part[j].yforce, part[j].zforce
// using direct sum
void shielded_freespace_particles_ds(int npart, PARTICLE* part, double alpha);
void shielded_freespace_pot_ds(int npart, PARTICLE* part, double alpha);


// calculates the electrostatic force on each particle, and assigns it to
// part[j].xforce, part[j].yforce, part[j].zforce using the treecode
void freespace_particles_tree(int npart, PARTICLE* part,int nterms,
			      double acc, TREE* tree, int numtree, int maxmemb);

// calculates the shielded electrostatic force on each particle,
// and assigns it to  part[j].xforce, part[j].yforce, part[j].zforce
// using the treecode
void shielded_freespace_particles_tree(int npart, PARTICLE* part, double alpha,
				       int nterms, double acc, TREE* tree,
				       int numtree, int maxmemb);




int box_calc_ds(int npart, PARTICLE* part, int numpan, PANEL* pan,
		double** panelarray, PARTICLE* bpart, double dlength,
		double dheight, double ddepth, double density);

int box_calc_tree(int npart, PARTICLE* part, int numpan, PANEL* pan,
		  double** panelarray, PARTICLE* bpart, double dlength,
		  double dheight, double ddepth, double density,
		  int nterms, double acc, TREE* tree, int numtree, int maxmemb);

// time integration

int vel_update(int npart, PARTICLE* part, double dt);
int loc_update(int npart, PARTICLE* part, double dt);



////////////////////////////
// Panel functions

// Panel initialization function
void panelinit(PANEL* panl, double x0, double x1, double y0, double y1,
	       double z0, double z1, double refpot, double xnrm,
	       double ynrm, double znrm, int type, int orientation);

// Fills the panel effects array using naive riemann integration
int panelconstantsreal(int numpan, PANEL* pan, double** panelarray);

// Fills the panel effects array using riemann integration: ill conditioned!
int badpanelconstantsreal(int numpan, PANEL* pan, double** panelarray);

// Solves for the panel strengths using direct summation
int panel_solve_ds(int numpan, int npart, PARTICLE* part, PANEL* pan,
		   double** panelarray);

// Solve for the panel strengths using the treecode.
int panel_solve_tree(int numpan, int npart, int nterms, double acc,
		     double** panelarray, PARTICLE* part, PANEL* pan,
		     TREE* tree);

// Create effective particles from the panel strengths
int effective_charge(int numpan, PANEL* pan, PARTICLE* bpart,
		     double dlength, double dheight, double ddepth,
		     double density, int npart);


/////////////////////////////////////
// Output particle properties
void output_particles(int npart, PARTICLE* part, string filename);

void output_potential(int npart, PARTICLE* part, double *pot, string filename);

void output_slice(int n, PARTICLE* part, double *pot, string headername);

void output_fields(int npart, PARTICLE* part,
		   double *Bx, double *By, double* Bz, string filename);

void plot_efield_ds(int nxnodes, int nynodes, int nznodes,
		    double dlength, double dheight, double ddepth,
		    int npart, PARTICLE* part);

void plot_efield_tree(int nxnodes, int nynodes, int nznodes, int npart,
		      int nterms, int maxmemb, double acc,
		      double dlength, double dheight, double ddepth,
		      PARTICLE* part, TREE* tree);


// Global variables
//double panelreg;  // regularization term for panels
//int bound_flg;   // type of boundary conditions
//BOUNDARY* boundaries; //array of boundary vertices and sides
//int npoints;  // number of integration points to use

// math functions
double sq(double x);

#endif //!defined(tree3D_h)
