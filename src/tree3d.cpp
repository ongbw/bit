/*---------------------------------------------------------------------------*/
//  tree3d.cpp                                                               //
//  contains function files for creating, deleting and modifying tree struct //
//  adding helmholtz recurrence relations
//  bwo (16/07/08)                                                           //
/*---------------------------------------------------------------------------*/

#include "tree3d.h"  // class files for the tree structure
#include "problem.h"
#define verbosity 0 // switch flag to 1 if you want output
 

// return a random double between 0 and 1
double RandL(void) {
  return (rand()/((double)RAND_MAX));
}

///////////////////////////////////////////////////////
// Function Name: treeempty
// Usage: deletes entire allocated tree
//
///////////////////////////////////////////////////////
// Inputs:
//         numtree - total number of tree nodes
//         tree    - tree to be deleted
//
///////////////////////////////////////////////////////
// Outputs: none
//
///////////////////////////////////////////////////////
// Functions Called: none   
//
///////////////////////////////////////////////////////
void treeempty(int numtree, TREE* tree)  {
  if (tree != NULL)  {
    for (int i=0; i<numtree; i++)  {
      if (tree[i].members != NULL)
        delete [] tree[i].members;
      tree[i].members=NULL;
      
      if (tree[i].moments != NULL)
        delete [] tree[i].moments;
      tree[i].moments=NULL;
    }
  }
  return;
}



///////////////////////////////////////////////////////
// Function Name: treeinit
// Usage: Initialize a tree cluster
//
///////////////////////////////////////////////////////
// Assumptions: global parameter NUMCHILD = 8 (oct tree)
//
///////////////////////////////////////////////////////
// Inputs:
//         tree1  - tree node we are initalizing
//         x0     - node min 
//         x1     - node max
//         y0     - node min 
//         y1     - node max
//         z0     - node min
//         z1     - node max
//
///////////////////////////////////////////////////////
// Outputs: (By Reference)
//        tree1->x0   -   node min
//        tree1->x1   -   node max
//        tree1->y0   -   node min
//        tree1->y1   -   node max 
//        tree1->z0   -   node min
//        tree1->z1   -   node max
//        tree1->xmid -   node midpoint x
//        tree1->ymid -   node midpoint y
//        tree1->zmid -   node midpoint z
//        tree1->aspectxy - aspect ratio of node
//        tree1->aspectxz - aspect ratio of node
//        tree1->aspectyz - aspect ratio of node
//        tree1->sqradius - squar radius of node
//
///////////////////////////////////////////////////////
// Functions Called:  none
//
///////////////////////////////////////////////////////
void treeinit(TREE* tree1, double x0, double x1, double y0, 
	      double y1, double z0, double z1)  {

  double x=0.5*(x1+x0);
  double y=0.5*(y1+y0);
  double z=0.5*(z1+z0);

  tree1->x0=x0;
  tree1->x1=x1;
  tree1->y0=y0;
  tree1->y1=y1;
  tree1->z0=z0;
  tree1->z1=z1;

  tree1->xmid=x;
  tree1->ymid=y;
  tree1->zmid=z;

  tree1->aspectxy=fabs((x1-x0)/(y1-y0));
  tree1->aspectxz=fabs((x1-x0)/(z1-z0));
  tree1->aspectyz=fabs((y1-y0)/(z1-z0));

  x=x1-x;
  y=y1-y;
  z=z1-z;
  tree1->sqradius=x*x+y*y+z*z;
  for (int i=0; i<NUMCHILD; i++)
    tree1->children[i]=-1;
}


///////////////////////////////////////////////////////
// Function Name:  treeshrink 
// Usage: shrink tree nodes to just bound particles
//
///////////////////////////////////////////////////////
// Inputs:
//         tree1 - node of oct tree
//         part  - particle list
//
///////////////////////////////////////////////////////
// Outputs: (By Reference)
//         tree1->x0
//         tree1->x1
//         tree1->y0
//         tree1->y1
//         tree1->z0
//         tree1->y1
//
///////////////////////////////////////////////////////
// Functions Called:
//         treeinit(tree1, minx, maxx, miny, maxy, minz, maxz)
//
///////////////////////////////////////////////////////
void treeshrink(TREE* tree1, PARTICLE* part)  {
  int j=tree1->members[0];
  double maxx=part[j].x;
  double maxy=part[j].y;
  double maxz=part[j].z;
  double minx=part[j].x;
  double miny=part[j].y;
  double minz=part[j].z;

  //loop over particles in tree node
  for (int i=0; i<tree1->nummemb; i++)  {
    j=tree1->members[i];
    if (part[j].x > maxx)
      maxx=part[j].x;
    if (part[j].x < minx)
      minx=part[j].x;

    if (part[j].y > maxy)
      maxy=part[j].y;
    if (part[j].y < miny)
      miny=part[j].y;

    if (part[j].z > maxz)
      maxz=part[j].z;
    if (part[j].z < minz)
      minz=part[j].z;
  }
  treeinit(tree1, minx, maxx, miny, maxy, minz, maxz);
}


///////////////////////////////////////////////////////
// Function Name: treedivide 
// Usage: Divide a tree cluster in x, y, or z direction. 
// The split variable points to the tree being split, 
// and the c1 and c2 variables point to two new trees, 
// if they are needed. split and c1 may point to the 
// same tree
//
///////////////////////////////////////////////////////
// Assumptions:
//         The list of particles mut be initalized
//         and there is not a list m-repeated particles
//         (Particles must be distinct for the most part)
//
///////////////////////////////////////////////////////
// Inputs:
//         dir   - direction node is to be split in
//                 0 - x
//                 1 - y
//                 2 - z
//         split - Tree node to be split 
//         part  - list of particles
//
///////////////////////////////////////////////////////
// Outputs: (By Refrence)
//         c1    - Child node 1  
//         c2    - Child node 2
//
///////////////////////////////////////////////////////
// Functions Called:
//         treeshrink(tree,part)
//
///////////////////////////////////////////////////////
int treedivide(int dir, TREE* split, TREE* c1, TREE* c2, PARTICLE* part)  {

  int i=0;
  int j=0;
  int k=0;
  int n=0;
  int ncreate=0;
  int index;
  int count[2];
  double midpoint=0.0;

  int* cmemb[2];
  TREE* point;

  // Get the number of particles that have to be assigned
  n=split->nummemb;
  cmemb[0]=new int[n];    // Allow enough space for all particles
  cmemb[1]=new int[n];


  // Divide in the x-direction

  if(dir == 0)  {
    count[0]=0;
    count[1]=0;
    midpoint=split->xmid;

    // Count and index particles in each side
    for (j=0; j<n; j++)  {
      index = split->members[j];
      if (part[index].x <= midpoint)  {
        cmemb[0][count[0]] = index;
        count[0]++;
      }
      else  {
        cmemb[1][count[1]] = index;
        count[1]++;
      }
    }

    point=c1;  

    for (j=0; j<2; j++)  {
      if (count[j] > 0)  {
        // Assign particles on this side to the tree
        point->nummemb=count[j];

        if (point->members != NULL)
          cout << "[treedivide] : Error!  Non-null pointer!" << endl;

        point->members=new int[count[j]];
        for (k=0; k<count[j]; k++)
          point->members[k]=cmemb[j][k];

        treeshrink(point, part);   
        ncreate++;
        point=c2;
      }
    }
  }

  // Divide in the y-direction
  else if(dir == 1) {
    count[0]=0;
    count[1]=0;
    midpoint=split->ymid;

    // Count and index particles in each side
    for (j=0; j<n; j++)  {
      index = split->members[j];
      if (part[index].y <= midpoint)  {
        cmemb[0][count[0]] = index;
        count[0]++;
      }
      else  {
        cmemb[1][count[1]] = index;
        count[1]++;
      }
    }
    point=c1;
    for (j=0; j<2; j++)  {

      if (count[j] > 0)  {
        // Assign particles on this side to the tree
        point->nummemb=count[j];

	if (point->members != NULL) 
          delete [] point->members;


        point->members=new int[count[j]];
        for (k=0; k<count[j]; k++)
          point->members[k]=cmemb[j][k];
        treeshrink(point, part);
        ncreate++;
        point=c2;
      }
    }
  }

  // divide in the z-direction
  else if(dir == 2) {
    count[0]=0;
    count[1]=0;
    midpoint=split->zmid;

    // Count and index particles in each side
    for (j=0; j<n; j++)  {
      index = split->members[j];
      if (part[index].z <= midpoint)  {
        cmemb[0][count[0]] = index;
        count[0]++;
      }
      else  {
        cmemb[1][count[1]] = index;
        count[1]++;
      }
    }
    point=c1;
    for (j=0; j<2; j++)  {
      if (count[j] > 0)  {
        // Assign particles on this side to the tree
        point->nummemb=count[j];

	if (point->members != NULL)
          delete [] point->members;

        point->members=new int[count[j]];
        for (k=0; k<count[j]; k++)
          point->members[k]=cmemb[j][k];
        treeshrink(point, part);
        ncreate++;
        point=c2;
      }
    }
  }

  // Clean up arrays
  delete [] cmemb[0];
  delete [] cmemb[1];

  return(ncreate);    // Return the number of new trees created
}



///////////////////////////////////////////////////////
// Function Name: treemake
// Usage: Recursively builds oct tree by sorting particles
//        about current nodes midpoint, then buliding 
//        childern nodes.  Subdivision is stopped when
//        the number of paricles in a cluster drops
//        below a fixed min.  
//        Particle sorting is done using Quick Sort
//
///////////////////////////////////////////////////////
// Assumptions:
//         The list of particles mut be initalized
//         and there is not a list m-repeated particles
//         (Particles must be distinct for the most part)
//
//         Aspect Ratio of Box Subdivision Hard coded
//         in this function. If subdivided box has a 
//         'bad' aspect ration, avoid subdivision!
//
///////////////////////////////////////////////////////
// Inputs:
//         npart   -  total number of particles 
//         maxmemb -  maximum number of members
//         part    -  complete list of particles (+ and -)
//
///////////////////////////////////////////////////////
// Outputs:(By Reference)
//         tree    - complete particle oct tree 
//         (not by ref)
//         numtree - total number of nodes in tree
///////////////////////////////////////////////////////
// Functions Called:
//         treeinit(&tree[i], minx, maxx, miny, maxy, minz, maxz);
//         treedivide(direction, &tree[i], &tree[numtree],
//                               &tree[numtree+1], part)
//
///////////////////////////////////////////////////////
int treemake(int npart, int nterms, int maxmemb, PARTICLE* part, TREE* tree)  {

  if (verbosity) {
    cout << "[treemake]: treemake function called" << endl;
    cout << "[treemake]: npart = " << npart << endl;
    cout << "[treemake]: maxmemb = " << maxmemb << endl;
  }

  int i=0;
  int j=0;
  int k=0;
  int xcreate=0;
  int ycreate=0;
  int zcreate=0;
  int nchild=0;
  int nchild2=0;
  int ind2=0;
  int numtree=1;
  //temporary variable for testing...
  int p=0;



  // Top level tree dimensions
  double minx=part[0].x;
  double maxx=part[0].x;
  double miny=part[0].y;
  double maxy=part[0].y;
  double minz=part[0].z;
  double maxz=part[0].z;
  double aspectlimit=1.0/sqrt(2.0);   // Limit on the tree aspect ratio

  // Initialize the top level tree
  if (verbosity) {
    cout << "[treemake]: initializing top level tree" << endl;
  }
  tree[0].nummemb=npart;
  tree[0].members=new int[npart];

  if (verbosity) {
    cout << "[treemake]: shrinking tree to fit members" << endl;
  }

  // Shrink the tree to fit the particles
  for (i=0; i<npart; i++)  {
    tree[0].members[i]=i;

    if (minx > part[i].x)
      minx=part[i].x;
    if (maxx < part[i].x)
      maxx=part[i].x;

    if (miny > part[i].y)
      miny=part[i].y;
    if (maxy < part[i].y)
      maxy=part[i].y;

    if (minz > part[i].z)
      minz=part[i].z;
    if (maxz < part[i].z)
      maxz=part[i].z;
  }
  treeinit(&tree[0], minx, maxx, miny, maxy, minz, maxz);


  i=-1;

  if (verbosity) {
    cout << "[treemake]: computing remaining trees" << endl;
  }
  // Compute the rest of the tree
  while (i < numtree-1)  {
    i++;
    nchild=0;
    nchild2=0;
   
    if (verbosity) {
      cout << "[treemake]: tree[" << i << "] requires subdivision." << endl;
    }
    if (tree[i].nummemb > maxmemb)  {
      // Division is needed
     
      if (verbosity) {
	cout << "[treemake]: splitting tree[" << i 
	     << "] in x-direction." << endl;
      }
      if ((tree[i].aspectxy > aspectlimit) && 
	  (tree[i].aspectxz > aspectlimit)) {
	// Divide in the x-direction
	
        xcreate=treedivide(0, &tree[i], &tree[numtree], 
                           &tree[numtree+1], part);      
	
	
        ind2=numtree+xcreate; 
	
	for (j = 0; j<xcreate; j++) {
          // Add reference to new lower level tree
          tree[i].children[nchild]=numtree+j;  
          nchild++;
	  
	  // now check for division in the y direction
          if (tree[numtree+j].aspectxy < aspectlimit && 
	      tree[numtree+j].aspectyz > aspectlimit && 
	      tree[numtree+j].nummemb > maxmemb)  {
	    // split the tree in the y direction
	    
	    ycreate=treedivide(1, &tree[numtree+j],&tree[numtree+j],
		  	       &tree[ind2], part);                      
	    
	    if (ycreate == 2)  {
	      // add reference to lower level tree
	      tree[i].children[nchild]=ind2;
	      nchild++;
	      ind2++;
	    }
	    
	    for (k=0; k<ycreate; k++) {
	      // now check for division in the z direction
	      if (tree[numtree+j+k*(nchild2+1)].aspectxz < aspectlimit &&
		  tree[numtree+j+k*(nchild2+1)].aspectyz < aspectlimit && 
		  tree[numtree+j+k*(nchild2+1)].nummemb > maxmemb)  {
  	        // split the tree in the z direction
		
		zcreate=treedivide(2, &tree[numtree+j+k*(nchild2+1)], 
				   &tree[numtree+j+k*(nchild2+1)],  
				   &tree[ind2], part);                      
		
		
	        if (zcreate == 2) {
		  // add reference to lower level tree
		  tree[i].children[nchild]=ind2;
		  nchild2++;
		  nchild++;
		  ind2++;
	        }
	      }
	    }
 	  } else {
	    // don't split in the y direction
	    // check for division in the z direction
	    if (tree[numtree+j].aspectxz < aspectlimit &&
		tree[numtree+j].aspectyz < aspectlimit && 
		tree[numtree+j].nummemb > maxmemb)  {
  	      // split the tree in the z direction
	      
	      zcreate=treedivide(2, &tree[numtree+j], 
				 &tree[numtree+j],  
				 &tree[ind2], part);
	      
	      if (zcreate == 2) {
		// add reference to lower level tree
		tree[i].children[nchild]=ind2;
		nchild++;
		ind2++;
	      }
	    }
	  } 
	} //for
      } else {
	// don't split in the x-direction (bad aspect ratio)
	if ((tree[i].aspectxy < aspectlimit) && 
	    (tree[i].aspectyz > aspectlimit)) {
	  // Divide in the y-direction
	  
	  ycreate=treedivide(1, &tree[i], &tree[numtree], 
			     &tree[numtree+1], part);      
	  
	  ind2=numtree+ycreate; 
	  
	  for (j = 0; j<ycreate; j++) {
	    // Add reference to new lower level tree
	    tree[i].children[nchild]=numtree+j;  
	    nchild++;
	    if (tree[numtree+j].aspectxz < aspectlimit &&
		tree[numtree+j].aspectyz < aspectlimit && 
		tree[numtree+j].nummemb > maxmemb)  {
	      // split the tree in the z direction
	      zcreate=treedivide(2, &tree[numtree+j], 
				 &tree[numtree+j],  
				 &tree[ind2], part);
	      
	      if (zcreate == 2) {
		// add reference to lower level tree
		tree[i].children[nchild]=ind2;
		nchild++;
		ind2++;
	      }	  
	    }
	  }
	} else {
	  // only split in the z direction
	  
	  
	  zcreate=treedivide(2, &tree[i], &tree[numtree],
			     &tree[numtree+1], part);
	  
	  for (j=0; j<zcreate; j++)  {
	    tree[i].children[nchild]=numtree+j;
	    nchild++;
	  }
	}
      }
      numtree+=nchild;    // Increment the total number of trees
    }
  }
  if (verbosity) {
    cout << "[treemake]: completed, number of trees = " << numtree << endl;
  }
  return(numtree);
}  



///////////////////////////////////////////////////////
// Function Name: bounding_box
// Usage: Finds the bound box around the outside
//        of all particles in the system
//
///////////////////////////////////////////////////////
// Assumptions:
//          none
//
///////////////////////////////////////////////////////
// Inputs:
//          npart  - number of particles
//          part   - list of all particles
//
///////////////////////////////////////////////////////
// Outputs: (By Reference)
//          Mx     - Max x
//          mx     - Min x
//          My     - Max y
//          my     - Min y
//          Mz     - Max z
//          mz     - Min z
//
///////////////////////////////////////////////////////
// Functions Called:
//          none
//
///////////////////////////////////////////////////////
int bounding_box(int npart, PARTICLE* part,double &Mx,
		 double &mx,double &My,double &my,double &Mz,double &mz){
  int i;
  // dimensions
  mx=part[0].x;
  Mx=part[0].x;
  my=part[0].y;
  My=part[0].y;
  mz=part[0].z;
  Mz=part[0].z;
  // fit box to particles
  for (i=0; i<npart; i++)  {
    if (mx > part[i].x)
      mx=part[i].x;
    if (Mx < part[i].x)
      Mx=part[i].x;

    if (my > part[i].y)
      my=part[i].y;
    if (My < part[i].y)
      My=part[i].y;

    if (mz > part[i].z)
      mz=part[i].z;
    if (Mz < part[i].z)
      Mz=part[i].z;
  }
  return 0;
}



///////////////////////////////////////////////////////
// Function Name: init_mesh
// Usage: initalize a mesh
//
///////////////////////////////////////////////////////
// Assumptions:
//         none
//
///////////////////////////////////////////////////////
// Inputs:
//         Nx - number of particles along x direction
//         Ny - number of particles along y direction
//         Nz - number of particles along z direction
//         dlength - length of unit cell in x
//         dheight - length of unit cell in y
//         ddepth  - length of unit cell in z
//         part    - list of all particles 
//
///////////////////////////////////////////////////////
// Outputs:
//         part    - location data
//
///////////////////////////////////////////////////////
// Functions Called: none
//
///////////////////////////////////////////////////////
void init_mesh(int Nx, int Ny, int Nz,
		  double dlength, double dheight, double ddepth, 
		  PARTICLE* part)  {
  double dx = 2*dlength/(Nx-1);
  double dy = 2*dheight/(Ny-1);
  double dz = 2*ddepth/(Nz-1);

  double weight = dx*dy*dz;

  int ind = 0;

  for (int i=0;i<Nx; i++) {
    for (int j=0;j<Ny; j++) {
      for (int k=0;k<Nz; k++) {
	part[ind].x = -dlength+i*dx;
	part[ind].y = -dheight+j*dy;
	part[ind].z = -ddepth+k*dz;
	part[ind].weight = weight;
	ind++;
      }
    }
  }
}


///////////////////////////////////////////////////////
// Function Name: boxparticles
// Usage: initalize particles in domain
//        -bound_box
//        -free_space
//        -periodic
//
///////////////////////////////////////////////////////
// Assumptions:
//         none
//
///////////////////////////////////////////////////////
// Inputs:
//         npart   - total number of macro particles
//         dlength - length of unit cell in x
//         dheight - length of unit cell in y
//         ddepth  - length of unit cell in z
//         part    - empty list of all particles 
//         density - plasma density
//
///////////////////////////////////////////////////////
// Outputs:
//         part    - all data
//
///////////////////////////////////////////////////////
// Functions Called:
//         RandL() - tree3d_v1.h
//         part_vol()
//         bounding_box()
//
///////////////////////////////////////////////////////
// 
void boxparticles(int npart, double dlength, double dheight, double ddepth, 
		  PARTICLE* part, double density)  {
  
  double dx,dy,dz;  // Spacing between particles for uniform spacing
  double weight;  // Particle weight  

  // Set the variables for uniformly-spaced particles
  weight = density * dlength * dheight *ddepth/ npart; 
  
  //srand((unsigned)time(0));
  srand(100);
  dx = dlength/2;
  dy = dheight/2;
  dz = ddepth/2;

  double temp;

  for (int i=0; i<npart; i++)  {
    part[i].charge = -1.0;
    part[i].weight = weight;

    temp = dx/2 + RandL()*dx;
    part[i].x=temp;
    temp = dy/2 + RandL()*dy;
    part[i].y=temp;
    temp = dz/2 + RandL()*dz;
    part[i].z = temp;

    // initial velocity set to zero
    // **MAYBE REQUIRES A BETTER CHOICE?***
    temp = (RandL()-0.5)*1e6;
    part[i].u=temp;
    temp = (RandL()-0.5)*1e6;
    part[i].v=temp;
    temp = (RandL()-0.5)*1e6;
    part[i].w=temp;

    part[i].tot_charge=part[i].charge*part[i].weight*elec;


  }
  return;
}


///////////////////////////////////////////////////////
// Function Name: part_vel 
// Usage: initialize particle velocity distribution
//        and macro particle charge
//
///////////////////////////////////////////////////////
// Assumptions:
//          none
//
///////////////////////////////////////////////////////
// Inputs:
//          part        - list of all particles
//          hot_or_cold - flag indicating hot or cold
//                        velocity distribution
//          s           - start particle index 
//          cnt         - end particle index
//
///////////////////////////////////////////////////////
// Outputs: (by refrence)
//          part.u
//          part.v
//          part.w
//          part.charge
//
///////////////////////////////////////////////////////
// Functions Called:
//          RandL()
//
///////////////////////////////////////////////////////

int part_vel(PARTICLE* part,int hot_or_cold,int s, int cnt){
  int i;
  double teref = 0.1;
  double mass = 9.1093897e-31;
  double v0 = sqrt(8*teref*elec/(pi*mass));     // Thermal velocity
  double beta = mass/(2*elec*teref);            // Maxwellian parameter
  double f=0.0;  // Fraction of Maxwellian
  double v=0.0;  // Particle velocity
  double theta;  // Velocity angle 
  double psi;    // Velocity angle  

  //Velocity Distribution
  if (hot_or_cold){
    // hot stream -- sample velocity from a Maxwellian
    for (i=s; i<cnt; i++)  {
      f=0.0;
      while (f < RandL())  {
        v=v0*(3.0-6.0*RandL());
        f=exp(-beta*v*v);
      }
      theta=2*pi*RandL();
      psi=pi*RandL();

      part[i].u=v*cos(theta)*sin(psi);
      part[i].v=v*sin(theta)*sin(psi);
      part[i].w=v*cos(psi);
      part[i].tot_charge=part[i].charge*part[i].weight*elec;
    }
  }
  else{  
    for (i=s; i<cnt; i++)  {
      part[i].u=0.0;
      part[i].v=0.0;
      part[i].w=0.0;
      part[i].tot_charge=part[i].charge*part[i].weight*elec;
    }
  }
  return 0;
}


///////////////////////////////////////////////////////
// Function Name: force_particle_ds 
// Usage: direct sum of force at particle index 
//        (using regularized kernel)
//
///////////////////////////////////////////////////////
// Assumptions:
//          none
//
///////////////////////////////////////////////////////
// Inputs:
//          ind    - index of particle
//          npart  - number of particles
//          part   - list of all particles
//
///////////////////////////////////////////////////////
// Outputs: (by reference)
//          force   - (xfor,yfor,zfor)
//
///////////////////////////////////////////////////////
// Functions Called:
//          none
//
///////////////////////////////////////////////////////
void force_particle_ds(int ind, int npart, PARTICLE* part, 
		       double *xfor, double *yfor, double *zfor) {

  double temp=0.0;
  double x_force=0.0;
  double y_force=0.0;
  double z_force=0.0;
  double dx;
  double dy;
  double dz;
  double rs;
  double r;

  for (int i=0; i<npart; i++)  {
    if (i != ind)  {
      // This skips the particle that is the input
      dx=part[i].x-part[ind].x;
      dy=part[i].y-part[ind].y;
      dz=part[i].z-part[ind].z;

      // regularization, needed for high order time convergence
      rs= dx*dx + dy*dy + dz*dz + DEL;
      r = sqrt (rs);

      temp=part[i].tot_charge/(r*rs);
      x_force-=dx*temp;
      y_force-=dy*temp;
      z_force-=dz*temp;
    } // endif not current index
  } //end loop over npart
 
  *xfor = x_force;
  *yfor = y_force;
  *zfor = z_force;
     

} // end force_particle_ds


///////////////////////////////////////////////////////
// Function Name: shielded_force_particle_ds 
// Usage: direct sum of screeened force at particle index 
//
///////////////////////////////////////////////////////
// Assumptions:
//          this code is specifically for a current in the z-direction
//
///////////////////////////////////////////////////////
// Inputs:
//          ind    - index of particle
//          npart  - number of particles
//          part   - list of all particles
//          alpha  - helmholtz parameter
//
///////////////////////////////////////////////////////
// Outputs: (by reference)
//          force   - (xfor,yfor,zfor)
//
///////////////////////////////////////////////////////
// Functions Called:
//          none
//
///////////////////////////////////////////////////////
void shielded_force_particle_ds(int ind, int npart, PARTICLE* part, 
				double *xfor, double *yfor, double*zfor,
				double alpha) {

  double temp=0.0;
  double x_force=0.0;
  double y_force=0.0;
  double z_force=0.0;
  double dx;
  double dy;
  double dz;
  double rs;
  double r;


  for (int i=0; i<npart; i++)  {
    if (i != ind)  {
      // This skips the particle that is the input
      dx=part[i].x-part[ind].x;
      dy=part[i].y-part[ind].y;
      dz=part[i].z-part[ind].z;

      // *** NOTE: REMOVING REGULARIZATION
      rs= dx*dx + dy*dy + dz*dz;
      r = sqrt (rs);

      temp=part[i].tot_charge*exp(-alpha*r)/rs*(alpha+1/r);
      x_force -= dy*temp;
      y_force += dx*temp;
      z_force -= 0.0;
    } // endif not current index
  } //end loop over npart

  *xfor = -x_force;
  *yfor = -y_force;
  *zfor = -z_force;
     
} // end shielded_force_particle_ds 



///////////////////////////////////////////////////////
// Function Name: shielded_pot_particle_ds 
// Usage: direct sum of screeened force at particle index 
//
///////////////////////////////////////////////////////
// Assumptions:
//          this code is specifically for a current in the z-direction
//
///////////////////////////////////////////////////////
// Inputs:
//          ind    - index of particle
//          npart  - number of particles
//          part   - list of all particles
//          alpha  - helmholtz parameter
//
///////////////////////////////////////////////////////
// Outputs: (by reference)
//          pot - potential
//
///////////////////////////////////////////////////////
// Functions Called:
//          none
//
///////////////////////////////////////////////////////
void shielded_pot_particle_ds(int ind, int npart, PARTICLE* part, 
			      double *pot,double alpha) {

  double temp=0.0;
  double dx;
  double dy;
  double dz;
  double rs;
  double r;

  if (verbosity){
    cout << "[shielded_pot_particle_ds]:" << endl;
    cout << "...location of particle[" << ind << "] = (" 
	 << part[ind].x << ", " 
	 << part[ind].y << ", " 
	 << part[ind].z << ")" << endl;
    cout << "...for particle[" << ind << "], contributions from" << endl;
  }

  for (int i=0; i<npart; i++)  {
    if (i != ind)  {
      if (verbosity) {
	cout << "...... particle[" << i << "] located at (" 
	     << part[i].x << ", " 
	     << part[i].y << ", " 
	     << part[i].z << ")" << endl;
      }
      // This skips the particle that is the input
      dx=part[i].x-part[ind].x;
      dy=part[i].y-part[ind].y;
      dz=part[i].z-part[ind].z;
      
      // *** NOTE: REMOVING REGULARIZATION
      rs= dx*dx + dy*dy + dz*dz;
      r = sqrt (rs);
      temp += part[i].tot_charge*exp(-alpha*r)/r;
      
      if (verbosity) {
	cout << "......... (dx,dy,dz) = ( " << dx << ", " << dy << ", " 
	     << dz << ")" << endl;
	cout << "......... part[" << i << "].tot_charge = " 
	     << part[i].tot_charge << endl;
	cout << "......... rs = " << rs << endl;
	cout << "......... r = " << r << endl;
	cout << "......... alpha = " << alpha << endl;
      }
    } // endif not current index
  } //end loop over npart
  *pot = -temp;
     
} // end shielded_pot_particle_ds 



///////////////////////////////////////////////////////
// Function Name: force_point_ds - 
// Usage: direct sum of force at a point (x,y,z)
//
///////////////////////////////////////////////////////
// Assumptions: uses regularized kernel
//
///////////////////////////////////////////////////////
// Inputs:
//          x     - x location
//          y     - y location
//          z     - z location
//          npart - number of particles
//          part  - list of all particles
//
///////////////////////////////////////////////////////
// Outputs: (by reference)
//          xfor  - x force
//          yfor  - y force
//          zfor  - z force
//
///////////////////////////////////////////////////////
// Functions Called:
//
///////////////////////////////////////////////////////
void force_point_ds(double x, double y, double z,
		    int npart, PARTICLE* part, 
		    double *xfor, double *yfor, double* zfor)  {
  int i;
  double temp=0.0;
  double x_force=0.0;
  double y_force=0.0;
  double z_force=0.0;
  double dx;
  double dy;
  double dz;
  double rs;
  double r;

  for (i=0; i<npart; i++)  {
    dx=part[i].x-x;
    dy=part[i].y-y;
    dz=part[i].z-z;

    rs= dx*dx + dy*dy + dz*dz + DEL;
    r = sqrt (rs);

    temp=part[i].tot_charge/(r*rs);
    x_force-=dx*temp;
    y_force-=dy*temp;
    z_force-=dz*temp;

  }

  *xfor = x_force;
  *yfor = y_force;
  *zfor = z_force;

} // force_point_ds


///////////////////////////////////////////////////////
// Function Name: shielded_force_point_ds - 
// Usage: direct sum of force at a point (x,y,z)
//
///////////////////////////////////////////////////////
// Assumptions: uses regularized kernel
//
///////////////////////////////////////////////////////
// Inputs:
//          x     - x location
//          y     - y location
//          z     - z location
//          npart - number of particles
//          part  - list of all particles
//          alpha - helmholtz parameter
//
///////////////////////////////////////////////////////
// Outputs: (by reference)
//          xfor  - x force
//          yfor  - y force
//          zfor  - z force
//
///////////////////////////////////////////////////////
// Functions Called:
//
///////////////////////////////////////////////////////
void shielded_force_point_ds(double x, double y, double z,
			     int npart, PARTICLE* part, 
			     double *xfor, double *yfor, double* zfor,
			     double alpha)  {
  int i;
  double temp=0.0;
  double x_force=0.0;
  double y_force=0.0;
  double z_force=0.0;
  double dx;
  double dy;
  double dz;
  double rs;
  double r;
  double expr;

  for (i=0; i<npart; i++)  {
    dx=part[i].x-x;
    dy=part[i].y-y;
    dz=part[i].z-z;

    // *** NOTE, removing regularization!!!
    rs= dx*dx + dy*dy + dz*dz;
    r = sqrt (rs);
    expr = exp(-alpha*r)/rs;

    temp=part[i].tot_charge*expr*(alpha + 1/r);
    x_force-=dy*temp;
    y_force= dx*temp;
    z_force = 0.0;
  }

  *xfor = -x_force;
  *yfor = -y_force;
  *zfor = -z_force;

} // end shielded_force_point_ds


///////////////////////////////////////////////////////
// Function Name: force_point_ds - 
// Usage: direct sum of force at a point (x,y,z)
//
///////////////////////////////////////////////////////
// Assumptions: uses regularized kernel
//
///////////////////////////////////////////////////////
// Inputs:
//          x     - x location
//          y     - y location
//          z     - z location
//          npart - number of particles
//          part  - list of all particles
//
///////////////////////////////////////////////////////
// Outputs: (by reference)
//          xfor  - x force
//          yfor  - y force
//          zfor  - z force
//
///////////////////////////////////////////////////////
// Functions Called:
//
///////////////////////////////////////////////////////
void force_point_ds(double x, double y, double z,
		    int npart, double* xfor, double* yfor,
		    double* zfor, PARTICLE* part)  {
  int i;
  double temp=0.0;
  double x_force=0.0;
  double y_force=0.0;
  double z_force=0.0;
  double dx;
  double dy;
  double dz;
  double rs;
  double r;

  for (i=0; i<npart; i++)  {
    dx=part[i].x-x;
    dy=part[i].y-y;
    dz=part[i].z-z;

    rs= dx*dx + dy*dy + dz*dz + DEL;
    r = sqrt (rs);

    temp=part[i].tot_charge/(r*rs);
    x_force-=dx*temp;
    y_force-=dy*temp;
    z_force-=dz*temp;

  }

  *xfor = x_force;
  *yfor = y_force;
  *zfor = z_force;

} // force_helmholtz_ds


///////////////////////////////////////////////////////
// Function Name: pot_point_ds - regularized (del)
// Usage: calculate potential at a point due to all 
//        particles, using direct sum.
//
///////////////////////////////////////////////////////
// Inputs:
//          x     - x location
//          y     - y location
//          z     - z location
//          npart - number of particles
//          part  - list of all particles
//
///////////////////////////////////////////////////////
// Outputs: (by reference)
//          pot - potential
//
///////////////////////////////////////////////////////
// Functions Called:
//
///////////////////////////////////////////////////////
void pot_point_ds(double x, double y, double z,
		  int npart, PARTICLE* part, double pot)  {

  double dx;
  double dy;
  double dz;
  double rs;
  double r;

  pot = 0.0;

  for (int i=0; i<npart; i++)  {
    dx=part[i].x-x;
    dy=part[i].y-y;
    dz=part[i].z-z;

    rs= dx*dx + dy*dy + dz*dz + DEL;
    r = sqrt (rs);

    pot-=0.5*inv2pi/r*part[i].tot_charge;
  }
}


///////////////////////////////////////////////////////
// Function Name: freespace_particles_ds
// Usage: Calculate force on each particle using direct summation 
//
///////////////////////////////////////////////////////
// Assumptions:
//           none
//
///////////////////////////////////////////////////////
// Inputs:
//           npart    - number of macro particles
//           part     - list of all particles 
//
///////////////////////////////////////////////////////
// Outputs: (by reference)
//           part.xforce
//           part.yforce
//           part.zforce
//           
//
///////////////////////////////////////////////////////
// Functions Called:
//           force_particle_ds()
//
///////////////////////////////////////////////////////
void freespace_particles_ds(int npart, PARTICLE* part)  {
  
  double xfor;
  double yfor;
  double zfor;
  
  // Compute the particle forces
  for (int i=0; i<npart; i++)  {
    if (verbosity) {
      cout << "[freespace_particle_ds]: calculating force on particle ["
	   << i << "]" << endl;
    }
    
    // Use direct sum to compute the electrostatic force from particles
    force_particle_ds(i, npart, part, &xfor, &yfor, &zfor);
    part[i].xforce=xfor*charge_constant*part[i].tot_charge;
    part[i].yforce=yfor*charge_constant*part[i].tot_charge;
    part[i].zforce=zfor*charge_constant*part[i].tot_charge;
  } // end particle loop
} // freespace_particles_ds



///////////////////////////////////////////////////////
// Function Name: shielded_freespace_particles_ds
// Usage: Calculate shielded force on each particle using direct summation 
//
///////////////////////////////////////////////////////
// Assumptions:
//           none
//
///////////////////////////////////////////////////////
// Inputs:
//           npart    - number of macro particles
//           part     - list of all particles 
//           alpha    - helmholtz parameter
//
///////////////////////////////////////////////////////
// Outputs: (by reference)
//           part.xforce
//           part.yforce
//           part.zforce
//           
//
///////////////////////////////////////////////////////
// Functions Called:
//           force_particle_ds()
//
///////////////////////////////////////////////////////
void shielded_freespace_particles_ds(int npart, PARTICLE* part, double alpha)  {
  
  double xfor;
  double yfor;
  double zfor;
  
  // Compute the particle forces
  for (int i=0; i<npart; i++)  {
    if (verbosity) {
      cout << "[shielded_freespace_particles_ds]: calculating force on particle ["
	   << i << "]" << endl;
    }
    
    // Use direct sum to compute the electrostatic force from particles
    shielded_force_particle_ds(i, npart, part, &xfor, &yfor, &zfor, alpha);
    part[i].xforce=xfor*charge_constant*part[i].tot_charge;
    part[i].yforce=yfor*charge_constant*part[i].tot_charge;
    part[i].zforce=zfor*charge_constant*part[i].tot_charge;
  } // end particle loop
} // end shielded_freespace_particles_ds


///////////////////////////////////////////////////////
// Function Name: shielded_freespace_pot_ds
// Usage: Calculate shielded force on each particle using direct summation 
//
///////////////////////////////////////////////////////
// Assumptions:
//           none
//
///////////////////////////////////////////////////////
// Inputs:
//           npart    - number of macro particles
//           part     - list of all particles 
//           alpha    - helmholtz parameter
//
///////////////////////////////////////////////////////
// Outputs: (by reference)
//           part.xforce
//           part.yforce
//           part.zforce
//           
//
///////////////////////////////////////////////////////
// Functions Called:
//           force_particle_ds()
//
///////////////////////////////////////////////////////
void shielded_freespace_pot_ds(int npart, PARTICLE* part, double alpha)  {
  
  double pot;
  
  // Compute the particle forces
  for (int i=0; i<npart; i++)  {
    if (verbosity) {
      cout << "[shielded_freespace_pot_ds]: calculating pot on particle ["
	   << i << "]" << endl;
    }
    
    // Use direct sum to compute the electrostatic force from particles
    shielded_pot_particle_ds(i, npart, part, &pot, alpha);
    part[i].xforce=pot*charge_constant*part[i].tot_charge;
  } // end particle loop
} // end shielded_freespace_pot_ds



///////////////////////////////////////////////////////
// Function Name: freespace_particles_tree
// Usage: Calculate force on each particle using the treecode
//
///////////////////////////////////////////////////////
// Assumptions: free space
//
///////////////////////////////////////////////////////
// Inputs:
//         npart    - number of macro particles
//         part     - list of all particles 
//         numtree - total number of tree nodes 
//         nterms  - number of terms in Taylor expantion
//         tree    - oct tree of particles, sorted
//         maxmemb - maximum number of particles per tree
//
///////////////////////////////////////////////////////
// Outputs: (by reference)
//           part[j].x_force
//           part[j].y_force
//           part[j].z_force
//           
//
///////////////////////////////////////////////////////
// Functions Called:
//           treemake, treemoments, treeforce, treeempty
//
///////////////////////////////////////////////////////
void freespace_particles_tree(int npart, PARTICLE* part, int nterms, 
			     double acc, TREE* tree, int numtree, 
			     int maxmemb)  {

  numtree=treemake(npart, nterms, maxmemb,  part, tree);

  //if (verbosity) {
    cout << "numtree = " << numtree << endl;
    //}
  for (int j = 0; j<numtree; j++) {
    if (verbosity) {
      cout << "tree[" << j << "]" << endl;
      cout << "... nummemb = " << tree[j].nummemb << endl;
      cout << "... x0 = " << tree[j].x0 << endl;
      cout << "... x1 = " << tree[j].x1 << endl;
      cout << "... y0 = " << tree[j].y0 << endl;
      cout << "... y1 = " << tree[j].y1 << endl;
      cout << "... z0 = " << tree[j].z0 << endl;
      cout << "... z1 = " << tree[j].z1 << endl;
      cout << "... xmid = " << tree[j].xmid << endl;
      cout << "... ymid = " << tree[j].ymid << endl;
      cout << "... zmid = " << tree[j].zmid << endl;
      for (int i=0; i< NUMCHILD; i++) {
	cout << "... numchild[" << i <<"] =" << tree[j].children[i] << endl;
      }
      for (int i=0; i<tree[j].nummemb; i++) {
	cout << "... member[" << i <<"] =" << tree[j].members[i] << endl;
      }
    }
  }

  if (verbosity) {
    cout << "calculating tree moments" << endl;
  }
  treemoments(numtree,nterms,tree,part);
  if (verbosity) {
    cout << "... finished calculating tree moments" << endl;
  }
  
  double xfor;
  double yfor;
  double zfor;
   
  // Compute the particle forces
  for (int i=0; i<npart; i++)  {
    if (verbosity) {
      cout << "...[freespace_particles_tree] : computing force for particle "
	   << i << endl;
    }
    xfor = 0.0;
    yfor = 0.0;
    zfor = 0.0;
    
    // Use the treecode to compute the electrostatic force from particles
    
    treeforce(0,i,nterms,part[i].x,part[i].y,part[i].z,
	      acc,&xfor,&yfor,&zfor,tree,part);
    
    part[i].xforce=xfor*charge_constant*part[i].tot_charge;
    part[i].yforce=yfor*charge_constant*part[i].tot_charge;
    part[i].zforce=zfor*charge_constant*part[i].tot_charge;
    
  } // end particle loop
  
  //empty tree
  treeempty(numtree,tree);
  
} //freespace_particles_tree


/////////////////////////////////////////////////////
// Function Name: shielded_freespace_particles_tree
// Usage: Calculate shielded force on each particle
//         using the treecode
//
///////////////////////////////////////////////////////
// Assumptions: free space
//
///////////////////////////////////////////////////////
// Inputs:
//         npart    - number of macro particles
//         part     - list of all particles 
//         alpha    - shielding parameter
//         numtree - total number of tree nodes 
//         nterms  - number of terms in Taylor expantion
//         tree    - oct tree of particles, sorted
//         maxmemb - maximum number of particles per tree
//
///////////////////////////////////////////////////////
// Outputs: (by reference)
//           part[j].x_force
//           part[j].y_force
//           part[j].z_force
//           
//
///////////////////////////////////////////////////////
// Functions Called:
//           treemake, treemoments, treeforce, treeempty
//
///////////////////////////////////////////////////////
void shielded_freespace_particles_tree(int npart, PARTICLE* part, 
				       double alpha, int nterms, 
				       double acc, TREE* tree, 
				       int numtree, int maxmemb)  {

  numtree=treemake(npart, nterms, maxmemb,  part, tree);

  if (verbosity) {
    cout << "numtree = " << numtree << endl;
  }
  for (int j = 0; j<numtree; j++) {
    if (verbosity) {
      cout << "tree[" << j << "]" << endl;
      cout << "... nummemb = " << tree[j].nummemb << endl;
      cout << "... x0 = " << tree[j].x0 << endl;
      cout << "... x1 = " << tree[j].x1 << endl;
      cout << "... y0 = " << tree[j].y0 << endl;
      cout << "... y1 = " << tree[j].y1 << endl;
      cout << "... z0 = " << tree[j].z0 << endl;
      cout << "... z1 = " << tree[j].z1 << endl;
      cout << "... xmid = " << tree[j].xmid << endl;
      cout << "... ymid = " << tree[j].ymid << endl;
      cout << "... zmid = " << tree[j].zmid << endl;
      for (int i=0; i< NUMCHILD; i++) {
	cout << "... numchild[" << i <<"] =" << tree[j].children[i] << endl;
      }
      for (int i=0; i<tree[j].nummemb; i++) {
	cout << "... member[" << i <<"] =" << tree[j].members[i] << endl;
      }
    }
  }

  if (verbosity) {
    cout << "calculating tree moments" << endl;
  }
  treemoments(numtree,nterms,tree,part);
  if (verbosity) {
    cout << "... finished calculating tree moments" << endl;
  }
  
  double xfor;
  double yfor;
  double zfor;
   
  // Compute the particle forces
  for (int i=0; i<npart; i++)  {
    if (verbosity) {
      cout << "...[freespace_particles_tree] : computing force for particle "
	   << i << endl;
    }
    xfor = 0.0;
    yfor = 0.0;
    zfor = 0.0;
    
    // Use the treecode to compute the electrostatic force from particles
    
    shielded_treeforce(0,i,nterms,part[i].x,part[i].y,part[i].z,alpha,
		       acc,&xfor,&yfor,&zfor,tree,part);
    
    part[i].xforce=xfor*charge_constant*part[i].tot_charge;
    part[i].yforce=yfor*charge_constant*part[i].tot_charge;
    part[i].zforce=zfor*charge_constant*part[i].tot_charge;
    
  } // end particle loop
  
  //empty tree
  treeempty(numtree,tree);
  
}//end shielded_freespace_particles_tree


///////////////////////////////////////////////////////
// Function Name: box_calc_ds
// Usage: Calculate force on each particle using direct summation 
//
///////////////////////////////////////////////////////
// Assumptions:
//           none
//
///////////////////////////////////////////////////////
// Inputs:
//           npart    - number of macro particles
//           part     - list of all particles 
//           bpart    - effective boundary particles
//           numpan   - number of panels
//           pan      - array of panels
//           panelarray - matrix for computing panel strengths
//           (dlength,dheight,ddepth) - dimensions of domain
//           density  - density of particles
//
///////////////////////////////////////////////////////
// Outputs: (by reference)
//           part.xforce
//           part.yforce
//           part.zforce
//           
//
///////////////////////////////////////////////////////
// Functions Called:
//           force_particle_ds()
//
///////////////////////////////////////////////////////
int box_calc_ds(int npart, PARTICLE* part, int numpan, PANEL* pan,	
		double** panelarray, PARTICLE* bpart, double dlength,
		double dheight, double ddepth, double density)  {
 
  double xfor;
  double yfor;
  double zfor;
  
  // compute panel strengths
  //panel_solve_ds(numpan,npart,part,pan,panelarray);
  //effective_charge(numpan,pan,bpart,dlength,dheight,ddepth,density,npart);

  // Compute the particle forces
  for (int i=0; i<npart; i++)  {

    // Use direct sum to compute the electrostatic force from particles
    force_particle_ds(i, npart, part, &xfor, &yfor, &zfor);
    part[i].xforce=xfor*charge_constant;
    part[i].yforce=yfor*charge_constant;
    part[i].zforce=zfor*charge_constant;
    
    // Use direct sum to compute the force from boundary
    //xfor=0.0;
    //yfor=0.0;
    //zfor=0.0;
    //force_point_ds(part[i].x, part[i].y, part[i].z, 
    //		  numpan, bpart,&xfor, &yfor, &zfor);
    //part[i].xforce+=xfor*charge_constant;
    //part[i].yforce+=yfor*charge_constant;
    //part[i].zforce+=zfor*charge_constant;
    
    // add in v x B effect, assume a 1T field
    //part[i].xforce+=part[i].v*charge_constant;
    //part[i].yforce+=-part[i].u*charge_constant;
    
    // Multiply by the particle weight
    part[i].xforce*=part[i].tot_charge;
    part[i].yforce*=part[i].tot_charge;
    part[i].zforce*=part[i].tot_charge;
  } // end particle loop
  return(0);
} // box_calc_ds


///////////////////////////////////////////////////////
// Function Name: box_calc_tree
// Usage: Calculate force on each particle using the treecode
//
///////////////////////////////////////////////////////
// Assumptions:
//           none
//
///////////////////////////////////////////////////////
// Inputs:
//           npart    - number of macro particles
//           part     - list of all particles 
//           bpart    - effective boundary particles
//           numpan   - number of panels
//           pan      - array of panels
//           panelarray - matrix for computing panel strengths
//           (dlength,dheight,ddepth) - dimensions of domain
//           density  - density of particles
//
///////////////////////////////////////////////////////
// Outputs: (by reference)
//           part.x_force
//           part.y_force
//           part.z_force
//           
//
///////////////////////////////////////////////////////
// Functions Called:
//           treemake, treemoments, treeforce, treeempty
//
///////////////////////////////////////////////////////
int box_calc_tree(int npart, PARTICLE* part, int numpan, PANEL* pan,	
		  double** panelarray, PARTICLE* bpart, double dlength,
		  double dheight, double ddepth, double density,
		  int nterms, double acc, TREE* tree, int numtree,
		  int maxmemb)  {

  int i;              // Index variable
  char output[128];
  int index=1;

  numtree=treemake(npart, nterms, maxmemb,  part, tree);
  treemoments(numtree,nterms,tree,part);
   
  // compute panel strengths using treecode
  //panel_solve_tree(numpan,npart,nterms,acc,panelarray,part,pan,tree);
  //effective_charge(numpan,pan,bpart,dlength,dheight,ddepth,density,npart);


  double xfor;
  double yfor;
  double zfor;

  //TREECODE

  // Compute the particle forces
  for (i=0; i<npart; i++)  {
    // Use the treecode to compute the electrostatic force from particles
    xfor=0.0;
    yfor=0.0;
    zfor=0.0;
    treeforce(0,-1,nterms,part[i].x,part[i].y,part[i].z,acc,
	      &xfor,&yfor,&zfor, tree, part);
    part[i].xforce+=xfor*charge_constant;
    part[i].yforce+=yfor*charge_constant;
    part[i].zforce+=zfor*charge_constant;
    
    // Use direct sum to compute the force from boundary
    //xfor=0.0;
    //yfor=0.0;
    //zfor=0.0;
    //force_point_ds(part[i].x, part[i].y, part[i].z, 
    //	       numpan, bpart,&xfor, &yfor, &zfor);
    //part[i].xforce+=xfor*charge_constant;
    //part[i].yforce+=yfor*charge_constant;
    //part[i].zforce+=zfor*charge_constant;
    
    // add in v x B effect, assume a 1T field
    //part[i].xforce+=part[i].v*charge_constant;
    //part[i].yforce+=-part[i].u*charge_constant;
    
    // Multiply by the particle weight
    part[i].xforce*=part[i].tot_charge;
    part[i].yforce*=part[i].tot_charge;
    part[i].zforce*=part[i].tot_charge;
  } // end particle loop

  //empty tree
  treeempty(numtree,tree);
  
  return(0);
} //box_calc_tree



///////////////////////////////////////////////////////
// Function Name: treemoments
// Usage: Compute moments for each cluster of particles
//
///////////////////////////////////////////////////////
// Assumptions:
//          none
//
///////////////////////////////////////////////////////
// Inputs:
//         numtree - total number of tree nodes 
//         nterms  - number of terms in Taylor expantion
//         tree    - oct tree of particles, sorted
//         part    - list of particles 
//
///////////////////////////////////////////////////////
// Outputs:(By Reference)
//         tree.moments - Creates and fills in moments       
//
///////////////////////////////////////////////////////
// Functions Called:
//          none
//
///////////////////////////////////////////////////////
void treemoments(int numtree, int nterms, TREE* tree, PARTICLE* part)  {
  int m;

  double dx;
  double dy;
  double dz;
  double mult;
  double xmul;
  double ymul;
  //double zmul;

  // Compute the factorial terms
  mult=1.0;

  if (verbosity) {
    cout << "... [treemoments]: number of trees = " << numtree << endl;
  }
  for(int ind=0; ind<numtree; ind++)  {
    //tree[ind].moments = NULL;
    
    if (verbosity) {
      cout << "... [treemoments]: allocating memory force tree[" << ind << "]" 
	   << endl;
    }

    //if (tree[ind].moments == NULL) {
    //  cout << "tree[" << ind << "] moments not initialized" << endl;
    //} else { 
    //  cout << "tree[" << ind << "] moments already initialized" << endl;
    // }
    if (tree[ind].moments == NULL){
      // allocate moments
      tree[ind].moments = new double**[nterms];
      for (int i=0; i<nterms; i++) {
	tree[ind].moments[i] = new double*[nterms];
	for (int j=0; j<nterms; j++) {
	  tree[ind].moments[i][j] = new double[nterms];
	}
      }
    } // end allocate moments
    
    if (verbosity) {
      cout << "... [treemoments]: setting moments to zero" << endl;    
    }
    // set all moments to zero
    for (int i=0; i<nterms; i++) {
      for (int j=0; j<nterms; j++) {
	for (int k=0; k<nterms; k++) {
	  tree[ind].moments[i][j][k] = 0.0;
	}
      }
    } // all moments set to zero
    
    
    
    for(int j1=0; j1<tree[ind].nummemb; j1++)  {
      int j;
      j=tree[ind].members[j1];
      dx=part[j].x-tree[ind].xmid;
      dy=part[j].y-tree[ind].ymid;
      dz=part[j].z-tree[ind].zmid;
      mult=part[j].tot_charge;
      ymul=mult;
      xmul=mult;
      
      for (int i = 0; i<nterms; i++) {
	for (int j = 0; j<nterms-i; j++) {
	  for (int k = 0; k<nterms-i-j; k++) {
	    tree[ind].moments[i][j][k]+=mult;
	    mult*=dz;
	  }
	  ymul*=dy;
	  mult=ymul;
	}
	xmul*=dx;
	ymul=xmul;
	mult=xmul;
      }
    } // finished sum over each particle in tree[ind]
  } // loop over all trees
  return;
}



///////////////////////////////////////////////////////
// Function Name: treetaylorforce
// Usage: Compute force for a cluster using far 
//        field taylor expansion
//
///////////////////////////////////////////////////////
// Assumptions:
//        none
//
///////////////////////////////////////////////////////
// Inputs:
//        nterms - number of terms in taylor expansion 
//        dx     - x_part - x_cluster_center
//        dy     - y_part - y_cluster_center
//        dz     - z_part - z_cluster_center
//        rs     - sqrt(dx^2+dy^2+dz^2+del)
//        moment - list of moments for given cluster
//
///////////////////////////////////////////////////////
// Outputs: (By reference)
//         force   - vector containing (fx,fy,fz);
//
///////////////////////////////////////////////////////
// Functions Called:
//         none
//
///////////////////////////////////////////////////////
void treetaylorforce(int nterms, double dx, double dy, double dz,
		     double rs, double *xfor, double *yfor, double *zfor,
		     double*** moment)  {

  if (verbosity) {
    cout << "[treetaylorforce] called" << endl;
  }

  double fac=1.0/rs;
  double sqfac=sqrt(fac);
  double ik;  // we will set this to 1/k
  double ii;  // we will set this to 1/i

  double derivs[MAXTERM][MAXTERM][MAXTERM];
  
  // zeroth order derivative
  derivs[0][0][0]=sqfac; 

  // first order derivatives
  derivs[1][0][0]=fac*dx*sqfac;       
  derivs[0][1][0]=fac*dy*sqfac;
  derivs[0][0][1]=fac*dz*sqfac;

  // recursion relation for derivatives in one direction
  for (int k=2; k<nterms+1;k++) {
    ik = 1.0/k;
    derivs[k][0][0] = fac*(2.0-ik)*derivs[k-1][0][0]*dx -
      fac*(1-ik)*derivs[k-2][0][0];

    derivs[0][k][0] = fac*(2.0-ik)*derivs[0][k-1][0]*dy -
      fac*(1-ik)*derivs[0][k-2][0];

    derivs[0][0][k] = fac*(2.0-ik)*derivs[0][0][k-1]*dz -
      fac*(1-ik)*derivs[0][0][k-2];
  }

  // derivatives for 
  // G(i,1,0), G(i,0,1), 
  // G(1,i,0), G(0,i,1),
  // G(1,0,i), G(0,1,i)


  derivs[1][1][0] = fac*dx*derivs[0][1][0] + 
    fac*2.0*dy*derivs[1][0][0];
 
  derivs[1][0][1] = fac*dx*derivs[0][0][1] + 
    fac*2.0*dz*derivs[1][0][0];

  derivs[0][1][1] = fac*dy*derivs[0][0][1] + 
    fac*2.0*dz*derivs[0][1][0];

  for (int k=2; k<nterms; k++) {
    derivs[1][0][k]=fac*(dx*derivs[0][0][k]+
			 2.0*dz*derivs[1][0][k-1]-
			 derivs[1][0][k-2]);
    
    derivs[0][1][k]=fac*(dy*derivs[0][0][k]+
			 2.0*dz*derivs[0][1][k-1]-
			 derivs[0][1][k-2]);
    
    derivs[0][k][1]=fac*(dz*derivs[0][k][0]+
			 2.0*dy*derivs[0][k-1][1]-
			 derivs[0][k-2][1]);
    
    derivs[1][k][0]=fac*(dx*derivs[0][k][0]+
			 2.0*dy*derivs[1][k-1][0]-
			 derivs[1][k-2][0]);
    
    derivs[k][1][0]=fac*(dy*derivs[k][0][0]+
			 2.0*dx*derivs[k-1][1][0]-
			 derivs[k-2][1][0]);
    
    derivs[k][0][1]=fac*(dz*derivs[k][0][0]+
			 2.0*dx*derivs[k-1][0][1]-
			 derivs[k-2][0][1]);
  }

  // derivatives for G(0,i,j), G(i,0,j), G(i,j,0), i,j >= 2

  for (int i=2; i<nterms-1; i++) {
    for (int j=2; j<nterms-i+1; j++) {
      ii = 1.0/i;
      derivs[i][j][0]=fac*(dx*(2.0-ii)*derivs[i-1][j][0]+
			   2.0*dy*derivs[i][j-1][0]-
			   (1-ii)*derivs[i-2][j][0]-
			   derivs[i][j-2][0]);
      
      derivs[i][0][j]=fac*(dx*(2.0-ii)*derivs[i-1][0][j]+
			   2.0*dz*derivs[i][0][j-1]-
			   (1-ii)*derivs[i-2][0][j]-
			   derivs[i][0][j-2]);
      
      derivs[0][i][j]=fac*(dy*(2.0-ii)*derivs[0][i-1][j]+
			   2.0*dz*derivs[0][i][j-1]-
			   (1-ii)*derivs[0][i-2][j]-
			   derivs[0][i][j-2]);   
    }
  }
  
  
  // 2 indices 1, other >= 1
  // deriv(1,1,1) is correct, but a little tricky!  
  //      b(1,1,1)=5.0*dz*fac*b(1,1,0)
  
  derivs[1][1][1]=fac*(dx*derivs[0][1][1] +
		      2.0*dy*derivs[1][0][1]+
		      2.0*dz*derivs[1][1][0]);
  
  for (int i=2; i<nterms-1; i++ ){
    derivs[1][1][i]=fac*(dx*derivs[0][1][i] + 
			 2.0*dy*derivs[1][0][i] +
			 2.0*dz*derivs[1][1][i-1] -
			 derivs[1][1][i-2]);

    derivs[1][i][1]=fac*(dx*derivs[0][i][1] +
			 2.0*dy*derivs[1][i-1][1] +
			 2.0*dz*derivs[1][i][0] -
			 derivs[1][i-2][1]);

    derivs[i][1][1]=fac*(dy*derivs[i][0][1] +
			 2.0*dx*derivs[i-1][1][1] +
			 2.0*dz*derivs[i][1][0] -
			 derivs[i-2][1][1]);
  }

  // 1 index 1, others >=2
  for (int i=2; i<nterms-2; i++) {
    for (int j=2; j<nterms-i+1; j++) {
      derivs[1][i][j]=fac*(dx*derivs[0][i][j] +
			   2.0*dy*derivs[1][i-1][j] +
			   2.0*dz*derivs[1][i][j-1] -
                           derivs[1][i-2][j] - 
			   derivs[1][i][j-2]);

      derivs[i][1][j]=fac*(dy*derivs[i][0][j] +
			   2.0*dx*derivs[i-1][1][j] +
			   2.0*dz*derivs[i][1][j-1] -
			   derivs[i-2][1][j] -
			   derivs[i][1][j-2]);

      derivs[i][j][1]=fac*(dz*derivs[i][j][0] +
			   2.0*dx*derivs[i-1][j][1] +
			   2.0*dy*derivs[i][j-1][1] -
                           derivs[i-2][j][1] -
			   derivs[i][j-2][1]);
    }
  }

  // all indices >=2
  for (int k=2; k<nterms-3; k++) {
    for (int j=2; j<nterms-1-k; j++) {
      for (int i=2; i<nterms-k-j+1; i++) {
	ii = 1.0/i;
	derivs[i][j][k]=fac*(2.0*dx*(1-0.5*ii)*derivs[i-1][j][k] +
			     2.0*dy*derivs[i][j-1][k] +
			     2.0*dz*derivs[i][j][k-1] - 
			     (1.0-ii)*derivs[i-2][j][k] - 
			     derivs[i][j-2][k] -
			     derivs[i][j][k-2]); 
      }
    }
  }

  // Add up the Taylor derivatives and moments to get the force
  for (int m=0; m<nterms; m++) {
    //    for (int k=0; k<nterms-m+1; k++)  {
    //for (int i=0; i<nterms-m-k+1; i++)  {
    for (int k=0; k<nterms-m; k++)  {
      for (int i=0; i<nterms-m-k; i++)  {
	*xfor+=moment[i][k][m]*derivs[i+1][k][m]*(i+1);
	*yfor+=moment[i][k][m]*derivs[i][k+1][m]*(k+1);
	*zfor+=moment[i][k][m]*derivs[i][k][m+1]*(m+1);

	// this gives the Greens function
	//*xfor+=moment[m*nterms*nterms+i*nterms+k]*derivs[i][k][m];
	//*yfor+=moment[m*nterms*nterms+i*nterms+k]*derivs[i][k][m];
	//*zfor+=moment[m*nterms*nterms+i*nterms+k]*derivs[i][k][m];
      }
    }
  }
  return;
} // end treetaylorforce



///////////////////////////////////////////////////////
// Function Name: shielded_treetaylorforce
// Usage: Compute shielded force for a cluster using far 
//        field taylor expansion
//
///////////////////////////////////////////////////////
// Assumptions:
//        none
//
///////////////////////////////////////////////////////
// Inputs:
//        nterms - number of terms in taylor expansion 
//        dx     - x_part - x_cluster_center
//        dy     - y_part - y_cluster_center
//        dz     - z_part - z_cluster_center
//        rs   - sqrt(dx^2+dy^2+dz^2+del)
//        moment - list of moments for given cluster
//        kappa  - screening parameter
//
///////////////////////////////////////////////////////
// Outputs: (By reference)
//         force   - vector containing (fx,fy,fz);
//
///////////////////////////////////////////////////////
// Functions Called:
//         none
//
///////////////////////////////////////////////////////
void shielded_treetaylorforce(int nterms, double dx, double dy, double dz,
			      double dist, double *xfor, double *yfor, 
			      double *zfor, double*** moment, double kappa)  {

  if (verbosity) {
    cout << "[shielded_treetaylorforce] called" << endl;
  }
  

  double ddx=2*dx;
  double ddy=2*dy;
  double ddz=2*dz;

  double kappax=kappa*dx;
  double kappay=kappa*dy;
  double kappaz=kappa*dz;
  
  double fac=1.0/dist;
  dist=sqrt(dist);


  double a[MAXTERM][MAXTERM][MAXTERM];
  double b[MAXTERM][MAXTERM][MAXTERM];


  // we should store this globally like pei jun
  double cf[MAXTERM];
  double cf1[MAXTERM];
  double cf2[MAXTERM];
  double cf3[MAXTERM];

  if (verbosity) {
    cout << "[shielded_treetaylorforce] : " 
	 << "computing vector of coefficients" << endl;
  }

  cf[0] = 1;
  for (int i=1;i<MAXTERM;i++) {
    cf[i] = 1.0+i;
    cf1[i] = 1.0/i;
    cf2[i] = 1.0-cf1[i]/2.0;
    cf3[i] = 1.0 - cf1[i];
  }


  if (verbosity) {
    cout << "[shielded_treetaylorforce] : " 
	 << "computing index 0,0,0" << endl;
  }
  
  // i=0, j=0, k=0
  b[0][0][0]=exp(-kappa*dist);
  a[0][0][0]=b[0][0][0]/dist;


  if (verbosity) {
    cout << "[shielded_treetaylorforce] : " 
	 << "computing two indices zero, other index greater than one" << endl;
  }

  // two indices are zero, an index is one
  b[1][0][0]=kappax*a[0][0][0];
  b[0][1][0]=kappay*a[0][0][0];
  b[0][0][1]=kappaz*a[0][0][0];

  a[1][0][0]=fac*dx*(a[0][0][0]+kappa*b[0][0][0]);
  a[0][1][0]=fac*dy*(a[0][0][0]+kappa*b[0][0][0]);
  a[0][0][1]=fac*dz*(a[0][0][0]+kappa*b[0][0][0]);

  // two indices are zero, the other index is greater than one
  for (int i=2; i<nterms+1; i++) {
    b[i][0][0]=cf1[i]*kappa*(dx*a[i-1][0][0]-a[i-2][0][0]);
    b[0][i][0]=cf1[i]*kappa*(dy*a[0][i-1][0]-a[0][i-2][0]);
    b[0][0][i]=cf1[i]*kappa*(dz*a[0][0][i-1]-a[0][0][i-2]);

    a[i][0][0]=fac*(ddx*cf2[i]*a[i-1][0][0]-cf3[i]*a[i-2][0][0] +
		  cf1[i]*kappa*(dx*b[i-1][0][0]-b[i-2][0][0]));
    a[0][i][0]=fac*(ddy*cf2[i]*a[0][i-1][0]-cf3[i]*a[0][i-2][0] +
		  cf1[i]*kappa*(dy*b[0][i-1][0]-b[0][i-2][0]));
    a[0][0][i]=fac*(ddz*cf2[i]*a[0][0][i-1]-cf3[i]*a[0][0][i-2] +
		  cf1[i]*kappa*(dz*b[0][0][i-1]-b[0][0][i-2]));
  }

  if (verbosity) {
    cout << "[shielded_treetaylorforce] : " 
	 << "computing one index zero, one index one, and other index > 1" << endl;
  }

  // 1 index 0, 1 index 1, other >=1
  b[1][1][0]=kappax*a[0][1][0];
  b[1][0][1]=kappax*a[0][0][1];
  b[0][1][1]=kappay*a[0][0][1];

  a[1][1][0]=fac*(dx*a[0][1][0]+ddy*a[1][0][0]+kappax*b[0][1][0]);
  a[1][0][1]=fac*(dx*a[0][0][1]+ddz*a[1][0][0]+kappax*b[0][0][1]);
  a[0][1][1]=fac*(dy*a[0][0][1]+ddz*a[0][1][0]+kappay*b[0][0][1]);

  for (int i = 2; i<nterms; i++) {
    b[1][0][i]=kappax*a[0][0][i];
    b[0][1][i]=kappay*a[0][0][i];
    b[0][i][1]=kappaz*a[0][i][0];
    b[1][i][0]=kappax*a[0][i][0];
    b[i][1][0]=kappay*a[i][0][0];
    b[i][0][1]=kappaz*a[i][0][0];

    a[1][0][i]=fac*(dx*a[0][0][i]+ddz*a[1][0][i-1]-a[1][0][i-2]+
		  kappax*b[0][0][i]); 
    a[0][1][i]=fac*(dy*a[0][0][i]+ddz*a[0][1][i-1]-a[0][1][i-2]+
		  kappay*b[0][0][i]);
    a[0][i][1]=fac*(dz*a[0][i][0]+ddy*a[0][i-1][1]-a[0][i-2][1]+
		  kappaz*b[0][i][0]);
    a[1][i][0]=fac*(dx*a[0][i][0]+ddy*a[1][i-1][0]-a[1][i-2][0]+
		  kappax*b[0][i][0]);
    a[i][1][0]=fac*(dy*a[i][0][0]+ddx*a[i-1][1][0]-a[i-2][1][0]+
		  kappay*b[i][0][0]);
    a[i][0][1]=fac*(dz*a[i][0][0]+ddx*a[i-1][0][1]-a[i-2][0][1]+
		  kappaz*b[i][0][0]);         
  }

  if (verbosity) {
    cout << "[shielded_treetaylorforce] : " 
	 << "computing one index zero, other indices > 2" << endl;
  }

  // 1 index 0, others >= 2

  for (int i=2; i<nterms-1; i++) {
    for (int j=2; j<nterms-i+1; j++) {
      b[i][j][0]=cf1[i]*kappa*(dx*a[i-1][j][0]-a[i-2][j][0]);
      b[i][0][j]=cf1[i]*kappa*(dx*a[i-1][0][j]-a[i-2][0][j]);
      b[0][i][j]=cf1[i]*kappa*(dy*a[0][i-1][j]-a[0][i-2][j]);

      a[i][j][0]=fac*(ddx*cf2[i]*a[i-1][j][0]+ddy*a[i][j-1][0]
                    -cf3[i]*a[i-2][j][0]-a[i][j-2][0]+
		    cf1[i]*kappa*(dx*b[i-1][j][0]-b[i-2][j][0]));
      a[i][0][j]=fac*(ddx*cf2[i]*a[i-1][0][j]+ddz*a[i][0][j-1]
                    -cf3[i]*a[i-2][0][j]-a[i][0][j-2]+
		    cf1[i]*kappa*(dx*b[i-1][0][j]-b[i-2][0][j]));
      a[0][i][j]=fac*(ddy*cf2[i]*a[0][i-1][j]+ddz*a[0][i][j-1]
                    -cf3[i]*a[0][i-2][j]-a[0][i][j-2]+
		    cf1[i]*kappa*(dy*b[0][i-1][j]-b[0][i-2][j]));
    }
  }

  if (verbosity) {
    cout << "[shielded_treetaylorforce] : " 
	 << "computing two indices one, other index greater than one" << endl;
  }

  //2 indices 1, other >= 1

  b[1][1][1]=kappax*a[0][1][1];
  a[1][1][1]=fac*(dx*a[0][1][1]+ddy*a[1][0][1]+ddz*a[1][1][0]+
		kappax*b[0][1][1]);

  for (int i=2; i<nterms-1; i++){
    b[1][1][i]=kappax*a[0][1][i];
    b[1][i][1]=kappax*a[0][i][1];
    b[i][1][1]=kappay*a[i][0][1];

    a[1][1][i]=fac*(dx*a[0][1][i]+ddy*a[1][0][i]+ddz*a[1][1][i-1]
		  -a[1][1][i-2]+kappax*b[0][1][i]);
    a[1][i][1]=fac*(dx*a[0][i][1]+ddy*a[1][i-1][1]+ddz*a[1][i][0]
		  -a[1][i-2][1]+kappax*b[0][i][1]);
    a[i][1][1]=fac*(dy*a[i][0][1]+ddx*a[i-1][1][1]+ddz*a[i][1][0]
		  -a[i-2][1][1]+kappay*b[i][0][1]);
  }

  if (verbosity) {
    cout << "[shielded_treetaylorforce] : " 
	 << "computing one index one, other two indices greater than one" << endl;
  }

  // 1 index 1, others >=2

  for (int i=2; i<nterms-2; i++) {
    for (int j=2; j<nterms-i+1; j++) {
      b[1][i][j]=kappax*a[0][i][j];
      b[i][1][j]=kappay*a[i][0][j];
      b[i][j][1]=kappaz*a[i][j][0];

      a[1][i][j]=fac*(dx*a[0][i][j]+ddy*a[1][i-1][j]+ddz*a[1][i][j-1]
		    -a[1][i-2][j]-a[1][i][j-2]+kappax*b[0][i][j]);
      a[i][1][j]=fac*(dy*a[i][0][j]+ddx*a[i-1][1][j]+ddz*a[i][1][j-1]
		    -a[i-2][1][j]-a[i][1][j-2]+kappay*b[i][0][j]);
      a[i][j][1]=fac*(dz*a[i][j][0]+ddx*a[i-1][j][1]+ddy*a[i][j-1][1]
		    -a[i-2][j][1]-a[i][j-2][1]+kappaz*b[i][j][0]);

    }
  }

  if (verbosity) {
    cout << "[shielded_treetaylorforce] : " 
	 << "computing all indices >= 1" << endl;
  }
  
  // all indices >=2

  for (int k=2; k<nterms-3; k++) {
    for (int j=2; j<nterms-1-k; j++) {
      for (int i=2; i<nterms-k-j+1; i++) {
	b[i][j][k]=cf1[i]*kappa*(dx*a[i-1][j][k]-a[i-2][j][k]);

	a[i][j][k]=fac*(ddx*cf2[i]*a[i-1][j][k]+ddy*a[i][j-1][k]
                      +ddz*a[i][j][k-1]-cf3[i]*a[i-2][j][k]
                      -a[i][j-2][k]-a[i][j][k-2]+
		      cf1[i]*kappa*(dx*b[i-1][j][k]-b[i-2][j][k]));
      }
    }
  }



  // Add up the Taylor derivatives and moments to get the force
  for (int m=0; m<nterms; m++) {
    for (int k=0; k<nterms-m; k++)  {
      for (int i=0; i<nterms-m-k; i++)  {
	//*xfor+=cf[i]*moment[i][k][m]*a[i+1][k][m];
	//*yfor+=cf[k]*moment[i][k][m]*a[i][k+1][m];
	//*zfor+=cf[m]*moment[i][k][m]*a[i][k][m+1];
	*xfor-=cf[k]*moment[i][k][m]*a[i][k+1][m];
	*yfor+=cf[i]*moment[i][k][m]*a[i+1][k][m];
	//*zfor+=cf[m]*moment[i][k][m]*a[i][k][m+1];
	*zfor = 0.0;
      }
    }
  }
  return;
} // end shielded_treetaylorforce


///////////////////////////////////////////////////////
// Function Name: shielded_treetaylorpot
// Usage: Compute shielded potential for a cluster using far 
//        field taylor expansion
//
///////////////////////////////////////////////////////
// Assumptions:
//        none
//
///////////////////////////////////////////////////////
// Inputs:
//        nterms - number of terms in taylor expansion 
//        dx     - x_part - x_cluster_center
//        dy     - y_part - y_cluster_center
//        dz     - z_part - z_cluster_center
//        rs   - sqrt(dx^2+dy^2+dz^2+del)
//        moment - list of moments for given cluster
//        kappa  - screening parameter
//
///////////////////////////////////////////////////////
// Outputs: (By reference)
//         pot - potential
//
///////////////////////////////////////////////////////
// Functions Called:
//         none
//
///////////////////////////////////////////////////////
void shielded_treetaylorpot(int nterms, double dx, double dy, double dz,
			    double dist, double *pot,
			    double*** moment, double kappa)  {

  if (verbosity) {
    cout << "[shielded_treetaylorpot] called" << endl;
  }
  

  double ddx=2*dx;
  double ddy=2*dy;
  double ddz=2*dz;

  double kappax=kappa*dx;
  double kappay=kappa*dy;
  double kappaz=kappa*dz;
  
  double fac=1.0/dist;
  dist=sqrt(dist);


  double a[MAXTERM][MAXTERM][MAXTERM];
  double b[MAXTERM][MAXTERM][MAXTERM];


  // we should store this globally like pei jun
  double cf[MAXTERM];
  double cf1[MAXTERM];
  double cf2[MAXTERM];
  double cf3[MAXTERM];

  if (verbosity) {
    cout << "[shielded_treetaylorpot] : " 
	 << "computing vector of coefficients" << endl;
  }

  cf[0] = 1;
  for (int i=1;i<MAXTERM;i++) {
    cf[i] = 1.0+i;
    cf1[i] = 1.0/i;
    cf2[i] = 1.0-cf1[i]/2.0;
    cf3[i] = 1.0 - cf1[i];
  }


  if (verbosity) {
    cout << "[shielded_treetaylorpot] : " 
	 << "computing index 0,0,0" << endl;
  }
  
  // i=0, j=0, k=0
  b[0][0][0]=exp(-kappa*dist);
  a[0][0][0]=b[0][0][0]/dist;


  if (verbosity) {
    cout << "[shielded_treetaylorpot] : " 
	 << "computing two indices zero, other index greater than one" << endl;
  }

  // two indices are zero, an index is one
  b[1][0][0]=kappax*a[0][0][0];
  b[0][1][0]=kappay*a[0][0][0];
  b[0][0][1]=kappaz*a[0][0][0];

  a[1][0][0]=fac*dx*(a[0][0][0]+kappa*b[0][0][0]);
  a[0][1][0]=fac*dy*(a[0][0][0]+kappa*b[0][0][0]);
  a[0][0][1]=fac*dz*(a[0][0][0]+kappa*b[0][0][0]);

  // two indices are zero, the other index is greater than one
  for (int i=2; i<nterms+1; i++) {
    b[i][0][0]=cf1[i]*kappa*(dx*a[i-1][0][0]-a[i-2][0][0]);
    b[0][i][0]=cf1[i]*kappa*(dy*a[0][i-1][0]-a[0][i-2][0]);
    b[0][0][i]=cf1[i]*kappa*(dz*a[0][0][i-1]-a[0][0][i-2]);

    a[i][0][0]=fac*(ddx*cf2[i]*a[i-1][0][0]-cf3[i]*a[i-2][0][0] +
		  cf1[i]*kappa*(dx*b[i-1][0][0]-b[i-2][0][0]));
    a[0][i][0]=fac*(ddy*cf2[i]*a[0][i-1][0]-cf3[i]*a[0][i-2][0] +
		  cf1[i]*kappa*(dy*b[0][i-1][0]-b[0][i-2][0]));
    a[0][0][i]=fac*(ddz*cf2[i]*a[0][0][i-1]-cf3[i]*a[0][0][i-2] +
		  cf1[i]*kappa*(dz*b[0][0][i-1]-b[0][0][i-2]));
  }

  if (verbosity) {
    cout << "[shielded_treetaylorpot] : " 
	 << "computing one index zero, one index one, and other index > 1" << endl;
  }

  // 1 index 0, 1 index 1, other >=1
  b[1][1][0]=kappax*a[0][1][0];
  b[1][0][1]=kappax*a[0][0][1];
  b[0][1][1]=kappay*a[0][0][1];

  a[1][1][0]=fac*(dx*a[0][1][0]+ddy*a[1][0][0]+kappax*b[0][1][0]);
  a[1][0][1]=fac*(dx*a[0][0][1]+ddz*a[1][0][0]+kappax*b[0][0][1]);
  a[0][1][1]=fac*(dy*a[0][0][1]+ddz*a[0][1][0]+kappay*b[0][0][1]);

  for (int i = 2; i<nterms; i++) {
    b[1][0][i]=kappax*a[0][0][i];
    b[0][1][i]=kappay*a[0][0][i];
    b[0][i][1]=kappaz*a[0][i][0];
    b[1][i][0]=kappax*a[0][i][0];
    b[i][1][0]=kappay*a[i][0][0];
    b[i][0][1]=kappaz*a[i][0][0];

    a[1][0][i]=fac*(dx*a[0][0][i]+ddz*a[1][0][i-1]-a[1][0][i-2]+
		  kappax*b[0][0][i]); 
    a[0][1][i]=fac*(dy*a[0][0][i]+ddz*a[0][1][i-1]-a[0][1][i-2]+
		  kappay*b[0][0][i]);
    a[0][i][1]=fac*(dz*a[0][i][0]+ddy*a[0][i-1][1]-a[0][i-2][1]+
		  kappaz*b[0][i][0]);
    a[1][i][0]=fac*(dx*a[0][i][0]+ddy*a[1][i-1][0]-a[1][i-2][0]+
		  kappax*b[0][i][0]);
    a[i][1][0]=fac*(dy*a[i][0][0]+ddx*a[i-1][1][0]-a[i-2][1][0]+
		  kappay*b[i][0][0]);
    a[i][0][1]=fac*(dz*a[i][0][0]+ddx*a[i-1][0][1]-a[i-2][0][1]+
		  kappaz*b[i][0][0]);         
  }

  if (verbosity) {
    cout << "[shielded_treetaylorpot] : " 
	 << "computing one index zero, other indices > 2" << endl;
  }

  // 1 index 0, others >= 2

  for (int i=2; i<nterms-1; i++) {
    for (int j=2; j<nterms-i+1; j++) {
      b[i][j][0]=cf1[i]*kappa*(dx*a[i-1][j][0]-a[i-2][j][0]);
      b[i][0][j]=cf1[i]*kappa*(dx*a[i-1][0][j]-a[i-2][0][j]);
      b[0][i][j]=cf1[i]*kappa*(dy*a[0][i-1][j]-a[0][i-2][j]);

      a[i][j][0]=fac*(ddx*cf2[i]*a[i-1][j][0]+ddy*a[i][j-1][0]
                    -cf3[i]*a[i-2][j][0]-a[i][j-2][0]+
		    cf1[i]*kappa*(dx*b[i-1][j][0]-b[i-2][j][0]));
      a[i][0][j]=fac*(ddx*cf2[i]*a[i-1][0][j]+ddz*a[i][0][j-1]
                    -cf3[i]*a[i-2][0][j]-a[i][0][j-2]+
		    cf1[i]*kappa*(dx*b[i-1][0][j]-b[i-2][0][j]));
      a[0][i][j]=fac*(ddy*cf2[i]*a[0][i-1][j]+ddz*a[0][i][j-1]
                    -cf3[i]*a[0][i-2][j]-a[0][i][j-2]+
		    cf1[i]*kappa*(dy*b[0][i-1][j]-b[0][i-2][j]));
    }
  }

  if (verbosity) {
    cout << "[shielded_treetaylorpot] : " 
	 << "computing two indices one, other index greater than one" << endl;
  }

  //2 indices 1, other >= 1

  b[1][1][1]=kappax*a[0][1][1];
  a[1][1][1]=fac*(dx*a[0][1][1]+ddy*a[1][0][1]+ddz*a[1][1][0]+
		kappax*b[0][1][1]);

  for (int i=2; i<nterms-1; i++){
    b[1][1][i]=kappax*a[0][1][i];
    b[1][i][1]=kappax*a[0][i][1];
    b[i][1][1]=kappay*a[i][0][1];

    a[1][1][i]=fac*(dx*a[0][1][i]+ddy*a[1][0][i]+ddz*a[1][1][i-1]
		  -a[1][1][i-2]+kappax*b[0][1][i]);
    a[1][i][1]=fac*(dx*a[0][i][1]+ddy*a[1][i-1][1]+ddz*a[1][i][0]
		  -a[1][i-2][1]+kappax*b[0][i][1]);
    a[i][1][1]=fac*(dy*a[i][0][1]+ddx*a[i-1][1][1]+ddz*a[i][1][0]
		  -a[i-2][1][1]+kappay*b[i][0][1]);
  }

  if (verbosity) {
    cout << "[shielded_treetaylorpot] : " 
	 << "computing one index one, other two indices greater than one" << endl;
  }

  // 1 index 1, others >=2

  for (int i=2; i<nterms-2; i++) {
    for (int j=2; j<nterms-i+1; j++) {
      b[1][i][j]=kappax*a[0][i][j];
      b[i][1][j]=kappay*a[i][0][j];
      b[i][j][1]=kappaz*a[i][j][0];

      a[1][i][j]=fac*(dx*a[0][i][j]+ddy*a[1][i-1][j]+ddz*a[1][i][j-1]
		    -a[1][i-2][j]-a[1][i][j-2]+kappax*b[0][i][j]);
      a[i][1][j]=fac*(dy*a[i][0][j]+ddx*a[i-1][1][j]+ddz*a[i][1][j-1]
		    -a[i-2][1][j]-a[i][1][j-2]+kappay*b[i][0][j]);
      a[i][j][1]=fac*(dz*a[i][j][0]+ddx*a[i-1][j][1]+ddy*a[i][j-1][1]
		    -a[i-2][j][1]-a[i][j-2][1]+kappaz*b[i][j][0]);

    }
  }

  if (verbosity) {
    cout << "[shielded_treetaylorpot] : " 
	 << "computing all indices >= 1" << endl;
  }
  
  // all indices >=2

  for (int k=2; k<nterms-3; k++) {
    for (int j=2; j<nterms-1-k; j++) {
      for (int i=2; i<nterms-k-j+1; i++) {
	b[i][j][k]=cf1[i]*kappa*(dx*a[i-1][j][k]-a[i-2][j][k]);

	a[i][j][k]=fac*(ddx*cf2[i]*a[i-1][j][k]+ddy*a[i][j-1][k]
                      +ddz*a[i][j][k-1]-cf3[i]*a[i-2][j][k]
                      -a[i][j-2][k]-a[i][j][k-2]+
		      cf1[i]*kappa*(dx*b[i-1][j][k]-b[i-2][j][k]));
      }
    }
  }



  // Add up the Taylor derivatives and moments to get the potential
  for (int m=0; m<nterms; m++) {
    for (int k=0; k<nterms-m; k++)  {
      for (int i=0; i<nterms-m-k; i++)  {
	*pot-=moment[i][k][m]*a[i][k][m];
      }
    }
  }
  return;
} // end shielded_treetaylorpot


///////////////////////////////////////////////////////
// Function Name: treeforce
// Usage: Compute force at location x,y,z due to 
//        all clusters of particles.  Algorithum
//        reverts to direct sum when clusters are 
//        very close to point x,y,z . 
//        (Recursive Function)
//
///////////////////////////////////////////////////////
// Assumptions:
//         eps, set in constant.h
//
///////////////////////////////////////////////////////
// Inputs:
//         ind    - index of tree cluster 
//                  ( when called index = 0 )
//         p      - particle index that we are computing 
//                  force on.  If -1, x,y is not a particle
//         nterms - number of terms in taylor expantion
//         x      - x location
//         y      - y location
//         z      - z location
//         acc    - acceptance criterion for using 
//                  cluster approximation
//         tree   - oct tree sorting particle locations
//         part   - list of all particles
//
///////////////////////////////////////////////////////
// Outputs: (By Reference)
//         force   - (fx,fy,fz)
//
///////////////////////////////////////////////////////
// Functions Called:
//         none
//
///////////////////////////////////////////////////////
double treeforce(int ind, int p, int nterms, 
		 double x, double y, double z,
		 double acc, double *xfor, double *yfor, double *zfor,
		 TREE* tree, PARTICLE* part)  {
  int go=1;
  int accept;

  double dx=x-tree[ind].xmid;
  double dy=y-tree[ind].ymid;
  double dz=z-tree[ind].zmid;
  double rs=dx*dx+dy*dy+dz*dz + DEL;  // Determine the square of the radius
  double r;

  // Determine whether or not to accept the cluster
  accept=(tree[ind].sqradius < acc*(rs-DEL));

  if (verbosity) {
    cout << "Checking tree " << ind << endl;
    cout << "cluster center = (" << tree[ind].xmid 
	 << ", " << tree[ind].ymid
	 << ", " << tree[ind].zmid << "]" << endl;
    cout << "cluster radius squared = " << tree[ind].sqradius << endl;
    cout << "distance squared from cluster center to r = " << rs <<endl;
    if (accept) {
      cout << "... tree accepted" << endl; 
    } else {
      cout << "... tree NOT accepted" << endl; 
    }
  }


  // On acceptance, use the Taylor expansion
  if (accept)  {
    if (verbosity) {
      cout << "... [treeforce]: using tree[" << ind 
	   << "] to approximate force for particle[" << p <<"]"  << endl;
    }
    treetaylorforce(nterms, dx, dy, dz, rs, xfor, yfor, zfor,
		    tree[ind].moments);
    if (verbosity) {
      cout << " ... [treeforce] ... finished" << endl;
    }
    return(0);
  }

  // Cluster was not accepted
  else  {
    if (verbosity){
      cout << "... [treeforce]: checking children clusters of tree[" <<
	ind << "]:" << endl;
      for (int i = 0; i<NUMCHILD; i++) {
	cout << "... [treeforce]: ... " << tree[ind].children[i] << endl;
      }
      cout << "... to approximate force for particle[" << p <<"]"  << endl;
    }
    // Check for child clusters
    for (int i=0; i<NUMCHILD; i++)  {
      if (tree[ind].children[i] != -1)  {
        go=0;
        // Call the recursion
        treeforce(tree[ind].children[i], p, nterms, x, y, z, acc, 
		  xfor, yfor, zfor, tree, part);
      }
    }

    // There are no child clusters, so do direct sum
    if (go)  {
      if (verbosity) {
	cout << "... [treeforce]: using direct sum from tree[" << ind 
	     << "] to approximate force for particle[" << p <<"]"  << endl;
      }
      int j;
      double temp=0.0;

      // Do direct summation on all particles in the cell
      if (p == -1)  {
        for (int j1=0; j1<tree[ind].nummemb; j1++)  {
          j=tree[ind].members[j1];
          dx=x-part[j].x;
          dy=y-part[j].y;
	  dz=z-part[j].z;
	  rs = dx*dx + dy*dy + dz*dz + DEL;
	  r = sqrt(rs);
          if (rs > eps)  {
            temp=part[j].tot_charge/(rs*r);
            *xfor+=dx*temp;
            *yfor+=dy*temp;
	    *zfor+=dz*temp;
          }
        }
      }
      else  {
        for (int j1=0; j1<tree[ind].nummemb; j1++)  {
          j=tree[ind].members[j1];
          if (j != p)  {
            dx=x-part[j].x;
            dy=y-part[j].y;
	    dz=z-part[j].z;
            rs=dx*dx+dy*dy+dz*dz + DEL;
	    r = sqrt(rs);
            //if (temp > eps)  {
	    temp=part[j].tot_charge/(rs*r);
	    *xfor+=dx*temp;
	    *yfor+=dy*temp;
	    *zfor+=dz*temp;
	    //} // if sufficiently far enough away from (x,y,z)
          } 
        }// end add contribution from particle tree[ind].member[j1]
      }
    }
  }

  return(0);
} // end treeforce


///////////////////////////////////////////////////////
// Function Name: shielded_treeforce
// Usage: Compute shielded force at location x,y,z due to 
//        all clusters of particles.  Algorithum
//        reverts to direct sum when clusters are 
//        very close to point x,y,z . 
//        (Recursive Function)
//
///////////////////////////////////////////////////////
// Assumptions:
//         eps, set in constant.h
//
///////////////////////////////////////////////////////
// Inputs:
//         ind    - index of tree cluster 
//                  ( when called index = 0 )
//         p      - particle index that we are computing 
//                  force on.  If -1, x,y is not a particle
//         nterms - number of terms in taylor expantion
//         x      - x location
//         y      - y location
//         z      - z location
//         alpha  - shielding parameter
//         acc    - acceptance criterion for using 
//                  cluster approximation
//         tree   - oct tree sorting particle locations
//         part   - list of all particles
//
///////////////////////////////////////////////////////
// Outputs: (By Reference)
//         force   - (fx,fy,fz)
//
///////////////////////////////////////////////////////
// Functions Called:
//         shielded_treetaylorforce
//
///////////////////////////////////////////////////////
double shielded_treeforce(int ind, int p, int nterms, 
			  double x, double y, double z, double alpha,
			  double acc, double *xfor, double *yfor, 
			  double *zfor, TREE* tree, PARTICLE* part)  {
  int go=1;
  int accept;

  double dx=x-tree[ind].xmid;
  double dy=y-tree[ind].ymid;
  double dz=z-tree[ind].zmid;
  // *** NOTE, removing regularization
  double rs=dx*dx+dy*dy+dz*dz;  // Determine the square of the radius
  double r;

  // Determine whether or not to accept the cluster
  accept=(tree[ind].sqradius < acc*(rs-DEL));

  if (verbosity) {
    cout << "[shielded_treeforce] : " << endl;
    cout << "Checking tree " << ind << endl;
    cout << "cluster center = (" << tree[ind].xmid 
	 << ", " << tree[ind].ymid
	 << ", " << tree[ind].zmid << "]" << endl;
    cout << "cluster radius squared = " << tree[ind].sqradius << endl;
    cout << "distance squared from cluster center to r = " << rs <<endl;
    if (accept) {
      cout << "... tree accepted" << endl; 
    } else {
      cout << "... tree NOT accepted" << endl; 
    }
  }


  // On acceptance, use the Taylor expansion
  if (accept)  {
    if (verbosity) {
      cout << "... [shieldedtreeforce]: using tree[" << ind 
	   << "] to approximate force for particle[" << p <<"]"  << endl;
      cout << "nterms = " << nterms << endl;
      cout << "dx = " << dx << endl;
      cout << "dy = " << dy << endl;
      cout << "dz = " << dz << endl;
      cout << "rs = " << rs << endl;
      cout << "alpha = " << alpha << endl;
      cout << "xfor = " << *xfor << endl;
      cout << "yfor = " << *yfor << endl;
      cout << "zfor = " << *zfor << endl;
      cout << "tree[" << ind << "].moments[0][0][0] = " 
	   << tree[ind].moments[0][0][0] << endl;
      cout << "tree[" << ind << "].moments[2][2][2] = " 
	   << tree[ind].moments[2][2][2] << endl;	
    }

    //treetaylorforce(nterms, dx, dy, dz, rs, xfor, yfor, zfor,
    //    tree[ind].moments);

    shielded_treetaylorforce(nterms, dx, dy, dz, rs, xfor, yfor, zfor,
    		     tree[ind].moments,alpha);

    if (verbosity) {
      cout << " ... [shieldedtreeforce] ... finished" << endl;
    }
    return(0);
  }

  // Cluster was not accepted
  else  {
    if (verbosity){
      cout << "... [shielded_treeforce]: checking children clusters of tree[" <<
	ind << "]:" << endl;
      for (int i = 0; i<NUMCHILD; i++) {
	cout << "... [shielded_treeforce]: ... " << tree[ind].children[i] << endl;
      }
      cout << "... to approximate force for particle[" << p <<"]"  << endl;
    }
    // Check for child clusters
    for (int i=0; i<NUMCHILD; i++)  {
      if (tree[ind].children[i] != -1)  {
        go=0;
        // Call the recursion
        shielded_treeforce(tree[ind].children[i],p,nterms,x,y,z,alpha,acc, 
			   xfor, yfor, zfor, tree, part);
      }
    }

    // There are no child clusters, so do direct sum
    if (go)  {
      if (verbosity) {
	cout << "... [shielded_treeforce]: using direct sum from tree[" << ind 
	     << "] to approximate force for particle[" << p <<"]"  << endl;
      }
      int j;
      double temp=0.0;
      double expr;

      // Do direct summation on all particles in the cell
      if (p == -1)  {
        for (int j1=0; j1<tree[ind].nummemb; j1++)  {
          j=tree[ind].members[j1];
          dx=x-part[j].x;
          dy=y-part[j].y;
	  dz=z-part[j].z;
	  // *** NOTE: removing regularization
	  rs = dx*dx + dy*dy + dz*dz;
	  r = sqrt(rs);
	  expr = exp(-alpha*r)/rs;
          
          if (rs > eps)  {
	    temp=part[j].tot_charge*expr*(alpha + 1/r);
            *xfor -= dy*temp;
            *yfor += dx*temp;
	    *zfor = 0.0;
          }
        }
      }
      else  {
        for (int j1=0; j1<tree[ind].nummemb; j1++)  {
          j=tree[ind].members[j1];
          if (j != p)  {
            dx=x-part[j].x;
            dy=y-part[j].y;
	    dz=z-part[j].z;
	  // *** NOTE: removing regularization
	    rs = dx*dx + dy*dy + dz*dz;
	    r = sqrt(rs);
          
            //if (temp > eps)  {
	    //temp=part[j].tot_charge*expr*(alpha + 1/r);
	    temp=part[j].tot_charge*exp(-alpha*r)/rs*(alpha + 1/r);
	    *xfor -=dy*temp;
	    *yfor +=dx*temp;
	    *zfor = 0.0;
	    //} // if sufficiently far enough away from (x,y,z)
          } 
        }// end add contribution from particle tree[ind].member[j1]
      }
    }
  }

  return(0);
} // end shielded_treeforce



///////////////////////////////////////////////////////
// Function Name: shielded_treepot
// Usage: Compute shielded potential at location x,y,z due to 
//        all clusters of particles.  Algorithum
//        reverts to direct sum when clusters are 
//        very close to point x,y,z . 
//        (Recursive Function)
//
///////////////////////////////////////////////////////
// Assumptions:
//
///////////////////////////////////////////////////////
// Inputs:
//         ind    - index of tree cluster 
//                  ( when called index = 0 )
//         p      - particle index that we are computing 
//                  force on.  If -1, x,y is not a particle
//         nterms - number of terms in taylor expantion
//         x      - x location
//         y      - y location
//         z      - z location
//         alpha  - shielding parameter
//         acc    - acceptance criterion for using 
//                  cluster approximation
//         tree   - oct tree sorting particle locations
//         part   - list of all particles
//
///////////////////////////////////////////////////////
// Outputs: (By Reference)
//         potential
//
///////////////////////////////////////////////////////
// Functions Called:
//         shielded_treetaylorpotential
//
///////////////////////////////////////////////////////
double shielded_treepot(int ind, int p, int nterms, 
			double x, double y, double z, double alpha,
			double acc, double *pot,
			TREE* tree, PARTICLE* part,
			int Nx, int Ny, int Nz,
			double dlength, double dheight, double ddepth)  {
  int go=1;
  int accept;

  if (verbosity) {
    cout << "[shielded_treepot]: input parameters" << endl;
    cout << "... checking tree[" << ind << "]" << endl;
    cout << "... checking particle[" << p << "] at location (" 
	 << x << ", " << y << ", " << z << ")" << endl;
    cout << "... alpha = " << alpha << endl;
    cout << "... acc = " << acc << endl;
  }
  double dx=x-tree[ind].xmid;
  double dy=y-tree[ind].ymid;
  double dz=z-tree[ind].zmid;
  // *** NOTE, removing regularization
  double rs=dx*dx+dy*dy+dz*dz;  // Determine the square of the radius
  double r;

  // Determine whether or not to accept the cluster
  accept=(tree[ind].sqradius < acc*(rs-DEL));

  if (verbosity) {
    cout << "[shielded_treepot] : " << endl;
    cout << "Checking tree " << ind << endl;
    cout << "cluster center = (" << tree[ind].xmid 
	 << ", " << tree[ind].ymid
	 << ", " << tree[ind].zmid << "]" << endl;
    cout << "cluster radius squared = " << tree[ind].sqradius << endl;
    cout << "distance squared from cluster center to r = " << rs <<endl;
    if (accept) {
      cout << "... tree accepted" << endl; 
    } else {
      cout << "... tree NOT accepted" << endl; 
    }
  }


  // On acceptance, use the Taylor expansion
  if (accept)  {
    if (verbosity) {
      cout << "... [shieldedtreepot]: using tree[" << ind 
	   << "] to approximate pot for particle[" << p <<"]"  << endl;
      cout << "nterms = " << nterms << endl;
      cout << "dx = " << dx << endl;
      cout << "dy = " << dy << endl;
      cout << "dz = " << dz << endl;
      cout << "rs = " << rs << endl;
      cout << "alpha = " << alpha << endl;
      cout << "potential = " << *pot << endl;
      cout << "tree[" << ind << "].moments[0][0][0] = " 
	   << tree[ind].moments[0][0][0] << endl;
      cout << "tree[" << ind << "].moments[2][2][2] = " 
	   << tree[ind].moments[2][2][2] << endl;	
    }

    shielded_treetaylorpot(nterms, dx, dy, dz, rs, pot,
			   tree[ind].moments,alpha);

    if (verbosity) {
      cout << " ... [shieldedtreepot] ... finished" << endl;
    }
    return(0);
  }

  // Cluster was not accepted
  else  {
    if (verbosity){
      cout << "... [shielded_treepot]: checking children clusters of tree[" <<
	ind << "]:" << endl;
      for (int i = 0; i<NUMCHILD; i++) {
	cout << "... [shielded_treepot]: ... " << tree[ind].children[i] << endl;
      }
      cout << "... to approximate potential for particle[" << p <<"]"  << endl;
    }
    // Check for child clusters
    for (int i=0; i<NUMCHILD; i++)  {
      if (tree[ind].children[i] != -1)  {
        go=0;
        // Call the recursion
        shielded_treepot(tree[ind].children[i],p,nterms,x,y,z,alpha,acc, 
			 pot, tree, part,Nx,Ny,Nz,dlength,dheight,ddepth);
      }
    }

    // There are no child clusters, so do direct sum
    if (go)  {
      if (verbosity) {
	cout << "... [shielded_treepot]: using direct sum from tree[" << ind 
	     << "] to approximate potential for particle[" << p <<"]"  << endl;
      }
      int j;
      double temp=0.0;
      double expr;

      // Do direct summation on all particles in the cell
      if (p == -1)  {
	for (int j1=0; j1<tree[ind].nummemb; j1++)  {
	  j=tree[ind].members[j1];
	  dx=x-part[j].x;
	  dy=y-part[j].y;
	  dz=z-part[j].z;
	  // *** NOTE: removing regularization
	  rs = dx*dx + dy*dy + dz*dz;
	  r = sqrt(rs);
	  
	  if (rs > eps)  {
	    *pot += part[j].tot_charge*exp(-alpha*r)/r;
	  }
	}
      }
      else  {
	if (verbosity) {
	  cout << "[shielded_treepot]:" << endl;
	  cout << "...location of particle[" << p << "] = (" 
	       << x << ", " << y << ", " << z << ")" << endl;
	  cout << "... potential for particle[" << p 
	       << "], contributions from" << endl;
	}
	//#pragma omp parallel 
	//{
	//#pragma omp for
	for (int j1=0; j1<tree[ind].nummemb; j1++)  {
	  j=tree[ind].members[j1];
	  if (j == p) {
	    double R;
	    R = min(2.0*dlength/Nx,2.0*dheight/Ny);
	    R = min(R,2.0*ddepth/Nz);
	    R = R/2.0;
	    double intG;
	    // alpha = 100, Nx = 100
	    //intG = -0.000432626092185;

	    // alpha = 100, Nx = 400
	    intG = -0.0000483811229725;

	    //intG =  4*pi*1.0/alpha/alpha*(-1+exp(-alpha*R)+alpha*R*exp(-alpha*R));
	    /*
	    cout << "alpha = " << alpha 
		 << "R = " << R 
		 << ", ind = " << j 
		 << ", intG = " << intG 
		 << ", rhs = " << part[j].tot_charge/part[j].weight << endl;
	    */
	    // analytically integrate out singularity
	    *pot += part[j].tot_charge*intG/part[j].weight;
	  } else {
	    if (verbosity) {
	      cout << "...... particle[" << j << "] located at (" 
		   << part[j].x << ", " 
		   << part[j].y << ", " 
		   << part[j].z << ")" << endl;
	    }
	    dx=x-part[j].x;
	    dy=y-part[j].y;
	    dz=z-part[j].z;
	    // *** NOTE: removing regularization
	    rs = dx*dx + dy*dy + dz*dz;
	    r = sqrt(rs);
	    if (verbosity) {
	      cout << "......... (dx,dy,dz) = ( " << dx << ", " << dy << ", " 
		   << dz << ")" << endl;
	      cout << "......... part[" << j << "].tot_charge = " 
		   << part[j].tot_charge << endl;
	      cout << "......... rs = " << rs << endl;
	      cout << "......... r = " << r << endl;
	      cout << "......... alpha = " << alpha << endl;
	    }
	    *pot -= part[j].tot_charge*exp(-alpha*r)/r;
	  } 
	}// end add contribution from particle tree[ind].member[j1]
	//} // end parallel for
      }
    }
  }
  return(0);
} // end shielded_treepot




///////////////////////////////////////////////////////
// Function Name: plot_efield_tree
// Usage: outputs the electric field on a mesh using
//        boundary integral corrected direct sum 
//
///////////////////////////////////////////////////////
// Assumptions:
//         none
//
///////////////////////////////////////////////////////
// Inputs:
//         nxnodes  - number of x mesh points
//         nynodes  - number of y mesh points
//         nz nodes - number of z mesh points
//         npart    - number of particles
//         nterms   - number of terms in taylor exp
//         maxmemb  - 
//         acc      - acceptance critron for taylor exp
//         dlength  - lenght of domain in x direction
//         dheight  - lenght of domain in y direction
//         ddepth   - length of domain in z direction
//         part     - list of particles
//         tree     - full oct tree of particles
//
///////////////////////////////////////////////////////
// Outputs:
//        Wries out a file in ascii Ex,Ey,Ez.
//
///////////////////////////////////////////////////////
// Functions Called:
//        treeforce
//
///////////////////////////////////////////////////////
//Note Updared for free space and pereodic
void plot_efield_tree(int nxnodes, int nynodes, int nznodes, int npart, 
		      int nterms, int maxmemb, double acc, 
		      double dlength, double dheight, double ddepth,
		      PARTICLE* part, TREE* tree)  {
  int i;
  int j;
  int k;
  int index=1;

  double dx=dlength/(nxnodes-1.0);
  double dy=dheight/(nynodes-1.0);
  double dz=ddepth/(nznodes-1.0);
  double x=0.0;
  double y=0.0;
  double z=0.0;
  double x_temp;
  double y_temp;
  double z_temp;

  double x_force=0.0;
  double y_force=0.0;
  double z_force=0.0;

  //double half=dlength*0.5;
  char output[128];
    
  sprintf(output,"tree-efield_%i.plt",index);
  ofstream file(output);

  //  x = 0.0;
  //  for (i=0; i<nxnodes; i++)  {
  //    y=0.0;
  //    for (j=0; j<nynodes; j++)  {
  //      z=0.0;
  x = 0.6;
  y = 0.2;
  z = 0.3;
  //    for (k=0; k<nznodes; k++)  {
  x_temp=0.0;
  y_temp=0.0;
  z_temp=0.0;

  // Use the treecode to compute the force
  treeforce(0, -1, nterms, x, y, z, acc, &x_temp, &y_temp, &z_temp, tree, part);
	
	x_force+=x_temp*charge_constant;
	y_force+=y_temp*charge_constant;
	z_force+=z_temp*charge_constant;
	file << setw(20) << setprecision(16) << std::scientific << endl;
	file << x << "\t" << y << "\t" << z << "\t" 
	     << x_force << "\t" << y_force << "\t" << z_force << endl;
	//	z+=dz;
	// }
	//      y+=dy;
	//    }
	//    x+=dx;
	//  }
  
  file.close();

  return;
}





///////////////////////////////////////////////////////
// Function Name: setboundaries
// Usage: initialize boundary of simulation
//
///////////////////////////////////////////////////////
// Assumptions:
//         domain is a rectangular box
//
///////////////////////////////////////////////////////
// Inputs:
//         xpanels  - number of yz panels
//         ypanels  - number of xz panels
//         zpanels  - number of xy panels
//         dlength  - length of x side
//         dheight  - length of y side
//         ddepth   - length of z size
//         pan      - empty panel array
//
///////////////////////////////////////////////////////
// Outputs: (by reference)
//         pan       - all data
//
///////////////////////////////////////////////////////
// Functions Called:
//         none
//
///////////////////////////////////////////////////////
// Global variables used
//
//         npoints - number of quadrature nodes to use
//
///////////////////////////////////////////////////////
// Set up boundaries in the box domain
// 
int setboundaries(int xpanels, int ypanels, int zpanels,
		  double dlength, double dheight, double ddepth,
                  PANEL* pan)  {

  int i=0;
  int j=0;
  int k=0;
  int dirichlet=0;
  int neumann=1;
  int totpan=0;
  int ind=0;
  
  int xytop=dirichlet;
  int xybot=dirichlet;
  int yztop=dirichlet;
  int yzbot=dirichlet;
  int xztop=dirichlet;
  int xzbot=dirichlet;
  
  
  // boundary conditions hard-coded here!!!
  double phixytopref=-50.0;
  double phixybotref=-50.0;
  double phixztopref=0.0;
  double phixzbotref=0.0;
  double phiyztopref=0.0;
  double phiyzbotref=0.0;
  double phixztopref2=-50.0;
  double phixzbotref2=-50.0;
  double phiyztopref2=-50.0;
  double phiyzbotref2=-50.0;
  
  
  double temp1=0.0;
  double temp2=0.0;
  double temp3=0.0;
  double dx;
  double dy;
  double dz;

  double* gpts = NULL;
  
  double* x;
  double halflength = dlength/2.0;
  x=new double[xpanels+1];
  temp1=dlength/(xpanels);
  x[0]=0.0-halflength;
  for (i=1; i<xpanels; i++) {
    x[i]=i*temp1-halflength;
  }
  x[xpanels]=halflength;  

  double* y;
  double halfheight = dheight/2.0;
  y=new double[ypanels+1];
  temp2=dheight/(ypanels);
  y[0]=0.0-halfheight;
  for (i=1; i<ypanels; i++) {
    y[i]=i*temp2-halfheight;
  }
  y[ypanels]=halfheight;

  double* z;
  double halfdepth = ddepth/2.0;
  z=new double[zpanels+1];
  temp3=ddepth/(zpanels);
  z[0]=0.0 -halfdepth;
  for (i=1; i<zpanels; i++) {
    z[i]=i*temp3 - halfdepth;
  }
  z[zpanels]=halfdepth;
  
  // Set up the bottom panels
  for (i=0; i<xpanels; i++)  {
    for (j=0; j<ypanels; j++) {
      ind=j+totpan;  // Increment the panel index
      // Initialize the panel properties
      panelinit(&pan[ind],x[i], x[i+1], y[j], y[j+1], 
		-halfdepth, -halfdepth, phixybotref, 0.0, 0.0, 1.0,
		xybot,1);
    }
    totpan+=ypanels;
  }
  
  // Set up the top panels
  for (i=0; i<xpanels; i++)  {
    for (j=0; j<ypanels; j++) {
      ind=j+totpan;  // Increment the panel index
      // Initialize the panel properties
      panelinit(&pan[ind],x[i], x[i+1], y[j], y[j+1], 
		halfdepth, halfdepth, phixytopref, 0.0, 0.0, -1.0,
		xytop,1);
    }
    totpan+=ypanels;
  }
  
  // Set up the left panels
  for (i=0; i<zpanels; i++)  {
    for (j=0; j<ypanels; j++) {
      ind=j+totpan;  // Increment the panel index
      // Initialize the panel properties
      if ((z[i]>4) && (z[i+1]<5)) { 
	panelinit(&pan[ind],-halflength , -halflength, y[j], y[j+1], 
		  z[i], z[i+1], phiyzbotref, 1.0, 0.0, 0.0,
		  yzbot,3);
      } else {
	panelinit(&pan[ind],-halflength , -halflength, y[j], y[j+1], 
		  z[i], z[i+1], phiyzbotref2, 1.0, 0.0, 0.0,
		  yzbot,3);
      }
    }
    totpan+=ypanels;
  }

  // Set up the right panels
  for (i=0; i<zpanels; i++)  {
    for (j=0; j<ypanels; j++) {
      ind=j+totpan;  // Increment the panel index
      // Initialize the panel properties
      if ((z[i]>4) && (z[i+1]<5)) { 
	panelinit(&pan[ind],halflength, halflength, y[j], y[j+1], 
		  z[i], z[i+1], phiyztopref, -1.0, 0.0, 0.0,
		  yztop,3);
      } else {
	panelinit(&pan[ind],halflength , halflength, y[j], y[j+1], 
		  z[i], z[i+1], phiyztopref2, -1.0, 0.0, 0.0,
		  yztop,3);
      }
    }
    totpan+=ypanels;
  }


  // Set up the front panels
  for (i=0; i<xpanels; i++)  {
    for (j=0; j<zpanels; j++) {
      ind=j+totpan;  // Increment the panel index
      // Initialize the panel properties
      if ((z[j]>4) && (z[j+1]<5)) { 
	panelinit(&pan[ind],x[i], x[i+1], -halfheight, -halfheight, 
		  z[j], z[j+1], phixzbotref, 0.0, 1.0, 0.0,
		  xzbot,2);
      } else {
	panelinit(&pan[ind],x[i], x[i+1], -halfheight, -halfheight, 
		  z[j], z[j+1], phixzbotref2, 0.0, 1.0, 0.0,
		  xzbot,2);
      }
    }
    totpan+=zpanels;
  }

  // Set up the back panels
  for (i=0; i<xpanels; i++)  {
    for (j=0; j<zpanels; j++) {
      ind=j+totpan;  // Increment the panel index
      // Initialize the panel properties
      if ((z[j]>4) && (z[j+1]<5)) { 
	panelinit(&pan[ind],x[i], x[i+1], halfheight, halfheight, 
		  z[j], z[j+1], phixztopref, 0.0, -1.0, 0.0,
		  xztop,2);
      } else {
	panelinit(&pan[ind],x[i], x[i+1], halfheight, halfheight, 
		  z[j], z[j+1], phixztopref2, 0.0, -1.0, 0.0,
		  xztop,2);
      }
    }
    totpan+=zpanels;
  }
  
  delete [] x;
  delete [] y;
  delete [] z;
  
  return(0);
}


///////////////////////////////////////////////////////
// Function Name: panelinit
// Usage: initalize a boundary panel by reference
//
///////////////////////////////////////////////////////
// Assumptions:
//         none
//
///////////////////////////////////////////////////////
// Inputs:
//         panl      - panel by ref
//         x0        - starting x-point
//         x1        - ending x-point
//         y0        - starting y-point
//         y1        - ending y-point
//         z0        - starting z-point
//         z1        - ending z-point
//         refpot    - potential on panel if dirichlet
//         xnrm      - unit norm of panel n = (xnrm,ynrm,znrm)
//         ynrm      - unit norm of panel n = (xnrm,ynrm,znrm)
//         znrm      - unit norm of panel n = (xnrm,ynrm,znrm)
//         type      - panel type: dirichlet or neumann
//
///////////////////////////////////////////////////////
// Outputs: (by reference)
//         panl->x0
//         panl->x1
//         panl->y0
//         panl->y1
//         panl->z0
//         panl->z1
//         panl->midx
//         panl->midy
//         panl->midz
//         panl->refpot
//         panl->xnrm
//         panl->ynrm
//         panl->znrm
//         panl->type
//         pan1->orientation
//
///////////////////////////////////////////////////////
// Functions Called:
//         none
//
///////////////////////////////////////////////////////
// Initialize a panel using the inputs
void panelinit(PANEL* panl, double x0, double x1, double y0, 
	       double y1, double z0, double z1, double refpot, 
	       double xnrm, double ynrm, double znrm, int type,
	       int orientation)  {

  panl->x0=x0;
  panl->x1=x1;
  panl->y0=y0;
  panl->y1=y1;
  panl->z0=z0;
  panl->z1=z1;
  panl->midx=0.5*(x1+x0);
  panl->midy=0.5*(y1+y0);
  panl->midz=0.5*(z1+z0);
  // assume it is dirichlet for now!
  panl->refpot=refpot;
  panl->xnrm=xnrm;
  panl->ynrm=ynrm;
  panl->znrm=znrm;
  panl->type=type;
  panl->orientation=orientation;
  switch (orientation){
  case 1: // xy panel
    panl->area = (x1-x0)*(y1-y0);
    break;
  case 2: // xz panel
    panl->area = (x1-x0)*(z1-z0);
    break;
  case 3: // yz panel
    panl->area = (z1-z0)*(y1-y0);
    break;
  } 
  // we will use the mid point of the panel for quadrature, so no need
  // for gauss points
}

///////////////////////////////////////////////////////
// Function Name: sq
// Usage: computes square of a real number
//
///////////////////////////////////////////////////////
// Assumptions:
double sq(double x){
  return x*x;
}


///////////////////////////////////////////////////////
// Function Name: panelconstantsreal
// Usage: Compute the panel effects array using Riemann
//        integration
//
///////////////////////////////////////////////////////
// Assumptions:
//          ??eps??
//
///////////////////////////////////////////////////////
// Inputs:
//          numpan     - number of panels
//          pan        - list of all panels
//
///////////////////////////////////////////////////////
// Outputs: (by refrence)
//          panelarray - quadrature on panels
//
///////////////////////////////////////////////////////
// Functions Called: 
//
///////////////////////////////////////////////////////
// 
int panelconstantsreal(int numpan, PANEL* pan, double** panelarray){

  int i;
  int j;
  int ind=0;
  double t;
  double r;

  //initialize panelarray
  for (i=0; i<numpan; i++) {
    for (j=0; j<numpan; j++) {
      panelarray[i][j] = 0.0;
    }
  }

  //Iterate through the panels to fill the array
  for (i=0; i<numpan; i++)  {
    for (j=0; j<numpan; j++)  {
      r = sqrt(sq(pan[i].midx-pan[j].midx) +
	       sq(pan[i].midy-pan[j].midy) + 
	       sq(pan[i].midz-pan[j].midz) +
	       DEL);
      panelarray[i][j] = -1/r/2*inv2pi*pan[j].area;
    }
  }
  return(0);
}

///////////////////////////////////////////////////////
// Function Name: badpanelconstantsreal
// Usage: Compute the panel effects array using Riemann
//        integration
//
///////////////////////////////////////////////////////
// Assumptions:
//          ??eps??
//
///////////////////////////////////////////////////////
// Inputs:
//          numpan     - number of panels
//          pan        - list of all panels
//
///////////////////////////////////////////////////////
// Outputs: (by refrence)
//          panelarray - quadrature on panels
//
///////////////////////////////////////////////////////
// Functions Called: 
//
///////////////////////////////////////////////////////
// 
int badpanelconstantsreal(int numpan, PANEL* pan, double** panelarray){

  int i;
  int j;
  int ind=0;
  double temp;

  //initialize panelarray
  for (i=0; i<numpan; i++) {
    for (j=0; j<numpan; j++) {
      panelarray[i][j] = 0.0;
    }
  }

  if (verbosity) {
    cout << "debugging panel array" << endl;
  }

  //Iterate through the panels to fill the array
  for (i=0; i<numpan; i++)  {
    if (verbosity) {
      cout << "i = " << i << endl ;
    }
    for (j=0; j<numpan; j++)  {
      if (verbosity) {
	cout << "j = " << j << ", orientation = "<< pan[j].orientation << endl;
      }
      temp = 0.0;
      switch (pan[j].orientation) {
      case 1:	// xy panels 
	if (verbosity) {
	  cout << "xy panels " << endl;
	  cout << pan[i].midx-pan[j].x0 << endl;
	  cout << pan[i].midy-pan[j].y0 << endl;
	  cout << pan[i].midz-pan[j].z0 << endl;
	  cout << " ..." << endl;
	  cout << pan[i].midx-pan[j].x1 << endl;
	  cout << pan[i].midy-pan[j].y0 << endl;
	  cout << pan[i].midz-pan[j].z0 << endl;
	  cout << " ..." << endl;
	  cout << pan[i].midx-pan[j].x1 << endl;
	  cout << pan[i].midy-pan[j].y1 << endl;
	  cout << pan[i].midz-pan[j].z0 << endl;
	  cout << " ..." << endl;
	  cout << pan[i].midx-pan[j].x0 << endl;
	  cout << pan[i].midy-pan[j].y1 << endl;
	  cout << pan[i].midz-pan[j].z0 << endl;
	}
	temp = (1.0/sqrt( sq(pan[i].midx-pan[j].x0) +
			  sq(pan[i].midy-pan[j].y0) +
			  sq(pan[i].midz-pan[j].z0)) +
		1.0/sqrt( sq(pan[i].midx-pan[j].x1) +
			  sq(pan[i].midy-pan[j].y0) +
			  sq(pan[i].midz-pan[j].z0)) +
		1.0/sqrt( sq(pan[i].midx-pan[j].x1) +
			  sq(pan[i].midy-pan[j].y1) +
			  sq(pan[i].midz-pan[j].z0)) +
		1.0/sqrt( sq(pan[i].midx-pan[j].x0) +
			  sq(pan[i].midy-pan[j].y1) +
			  sq(pan[i].midz-pan[j].z0)));
	break;
      case 2: // xz panels
	if (verbosity) {
	  cout << "xz panels " << endl;
	  cout << pan[i].midx-pan[j].x0 << endl;
	  cout << pan[i].midy-pan[j].y0 << endl;
	  cout << pan[i].midz-pan[j].z0 << endl;
	  cout << " ..." << endl;
	  cout << pan[i].midx-pan[j].x1 << endl;
	  cout << pan[i].midy-pan[j].y0 << endl;
	  cout << pan[i].midz-pan[j].z0 << endl;
	  cout << " ..." << endl;
	  cout << pan[i].midx-pan[j].x1 << endl;
	  cout << pan[i].midy-pan[j].y0 << endl;
	  cout << pan[i].midz-pan[j].z1 << endl;
	  cout << " ..." << endl;
	  cout << pan[i].midx-pan[j].x0 << endl;
	  cout << pan[i].midy-pan[j].y0 << endl;
	  cout << pan[i].midz-pan[j].z1 << endl;
	}
	temp = (1.0/sqrt( sq(pan[i].midx-pan[j].x0) +
			  sq(pan[i].midy-pan[j].y0) +
			  sq(pan[i].midz-pan[j].z0)) +
		1.0/sqrt( sq(pan[i].midx-pan[j].x1) +
			  sq(pan[i].midy-pan[j].y0) +
			  sq(pan[i].midz-pan[j].z0)) +
		1.0/sqrt( sq(pan[i].midx-pan[j].x1) +
			  sq(pan[i].midy-pan[j].y0) +
			  sq(pan[i].midz-pan[j].z1)) +
		1.0/sqrt( sq(pan[i].midx-pan[j].x0) +
			  sq(pan[i].midy-pan[j].y0) +
			  sq(pan[i].midz-pan[j].z1)));
	break;
      case 3: // yz panels
	if (verbosity) {
	  cout << "yz panels " << endl;
	  cout << pan[i].midx-pan[j].x0 << endl;
	  cout << pan[i].midy-pan[j].y0 << endl;
	  cout << pan[i].midz-pan[j].z0 << endl;
	  cout << " ..." << endl;
	  cout << pan[i].midx-pan[j].x0 << endl;
	  cout << pan[i].midy-pan[j].y0 << endl;
	  cout << pan[i].midz-pan[j].z1 << endl;
	  cout << " ..." << endl;
	  cout << pan[i].midx-pan[j].x0 << endl;
	  cout << pan[i].midy-pan[j].y1 << endl;
	  cout << pan[i].midz-pan[j].z1 << endl;
	  cout << " ..." << endl;
	  cout << pan[i].midx-pan[j].x0 << endl;
	  cout << pan[i].midy-pan[j].y1 << endl;
	  cout << pan[i].midz-pan[j].z0 << endl;
	  cout << " ..." << endl;
	}
	temp = (1.0/sqrt( sq(pan[i].midx-pan[j].x0) +
			  sq(pan[i].midy-pan[j].y0) +
			  sq(pan[i].midz-pan[j].z0)) +
		1.0/sqrt( sq(pan[i].midx-pan[j].x0) +
			  sq(pan[i].midy-pan[j].y0) +
			  sq(pan[i].midz-pan[j].z1)) +
		1.0/sqrt( sq(pan[i].midx-pan[j].x0) +
			  sq(pan[i].midy-pan[j].y1) +
			  sq(pan[i].midz-pan[j].z1)) +
		1.0/sqrt( sq(pan[i].midx-pan[j].x0) +
			  sq(pan[i].midy-pan[j].y1) +
			  sq(pan[i].midz-pan[j].z0)));
	break;
      }
      panelarray[i][j] = inv2pi*temp/4.0*pan[j].area;
    }
  }
  if (verbosity){
    // output matrix and vectors 
    cout << "panelarray = [" <<endl;
    for (i=0; i<numpan; i++) {
      for (j=0; j<numpan; j++) {
	cout << panelarray[i][j] << ", ";
      }
      cout << endl;
    }
    cout << "]";
  }

  return(0);
}




///////////////////////////////////////////////////////
// Function Name: panel_solve_ds
// Usage: Solve for the panel strengths using direct sum
//
///////////////////////////////////////////////////////
// Assumptions:
//          none
///////////////////////////////////////////////////////
// Inputs:  
//          numpan     - number of panels
//          npart      - number of particles
//          part       - list of all particles
//          pan        - list of all panels
//          panelarray - quadrature on panels
//
///////////////////////////////////////////////////////
// Outputs: (by reference)
//          pan[i].str -- where sigma gets stored 
//
///////////////////////////////////////////////////////
// Functions Called:
//          gmres() - using solver.h
///////////////////////////////////////////////////////
// 
/*
int panel_solve_ds(int numpan, int npart, PARTICLE* part, 
		   PANEL* pan, double** panelarray)  {

  typedef double *pdoub;

  pdoub b;
  pdoub sigma;
  
  b = new double[numpan];
  sigma = new double[numpan];

  double temp;

  for (int i=0; i<numpan; i++)  {
    b[i]=pan[i].refpot/2.0;
    sigma[i]=0.0;
  }


 // Get the particle contributions at each panel center
  for (int i=0; i<numpan; i++)  {
    pot_point_ds(pan[i].midx,pan[i].midy,pan[i].midz,
		 npart,part,temp);
    b[i] += temp;
  }
  

  // now do matrix solve for constants sigma, and store to pan.str
  solver A_M;
  A_M.makeFram(numpan);
  A_M.gmres(panelarray,b,sigma,numpan);

  for (init i=0; i<numpan; i++){
    pan[i].str = sigma[i];
  }

  // Clean up arrays
  delete [] b;
  delete [] sigma;
}
*/

/*
///////////////////////////////////////////////////////
// Function Name: panelsolvetree  
// Usage: Solve for the panel strengths using the 
//        treecode and matrix inversion
//
///////////////////////////////////////////////////////
// Assumptions:
//          none
//         
///////////////////////////////////////////////////////
// Inputs:  
//          numpan     - number of panles
//          npart      - number of paricles
//          nterms     - number of terms in taylor exp
//          acc        - acceptance crterion for 
//                       aylor approximation
//          panelarray - 
//          part       - list of all particles
//          pan        - list of all panles
//          tree       - enter quad tree
//          invert     - flg for type of matrix invertion
//
///////////////////////////////////////////////////////
// Outputs: (by refrence)
//          pan[i].str 
//
///////////////////////////////////////////////////////
// Functions Called:
//          treeforce()
//          gmresclean() - feval.cpp
//
///////////////////////////////////////////////////////
// 
int panel_solve_tree(int numpan, int npart, int nterms, double acc,
                   double** panelarray, PARTICLE* part, PANEL* pan,
                   TREE* tree)  {

  int i;
  int j;
  int ind;
  int restult;

  typedef double *pdoub;

  pdoub b;
  pdoub sigma;
  
  b = new double[numpan];
  sigma = new double[numpan];

  double temp;

  for (i=0; i<numpan; i++)  {
    b[i]=pan[i].refpot;
    sigma[i]=0.0;
  }

  // Get the particle contributions at each panel center
  //for (i=0; i<numpan; i++)  {
  //  potreal_point(pan[i].midx,pan[i].midy,pan[i].midz,
  //		  npart,temp,part);
  //b[i] += temp;
  //}


  double x;
  double y;
  double z;
  double rs;
  double r;

  //   // Get the particle contributions at each panel center
     for (i=0; i<numpan; i++)  {
       for (j=0; j<npart; j++)  {
         x=pan[i].midx-part[j].x;
         y=pan[i].midy-part[j].y;
         z=pan[i].midz-part[j].z;
         rs=x*x+y*y+z*z;
         r = sqrt(rs);
         if (rs==0.0) {
     	     cout << "Error: Particle has zero distance from panel center : "
   	     << j << endl;
         }
         // what sign is this supposed to be?
         //b[i]+=charge_constant*0.5/r*part[j].tot_charge;
         b[i]+=0.5/r*part[j].tot_charge;
       }
     }

  if (verbosity){
    // output matrix and vectors 
    cout << "A = [" <<endl;
    for (i=0; i<numpan; i++) {
      for (j=0; j<numpan; j++) {
	cout << panelarray[i][j] << ", ";
      }
      cout << endl;
    }
    cout << "]" << endl;
    cout << "b = [";
    for (i=0; i<numpan; i++) {
      cout << b[i] << endl;
    }
    cout << "]";
  }

  // now do matrix solve for constants sigma, and store to pan.str(?)
  solver A_M;
  A_M.makeFram(numpan);
  A_M.gmres(panelarray,b,sigma,numpan);

  for (i=0; i<numpan; i++){
    pan[i].str = sigma[i];
  }

  //   // Clean up arrays
  //   delete [] phis;
  //   delete [] gammas;

  return(0);
}
*/

///////////////////////////////////////////////////////
// Function Name: treetaylorpotential
// Usage: Compute force for a cluster using far 
//        field taylor expansion
//
///////////////////////////////////////////////////////
// Assumptions:
//        none
//
///////////////////////////////////////////////////////
// Inputs:
//        nterms - number of terms in taylor expansion 
//        dx     - x_part - x_cluster_center
//        dy     - y_part - y_cluster_center
//        dz     - z_part - z_cluster_center
//        rs     - sqrt(dx^2+dy^2+dz^2+del)
//        moment - list of moments for given cluster
//
///////////////////////////////////////////////////////
// Outputs: (By Reference)
//         pot - potential at (x,y,z) due to cluster
//
///////////////////////////////////////////////////////
// Functions Called:
//         none
//
///////////////////////////////////////////////////////
void treetaylorpotential(int nterms, double dx, double dy, double dz,
			 double rs, double* pot, double*** moment)  {

  int i;
  int j;
  int k;
  int k1;
  int k2;
  int m;

  double fac=1.0/rs;
  double sqfac=sqrt(fac);
  double ik;  // we will set this to 1/k
  double ii;  // we will set this to 1/i

  double derivs[MAXTERM][MAXTERM][MAXTERM];
  
  // zeroth order derivative
  derivs[0][0][0]=sqfac; 

  // first order derivatives
  derivs[1][0][0]=fac*dx*sqfac;       
  derivs[0][1][0]=fac*dy*sqfac;
  derivs[0][0][1]=fac*dz*sqfac;

  // recursion relation for derivatives in one direction
  for (k=2; k<nterms+1;k++) {
    ik = 1.0/k;
    derivs[k][0][0] = fac*(2.0-ik)*derivs[k-1][0][0]*dx -
      fac*(1-ik)*derivs[k-2][0][0];

    derivs[0][k][0] = fac*(2.0-ik)*derivs[0][k-1][0]*dy -
      fac*(1-ik)*derivs[0][k-2][0];

    derivs[0][0][k] = fac*(2.0-ik)*derivs[0][0][k-1]*dz -
      fac*(1-ik)*derivs[0][0][k-2];
  }

  // derivatives for 
  // G(i,1,0), G(i,0,1), 
  // G(1,i,0), G(0,i,1),
  // G(1,0,i), G(0,1,i)


  derivs[1][1][0] = fac*dx*derivs[0][1][0] + 
    fac*2.0*dy*derivs[1][0][0];
 
  derivs[1][0][1] = fac*dx*derivs[0][0][1] + 
    fac*2.0*dz*derivs[1][0][0];

  derivs[0][1][1] = fac*dy*derivs[0][0][1] + 
    fac*2.0*dz*derivs[0][1][0];

  for (k=2; k<nterms; k++) {
    derivs[1][0][k]=fac*(dx*derivs[0][0][k]+
			 2.0*dz*derivs[1][0][k-1]-
			 derivs[1][0][k-2]);
    
    derivs[0][1][k]=fac*(dy*derivs[0][0][k]+
			 2.0*dz*derivs[0][1][k-1]-
			 derivs[0][1][k-2]);
    
    derivs[0][k][1]=fac*(dz*derivs[0][k][0]+
			 2.0*dy*derivs[0][k-1][1]-
			 derivs[0][k-2][1]);
    
    derivs[1][k][0]=fac*(dx*derivs[0][k][0]+
			 2.0*dy*derivs[1][k-1][0]-
			 derivs[1][k-2][0]);
    
    derivs[k][1][0]=fac*(dy*derivs[k][0][0]+
			 2.0*dx*derivs[k-1][1][0]-
			 derivs[k-2][1][0]);
    
    derivs[k][0][1]=fac*(dz*derivs[k][0][0]+
			 2.0*dx*derivs[k-1][0][1]-
			 derivs[k-2][0][1]);
  }

  // derivatives for G(0,i,j), G(i,0,j), G(i,j,0), i,j >= 2

  for (i=2; i<nterms-1; i++) {
    for (j=2; j<nterms-i+1; j++) {
      ii = 1.0/i;
      derivs[i][j][0]=fac*(dx*(2.0-ii)*derivs[i-1][j][0]+
			   2.0*dy*derivs[i][j-1][0]-
			   (1-ii)*derivs[i-2][j][0]-
			   derivs[i][j-2][0]);
      
      derivs[i][0][j]=fac*(dx*(2.0-ii)*derivs[i-1][0][j]+
			   2.0*dz*derivs[i][0][j-1]-
			   (1-ii)*derivs[i-2][0][j]-
			   derivs[i][0][j-2]);
      
      derivs[0][i][j]=fac*(dy*(2.0-ii)*derivs[0][i-1][j]+
			   2.0*dz*derivs[0][i][j-1]-
			   (1-ii)*derivs[0][i-2][j]-
			   derivs[0][i][j-2]);   
    }
  }
  
  
  // 2 indices 1, other >= 1
  // deriv(1,1,1) is correct, but a little tricky!  
  //      b(1,1,1)=5.0*dz*fac*b(1,1,0)
  
  derivs[1][1][1]=fac*(dx*derivs[0][1][1] +
		      2.0*dy*derivs[1][0][1]+
		      2.0*dz*derivs[1][1][0]);
  
  for (i=2; i<nterms-1; i++ ){
    derivs[1][1][i]=fac*(dx*derivs[0][1][i] + 
			 2.0*dy*derivs[1][0][i] +
			 2.0*dz*derivs[1][1][i-1] -
			 derivs[1][1][i-2]);

    derivs[1][i][1]=fac*(dx*derivs[0][i][1] +
			 2.0*dy*derivs[1][i-1][1] +
			 2.0*dz*derivs[1][i][0] -
			 derivs[1][i-2][1]);

    derivs[i][1][1]=fac*(dy*derivs[i][0][1] +
			 2.0*dx*derivs[i-1][1][1] +
			 2.0*dz*derivs[i][1][0] -
			 derivs[i-2][1][1]);
  }

  // 1 index 1, others >=2
  for (i=2; i<nterms-2; i++) {
    for (j=2; j<nterms-i+1; j++) {
      derivs[1][i][j]=fac*(dx*derivs[0][i][j] +
			   2.0*dy*derivs[1][i-1][j] +
			   2.0*dz*derivs[1][i][j-1] -
                           derivs[1][i-2][j] - 
			   derivs[1][i][j-2]);

      derivs[i][1][j]=fac*(dy*derivs[i][0][j] +
			   2.0*dx*derivs[i-1][1][j] +
			   2.0*dz*derivs[i][1][j-1] -
			   derivs[i-2][1][j] -
			   derivs[i][1][j-2]);

      derivs[i][j][1]=fac*(dz*derivs[i][j][0] +
			   2.0*dx*derivs[i-1][j][1] +
			   2.0*dy*derivs[i][j-1][1] -
                           derivs[i-2][j][1] -
			   derivs[i][j-2][1]);
    }
  }

  // all indices >=2
  for (k=2; k<nterms-3; k++) {
    for (j=2; j<nterms-1-k; j++) {
      for (i=2; i<nterms-k-j+1; i++) {
	ii = 1.0/i;
	derivs[i][j][k]=fac*(2.0*dx*(1-0.5*ii)*derivs[i-1][j][k] +
			     2.0*dy*derivs[i][j-1][k] +
			     2.0*dz*derivs[i][j][k-1] - 
			     (1.0-ii)*derivs[i-2][j][k] - 
			     derivs[i][j-2][k] -
			     derivs[i][j][k-2]); 
      }
    }
  }

  // Add up the Taylor derivatives and moments to get the force
  for (m=0; m<nterms; m++) {
    for (k=0; k<nterms-m+1; k++)  {
      for (i=0; i<nterms-m-k+1; i++)  {
	*pot+=moment[i][k][m]*derivs[i][k][m];
      }
    }
  }
  return;
}


///////////////////////////////////////////////////////
// Function Name: treepotential
// Usage: Compute the potential at location x,y,z due to 
//        all clusters of particles.  Algorithum
//        reverts to direct sum when clusters are 
//        very close to point x,y,z . 
//        (Recursive Function)
//
///////////////////////////////////////////////////////
// Assumptions:
//         eps, set in constant.h
//
///////////////////////////////////////////////////////
// Inputs:
//         ind    - index of tree cluster 
//                  ( when called index = 0 )
//         p      - particle index that we are computing 
//                  force on.  If -1, x,y is not a particle
//         nterms - number of terms in taylor expantion
//         x      - x location
//         y      - y location
//         z      - z location
//         acc    - acceptance criterion for using 
//                  cluster approximation
//         tree   - oct tree sorting particle locations
//         part   - list of all particles
//
///////////////////////////////////////////////////////
// Outputs: (By Reference)
//         pot - potential at (x,y,z)
//
///////////////////////////////////////////////////////
// Functions Called:
//         none
//
///////////////////////////////////////////////////////
double treepotential(int ind, int p, int nterms, double x, double y, double z,
		     double acc, double* pot, TREE* tree, PARTICLE* part)  {
  int go=1;
  int accept;

  double dx=x-tree[ind].xmid;
  double dy=y-tree[ind].ymid;
  double dz=z-tree[ind].zmid;
  double rs=dx*dx+dy*dy+dz*dz + DEL;  // Determine the square of the radius
  double r;

  // Determine whether or not to accept the cluster
  accept=(tree[ind].sqradius < acc*(rs-DEL));

  if (verbosity) {
    cout << "Checking tree " << ind << endl;
    cout << "cluster center = (" << tree[ind].xmid 
	 << ", " << tree[ind].ymid
	 << ", " << tree[ind].zmid << "]" << endl;
    cout << "cluster radius squared = " << tree[ind].sqradius << endl;
    cout << "distance squared from cluster center to r = " << rs <<endl;
    if (accept) {
      cout << "... tree accepted" << endl; 
    } else {
      cout << "... tree NOT accepted" << endl; 
    }
  }


  // On acceptance, use the Taylor expansion
  if (accept)  {
    if (verbosity) {
      cout << "... using tree approximation" << endl;
    }
    treetaylorpotential(nterms, dx, dy, dz, rs, pot,
			tree[ind].moments);
    return(0);
  }

  // Cluster was not accepted
  else  {
    if (verbosity){
      cout << "... children indices: " << endl;
      for (int i = 0; i<NUMCHILD; i++) {
	cout << "...... " << tree[ind].children[i] << endl;
      }
    }
    // Check for child clusters
    for (int i=0; i<NUMCHILD; i++)  {
      if (tree[ind].children[i] != -1)  {
        go=0;
        // Call the recursion
        treepotential(tree[ind].children[i], p, nterms, x, y, z, acc, 
		      pot, tree, part);
      }
    }

    // There are no child clusters, so do direct sum
    if (go)  {
      if (verbosity) {
	cout << "using direct sum ... " << endl;
      }
      int j;
      double temp=0.0;
      // Do direct summation on all particles in the cell
      if (p == -1)  {
        for (int j1=0; j1<tree[ind].nummemb; j1++)  {
          j=tree[ind].members[j1];
          dx=x-part[j].x;
          dy=y-part[j].y;
	  dz=z-part[j].z;
          temp=dx*dx+dy*dy+dz*dz+DEL;
	  r = sqrt(temp);
          if (temp > eps)  {
            temp=part[j].tot_charge/(r);
            *pot+=temp;
          }
        }
      }
      else  {
        for (int j1=0; j1<tree[ind].nummemb; j1++)  {
          j=tree[ind].members[j1];
          if (j != p)  {
            dx=x-part[j].x;
            dy=y-part[j].y;
	    dz=z-part[j].z;
            temp=dx*dx+dy*dy+dz*dz;
	    r = sqrt(temp);
            if (temp > eps)  {
              temp=part[j].tot_charge/(r);
              *pot+=temp;
            }
          }
        }
      }
    }
  }

  return(0);
}


///////////////////////////////////////////////////////
// Function Name: effective_charge
// Usage: initialize effective charges due to boundary
//        conditions
//
///////////////////////////////////////////////////////
// Assumptions:
//          none
///////////////////////////////////////////////////////
// Inputs:  
//          numpan     - number of panels
//          pan        - list of all panels
//
///////////////////////////////////////////////////////
// Outputs: (by reference)
//          bpart -- list of effective particles 
//
///////////////////////////////////////////////////////
/*
void effective_charge(int numpan, PANEL* pan, PARTICLE* bpart,
		      double dlength, double dheight, double ddepth,
		      double density, int npart) {

  for (int i=0; i<numpan; i++ ){
    bpart[i].x = pan[i].midx;
    bpart[i].y = pan[i].midy;
    bpart[i].z = pan[i].midz;
    bpart[i].charge = pan[i].str*pan[i].area;
    bpart[i].tot_charge = bpart[i].charge;
    bpart[i].weight = 1.0;
  }

}
*/

///////////////////////////////////////////////////////
// Function Name: vel_update
// Usage: time integation, part of forward euler
//
///////////////////////////////////////////////////////
// Assumptions:
//          none
//
///////////////////////////////////////////////////////
// Inputs:
//          npart  - number of macro particles
//          part   - particle array
//          dt     - time step
//
///////////////////////////////////////////////////////
// Outputs:(by reference)
//          part.u
//          part.v
//          part.w
//
///////////////////////////////////////////////////////
// Functions Called:
//          none
//
///////////////////////////////////////////////////////
//Update Velocity
int vel_update(int npart, PARTICLE*  part, double dt){
  int i;
  
  for (i=0; i<npart; i++){
    part[i].u+=dt*part[i].xforce* invmass;
    part[i].v+=dt*part[i].yforce* invmass;
    part[i].w+=dt*part[i].zforce* invmass;
  }
  return(0);
}


///////////////////////////////////////////////////////
// Function Name: loc_update
// Usage: time integation
//
///////////////////////////////////////////////////////
// Assumptions:
//          none
//
///////////////////////////////////////////////////////
// Inputs:
//          npart  - number of macro particles
//          part   - particle array
//          dt     - time step
//
///////////////////////////////////////////////////////
// Outputs:(by reference)
//          part.x
//          part.y
//          part.z
//
///////////////////////////////////////////////////////
// Functions Called:
//          Poor_Mans_Perodic()
//
///////////////////////////////////////////////////////
//Update Possition
int loc_update(int npart, PARTICLE*  part, double dt){
  int i;
  for (i=0;i<npart;i++){
    part[i].x+=dt*part[i].u;
    part[i].y+=dt*part[i].v;
    part[i].z+=dt*part[i].w;
  }
  return(0);
}




///////////////////////////////////////////////////////
// Function Name: output_potential
// Usage: Write an ascii file with all the 
//        particle locations and potentials
//
///////////////////////////////////////////////////////
// Assumptions:
//         none
//
///////////////////////////////////////////////////////
// Inputs:
//         npart    - Number of Particles
//         part     - list of all particles
//         filename - Name of Output file
//         potential - vector of potentials
//
///////////////////////////////////////////////////////
// Outputs:
//         none
//
///////////////////////////////////////////////////////
// Functions Called:
//         none
//
///////////////////////////////////////////////////////
void output_potential(int npart, PARTICLE* part, double* pot, string filename)  {
  int i;

  //generate three files here -- xy slice, xz slice, yz slice

  // Output the particle properties to file
  ofstream pfor(filename.c_str());


  //pfor << "i\tx\ty\tz\tu\tv\tw" << endl;
  for (i=0; i<npart; i++)  {
    pfor.precision(14);
    pfor << i << "\t" 
	 << part[i].x << "\t" << part[i].y << "\t" << part[i].z 
	 << "\t" << pot[i] << endl;
  }
  pfor.close();

  return;
} // end output potential


///////////////////////////////////////////////////////
// Function Name: output_fields
// Usage: Write an ascii file with all the 
//        particle locations and fields
//
///////////////////////////////////////////////////////
// Assumptions:
//         none
//
///////////////////////////////////////////////////////
// Inputs:
//         npart    - Number of Particles
//         part     - list of all particles
//         filename - Name of Output file
//         [Bx,By,Bz] - vector of magnetic
//
///////////////////////////////////////////////////////
// Outputs:
//         none
//
///////////////////////////////////////////////////////
// Functions Called:
//         none
//
///////////////////////////////////////////////////////
void output_fields(int npart, PARTICLE* part, 
		   double* Bx, double* By, double* Bz, string filename)  {
  int i;

  // Output the particle properties to file
  ofstream pfor(filename.c_str());
  pfor << "i\tx\ty\tz\tBx\tBy\tBz" << endl;
  for (i=0; i<npart; i++)  {
    pfor.precision(14);
    pfor << i << "\t" 
	 << part[i].x << "\t" << part[i].y << "\t" << part[i].z 
	 << "\t" << Bx[i] << "\t" << By[i] << "\t" << Bz[i] << endl;
  }
  pfor.close();

  return;
} // end output_fields

///////////////////////////////////////////////////////
// Function Name: output_slice
// Usage: Write an ascii file with 
//        particle locations and potentials at slices
//
///////////////////////////////////////////////////////
// Assumptions:
//         none
//
///////////////////////////////////////////////////////
// Inputs:
//         npart    - Number of Particles
//         part     - list of all particles
//         headername - header name of Output file
//         potential - vector of potentials
//
///////////////////////////////////////////////////////
// Outputs:
//         none
//
///////////////////////////////////////////////////////
// Functions Called:
//         none
//
///////////////////////////////////////////////////////
void output_slice(int N, PARTICLE* part, double* pot, string headername)  {


  string filename;
  filename = headername + "_xy.dat";
  // Output the particle properties to file
  ofstream pfor(filename.c_str());

  int i = N/2;
  int ind;
  for (int j=0; j<N; j++) {
    for (int k=0; k<N; k++) {
      ind = i*N*N + j*(N) + k;
      pfor.precision(14);
      pfor << pot[ind] << endl;
    }
  }
  pfor.close();

  filename = headername + "_xz.dat";
  // Output the particle properties to file
  ofstream pfor_xz(filename.c_str());

  int j = N/2;
  for (int i=0; i<N; i++) {
    for (int k=0; k<N; k++) {
      ind = i*N*N + j*(N) + k;
      pfor_xz.precision(14);
      pfor_xz << pot[ind] << endl;
    }
  }
  pfor_xz.close();

  filename = headername + "_yz.dat";
  // Output the particle properties to file
  ofstream pfor_yz(filename.c_str());

  int k = N/2;
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
      ind = i*N*N + j*(N) + k;
      pfor_yz.precision(14);
      pfor_yz << pot[ind] << endl;
    }
  }
  pfor_yz.close();

} // end output_slice


///////////////////////////////////////////////////////
// Function Name: output_particles
// Usage: Write an ascii file with all the 
//        particle locations and fields
//
///////////////////////////////////////////////////////
// Assumptions:
//         none
//
///////////////////////////////////////////////////////
// Inputs:
//         npart    - Number of Particles
//         part     - list of all particles
//         filename - Name of Output file
//
///////////////////////////////////////////////////////
// Outputs:
//         none
//
///////////////////////////////////////////////////////
// Functions Called:
//         none
//
///////////////////////////////////////////////////////
void output_particles(int npart, PARTICLE* part, string filename)  {
  int i;

  // Output the particle properties to file
  //ofstream pfor(filename.c_str());

  ofstream pfor;
  pfor.open (filename.c_str(),  ios::app);

  //pfor << "i\tx\ty\tz\tu\tv\tw" << endl;
  for (i=0; i<npart; i++)  {
    pfor.precision(14);
    pfor << part[i].x << "\t" << part[i].y << "\t" << part[i].z 
      //<< endl;
    //pfor << "\t" << part[i].u << "\t" << part[i].v << "\t" 
    	 << "\t" << part[i].w << endl;
  }
  pfor.close();

  return;
}




///////////////////////////////////////////////////////
// Function Name: plot_efield_ds
// Usage: outputs the electic field on a mesh using
//        direct summation
//
///////////////////////////////////////////////////////
// Assumptions:
//         none
//
///////////////////////////////////////////////////////
// Inputs:
//         nxnodes  - number of x mesh points
//         nynodes  - number of y mesh points
//         nznodes  - number of z mesh points
//         dlength  - length of domain in x direction
//         dheight  - length of domain in y direction
//         ddepth   - length of domain in z direction 
//         npart    - number of particles
//         part     - list of particles
//
///////////////////////////////////////////////////////
// Outputs:
//        Writes out a file in ascii Ex, Ey, Ez.
//
///////////////////////////////////////////////////////
// Functions Called:
//        force_point_ds(x, y, z npart, &x_temp, 
//                       &y_temp, &z_temp, part);
///////////////////////////////////////////////////////

void plot_efield_ds(int nxnodes, int nynodes, int nznodes,
		    double dlength, double dheight, double ddepth,
		    int npart, PARTICLE* part)  {
  int i;
  int j;
  int k;
  int index=1;

  double dx=dlength/(nxnodes-1.0);
  double dy=dheight/(nynodes-1.0);
  double dz=ddepth/(nznodes-1.0);
  double x=0.0;
  double y=0.0;
  double z=0.0;
  double x_temp;
  double y_temp;
  double z_temp;
  double x_force=0.0;
  double y_force=0.0;
  double z_force=0.0;

  //  double half=dlength*0.5;
  char output[128];
    
  sprintf(output,"ds-efield_%i.plt",index);
  ofstream file(output);

  //file << "title='Direct Sum Electric Field'\nvariables=x, y, z, E_x, E_y, E_z\n";
  file << "zone t=main, i=" << nxnodes << ", j=" 
       << nynodes << ", k=" << nznodes << ", f=point\n";
  file << "particle locations:\n";
  for (i=0; i<npart; i++) {
    file << "particle[" << i << "] at (" 
	 << part[i].x << "," 
	 << part[i].y << "," 
	 << part[i].z << ")\n";
  }
  
  //x=0.0;
  //  for (i=0; i<nxnodes; i++)  {
  //  y=0.0;
  //  for (j=0; j<nynodes; j++)  {
  //    z=0.0;
  //    for (k=0; k<nznodes; k++)  {
  x = 0.6;
  y = 0.2;
  z = 0.3;

  	x_temp=0.0;
  	y_temp=0.0;
  	z_temp=0.0;
	
	// Use direct sum to compute the force 
	force_point_ds(x,y,z,npart,part,&x_temp,&y_temp,&z_temp);
  
	x_force+=x_temp*charge_constant;
	y_force+=y_temp*charge_constant;
	z_force+=z_temp*charge_constant;
  
	file << setw(20) << setprecision(16) <<std::scientific << endl;
	file << x << "\t" << y << "\t" << z << "\t" 
	     << x_force << "\t" << y_force << "\t" << z_force << endl;
	//	z+=dz;
	//}
	//y+=dy;
	//    }
	//    x+=dx;
	//  }
  
  file.close();

  return;
}
