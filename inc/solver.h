//-------------------------------------------------------
//Title:        Matrix Solver (GMRES) 
//purpos:       Solve Ax=B
//By:           Jerry Emhoff
//Rewriten By:  Andrew J. Christlieb
//Date:         8/26/03-8/26/03
//Contact: Andrew J. Christlieb, Ph. D.
//         Assistant Professor
//         Department of Mathematics
//         University of Michigan
//         2470 East Hall
//         Ann Arbor, MI 48109.
//         
//         Office: 4851 East Hall
//         E-mail: christli@umich.edu
//         Tel:    (734) 763-5725
//------------------------------------------------------
//
// Modified by bwo, 05-19-2009

#ifndef solveh
#define solveh 1

#include <cmath>
#include <iostream>
using namespace std;

typedef double *pdoub;  
typedef pdoub  *ppdoub; 

class solver{
  int MAXPAN;
  ppdoub q;
  pdoub  v;
  ppdoub h;
  pdoub  y;

  ppdoub q2;
  ppdoub v2;
  ppdoub r2;

  double bnorm;
  //
  public:
  solver(){MAXPAN=0;q='\0';v='\0';h='\0';y='\0';
                      q2='\0';v2='\0';r2='\0';
                      bnorm=1;};
 
  void makeFram(int pan){
    int i,j;

    if(pan>MAXPAN){
      clear_all();
      MAXPAN=pan;
    };

    //make main arrays
    j=pan+1;
    q=new pdoub[pan];
    for(i=0;i<pan;i++){
      q[i]=new double[j];
    };
    v=new double[pan];
    h=new pdoub[pan+2];
    for(i=0;i<pan+2;i++){
      h[i]=new double[j];
    };
    y=new double[j];

    //make tmp arrays:
    q2=new pdoub[pan+2];
    for(i=0;i<pan+2;i++){
      q2[i]=new double[j];
    };
    v2=new pdoub[pan+2];
    for(i=0;i<pan+2;i++){
      v2[i]=new double[j];
    };
    r2=new pdoub[pan+2];
    for(i=0;i<pan+2;i++){
      r2[i]=new double[j];
    };
  };

  ~solver(){
    clear_all();
    //cout << "Cleaing up memory used in matrix solver.\n";
  };

  void clear_all(void){
    clear_q();
    clear_v();
    clear_h();
    clear_y();
    clear_q2();
    clear_v2();
    clear_r2();
  };  
  void clear_q(void){
    int i;
    if(q!='\0'){
      for(i=0;i<MAXPAN;i++){
        delete[] q[i];
      };
    };
    delete[] q;
  };
  void clear_v(void){delete[] v;};
  void clear_h(void){
    int i;
    if(h!='\0'){
      for(i=0;i<MAXPAN+2;i++){
        delete[] h[i];
      };
    };
    delete[] h;
  };
  void clear_y(void){delete[] y;};

  void clear_q2(void){
    int i;
    if(q2!='\0'){
      for(i=0;i<MAXPAN+2;i++){
        delete[] q2[i];
      };
    };
    delete[] q2;
  };
  void clear_v2(void){
    int i;
    if(v2!='\0'){
      for(i=0;i<MAXPAN+2;i++){
        delete[] v2[i];
      };
    };
    delete[] v2;
  };
  void clear_r2(void){
    int i;
    if(r2!='\0'){
      for(i=0;i<MAXPAN+2;i++){
        delete[] r2[i];
      };
    };
    delete[] r2;
  };
 
  //-------------------------------
  //MATRIX SOLVER:
  //  GMRES algorithm 
  //  Solve Ax=b
  //-------------------------------
  //A=a[up to MAXPAN][up to MAXPAN]
  //b=b[up to MAXPAN]
  //x=x[up to MAXPAN]
  //num= max counter (N-1).
  int gmres(double **a, double *b, double *x, int num)  { 
     
    int i=0; 
    int j=0; 
    int k; 
    int n=0; 
    int go=1; 
    int n1; 

    if(num>MAXPAN-1){makeFram(num+1);};
 
    double temp=0.0; 
    double error=0.0; 
    double maxerror=0.0; 
    double errlim=1e-5; 
 
    //--------------------------
    //zero h matrix
    //--------------------------
    for (i=0; i<num+2; i++)  { 
      for(j=0; j<num+1;j++) {h[i][j]=0;}; 
    }; 
 
    //--------------------------
    //set bnorm
    //--------------------------
    bnorm=0.0;
    for (i=0; i<num; i++) {bnorm+=b[i]*b[i];}; 
    bnorm=sqrt(bnorm); 

    //--------------------------
    // q(0)=b/||b|| 
    //--------------------------
    for (i=0; i<num; i++) {q[i][0]=b[i]/bnorm;}; 

    while (go) { 
        n1=n+1;
 
        //--------------------------
        //v=Aq(n) 
        //--------------------------
        for (i=0; i<num; i++)  { 
            v[i]=0.0; 
            for (j=0; j<num; j++) {v[i]+=a[i][j]*q[j][n];}; 
        } 
 
        for (j=0; j<n1; j++)  { 

            //h(j,n)=q(j)*v 
            h[j][n]=0.0; 
            for (i=0; i<num; i++) {h[j][n]+=q[i][j]*v[i];}; 

            //v=v-h(j,n)*q(j) 
            for (i=0; i<num; i++) {v[i]-=h[j][n]*q[i][j];}; 
        }; 

        //--------------------------
        // h(n+1,n)=||v|| 
        //--------------------------
        temp=0.0; 
        for (i=0; i<num; i++) {temp+=v[i]*v[i];}; 
        h[n1][n]=sqrt(temp); 

        if ((n1 < num) && (h[n1][n] < 1e-14))  { 
	  cout << h[n1][n] << endl;
	  cout << "num = " << num;
	  cout << "ERROR: Breakdown!  Quitting..." << n << endl; 
	  go=0; 
        }; 
 
        //--------------------------
        // q(n+1)=v/h(n+1,n) 
        //--------------------------
        for (i=0; i<num; i++) { q[i][n1]=v[i]/h[n1][n];}; 
 
        //---------------------------------------------------
        //End of pure Arnoldi iteration 
        //Minimize ||(H*y-||b||e1)|| (Least squares problem) 
        //---------------------------------------------------
        minimize(n1); 
 
        if (n1==num) {go=0;}; 
 
        //--------------------------
        // x=Q*y 
        //--------------------------
        if (n1%2 || !go)  { 
            maxerror=0.0; 
            for (i=0; i<num; i++)  { 
                temp=x[i]; 
                x[i]=0; 
                for (j=0; j<n1; j++)  {x[i]+=q[i][j]*y[j];}; 
 
                error=fabs((temp-x[i])/temp); 
                if (error > maxerror){maxerror=error;}; 
            }; 

            if (maxerror < errlim && n>1) {go=0;}; 
        }; 
        n++; 
    }; 
 
    return(0); 
  }; 
 
  //--------------------------------------------
  //Assume:
  //computes updated y[MAXPAN+1]
  //--------------------------------------------
  int minimize (int n)  { 
 
    int i=0; 
    int j=0; 
    int k=0; 
    int ind=0; 
    int m=n+1; 
 
    //------------------------------------
    //copy needed data from h[i][j] to tmp
    for(j=0; j<n; j++)  { 
        for (i=0; i<j+2; i++) {v2[i][j]=h[i][j];}; 
    }; 
 
    for (i=0; i<n; i++)  { 
        ind=i+2; 
        r2[i][i]=0; 
        for (k=0; k<ind; k++) {r2[i][i]+=v2[k][i]*v2[k][i];}; 

        //Calculate rjj 
        r2[i][i]=sqrt(r2[i][i]); 

        //Calculate q 
        for (k=0; k<ind; k++) {q2[k][i]=v2[k][i]/r2[i][i];}; 
 
        for (j=i+1; j<n; j++)  { 
            r2[i][j]=0.0; 

            //Do the dot product... 
            for(k=0; k<ind; k++) {r2[i][j]+=q2[k][i]*v2[k][j];}; 

            //Calculate v 
            for(k=0; k<ind; k++) {v2[k][j]-=r2[i][j]*q2[k][i];}; 
        }; 
    }; 
 
    //Back substitution 
    for (i=n-1; i>-1; i--)  { 
        y[i]=bnorm*q2[0][i]; 
        for (j=i+1; j<n; j++) { y[i]-=r2[i][j]*y[j];};
        y[i]/=r2[i][i]; 
    }; 
 
    return(0); 
  }; 


};

#endif
