#ifndef _UTILS_H
#define _UTILS_H


//SYS LIBRARIES
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

//ROOT LIBRARIES
#include <TLorentzVector.h>

#include "constants.h"

using namespace std;



// get random number in the interval [min, max].
double inline randomIntv(double min, double max){
  return min + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(max-min)));
};
// get random number in the interval [min, max].
double inline random(double min, double max){
  return min + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(max-min)));
};


// print TLorentzVector
void prnLV(string s, TLorentzVector v) {
   cout << s << v.E() <<" "<< v.Px() <<" "<< v.Py() <<" "<< v.Pz() << endl;
}


// Legendre polinoms
double PLeg(int l, double ct) {
  if     (l==0) { return 1; }
  else if(l==1) { return ct; }
  else if(l==2) { return 0.5*(3*ct*ct - 1.); }
  else if(l==3) { return 0.5*(5*ct*ct*ct - 3*ct); }
  else if(l==4) { return 0.125*(35*ct*ct*ct*ct - 30*ct*ct + 3.); }
  else if(l==5) { return 0.125*(63*ct*ct*ct*ct*ct - 70*ct*ct*ct + 15*ct); }
  else if(l==6) { return 0.0625*(231*ct*ct*ct*ct*ct*ct - 315*ct*ct*ct*ct + 105*ct*ct - 5.); }  
  else {
    cerr << " PL : Wrong L! Stop." << endl;
    exit(1); 
  }
}






double getMax(int n, double a[]) {
  double amax=a[0];
  for(int i=0; i<n; i++) if(a[i]>amax) amax=a[i]; 
  return amax;
}
double getMax(int n, double a[], double b[]) {
  double amax=a[0];
  for(int i=0; i<n; i++) {
    if(a[i]>amax) amax=a[i]; 
    if(b[i]>amax) amax=b[i];
  } 
  return amax;
}
double getMin(int n, double a[]) {
  double amin=a[0];
  for(int i=0; i<n; i++) if(a[i]<amin) amin=a[i]; 
  return amin;
}



double getLoopVal(int N, double min, double max, int i) {
  double val = min + (max-min) *i /(N-1.);
  return val;
}



#endif

