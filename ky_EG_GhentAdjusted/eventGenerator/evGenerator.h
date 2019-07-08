#ifndef _EVGENERATOR_H
#define _EVGENERATOR_H

#include "string.h"
#include <time.h>

#include <TLorentzVector.h>

#include "../generalUtils/constants.h"
#include "../generalUtils/kinematics.h"
#include "sigmaKY.h"




class evGenerator {

  string type;
  double m1;
  double m2;
  double m12;
  double m22;
  double Ebeam;
  double Q2min, Q2max;
  double Wmin, Wmax;
  double d5sigmaMax, d5sigmaMax2;
  double intKLambda, intKSigma;
  sigmaKY  *model;
  sigmaKY  *model2;
  
  int nEvent;

public:

evGenerator(string dataPath, string t, double E, 
            double q2min, double q2max,
	    double wmin,  double wmax)
  {
  type = t;
  Ebeam = E;
  Q2min = q2min;      
  Q2max = q2max;
  Wmin = wmin;
  Wmax = wmax;
  nEvent = 0;
    
  if(type == "KLambda") {
    m1 = massKaon; 
    m2 = massLambda;
  } else if(type == "KSigma") {
    m1 = massKaon; 
    m2 = massSigma0;
  } else if(type == "KLambdaKSigma") {
    m1  = massKaon; 
    m2  = massLambda;
    m12 = massKaon; 
    m22 = massSigma0;
  } else {
    cerr << " Error! evGenerator::evGenerator Wrong reaction type: " <<type<< endl;
    cerr << " It must be KLambda, KSigma or KLambdaKSigma";  
    exit(1);
  }  
  
  // read data file
  if(type == "KLambda" || type == "KSigma") {
    model = new sigmaKY(dataPath,type);
    model->use_costeta();
  } else {
    model  = new sigmaKY(dataPath,"KLambda");
    model-> use_costeta();
    model2 = new sigmaKY(dataPath,"KSigma");
    model2->use_costeta();
  }
  
  
  // Find maximum of the cross section
  if(type == "KLambda" || type == "KSigma") {
    d5sigmaMax = model->d5sigma_max(Ebeam, Q2min, Q2max, Wmin, Wmax);
  } else {
    d5sigmaMax  = model->d5sigma_max(Ebeam, Q2min, Q2max, Wmin, Wmax);
    d5sigmaMax2 = model2->d5sigma_max(Ebeam, Q2min, Q2max, Wmin, Wmax);
    
    int nQ2=50, nW=50;
    double intKLambda=0, intKSigma=0;
    for(int iQ2=0; iQ2<nQ2; iQ2++) {
    for(int iW=0;  iW<nW;  iW++) {
      double Q2 = Q2min + (iQ2)*((Q2max-Q2min)/(nQ2-1));
      double W = Wmin + (iW)*((Wmax-Wmin)/(nW-1));
      intKLambda += model->sigma_int(Ebeam, Q2, W);
      intKSigma += model2->sigma_int(Ebeam, Q2, W);
    }
    }
    cout << " SSS " << intKLambda <<" "<< intKSigma << endl;
  
  }

  // initialize random seed:
  srand (time(NULL));
 
};

~evGenerator(){
    delete model;
    if(type == "KSigmaKLambda") delete model2;
};




		 
		 
// all parameters are output parameters
// 4-momenta of particles in the final state in LAB frame:
//    Pefin, PK, PL -> electron, Kaon, Lambda/Sigma
// 4 momenta of 

void getEvent(double &Q2, double &W, 
              TLorentzVector &Pefin, TLorentzVector &PK, TLorentzVector &PY,
	      TLorentzVector &Ppfin, TLorentzVector &Ppim, TLorentzVector &Pgam) {

  int nTry=0;
  nEvent=0;


  while(true) {
 
    Q2 = randomIntv(Q2min, Q2max);
    W  = randomIntv(Wmin, Wmax);
    double cosThetaK = randomIntv(-0.999999,0.999999);
    double thetaK = acos(cosThetaK);
    double phiK = randomIntv(0.0, 2.*constantPi);

    double d5sigma;
    //if(Q2<3 && W<2) {
    //  d5sigma = 0;
    //} else {
      d5sigma = model->d5sigma(Ebeam, Q2, W, thetaK, phiK);
    //}
    nTry++;
        
    //int 
    randomIntv(0.,intKLambda+intKSigma);
    
    
    
        
    if(randomIntv(0.,1.) < d5sigma/d5sigmaMax) {
       nEvent++;
       
       //cout << "   NNN  " << nTry << " " << nEvent << endl;

            
       double W2 = W*W;
       double omega = getomega(Q2,W);
       if(omega <=0. ) continue;
       double Ee = Ebeam - omega;
       if(Ee<massElectron) continue;
       double arg = 1. - Q2/(2.*Ebeam*Ee);
       if( fabs(arg) > 1. ) continue;
       double theta = acos(arg);
       double phi = randomIntv(0.0, 2.*constantPi);
       double gamma = getGamma(Ebeam,Q2,W);
       
       // 4-momentum of final electron in LAB frame
       Pefin.SetXYZT(1.,1.,1.,Ee);
       Pefin.SetRho(sqrt(Ee*Ee - massElectron2)); 
       Pefin.SetTheta(theta);
       Pefin.SetPhi(phi); 
     
       cms2lab(W, Q2, phi, Ebeam, thetaK, phiK, m1, m2, PK, PY);     


       // Decay of Lambda into proton and pi minus. 
       if(type == "KLambda" ) {
	 getLdecayProd(PY, Ppfin, Ppim);
       } else {
       // Decay of Sigma into proton and pi minus and gamma. 
         TLorentzVector PL;
         getSdecayProd(PY, PL, Ppfin, Ppim, Pgam);
       }
       
       //cout << " PY " << PY.E() <<" "<< PY.Px() <<" "<< PY.Py() <<" "<< PY.Pz() << endl; 
       //cout << " PS " << (Ppfin+Ppim).E() <<" "<< (Ppfin+Ppim).Px() <<" "<< (Ppfin+Ppim).Py() <<" "<< (Ppfin+Ppim).Pz() << endl; 
       //cout << " PK " << (PK).E() <<" "<< (PK).Px() <<" "<< (PK).Py() <<" "<< (PK).Pz() << endl; 
       //cout << " Pp " << (Ppfin).E() <<" "<< (Ppfin).Px() <<" "<< (Ppfin).Py() <<" "<< (Ppfin).Pz() << endl; 
       //cout << " Pi " << (Ppim).E()  <<" "<< (Ppim).Px()  <<" "<< (Ppim).Py()  <<" "<< (Ppim).Pz()  << endl; 
       
       return;
    
    }

  
  }

}; //end getEvent(...)

	
};






#endif 
