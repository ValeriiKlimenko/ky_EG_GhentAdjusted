#ifndef _KINEMATICS_H
#define _KINEMATICS_H


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


// calculate energy of virtual photon in LAB frame (omega)
// from Q2 and W.
// (It's frequently used.)
double inline getomega(double Q2, double W)
{
  double om = (W*W + Q2 - massProton2) / (2.*massProton);
  if (om <= 0.) return 0.;
  return om;
}


// calculate energy of virtual photon in CMS frame 
// from Q2 and W.
double inline getomegastar(double Q2, double W)
{
  double omstar = (W*W - Q2 - massProton2) / (2.*W);
  if (omstar <= 0.) return 0.;
  return omstar;
}



double getKH(double Q2, double W) {
   getomega(Q2, W) - (Q2/2./massProton);
}



// longitudinal polarization parameter
double getEpsilon(double Ebeam, double Q2, double W)
{
  double omega = getomega(Q2, W);  
  return 1./
            (1. + ( 2.*(Q2 + omega*omega) )
	         /( 4.*Ebeam*(Ebeam-omega) - Q2 ) );
}  

// 
double getEpsilonL(double Ebeam, double Q2, double W)
{
  double omega =  getomega(Q2, W);
  double eps = getEpsilon(Ebeam, Q2, W);
  return (sqrt(Q2)/omega)*eps;
}  



// calculate virtual photon flux (Gamma)
double getGamma(double Ebeam, double Q2, double W)
{
  double mp=massProton;
  double mp2 = massProton2;
  double W2 = W*W;
  
  double omega = getomega(Q2,W);;
  if(omega<0. || omega>Ebeam) return 0.;
  double eps2 = Ebeam - omega;
  return    ( constantAlpha/(2.*constantPi2) )
          * (eps2/Ebeam)
          * (W2-mp2)/(2.*mp*Q2)
	  * (1./(1.-getEpsilon(Ebeam, Q2, W)));

}






// calculate 4-momenta of Kaon (PK) and Lambda in LAB frame
// from:  
// W, Q2,  Phi of electron in LAB frame, Ebeam.
// thetaK, phiK are angles of Kaon in CM frame, 
// mK, mL are masses of Kaon and Lambda
// output is PK and PL
void cms2lab(double W, double Q2, double phi, double Ebeam, 
	     double thetaK, double phiK, double mK, double mL,
	     TLorentzVector &PK, TLorentzVector &PL) 
{

  double mp=massProton;
  double pi=constantPi;
  double mp2 = massProton2;
  double W2 = W*W;
  double thetaL = pi-thetaK;
  double phiL = phiK+pi; if(phiL>2.*pi) phiL=phiL-2.*pi;

  //Energy of virtual photon is LAB frame
  double Egamma_lab = (W2 + Q2 - mp2)/(2.*mp);
  //Energy of virtual photon is CMS-frame
  double Egamma_cms = (W2 - Q2 - mp2)/(2.*W);

   //4-momenta of virt. photon and proton in  in CMS
  TLorentzVector Pgamma_cms,Peini_cms,Pp_cms;
  double PgammaZ_cms = sqrt(Egamma_cms*Egamma_cms+Q2);
  Pgamma_cms.SetXYZT(0,0,PgammaZ_cms,Egamma_cms); 
  Pp_cms.SetXYZT(0,0,PgammaZ_cms,sqrt(mp2+PgammaZ_cms*PgammaZ_cms));

  //4-momenta of K and L in CMS. Set energy components 
  PK.SetXYZT(1.,1.,1.,(W*W+mK*mK-mL*mL)/2./W);
  PL.SetXYZT(1.,1.,1.,(W*W+mL*mL-mK*mK)/2./W);
  //Set momentum muduli and angles
  PK.SetRho(sqrt(PK.E()*PK.E()-mK*mK));
  PL.SetRho(sqrt(PL.E()*PL.E()-mL*mL));
  PK.SetTheta(thetaK);
  PK.SetPhi(phiK);
  PL.SetTheta(thetaL);
  PL.SetPhi(phiL);
  PK.Boost(0,0,sqrt(Egamma_lab*Egamma_lab+Q2)/(Egamma_lab+mp));
  PL.Boost(0,0,sqrt(Egamma_lab*Egamma_lab+Q2)/(Egamma_lab+mp));
  // Now we are in the LAB-frame (let's mark it frame 2) with 
  // Z2 - along the virtual photon, X2 - in the ee'-PLane and 
  // Y2 - perpendicular to the e'-PLane.
  // The only thing left is to rotate our coordinate axes into usual 
  // LAB-frame-axes-orientation (lets call it frame 1): 
  // Z1-axis along the beam, 
  // Y1 - up, X1 - to the south.
  // For that purpose we will find coordinates of unit vectors along 
  // the axes of frame 2 in frame 1.
  TVector3 Pefin,Pgamma,PX,PY; 
  // We define 3-vector of scattered electron in LAB-frame 1 
  // theta is from the formula for Q2 in the LAB frame: 
  // Q2 = 2 E_beam E_e' (1-cos(theta_e'))
  Pefin.SetXYZ(0.,0.,1.);
  Pefin.SetTheta(acos( 1. - Q2/(2.*Ebeam*(Ebeam-Egamma_lab))) );
  Pefin.SetPhi(phi); 

  // 3-vector of virtual photon in LAB-frame
  Pgamma.SetXYZ(0.,0.,1.);
  Pgamma.SetTheta(acos( (Q2+2.*Ebeam*Egamma_lab)
                        /(2.*Ebeam*(sqrt(Q2+Egamma_lab*Egamma_lab))) ) );
  if (phi < pi)  Pgamma.SetPhi(phi + pi);
  else           Pgamma.SetPhi(phi - pi);
  //Pgamma is along Z2 and now it is written in the frame 1

  // PY is along Y2 (perpendicular to the scattering PLane)
  PY = (Pgamma.Cross(Pefin)).Unit();
  // PX is along X2 (in scattering PLane)
  PX = (PY.Cross(Pgamma)).Unit();

  // Now unit vectors of frame 2 are written in the frame 1. 
  // We define rotation of frame 2 -> frame 1 via the Euler Angles  	   
  TRotation vrot;
  vrot.SetXEulerAngles(atan2(PX.Z(),PY.Z()), acos(Pgamma.Z()), 
                       atan2(Pgamma.X(),-1.*Pgamma.Y()));
  	   
  // transfom our 4-momenta of final hadrons according to axes transformation
  PK.Transform(vrot);
  PL.Transform(vrot);

  return;
} //end cms2lab(...)




// calculate Theta and Phi of the Kaon in CMS frame
// from:  
// Q2, Ebeam and 4-momenta of electon (Pefin), Kaon (PK), and Lambda (PL)
// in LAB frame
void lab2cms(double Q2, double Ebeam, 
     TLorentzVector Pefin, TLorentzVector PK, TLorentzVector PL, 
     double &theta_hadr, double &phi_hadr)
{

  
  double mp=massProton;

  double a_gamma, b_gamma, a_beta,b_beta, beta;
  TVector3 Vect3_gamma, Vect3_beta, V3_anti_z(0,0,-1);
  TLorentzVector Peini, Pgamma; 
  TRotation rot; 
 
  Peini.SetXYZT(0.,0.,Ebeam, Ebeam);

  Pgamma = Peini - Pefin;

  TVector3 uz = Pgamma.Vect().Unit();
  TVector3 ux = (Peini.Vect().Cross(Pefin.Vect())).Unit();
  ux.Rotate(3.*constantPi/2,uz);
  rot.SetZAxis(uz,ux).Invert();
  PK.Transform(rot);
  PL.Transform(rot);
  Pgamma.Transform(rot);

  beta = sqrt(Pgamma.E()*Pgamma.E()+Q2)/(Pgamma.E()+mp);
 
  PK.Boost(0.,0.,-beta);
  PL.Boost(0.,0.,-beta);
  Pgamma.Boost(0.,0.,-beta);

  theta_hadr = PK.Theta();

  if (PK.Phi()>0) phi_hadr = PK.Phi();
  else            phi_hadr = PK.Phi()+2.*constantPi;

  return;
 
 };




void getLdecayProd(TLorentzVector &PL,
                   TLorentzVector &Ppfin, TLorentzVector &Ppim) {

   //4-momenta Ppfin and Ppim in PL rest frame
   //Set energy components
   double mL=massLambda;
   Ppfin.SetXYZT(1.,1.,1.,(mL*mL+massProton2-massPion2)/2./mL);
   Ppim.SetXYZT (1.,1.,1.,(mL*mL+massPion2-massProton2)/2./mL);
   
   //Set momentum muduli and angles
   Ppfin.SetRho(sqrt(Ppfin.E()*Ppfin.E()-massProton2));
   Ppim.SetRho (sqrt(Ppim.E()*Ppim.E()  -massPion2));
   double cost      = randomIntv(-1., 1.);
   double thetapfin = acos(cost);
   double phipfin   = randomIntv(0.0, 2.*constantPi);
   double thetapim = constantPi-thetapfin;
   double phipim   = phipfin+constantPi; if(phipim>2.*constantPi) phipim=phipim-2.*constantPi;
   Ppfin.SetTheta(thetapfin);
   Ppfin.SetPhi(phipfin);
   Ppim.SetTheta(thetapim);
   Ppim.SetPhi(phipim);
        
   //   cout << " SSSS " <<  Ppfin.E() + Ppim.E() << " " << endl;
   //	cout << " EECMS " <<  Ppfin.E() << " " <<  Ppim.E() << " " << endl;
   //	cout << " PPCMS " <<  Ppfin.Vect().Mag() << " " <<  Ppim.Vect().Mag() << " " << endl;
   //	cout << " MMCMS " <<  Ppfin.M() << " " <<  Ppim.Mag() << " " << endl;
	 
	 
   // boost to the PL lab frame.
   double vx = PL.Px()/PL.E();
   double vy = PL.Py()/PL.E();
   double vz = PL.Pz()/PL.E();
   Ppfin.Boost(vx,vy,vz);
   Ppim.Boost (vx,vy,vz);
	//cout << " PPLAB " <<  Ppfin.Vect().Mag() << " " <<  Ppim.Vect().Mag() << " " << endl;

 	//cout << " PPLAB " <<  Ppfin.Vect().Z() << " " <<  Ppim.Vect().Z() << " " << endl;
  


   //Ppfin.SetXYZT(0.1,0.1,0.1,(mL*mL+massProton2-massPion2)/2./mL);
   //Ppim.SetXYZT (0.1,0.1,0.1,(mL*mL+massPion2-massProton2)/2./mL);
   //Ppfin.Boost(vx,vy,vz);
   //Ppim.Boost (vx,vy,vz);
 
   //cout << " PPBOOST " <<  Ppfin.Vect().Mag() << " " <<  Ppim.Vect().Mag() << " " << endl;

}


void getSdecayProd(TLorentzVector &PS, TLorentzVector &PL,
                   TLorentzVector &Ppfin, TLorentzVector &Ppim, TLorentzVector &Pgam) {

   //4-momenta PL and Pgam in the PS rest frame
   //Set energy components 
   double mS = massSigma0;
   Pgam.SetXYZT(1.,1.,1.,(mS*mS+0*0        -massLambda2)/2./mS);
   PL.SetXYZT  (1.,1.,1.,(mS*mS+massLambda2-0*0        )/2./mS);
   
   //Set momentum muduli and angles
   Pgam.SetRho(sqrt(Pgam.E()*Pgam.E()-0.));
   PL.SetRho  (sqrt(PL.E()*PL.E()    -massLambda2));
   double cost = randomIntv(-1., 1.);
   double thetagam = acos(cost);
   double phigam   = randomIntv(0.0, 2.*constantPi);
   double thetaL   = constantPi-thetagam;
   double phiL     = phigam+constantPi; if(phiL>2.*constantPi) phiL=phiL-2.*constantPi;
   Pgam.SetTheta(thetagam);
   Pgam.SetPhi(phigam);
   PL.SetTheta(thetaL);
   PL.SetPhi(phiL);
    
   // boost to the PL lab frame.
   double vx = PS.Px()/PS.E();
   double vy = PS.Py()/PS.E();
   double vz = PS.Pz()/PS.E();
   Pgam.Boost(vx,vy,vz);
   PL.Boost  (vx,vy,vz);

   //prnLV(" PS         ",   PS);
   //prnLV(" Pg+PL      ", (Pgam+PL));

   getLdecayProd(PL, Ppfin, Ppim);

   //prnLV(" Pg+Ppr+Ppi ", (Pgam+Ppfin+Ppim)); 
   //cout <<" MMM "<< PS.M() <<" "<<PL.M()<<" "<<Pgam.M()<<" "<<Ppfin.M()<<" "<<Ppim.M()<< endl; 
   
 
}





// Calculate phi in sector frame and sector number 
// from phi in lab frame
// Angles are in degrees!!!
double getPhiSectDeg(double phi, int &sect) {

  double phiSect;

  if     ( phi >= (-30)  && phi <= (30)  ) { phiSect = phi;     sect=1; }
  else if( phi >= (30)   && phi <= (90)  ) { phiSect = phi-60;  sect=2; }
  else if( phi >= (90)   && phi <= (150) ) { phiSect = phi-120; sect=3; }
  else if( phi >= (150)  && phi <= (180) ) { phiSect = phi-180; sect=4; }
  else if( phi >= (-90)  && phi <= (-30) ) { phiSect = phi+60;  sect=6; }
  else if( phi >= (-150) && phi <= (-90) ) { phiSect = phi+120; sect=5; }
  else if( phi >= (-180) && phi <= (-150)) { phiSect = phi+180; sect=4; }
  else { 
    cerr << " Error! getPhiSectDeg(...). wrong phi = " << phi << endl;
    exit(1);
  } 
  if( phiSect < (-30) || phiSect > (30) ) {
    cerr << " Error! getPhiSectDeg(...). wrong phiSect = " << phiSect << endl;
    exit(1);
  } 

  return phiSect;

}


#endif

