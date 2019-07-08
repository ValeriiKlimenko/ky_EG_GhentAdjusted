#ifndef _SIGMAKY_H
#define _SIGMAKY_H

#include "string.h"

#include "../generalUtils/utils.h"
#include "../generalUtils/constants.h"
#include "../generalUtils/kinematics.h"


class sigmaKY {
	int varsnum;
	int obsnum;
	char vars[7][20];
	char vardims[7][20];
	int count_Q2;
	int count_W;
	int count_cos_t;
	char eps_or_Eb;
	char teta_or_cost;
	char rads_or_degrees;
	double beam_energy;
	double step_Q2;
	double step_W;
	double step_cos_t;
	
	int maxIQ2 = 61;
	int maxIW = 35;
	int maxIcos_t = 21;	
	double Q2[61];	
	double W[35];	
	double cos_t[21];
	double L[61][35][21];
	double T[61][35][21];
	double LT[61][35][21];
	double TT[61][35][21];
	double SINT[61][35];
	double SINTint[61];	
	double SINText[61];		
	double sintQ2[856];
	double sintXS[856];
	double factQ2[61][35];
	
	double deltaSS;
	string correctQ2;
	
	

	double dsigmaXX_dphi(double XX[61][35][21], int iQ2, double _Q2, int iW, double _W);					//complete
	double dsigmaXX_domega(double XX[61][35][21], int iQ2, double _Q2, int iW, double _W, int icos_t, double _cos_t);					//complete
	int srch_Q2(double value);
	int srch_W(double value);
	int srch_cos_t(double value);
	double fxxx(int i, int j, int k, double eps, double phi);
	double fxxx(int i, int j, int k, double eps);
public:
	sigmaKY();
	sigmaKY(string dataPath, string path, string, double, double, double);
	//~sigmaKY(){ delete hres; }

			//switches
	void use_beam_energy()  {eps_or_Eb = 1;}
	void use_epsilon()	{eps_or_Eb = 0;}
	void use_theta()	{teta_or_cost = 1;}
	void use_costeta()	{teta_or_cost = 0;}
	void use_radians()	{rads_or_degrees = 1;}
	void use_degrees()	{rads_or_degrees = 0;}
	double StepW() {return step_W;}
	double StepQ2() {return step_Q2;}
	double MinW() {return W[0];}
	double MinQ2() {return Q2[0];}
	double MaxW() {return W[count_W-1];}
	double MaxQ2() {return Q2[count_Q2-1];}
	double BeamEnergy() {return beam_energy;}
	int get_count_cos_t()	{return count_cos_t;}
	double *get_cos_t() {return cos_t;}
	void setBeamEnergy(double _beam_energy) {beam_energy = _beam_energy;}

        double d2sigma(double Ebeam, double Q2, double W, 
	               double thetaK, double phiK);

	double d5sigma(double beam_energy, double _Q2, double _W, 
	               double thetaK, double phiK);

        double d5sigma_max(double Ebeam, double Q2min, double Q2max,
                           double Wmin,  double Wmax  );

	double dsigma_dcos(double _beam_energy, double _Q2, double _W, double teta);					//complete

        double sigma_int(double Ebeam, double Q2, double W);

	double *dsigma_dcos(double _beam_energy, double _Q2, double _W);								//complete

	double dsigma_dphi(double _beam_enegry, double _Q2, double _W, double phi);					//complete but perhaps incorrectly

	double *dsigma_dphi(double _beam_enegry, double _Q2, double _W);								//complete

	double dsigma_domega(double _beam_energy, double _Q2, double _W, double teta, double phi);	//complete

	double *dsigma_domega(double _beam_energy, double _Q2, double _W, double teta);				//complete

	double sigma(double _beam_energy, double _Q2, double _W);

	double **sigma_grid(double _beam_energy);
};





sigmaKY::sigmaKY(string dataPath, string name, string addr="", 
                 double a12=0., double a32=0., double s12=0.)
{

	correctQ2 = "no";

        string KLfName = "work_K+Lambda0_strange_calc_RPR-2011.txt";
        string KSfName = "work_K+Sigma0_strange_calc_RPR-2007.txt";
	string KLint   = "KL_int.txt";
	string KSint   = "KS_int.txt";


        ifstream input;
        ifstream input_int;
	
        if(name=="KLambda") {
	  input.open(dataPath+"/"+KLfName);
	  input_int.open(dataPath+"/"+KLint);
	} else if(name=="KSigma") {
	  input.open(dataPath+"/"+KSfName);
	  input_int.open(dataPath+"/"+KSint);
	} else {
	  input.open(name.c_str());
	}
	

	this->use_beam_energy();
	this->use_costeta();
	this->use_radians();

	varsnum = 7;
	obsnum = 44835;
	count_Q2 = 61;
	count_W = 35;
	count_cos_t = 21;


	if (!input)	{
		std::cout << "Unable to open file " << name << '\n';
		exit(1);
	}
	
	
	double _Q2, _W, _cos_t;
	char junk[20];
	for(int i = 0; i < varsnum; i++) { input >> vars[i]; }
	for(int i = 0; i < varsnum; i++) input >> vardims[i];
	for(int i = 0; i < varsnum; i++) input >> junk;
	
	for(int i = 0; i < count_Q2; i++)	{
		for(int j = 0; j < count_W; j++)	{
			for(int k = 0; k < count_cos_t; k++)	{
				input >> _Q2 >> _W >> _cos_t >> L[i][j][k] >> T[i][j][k] >> LT[i][j][k] >> TT[i][j][k];
				if (i == 0 && j ==0) cos_t[k] = _cos_t;
			}
			if(i == 0) W[j] = _W;
		}
		Q2[i] = _Q2;
	}	
	input.close();
	
	
	string tmp;
	getline(input_int, tmp);
	getline(input_int, tmp);
	for(int i=0; i<856; i++) {
	  input_int >> sintQ2[i] >> sintXS[i];
	}
	input_int.close();

	step_Q2 = (Q2[count_Q2-1] - Q2[0]) / (double)(count_Q2-1);
	step_W = (W[count_W-1] - W[0]) / (double)(count_W-1);
	step_cos_t = (cos_t[count_cos_t-1] - cos_t[0]) / (double)(count_cos_t-1);
		



	// calculate integrated cross section
	double sint1 = sigma_int(11, 4, 2);
	for(int iQ2=0; iQ2 < maxIQ2; iQ2++) {
	SINTint[iQ2]=0;
	//if(Q2[iQ2] == 3.6) cout << " " << iQ2 << " "<< Q2[iQ2] << endl;
	for(int iW=0; iW < maxIW; iW++) {
	
	   //cout << " ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZz " << iQ2 <<" "<<iW <<" " << endl;
	
	
           SINT[iQ2][iW] = sigma_int(8, Q2[iQ2], W[iW]);
	   
	   //cout << " ZZZ " << iQ2 <<" "<<iW <<" "<< SINT[iQ2][iW] << endl;
	   
	   SINTint[iQ2] += SINT[iQ2][iW];	   
	}
	   //cout << " NNN " << iQ2 <<" "<< SINTint[iQ2] << endl;

	
	
	}	
	for(int iQ2=0; iQ2 < maxIQ2; iQ2++) {
	  for(int iext=0; iext<(856-1); iext++) {
	  if(Q2[iQ2] >= sintQ2[iext] && Q2[iQ2] <= sintQ2[iext+1]) {
	    SINText[iQ2] = sintXS[iext] 
	                 +  (sintQ2[iext+1]-sintQ2[iext])
			   *(sintXS[iext+1]-sintXS[iext])/(sintQ2[iext+1]-sintQ2[iext]);
	    //cout << " SSS+- " << sintXS[iext] <<" "<< sintXS[iext+1] << " " <<  SINText[iQ2] << endl;
	    //cout << " SSS " << iQ2 <<" "<<  Q2[iQ2] <<" "<< sintQ2[iext] <<" "<< sintQ2[iext+1] 
	    //     <<" "<< SINText[iQ2] << endl;
	  }
	  }
	  
	}
	
	
	for(int iQ2=0; iQ2 < maxIQ2; iQ2++) {
        if(iQ2 >= 18) {
          double sig0    = SINTint[18];
          double sigext0 = SINText[18];
	  double sig     = SINTint[iQ2];
          double sigext  = SINText[iQ2];
  	  double fact = 1.;
	   
	    //cout << " QQQQ2 " << iQ2 <<" "<< Q2[iQ2] << " | "
	    //     << " "<< sigext0 <<" "<< sigext << " " << sig0 << " "<< sig << " | "
	    //  <<" "<< sigext/sigext0 <<" "<< sig/sig0 << endl;

	  if(sigext>0. && sigext0>0. && sig0>0. && sig>0.) {
	    	  
	    //cout << " QQ2 " << iQ2 <<" "<< Q2[iQ2] << " | "
	    //     << " "<< sigext0 <<" "<< sigext << " " << sig0 << " "<< sig << " | "
            //     <<" "<< sigext/sigext0 <<" "<< sig/sig0 << endl;
	    fact = (sigext/sigext0)*(sig0/sig);
	    
	  }
	    //cout << " FFF " <<  fact << endl;
	  	
	}
	}
	
	
	

	for(int iQ2=0; iQ2 < maxIQ2; iQ2++) {
	for(int iW=0; iW < maxIW; iW++) {
  	 factQ2[iQ2][iW] = 1.;
        if(iQ2 >= 18) {
          double sig0    = SINT[18][iW];
          double sigext0 = SINText[18];
	  double sig     = SINT[iQ2][iW];
          double sigext  = SINText[iQ2];
  	  factQ2[iQ2][iW] = 1.;

	  if(sigext>0. && sigext0>0. && sig0>0. && sig>0.) {
	    factQ2[iQ2][iW] = (sigext/sigext0)*(sig0/sig);
	    	  
	    //cout << " Q2W " << iQ2 <<" "<< iW << " " << Q2[iQ2] << " "<< W[iW]  << " | "
	    //     << " "<< sigext0 <<" "<< sigext << " " << sig0 << " "<< sig << " | "
		// <<" "<< sigext/sigext0 <<" "<< sig/sig0 << endl;
	    
	  }
	  	  
	  //if(factQ2[iQ2][iW] > 1.1) cout << " RRRRR " << endl;
	  //cout << " FFF " << factQ2[iQ2][iW] << endl;
	  	
	}
	}
	}
		
		
	
	
	/*
	double sint2 = 0;
	double dcos = 0.05;
	double dphi = 0.05;	
	for(double cost=-1; cost<1; cost+=dcos) {
	for(double phi=0; phi<constantPi*2; phi+=dphi) {
	   sint2 += dsigma_domega(11, 4, 2, cost, phi)*dcos*dphi;
	
	}
	}	
	*/	
	//cout << "sint1 " << sint1 <<" " << endl;
	
	
	correctQ2 = "yes";
	
	
	
}



int sigmaKY::srch_Q2(double value)
{
	if (value <= Q2[count_Q2-1] && value >= Q2[0])	{
		int pre_index = (int)((value - Q2[0]) / step_Q2);
		int index = pre_index;
		if (Q2[pre_index] <= value && Q2[pre_index+1] > value) index = pre_index;
			else if (Q2[pre_index] > value) index =  pre_index-1;
				else if (Q2[pre_index+1] <= value) index = pre_index+1;
		if(index < 0) index=0;
		if(index >= maxIQ2) index=maxIQ2-1;			
		return index;
				
	} else {
		std::cout << "Select the value of Q^2 in the interval [" << Q2[0] << ", " << Q2[count_Q2-1] << "] GeV^2\n";
		cout << "Q2 " << value << endl;
		exit(1);
	}
}
int sigmaKY::srch_W(double value)
{
        if(value<W[0]) return -1;
	if (value <= W[count_W-1] && value >= W[0])	{
		int pre_index = (int)((value - W[0]) / step_W);
		int index = pre_index;
		if (W[pre_index] <= value && W[pre_index+1] > value) index = pre_index;
			else if (W[pre_index] > value) index =  pre_index-1;
				else if (W[pre_index+1] <= value) index = pre_index+1;
		if(index < 0) index=0;
		if(index >= maxIW) index=maxIW-1;			
		return index;
		
				
	} else {
		std::cout << "Select the value of W in the interval [" << W[0] << ", " << W[count_W-1] << "] GeV\n";
		exit(1);
	}
}
int sigmaKY::srch_cos_t(double value)
{
	if (value <= cos_t[count_cos_t-1] && value >= cos_t[0])	{
		int pre_index = (int)((value - cos_t[0]) / step_cos_t);
		int index = pre_index;
		if (cos_t[pre_index] <= value && cos_t[pre_index+1] > value) index = pre_index;
			else if (cos_t[pre_index] > value) index =  pre_index-1;
				else if (cos_t[pre_index+1] <= value) index = pre_index+1;
		if(index < 0) index=0;
		if(index >= maxIcos_t) index=maxIcos_t-1;			
		return index;
	} else {
		std::cout << "Select the value of cos(teta) in the interval [-1, 1]\n";
		exit(1);
	}
}




double sigmaKY::fxxx(int i, int j, int k, double eps, double phi)	{

        //cout << " TTT1 " << T[i][j][k] <<" "<<  L[i][j][k] << " " << TT[i][j][k] << " " << LT[i][j][k] << endl;


	return T[i][j][k] + eps * L[i][j][k] + eps * TT[i][j][k] * cos(2*phi) + sqrt(eps*(1+eps)) * LT[i][j][k] * cos(phi);


}
double sigmaKY::fxxx(int i, int j, int k, double eps)	{


        //cout << " III " << i <<" "<<  j <<" "<< k << endl;
        //cout << " TTT2 " << T[i][j][k] <<" "<<  L[i][j][k] << endl;

	return T[i][j][k] + eps * L[i][j][k];

}




double sigmaKY::dsigma_domega(double _beam_energy, double _Q2, double _W, double teta, double phi)
{
	double _cos_t;
	if (teta_or_cost) _cos_t = cos(teta);
		else _cos_t = teta;
		
	int iQ2 = this->srch_Q2(_Q2);
	int iW = this->srch_W(_W);
	if(iW==-1) return 0.;
	int icos_t = this->srch_cos_t(_cos_t);
	double nu = (_Q2 + _W*_W - massProton2) / (2 * massProton);
	double eps;
	if (eps_or_Eb) eps = 1 / (1 + 2*((_Q2 + nu*nu) / (4*(_beam_energy - nu)*_beam_energy - _Q2)));
	else eps = _beam_energy;

	double _const = step_Q2 * step_W * step_cos_t;
	double interpolated_value = 0;

	interpolated_value += this->fxxx(iQ2,iW,icos_t,eps,phi)/_const * (Q2[iQ2+1] - _Q2) * (W[iW+1] - _W) * (cos_t[icos_t+1] - _cos_t);
	interpolated_value += this->fxxx(iQ2,iW,icos_t+1,eps,phi)/_const * (Q2[iQ2+1] - _Q2) * (W[iW+1] - _W) * (_cos_t - cos_t[icos_t]);
	interpolated_value += this->fxxx(iQ2,iW+1,icos_t,eps,phi)/_const * (Q2[iQ2+1] - _Q2) * (_W - W[iW]) * (cos_t[icos_t+1] - _cos_t);
	interpolated_value += this->fxxx(iQ2,iW+1,icos_t+1,eps,phi)/_const * (Q2[iQ2+1] - _Q2) * (_W - W[iW]) * (_cos_t - cos_t[icos_t]);
	interpolated_value += this->fxxx(iQ2+1,iW,icos_t,eps,phi)/_const * (_Q2 - Q2[iQ2]) * (W[iW+1] - _W) * (cos_t[icos_t+1] - _cos_t);
	interpolated_value += this->fxxx(iQ2+1,iW,icos_t+1,eps,phi)/_const * (_Q2 - Q2[iQ2]) * (W[iW+1] - _W) * (_cos_t - cos_t[icos_t]);
	interpolated_value += this->fxxx(iQ2+1,iW+1,icos_t,eps,phi)/_const * (_Q2 - Q2[iQ2]) * (_W - W[iW]) * (cos_t[icos_t+1] - _cos_t);
	interpolated_value += this->fxxx(iQ2+1,iW+1,icos_t+1,eps,phi)/_const * (_Q2 - Q2[iQ2]) * (_W - W[iW]) * (_cos_t - cos_t[icos_t]);


        if(correctQ2=="yes") {
        if(iQ2 >= 18) {
          double sig0    = SINTint[18]; //SINT[18][iW];
          double sigext0 = SINText[18];
	  double sig     = SINTint[iQ2]; //SINT[iQ2][iW];
          double sigext  = SINText[iQ2];
  	  double fact = 1.;

	  if(sigext>0. && sigext0>0. && sig0>0. && sig>0.) {
	    fact = (sigext/sigext0)*(sig0/sig);
	  }
	  //if(sig<1.E-5) fact=1.;
	  
	  if(fact!=1) {
	  //cout << " iQ2 iW " << iQ2 <<" "<< iW << " " << sig <<" "<< sigext << endl;
	  //cout << " ssss " << sig0 <<" "<< sigext0 << " " << endl;
	  //cout << " ssss " << sig <<" "<< sigext << " " << endl;
	  //cout << " Q2 FAC " << _Q2 <<" "<< fact << " " << endl;
	  }

	  interpolated_value *= fact;  
        
	}
        }

	return interpolated_value;
}



double sigmaKY::d5sigma(double Ebeam, double Q2, double W, double thetaK, double phiK)
{
  double mp=massProton;
  double mp2 = mp*mp;
  double pi=constantPi;
  double alpha=constantAlpha;
  double W2 = W*W;
  
  double omega = (W2 + Q2 -mp2)/(2.*mp);
  if(omega<0. || omega>Ebeam) return 0.;
  double eps2 = Ebeam - omega;
  double k2 = Q2 + omega*omega;
  if(k2<0) return 0.;
  double arg = 1. - Q2/(2.*Ebeam*eps2);
  if(arg>1. || arg<-1.) return 0.;
  double theta = 2.*acos(arg);
  double epsilon = 1./( 1 + (2.*k2/Q2)*tan(theta/2.) );
  double gamma = getGamma(Ebeam, Q2, W);
    
  double d5sig = gamma *
                 this->dsigma_domega(Ebeam, Q2, W, cos(thetaK), phiK);
  
  return d5sig;
  
  

}



double sigmaKY::d5sigma_max(double Ebeam, double Q2min, double Q2max,
                                       double Wmin,  double Wmax  )
{

        double pi=constantPi;

	// Find maximum of the cross section
        int nQ2=100;
        int nW=20;
        int nCosThetaK=20;
        int nPhiK=20;
	double d5sigmaMax=0.;
        for(int iQ2=0; iQ2<nQ2; iQ2++) {
	double Q2 = Q2min + (Q2max-Q2min)*iQ2/(nQ2-1);
        for(int iW=0;   iW<nW;  iW++) {
	double W =  Wmin +  (Wmax-Wmin)*iW/(nW-1);
	for(int iCosThK=0; iCosThK<20; iCosThK++) {
	double cosThetaK = -0.9999 + (0.9999-(-0.9999))*iCosThK/(nCosThetaK-1);
	double thetaK = acos(cosThetaK);
        for(int iPhiK=0; iPhiK<nPhiK;  iPhiK++) {
	double phiK = 0. + (2.*pi-0.)*iPhiK/(nPhiK-1);
          double d5sig = d5sigma(Ebeam, Q2, W, thetaK, phiK);

          if(d5sig>d5sigmaMax) d5sigmaMax=d5sig;
    
	}
        }
	}
        }
        return d5sigmaMax;

}







double sigmaKY::dsigma_dcos(double _beam_energy, double _Q2, double _W, double teta)
{
	double _cos_t;
	if (teta_or_cost) _cos_t = cos(teta);
        	else _cos_t = teta;
	int iQ2 = this->srch_Q2(_Q2);
	int iW = this->srch_W(_W);
	int icos_t = this->srch_cos_t(_cos_t);
	if(iW==-1) return 0.;
	double nu = (_Q2 + _W*_W - massProton2) / (2 * massProton);
	double eps;
	if (eps_or_Eb) eps = 1 / (1 + 2*((_Q2 + nu*nu) / (4*(_beam_energy - nu)*_beam_energy - _Q2)));
		else eps = _beam_energy;
	//cout << " dddddddddddddddddddd " << iQ2 <<" "<< iW <<" "<< icos_t <<" "<< _W << endl;
	double _const = (Q2[iQ2+1] - Q2[iQ2]) * (W[iW+1] - W[iW]) * (cos_t[icos_t+1] - cos_t[icos_t]);
	
	//cout << " ddd " <<  (Q2[iQ2+1] - Q2[iQ2]) <<" "<< (W[iW+1] - W[iW]) <<" "<< (cos_t[icos_t+1] - cos_t[icos_t]) << endl;
	//cout << " ddd " <<  iW <<" "<< W[iW+1] <<" "<< W[iW]  << endl;

	
	double interpolated_value = 0;




	interpolated_value += this->fxxx(iQ2,iW,icos_t,eps)/_const * (Q2[iQ2+1] - _Q2) * (W[iW+1] - _W) * (cos_t[icos_t+1] - _cos_t);
	interpolated_value += this->fxxx(iQ2,iW,icos_t+1,eps)/_const * (Q2[iQ2+1] - _Q2) * (W[iW+1] - _W) * (_cos_t - cos_t[icos_t]);
	interpolated_value += this->fxxx(iQ2,iW+1,icos_t,eps)/_const * (Q2[iQ2+1] - _Q2) * (_W - W[iW]) * (cos_t[icos_t+1] - _cos_t);
	interpolated_value += this->fxxx(iQ2,iW+1,icos_t+1,eps)/_const * (Q2[iQ2+1] - _Q2) * (_W - W[iW]) * (_cos_t - cos_t[icos_t]);
	interpolated_value += this->fxxx(iQ2+1,iW,icos_t,eps)/_const * (_Q2 - Q2[iQ2]) * (W[iW+1] - _W) * (cos_t[icos_t+1] - _cos_t);
	interpolated_value += this->fxxx(iQ2+1,iW,icos_t+1,eps)/_const * (_Q2 - Q2[iQ2]) * (W[iW+1] - _W) * (_cos_t - cos_t[icos_t]);
	interpolated_value += this->fxxx(iQ2+1,iW+1,icos_t,eps)/_const * (_Q2 - Q2[iQ2]) * (_W - W[iW]) * (cos_t[icos_t+1] - _cos_t);
	interpolated_value += this->fxxx(iQ2+1,iW+1,icos_t+1,eps)/_const * (_Q2 - Q2[iQ2]) * (_W - W[iW]) * (_cos_t - cos_t[icos_t]);


	/*
	//cout << " uuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuu " << interpolated_value << endl;
	//cout << " isv1 " << this->fxxx(iQ2,iW,icos_t,eps)/_const * (Q2[iQ2+1] - _Q2) * (W[iW+1] - _W) * (cos_t[icos_t+1] - _cos_t) << endl;
	//cout << " vvv " <<  this->fxxx(iQ2,iW,icos_t,eps) << " " << _const << endl;
	//cout << " vvv " << (Q2[iQ2+1] - _Q2) <<" "<< (W[iW+1] - _W) << " "<< (cos_t[icos_t+1] - _cos_t) <<" " << endl;
	//cout << " vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv " <<  endl;
	
	cout << " isv2 " << this->fxxx(iQ2,iW,icos_t+1,eps)/_const * (Q2[iQ2+1] - _Q2) * (W[iW+1] - _W) * (_cos_t - cos_t[icos_t]) << endl;
	cout << " isv3 " << this->fxxx(iQ2,iW+1,icos_t,eps)/_const * (Q2[iQ2+1] - _Q2) * (_W - W[iW]) * (cos_t[icos_t+1] - _cos_t) << endl;
	cout << " isv4 " << this->fxxx(iQ2,iW+1,icos_t+1,eps)/_const * (Q2[iQ2+1] - _Q2) * (_W - W[iW]) * (_cos_t - cos_t[icos_t]) << endl;
	cout << " isv5 " << this->fxxx(iQ2+1,iW,icos_t,eps)/_const * (_Q2 - Q2[iQ2]) * (W[iW+1] - _W) * (cos_t[icos_t+1] - _cos_t) << endl;
	cout << " isv6 " << this->fxxx(iQ2+1,iW,icos_t+1,eps)/_const * (_Q2 - Q2[iQ2]) * (W[iW+1] - _W) * (_cos_t - cos_t[icos_t]) << endl;
	cout << " isv7 " << this->fxxx(iQ2+1,iW+1,icos_t,eps)/_const * (_Q2 - Q2[iQ2]) * (_W - W[iW]) * (cos_t[icos_t+1] - _cos_t) << endl;
	cout << " isv8 " << this->fxxx(iQ2+1,iW+1,icos_t+1,eps)/_const * (_Q2 - Q2[iQ2]) * (_W - W[iW]) * (_cos_t - cos_t[icos_t]) << endl;


	//cout << " ddddd " << interpolated_value << endl;
	*/

	return 2 * M_PI * interpolated_value; // * (389.5*(sigt + eps*sigl)); // * hres->dsigma_dcos(getEbeamEps(eps, _Q2, _W), _Q2, _W, teta);
}


double sigmaKY::sigma_int(double Ebeam, double Q2, double W)
{
   double sig = 0.;
   int np = 41;
   double dC=2./(np-1.);
   for(int iC=0; iC<(np-1); iC++) {
     double ct = -1. + dC/2 + iC*dC;
     
     sig += dsigma_dcos(Ebeam, Q2, W, ct)*dC;
   
     //cout << " DDDDD " << dsigma_dcos(Ebeam, Q2, W, ct)*dC << endl;
   
   }
   return sig;   
}







#endif 
