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

#include "../generalUtils/utils.h"
#include "../generalUtils/constants.h"
#include "../generalUtils/kinematics.h"
#include "evGenerator.h"


using namespace std;


int main(int argc, char *argv[])
{

        int channel;
	string channelName, outputFileName, decMode, dataPath;
	double Ebeam, Q2min, Q2max, Wmin, Wmax;
	int nEvents;
        double jr, mr, gr, a12, a32, s12, onlyres;

	cout << "\nEvent generator started. " <<  endl;
	
	if( argc != 2 ) {
	  cerr << " Usage:  ./eg_ky <configuration file name>" << endl;
	  cerr << " STOP!" << endl;
	  return 1;
	}	
	  

	// read config. file
	ifstream input;
	string fName = argv[1];
	input.open(fName.c_str());
	if(!input.good()) {
	  cerr << " eg_ky.cpp: can not open file: " << fName.c_str() << endl;
	  cerr << " STOP!" << endl;
	  return 1;
	}
	std::getline(input,channelName); 
	cout << " Configuration file is " << fName << endl;
	cout << " Channel is " << channelName << endl;
	input >> Ebeam >> Q2min >> Q2max >> Wmin >> Wmax >> nEvents;
	cout << " Ebeam is " << Ebeam << endl;
	cout << " Q2min is " << Q2min << endl;
	cout << " Q2max is " << Q2max << endl;
	cout << " Wmin is " << Wmin << endl;
	cout << " Wmax is " << Wmax << endl;
	cout << " nEvents is " << nEvents << endl;
	input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	std::getline(input,decMode); decMode.erase(remove(decMode.begin(),decMode.end(),' '),decMode.end());
	cout << " decMode is " << decMode << endl;
	std::getline(input,outputFileName); outputFileName.erase(remove(outputFileName.begin(),outputFileName.end(),' '),outputFileName.end());
	cout << " outputFileName is " << outputFileName << endl;
	std::getline(input,dataPath); dataPath.erase(remove(dataPath.begin(),dataPath.end(),' '),dataPath.end());
	cout << " dataPath is " << dataPath << endl;
	input.close();


        // PARAMETERS
	//channelName = "KLambda"; //"KLambda" or "KSigma"
	//Ebeam = 6.6;
	//Q2min = 0.05;
	//Q2max = 2.0;
	//Wmin = 1.5;
	//Wmax = 3.0;
        //nEvents = 10000;
	//outputFileName = "output.lund";
	//dataPath = "./data";
	
	

	// initilize event generator
	evGenerator eg(dataPath, channelName, Ebeam,  Q2min, Q2max, Wmin, Wmax);
        cout << " evGenerator has been initialized\n" << endl;
	     
 
        // 4-momenta of final electron, Kaon and Lambda, proton, pi- in LAB frame
	TLorentzVector Peini( TVector3(0.,0.,Ebeam), sqrt(Ebeam*Ebeam+0.005*0.005) );
	TLorentzVector Ppini( TVector3(0.,0.,0.),    0.9383);
	TLorentzVector Pefin, PK, PL, Ppfin, Ppim, Pgam;
	double Q2,W;
	// output  
	ofstream output(outputFileName.c_str()); 
	// Loop through events
        for (int i=0; i<nEvents; i++) {	
		
		
	  // get event. i.e. get 4-momenta of final state particle. 
	  // Values of Q2 and W are also returned.
	  // for K+Sigma events Pgam is a decay gamma, 
          eg.getEvent(Q2, W, Pefin, PK, PL, Ppfin, Ppim, Pgam);


	  int nParticles = 0;
	  if(decMode=="nodec") 
	     { nParticles = 3; } 
	  else if(decMode=="dec") 
	     { 
	        if(channelName=="KLambda") nParticles = 4; 
		else if(channelName=="KSigma") nParticles = 5;  }
	  else {
	    cerr << " eg_ky.cpp: Wrong decay mode: " << decMode << endl;
	    cerr << " STOP!" << endl;
	    return 1;
	  }


          // OUTPUT in LUND format
	  // header
	  output << nParticles << " 1 1 0 0 0 0 "
	         <<" "<< W <<" "<< Q2 <<" "<< getomega(Q2, W) 
		 << endl;
	  // electron
	  output 
	    << "1 -1 1 " << lundIdElectron << " 0 0 "
	    << Pefin.Px() <<" "<< Pefin.Py() <<" "<< Pefin.Pz()  <<" "<< Pefin.E() <<" 0.0005"
	    <<" 0 0 0 "
	    << endl;
	  // Kaon+
	  output 
	    << "2 1 1 " << lundIdKaonPlus << " 0 0 " << PK.Px() <<" "<< PK.Py() <<" "<< PK.Pz() 
	    <<" "<< PK.E() <<" 0.4936"
	    <<" 0 0 0 "
	    << endl;
	  if(channelName=="KLambda"){
	    // Lambda
	    if(decMode=="nodec") {
	       output 
	   	 << "3 0 1 " << lundIdLambda << " 0 0 "
	   	 << PL.Px() <<" "<< PL.Py() <<" "<< PL.Pz() 
	   	 <<" "<< PL.E() <<" 1.115"
	   	 <<" 0 0 0 "
	   	 << endl;
	    } else {
	      output 
	   	<< "3 1 1 " << lundIdProton << " 0 0 "
	   	<< Ppfin.Px() <<" "<< Ppfin.Py() <<" "<< Ppfin.Pz() 
	   	<<" "<< Ppfin.E() <<" 0.9383"
	   	<<" 0 0 0 "
	   	<< endl;
	      output 
	   	<< "4 -1 1 " << lundIdPiMinus << " 0 0 "
	   	<< Ppim.Px() <<" "<< Ppim.Py() <<" "<< Ppim.Pz() 
	   	<<" "<< Ppim.E() <<" 0.1396"
	   	<<" 0 0 0 "
	   	<< endl;
            }
	  } else if(channelName=="KSigma") {
	  // Sigma zero
	    if(decMode=="nodec") {
	       output 
	   	 << "3 0 1 " << lundIdSigmaZero << " 0 0 "
	   	 << PL.Px() <<" "<< PL.Py() <<" "<< PL.Pz() 
	   	 <<" "<< PL.E() <<" 1.115"
	   	 <<" 0 0 0 "
	   	 << endl;
	    } else {
	      output 
	   	<< "3 1 1 " << lundIdProton << " 0 0 "
	   	<< Ppfin.Px() <<" "<< Ppfin.Py() <<" "<< Ppfin.Pz() 
	   	<<" "<< Ppfin.E() <<" 0.9383"
	   	<<" 0 0 0 "
	   	<< endl;
	      output 
	   	<< "4 -1 1 " << lundIdPiMinus << " 0 0 "
	   	<< Ppim.Px() <<" "<< Ppim.Py() <<" "<< Ppim.Pz() 
	   	<<" "<< Ppim.E() <<" 0.1396"
	   	<<" 0 0 0 "
	   	<< endl;
	      output 
	   	<< "5  0 1 " << lundIdGamma << " 0 0 "
	   	<< Pgam.Px() <<" "<< Pgam.Py() <<" "<< Pgam.Pz() 
	   	<<" "<< Pgam.E() <<" 0."
	   	<<" 0 0 0 "
	   	<< endl;
            }
	  } else {
	    cerr << " eg_ky.cpp: Wrong channel name " << channel << endl;
	    cerr << " STOP!" << endl;
	  }
	  
	  // test 
	  //TLorentzVector PfinSum = Pefin + PK + Ppfin + Ppim;
	  //cout << " Pfin " << PfinSum.E() <<" "<< PfinSum.Px() <<" "<< PfinSum.Py() <<" "<< PfinSum.Pz()
	  //     << endl; 
	  //cout << " miss 0 "  << (Peini + Ppini - Pefin - PK - Ppfin - Ppim).M() << endl;
	  //cout << " miss e "  << (Peini + Ppini - PK - Ppfin - Ppim).M() << endl;
	  //cout << " miss p "  << (Peini + Ppini - Pefin - PK - Ppim).M() << endl;
	  //cout << " miss K "  << (Peini + Ppini - Pefin - Ppfin - Ppim).M() << endl;
	  //cout << " miss pi " << (Peini + Ppini - Pefin - PK - Ppfin).M() << endl;
	  //cout << " Mass L "  << (Ppfin + Ppim).M() << endl;
	  
	  
	  if( int((i+1)/1000)*1000 == (i+1) | (i+1)==1) {
	    cout << " Event # " << i+1 << endl; 
	  }
	}
	output.close();
	

	return 0;
	
	
}
