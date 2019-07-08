
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

//ROOT LIBs
#include "TLorentzVector.h"
#include "TApplication.h"
#include "TRint.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TSpline.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TPad.h"

#include "../generalUtils/utils.h"
#include "../generalUtils/constants.h"
#include "../generalUtils/kinematics.h"
#include "sigmaKY.h"
#include "evGenerator.h"


using namespace std;




int main(int argc, char **argv)
{

	TRint* app = new TRint ("Rint", &argc, argv);

	
	double Ebeam = 6.6;
	double Q2min = 1.00;
	double Q2max = 1.0001;
	double Wmin = 2.00;
	double Wmax = 2.001;
        int nEventMax = 100000;


	double Q2avr = 0.5*(Q2max+Q2min);
	double Wavr  = 0.5*(Wmax+Wmin);
	double epsAvr = getEpsilon(Ebeam, Q2avr, Wavr);
	double Q2,W,cosThetaK,thetaK,phiK,d5sigma;
         

	// initilize event generator
        evGenerator eg("./data","KLambda", Ebeam,  Q2min, Q2max, Wmin, Wmax);
		
	//histogram for events
        TH1F *hTh = new TH1F("qq","#gamma_{v} p #rightarrow #Lambda^{0} K^{+}, Q^{2}=0.65, W=1.85, #varepsilon=0.54",240,-1.,1.);
	hTh->GetXaxis()->SetTitle("cos(#theta)");
	hTh->GetYaxis()->SetTitle("arb. units.");

        // 4-momenta of final electron, Kaon and Lambda LAB frame
	TLorentzVector Pefin, PK, PY, Ppfin, Ppim, Pgam;
	// Loop through events
        for (int i=0; i<nEventMax; i++) {	
	   // get event (4-momenta)
	   eg.getEvent(Q2, W, Pefin, PK, PY, Ppfin, Ppim, Pgam);
	   // virtual photon flux
	   double gamma = getGamma(Ebeam,Q2,W);
	   // calculate angles of Kaon in CMS from  4-momenta of particles in LAB
	   lab2cms(Q2, Ebeam, Pefin, PK, PY, thetaK, phiK);
	   // fill histograms of dSigma/dCos(ThetaK)
           hTh->Fill(cos(thetaK),1./gamma);
	   
	   
	
	}
		
		
	// initilize model cross section 
	sigmaKY xsec("./data","KLambda");
	// Fill histogram for the model crtoss section 
        TH1F *hM = new TH1F("h","h",20,-1.,1.);
	hM->GetXaxis()->SetTitle("cos(#theta)");
	hM->GetYaxis()->SetTitle("Nevents");
	for(int i=1; i<=hM->GetXaxis()->GetNbins(); i++) {
	  double cost = hM->GetBinCenter(i);
	  double xsecVal = xsec.dsigma_dcos(Ebeam, Q2avr, Wavr, cost);
          hM->SetBinContent(i, xsecVal );
	}
		
		
        gStyle->SetOptStat(0);
	TCanvas *can = new TCanvas("c", "c" , 800, 600);
	can->Divide(1,1);
	can->cd(1);
	
	
	
	hTh->SetLineColor(1);
	hTh->SetLineStyle(1);
	hTh->SetLineWidth(2);

	hM->SetLineColor(2);
	hM->SetLineStyle(1);
	hM->SetLineWidth(3);	
	
	hTh->Scale( 100./hTh->Integral("width"));
	hM->Scale( 100./hM->Integral("width"));
	
	hTh->SetLineColor(1);
	hTh->SetLineStyle(1);
	hTh->SetLineWidth(2);
        hTh->DrawCopy("");
	hM->DrawCopy("LSAME");
	can->Print( "test.eps","eps" );
        can->WaitPrimitive();



	return 0;
}
