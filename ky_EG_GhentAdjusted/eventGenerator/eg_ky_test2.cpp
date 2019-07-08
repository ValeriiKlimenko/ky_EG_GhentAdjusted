
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
	
	double Ebeam = 11;
	double Q2min = 1.0;
	double Q2max = 11.999;
	double Wmin = 2.5;
	double Wmax = 5.0;
        int nEventMax = 100000;

	double Q2avr = 0.5*(Q2max+Q2min);
	double Wavr  = 0.5*(Wmax+Wmin);
	double epsAvr = getEpsilon(Ebeam, Q2avr, Wavr);
	double Q2,W,cosThetaK,thetaK,phiK,d5sigma;
         
	//string chname = "KLambda";
	string chname = "KSigma";

	// initilize event generator
        evGenerator eg("./data", chname.c_str(), Ebeam, Q2min, Q2max, Wmin, Wmax);
	
	//histogram for events
        string sn = "#gamma_{v} p #rightarrow K^{+}"; 
	if(chname=="KLambda") sn += " #Lambda"; else sn += " #Sigma^{0}";
	TH1F *hQ2 = new TH1F("Q2",sn.c_str(),50,  0.05,12);
	hQ2->GetXaxis()->SetTitle("Q^{2}, GeV^{2}");
	hQ2->GetYaxis()->SetTitle("N_{ev}");
        TH1F *hW = new TH1F("W",sn.c_str(),50, 1.5,4.5);
	hW->GetXaxis()->SetTitle("W, GeV");
	hW->GetYaxis()->SetTitle("N_{ev}");
        TH2F *hWQ2 = new TH2F("WQ2",sn.c_str(), 50,1.5,4.5,  50,0.05,12);
	hWQ2->GetXaxis()->SetTitle("W, GeV");
	hWQ2->GetYaxis()->SetTitle("Q^{2}, GeV^{2}");

        // 4-momenta of final electron, Kaon and Lambda LAB frame
	TLorentzVector Pefin, PK, PY, Ppfin, Ppim, Pgam;
	// Loop through events
	int nEvents = 0;
        for (int i=0; i<nEventMax; i++) {
	    nEvents++;
	   if((nEvents/10000)*10000 == nEvents ) {
	     cout << " Event # " << nEvents << endl;
	   }	
	   eg.getEvent(Q2, W, Pefin, PK, PY, Ppfin, Ppim, Pgam);
	   double gamma = getGamma(Ebeam,Q2,W);
	   lab2cms(Q2, Ebeam, Pefin, PK, PY, thetaK, phiK);
           hQ2->Fill(Q2,1./gamma);
           hW->Fill(W,1./gamma);
	   hWQ2->Fill(W,Q2,1./gamma);	   
	}
		
	
	/*	
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
	*/
	
	
	
	
		
        gStyle->SetOptStat(0);
	TCanvas *can = new TCanvas("c", "c" , 1000, 500);
	can->Divide(1,1);
	
	hQ2->SetLineColor(1);
	hQ2->SetLineStyle(1);
	hQ2->SetLineWidth(2);
	hQ2->SetMinimum(0);


	hW->SetLineColor(1);
	hW->SetLineStyle(1);
	hW->SetLineWidth(2);	
	hW->SetMinimum(0);
	
	
	
	hWQ2->SetLabelSize(0.05,"X");
	hWQ2->SetTitleSize(0.05,"X");
	hWQ2->SetLabelSize(0.05,"Y");
	hWQ2->SetTitleSize(0.05,"Y");
	hW->SetLabelSize(0.05,"X");
	hW->SetTitleSize(0.05,"X");
	hW->SetLabelSize(0.05,"Y");
	hW->SetTitleSize(0.05,"Y");
	hQ2->SetLabelSize(0.05,"X");
	hQ2->SetTitleSize(0.05,"X");
	hQ2->SetLabelSize(0.05,"Y");
	hQ2->SetTitleSize(0.05,"Y");
	
	
	can->cd(1); hWQ2->DrawCopy("COLZ");
	can->Print( "q2w.eps","eps" );
        can->WaitPrimitive();

	can->cd(1); hQ2->DrawCopy("");
	can->Print( "q2.eps","eps" );
        can->WaitPrimitive();

	can->cd(1);
	hW->DrawCopy("");
	can->Print( "w.eps","eps" );
        can->WaitPrimitive();



	return 0;
}
