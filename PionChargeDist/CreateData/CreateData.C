#include<cstdio>
#include<TROOT.h>
#include<TRint.h>
#include<stdlib.h>
#include<TCanvas.h>
#include<TStyle.h>
#include<TNtuple.h>
#include<TGraph.h>
#include<TGraphErrors.h>
#include<TLegend.h>
#include<TMultiGraph.h>
#include<TFrame.h>
#include<TPaveText.h>
#include<TFile.h>
#include<TH1.h>
#include<TF1.h>
#include<TMath.h>
#include<Riostream.h>
//#include<iostream>
#include<TImage.h>
#include<TRandom.h>

// MARCO FUNCTIONS
#include "../../MLibraries/MConfigGraphs.h"

// Return Form Factor relative to a monopole:
// F(Q2) = 1.0 / (1.0 + Q2/G2)
Double_t monopole(Float_t Q2, Float_t G2 = 0.507205) {
	return (1.0 / (1.0+(Q2/G2)));
}

// Return Form Factor relative to a dipole:
// F(Q2) = 1.0 / (1.0 + Q2/G2)^2
Double_t dipole(Float_t Q2, Float_t G2 = 0.71) {
	return (1.0 / pow(1.0+(Q2/G2),2));
}

int CreateData(Float_t Q2min=0.0, Float_t Q2max=1000.0, Int_t npts=5000, TString outfile="out.txt") {
	Int_t i;
	Float_t Q2, FQ2;
	ofstream fout;
	fout.open(outfile);
	if(!fout.is_open()) {
		cout << endl << "<Marco> Error opening " << outfile << endl
		     << "        Aborting program execution." << endl << endl;
		return 1;
	}
	for(i=0;i<npts;i++) {
		Q2 = Q2min + i*(Q2max-Q2min)/((double)(npts-1));
		//FQ2 = monopole(Q2);
		FQ2 = dipole(Q2);
		// Q2 F(Q2) sF
		fout << Form("%.8f\t%.8f\t%.8f\n", Q2, FQ2, 0.05*FQ2);
	}
	fout.close();
	return 0;
}

void Instructions(void) {
	cout << endl << endl << endl
	     << "Usage:" << endl
	     << "  ./CreateData <Q2min> <Q2max> <npts> <outfile>" << endl << endl
	     << "where:" << endl
	     << "  - <Q2min> is the minimum value of Q2 to be created" << endl
	     << "  - <Q2max> is the maximum value of Q2 to be created" << endl
	     << "  - <npts>  is the number of points in Q2 to be created" << endl
	     << "  - <outfile> is the filename of the output" << endl << endl
	     << "by Marco Antonio Pannunzio Carmignotto, 18 of August, 2013" << endl << endl; 

	return;
}

int main(int argc, char *argv[]) {
	if(argc==1) Instructions();
	return CreateData();
}
