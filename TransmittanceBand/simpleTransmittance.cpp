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
#include "../MLibraries/MConfigGraphs.h"
#include "../MLibraries/MMultiFits.h"

// Open a csv file and read lines for wavelengths specified by wavelength.
TNtuple *openCSV(TNtuple wavelength, TString filename) {
	Int_t i;
	TNtuple *Tdata = new TNtuple("Tdata",filename.Data(),"w:T");
	Tdata->ReadFile(filename.Data());
	
	
}

int simpleTransmittance(void) {
}

int main(int argc, char *argv[]) {
        return simpleTransmittance(void);
}
