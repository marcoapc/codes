#include "../MLibraries/MallIncludes.h"
#include "../MLibraries/MConfigGraphs.h"

int ReadCalib(TString inFile="",TString outFile="") {

	Int_t i;

	TNtuple *fin = new TNtuple("fin","Input Data","x");
	fin->ReadFile(inFile);
	Double_t npts = fin->GetEntries();
	cout << "<Marco> Number of entries: " << npts << endl;

	Double_t range[2];
	range[0] = fin->GetMinimum("x");
	range[1] = fin->GetMaximum("x");

	TH1F *h = new TH1F("h","Calibration",((int)(range[1]-range[0]+1)),((double)(range[0])),((double)(range[1])));

	for(i=0; i<npts; i++) h->Fill(fin->GetEntry(i));

	TCanvas *c1 = new TCanvas("c1","Data",600,600);
	MConfigHist(h,inFile.Data(),"Channel number","counts");
	h->Draw();
}
