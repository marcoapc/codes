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
#include<TTree.h>
#include<TBranch.h>
#include<TLeafF.h>
#include<TBranchRef.h>

// MARCO FUNCTIONS
#include "MConfigGraphs.h"
#include "./MMultiFits.h"

int M_plotvars(TString runN, TString var1, TString var2="") {

	// Opening file
	TString filename = "hms" + runN + ".root";
	TFile *f = new TFile(filename);

	// Cheching file contents
	f->Print();
	//h9010->Print();

	// Getting tree from file
	TTree *t1 = (TTree*)f->Get("h9010");
	Int_t nentries = (Int_t)t1->GetEntries();

	if(var2.Length() == 0) {
		// Creating selectable histogram
		TString descript = "";//"Tracked Total shower energy of chosen track";
		Float_t v1;
		t1->SetBranchAddress(var1,&v1);
		cout << "Marco: " << v1.GetMaximum() << endl;
		TH1F *h_v1 = new TH1F("h_v1",var1 + " - no cut (run " + runN + ") - " + descript,200,0.,2500.);
		for(Int_t i=0; i<nentries; i++) {
			t1->GetEntry(i);
			h_v1->Fill(v1);
		}
		TCanvas *c4 = new TCanvas("c4","c4",600,600);
		c4->SetLogy();
		c4->SetGrid();
		c4->GetFrame()->SetBorderSize(10);
		h_v1->Draw();
		h_v1->GetXaxis()->SetTitle(var1);
		h_v1->GetYaxis()->SetTitle("counts");
		pad2png(c4,var1 + ".png");
	}
/*
	// Calo vs Cherenkov
	Float_t c6ycutmax = 48.0, c6ycutmin = 0.4, //48.0 0.4
                c6xcutmax = 2.3, c6xcutmin = 0.35;
	TH2F *h_caloCher = new TH2F("h_caloCher","Calorimeter vs Gas Cherenkov signal",200,0.0,2.5,200,0.0,50.0);
	Float_t hsshtrk;
	t1->SetBranchAddress("hsshtrk",&hsshtrk);
	for(Int_t i=0; i<nentries; i++) {
		t1->GetEntry(i);
		h_caloCher->Fill(hsshtrk,hcer_npe);
	}
	TCanvas *c6 = new TCanvas("c6","c6",600,600);
	c6->SetGrid();
	c6->GetFrame()->SetBorderSize(10);
	h_caloCher->Draw("colz");
	h_caloCher->GetXaxis()->SetTitle("hsshtrk - calo");
	h_caloCher->GetYaxis()->SetTitle("hcer_npe - Gas Cherenkov");
	TLine *c6l[4];
	c6l[0] = new TLine(c6xcutmin,c6ycutmin,c6xcutmax,c6ycutmin);
	c6l[1] = new TLine(c6xcutmax,c6ycutmin,c6xcutmax,c6ycutmax);
	c6l[2] = new TLine(c6xcutmax,c6ycutmax,c6xcutmin,c6ycutmax);
	c6l[3] = new TLine(c6xcutmin,c6ycutmax,c6xcutmin,c6ycutmin);
	for(Int_t i=0; i<4; i++){
		c6l[i]->SetLineColor(1);
		//c6l[i]->SetLineWidth(2);
		c6l[i]->Draw();
	}
	pad2png(c6,"Cher_calo.png");
*/

	//f->Close();
	return 0;
}

int main(int argc, char *argv[]) {
	TString inData = argv[1];

	return M_plotvars();
}
