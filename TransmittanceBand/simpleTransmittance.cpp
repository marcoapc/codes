#include<cstdio>
#include<string.h>
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
#include "../MLibraries/McommentNlines.h"

// Open a csv file and read lines for wavelengths specified by wavelength.
TNtuple *openCSV(TBranch *wavelength, TString filename) {
	Int_t i;
	TNtuple *Tdata = new TNtuple("Tdata",filename.Data(),"w:T");
	Tdata->ReadFile(commentNlines(filename));
	
	
}

int simpleTransmittance(void) {
/*	Float_t count;

	Float_t wl_range[2] = {200.0, 900.0};
	Float_t wl_steps = 10.0;

	TTree *tr = new TTree("tr","Transmittance");
	Float_t wvl;
	TBranch *wl = tr->Branch("wl",&wvl,"wl/F");
	for(wvl=wl_range[0]; wvl<=wl_range[1]; wvl+=wl_steps) wl->Fill();
	tr->SetEntries(wl->GetEntries());


	TString filename = "/media/marco/Windows/Users/Marco/Documents/Marco_CUA/KaonDetector/Aerogel_Characteristics/TransmittanceData/transmittancesp30blast1samples/teste.csv";*/

	Float_t linewidth = 2.0;

	TCanvas *cf = new TCanvas("cf","Transmittance Blast",1200,600);
	MConfigCanvas(cf);
	
	TNtuple *Tdata1 = new TNtuple("Tdata1","Data1","w:T");
	Tdata1->ReadFile(commentNlines("/home/marco/Documents/Marco_CUA/KaonDetector/Aerogel_Characteristics/TransmittanceData/transmittancesp30blast1samples/Sample297.Sample.Raw.csv"));
	Tdata1->Draw("w:T","","goff");
	TGraph *gr1 = new TGraph(Tdata1->GetSelectedRows(),Tdata1->GetV1(),Tdata1->GetV2());
	gr1->SetLineColor(2);
	gr1->SetLineWidth(linewidth);
	gr1->SetTitle("Transmittance SP30 from 1st Blast detector");
	gr1->GetXaxis()->SetTitle("Wavelength / nm");
	gr1->GetYaxis()->SetTitle("Transmittance / %");
	gr1->Draw("al");

	TNtuple *Tdata2 = new TNtuple("Tdata2","Data2","w:T");
	Tdata2->ReadFile(commentNlines("Sample299.Sample.Raw.csv"));
	Tdata2->Draw("w:T","","goff");
	TGraph *gr2 = new TGraph(Tdata2->GetSelectedRows(),Tdata2->GetV1(),Tdata2->GetV2());
	gr2->SetLineColor(3);
	gr2->SetLineWidth(linewidth);
	gr2->Draw("same");

	TNtuple *Tdata3 = new TNtuple("Tdata3","Data3","w:T");
	Tdata3->ReadFile(commentNlines("Sample301.Sample.Raw.csv"));
	Tdata3->Draw("w:T","","goff");
	TGraph *gr3 = new TGraph(Tdata3->GetSelectedRows(),Tdata3->GetV1(),Tdata3->GetV2());
	gr3->SetLineColor(4);
	gr3->SetLineWidth(linewidth);
	gr3->Draw("same");

	TNtuple *Tdata4 = new TNtuple("Tdata4","Data4","w:T");
	Tdata4->ReadFile(commentNlines("Sample303.Sample.Raw.csv"));
	Tdata4->Draw("w:T","","goff");
	TGraph *gr4 = new TGraph(Tdata4->GetSelectedRows(),Tdata4->GetV1(),Tdata4->GetV2());
	gr4->SetLineColor(5);
	gr4->SetLineWidth(linewidth);
	gr4->Draw("same");

	TNtuple *Tdata5 = new TNtuple("Tdata5","Data5","w:T");
	Tdata5->ReadFile(commentNlines("Sample305.Sample.Raw.csv"));
	Tdata5->Draw("w:T","","goff");
	TGraph *gr5 = new TGraph(Tdata5->GetSelectedRows(),Tdata5->GetV1(),Tdata5->GetV2());
	gr5->SetLineColor(6);
	gr5->SetLineWidth(linewidth);
	gr5->Draw("same");

	TNtuple *Tdata6 = new TNtuple("Tdata6","Data6","w:T");
	Tdata6->ReadFile(commentNlines("Sample302.Sample.Raw.csv"));
	Tdata6->Draw("w:T","","goff");
	TGraph *gr6 = new TGraph(Tdata6->GetSelectedRows(),Tdata6->GetV1(),Tdata6->GetV2());
	gr6->SetLineColor(7);
	gr6->SetLineWidth(linewidth);
	gr6->Draw("same");

	TLegend *leg = new TLegend(0.14,0.68,0.4,0.87);
	leg->SetTextFont(50);
	leg->SetTextSize(0.025);
	leg->AddEntry(gr1,"Top layer, touching detector's wall","l");
	leg->AddEntry(gr2,"Top layer, touching detector's wall","l");
	leg->AddEntry(gr3,"Top layer, far from the wall","l");
	leg->AddEntry(gr4,"2nd layer, far from the wall","l");
	leg->AddEntry(gr6,"2nd layer, far from the wall","l");
	leg->AddEntry(gr5,"Top layer, center","l");
	leg->Draw();

	pad2png(cf,"TransmittanceBlast.png");

	return 0;
}

int main(int argc, char *argv[]) {
	int u;
	//commentNlines("teste.txt");
	//commentNlines("/media/marco/Windows/Users/Marco/Documents/Marco_CUA/KaonDetector/Aerogel_Characteristics/TransmittanceData/transmittancesp30blast1samples/teste.csv");
	//return 0;
        return simpleTransmittance();
}
