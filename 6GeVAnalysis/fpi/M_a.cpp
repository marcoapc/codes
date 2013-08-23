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

int M_a() {
	// RUN
	TString runN = "47031";	

	// Opening file
	TString filename = "hms" + runN + ".root";
	TFile *f = new TFile(filename);

	// Cheching file contents
	f->Print();
	//h9010->Print();

	// Getting tree from file
	TTree *t1 = (TTree*)f->Get("h9010");

	// Creating Q2 histogram
	Float_t q2;
	Int_t nentries = (Int_t)t1->GetEntries();
	t1->SetBranchAddress("Q2",&q2);
	TH1F *hQ2 = new TH1F("hQ2","Q2 - no cut (run " + runN + ")",100,0.5,2.0);
	for(Int_t i=0; i<nentries; i++) {
		t1->GetEntry(i);
		hQ2->Fill(q2);
	}
	/*TCanvas *c1 = new TCanvas("c1","c1",600,600);
	c1->SetGrid();
	c1->GetFrame()->SetBorderSize(10);
	hQ2->Draw();
	hQ2->GetXaxis()->SetTitle("X title");
	hQ2->GetYaxis()->SetTitle("Ytitle");
	pad2png(c1,"Q2_hist.png");
	*/

	// Creating W histogram
	Float_t w;
	t1->SetBranchAddress("W",&w);
	Int_t w_nbin=200;
	Double_t w_min=0.5, w_max=1.8;
	TH1F *hw = new TH1F("hw","W - no cut (run " + runN + ")",w_nbin,w_min,w_max);
	for(Int_t i=0; i<nentries; i++) {
		t1->GetEntry(i);
		hw->Fill(w);
	}
	//Plots
	TCanvas *c2 = new TCanvas("c2","c2",600,600);
	c2->SetGrid();
	c2->GetFrame()->SetBorderSize(10);
	hw->Draw();
	hw->GetXaxis()->SetTitle("W");
	hw->GetYaxis()->SetTitle("counts");
	hw->GetYaxis()->SetTitleOffset(1.3);
/*	//Fits
	TF1 *f1 = new TF1("f1","gaus",0.93,0.95);
	TF1 *f2 = new TF1("f2","gaus",1.17,1.27);
	Double_t par1[6];
	hw->Fit(f1,"R");
	hw->Fit(f2,"R+");
	f1->GetParameters(&par1[0]);
	f2->GetParameters(&par1[3]);
	cout << "Marco:\n  Dist peaks = " << (par1[4]-par1[1]) << endl;
*/
	//Save
	pad2png(c2,"W_hist.png");


	// Creating hcer_npe histogram
	Float_t hcer_npe;
	t1->SetBranchAddress("hcer_npe",&hcer_npe);
	TH1F *h_hcer_npe = new TH1F("h_hcer_npe","hcer_npe no cut (run " + runN + ")",100,0.,40.);
	for(Int_t i=0; i<nentries; i++) {
		t1->GetEntry(i);
		h_hcer_npe->Fill(hcer_npe);
	}
/*	TCanvas *c3 = new TCanvas("c3","c3",600,600);
	c3->SetGrid();
	c3->GetFrame()->SetBorderSize(10);
	h_hcer_npe->Draw();
	h_hcer_npe->GetXaxis()->SetTitle("hcer_npe");
	h_hcer_npe->GetYaxis()->SetTitle("counts");
	pad2png(c3,"hcer_npe.png");
*/

	// Creating selectable histogram
	TString v_plot = "hstof"; //
	TString descript = "";//"Tracked Total shower energy of chosen track";
	Float_t v1;
	t1->SetBranchAddress(v_plot,&v1);
	TH1F *h_v1 = new TH1F("h_v1",v_plot + " - no cut (run " + runN + ") - " + descript,100,0.,2500.);
	for(Int_t i=0; i<nentries; i++) {
		t1->GetEntry(i);
		h_v1->Fill(v1);
	}
	TCanvas *c4 = new TCanvas("c4","c4",600,600);
	c4->SetLogy();
	c4->SetGrid();
	c4->GetFrame()->SetBorderSize(10);
	h_v1->Draw();
	h_v1->GetXaxis()->SetTitle(v_plot);
	h_v1->GetYaxis()->SetTitle("counts");
	pad2png(c4,v_plot + ".png");

	// Trying some 3D histogram
	TH2F *hangle = new TH2F("hangle","Angular Dist",200,180.0,360.0,200,0.0,100.0);
	Float_t hstheta, hsphi;
	t1->SetBranchAddress("hstheta",&hstheta);
	t1->SetBranchAddress("hsphi",&hsphi);
	for(Int_t i=0; i<nentries; i++) {
		t1->GetEntry(i);
		hangle->Fill(hsphi*180.0/3.1415,hstheta*180.0/3.1415);
	}
/*	TCanvas *c5 = new TCanvas("c5","c5",600,600);
	c5->SetGrid();
	c5->GetFrame()->SetBorderSize(10);
	hangle->Draw("colz");
	hangle->GetXaxis()->SetTitle("hsphi / deg");
	hangle->GetYaxis()->SetTitle("hstheta / deg");
	pad2png(c5,"AngularDist.png");
*/

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

	// TOF vs calo
//	Float_t c7ycutmax = 48.0, c7ycutmin = 0.4, //48.0 0.4
//              c7xcutmax = 2.3, c7xcutmin = 0.35;
	TH2F *h_caloTOF = new TH2F("h_caloTOF","Calorimeter vs TOF signal",200,0.0,2.5,200,0.0,1500.0);
	Float_t hstof;
	t1->SetBranchAddress("hstof",&hstof);
	for(Int_t i=0; i<nentries; i++) {
		t1->GetEntry(i);
		h_caloTOF->Fill(hsshtrk,hstof);
	}
	TCanvas *c7 = new TCanvas("c7","c7",600,600);
	c7->SetGrid();
	c7->GetFrame()->SetBorderSize(10);
	h_caloTOF->Draw("colz");
	h_caloTOF->GetXaxis()->SetTitle("hsshtrk - calo");
	h_caloTOF->GetYaxis()->SetTitle("hstof - Time of Flight");
/*	TLine *c7l[4];
	c7l[0] = new TLine(c7xcutmin,c7ycutmin,c7xcutmax,c7ycutmin);
	c7l[1] = new TLine(c7xcutmax,c7ycutmin,c7xcutmax,c7ycutmax);
	c7l[2] = new TLine(c7xcutmax,c7ycutmax,c7xcutmin,c7ycutmax);
	c7l[3] = new TLine(c7xcutmin,c7ycutmax,c7xcutmin,c7ycutmin);
	for(Int_t i=0; i<4; i++){
		c7l[i]->SetLineColor(1);
		//c7l[i]->SetLineWidth(2);
		c7l[i]->Draw();
	}*/
	pad2png(c7,"TOF_calo.png");

	// hse vs W
//	Float_t c7ycutmax = 48.0, c7ycutmin = 0.4, //48.0 0.4
//              c7xcutmax = 2.3, c7xcutmin = 0.35;
	TH2F *h_hseW = new TH2F("h_hseW","W vs Lab Total Energy",w_nbin,w_min,w_max,200,0.0,5.0);
	Float_t hse;
	t1->SetBranchAddress("hse",&hse);
	for(Int_t i=0; i<nentries; i++) {
		t1->GetEntry(i);
		h_hseW->Fill(w,hse);
	}
	TCanvas *c8 = new TCanvas("c8","c8",600,600);
	c8->SetGrid();
	c8->GetFrame()->SetBorderSize(10);
	h_hseW->Draw("colz");
	h_hseW->GetXaxis()->SetTitle("W");
	h_hseW->GetYaxis()->SetTitle("hse - Lab Total Energy");
/*	TLine *c7l[4];
	c7l[0] = new TLine(c7xcutmin,c7ycutmin,c7xcutmax,c7ycutmin);
	c7l[1] = new TLine(c7xcutmax,c7ycutmin,c7xcutmax,c7ycutmax);
	c7l[2] = new TLine(c7xcutmax,c7ycutmax,c7xcutmin,c7ycutmax);
	c7l[3] = new TLine(c7xcutmin,c7ycutmax,c7xcutmin,c7ycutmin);
	for(Int_t i=0; i<4; i++){
		c7l[i]->SetLineColor(1);
		//c7l[i]->SetLineWidth(2);
		c7l[i]->Draw();
	}*/
	pad2png(c8,"W_hse.png");

	// Updating W histogram with gas cherenkov and calo cut
	c2->cd();
	//t1->Draw("W>>hw2",Form("hsshtrk >= %f && hsshtrk <= %f && hcer_npe >= %f && hcer_npe <= %f",c6xcutmin,c6xcutmax,c6ycutmin,c6ycutmax),"goff");
	//hw2->Draw("same");
	TH1F *hw2 = new TH1F("hw2","w cutted",w_nbin,w_min,w_max);
	for(Int_t i=0; i<nentries; i++) {
		t1->GetEntry(i);
		//if(i==0) cout << hsshtrk << " \t" << hcer_npe << endl;
		if(hsshtrk >= c6xcutmin && hsshtrk <= c6xcutmax && hcer_npe >= c6ycutmin && hcer_npe <= c6ycutmax)
			hw2->Fill(w);
	}
	hw2->SetFillColor(3);
	hw2->SetLineColor(3);
	hw2->SetFillStyle(3001);
	hw2->Draw("same");
//,hsshtrk >= c6xcutmin && hsshtrk <= c6xcutmax && hcer_npe >= c6ycutmin && hcer_npe <= c6ycutmax");
	pad2png(c2,"W_hist.png");

	//f->Close();
	return 0;
}

int main(int argc, char *argv[]) {
	TString inData = argv[1];

	return M_a02();
}
