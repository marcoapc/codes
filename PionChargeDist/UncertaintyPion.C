/*#include<cstdio>
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
#include<TLatex.h>*/

// MARCO FUNCTIONS
#include "../MLibraries/MallIncludes.h"
#include "../MLibraries/MConfigGraphs.h"
#include "../MLibraries/MMultiFits.h"
#include "CreateErrorBand.C"

//TF1 *MakeGe(TNtuple *pdata, Double_t *range);
//TF1 *MakeGm(TNtuple *pdata, Double_t *range);
void PlotFF(TNtuple *FF, Double_t R1, Double_t *range);
Double_t CalcR1(Double_t *range, TNtuple *FF, Double_t MaxRangeForFit, Int_t Mgraph=0);
TF1 *FF1fit(TNtuple *FF, Double_t *range, Int_t q2f=1);
//Double_t CalcR2(TNtuple *pdata, Double_t *range, Double_t *rangeSmallQ2, TNtuple *FF, Double_t FracRangeForFit);

// Compute the n^th zero (n = 1, 2, ...) of the Bessel Function J0
Double_t ZeroBesselJ0(Int_t n, Double_t min=0.0) {
	Double_t max, val;
	if(min==0.0) max = 3.93;
	else max = min + TMath::Pi();

	TF1 *J0 = new TF1("J0","TMath::BesselJ0(x)",min,max);
	
	if(n==1) val = J0->GetX(0.0,min,max);
	else {
		Double_t tmp;
		tmp = J0->GetX(0.0,min,max);
		val = ZeroBesselJ0(n-1,max);
	}
	delete J0;
	return val;
}

int UncertaintyPion(TString inData) {
	// Variables
	char line[500];
	Double_t mu = 2.792782;
	Int_t i, j, k;
	TString outFile = "pionUncert"; // "n1 + .root" will be added later

	// Important information/config
	
/*	// Opening data input file
	TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
	dir.ReplaceAll("ChargeDist.C","");
	dir.ReplaceAll("/./","/");
	cout << "<Marco> current folder: " << dir.Data() << endl;
	ifstream in;
	in.open(Form("%s%s",dir.Data(),inData.Data()));
	if(in.fail()) {
		cout << "<Marco> Error opening file." << endl;
		cout << "\t" << dir.Data() << inData.Data() << endl;
		return 1;
	}*/

	// Reading data
	TNtuple *FF = new TNtuple("FF","PionData","Q2:F1:sF1");
	FF->ReadFile(inData);
	cout << "<Marco> Number of data points: " << FF->GetEntries() << endl;
	Double_t range[2];
	range[0] = FF->GetMinimum("Q2");
	range[1] = FF->GetMaximum("Q2");
	// pdata->Scan(); // to show data

// APAGAR!!!! - Usar somente para mostrar maior range (CAUTION: Q2 range used to evaluate charge distribution will be changed too)
range[1] = 30.0;

	//TF1 *ge_paper = MakeGe(pdata, range);
	//TF1 *gm_paper = MakeGm(pdata, range);

	/*********************
	 *   FORM FACTORS    *
	 * F1(Q2) AND F2(Q2) *
	 *********************/
/*
	// Creating NTuples with F1(Q2) and F2(Q2)
	const unsigned int Npts = 500;
	TNtuple *FF = new TNtuple("FF","Form Factor","Q2:Ge:Gm:F1:F2");
	Double_t Q2;
	Float_t a[4];
	Float_t tau;
	range[0]=0.05; range[1]=30.0;
	for(Int_t i=0; i<Npts; i++) {
		Q2 = range[0] + (((double)i)/(((double)Npts)-1.0))*(range[1]-range[0]); // Q2
		tau = Q2/(4.0*pow(0.93827,2.0)); // Proton mass here!
		a[0] = ge_paper->Eval(Q2);       // Ge
		a[1] = mu*gm_paper->Eval(Q2);    // Gm
		//a[0] = ge_pol->EvalPar(Q2);       // Ge
		//a[1] = gm_pol->EvalPar(Q2);       // Gm
		a[2] = (a[0] + tau*a[1])/(1.0+tau); // F1
		a[3] = (a[1] - a[0])/(1.0+tau);     // F2
		FF->Fill(Q2,a[0],a[1],a[2],a[3]);
	}
	cout << "<Marco> Last: " << Q2 << a[3] << endl;
	//FF->Scan();
*/	
	// Fitting F x Q2 for Q2->0 (analysis of R)
	///// CONFIGURE HERE! //////
	Double_t MaxRangeForFit = 0.1; // 100*x% of data range - can change as one wants! - using 0.1

	// Calculating R1 and R2 in fm, from Q2->0 data
	Double_t R1 = CalcR1(range, FF, MaxRangeForFit);

	// Plotting F1
	PlotFF(FF,R1,range);

	/*************************************
	 *************************************
	 ** Calculating charge distribution **
	 *************************************
	 *************************************/

	/******************
	 * ro_ch, from F1 *
	 ******************/
	// upper limit of summation
	Int_t n1;
	Double_t Ecut_F1 = range[1];

	R1 = R1*5.068;	// CONVERT R1 TO GeV^-1!!!!!!
	for(n1=1; (ZeroBesselJ0(n1)/R1) <= TMath::Sqrt(Ecut_F1); n1++);
	n1--; // The last term didn't match the condition, so subtr. one

	// Output file
	TFile *f = new TFile(outFile+Form("%d.root",n1),"RECREATE");

//APAGAR!!!!!! Usar somente se quiser especificar o numero de termos (comentar as linhas <>, do ciclo
//n1=3;


	cout << "<Marco> Considering " << n1 << " terms for ro_ch (F1)." << endl;

	// Doesn't make any difference, just create plot (don't use Ffitted)
	TF1 *Ffitted = FF1fit(FF, range, 1);

	/////////////////////////////////
	// Making several fits to data //
	/////////////////////////////////
	const int nfits = 20; // <- NUMBER OF FITS! (200)
	TRandom *r0 = new TRandom();
	TCanvas *cf = new TCanvas("cf","MultiFits",600,600);
	MConfigCanvas(cf);
	FF->Draw("Q2:F1:sF1","","goff");
	TGraphErrors *origData = new TGraphErrors(FF->GetSelectedRows(),FF->GetV1(),FF->GetV2(),0,FF->GetV3());
	MConfigPoints(origData,2);
	MConfigAxis(origData,"","Q^{2} (GeV^{2})","F_{#pi}(Q^{2})");
	origData->Draw("ap");
	
	//TGraphErrors *newData;
	TGraph *newData;
	TF1 *fits[nfits];
	const int nparFit=3;
	Double_t par0[nparFit] = {0.5, 1.2, 1.1};
	//Double_t par0[nparFit] = {0.25, 0.2, 2.75};
	TH1F *parFit[nparFit];
	parFit[0] = new TH1F("par0","MultiFit par[0]",50,0.0,1.0);
	parFit[1] = new TH1F("par1","MultiFit par[1]",50,0.0,2.0);
	parFit[2] = new TH1F("par2","MultiFit par[2]",50,0.0,2.0);
	for(i=0; i<nfits; i++) {
		newData=MakeNewPoints(FF,r0);
		fits[i]=ExecuteFit(newData,cf,par0,1); //last 1 means quiet fit, don't show fitted paramenters. If 0, will show fitted parmts.
		for(j=0; j<nparFit; j++)
			parFit[j]->Fill(fits[i]->GetParameter(j));
		delete newData;
	}
	origData->Draw("psame");
	pad2png(cf,"MultiFits.png");
	cf->Print("MultiFits.eps");

	// New plot of fittings, with Q2xF
	TCanvas *cf2 = new TCanvas("cf2","MultiFits2",600,600);
	MConfigCanvas(cf2);
	FF->Draw("Q2:(Q2*F1):(Q2*sF1)","","goff");
	TGraphErrors *origData2 = new TGraphErrors(FF->GetSelectedRows(),FF->GetV1(),FF->GetV2(),0,FF->GetV3());
	MConfigPoints(origData2,2);
	MConfigAxis(origData2,"","Q^{2} (GeV^{2})","Q^{2}.F_{#pi}(Q^{2})");
	//MConfigAxis(origData2,"Multi Fits of F_{#pi}","Q^{2} (GeV^{2})","Q^{2}.F_{#pi}(Q^{2})");
	origData2->GetXaxis()->SetLimits(0.0,1.1*range[1]);
	origData2->GetYaxis()->SetRangeUser(0.0,1.1);
	origData2->Draw("ap");
	TF1 *fitsQ2[nfits];

	for(i=0; i<nfits; i++) {
		fitsQ2[i] = new TF1("fitsQ2i",MfitFQ2,range[0],range[1],3); //care with npar=3
		fitsQ2[i]->SetParameters(fits[i]->GetParameters());
		MConfigLines(fitsQ2[i],13,1.0);
		fitsQ2[i]->Draw("same");
	}
	origData2->Draw("psame");
	// Ideal monopole
	TF1 *monopole = new TF1("monopole","x*1.0/(1.0+[0]*x)",0.0,range[1]);
	monopole->SetLineColor(3);
	monopole->SetLineStyle(1);
	monopole->SetParameter(0,pow(R1/5.068/0.197327,2.0)/6.0);
	monopole->Draw("same");
	// Ideal dipole
	TF1 *dipole = new TF1("dipole","x*1.0/pow(1.0+[0]*x,2.0)",0.0,range[1]);
	dipole->SetLineColor(6);
	dipole->SetLineStyle(1);
	dipole->SetParameter(0,pow(R1/5.068/0.197327,2.0)/12.0);
	dipole->Draw("same");

	//Drawing JLAB 12 data, if wanted
	Int_t drawJlab12 = 1;
	TGraphErrors *jlab12;
	if(drawJlab12) {
		Double_t xjlab12[7] = {0.3,1.6,2.45,3.5,4.5,5.25,6.0};
		Double_t yjlab12[7] = {0.5,0.5,0.5,0.5,0.5,0.5,0.5};
		Double_t ejlab12[7] = {0.026,0.018,0.015,0.020,0.0215,0.025,0.033};
		jlab12 = new TGraphErrors(7,xjlab12,yjlab12,0,ejlab12);
		jlab12->SetMarkerColor(4);
		jlab12->SetMarkerStyle(33);
		jlab12->SetMarkerSize(1);
		jlab12->Draw("psame");
	}

	//Drawing EIC data, if wanted
	Int_t drawEIC = 1;
	TGraphErrors *eic5, *eic10, *eic15;
	if(drawEIC) {
		Double_t xeic5[6] = {5.0,6.0,7.5,10.0,12.5,15.0};
		Double_t yeic5[6] = {0.6,0.6,0.6,0.6,0.6,0.6};
		Double_t eeic5[6] = {0.05,0.05,0.05,0.04,0.04,0.04};
		eic5 = new TGraphErrors(6,xeic5,yeic5,0,eeic5);
		eic5->SetMarkerColor(6);
		eic5->SetMarkerStyle(33);
		eic5->SetMarkerSize(1);
		eic5->Draw("psame");
	}
	if(drawEIC) {
		Double_t xeic10[5] = {10.0,12.5,15.0,17.5,20.0};
		Double_t yeic10[5] = {0.7,0.7,0.7,0.7,0.7};
		Double_t eeic10[5] = {0.04,0.03,0.04,0.04,0.04};
		eic10 = new TGraphErrors(5,xeic10,yeic10,0,eeic10);
		eic10->SetMarkerColor(8);
		eic10->SetMarkerStyle(33);
		eic10->SetMarkerSize(1);
		eic10->Draw("psame");
	}
	if(drawEIC) {
		Double_t xeic15[4] = {15.0,17.5,20.0,25.0};
		Double_t yeic15[4] = {0.8,0.8,0.8,0.8};
		Double_t eeic15[4] = {0.03,0.03,0.03,0.03};
		eic15 = new TGraphErrors(4,xeic15,yeic15,0,eeic15);
		eic15->SetMarkerColor(28);
		eic15->SetMarkerStyle(33);
		eic15->SetMarkerSize(1);
		eic15->Draw("psame");
	}


	//Drawing bands, if wanted
	Int_t drawBands = 0; //1-YES, 0-NO
	TGraph *grshade1 = new TGraph(5);
	TGraph *grshadec1 = new TGraph(5);
	if(drawBands) {
		Double_t shade1xlim[2] = {12.0, 15.0};
		Double_t shade1ylim[2] = {0.2, 0.4};
		grshade1->SetPoint(0,shade1xlim[0],shade1ylim[0]);
		grshade1->SetPoint(1,shade1xlim[1],shade1ylim[0]);
		grshade1->SetPoint(2,shade1xlim[1],shade1ylim[1]);
		grshade1->SetPoint(3,shade1xlim[0],shade1ylim[1]);
		grshade1->SetPoint(4,shade1xlim[0],shade1ylim[0]);
		grshadec1->SetPoint(0,shade1xlim[0],shade1ylim[0]);
		grshadec1->SetPoint(1,shade1xlim[1],shade1ylim[0]);
		grshadec1->SetPoint(2,shade1xlim[1],shade1ylim[1]);
		grshadec1->SetPoint(3,shade1xlim[0],shade1ylim[1]);
		grshadec1->SetPoint(4,shade1xlim[0],shade1ylim[0]);
		grshade1->SetFillStyle(3144);
		grshade1->SetFillColor(2);
		//grshadec1->SetLineColor(1);
		//grshadec1->SetLineWidth(1);
		grshade1->Draw("fsame");
		grshadec1->Draw("lsame");
	}

	// Legend 
	TLegend *leg = new TLegend(0.43,0.72,0.95,0.98);
	//leg->SetTextFont(50);
	leg->SetTextSize(0.025);
	leg->AddEntry(origData2,"Data points","p");
	leg->AddEntry(fitsQ2[0],"Evaluated fits","l");
	leg->AddEntry(monopole,Form("Monopole (R = %.3f fm)",sqrt(6.0*pow(0.197327,2.0)*monopole->GetParameter(0))));
	leg->AddEntry(dipole,Form("Dipole (R = %.3f fm)",sqrt(12.0*pow(0.197327,2.0)*dipole->GetParameter(0))));
	if(drawBands) {
		leg->AddEntry(grshade1,"Future experiment - JLab @ 12GeV","f");
	}
	if(drawJlab12) leg->AddEntry(jlab12,"Future experiment - JLab @ 12GeV","p");
	if(drawEIC) {
		leg->AddEntry(eic5,"Future experiment - EIC E_{p}=5GeV","p");
		leg->AddEntry(eic10,"Future experiment - EIC E_{p}=10GeV","p");
		leg->AddEntry(eic15,"Future experiment - EIC E_{p}=15GeV","p");
	}
	leg->Draw();
	// Saving figure
	pad2png(cf2,"MultiFitsQ2.png");	
	cf2->Print("MultiFitsQ2.eps");

	// Histograms of the fitted parameters
	TCanvas *pp[nparFit];
	for(j=0; j<nparFit; j++) {
		pp[j] = new TCanvas(Form("pp[%d]",j),Form("par[%d] from MultiFit",j),600,600);
		MConfigCanvas(pp[j]);
		MConfigHist(parFit[j],Form("par[%d]",j),"Parameter value","counts");
		parFit[j]->Draw();
		pad2png(pp[j],Form("MultiFit_par[%d].png",j));
		pp[j]->Print(Form("MultiFit_par%d.eps",j));
	}

	////////////////////////////////////////
	// EVALUATION OF CHARGE DISTRIBUTION! //
	////////////////////////////////////////

// IF ONE WANTS TO RUN A SERIES OF CALCULATION, UNCOMMENT THE <> LINES

	//Writting exact solutions (for ideal complete Q2 data), from analytical solution
	//MONOPOLE
	TF1 *solMonopole = new TF1("solMonopole","(1.0/(2.0*TMath::Pi()))*(6.0/pow([0],2))*TMath::BesselK0(sqrt(6.0)*x/[0])",0.0,2.0); //[0] is R1 in units of b
	Double_t par_solMonopole[1];
	par_solMonopole[0] = R1/5.068; //in fm
	solMonopole->SetParameters(par_solMonopole);
	MConfigLines(solMonopole,2);
	//DIPOLE
	TF1 *solDipole = new TF1("solDipole","(6.0*sqrt(3.0)*x/(TMath::Pi()*pow([0],3)))*TMath::BesselK1(sqrt(12.0)*x/[0])",0.0,2.0); //[0] is R1 in units of b
	Double_t par_solDipole[1];
	par_solDipole[0] = R1/5.068; //in fm
	solDipole->SetParameters(par_solDipole);
	MConfigLines(solDipole,3);
	//MONOPOLE + DIPOLE complete
	TF1 *solMonDi = new TF1("solMonDi","[0]*(1.0/(2.0*TMath::Pi()*[1]))*TMath::BesselK0(x/sqrt([1])) + (1.0-[0])*(x/(4.0*TMath::Pi()*pow([2],3.0/2.0)))*TMath::BesselK1(x/sqrt([2]))",0.0,2.0);
	Double_t par_solMonDi[3];
	par_solMonDi[0] = parFit[0]->GetMean(1);
	par_solMonDi[1] = parFit[1]->GetMean(1)*0.197327*0.197327;
	par_solMonDi[2] = parFit[2]->GetMean(1)*0.197327*0.197327;
	for(i=0;i<3;i++) 
		cout << "parFit[" << i << "]_mean = " << parFit[i]->GetMean(1) << endl;
	solMonDi->SetParameters(par_solMonDi);
	MConfigLines(solMonDi,4);

	//Creating date
	TDatime da;
	da.GetDate();
	TString da2 = da.AsSQLString();

	// Creating the NTuple and evaluating the sum (eq 5)
	cout << "<Marco> Evaluating charge distribution..." << endl;
	Float_t rangeB[2] = {0.0, 0.8};
	Int_t nb = 100; // number of points (in b) to be calculated
	Float_t ba, Xn, tauN, F1, F2;
	Float_t Q2n, sum, bb;
	Double_t xlim[2] = {0.0, 1.0};
	Double_t ylim[2] = {0.0, 5.0};
//<>
//for(n1=1;n1<=30;n1++) {
//<>
	TTree *tr = new TTree("tr","Tree_ro1");
	//TNtuple *ro1 = new TNtuple("ro1","ro_ch(b)","b:ro");
	TBranch *b = tr->Branch("b",&ba,"b/F");
	// Filling b
	for(ba=rangeB[0]; ba<=rangeB[1]; ba += (rangeB[1]-rangeB[0])/((double)(nb-1))) { b->Fill(); }
	tr->SetEntries(b->GetEntries());
	// Running the charge calculation for each fit of F(Q^2)
	TBranch *ros[nfits];
	for(j=0; j<nfits; j++) {
		cout << "Calculating charge dist. for fit " << j+1 << "/" << nfits << endl;
		ros[j] = tr->Branch(Form("ro%d",j),&sum,Form("ro%d/F",j));
		// Running the calculation for each b
		if(j==0) { cout << endl << "n1=" << n1 << endl << "    invoqued Q2 = "; }
		for(bb=rangeB[0]; bb<=rangeB[1]; bb += (rangeB[1]-rangeB[0])/((double)(nb-1))) {
			sum=0.0;
			if(bb<=(R1/5.068)) {
				for(i=1; i<=n1; i++) {
					// R1 in GeV^-1
					Xn = ZeroBesselJ0(i);
					Q2n = pow(Xn/R1,2.0);
					if(j==0 && bb==rangeB[0]) { cout << Q2n << " "; };
					//tauN = Q2n/(4.0*pow(0.93827,2.0));
					F1=fits[j]->Eval(Q2n);
					sum += pow(TMath::BesselJ1(Xn),-2.0)*F1*TMath::BesselJ0(Xn*bb/(R1/5.068));
				}
				sum /= (TMath::Pi()*pow(R1/5.068,2.0));
				//cout << "<Marco> bb = " << bb << "\t ro = " << sum << endl;
			}
			else sum=0.0;
			ros[j]->Fill();
		}
		if(j==0) { cout << endl; };
	}
	//If want to save data, uncomment here
	f->Write();
	
	// Plotting
	cout << "Plotting charge distribution evaluation..." << endl;
	TCanvas *c5 = new TCanvas("c5","ChargeDistributionF1",600,600);
	c5->SetGrid();
	c5->GetFrame()->SetBorderSize(10);
	TGraph *gr5[nfits];
	for(i=0; i<nfits; i++) {
		tr->Draw(Form("b:ro%d",i),"","goff"); // Apply data selection here
		gr5[i] = new TGraph(tr->GetSelectedRows(),tr->GetV1(),tr->GetV2());
		//gr5[i]->SetMarkerColor(2);
		//gr5[i]->SetMarkerStyle(8);
		//gr5[i]->SetMarkerSize(1);
		//gStyle->SetEndErrorSize(3);
		//gr5[i]->Draw("AP");
		MConfigGraphLines(gr5[i],13,0.2);
//<>
MConfigGraphLines(gr5[i],13,2);
//<>

		if(i==0) {
			//gr5[i]->SetTitle("Pion Analysis - #rho_{ch} (from F_{#pi})");
			gr5[i]->GetXaxis()->SetTitle("b (fm)");
			gr5[i]->GetYaxis()->SetTitle("#rho_{ch}(b) (fm^{-2})");
			gr5[i]->GetXaxis()->CenterTitle();
			gr5[i]->GetXaxis()->SetLimits(xlim[0],xlim[1]);
			gr5[i]->GetYaxis()->CenterTitle();
			gr5[i]->GetYaxis()->SetTitleOffset(1.3);
			gr5[i]->GetYaxis()->SetRangeUser(ylim[0],ylim[1]);
			gr5[i]->Draw("al");
		}
		else gr5[i]->Draw("same");
	}
	//Including analytical solutions
	solMonopole->Draw("lsame");
	solDipole->Draw("lsame");
	solMonDi->Draw("lsame");
	//Creating Legend
	TLegend *legb = new TLegend(0.25,0.72,0.85,0.87);
	//legb->SetTextFont(50);
	legb->SetTextSize(0.025);
	legb->AddEntry(gr5[0],"Fitted (monopole + dipole) functions","l");
	legb->AddEntry(solMonDi,"Complete (monopole + dipole) in Q^{2}","l");
	legb->AddEntry(solMonopole,"Monopole - analytical solution","l");
	legb->AddEntry(solDipole,"Dipole - analytical solution","l");
	legb->Draw();
	//Writting on plot
//<>
	TLatex l;
	l.SetTextSize(0.04);
	l.DrawLatex(xlim[0]+(xlim[1]-xlim[0])*0.75,ylim[0]+(ylim[1]-ylim[0])*0.63,Form("#color[4]{n = %d}",n1));
	l.DrawLatex(xlim[0]+(xlim[1]-xlim[0])*0.68,ylim[0]+(ylim[1]-ylim[0])*0.53,Form("#color[2]{Q^{2} = %.1f GeV^{2}}",pow(ZeroBesselJ0(n1)/R1,2.0)));
	TLatex l2;
	l2.SetTextSize(0.02);
	l2.DrawLatex(xlim[0]+(xlim[1]-xlim[0])*0.68,ylim[0]+(ylim[1]-ylim[0])*0.97,"MAPC @ " + da2);
//<>
	// Saving Ge Canvas
	pad2png(c5,Form("ro_ch-F1_n%d.png",n1));
	c5->Print(Form("ro_ch-F1_n%d.eps",n1));
//<>
/*for(j=0;j<nfits;j++) delete ros[j], gr5[j];
delete b, c5, legb, l, l2, tr;
}*/
//<>

	// To create error band
	CreateErrorBand(tr,nfits,solMonopole);

	// Concluding code
	f->Close();
	cout << "<Marco> Done!" << endl;
	return 0;
}

// Function to plot FF
void PlotFF(TNtuple *FF, Double_t R1, Double_t *range) {

	//To use my R1, comment following line. Otherwise using PDG radius
	R1 = 0.672; //fm

	// First simple fits of Q^4F x Q^2
	TCanvas *FF1 = new TCanvas("FF1", "Q2F1xQ2",600,600);
	FF1->SetGrid();
	FF->Draw("Q2:(F1*Q2):(sF1*Q2)","","goff");
	TMultiGraph *FF1m = new TMultiGraph();
	//TGraphErrors *FF1gr[1];
	TGraphErrors *FF1gr[1];
	FF1gr[0] = new TGraphErrors(FF->GetSelectedRows(),FF->GetV1(),FF->GetV2(),0,FF->GetV3());
	//FF1gr[0] = new TGraphErrors(FF->GetSelectedRows(),FF->GetV1(),FF->GetV2(),0,FF->GetV3());
	FF1gr[0]->SetMarkerColor(4);
	FF1gr[0]->SetMarkerStyle(8);
	FF1gr[0]->SetMarkerSize(1);
	//FF1gr[1] = new TGraph(FF->GetSelectedRows(),FF->GetV1(),FF->GetV2());
	//FF1gr[1]->SetLineColor(2);
	//FF1gr[1]->SetLineWidth(2);
	//FF1->Clear();
	FF1m->Add(FF1gr[0],"p");
	//FF1m->Add(FF1gr[1],"c");
	//FF1m->SetTitle("Q^{2}.F_{#pi} x Q^{2}");
	FF1m->Draw("al");

	TF1 *monopole = new TF1("monopole","x*1.0/(1.0+[0]*x)",0.0,range[1]);
	monopole->SetLineColor(3);
	monopole->SetLineStyle(1);
	monopole->SetParameter(0,pow(R1/0.197327,2.0)/6.0);
	monopole->Draw("same");

	TF1 *dipole = new TF1("dipole","x*1.0/pow(1.0+[0]*x,2.0)",0.0,range[1]);
	dipole->SetLineColor(6);
	dipole->SetLineStyle(1);
	dipole->SetParameter(0,pow(R1/0.197327,2.0)/12.0);
	dipole->Draw("same");

	FF1m->GetXaxis()->SetTitle("Q^{2} (GeV^{2})");
	FF1m->GetYaxis()->SetTitle("Q^{2}.F_{#pi}");
	FF1m->GetXaxis()->SetLimits(0.0,1.1*range[1]);
	FF1m->GetXaxis()->CenterTitle();
	FF1m->GetYaxis()->CenterTitle();
	FF1m->GetYaxis()->SetTitleOffset(1.5);
	
	// Legend 
	TLegend *leg = new TLegend(0.14,0.72,0.5,0.87);
	//leg->SetTextFont(50);
	leg->SetTextSize(0.025);
	leg->AddEntry(FF1gr[0],"Data points","lp");
	leg->AddEntry(monopole,Form("Monopole (R = %.3f fm)",sqrt(6.0*pow(0.197327,2.0)*monopole->GetParameter(0))));
	leg->AddEntry(dipole,Form("Dipole (R = %.3f fm)",sqrt(12.0*pow(0.197327,2.0)*dipole->GetParameter(0))));
	leg->Draw();

	pad2png(FF1,"Pion_Q2F1.png");
	FF1->Print("Pion_Q2F1.eps");
	//delete FF1;

	return;
}

//Function to calculate F1 and plot it
Double_t CalcR1(Double_t *range, TNtuple *FF, Double_t MaxRangeForFit, Int_t Mgraph) {
	Double_t rangeSmallQ2[2];
	rangeSmallQ2[0]=range[0];
	rangeSmallQ2[1]=MaxRangeForFit;

	// For F1
	TCanvas *c3 = new TCanvas("c3","Pion_F1_SmallQ2_Fit",600,600);
	//c3->SetLogy();
	c3->SetGrid();
	c3->GetFrame()->SetBorderSize(10);
	// Preparing plots
	FF->Draw("Q2:F1:sF1","","goff");
	TGraphErrors *gr3 = new TGraphErrors(FF->GetSelectedRows(),FF->GetV1(),FF->GetV2(),0,FF->GetV3());
	TString opt = "AP";
	gr3->SetMarkerColor(1);
	gr3->SetMarkerStyle(8);
	gr3->SetMarkerSize(1);
	gStyle->SetEndErrorSize(3);
	gr3->SetTitle("F_{#pi} x Q^{2} - analysis of R1"); 
	gr3->Draw(opt.Data());
	gr3->GetXaxis()->SetTitle("Q^{2} (GeV^{2})");
	gr3->GetYaxis()->SetTitle("F_{#pi} (Q^{2})");
	gr3->GetXaxis()->CenterTitle();
	gr3->GetYaxis()->CenterTitle();
	gr3->GetYaxis()->SetTitleOffset(1.3);
	gr3->GetXaxis()->SetLimits(0.0,0.15);
	gr3->Draw(opt.Data());
	gr3->GetYaxis()->SetRangeUser(0.5,1.0);
	//gr3->GetYaxis()->SetNdivisions(510,kFALSE);
	TF1 *polfit = new TF1("polfit","pol1",rangeSmallQ2[0],rangeSmallQ2[1]);
	polfit->SetLineColor(6);
	polfit->SetLineStyle(1);
	TF1 *dF1dQ2 = new TF1("dF1dQ2","expo",rangeSmallQ2[0],rangeSmallQ2[1]);
	dF1dQ2->SetLineColor(2);
	dF1dQ2->SetLineStyle(1);
	TF1 *monopole = new TF1("monopole","1.0/(1.0+[0]*x)",rangeSmallQ2[0],rangeSmallQ2[1]);
	monopole->SetLineColor(3);
	monopole->SetLineStyle(2);
	TF1 *monopole2 = new TF1("monopole2","[0]/(1.0+[1]*x)",rangeSmallQ2[0],rangeSmallQ2[1]);
	monopole2->SetLineColor(7);
	monopole2->SetLineStyle(2);
	TF1 *dipole = new TF1("dipole","1.0/((1.0+[0]*x)**2)",rangeSmallQ2[0],rangeSmallQ2[1]);
	dipole->SetLineColor(4);
	dipole->SetLineStyle(3);
	cout << "<Marco> Fitting polfit: a+b.x" << endl;
	gr3->Fit("polfit","R0");
	cout << "<Marco> Fitting expo: exp(a+b.x)" << endl;
	gr3->Fit("dF1dQ2","R0");
	cout << "<Marco> Fitting monopole: 1/(1+a.x)" << endl;
	gr3->Fit("monopole","R0");
	cout << "<Marco> Fitting monopole2: a/(1+b.x)" << endl;
	gr3->Fit("monopole2","R0");
	cout << "<Marco> Fitting dipole: 1/(1+a.x)^2" << endl;
	gr3->Fit("dipole","R0");
	polfit->SetRange(range[0],range[1]);
	dF1dQ2->SetRange(range[0],range[1]);
	monopole->SetRange(range[0],range[1]);
	monopole2->SetRange(range[0],range[1]);
	dipole->SetRange(range[0],range[1]);
	polfit->Draw("same");
	dF1dQ2->Draw("same");
	monopole->Draw("same");
	monopole2->Draw("same");
	dipole->Draw("same");

	// Calculating R1
	Double_t R1_poli, R1_mono, R1_mono2, R1_di, R1;
	R1_poli = sqrt(-6.0*pow(0.197327,2)*polfit->GetParameter(1));
	R1_mono = sqrt(6.0*pow(0.197327,2.0)*monopole->GetParameter(0));
	R1_mono2 = sqrt(6.0*monopole2->GetParameter(0)*pow(0.197327,2.0)*monopole2->GetParameter(1));
	R1_di    = sqrt(12.0*pow(0.197327,2.0)*dipole->GetParameter(0));

	R1 = R1_mono;
	// 1 fm = 5.068 GeV^{-1}
	cout << "<Marco> Evaluating R1 = Sqrt(<r^2>)" << endl
             << "<Marco> polfit = " << polfit->GetParameter(0) << " + Q^2 * " << polfit->GetParameter(1) << endl
	     << "        #sqrt{<r^2>} = " << R1_poli << " fm" << endl
	     << "<Marco> dF1dQ2 = " << dF1dQ2->GetParameter(1) << endl
	     << "<Marco> Monopole = " << R1_mono << endl
	     << "<Marco> Monopole2 = " << R1_mono2 << endl
	     << "<Marco> Dipole = " << R1_di << endl << endl 
	     << "<Marco> VERY IMPORTANT:" << endl
             << "<Marco> Using R1 = " << R1 << endl << endl;
	
	// Legend F1
	//TLegend *legF1 = new TLegend(0.3,0.8,0.95,0.92);
	TLegend *legF1 = new TLegend(0.14,0.14,0.7,0.50);
	//legF1->SetTextFont(50);
	legF1->SetTextSize(0.025);
	legF1->SetHeader(Form("Fitted Q^{2} range: [%.2f, %.2f] GeV^{2}",rangeSmallQ2[0],rangeSmallQ2[1]));
	//if(PointsOrFits==1) legF1->AddEntry(gr3,"Fitted data");
	//else 
	legF1->AddEntry(gr3,"Data points","p");
	legF1->AddEntry(polfit,Form("Linear polynomium \t#sqrt{<r^{2}>} = %.3f fm",R1_poli));
	legF1->AddEntry(dF1dQ2,Form("exp(a+b.Q^{2})"));// \t#sqrt{<r^{2}>} = %.3f fm",R1_exp));
	legF1->AddEntry(monopole,Form("#frac{1}{1+A.Q^{2}} \t#sqrt{<r^{2}>} = %.3f fm",R1_mono));
	legF1->AddEntry(monopole2,Form("#frac{A}{1+B.Q^{2}} \t#sqrt{<r^{2}>} = %.3f fm",R1_mono2));
	legF1->AddEntry(dipole,Form("#frac{1}{(1+A.Q^{2})^{2}} \t#sqrt{<r^{2}>} = %.3f fm",R1_di));
	legF1->Draw();
/*
	TPaveText *ptF1 = new TPaveText(0.14,0.14,0.5,0.4,"NDC");
	ptF1->AddText(" ");
	ptF1->AddText("Fitted function: exp(a*Q^{2} + b), Q^{2} #rightarrow 0");
	ptF1->AddText("where a = #frac{d log(F)}{d Q^{2}}");
	ptF1->AddText("so, R1 = 5 * #sqrt{|<b^{2}>|} = 5 * #sqrt{-4*a}");
	ptF1->AddText(Form("R1 = %.3f fm",R1));
	ptF1->SetLabel("R evaluation");
	ptF1->Draw();
*/

	// Saving F1 Canvas
	pad2png(c3,"Pion_F1_Rfit.png");
	c3->Print("Pion_F1_Rfit.eps");
	//delete c3;

	return R1;
}

TF1 *FF1fit(TNtuple *FF, Double_t *range, Int_t q2f) {

	//Int_t q2f=1; //q2f=0 plot (F)x(Q2). q2f=1 plot (Q2.F)x(Q2)

	// Functions F1 and F2
	TCanvas *fF1 = new TCanvas("fF1","Fit Fpi",600,600);
	fF1->SetGrid();
	fF1->GetFrame()->SetBorderSize(10);
	FF->Draw("Q2:F1:sF1","","goff");
	TGraphErrors *pts = new TGraphErrors(FF->GetSelectedRows(),FF->GetV1(),FF->GetV2(),0,FF->GetV3());
	FF->Draw("Q2:(Q2*F1):(Q2*sF1)","","goff");
	TGraphErrors *ptsPlot = new TGraphErrors(FF->GetSelectedRows(),FF->GetV1(),FF->GetV2(),0,FF->GetV3());
	pts->SetTitle("Pion data - #rho_{ch} (from F1)");
	ptsPlot->SetTitle("Pion data - #rho_{ch} (from F1)");
	pts->GetXaxis()->SetTitle("Q^{2} (GeV^{2})");
	ptsPlot->GetXaxis()->SetTitle("Q^{2} (GeV^{2})");
	pts->GetYaxis()->SetTitle("F_{#pi} (GeV^{2})");
	ptsPlot->GetYaxis()->SetTitle("Q^{2}.F_{#pi} (GeV^{2})");
	pts->GetYaxis()->SetTitleOffset(1.4);
	ptsPlot->GetYaxis()->SetTitleOffset(1.4);
	pts->GetYaxis()->SetRangeUser(0.0,1.0);
	ptsPlot->GetYaxis()->SetRangeUser(0.0,1.0);
	pts->SetMarkerColor(1);
	ptsPlot->SetMarkerColor(1);
	pts->SetMarkerStyle(8);
	ptsPlot->SetMarkerStyle(8);
	pts->SetMarkerSize(1);
	ptsPlot->SetMarkerSize(1);
	gStyle->SetEndErrorSize(3);
	//FF->Draw("Q2:(Q2*F1)");
	if(q2f) ptsPlot->Draw("ap");
	else pts->Draw("ap");
	// Legend Ge
	//TLegend *legfF = new TLegend(0.33,0.7,0.95,0.9);
	TLegend *legfF = new TLegend(0.16,0.67,0.78,0.87);
	legfF->AddEntry(ptsPlot,"Experimental data","p");
	//legfF->SetTextFont(72);
	legfF->SetTextSize(0.025);
	TF1 *fitF1[2];
	TF1 *fitF1Plot[2];
	//start loop here	
	fitF1[0] = new TF1("fitF1_0","1.0/(1.0+x*[0])",range[0],range[1]);
	fitF1[1] = new TF1("fitF1_1","1.0/((1.0+x*[0])**2)",range[0],range[1]);
	fitF1[2] = new TF1("fitF1_2","[0]/(1.0+x*[1]) + (1.0-[0])/((1.0+x*[2])**2)",range[0],range[1]);
	fitF1[2]->SetParLimits(0,0.0,1.0);
	//fitF1[2]->SetParLimits(1,0.0,10.0);
	//fitF1[2]->SetParLimits(2,0.0,10.0);
	//Double_t par_fitF1[3] = {0.25, 0.2, 0.2};
	Double_t par_fitF1[3] = {0.25, 0.2, 2.75};
	fitF1[2]->SetParameters(par_fitF1);

	TF1 **Plot;
	if(q2f) Plot = fitF1Plot;
	else Plot = fitF1;

	pts->Fit("fitF1_0","R0");
	pts->Fit("fitF1_1","R0");
	pts->Fit("fitF1_2","R0");
	fitF1[1]->GetParameters(par_fitF1);
	fitF1Plot[0] = new TF1("fitF1P_0","x*fitF1_0",range[0],range[1]);
	Plot[0]->SetLineColor(2);
	Plot[0]->SetLineStyle(3);
	legfF->AddEntry(Plot[0],"Monopole #rightarrow #frac{1}{(1+A*Q^{2})}","l");
	fitF1Plot[1] = new TF1("fitF1P_1","x*fitF1_1",range[0],range[1]);
	Plot[1]->SetLineColor(3);
	Plot[1]->SetLineStyle(2);
	legfF->AddEntry(Plot[1],"Dipole #rightarrow #frac{1}{(1+A*Q^{2})^{2}}","l");
	fitF1Plot[2] = new TF1("fitF1P_2","x*fitF1_2",range[0],range[1]);
	Plot[2]->SetLineColor(4);
	Plot[2]->SetLineStyle(1);
	legfF->AddEntry(Plot[2],"Combination #rightarrow A.#frac{1}{(1+B*Q^{2})} + (1-A).#frac{1}{(1+C*Q^{2})^{2}}");

	Plot[0]->Draw("acsame");
	Plot[1]->Draw("acsame");
	Plot[2]->Draw("acsame");
	//end loop
	legfF->Draw();
	if(q2f)	{
		pad2png(fF1,"Pion_FF_Fitting_q2f1.png");
		fF1->Print("Pion_FF_Fitting_q2f1.eps");
	}
	else {
		pad2png(fF1,"Pion_FF_Fitting.png");
		fF1->Print("Pion_FF_Fitting.eps");
	}

	return fitF1[2]; //choose here which function to return
}

int main(int argc, char *argv[]) {
	TString inData = argv[1];

	return UncertaintyPion(inData);
}
