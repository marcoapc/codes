#include<TRandom.h>

// Function to be fitted (can change the function, not the input arguments!). Change number of parameters in the function creation if it changes here.
Double_t MfitF(Double_t *x, Double_t *par) {
	//// y = [0] * 1.0/(1.0+[1]*x) + (1-[0]) * 1.0/pow(1.0+[2]*x,2.0)
	////return (par[0] + x[0]*par[1] + x[0]*x[0]*par[2]);
	return (par[0] * 1.0/(1.0+par[1]*x[0]) + (1.0-par[0]) * 1.0/pow(1.0+par[2]*x[0],2.0));
	//return (par[0]/pow(1.0+par[1]*x[0],2));
}

Double_t MfitFQ2(Double_t *x, Double_t *par) {
	return (x[0]*MfitF(x,par));
}

TF1 *CreateTF(Double_t xmin, Double_t xmax, Int_t npar) {
	TF1 *f = new TF1("f",MfitFQ2,xmin,xmax,npar);
	return f;
}

// Function that receives a TNtuple "in" with columns "x", "y" and "sy" and returns a pointer to a created TGraph "out". The pointsof out are created for each (x,y) point of in, in the position (x,y+). If one wants to also use errorbars in the fitting, change TGraph to TGraphErrors.
//TGraphErrors *MakeNewPoints(TNtuple *in, TRandom *r0) {
TGraph *MakeNewPoints(TNtuple *in, TRandom *r0) {
	Int_t i, nEntries;
	Float_t x, y, sy;
	Double_t *xnew, *ynew, *synew;

	// Preparing x, y, sy
	//in->SetBranchAddress("x",&x);
	//in->SetBranchAddress("y",&y);
	//in->SetBranchAddress("sy",&sy);
	// Preparing x, y, sy
	in->SetBranchAddress("Q2",&x);
	in->SetBranchAddress("F1",&y);
	in->SetBranchAddress("sF1",&sy);

	// Reading number of elements in "in"
	nEntries = in->GetEntries();

	// Allocating memory for point
	xnew = (Double_t*)TStorage::Alloc(nEntries*sizeof(Double_t));
	ynew = (Double_t*)TStorage::Alloc(nEntries*sizeof(Double_t));
	synew = (Double_t*)TStorage::Alloc(nEntries*sizeof(Double_t));

	// Creating points
	for(i=0; i<nEntries; i++) {
		in->GetEntry(i);
		xnew[i]  = x;
		ynew[i]  = r0->Gaus(y,sy);
		synew[i] = sy;
		//cout << "  i=" << i << " -> x=" << xnew[i] << " \t y=" << y << " \tsy=" << sy << " \tynew=" << ynew[i] << endl;
	}
		
	//TGraphErrors *out = new TGraphErrors(nEntries,xnew,ynew,0,synew);
	TGraph *out = new TGraph(nEntries,xnew,ynew);
	MConfigPoints(out,3);

	return out;
}

// Function that receive a TGraph "gr" and a TCanvas "c", fit the curve fiven by the function MfitF (showing it in c). Return a pointer to a vector with the fitted function.
TF1 *ExecuteFit(TGraph *gr, TCanvas *c, Double_t *par0, Int_t quiet=0) {
	Int_t i;
	Double_t xmin, xmax;
	//Double_t *par;
	Float_t nEntries;

	const int npar = 3;

	nEntries = gr->GetN();
	xmin = TMath::MinElement(nEntries,gr->GetX());
	xmax = TMath::MaxElement(nEntries,gr->GetX());
	//cout << "xmin=" << xmin << " \txmax=" << xmax << endl;

	TF1 *f = new TF1("f",MfitF,xmin,xmax,npar);
	f->SetParLimits(0,0.0,1.0);
	f->SetParLimits(1,0.0,20.0);
	f->SetParLimits(2,0.0,20.0);
	MConfigLines(f, 13, 0.2);

	f->SetParameters(par0);
	if(quiet==0) gr->Fit("f","R0");
	else gr->Fit("f","QR0");
	f->Draw("same");

	//delete f;
	return f;
}
