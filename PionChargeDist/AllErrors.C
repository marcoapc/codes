#include "../MLibraries/MallIncludes.h"
#include "../MLibraries/MConfigGraphs.h"

const int maxnumfiles = 10;

//int AllErrors(TTree *tr, Int_t n, TF1 *model, TString var="ro") {
int AllErrors(TString *inData, Int_t nfiles) {

	// Preparing variables
	Int_t i,j,k;
	TBranch *br;
	Float_t T;
	Int_t npts;
	Int_t n; // Number of branches ro%d (number of fits)
	Double_t *r, *a2, *a, *Tstd, *Tmean, *Tbandmin, *Tbandmax, *Tbandshape;
	TString var = "ro"; // core name of the branches with fit data. Example: ro1, ro2, ro3, ... must have var = "ro"
	TGraph *grmean[maxnumfiles];
	TGraph *grmin[maxnumfiles];
	TGraph *grmax[maxnumfiles];
	TGraph *grshade[maxnumfiles];

	// Preparing for plots
	TCanvas *c1 = new TCanvas("c1","c1",600,600);
	MConfigCanvas(c1);
	TLegend *leg = new TLegend(0.43,0.72,0.95,0.98);
	leg->SetTextSize(0.025);

	// Exact solution for monopole
	TF1 *solMonopole = new TF1("solMonopole","(1.0/(2.0*TMath::Pi()))*(6.0/pow([0],2))*TMath::BesselK0(sqrt(6.0)*x/[0])",0.0,2.0); //[0] is R1 in units of b
        Double_t par_solMonopole[1];
        par_solMonopole[0] = 0.672; //in fm
        solMonopole->SetParameters(par_solMonopole);
        MConfigLines(solMonopole,2);

	// Opening files and getting trees
	TFile *f[maxnumfiles];
	TTree *t[maxnumfiles];
	for(k=0;k<nfiles;k++) {
		cout << endl << "<Marco> Opening file: " << inData[k] << endl;
		f[k] = new TFile(inData[k]);
		
		//Getting tree name
		TIter nextkey(f[k]->GetListOfKeys());
		TKey *key, *oldkey=0;
		while((key = (TKey*)nextkey())) {
			TObject *obj = key->ReadObj();
			if(obj->IsA()->InheritsFrom(TTree::Class())) {
				t[k] = (TTree*)obj;
				break;
			}
		}

		cout << "\tGot the tree:" << endl;
		t[k]->ls();

		//////////////////////////////
		// Constructing Error Bands //
		//////////////////////////////
	
		// Getting radius information (b)
		br = t[k]->GetBranch("b");
		br->SetAddress(&T);
		npts = br->GetEntries();
		
		// Number of fits (from the number of branches minus 1)
		n = t[k]->GetNbranches() - 1; // 1 from the branch "b"
		cout << Form("\tFound %d fits",n) << endl;
	
		// Preparing variables to calculate stdev and mean
		r = (Double_t *)malloc(npts*sizeof(Double_t));
		a2 = (Double_t *)malloc(npts*sizeof(Double_t));
		a = (Double_t *)malloc(npts*sizeof(Double_t));
		Tstd = (Double_t *)malloc(npts*sizeof(Double_t));
		Tmean = (Double_t *)malloc(npts*sizeof(Double_t));
		Tbandmin = (Double_t *)malloc(npts*sizeof(Double_t));
		Tbandmax = (Double_t *)malloc(npts*sizeof(Double_t));
		Tbandshape = (Double_t *)malloc(npts*sizeof(Double_t));

		// Zeroing some incremental variables
		for(j=0;j<npts;j++) {
			a2[j]=0.0;
			a[j]=0.0;
		}

		// Reading radius
		for(j=0;j<npts;j++) {
			br->GetEntry(j);
			r[j] = T;
		}

		// reading fits - i runs over all "n" fits
		for(i=0;i<n;i++) {
			br = t[k]->GetBranch(var+Form("%d",i));
			br->SetAddress(&T);
			// j runs over the radius (b)
			for(j=0;j<npts;j++) {
				br->GetEntry(j);
				a2[j] += T*T;
				a[j] += T;
			}
		}

		// Calculating mean and std for band construction
		for(j=0;j<npts;j++) {
			Tmean[j] = a[j]/n; //for computing Tstd
			Tstd[j] = sqrt((a2[j]-2.0*a[j]*Tmean[j]+n*Tmean[j]*Tmean[j])/(n-1.0));
			//Tmean[j] = abs(solMonopole->Eval(r[j])-Tmean[j]);

			Tbandmin[j] = Tmean[j];
			Tbandmax[j] = Tmean[j]+Tstd[j];
			//Tbandmin[j] = Tmean[j]-Tstd[j];
			//Tbandmax[j] = Tmean[j]+Tstd[j];
		}

		// Plotting
		grmean[k] = new TGraph(npts,r,Tmean);
		grmin[k]  = new TGraph(npts,r,Tbandmin);
		grmax[k]  = new TGraph(npts,r,Tbandmax);
		grshade[k] = new TGraph(2*npts);
		for(j=0;j<npts;j++) {
			grshade[k]->SetPoint(j,r[j],Tbandmax[j]);
			grshade[k]->SetPoint(npts+j,r[npts-j-1],Tbandmin[npts-j-1]);
		}
		grshade[k]->SetFillStyle(3144);
		grshade[k]->SetFillColor(k+1);
		if(k==0) {
			MConfigAxis(grmax[k],"","Radius / fm","Charge Distribution / fm^{-2}");
			//grmin->GetXaxis()->SetLimits(0.0,2.0);
			//grmin->GetYaxis()->SetRangeUser(0.0,5.0);
			grmax[k]->Draw("al");
		}
		else grmax[k]->Draw("lsame");
		grmean[k]->SetLineColor(2);
		grmean[k]->SetLineWidth(2);
		grmin[k]->Draw("lsame");
		grshade[k]->Draw("fsame");
		grmean[k]->Draw("lsame");

/*
	// Creating error graph
	const int ngraphs = 4;
	Int_t n1[ngraphs] = {3,3,5,6}; // terms in the expansion for each type
	TCanvas *c2 = new TCanvas("c2","c2",600,600);
	Double_t *err[ngraphs];
	for(i=0;i<ngraphs;i++) {
		Double_t *err[i] = (Double_t *)malloc(npts*sizeof(Double_t));
		for(j=0;j<npts;j++) 
			err[i][j] = model->Eval(r[j]) - Tmean[r]
*/
		// Free memory
		free(r);
		free(a2);
		free(a);
		free(Tstd);
		free(Tmean);
		free(Tbandmin);
		free(Tbandmax);
		free(Tbandshape);
	}
	
	// Saving canvas
	leg->AddEntry(grshade[0],"Present data fit","f");
	leg->AddEntry(grshade[1],"Prediction for JLab 12GeV","f");
	leg->Draw();
	solMonopole->Draw("lsame");
	pad2png(c1,"DistBand.png");
	c1->Print("DistBand.eps");

	return 0;
}

int main(int argc, char *argv[]) {
	if(argc>maxnumfiles) {
		cout << "MError: Max number of files is " << maxnumfiles << ". Cannot work with " << argc << " files." << endl;
		return 1;
	}

	TString inData[maxnumfiles];
	for(Int_t i=1; i<argc; i++) inData[i-1] = argv[i];

	return AllErrors(inData,argc-1);
}
