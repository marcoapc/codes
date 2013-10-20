// - n is the number of fits

int CreateErrorBand(TTree *tr, Int_t n, TF1 *model, TString var="ro") {

	// Preparing variables
	Int_t i,j;
	TBranch *br;
	Float_t T;

	// Getting radius information (b)
	br = tr->GetBranch("b");
	br->SetAddress(&T);
	Int_t npts = br->GetEntries();

	// Preparing variables to calculate stdev and mean
	Double_t *r = (Double_t *)malloc(npts*sizeof(Double_t));
	Double_t *a2 = (Double_t *)malloc(npts*sizeof(Double_t));
	Double_t *a = (Double_t *)malloc(npts*sizeof(Double_t));
	Double_t *Tstd = (Double_t *)malloc(npts*sizeof(Double_t));
	Double_t *Tmean = (Double_t *)malloc(npts*sizeof(Double_t));
	Double_t *Tbandmin = (Double_t *)malloc(npts*sizeof(Double_t));
	Double_t *Tbandmax = (Double_t *)malloc(npts*sizeof(Double_t));
	Double_t *Tbandshape = (Double_t *)malloc(npts*sizeof(Double_t));

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
		br = tr->GetBranch(var+Form("%d",i));
		br->SetAddress(&T);
		// j runs over the radius (b)
		for(j=0;j<npts;j++) {
			br->GetEntry(j);
			//if(j==27) cout << "i=" << i << "\tT=" << T << endl;
			a2[j] += T*T;
			a[j] += T;
		}
	}

	// Calculating mean and std for band construction
	for(j=0;j<npts;j++) {
		Tmean[j] = a[j]/n;
		Tstd[j] = sqrt((a2[j]-2.0*a[j]*Tmean[j]+n*Tmean[j]*Tmean[j])/(n-1.0));
		//cout << "j=" << j << "\tr=" << r[j] << "\tTmean=" << Tmean[j] << "\tTstd=" << Tstd[j] << endl;

		Tbandmin[j] = Tmean[j]-Tstd[j];
		Tbandmax[j] = Tmean[j]+Tstd[j];
	}

	// Plotting
	TCanvas *c1 = new TCanvas("c1","c1",600,600);
	MConfigCanvas(c1);
	TGraph *grmean = new TGraph(npts,r,Tmean);
	TGraph *grmin  = new TGraph(npts,r,Tbandmin);
	TGraph *grmax  = new TGraph(npts,r,Tbandmax);
	TGraph *grshade = new TGraph(2*npts);
	for(j=0;j<npts;j++) {
		grshade->SetPoint(j,r[j],Tbandmax[j]);
		grshade->SetPoint(npts+j,r[npts-j-1],Tbandmin[npts-j-1]);
	}
	grshade->SetFillStyle(3144);
	grshade->SetFillColor(2);
	grmin->SetTitle("");
	grmin->GetXaxis()->SetTitle("Radius / fm");
	//grmin->GetXaxis()->SetLimits(0.0,2.0);
	grmin->GetYaxis()->SetTitle("Charge Distribution / fm^{-2}");
	grmin->GetYaxis()->SetTitleOffset(1.3);
	//grmin->GetYaxis()->SetRangeUser(0.0,5.0);
	grmean->SetLineColor(2);
	grmean->SetLineWidth(2);
	grmin->Draw("al");
	grmax->Draw("lsame");
	grshade->Draw("fsame");
	grmean->Draw("lsame");

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

	pad2png(c1,"DistBand.png");
	c1->Print("DistBand.eps");
*/
	return 0;
}
		
