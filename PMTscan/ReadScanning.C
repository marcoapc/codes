#include<cstdio>
#include<stdlib.h>
#include<TCanvas.h>
#include<TStyle.h>
#include<TH1.h>
#include<TF1.h>
#include<TMath.h>
#include<Riostream.h>

int ReadScanning(TString inFile, TString outFile, Double_t gain=0.0, Int_t nlin=25, Int_t ncol=25, Int_t nptos=500, Int_t GraphOpt=0) {
        /* Running command example in root:
	.x ReadScanning.C("600.txt","600.root",1,100,100,1000)
        */

	// Variables
	Int_t runDir=-1;
	Int_t xlim[2];
        Double_t zlim[2];
	char nzi[5], nzj[5], NameOutput[50];

	/* IMPORTANT INFORMATION */
	//Int_t refCol=2, sigCol=1; // To analyze Reference PMT, uncomment this line and comment the next one.
	Int_t refCol=0, sigCol=1; // To analyze Scanned PMT, uncomment this line and comment the previous one.
	Int_t zigzag=0; Int_t AnalyzePED=0;
	Int_t useZLim=0; zlim[0]=0; zlim[1]=6;
	Int_t DeltaPosPeak=1;
	xlim[0]=1; xlim[1]=150;

	// Getting current directory
	TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
	dir.ReplaceAll("ReadScanning.C","");
	dir.ReplaceAll("/./","/");
	printf("Current Directory: %s\n",dir.Data());

	// Opening input file
	ifstream in;
	in.open(Form("%s%s",dir.Data(),inFile.Data()));
	if(in.fail()) {
		printf(" Marco: Error opening file.\n");
		printf("        %s%s\n",dir.Data(),inFile.Data());
		return 1;
	}

	// Creating output file
	TFile froot(outFile.Data(),"recreate");
	TTree *T = new TTree("T","PMT_scanning");

	// Preparing histograms
	Int_t nbins=1024;
	TH1F *h = new TH1F("h","",nbins,1,nbins);
	h->SetFillColor(41);
	Double_t sizeX = 4.5;
	Double_t sizeY = 4.5;
	inFile.ReplaceAll(".txt","");
	TH2F *pos = new TH2F("pos",inFile.Data(),ncol,1,ncol,nlin,1,nlin);
	T->Branch("pos","TH2F",&pos,32000,0);

	// Preparing fits (pedestal, SEP)
	//TF1 *pedestal = new TF1("pedestal","gaus(0)",0,0);
	Double_t par[3];
	TF1 *fSEP = new TF1("fSEP","gaus(0)",0,0);
	fSEP->SetNpx(200);
	Double_t parSEP[3];
	TF1 *f = new TF1("f","gaus(0)",1,nbins);
	f->SetNpx(5000);
	Double_t parF[3];
	f->SetLineColor(3);

	// Preparing Canvas
	Double_t Ww = 600;
	Double_t Wh = 600;

	// Preparing to read data
	Int_t sig, ref, flagFirst=1, ntot=0;
	Int_t i, j, n=0, nsomados=0;
	i=nlin; j=ncol;
	Double_t somaA=0.0;

	// Reading 1st set of data, to analyze SEP and pedestal
	while(ntot<nptos) {
		if(refCol==1 && sigCol==2) in >> ref >> sig;
		else if(refCol==0 && sigCol==1) in >> sig;
		else if(refCol==2 && sigCol==1) in >> sig >> ref;
		else {
			printf(" Marco: problems with refCol refSig.\n");
			return 1;
		}
		//printf("ref: %d\tsig:%d\n",ref,sig);
		//return 0;
		printf("\r%5d x %5d - %d     ",i,j,n+1);
		if(in.eof()) break;
		if(in.fail()) {
			printf(" Marco: Fail reading input file.\n");
			return 1;
		}
		n++; ntot++;
	
		h->Fill(((double)sig));
	}
	
	// Analyzing pedestal	
	par[1]=0.0; // pedestal center
	if(AnalyzePED==1) {
		Double_t soma=0.0;
		for(int k1=h->GetMaximumBin()-DeltaPosPeak;k1 <= h->GetMaximumBin()+DeltaPosPeak;k1++) {
			par[1]+=h->GetBinContent(k1)*k1;
			soma+=h->GetBinContent(k1);
		}
		par[1]/=soma;
		/* For fitting a gaussian 
		par[0]=20*h->GetMaximum();
		par[1]=h->GetMaximumBin();
		par[2]=0.03;
		pedestal->SetParameters(&par[0]);
		Int_t DeltaPosPeak=3;
		pedestal->SetParLimits(0,0.0,10000);
		pedestal->SetParLimits(1,par[1]-DeltaPosPeak,par[1]+DeltaPosPeak);
		pedestal->SetParLimits(2,0.0,0.3);
		pedestal->SetRange(par[1]-DeltaPosPeak,par[1]+DeltaPosPeak);
		h->Fit("pedestal","MNB");
		pedestal->SetNpx(100);
		pedestal->FixParameter(1,pedestal->GetParameter(1));
		h->Fit("pedestal","MRB");
		//printf("A_0: %f\nx_0: %f\n",par[0],par[1]);
		pedestal->GetParameters(par);
		*/
		//h->Draw();
		//pedestal->Draw("same");
		printf("\nPedestal center: %f\n",par[1]);
		// preparing canvas, TH and other variables
		if(GraphOpt==1) {
			TCanvas *c1 = new TCanvas("c1","",Ww,Wh);
			h->GetXaxis()->SetRangeUser(xlim[0],xlim[1]);
			h->GetYaxis()->SetRangeUser(0,nptos/2.0);
			// Duplicating original histogram
			TH1F *horig = (TH1F*)h->Clone("horig");
		}
	}

	// Removing pedestal for SEP analyzis
	for(int k2=((int)par[1])-DeltaPosPeak;k2 <= ((int)par[1])+DeltaPosPeak;k2++) h->SetBinContent(k2,0.0);
	// Fitting SEP if there is no input calib
	if(gain==0.0) {
		// fitting a gaussian
		parSEP[0]=h->GetMaximum();
		parSEP[1]=h->GetMaximumBin();
		parSEP[2]=3.0;
		fSEP->SetParameters(&parSEP[0]);
		fSEP->SetParLimits(0,0,20*nptos);
		fSEP->SetParLimits(2,0,20);
		fSEP->SetRange(parSEP[1]-10.0,parSEP[1]+10.0);
		h->Fit("fSEP","MNB");
		fSEP->GetParameters(parSEP);
		//h->Draw("same");
		printf("SEP center: %f\n",parSEP[1]);
		gain=(parSEP[1]-par[1]); // FIX IT!!!!!
		printf("  Marco: WARNING! Check the gain equation (attenuation, ADC features, ...)");
	}
	
	// Saving data on file
	if(h->GetMean() > par[1]) pos->SetBinContent(i,j,(h->GetMean()-par[1])/gain);
	printf("\n%d x %d       \n",i,j);
	T->Fill();
			
	// Drawing
	if(GraphOpt==1) { 
		horig->Draw();
		if(n!=nptos) f->Draw("same");
		else fSEP->Draw("same");
		// Saving image
		switch(((int)TMath::Log10(i))+1) {
			case 1:
				sprintf(nzi,"00");
				break;
			case 2:
				sprintf(nzi,"0");
				break;
			default:
				sprintf(nzi,"");
				break;
		}
		switch(((int)TMath::Log10(j))+1) {
			case 1:
				sprintf(nzj,"00");
				break;
			case 2:
				sprintf(nzj,"0");
				break;
			default:
				sprintf(nzj,"");
				break;
		}
		sprintf(NameOutput,"hist_%s%dx%s%d.png",nzi,i,nzj,j);
		c1->SaveAs(NameOutput);
	}

	if(GraphOpt==1) {
		// Hist for the next position
		delete h;
		TH1F *h = new TH1F("h","",nbins,1,nbins);
		h->SetFillColor(41);
	}

	// updating position
	if(zigzag==1) {
		if(runDir==-1 && j!=1) j--;
		else if(runDir==1 && j!=ncol) j++;
		else if(runDir==-1 && j==1) {i--; runDir=1;}
		else if(runDir==1 && j==ncol) {i--; runDir=-1;}
	}
	else {
		if(runDir==-1 && j!=1) j--;
		else if(runDir==1 && j!=ncol) j++;
		else if(runDir==-1 && j==1) {i--; j=ncol;}
		else if(runDir==1 && j==ncol) {i--; j=1;}
	}

	somaA=0.0;
	nsomados=0;
	n=0;

	// Reading all the next data
	while(i>=1 && j>=1) {
		if(refCol==1 && sigCol==2) in >> ref >> sig;
		else if(refCol==0 && sigCol==1) in >> sig;
		else if(refCol==2 && sigCol==1) in >> sig >> ref;
		else {
			printf(" Marco: problems with refCol refSig.\n");
			return 1;
		}
		//printf("\r%5d x %5d - %d     ",i,j,n+1);
		//printf("\r%d",n+1);
		n++; ntot++;
	
		// Inserting readed point into histogram
		if(!in.eof()) {// && sig>par[1]+DeltaPosPeak) { // remove pedestal from mean
			somaA += (double)sig;
			nsomados++;
		}
	
		if(n==nptos || in.eof()) {
			if(nsomados!=0) pos->SetBinContent(i,j,(somaA/nsomados-par[1])/gain);
			else pos->SetBinContent(i,j,0.0);
			// updating position
			if(zigzag==1) {
				if(runDir==-1 && j!=1) j--;
				else if(runDir==1 && j!=ncol) j++;
				else if(runDir==-1 && j==1) {i--; runDir=1;}
				else if(runDir==1 && j==ncol) {i--; runDir=-1;}
			}
			else {
				if(runDir==-1 && j!=1) j--;
				else if(runDir==1 && j!=ncol) j++;
				else if(runDir==-1 && j==1) {i--; j=ncol;}
				else if(runDir==1 && j==ncol) {i--; j=1;}
			}
	
			printf("\r%3d x %3d       ",i,j);
			n=0; nsomados=0; somaA=0.0;
		}
		if(in.eof()) break;
		if(in.fail()) {
			printf(" Marco: Fail reading input file.\n");
			return 1;
		}
	}

	// Saving output file
	T->Print();
	froot.Write();

	// Final plot
	inFile.ReplaceAll(".txt","");
	TCanvas *c2 = new TCanvas("c2","c2",Ww,Wh);
	gStyle->SetPalette(1);
	pos.SetXTitle("Horizontal position (a.u.)");
	pos.SetYTitle("Vertical position (a.u.)");
	if(useZLim==1) {
		pos->SetMinimum(zlim[0]);
		pos->SetMaximum(zlim[1]);
	}
	//pos.SetAxisRange(0,pos->GetMaximum(),"Z");
	gStyle->SetPadTickY(1);
	gStyle->SetPadTickX(1);
	//TPaletteAxis *palette = (TPaletteAxis*)pos->GetListOfFunctions()->FindObject("palette");
	//palette->SetName("Number photo-electrons");
	//pos->GetZaxis(50,0.0,100.0);
	//pos->Draw("lego2 0");
	pos->Draw("COLZ");
	gStyle->SetOptStat(0);

	// Saving
	outFile.ReplaceAll(".root",".png");
	c2->SaveAs(outFile.Data());

	printf("n=%d\n",n);

	in.close();
	return 0;
}
