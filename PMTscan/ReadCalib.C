#include "../MLibraries/MallIncludes.h"
#include "../MLibraries/MConfigGraphs.h"

int ReadCalib(TString inFile="", Double_t HV=0.0, Double_t xmin=40.0, Double_t xmax=0.0, Double_t ymax=0.0, TString PMT="") {

	Int_t i;
	if(PMT.IsNull()) PMT = "unknown PMT";

	TNtuple *fin = new TNtuple("fin","Input Data","x");
	fin->ReadFile(inFile);

	//Creating Title
	TString title = inFile(inFile.Last('/')+1,inFile.Last('.')-inFile.Last('/')-1);
	cout << title << endl;

	Double_t npts = fin->GetEntries();
	cout << "<Marco> Number of entries: " << npts << endl;

	Double_t rangeX[2];
	rangeX[0] = xmin; //0.0; //fin->GetMinimum("x");
	if(xmax<=xmin) rangeX[1] = fin->GetMaximum("x");
	else rangeX[1] = xmax;

	TH1F *h = new TH1F(PMT,"Calibration",1101,0,1100);

	TString name = "h1";

	TCanvas *c1 = new TCanvas("c1","Data",600,600);
	fin->Draw("x>>"+PMT,"","goff");
	MConfigHist(h,"run " + title,"Channel number","counts");
	//c1->SetLogy();
	if(ymax!=0.0) h->GetYaxis()->SetRangeUser(0.0,ymax);
	else ymax=h->GetBinContent(h->GetMaximumBin());
	h->GetXaxis()->SetRangeUser(rangeX[0],rangeX[1]);
	h->Draw();

	//Fitting pedestal
	Double_t deltaPed = 5.0;
	TF1 *fped = new TF1("fped","gaus(0)",h->GetMaximumBin()-deltaPed,h->GetMaximumBin()+deltaPed);
	Double_t parPed[3];
	parPed[0] = h->GetBinContent(h->GetMaximumBin())*deltaPed;
	parPed[1] = h->GetMaximumBin();
	parPed[2] = 0.2;
	//cout << "0 = " << parPed[0] << endl << "1 = " << parPed[1] << endl;
	fped->SetParameters(parPed);
	MConfigLines(fped,2);
	h->Fit("fped","MR","same");

	//Fitting SEP
	Double_t deltaSEP = 20.0;
	Double_t firstCh = fped->GetParameter(1) + 20.0*fped->GetParameter(2);
	TH1F *hnoPed = new TH1F("hnoPed","Calibration",1101,0,1100);
	fin->Draw("x>>hnoPed",Form("x>%f",firstCh),"goff");
	hnoPed->SetFillColor(11);
	hnoPed->Draw("same");
	TF1 *fSEP = new TF1("fSEP","gaus(0)",hnoPed->GetMaximumBin()-deltaSEP,hnoPed->GetMaximumBin()+deltaSEP);
	Double_t parSEP[3];
	parSEP[0] = h->GetBinContent(h->GetMaximumBin())*deltaSEP/2.0;
	parSEP[1] = h->GetMaximumBin();
	parSEP[2] = 5.0;
	cout << "0 = " << parSEP[0] << endl << "1 = " << parSEP[1] << endl;
	fSEP->SetParameters(parSEP);
	fSEP->SetParLimits(1,hnoPed->GetMaximumBin()-deltaSEP,hnoPed->GetMaximumBin()+deltaSEP);
	fSEP->SetParLimits(2,1.0,deltaSEP);
	MConfigLines(fSEP,3);
	hnoPed->Fit("fSEP","MR","same");

	TLatex l;
	l.SetTextSize(0.04);
	if(HV!=0.0) l.DrawLatex(rangeX[0]+0.6*(rangeX[1]-rangeX[0]),0.85*ymax,Form("%.2f kV",HV));
	l.DrawLatex(rangeX[0]+0.7*(rangeX[1]-rangeX[0]),0.65*ymax,Form("Ped = %.3f",fped->GetParameter(1)));
	l.DrawLatex(rangeX[0]+0.7*(rangeX[1]-rangeX[0]),0.6*ymax,Form("SEP = %.3f",fSEP->GetParameter(1)));
	l.DrawLatex(rangeX[0]+0.7*(rangeX[1]-rangeX[0]),0.53*ymax,Form("#Delta = %.3f",fSEP->GetParameter(1)-fped->GetParameter(1)));
	
	TDatime da;
	da.GetDate();
	TString da2 = da.AsSQLString();
	cout << da2 << endl;
	TLatex lm;
	lm.SetTextSize(0.02);
	lm.DrawLatex(rangeX[0]+0.68*(rangeX[1]-rangeX[0]),0.01*ymax,"MAPC @ " + da2);

	pad2png(c1,inFile(0,inFile.Last('.')) + ".png");

	return 0;
}
