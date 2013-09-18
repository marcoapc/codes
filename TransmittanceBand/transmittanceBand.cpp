// MARCO FUNCTIONS TO INCLUDE
#include "../MLibraries/MallIncludes.h"
#include "../MLibraries/MConfigGraphs.h"
#include "../MLibraries/MMultiFits.h"
#include "../MLibraries/McommentNlines.h"

//#include "../MLibraries/TestList.h"

void Instructions(void);

int transmittanceBand(TString maindir) {
	Int_t i,j,k;
	Float_t count;
	Float_t wl_range[2] = {200.0, 900.0};
	Float_t wl_steps = 10.0;
	
	// For calculating stdev and mean
	Int_t wl_size = 1+((int)(wl_range[1]-wl_range[0])/(wl_steps));
	cout << "Number of walelengths: " << wl_size << endl << endl;
	Double_t *w = (Double_t *)malloc(wl_size*sizeof(Double_t));
	Double_t *a2 = (Double_t *)malloc(wl_size*sizeof(Double_t));
	Double_t *a = (Double_t *)malloc(wl_size*sizeof(Double_t));
	Double_t *Tstd = (Double_t *)malloc(wl_size*sizeof(Double_t));
	Double_t *Tmean = (Double_t *)malloc(wl_size*sizeof(Double_t));
	Double_t *Tbandmin = (Double_t *)malloc(wl_size*sizeof(Double_t));
	Double_t *Tbandmax = (Double_t *)malloc(wl_size*sizeof(Double_t));
	Double_t *Tbandshape = (Double_t *)malloc(wl_size*sizeof(Double_t));
	for(k=0;k<wl_size;k++) {
		a2[k]=0.0;
		a[k]=0.0;
	}

	// Creating Tree with all the files
	TTree *tr = new TTree("tr","Transmittance");
	Float_t wvl;
	TBranch *wl = tr->Branch("wl",&wvl,"wl/F");
	for(wvl=wl_range[0],i=0; wvl<=wl_range[1]; wvl+=wl_steps) {
		wl->Fill();
		w[i++]=wvl;
	}
	tr->SetEntries(wl->GetEntries());

	// Finding csv files in maindir
	Int_t ind=0;
	Float_t T;
	Float_t *readRow;
	const int MaxNumFiles = 500;
	TBranch *br[MaxNumFiles];
	if(!maindir.EndsWith("/")) maindir += '/';
	TSystemDirectory dir(maindir.Data(), maindir.Data());
	TList *tempfiles = dir.GetListOfFiles();
	if(tempfiles) {
		TSystemFile *f;
		TString fname;
		TIter next(tempfiles);
		while((f=(TSystemFile*)next())) {
			fname = f->GetName();
			if(!f->IsDirectory() && fname.EndsWith(".csv")) {
				//DO SOMETHING HERE!
	cout << fname.Data() << endl;
	br[ind] = tr->Branch(Form("f%d",ind),&T,Form("T%d/F",ind));
	TNtuple *Tdata = new TNtuple("Tdata","Tdata","w:Td");
	Tdata->ReadFile(commentNlines(maindir+fname));
	cout << "Entries: " << Tdata->GetEntries() << endl;
	for(i=Tdata->GetEntries(),wvl=wl_range[0],k=0; i>=0 && wvl<=wl_range[1];) {
		Tdata->GetEntry(i);
		readRow = Tdata->GetArgs();
		if(readRow[0]==wvl) {
			T = readRow[1];
			br[ind]->Fill();
			a2[k] += T*T;
			a[k] += T;
			//cout << "Fill: wvl=" << wvl << "\tT=" << T << endl;
			wvl+=wl_steps;
			i--;
			k++;
		}
		else if(readRow[0]>wvl) {
			T = 0.0;
			br[ind]->Fill();
			a2[k] += T*T;
			a[k] += T;
			//cout << "Fill: wvl=" << wvl << "\tT=" << T << endl;
			wvl+=wl_steps;
			k++;
		}
		else { i--; }
	}
	ind++;
	Tdata->Delete();
	//break;
				//END DOING SOMETHING
			}
		}
	}

	// Plot all used data
	Float_t linewidth = 1.0;
	TCanvas *cf = new TCanvas("cf","Transmittance Blast",1200,600);
	MConfigCanvas(cf);
	TGraph *gr[MaxNumFiles];
	for(i=0;i<ind;i++) {
		tr->Draw(Form("wl:f%d",i),"","goff");
		gr[i] = new TGraph(tr->GetSelectedRows(),tr->GetV1(),tr->GetV2());
		gr[i]->SetLineColor(11);
		gr[i]->SetLineWidth(linewidth);
		if(i!=0) gr[i]->Draw("lsame");
		else {
			gr[i]->SetTitle("All Data Used");
			gr[i]->GetXaxis()->SetTitle("Wavelength / nm");
			gr[i]->GetYaxis()->SetTitle("Transmittance / %");
			gr[i]->Draw("al");
		}
	}
	pad2png(cf,"allData.png");
	
	// Calculating mean and stdev for band construction
	//TLeaf *lf = wl->GetLeaf("wl");
	for(k=0;k<wl_size;k++) {
		Tmean[k] = a[k]/ind;
		Tstd[k] = sqrt((a2[k]-2.0*a[k]*Tmean[k]+ind*Tmean[k]*Tmean[k])/(ind-1.0));
		
		Tbandmin[k] = Tmean[k]-Tstd[k];
		Tbandmax[k] = Tmean[k]+Tstd[k];
		// To print on screen Tmean and Tstd (not working)
		//tr->GetEntry(k);
		//cout << Form("wl=%.0f\tMean=%.2f\tStd=%.2f",wl->GetValue(),Tmean[k],Tstd[k]) << endl;
	}

	TCanvas *c2 = new TCanvas("c2","Band",600,600);
	TGraph *grmin = new TGraph(wl_size,w,Tbandmin);
	TGraph *grmax = new TGraph(wl_size,w,Tbandmax);
	TGraph *grshade = new TGraph(2*wl_size);
	for(k=0;k<wl_size;k++) {
		cout << "wl=" << w[k] << "\tTmin=" << Tbandmin[k] << "\tTmax=" << Tbandmax[k] << endl;
		grshade->SetPoint(k,w[k],Tbandmax[k]);
		grshade->SetPoint(wl_size+k,w[wl_size-k-1],Tbandmin[wl_size-k-1]);
	}
	grshade->SetFillStyle(3013);
	grshade->SetFillColor(16);
	grmin->SetTitle("");
	grmin->GetXaxis()->SetTitle("Wavelength / nm");
	grmin->GetXaxis()->SetLimits(200.0,900.0);
	grmin->GetYaxis()->SetTitle("Transmittance / %");
	grmin->GetYaxis()->SetTitleOffset(1.3);
	grmin->GetYaxis()->SetRangeUser(0.0,100.0);
	grmin->Draw("al");
	grmax->Draw("lsame");
	grshade->Draw("fsame");

	TLegend *leg = new TLegend(0.70,0.14,0.87,0.2);
	leg->SetTextFont(50);
	leg->SetTextSize(0.025);
	leg->AddEntry(grshade,"SP30","f");
	leg->Draw();

	pad2png(c2,"band.png");
	
	
/*	TNtuple *Tdata1 = new TNtuple("Tdata1","Data1","w:T");
	Tdata1->ReadFile(commentNlines("Sample297.Sample.Raw.csv"));
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
*/
	//pad2png(cf,"TransmittanceBlast.png");

	return 0;
}

int main(int argc, char *argv[]) {
	if(argc!=2) {
		Instructions();
		return 0;
	}
        return transmittanceBand(argv[1]);
}

void Instructions(void) {
	cout << endl << endl << "Usage:" << endl
	     << "    ./transmittanceBand <mainFolder>" << endl << endl
	     << "where" << endl << "  <mainFodler> is the folder with csv files." << endl
             << endl << "by Marco Antonio Pannunzio Carmignotto" << "   09/11/2013" << endl << endl;
	return;
}
