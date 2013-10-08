// Function to save Canvas in an PNG file
void pad2png(TCanvas *c, TString name) {
        TImage *img = TImage::Create();
        img->FromPad(c);
        img->WriteImage(name.Data());
        delete img;
        return;
}

// Function to configure a TGraph or TGraphErrors "gr" to be plotted as points
void MConfigPoints(TGraph *gr, Int_t color=1, Int_t size=1, Int_t style=8, Int_t errorSize=3) {
        gr->SetMarkerSize(size);
        gr->SetMarkerStyle(style);
        gr->SetMarkerColor(color);
        gStyle->SetEndErrorSize(errorSize);

        return;
}

// Function to configure a TGraph or TGraphErrors "gr" to be plotted as lines
void MConfigGraphLines(TGraph *f, Int_t color=1, Int_t width=2, Int_t style=1) {
	f->SetLineColor(color);
	f->SetLineWidth(width);
	f->SetLineStyle(style);

	return;
}

// Function to config a histogram TH1F
void MConfigHist(TH1F *h, TString title, TString xlabel, TString ylabel, Int_t color=4, Int_t fill=1) {
	h->SetTitle(title.Data());
	h->SetXTitle(xlabel.Data());
	h->SetYTitle(ylabel.Data());
	h->GetXaxis()->CenterTitle();
	h->GetYaxis()->CenterTitle();
	h->GetYaxis()->SetTitleOffset(1.3);
	if(fill) h->SetFillColor(color);
	else h->SetLineColor(color);
	
	return;
}

// Function to configure a TF1 "f" to be ploted with configured parameters
void MConfigLines(TF1 *f, Int_t color=1, Int_t width=2, Int_t style=1) {
	f->SetLineColor(color);
	f->SetLineWidth(width);
	f->SetLineStyle(style);

	return;
}

void MConfigCanvas(TCanvas *c1, Int_t islogy=0) {
	c1->SetGrid();
        c1->GetFrame()->SetBorderSize(10);
	if(islogy) c1->SetLogy();

	return;
}

void MConfigAxis(TGraph *gr, TString title, TString xlabel, TString ylabel) {
	if(!title.Length()) gr->SetTitle(title.Data());
	gr->GetXaxis()->SetTitle(xlabel.Data());
	gr->GetXaxis()->CenterTitle();
	gr->GetYaxis()->CenterTitle();
	gr->GetYaxis()->SetTitle(ylabel.Data());
	gr->GetYaxis()->SetTitleOffset(1.3);

	return;
}
