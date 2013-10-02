#include<cstdio>
#include<stdlib.h>
#include<Riostream.h>

int ConvertADCFiles(TString inFile="", TString outFile="", Int_t ch0=-2, Int_t ch1=-2, Int_t ch2=-2, Int_t ch3=-2, Int_t ch4=-2, Int_t ch5=-2, Int_t ch6=-2, Int_t ch7=-2, Int_t ch8=-2, Int_t ch9=-2, Int_t ch10=-2, Int_t ch11=-2) {

	// Variables
	Int_t ch[12];
	Int_t Nch;
	Int_t i;
	int l[12];
	char line[100];
	stringstream ss, s2;

	// Printing instructions
	if(inFile.IsNull() || outFile.IsNull() || ch0==-2) {
		Instructions();
		return 0;
	}

	// Constructing channels to be saved vector
	//ch[0]=0; Nch=1;
	if(ch0==-1) {for(i=0; i<12; i++) ch[i]=i; Nch=12;}
	else {
		if(ch0>=0) {ch[0]=ch0; Nch=1;}
		if(ch1>=0) {ch[1]=ch1; Nch=2;}
		if(ch2>=0) {ch[2]=ch2; Nch=3;}
		if(ch3>=0) {ch[3]=ch3; Nch=4;}
		if(ch4>=0) {ch[4]=ch4; Nch=5;}
		if(ch5>=0) {ch[5]=ch5; Nch=6;}
		if(ch6>=0) {ch[6]=ch6; Nch=7;}
		if(ch7>=0) {ch[7]=ch7; Nch=8;}
		if(ch8>=0) {ch[8]=ch8; Nch=9;}
		if(ch9>=0) {ch[9]=ch9; Nch=10;}
		if(ch10>=0) {ch[10]=ch10; Nch=11;}
		if(ch11>=0) {ch[11]=ch11; Nch=12;}
	}

	// Getting current directory
	TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
	dir.ReplaceAll("ConvertADCFiles.C","");
	dir.ReplaceAll("/./","/");
	printf("Current Directory: %s\n",dir.Data());

	// Opening input file
	ifstream in;
	in.open(Form("%s%s",dir.Data(),inFile.Data()));
	if(in.fail()) {
		printf(" Marco: Error opening input file.\n");
		printf("        %s%s\n",dir.Data(),inFile.Data());
		return 1;
	}
	
	// Opening output file
	ofstream out;
	out.open(Form("%s%s",dir.Data(),outFile.Data()));
	if(out.fail()) {
		printf(" Marco: Error opening output file.\n");
		printf("        %s%s\n",dir.Data(),outFile.Data());
		return 1;
	}

	while(1) {
		in.getline(line,100);
		if(in.eof()) break;
		if(in.fail()) {
			printf(" Marco: Fail reading input file.\n");
			return 1;
		}
		ss << line;
		for(i=0;i<3;i++) ss.getline(line, 20, ',');
		for(i=0;i<12;i++) {
			ss.getline(line, 20, ',');
			sscanf(line,"%x",&l[i]);
		}
		
		for(i=0;i<Nch;i++) out << l[ch[i]] << " ";
		out << endl; 
	}
	
	in.close();
	out.close();
	return 0;
}
void Instructions(void) {
        printf("\nInstructions:\n");
        printf("  .\\ConvertADCFiles <input> <output> <ch0> <ch1> ...\n");
        printf("      - <input> - input file name.\n");
        printf("      - <output> - output file name.\n");
        printf("      - <ch0>    - First channel to be saved.\n");
        printf("      - <ch1>    - Second channel to be saved.\n");
        printf("      - ...      - Sequence of numbers of all the channels.\n");
        printf("      * Channels number in interval [0, 12].\n");
	printf("      * Special usage: if <ch1> is -1, all channels will be saved.\n\n");
        printf(" by Marco Antonio Pannunzio Carmignotto - 04/20/2012\n\n");

        return;
}
