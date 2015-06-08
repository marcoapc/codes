Int_t countDig(Float_t num);
Float_t twoDig(Float_t num);

// The order of the columns should be:
// x:y1:sy1:y2:sy2: ...
// To comment a line, insert a '#' as the FIRST character of that line
int csv2latex(TString filename = "filetest.csv") {

	ifstream in;
	in.open(filename.Data());
	if(in.fail()) {
		cout << "Error opening file: " << filename << endl;
		return 1;
	}

	ofstream out;
	TString outfile = filename(0,filename.Last('.'))+".tex";
	out.open(outfile.Data());
	if(out.fail()) {
		cout << "Error creating output file: " << outfile << endl;
		return 1;
	}

	// Date
	TDatime da;
	da.GetDate();
	out << "\% Created by csv2latex at " << da.AsSQLString() << endl;

	char line[500];
	stringstream ss;
	Int_t n=0;
	Float_t tmp1, tmp2;
	TString format;
	Int_t ndig;
	Bool_t stop;
	while(1) {
		in.getline(line,500);
		if(in.eof() || in.fail()) break;
		ss.clear();
		ss << line;
		if(line[0] != '#') {
			// Reading x
			ss.getline(line, 20, ',');
			if(ss.fail()) break;
			sscanf(line,"%f",&tmp1);
			out << Form("%.0f ",tmp1);
			// Reading all y:sy pairs
			stop=0;
			while(!stop) {
				ss.getline(line, 20, ',');
				if(ss.fail()) stop=1;
				sscanf(line,"%f",&tmp1);
				ss.getline(line, 20, ',');
				if(ss.fail()) stop=1;
				if(!stop) {
					sscanf(line,"%f",&tmp2);
					ndig = countDig(tmp2);
					format = Form("& %c.%df (%c.0f) ",'%',ndig,'%');
					out << Form(format.Data(),tmp1,twoDig(tmp2));
				}
			}
		}
		out << "\\\\" << endl;
	}

	in.close();
	out.close();

	return 0;
}

Int_t countDig(Float_t num) {
	Int_t n=0;
	while(num<10.0) {
		num*=10.0;
		n++;
	}
	return n;
}

Float_t twoDig(Float_t num) {
	while(num<10.0) num*=10.0;
	return num;
}
