// comment the first N lines of a file, if it is not commented yet, with the comment char comChar
TString commentNlines(TString filename, Int_t nlines=1, TString comChar="#") {
	Int_t i;
	ifstream arq;
	ofstream arqout;
	const int maxLineSize = 500;
	char line[maxLineSize];
	TString fileout = ".temp_" + filename(filename.Last('/')+1,filename.Length()-filename.Last('/'));
	
	// Open files
	arq.open(filename.Data());
	arqout.open(fileout.Data());
	if(!arq) cout << "Error opening input file: " << filename.Data() << endl;
	if(!arqout) cout << "Error opening output file: " << fileout.Data() << endl;
	if(!arq || !arqout) {
		arq.close();
		arqout.close();
		TString a = "";
		return a;
	}

	// Check if the first N character of the file is equal to comChar
	TString l;
	i=0;
	while(arq.getline(line,maxLineSize)) {
		l=line;
		if((i++)<nlines && l.CompareTo(comChar)!=1)  // there is not a comChar in the 1st pos
			arqout << "# " + l << endl;
		else
			arqout << l << endl;
	}

	//arqtemp = fopen("Mttt.txt","w");
	arqout.close();
	arq.close();

	return fileout;
}
