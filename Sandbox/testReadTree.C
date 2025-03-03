#include "../Tools/BasicHeaders.h"

void testReadTree() {


	TFile* file = TFile::Open("DatasetToTree.root", "READ");

	TFile* infile = (TFile*) file->Get("DatasetToTree.root"); 

	// infile->cd();
	
	infile->Print();

}