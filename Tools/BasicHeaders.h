#ifndef basic_headers
#define basic_headers

//basic headers
//c headers
#include <iostream>
#include <fstream>
#include <cmath>

//root headers
#include "TAxis.h"
#include "TCanvas.h"
#include "TClonesArray.h"
#include "TEfficiency.h"
#include "TF1.h"
#include "TFile.h"
#include "TH2.h"
#include "TH3.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLorentzVector.h"
#include "TNtuple.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TRandom3.h"
#include "TROOT.h"
#include "TRotation.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "TVector3.h"

//roofit headers
#include "RooAbsCategory.h"
#include "RooAbsReal.h"
#include "RooAbsRealLValue.h"
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooArgSet.h"
#include "RooCategoryProxy.h"
#include "RooChebychev.h"
#include "RooConstVar.h"
#include "RooCrystalBall.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooHist.h"
#include "RooHypatia2.h"
#include "RooPlot.h"
#include "RooRealProxy.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"

TFile* openFile(const char* fileName){

	TFile* file = TFile::Open(fileName, "READ");
	if (!file) {
		cout << fileName << " not found. Check the directory of the file." << endl;
		exit(1);
	}

	return file;
}

#endif