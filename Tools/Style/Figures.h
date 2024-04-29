#include "TH1.h"
#include "TH2.h"

#include "RooDataSet.h"
#include "RooRealVar.h"

TH2* TH2fromRooDataSet(RooDataSet data, const char* histoName, RooRealVar var1, Int_t nVar1Bins, Float_t var1Min, Float_t var1Max, RooRealVar var2, Int_t nVar2Bins, Float_t var2Min, Float_t var2Max) {
	TH2* histo = dynamic_cast<TH2*>(data.createHistogram(histoName, var1, RooFit::Binning(nVar1Bins, var1Min, var1Max), RooFit::YVar(var2, RooFit::Binning(nVar2Bins, var2Min, var2Max))));

	histo->SetTitle(" ");
	histo->GetXaxis()->CenterTitle();

	histo->GetYaxis()->CenterTitle();

	histo->GetZaxis()->SetMaxDigits(3);
	return histo;
}
