#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/FitShortcuts.h"
#include "../Tools/Style/Legends.h"

#include "../Tools/Style/FitDistributions.h"

#include "../Tools/RooFitPDFs/CosThetaPolarizationPDF.h"

// fit the efficiency as a function of cos theta from the reco MC dataset

// code based on the RooFit tutorials https://root.cern/doc/master/rf701__efficiencyfit_8C.html and https://root.cern/doc/master/rf703__effpdfprod_8C.html

// reduce the whole dataset (N dimensions) to (cos theta, phi) and the selection category, allowing cutting off along phi
RooDataSet* ReducedMCDataset(RooDataSet* allDataset, RooWorkspace& wspace, Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS", Int_t phiMin = -180, Int_t phiMax = 180) {
	if (allDataset == nullptr) {
		cerr << "Null RooDataSet provided to the reducer method!!" << endl;
		return nullptr;
	}

	const char* kinematicCut = Form("(centrality >= %d && centrality < %d) && (rapidity > %f && rapidity < %f) && (pt > %d && pt < %d) && (phi%s > %d && phi%s < %d)", 2 * gCentralityBinMin, 2 * gCentralityBinMax, gRapidityMin, gRapidityMax, ptMin, ptMax, refFrameName, phiMin, refFrameName, phiMax);

	RooDataSet* reducedDataset = (RooDataSet*)allDataset->reduce(RooArgSet(*(wspace.cat("recoCategory")), *(wspace.var(CosThetaVarName(refFrameName))), *(wspace.var(PhiVarName(refFrameName)))), kinematicCut);

	wspace.import(*reducedDataset, RooFit::Rename(Form("(cos theta, phi) %s reduced dataset", refFrameName)));

	return reducedDataset;
}

void efficiencyFit_1D(Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS", Int_t phiMin = -180, Int_t phiMax = 180, Int_t iState = 1) {
	writeExtraText = true; // if extra text
	extraText = "      Internal";

	Int_t nCosThetaBins = 50;
	Float_t cosThetaMin = -1, cosThetaMax = 1;

	// set up the data
	using namespace RooFit;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	const char* filename = Form("../Files/Y%dSReconstructedMCWeightedDataset.root", iState);

	TFile* f = TFile::Open(filename, "READ");
	if (!f) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
		return;
	}

	cout << "File " << filename << " opened" << endl;

	RooDataSet* allDataset = (RooDataSet*)f->Get("MCdataset");

	// import the dataset to a workspace
	RooWorkspace wspace("workspace");
	wspace.import(*allDataset);

	auto* data = ReducedMCDataset(allDataset, wspace, ptMin, ptMax, refFrameName, phiMin, phiMax);

	RooRealVar normFactor("normFactor", "normFactor", 0, 1);
	RooRealVar lambdaTheta("lambdaTheta", "lambdaTheta", -2.0, 2.0);

	RooRealVar cosTheta = *wspace.var(CosThetaVarName(refFrameName));

	RooCategory recoCat = *wspace.cat("recoCategory");

	// define efficiency objects

	RooFormulaVar effFunc("effFunc", "@0*(1+@1*@2*@2)", RooArgList(normFactor, lambdaTheta, cosTheta));

	RooEfficiency effPDF("effPDF", "effPDF", effFunc, recoCat, "selected");

	//auto cosThetaPDF = CosThetaPolarizationPDF("cosThetaPDF", " ", cosTheta, lambdaTheta);
	// Fit conditional efficiency pdf to data

	auto* fitResult = effPDF.fitTo(*data, ConditionalObservables(cosTheta), Save(), PrintLevel(+1), Range(cosThetaMin, cosThetaMax), AsymptoticError(DoAsymptoticError));

	fitResult->Print("v");

	// draw

	TCanvas* canvas = new TCanvas("canvas", "canvas", 650, 600);

	RooPlot* frame = cosTheta.frame(Title(" "), Bins(nCosThetaBins), Range(cosThetaMin, cosThetaMax));
	frame->SetYTitle("Efficiency");

	data->plotOn(frame, Efficiency(recoCat), DrawOption("P0Z"), Name("data"));

	effPDF.plotOn(frame, LineColor(kRed + 1), Name("polaResult"));

	//frame->GetYaxis()->SetRangeUser(0, 1000);
	frame->SetMaximum(1);

	frame->Draw();

	// cosmetics

	TLegend legend(.22, .88, .5, .65);
	legend.SetTextSize(.05);
	legend.SetHeader(Form("%s, %s", CentralityRangeText(gCentralityBinMin, gCentralityBinMax), DimuonPtRangeText(ptMin, ptMax)));
	legend.AddEntry(frame->findObject("data"), Form("selected #varUpsilon(%dS) MC candidates", iState), "EP");
	legend.AddEntry(frame->findObject("polaResult"), Form("distribution fit: #lambda_{#theta} = %.3f #pm %.3f", lambdaTheta.getVal(), lambdaTheta.getError()), "l");

	legend.DrawClone();

	gPad->Update();

	CMS_lumi(canvas, Form("#varUpsilon(%dS) Hydjet-embedded MC", iState));

	canvas->SaveAs(Form("EfficiencyMaps/%dS/Efficiency1D_CosTheta%s_cent%dto%d_pt%dto%dGeV_phi%dto%d.png", iState, refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, phiMin, phiMax), "RECREATE");
}
