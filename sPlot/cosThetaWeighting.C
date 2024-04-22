#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/Datasets/RooDataSetHelpers.h"
#include "../Tools/Datasets/SPlotHelpers.h"

#include "../Tools/FitShortcuts.h"
#include "../Tools/Style/Legends.h"

#include "../Tools/RooFitPDFs/InvariantMassModels.h"
#include "../Tools/Style/FitDistributions.h"

void cosThetaWeighting(Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS", const char* filename = "../Files/WeightedUpsilonSkimmedDataset.root") {
	TFile* f = TFile::Open(filename, "READ");
	if (!f) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
		return;
	}

	cout << "File " << filename << " opened" << endl;

	writeExtraText = true; // if extra text
	extraText = "      Internal";

	Int_t nCosThetaBins = 20;
	Float_t cosThetaMin = -1, cosThetaMax = 1;

	/// Set up the data
	using namespace RooFit;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	RooWorkspace wspace = SetUpWorkspace(filename, refFrameName);

	auto* data = InvMassCosThetaPhiDataset(wspace, ptMin, ptMax, refFrameName);

	// read variables in the reduced dataset in the workspace
	RooRealVar* invMass = wspace.var("mass");

	RooRealVar* cosThetaCS = wspace.var(CosThetaVarName(refFrameName));

	Long64_t nEntries = data->sumEntries();

	/// Invariant mass model

	// signal: one double-sided Crystal Ball PDF (symmetric Gaussian core) per Y resonance

	const char* signalShapeName = "SymDSCB";

	// background
	int order = 2;
	const char* bkgShapeName = Form("ChebychevOrder%d", order);

	auto* invMassModel = MassFitModel(wspace, signalShapeName, bkgShapeName, ptMin, ptMax, nEntries);

	auto* fitResult = invMassModel->fitTo(*data, Save(), Extended(kTRUE), PrintLevel(-1), NumCPU(NCPUs), Range(MassBinMin, MassBinMax), AsymptoticError(DoAsymptoticError), SumW2Error(!DoAsymptoticError));

	fitResult->Print("v");

	/// Draw the invariant mass distribution, to check the fit
	TCanvas* massCanvas = DrawMassFitDistributions(wspace, data, fitResult->floatParsFinal().getSize(), ptMin, ptMax);

	massCanvas->SaveAs(Form("1D/invMassFit_%s_cent%dto%d_pt%dto%dGeV.png", bkgShapeName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax), "RECREATE");

	/// SPlot time!

	SPlot sData = CreateSPlot(wspace, data, invMassModel);

	/// Draw the cos theta distribution with and without sWeights

	TCanvas* canvas = new TCanvas("canvas", "canvas", 650, 600);

	RooPlot* frame = cosThetaCS->frame(Title(" "), Bins(nCosThetaBins), Range(cosThetaMin, cosThetaMax));
	data->plotOn(frame, DrawOption("P0Z"), Name("data"), DataError(RooAbsData::SumW2));

	// create weighted data sets
	RooDataSet data_weight1S = GetSWeightedDataset(data, "1S");

	data_weight1S.plotOn(frame, DrawOption("P0Z"), MarkerColor(gColor1S), DataError(RooAbsData::SumW2), Name("data1S"));

	frame->GetYaxis()->SetMaxDigits(3);

	frame->Draw();

	gPad->RedrawAxis();

	frame->SetMaximum(nEntries / 30);

	TLatex text;
	text.SetTextAlign(22);
	text.SetTextSize(0.05);
	text.DrawLatexNDC(.55, .85, Form("%s, %s", CentralityRangeText(gCentralityBinMin, gCentralityBinMax), DimuonPtRangeText(ptMin, ptMax)));
	//text.DrawLatexNDC(.48, .8, Form("%.0f < m_{#mu#mu} < %.0f GeV/c^{2}", massMin, massMax));

	TLegend legend(.3, .8, .5, .65);
	legend.SetTextSize(.05);

	legend.AddEntry(frame->findObject("data"), "corrected dimuon events", "lp");
	legend.AddEntry(frame->findObject("data1S"), "with #varUpsilon(1S) yield sWeights", "lp");
	legend.DrawClone();

	gPad->Update();

	//CMS_lumi(canvas, gCMSLumiText);
	canvas->SaveAs(Form("1D/cosThetaCS_cent%dto%d_pt%dto%dGeV.png", gCentralityBinMin, gCentralityBinMax, ptMin, ptMax), "RECREATE");
}
