#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/FitShortcuts.h"
#include "../Tools/Style/Legends.h"
#include "../Tools/Graphs.h"

#include "../Tools/RooFitPDFs/InvariantMassModels.h"
#include "../Tools/RooFitPDFs/ErrorFuncTimesExp.h"
#include "../Tools/Style/FitDistributions.h"

void invariantMass(Int_t ptMin = 0, Int_t ptMax = gPtMax, const char* refFrameName = "HX", Float_t cosThetaMin = -1, Float_t cosThetaMax = 1, Int_t phiMin = -180, Int_t phiMax = 180, Float_t massMin = MassBinMin, Float_t massMax = MassBinMax) {
	/// Set up the data
	using namespace RooFit;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	const char* filename = Form("../Files/UpsilonSkimmedDataset%s.root", gMuonAccName);

	RooWorkspace wspace = SetUpWorkspace(filename);

	auto allDataset = (RooDataSet*)wspace.data("dataset");

	RooDataSet data = ReducedMassDataset(allDataset, wspace, ptMin, ptMax, refFrameName, massMin, massMax, cosThetaMin, cosThetaMax, phiMin, phiMax);

	RooRealVar mass = *wspace.var("mass");
	mass.setRange("MassFitRange", massMin, massMax);

	Long64_t nEntries = data.sumEntries();

	/// Invariant mass model

	const char* signalShapeName = "SymDSCB";
	const char* bkgShapeName = "ExpTimesErr";

	const char* fitModelName = GetFitModelName(signalShapeName, ptMin, ptMax, refFrameName, cosThetaMin, cosThetaMax, phiMin, phiMax);

	BuildInvariantMassModel(wspace, signalShapeName, bkgShapeName, fitModelName, nEntries, false);

	auto* fitResult = RawInvariantMassFit(wspace, data);

	auto maxYield = wspace.var("yield1S")->getVal();

	TCanvas* canvas = MyCanvas("canvas", 650, 600);

	RooPlot* frame = InvariantMassRooPlot(wspace, data);
	frame->GetXaxis()->SetTitle("#it{m}_{#mu^{#plus}#mu^{#font[122]{\55}}} (GeV/#it{c}^{ 2})");
	frame->Draw();

	frame->SetMaximum(frame->GetMaximum() * 1.1);

	// legend and text

	auto legend = LegendBlock(0.2, 0.8, 6);
	legend->AddEntry("data", "Data", "ep");
	legend->AddEntry("total", "Total fit", "l");
	legend->AddEntry("1S", "#varUpsilon(1S)", "l");
	legend->AddEntry("2S", "#varUpsilon(2S)", "l");
	legend->AddEntry("3S", "#varUpsilon(3S)", "l");
	legend->AddEntry("background", "Background", "l");

	legend->Draw();

	TPaveText* kinematics = new TPaveText(0.6, 0.9, 0.95, 0.5, "NDCNB");
	// TPaveText* text = new TPaveText(0.65, 0.90, 0.95, 0.60, "NDCNB");
	kinematics->SetTextSize(0.045);
	kinematics->SetFillColor(4000);
	kinematics->SetBorderSize(0);
	kinematics->AddText(CentralityRangeText(gCentralityBinMin, gCentralityBinMax));
	kinematics->AddText(DimuonRapidityRangeText(gRapidityMin, gRapidityMax));
	kinematics->AddText(Form("%d < %s < %d %s", ptMin, gDimuonPtVarTitle, ptMax, gPtUnit));
	kinematics->AddText(" ");
	kinematics->AddText("Helicity frame");
	kinematics->AddText(Form("%.2f < cos #theta < %.2f", cosThetaMin, cosThetaMax));
	kinematics->AddText(Form("%d < |%s| < %d %s", phiMin, gPhiSymbol, phiMax, gPhiUnit));

	kinematics->SetAllWith("", "align", 32);
	kinematics->Draw();

	canvas->Modified();
	canvas->Update();

	CMS_lumi(canvas, gCMSLumiText, 11);

	const char* totalFitModelName = GetTotalFitModelName(bkgShapeName, signalShapeName, ptMin, ptMax, refFrameName, cosThetaMin, cosThetaMax, phiMin, phiMax);

	SaveMyCanvas(canvas, totalFitModelName);
}
