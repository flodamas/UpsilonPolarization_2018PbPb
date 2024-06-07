#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/Datasets/RooDataSetHelpers.h"
#include "SPlotHelpers.h"

#include "../Polarization/PolarFitHelpers.h"

#include "../Tools/Style/Legends.h"


// compare the sPlot and raw yield extraction methods
void compareRawCosThetaDistrib_customizedFits(Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS", const Int_t nCosThetaBins = 10, Int_t phiMin = -180, Int_t phiMax = 180, Int_t iState = gUpsilonState, const char* filename = "../Files/UpsilonSkimmedDataset.root") { //possible refFrame names: CS or HX
	
	writeExtraText = true; // if extra text
	extraText = "      Internal";

	using namespace RooFit;
	using namespace RooStats;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
	
	/// Set up the data
	RooWorkspace wspace = SetUpWorkspace(filename);

	// reduce the whole dataset (N dimensions) to (invariant mass, cos theta, phi)
	auto data = InvMassCosThetaPhiDataset(wspace, ptMin, ptMax);

	Long64_t nEntries = data.sumEntries();

	/// Bin edges and width 
	// Set the bin edges along the cosTheta axis depending on the number of bins 
	// (The bin edges are pre-defined in the function, so need to modify them if different bin edges are required)
	vector<Double_t> cosThetaBinEdges = setCosThetaBinEdges(nCosThetaBins);

	// bin width
	Double_t cosThetaStep = (cosThetaBinEdges[nCosThetaBins] - cosThetaBinEdges[0]) / nCosThetaBins;

	/// Set up the variables
	RooRealVar invMass = *wspace.var("mass");

	RooRealVar cosTheta = *wspace.var(CosThetaVarName(refFrameName));

	cosTheta.setRange(cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins]);
	cosTheta.setBins(nCosThetaBins);

	RooRealVar* yield1S = new RooRealVar("yield1S", "", 1000);
	RooRealVar* yield2S = new RooRealVar("yield2S", "", 100);
	RooRealVar* yield3S = new RooRealVar("yield3S", "", 10);

	/// Assign signal and background shape name to read the file for the yield extraction results
	const char* signalShapeName = "SymDSCB";

	// background shape array: ChebychevOrderN or ExpTimesErr	
	// const char* bkgShapeName[] = {
	//   // "ChebychevOrder1",
	//   "ChebychevOrder2",
	//   "ChebychevOrder2",
	//   "ChebychevOrder2",
	//   "ChebychevOrder2",
	//   "ChebychevOrder2",

	//   "ChebychevOrder2",
	//   "ChebychevOrder2",
	//   "ChebychevOrder2",
	//   "ChebychevOrder2",
	//   "ChebychevOrder2",

	//   "ChebychevOrder2",
	//   "ChebychevOrder2",
	//   "ChebychevOrder2",
	//   "ChebychevOrder2",
	//   "ChebychevOrder2",

	//   "ChebychevOrder2",
	//   // "ChebychevOrder1"
	// };

	const char* bkgShapeName[] = {
	  "ExpTimesErr",
	  "ExpTimesErr",
	  "ExpTimesErr",
	  "ExpTimesErr",
	  "ExpTimesErr",
	  "ExpTimesErr",
	  "ExpTimesErr",
	  "ExpTimesErr",
	  "ExpTimesErr",
	  "ExpTimesErr",
	  "ExpTimesErr",
	  "ExpTimesErr",
	  "ExpTimesErr",
	  "ExpTimesErr",
	  "ExpTimesErr",
	  "ExpTimesErr"
	};

	// a hisogram for raw yields before applying corrections
	TH1D* yieldHist = new TH1D("yieldHist", " ", nCosThetaBins, cosThetaBinEdges.data());

	auto invMassModel = MassFitModel(wspace, signalShapeName, bkgShapeName[0], ptMin, ptMax, nEntries);

	/// SPlot time!
	auto* sData = SWeightedDataset(wspace, ptMin, ptMax, signalShapeName, bkgShapeName[0]);

	// create weighted data sets
	RooDataSet data_weight1S = GetSpeciesSWeightedDataset(sData, "1S");

	Bool_t isCSframe = (strcmp(refFrameName, "CS") == 0) ? kTRUE : kFALSE;

	Float_t maxYield = 0;

	TCanvas* massCanvas = 0;

	/// apply weights and errors to each costheta bin
	for (Int_t iCosTheta = 0; iCosTheta < nCosThetaBins; iCosTheta++) {

		// get yields and their uncertainties
		const char* fitModelName = GetFitModelName(signalShapeName, ptMin, ptMax, isCSframe, cosThetaBinEdges[iCosTheta], cosThetaBinEdges[iCosTheta + 1], phiMin, phiMax);

		RooArgSet signalYields = GetSignalYields(yield1S, yield2S, yield3S, Form("RawData_%s", bkgShapeName[iCosTheta]), fitModelName);

		double yield1SVal = yield1S->getVal();

		double yield1SUnc = yield1S->getError();

		yieldHist->SetBinContent(iCosTheta + 1, yield1SVal);
		
		yieldHist->SetBinError(iCosTheta + 1, yield1SUnc);

		if ((yield1S->getVal()) > maxYield) maxYield = yield1S->getVal();
	}

	RooDataHist yieldRooDataHist("yieldRooDataHist", " ", cosTheta, Import(*yieldHist));

	/// compare raw yield and sPlot
	TCanvas* compCanvas = new TCanvas("compCanvas", "", 650, 600);

	RooPlot* frame = cosTheta.frame(Title(" "), Range(cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins]));

	yieldRooDataHist.plotOn(frame, Binning(nCosThetaBins), DrawOption("P0Z"), MarkerColor(kAzure + 2), Name("rawYield"));

	data_weight1S.plotOn(frame, Binning(nCosThetaBins), DrawOption("P0Z"), MarkerColor(gColor1S), Name("sData"));

	frame->GetYaxis()->SetMaxDigits(3);
	frame->SetMaximum(maxYield * 2.5);

	frame->Draw();

	gPad->RedrawAxis();

	TLegend legend(.22, .85, .45, .65);
	legend.SetTextSize(.05);
	legend.SetHeader(Form("%d < p_{T}^{#mu#mu} < %d GeV/c, %d#circ < #varphi_{%s} < %d#circ", ptMin, ptMax, phiMin, refFrameName, phiMax));
	legend.AddEntry(frame->findObject("rawYield"), "#varUpsilon(1S) raw yield", "lp");
	legend.AddEntry(frame->findObject("sData"), "sWeighted data", "lp");

	legend.DrawClone();

	gPad->Update();

	CMS_lumi(compCanvas, gCMSLumiText, 10);

	gSystem->mkdir("1D", kTRUE);
	compCanvas->SaveAs(Form("1D/compareRawCosTheta%s_cent%dto%d_pt%dto%dGeV_phi%dto%d_costheta%.1fto%.1f.png", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, phiMin, phiMax, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins]), "RECREATE");
}
