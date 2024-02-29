#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/FitShortcuts.h"
#include "../Tools/Style/Legends.h"

#include "../Tools/RooFitPDFs/InvariantMassModels.h"
#include "../Tools/Style/FitDistributions.h"

#include "../Tools/RooFitPDFs/CosThetaPolarizationPDF.h"
#include "../Tools/RooFitPDFs/cosThetaPolarFunc.h"

void correctedYield_1D_customizedFits(Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS"/*, const Int_t nCosThetaBins = 10, Float_t cosThetaMin = -1, Float_t cosThetaMax = 1.*/, Int_t phiMin = 0, Int_t phiMax = 180) {
	writeExtraText = true; // if extra text
	extraText = "      Internal";

	using namespace RooFit;
	using namespace RooStats;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	// /// Set up the variables
	// RooRealVar cosTheta("cosTheta", "", -1, 1);

	RooRealVar* yield1S = new RooRealVar("yield1S", "", 1000);
	RooRealVar* yield2S = new RooRealVar("yield2S", "", 100);
	RooRealVar* yield3S = new RooRealVar("yield3S", "", 10);

	/// Bin width
	const Int_t nCosThetaBins = 7;
	Float_t cosThetaBinEdges[nCosThetaBins+1] = {-0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7}; 
	// Float_t cosThetaBinEdges[nCosThetaBins+1] = {-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8}; 

	/// Set up the variables
	RooRealVar cosTheta("cosTheta", "", cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins]);
	cosTheta.setBins(nCosThetaBins);

	/// Assign signal and background shape name to read the file for the yield extraction results
	const char* signalShapeName = "SymDSCB";

	// background shape array: ChebychevOrderN or ExpTimesErr
	const char* bkgShapeName[] = {
		"ChebychevOrder2", 
		"ChebychevOrder2",
		"ChebychevOrder2",
		"ChebychevOrder2",
		"ChebychevOrder2",
		"ChebychevOrder2",
		"ChebychevOrder2",
		"ChebychevOrder2",
		"ChebychevOrder2",
		// "ChebychevOrder2"
	};

	/// "Standard" procedure: extract the yields per bin
	cout << "?" << endl;
	TH1D* standardCorrectedHist = new TH1D("standardCorrectedHist", " ", nCosThetaBins, cosThetaBinEdges);
	// TH1D* standardCorrectedHist = new TH1D("standardCorrectedHist", " ", nCosThetaBins, cosThetaMin, cosThetaMax);

	// Float_t cosThetaStep = ((cosThetaMax - cosThetaMin) / nCosThetaBins);

	TCanvas* massCanvas = 0;
	Bool_t isCSframe = (strcmp(refFrameName, "CS") == 0) ? kTRUE : kFALSE;

	Float_t maxYield = 0;

	for (Int_t iCosTheta = 0; iCosTheta < nCosThetaBins; iCosTheta++) {
		// Float_t cosThetaVal = cosThetaBinEdges[iCosTheta];

		const char* fitModelName = GetFitModelName(signalShapeName, ptMin, ptMax, isCSframe, cosThetaBinEdges[iCosTheta], cosThetaBinEdges[iCosTheta+1], phiMin, phiMax);

		RooArgSet signalYields = GetSignalYields(yield1S, yield2S, yield3S, bkgShapeName[iCosTheta], fitModelName);

		standardCorrectedHist->SetBinContent(iCosTheta + 1, yield1S->getVal());
		standardCorrectedHist->SetBinError(iCosTheta + 1, yield1S->getError());

		if (yield1S->getVal() > maxYield) maxYield = yield1S->getVal();
	}
	cout << yield1S->getError() << endl;

	RooDataHist correctedHist("correctedHist", " ", cosTheta, standardCorrectedHist);

	/// Draw the cos theta distributions

	TCanvas* canvas = new TCanvas("canvas", "canvas", 650, 600);

	RooPlot* frame = cosTheta.frame(Title(" "), Range(cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins]));
	frame->SetXTitle(Form("cos #theta_{%s}", refFrameName));

	correctedHist.plotOn(frame, DrawOption("P0Z"), MarkerColor(kAzure + 2), Name("dataPoints"));

	frame->GetYaxis()->SetRangeUser(0, 2 * maxYield);
	frame->GetYaxis()->SetMaxDigits(3);

	frame->Draw();

	gPad->RedrawAxis();

	frame->SetMaximum(2 * maxYield);

	/// Polarization fit
	// with RooFit
	cout << endl
	     << "Distribution fit for polarization paramaters extraction" << endl;

	RooRealVar lambdaTheta("lambdaTheta", "lambdaTheta", -1, 1);
	auto cosThetaPDF_1S = CosThetaPolarizationPDF("cosThetaPDF_1S", " ", cosTheta, lambdaTheta);

   	enableBinIntegrator(cosThetaPDF_1S, cosTheta.numBins());

	auto* polarizationFitResult = cosThetaPDF_1S.fitTo(correctedHist, Save(), Extended(kTRUE)/*, PrintLevel(+1)*/, NumCPU(NCPUs), Range(cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins]), SumW2Error(kFALSE)/*, AsymptoticError(DoAsymptoticError)*/);

	polarizationFitResult->Print("v");

	cosThetaPDF_1S.plotOn(frame, LineColor(kRed + 1), Name("polaResult"));

	frame->Draw();

	// cosmetics

	TLegend legend(.22, .88, .5, .68);
	legend.SetTextSize(.05);
	legend.SetHeader(Form("centrality %d-%d%%, %d < p_{T}^{#mu#mu} < %d GeV/c", gCentralityBinMin, gCentralityBinMax, ptMin, ptMax));
	legend.AddEntry(frame->findObject("dataPoints"), "#varUpsilon(1S) corrected yield", "lp");
	legend.AddEntry(frame->findObject("polaResult"), Form("distribution fit: #lambda_{#theta} = %.2f #pm %.2f", lambdaTheta.getVal(), lambdaTheta.getError()), "l");

	legend.DrawClone();

	gPad->Update();

	//CMS_lumi(canvas, gCMSLumiText);

	gSystem->mkdir("DistributionFits/1D", kTRUE);
	canvas->SaveAs(Form("DistributionFits/1D/compareCorrectedCosTheta%s_cent%dto%d_pt%dto%dGeV_phi%dto%d.png", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, phiMin, phiMax), "RECREATE");
}