#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/FitShortcuts.h"
#include "../Tools/Style/Legends.h"

#include "../Tools/RooFitPDFs/InvariantMassModels.h"
#include "../Tools/Style/FitDistributions.h"

#include "../Tools/RooFitPDFs/CosThetaPolarizationPDF.h"
#include "../Tools/RooFitPDFs/cosThetaPolarFunc.h"

#include "../ReferenceFrameTransform/Transformations.h"

TEfficiency* rebinTEff3DMap(TEfficiency* TEff3DMap, Int_t phiMin = -180, Int_t phiMax = 180, Int_t ptMin = 0, Int_t ptMax = 30, Int_t nCosThetaBins = 10, Double_t* cosThetaBinEdges = nullptr) {
	// rebin efficiency maps based on costheta, phi, and pT selection 
	TH3D* hPassed = (TH3D*) TEff3DMap->GetPassedHistogram();
	TH3D* hTotal = (TH3D*) TEff3DMap->GetTotalHistogram();

	Int_t iPhiMin = hPassed->GetYaxis()->FindBin(phiMin);
	Int_t iPhiMax = hPassed->GetYaxis()->FindBin(phiMax);

	Int_t iPtMin = hPassed->GetZaxis()->FindBin(ptMin);
	Int_t iPtMax = hPassed->GetZaxis()->FindBin(ptMax);

	TH1D* hPassedCosTheta = (TH1D*) hPassed->ProjectionX("hPassedCosTheta", iPhiMin, iPhiMax-1, iPtMin, iPtMax-1, "eo");
	TH1D* hTotalCosTheta = (TH1D*) hTotal->ProjectionX("hTotalCosTheta", iPhiMin, iPhiMax-1, iPtMin, iPtMax-1, "eo");

	TH1D* hPassedCosTheta_Rebin = (TH1D*) hPassedCosTheta->Rebin(nCosThetaBins, "hPassedCosTheta_Rebin", cosThetaBinEdges);
	TH1D* hTotalCosTheta_Rebin = (TH1D*) hTotalCosTheta->Rebin(nCosThetaBins, "hTotalCosTheta_Rebin", cosThetaBinEdges);

	TEfficiency* TEffMapCosTheta = new TEfficiency("TEffMapCosTheta", "cos #theta_{CS}; efficiency", nCosThetaBins, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins]);
	
	TEffMapCosTheta->SetPassedHistogram(*hPassedCosTheta_Rebin, "f");
	TEffMapCosTheta->SetTotalHistogram(*hTotalCosTheta_Rebin, "f"); 

	return TEffMapCosTheta;
}

TH1D* rebin3DUnc(TH3D* systEff, Int_t phiMin = -180, Int_t phiMax = 180, Int_t ptMin = 0, Int_t ptMax = 30, Int_t nCosThetaBins = 10, Double_t* cosThetaBinEdges = nullptr){
	// rebin efficiency maps based on costheta, phi, and pT selection 
	// uncertainty sum is sqrt(pow(unc1, 2) + pow(unc2, 2)), so fold it manually
	TH1D* h1DUnc = new TH1D("h1DUnc", "cos #theta_{CS}; uncertainty", nCosThetaBins, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins]);

	Int_t iPhiMin = systEff->GetYaxis()->FindBin(phiMin);
	Int_t iPhiMax = systEff->GetYaxis()->FindBin(phiMax);

	Int_t iPtMin = systEff->GetZaxis()->FindBin(ptMin);
	Int_t iPtMax = systEff->GetZaxis()->FindBin(ptMax);

	for (int iCosThetaBin=1; iCosThetaBin<=nCosThetaBins; iCosThetaBin++) {

		Double_t phiSumSystEff = 0;

		// sum uncertainties along the phi axis
		for (int iPhiBin=iPhiMin; iPhiBin<=iPhiMax; iPhiBin++) {

			Double_t ptSumSystEff = 0;

			// sum uncertainties along the pt axis
			for (int iPtBin=iPtMin; iPtBin<=iPtMax; iPtBin++) {

				ptSumSystEff = TMath::Hypot(ptSumSystEff, systEff->GetBinContent(iCosThetaBin, iPhiBin, iPtBin));
			}

			phiSumSystEff = TMath::Hypot(phiSumSystEff, ptSumSystEff);
		}

		h1DUnc->SetBinContent(iCosThetaBin, phiSumSystEff);
	}	

	return h1DUnc;
}

void rawYield_1D_customizedFits(Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS" /*, const Int_t nCosThetaBins = 10, Float_t cosThetaMin = -1, Float_t cosThetaMax = 1.*/, Int_t phiMin = -180, Int_t phiMax = 180) {
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
	const Int_t nCosThetaBins = 8;
	// Double_t cosThetaBinEdges[nCosThetaBins+1] = {-0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7}; 
	Double_t cosThetaBinEdges[nCosThetaBins + 1] = {-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8}; 

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
	TH1D* standardCorrectedHist = new TH1D("standardCorrectedHist", " ", nCosThetaBins, cosThetaBinEdges);
	// TH1D* standardCorrectedHist = new TH1D("standardCorrectedHist", " ", nCosThetaBins, cosThetaMin, cosThetaMax);

	// Float_t cosThetaStep = ((cosThetaMax - cosThetaMin) / nCosThetaBins);

	// acceptance maps
	TFile* acceptanceFile = TFile::Open("../MonteCarlo/AcceptanceMaps/1S/AcceptanceResults.root", "READ");
	if (!acceptanceFile) {
		cout << "Acceptance file not found. Check the directory of the file." << endl;
		return;
	}

	auto* accMap = (TEfficiency*)acceptanceFile->Get(Form("AccMatrix%s", refFrameName));

	// rebin acceptance maps based on costheta, phi, and pT selection 
	TEfficiency* accMapCosTheta = rebinTEff3DMap(accMap, phiMin, phiMax, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges);

	// efficiency maps
	TFile* efficiencyFile = TFile::Open("../MonteCarlo/EfficiencyMaps/1S/EfficiencyResults.root", "READ");
	if (!efficiencyFile) {
		cout << "Efficiency file not found. Check the directory of the file." << endl;
		return;
	}

	auto* effMap = (TEfficiency*)efficiencyFile->Get(Form("NominalEff_%s", refFrameName));
	auto* systEff = (TH3D*)efficiencyFile->Get(Form("RelatSystEff_%s", refFrameName));

	// rebin efficiency maps based on costheta, phi, and pT selection 
	TEfficiency* effMapCosTheta = rebinTEff3DMap(effMap, phiMin, phiMax, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges);

	// rebin uncertainty map based on costheta, phi, and pT selection 
	TH1D* systEffCosTheta = rebin3DUnc(systEff, phiMin, phiMax, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges);

	Double_t errorWeightLow = 0, errorWeightHigh = 0;

	TCanvas* massCanvas = 0;
	Bool_t isCSframe = (strcmp(refFrameName, "CS") == 0) ? kTRUE : kFALSE;

	Float_t maxYield = 0;

	for (Int_t iCosTheta = 0; iCosTheta < nCosThetaBins; iCosTheta++) {
		Double_t weight = 0;
		// Float_t cosThetaVal = cosThetaBinEdges[iCosTheta];

		// get the corresponding weights
		double acceptance = accMapCosTheta->GetEfficiency(iCosTheta + 1);
		double efficiency = effMapCosTheta->GetEfficiency(iCosTheta + 1);

		weight = 1. / (acceptance * efficiency);
		cout << "acceptance: " << acceptance << endl;
		cout << "efficiency: " << efficiency << endl;
		cout << "weight: " << weight << endl;

		// propagate both scale factor uncertainties and efficiency stat errors to the weight
		errorWeightLow = weight * TMath::Hypot(TMath::Hypot(systEff->GetBinContent(iCosTheta + 1), effMapCosTheta->GetEfficiencyErrorLow(iCosTheta + 1) / efficiency), accMapCosTheta->GetEfficiencyErrorLow(iCosTheta + 1) / acceptance);

		errorWeightHigh = weight * TMath::Hypot(TMath::Hypot(systEff->GetBinContent(iCosTheta + 1), effMapCosTheta->GetEfficiencyErrorUp(iCosTheta + 1) / efficiency), accMapCosTheta->GetEfficiencyErrorUp(iCosTheta + 1) / acceptance);

		const char* fitModelName = GetFitModelName(signalShapeName, ptMin, ptMax, isCSframe, cosThetaBinEdges[iCosTheta], cosThetaBinEdges[iCosTheta + 1], phiMin, phiMax);

		RooArgSet signalYields = GetSignalYields(yield1S, yield2S, yield3S, Form("RawData_%s",bkgShapeName[iCosTheta]), fitModelName);

		standardCorrectedHist->SetBinContent(iCosTheta + 1, (yield1S->getVal()) * weight);
		standardCorrectedHist->SetBinError(iCosTheta + 1, TMath::Hypot(yield1S->getError(), errorWeightLow)); //this one is wrong, will fix it 

		if ((yield1S->getVal()) * weight > maxYield) maxYield = (yield1S->getVal()) * weight;
	}
	cout << "Error: " << yield1S->getError() << endl;

	RooDataHist correctedHist("correctedHist", " ", cosTheta, standardCorrectedHist);

	/// Draw the cos theta distributions

	TCanvas* canvas = new TCanvas("canvas", "canvas", 650, 600);

	RooPlot* frame = cosTheta.frame(Title(" "), Range(cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins]));

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

	auto* polarizationFitResult = cosThetaPDF_1S.fitTo(correctedHist, Save(), Extended(kTRUE) /*, PrintLevel(+1)*/, NumCPU(NCPUs), Range(cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins]), SumW2Error(kFALSE) /*, AsymptoticError(DoAsymptoticError)*/);

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
	canvas->SaveAs(Form("DistributionFits/1D/compareRawYieldCosTheta%s_cent%dto%d_pt%dto%dGeV_phi%dto%d.png", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, phiMin, phiMax), "RECREATE");
}