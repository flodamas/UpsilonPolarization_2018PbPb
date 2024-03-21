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
	
	/// rebin efficiency maps based on costheta, phi, and pT selection 
	
	// extract the numerator and the denominator from the 3D TEfficiency Map
	TH3D* hPassed = (TH3D*) TEff3DMap->GetPassedHistogram();
	TH3D* hTotal = (TH3D*) TEff3DMap->GetTotalHistogram();

	// obtain the bin numbers of the boundaries on phi and pt
	Int_t iPhiMin = hPassed->GetYaxis()->FindBin(phiMin);
	Int_t iPhiMax = hPassed->GetYaxis()->FindBin(phiMax);

	Int_t iPtMin = hPassed->GetZaxis()->FindBin(ptMin);
	Int_t iPtMax = hPassed->GetZaxis()->FindBin(ptMax);

	// obtain the projection histogram along the costheta axis within boundaries of phi and pt 
	// (option e: calculate errors, o: only bins inside the selected range will be filled)
	TH1D* hPassedCosTheta = (TH1D*) hPassed->ProjectionX("hPassedCosTheta", iPhiMin, iPhiMax-1, iPtMin, iPtMax-1, "eo"); 
	TH1D* hTotalCosTheta = (TH1D*) hTotal->ProjectionX("hTotalCosTheta", iPhiMin, iPhiMax-1, iPtMin, iPtMax-1, "eo");

	// rebin the projection histogram (default: 20 bins from -1 to 1)
	TH1D* hPassedCosTheta_Rebin = (TH1D*) hPassedCosTheta->Rebin(nCosThetaBins, "hPassedCosTheta_Rebin", cosThetaBinEdges);
	TH1D* hTotalCosTheta_Rebin = (TH1D*) hTotalCosTheta->Rebin(nCosThetaBins, "hTotalCosTheta_Rebin", cosThetaBinEdges);

	// define TEfficiency using the final numerator and denominator
	TEfficiency* TEffMapCosTheta = new TEfficiency("TEffMapCosTheta", "cos #theta_{CS}; efficiency", nCosThetaBins, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins]);
	
	TEffMapCosTheta->SetPassedHistogram(*hPassedCosTheta_Rebin, "f");
	TEffMapCosTheta->SetTotalHistogram(*hTotalCosTheta_Rebin, "f"); 

	return TEffMapCosTheta;
}

TH1D* rebin3DUnc(TH3D* systEff, Int_t phiMin = -180, Int_t phiMax = 180, Int_t ptMin = 0, Int_t ptMax = 30, Int_t nCosThetaBins = 10, Double_t* cosThetaBinEdges = nullptr){
	
	/// rebin efficiency maps based on costheta, phi, and pT selection 
	// uncertainty addition is sqrt(pow(unc1, 2) + pow(unc2, 2)), so fold it manually
	
	TH1D* h1DSystEff = new TH1D("h1DSystEff", ";cos #theta_{CS};Relative systematic uncertainty", nCosThetaBins, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins]);

	Int_t iPhiMin = systEff->GetYaxis()->FindBin(phiMin);
	Int_t iPhiMax = systEff->GetYaxis()->FindBin(phiMax);

	Int_t iPtMin = systEff->GetZaxis()->FindBin(ptMin);
	Int_t iPtMax = systEff->GetZaxis()->FindBin(ptMax);

	for (int iCosTheta=1; iCosTheta<=nCosThetaBins; iCosTheta++) {

		Double_t phiSumSystEff = 0;

		// sum uncertainties along the phi axis
		for (int iPhi=iPhiMin; iPhi<=iPhiMax; iPhi++) {

			Double_t ptSumSystEff = 0;

			// sum uncertainties along the pt axis
			for (int iPt=iPtMin; iPt<=iPtMax; iPt++) {

				ptSumSystEff = TMath::Hypot(ptSumSystEff, systEff->GetBinContent(iCosTheta, iPhi, iPt));
			}

			phiSumSystEff = TMath::Hypot(phiSumSystEff, ptSumSystEff);
		}

		h1DSystEff->SetBinContent(iCosTheta, phiSumSystEff);
	}	

	return h1DSystEff;
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

	Bool_t isCSframe = (strcmp(refFrameName, "CS") == 0) ? kTRUE : kFALSE;

	Double_t errorWeightLow = 0, errorWeightHigh = 0;
	
	Float_t maxYield = 0;

	TCanvas* massCanvas = 0;

	// define arrays for TGraph (to draw AsymmError)
	double finalDataPoints[nCosThetaBins]; 
	double finalErrHigh[nCosThetaBins]; 
	double finalErrLow[nCosThetaBins];
	double cosThetaBinCenter[nCosThetaBins];

	/// apply weights and errors to each costheta bin
	for (Int_t iCosTheta = 0; iCosTheta < nCosThetaBins; iCosTheta++) {
		Double_t weight = 0;

		// get the corresponding weights
		double acceptance = accMapCosTheta->GetEfficiency(iCosTheta + 1);
		double efficiency = effMapCosTheta->GetEfficiency(iCosTheta + 1);

		weight = 1. / (acceptance * efficiency);
		cout << "acceptance: " << acceptance << endl;
		cout << "efficiency: " << efficiency << endl;
		cout << "weight: " << weight << endl;

		// propagate both scale factor uncertainties and efficiency stat errors to the weight
		double relSystUnc = systEff->GetBinContent(iCosTheta + 1);

		double relEffUncHigh = effMapCosTheta->GetEfficiencyErrorUp(iCosTheta + 1) / efficiency;
		double relEffUncLow = effMapCosTheta->GetEfficiencyErrorLow(iCosTheta + 1) / efficiency;

		double relAccUncHigh = accMapCosTheta->GetEfficiencyErrorUp(iCosTheta + 1) / acceptance;
		double relAccUncLow = accMapCosTheta->GetEfficiencyErrorLow(iCosTheta + 1) / acceptance;

		errorWeightHigh = weight * TMath::Hypot(TMath::Hypot(relSystUnc, relEffUncHigh), relAccUncHigh);

		errorWeightLow = weight * TMath::Hypot(TMath::Hypot(relSystUnc, relEffUncLow), relAccUncLow);

		const char* fitModelName = GetFitModelName(signalShapeName, ptMin, ptMax, isCSframe, cosThetaBinEdges[iCosTheta], cosThetaBinEdges[iCosTheta + 1], phiMin, phiMax);

		RooArgSet signalYields = GetSignalYields(yield1S, yield2S, yield3S, Form("RawData_%s",bkgShapeName[iCosTheta]), fitModelName);

		double yield1SVal = yield1S->getVal();

		double yield1SErr = yield1S->getError();

		standardCorrectedHist->SetBinContent(iCosTheta + 1, yield1SVal * weight);

		// standardCorrectedHist->SetBinError(iCosTheta + 1, yield1SErr * weight); 
		// standardCorrectedHist->SetBinError(iCosTheta + 1, yield1SVal * weight * TMath::Hypot(yield1SErr / yield1SVal, relAccUncHigh)); 
		// standardCorrectedHist->SetBinError(iCosTheta + 1, yield1SVal * weight * TMath::Hypot(yield1SErr / yield1SVal, relEffUncHigh));
		// standardCorrectedHist->SetBinError(iCosTheta + 1, yield1SVal * TMath::Hypot(weight * yield1SErr / yield1SVal, errorWeightLow)); 	

		/// fill arrays for TGraphAsymmError
		cosThetaBinCenter[iCosTheta] = (cosThetaBinEdges[iCosTheta] + cosThetaBinEdges[iCosTheta + 1]) / 2.;

		finalDataPoints[iCosTheta] = yield1SVal * weight;

		finalErrHigh[iCosTheta] = yield1SVal * TMath::Hypot(weight * yield1SErr / yield1SVal, errorWeightHigh);
		finalErrLow[iCosTheta] = yield1SVal * TMath::Hypot(weight * yield1SErr / yield1SVal, errorWeightLow);

		standardCorrectedHist->SetBinError(iCosTheta + 1, finalErrHigh[iCosTheta]); 	

		if ((yield1S->getVal()) * weight > maxYield) maxYield = (yield1S->getVal()) * weight;
		cout << "Error: " << finalErrHigh[iCosTheta] << endl;
	}
	// cout << "Error: " << yield1S->getError() << endl;

	/// TGraphAsymmErrors - comment out for now

	// TGraphAsymmErrors *correctedGraph = new TGraphAsymmErrors(nCosThetaBins, cosThetaBinCenter, finalDataPoints, 0, 0, finalErrLow, finalErrHigh);

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

	// cout << endl
	//      << "Distribution fit for polarization paramaters extraction" << endl;

	// RooRealVar lambdaTheta("lambdaTheta", "lambdaTheta", -1, 1);
	// auto cosThetaPDF_1S = CosThetaPolarizationPDF("cosThetaPDF_1S", " ", cosTheta, lambdaTheta);

   	// enableBinIntegrator(cosThetaPDF_1S, cosTheta.numBins());

	// auto* polarizationFitResult = cosThetaPDF_1S.chi2FitTo(correctedHist, Save(), Extended(kTRUE) /*, PrintLevel(+1)*/, NumCPU(NCPUs), Range(cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins]), SumW2Error(kFALSE) /*, AsymptoticError(DoAsymptoticError)*/);

	// polarizationFitResult->Print("v");

	// cosThetaPDF_1S.plotOn(frame, LineColor(kRed + 1), Name("polaResult"));

	// frame->Draw();

	// // cosmetics

	// TLegend legend(.22, .88, .5, .68);
	// legend.SetTextSize(.05);
	// legend.SetHeader(Form("centrality %d-%d%%, %d < p_{T}^{#mu#mu} < %d GeV/c", gCentralityBinMin, gCentralityBinMax, ptMin, ptMax));
	// legend.AddEntry(frame->findObject("dataPoints"), "#varUpsilon(1S) corrected yield", "lp");
	// legend.AddEntry(frame->findObject("polaResult"), Form("distribution fit: #lambda_{#theta} = %.2f #pm %.2f", lambdaTheta.getVal(), lambdaTheta.getError()), "l");

	// legend.DrawClone();

	// gPad->Update();

	// //CMS_lumi(canvas, gCMSLumiText);

	// gSystem->mkdir("DistributionFits/1D", kTRUE);
	// canvas->SaveAs(Form("DistributionFits/1D/compareRawYieldCosTheta%s_cent%dto%d_pt%dto%dGeV_phi%dto%d.png", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, phiMin, phiMax), "RECREATE");

	// with Root Fit function
	
	TCanvas* canvas2 = new TCanvas("canvas2", "canvas2", 650, 600);

	TF1* PolarFunc = cosThetaPolarFunc(maxYield);

	TFitResultPtr fitResults = standardCorrectedHist->Fit("PolarFunc", "ESV", "", cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins]); //L:log likelihood fit (default: chi2 method), E: NINOS

	cout << "Error of hist (bin1): "<< standardCorrectedHist->GetBinError(1) << endl;
	// Fit results

	double chi2 = fitResults->Chi2();;
	double lambdaVal = fitResults->Parameter(0);
	double lambdaErr = fitResults->ParError(0);	

	gStyle->SetOptFit(1011);

	// cosmetics

	standardCorrectedHist->SetMarkerStyle(20);
    standardCorrectedHist->SetMarkerSize(1);
    standardCorrectedHist->SetMarkerColor(kAzure + 2);

	standardCorrectedHist->GetYaxis()->SetRangeUser(0, 2 * maxYield);
	standardCorrectedHist->GetYaxis()->SetMaxDigits(3);

	standardCorrectedHist->SetXTitle(Form("cos #theta_{%s}", refFrameName));
	// standardCorrectedHist->SetYTitle(Form("Events / ( %0.1f )", cosThetaStep));

	TLegend legend2(.22, .88, .5, .68);
	legend2.SetTextSize(.05);
	legend2.SetHeader(Form("centrality %d-%d%%, %d < p_{T}^{#mu#mu} < %d GeV/c", gCentralityBinMin, gCentralityBinMax, ptMin, ptMax));
	legend2.AddEntry(standardCorrectedHist, "#varUpsilon(1S) corrected yield", "lp");
	legend2.AddEntry(PolarFunc, Form("distribution fit: #lambda_{#theta} = %.2f #pm %.2f", lambdaVal, lambdaErr), "l");

	legend2.DrawClone();

	gPad->Update();

	canvas2->SaveAs(Form("DistributionFits/1D/ROOTFIT_compareCorrectedCosTheta%s_cent%dto%d_pt%dto%dGeV_phi%dto%d.png", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, phiMin, phiMax), "RECREATE");
}