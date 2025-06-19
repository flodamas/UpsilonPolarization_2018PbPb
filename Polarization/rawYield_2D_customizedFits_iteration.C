#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/Parameters/PhaseSpace.h"

#include "../MonteCarlo/AccEffHelpers.h"

#include "../Polarization/PolarFitHelpers.h"

#include "../Tools/FitShortcuts.h"
#include "../Tools/Style/Legends.h"

#include "../Tools/Style/FitDistributions.h"

#include "../Tools/Parameters/CentralityValues.h"
#include "../Tools/Parameters/EfficiencyWeights.h"
#include "../Tools/Parameters/MuonScaleFactors.h"

#include "../ReferenceFrameTransform/Transformations.h"

#include "../Tools/RooFitPDFs/PolarFunc.h"

#include "../Tools/RooFitPDFs/CosThetaPolarizationPDF.h"

#include "../MonteCarlo/acceptanceMap_noGenFilter.C"
#include "../MonteCarlo/weightedEfficiencyMaps.C"

using namespace RooFit;

TEfficiency* getAcceptance3DMap(const char* refFrameName, Double_t lambdaTheta = 0, Double_t lambdaPhi = 0, Double_t lambdaThetaPhi = 0) {
	/// get acceptance and efficiency in 1D
	TString nominalMapName = NominalTEfficiency3DName(refFrameName, lambdaTheta, lambdaPhi, lambdaThetaPhi);

	TString fileName = Form("%s/AcceptanceResults%s.root", AcceptanceResultsPath(gMuonAccName.Data()), "_fullPhi");

	// get acceptance maps
	// TFile* acceptanceFile = openFile(fileName.Data());
	TFile* acceptanceFile = TFile::Open(fileName, "READ");

	if (!acceptanceFile) {
		std::cerr << "Error: acceptanceFile is null." << std::endl;
		acceptanceFile->Close();
		delete acceptanceFile;

		// acceptanceMap_noGenFilter(0, 30, gUpsilonState, lambdaTheta, lambdaPhi, lambdaThetaPhi, isPhiFolded, "MuonUpsilonTriggerAcc");

		return nullptr;
	}

	acceptanceFile = openFile(fileName.Data());

	/// get acceptance maps
	auto* accMap = (TEfficiency*)acceptanceFile->Get(nominalMapName.Data());

	/// check if the acceptance map is null and create it if it doesn't exist
	if (!accMap) {
		std::cerr << "Error: accMap is null." << std::endl;
		acceptanceFile->Close();
		delete acceptanceFile;
		delete accMap;

		// create the acceptance map
		cout << Form("Acceptance map not found. Creating a new one with (lambdaTheta, phi, thetaPhi) = (%.2f, %.2f, %.2f)....", lambdaTheta, lambdaPhi, lambdaThetaPhi) << endl;

		acceptanceMap_noGenFilter(0, 30, kFALSE, "UpsilonTriggerThresholds", lambdaTheta, lambdaPhi, lambdaThetaPhi, gUpsilonState);

        // force ROOT to forget cached file
        gROOT->GetListOfFiles()->Remove(gROOT->GetFile(fileName.Data()));
        delete TFile::Open(fileName.Data());  // flush disk buffers once

		TFile* newAcceptanceFile = openFile(fileName.Data());
        // newAcceptanceFile->cd();
        // newAcceptanceFile->ls();  // now you should see your object

        cout << "nominal map name: " << nominalMapName << endl;
		auto* newAccMap = (TEfficiency*)newAcceptanceFile->Get(nominalMapName.Data());
		cout << "Acceptance map is loaded." << endl;

		/// check if the acceptance map is still null after creating it
		if (!newAccMap) {
			std::cerr << "Error: accMap is still null after creating it." << std::endl;
			// accMap->ls();
			return nullptr;
		}

		return newAccMap;
	}

	return accMap;
}

TEfficiency* getEfficiency3DMap(const char* refFrameName, Double_t lambdaTheta = 0, Double_t lambdaPhi = 0, Double_t lambdaThetaPhi = 0) {
	/// get acceptance and efficiency in 1D
	TString nominalMapName = NominalTEfficiency3DName(refFrameName, lambdaTheta, lambdaPhi, lambdaThetaPhi);

	TString fileName = Form("%s/EfficiencyResults%s.root", EfficiencyResultsPath(gMuonAccName.Data()), "_fullPhi");

	// get efficiency maps
	TFile* efficiencyFile = openFile(fileName.Data());

	if (!efficiencyFile) {
		std::cerr << "Error: efficiencyFile is null." << std::endl;

		return nullptr;
	}

	auto* effMap = (TEfficiency*)efficiencyFile->Get(nominalMapName.Data());

	if (!effMap) {
		std::cerr << "Error: effMap is null." << std::endl;
		weightedEfficiencyMaps(0, 30, "UpsilonTriggerThresholds", lambdaTheta, lambdaPhi, lambdaThetaPhi, kFALSE, gUpsilonState);

        // force ROOT to forget cached file
        gROOT->GetListOfFiles()->Remove(gROOT->GetFile(fileName.Data()));
        delete TFile::Open(fileName.Data());  // flush disk buffers once

		efficiencyFile = openFile(fileName.Data());
        efficiencyFile->cd();
        efficiencyFile->ls();  // now you should see your object

        cout << "nominal map name: " << nominalMapName.Data() << endl;
		effMap = (TEfficiency*)efficiencyFile->Get(nominalMapName.Data());
		cout << "Efficiency map is loaded." << endl;

		if (!effMap) {
			std::cerr << "Error: effMap is still null after creating it." << std::endl;
			// effMap->ls();
			return nullptr;
		}
	}

	return effMap;
}

TH3D* getSysEff3DMap(const char* refFrameName, Double_t lambdaTheta = 0, Double_t lambdaPhi = 0, Double_t lambdaThetaPhi = 0) {
	/// get acceptance and efficiency in 1D
	const char* nominalMapName = SystTEfficiency3DName(refFrameName, lambdaTheta, lambdaPhi, lambdaThetaPhi);

	TString fileName = Form("%s/EfficiencyResults%s.root", EfficiencyResultsPath(gMuonAccName.Data()), "_fullPhi");

	// get relative systematic uncertainty of efficiency
	TFile* efficiencyFile = openFile(fileName.Data());

	if (!efficiencyFile) {
		std::cerr << "Error: efficiencyFile is null." << std::endl;
		return nullptr;
	}

	auto* systEff = (TH3D*)efficiencyFile->Get(nominalMapName);

	if (!systEff) {
		std::cerr << "Error: systEff is null." << std::endl;
		return nullptr;
	}

	return systEff;
}

TLine* drawLine(double x1, double y1, double x2, double y2) {
	TLine* line = new TLine(x1, y1, x2, y2);
	line->SetLineStyle(2);
	line->SetLineColor(kBlack);

	line->Draw("SAME");
	return line;
}

TGraphAsymmErrors* getConfidenceIntervalBands(int nCosThetaBins = 5, const vector<Double_t>& cosThetaBinEdges = {}, int phi = 30, TF2* polarFunc2D = nullptr, TFitResultPtr fitResults = nullptr, double confidenceLevel = 0.6827) {
	/// Get the confidence interval bands for a histogram
	/// Get nominal parameter values
	double n = fitResults->Parameter(0);
	double lambdaTheta = fitResults->Parameter(1);
	double lambdaPhi = fitResults->Parameter(2);
	double lambdaThetaPhi = fitResults->Parameter(3);

	polarFunc2D->SetParameters(n, lambdaTheta, lambdaPhi, lambdaThetaPhi);

	TMatrixDSym cov = fitResults->GetCovarianceMatrix();

	TGraphAsymmErrors* band = new TGraphAsymmErrors();

	const int nPoints = 1000; // Number of points for the band
	double step = (cosThetaBinEdges[nCosThetaBins] - cosThetaBinEdges[0]) / nPoints;	 
	
	for (int i = 0; i < nPoints; ++i) {
		double x = cosThetaBinEdges[0] + i * step;

		/// Set phi value
		// double phiRad = phi * TMath::Pi() / 180.0; // Convert phi to radians

		/// Evaluate the function 
		double y = polarFunc2D->Eval(x, phi);

		// Compute derivatives
		double dfdn = 1.0 / (3 + lambdaTheta) * (1 + lambdaTheta * x * x + lambdaPhi * pow(sin(acos(x)), 2) * cos(2 * phi) + lambdaThetaPhi * sin(2 * acos(x)) * cos(phi));
		double dfdlambdaTheta = -n / pow(3 + lambdaTheta, 2) * (1 + lambdaTheta * x * x + lambdaPhi * pow(sin(acos(x)), 2) * cos(2 * phi) + lambdaThetaPhi * sin(2 * acos(x)) * cos(phi)) + n / (3 + lambdaTheta) * (x * x);
		double dfdlambdaPhi = n / (3 + lambdaTheta) * (pow(sin(acos(x)), 2) * cos(2 * phi));
		double dfdlambdaThetaPhi = n / (3 + lambdaTheta) * (sin(2 * acos(x)) * cos(phi));

		TVectorD grad(4);
		grad[0] = dfdn;
		grad[1] = dfdlambdaTheta;
		grad[2] = dfdlambdaPhi;
		grad[3] = dfdlambdaThetaPhi;

		// Error^2 = grad^T * Cov * grad
		
		double err2 = cov.Similarity(grad);
		double z = TMath::NormQuantile(0.5 * (1 + confidenceLevel));  // e.g., 1.96 for 95%
		double err = z * sqrt(err2);

		// double err = sqrt(err2);

		band->SetPoint(i, x, y);
		band->SetPointError(i, 0, 0, err, err); // symmetric
	}

	return band;
}

/// fit the corrected histo and extract polarization parameters
RooArgSet extractPolarParam(TH2D* correctedHist, std::vector<TH1D*>& correctedHist1DCosTheta, TString refFrameName = "CS",
                            Int_t ptMin = 0, Int_t ptMax = 30,
                            Int_t nCosThetaBins = 5, const vector<Double_t>& cosThetaBinEdges = {},
                            Int_t nPhiBins = 6, const vector<Double_t>& phiBinEdges = {},
                            Bool_t isPhiFolded = kFALSE, Bool_t applyEff = kTRUE) {
	
	writeExtraText = true; // if extra text
	extraText = "       Internal";
	// extraText = "       Simulation Preliminary";

	TH2D* fitHist = (TH2D*)correctedHist->Clone("fitHist");

	Float_t maxYield = fitHist->GetMaximum();

	TF2* polarFunc2D = getGeneralPolarFunc(maxYield);

	/// Polarization fit
	/// with Root Fit function (E: minos on, S: save results, V: Verbose, I: integral, M: imporve algorithm, R: specificed range)
	TFitResultPtr fitResults = fitHist->Fit("generalPolarFunc", "ESVIMR"); // chi2 fit to the integrated bin

	/// Fit results
	double chi2 = fitResults->Chi2();
	double nDOF = nCosThetaBins * nPhiBins - polarFunc2D->GetNpar();

	double pValue = TMath::Prob(chi2, nDOF); // https://root.cern.ch/doc/master/namespaceTMath.html#a3aba6abaf2bc605680c7be8fc5fd98aa

	double normVal = fitResults->Parameter(0);
	double normErr = fitResults->ParError(0);

	double lambdaThetaVal = fitResults->Parameter(1);
	double lambdaThetaErr = fitResults->ParError(1);

	double lambdaPhiVal = fitResults->Parameter(2);
	double lambdaPhiErr = fitResults->ParError(2);

	double lambdaThetaPhiVal = fitResults->Parameter(3);
	double lambdaThetaPhiErr = fitResults->ParError(3);

	double lambdaTildeVal = (lambdaThetaVal + 3. * lambdaPhiVal) / (1. - lambdaPhiVal);
	double lambdaTildeErr = TMath::Hypot(1. / (1. - lambdaPhiVal) * lambdaPhiErr, (3. - lambdaThetaVal - 6. * lambdaPhiVal) / TMath::Power((1. - lambdaPhiVal), 2) * lambdaPhiErr);

	RooRealVar* lambdaTheta = new RooRealVar("lambdaTheta", "lambdaTheta", -3., 3.);
	RooRealVar* lambdaPhi = new RooRealVar("lambdaPhi", "lambdaPhi", -3., 3.);
	RooRealVar* lambdaThetaPhi = new RooRealVar("lambdaThetaPhi", "lambdaThetaPhi", -3., 3.);
	RooRealVar* lambdaTilde = new RooRealVar("lambdaTilde", "lambdaTilde", -3., 3.);

	lambdaTheta->setVal(lambdaThetaVal);
	lambdaTheta->setError(lambdaThetaErr);

	lambdaPhi->setVal(lambdaPhiVal);
	lambdaPhi->setError(lambdaPhiErr);

	lambdaThetaPhi->setVal(lambdaThetaPhiVal);
	lambdaThetaPhi->setError(lambdaThetaPhiErr);

	lambdaTilde->setVal(lambdaTildeVal);
	lambdaTilde->setError(lambdaTildeErr);

	/// store fit results
	RooArgSet fittedParams(*lambdaTheta, *lambdaPhi, *lambdaThetaPhi, *lambdaTilde);

	/// draw the fit results
	TCanvas* mc2DCanvas = draw2DMap(fitHist, refFrameName.Data(), nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kTRUE, kFALSE, 1, isPhiFolded);

	fitHist->GetZaxis()->SetRangeUser(0, maxYield * 2);
	cout << "maxYield: " << maxYield << endl;

	fitHist->GetZaxis()->SetTitle("#varUpsilon(1S) weighted events");
	fitHist->GetZaxis()->SetTitleSize(0.055);
	fitHist->GetZaxis()->SetTitleOffset(1.6);
	fitHist->GetZaxis()->SetLabelSize(0.044);

	fitHist->GetYaxis()->SetTitleSize(0.055);
	fitHist->GetYaxis()->SetTitleOffset(1.6);
	fitHist->GetYaxis()->SetLabelSize(0.044);
	fitHist->GetYaxis()->SetLabelOffset(0);

	fitHist->GetXaxis()->SetTitleSize(0.055);
	fitHist->GetXaxis()->SetTitleOffset(1.2);
	fitHist->GetXaxis()->SetLabelSize(0.044);

	// polarFunc2D->Draw("SURFACE SAME");
	TPaveText* kinematicsText = new TPaveText(0.28, 0.89, 0.85, 0.93, "NDCNB"); //0.375,0.8521739,0.9459877,0.9043478
	kinematicsText->SetFillColor(4000);
	kinematicsText->SetBorderSize(0);
	kinematicsText->AddText(Form("%s, %s", DimuonRapidityRangeText(gRapidityMin, gRapidityMax), DimuonPtRangeText(ptMin, ptMax)));
	kinematicsText->SetAllWith("", "align", 12);
	kinematicsText->Draw("SAME");

	/// add results to the plot
	TLegend legend2(.16, .78, .46, .89); // 0.159292,0.7791304,0.460177,0.8886957
	legend2.SetTextSize(.045);
	legend2.SetFillColor(0);
	legend2.SetFillStyle(1001);

	legend2.AddEntry(fitHist, "#varUpsilon(1S) corrected yield", "lp");
	// legend2.AddEntry(polarFunc2D, Form("fit: #lambda_{#theta}  = %.2f #pm %.2f    #lambda_{#varphi} = %.2f #pm %.2f", lambdaThetaVal, lambdaThetaErr, lambdaPhiVal, lambdaPhiErr), "l");
	// legend2.AddEntry((TObject*)0, Form("     #lambda_{#theta#varphi} = %.2f #pm %.2f  #tilde{#lambda}  = %.2f #pm %.2f", lambdaThetaPhiVal, lambdaThetaPhiErr, lambdaTildeVal, lambdaTildeErr), "");
	legend2.AddEntry(polarFunc2D, "fit: ", "l");

	legend2.DrawClone();

	TPaveText* resultTextRight = new TPaveText(0.28, 0.72, 0.65, 0.83, "NDC"); // Adjust coordinates
	resultTextRight->SetFillColor(0);                                          // White background
	resultTextRight->SetFillStyle(1001);                                       // Solid fill
	resultTextRight->SetBorderSize(0);                                         // Optional: Thin border
	resultTextRight->SetTextSize(0.045);
	resultTextRight->SetTextAlign(12); // Align text left
	resultTextRight->SetMargin(0.03);
	resultTextRight->AddText(Form("#lambda_{#theta}  = %.2f #pm %.2f ", lambdaThetaVal, lambdaThetaErr));
	resultTextRight->AddText(Form("#lambda_{#theta#varphi} = %.2f #pm %.2f", lambdaThetaPhiVal, lambdaThetaPhiErr));
	resultTextRight->Draw();

	TPaveText* resultTextLeft = new TPaveText(0.57, 0.74, 0.83, 0.825, "NDC"); // Adjust coordinates
	resultTextLeft->SetFillColor(0);                                           // White background
	resultTextLeft->SetFillStyle(1001);                                        // Solid fill
	resultTextLeft->SetBorderSize(0);                                          // Optional: Thin border
	resultTextLeft->SetTextSize(0.045);
	resultTextLeft->SetTextAlign(12); // Align text left
	resultTextLeft->SetMargin(0.03);
	resultTextLeft->AddText(Form("#lambda_{#varphi} = %.2f #pm %.2f", lambdaPhiVal, lambdaPhiErr));
	resultTextLeft->AddText(Form("#tilde{#lambda}  = %.2f #pm %.2f", lambdaTildeVal, lambdaTildeErr));
	resultTextLeft->Draw();

	TLatex textChi2;
	textChi2.SetTextAlign(22);
	textChi2.SetTextSize(0.04);
	textChi2.DrawLatexNDC(0.72, 0.044, Form("#chi^{2} / n_{dof} = %.2f, p-value = %.2f", chi2 / nDOF, pValue));

	// gPad->RedrawAxis();

	gPad->Update();

	/// draw 1D histograms along the cosTheta axis with the fit
	TCanvas* corrected1DFitCanvas = new TCanvas("corrected1DFitCanvas", "corrected1DFitCanvas", 1500, 500);
	corrected1DFitCanvas->Divide(3, 1);

	for (Int_t iPhi = 0; iPhi < nPhiBins; iPhi++) {
		/// cosmetics for 1D histograms along the cosTheta axis
		corrected1DFitCanvas->cd(iPhi + 1);
		
		correctedHist1DCosTheta[iPhi]->GetXaxis()->SetTitle(CosThetaVarTitle(refFrameName.Data()));
		correctedHist1DCosTheta[iPhi]->GetXaxis()->CenterTitle();
		correctedHist1DCosTheta[iPhi]->GetXaxis()->SetNdivisions(-nCosThetaBins);

		correctedHist1DCosTheta[iPhi]->GetYaxis()->SetTitle("#varUpsilon(1S) weighted events");
		correctedHist1DCosTheta[iPhi]->GetYaxis()->SetRangeUser(0, correctedHist->GetMaximum() * 1.5);

		correctedHist1DCosTheta[iPhi]->SetMarkerStyle(20);
		correctedHist1DCosTheta[iPhi]->SetMarkerSize(1.2);
		correctedHist1DCosTheta[iPhi]->SetLineWidth(2);
		
		correctedHist1DCosTheta[iPhi]->Draw("E");

		/// slice of the 2D fit along the cosTheta axis
		TF1* polarFuncCosTheta = new TF1("polarFuncCosTheta", "[0] / (3 + [1]) * (1 + [1] * x * x + [2] * TMath::Sin(TMath::ACos(x)) * TMath::Sin(TMath::ACos(x)) * TMath::Cos(2. * [4] * pi / 180.) + [3] * TMath::Sin(2 * TMath::ACos(x)) * TMath::Cos([4] * pi / 180.))", -1, 1);
		polarFuncCosTheta->SetParameter(0, normVal);
		polarFuncCosTheta->SetParameter(1, lambdaThetaVal);
		polarFuncCosTheta->SetParameter(2, lambdaPhiVal);
		polarFuncCosTheta->SetParameter(3, lambdaThetaPhiVal);
		
		double phiBinCenter = phiBinEdges[iPhi] + (phiBinEdges[iPhi + 1] - phiBinEdges[iPhi]) / 2.;
		polarFuncCosTheta->SetParameter(4, phiBinCenter); // set the phi value for the fit

		polarFuncCosTheta->SetLineColor(kRed);
		polarFuncCosTheta->SetLineWidth(2);
		polarFuncCosTheta->SetNpx(1000);

		/// 1-, 2-, 3-Sigma band
		TGraphAsymmErrors* errorBand1sigma = getConfidenceIntervalBands(nCosThetaBins, cosThetaBinEdges, phiBinCenter, polarFunc2D, fitResults, 0.6827);
		TGraphAsymmErrors* errorBand2sigma = getConfidenceIntervalBands(nCosThetaBins, cosThetaBinEdges, phiBinCenter, polarFunc2D, fitResults, 0.9545);
		TGraphAsymmErrors* errorBand3sigma = getConfidenceIntervalBands(nCosThetaBins, cosThetaBinEdges, phiBinCenter, polarFunc2D, fitResults, 0.9973);

		// new TH1D("errorBand", "", 1000, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins]);
		// (TVirtualFitter::GetFitter())->GetConfidenceIntervals(errorBand, 0.68);
		
		gStyle->SetOptFit(1011);

		/// cosmetics
		errorBand1sigma->SetFillColorAlpha(kRed - 9, 0.2);
		errorBand1sigma->SetMarkerSize(0);

		errorBand2sigma->SetFillColorAlpha(kOrange-3, 0.2);
		errorBand2sigma->SetMarkerSize(0);

		errorBand3sigma->SetFillColorAlpha(kAzure+1, 0.15);
		errorBand3sigma->SetMarkerSize(0);

		errorBand3sigma->Draw("E3 SAME");	
		errorBand2sigma->Draw("E3 SAME");
		errorBand1sigma->Draw("E3 SAME");

		polarFuncCosTheta->Draw("SAME");
				
		// /// contour plot
		// // (ref: https://root-forum.cern.ch/t/roofit-minos-errors-for-2-parameters-of-interest/16157)
		
		// // set the confidence level
		// gMinuit->SetErrorDef(2.30); // 1 sigma corresponds to delchi2 = 2.30 
		// TGraph* contourPlot1 = (TGraph*)gMinuit->Contour(1000, 1, 0); // Contour(number of points, lambda_theta, normalization factor)

		// gMinuit->SetErrorDef(6.18); // 2 sigma corresponds to delchi2 = 6.18
		// TGraph* contourPlot2 = (TGraph*)gMinuit->Contour(1000, 1, 0);	

		// TCanvas* contourCanvas = drawContourPlots(ptMin, ptMax, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], refFrameName, contourPlot1, contourPlot2);


		/// legend
		TLatex* text1D = new TLatex();
		text1D->SetTextAlign(22);
		text1D->SetTextSize(0.05);
		/// Put texts inside the plot
		text1D->DrawLatexNDC(.56, .84, Form("%s, %d < p_{T} < %d GeV/c, %d < |#varphi| < %d #circ", refFrameName.Data(), ptMin, ptMax, (int)phiBinEdges[iPhi], (int)phiBinEdges[iPhi + 1]));
		// text1D->DrawLatexNDC(.55, .72, Form("(#lambda_{#theta}, #lambda_{#varphi}, #lambda_{#theta#varphi}) = (%.2f#pm%.2f, %.2f#pm%.2f, %.2f#pm%.2f)", lambdaThetaVal, lambdaThetaErr, lambdaPhiVal, lambdaPhiErr, lambdaThetaPhiVal, lambdaThetaPhiErr));
		
		TLegend* legend1D = new TLegend(.17, .25, .48, .35);
		legend1D->SetTextSize(0.04);
		legend1D->SetFillColor(0);
		legend1D->SetFillStyle(1001);
		legend1D->SetBorderSize(0);
		legend1D->AddEntry(correctedHist1DCosTheta[iPhi], "#varUpsilon(1S) corrected yield", "pe");
		legend1D->AddEntry(polarFuncCosTheta, "fit:", "l");
		legend1D->Draw("SAME");

		TPaveText* resultTextRight1D = new TPaveText(0.29, 0.185, 0.66, 0.295, "NDC"); // Adjust coordinates
		resultTextRight1D->SetFillColor(0);                                          // White background
		resultTextRight1D->SetFillStyle(1001);                                       // Solid fill
		resultTextRight1D->SetBorderSize(0);                                         // Optional: Thin border
		resultTextRight1D->SetTextSize(0.045);
		resultTextRight1D->SetTextAlign(12); // Align text left
		resultTextRight1D->SetMargin(0.03);
		resultTextRight1D->AddText(Form("#lambda_{#theta}  = %.2f #pm %.2f ", lambdaThetaVal, lambdaThetaErr));
		resultTextRight1D->AddText(Form("#lambda_{#theta#varphi} = %.2f #pm %.2f", lambdaThetaPhiVal, lambdaThetaPhiErr));
		resultTextRight1D->Draw();

		TPaveText* resultTextLeft1D = new TPaveText(0.61, 0.20, 0.88, 0.295, "NDC"); // Adjust coordinates
		resultTextLeft1D->SetFillColor(0);                                           // White background
		resultTextLeft1D->SetFillStyle(1001);                                        // Solid fill
		resultTextLeft1D->SetBorderSize(0);                                          // Optional: Thin border
		resultTextLeft1D->SetTextSize(0.045);
		resultTextLeft1D->SetTextAlign(12); // Align text Left1D
		resultTextLeft1D->SetMargin(0.03);
		resultTextLeft1D->AddText(Form("#lambda_{#varphi} = %.2f #pm %.2f", lambdaPhiVal, lambdaPhiErr));
		resultTextLeft1D->AddText(Form("#tilde{#lambda}  = %.2f #pm %.2f", lambdaTildeVal, lambdaTildeErr));
		resultTextLeft1D->Draw();

		CMS_lumi((TPad*)gPad, gCMSLumiText);	
		
		gPad->RedrawAxis();

		gPad->Update();
	}

	gSystem->mkdir(Form("closureTest/%s", applyEff ? "HydjetAccEff": "PythiaAcc"), kTRUE);
	if (isPhiFolded)
		mc2DCanvas->SaveAs(Form("closureTest/%s/2Dfit_%s_pt%dto%d_cosTheta%.2fto%.2f_absphi%dto%d_lambdaTheta%.2f_Phi%.2f_ThetaPhi%.2f.png", applyEff ? "HydjetAccEff": "PythiaAcc", refFrameName.Data(), ptMin, ptMax, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], (int)phiBinEdges[0], (int)phiBinEdges[nPhiBins], lambdaThetaVal, lambdaPhiVal, lambdaThetaPhiVal));
	else
		mc2DCanvas->SaveAs(Form("closureTest/%s/2Dfit_%s_pt%dto%d_cosTheta%.2fto%.2f_phi%dto%d_lambdaTheta%.2f_Phi%.2f_ThetaPhi%.2f.png", applyEff ? "HydjetAccEff": "PythiaAcc", refFrameName.Data(), ptMin, ptMax, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], (int)phiBinEdges[0], (int)phiBinEdges[nPhiBins], lambdaThetaVal, lambdaPhiVal, lambdaThetaPhiVal));

	return fittedParams;
}

/// correct the polarized MC with weights (1 / (acceptance and/or efficiency))
void correctMC2DHist(TH2D* polarizedHist, TH2D* correctedHist, std::vector<TH1D*>& rawYieldHist1DCosTheta, std::vector<TH1D*>& rawYieldHist1DPhi, std::vector<TH1D*>& correctedHist1DCosTheta, std::vector<TH1D*>& correctedHist1DPhi, TString refFrameName = "CS",
                     RooArgSet polarParams = RooArgSet(),
                     Int_t ptMin = 0, Int_t ptMax = 30,
                     Int_t nCosThetaBins = 20, const vector<Double_t>& cosThetaBinEdges = {},
                     Int_t nPhiBins = 18, const vector<Double_t>& phiBinEdges = {},
                     Bool_t drawPlot = kFALSE, Bool_t isPhiFolded = kTRUE, Bool_t applyAcc = kTRUE, Bool_t applyEff = kFALSE,
					 float lambdaTheta0 = 0, float lambdaPhi0 = 0, float lambdaThetaPhi0 = 0, 
					 int nLegendRows = 1) {
	writeExtraText = true; // if extra text
	extraText = "       Internal";
	// extraText = "       Simulation Preliminary";

	/// initialize variables used in the loop below
	Double_t totalRelUncHigh = 0, totalRelUncLow = 0;
	Double_t totalUncHigh = 0, totalUncLow = 0;

	Float_t maxYield = 0;

	RooRealVar* lambdaThetaVar = (RooRealVar*)polarParams.find("lambdaTheta");
	RooRealVar* lambdaPhiVar = (RooRealVar*)polarParams.find("lambdaPhi");
	RooRealVar* lambdaThetaPhiVar = (RooRealVar*)polarParams.find("lambdaThetaPhi");

	/// set the polarization parameters for acceptance and efficiency MC samples to the input polarization parameters
	float lambdaTheta = lambdaThetaVar->getVal();
	float lambdaPhi = lambdaPhiVar->getVal();
	float lambdaThetaPhi = lambdaThetaPhiVar->getVal();

	/// get acceptance and efficiency 3D map
	TEfficiency* accMap;
	TEfficiency* effMap;
	TH3D* systEff;
	
	if (!applyEff) {
		/// apply underlying polarization parameters to only acceptance
		accMap = getAcceptance3DMap(refFrameName.Data(), lambdaTheta, lambdaPhi, lambdaThetaPhi);
		// accMap = getAcceptance3DMap(refFrameName.Data(), 1, 0, 0);
    	effMap = getEfficiency3DMap(refFrameName.Data(), 0, 0, 0);
		systEff = getSysEff3DMap(refFrameName.Data(), 0, 0, 0);
	}
	
	else {
		/// apply underlying polarization parameters to only efficiency
		// accMap = getAcceptance3DMap(refFrameName.Data(), 0, 0, 0);
		// effMap = getEfficiency3DMap(refFrameName.Data(), lambdaTheta, lambdaPhi, lambdaThetaPhi);
    	// systEff = getSysEff3DMap(refFrameName.Data(), lambdaTheta, lambdaPhi, lambdaThetaPhi);

		/// apply underlying polarization parameters to both acceptance and efficiency maps
		accMap = getAcceptance3DMap(refFrameName.Data(), lambdaTheta, lambdaPhi, lambdaThetaPhi);
    	effMap = getEfficiency3DMap(refFrameName.Data(), lambdaTheta, lambdaPhi, lambdaThetaPhi);
		systEff = getSysEff3DMap(refFrameName.Data(), lambdaTheta, lambdaPhi, lambdaThetaPhi);
	}

	/// rebin acceptance and efficiency, efficiency systematic Uncertainty
	TEfficiency* accMapCosThetaPhi = rebinTEff3DMap(accMap, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, isPhiFolded);
	TEfficiency* effMapCosThetaPhi = rebinTEff3DMap(effMap, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, isPhiFolded);
	TH2D* systEffCosThetaPhi = rebinRel3DUncMap(effMap, systEff, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, isPhiFolded);

	TH2D* hTotalCosThetaPhiAcc = (TH2D*)(accMapCosThetaPhi->GetTotalHistogram())->Clone("hTotalCosThetaPhiAcc");
	hTotalCosThetaPhiAcc->SetName("hTotalCosThetaPhiAcc");

	TH2D* hPassedCosThetaPhiAcc = (TH2D*)accMapCosThetaPhi->GetPassedHistogram();
	hPassedCosThetaPhiAcc->SetName("hPassedCosThetaPhiAcc");

	TH2D* hTotalCosThetaPhiEff = (TH2D*)(effMapCosThetaPhi->GetTotalHistogram())->Clone("hTotalCosThetaPhiEff");
	hTotalCosThetaPhiEff->SetName("hTotalCosThetaPhiEff");

	TH2D* hPassedCosThetaPhiEff = (TH2D*)effMapCosThetaPhi->GetPassedHistogram();
	hPassedCosThetaPhiEff->SetName("hPassedCosThetaPhiEff");

	TH2D* hRatioCosThetaPhi = (TH2D*)hTotalCosThetaPhiAcc->Clone("hRatioCosThetaPhi");

	hRatioCosThetaPhi->Divide(hPassedCosThetaPhiAcc, hTotalCosThetaPhiEff, 1, 1);
	
	/// draw acceptance and efficiency graph for check
	TCanvas* accCanvas = nullptr;
	TCanvas* effCanvas = nullptr;
	
	TCanvas* accTotalCanvas = nullptr;
	TCanvas* accPassedCanvas = nullptr;
	TCanvas* effTotalCanvas = nullptr;
	TCanvas* effPassedCanvas = nullptr;

	TCanvas* ratioCanvas = nullptr;
	
	/// draw 2D acc and eff maps
	if (!applyEff) {
		accCanvas = DrawEfficiency2DHist(accMapCosThetaPhi, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, gUpsilonState, kTRUE, kTRUE, kFALSE, "_TriggerAcc", isPhiFolded, kTRUE, lambdaTheta, lambdaPhi, lambdaThetaPhi, nLegendRows);
		effCanvas = DrawEfficiency2DHist(effMapCosThetaPhi, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, gUpsilonState, kFALSE, kFALSE, kFALSE, "_TriggerAcc", isPhiFolded, kTRUE, 0, 0, 0, nLegendRows);
	
		accTotalCanvas = draw2DMap(hTotalCosThetaPhiAcc, refFrameName.Data(), nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kFALSE, 1, isPhiFolded, nLegendRows);
        hTotalCosThetaPhiAcc->GetZaxis()->SetTitle("#varUpsilon(1S) acceptance denominator");
		display2DMapContents(hTotalCosThetaPhiAcc, nCosThetaBins, nPhiBins, kFALSE, 0.04, kBlack);
         
        accPassedCanvas = draw2DMap(hPassedCosThetaPhiAcc, refFrameName.Data(), nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kFALSE, 1, isPhiFolded, nLegendRows);
        hPassedCosThetaPhiAcc->GetZaxis()->SetTitle("#varUpsilon(1S) acceptance numerator");
		display2DMapContents(hPassedCosThetaPhiAcc, nCosThetaBins, nPhiBins, kFALSE);
        
	}

	else {
		accCanvas = DrawEfficiency2DHist(accMapCosThetaPhi, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, gUpsilonState, kTRUE, kTRUE, kFALSE, "_TriggerAcc", isPhiFolded, kTRUE, lambdaTheta, lambdaPhi, lambdaThetaPhi, nLegendRows);
		effCanvas = DrawEfficiency2DHist(effMapCosThetaPhi, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, gUpsilonState, kFALSE, kTRUE, kFALSE, "_TriggerAcc", isPhiFolded, kTRUE, lambdaTheta, lambdaPhi, lambdaThetaPhi, nLegendRows);
        
		ratioCanvas = draw2DMap(hRatioCosThetaPhi, refFrameName.Data(), nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kFALSE, 1, isPhiFolded, nLegendRows);
		display2DMapContents(hRatioCosThetaPhi, nCosThetaBins, nPhiBins, kFALSE, 0.04, kBlack, 4);

		accTotalCanvas = draw2DMap(hTotalCosThetaPhiAcc, refFrameName.Data(), nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kFALSE, 1, isPhiFolded, nLegendRows);
		hTotalCosThetaPhiAcc->GetZaxis()->SetTitle("#varUpsilon(1S) acceptance denominator");
		display2DMapContents(hTotalCosThetaPhiAcc, nCosThetaBins, nPhiBins, kFALSE);
       
        accPassedCanvas = draw2DMap(hPassedCosThetaPhiAcc, refFrameName.Data(), nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kFALSE, 1, isPhiFolded, nLegendRows);
        hPassedCosThetaPhiAcc->GetZaxis()->SetTitle("#varUpsilon(1S) acceptance numerator");
		display2DMapContents(hPassedCosThetaPhiAcc, nCosThetaBins, nPhiBins, kFALSE);
        
        effTotalCanvas = draw2DMap(hTotalCosThetaPhiEff, refFrameName.Data(), nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kFALSE, 1, isPhiFolded, nLegendRows);
        hTotalCosThetaPhiEff->GetZaxis()->SetTitle("#varUpsilon(1S) efficiency denominator");
		display2DMapContents(hTotalCosThetaPhiEff, nCosThetaBins, nPhiBins, kFALSE);

        effPassedCanvas = draw2DMap(hPassedCosThetaPhiEff, refFrameName.Data(), nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kFALSE, 1, isPhiFolded, nLegendRows);
        hPassedCosThetaPhiEff->GetZaxis()->SetTitle("#varUpsilon(1S) efficiency numerator");
		display2DMapContents(hPassedCosThetaPhiEff, nCosThetaBins, nPhiBins, kFALSE);
	}
	
	TCanvas* dummyCanvas = new TCanvas("dummyCanvas", "dummyCanvas", 600, 600); // this is dummy canvas due to the overwrite of the fitted histogram

	/// apply acc x eff correction weights and errors to each costheta bin
	for (Int_t iPhi = 0; iPhi < nPhiBins; iPhi++) {
		correctedHist1DCosTheta[iPhi] = new TH1D(Form("correctedHist1DCosTheta_phi%dto%d", (int)phiBinEdges[iPhi], (int)phiBinEdges[iPhi + 1]), "", nCosThetaBins, cosThetaBinEdges.data());
		rawYieldHist1DCosTheta[iPhi] = new TH1D(Form("rawYieldHist1DCosTheta_phi%dto%d", (int)phiBinEdges[iPhi], (int)phiBinEdges[iPhi + 1]), "", nCosThetaBins, cosThetaBinEdges.data());

		for (Int_t iCosTheta = 0; iCosTheta < nCosThetaBins; iCosTheta++) {
			if (iPhi == 0) {
				rawYieldHist1DPhi[iCosTheta] = new TH1D(Form("rawYieldHist1DPhi_cosTheta%dto%d", (int)cosThetaBinEdges[iCosTheta], (int)cosThetaBinEdges[iCosTheta + 1]), "", nPhiBins, phiBinEdges.data());
				correctedHist1DPhi[iCosTheta] = new TH1D(Form("correctedHist1DPhi_cosTheta%dto%d", (int)cosThetaBinEdges[iCosTheta], (int)cosThetaBinEdges[iCosTheta + 1]), "", nPhiBins, phiBinEdges.data());
			}

			Double_t weight = 1;

			// get the global bin number of Efficiency
			Double_t binCenterCosTheta = hTotalCosThetaPhiAcc->GetXaxis()->GetBinCenter(iCosTheta + 1);
			Double_t binCenterPhi = hTotalCosThetaPhiAcc->GetYaxis()->GetBinCenter(iPhi + 1);

			Int_t iGlobalBin = hTotalCosThetaPhiAcc->FindFixBin(binCenterCosTheta, binCenterPhi);

			// get the corresponding weights
			double acceptance = accMapCosThetaPhi->GetEfficiency(iGlobalBin);
			double efficiency = effMapCosThetaPhi->GetEfficiency(iGlobalBin);
			double residual = hRatioCosThetaPhi->GetBinContent(iCosTheta + 1, iPhi + 1);
			// double efficiency = 1;

			cout << "acceptance: " << acceptance << endl;
			cout << "efficiency: " << efficiency << endl;

			double relEffUncHigh = 0;
			double relEffUncLow = 0;

			double relAccUncHigh = 0;
			double relAccUncLow = 0;

			/// calculate weight
			/// apply both eff and acc
			if (applyAcc && applyEff) {
				if (acceptance == 0 || efficiency == 0) {
					weight = 0.;
					relEffUncHigh = 0.;
					relAccUncHigh = 0.;
				}

				else {
					// weight = 1. / (acceptance * efficiency);
					weight = 1. / (acceptance * efficiency) * residual;

					relEffUncHigh = effMapCosThetaPhi->GetEfficiencyErrorUp(iGlobalBin) / efficiency;
					relAccUncHigh = accMapCosThetaPhi->GetEfficiencyErrorUp(iGlobalBin) / acceptance;
				}
			}

			/// apply only eff
			else if (!applyAcc && applyEff) {
				if (efficiency == 0) {
					weight = 0.;
					relEffUncHigh = 0.;
				}
				else {
					weight = 1. / (efficiency);
					relEffUncHigh = effMapCosThetaPhi->GetEfficiencyErrorUp(iGlobalBin) / efficiency;
				}
			}

			/// apply only acc
			else if (applyAcc && !applyEff) {
				if (acceptance == 0) {
					weight = 0.;
					relAccUncHigh = 0.;
				}
				else {
					weight = 1. / (acceptance);
					relAccUncHigh = accMapCosThetaPhi->GetEfficiencyErrorUp(iGlobalBin) / acceptance;
				}
			}

			/// no detector effects
			else
				weight = 1.;

			// weightMap->SetBinContent(iCosTheta + 1, iPhi + 1, weight);

			// propagate both scale factor uncertainties and efficiency stat errors to the weight
			// double relSystUnc = systEff->GetBinContent(iCosTheta + 1, iPhi + 1) / efficiency;

			// totalRelUncHigh = TMath::Hypot(TMath::Hypot(relSystUnc, relEffUncHigh), relAccUncHigh);
			// totalRelUncLow = TMath::Hypot(TMath::Hypot(relSystUnc, relEffUncLow), relAccUncLow);

			// totalRelUncHigh = TMath::Hypot(relEffUncHigh, relAccUncHigh);
			// totalRelUncLow = TMath::Hypot(relEffUncLow, relAccUncLow);

			// totalUncHigh = totalRelUncHigh * efficiency * acceptance;
			// totalUncLow = totalRelUncLow * efficiency * acceptance;

			// full relative uncertainty on weight
			double relWeightUnc = TMath::Hypot(relAccUncHigh, relEffUncHigh);
			double deltaW = weight * relWeightUnc;

			double recoMCVal = polarizedHist->GetBinContent(iCosTheta + 1, iPhi + 1);

			double recoMCUnc = polarizedHist->GetBinError(iCosTheta + 1, iPhi + 1);

			double residualVal = hRatioCosThetaPhi->GetBinContent(iCosTheta + 1, iPhi + 1);
			double residualUnc = hRatioCosThetaPhi->GetBinError(iCosTheta + 1, iPhi + 1);

			// cout << "recoMCVal: " << recoMCVal << endl;
			// cout << "recoMCUnc: " << recoMCUnc << endl;

			// now propagate both sources of uncertainty:
			double totalError = TMath::Sqrt(pow(weight * recoMCUnc, 2) + pow(recoMCVal * deltaW, 2)); // error from acc × eff)

			// cout << "weight: " << weight << endl;
			// cout << "relEffUncHigh: " << relEffUncHigh << endl;
			// cout << "relEffUncLow: " << relEffUncLow << endl;
			// cout << "relAccUncHigh: " << relAccUncHigh << endl;
			// cout << "relAccUncLow: " << relAccUncLow << endl;
			// cout << "relWeightUnc: " << relWeightUnc << endl;
			// cout << "deltaW: " << deltaW << endl;
			// cout << "residualVal: " << residualVal << endl;
			// cout << "residualUnc: " << residualUnc << endl;
			// cout << "recoMCUnc: " << recoMCUnc << endl;
			// cout << "recoMCVal: " << recoMCVal << endl;

			// cout << "totalError: " << totalError << endl;
			/// set the bin contents reflecting weights
			/// yield with acceptance x efficiency correction
			correctedHist->SetBinContent(iCosTheta + 1, iPhi + 1, recoMCVal * weight);
			// correctedHist->SetBinContent(iCosTheta + 1, iPhi + 1, residualVal);
			
			rawYieldHist1DCosTheta[iPhi]->SetBinContent(iCosTheta + 1, recoMCVal);
			rawYieldHist1DPhi[iCosTheta]->SetBinContent(iPhi + 1, recoMCVal);

			correctedHist1DCosTheta[iPhi]->SetBinContent(iCosTheta + 1, recoMCVal * weight);
			correctedHist1DPhi[iCosTheta]->SetBinContent(iPhi + 1, recoMCVal * weight);

			// standardCorrectedMap->SetBinError(iCosTheta + 1, iPhi + 1, TMath::Hypot(yield1SUnc / yield1SVal, totalRelUncHigh) * yield1SVal * weight);

			// correctedHist->SetBinError(iCosTheta + 1, iPhi + 1, recoMCUnc * weight);
			correctedHist->SetBinError(iCosTheta + 1, iPhi + 1, totalError);
			// correctedHist->SetBinError(iCosTheta + 1, iPhi + 1, TMath::Hypot(recoMCUnc / recoMCVal, totalRelUncHigh));
			// correctedHist->SetBinError(iCosTheta + 1, iPhi + 1, residualUnc);
			
			rawYieldHist1DCosTheta[iPhi]->SetBinError(iCosTheta + 1, recoMCUnc);
			rawYieldHist1DPhi[iCosTheta]->SetBinError(iPhi + 1, recoMCUnc);

			correctedHist1DCosTheta[iPhi]->SetBinError(iCosTheta + 1, totalError);
			correctedHist1DPhi[iCosTheta]->SetBinError(iPhi + 1, totalError);

			cout << "correctedHist->GetBinContent(" << iCosTheta + 1 << ", " << iPhi + 1 << ") = " << correctedHist->GetBinContent(iCosTheta + 1, iPhi + 1) << endl;
			cout << "correctedHist->GetBinError(" << iCosTheta + 1 << ", " << iPhi + 1 << ") = " << correctedHist->GetBinError(iCosTheta + 1, iPhi + 1) << endl;
			// /// fill uncertainty histograms
			// relSystEffCosThetaPhi->SetBinContent(iCosTheta +1, iPhi + 1, relSystUnc);
			// statHighEffCosThetaPhi->SetBinContent(iCosTheta + 1, iPhi + 1, relEffUncHigh);
			// statLowEffCosThetaPhi->SetBinContent(iCosTheta + 1, iPhi + 1, relEffUncLow);
			// statHighAccCosThetaPhi->SetBinContent(iCosTheta + 1, iPhi + 1, relAccUncHigh);
			// statLowAccCosThetaPhi->SetBinContent(iCosTheta + 1, iPhi + 1, relAccUncLow);
			// yield1SUncCosThetaPhi->SetBinContent(iCosTheta + 1, iPhi + 1, yield1SUnc / yield1SVal);
			// totalRelUncCosThetaPhi->SetBinContent(iCosTheta + 1, iPhi + 1, TMath::Hypot(yield1SUnc / yield1SVal, totalRelUncHigh));

			// if ((yield1S->getVal()) * weight > maxYield) maxYield = (yield1S->getVal()) * weight;
		}
	}

	if (drawPlot) {
		TCanvas* correctedCanvas = draw2DMap(correctedHist, refFrameName.Data(), nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kFALSE, 1, kFALSE, nLegendRows);
		display2DMapContents(correctedHist, nCosThetaBins, nPhiBins, kFALSE, 0.035, kWhite, 0);

		/// Styles of the texts in the plot
		TLatex* legend1 = new TLatex();
		legend1->SetTextAlign(22);
		legend1->SetTextSize(0.05);

		/// Put texts inside the plot
		legend1->DrawLatexNDC(.50, .88, "Corrected MC");
		legend1->DrawLatexNDC(.48, .80, Form("Input: #lambda_{#theta} = %.2f, #lambda_{#varphi} = %.2f, #lambda_{#theta#varphi} = %.2f", lambdaTheta0, lambdaPhi0, lambdaThetaPhi0));

		gPad->RedrawAxis();

		gPad->Update();
	
		/// draw 1D histograms along the cosTheta axis with the fit
		TCanvas* rawYield1DCanvas = new TCanvas("rawYield1DCanvas", "rawYield1DCanvas", 1500, 500);
		rawYield1DCanvas->Divide(3, 1, 0.001, 0.001);

		TCanvas* corrected1DCanvas = new TCanvas("corrected1DCanvas", "corrected1DCanvas", 1500, 500);
		corrected1DCanvas->Divide(3, 1, 0.001, 0.001);

		for (Int_t iPhi = 0; iPhi < nPhiBins; iPhi++) {
			/// cosmetics for 1D histograms along the cosTheta axis
			rawYield1DCanvas->cd(iPhi + 1);
			gPad->SetLeftMargin(0.17);

			rawYieldHist1DCosTheta[iPhi]->GetXaxis()->SetTitle(CosThetaVarTitle(refFrameName.Data()));
			rawYieldHist1DCosTheta[iPhi]->GetXaxis()->CenterTitle();
			rawYieldHist1DCosTheta[iPhi]->GetXaxis()->SetNdivisions(-nCosThetaBins);

			rawYieldHist1DCosTheta[iPhi]->GetYaxis()->SetTitle("#varUpsilon(1S) raw events");
			rawYieldHist1DCosTheta[iPhi]->GetYaxis()->SetRangeUser(0, rawYieldHist1DCosTheta[iPhi]->GetMaximum() * 1.5);
			rawYieldHist1DCosTheta[iPhi]->GetYaxis()->SetTitleOffset(1.4);

			rawYieldHist1DCosTheta[iPhi]->SetMarkerStyle(20);
			rawYieldHist1DCosTheta[iPhi]->SetMarkerSize(1.2);
			rawYieldHist1DCosTheta[iPhi]->SetLineWidth(2);
			
			rawYieldHist1DCosTheta[iPhi]->Draw("E");

			/// legend
			TLatex* text1D = new TLatex();
			text1D->SetTextAlign(22);
			text1D->SetTextSize(0.05);
			/// Put texts inside the plot
			text1D->DrawLatexNDC(.56, .86, Form("%s, %d < p_{T} < %d GeV/c", refFrameName.Data(), ptMin, ptMax));
			text1D->DrawLatexNDC(.56, .80, Form("%d < |#varphi| < %d #circ", (int)phiBinEdges[iPhi], (int)phiBinEdges[iPhi + 1]));

			CMS_lumi((TPad*)gPad, gCMSLumiText);	
			
			gPad->RedrawAxis();

			gPad->Update();

			corrected1DCanvas->cd(iPhi + 1);
			gPad->SetLeftMargin(0.17);

			correctedHist1DCosTheta[iPhi]->GetXaxis()->SetTitle(CosThetaVarTitle(refFrameName.Data()));
			correctedHist1DCosTheta[iPhi]->GetXaxis()->CenterTitle();
			correctedHist1DCosTheta[iPhi]->GetXaxis()->SetNdivisions(-nCosThetaBins);

			correctedHist1DCosTheta[iPhi]->GetYaxis()->SetTitle("#varUpsilon(1S) weighted events");
			correctedHist1DCosTheta[iPhi]->GetYaxis()->SetRangeUser(0, correctedHist1DCosTheta[iPhi]->GetMaximum() * 1.5);
			correctedHist1DCosTheta[iPhi]->GetYaxis()->SetTitleOffset(1.4);
			
			correctedHist1DCosTheta[iPhi]->SetMarkerStyle(20);
			correctedHist1DCosTheta[iPhi]->SetMarkerSize(1.2);
			correctedHist1DCosTheta[iPhi]->SetLineWidth(2);

			correctedHist1DCosTheta[iPhi]->Draw("E");

			/// legend
			TLatex* text1Dcorrected = new TLatex();
			text1Dcorrected->SetTextAlign(22);
			text1Dcorrected->SetTextSize(0.05);
			/// Put texts inside the plot
			text1Dcorrected->DrawLatexNDC(.56, .86, Form("%s, %d < p_{T} < %d GeV/c", refFrameName.Data(), ptMin, ptMax));
			text1Dcorrected->DrawLatexNDC(.56, .80, Form("%d < |#varphi| < %d #circ", (int)phiBinEdges[iPhi], (int)phiBinEdges[iPhi + 1]));
		
			CMS_lumi((TPad*)gPad, gCMSLumiText);	
		}

		/// draw 1D histograms along the cosTheta axis with the fit
		TCanvas* rawYield1DPhiCanvas = new TCanvas("rawYield1DPhiCanvas", "rawYield1DPhiCanvas", 1500, 300);
		rawYield1DPhiCanvas->Divide(5, 1, 0.002, 0.001); // last two arguments: Reduce horizontal and vertical spacing

		TCanvas* corrected1DPhiCanvas = new TCanvas("corrected1DPhiCanvas", "corrected1DPhiCanvas", 1500, 300);
		corrected1DPhiCanvas->Divide(5, 1, 0.002, 0.001); // last two arguments: Reduce horizontal and vertical spacing

		for (Int_t iCosTheta = 0; iCosTheta < nCosThetaBins; iCosTheta++) {
			/// cosmetics for 1D histograms along the cosTheta axis
			rawYield1DPhiCanvas->cd(iCosTheta + 1);
			gPad->SetLeftMargin(0.17);

			rawYieldHist1DPhi[iCosTheta]->GetXaxis()->SetTitle(AbsPhiAxisTitle(refFrameName.Data()));
			rawYieldHist1DPhi[iCosTheta]->GetXaxis()->CenterTitle();
			rawYieldHist1DPhi[iCosTheta]->GetXaxis()->SetNdivisions(-nPhiBins);

			rawYieldHist1DPhi[iCosTheta]->GetYaxis()->SetTitle("#varUpsilon(1S) raw events");
			rawYieldHist1DPhi[iCosTheta]->GetYaxis()->SetRangeUser(0, rawYieldHist1DPhi[iCosTheta]->GetMaximum() * 1.5);
			rawYieldHist1DPhi[iCosTheta]->GetYaxis()->SetTitleOffset(1.4);

			rawYieldHist1DPhi[iCosTheta]->SetMarkerStyle(20);
			rawYieldHist1DPhi[iCosTheta]->SetMarkerSize(1.2);
			rawYieldHist1DPhi[iCosTheta]->SetLineWidth(2);
			
			rawYieldHist1DPhi[iCosTheta]->Draw("E");

			/// legend
			TLatex* text1D = new TLatex();
			text1D->SetTextAlign(22);
			text1D->SetTextSize(0.05);
			/// Put texts inside the plot
			text1D->DrawLatexNDC(.56, .84, Form("%s, %d < p_{T} < %d GeV/c", refFrameName.Data(), ptMin, ptMax));
			text1D->DrawLatexNDC(.56, .79, Form("%0.2f < cos#theta < %0.2f", cosThetaBinEdges[iCosTheta], cosThetaBinEdges[iCosTheta + 1]));

			CMS_lumi((TPad*)gPad, gCMSLumiText);	
			
			gPad->RedrawAxis();

			gPad->Update();

			corrected1DPhiCanvas->cd(iCosTheta + 1);
			gPad->SetLeftMargin(0.17);

			correctedHist1DPhi[iCosTheta]->GetXaxis()->SetTitle(AbsPhiAxisTitle(refFrameName.Data()));
			correctedHist1DPhi[iCosTheta]->GetXaxis()->CenterTitle();
			correctedHist1DPhi[iCosTheta]->GetXaxis()->SetNdivisions(-nPhiBins);

			correctedHist1DPhi[iCosTheta]->GetYaxis()->SetTitle("#varUpsilon(1S) weighted events");
			correctedHist1DPhi[iCosTheta]->GetYaxis()->SetRangeUser(0, correctedHist1DPhi[iCosTheta]->GetMaximum() * 1.5);
			correctedHist1DPhi[iCosTheta]->GetYaxis()->SetTitleOffset(1.4);
			
			correctedHist1DPhi[iCosTheta]->SetMarkerStyle(20);
			correctedHist1DPhi[iCosTheta]->SetMarkerSize(1.2);
			correctedHist1DPhi[iCosTheta]->SetLineWidth(2);

			correctedHist1DPhi[iCosTheta]->Draw("E");

			/// legend
			TLatex* text1Dcorrected = new TLatex();
			text1Dcorrected->SetTextAlign(22);
			text1Dcorrected->SetTextSize(0.05);
			/// Put texts inside the plot
			text1Dcorrected->DrawLatexNDC(.56, .84, Form("%s, %d < p_{T} < %d GeV/c", refFrameName.Data(), ptMin, ptMax));
			text1Dcorrected->DrawLatexNDC(.56, .79, Form("%0.2f < cos#theta < %0.2f", cosThetaBinEdges[iCosTheta], cosThetaBinEdges[iCosTheta + 1]));

			CMS_lumi((TPad*)gPad, gCMSLumiText);
			
			gPad->RedrawAxis();
			
			gPad->Update();
		}

		gSystem->mkdir(Form("closureTest/%s", applyEff ? "HydjetAccEff": "PythiaAcc"), kTRUE);
		correctedCanvas->SaveAs(Form("closureTest/%s/2Dcorrected_%s_pt%dto%d_cosTheta%.2fto%.2f_phi%dto%d_lambdaTheta%.2f_phi%.2f_ThetaPhi%.2f.pdf", applyEff ? "HydjetAccEff": "PythiaAcc", refFrameName.Data(), ptMin, ptMax, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], (int)phiBinEdges[0], (int)phiBinEdges[nPhiBins], lambdaTheta, lambdaPhi, lambdaThetaPhi));
		correctedCanvas->SaveAs(Form("closureTest/%s/2Dcorrected_%s_pt%dto%d_cosTheta%.2fto%.2f_phi%dto%d_lambdaTheta%.2f_phi%.2f_ThetaPhi%.2f.png", applyEff ? "HydjetAccEff": "PythiaAcc", refFrameName.Data(), ptMin, ptMax, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], (int)phiBinEdges[0], (int)phiBinEdges[nPhiBins], lambdaTheta, lambdaPhi, lambdaThetaPhi));
	}

	gPad->RedrawAxis();

	/// save 2D acc and eff plots
	gSystem->mkdir(Form("closureTest/%s", applyEff ? "HydjetAccEff": "PythiaAcc"), kTRUE);
	if (isPhiFolded) {
		/// png
		accCanvas->SaveAs(Form("closureTest/%s/2Dacc_%s_pt%dto%d_cosTheta%.2fto%.2f_absphi%dto%d_lambdaTheta%.2f_phi%.2f_ThetaPhi%.2f.png", applyEff ? "HydjetAccEff": "PythiaAcc", refFrameName.Data(), ptMin, ptMax, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], (int)phiBinEdges[0], (int)phiBinEdges[nPhiBins], lambdaTheta, lambdaPhi, lambdaThetaPhi));
		effCanvas->SaveAs(Form("closureTest/%s/2Deff_%s_pt%dto%d_cosTheta%.2fto%.2f_absphi%dto%d_lambdaTheta%.2f_phi%.2f_ThetaPhi%.2f.png", applyEff ? "HydjetAccEff": "PythiaAcc", refFrameName.Data(), ptMin, ptMax, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], (int)phiBinEdges[0], (int)phiBinEdges[nPhiBins], lambdaTheta, lambdaPhi, lambdaThetaPhi));
	} else {
		/// pdf
		accCanvas->SaveAs(Form("closureTest/%s/2Dacc_%s_pt%dto%d_cosTheta%.2fto%.2f_phi%dto%d_lambdaTheta%.2f_phi%.2f_ThetaPhi%.2f.pdf", applyEff ? "HydjetAccEff": "PythiaAcc", refFrameName.Data(), ptMin, ptMax, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], (int)phiBinEdges[0], (int)phiBinEdges[nPhiBins], lambdaTheta, lambdaPhi, lambdaThetaPhi));
		effCanvas->SaveAs(Form("closureTest/%s/2Deff_%s_pt%dto%d_cosTheta%.2fto%.2f_phi%dto%d_lambdaTheta%.2f_phi%.2f_ThetaPhi%.2f.pdf", applyEff ? "HydjetAccEff": "PythiaAcc", refFrameName.Data(), ptMin, ptMax, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], (int)phiBinEdges[0], (int)phiBinEdges[nPhiBins], lambdaTheta, lambdaPhi, lambdaThetaPhi));
		///png
		accCanvas->SaveAs(Form("closureTest/%s/2Dacc_%s_pt%dto%d_cosTheta%.2fto%.2f_phi%dto%d_lambdaTheta%.2f_phi%.2f_ThetaPhi%.2f.png", applyEff ? "HydjetAccEff": "PythiaAcc", refFrameName.Data(), ptMin, ptMax, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], (int)phiBinEdges[0], (int)phiBinEdges[nPhiBins], lambdaTheta, lambdaPhi, lambdaThetaPhi));
		effCanvas->SaveAs(Form("closureTest/%s/2Deff_%s_pt%dto%d_cosTheta%.2fto%.2f_phi%dto%d_lambdaTheta%.2f_phi%.2f_ThetaPhi%.2f.png", applyEff ? "HydjetAccEff": "PythiaAcc", refFrameName.Data(), ptMin, ptMax, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], (int)phiBinEdges[0], (int)phiBinEdges[nPhiBins], lambdaTheta, lambdaPhi, lambdaThetaPhi));
	}

	return;
}

void getYieldHist(TH2D* angDistHist2D, TString refFrameName = "CS", TString muonAccName = "UpsilonTriggerThresholds",
                  Float_t lambdaTheta = 0, Float_t lambdaPhi = 0, Float_t lambdaThetaPhi = 0,
                  Int_t ptMin = 0, Int_t ptMax = 30,
                  Int_t nCosThetaBins = 5, const vector<Double_t>& cosThetaBinEdges = {},
                  Int_t nPhiBins = 6, const vector<Double_t>& phiBinEdges = {},
                  Bool_t isPhiFolded = kFALSE, const char* signalShapeName = "SymDSCB", const char* defaultBkgShapeName = "ExpTimesErr", int nLegendRows = 1) {
	writeExtraText = true; // if extra text
	extraText = "       Internal";
	// extraText = "       Simulation Preliminary";

	// bin width
	Double_t cosThetaStep = (cosThetaBinEdges[nCosThetaBins] - cosThetaBinEdges[0]) / nCosThetaBins;

	Double_t phiStep = (phiBinEdges[nPhiBins] - phiBinEdges[0]) / nPhiBins;

	RooRealVar* yield1S = new RooRealVar("yield1S", "", 1000);
	RooRealVar* yield2S = new RooRealVar("yield2S", "", 100);
	RooRealVar* yield3S = new RooRealVar("yield3S", "", 10);

	/// Assign signal and background shape name to read the file for the yield extraction results
	// const char* signalShapeName = "SymDSCB";
	// const char* signalShapeName = "Johnson";

	// background shape array: ChebychevOrderN or ExpTimesErr

	const Int_t nCosThetaBinsMax = 20;
	const Int_t nPhiBinsMax = 10;

	std::string bkgShapeName[nCosThetaBinsMax][nPhiBinsMax];

	// // fill the background shape array with ChebychevOrder2
	if (strcmp(defaultBkgShapeName, "Chebychev") == 0) {
		std::fill(&bkgShapeName[0][0], &bkgShapeName[0][0] + nCosThetaBinsMax * nPhiBinsMax, "ChebychevOrder2");

		// // exceptions
		// // bkgShapeName[1][1] = "ChebychevOrder1";
		// // bkgShapeName[1][3] = "ChebychevOrder1";
		// // bkgShapeName[1][4] = "ChebychevOrder1";
		// // bkgShapeName[1][5] = "ChebychevOrder1";
	}

	if (strcmp(defaultBkgShapeName, "ExpTimesErr") == 0) {
		// fill the background shape array with ExpTimesErr
		std::fill(&bkgShapeName[0][0], &bkgShapeName[0][0] + nCosThetaBinsMax * nPhiBinsMax, "ExpTimesErr");

		// // exceptions
		// // HX, 2 < pT < 6 GeV/c (-180 to 180)
		// bkgShapeName[0][2] = "ChebychevOrder1";
		// bkgShapeName[0][3] = "ChebychevOrder1";

		// bkgShapeName[3][2] = "ChebychevOrder2";
		// bkgShapeName[3][3] = "ChebychevOrder2";

		// bkgShapeName[3][0] = "ChebychevOrder2";
		// bkgShapeName[3][5] = "ChebychevOrder2";
	}

	Float_t maxYield = 0;

	TCanvas* massCanvas = 0;

	for (Int_t iCosTheta = 0; iCosTheta < nCosThetaBins; iCosTheta++) {
        for (Int_t iPhi = 0; iPhi < nPhiBins; iPhi++) {
            
			/// get bin center values
			Double_t binCenterCosTheta = angDistHist2D->GetXaxis()->GetBinCenter(iCosTheta + 1);
			Double_t binCenterPhi = angDistHist2D->GetYaxis()->GetBinCenter(iPhi + 1);

			/// Get yields and their uncertainties
			Int_t absiPhi = angDistHist2D->GetYaxis()->FindBin(fabs(binCenterPhi)) - 1;
			cout << "absiPhi: " << absiPhi << endl;

			const char* fitModelName = GetFitModelName(signalShapeName, ptMin, ptMax, refFrameName.Data(), cosThetaBinEdges[iCosTheta], cosThetaBinEdges[iCosTheta + 1], (Int_t)phiBinEdges[absiPhi], (Int_t)phiBinEdges[absiPhi + 1]);

			RooArgSet signalYields = GetSignalYields(yield1S, yield2S, yield3S, Form("%s", bkgShapeName[iCosTheta][absiPhi].c_str()), fitModelName, Form("%s", muonAccName.Data()));

			cout << "signalYields: " << signalYields << endl;

			/// get the yield value
			double yield1SVal = (yield1S->getVal());

			double yield1SUnc = (yield1S->getError());

			/// Set the bin contents reflecting weights
			/// only raw yield itself before correction
			angDistHist2D->SetBinContent(iCosTheta + 1, iPhi + 1, yield1SVal);

			/// Set the bin error of the raw yield
			// yieldMap->SetBinError(iCosTheta + 1, yield1SErr * weight);
			// yieldMap->SetBinError(iCosTheta + 1, yield1SVal * weight * TMath::Hypot(yield1SErr / yield1SVal, relAccUncHigh));
			// yieldMap->SetBinError(iCosTheta + 1, yield1SVal * weight * TMath::Hypot(yield1SErr / yield1SVal, relEffUncHigh));

			angDistHist2D->SetBinError(iCosTheta + 1, iPhi + 1, yield1SUnc);

			// if ((yield1S->getVal()) * weight > maxYield) maxYield = (yield1S->getVal()) * weight / 2.;
			// if ((yield1S->getVal()) * weight > maxYield) maxYield = (yield1S->getVal()) * weight;
        }
    }

    /// draw signal yield histogram
	TCanvas* yieldCanvas = draw2DMap(angDistHist2D, refFrameName.Data(), nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kFALSE, 1, isPhiFolded, nLegendRows);
	display2DMapContents(angDistHist2D, nCosThetaBins, nPhiBins, kFALSE);
	
	/// Styles of the texts in the plot
	TLatex* legend1 = new TLatex();
	legend1->SetTextAlign(22);
	legend1->SetTextSize(0.05);

	/// Put texts inside the plot
	legend1->DrawLatexNDC(.50, .88, "Reconstructed MC");
	legend1->DrawLatexNDC(.50, .80, Form("#lambda_{#theta} = %.2f, #lambda_{#varphi} = %.2f, #lambda_{#theta#varphi} = %.2f", lambdaTheta, lambdaPhi, lambdaThetaPhi));

	gPad->RedrawAxis();

	gPad->Update();

	// save the plot
	gSystem->mkdir("YieldMaps/", kTRUE);
    if (isPhiFolded)
	    yieldCanvas->SaveAs(Form("YieldMaps/Ups1SYield_%s_pt%dto%d_cosTheta%.2fto%.2f_absphi%dto%d_LambdaTheta%.2f_Phi%.2f_ThetaPhi%.2f.png", refFrameName.Data(), ptMin, ptMax, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], (int)phiBinEdges[0], (int)phiBinEdges[nPhiBins], lambdaTheta, lambdaPhi, lambdaThetaPhi), "RECREATE");
    else
       yieldCanvas->SaveAs(Form("YieldMaps/Ups1SYield_%s_pt%dto%d_cosTheta%.2fto%.2f_phi%dto%d_LambdaTheta%.2f_Phi%.2f_ThetaPhi%.2f.png", refFrameName.Data(), ptMin, ptMax, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], (int)phiBinEdges[0], (int)phiBinEdges[nPhiBins], lambdaTheta, lambdaPhi, lambdaThetaPhi), "RECREATE");
	
    return;
}

void rawYield_2D_customizedFits_iteration(TString refFrameName = "CS", TString muonAccName = "UpsilonTriggerThresholds",
                                          Int_t ptMin = 2, Int_t ptMax = 6,
                                          const Int_t nCosThetaBins = 5, Double_t cosThetaMin = -0.7, Double_t cosThetaMax = 0.7,
                                          const Int_t nPhiBins = 6, Int_t phiMin = -180, Int_t phiMax = 180,
                                          Bool_t isPhiFolded = kFALSE, Bool_t applyAcc = kTRUE, Bool_t applyEff = kFALSE, int totNItrs = 0, int nLegendRows = 1) {

	/// set bin edges and width
	vector<Double_t> cosThetaBinEdges = setCosThetaBinEdges(nCosThetaBins, cosThetaMin, cosThetaMax);

	vector<Double_t> phiBinEdges = setPhiBinEdges(nPhiBins, phiMin, phiMax);

	/// bin width
	Double_t cosThetaStep = (cosThetaBinEdges[nCosThetaBins] - cosThetaBinEdges[0]) / nCosThetaBins;

	Double_t phiStep = (phiBinEdges[nPhiBins] - phiBinEdges[0]) / nPhiBins;

	/// create 2D angular distribution histogram
	TH2D* angDistHist2D = new TH2D("angDistHist2D", "; cos #theta; #varphi (#circ); Number of generated #varUpsilon(1S) events", nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);
	
	/// mean nothing, just for the file name
	Float_t lambdaTheta0 = 0, lambdaPhi0 = 0, lambdaThetaPhi0 = 0;

	/// define polarization parameters that will be used for correction (acceptance and efficiency)
	RooRealVar* lambdaTheta = new RooRealVar("lambdaTheta", "lambdaTheta", -3., 3.);
	RooRealVar* lambdaPhi = new RooRealVar("lambdaPhi", "lambdaPhi", -3., 3.);
	RooRealVar* lambdaThetaPhi = new RooRealVar("lambdaThetaPhi", "lambdaThetaPhi", -3., 3.);
	RooRealVar* lambdaTilde = new RooRealVar("lambdaTilde", "lambdaTilde", -3., 3.);

	/// always start with null polarization assumption in acceptance and efficiency
	lambdaTheta->setVal(0);
	lambdaPhi->setVal(0);
	lambdaThetaPhi->setVal(0);
	lambdaTilde->setVal(0);

	/// place holder of the polarization parameters before and after iterations
	vector<double> lambdaThetaArr, lambdaPhiArr, lambdaThetaPhiArr, lambdaTildeArr, numItrArr;
	vector<double> lambdaThetaErrArr, lambdaPhiErrArr, lambdaThetaPhiErrArr, lambdaTildeErrArr, xErrArr;

	/// set the initial values of the polarization parameters for acceptance and efficiency correction
	lambdaThetaArr.push_back(lambdaTheta->getVal());
	lambdaPhiArr.push_back(lambdaPhi->getVal());
	lambdaThetaPhiArr.push_back(lambdaThetaPhi->getVal());
	lambdaTildeArr.push_back(lambdaTilde->getVal());
	numItrArr.push_back(-1);

	lambdaThetaErrArr.push_back(lambdaTheta->getError());
	lambdaPhiErrArr.push_back(lambdaPhi->getError());
	lambdaThetaPhiErrArr.push_back(lambdaThetaPhi->getError());
	lambdaTildeErrArr.push_back(lambdaTilde->getError());
	
	xErrArr.push_back(0.);

	/// store lambda parameters
	RooArgSet polarParams(*lambdaTheta, *lambdaPhi, *lambdaThetaPhi, *lambdaTilde);

    getYieldHist(angDistHist2D, refFrameName, muonAccName, lambdaTheta0, lambdaPhi0, lambdaThetaPhi0, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, isPhiFolded, "SymDSCB", "ExpTimesErr", nLegendRows);

	/// define histograms
	TH2D* correctedHist = new TH2D("correctedHist", "; cos #theta; #varphi (#circ); Number of generated #varUpsilon(1S) events", nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);
	std::vector<TH1D*> correctedHist1DCosTheta(nPhiBins);
	std::vector<TH1D*> correctedHist1DPhi(nCosThetaBins);

	std::vector<TH1D*> rawYieldHist1DCosTheta(nPhiBins);
	std::vector<TH1D*> rawYieldHist1DPhi(nCosThetaBins);
	
	TH2D* fittedHist = new TH2D("fittedHist", "; cos #theta; #varphi (#circ); Number of generated #varUpsilon(1S) events", nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);

	/// loop over the number of iterations
	int nItrs = 0; /// number of iterations

	for (int iItr = 0; iItr <= totNItrs; iItr++) {
		/// correct the MC
		correctMC2DHist(angDistHist2D, correctedHist, rawYieldHist1DCosTheta, rawYieldHist1DPhi, correctedHist1DCosTheta, correctedHist1DPhi, refFrameName, polarParams, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kTRUE, isPhiFolded, applyAcc, applyEff, lambdaTheta0, lambdaPhi0, lambdaThetaPhi0, nLegendRows);

		/// due to drawing issue, separate correctedHist and fittedHist
		correctMC2DHist(angDistHist2D, fittedHist, rawYieldHist1DCosTheta, rawYieldHist1DPhi, correctedHist1DCosTheta, correctedHist1DPhi, refFrameName, polarParams, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, isPhiFolded, applyAcc, applyEff, lambdaTheta0, lambdaPhi0, lambdaThetaPhi0, nLegendRows);

		/// extract parameters from the corrected hist
		polarParams = extractPolarParam(fittedHist, correctedHist1DCosTheta, refFrameName, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, isPhiFolded, applyEff);
		cout << "extracted polarization parameters: " << endl;
		polarParams.Print("v");

		/// update the polarization parameters
		lambdaThetaArr.push_back(lambdaTheta->getVal());
		lambdaPhiArr.push_back(lambdaPhi->getVal());
		lambdaThetaPhiArr.push_back(lambdaThetaPhi->getVal());
		lambdaTildeArr.push_back(lambdaTilde->getVal());

		numItrArr.push_back(iItr); // x-axis values for the graph

		/// update the errors of the polarization parameters
		lambdaThetaErrArr.push_back(lambdaTheta->getError());
		lambdaPhiErrArr.push_back(lambdaPhi->getError());
		lambdaThetaPhiErrArr.push_back(lambdaThetaPhi->getError());
		lambdaTildeErrArr.push_back(lambdaTilde->getError());
	
		xErrArr.push_back(0.); // x-axis error for the graph

		// cout << "lambdaTheta arr size: " << lambdaThetaArr.size() << endl;
		// cout << lambdaThetaArr[0] << ", " << lambdaThetaArr[1] << endl;
		// cout << "lambdaPhi arr size: " << lambdaPhiArr.size() << endl;
		// cout << lambdaPhiArr[0] << ", " << lambdaPhiArr[1] << endl;
		// cout << "lambdaThetaPhi arr size: " << lambdaThetaPhiArr.size() << endl;
		// cout << lambdaThetaPhiArr[0] << ", " << lambdaThetaPhiArr[1] << endl;
		// cout << "lambdaTilde arr size: " << lambdaTildeArr.size() << endl;
		// cout << lambdaTildeArr[0] << ", " << lambdaTildeArr[1] << endl;
		// cout << "numItr arr size: " << numItrArr.size() << endl;
		// cout << numItrArr[0] << ", " << numItrArr[1] << endl;

		nItrs++;
		
		/// if the difference between the new paramter and the previous one is less than 0.01, stop the iterative procedure
		if ((fabs(lambdaThetaArr[iItr] - lambdaThetaArr[iItr + 1]) < 0.01) && (fabs(lambdaPhiArr[iItr] - lambdaPhiArr[iItr + 1]) < 0.01) && (fabs(lambdaThetaPhiArr[iItr] - lambdaThetaPhiArr[iItr + 1]) < 0.01)) {
			// cout << "converged after " << iItr << " iterations" << endl;
			// cout << "lambdaTheta diff: " << lambdaThetaArr[iItr] << ", " << lambdaThetaArr[iItr + 1] << ", " << fabs(lambdaThetaArr[iItr + 1] - lambdaThetaArr[iItr]) << endl;
			// cout << "lambdaPhi diff: " << lambdaPhiArr[iItr] << ", " << fabs(lambdaPhiArr[iItr + 1] - lambdaPhiArr[iItr]) << endl;
			// cout << "lambdaThetaPhi diff: " << lambdaThetaPhiArr[iItr] << ", " << fabs(lambdaThetaPhiArr[iItr + 1] - lambdaThetaPhiArr[iItr]) << endl;

			// break;
		}
	}

	/// save the last iteration results
	const char* fitModelName = GetFitModelName("SymDSCB", ptMin, ptMax, refFrameName.Data(), cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], (Int_t)phiBinEdges[0], (Int_t)phiBinEdges[nPhiBins]);

	SavePolarizationFitParameters(&polarParams, "iteration", fitModelName, muonAccName);

	/// draw plots for polarization parameters vs # of Iterations
	TCanvas* iterCanvas = new TCanvas("iterCanvas", "iterCanvas", 650, 600);
	iterCanvas->SetRightMargin(0.05);
	iterCanvas->cd();

	TGraphErrors* lambdaThetaGraph = new TGraphErrors(lambdaThetaArr.size(), numItrArr.data(), lambdaThetaArr.data(), xErrArr.data(), lambdaThetaErrArr.data());
	lambdaThetaGraph->SetTitle("; Number of Iterations; Polarization Parameter Values");
	lambdaThetaGraph->GetXaxis()->CenterTitle();
	lambdaThetaGraph->GetYaxis()->CenterTitle();
	lambdaThetaGraph->SetMarkerStyle(20);
	lambdaThetaGraph->SetMarkerColor((TColor::GetColor("#FF1F5B")));
	lambdaThetaGraph->SetLineColor((TColor::GetColor("#FF1F5B")));
	lambdaThetaGraph->SetLineWidth(2);
	lambdaThetaGraph->SetMarkerSize(1);
	lambdaThetaGraph->SetMinimum(-1.2);
	lambdaThetaGraph->SetMaximum(1.5);

	gPad->Update(); // ensures the histogram behind the graph is created

	if ((nItrs - 1) == 0) lambdaThetaGraph->GetXaxis()->SetNdivisions(nItrs + 1);
	else lambdaThetaGraph->GetXaxis()->SetNdivisions(nItrs + 2);
	lambdaThetaGraph->Draw("APL");

	TGraphErrors* lambdaPhiGraph = new TGraphErrors(lambdaPhiArr.size(), numItrArr.data(), lambdaPhiArr.data(), xErrArr.data(), lambdaPhiErrArr.data());
	lambdaPhiGraph->SetMarkerStyle(20);
	lambdaPhiGraph->SetMarkerColor(TColor::GetColor("#009ADE"));
	lambdaPhiGraph->SetLineColor(TColor::GetColor("#009ADE"));
	lambdaPhiGraph->SetLineWidth(2);
	lambdaPhiGraph->SetMarkerSize(1);
	lambdaPhiGraph->Draw("PL same");

	TGraphErrors* lambdaThetaPhiGraph = new TGraphErrors(lambdaThetaPhiArr.size(), numItrArr.data(), lambdaThetaPhiArr.data(), xErrArr.data(), lambdaThetaPhiErrArr.data());
	lambdaThetaPhiGraph->SetMarkerStyle(20);
	lambdaThetaPhiGraph->SetMarkerColor(TColor::GetColor("#bababa"));
	lambdaThetaPhiGraph->SetLineColor(TColor::GetColor("#bababa"));
	lambdaThetaPhiGraph->SetLineWidth(2);
	lambdaThetaPhiGraph->SetMarkerSize(1);
	lambdaThetaPhiGraph->Draw("PL same");

	if ((nItrs - 1) == 0) {
		drawLine(-1.1, 0, 0, 0);
		drawLine(-1.1, lambdaTheta0, 0, lambdaTheta0);
		drawLine(-1.1, lambdaPhi0, 0, lambdaPhi0);
		drawLine(-1.1, lambdaThetaPhi0, 0, lambdaThetaPhi0);
	}

    else if ((nItrs - 1) == 2) {
		cout << "nItrs: " << nItrs << endl;
		drawLine(-1.3, 0, 2.3, 0);
		drawLine(-1.3, lambdaTheta0, 2.3, lambdaTheta0);
		drawLine(-1.3, lambdaPhi0, 2.3, lambdaPhi0);
		drawLine(-1.3, lambdaThetaPhi0, 2.3, lambdaThetaPhi0);
    }

	else {
		cout << "nItrs: " << nItrs << endl;
		drawLine(-1.575, 0, (nItrs - 1) * 1.12, 0);
		drawLine(-1.575, lambdaTheta0, (nItrs - 1) * 1.12, lambdaTheta0);
		drawLine(-1.575, lambdaPhi0, (nItrs - 1) * 1.12, lambdaPhi0);
		drawLine(-1.575, lambdaThetaPhi0, (nItrs - 1) * 1.12, lambdaThetaPhi0);
	}

	TLegend* legend = new TLegend(0.16, 0.78, 0.46, 0.91);
	legend->SetBorderSize(0);
	legend->SetFillColor(0);
	legend->SetTextSize(0.05);

	legend->AddEntry(lambdaThetaGraph, "#lambda_{#theta}", "p");
	legend->AddEntry(lambdaPhiGraph, "#lambda_{#varphi}", "p");
	legend->AddEntry(lambdaThetaPhiGraph, "#lambda_{#theta#varphi}", "p");
	legend->Draw();

	TLatex* legend1 = new TLatex();
	legend1->SetTextAlign(23);
	legend1->SetTextSize(0.04);
	// legend1->DrawLatexNDC(.65, .84, "data");
	// legend1->DrawLatexNDC(.65, .84, Form("Input: (#lambda_{#theta},  #lambda_{#varphi}, #lambda_{#theta#varphi}) = (%.2f, %.2f, %.2f)", lambdaTheta0, lambdaPhi0, lambdaThetaPhi0));

	Bool_t isCSframe;

	if (refFrameName == (TString) "CS")
		isCSframe = true;
	else
		isCSframe = false;

	TPaveText* text = new TPaveText(0.38, 0.85, 0.95, 0.90, "NDCNB");
	text->SetFillColor(4000);
	text->SetBorderSize(0);
	text->AddText(Form("%s, %s", isCSframe ? "Collins-Soper frame" : "Helicity frame", DimuonPtRangeText(ptMin, ptMax)));
	// text->AddText(DimuonPtRangeText(ptMin, ptMax));
	text->SetAllWith("", "align", 32);
	text->Draw("same");

	gPad->Update();   // important to get pad and axis ranges

	// Get y-axis range (used for placing the box slightly below the axis)
	double yMin = gPad->GetUymin();
	double yMax = gPad->GetUymax();
	double yRange = yMax - yMin;
	
	// Define box size and position around x = -1
	double box_x1 = -1.21;
	double box_x2 = -0.81;
	double box_y1 = yMin - 0.06 * yRange - 0.1;  // slightly below visible axis
	double box_y2 = yMin + 0.01 * yRange;  // near the bottom of plot
	
	// Draw box to cover the default label
	TBox* cover = new TBox(box_x1, box_y1, box_x2, box_y2);
	cover->SetFillColor(kWhite);
	cover->SetLineColor(kWhite);
	cover->SetFillStyle(1001);   // solid fill
	cover->Draw("same");
	// cout << box_x1 << ", " << box_x2 << ", " << box_y1 << ", " << box_y2 << endl;

	// Add custom label at x = -1
	TLatex* t = new TLatex(-1, -1.37, "Initial");
	t->SetTextAlign(21);  // center align
	t->SetTextSize(0.047);
	t->Draw("same");

	// TPaveText* refFrameText = RefFrameText(isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax, 0.61, 0.34, 0.97, 0.56);
	// TPaveText* kinematicsText = KinematicsText(gCentralityBinMin, gCentralityBinMax, ptMin,  ptMax, 0.57, 0.95, 0.95, 0.65);

	// TPaveText* refFrameText = new TPaveText( 0.61, 0.34, 0.97, 0.56, "NDCNB");
	// refFrameText->SetFillColor(4000);
	// refFrameText->SetBorderSize(0);
	// refFrameText->AddText(isCSframe ? "Collins-Soper frame" : "Helicity frame");

	// refFrameText->Draw("same");

	// TPaveText* kinematicsText = new TPaveText(0.57, 0.95, 0.95, 0.65, "NDCNB");
	// kinematicsText->SetFillColor(4000);
	// kinematicsText->SetBorderSize(0);
	// // kinematicsText->AddText(CentralityRangeText(gCentralityBinMin, gCentralityBinMax));
	// // kinematicsText->AddText(gMuonPtCutText);
	// // kinematicsText->AddText(DimuonRapidityRangeText(gRapidityMin, gRapidityMax));
	// kinematicsText->AddText(DimuonPtRangeText(ptMin, ptMax));

	// kinematicsText->SetAllWith("", "align", 32);
	// kinematicsText->Draw("same");
	
	CMS_lumi(iterCanvas, gCMSLumiText);	
	// CMS_lumi(iterCanvas, "#varUpsilon(1S) Pythia 8 (5.02 TeV)");

	gPad->RedrawAxis();

	iterCanvas->SaveAs(Form("ClosureTest/%s/Iterations_totItr%d_%s_pt%dto%d_cosTheta%.2fto%.2f_phi%dto%d_LambdaTheta%.2f_Phi%.2f_ThetaPhi%.2f.png", applyEff ? "HydjetAccEff": "PythiaAcc", totNItrs, refFrameName.Data(), ptMin, ptMax, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], (int)phiBinEdges[0], (int)phiBinEdges[nPhiBins], lambdaTheta0, lambdaPhi0, lambdaThetaPhi0), "RECREATE");

	return;
}