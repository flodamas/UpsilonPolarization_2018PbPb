#include "../Tools/Style/tdrStyle.C"
#include "../Tools/Style/CMS_lumi.C"

#include "../AnalysisParameters.h"

// extraction of the polarization paramaters from the fit of the (cos theta, phi) yield distribution

void fit2Dmap(Int_t ptMin = 0, Int_t ptMax = 30, Bool_t isCSframe = kTRUE) {
	// get the raw yield distribution
	TFile* rawYieldFile = TFile::Open("../SignalExtraction/YieldResults/NominalFit_2DAnalysis.root", "READ");
	if (!rawYieldFile) {
		cout << "Raw yield results file not found. Check the directory of the file." << endl;
		return;
	}
	auto* RawYieldMap = (TH2D*)rawYieldFile->Get(Form("RawYield1S_%s_pt%dto%d", (isCSframe) ? "CS" : "HX", ptMin, ptMax));

	// get the acceptance map
	TFile* acceptanceFile = TFile::Open("../MonteCarlo/AcceptanceMaps/1S/AcceptanceResults.root", "READ");
	if (!acceptanceFile) {
		cout << "Acceptance file not found. Check the directory of the file." << endl;
		return;
	}
	auto* AccMap = (TEfficiency*)acceptanceFile->Get(Form("Analysis%s_pt%dto%d", (isCSframe) ? "CS" : "HX", ptMin, ptMax));

	// get the efficiency map
	TFile* efficiencyFile = TFile::Open("../MonteCarlo/EfficiencyMaps/1S/EfficiencyResults.root", "READ");
	if (!efficiencyFile) {
		cout << "Efficiency file not found. Check the directory of the file." << endl;
		return;
	}
	auto* NomEffMap = (TEfficiency*)efficiencyFile->Get(Form("NominalEff_%s_pt%dto%d", (isCSframe) ? "CS" : "HX", ptMin, ptMax));

	TH2D* NormYieldMap = new TH2D("NormYieldMap", ";cos #theta_{CS}; #varphi_{CS} (#circ);d^{2}N_{#varUpsilon(1S)} / d cos #theta d#varphi", NCosThetaBinsCS, CosThetaBinningCS, NPhiBinsCS, PhiBinningCS);

	// WATCH OUT: the index of an histogram starts from 1 while the binning defined in AnalysisParameters.h starts from 0
	// see https://root.cern.ch/doc/master/classTH1.html#a641262682144d465d7e2bc6101a04bf6 for the numbering convention

	for (int iCosTheta = 0; iCosTheta < NCosThetaBinsCS; iCosTheta++) {
		for (int iPhi = 0; iPhi < NPhiBinsCS; iPhi++) {
			int globalBin = RawYieldMap->GetBin(iCosTheta + 1, iPhi + 1); // so that we make sure there is no conflict between the bin numberings

			double normYield = 0;

			double rawYield = RawYieldMap->GetBinContent(globalBin);
			double accValue = AccMap->GetEfficiency(globalBin);
			double effValue = NomEffMap->GetEfficiency(globalBin);

			// make sure the acceptance AND the efficiency are defined for this bin
			if (accValue == 0 || effValue == 0)
				normYield = 0;
			else
				normYield = rawYield / (accValue * effValue * fabs(CosThetaBinningCS[iCosTheta] - CosThetaBinningCS[iCosTheta + 1]) * fabs(PhiBinningCS[iPhi] - PhiBinningCS[iPhi + 1]));

			NormYieldMap->SetBinContent(globalBin, normYield);

			NormYieldMap->SetBinError(globalBin, RawYieldMap->GetBinError(globalBin) * normYield / rawYield); // scale the raw yield error
		}
	}

	// draw the corrected yield distribution

	gStyle->SetPadRightMargin(0.18);
	gStyle->SetPalette(kRainBow);
	gStyle->SetNumberContours(256);

	auto* canvas = new TCanvas("canvas", "", 700, 600);

	NormYieldMap->Draw("COLZ");
	NormYieldMap->GetZaxis()->SetMaxDigits(3);
	NormYieldMap->GetXaxis()->CenterTitle();
	NormYieldMap->GetXaxis()->SetRangeUser(-1, 1);
	NormYieldMap->GetYaxis()->CenterTitle();
	NormYieldMap->GetYaxis()->SetRangeUser(-190, 240);

	const char* legendText = Form("%d < p_{T}^{#mu#mu} < %d GeV, Collins-Soper frame", ptMin, ptMax);

	TLatex* legend = new TLatex(.5, .88, legendText);
	legend->SetTextAlign(22);
	//legend->SetTextSize(0.045);
	legend->DrawLatexNDC(.5, .88, legendText);

	CMS_lumi(canvas, "2018 PbPb data");

	/* time to fit the distribution!!
	 very basic for now, should switch to something more fancy with RooFit later...

	 parameter 0: normalization factor
	 parameter 1: lambda theta
	 parameter 2: lambda phi
	 parameter 3: lambda theta-phi

	 variable x[0]: cos theta (x-axis of the fitted histogram)
	 variable x[1]: phi (y axis) -> to be converted to radians!!!
*/
	TF2* fitFunc = new TF2("fitFunc", "([0]/(3+[1])) * (1 + [1]*x[0]*x[0] + [2]*(1-x[0]*x[0])*(cos(2*TMath::Pi()*x[1]/180)) + [3]*sin(2*acos(x[0]))*cos(TMath::Pi()*x[1]/180))", CosThetaBinningCS[0], CosThetaBinningCS[NCosThetaBinsCS], PhiBinningCS[0], PhiBinningCS[NPhiBinsCS]); // function only defined in the region of the measurement

	fitFunc->SetParName(0, "normalization factor");
	fitFunc->SetParName(1, "lambda theta");
	fitFunc->SetParName(2, "lambda phi");
	fitFunc->SetParName(3, "lambda theta-phi");

	fitFunc->SetParameters(40, 0, 0, 0); // initialize parameters values
	fitFunc->SetParLimits(1, -1.5, 1.5);
	fitFunc->SetParLimits(2, -1.5, 1.5);
	fitFunc->SetParLimits(3, -1., 1.);

	auto fitResult = NormYieldMap->Fit(fitFunc, "QIEMSRN"); // "R" option: only fit in the defined region
	fitResult->Print("v");
	//fitFunc->Draw("SAME");

	TLatex* fitLegend = new TLatex(.5, .88, " ");
	fitLegend->SetTextAlign(22);
	fitLegend->SetTextSize(0.042);
	fitLegend->DrawLatexNDC(.5, .19, Form("2D fit result: #lambda_{#theta} = %.2f #pm %.2f, #lambda_{#varphi} = %.2f #pm %.2f", fitFunc->GetParameter(1), fitFunc->GetParError(1), fitFunc->GetParameter(2), fitFunc->GetParError(2)));

	canvas->SaveAs(Form("DistributionFits/2D/%s_pt%dto%d.png", (isCSframe) ? "CS" : "HX", ptMin, ptMax), "RECREATE");
}
