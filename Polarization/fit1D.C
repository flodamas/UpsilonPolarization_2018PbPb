#include "../AnalysisParameters.h"

// extraction of the polarization paramaters from the fit of the (cos theta, phi) yield distribution

void fit1D(Int_t ptMin = 16, Int_t ptMax = 30, Bool_t isCSframe = kTRUE) {
	// get the raw yield distribution
	TFile* rawYieldFile = TFile::Open("../SignalExtraction/YieldResults/NominalFit_1DAnalysis.root", "READ");
	if (!rawYieldFile) {
		cout << "Raw yield results file not found. Check the directory of the file." << endl;
		return;
	}
	auto* RawYieldMap = (TH1D*)rawYieldFile->Get(Form("RawYield1S_%s_pt%dto%d", (isCSframe) ? "CS" : "HX", ptMin, ptMax));

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

	TH1D* NormYieldMap = new TH1D("NormYieldMap", ";cos #theta_{CS};dN_{#varUpsilon(1S)} / d cos #theta", NCosThetaBinsCS, CosThetaBinningCS);

	// WATCH OUT: the index of an histogram starts from 1 while the binning defined in AnalysisParameters.h starts from 0
	// see https://root.cern.ch/doc/master/classTH1.html#a641262682144d465d7e2bc6101a04bf6 for the numbering convention

	for (int iCosTheta = 0; iCosTheta < NCosThetaBinsCS; iCosTheta++) {
		int globalBin = RawYieldMap->GetBin(iCosTheta + 1); // so that we make sure there is no conflict between the bin numberings

		double normYield = 0;

		double rawYield = RawYieldMap->GetBinContent(globalBin);
		/*
			double accValue = AccMap->GetEfficiency(globalBin);
			double effValue = NomEffMap->GetEfficiency(globalBin);

			// make sure the acceptance AND the efficiency are defined for this bin
			if (accValue == 0 || effValue == 0)
				normYield = 0;
			else
				normYield = rawYield / (accValue * effValue * fabs(CosThetaBinningCS[iCosTheta] - CosThetaBinningCS[iCosTheta + 1]));
*/
		normYield = rawYield;
		NormYieldMap->SetBinContent(globalBin, normYield);

		NormYieldMap->SetBinError(globalBin, RawYieldMap->GetBinError(globalBin) * normYield / rawYield); // scale the raw yield error
	}

	// draw the corrected yield distribution

	auto* canvas = new TCanvas("canvas", "", 700, 600);

	NormYieldMap->GetXaxis()->CenterTitle();
	NormYieldMap->GetXaxis()->SetRangeUser(-1, 1);
	NormYieldMap->GetYaxis()->CenterTitle();
	NormYieldMap->Draw("APZ");

	const char* legendText = Form("%d < p_{T}^{#mu#mu} < %d GeV, Collins-Soper frame", ptMin, ptMax);

	TLatex* legend = new TLatex(.5, .88, legendText);
	legend->SetTextAlign(22);
	//legend->SetTextSize(0.045);
	legend->DrawLatexNDC(.5, .88, legendText);

	CMS_lumi(canvas, lumi_PbPb);

	/* time to fit the distribution!!
	 very basic for now, should switch to something more fancy with RooFit later...

	 parameter 0: normalization factor
	 parameter 1: lambda theta

	 variable x: cos theta
*/
	TF1* fitFunc = new TF1("fitFunc", "([0]/(3+[1])) * (1 + [1]*x*x)", CosThetaBinningCS[0], CosThetaBinningCS[NCosThetaBinsCS]); // function only defined in the region of the measurement

	fitFunc->SetParName(0, "normalization factor");
	fitFunc->SetParName(1, "lambda theta");

	//	fitFunc->SetParameters(40, 0); // initialize parameters values
	//	fitFunc->SetParLimits(1, -1.5, 1.5);

	auto fitResult = NormYieldMap->Fit(fitFunc, "QIEMSR"); // "R" option: only fit in the defined region
	fitResult->Print("v");

	TLatex* fitLegend = new TLatex(.5, .88, " ");
	fitLegend->SetTextAlign(22);
	fitLegend->SetTextSize(0.042);
	//	fitLegend->DrawLatexNDC(.5, .19, Form("2D fit result: #lambda_{#theta} = %.2f #pm %.2f, #lambda_{#varphi} = %.2f #pm %.2f", fitFunc->GetParameter(1), fitFunc->GetParError(1), fitFunc->GetParameter(2), fitFunc->GetParError(2)));

	canvas->SaveAs(Form("DistributionFits/1D/%s_pt%dto%d.png", (isCSframe) ? "CS" : "HX", ptMin, ptMax), "RECREATE");
}
