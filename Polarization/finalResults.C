#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/FitShortcuts.h"
#include "../Tools/Style/Legends.h"

std::vector<std::vector<std::vector<double>>> readSystematicUncertainties() {
    std::ifstream inFile("../SystematicUncertainties/systematic_errors_total.txt");
    std::vector<std::vector<std::vector<double>>> sysError(2);
	std::string line;

    if (!inFile) {
        std::cerr << "Error: Could not open file!" << std::endl;
        return sysError;
    }

    /// Skip the header lines
    std::getline(inFile, line); // Skip the first header line

    /// Read the remaining lines
    while (std::getline(inFile, line)) {

        std::istringstream iss(line);
        std::string refFrameName, ptRange;
        double lambdaTheta, lambdaPhi, lambdaThetaPhi, lambdaTilde;

        // Parse the line
        iss >> refFrameName >> ptRange >> lambdaTheta >> lambdaPhi >> lambdaThetaPhi >> lambdaTilde;

		cout << "refFrameName: " << refFrameName << endl;
		cout << (refFrameName == "CS") << endl;
		cout << (refFrameName == "HX") << endl;

        // Save the data to the 2D vector
        if (refFrameName == "CS") sysError[0].push_back({lambdaTheta, lambdaPhi, lambdaThetaPhi, lambdaTilde});
		else sysError[1].push_back({lambdaTheta, lambdaPhi, lambdaThetaPhi, lambdaTilde});
    }

    inFile.close();
    return sysError;
}

void createSysErrorBox(TGraphErrors* gSys, const char* lineColor = "#A0B1BA") {
	int nPoints = gSys->GetN();

	for (int i = 1; i < nPoints; ++i) {
		double x, y;
		gSys->GetPoint(i, x, y);
		double ex = 0.5;                         // x half-width, as set in SetPointError
		double eyHigh = gSys->GetErrorYhigh(i);    // upper y error
		double eyLow  = gSys->GetErrorYlow(i);     // lower y error

		TBox *box = new TBox(x - ex, y - eyLow, x + ex, y + eyHigh);
		box->SetLineColorAlpha(TColor::GetColor(lineColor), 0.6);                // Outline color
		box->SetLineWidth(1);
		box->Draw("same");
	}
}

void finalResults(const char* bkgShapeName = "ExpTimesErr", const Int_t nCosThetaBins = 5, Double_t cosThetaMin = -0.7, Double_t cosThetaMax = 0.7, const Int_t nPhiBins = 3, Int_t phiMin = 0, Int_t phiMax = 180) {
	writeExtraText = true; // if extra text
	extraText = "       Internal";

	const int NRefFrame = 2;

	RooRealVar* lambdaTheta = new RooRealVar("lambdaTheta", "lambdaTheta", 0);
	RooRealVar* lambdaPhi = new RooRealVar("lambdaPhi", "lambdaPhi", 0);
	RooRealVar* lambdaThetaPhi = new RooRealVar("lambdaThetaPhi", "lambdaThetaPhi", 0);
	RooRealVar* lambdaTilde = new RooRealVar("lambdaTilde", "lambdaTilde", 0);

	TH1D* lambdaThetaHist[NRefFrame]; 
	TH1D* lambdaPhiHist[NRefFrame]; 
	TH1D* lambdaThetaPhiHist[NRefFrame]; 
	TH1D* lambdaTildeHist[NRefFrame];

	const char* CSColor = "#009ADE";
	const char* HXColor = "#FF1F5B";

	const char* color[NRefFrame] = {CSColor, HXColor};

	const char* signalShapeName = "SymDSCB";
	const char* methodName = "rootFit";

	std::vector<std::vector<std::vector<double>>> totalSysError = readSystematicUncertainties();

	TGraphErrors *gSysLambTheta[NRefFrame]; 
	TGraphErrors *gSysLambPhi[NRefFrame]; 
	TGraphErrors *gSysLambThetaPhi[NRefFrame]; 
	TGraphErrors *gSysLambTilde[NRefFrame];

	const char* refFrameName = "CS";

	/// Loop over the reference frames
	for (Int_t iRefFrame = 0; iRefFrame < NRefFrame; iRefFrame++) {
			
		if (iRefFrame == 0) {
			refFrameName = "CS";
		} else {
			refFrameName = "HX";
		}

		/// generate the polarization parameter histograms
		lambdaThetaHist[iRefFrame] = new TH1D(Form("lambdaTildeHist%s", "refFrameName"), "; p_{T} (GeV/c); #lambda_{#theta}", NPtBins, gPtBinning);
		lambdaPhiHist[iRefFrame] = new TH1D(Form("lambdaTildeHist%s", "refFrameName"), "; p_{T} (GeV/c); #lambda_{#varphi}", NPtBins, gPtBinning);
		lambdaThetaPhiHist[iRefFrame] = new TH1D(Form("lambdaTildeHist%s", "refFrameName"), "; p_{T} (GeV/c); #lambda_{#theta#varphi}", NPtBins, gPtBinning);
		lambdaTildeHist[iRefFrame] = new TH1D(Form("lambdaTildeHist%s", "refFrameName"), "; p_{T} (GeV/c); #tilde{#lambda}", NPtBins, gPtBinning);
			
		gSysLambTheta[iRefFrame] = new TGraphErrors(NPtBins);
		gSysLambPhi[iRefFrame] = new TGraphErrors(NPtBins);
		gSysLambThetaPhi[iRefFrame] = new TGraphErrors(NPtBins);
		gSysLambTilde[iRefFrame] = new TGraphErrors(NPtBins);
        
		/// fill the polarization parameter histograms
		// for (Int_t ibin = 1; ibin <= NPtBins; ibin++) { // start from 0 < pT < 2 GeV bin
		for (Int_t ibin = 2; ibin <= NPtBins; ibin++) { // start from 2 < pT < 6 GeV bin
			
			/// get polarization parameters
			const char* fitModelName = GetFitModelName(signalShapeName, gPtBinning[ibin - 1], gPtBinning[ibin], refFrameName, cosThetaMin, cosThetaMax, phiMin, phiMax);

			RooArgSet polarParams = GetPolarParams(lambdaTheta, lambdaPhi, lambdaThetaPhi, lambdaTilde, methodName, fitModelName, bkgShapeName);

			double lambdaThetaVal = lambdaTheta->getVal();
			double lambdaThetaUnc = lambdaTheta->getError();

			double lambdaPhiVal = lambdaPhi->getVal();
			double lambdaPhiUnc = lambdaPhi->getError();

			double lambdaThetaPhiVal = lambdaThetaPhi->getVal();
			double lambdaThetaPhiUnc = lambdaThetaPhi->getError();

			double lambdaTildeVal = lambdaTilde->getVal();
			double lambdaTildeUnc = lambdaTilde->getError();

			cout << "lambdaThetaVal: " << lambdaThetaVal << endl;

			/// Fill the histograms
			lambdaThetaHist[iRefFrame]->SetBinContent(ibin, lambdaThetaVal);
			lambdaThetaHist[iRefFrame]->SetBinError(ibin, lambdaThetaUnc);

			lambdaPhiHist[iRefFrame]->SetBinContent(ibin, lambdaPhiVal);
			lambdaPhiHist[iRefFrame]->SetBinError(ibin, lambdaPhiUnc);

			lambdaThetaPhiHist[iRefFrame]->SetBinContent(ibin, lambdaThetaPhiVal);
			lambdaThetaPhiHist[iRefFrame]->SetBinError(ibin, lambdaThetaPhiUnc);

			lambdaTildeHist[iRefFrame]->SetBinContent(ibin, lambdaTildeVal);
			lambdaTildeHist[iRefFrame]->SetBinError(ibin, lambdaTildeUnc);

			/// Fill the systematic uncertainties
			double x = (gPtBinning[ibin - 1] + gPtBinning[ibin]) / 2;
			double y = lambdaThetaVal;
			double xError = 0.5;

			gSysLambTheta[iRefFrame]->SetPoint(ibin - 1, x, lambdaThetaVal);
			gSysLambTheta[iRefFrame]->SetPointError(ibin - 1, xError, totalSysError[iRefFrame][ibin - 1][0]);

			gSysLambPhi[iRefFrame]->SetPoint(ibin - 1, x, lambdaPhiVal);
			gSysLambPhi[iRefFrame]->SetPointError(ibin - 1, xError, totalSysError[iRefFrame][ibin - 1][1]);

			gSysLambThetaPhi[iRefFrame]->SetPoint(ibin - 1, x, lambdaThetaPhiVal);
			gSysLambThetaPhi[iRefFrame]->SetPointError(ibin - 1, xError, totalSysError[iRefFrame][ibin - 1][2]);

			gSysLambTilde[iRefFrame]->SetPoint(ibin - 1, x, lambdaTildeVal);
			gSysLambTilde[iRefFrame]->SetPointError(ibin - 1, xError, totalSysError[iRefFrame][ibin - 1][3]);
		
		}
	}

	/// Draw the polarization parameters histograms
	TCanvas* polarParamsCanvas = new TCanvas("polarParamsCanvas", "polarParamsCanvas", 900, 800);

	TPad* pad1[NRefFrame]; 
	TPad* pad2[NRefFrame]; 
	TPad* pad3[NRefFrame]; 

	TLine* zeroLine = new TLine(gPtBinning[0], 0, gPtBinning[NPtBins], 0);
	zeroLine->SetLineStyle(kDashed);

	TLatex text;
	text.SetTextAlign(22);
	text.SetTextSize(0.08);

	CMS_lumi(polarParamsCanvas, gCMSLumiText);

	for (Int_t iRefFrame = 0; iRefFrame < NRefFrame; iRefFrame++) {
			
		if (iRefFrame == 0) {
			refFrameName = "CS";
		} else {
			refFrameName = "HX";
		}

		float xPadMin = 0.5 * iRefFrame;
		float xPadMax = 0.5 * (iRefFrame + 1);

		/// Draw lambdaTheta pad
		polarParamsCanvas->cd();

		pad1[iRefFrame] = new TPad(Form("pad1%s", refFrameName), "pad1", xPadMin, 0.672, xPadMax, 1.0);
		
		pad1[iRefFrame]->SetTopMargin(0.12);
		pad1[iRefFrame]->SetBottomMargin(0.0);
		if (iRefFrame == 0) {pad1[iRefFrame]->SetLeftMargin(0.17); pad1[iRefFrame]->SetRightMargin(0.0);}
		else {pad1[iRefFrame]->SetLeftMargin(0.0); pad1[iRefFrame]->SetRightMargin(0.17);}
		pad1[iRefFrame]->Draw();
		pad1[iRefFrame]->cd();
		
		/// cosmetics
		lambdaThetaHist[iRefFrame]->Draw("PL");
		lambdaThetaHist[iRefFrame]->SetMarkerStyle(20);
		lambdaThetaHist[iRefFrame]->SetMarkerSize(1.4);
		lambdaThetaHist[iRefFrame]->SetMarkerColor(TColor::GetColor(color[iRefFrame]));
		
		lambdaThetaHist[iRefFrame]->SetLineColor(TColor::GetColor(color[iRefFrame]));
		lambdaThetaHist[iRefFrame]->SetLineWidth(3);

		lambdaThetaHist[iRefFrame]->GetXaxis()->CenterTitle();
		lambdaThetaHist[iRefFrame]->GetXaxis()->SetLabelSize(0);
		lambdaThetaHist[iRefFrame]->GetXaxis()->SetTitleSize(0);

		lambdaThetaHist[iRefFrame]->GetYaxis()->SetRangeUser(-1.2, 1.2);
		lambdaThetaHist[iRefFrame]->GetYaxis()->CenterTitle();

		if (iRefFrame == 0) {
			lambdaThetaHist[iRefFrame]->GetYaxis()->SetTitleSize(0.105);
			lambdaThetaHist[iRefFrame]->GetYaxis()->SetLabelSize(0.085);
		}
		else {
			lambdaThetaHist[iRefFrame]->GetYaxis()->SetTitleSize(0);
			lambdaThetaHist[iRefFrame]->GetYaxis()->SetLabelSize(0);
		}

		lambdaThetaHist[iRefFrame]->GetYaxis()->SetTitleOffset(0.75);

		/// Draw systematic error band
		gSysLambTheta[iRefFrame]->SetFillColorAlpha(TColor::GetColor(color[iRefFrame]), 0.3);

		gSysLambTheta[iRefFrame]->Draw("E2 SAME");
		
		createSysErrorBox(gSysLambTheta[iRefFrame], color[iRefFrame]);

		TLatex refFrameText;
		refFrameText.SetTextAlign(22);
		refFrameText.SetTextSize(0.08);
		if (iRefFrame == 0){
			refFrameText.DrawLatexNDC(0.82, 0.8, "Collins-Soper"); 
			text.DrawLatexNDC(.31, .8, "#varUpsilon(1S) #rightarrow #mu^{+}#mu^{-}");
		}
		else
			refFrameText.DrawLatexNDC(0.71, 0.8, "Helicity");

		CMS_lumi(pad1[iRefFrame], gCMSLumiText);

		zeroLine->Draw("SAME");

		/// Draw lambdaPhi pad
		polarParamsCanvas->cd();

		pad2[iRefFrame] = new TPad(Form("pad2%s", refFrameName), "pad2", xPadMin, 0.384, xPadMax, 0.672);

		pad2[iRefFrame]->SetTopMargin(0.0);
		pad2[iRefFrame]->SetBottomMargin(0.0);
		if (iRefFrame == 0) {pad2[iRefFrame]->SetLeftMargin(0.17); pad2[iRefFrame]->SetRightMargin(0.0);}
		else {pad2[iRefFrame]->SetLeftMargin(0.0); pad2[iRefFrame]->SetRightMargin(0.17);}
		pad2[iRefFrame]->Draw();
		pad2[iRefFrame]->cd();

		lambdaPhiHist[iRefFrame]->Draw();
		/// cosmetics
		lambdaPhiHist[iRefFrame]->SetMarkerStyle(20);
		lambdaPhiHist[iRefFrame]->SetMarkerSize(1.4);
		lambdaPhiHist[iRefFrame]->SetMarkerColor(TColor::GetColor(color[iRefFrame]));
		
		lambdaPhiHist[iRefFrame]->SetLineColor(TColor::GetColor(color[iRefFrame]));
		lambdaPhiHist[iRefFrame]->SetLineWidth(3);

		lambdaPhiHist[iRefFrame]->GetXaxis()->CenterTitle();

		lambdaPhiHist[iRefFrame]->GetYaxis()->SetRangeUser(-1.2, 1.2);
		lambdaPhiHist[iRefFrame]->GetYaxis()->CenterTitle();

		if (iRefFrame == 0) {
			lambdaPhiHist[iRefFrame]->GetYaxis()->SetTitleSize(0.12);
			lambdaPhiHist[iRefFrame]->GetYaxis()->SetLabelSize(0.10);
		}

		else {
			lambdaPhiHist[iRefFrame]->GetYaxis()->SetTitleSize(0);
			lambdaPhiHist[iRefFrame]->GetYaxis()->SetLabelSize(0);
		}

		lambdaPhiHist[iRefFrame]->GetYaxis()->SetTitleOffset(0.65);

		gSysLambPhi[iRefFrame]->SetFillColorAlpha(TColor::GetColor(color[iRefFrame]), 0.3);

		gSysLambPhi[iRefFrame]->Draw("E2 SAME");
		
		createSysErrorBox(gSysLambPhi[iRefFrame], color[iRefFrame]);

		TLatex kinematicText;
		kinematicText.SetTextAlign(13);
		kinematicText.SetTextSize(0.085);
		if (iRefFrame == 0) kinematicText.DrawLatexNDC(.29, .85, Form("%s, %s", CentralityRangeText(), DimuonRapidityRangeText(0, 2.4)));
		else kinematicText.DrawLatexNDC(.14, .85, Form("%s, %s", CentralityRangeText(), DimuonRapidityRangeText(0, 2.4)));

		zeroLine->Draw("SAME");

		/// Draw lambdaThetaPhi pad
		polarParamsCanvas->cd();

		pad3[iRefFrame] = new TPad(Form("pad3%s", refFrameName), "pad3", xPadMin, 0.00, xPadMax, 0.384);

		pad3[iRefFrame]->SetTopMargin(0.0);
		pad3[iRefFrame]->SetBottomMargin(0.25);
		if (iRefFrame == 0) {pad3[iRefFrame]->SetLeftMargin(0.17); pad3[iRefFrame]->SetRightMargin(0.0);}
		else {pad3[iRefFrame]->SetLeftMargin(0.0); pad3[iRefFrame]->SetRightMargin(0.17);}
		pad3[iRefFrame]->Draw();
		pad3[iRefFrame]->cd();

		lambdaThetaPhiHist[iRefFrame]->Draw();
		
		/// cosmetics
		lambdaThetaPhiHist[iRefFrame]->SetMarkerStyle(20);
		lambdaThetaPhiHist[iRefFrame]->SetMarkerSize(1.4);
		lambdaThetaPhiHist[iRefFrame]->SetMarkerColor(TColor::GetColor(color[iRefFrame]));
		
		lambdaThetaPhiHist[iRefFrame]->SetLineColor(TColor::GetColor(color[iRefFrame]));
		lambdaThetaPhiHist[iRefFrame]->SetLineWidth(3);

		lambdaThetaPhiHist[iRefFrame]->GetXaxis()->CenterTitle();
		lambdaThetaPhiHist[iRefFrame]->GetXaxis()->SetTitleSize(0.08);
		lambdaThetaPhiHist[iRefFrame]->GetXaxis()->SetTitleOffset(1.2);

		lambdaThetaPhiHist[iRefFrame]->GetXaxis()->SetLabelSize(0.075);

		lambdaThetaPhiHist[iRefFrame]->GetYaxis()->SetRangeUser(-1.2, 1.2);
		lambdaThetaPhiHist[iRefFrame]->GetYaxis()->CenterTitle();
		
		if (iRefFrame == 0) {
			lambdaThetaPhiHist[iRefFrame]->GetYaxis()->SetTitleSize(0.09);
			lambdaThetaPhiHist[iRefFrame]->GetYaxis()->SetLabelSize(0.075);
		}

		else {
			lambdaThetaPhiHist[iRefFrame]->GetYaxis()->SetTitleSize(0);
			lambdaThetaPhiHist[iRefFrame]->GetYaxis()->SetLabelSize(0);
		}

		lambdaThetaPhiHist[iRefFrame]->GetYaxis()->SetTitleOffset(0.85);

		gSysLambThetaPhi[iRefFrame]->SetFillColorAlpha(TColor::GetColor(color[iRefFrame]), 0.3);

		gSysLambThetaPhi[iRefFrame]->Draw("E2 SAME");
		
		createSysErrorBox(gSysLambThetaPhi[iRefFrame], color[iRefFrame]);

		zeroLine->Draw("SAME");
	}

	/// draw lambdaTilde histogram
	const char* fitModelName2 = GetFitModelName(signalShapeName, gPtBinning[1], gPtBinning[NPtBins], "CSHX", cosThetaMin, cosThetaMax, phiMin, phiMax);

	gSystem->mkdir("ParametersResults", kTRUE);
	polarParamsCanvas->SaveAs(Form("ParametersResults/polarParams_%s_%s_%s.png", methodName, fitModelName2, bkgShapeName), "RECREATE");

	TCanvas* invPolarParamsCanvas = new TCanvas("invPolarParamsCanvas", "invPolarParamsCanvas", 500, 500);

	lambdaTildeHist[0]->Draw();
	
	lambdaTildeHist[0]->SetMarkerStyle(20);
	lambdaTildeHist[0]->SetMarkerSize(1.4);
	lambdaTildeHist[0]->SetMarkerColor(TColor::GetColor(CSColor));
	
	lambdaTildeHist[0]->SetLineColor(TColor::GetColor(CSColor));
	lambdaTildeHist[0]->SetLineWidth(3);

	lambdaTildeHist[1]->Draw("same");

	lambdaTildeHist[1]->SetMarkerStyle(20);
	lambdaTildeHist[1]->SetMarkerSize(1.4);
	lambdaTildeHist[1]->SetMarkerColor(TColor::GetColor(HXColor));

	lambdaTildeHist[1]->SetLineColor(TColor::GetColor(HXColor));
	lambdaTildeHist[1]->SetLineWidth(3);

	// cosmetics
	lambdaTildeHist[0]->GetXaxis()->CenterTitle();

	lambdaTildeHist[0]->GetYaxis()->SetRangeUser(-1.4, 1.4);
	lambdaTildeHist[0]->GetYaxis()->CenterTitle();

	gSysLambTilde[0]->SetFillColorAlpha(TColor::GetColor(CSColor), 0.3);

	gSysLambTilde[0]->Draw("E2 SAME");

	createSysErrorBox(gSysLambTilde[0], CSColor);

	gSysLambTilde[1]->SetFillColorAlpha(TColor::GetColor(HXColor), 0.3);

	gSysLambTilde[1]->Draw("E2 SAME");

	createSysErrorBox(gSysLambTilde[1], HXColor);

	zeroLine->Draw("SAME");

	CMS_lumi(invPolarParamsCanvas, gCMSLumiText);

	text.SetTextSize(0.058);
	text.DrawLatexNDC(.32, .85, "#varUpsilon(1S) #rightarrow #mu^{+}#mu^{-}");

	TLegend* legend = new TLegend(0.6, 0.77, 0.9, 0.9);
	legend->SetBorderSize(0);
	legend->SetFillStyle(0);
	legend->SetTextSize(0.05);
	legend->AddEntry(lambdaTildeHist[0], "Collins-Soper", "lp");
	legend->AddEntry(lambdaTildeHist[1], "Helicity", "lp");
	legend->Draw("SAME");

	invPolarParamsCanvas->SaveAs(Form("ParametersResults/invariantParameter_%s_%s_%s.png", methodName, fitModelName2, bkgShapeName), "RECREATE");

    /// output file for the polarization parameter table
    const char* outFileName = "polarizationParameters.txt";
    std::ofstream outFile(outFileName);

    /// write table headers for the polarization parameter table
    outFile << "\\begin{table}[htpb]" << std::endl;
    outFile << "    \\centering" << std::endl;
    // outFile << "    \\resizebox{\\textwidth}{!}{" << std::endl;
    outFile << "        \\begin{tabular}{|c|c|c|c|}" << std::endl;
    outFile << "        \\hline" << std::endl;
    outFile << "        & \\pt (\\GeVc) & Collins-Soper & Helicity \\\\ " << std::endl;
    outFile << "        \\hline" << std::endl;
    outFile << "        \\hline" << std::endl;

	/// loop over the polarization parameter names
	for (Int_t index = 0; index < 3; index++){
		if (index == 0) outFile << "	\\multirow{3}*{\\begin{tabular}[l]{@{}c@{}} $\\lambda_\\theta$ \\end{tabular}}" << std::endl;	
		if (index == 1) outFile << "	\\multirow{3}*{\\begin{tabular}[l]{@{}c@{}} $\\lambda_\\varphi$ \\end{tabular}}" << std::endl;
		if (index == 2) outFile << "	\\multirow{3}*{\\begin{tabular}[l]{@{}c@{}} $\\lambda_{\\theta\\varphi}$ \\end{tabular}}" << std::endl;
		
		/// loop over the pt bins
		for (Int_t ibin = 1; ibin < NPtBins; ibin++) {

			/// set the pt bin range (second column of the table)
			outFile << std::fixed << std::setprecision(0);
			outFile << "        & $" << gPtBinning[ibin] << " < \\pt < " << gPtBinning[ibin + 1] << "$"; 

			/// loop over the reference frames
			for (Int_t iRefFrame = 0; iRefFrame < NRefFrame; iRefFrame++) {

				if (iRefFrame == 0) {
					refFrameName = "CS";
				} else {
					refFrameName = "HX";
				}
				
				outFile << std::fixed << std::setprecision(3);
				if (index == 0) outFile << " & $" << lambdaThetaHist[iRefFrame]->GetBinContent(ibin + 1) << " \\pm " << lambdaThetaHist[iRefFrame]->GetBinError(ibin + 1) << " \\pm " << gSysLambTheta[iRefFrame]->GetErrorY(ibin) << "$"; 
				else if (index == 1) outFile << " & $" << lambdaPhiHist[iRefFrame]->GetBinContent(ibin + 1) << " \\pm " << lambdaPhiHist[iRefFrame]->GetBinError(ibin + 1) << " \\pm " << gSysLambPhi[iRefFrame]->GetErrorY(ibin) << "$"; 
				else if (index == 2) outFile << " & $" << lambdaThetaPhiHist[iRefFrame]->GetBinContent(ibin + 1) << " \\pm " << lambdaThetaPhiHist[iRefFrame]->GetBinError(ibin + 1) << " \\pm " << gSysLambThetaPhi[iRefFrame]->GetErrorY(ibin) << "$";
			
			}
			outFile << " \\\\" << std::endl;
		}

		outFile << "        \\hline" << std::endl;
	}

	/// write the table footer
	outFile << "        \\end{tabular}" << std::endl;
	// outFile << "    }" << std::endl;
	outFile << "    \\caption{Polarization parameters measured for each \pt range for both helicity and Collins-Soper frames. The first uncertainty is statistical while the second is systematic.}" << std::endl;
	outFile << "    \\label{tab:pola_results}" << std::endl;
	outFile << "\\end{table}" << std::endl;
	
	/// close the file
	outFile.close();

	return;
}
