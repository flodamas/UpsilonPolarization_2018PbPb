#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/FitShortcuts.h"
#include "../Tools/Style/Legends.h"

void createSysErrorBox(TGraphErrors* gSys, const char* lineColor = "#A0B1BA") {
	int nPoints = gSys->GetN();

	for (int i = 1; i < nPoints; ++i) {
		double x, y;
		gSys->GetPoint(i, x, y);

		double ex = 0.5;                         // x half-width, as set in SetPointError
		double eyHigh = gSys->GetErrorYhigh(i);    // upper y error
		double eyLow  = gSys->GetErrorYlow(i);     // lower y error

		TBox *box = new TBox(x - ex, y - eyLow, x + ex, y + eyHigh);
		box->SetLineColorAlpha(TColor::GetColor(lineColor), 0.9);                // Outline color
		box->SetLineWidth(1);
		box->Draw("same");
	}
}

void createSysErrorBox(TGraphAsymmErrors* gSys, const char* lineColor = "#A0B1BA") {
	int nPoints = gSys->GetN();

	for (int i = 0; i < nPoints; ++i) {
		double x, y;
		gSys->GetPoint(i, x, y);

		double ex = 0.5;                         // x half-width, as set in SetPointError
		double eyHigh = gSys->GetErrorYhigh(i);    // upper y error
		double eyLow  = gSys->GetErrorYlow(i);     // lower y error

		TBox *box = new TBox(x - ex, y - eyLow, x + ex, y + eyHigh);
		box->SetLineColorAlpha(TColor::GetColor(lineColor), 0.9);                // Outline color
		box->SetLineWidth(1);
		box->Draw("same");
	}
}

void createLegendBox(double x, double y, const char* color = "#A0B1BA", Bool_t QMPoster = kFALSE) {
	float dx; // x half-width, as set in SetPointError
	float dy; // y half-width, as set in SetPointError
	
	if (QMPoster) {
		dx = 0.2;
		dy = 0.05;
	}
	else {
		dx = 0.5;
		dy = 0.08;
	}

	TBox *fillbox = new TBox(x - dx, y - dy, x + dx, y + dy);
	fillbox->SetFillColorAlpha(TColor::GetColor(color), 0.3);                // Outline color
	fillbox->SetFillStyle(1001);
	fillbox->Draw("same");

	TBox *linebox = new TBox(x - dx, y - dy, x + dx, y + dy);
	linebox->SetLineColorAlpha(TColor::GetColor(color), 0.9);  // Outline color
	linebox->SetLineWidth(1);
	linebox->Draw("same");

	gPad->Update();
}

/// read the pp data from the BPH_11_023_SupplementalMaterial.txt file
std::vector<std::vector<std::vector<TGraphAsymmErrors*>>> readppData(std::vector<std::vector<std::vector<TGraphAsymmErrors*>>> &gppDataTotErr, Bool_t QMPoster = kFALSE) {

	/// Define the variables
	TString label1, polarParamName, refFrameName, label4, label5, label6, label7, label8, label9, label10, label11, label12, label13;

	Int_t polarParamIdx, refFrameIdx;

	Bool_t isY1S = kFALSE;

	/// Define the vector to hold the histograms
	std::vector<std::vector<std::vector<TGraphAsymmErrors*>>> ppDataHists(4, std::vector<std::vector<TGraphAsymmErrors*>>(2, std::vector<TGraphAsymmErrors*>(2, nullptr))); // [polarParamIdx][refFrameIdx][yBinIdx]

	/// pt binning that pp analysis used
	int nPtBins;
	std::vector<double> pTBinning;

	// const int nPtBins = 5;
	// double pTBinning[nPtBins + 1] = {10, 12, 16, 20, 30, 50};
	// const int nPtBins = 2; // for QM poster
	// double pTBinning[nPtBins + 1] = {12, 16, 20}; // for QM poster
	if (QMPoster) {
		nPtBins = 2;
		pTBinning = {12, 16, 20};
	}
	else {
		nPtBins = 5;
		pTBinning = {10, 12, 16, 20, 30, 50};
	}
	
	int ptBinIdx[2] = {0, 0};

	/// Open the txt file
	std::ifstream file("BPH_11_023_SupplementalMaterial.txt");
	if (!file) {
		cout << "pp data file not found...." << endl;
		return {};
	}

	cout << "pp data file (BPH_11_023) opened..." << endl;

	/// Read the txt data line by line 
	std::string line;
	
    while (std::getline(file, line)) {

		if (line[0] == '%') continue; // skip the comments
		else if (line.empty()) { // skip the empty lines
			ptBinIdx[0] = 0;
			ptBinIdx[1] = 0;
			cout << "" << endl;
			continue;  
		}

		else if (line[0] == 'T') { // read the table title
			/// Read the table title
			std::istringstream iss(line);
			iss >> label1 >> polarParamName >> refFrameName >> label4 >> label5 ;

			/// switch polarParamIdx
			if (polarParamName == (TString) "Lambda-Theta,") polarParamIdx = 0;
			else if (polarParamName == (TString) "Lambda-Phi,") polarParamIdx = 1;
			else if (polarParamName == (TString) "Lambda-Theta-Phi,") polarParamIdx = 2;
			else if (polarParamName == (TString) "Lambda-Tilde,") polarParamIdx = 3;
			else {
				polarParamIdx = -1;
				continue;
			}

			/// switch refFrameIdx
			if (refFrameName == (TString) "CS") refFrameIdx = 0;
			else if (refFrameName == (TString) "HX") refFrameIdx = 1;
			else {
				refFrameIdx = -1; 
				continue;
			}

			/// create the histograms if the values are for Y1S
			if (label5 == (TString) "Y(1S):") {
				isY1S = kTRUE;
				ppDataHists[polarParamIdx][refFrameIdx][0] = new TGraphAsymmErrors(nPtBins);
				ppDataHists[polarParamIdx][refFrameIdx][1] = new TGraphAsymmErrors(nPtBins);
				// ppDataHists[polarParamIdx][refFrameIdx][0] = new TH1D(Form("ppData_%s_%s_y0to0.6", polarParamName.Data(), refFrameName.Data()), "", nPtBins, pTBinning);
				// ppDataHists[polarParamIdx][refFrameIdx][1] = new TH1D(Form("ppData_%s_%s_y0.6to1.2", polarParamName.Data(), refFrameName.Data()), "", nPtBins, pTBinning);

				gppDataTotErr[polarParamIdx][refFrameIdx][0] = new TGraphAsymmErrors(nPtBins);
				gppDataTotErr[polarParamIdx][refFrameIdx][1] = new TGraphAsymmErrors(nPtBins);
			}
			else {
				isY1S = kFALSE;
				continue;
			}

			/// Print the polarParamName and refFrameName
			cout << polarParamName << " " << refFrameName << endl;
			continue;
		}

		else if (line[0] == 'p') continue; // skip the header line
		else if (!isY1S) continue; // skip the lines that are not for Y1S
		else if (polarParamIdx == -1) continue; // skip the lines if the polarization parameter doesn't exist in the line
		else if (refFrameIdx == -1) continue; // skip the lines if the reference frame doesn't exist in the line

		/// Define variables to hold data from the txt
		Float_t pT_min, pT_max, y_min, y_max; // kinematics
		Float_t Lambda_Theta, Lambda_Phi, Lambda_Theta_Phi, Lambda_Tilde, Lambda; // polarization paramters
		Float_t TU_68_3m, TU_68_3p, TU_95_5m, TU_95_5p, TU_99_7m, TU_99_7p; // Total uncertainties down (m) and up (p)
		Float_t SU_68_3m, SU_68_3p; // Statistical uncertainties down (m) and up (p)

		std::stringstream ss(line); // convert the line to stringstream
		
		/// Read the data from the stringstream
		ss >> pT_min >> pT_max >> y_min >> y_max >> Lambda >> TU_68_3m >> TU_68_3p >> TU_95_5m >> TU_95_5p >> TU_99_7m >> TU_99_7p >> SU_68_3m >> SU_68_3p;
		cout << "pT_min: " << pT_min << ", pT_max: " << pT_max << ", y_min: " << y_min << ", y_max: " << y_max << ", " << polarParamName.Data() << ": " << Lambda << ", TU_68_3m: " << TU_68_3m << ", TU_68_3p: " << TU_68_3p << ", TU_95_5m: " << TU_95_5m << ", TU_95_5p: " << TU_95_5p << ", TU_99_7m: " << TU_99_7m << ", TU_99_7p: " << TU_99_7p << ", SU_68_3m: " << SU_68_3m << ", SU_68_3p: " << SU_68_3p << endl;

		if (QMPoster) {
			if (pT_min == 10) continue; // skip the lines that have pT_min = 0 for QM poster
			else if (pT_min == 20) continue; // skip the lines that have pT_min = 20 for QM poster
			else if (pT_min == 30) continue; // skip the lines that have pT_min = 30 for QM poster
		}

		/// obtain the mid-point of the pT bin (this has to be changed using the pT sectrum)
		Float_t pT = (pT_min + pT_max) / 2;

		/// Fill the histograms
		if (y_min == 0) { // y0to0.6
			/// fill data with statistical uncertainties
			ppDataHists[polarParamIdx][refFrameIdx][0]->SetPoint(ptBinIdx[0], pT, Lambda);
			ppDataHists[polarParamIdx][refFrameIdx][0]->SetPointError(ptBinIdx[0], pT-pT_min, pT_max-pT, -SU_68_3m, SU_68_3p);
			/// fill data with total uncertainties
			gppDataTotErr[polarParamIdx][refFrameIdx][0]->SetPoint(ptBinIdx[0], pT, Lambda);
			gppDataTotErr[polarParamIdx][refFrameIdx][0]->SetPointError(ptBinIdx[0], 0.5, 0.5, -TU_68_3m, TU_68_3p);
			ptBinIdx[0]++;
		}
		else {	// y0.6to1.2
			/// fill data with statistical uncertainties
			ppDataHists[polarParamIdx][refFrameIdx][1]->SetPoint(ptBinIdx[1], pT, Lambda);
			ppDataHists[polarParamIdx][refFrameIdx][1]->SetPointError(ptBinIdx[1], pT-pT_min, pT_max-pT, -SU_68_3m, SU_68_3p);
			/// fill data with total uncertainties
			gppDataTotErr[polarParamIdx][refFrameIdx][1]->SetPoint(ptBinIdx[1], pT, Lambda);
			gppDataTotErr[polarParamIdx][refFrameIdx][1]->SetPointError(ptBinIdx[1], 0.5, 0.5, -TU_68_3m, TU_68_3p);
			ptBinIdx[1]++;
		}
		
		cout << "ptBinIdx: " << ptBinIdx << endl;
	}

	return ppDataHists;
}

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

		// cout << "refFrameName: " << refFrameName << endl;
		// cout << (refFrameName == "CS") << endl;
		// cout << (refFrameName == "HX") << endl;

        // Save the data to the 2D vector
        if (refFrameName == "CS") sysError[0].push_back({lambdaTheta, lambdaPhi, lambdaThetaPhi, lambdaTilde});
		else sysError[1].push_back({lambdaTheta, lambdaPhi, lambdaThetaPhi, lambdaTilde});
    }

    inFile.close();
    return sysError;
}

void finalResults(Bool_t QMPoster = kFALSE, const char* bkgShapeName = "ExpTimesErr", const Int_t nCosThetaBins = 5, Double_t cosThetaMin = -0.7, Double_t cosThetaMax = 0.7, const Int_t nPhiBins = 3, Int_t phiMin = 0, Int_t phiMax = 180) {
	writeExtraText = true; // if extra text
	extraText = "       Preliminary";

	/// define the variables
	const int NRefFrame = 2;
	const int NLambParams = 4;

	const char* signalShapeName = "SymDSCB";
	const char* methodName = "rootFit";

	/// set the colors
	const char* CSColor = "#009ADE"; // blue
	const char* HXColor = "#FF1F5B"; // red

	const char* ppDataRap1Color = "#bababa"; // gold
	const char* ppDataRap2Color = "#ff9d3a"; // grey

	const char* color[NRefFrame] = {CSColor, HXColor};
	const char* ppDataColor[2] = {ppDataRap1Color, ppDataRap2Color};

	/// define the variables as a place holder of the polarization parameters
	RooRealVar* lambdaTheta = new RooRealVar("lambdaTheta", "lambdaTheta", 0);
	RooRealVar* lambdaPhi = new RooRealVar("lambdaPhi", "lambdaPhi", 0);
	RooRealVar* lambdaThetaPhi = new RooRealVar("lambdaThetaPhi", "lambdaThetaPhi", 0);
	RooRealVar* lambdaTilde = new RooRealVar("lambdaTilde", "lambdaTilde", 0);

	/// define the variables for plot titles and names
	const char* lambParamName[NLambParams] = {"lambdaTheta", "lambdaPhi", "lambdaThetaPhi", "lambdaTilde"};
	const char* lambParamTitle[NLambParams] = {"#lambda_{#theta}", "#lambda_{#varphi}", "#lambda_{#theta#varphi}", "#tilde{#lambda}"};
	const char* refFrameName[NRefFrame] = {"CS", "HX"};

	/// define the histograms
	TH1D* lambdaHist[NLambParams][NRefFrame]; // store the polarization parameters with statistical uncertainties
	TGraphErrors *gSysLamb[NLambParams][NRefFrame]; // store the systematic uncertainties

	/// read the systematic uncertainties
	std::vector<std::vector<std::vector<double>>> totalSysError = readSystematicUncertainties();

	/// read the pp data
	std::vector<std::vector<std::vector<TGraphAsymmErrors*>>> gppDataTotErr(4, std::vector<std::vector<TGraphAsymmErrors*>>(2, std::vector<TGraphAsymmErrors*>(2, nullptr))); // place holder for the pp data with total uncertainties

	std::vector<std::vector<std::vector<TGraphAsymmErrors*>>> ppDataHists = readppData(gppDataTotErr, QMPoster); // store the pp data with statistical uncertainties (CL 68.3%)

	/// Loop over the polarization parameters to fill the histograms
	for (Int_t iLambParam = 0; iLambParam < NLambParams; iLambParam++) {

		/// Loop over the reference frames
		for (Int_t iRefFrame = 0; iRefFrame < NRefFrame; iRefFrame++) {

			/// generate the polarization parameter histograms
			lambdaHist[iLambParam][iRefFrame] = new TH1D(Form("%sHist%s", lambParamName[iLambParam], refFrameName[iRefFrame]), Form("; p_{T} (GeV/c); %s", lambParamTitle[iLambParam]), NPtBins, gPtBinning);
			
			gSysLamb[iLambParam][iRefFrame] = new TGraphErrors(NPtBins);
			
			/// fill the polarization parameter histograms
			// for (Int_t ibin = 1; ibin <= NPtBins; ibin++) { // start from 0 < pT < 2 GeV bin
			for (Int_t ibin = 2; ibin <= NPtBins; ibin++) { // start from 2 < pT < 6 GeV bin
				
				/// get polarization parameters
				const char* fitModelName = GetFitModelName(signalShapeName, gPtBinning[ibin - 1], gPtBinning[ibin], refFrameName[iRefFrame], cosThetaMin, cosThetaMax, phiMin, phiMax);

				RooArgSet polarParams = GetPolarParams(lambdaTheta, lambdaPhi, lambdaThetaPhi, lambdaTilde, methodName, fitModelName, bkgShapeName);

				double lambdaVal, lambdaUnc; 

				/// Fill the histograms
				if (iLambParam == 0) {
					lambdaVal = lambdaTheta->getVal();
					lambdaUnc = lambdaTheta->getError();
				}
				else if (iLambParam == 1) {
					lambdaVal = lambdaPhi->getVal();
					lambdaUnc = lambdaPhi->getError();
				}
				else if (iLambParam == 2) {
					lambdaVal = lambdaThetaPhi->getVal();
					lambdaUnc = lambdaThetaPhi->getError();
				}
				else if (iLambParam == 3) {
					lambdaVal = lambdaTilde->getVal();
					lambdaUnc = lambdaTilde->getError();
				}

				cout << lambParamName[iLambParam] << "Val: " << lambdaVal << endl;

				lambdaHist[iLambParam][iRefFrame]->SetBinContent(ibin, lambdaVal);
				lambdaHist[iLambParam][iRefFrame]->SetBinError(ibin, lambdaUnc);

				/// Fill the systematic uncertainties
				double x = (gPtBinning[ibin - 1] + gPtBinning[ibin]) / 2;
				double y = lambdaVal;
				double xError = 0.5;

				gSysLamb[iLambParam][iRefFrame]->SetPoint(ibin - 1, x, y);
				gSysLamb[iLambParam][iRefFrame]->SetPointError(ibin - 1, xError, totalSysError[iRefFrame][ibin - 1][iLambParam]);
			}
		}
	}

	/// Draw the polarization parameters histograms
	TCanvas* polarParamsCanvas = new TCanvas("polarParamsCanvas", "polarParamsCanvas", 900, 800);

	TPad* pad[NLambParams - 1][NRefFrame]; 

	TLine* zeroLine;
	
	if (QMPoster) zeroLine = new TLine(gPtBinning[3], 0, gPtBinning[NPtBins], 0); // for QM poster
	else zeroLine = new TLine(gPtBinning[0], 0, gPtBinning[NPtBins], 0);
	
	zeroLine->SetLineStyle(kDashed);

	TLatex text;
	text.SetTextAlign(22);
	text.SetTextSize(0.08);

	/// loop for cosmetics and draw the histograms
	for (Int_t iLambParam = 0; iLambParam < NLambParams - 1; iLambParam++) {
		/// Loop over the reference frames
		for (Int_t iRefFrame = 0; iRefFrame < NRefFrame; iRefFrame++) {

			float xPadMin = 0.5 * iRefFrame;
			float xPadMax = 0.5 * (iRefFrame + 1);

			/// Draw lambdaTheta pad
			polarParamsCanvas->cd();

			/// top pad for lambdaTheta
			if (iLambParam == 0) {
				pad[iLambParam][iRefFrame] = new TPad(Form("pad1%s", refFrameName[iRefFrame]), "pad1", xPadMin, 0.672, xPadMax, 1.0);
				pad[iLambParam][iRefFrame]->SetTopMargin(0.12);
				pad[iLambParam][iRefFrame]->SetBottomMargin(0.0);

				lambdaHist[iLambParam][iRefFrame]->GetXaxis()->SetLabelSize(0);
				lambdaHist[iLambParam][iRefFrame]->GetXaxis()->SetTitleSize(0);
				lambdaHist[iLambParam][iRefFrame]->GetYaxis()->SetTitleOffset(0.75);
			}

			else if (iLambParam == 1) {
				pad[iLambParam][iRefFrame] = new TPad(Form("pad2%s", refFrameName[iRefFrame]), "pad2", xPadMin, 0.384, xPadMax, 0.672);
				pad[iLambParam][iRefFrame]->SetTopMargin(0.0);
				pad[iLambParam][iRefFrame]->SetBottomMargin(0.0);

				lambdaHist[iLambParam][iRefFrame]->GetYaxis()->SetTitleOffset(0.65);
			}

			else {
				pad[iLambParam][iRefFrame] = new TPad(Form("pad3%s", refFrameName[iRefFrame]), "pad3", xPadMin, 0.00, xPadMax, 0.384);
				pad[iLambParam][iRefFrame]->SetTopMargin(0.0);
				pad[iLambParam][iRefFrame]->SetBottomMargin(0.25);

				lambdaHist[iLambParam][iRefFrame]->GetXaxis()->SetTitleSize(0.08);
				lambdaHist[iLambParam][iRefFrame]->GetXaxis()->SetTitleOffset(1.2);
				lambdaHist[iLambParam][iRefFrame]->GetXaxis()->SetLabelSize(0.075);

				lambdaHist[iLambParam][iRefFrame]->GetYaxis()->SetTitleOffset(0.85);
			}

			if (iRefFrame == 0) {
				pad[iLambParam][iRefFrame]->SetLeftMargin(0.17); 
				pad[iLambParam][iRefFrame]->SetRightMargin(0.0);
				ppDataHists[iLambParam][iRefFrame][0]->SetMarkerStyle(33);
				ppDataHists[iLambParam][iRefFrame][1]->SetMarkerStyle(33);
			}
			else {
				pad[iLambParam][iRefFrame]->SetLeftMargin(0.0); 
				pad[iLambParam][iRefFrame]->SetRightMargin(0.17);
				ppDataHists[iLambParam][iRefFrame][0]->SetMarkerStyle(27);
				ppDataHists[iLambParam][iRefFrame][1]->SetMarkerStyle(27);
			}

			pad[iLambParam][iRefFrame]->Draw();
			pad[iLambParam][iRefFrame]->cd();

			lambdaHist[iLambParam][iRefFrame]->Draw("PL");
			ppDataHists[iLambParam][iRefFrame][0]->Draw("SAME PZ");
			ppDataHists[iLambParam][iRefFrame][1]->Draw("SAME PZ");

			/// cosmetics
			lambdaHist[iLambParam][iRefFrame]->SetMarkerStyle(20);
			lambdaHist[iLambParam][iRefFrame]->SetMarkerSize(1.1);
			lambdaHist[iLambParam][iRefFrame]->SetMarkerColor(TColor::GetColor(color[iRefFrame]));
			
			lambdaHist[iLambParam][iRefFrame]->SetLineColor(TColor::GetColor(color[iRefFrame]));
			lambdaHist[iLambParam][iRefFrame]->SetLineWidth(3);

			lambdaHist[iLambParam][iRefFrame]->GetXaxis()->CenterTitle();
			if (QMPoster) lambdaHist[iLambParam][iRefFrame]->GetXaxis()->SetRangeUser(12, 20); // this range is for QM poster
			else lambdaHist[iLambParam][iRefFrame]->GetXaxis()->SetRangeUser(0, 20);

			if (QMPoster) {
				lambdaHist[iLambParam][iRefFrame]->GetYaxis()->SetRangeUser(-0.6, 0.6); // this range is for QM poster
				lambdaHist[iLambParam][iRefFrame]->GetYaxis()->SetNdivisions(406); // for QM poster
			}
			else lambdaHist[iLambParam][iRefFrame]->GetYaxis()->SetRangeUser(-1.2, 1.2);
			lambdaHist[iLambParam][iRefFrame]->GetYaxis()->CenterTitle();

			// ppDataHists[iLambParam][iRefFrame][0]->SetMarkerStyle(33);
			ppDataHists[iLambParam][iRefFrame][0]->SetMarkerSize(1.6);
			ppDataHists[iLambParam][iRefFrame][0]->SetMarkerColor(TColor::GetColor(ppDataColor[0]));
			ppDataHists[iLambParam][iRefFrame][0]->SetLineColor(TColor::GetColor(ppDataColor[0]));
			ppDataHists[iLambParam][iRefFrame][0]->SetLineWidth(2);

			// ppDataHists[iLambParam][iRefFrame][1]->SetMarkerStyle(27);
			ppDataHists[iLambParam][iRefFrame][1]->SetMarkerSize(1.6);
			ppDataHists[iLambParam][iRefFrame][1]->SetMarkerColor(TColor::GetColor(ppDataColor[1]));
			ppDataHists[iLambParam][iRefFrame][1]->SetLineColor(TColor::GetColor(ppDataColor[1]));
			// ppDataHists[iLambParam][iRefFrame][1]->SetLineStyle(kDashed);
			ppDataHists[iLambParam][iRefFrame][1]->SetLineWidth(2);

			if (iRefFrame == 0) {
				if (iLambParam == 0) {
				lambdaHist[iLambParam][iRefFrame]->GetYaxis()->SetTitleSize(0.105);
				lambdaHist[iLambParam][iRefFrame]->GetYaxis()->SetLabelSize(0.085);
				}
				else if (iLambParam == 1) {
				lambdaHist[iLambParam][iRefFrame]->GetYaxis()->SetTitleSize(0.12);
				lambdaHist[iLambParam][iRefFrame]->GetYaxis()->SetLabelSize(0.10);	
				}
				else {
				lambdaHist[iLambParam][iRefFrame]->GetYaxis()->SetTitleSize(0.09);
				lambdaHist[iLambParam][iRefFrame]->GetYaxis()->SetLabelSize(0.075);
				}
			}
			else {
				lambdaHist[iLambParam][iRefFrame]->GetYaxis()->SetTitleSize(0);
				lambdaHist[iLambParam][iRefFrame]->GetYaxis()->SetLabelSize(0);
			}

			/// Draw systematic error band
			gSysLamb[iLambParam][iRefFrame]->SetFillColorAlpha(TColor::GetColor(color[iRefFrame]), 0.3);
			gSysLamb[iLambParam][iRefFrame]->Draw("E2 SAME");
		
			gppDataTotErr[iLambParam][iRefFrame][0]->SetFillColorAlpha(TColor::GetColor(ppDataColor[0]), 0.3);
			gppDataTotErr[iLambParam][iRefFrame][0]->Draw("E2 SAME");

			gppDataTotErr[iLambParam][iRefFrame][1]->SetFillColorAlpha(TColor::GetColor(ppDataColor[1]), 0.3);
			gppDataTotErr[iLambParam][iRefFrame][1]->Draw("E2 SAME");

			createSysErrorBox(gSysLamb[iLambParam][iRefFrame], color[iRefFrame]);
			createSysErrorBox(gppDataTotErr[iLambParam][iRefFrame][0], ppDataColor[0]);
			createSysErrorBox(gppDataTotErr[iLambParam][iRefFrame][1], ppDataColor[1]);

			gPad->Modified();
			gPad->Update();

			if (iLambParam == 0) {
				TLatex refFrameText;
				refFrameText.SetTextAlign(22);
				refFrameText.SetTextSize(0.08);

				/// draw the box around the legend symbols
				if (QMPoster) {
					createLegendBox(12.83, -0.477, color[iRefFrame], kTRUE);					
				}
				else {
					createLegendBox(2.07, -0.95, color[iRefFrame]);
				}

				if (iRefFrame == 0){
					refFrameText.DrawLatexNDC(0.82, 0.8, "Collins-Soper"); 
					text.DrawLatexNDC(.33, .8, "#varUpsilon(1S) #rightarrow #mu^{+}#mu^{-}");
					CMS_lumi(pad[iLambParam][iRefFrame], "");
	
					/// legend
					TLegend *legend = new TLegend(0.20, 0.02, 0.65, 0.16);
					legend->SetBorderSize(0);
					legend->SetFillStyle(0);
					legend->SetTextSize(0.06);
					legend->SetEntrySeparation(0.);
					legend->AddEntry(lambdaHist[iLambParam][iRefFrame], "CMS, PbPb, #sqrt{S_{NN}} = 5.02 TeV, 0 < |y| < 2.4, CS", "lep");
					legend->Draw("SAME");
				}
				else {
					refFrameText.DrawLatexNDC(0.71, 0.8, "Helicity");
					TLatex *tex = new TLatex(0.83,0.916,"PbPb 1.61 nb^{#minus1} (5.02 TeV)");
					tex->SetNDC();
					tex->SetTextAlign(31);
					tex->SetTextFont(42);
					tex->SetTextSize(0.06);
					tex->SetLineWidth(2);
					tex->Draw();

					TLegend *legend = new TLegend(0.03, 0.02, 0.46, 0.16);
					legend->SetBorderSize(0);
					legend->SetFillStyle(0);
					legend->SetTextAlign(12);
					legend->SetTextSize(0.06);
					legend->AddEntry(lambdaHist[iLambParam][iRefFrame], "CMS, PbPb, #sqrt{S_{NN}} = 5.02 TeV, 0 < |y| < 2.4, HX", "lep");
					legend->Draw("SAME");
				}
			}

			else if (iLambParam == 1) {
				if (QMPoster) {
					createLegendBox(12.83, 0.475, ppDataColor[0], kTRUE);
					createLegendBox(12.83, 0.3, ppDataColor[1], kTRUE);						
				}
				else {
					createLegendBox(2.07, 0.95, ppDataColor[0]);
					createLegendBox(2.07, 0.6, ppDataColor[1]);
				}

				if (iRefFrame == 0) {
					/// legend
					TLegend *legend = new TLegend(0.2, 0.68, 0.65, 0.97);
					legend->SetBorderSize(0);
					legend->SetFillStyle(0);
					legend->SetTextSize(0.07);
					legend->AddEntry(ppDataHists[iLambParam][iRefFrame][0], "CMS, pp, #sqrt{S} = 7 TeV, 0 < |y| < 0.6, CS", "lep");
					legend->AddEntry(ppDataHists[iLambParam][iRefFrame][1], "CMS, pp, #sqrt{S} = 7 TeV, 0.6 < |y| < 1.2, CS", "lep");
					legend->Draw("SAME");
				}
				else {
					/// legend
					TLegend *legend = new TLegend(0.03, 0.68, 0.46, 0.97);
					legend->SetBorderSize(0);
					legend->SetFillStyle(0);
					legend->SetTextAlign(12);
					legend->SetTextSize(0.07);
					legend->AddEntry(ppDataHists[iLambParam][iRefFrame][0], "CMS, pp, #sqrt{S} = 7 TeV, 0 < |y| < 0.6, HX", "lep");
					legend->AddEntry(ppDataHists[iLambParam][iRefFrame][1], "CMS, pp, #sqrt{S} = 7 TeV, 0.6 < |y| < 1.2, HX", "lep");
					legend->Draw("SAME");
				}
			}

			else if (iLambParam == 2) {
				if (iRefFrame == 0) {
					
					if (QMPoster) {
						/// white background
						TPaveText *paveText = new TPaveText(0.96, 0.01, 1, 0.24, "brNDC");
						paveText->SetFillColor(kWhite);  // White background
						paveText->SetFillStyle(1001);
						paveText->SetTextFont(42);
						paveText->SetBorderSize(0);
						paveText->SetMargin(0);
						paveText->Draw();
						
						TLatex *latex = new TLatex();
						latex->SetNDC();  // Use normalized coordinates
						latex->SetTextSize(0.075);
						latex->DrawLatex(0.975, 0.18, "1"); 
					}
					else {/// white background
						TPaveText *paveText = new TPaveText(0.96, 0.01, 1, 0.24, "brNDC");
						paveText->SetFillColor(kWhite);  // White background
						paveText->SetFillStyle(1001);
						paveText->SetTextFont(42);
						paveText->SetBorderSize(0);
						paveText->SetMargin(0);
						paveText->Draw();
						
						TLatex *latex = new TLatex();
						latex->SetNDC();  // Use normalized coordinates
						latex->SetTextSize(0.075);
						latex->DrawLatex(0.987, 0.18, "0"); 
					}
						gPad->Update();
				}
			}
			zeroLine->Draw("SAME");
		}
	}

	/// save the canvas
	const char* fitModelName2 = GetFitModelName(signalShapeName, gPtBinning[1], gPtBinning[NPtBins], "CSHX", cosThetaMin, cosThetaMax, phiMin, phiMax);
	
	gSystem->mkdir("ParametersResults", kTRUE);
	if (QMPoster) polarParamsCanvas->SaveAs(Form("ParametersResults/polarParams_%s_%s_%s_QMPoster.png", methodName, fitModelName2, bkgShapeName), "RECREATE");
	else polarParamsCanvas->SaveAs(Form("ParametersResults/polarParams_%s_%s_%s.png", methodName, fitModelName2, bkgShapeName), "RECREATE");

	/// draw lambdaTilde histogram
	TCanvas* invPolarParamsCanvas = new TCanvas("invPolarParamsCanvas", "invPolarParamsCanvas", 500, 500);

	lambdaHist[3][0]->Draw();
	
	lambdaHist[3][0]->SetMarkerStyle(20);
	lambdaHist[3][0]->SetMarkerSize(1.4);
	lambdaHist[3][0]->SetMarkerColor(TColor::GetColor(CSColor));
	
	lambdaHist[3][0]->SetLineColor(TColor::GetColor(CSColor));
	lambdaHist[3][0]->SetLineWidth(3);

	lambdaHist[3][1]->Draw("same");

	lambdaHist[3][1]->SetMarkerStyle(20);
	lambdaHist[3][1]->SetMarkerSize(1.4);
	lambdaHist[3][1]->SetMarkerColor(TColor::GetColor(HXColor));

	lambdaHist[3][1]->SetLineColor(TColor::GetColor(HXColor));
	lambdaHist[3][1]->SetLineWidth(3);

	// cosmetics
	lambdaHist[3][0]->GetXaxis()->CenterTitle();

	lambdaHist[3][0]->GetYaxis()->SetRangeUser(-1.4, 1.4);
	lambdaHist[3][0]->GetYaxis()->CenterTitle();

	gSysLamb[3][0]->SetFillColorAlpha(TColor::GetColor(CSColor), 0.3);

	gSysLamb[3][0]->Draw("E2 SAME");

	createSysErrorBox(gSysLamb[3][0], CSColor);

	gSysLamb[3][1]->SetFillColorAlpha(TColor::GetColor(HXColor), 0.3);

	gSysLamb[3][1]->Draw("E2 SAME");

	createSysErrorBox(gSysLamb[3][1], HXColor);

	zeroLine->Draw("SAME");

	CMS_lumi(invPolarParamsCanvas, gCMSLumiText);

	text.SetTextSize(0.058);
	text.DrawLatexNDC(.32, .85, "#varUpsilon(1S) #rightarrow #mu^{+}#mu^{-}");

	TLegend* legend = new TLegend(0.6, 0.77, 0.9, 0.9);
	legend->SetBorderSize(0);
	legend->SetFillStyle(0);
	legend->SetTextSize(0.05);
	legend->AddEntry(lambdaHist[3][0], "Collins-Soper", "lp");
	legend->AddEntry(lambdaHist[3][1], "Helicity", "lp");
	legend->Draw("SAME");

	invPolarParamsCanvas->SaveAs(Form("ParametersResults/invariantParameter_%s_%s_%s.png", methodName, fitModelName2, bkgShapeName), "RECREATE");

    /// output file for the polarization parameter table
    const char* outFileName = "polarizationParameters_table.txt";
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
				
				outFile << std::fixed << std::setprecision(3);
				outFile << " & $" << lambdaHist[index][iRefFrame]->GetBinContent(ibin + 1) << " \\pm " << lambdaHist[index][iRefFrame]->GetBinError(ibin + 1) << " \\pm " << gSysLamb[index][iRefFrame]->GetErrorY(ibin) << "$"; 
			}
			outFile << " \\\\" << std::endl;
		}

		outFile << "        \\hline" << std::endl;
	}

	/// write the table footer
	outFile << "        \\end{tabular}" << std::endl;
	// outFile << "    }" << std::endl;
	outFile << "    \\caption{Polarization parameters measured for each \\pt range for both helicity and Collins-Soper frames. The first uncertainty is statistical while the second is systematic.}" << std::endl;
	outFile << "    \\label{tab:pola_results}" << std::endl;
	outFile << "\\end{table}" << std::endl;
	
	/// close the file
	outFile.close();

	return;
}
