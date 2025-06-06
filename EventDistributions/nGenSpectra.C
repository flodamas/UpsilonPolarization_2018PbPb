#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

void makePtSpectra_pythiaSignal(int iState = 1) {
	// Read GenOnly Nofilter file with polarization weights
	const char* filename = Form("../Files/OniaTree_Y%dS_GENONLY_NoFilter.root", iState);
    // const char* filename = Form("../Files/Oniatree_Upsilon%dS_5p02TeV_TuneCP5_141X_GenOnly_merge.root", iState);

	TFile* file = TFile::Open(filename, "READ");
	if (!file) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
		return;
	}

	cout << "File " << filename << " opened" << endl;

	TTree* OniaTree = (TTree*)file->Get("hionia/myTree");

	/// OniaTree variables, quite old version (since this genonly file from 2015)
	Int_t Gen_QQ_size;

	TClonesArray* Gen_QQ_4mom = nullptr;
	TClonesArray* Gen_QQ_mumi_4mom = nullptr;
	TClonesArray* Gen_QQ_mupl_4mom = nullptr;

	OniaTree->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size);
	OniaTree->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom);

	OniaTree->SetBranchAddress("Gen_QQ_mumi_4mom", &Gen_QQ_mumi_4mom);
	OniaTree->SetBranchAddress("Gen_QQ_mupl_4mom", &Gen_QQ_mupl_4mom);

	TLorentzVector* gen_QQ_LV = new TLorentzVector();
	TLorentzVector* gen_mumi_LV = new TLorentzVector();
	TLorentzVector* gen_mupl_LV = new TLorentzVector();

    Long64_t totEntries = OniaTree->GetEntries(); 
    
    // const int nBins = 6;
    // float ptBinning[nBins + 1] = {0, 2, 4, 6, 9, 12, 30};

    double ptBinning[] = {0.0, 2.0, 6.0, 12.0, 30.0};
    int nBins = sizeof(ptBinning)/sizeof(ptBinning[0]) - 1;

    TH1D* hGenPt = new TH1D("hGenPt", "", nBins, ptBinning);
    hGenPt->Sumw2();

 	// Loop over the events
    for (Long64_t iEvent = 0; iEvent < (totEntries); iEvent++) {
		if (iEvent % 10000 == 0) {
			cout << Form("\rProcessing event %lld / %lld (%.0f%%)", iEvent, totEntries, 100. * iEvent / totEntries) << flush;
		}

		OniaTree->GetEntry(iEvent);

		// loop over all gen upsilons
		for (int iGen = 0; iGen < Gen_QQ_size; iGen++) {
            gen_QQ_LV = (TLorentzVector*)Gen_QQ_4mom->At(iGen);

            if (fabs(gen_QQ_LV->Rapidity()) < gRapidityMin || fabs(gen_QQ_LV->Rapidity()) > gRapidityMax) continue;
            hGenPt->Fill(gen_QQ_LV->Pt());

        }
    }

    TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
    hGenPt->SetTitle("");
    hGenPt->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    hGenPt->GetYaxis()->SetTitle("Counts");
    hGenPt->SetLineColor(kBlue);
    hGenPt->SetMarkerStyle(20);
    hGenPt->SetMarkerColor(kBlue);
    hGenPt->SetMarkerSize(1.0);
    hGenPt->Draw("E1");

    // Save the histogram to a file
    TFile* outFile = new TFile("GenPtSpectra_pythiaSignal.root", "RECREATE");
    hGenPt->Write();
    outFile->Close();
}

void makePtSpectra_pythiaSignal_141X(int iState = 1) {
	// Read GenOnly Nofilter file with polarization weights
	// const char* filename = Form("../Files/OniaTree_Y%dS_GENONLY_NoFilter.root", iState);
    const char* filename = Form("../Files/Oniatree_Upsilon%dS_5p02TeV_TuneCP5_141X_GenOnly_merge.root", iState);

	TFile* file = TFile::Open(filename, "READ");
	if (!file) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
		return;
	}

	cout << "File " << filename << " opened" << endl;

	TTree* OniaTree = (TTree*)file->Get("hionia/myTree");

	/// OniaTree variables, quite old version (since this genonly file from 2015)
	// Int_t Gen_QQ_size;
    
    Short_t Gen_QQ_size;
	TClonesArray* Gen_QQ_4mom = nullptr;
    
    TClonesArray* Gen_mu_4mom = nullptr;

    TClonesArray* CloneArr_QQ = nullptr;
    TClonesArray* CloneArr_mu = nullptr;

    Short_t Gen_QQ_mupl_idx[1000];
    Short_t Gen_QQ_mumi_idx[1000];

	OniaTree->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size);
	OniaTree->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom);

	OniaTree->SetBranchAddress("Gen_mu_4mom", &Gen_mu_4mom);

    OniaTree->SetBranchAddress("Gen_QQ_mumi_idx", &Gen_QQ_mumi_idx);
    OniaTree->SetBranchAddress("Gen_QQ_mupl_idx", &Gen_QQ_mupl_idx);

	TLorentzVector* gen_QQ_LV = new TLorentzVector();
	TLorentzVector* gen_mumi_LV = new TLorentzVector();
	TLorentzVector* gen_mupl_LV = new TLorentzVector();

    Long64_t totEntries = OniaTree->GetEntries(); 
    
    // const int nBins = 6;
    // float ptBinning[nBins + 1] = {0, 2, 4, 6, 9, 12, 30};

    double ptBinning[] = {0.0, 2.0, 6.0, 12.0, 30.0};
    int nBins = sizeof(ptBinning)/sizeof(ptBinning[0]) - 1;

    TH1D* hGenPt = new TH1D("hGenPt", "", nBins, ptBinning);
    hGenPt->Sumw2();

 	// Loop over the events
    for (Long64_t iEvent = 0; iEvent < (totEntries); iEvent++) {
		if (iEvent % 10000 == 0) {
			cout << Form("\rProcessing event %lld / %lld (%.0f%%)", iEvent, totEntries, 100. * iEvent / totEntries) << flush;
		}

		OniaTree->GetEntry(iEvent);

		// loop over all gen upsilons
		for (int iGen = 0; iGen < Gen_QQ_size; iGen++) {
            gen_QQ_LV = (TLorentzVector*)Gen_QQ_4mom->At(iGen);

            if (fabs(gen_QQ_LV->Rapidity()) < gRapidityMin || fabs(gen_QQ_LV->Rapidity()) > gRapidityMax) continue;

            hGenPt->Fill(gen_QQ_LV->Pt());

        }
    }

    TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
    hGenPt->SetTitle("");
    hGenPt->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    hGenPt->GetYaxis()->SetTitle("Counts");
    hGenPt->SetLineColor(kBlue);
    hGenPt->SetMarkerStyle(20);
    hGenPt->SetMarkerColor(kBlue);
    hGenPt->SetMarkerSize(1.0);
    hGenPt->Draw("E1");

    // Save the histogram to a file
    TFile* outFile = new TFile("GenPtSpectra_pythiaSignal_141X.root", "RECREATE");
    hGenPt->Write();
    outFile->Close();
}

TH1D* convertGraphToHist(TGraphErrors* hppXSecGraph, double& totalXS) {

    double ptEdges[] = {0.0, 2.0, 4.0, 6.0, 9.0, 12.0, 30.0};
    int nBins = sizeof(ptEdges)/sizeof(ptEdges[0]) - 1;

    double ptEdges_combined[] = {0.0, 2.0, 6.0, 12.0, 30.0};
    int nBins_combined = sizeof(ptEdges_combined)/sizeof(ptEdges_combined[0]) - 1;

    TH1D* hGraphAsHist = new TH1D("hGraphAsHist", "", nBins_combined, ptEdges_combined);
    hGraphAsHist->Sumw2();

    for (int i = 0; i < nBins_combined; ++i) {
        double x, y;
        double yErr;
        double ySum = 0, yErrSum = 0;
        
        int n = 0;

        cout << "i: " << i << ", n: " << n << endl;
        cout << "ptEdges[i + n + 1]: " << ptEdges[i + n + 1] << ", ptEdges_combined[i + 1]: " << ptEdges_combined[i + 1] << endl;
        // if (ptEdges[i + n + 1] == ptEdges_combined[i + 1]) n++;
        while (true) {
            if (n == 0) {
                hppXSecGraph->GetPoint(i + n, x, y);
                ySum = y;
                yErr = hppXSecGraph->GetErrorY(i + n);
                yErrSum = yErr;
            }
            else {
                hppXSecGraph->GetPoint(i + n, x, y);
                yErr = hppXSecGraph->GetErrorY(i + n);
                ySum = (ySum / pow(yErrSum, 2) + y / pow(yErr, 2)) /  (1. / pow(yErrSum, 2) + 1. / pow(yErr, 2));
                yErrSum = sqrt(1. / (1. / pow(yErrSum, 2) + 1. / pow(yErr, 2)));
            }
            
            double width = ptEdges[i + n + 1] - ptEdges[i + n];
            totalXS += y * width;

            cout << "ptEdges[i + n + 1]: " << ptEdges[i + n + 1] << ", ptEdges_combined[i + 1]: " << ptEdges_combined[i + 1] << endl;
            cout << "i: " << i << ", n: " << n << ", x: " << x << ", y: " << y << ", yErr: " << yErr << endl;
            cout << "ySum: " << ySum << ", yErrSum: " << yErrSum << endl;
            cout << "totalXS: " << totalXS << endl;

            if (ptEdges[i + n + 1] == ptEdges_combined[i + 1]) break;
            n++;
        }

        cout << "ptEdges[i + n + 1]: " << ptEdges[i + n + 1] << ", ptEdges_combined[i + 1]: " << ptEdges_combined[i + 1] << endl;
        cout << "ySum: " << ySum << ", yErrSum: " << yErrSum << endl;
        cout << " " << endl;

        hGraphAsHist->SetBinContent(i+1, ySum);
        hGraphAsHist->SetBinError(i+1, yErrSum);

    }

    return hGraphAsHist;
}

void errorPropagation(TH1D* hratio, TH1D* hnum, TH1D* hden) {
    // Error propagation for the ratio
    int nBins = hratio->GetNbinsX();

    for (int i = 1; i <= nBins; ++i) {
        double numContent = hnum->GetBinContent(i);
        double numError = hnum->GetBinError(i);
        double denContent = hden->GetBinContent(i);
        double denError = hden->GetBinError(i);

        double ratio = (denContent !=0) ? numContent / denContent : 0;
        double ratioError = 0;
        
        if (numContent != 0 && denContent != 0) {
            cout << "i: " << i << ", ratio: " << ratio << endl;
            cout << "numContent: " << numContent << ", numError: " << numError << ", denContent: " << denContent << ", denError: " << denError << endl;
            ratioError = ratio * std::sqrt(std::pow(numError / numContent, 2) + std::pow(denError / denContent, 2));
            // ratioError = std::sqrt(std::pow(numError / numContent, 2) + std::pow(denError / denContent, 2));
            cout << "ratioError: " << ratioError << endl;
        }
        else {
            ratioError = 0;
        }

        hratio->SetBinError(i, ratioError);
    }

}

void nGenSpectra() {
	writeExtraText = true; // if extra text
	extraText = "       Internal";

	// const char* filename = "GenPtSpectra_pythiaSignal.root";
	const char* filename = "GenPtSpectra_pythiaSignal_141X.root";

	TFile* file = TFile::Open(filename, "READ");
	if (!file) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
		return;
	}

	cout << "File " << filename << " opened" << endl;

    TH1D* hGenPt = (TH1D*)file->Get("hGenPt");
    if (!hGenPt) {
        std::cerr << "Error: Could not find histogram 'hGenPt'!" << std::endl;
        file->Close();
        return;
    }

	TTree* OniaTree = (TTree*)file->Get("hionia/myTree");

    /// open the published HEPData file (https://www.hepdata.net/record/ins1674529?version=2)
	TFile* HEPDataFile = TFile::Open("./HEPData-ins1674529-v2-Table_1.root", "READ");
    if (!HEPDataFile || HEPDataFile->IsZombie()) {
        std::cerr << "Error: Could not open file!" << std::endl;
        return;
    }

    TDirectory* dir = (TDirectory*)HEPDataFile->Get("Table 1");
    if (!dir) {
        std::cerr << "Error: Could not find directory 'Table 1'!" << std::endl;
        HEPDataFile->Close();
        return;
    }

    TGraphErrors* hppXSecGraph = (TGraphErrors*)dir->Get("Graph1D_y1");
	if (!hppXSecGraph) {
		std::cerr << "Error: Could not find graph 'Graph1D_y1'!" << std::endl;
		HEPDataFile->Close();
		return;
	}

	gStyle->SetPadLeftMargin(.18);

	auto canvas = new TCanvas("canvas", "", 600, 650);
	TPad* pad1 = new TPad("pad1", "pad1", 0, 0.4, 1, 1.0);
	pad1->SetBottomMargin(0.04);
	pad1->SetTopMargin(0.1);
	pad1->Draw();
	pad1->cd();

    hGenPt->SetTitle("");
    hGenPt->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    hGenPt->GetYaxis()->SetTitle("d#sigma^{2} / dydp_{T} (nb/GeV/c)");
    hGenPt->GetXaxis()->SetLabelSize(0);
    hGenPt->SetLineColor(kBlue);
    hGenPt->SetMarkerStyle(20);
    hGenPt->SetMarkerColor(kBlue);
    hGenPt->SetMarkerSize(1.0);
    hGenPt->Sumw2();
    
    double totalXS = 0;
    TH1D* hppXSecHist = convertGraphToHist(hppXSecGraph, totalXS);

    hGenPt->Scale(totalXS / hGenPt->Integral("width"));

    hGenPt->GetYaxis()->SetRangeUser(0, std::max(hppXSecHist->GetMaximum(), hGenPt->GetMaximum()) * 1.2);
    hGenPt->Draw("E1");

	hppXSecHist->Draw("SAME P");

	TPaveText* header = new TPaveText(.3, .85, .9, .75, "NDCNB");
	header->SetFillColor(4000);
	header->SetBorderSize(0);
	header->SetTextSize(.07);
	header->AddText(Form("%s, %s", CentralityRangeText(), DimuonRapidityRangeText(0., 2.4)));
	header->SetAllWith("", "align", 12);
	header->Draw();

	TLegend* legend = new TLegend(.5, .7, .85, .55);
	legend->SetTextSize(.07);
	legend->AddEntry(hppXSecHist, "Published pp XS", "lep");
    legend->AddEntry(hGenPt, "Gen MC", "lep");

	legend->Draw();

	// ratio plot
	canvas->cd();

	TPad* pad2 = new TPad("bottomPad", "bottomPad", 0, 0.0, 1, .4);
	pad2->SetTopMargin(0.025);
	pad2->SetBottomMargin(0.3);
	pad2->SetTicks(1, 1);
	pad2->Draw();
	pad2->cd();

    TH1D* hRatio = (TH1D*)hppXSecHist->Clone("hRatio");
    hRatio->Divide(hGenPt);

    hRatio->SetTitle(" ");
    hRatio->GetYaxis()->SetTitleOffset(0.65);
    hRatio->GetYaxis()->SetTitle("Data / MC");
    hRatio->GetYaxis()->SetTitleSize(0.13);
    hRatio->GetYaxis()->SetLabelSize(0.11);
    hRatio->GetYaxis()->CenterTitle();
    hRatio->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    hRatio->GetXaxis()->SetLabelSize(0.11);
    hRatio->GetXaxis()->SetTitleSize(0.13);
    hRatio->GetXaxis()->SetTickSize(0.06);
    hRatio->GetYaxis()->SetNdivisions(505);
    hRatio->Draw("E1");
    // hRatio->SetMinimum(0.1);
    // hRatio->SetMaximum(2.3);
    hRatio->SetMarkerStyle(20);

    errorPropagation(hRatio, hppXSecHist, hGenPt);

	TF1* fitFunc = new TF1("fitFunc", "[0]/([1] + x)", gPtFineBinning[0], gPtFineBinning[NPtFineBins]);
    // TF1* fitFunc = new TF1("fitFunc", "[0]/([1] + x)", 0, 30);
    // TF1* fitFunc = new TF1("fitFunc", "([0]  + [1]*x*x) / ( x - [2])^3", gPtFineBinning[0], gPtFineBinning[NPtFineBins]);
    // fitFunc->SetParameters(11.0, 6.0); 
    // fitFunc->FixParameter(0, 11.0);
    // fitFunc->FixParameter(1, 5.8);

    // TF1* fitFunc = new TF1("fitFunc", "[0]/pow(1 + (x / [1]), [2])", 0, 30);
    // fitFunc->SetParameters(2.5, 2.5, 2.2);  // A, B, C
    
    // TF1* fitFunc = new TF1("fitFunc", "[0]/pow(([1] + x), [2])", gPtFineBinning[0], gPtFineBinning[NPtFineBins]);
    // fitFunc->SetParameters(10.0, 0.1, 5.0); 

    // TF1* fitFunc = new TF1("fitFunc", "([0] + [1] * pow(x, 2.))/pow((x - [2]), 3)", gPtFineBinning[0], gPtFineBinning[NPtFineBins]);

    // fitFunc->SetParLimits(0, 1.0, 5.0);  // A
    // fitFunc->SetParLimits(1, 0.1, 10.0); // B
    // fitFunc->SetParLimits(2, 0.5, 5.0);  // C

    // Optionally set limits to keep it physically reasonable
    // fitFunc->SetParLimits(0, 7, 20);     // A
    // fitFunc->SetParLimits(1, 4, 10);   // B
    // fitFunc->SetParLimits(2, 0.1, 5.0);   // C

    auto fitResult = hRatio->Fit(fitFunc, "QEMSR", "", gPtFineBinning[0], gPtFineBinning[NPtFineBins]);
    fitResult->Print("v");

  	// legend with fit result info

	TLegend* fitLegend = new TLegend(.4, .9, .7, .65);
	fitLegend->SetTextSize(.09);
	fitLegend->AddEntry(fitFunc, Form("#frac{A}{(p_{T} + B)}  fit (#chi^{2} / n_{dof} = %.1f / %d)", fitResult->Chi2(), fitResult->Ndf()), "l");
		// fitLegend->AddEntry(fitFunc, Form("#frac{A + B p_{T}^{2}}{(p_{T} - C)^{3}}  fit (#chi^{2} / n_{dof} = %.1f / %d)", fitResult->Chi2(), fitResult->Ndf()), "l");
    // fitLegend->AddEntry(fitFunc, Form("#frac{A}{(p_{T} + B)^{C}}  fit (#chi^{2} / n_{dof} = %.1f / %d)", fitResult->Chi2(), fitResult->Ndf()), "l");
	
    fitLegend->Draw();

	TPaveText* fitHeader = new TPaveText(.65, .7, .9, .45, "NDCNB");
	fitHeader->SetFillColor(4000);
	fitHeader->SetBorderSize(0);
	fitHeader->SetTextSize(.09);
	fitHeader->AddText(Form("A = %.1f #pm %.1f", fitFunc->GetParameter(0), fitFunc->GetParError(0)));
	fitHeader->AddText(Form("B = %.1f #pm %.1f", fitFunc->GetParameter(1), fitFunc->GetParError(1)));
	// fitHeader->AddText(Form("C = %.1f #pm %.1f", fitFunc->GetParameter(2), fitFunc->GetParError(2)));

	fitHeader->SetAllWith("", "align", 12);
	fitHeader->Draw();  

	canvas->Modified();
	canvas->Update();
	canvas->cd();

	pad1->Draw();

	pad2->Draw();

	CMS_lumi(pad1, gCMSLumiText);

}