#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/Parameters/PhaseSpace.h"

#include "../MonteCarlo/AccEffHelpers.h"

//#include "../Tools/Style/FitDistributions.h"
#include "../Tools/Style/Legends.h"

#include "../Tools/Parameters/CentralityValues.h"
#include "../Tools/Parameters/EfficiencyWeights.h"
#include "../Tools/Parameters/MuonScaleFactors.h"

#include "../ReferenceFrameTransform/Transformations.h"

/// skim the GenOnly file to get the pT spectra of the generated upsilons
void makePtSpectra_pythiaSignal(int iState = 1, TString muonAccName = "UpsilonTriggerThresholds") {
	// Read GenOnly Nofilter file with polarization weights
	const char* filename = Form("../Files/OniaTree_Y%dS_GENONLY_NoFilter.root", iState);

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

	Bool_t withinAcceptance;
    
	Long64_t totEntries = OniaTree->GetEntries(); 
    
    // const int nBins = 6;
    // float ptBinning[nBins + 1] = {0, 2, 4, 6, 9, 12, 30};

    // double ptBinning[] = {0.0, 2.0, 6.0, 12.0, 30.0};
    double ptBinning[] = {0.0, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 12.0, 16.0, 20.0, 30.0};
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

			// single-muon acceptance cuts
			gen_mumi_LV = (TLorentzVector*)Gen_QQ_mumi_4mom->At(iGen);
			gen_mupl_LV = (TLorentzVector*)Gen_QQ_mupl_4mom->At(iGen);

			withinAcceptance = MuonKinematicsWithinLimits(*gen_mupl_LV, muonAccName) && MuonKinematicsWithinLimits(*gen_mumi_LV, muonAccName);
			if (!withinAcceptance) continue;

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
    TFile* outFile = new TFile(Form("GenPtSpectra_pythiaSignal_finerBin_%s.root", muonAccName.Data()), "RECREATE");
    hGenPt->Write();
    outFile->Close();
}

/// skim the GenOnly file to get the pT spectra of the generated upsilons
void makePtSpectra_pythiaSignal_141X(int iState = 1, TString muonAccName = "UpsilonTriggerThresholds") {
	// Read GenOnly Nofilter file with polarization weights
    // const char* filename = Form("../Files/Oniatree_Upsilon%dS_5p02TeV_TuneCP5_141X_GenOnly_merge.root", iState);
	// const char* filename = Form("../Files/Oniatree_Upsilon%dS_5p02TeV_TuneCP5_141X_GenOnly_combined.root", iState);
	const char* filename = Form("../Files/Oniatree_Upsilon%dS_5p02TeV_TuneCP5_141X_GenOnly_50M.root", iState);

	TFile* file = TFile::Open(filename, "READ");
	if (!file) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
		return;
	}

	cout << "File " << filename << " opened" << endl;

	TTree* OniaTree = (TTree*)file->Get("hionia/myTree");

	/// OniaTree variables, quite old version (since this genonly file from 2015)
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

    // double ptBinning[] = {0.0, 2.0, 6.0, 12.0, 30.0};
    double ptBinning[] = {0.0, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 12.0, 16.0, 20.0, 30.0};
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
            
			// single-muon acceptance
			// positive muon first
			gen_mupl_LV = (TLorentzVector*)Gen_mu_4mom->At(Gen_QQ_mupl_idx[iGen]);

			if (!MuonKinematicsWithinLimits(*gen_mupl_LV, muonAccName)) continue;

			// then negative muon
			gen_mumi_LV = (TLorentzVector*)Gen_mu_4mom->At(Gen_QQ_mumi_idx[iGen]);

			if (!MuonKinematicsWithinLimits(*gen_mumi_LV, muonAccName)) continue;

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
    TFile* outFile = new TFile(Form("GenPtSpectra_pythiaSignal_141X_finerBin_%s.root", muonAccName.Data()), "RECREATE");
    hGenPt->Write();
    outFile->Close();
}

/// skim the Hydjet Gen file to get the pT spectra of the generated upsilons
void makePtSpectra_HydjetGen(int iState = 1, TString muonAccName = "UpsilonTriggerThresholds") {
    /// Read Hydjet Gen file
	// const char* filename = Form("../Files/OniaTree_Y%dS_pThat2_HydjetDrumMB_miniAOD.root", iState);
	const char* filename = Form("../Files/Oniatree_Upsilon%dS_5p02TeV_TuneCP5_141X_GenOnly_Hydjet_merge.root", iState);
	TFile* file = TFile::Open(filename, "READ");
	if (!file) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
		return;
	}

	cout << "File " << filename << " opened" << endl;

	writeExtraText = true; // if extra text
	extraText = "       Internal";

	TTree* OniaTree = (TTree*)file->Get("hionia/myTree");

    /// OniaTree variables

	Float_t Gen_weight;

	TClonesArray* Gen_QQ_4mom = nullptr;
	TClonesArray* Gen_mu_4mom = nullptr;

	Float_t zVtx;
	ULong64_t HLTriggers;
	ULong64_t Reco_QQ_trig[1000];
	Int_t Centrality;
	Float_t HFmean;

	TClonesArray* CloneArr_QQ = nullptr;
	TClonesArray* CloneArr_mu = nullptr;
	Short_t Reco_QQ_size;
	Short_t Reco_QQ_sign[1000];
	Short_t Reco_QQ_mupl_idx[1000];
	Short_t Reco_QQ_mumi_idx[1000];
	Short_t Gen_QQ_whichRec[1000];

	Short_t Gen_QQ_mupl_idx[1000];
	Short_t Gen_QQ_mumi_idx[1000];

	ULong64_t Reco_mu_trig[1000];
	Int_t Reco_mu_SelectionType[1000];
	//(parameters for quality cuts)
	Float_t Reco_QQ_VtxProb[1000];
	Int_t Reco_mu_nPixWMea[1000];
	Int_t Reco_mu_nTrkWMea[1000];
	Float_t Reco_mu_dxy[1000];
	Float_t Reco_mu_dz[1000];

	Short_t Gen_QQ_size;

	// event variables
	OniaTree->SetBranchAddress("Gen_weight", &Gen_weight);
	OniaTree->SetBranchAddress("Centrality", &Centrality);
	OniaTree->SetBranchAddress("SumET_HF", &HFmean);
	OniaTree->SetBranchAddress("HLTriggers", &HLTriggers);
	OniaTree->SetBranchAddress("zVtx", &zVtx);

	// gen-level variables
	OniaTree->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size);
	OniaTree->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom);
	OniaTree->SetBranchAddress("Gen_QQ_whichRec", Gen_QQ_whichRec);

	OniaTree->SetBranchAddress("Gen_mu_4mom", &Gen_mu_4mom);
	OniaTree->SetBranchAddress("Gen_QQ_mupl_idx", Gen_QQ_mupl_idx);
	OniaTree->SetBranchAddress("Gen_QQ_mumi_idx", Gen_QQ_mumi_idx);

	OniaTree->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig);
	OniaTree->SetBranchAddress("Reco_QQ_4mom", &CloneArr_QQ);
	OniaTree->SetBranchAddress("Reco_mu_4mom", &CloneArr_mu);
	OniaTree->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size);
	OniaTree->SetBranchAddress("Reco_QQ_sign", &Reco_QQ_sign);
	OniaTree->SetBranchAddress("Reco_QQ_mupl_idx", Reco_QQ_mupl_idx);
	OniaTree->SetBranchAddress("Reco_QQ_mumi_idx", Reco_QQ_mumi_idx);

	OniaTree->SetBranchAddress("Reco_mu_trig", Reco_mu_trig);
	OniaTree->SetBranchAddress("Reco_mu_SelectionType", Reco_mu_SelectionType);

	OniaTree->SetBranchAddress("Reco_QQ_VtxProb", &Reco_QQ_VtxProb);
	OniaTree->SetBranchAddress("Reco_mu_nPixWMea", &Reco_mu_nPixWMea);
	OniaTree->SetBranchAddress("Reco_mu_nTrkWMea", &Reco_mu_nTrkWMea);
	OniaTree->SetBranchAddress("Reco_mu_dxy", &Reco_mu_dxy);
	OniaTree->SetBranchAddress("Reco_mu_dz", &Reco_mu_dz);

	// loop variables
	TLorentzVector* genLorentzVector = new TLorentzVector();
	TLorentzVector* gen_QQ_LV = new TLorentzVector();
	TLorentzVector* gen_mupl_LV = new TLorentzVector();
	TLorentzVector* gen_mumi_LV = new TLorentzVector();
	TLorentzVector* recoLorentzVector = new TLorentzVector();

    Int_t hiBin;
	
    Long64_t totEntries = OniaTree->GetEntries();

    double ptBinning[] = {0.0, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 12.0, 16.0, 20.0, 30.0};
    // double ptBinning[] = {0.0, 2.0, 6.0, 12.0, 30.0};
    int nBins = sizeof(ptBinning)/sizeof(ptBinning[0]) - 1;

    TH1D* hGenPtHydjet = new TH1D("hGenPtHydjet", "", nBins, ptBinning);
    hGenPtHydjet->Sumw2();

	for (Long64_t iEvent = 0; iEvent < totEntries; iEvent++) {
		if (iEvent % 10000 == 0) {
			cout << Form("\rProcessing event %lld / %lld (%.0f%%)", iEvent, totEntries, 100. * iEvent / totEntries) << flush;
		}
		// if (iEvent > 100) break; // for testing

		OniaTree->GetEntry(iEvent);

        hiBin = GetHiBinFromhiHF(HFmean);

		// loop over all gen upsilons
		for (int iGen = 0; iGen < Gen_QQ_size; iGen++) {
            gen_QQ_LV = (TLorentzVector*)Gen_QQ_4mom->At(iGen);
            double gen_QQ_pt = gen_QQ_LV->Pt();
            
            if (fabs(gen_QQ_LV->Rapidity()) < gRapidityMin || fabs(gen_QQ_LV->Rapidity()) > gRapidityMax) continue;
           
						// single-muon acceptance
			// positive muon first
			gen_mupl_LV = (TLorentzVector*)Gen_mu_4mom->At(Gen_QQ_mupl_idx[iGen]);

			if (!MuonKinematicsWithinLimits(*gen_mupl_LV, muonAccName)) continue;

			// then negative muon
			gen_mumi_LV = (TLorentzVector*)Gen_mu_4mom->At(Gen_QQ_mumi_idx[iGen]);

			if (!MuonKinematicsWithinLimits(*gen_mumi_LV, muonAccName)) continue;

			// hGenPtHydjet->Fill(gen_QQ_pt, Gen_weight * FindNcoll(hiBin));
			// hGenPtHydjet->Fill(gen_QQ_pt, Gen_weight);
			hGenPtHydjet->Fill(gen_QQ_pt);
            // hGenPtHydjet->Fill(gen_QQ_pt, FindNcoll(hiBin));
        }
    }

    TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
    hGenPtHydjet->SetTitle("");
    hGenPtHydjet->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    hGenPtHydjet->GetYaxis()->SetTitle("Counts");
    hGenPtHydjet->SetLineColor(kBlue);
    hGenPtHydjet->SetMarkerStyle(20);
    hGenPtHydjet->SetMarkerColor(kBlue);
    hGenPtHydjet->SetMarkerSize(1.0);
    hGenPtHydjet->Draw("E1");

    // Save the histogram to a file
    // TFile* outFile = new TFile(Form("GenPtSpectra_HydjetGen_finerBin_noGenWeight_%s.root", muonAccName.Data()), "RECREATE");
	// TFile* outFile = new TFile(Form("GenPtSpectra_HydjetGen_finerBin_noNColl_%s.root", muonAccName.Data()), "RECREATE");
	TFile* outFile = new TFile(Form("GenPtSpectra_HydjetGen_finerBin_noNColl_noGenWeight_141X_%s.root", muonAccName.Data()), "RECREATE");
    hGenPtHydjet->Write();
    outFile->Close();

}

void nGenSpectra_MC(TString muonAccName = "UpsilonTriggerThresholds") {
	writeExtraText = true; // if extra text
	extraText = "       Internal";

    /// open the Pythia signal gen pt spectra file
	const char* pythiaFilename = Form("GenPtSpectra_pythiaSignal_141X_finerBin_%s.root", muonAccName.Data());

	TFile* pythiaFile = TFile::Open(pythiaFilename, "READ");
	if (!pythiaFile) {
		cout << "File " << pythiaFilename << " not found. Check the directory of the file." << endl;
		return;
	}

	cout << "File " << pythiaFilename << " opened" << endl;

    TH1D* hGenPt = (TH1D*)pythiaFile->Get("hGenPt");
    if (!hGenPt) {
        std::cerr << "Error: Could not find histogram 'hGenPt'!" << std::endl;
        pythiaFile->Close();
        return;
    }

	/// open the Pythia embedded Hydjet GEN pt spectra file with gen weights
	const char* hydjetFilename = Form("GenPtSpectra_HydjetGen_finerBin_noNColl_%s.root", muonAccName.Data());
	// const char* hydjetFilename = Form("GenPtSpectra_HydjetGen_finerBin_noNColl_141X_%s.root", muonAccName.Data());

	TFile* hydjetFile = TFile::Open(hydjetFilename, "READ");
	if (!hydjetFile) {
		cout << "File " << hydjetFilename << " not found. Check the directory of the file." << endl;
		return;
	}

	cout << "File " << hydjetFilename << " opened" << endl;

    // TH1D* hGenPt = (TH1D*)file->Get("hGenPt");
	TH1D* hGenPtHydjet_wWeight = (TH1D*)hydjetFile->Get("hGenPtHydjet");
    if (!hGenPtHydjet_wWeight) {
        std::cerr << "Error: Could not find histogram 'hGenPt'!" << std::endl;
        hydjetFile->Close();
        return;
    }

    // /// open the Hydjet gen pt spectra file
    // const char* hydjetFilename2 = Form("GenPtSpectra_HydjetGen_finerBin_noGenWeight_%s.root", muonAccName.Data());
	// // const char* hydjetFilename2 = Form("GenPtSpectra_HydjetGen_finerBin_noNColl_noGenWeight_141X_%s.root", muonAccName.Data());

    // TFile* hydjetFile2 = TFile::Open(hydjetFilename2, "READ");

    // if (!hydjetFile2 || hydjetFile2->IsZombie()) {
    //     std::cerr << "Error: Could not open file!" << std::endl;
    //     return;
    // }

    // cout << "File " << hydjetFilename2 << " opened" << endl;

    // TH1D* hGenPtHydjet_woWeight = (TH1D*)hydjetFile2->Get("hGenPtHydjet");

    // if (!hGenPtHydjet_woWeight) {
    //     std::cerr << "Error: Could not find histogram 'hGenPtHydjet'!" << std::endl;
    //     hydjetFile2->Close();
    //     return;
    // }

	gStyle->SetPadLeftMargin(.18);

	auto canvas = new TCanvas("canvas", "", 600, 650);
	TPad* pad1 = new TPad("pad1", "pad1", 0, 0.4, 1, 1.0);
	pad1->SetBottomMargin(0.04);
	pad1->SetTopMargin(0.1);
	pad1->Draw();
	pad1->cd();

    hGenPt->SetTitle("");
    hGenPt->GetXaxis()->SetTitle(gPtAxisTitle);
    // hGenPt->GetYaxis()->SetTitle("d#sigma^{2} / dydp_{T} (nb/GeV/c)");
	hGenPt->GetYaxis()->SetTitle(Form("dN(#varUpsilon(1S)) / d%s (%s)^{-1}", gPtVarName, gPtUnit));
    hGenPt->GetXaxis()->SetLabelSize(0);
    hGenPt->SetLineColor(kBlue);
    hGenPt->SetMarkerStyle(20);
    hGenPt->SetMarkerColor(kBlue);
    hGenPt->SetMarkerSize(1.0);
    hGenPt->Sumw2();
    
    hGenPt->Scale(1 / hGenPt->Integral());

    hGenPt->Draw("E1");

    hGenPtHydjet_wWeight->Scale(1 / hGenPtHydjet_wWeight->Integral());
    hGenPtHydjet_wWeight->SetLineColor(kBlack);
    hGenPtHydjet_wWeight->SetMarkerStyle(22);
    hGenPtHydjet_wWeight->SetMarkerColor(kBlack);
    hGenPtHydjet_wWeight->Draw("E1 same");

    // hGenPtHydjet_woWeight->Scale(1 / hGenPtHydjet_woWeight->Integral());
    // hGenPtHydjet_woWeight->SetLineColor(kRed);
    // hGenPtHydjet_woWeight->SetMarkerStyle(23);
    // hGenPtHydjet_woWeight->SetMarkerColor(kRed);
    // hGenPtHydjet_woWeight->Draw("E1 same");

    hGenPt->GetYaxis()->SetRangeUser(0, std::max(hGenPtHydjet_wWeight->GetMaximum(), hGenPt->GetMaximum()) * 1.2);

	TPaveText* header = new TPaveText(.3, .85, .9, .75, "NDCNB");
	header->SetFillColor(4000);
	header->SetBorderSize(0);
	header->SetTextSize(.07);
	header->AddText(Form("%s, %s", CentralityRangeText(), DimuonRapidityRangeText(0., 2.4)));
	header->SetAllWith("", "align", 12);
	header->Draw();

	TLegend* legend = new TLegend(.5, .7, .85, .55);
	legend->SetTextSize(.05);
	// legend->AddEntry(hGenPtHydjet, "Hydjet Gen MC", "lep");
    legend->AddEntry(hGenPt, "Pythia Signal Gen MC", "lep");
	// legend->AddEntry(hGenPtHydjet_woWeight, "Hydjet w/o GenWeight", "lep");
    legend->AddEntry(hGenPtHydjet_wWeight, "Hydjet with GenWeight", "lep");
	
	legend->Draw();

	// ratio plot
	canvas->cd();

	TPad* pad2 = new TPad("bottomPad", "bottomPad", 0, 0.0, 1, .4);
	pad2->SetTopMargin(0.025);
	pad2->SetBottomMargin(0.3);
	pad2->SetTicks(1, 1);
	pad2->Draw();
	pad2->cd();

    TH1D* hRatio = (TH1D*)hGenPt->Clone("hRatio");
    hRatio->Divide(hGenPtHydjet_wWeight);

    hRatio->SetTitle(" ");
    hRatio->GetYaxis()->SetTitleOffset(0.65);
    hRatio->GetYaxis()->SetTitle("Pythia / Hydjet");
    hRatio->GetYaxis()->SetTitleSize(0.11);
    hRatio->GetYaxis()->SetLabelSize(0.11);
    hRatio->GetYaxis()->CenterTitle();
    hRatio->GetXaxis()->SetTitle(gPtAxisTitle);
    hRatio->GetXaxis()->SetLabelSize(0.11);
    hRatio->GetXaxis()->SetTitleSize(0.13);
    hRatio->GetXaxis()->SetTickSize(0.06);
    hRatio->GetYaxis()->SetNdivisions(505);
    hRatio->SetLineColor(kBlack);
    hRatio->GetYaxis()->SetRangeUser(0., 2.3);   

    hRatio->Draw("E1");
    // hRatio->SetMinimum(0.1);
    // hRatio->SetMaximum(2.3);
    hRatio->SetMarkerStyle(20); 
    // errorPropagation(hRatio, hppXSecHist, hGenPt);

	// TF1* fitFunc = new TF1("fitFunc", "[0]/([1] + x)", gPtFineBinning[0], gPtFineBinning[NPtFineBins]);
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
    
    // // TF1* fitFunc = new TF1("fitFunc", "[0] + [1] / pow((x + [2]), [3])", gPtFineBinning[0], gPtFineBinning[NPtFineBins]);
   
    // // TF1* fitFunc = new TF1("fitFunc", "[0] * exp(-x/[1])/(x + [2])", gPtFineBinning[0], gPtFineBinning[NPtFineBins]);

    // // fitFunc->SetParLimits(0, 1.0, 5.0);  // A
    // // fitFunc->SetParLimits(1, 0.1, 10.0); // B
    // // fitFunc->SetParLimits(2, 0.5, 5.0);  // C

    // // Optionally set limits to keep it physically reasonable
    // // fitFunc->SetParLimits(0, 7, 20);     // A
    // // fitFunc->SetParLimits(1, 4, 10);   // B
    // // fitFunc->SetParLimits(2, 0.1, 5.0);   // C

    // auto fitResult = hRatio->Fit(fitFunc, "QEMSR", "", gPtFineBinning[0], gPtFineBinning[NPtFineBins]);
    // fitResult->Print("v");

  	// // legend with fit result info

	// TLegend* fitLegend = new TLegend(.4, .9, .7, .65);
	// fitLegend->SetTextSize(.09);
	// fitLegend->AddEntry(fitFunc, Form("#frac{A}{(p_{T} + B)}  fit (#chi^{2} / n_{dof} = %.1f / %d)", fitResult->Chi2(), fitResult->Ndf()), "l");
	// 	// fitLegend->AddEntry(fitFunc, Form("#frac{A + B p_{T}^{2}}{(p_{T} - C)^{3}}  fit (#chi^{2} / n_{dof} = %.1f / %d)", fitResult->Chi2(), fitResult->Ndf()), "l");
    // // fitLegend->AddEntry(fitFunc, Form("#frac{A}{(p_{T} + B)^{C}}  fit (#chi^{2} / n_{dof} = %.1f / %d)", fitResult->Chi2(), fitResult->Ndf()), "l");
	
    // fitLegend->Draw();

	// TPaveText* fitHeader = new TPaveText(.65, .7, .9, .45, "NDCNB");
	// fitHeader->SetFillColor(4000);
	// fitHeader->SetBorderSize(0);
	// fitHeader->SetTextSize(.09);
	// fitHeader->AddText(Form("A = %.1f #pm %.1f", fitFunc->GetParameter(0), fitFunc->GetParError(0)));
	// fitHeader->AddText(Form("B = %.1f #pm %.1f", fitFunc->GetParameter(1), fitFunc->GetParError(1)));
	// // fitHeader->AddText(Form("C = %.1f #pm %.1f", fitFunc->GetParameter(2), fitFunc->GetParError(2)));

	// fitHeader->SetAllWith("", "align", 12);
	// fitHeader->Draw();  

	canvas->Modified();
	canvas->Update();
	canvas->cd();

	pad1->Draw();

	pad2->Draw();

	CMS_lumi(pad1, gCMSLumiText);

}