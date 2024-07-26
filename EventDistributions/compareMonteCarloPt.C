#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../MonteCarlo/AccEffHelpers.h"

#include "../Tools/Style/Legends.h"

#include "../ReferenceFrameTransform/Transformations.h"

TEfficiency* getNoFilterGenPt(Int_t iState = 1){

	// Read GenOnly Nofilter file
	const char* filename = Form("../Files/OniaTree_Y%dS_GENONLY_NoFilter.root", iState);	

	TFile* file = TFile::Open(filename, "READ");
	if (!file) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
		exit(1);
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

	TEfficiency* hAccPt = new TEfficiency("", ";p^{#varUpsilon}_{T} (GeV/c);Acceptance of #varUpsilon(1S)", 50, 0, 50);

	TLorentzVector* gen_QQ_LV = new TLorentzVector();
	TLorentzVector* gen_mumi_LV = new TLorentzVector();
	TLorentzVector* gen_mupl_LV = new TLorentzVector();
	
	Bool_t withinAcceptance;

	Long64_t totEntries = OniaTree->GetEntries();

	// Loop over the events
	for (Long64_t iEvent = 0; iEvent < (totEntries); iEvent++) {
		if (iEvent % 10000 == 0) {
			cout << Form("\rProcessing event %lld / %lld (%.0f%%)", iEvent, totEntries, 100. * iEvent / totEntries) << flush;
		}

		OniaTree->GetEntry(iEvent);

		// if (iEvent > 50000) break;

		// loop over all gen upsilons
		for (int iGen = 0; iGen < Gen_QQ_size; iGen++) {
			gen_QQ_LV = (TLorentzVector*)Gen_QQ_4mom->At(iGen);

			if (fabs(gen_QQ_LV->Rapidity()) < gRapidityMin || fabs(gen_QQ_LV->Rapidity()) > gRapidityMax) continue; // upsilon within fiducial region

			// single-muon acceptance cuts
			gen_mumi_LV = (TLorentzVector*)Gen_QQ_mumi_4mom->At(iGen);
			gen_mupl_LV = (TLorentzVector*)Gen_QQ_mupl_4mom->At(iGen);

			withinAcceptance = (fabs(gen_mupl_LV->Eta()) < 2.4) && (gen_mupl_LV->Pt() > gMuonPtCut) && (fabs(gen_mumi_LV->Eta()) < 2.4) && (gen_mumi_LV->Pt() > gMuonPtCut);

			hAccPt->Fill(withinAcceptance, gen_QQ_LV->Pt());
		}

	}
	
	return hAccPt;
}

TEfficiency* getHydjetGenPt(Int_t iState = 1){

	// Read Gen from Hydjet MC
	const char* filename = Form("../Files/OniaTree_Y%dS_pThat2_HydjetDrumMB_miniAOD.root", iState);	

	TFile* file = TFile::Open(filename, "READ");
	if (!file) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
		exit(1);
	}

	cout << "File " << filename << " opened" << endl;

	TTree* OniaTree = (TTree*)file->Get("hionia/myTree");

	/// OniaTree variables

	Float_t Gen_weight;

	TClonesArray* Gen_QQ_4mom = nullptr;
	TClonesArray* Gen_mu_4mom = nullptr;

	Short_t Gen_QQ_whichRec[1000];

	Short_t Gen_QQ_mupl_idx[1000];
	Short_t Gen_QQ_mumi_idx[1000];

	Short_t Gen_QQ_size;

	// event variables
	OniaTree->SetBranchAddress("Gen_weight", &Gen_weight);

	// gen-level variables
	OniaTree->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size);
	OniaTree->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom);
	OniaTree->SetBranchAddress("Gen_QQ_whichRec", Gen_QQ_whichRec);

	OniaTree->SetBranchAddress("Gen_mu_4mom", &Gen_mu_4mom);
	OniaTree->SetBranchAddress("Gen_QQ_mupl_idx", Gen_QQ_mupl_idx);
	OniaTree->SetBranchAddress("Gen_QQ_mumi_idx", Gen_QQ_mumi_idx);

	TEfficiency* hAccPt = new TEfficiency("", ";p^{#varUpsilon}_{T} (GeV/c);Acceptance of #varUpsilon(1S)", 50, 0, 50);

	// loop variables
	TLorentzVector* gen_QQ_LV = new TLorentzVector();
	TLorentzVector* gen_mumi_LV = new TLorentzVector();
	TLorentzVector* gen_mupl_LV = new TLorentzVector();

	Bool_t withinAcceptance;

	Long64_t totEntries = OniaTree->GetEntries();

	// Loop over the events
	for (Long64_t iEvent = 0; iEvent < (totEntries); iEvent++) {
		if (iEvent % 10000 == 0) {
			cout << Form("\rProcessing event %lld / %lld (%.0f%%)", iEvent, totEntries, 100. * iEvent / totEntries) << flush;
		}

		// if (iEvent > 50000) break;

		OniaTree->GetEntry(iEvent);

		// loop over all gen upsilons
		for (int iGen = 0; iGen < Gen_QQ_size; iGen++) {
			gen_QQ_LV = (TLorentzVector*)Gen_QQ_4mom->At(iGen);

			if (fabs(gen_QQ_LV->Rapidity()) < gRapidityMin || fabs(gen_QQ_LV->Rapidity()) > gRapidityMax) continue; // upsilon within fiducial region

			// single-muon acceptance cuts
			gen_mupl_LV = (TLorentzVector*)Gen_mu_4mom->At(Gen_QQ_mupl_idx[iGen]);
			gen_mumi_LV = (TLorentzVector*)Gen_mu_4mom->At(Gen_QQ_mumi_idx[iGen]);

			withinAcceptance = (fabs(gen_mupl_LV->Eta()) < 2.4) && (gen_mupl_LV->Pt() > gMuonPtCut) && (fabs(gen_mumi_LV->Eta()) < 2.4) && (gen_mumi_LV->Pt() > gMuonPtCut);

			hAccPt->Fill(withinAcceptance, gen_QQ_LV->Pt());
		}
	}

	return hAccPt;
}


void compareMonteCarloPt(Int_t iState = 1) {

	// Put the text, "CMS Internal", on the right top of the plot
	writeExtraText = true; // if extra text
	extraText = "       Internal";

	TEfficiency* hAccNoFilterGenPt = getNoFilterGenPt(iState);

	TEfficiency* hAccHydjetGenPt = getHydjetGenPt(iState);

	auto canvas = new TCanvas("canvas", "", 650, 650);

	canvas->SetLeftMargin(0.2);
	canvas->SetRightMargin(0.05);

	TH1D* hDummy = new TH1D("hDummy", ";p^{#varUpsilon}_{T} (GeV/c);Acceptance of #varUpsilon(1S)", 50, 0, 50);

	hDummy->GetXaxis()->SetTitleSize(0.075);

	hDummy->GetYaxis()->SetTitleOffset(1.2);
	hDummy->GetYaxis()->SetTitleSize(0.075);
	hDummy->GetYaxis()->SetLabelSize(0.065);

	hDummy->GetYaxis()->SetRangeUser(0, 1);

	hAccNoFilterGenPt->SetLineColor(kRed + 1);
	hAccNoFilterGenPt->SetMarkerStyle(8);
	hAccNoFilterGenPt->SetMarkerColor(kRed + 1);
	hAccNoFilterGenPt->SetLineWidth(2);

	hAccHydjetGenPt->SetLineColor(kAzure + 1);
	hAccHydjetGenPt->SetMarkerStyle(8);
	hAccHydjetGenPt->SetMarkerColor(kAzure + 1);
	hAccHydjetGenPt->SetLineWidth(2);

	hDummy->Draw();

	hAccNoFilterGenPt->Draw("P SAME");

	hAccHydjetGenPt->Draw("P SAME");

	TLegend* legend = new TLegend(.22, .9, .52, .8);
	legend->SetTextSize(.05);
	legend->AddEntry(hAccHydjetGenPt, "Hydjet Gen", "lep");
	legend->AddEntry(hAccNoFilterGenPt, "No Filter Gen", "lep");

	legend->Draw();

	gSystem->mkdir("plots", kTRUE);
	canvas->SaveAs("plots/MCPtAcceptance.png", "RECREATE");
}