#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../MonteCarlo/AccEffHelpers.h"

#include "../Tools/Style/Legends.h"

#include "../ReferenceFrameTransform/Transformations.h"

TEfficiency* getNoFilterGenRapidity(Int_t iState = 1){

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

	TEfficiency* hAccRapidity = new TEfficiency("", ";y^{#varUpsilon};Acceptance of #varUpsilon(1S)", 50, -3.5, 3.5);

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

			hAccRapidity->Fill(withinAcceptance, gen_QQ_LV->Rapidity());
		}

	}
	
	return hAccRapidity;
}

TEfficiency* getHydjetGenRapidity(Int_t iState = 1){

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

	TEfficiency* hAccRapidity = new TEfficiency("", ";y^{#varUpsilon};Acceptance of #varUpsilon(1S)", 50, -3.5, 3.5);

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

			hAccRapidity->Fill(withinAcceptance, gen_QQ_LV->Rapidity());
		}
	}

	return hAccRapidity;
}


void compareMonteCarloRapidity(Int_t iState = 1) {

	// Put the text, "CMS Internal", on the right top of the plot
	writeExtraText = true; // if extra text
	extraText = "       Internal";

	TEfficiency* hAccNoFilterGenRapidity = getNoFilterGenRapidity(iState);

	TEfficiency* hAccHydjetGenRapidity = getHydjetGenRapidity(iState);

	TH1D* hPassedNoFilterGenRapidity = (TH1D*) hAccNoFilterGenRapidity->GetPassedHistogram();
	TH1D* hTotalNoFilterGenRapidity = (TH1D*) hAccNoFilterGenRapidity->GetTotalHistogram();

	TH1D* hPassedHydjetGenRapidity = (TH1D*) hAccHydjetGenRapidity->GetPassedHistogram();
	TH1D* hTotalHydjetGenRapidity = (TH1D*) hAccHydjetGenRapidity->GetTotalHistogram();

	TH1D* hPassedNoFilterGenRapidityCopy = (TH1D*) hPassedNoFilterGenRapidity->Clone();
	TH1D* hTotalNoFilterGenRapidityCopy = (TH1D*) hTotalNoFilterGenRapidity->Clone();

	TH1D* hPassedHydjetGenRapidityCopy = (TH1D*) hPassedHydjetGenRapidity->Clone();
	TH1D* hTotalHydjetGenRapidityCopy = (TH1D*) hTotalHydjetGenRapidity->Clone();

	auto canvas = new TCanvas("canvas", "", 1200, 600);

	canvas->Divide(2);
	canvas->cd(1);

	gPad->SetLeftMargin(0.2);
	gPad->SetRightMargin(0.015);
	gPad->SetTopMargin(0.07);

	TH1D* hAccDummy = new TH1D("hAccDummy", ";y^{#varUpsilon};Acceptance of #varUpsilon(1S)", 50, -3., 3.);

	hAccDummy->GetXaxis()->SetTitleSize(0.075);

	hAccDummy->GetYaxis()->SetTitleOffset(1.2);
	hAccDummy->GetYaxis()->SetTitleSize(0.075);
	hAccDummy->GetYaxis()->SetLabelSize(0.065);

	hAccDummy->GetYaxis()->SetRangeUser(0, 1);

	hAccNoFilterGenRapidity->SetLineColor(kRed + 1);
	hAccNoFilterGenRapidity->SetMarkerStyle(8);
	hAccNoFilterGenRapidity->SetMarkerColor(kRed + 1);
	hAccNoFilterGenRapidity->SetLineWidth(2);

	hAccHydjetGenRapidity->SetLineColor(kAzure + 1);
	hAccHydjetGenRapidity->SetMarkerStyle(8);
	hAccHydjetGenRapidity->SetMarkerColor(kAzure + 1);
	hAccHydjetGenRapidity->SetLineWidth(2);

	hAccDummy->Draw();

	hAccNoFilterGenRapidity->Draw("P SAME");

	hAccHydjetGenRapidity->Draw("P SAME");

	TLegend* legend = new TLegend(.22, .9, .52, .8);
	legend->SetTextSize(.05);
	legend->AddEntry(hAccHydjetGenRapidity, "Hydjet Gen", "lep");
	legend->AddEntry(hAccNoFilterGenRapidity, "No Filter Gen", "lep");

	legend->Draw();

	canvas->cd(2);

	gPad->SetLeftMargin(0.2);
	gPad->SetRightMargin(0.015);
	gPad->SetTopMargin(0.07);

	hPassedNoFilterGenRapidityCopy->Scale(1. / hPassedNoFilterGenRapidity->Integral(), "Width");
	hTotalNoFilterGenRapidityCopy->Scale(1. / hTotalNoFilterGenRapidity->Integral(), "Width");

	hPassedHydjetGenRapidityCopy->Scale(1. / hPassedHydjetGenRapidity->Integral(), "Width");
	hTotalHydjetGenRapidityCopy->Scale(1. / hTotalHydjetGenRapidity->Integral(), "Width");


	TH1D* hDummy = new TH1D("hDummy", ";y^{#varUpsilon};dN(#varUpsilon(1S)) / y", 50, -3., 3.);

	hDummy->GetXaxis()->SetTitleSize(0.075);

	hDummy->GetYaxis()->SetTitleOffset(1.1);
	hDummy->GetYaxis()->SetTitleSize(0.075);
	hDummy->GetYaxis()->SetLabelSize(0.065);

	hDummy->GetYaxis()->SetRangeUser(0, hTotalHydjetGenRapidityCopy->GetMaximum() * 1.6);

	hDummy->GetYaxis()->SetMaxDigits(2);

	hTotalNoFilterGenRapidityCopy->SetLineColor(kRed + 1);
	hTotalNoFilterGenRapidityCopy->SetMarkerStyle(8);
	hTotalNoFilterGenRapidityCopy->SetMarkerColor(kRed + 1);
	hTotalNoFilterGenRapidityCopy->SetLineWidth(2);

	hPassedNoFilterGenRapidityCopy->SetLineColor(kPink - 4);
	hPassedNoFilterGenRapidityCopy->SetMarkerStyle(8);
	hPassedNoFilterGenRapidityCopy->SetMarkerColor(kPink - 4);
	hPassedNoFilterGenRapidityCopy->SetLineWidth(2);

	hTotalHydjetGenRapidityCopy->SetLineColor(kAzure + 1);
	hTotalHydjetGenRapidityCopy->SetMarkerStyle(8);
	hTotalHydjetGenRapidityCopy->SetMarkerColor(kAzure + 1);
	hTotalHydjetGenRapidityCopy->SetLineWidth(2);

	hPassedHydjetGenRapidityCopy->SetLineColor(kCyan - 3);
	hPassedHydjetGenRapidityCopy->SetMarkerStyle(8);
	hPassedHydjetGenRapidityCopy->SetMarkerColor(kCyan - 3);
	hPassedHydjetGenRapidityCopy->SetLineWidth(2);

	hDummy->Draw();

	hTotalNoFilterGenRapidityCopy->Draw("PE SAME");

	hPassedNoFilterGenRapidityCopy->Draw("PE SAME");

	hTotalHydjetGenRapidityCopy->Draw("PE SAME");

	hPassedHydjetGenRapidityCopy->Draw("PE SAME");

	TLegend* legend2 = new TLegend(.42, .9, .72, .66);
	legend2->SetTextSize(.05);
	
	legend2->AddEntry(hTotalHydjetGenRapidityCopy, "Hydjet denominator", "lep");
	legend2->AddEntry(hPassedHydjetGenRapidityCopy, "Hydjet numerator", "lep");

	legend2->AddEntry(hTotalNoFilterGenRapidityCopy, "No Filter denominator", "lep");
	legend2->AddEntry(hPassedNoFilterGenRapidityCopy, "No Filter numerator", "lep");

	legend2->Draw();

	gSystem->mkdir("plots", kTRUE);
	canvas->SaveAs("plots/MCRapidityAcceptance.png", "RECREATE");
}