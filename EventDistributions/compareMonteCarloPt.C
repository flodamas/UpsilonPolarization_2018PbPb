#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../MonteCarlo/AccEffHelpers.h"

#include "../Tools/Style/Legends.h"
#include "../Tools/Parameters/CentralityValues.h"

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

	TEfficiency* hAccPt = new TEfficiency("", ";p^{#varUpsilon}_{T} (GeV/c);Acceptance of #varUpsilon(1S)", NPtFineBins, gPtFineBinning);

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

	Int_t Centrality;

	Short_t Gen_QQ_whichRec[1000];

	Short_t Gen_QQ_mupl_idx[1000];
	Short_t Gen_QQ_mumi_idx[1000];

	Short_t Gen_QQ_size;
	Short_t Reco_QQ_sign[1000];

	// event variables
	OniaTree->SetBranchAddress("Gen_weight", &Gen_weight);
	OniaTree->SetBranchAddress("Centrality", &Centrality);

	// gen-level variables
	OniaTree->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size);
	OniaTree->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom);
	OniaTree->SetBranchAddress("Gen_QQ_whichRec", Gen_QQ_whichRec);

	OniaTree->SetBranchAddress("Gen_mu_4mom", &Gen_mu_4mom);
	OniaTree->SetBranchAddress("Gen_QQ_mupl_idx", Gen_QQ_mupl_idx);
	OniaTree->SetBranchAddress("Gen_QQ_mumi_idx", Gen_QQ_mumi_idx);

	OniaTree->SetBranchAddress("Reco_QQ_sign", &Reco_QQ_sign);

	TEfficiency* hAccPt = new TEfficiency("", ";p^{#varUpsilon}_{T} (GeV/c);Acceptance of #varUpsilon(1S)", NPtFineBins, gPtFineBinning);

	// loop variables
	TLorentzVector* gen_QQ_LV = new TLorentzVector();
	TLorentzVector* gen_mumi_LV = new TLorentzVector();
	TLorentzVector* gen_mupl_LV = new TLorentzVector();

	Bool_t withinAcceptance;

	Long64_t totEntries = OniaTree->GetEntries();

	double eventWeight;

	// Loop over the events
	for (Long64_t iEvent = 0; iEvent < (totEntries); iEvent++) {
		if (iEvent % 10000 == 0) {
			cout << Form("\rProcessing event %lld / %lld (%.0f%%)", iEvent, totEntries, 100. * iEvent / totEntries) << flush;
		}

		// if (iEvent > 50000) break;

		OniaTree->GetEntry(iEvent);

		if (Centrality >= 2 * gCentralityBinMax) continue;

		eventWeight = (Gen_weight * FindNcoll(Centrality));

		// loop over all gen upsilons
		for (int iGen = 0; iGen < Gen_QQ_size; iGen++) {
			gen_QQ_LV = (TLorentzVector*)Gen_QQ_4mom->At(iGen);

			if (fabs(gen_QQ_LV->Rapidity()) < gRapidityMin || fabs(gen_QQ_LV->Rapidity()) > gRapidityMax) continue; // upsilon within fiducial region

			// go to reco level
			Int_t iReco = Gen_QQ_whichRec[iGen];

			if (Reco_QQ_sign[iReco] != 0) continue; // only opposite-sign muon pairs

			// single-muon acceptance cuts
			gen_mupl_LV = (TLorentzVector*)Gen_mu_4mom->At(Gen_QQ_mupl_idx[iGen]);
			gen_mumi_LV = (TLorentzVector*)Gen_mu_4mom->At(Gen_QQ_mumi_idx[iGen]);

			withinAcceptance = (fabs(gen_mupl_LV->Eta()) < 2.4) && (gen_mupl_LV->Pt() > gMuonPtCut) && (fabs(gen_mumi_LV->Eta()) < 2.4) && (gen_mumi_LV->Pt() > gMuonPtCut);

			hAccPt->FillWeighted(withinAcceptance, eventWeight, gen_QQ_LV->Pt());
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

	TH1D* hPassedNoFilterGenPt = (TH1D*) hAccNoFilterGenPt->GetPassedHistogram();
	TH1D* hTotalNoFilterGenPt = (TH1D*) hAccNoFilterGenPt->GetTotalHistogram();

	TH1D* hPassedHydjetGenPt = (TH1D*) hAccHydjetGenPt->GetPassedHistogram();
	TH1D* hTotalHydjetGenPt = (TH1D*) hAccHydjetGenPt->GetTotalHistogram();

	TH1D* hPassedNoFilterGenPtCopy = (TH1D*) hPassedNoFilterGenPt->Clone();
	TH1D* hTotalNoFilterGenPtCopy = (TH1D*) hTotalNoFilterGenPt->Clone();

	TH1D* hPassedHydjetGenPtCopy = (TH1D*) hPassedHydjetGenPt->Clone();
	TH1D* hTotalHydjetGenPtCopy = (TH1D*) hTotalHydjetGenPt->Clone();

	auto canvas = new TCanvas("canvas", "", 1200, 700);

	canvas->Divide(2);
	canvas->cd(1);

	gPad->SetLeftMargin(0.21);
	gPad->SetRightMargin(0.03);
	gPad->SetTopMargin(0.07);

	TH1D* hAccDummy = new TH1D("hAccDummy", ";p^{#varUpsilon}_{T} (GeV/c);Acceptance of #varUpsilon(1S)", NPtFineBins, gPtFineBinning);

	hAccDummy->GetXaxis()->SetTitleSize(0.075);

	hAccDummy->GetYaxis()->SetTitleOffset(1.2);
	hAccDummy->GetYaxis()->SetTitleSize(0.075);
	hAccDummy->GetYaxis()->SetLabelSize(0.065);

	hAccDummy->GetYaxis()->SetRangeUser(0, 1);

	hAccNoFilterGenPt->SetLineColor(kRed + 1);
	hAccNoFilterGenPt->SetMarkerStyle(8);
	hAccNoFilterGenPt->SetMarkerColor(kRed + 1);
	hAccNoFilterGenPt->SetLineWidth(2);

	hAccHydjetGenPt->SetLineColor(kAzure + 1);
	hAccHydjetGenPt->SetMarkerStyle(8);
	hAccHydjetGenPt->SetMarkerColor(kAzure + 1);
	hAccHydjetGenPt->SetLineWidth(2);

	hAccDummy->Draw();

	hAccNoFilterGenPt->Draw("P SAME");

	hAccHydjetGenPt->Draw("P SAME");

	TLegend* legend = new TLegend(.22, .9, .52, .8);
	legend->SetTextSize(.05);
	legend->AddEntry(hAccHydjetGenPt, "Hydjet Gen", "lep");
	legend->AddEntry(hAccNoFilterGenPt, "No Filter Gen", "lep");

	legend->Draw();

	canvas->cd(2);

	gPad->SetLeftMargin(0.1);
	gPad->SetRightMargin(0.03);
	gPad->SetTopMargin(0.07);

	TPad* pad1 = new TPad("pad1", "pad1", 0., 0.3, 1, 1.0);
	pad1->SetBottomMargin(0.04);
	pad1->SetLeftMargin(0.2);
	pad1->SetTopMargin(0.1);
	pad1->Draw();
	pad1->cd();

	double NoFilterSum = hTotalNoFilterGenPt->Integral();
	double HydjetSum = hTotalHydjetGenPt->Integral();

	hPassedNoFilterGenPtCopy->Scale(1. / hPassedNoFilterGenPt->Integral(), "Width");
	hTotalNoFilterGenPtCopy->Scale(1. / hTotalNoFilterGenPt->Integral(), "Width");

	hPassedHydjetGenPtCopy->Scale(1. / hPassedHydjetGenPt->Integral(), "Width");
	hTotalHydjetGenPtCopy->Scale(1. / hTotalHydjetGenPt->Integral(), "Width");

	TH1D* hPassedRatio = (TH1D*) hPassedHydjetGenPtCopy->Clone();
	hPassedRatio->Divide(hPassedNoFilterGenPtCopy);

	TH1D* hDummy = new TH1D("hDummy", ";p^{#varUpsilon}_{T} (GeV/c);dN(#varUpsilon(1S)) / dp_{T} (GeV/c)^{-1}", NPtFineBins, gPtFineBinning);

	hDummy->GetXaxis()->SetTitleSize(0.077);
	hDummy->GetXaxis()->SetLabelSize(0.);

	hDummy->GetYaxis()->SetTitleOffset(1.2);
	hDummy->GetYaxis()->SetTitleSize(0.075);
	hDummy->GetYaxis()->SetLabelSize(0.065);

	hDummy->GetYaxis()->SetRangeUser(0, hTotalNoFilterGenPtCopy->GetMaximum() * 1.2);

	// hDummy->GetYaxis()->SetMaxDigits(2);

	hTotalNoFilterGenPtCopy->SetLineColor(kRed + 1);
	hTotalNoFilterGenPtCopy->SetMarkerStyle(8);
	hTotalNoFilterGenPtCopy->SetMarkerColor(kRed + 1);
	hTotalNoFilterGenPtCopy->SetLineWidth(2);

	hPassedNoFilterGenPtCopy->SetLineColor(kPink - 4);
	hPassedNoFilterGenPtCopy->SetMarkerStyle(8);
	hPassedNoFilterGenPtCopy->SetMarkerColor(kPink - 4);
	hPassedNoFilterGenPtCopy->SetLineWidth(2);

	hTotalHydjetGenPtCopy->SetLineColor(kAzure + 1);
	hTotalHydjetGenPtCopy->SetMarkerStyle(8);
	hTotalHydjetGenPtCopy->SetMarkerColor(kAzure + 1);
	hTotalHydjetGenPtCopy->SetLineWidth(2);

	hPassedHydjetGenPtCopy->SetLineColor(kCyan - 3);
	hPassedHydjetGenPtCopy->SetMarkerStyle(8);
	hPassedHydjetGenPtCopy->SetMarkerColor(kCyan - 3);
	hPassedHydjetGenPtCopy->SetLineWidth(2);

	hDummy->Draw();

	hTotalNoFilterGenPtCopy->Draw("PE SAME");

	hPassedNoFilterGenPtCopy->Draw("PE SAME");

	hTotalHydjetGenPtCopy->Draw("PE SAME");

	hPassedHydjetGenPtCopy->Draw("PE SAME");

	TLegend* legend2 = new TLegend(.5, .85, .8, .55);
	legend2->SetTextSize(.05);
	
	legend2->AddEntry(hTotalHydjetGenPtCopy, "Hydjet denominator", "lep");
	legend2->AddEntry(hPassedHydjetGenPtCopy, "Hydjet numerator", "lep");
	
	legend2->AddEntry(hTotalNoFilterGenPtCopy, "No Filter denominator", "lep");
	legend2->AddEntry(hPassedNoFilterGenPtCopy, "No Filter numerator", "lep");

	legend2->Draw();

	// ratio plot
	canvas->cd(2);

	TPad* pad2 = new TPad("bottomPad", "bottomPad", 0, 0.0, 1, .31);
	pad2->SetTopMargin(0.025);
	pad2->SetBottomMargin(0.3);
	pad2->SetLeftMargin(0.2);
	pad2->SetTicks(1, 1);
	pad2->Draw();
	pad2->cd();

	hPassedRatio->SetTitle(" ");
	hPassedRatio->GetYaxis()->SetTitleOffset(0.5);
	hPassedRatio->GetYaxis()->SetTitle("#splitline{Hydjet /}{NoFilter}");
	hPassedRatio->GetYaxis()->SetTitleSize(0.15);
	hPassedRatio->GetYaxis()->SetLabelSize(0.15);
	hPassedRatio->GetYaxis()->CenterTitle();

	hPassedRatio->GetXaxis()->SetTitle("p^{#varUpsilon}_{T} (GeV/c)");
	hPassedRatio->GetXaxis()->SetLabelSize(0.15);
	hPassedRatio->GetXaxis()->SetTitleSize(0.15);
	hPassedRatio->GetXaxis()->SetTickSize(0.06);

	hPassedRatio->GetYaxis()->SetNdivisions(505);
	
	hPassedRatio->Draw("PZ");

	hPassedRatio->SetMinimum(0.);
	hPassedRatio->SetMaximum(3.);

	gSystem->mkdir("plots", kTRUE);
	canvas->SaveAs("plots/MCPtAcceptance.png", "RECREATE");
}