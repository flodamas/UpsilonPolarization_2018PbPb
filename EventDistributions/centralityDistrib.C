#include "../AnalysisParameters.h"

#include "../Tools/Parameters/CentralityValues.h"

#include "../Tools/Parameters/PhaseSpace.h"

// plot the event centrality distribution for triggered data, reco MC with and w/o Ncoll weighting to show the importance of reweighting the Hydjet MC sample

TH1D* CentDistribution(const char* inputFileName = "OniaTree_Y1S_pThat2_HydjetDrumMB_miniAOD.root", bool isMC = false, bool doNcollWeighting = false) {
	TFile* infile = TFile::Open(Form("../Files/%s", inputFileName), "READ");
	TTree* OniaTree = (TTree*)infile->Get("hionia/myTree");

	/// OniaTree variables
	Int_t Centrality;
	ULong64_t HLTriggers;
	Float_t zVtx;

	Short_t Gen_QQ_size;
	Float_t Gen_weight;

	Short_t Reco_QQ_size;
	TClonesArray* CloneArr_QQ = nullptr;
	Float_t Reco_QQ_VtxProb[1000];
	Short_t Reco_QQ_mupl_idx[1000];
	Short_t Reco_QQ_mumi_idx[1000];
	ULong64_t Reco_QQ_trig[1000];

	TClonesArray* CloneArr_mu = nullptr;
	Int_t Reco_mu_SelectionType[1000];
	Int_t Reco_mu_nPixWMea[1000];
	Int_t Reco_mu_nTrkWMea[1000];
	Float_t Reco_mu_dxy[1000];
	Float_t Reco_mu_dz[1000];

	OniaTree->SetBranchAddress("Centrality", &Centrality);
	OniaTree->SetBranchAddress("HLTriggers", &HLTriggers);
	OniaTree->SetBranchAddress("zVtx", &zVtx);

	if (isMC) {
		OniaTree->SetBranchAddress("Gen_weight", &Gen_weight);
		OniaTree->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size);
	}

	OniaTree->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size);
	OniaTree->SetBranchAddress("Reco_QQ_4mom", &CloneArr_QQ);
	OniaTree->SetBranchAddress("Reco_QQ_VtxProb", &Reco_QQ_VtxProb);
	OniaTree->SetBranchAddress("Reco_QQ_mupl_idx", Reco_QQ_mupl_idx);
	OniaTree->SetBranchAddress("Reco_QQ_mumi_idx", Reco_QQ_mumi_idx);
	OniaTree->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig);

	OniaTree->SetBranchAddress("Reco_mu_4mom", &CloneArr_mu);

	OniaTree->SetBranchAddress("Reco_mu_SelectionType", Reco_mu_SelectionType);
	OniaTree->SetBranchAddress("Reco_mu_nPixWMea", &Reco_mu_nPixWMea);
	OniaTree->SetBranchAddress("Reco_mu_nTrkWMea", &Reco_mu_nTrkWMea);
	OniaTree->SetBranchAddress("Reco_mu_dxy", &Reco_mu_dxy);
	OniaTree->SetBranchAddress("Reco_mu_dz", &Reco_mu_dz);

	auto histo = new TH1D("histo", ";centrality (in %);Self-normalized distribution", 180, 0, 90);
	histo->SetDirectory(0);

	Long64_t totEntries = OniaTree->GetEntries();

	bool skipEvent = false;

	for (Long64_t iEvent = 0; iEvent < (totEntries); iEvent++) {
		if (iEvent % 10000 == 0) {
			cout << Form("\rProcessing event %lld / %lld (%.0f%%)", iEvent, totEntries, 100. * iEvent / totEntries) << flush;
		}

		OniaTree->GetEntry(iEvent);

		if (!((HLTriggers & (ULong64_t)(1 << (gUpsilonHLTBit - 1))) == (ULong64_t)(1 << (gUpsilonHLTBit - 1)))) continue; // must fire the upsilon HLT path

		//		if (fabs(zVtx) > 20) continue;

		float weight = 1;

		if (isMC) {
			if (doNcollWeighting) weight *= FindNcoll(Centrality);

			for (int iGen = 0; iGen < Gen_QQ_size; iGen++) {
				weight *= Gen_weight; // one gen upsilon per event in principle
			}
		}

		// test: only counting events containing at least one potential upsilon candidate
		skipEvent = true;

		// loop over reconstructed dimuon candidates
		for (int iQQ = 0; iQQ < Reco_QQ_size; iQQ++) {
			// dimuon quality conditions

			if (!((Reco_QQ_trig[iQQ] & (ULong64_t)(1 << (gUpsilonHLTBit - 1))) == (ULong64_t)(1 << (gUpsilonHLTBit - 1)))) continue;

			if (Reco_QQ_VtxProb[iQQ] < 0.05) continue;

			TLorentzVector* Reco_QQ_4mom = (TLorentzVector*)CloneArr_QQ->At(iQQ);

			if (fabs(Reco_QQ_4mom->Eta()) > 2.4) continue;

			if (Reco_QQ_4mom->M() < 9 || Reco_QQ_4mom->M() > 10) continue;

			// muon quality criteria

			int iMuPlus = Reco_QQ_mupl_idx[iQQ];
			int iMuMinus = Reco_QQ_mumi_idx[iQQ];

			if (!(((Reco_mu_SelectionType[iMuPlus] & 2) && (Reco_mu_SelectionType[iMuPlus] & 8)) && ((Reco_mu_SelectionType[iMuMinus] & 2) && (Reco_mu_SelectionType[iMuMinus] & 8)))) continue; // tracker and global muons

			if (!(((Reco_mu_nTrkWMea[iMuPlus] > 5) && (Reco_mu_nPixWMea[iMuPlus] > 0) && (fabs(Reco_mu_dxy[iMuPlus]) < 0.3) && (fabs(Reco_mu_dz[iMuPlus]) < 20.)) && ((Reco_mu_nTrkWMea[iMuPlus] > 5) && (Reco_mu_nPixWMea[iMuPlus] > 0) && (fabs(Reco_mu_dxy[iMuPlus]) < 0.3) && (fabs(Reco_mu_dz[iMuPlus]) < 20.)))) continue; // hybrid-soft

			TLorentzVector* Reco_mupl_4mom = (TLorentzVector*)CloneArr_mu->At(iMuPlus);

			if (!MuonUpsilonTriggerAcc(*Reco_mupl_4mom)) continue;

			TLorentzVector* Reco_mumi_4mom = (TLorentzVector*)CloneArr_mu->At(iMuMinus);

			if (!MuonUpsilonTriggerAcc(*Reco_mumi_4mom)) continue;

			// if we arrive, it means that the reco'ed dimuon is a good upsilon candidate
			skipEvent = false;
			break; // no need to verify on the other dimuons
		}

		if (skipEvent) continue;

		histo->Fill(Centrality / 2., weight);
	}

	infile->Close();

	histo->Scale(1. / histo->Integral());

	return histo;
}

void centralityDistrib() {
	writeExtraText = true; // if extra text
	extraText = "       Internal";

	gStyle->SetPadRightMargin(.03);

	auto canvas = new TCanvas("canvas", "", 600, 650);
	TPad* pad1 = new TPad("pad1", "pad1", 0, 0.4, 1, 1.0);
	pad1->SetBottomMargin(0.04);
	pad1->SetTopMargin(0.1);
	pad1->Draw();
	pad1->cd();

	auto distribData = CentDistribution("OniaTree_miniAOD_PbPbPrompt_112X_DATA_ep.root");

	distribData->SetMarkerColor(kAzure + 1);
	distribData->SetLineWidth(2);
	distribData->SetLineColor(kAzure + 1);

	distribData->GetXaxis()->SetLabelOffset(1);

	distribData->Draw("PZ");

	auto distribMCnoNcoll = CentDistribution("OniaTree_Y1S_pThat2_HydjetDrumMB_miniAOD.root", true, false);

	distribMCnoNcoll->SetLineColor(kRed + 1);
	distribMCnoNcoll->SetMarkerColor(kRed + 1);
	distribMCnoNcoll->SetLineWidth(2);

	distribMCnoNcoll->Draw("SAME PZ");

	auto distribMCwithNcoll = CentDistribution("OniaTree_Y1S_pThat2_HydjetDrumMB_miniAOD.root", true, true);

	distribMCwithNcoll->SetLineColor(kOrange + 1);
	distribMCwithNcoll->SetMarkerColor(kOrange + 1);
	distribMCwithNcoll->SetLineWidth(2);

	distribMCwithNcoll->Draw("SAME PZ");

	// few cosmetics

	TLegend* legend = new TLegend(.35, .8, .65, .5);
	legend->SetTextSize(.06);
	legend->AddEntry(distribData, "#varUpsilon(1S) data candidates", "lep");
	legend->AddEntry((TObject*)0, "Hydjet-embedded MC events", "");
	legend->AddEntry(distribMCnoNcoll, "w/o N_{coll} weighting", "lep");
	legend->AddEntry(distribMCwithNcoll, "with N_{coll} weighting", "lep");

	legend->Draw();

	// ratio plot
	canvas->cd();

	TPad* pad2 = new TPad("bottomPad", "bottomPad", 0, 0.0, 1, .4);
	pad2->SetTopMargin(0.03);
	pad2->SetBottomMargin(0.25);
	pad2->SetTicks(1, 1);
	pad2->Draw();
	pad2->cd();

	TH1D* ratioPlot = (TH1D*)distribData->Clone("ratioPlot");

	ratioPlot->Divide(distribMCwithNcoll);

	ratioPlot->SetTitle(" ");
	ratioPlot->GetYaxis()->SetTitleOffset(0.63);
	ratioPlot->GetYaxis()->SetTitle("Data / N_{coll} weighted MC");
	ratioPlot->GetYaxis()->SetTitleSize(0.11);
	ratioPlot->GetYaxis()->SetLabelSize(0.09);
	ratioPlot->GetYaxis()->CenterTitle();

	ratioPlot->GetXaxis()->SetLabelSize(0.09);
	ratioPlot->GetXaxis()->SetTitleSize(0.11);
	ratioPlot->GetXaxis()->SetTickSize(0.05);
	ratioPlot->GetXaxis()->SetLabelOffset(0.01);

	//ratioPlot->GetYaxis()->SetNdivisions(505);

	ratioPlot->SetMarkerColor(kMagenta + 2);
	ratioPlot->SetLineColor(kMagenta + 2);

	ratioPlot->Draw("SAME PZ");

	ratioPlot->SetMinimum(0.);
	ratioPlot->SetMaximum(1.92);

	canvas->Modified();
	canvas->Update();
	canvas->cd();

	pad1->Draw();

	//pad2->Draw();

	CMS_lumi(canvas, gCMSLumiText);

	canvas->SaveAs("plots/centralityDistrib.pdf", "RECREATE");
}
