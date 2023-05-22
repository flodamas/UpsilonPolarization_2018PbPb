#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TRotation.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"

#include <iostream>
#include <fstream>
#include <cmath>

#include "../Tools/Style/tdrStyle.C"
#include "../Tools/Style/CMS_lumi.C"

#include "../Tools/Style/FitDistributions.h"
#include "../Tools/Style/Legends.h"

#include "../Tools/Parameters/AnalysisParameters.h"

#include "../Tools/Parameters/CentralityValues.h"

#include "../ReferenceFrameTransform/Transformations.h"

void recoDistribution(Int_t ptMin = 0, Int_t ptMax = 30) {
	const char* filename = "../Files/OniaTree_Y1S_miniAOD_HydjetEmbeddedMC.root";
	TFile* file = TFile::Open(filename, "READ");
	if (!file) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
		return;
	}

	cout << "File " << filename << " opened" << endl;

	writeExtraText = true; // if extra text
	extraText = "       Internal";

	Int_t centMin = 0, centMax = 90;

	TTree* OniaTree = (TTree*)file->Get("hionia/myTree");

	/// OniaTree variables

	Float_t SumET_HF;

	TClonesArray* Reco_QQ_4mom = nullptr;

	ULong64_t HLTriggers;
	ULong64_t Reco_QQ_trig[1000];
	Int_t Centrality;
	TClonesArray* CloneArr_QQ = nullptr;
	TClonesArray* CloneArr_mu = nullptr;
	Short_t Reco_QQ_size;
	Short_t Reco_QQ_sign[1000];
	Short_t Reco_QQ_mupl_idx[1000];
	Short_t Reco_QQ_mumi_idx[1000];
	Short_t Reco_QQ_whichGen[1000];

	ULong64_t Reco_mu_trig[1000];
	Int_t Reco_mu_SelectionType[1000];
	//(parameters for quality cuts)
	Float_t Reco_QQ_VtxProb[1000];
	Int_t Reco_mu_nPixWMea[1000];
	Int_t Reco_mu_nTrkWMea[1000];
	Float_t Reco_mu_dxy[1000];
	Float_t Reco_mu_dz[1000];

	// event variables
	OniaTree->SetBranchAddress("SumET_HF", &SumET_HF);
	OniaTree->SetBranchAddress("HLTriggers", &HLTriggers);
	OniaTree->SetBranchAddress("Centrality", &Centrality);

	OniaTree->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig);
	OniaTree->SetBranchAddress("Reco_QQ_whichGen", Reco_QQ_whichGen);
	OniaTree->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom);
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

	// (cos theta, phi) 2D distribution maps for CS and HX frames

	Int_t nCosThetaBins = 20;
	Float_t cosThetaMin = -1, cosThetaMax = 1;

	//	Double_t cosThetaBinning[] = {-1, -0.8, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1};
	//	nCosThetaBins = sizeof(cosThetaBinning) / sizeof(Double_t) - 1;

	Int_t nPhiBins = 21;
	Float_t phiMin = -180, phiMax = 240;

	//Double_t phiBinning[] = {-TMath::Pi(), -2.5, -2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2, 2.5, 180};
	//	nPhiBins = sizeof(cosThetaBinning) / sizeof(Double_t) - 1;

	TH2F* hCS = new TH2F("hCS", ";cos #theta_{CS}; #varphi_{CS} (#circ);# reco QQ", nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);
	TH2F* hHX = new TH2F("hHX", ";cos #theta_{HX}; #varphi_{HX} (#circ);# reco QQ", nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);

	TLorentzVector* recoLV = new TLorentzVector();

	double weight, muTrkSF, muIdSF, muTrigSF;
	Bool_t allGood, firesTrigger, isRecoMatched, goodVertexProba, withinAcceptance, trackerAndGlobalMuons, hybridSoftMuons;

	Long64_t totEntries = OniaTree->GetEntries();

	for (Long64_t iEvent = 0; iEvent < (totEntries); iEvent++) {
		if (iEvent % 10000 == 0) {
			cout << Form("\rProcessing event %lld / %lld (%.0f%%)", iEvent, totEntries, 100. * iEvent / totEntries) << flush;
		}

		OniaTree->GetEntry(iEvent);
		//recoLV->Clear();

		// event selection

		if (Centrality >= 2 * 90) continue; // discard events with centrality >= 90% in 2018 data

		firesTrigger = ((HLTriggers & (ULong64_t)(1 << (UpsilonHLTBit - 1))) == (ULong64_t)(1 << (UpsilonHLTBit - 1)));

		if (!firesTrigger) continue;

		// get N_coll from HF
		Int_t hiBin = GetHiBinFromhiHF(SumET_HF);

		weight = FindNcoll(hiBin);

		// loop over all reco QQ candidates
		for (int iReco = 0; iReco < Reco_QQ_size; iReco++) {
			recoLV = (TLorentzVector*)Reco_QQ_4mom->At(iReco);

			if (recoLV->Pt() < ptMin || recoLV->Pt() > ptMax) continue; // pt bin of interest

			if (fabs(recoLV->Rapidity()) > 2.4) continue; // upsilon within acceptance

			if (Reco_QQ_whichGen[iReco] < 0) continue; //matching gen

			if (Reco_QQ_sign[iReco] != 0) continue; // only opposite-sign muon pairs

			if (Reco_QQ_VtxProb[iReco] < 0.01) continue;

			/// single-muon selection criteria
			int iMuPlus = Reco_QQ_mupl_idx[iReco];
			int iMuMinus = Reco_QQ_mumi_idx[iReco];

			// global AND tracker muons
			trackerAndGlobalMuons = (Reco_mu_SelectionType[iMuPlus] & 2) && (Reco_mu_SelectionType[iMuPlus] & 8) && (Reco_mu_SelectionType[iMuMinus] & 2) && (Reco_mu_SelectionType[iMuMinus] & 8);

			if (!trackerAndGlobalMuons) continue;

			// passing hybrid-soft Id
			hybridSoftMuons = (Reco_mu_nTrkWMea[iMuPlus] > 5) && (Reco_mu_nPixWMea[iMuPlus] > 0) && (fabs(Reco_mu_dxy[iMuPlus]) < 0.3) && (fabs(Reco_mu_dz[iMuPlus]) < 20.) && (Reco_mu_nTrkWMea[iMuMinus] > 5) && (Reco_mu_nPixWMea[iMuMinus] > 0) && (fabs(Reco_mu_dxy[iMuMinus]) < 0.3) && (fabs(Reco_mu_dz[iMuMinus]) < 20.);

			if (!hybridSoftMuons) continue;

			// acceptance
			TLorentzVector* Reco_mupl_4mom = (TLorentzVector*)CloneArr_mu->At(iMuPlus);

			TLorentzVector* Reco_mumi_4mom = (TLorentzVector*)CloneArr_mu->At(iMuMinus);

			withinAcceptance = (1.2 < fabs(Reco_mupl_4mom->Eta()) < 2.4) && (Reco_mupl_4mom->Pt() > 3.) && (1.2 < fabs(Reco_mumi_4mom->Eta()) < 2.4) && (Reco_mumi_4mom->Pt() > 3.);

			if (!withinAcceptance) continue;

			// get positive muon's coordinates in the studied reference frames
			TVector3 muPlus_CS = MuPlusVector_CollinsSoper(*recoLV, *Reco_mupl_4mom);

			TVector3 muPlus_HX = MuPlusVector_Helicity(*recoLV, *Reco_mupl_4mom);

			hCS->Fill(muPlus_CS.CosTheta(), muPlus_CS.Phi() * 180 / TMath::Pi(), weight);
			hHX->Fill(muPlus_HX.CosTheta(), muPlus_HX.Phi() * 180 / TMath::Pi(), weight);
		}
	}

	gStyle->SetPadLeftMargin(.15);
	//gStyle->SetTitleYOffset(.9);
	gStyle->SetPadRightMargin(0.18);
	gStyle->SetPalette(kRainBow);

	TLatex* legend = new TLatex();
	legend->SetTextAlign(22);
	legend->SetTextSize(0.05);

	auto* canvasCS = new TCanvas("canvasCS", "", 700, 600);
	hCS->Draw("COLZ");

	CMS_lumi(canvasCS, "#varUpsilon(1S) Hydjet-embedded PbPb MC");

	legend->DrawLatexNDC(.5, .88, Form("centrality %d-%d%%, %d < p_{T}^{#mu#mu} < %d GeV", centMin, centMax, ptMin, ptMax));

	gPad->Update();

	hCS->GetYaxis()->SetRangeUser(-190, 240);
	//	hCS->GetPaintedHistogram()->GetZaxis()->SetRangeUser(0, 0.9);

	canvasCS->SaveAs(Form("RecoMaps/CS_cent%dto%d_pt%dto%dGeV.png", centMin, centMax, ptMin, ptMax), "RECREATE");

	auto* canvasHX = new TCanvas("canvasHX", "", 700, 600);
	hHX->Draw("COLZ");

	CMS_lumi(canvasHX, "#varUpsilon(1S) Hydjet-embedded PbPb MC");

	legend->DrawLatexNDC(.5, .88, Form("centrality %d-%d%%, %d < p_{T}^{#mu#mu} < %d GeV", centMin, centMax, ptMin, ptMax));

	gPad->Update();

	hHX->GetYaxis()->SetRangeUser(-190, 240);
	//	hHX->GetPaintedHistogram()->GetZaxis()->SetRangeUser(0, 0.9);

	canvasHX->SaveAs(Form("RecoMaps/HX_cent%dto%d_pt%dto%dGeV.png", centMin, centMax, ptMin, ptMax), "RECREATE");
}
