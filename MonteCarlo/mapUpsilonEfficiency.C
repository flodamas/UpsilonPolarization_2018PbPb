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

#include "../AnalysisParameters.h"

#include "../Tools/Parameters/CentralityValues.h"

#include "../Tools/Parameters/MuonScaleFactors.h"

#include "../ReferenceFrameTransform/Transformations.h"

// helping function to compute the product of the muon trigger efficiency scale factors when the two muons pass the L3 filter
double DimuonL3TriggerWeight(float pt_mupl, float eta_mupl, float pt_mumi, float eta_mumi, int indexSF = 0) {
	double T1_ = tnp_weight_trg_pbpb_mc(pt_mupl, eta_mupl, 3, indexSF);
	double T2_ = tnp_weight_trg_pbpb_mc(pt_mumi, eta_mumi, 3, indexSF);
	double T1 = tnp_weight_trg_pbpb_mc(pt_mupl, eta_mupl, 2, indexSF);
	double T2 = tnp_weight_trg_pbpb_mc(pt_mumi, eta_mumi, 2, indexSF);
	double den_ = T1_ * T2 + (T1 - T1_) * T2_;
	double num_ = T1_ * tnp_weight_trg_pbpb(pt_mupl, eta_mupl, 3, indexSF) * T2 * tnp_weight_trg_pbpb(pt_mumi, eta_mumi, 2, indexSF) + (T1 * tnp_weight_trg_pbpb(pt_mupl, eta_mupl, 2, indexSF) - T1_ * tnp_weight_trg_pbpb(pt_mupl, eta_mupl, 3, indexSF)) * T2_ * tnp_weight_trg_pbpb(pt_mumi, eta_mumi, 3, indexSF);

	if (den_ <= 0 || num_ <= 0) {
		cout << "ERROR wrong calculation" << endl;
		return 0;
	}

	return num_ / den_;
}

// return the 2D map of the relative systematic efficiency uncertainty
TH2D* RelSystEffHist(TEfficiency* hNominal, TEfficiency* hTrk_systUp, TEfficiency* hTrk_systDown, TEfficiency* hMuId_systUp, TEfficiency* hMuId_systDown, TEfficiency* hTrig_systUp, TEfficiency* hTrig_systDown, TEfficiency* hTrk_statUp, TEfficiency* hTrk_statDown, TEfficiency* hMuId_statUp, TEfficiency* hMuId_statDown, TEfficiency* hTrig_statUp, TEfficiency* hTrig_statDown) {
	// clone one of the input efficiency map, to make sure we get the correct binning
	TH2D* hTotalSyst = (TH2D*)hNominal->GetCopyPassedHisto(); // will rename it outside this function

	// useful loop variables
	int globalBin;
	double nominalEff;
	double trk_systEff, muId_systEff, trig_systEff, systEffSquared;
	double trk_statEff, muId_statEff, trig_statEff, statEffSquared;
	double finalSystEff;

	for (int iCosThetaBin = 1; iCosThetaBin <= hTotalSyst->GetNbinsX(); iCosThetaBin++) {
		for (int iPhiBin = 1; iPhiBin <= hTotalSyst->GetNbinsY(); iPhiBin++) {
			globalBin = hNominal->GetGlobalBin(iCosThetaBin, iPhiBin); // to run though the efficiency maps

			nominalEff = hNominal->GetEfficiency(globalBin);

			// compute the systematic uncertainties associated to the variation of the muon SF uncertainties

			// instructions can be found here: https://twiki.cern.ch/twiki/pub/CMS/HIMuonTagProbe/TnpHeaderFile.pdf#page=5

			// from muon SF systematics
			trk_systEff = max(abs(hTrk_systUp->GetEfficiency(globalBin) - nominalEff), abs(hTrk_systDown->GetEfficiency(globalBin) - nominalEff));

			muId_systEff = max(abs(hMuId_systUp->GetEfficiency(globalBin) - nominalEff), abs(hMuId_systDown->GetEfficiency(globalBin) - nominalEff));

			trig_systEff = max(abs(hTrig_systUp->GetEfficiency(globalBin) - nominalEff), abs(hTrig_systDown->GetEfficiency(globalBin) - nominalEff));

			systEffSquared = trk_systEff * trk_systEff + muId_systEff * muId_systEff + trig_systEff * trig_systEff;

			// from muon SF statistical uncertainties
			trk_statEff = max(abs(hTrk_statUp->GetEfficiency(globalBin) - nominalEff), abs(hTrk_statDown->GetEfficiency(globalBin) - nominalEff));

			muId_statEff = max(abs(hMuId_statUp->GetEfficiency(globalBin) - nominalEff), abs(hMuId_statDown->GetEfficiency(globalBin) - nominalEff));

			trig_statEff = max(abs(hTrig_statUp->GetEfficiency(globalBin) - nominalEff), abs(hTrig_statDown->GetEfficiency(globalBin) - nominalEff));

			statEffSquared = trk_statEff * trk_statEff + muId_statEff * muId_statEff + trig_statEff * trig_statEff;

			finalSystEff = sqrt(systEffSquared + statEffSquared);

			// store the RELATIVE SYSTEMATIC UNCERTAINTIES (not in % though) with respect to the nominal efficiency, more useful at the end
			hTotalSyst->SetBinContent(iCosThetaBin, iPhiBin, (nominalEff == 0) ? 0 : finalSystEff / nominalEff);

			//cout << "Efficiency in global bin " << globalBin << " = " << nominalEff << " +/- " << finalSystEff << endl;
		}
	}

	return hTotalSyst;
}

void mapUpsilonEfficiency(Int_t iState = 1, Int_t ptMin = 0, Int_t ptMax = 30) {
	const char* filename = Form("../Files/OniaTree_Y%dS_pThat2_HydjetDrumMB_miniAOD.root", iState);
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

	ULong64_t HLTriggers;
	ULong64_t Reco_QQ_trig[1000];
	Int_t Centrality;
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
	OniaTree->SetBranchAddress("HLTriggers", &HLTriggers);

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

	/// (cos theta, phi) 2D distribution maps for CS and HX frames

	Int_t nCosThetaBins = 20;
	Float_t cosThetaMin = -1, cosThetaMax = 1;

	//	Double_t cosThetaBinning[] = {-1, -0.8, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1};
	//	nCosThetaBins = sizeof(cosThetaBinning) / sizeof(Double_t) - 1;

	Int_t nPhiBins = 21;
	Float_t phiMin = -180, phiMax = 240;

	//Double_t phiBinning[] = {-TMath::Pi(), -2.5, -2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2, 2.5, 180};
	//	nPhiBins = sizeof(cosThetaBinning) / sizeof(Double_t) - 1;

	const char* titleCS = Form(";cos #theta_{CS}; #varphi_{CS} (#circ);#varUpsilon(%dS) total efficiency", iState);

	const char* titleHX = Form(";cos #theta_{HX}; #varphi_{HX} (#circ);#varUpsilon(%dS) total efficiency", iState);

	// we want to estimate the uncertainties from scale factors at the same time
	// instructions can be found here: https://twiki.cern.ch/twiki/pub/CMS/HIMuonTagProbe/TnpHeaderFile.pdf#page=5

	// for a given type of muon scale factor (i.e., tracking, muId, trigger) we need to compute the efficiency with up and down variations, for stat and syst uncertainties

	int indexNominal = 0;
	int indexSystUp = -1, indexSystDown = -2;
	int indexStatUp = +1, indexStatDown = +2;

	double dimuWeight_nominal = 0;
	double dimuWeight_trk_systUp, dimuWeight_trk_systDown, dimuWeight_trk_statUp, dimuWeight_trk_statDown;
	double dimuWeight_muId_systUp, dimuWeight_muId_systDown, dimuWeight_muId_statUp, dimuWeight_muId_statDown;
	double dimuWeight_trig_systUp, dimuWeight_trig_systDown, dimuWeight_trig_statUp, dimuWeight_trig_statDown;

	// Collins-Soper
	TEfficiency* hCS_nominal = new TEfficiency(Form("NominalEff_CS_pt%dto%d", ptMin, ptMax), titleCS, nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);

	TEfficiency* hCS_trk_systUp = new TEfficiency("hCS_trk_systUp", titleCS, nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);
	TEfficiency* hCS_trk_systDown = new TEfficiency("hCS_trk_systDown", titleCS, nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);
	TEfficiency* hCS_trk_statUp = new TEfficiency("hCS_trk_statUp", titleCS, nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);
	TEfficiency* hCS_trk_statDown = new TEfficiency("hCS_trk_statDown", titleCS, nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);

	TEfficiency* hCS_muId_systUp = new TEfficiency("hCS_muId_systUp", titleCS, nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);
	TEfficiency* hCS_muId_systDown = new TEfficiency("hCS_muId_systDown", titleCS, nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);
	TEfficiency* hCS_muId_statUp = new TEfficiency("hCS_muId_statUp", titleCS, nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);
	TEfficiency* hCS_muId_statDown = new TEfficiency("hCS_muId_statDown", titleCS, nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);

	TEfficiency* hCS_trig_systUp = new TEfficiency("hCS_trig_systUp", titleCS, nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);
	TEfficiency* hCS_trig_systDown = new TEfficiency("hCS_trig_systDown", titleCS, nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);
	TEfficiency* hCS_trig_statUp = new TEfficiency("hCS_trig_statUp", titleCS, nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);
	TEfficiency* hCS_trig_statDown = new TEfficiency("hCS_trig_statDown", titleCS, nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);

	// Helicity
	TEfficiency* hHX_nominal = new TEfficiency(Form("NominalEff_HX_pt%dto%d", ptMin, ptMax), titleHX, nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);

	TEfficiency* hHX_trk_systUp = new TEfficiency("hHX_trk_systUp", titleCS, nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);
	TEfficiency* hHX_trk_systDown = new TEfficiency("hHX_trk_systDown", titleCS, nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);
	TEfficiency* hHX_trk_statUp = new TEfficiency("hHX_trk_statUp", titleCS, nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);
	TEfficiency* hHX_trk_statDown = new TEfficiency("hHX_trk_statDown", titleCS, nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);

	TEfficiency* hHX_muId_systUp = new TEfficiency("hHX_muId_systUp", titleCS, nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);
	TEfficiency* hHX_muId_systDown = new TEfficiency("hHX_muId_systDown", titleCS, nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);
	TEfficiency* hHX_muId_statUp = new TEfficiency("hHX_muId_statUp", titleCS, nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);
	TEfficiency* hHX_muId_statDown = new TEfficiency("hHX_muId_statDown", titleCS, nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);

	TEfficiency* hHX_trig_systUp = new TEfficiency("hHX_trig_systUp", titleCS, nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);
	TEfficiency* hHX_trig_systDown = new TEfficiency("hHX_trig_systDown", titleCS, nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);
	TEfficiency* hHX_trig_statUp = new TEfficiency("hHX_trig_statUp", titleCS, nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);
	TEfficiency* hHX_trig_statDown = new TEfficiency("hHX_trig_statDown", titleCS, nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);

	// loop variables
	TLorentzVector* genLorentzVector = new TLorentzVector();

	double eventWeight, totalWeight;
	double dimuTrigWeight_nominal = -1, dimuTrigWeight_systUp = -1, dimuTrigWeight_systDown = -1, dimuTrigWeight_statUp = -1, dimuTrigWeight_statDown = -1;

	Bool_t allGood, firesTrigger, isRecoMatched, dimuonMatching, goodVertexProba, passHLTFilterMuons, trackerAndGlobalMuons, hybridSoftMuons;

	Long64_t totEntries = OniaTree->GetEntries();

	for (Long64_t iEvent = 0; iEvent < totEntries; iEvent++) {
		if (iEvent % 10000 == 0) {
			cout << Form("\rProcessing event %lld / %lld (%.0f%%)", iEvent, totEntries, 100. * iEvent / totEntries) << flush;
		}

		OniaTree->GetEntry(iEvent);
		//genLorentzVector->Clear();

		// event selection

		if (Centrality >= 2 * gCentralityBinMax) continue; // discard events with centrality >= 90% in 2018 data

		firesTrigger = ((HLTriggers & (ULong64_t)(1 << (gUpsilonHLTBit - 1))) == (ULong64_t)(1 << (gUpsilonHLTBit - 1)));

		eventWeight = Gen_weight * FindNcoll(Centrality);

		// loop over all gen upsilons
		for (int iGen = 0; iGen < Gen_QQ_size; iGen++) {
			genLorentzVector = (TLorentzVector*)Gen_QQ_4mom->At(iGen);

			// fiducial region
			if (genLorentzVector->Pt() < ptMin || genLorentzVector->Pt() > ptMax) continue; // pt bin of interest

			if (fabs(genLorentzVector->Rapidity()) < gRapidityMin || fabs(genLorentzVector->Rapidity()) > gRapidityMax) continue;

			// single-muon acceptance

			// positive muon first
			genLorentzVector = (TLorentzVector*)Gen_mu_4mom->At(Gen_QQ_mupl_idx[iGen]);

			if (genLorentzVector->Pt() < 3.5) continue;

			if (fabs(genLorentzVector->Eta()) > 2.4) continue;

			// then negative muon
			genLorentzVector = (TLorentzVector*)Gen_mu_4mom->At(Gen_QQ_mumi_idx[iGen]);

			if (genLorentzVector->Pt() < 3.5) continue;

			if (fabs(genLorentzVector->Eta()) > 2.4) continue;

			// go to reco level
			Int_t iReco = Gen_QQ_whichRec[iGen];

			if (Reco_QQ_sign[iReco] != 0) continue; // only opposite-sign muon pairs

			/// all the reconstructed upsilons must pass the conditions below!

			isRecoMatched = iReco > -1;

			if (isRecoMatched) {
				dimuonMatching = (Reco_QQ_trig[iReco] & (ULong64_t)(1 << (gUpsilonHLTBit - 1))) == (ULong64_t)(1 << (gUpsilonHLTBit - 1));

				goodVertexProba = Reco_QQ_VtxProb[iReco] > 0.01;

				/// single-muon selection criteria
				int iMuPlus = Reco_QQ_mupl_idx[iReco];
				int iMuMinus = Reco_QQ_mumi_idx[iReco];

				// HLT filters
				bool mupl_L2Filter = ((Reco_mu_trig[iMuPlus] & ((ULong64_t)pow(2, gL2FilterBit))) == ((ULong64_t)pow(2, gL2FilterBit)));
				bool mupl_L3Filter = ((Reco_mu_trig[iMuPlus] & ((ULong64_t)pow(2, gL3FilterBit))) == ((ULong64_t)pow(2, gL3FilterBit)));
				bool mumi_L2Filter = ((Reco_mu_trig[iMuMinus] & ((ULong64_t)pow(2, gL2FilterBit))) == ((ULong64_t)pow(2, gL2FilterBit)));
				bool mumi_L3Filter = ((Reco_mu_trig[iMuMinus] & ((ULong64_t)pow(2, gL3FilterBit))) == ((ULong64_t)pow(2, gL3FilterBit)));

				passHLTFilterMuons = (mupl_L2Filter && mumi_L3Filter) || (mupl_L3Filter && mumi_L2Filter) || (mupl_L3Filter && mumi_L3Filter);

				// global AND tracker muons
				trackerAndGlobalMuons = (Reco_mu_SelectionType[iMuPlus] & 2) && (Reco_mu_SelectionType[iMuPlus] & 8) && (Reco_mu_SelectionType[iMuMinus] & 2) && (Reco_mu_SelectionType[iMuMinus] & 8);

				// passing hybrid-soft Id
				hybridSoftMuons = (Reco_mu_nTrkWMea[iMuPlus] > 5) && (Reco_mu_nPixWMea[iMuPlus] > 0) && (fabs(Reco_mu_dxy[iMuPlus]) < 0.3) && (fabs(Reco_mu_dz[iMuPlus]) < 20.) && (Reco_mu_nTrkWMea[iMuMinus] > 5) && (Reco_mu_nPixWMea[iMuMinus] > 0) && (fabs(Reco_mu_dxy[iMuMinus]) < 0.3) && (fabs(Reco_mu_dz[iMuMinus]) < 20.);

				allGood = firesTrigger && isRecoMatched && dimuonMatching && goodVertexProba && passHLTFilterMuons && trackerAndGlobalMuons && hybridSoftMuons;

				// get muon coordinates

				TLorentzVector* Reco_mupl_LV = (TLorentzVector*)CloneArr_mu->At(iMuPlus);
				double Reco_mupl_eta = Reco_mupl_LV->Eta();
				double Reco_mupl_pt = Reco_mupl_LV->Pt();

				TLorentzVector* Reco_mumi_LV = (TLorentzVector*)CloneArr_mu->At(iMuMinus);
				double Reco_mumi_eta = Reco_mumi_LV->Eta();
				double Reco_mumi_pt = Reco_mumi_LV->Pt();

				TVector3 muPlus_CS = MuPlusVector_CollinsSoper(*(TLorentzVector*)CloneArr_QQ->At(iReco), *Reco_mupl_LV);
				double cosThetaCS = muPlus_CS.CosTheta();
				double phiCS = muPlus_CS.Phi() * 180 / TMath::Pi();

				TVector3 muPlus_HX = MuPlusVector_Helicity(*(TLorentzVector*)CloneArr_QQ->At(iReco), *Reco_mupl_LV);
				double cosThetaHX = muPlus_HX.CosTheta();
				double phiHX = muPlus_HX.Phi() * 180 / TMath::Pi();

				/// muon scale factors

				// muon trigger SF is tricky, need to know which muon passed which trigger filter
				bool mupl_isL2 = (mupl_L2Filter && !mupl_L3Filter);
				bool mupl_isL3 = (mupl_L2Filter && mupl_L3Filter);
				bool mumi_isL2 = (mumi_L2Filter && !mumi_L3Filter);
				bool mumi_isL3 = (mumi_L2Filter && mumi_L3Filter);

				if (mupl_isL2 && mumi_isL3) {
					dimuTrigWeight_nominal = tnp_weight_trg_pbpb(Reco_mupl_pt, Reco_mupl_eta, 2, indexNominal) * tnp_weight_trg_pbpb(Reco_mumi_pt, Reco_mumi_eta, 3, indexNominal);

					dimuTrigWeight_systUp = tnp_weight_trg_pbpb(Reco_mupl_pt, Reco_mupl_eta, 2, indexSystUp) * tnp_weight_trg_pbpb(Reco_mumi_pt, Reco_mumi_eta, 3, indexSystUp);

					dimuTrigWeight_systDown = tnp_weight_trg_pbpb(Reco_mupl_pt, Reco_mupl_eta, 2, indexSystDown) * tnp_weight_trg_pbpb(Reco_mumi_pt, Reco_mumi_eta, 3, indexSystDown);

					dimuTrigWeight_statUp = tnp_weight_trg_pbpb(Reco_mupl_pt, Reco_mupl_eta, 2, indexStatUp) * tnp_weight_trg_pbpb(Reco_mumi_pt, Reco_mumi_eta, 3, indexStatUp);

					dimuTrigWeight_statDown = tnp_weight_trg_pbpb(Reco_mupl_pt, Reco_mupl_eta, 2, indexStatDown) * tnp_weight_trg_pbpb(Reco_mumi_pt, Reco_mumi_eta, 3, indexStatDown);

				}

				else if (mupl_isL3 && mumi_isL2) {
					dimuTrigWeight_nominal = tnp_weight_trg_pbpb(Reco_mupl_pt, Reco_mupl_eta, 3, indexNominal) * tnp_weight_trg_pbpb(Reco_mumi_pt, Reco_mumi_eta, 2, indexNominal);

					dimuTrigWeight_systUp = tnp_weight_trg_pbpb(Reco_mupl_pt, Reco_mupl_eta, 3, indexSystUp) * tnp_weight_trg_pbpb(Reco_mumi_pt, Reco_mumi_eta, 2, indexSystUp);

					dimuTrigWeight_systDown = tnp_weight_trg_pbpb(Reco_mupl_pt, Reco_mupl_eta, 3, indexSystDown) * tnp_weight_trg_pbpb(Reco_mumi_pt, Reco_mumi_eta, 2, indexSystDown);

					dimuTrigWeight_statUp = tnp_weight_trg_pbpb(Reco_mupl_pt, Reco_mupl_eta, 3, indexStatUp) * tnp_weight_trg_pbpb(Reco_mumi_pt, Reco_mumi_eta, 2, indexStatUp);

					dimuTrigWeight_statDown = tnp_weight_trg_pbpb(Reco_mupl_pt, Reco_mupl_eta, 3, indexStatDown) * tnp_weight_trg_pbpb(Reco_mumi_pt, Reco_mumi_eta, 2, indexStatDown);

				}

				else if (mupl_isL3 && mumi_isL3) {
					dimuTrigWeight_nominal = DimuonL3TriggerWeight(Reco_mupl_pt, Reco_mupl_eta, Reco_mumi_pt, Reco_mumi_eta, indexNominal);

					dimuTrigWeight_systUp = DimuonL3TriggerWeight(Reco_mupl_pt, Reco_mupl_eta, Reco_mumi_pt, Reco_mumi_eta, indexSystUp);

					dimuTrigWeight_systDown = DimuonL3TriggerWeight(Reco_mupl_pt, Reco_mupl_eta, Reco_mumi_pt, Reco_mumi_eta, indexSystDown);

					dimuTrigWeight_statUp = DimuonL3TriggerWeight(Reco_mupl_pt, Reco_mupl_eta, Reco_mumi_pt, Reco_mumi_eta, indexStatUp);

					dimuTrigWeight_statDown = DimuonL3TriggerWeight(Reco_mupl_pt, Reco_mupl_eta, Reco_mumi_pt, Reco_mumi_eta, indexStatDown);
				}

				// dimuon efficiency weight = product of the total scale factors
				dimuWeight_nominal = tnp_weight_trk_pbpb(Reco_mupl_eta, indexNominal) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexNominal) * dimuTrigWeight_nominal;

				totalWeight = eventWeight * dimuWeight_nominal;
				hCS_nominal->FillWeighted(allGood, totalWeight, cosThetaCS, phiCS);
				hHX_nominal->FillWeighted(allGood, totalWeight, cosThetaHX, phiHX);

				/// variations for muon tracking SF (keeping the nominal efficiency for muon Id and trigger)

				// tracking, syst up
				dimuWeight_trk_systUp = tnp_weight_trk_pbpb(Reco_mupl_eta, indexSystUp) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexSystUp) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexNominal) * dimuTrigWeight_nominal;

				totalWeight = eventWeight * dimuWeight_trk_systUp;
				hCS_trk_systUp->FillWeighted(allGood, totalWeight, cosThetaCS, phiCS);
				hHX_trk_systUp->FillWeighted(allGood, totalWeight, cosThetaHX, phiHX);

				// tracking, syst down
				dimuWeight_trk_systDown = tnp_weight_trk_pbpb(Reco_mupl_eta, indexSystDown) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexSystDown) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexNominal) * dimuTrigWeight_nominal;

				totalWeight = eventWeight * dimuWeight_trk_systDown;
				hCS_trk_systDown->FillWeighted(allGood, totalWeight, cosThetaCS, phiCS);
				hHX_trk_systDown->FillWeighted(allGood, totalWeight, cosThetaHX, phiHX);

				// tracking, stat up
				dimuWeight_trk_statUp = tnp_weight_trk_pbpb(Reco_mupl_eta, indexStatUp) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexStatUp) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexNominal) * dimuTrigWeight_nominal;

				totalWeight = eventWeight * dimuWeight_trk_statUp;
				hCS_trk_statUp->FillWeighted(allGood, totalWeight, cosThetaCS, phiCS);
				hHX_trk_statUp->FillWeighted(allGood, totalWeight, cosThetaHX, phiHX);

				// tracking, stat down
				dimuWeight_trk_statDown = tnp_weight_trk_pbpb(Reco_mupl_eta, indexStatDown) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexStatDown) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexNominal) * dimuTrigWeight_nominal;

				totalWeight = eventWeight * dimuWeight_trk_statDown;
				hCS_trk_statDown->FillWeighted(allGood, totalWeight, cosThetaCS, phiCS);
				hHX_trk_statDown->FillWeighted(allGood, totalWeight, cosThetaHX, phiHX);

				/// variations for muon Id SF (keeping the nominal efficiency for tracking and trigger)

				// Id, syst up
				dimuWeight_muId_systUp = tnp_weight_trk_pbpb(Reco_mupl_eta, indexNominal) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexSystUp) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexSystUp) * dimuTrigWeight_nominal;

				totalWeight = eventWeight * dimuWeight_muId_systUp;
				hCS_muId_systUp->FillWeighted(allGood, totalWeight, cosThetaCS, phiCS);
				hHX_muId_systUp->FillWeighted(allGood, totalWeight, cosThetaHX, phiHX);

				// Id, syst down
				dimuWeight_muId_systDown = tnp_weight_trk_pbpb(Reco_mupl_eta, indexNominal) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexSystDown) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexSystDown) * dimuTrigWeight_nominal;

				totalWeight = eventWeight * dimuWeight_muId_systDown;
				hCS_muId_systDown->FillWeighted(allGood, totalWeight, cosThetaCS, phiCS);
				hHX_muId_systDown->FillWeighted(allGood, totalWeight, cosThetaHX, phiHX);

				// Id, stat up
				dimuWeight_muId_statUp = tnp_weight_trk_pbpb(Reco_mupl_eta, indexNominal) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexStatUp) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexStatUp) * dimuTrigWeight_nominal;

				totalWeight = eventWeight * dimuWeight_muId_statUp;
				hCS_muId_statUp->FillWeighted(allGood, totalWeight, cosThetaCS, phiCS);
				hHX_muId_statUp->FillWeighted(allGood, totalWeight, cosThetaHX, phiHX);

				// Id, stat down
				dimuWeight_muId_statDown = tnp_weight_trk_pbpb(Reco_mupl_eta, indexNominal) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexStatDown) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexStatDown) * dimuTrigWeight_nominal;

				totalWeight = eventWeight * dimuWeight_muId_statDown;
				hCS_muId_statDown->FillWeighted(allGood, totalWeight, cosThetaCS, phiCS);
				hHX_muId_statDown->FillWeighted(allGood, totalWeight, cosThetaHX, phiHX);

				/// variations for trigger SF (keeping the nominal efficiency for tracking and muon Id)

				// trigger, syst up
				dimuWeight_trig_systUp = tnp_weight_trk_pbpb(Reco_mupl_eta, indexNominal) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexNominal) * dimuTrigWeight_systUp;

				totalWeight = eventWeight * dimuWeight_trig_systUp;
				hCS_trig_systUp->FillWeighted(allGood, totalWeight, cosThetaCS, phiCS);
				hHX_trig_systUp->FillWeighted(allGood, totalWeight, cosThetaHX, phiHX);

				// trigger, syst down
				dimuWeight_trig_systDown = tnp_weight_trk_pbpb(Reco_mupl_eta, indexNominal) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexNominal) * dimuTrigWeight_systDown;

				totalWeight = eventWeight * dimuWeight_trig_systDown;
				hCS_trig_systDown->FillWeighted(allGood, totalWeight, cosThetaCS, phiCS);
				hHX_trig_systDown->FillWeighted(allGood, totalWeight, cosThetaHX, phiHX);

				// trigger, stat up
				dimuWeight_trig_statUp = tnp_weight_trk_pbpb(Reco_mupl_eta, indexNominal) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexNominal) * dimuTrigWeight_statUp;

				totalWeight = eventWeight * dimuWeight_trig_statUp;
				hCS_trig_statUp->FillWeighted(allGood, totalWeight, cosThetaCS, phiCS);
				hHX_trig_statUp->FillWeighted(allGood, totalWeight, cosThetaHX, phiHX);

				// trigger, stat down
				dimuWeight_trig_statDown = tnp_weight_trk_pbpb(Reco_mupl_eta, indexNominal) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexNominal) * dimuTrigWeight_statDown;

				totalWeight = eventWeight * dimuWeight_trig_statDown;
				hCS_trig_statDown->FillWeighted(allGood, totalWeight, cosThetaCS, phiCS);
				hHX_trig_statDown->FillWeighted(allGood, totalWeight, cosThetaHX, phiHX);
			}
		} // end of gen upsilon loop
	}

	cout << endl;

	/// display the nominal results

	gStyle->SetPadLeftMargin(.15);
	gStyle->SetPadRightMargin(0.18);
	gStyle->SetPalette(kRainBow);
	gStyle->SetNumberContours(256);

	const char* legendText = Form("cent. %d-%d%%, |y^{#mu#mu}| < 2.4, %d < p_{T}^{#mu#mu} < %d GeV", gCentralityBinMin, gCentralityBinMax, ptMin, ptMax);

	TLatex* legend = new TLatex(.5, .88, legendText);
	legend->SetTextAlign(22);
	legend->SetTextSize(0.042);

	// Collins-Soper

	auto* canvasCS = new TCanvas("canvasCS", "", 700, 600);
	hCS_nominal->Draw("COLZ");

	CMS_lumi(canvasCS, Form("#varUpsilon(%dS) Hydjet-embedded MC", iState));

	legend->DrawLatexNDC(.5, .88, legendText);

	gPad->Update();

	hCS_nominal->GetPaintedHistogram()->GetYaxis()->SetRangeUser(-190, 240);
	hCS_nominal->GetPaintedHistogram()->GetZaxis()->SetRangeUser(0, 1);

	canvasCS->SaveAs(Form("EfficiencyMaps/%dS/NominalEff_CS_pt%dto%d.png", iState, ptMin, ptMax), "RECREATE");

	// Helicity
	auto* canvasHX = new TCanvas("canvasHX", "", 700, 600);
	hHX_nominal->Draw("COLZ");

	CMS_lumi(canvasHX, Form("#varUpsilon(%dS) Hydjet-embedded MC", iState));

	legend->DrawLatexNDC(.5, .88, legendText);

	gPad->Update();

	hHX_nominal->GetPaintedHistogram()->GetYaxis()->SetRangeUser(-190, 240);
	hHX_nominal->GetPaintedHistogram()->GetZaxis()->SetRangeUser(0, 1);

	canvasHX->SaveAs(Form("EfficiencyMaps/%dS/NominalEff_HX_pt%dto%d.png", iState, ptMin, ptMax), "RECREATE");

	/// compute the systematics in this macro since we have all the ingredients for that
	// instructions can be found here: https://twiki.cern.ch/twiki/pub/CMS/HIMuonTagProbe/TnpHeaderFile.pdf#page=5

	// store the RELATIVE SYSTEMATIC UNCERTAINTIES (not in % though) with respect to the nominal efficiency, more useful at the end

	gStyle->SetPadRightMargin(0.19);
	gStyle->SetTitleOffset(1.3, "Z");

	// Collins-Soper
	auto* hSystCS = RelSystEffHist(hCS_nominal, hCS_trk_systUp, hCS_trk_systDown, hCS_muId_systUp, hCS_muId_systDown, hCS_trig_systUp, hCS_trig_systDown, hCS_trk_statUp, hCS_trk_statDown, hCS_muId_statUp, hCS_muId_statDown, hCS_trig_statUp, hCS_trig_statDown);
	hSystCS->SetName(Form("RelatSystEff_CS_pt%dto%d", ptMin, ptMax));
	hSystCS->SetTitle(";cos #theta_{CS}; #varphi_{CS} (#circ);Relative efficiency uncertainty");

	auto* canvasCSsyst = new TCanvas("canvasCSsyst", "", 700, 600);
	hSystCS->Draw("COLZ");

	CMS_lumi(canvasCSsyst, Form("#varUpsilon(%dS) Hydjet-embedded MC", iState));

	legend->DrawLatexNDC(.5, .88, legendText);

	gPad->Update();

	hSystCS->GetYaxis()->SetRangeUser(-190, 240);

	canvasCSsyst->SaveAs(Form("EfficiencyMaps/%dS/RelSystEff_CS_pt%dto%d.png", iState, ptMin, ptMax), "RECREATE");

	// Helicity
	auto* hSystHX = RelSystEffHist(hHX_nominal, hHX_trk_systUp, hHX_trk_systDown, hHX_muId_systUp, hHX_muId_systDown, hHX_trig_systUp, hHX_trig_systDown, hHX_trk_statUp, hHX_trk_statDown, hHX_muId_statUp, hHX_muId_statDown, hHX_trig_statUp, hHX_trig_statDown);
	hSystHX->SetName(Form("RelatSystEff_HX_pt%dto%d", ptMin, ptMax));
	hSystHX->SetTitle(";cos #theta_{HX}; #varphi_{HX} (#circ);Relative efficiency uncertainty");

	auto* canvasHXsyst = new TCanvas("canvasHXsyst", "", 700, 600);
	hSystHX->Draw("COLZ");

	CMS_lumi(canvasHXsyst, Form("#varUpsilon(%dS) Hydjet-embedded MC", iState));

	legend->DrawLatexNDC(.5, .88, legendText);

	gPad->Update();

	hSystHX->GetYaxis()->SetRangeUser(-190, 240);

	canvasHXsyst->SaveAs(Form("EfficiencyMaps/%dS/RelSystEff_HX_pt%dto%d.png", iState, ptMin, ptMax), "RECREATE");

	/// save the nominal efficiency results and the corresponding systematics in a file for later usage
	const char* outputFileName = Form("EfficiencyMaps/%dS/EfficiencyResults.root", iState);
	TFile outputFile(outputFileName, "RECREATE");

	hCS_nominal->Write();
	hHX_nominal->Write();
	hSystCS->Write();
	hSystHX->Write();
	outputFile.Close();

	cout << endl
	     << "Nominal efficiency and corresponding systematic uncertainty maps saved in " << outputFileName << endl;
}
