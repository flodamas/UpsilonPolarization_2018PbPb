#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/Parameters/PhaseSpace.h"

#include "AccEffHelpers.h"
#include "../Polarization/PolarFitHelpers.h"

#include "../Tools/Style/Legends.h"

#include "../Tools/Parameters/CentralityValues.h"
#include "../Tools/Parameters/EfficiencyWeights.h"
#include "../Tools/Parameters/MuonScaleFactors.h"

#include "../ReferenceFrameTransform/Transformations.h"

void GenRecoResolution(Int_t iState = gUpsilonState, Float_t lambdaTheta = 0, Float_t lambdaPhi = 0, Float_t lambdaThetaPhi = 0, Bool_t isPhiFolded = kFALSE, TString muonAccName = "UpsilonTriggerThresholds") { // accName = "MuonSimpleAcc", "MuonWithin2018PbPbAcc", or "MuonUpsilonTriggerAcc"

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

    	// for a given type of muon scale factor (i.e., tracking, muId, trigger) we need to compute the efficiency with up and down variations, for stat and syst uncertainties

	int indexNominal = 0;
	int indexSystUp = -1, indexSystDown = -2;
	int indexStatUp = +1, indexStatDown = +2;

	double dimuWeight_nominal = 0;
	double dimuWeight_trk_systUp, dimuWeight_trk_systDown, dimuWeight_trk_statUp, dimuWeight_trk_statDown;
	double dimuWeight_muId_systUp, dimuWeight_muId_systDown, dimuWeight_muId_statUp, dimuWeight_muId_statDown;
	double dimuWeight_trig_systUp, dimuWeight_trig_systDown, dimuWeight_trig_statUp, dimuWeight_trig_statDown;
	double dimuWeight_total_systUp, dimuWeight_total_systDown, dimuWeight_total_statUp, dimuWeight_total_statDown;

	// loop variables
	TLorentzVector* genLorentzVector = new TLorentzVector();
	TLorentzVector* gen_QQ_LV = new TLorentzVector();
	TLorentzVector* gen_mupl_LV = new TLorentzVector();
	TLorentzVector* gen_mumi_LV = new TLorentzVector();
	TLorentzVector* recoLorentzVector = new TLorentzVector();

	double eventWeight, dimuonPtWeight, totalWeightCS, totalWeightHX, totalWeightLab;
	double dimuTrigWeight_nominal = -1, 
		   dimuTrigWeight_systUp = -1, dimuTrigWeight_systDown = -1, 
		   dimuTrigWeight_statUp = -1, dimuTrigWeight_statDown = -1;

	Float_t weightCS = 0, weightHX = 0;

	Bool_t allGood, firesTrigger, isRecoMatched, dimuonMatching, goodVertexProba, passHLTFilterMuons, trackerAndGlobalMuons, hybridSoftMuons;

	Int_t hiBin;

	Int_t ptBinIdx = 0;

	/// define the 2D histograms
	TH2D* genCosThetaPhiHistCS[NPtBins];
	TH2D* genCosThetaPhiHistHX[NPtBins];
	
	TH2D* recoCosThetaPhiHistCS[NPtBins];
	TH2D* recoCosThetaPhiHistHX[NPtBins];

	TH2D* ratioCosThetaPhiHistCS[NPtBins];
	TH2D* ratioCosThetaPhiHistHX[NPtBins];

	for (int iPtBin = 0; iPtBin < NPtBins; iPtBin++) {
		genCosThetaPhiHistCS[iPtBin] = new TH2D(Form("genCosThetaPhiHistCS%d", iPtBin), "", NCosThetaBins, gCosThetaBinning, NFullPhiBins, gFullPhiBinning);
		genCosThetaPhiHistHX[iPtBin] = new TH2D(Form("genCosThetaPhiHistHX%d", iPtBin), "", NCosThetaBins, gCosThetaBinning, NFullPhiBins, gFullPhiBinning);
		recoCosThetaPhiHistCS[iPtBin] = new TH2D(Form("recoCosThetaPhiHistCS%d", iPtBin), "", NCosThetaBins, gCosThetaBinning, NFullPhiBins, gFullPhiBinning);
		recoCosThetaPhiHistHX[iPtBin] = new TH2D(Form("recoCosThetaPhiHistHX%d", iPtBin), "", NCosThetaBins, gCosThetaBinning, NFullPhiBins, gFullPhiBinning);
		ratioCosThetaPhiHistCS[iPtBin] = new TH2D(Form("ratioCosThetaPhiHistCS%d", iPtBin), "", NCosThetaBins, gCosThetaBinning, NFullPhiBins, gFullPhiBinning);
		ratioCosThetaPhiHistHX[iPtBin] = new TH2D(Form("ratioCosThetaPhiHistHX%d", iPtBin), "", NCosThetaBins, gCosThetaBinning, NFullPhiBins, gFullPhiBinning);
	}

	Long64_t totEntries = OniaTree->GetEntries();

    for (Long64_t iEvent = 0; iEvent < totEntries; iEvent++) {
		if (iEvent % 10000 == 0) {
			cout << Form("\rProcessing event %lld / %lld (%.0f%%)", iEvent, totEntries, 100. * iEvent / totEntries) << flush;
		}

		OniaTree->GetEntry(iEvent);
		//genLorentzVector->Clear();

		// if (iEvent > 100) break; // for testing

		// event selection

		if (Centrality >= 2 * gCentralityBinMax) continue; // discard events with centrality >= 90% in 2018 data

		firesTrigger = ((HLTriggers & (ULong64_t)(1 << (gUpsilonHLTBit - 1))) == (ULong64_t)(1 << (gUpsilonHLTBit - 1)));

		hiBin = GetHiBinFromhiHF(HFmean);

		eventWeight = Gen_weight * FindNcoll(Centrality); // * Get_zPV_weight(zVtx);

		// loop over all gen upsilons
		for (int iGen = 0; iGen < Gen_QQ_size; iGen++) {   
			gen_QQ_LV = (TLorentzVector*)Gen_QQ_4mom->At(iGen);
			// fiducial region
			//if (genLorentzVector->Pt() < ptMin || genLorentzVector->Pt() > ptMax) continue; // pt bin of interest

			if (fabs(gen_QQ_LV->Rapidity()) < gRapidityMin || fabs(gen_QQ_LV->Rapidity()) > gRapidityMax) continue;

			if (gen_QQ_LV->Pt() > gPtMax) continue;

			// single-muon acceptance

			// positive muon first
			gen_mupl_LV = (TLorentzVector*)Gen_mu_4mom->At(Gen_QQ_mupl_idx[iGen]);

			if (!MuonKinematicsWithinLimits(*gen_mupl_LV, muonAccName)) continue;

			// then negative muon
			gen_mumi_LV = (TLorentzVector*)Gen_mu_4mom->At(Gen_QQ_mumi_idx[iGen]);

			if (!MuonKinematicsWithinLimits(*gen_mumi_LV, muonAccName)) continue;

			// go to reco level
			Int_t iReco = Gen_QQ_whichRec[iGen];

			if (Reco_QQ_sign[iReco] != 0) continue; // only opposite-sign muon pairs

			/// all the reconstructed upsilons must pass the conditions below!

			isRecoMatched = iReco > -1;

			if (isRecoMatched) {
				recoLorentzVector = (TLorentzVector*)CloneArr_QQ->At(iReco);
				double reco_QQ_pt = recoLorentzVector->Pt();

				dimuonPtWeight = Get_RecoPtWeight(recoLorentzVector->Rapidity(), reco_QQ_pt);

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

				/// numerator for the efficiency
				allGood = firesTrigger && isRecoMatched && dimuonMatching && goodVertexProba && passHLTFilterMuons && trackerAndGlobalMuons && hybridSoftMuons;

				// get muon coordinates
				TLorentzVector* Reco_mupl_LV = (TLorentzVector*)CloneArr_mu->At(iMuPlus);
				double Reco_mupl_eta = Reco_mupl_LV->Eta();
				double Reco_mupl_pt = Reco_mupl_LV->Pt();

				TLorentzVector* Reco_mumi_LV = (TLorentzVector*)CloneArr_mu->At(iMuMinus);
				double Reco_mumi_eta = Reco_mumi_LV->Eta();
				double Reco_mumi_pt = Reco_mumi_LV->Pt();

				// double cosThetaLab_gen = gen_mupl_LV->CosTheta();
				double cosThetaLab = Reco_mupl_LV->CosTheta();

				// double phiLab_gen = 0;
				double phiLab = 0;

				if (isPhiFolded == kTRUE) {
					// phiLab_gen = fabs(gen_mupl_LV->Phi() * 180 / TMath::Pi());
					phiLab = fabs(Reco_mupl_LV->Phi() * 180 / TMath::Pi());
				} else {
					// phiLab_gen = gen_mupl_LV->Phi() * 180 / TMath::Pi();
					phiLab = Reco_mupl_LV->Phi() * 180 / TMath::Pi();
				}

				TVector3 muPlus_CS_gen = MuPlusVector_CollinsSoper(*gen_QQ_LV, *gen_mupl_LV);
				TVector3 muPlus_CS = MuPlusVector_CollinsSoper(*recoLorentzVector, *Reco_mupl_LV);

				// double cosThetaCS_gen = muPlus_CS_gen.CosTheta();
				double cosThetaCS = muPlus_CS.CosTheta();
				// double phiCS_gen = 0;
				double phiCS = 0;

				if (isPhiFolded == kTRUE) {
					// phiCS_gen = fabs(muPlus_CS_gen.Phi() * 180 / TMath::Pi());
					phiCS = fabs(muPlus_CS.Phi() * 180 / TMath::Pi());
				} else {
					// phiCS_gen = muPlus_CS_gen.Phi() * 180 / TMath::Pi();
					phiCS = muPlus_CS.Phi() * 180 / TMath::Pi();
				}

				TVector3 muPlus_HX_gen = MuPlusVector_Helicity(*gen_QQ_LV, *gen_mupl_LV);
				TVector3 muPlus_HX = MuPlusVector_Helicity(*recoLorentzVector, *Reco_mupl_LV);

				// double cosThetaHX_gen = muPlus_HX_gen.CosTheta();
				double cosThetaHX = muPlus_HX.CosTheta();
				// double phiHX_gen = 0;
				double phiHX = 0;

				if (isPhiFolded == kTRUE) {
					// phiHX_gen = fabs(muPlus_HX_gen.Phi() * 180 / TMath::Pi());
					phiHX = fabs(muPlus_HX.Phi() * 180 / TMath::Pi());
				} else {
					phiHX = muPlus_HX.Phi() * 180 / TMath::Pi();
					// phiHX_gen = muPlus_HX_gen.Phi() * 180 / TMath::Pi();
				}

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

				/// dimuon efficiency weight = product of the total scale factors
				dimuWeight_nominal = tnp_weight_trk_pbpb(Reco_mupl_eta, indexNominal) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexNominal) * dimuTrigWeight_nominal;

				if (isPhiFolded == kTRUE) {
					weightCS = 1 + lambdaTheta * TMath::Power(muPlus_CS_gen.CosTheta(), 2) + lambdaPhi * TMath::Power(std::sin(muPlus_CS_gen.Theta()), 2) * std::cos(2 * fabs(muPlus_CS_gen.Phi())) + lambdaThetaPhi * std::sin(2 * muPlus_CS_gen.Theta()) * std::cos(fabs(muPlus_CS_gen.Phi()));
					weightHX = 1 + lambdaTheta * TMath::Power(muPlus_HX_gen.CosTheta(), 2) + lambdaPhi * TMath::Power(std::sin(muPlus_HX_gen.Theta()), 2) * std::cos(2 * fabs(muPlus_HX_gen.Phi())) + lambdaThetaPhi * std::sin(2 * muPlus_HX_gen.Theta()) * std::cos(fabs(muPlus_HX_gen.Phi()));
				}

				else {
					weightCS = 1 + lambdaTheta * TMath::Power(muPlus_CS_gen.CosTheta(), 2) + lambdaPhi * TMath::Power(std::sin(muPlus_CS_gen.Theta()), 2) * std::cos(2 * muPlus_CS_gen.Phi()) + lambdaThetaPhi * std::sin(2 * muPlus_CS_gen.Theta()) * std::cos(muPlus_CS_gen.Phi());
					weightHX = 1 + lambdaTheta * TMath::Power(muPlus_HX_gen.CosTheta(), 2) + lambdaPhi * TMath::Power(std::sin(muPlus_HX_gen.Theta()), 2) * std::cos(2 * muPlus_HX_gen.Phi()) + lambdaThetaPhi * std::sin(2 * muPlus_HX_gen.Theta()) * std::cos(muPlus_HX_gen.Phi());
				}

				/// total weight
				totalWeightLab = eventWeight * dimuonPtWeight * dimuWeight_nominal;
				totalWeightCS = eventWeight * dimuonPtWeight * dimuWeight_nominal * weightCS;
				totalWeightHX = eventWeight * dimuonPtWeight * dimuWeight_nominal * weightHX;

				/// set the pt bin index
				if (reco_QQ_pt > gPtBinning[0] && reco_QQ_pt <= gPtBinning[1]) ptBinIdx = 0;
				else if (reco_QQ_pt > gPtBinning[1] && reco_QQ_pt <= gPtBinning[2]) ptBinIdx = 1;
				else if (reco_QQ_pt > gPtBinning[2] && reco_QQ_pt <= gPtBinning[3]) ptBinIdx = 2;
				else if (reco_QQ_pt > gPtBinning[3] && reco_QQ_pt <= gPtBinning[4]) ptBinIdx = 3;
				else continue;

				/// fill histograms
				genCosThetaPhiHistCS[ptBinIdx]->Fill(muPlus_CS_gen.CosTheta(), muPlus_CS_gen.Phi() * 180 / TMath::Pi());
				genCosThetaPhiHistHX[ptBinIdx]->Fill(muPlus_HX_gen.CosTheta(), muPlus_HX_gen.Phi() * 180 / TMath::Pi());
				
				recoCosThetaPhiHistCS[ptBinIdx]->Fill(cosThetaCS, phiCS);
				recoCosThetaPhiHistHX[ptBinIdx]->Fill(cosThetaHX, phiHX);
			}
		}
	}

	/// save the histograms
	gSystem->mkdir("GenRecoResolution", kTRUE);
	TFile *outputFile = new TFile("GenRecoResolution/GenRecoResolution.root", "RECREATE");
	outputFile->cd();
	
	/// write histograms to file
	for (int iPtBin = 0; iPtBin < NPtBins; iPtBin++) {
		ratioCosThetaPhiHistCS[iPtBin]->Divide(recoCosThetaPhiHistCS[iPtBin], genCosThetaPhiHistCS[iPtBin]);
		ratioCosThetaPhiHistHX[iPtBin]->Divide(recoCosThetaPhiHistHX[iPtBin], genCosThetaPhiHistHX[iPtBin]);

		genCosThetaPhiHistCS[iPtBin]->Write();
		genCosThetaPhiHistHX[iPtBin]->Write();
		recoCosThetaPhiHistCS[iPtBin]->Write();
		recoCosThetaPhiHistHX[iPtBin]->Write();
		ratioCosThetaPhiHistCS[iPtBin]->Write();
		ratioCosThetaPhiHistHX[iPtBin]->Write();
	}

	outputFile->Close();

	// draw an example histograms in HX frame (2 < pt < 6 GeV/c)
	TCanvas* genCanvasHX = new TCanvas("genCanvasHX", "genCanvasHX", 800, 600);
	genCosThetaPhiHistHX[1]->Draw("COLZ");

	TCanvas* recoCanvasHX = new TCanvas("recoCanvasHX", "recoCanvasHX", 800, 600);
	recoCosThetaPhiHistHX[1]->Draw("COLZ");

	TCanvas* ratioCanvasHX = new TCanvas("ratioCanvasHX", "ratioCanvasHX", 800, 600);
	ratioCosThetaPhiHistHX[1]->Divide(recoCosThetaPhiHistHX[1], genCosThetaPhiHistHX[1]);
	ratioCosThetaPhiHistHX[1]->Draw("COLZ");

	return;
}

void draw2DGenRecoResolution() {

	std::vector<Double_t> cosThetaEdges = setCosThetaBinEdges(NCosThetaBins, gCosThetaMin, gCosThetaMax);

	std::vector<Double_t> phiEdges = setPhiBinEdges(NFullPhiBins, gFullPhiMin, gFullPhiMax);

	TFile* file = TFile::Open("./GenRecoResolution/GenRecoResolution.root", "READ");
	if (!file) {
		cout << "File GenRecoResolution.root not found. Check the directory of the file." << endl;
		return;
	}

	cout << "File GenRecoResolution.root opened" << endl;
	
	// define the 2D histograms
	TH2D* genCosThetaPhiHistCS[NPtBins];
	TH2D* genCosThetaPhiHistHX[NPtBins];
	TH2D* recoCosThetaPhiHistCS[NPtBins];
	TH2D* recoCosThetaPhiHistHX[NPtBins];
	TH2D* ratioCosThetaPhiHistCS[NPtBins];
	TH2D* ratioCosThetaPhiHistHX[NPtBins];

	TCanvas* ratioCosThetaPhiCanvasCS[NPtBins];
	TCanvas* ratioCosThetaPhiCanvasHX[NPtBins];
	
	for (int iPtBin = 0; iPtBin < NPtBins; iPtBin++) {
		/// get the histograms
		genCosThetaPhiHistCS[iPtBin] = (TH2D*)file->Get(Form("genCosThetaPhiHistCS%d", iPtBin));
		genCosThetaPhiHistHX[iPtBin] = (TH2D*)file->Get(Form("genCosThetaPhiHistHX%d", iPtBin));
		recoCosThetaPhiHistCS[iPtBin] = (TH2D*)file->Get(Form("recoCosThetaPhiHistCS%d", iPtBin));
		recoCosThetaPhiHistHX[iPtBin] = (TH2D*)file->Get(Form("recoCosThetaPhiHistHX%d", iPtBin));
		ratioCosThetaPhiHistCS[iPtBin] = (TH2D*)file->Get(Form("ratioCosThetaPhiHistCS%d", iPtBin));
		ratioCosThetaPhiHistHX[iPtBin] = (TH2D*)file->Get(Form("ratioCosThetaPhiHistHX%d", iPtBin));

		/// draw ratio histograms CS
		ratioCosThetaPhiCanvasCS[iPtBin] = draw2DMap(ratioCosThetaPhiHistCS[iPtBin], "CS", NCosThetaBins, cosThetaEdges, NFullPhiBins, phiEdges, kFALSE, kFALSE, 1, kFALSE);
		display2DMapContents(ratioCosThetaPhiHistCS[iPtBin], NCosThetaBins, NFullPhiBins, kFALSE, 0.04, kWhite, 2);

		ratioCosThetaPhiHistCS[iPtBin]->GetZaxis()->SetTitleOffset(1.6);
		ratioCosThetaPhiHistCS[iPtBin]->GetZaxis()->SetTitleSize(0.05);
		ratioCosThetaPhiHistCS[iPtBin]->GetZaxis()->SetTitle("#varUpsilon(1S) N_{RECO} / N_{GEN}");

		// ratioCosThetaPhiCanvasCS[iPtBin]->SetRightMargin(0.20);
		// gStyle->SetPadRightMargin(0.2);

		TPaveText* kinematicsTextCS = KinematicsText_v2(gCentralityBinMin, gCentralityBinMax, (int)gPtBinning[iPtBin], (int)gPtBinning[iPtBin + 1]);  // (HX, 2to6, -0.42to-0.14,
		kinematicsTextCS->Draw();

		/// save the histograms
		ratioCosThetaPhiCanvasCS[iPtBin]->SaveAs(Form("./GenRecoResolution/ratioCosThetaPhiCS_pt%dto%d.png", (int)gPtBinning[iPtBin], (int)gPtBinning[iPtBin + 1]));

		/// draw ratio histograms HX
		ratioCosThetaPhiCanvasHX[iPtBin] = draw2DMap(ratioCosThetaPhiHistHX[iPtBin], "HX", NCosThetaBins, cosThetaEdges, NFullPhiBins, phiEdges, kFALSE, kFALSE, 1, kFALSE);
		display2DMapContents(ratioCosThetaPhiHistHX[iPtBin], NCosThetaBins, NFullPhiBins, kFALSE, 0.04, kWhite, 2);
		
		ratioCosThetaPhiHistHX[iPtBin]->GetZaxis()->SetTitleOffset(1.8);
		ratioCosThetaPhiHistHX[iPtBin]->GetZaxis()->SetTitleSize(0.05);
		ratioCosThetaPhiHistHX[iPtBin]->GetZaxis()->SetTitle("#varUpsilon(1S) N_{RECO} / N_{GEN}");
		
		// ratioCosThetaPhiCanvasHX[iPtBin]->SetRightMargin(0.20);

		// gStyle->SetPadRightMargin(0.2);

		TPaveText* kinematicsTextHX = KinematicsText_v2(gCentralityBinMin, gCentralityBinMax, (int)gPtBinning[iPtBin], (int)gPtBinning[iPtBin + 1]);  // (HX, 2to6, -0.42to-0.14,
		kinematicsTextHX->Draw();

		/// save the histograms
		ratioCosThetaPhiCanvasHX[iPtBin]->SaveAs(Form("./GenRecoResolution/ratioCosThetaPhiHX_pt%dto%d.png", (int)gPtBinning[iPtBin], (int)gPtBinning[iPtBin + 1]));
	}

	return;
}