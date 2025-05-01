#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/Parameters/PhaseSpace.h"

#include "AccEffHelpers.h"

#include "../Tools/Style/Legends.h"

#include "../Tools/Parameters/CentralityValues.h"
#include "../Tools/Parameters/EfficiencyWeights.h"
#include "../Tools/Parameters/MuonScaleFactors.h"

#include "../ReferenceFrameTransform/Transformations.h"

void GenRecoResolution(Int_t ptMin = 0, Int_t ptMax = 2, Int_t iState = gUpsilonState, Float_t lambdaTheta = 0, Float_t lambdaPhi = 0, Float_t lambdaThetaPhi = 0, Bool_t isPhiFolded = kFALSE, TString muonAccName = "UpsilonTriggerThresholds") { // accName = "MuonSimpleAcc", "MuonWithin2018PbPbAcc", or "MuonUpsilonTriggerAcc"

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

	// TString MuonAccName = "";

	Bool_t allGood, firesTrigger, isRecoMatched, dimuonMatching, goodVertexProba, passHLTFilterMuons, trackerAndGlobalMuons, hybridSoftMuons;

	Int_t hiBin;

	TH2D* genCosThetaPhiHistCS = new TH2D("genCosThetaPhiHistCS", "genCosThetaPhiHistCS", NCosThetaBins, gCosThetaBinning, NFullPhiBins, gFullPhiBinning);
	TH2D* genCosThetaPhiHistHX = new TH2D("genCosThetaPhiHistHX", "genCosThetaPhiHistHX", NCosThetaBins, gCosThetaBinning, NFullPhiBins, gFullPhiBinning);
	
	TH2D* recoCosThetaPhiHistCS = new TH2D("recoCosThetaPhiHistCS", "recoCosThetaPhiHistCS", NCosThetaBins, gCosThetaBinning, NFullPhiBins, gFullPhiBinning);
	TH2D* recoCosThetaPhiHistHX = new TH2D("recoCosThetaPhiHistHX", "recoCosThetaPhiHistHX", NCosThetaBins, gCosThetaBinning, NFullPhiBins, gFullPhiBinning);

	TH2D* ratioCosThetaPhiHistCS = new TH2D("ratioCosThetaPhiHistCS", "ratioCosThetaPhiHistCS", NCosThetaBins, gCosThetaBinning, NFullPhiBins, gFullPhiBinning);
	TH2D* ratioCosThetaPhiHistHX = new TH2D("ratioCosThetaPhiHistHX", "ratioCosThetaPhiHistHX", NCosThetaBins, gCosThetaBinning, NFullPhiBins, gFullPhiBinning);

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

				// dimuon efficiency weight = product of the total scale factors
				dimuWeight_nominal = tnp_weight_trk_pbpb(Reco_mupl_eta, indexNominal) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexNominal) * dimuTrigWeight_nominal;

				if (isPhiFolded == kTRUE) {
					weightCS = 1 + lambdaTheta * TMath::Power(muPlus_CS_gen.CosTheta(), 2) + lambdaPhi * TMath::Power(std::sin(muPlus_CS_gen.Theta()), 2) * std::cos(2 * fabs(muPlus_CS_gen.Phi())) + lambdaThetaPhi * std::sin(2 * muPlus_CS_gen.Theta()) * std::cos(fabs(muPlus_CS_gen.Phi()));
					weightHX = 1 + lambdaTheta * TMath::Power(muPlus_HX_gen.CosTheta(), 2) + lambdaPhi * TMath::Power(std::sin(muPlus_HX_gen.Theta()), 2) * std::cos(2 * fabs(muPlus_HX_gen.Phi())) + lambdaThetaPhi * std::sin(2 * muPlus_HX_gen.Theta()) * std::cos(fabs(muPlus_HX_gen.Phi()));
				}

				else {
					weightCS = 1 + lambdaTheta * TMath::Power(muPlus_CS_gen.CosTheta(), 2) + lambdaPhi * TMath::Power(std::sin(muPlus_CS_gen.Theta()), 2) * std::cos(2 * muPlus_CS_gen.Phi()) + lambdaThetaPhi * std::sin(2 * muPlus_CS_gen.Theta()) * std::cos(muPlus_CS_gen.Phi());
					weightHX = 1 + lambdaTheta * TMath::Power(muPlus_HX_gen.CosTheta(), 2) + lambdaPhi * TMath::Power(std::sin(muPlus_HX_gen.Theta()), 2) * std::cos(2 * muPlus_HX_gen.Phi()) + lambdaThetaPhi * std::sin(2 * muPlus_HX_gen.Theta()) * std::cos(muPlus_HX_gen.Phi());
				}

				// total weight
				totalWeightLab = eventWeight * dimuonPtWeight * dimuWeight_nominal;
				totalWeightCS = eventWeight * dimuonPtWeight * dimuWeight_nominal * weightCS;
				totalWeightHX = eventWeight * dimuonPtWeight * dimuWeight_nominal * weightHX;

				// fill histograms
				if (reco_QQ_pt > 2 && reco_QQ_pt < 6) {
					genCosThetaPhiHistCS->Fill(muPlus_CS_gen.CosTheta(), muPlus_CS_gen.Phi() * 180 / TMath::Pi());
					genCosThetaPhiHistHX->Fill(muPlus_HX_gen.CosTheta(), muPlus_HX_gen.Phi() * 180 / TMath::Pi());
					recoCosThetaPhiHistCS->Fill(cosThetaCS, phiCS);
					recoCosThetaPhiHistHX->Fill(cosThetaHX, phiHX);

					// cout << "gen CS CosTheta: " << muPlus_CS_gen.CosTheta() << ", gen CS Phi: " << muPlus_CS_gen.Phi() * 180 / TMath::Pi() << endl;
					// cout << "gen HX CosTheta: " << muPlus_HX_gen.CosTheta() << ", gen HX Phi: " << muPlus_HX_gen.Phi() * 180 / TMath::Pi() << endl;
					// cout << cosThetaCS << " " << phiCS << " " << cosThetaHX << " " << phiHX << endl;
				}
				
			}
		}
	}

	// draw histograms
	TCanvas* genCanvasHX = new TCanvas("genCanvasHX", "genCanvasHX", 800, 600);
	genCosThetaPhiHistHX->Draw("COLZ");

	TCanvas* recoCanvasHX = new TCanvas("recoCanvasHX", "recoCanvasHX", 800, 600);
	recoCosThetaPhiHistHX->Draw("COLZ");

	TCanvas* ratioCanvasHX = new TCanvas("ratioCanvasHX", "ratioCanvasHX", 800, 600);
	ratioCosThetaPhiHistHX->Divide(recoCosThetaPhiHistHX, genCosThetaPhiHistHX);
	ratioCosThetaPhiHistHX->Draw("COLZ");
}