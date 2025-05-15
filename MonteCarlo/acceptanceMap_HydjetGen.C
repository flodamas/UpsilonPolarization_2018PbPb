#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/Parameters/PhaseSpace.h"

#include "AccEffHelpers.h"

#include "../Polarization/PolarFitHelpers.h"

//#include "../Tools/Style/FitDistributions.h"
#include "../Tools/Style/Legends.h"

#include "../Tools/Parameters/CentralityValues.h"
#include "../Tools/Parameters/EfficiencyWeights.h"
#include "../Tools/Parameters/MuonScaleFactors.h"

#include "../ReferenceFrameTransform/Transformations.h"

void acceptanceMap_HydjetGen(Bool_t isPhiFolded = kTRUE, TString muonAccName = gMuonAccName, Float_t lambdaTheta = 0, Float_t lambdaPhi = 0, Float_t lambdaThetaPhi = 0, Int_t iState = 1) {
    // Read HydjetGen file with polarization weights
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

    // reco-level variables
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

	// (cos theta, phi, pT) 3D maps for final acceptance correction, variable size binning for the stats
	TEfficiency* accMatrixLab = TEfficiency3D(NominalTEfficiency3DName("Lab", lambdaTheta, lambdaPhi, lambdaThetaPhi), "Lab", iState, isPhiFolded);
	TEfficiency* accMatrixCS = TEfficiency3D(NominalTEfficiency3DName("CS", lambdaTheta, lambdaPhi, lambdaThetaPhi), "CS", iState, isPhiFolded);
	TEfficiency* accMatrixHX = TEfficiency3D(NominalTEfficiency3DName("HX", lambdaTheta, lambdaPhi, lambdaThetaPhi), "HX", iState, isPhiFolded);

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
	double dimuWeight_total_systUp, dimuWeight_total_systDown, dimuWeight_total_statUp, dimuWeight_total_statDown;

	// loop variables
	TLorentzVector* genLorentzVector = new TLorentzVector();
	TLorentzVector* gen_QQ_LV = new TLorentzVector();
	TLorentzVector* gen_mupl_LV = new TLorentzVector();
	TLorentzVector* gen_mumi_LV = new TLorentzVector();
	TLorentzVector* recoLorentzVector = new TLorentzVector();

    Bool_t withinAcceptance;

	Double_t cosThetaCS_gen, phiCS_gen, cosThetaHX_gen, phiHX_gen, cosThetaLab_gen, phiLab_gen;
	Double_t cosThetaCS_reco, phiCS_reco, cosThetaHX_reco, phiHX_reco, cosThetaLab_reco, phiLab_reco;

    double eventWeight, dimuonPtWeight, totalWeightCS, totalWeightHX, totalWeightLab;
	double dimuTrigWeight_nominal = -1,
	       dimuTrigWeight_systUp = -1, dimuTrigWeight_systDown = -1,
	       dimuTrigWeight_statUp = -1, dimuTrigWeight_statDown = -1;

	Float_t weightCS = 0, weightHX = 0;

	Bool_t allGood, firesTrigger, isRecoMatched, dimuonMatching, goodVertexProba, passHLTFilterMuons, trackerAndGlobalMuons, hybridSoftMuons;

	Int_t hiBin;

	Long64_t totEntries = OniaTree->GetEntries();
	
    double counter = 0;

	// Loop over the events
	for (Long64_t iEvent = 0; iEvent < (totEntries); iEvent++) {
		if (iEvent % 10000 == 0) {
			cout << Form("\rProcessing event %lld / %lld (%.0f%%)", iEvent, totEntries, 100. * iEvent / totEntries) << flush;
		}

		// if (iEvent > 100) break; // for testing purposes

        OniaTree->GetEntry(iEvent);

		// if (Centrality >= 2 * gCentralityBinMax) continue; // discard events with centrality >= 90% in 2018 data

		// firesTrigger = ((HLTriggers & (ULong64_t)(1 << (gUpsilonHLTBit - 1))) == (ULong64_t)(1 << (gUpsilonHLTBit - 1)));

		hiBin = GetHiBinFromhiHF(HFmean);

		eventWeight = Gen_weight * FindNcoll(hiBin); // * Get_zPV_weight(zVtx);

        // loop over all gen upsilons
		for (int iGen = 0; iGen < Gen_QQ_size; iGen++) {

            gen_QQ_LV = (TLorentzVector*)Gen_QQ_4mom->At(iGen);
   
            if (fabs(gen_QQ_LV->Rapidity()) < gRapidityMin || fabs(gen_QQ_LV->Rapidity()) > gRapidityMax) continue; // upsilon within fiducial region

            // if (gen_QQ_LV->Pt() > gPtMax) continue;

			// go to reco level (apply reco level weights to Hydjet GEN acceptance)
			Int_t iReco = Gen_QQ_whichRec[iGen];

			// if (Reco_QQ_sign[iReco] != 0) continue; // only opposite-sign muon pairs

            isRecoMatched = iReco > -1;
            
            if (isRecoMatched) {

                recoLorentzVector = (TLorentzVector*)CloneArr_QQ->At(iReco);
				double reco_QQ_pt = recoLorentzVector->Pt();

				dimuonPtWeight = Get_RecoPtWeight(recoLorentzVector->Rapidity(), reco_QQ_pt);

				// dimuonMatching = (Reco_QQ_trig[iReco] & (ULong64_t)(1 << (gUpsilonHLTBit - 1))) == (ULong64_t)(1 << (gUpsilonHLTBit - 1));

				// goodVertexProba = Reco_QQ_VtxProb[iReco] > 0.01;

				/// single-muon selection criteria
				int iMuPlus = Reco_QQ_mupl_idx[iReco];
				int iMuMinus = Reco_QQ_mumi_idx[iReco];

				// HLT filters
				bool mupl_L2Filter = ((Reco_mu_trig[iMuPlus] & ((ULong64_t)pow(2, gL2FilterBit))) == ((ULong64_t)pow(2, gL2FilterBit)));
				bool mupl_L3Filter = ((Reco_mu_trig[iMuPlus] & ((ULong64_t)pow(2, gL3FilterBit))) == ((ULong64_t)pow(2, gL3FilterBit)));
				bool mumi_L2Filter = ((Reco_mu_trig[iMuMinus] & ((ULong64_t)pow(2, gL2FilterBit))) == ((ULong64_t)pow(2, gL2FilterBit)));
				bool mumi_L3Filter = ((Reco_mu_trig[iMuMinus] & ((ULong64_t)pow(2, gL3FilterBit))) == ((ULong64_t)pow(2, gL3FilterBit)));

				// passHLTFilterMuons = (mupl_L2Filter && mumi_L3Filter) || (mupl_L3Filter && mumi_L2Filter) || (mupl_L3Filter && mumi_L3Filter);

				// // global AND tracker muons
				// trackerAndGlobalMuons = (Reco_mu_SelectionType[iMuPlus] & 2) && (Reco_mu_SelectionType[iMuPlus] & 8) && (Reco_mu_SelectionType[iMuMinus] & 2) && (Reco_mu_SelectionType[iMuMinus] & 8);

				// // passing hybrid-soft Id
				// hybridSoftMuons = (Reco_mu_nTrkWMea[iMuPlus] > 5) && (Reco_mu_nPixWMea[iMuPlus] > 0) && (fabs(Reco_mu_dxy[iMuPlus]) < 0.3) && (fabs(Reco_mu_dz[iMuPlus]) < 20.) && (Reco_mu_nTrkWMea[iMuMinus] > 5) && (Reco_mu_nPixWMea[iMuMinus] > 0) && (fabs(Reco_mu_dxy[iMuMinus]) < 0.3) && (fabs(Reco_mu_dz[iMuMinus]) < 20.);

				// /// numerator for the efficiency
				// allGood = firesTrigger && isRecoMatched && dimuonMatching && goodVertexProba && passHLTFilterMuons && trackerAndGlobalMuons && hybridSoftMuons;
				
                // get reco muon coordinates
				TLorentzVector* Reco_mupl_LV = (TLorentzVector*)CloneArr_mu->At(iMuPlus);
				double Reco_mupl_eta = Reco_mupl_LV->Eta();
				double Reco_mupl_pt = Reco_mupl_LV->Pt();

				TLorentzVector* Reco_mumi_LV = (TLorentzVector*)CloneArr_mu->At(iMuMinus);
				double Reco_mumi_eta = Reco_mumi_LV->Eta();
				double Reco_mumi_pt = Reco_mumi_LV->Pt();
                
				TVector3 Reco_mupl_vecLab = Reco_mupl_LV->Vect();
				TVector3 Reco_mupl_vecCS = MuPlusVector_CollinsSoper(*recoLorentzVector, *Reco_mupl_LV);
				TVector3 Reco_mupl_vecHX = MuPlusVector_Helicity(*recoLorentzVector, *Reco_mupl_LV);

				/// cosTheta and phi
				/// lab
				cosThetaLab_reco = Reco_mupl_vecLab.CosTheta();
				phiLab_reco = Reco_mupl_vecLab.Phi() * 180 / TMath::Pi();
				if (isPhiFolded == kTRUE) phiLab_reco = fabs(phiLab_reco);
				
				/// CS
				cosThetaCS_reco = Reco_mupl_vecCS.CosTheta();
				phiCS_reco = Reco_mupl_vecCS.Phi() * 180 / TMath::Pi();
				if (isPhiFolded == kTRUE) phiCS_reco = fabs(phiCS_reco);
				
				/// HX
				cosThetaHX_reco = Reco_mupl_vecHX.CosTheta();
				phiHX_reco = Reco_mupl_vecHX.Phi() * 180 / TMath::Pi();
				if (isPhiFolded == kTRUE) phiHX_reco = fabs(phiHX_reco);

				/// get gen muon coordinates
				/// LV
				gen_mupl_LV = (TLorentzVector*)Gen_mu_4mom->At(Gen_QQ_mupl_idx[iGen]);
				gen_mumi_LV = (TLorentzVector*)Gen_mu_4mom->At(Gen_QQ_mumi_idx[iGen]);

				/// vector
                TVector3 Gen_mupl_vecLab = gen_mupl_LV->Vect();
				TVector3 Gen_mupl_vecCS = MuPlusVector_CollinsSoper(*gen_QQ_LV, *gen_mupl_LV);
				TVector3 Gen_mupl_vecHX = MuPlusVector_Helicity(*gen_QQ_LV, *gen_mupl_LV);

				/// cosTheta and phi
				/// lab
                cosThetaLab_gen = Gen_mupl_vecLab.CosTheta();
                phiLab_gen = Gen_mupl_vecLab.Phi() * 180 / TMath::Pi();
                if (isPhiFolded == kTRUE) phiLab_gen = fabs(phiLab_gen);
				
				/// CS
                cosThetaCS_gen = Gen_mupl_vecCS.CosTheta();
                phiCS_gen = Gen_mupl_vecCS.Phi() * 180 / TMath::Pi();
				if (isPhiFolded == kTRUE) phiCS_gen = fabs(phiCS_gen);
                // cout << "cosThetaCS_gen: " << cosThetaCS_gen << endl;
                // cout << "phiCS_gen: " << phiCS_gen << endl;
				
				/// HX
                cosThetaHX_gen = Gen_mupl_vecHX.CosTheta();
                phiHX_gen = Gen_mupl_vecHX.Phi() * 180 / TMath::Pi();
				if (isPhiFolded == kTRUE) phiHX_gen = fabs(phiHX_gen);

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

				else {
					dimuTrigWeight_nominal = 1;
					dimuTrigWeight_systUp = 1;
					dimuTrigWeight_systDown = 1;
					dimuTrigWeight_statUp = 1;
					dimuTrigWeight_statDown = 1;
				}
				// dimuon efficiency weight = product of the total scale factors
				dimuWeight_nominal = tnp_weight_trk_pbpb(Reco_mupl_eta, indexNominal) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexNominal) * dimuTrigWeight_nominal;

				// single-muon acceptance
                withinAcceptance = MuonKinematicsWithinLimits(*gen_mupl_LV, muonAccName) && MuonKinematicsWithinLimits(*gen_mumi_LV, muonAccName);
    
                // totalWeightLab = eventWeight * dimuonPtWeight * dimuWeight_nominal;
				totalWeightLab = 1;
                accMatrixLab->FillWeighted(withinAcceptance, totalWeightLab, cosThetaLab_reco, phiLab_reco, reco_QQ_pt);

                // Reference frame transformations

                if (isPhiFolded == kTRUE)
                    weightCS = 1 + lambdaTheta * TMath::Power(Gen_mupl_vecCS.CosTheta(), 2) + lambdaPhi * TMath::Power(std::sin(Gen_mupl_vecCS.Theta()), 2) * std::cos(2 * fabs(Gen_mupl_vecCS.Phi())) + lambdaThetaPhi * std::sin(2 * Gen_mupl_vecCS.Theta()) * std::cos(fabs(Gen_mupl_vecCS.Phi()));
                else
                    weightCS = 1 + lambdaTheta * TMath::Power(Gen_mupl_vecCS.CosTheta(), 2) + lambdaPhi * TMath::Power(std::sin(Gen_mupl_vecCS.Theta()), 2) * std::cos(2 * Gen_mupl_vecCS.Phi()) + lambdaThetaPhi * std::sin(2 * Gen_mupl_vecCS.Theta()) * std::cos(Gen_mupl_vecCS.Phi());

                // total weight
				// totalWeightCS = eventWeight * dimuonPtWeight * dimuWeight_nominal * weightCS;
				totalWeightCS = 1;

                // cout << "weightCS: " << weightCS << endl;
                // cout << "cosThetaCS_gen: " << cosThetaCS_gen << endl;
                // cout << "phiCS_gen: " << phiCS_gen << endl;

                accMatrixCS->FillWeighted(withinAcceptance, totalWeightCS, cosThetaCS_reco, phiCS_reco, reco_QQ_pt);

                if (isPhiFolded == kTRUE)
                    weightHX = 1 + lambdaTheta * TMath::Power(Gen_mupl_vecHX.CosTheta(), 2) + lambdaPhi * TMath::Power(std::sin(Gen_mupl_vecHX.Theta()), 2) * std::cos(2 * fabs(Gen_mupl_vecHX.Phi())) + lambdaThetaPhi * std::sin(2 * Gen_mupl_vecHX.Theta()) * std::cos(fabs(Gen_mupl_vecHX.Phi()));
                else
                    weightHX = 1 + lambdaTheta * TMath::Power(Gen_mupl_vecHX.CosTheta(), 2) + lambdaPhi * TMath::Power(std::sin(Gen_mupl_vecHX.Theta()), 2) * std::cos(2 * Gen_mupl_vecHX.Phi()) + lambdaThetaPhi * std::sin(2 * Gen_mupl_vecHX.Theta()) * std::cos(Gen_mupl_vecHX.Phi());
                
                // totalWeightHX = eventWeight * dimuonPtWeight * dimuWeight_nominal * weightHX;
				totalWeightHX = 1;

                // cout << "weightHX: " << weightHX << endl;
                // cout << "cosThetaHX_gen: " << cosThetaHX_gen << endl;
                // cout << "phiHX_gen: " << phiHX_gen << endl;

                accMatrixHX->FillWeighted(withinAcceptance, totalWeightHX, cosThetaHX_reco, phiHX_reco, reco_QQ_pt);
				
				// cout << "Event: " << iEvent << endl;
				// cout << "withinAcceptance: " << withinAcceptance << endl;
				// cout << "weightHX: " << totalWeightHX << endl;
				// cout << "cosThetaHX: " << cosThetaHX_reco << endl;
				// cout << "phiHX: " << phiHX_reco << endl;
				// cout << "reco_QQ_pt: " << reco_QQ_pt << endl;

				// cout << "eventWeight: " << eventWeight << endl;
				// cout << "dimuonPtWeight: " << dimuonPtWeight << endl;
				// cout << "dimuWeight_nominal: " << dimuWeight_nominal << endl;
				// cout << "weightHX: " << weightHX << endl;

				// cout << "Reco_mupl_pt: " << Reco_mupl_pt << endl;
				// cout << "Reco_mupl_eta: " << Reco_mupl_eta << endl;
				// cout << "Reco_mumi_pt: " << Reco_mumi_pt << endl;
				// cout << "Reco_mumi_eta: " << Reco_mumi_eta << endl;

				// cout << "tnp_weight_trk_pbpb(Reco_mupl_eta, indexNominal): " << tnp_weight_trk_pbpb(Reco_mupl_eta, indexNominal) << endl;
				// cout << "tnp_weight_trk_pbpb(Reco_mumi_eta, indexNominal): " << tnp_weight_trk_pbpb(Reco_mumi_eta, indexNominal) << endl;
				// cout << "tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexNominal): " << tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexNominal) << endl;
				// cout << "tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexNominal): " << tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexNominal) << endl;
				// cout << "dimuTrigWeight_nominal: " << dimuTrigWeight_nominal << endl;

				// cout << "mupl_isL2: " << mupl_isL2 << endl;
				// cout << "mupl_isL3: " << mupl_isL3 << endl;
				// cout << "mumi_isL2: " << mumi_isL2 << endl;
				// cout << "mumi_isL3: " << mumi_isL3 << endl;
				// cout << "mupl_L2Filter: " << mupl_L2Filter << endl;
				// cout << "mupl_L3Filter: " << mupl_L3Filter << endl;
				// cout << "mumi_L2Filter: " << mumi_L2Filter << endl;
				// cout << "mumi_L3Filter: " << mumi_L3Filter << endl;

				// cout << "" << endl;

			}
        }
    }
	// Set the plot styles
	gStyle->SetPadLeftMargin(.15);
	//gStyle->SetTitleYOffset(.9);
	gStyle->SetPadRightMargin(0.18);
	SetColorPalette(gAcceptanceColorPaletteName);

 	/// save the results in a file for later usage
    const char* path = AcceptanceResultsPath(muonAccName.Data());

    gSystem->mkdir(path, kTRUE);
    // const char* outputFileName = Form("%s/AcceptanceResults_Hydjet%s.root", path, isPhiFolded ? "" : "_fullPhi");
	const char* outputFileName = Form("%s/AcceptanceResults_noWeights_Hydjet%s.root", path, isPhiFolded ? "" : "_fullPhi");
 
    TFile outputFile(outputFileName, "UPDATE");
 
    //accMatrixLab->Write();
    accMatrixCS->Write();
    accMatrixHX->Write();
     
    outputFile.Write();
    outputFile.Close();
 
    if (BeVerbose) cout << "\nAcceptance maps saved in " << outputFileName << endl;     
}

TEfficiency* getAcceptance3DMap(TString fileName = "", const char* refFrameName = "HX", Double_t lambdaTheta = 0, Double_t lambdaPhi = 0, Double_t lambdaThetaPhi = 0, Bool_t isPhiFolded = kFALSE) {
	/// get acceptance and efficiency in 1D
	TString nominalMapName = NominalTEfficiency3DName(refFrameName, lambdaTheta, lambdaPhi, lambdaThetaPhi);

	// TString fileName = Form("%s/AcceptanceResults_Hydjet%s.root", AcceptanceResultsPath(gMuonAccName), isPhiFolded ? "" : "_fullPhi");

	// get acceptance maps
	// TFile* acceptanceFile = openFile(fileName.Data());
	TFile* acceptanceFile = TFile::Open(fileName, "READ");

	if (!acceptanceFile) {
		std::cerr << "Error: acceptanceFile is null." << std::endl;
		acceptanceFile->Close();
		delete acceptanceFile;

		// acceptanceMap_noGenFilter(0, 30, gUpsilonState, lambdaTheta, lambdaPhi, lambdaThetaPhi, isPhiFolded, "MuonUpsilonTriggerAcc");

		return nullptr;
	}

	acceptanceFile = openFile(fileName.Data());

	auto* accMap = (TEfficiency*)acceptanceFile->Get(nominalMapName.Data());

	if (!accMap) {
		std::cerr << "Error: accMap is null." << std::endl;
		acceptanceFile->Close();
		delete acceptanceFile;
		delete accMap;

		// create the acceptance map
		cout << Form("Acceptance map not found. Creating a new one with (lambdaTheta, phi, thetaPhi) = (%.2f, %.2f, %.2f)....", lambdaTheta, lambdaPhi, lambdaThetaPhi) << endl;

		acceptanceMap_HydjetGen(isPhiFolded, "UpsilonTriggerThresholds", lambdaTheta, lambdaPhi, lambdaThetaPhi, gUpsilonState);

        // // Force TFile closure before reopening — if that function didn't already do it
        // gSystem->ProcessEvents();  // Sometimes helps
        // gSystem->Sleep(300);       // Give ROOT some time if on network drive

        // force ROOT to forget cached file
        gROOT->GetListOfFiles()->Remove(gROOT->GetFile(fileName.Data()));
        delete TFile::Open(fileName.Data());  // flush disk buffers once

		TFile* newAcceptanceFile = openFile(fileName.Data());
        // newAcceptanceFile->cd();
        // newAcceptanceFile->ls();  // now you should see your object

        cout << "nominal map name: " << nominalMapName << endl;
		auto* newAccMap = (TEfficiency*)newAcceptanceFile->Get(nominalMapName.Data());
		cout << "Acceptance map is loaded." << endl;
		if (!newAccMap) {
			std::cerr << "Error: accMap is still null after creating it." << std::endl;
			// accMap->ls();
			return nullptr;
		}

		return newAccMap;
	}

	return accMap;
}

void drawAcceptanceMap(TString refFrameName = "HX",
                      Float_t lambdaTheta = 0, Float_t lambdaPhi = 0, Float_t lambdaThetaPhi = 0,
                      Int_t ptMin = 0, Int_t ptMax = 30,
                      const Int_t nCosThetaBins = 5, Double_t cosThetaMin = -0.7, Double_t cosThetaMax = 0.7,
                      const Int_t nPhiBins = 6, Int_t phiMin = -180, Int_t phiMax = 180,
                      Bool_t isPhiFolded = kTRUE) {

	writeExtraText = true; // if extra text
	extraText = "       Internal";
	// extraText = "       Simulation Preliminary";

	/// set bin edges and width
	vector<Double_t> cosThetaBinEdges = setCosThetaBinEdges(nCosThetaBins, cosThetaMin, cosThetaMax);

	vector<Double_t> phiBinEdges = setPhiBinEdges(nPhiBins, phiMin, phiMax);

	/// bin width
	Double_t cosThetaStep = (cosThetaBinEdges[nCosThetaBins] - cosThetaBinEdges[0]) / nCosThetaBins;

	Double_t phiStep = (phiBinEdges[nPhiBins] - phiBinEdges[0]) / nPhiBins;

	/// get acceptance 3D map
	TEfficiency* accMap;
	TEfficiency* accMap_flat;

	TString fileName = Form("%s/AcceptanceResults_noWeights_Hydjet%s.root", AcceptanceResultsPath(gMuonAccName), isPhiFolded ? "" : "_fullPhi");
	// TString fileName = Form("%s/AcceptanceResults_Hydjet%s.root", AcceptanceResultsPath(gMuonAccName), isPhiFolded ? "" : "_fullPhi");
	TString fileName_flat = Form("%s/AcceptanceResults%s.root", AcceptanceResultsPath(gMuonAccName), isPhiFolded ? "" : "_fullPhi");

    accMap = getAcceptance3DMap(fileName, refFrameName.Data(), lambdaTheta, lambdaPhi, lambdaThetaPhi, isPhiFolded);
	accMap_flat = getAcceptance3DMap(fileName_flat, refFrameName.Data(), lambdaTheta, lambdaPhi, lambdaThetaPhi, isPhiFolded);
	
    TEfficiency* accMapCosThetaPhi = rebinTEff3DMap(accMap, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges);
	TEfficiency* accMapCosThetaPhi_flat = rebinTEff3DMap(accMap_flat, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges);
	
    TH2D* hTotalCosThetaPhi = (TH2D*)accMapCosThetaPhi->GetTotalHistogram();
    // hTotalCosThetaPhi->Scale(1. / hTotalCosThetaPhi->Integral(0, nCosThetaBins + 1, 0, nPhiBins + 1));
    
	TH2D* hPassedCosThetaPhi = (TH2D*)accMapCosThetaPhi->GetPassedHistogram();

	TH2D* hTotalCosThetaPhi_flat = (TH2D*)accMapCosThetaPhi_flat->GetTotalHistogram();
	TH2D* hPassedCosThetaPhi_flat = (TH2D*)accMapCosThetaPhi_flat->GetPassedHistogram();

	TH2D* hRatioCosThetaPhi = (TH2D*)hTotalCosThetaPhi->Clone("hRatioCosThetaPhi");
	hRatioCosThetaPhi->Divide(hPassedCosThetaPhi_flat, hPassedCosThetaPhi, 1, 1, "B");

    /// draw acceptance graph for check
	TCanvas* accCanvas = nullptr;
	TCanvas* accTotalCanvas = nullptr;
    TCanvas* accPassedCanvas = nullptr;
	TCanvas* accRatioCanvas = nullptr;

    accCanvas = DrawEfficiency2DHist(accMapCosThetaPhi, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, gUpsilonState, kTRUE, kTRUE, kFALSE, "_TriggerAcc", isPhiFolded, kTRUE, lambdaTheta, lambdaPhi, lambdaThetaPhi);
   
    accTotalCanvas = draw2DMap(hTotalCosThetaPhi, refFrameName.Data(), nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kFALSE, 1, isPhiFolded);
    display2DMapContents(hTotalCosThetaPhi, nCosThetaBins, nPhiBins, kFALSE);

    accPassedCanvas = draw2DMap(hPassedCosThetaPhi, refFrameName.Data(), nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kFALSE, 1, isPhiFolded);
    display2DMapContents(hPassedCosThetaPhi, nCosThetaBins, nPhiBins, kFALSE);

	accRatioCanvas = draw2DMap(hRatioCosThetaPhi, refFrameName.Data(), nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kFALSE, 1, isPhiFolded);
	display2DMapContents(hRatioCosThetaPhi, nCosThetaBins, nPhiBins, kFALSE);

	return;
}