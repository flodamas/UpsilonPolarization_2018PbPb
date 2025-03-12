#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/Parameters/PhaseSpace.h"

#include "../MonteCarlo/AccEffHelpers.h"

#include "../Polarization/PolarFitHelpers.h"

//#include "../Tools/Style/FitDistributions.h"
#include "../Tools/Style/Legends.h"

#include "../Tools/Parameters/CentralityValues.h"
#include "../Tools/Parameters/EfficiencyWeights.h"
#include "../Tools/Parameters/MuonScaleFactors.h"

#include "../ReferenceFrameTransform/Transformations.h"

void averageSF_MC(Int_t iState = gUpsilonState, Bool_t isPhiFolded = kFALSE, TString accName = "MuonUpsilonTriggerAcc") { // accName = "MuonSimpleAcc", "MuonWithin2018PbPbAcc", or "MuonUpsilonTriggerAcc")

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

    // Create a 3D histogram to store sums of SF and counts per bin
    TH3D* hSF_sum[2][4]; 
    hSF_sum[0][0] = new TH3D("hSF_sumCS_total", "Sum of SF in CS", NPtBins, gPtBinning, NCosThetaBins, gCosThetaBinning, NPhiBins, gPhiBinning);
    hSF_sum[1][0] = new TH3D("hSF_sumHX_total", "Sum of SF in HX", NPtBins, gPtBinning, NCosThetaBins, gCosThetaBinning, NPhiBins, gPhiBinning);
    hSF_sum[0][1] = new TH3D("hSF_sumCS_trk", "Sum of SF in CS", NPtBins, gPtBinning, NCosThetaBins, gCosThetaBinning, NPhiBins, gPhiBinning);
    hSF_sum[1][1] = new TH3D("hSF_sumHX_trk", "Sum of SF in HX", NPtBins, gPtBinning, NCosThetaBins, gCosThetaBinning, NPhiBins, gPhiBinning);
	hSF_sum[0][2] = new TH3D("hSF_sumCS_trig", "Sum of SF in CS", NPtBins, gPtBinning, NCosThetaBins, gCosThetaBinning, NPhiBins, gPhiBinning);
    hSF_sum[1][2] = new TH3D("hSF_sumHX_trig", "Sum of SF in HX", NPtBins, gPtBinning, NCosThetaBins, gCosThetaBinning, NPhiBins, gPhiBinning);
	hSF_sum[0][3] = new TH3D("hSF_sumCS_muID", "Sum of SF in CS", NPtBins, gPtBinning, NCosThetaBins, gCosThetaBinning, NPhiBins, gPhiBinning);
    hSF_sum[1][3] = new TH3D("hSF_sumHX_muID", "Sum of SF in HX", NPtBins, gPtBinning, NCosThetaBins, gCosThetaBinning, NPhiBins, gPhiBinning);      

    TH3D* hSF_count[2][4]; 
    hSF_count[0][0] = new TH3D("hSF_countCS_total", "Count of events", NPtBins, gPtBinning, NCosThetaBins, gCosThetaBinning, NPhiBins, gPhiBinning);
    hSF_count[1][0] = new TH3D("hSF_countHX_total", "Count of events", NPtBins, gPtBinning, NCosThetaBins, gCosThetaBinning, NPhiBins, gPhiBinning);
    hSF_count[0][1] = new TH3D("hSF_countCS_trk", "Count of events", NPtBins, gPtBinning, NCosThetaBins, gCosThetaBinning, NPhiBins, gPhiBinning);
    hSF_count[1][1] = new TH3D("hSF_countHX_trk", "Count of events", NPtBins, gPtBinning, NCosThetaBins, gCosThetaBinning, NPhiBins, gPhiBinning);
    hSF_count[0][2] = new TH3D("hSF_countCS_trig", "Count of events", NPtBins, gPtBinning, NCosThetaBins, gCosThetaBinning, NPhiBins, gPhiBinning);
    hSF_count[1][2] = new TH3D("hSF_countHX_trig", "Count of events", NPtBins, gPtBinning, NCosThetaBins, gCosThetaBinning, NPhiBins, gPhiBinning);
    hSF_count[0][3] = new TH3D("hSF_countCS_muID", "Count of events", NPtBins, gPtBinning, NCosThetaBins, gCosThetaBinning, NPhiBins, gPhiBinning);
    hSF_count[1][3] = new TH3D("hSF_countHX_muID", "Count of events", NPtBins, gPtBinning, NCosThetaBins, gCosThetaBinning, NPhiBins, gPhiBinning);

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
	double dimuWeight_total_systUp, dimuWeight_total_systDown;

	// loop variables
	TLorentzVector* genLorentzVector = new TLorentzVector();
	TLorentzVector* recoLorentzVector = new TLorentzVector();

	double eventWeight, dimuonPtWeight, totalWeightCS, totalWeightHX;
	double dimuTrigWeight_nominal = -1, dimuTrigWeight_systUp = -1, dimuTrigWeight_systDown = -1, dimuTrigWeight_statUp = -1, dimuTrigWeight_statDown = -1;
	double dimuTrkWeight_nominal = -1, dimuMuIdWeight_nominal = -1;

	Float_t weightCS = 0, weightHX = 0;

	TString MuonAccName = "";

	Bool_t allGood, firesTrigger, isRecoMatched, dimuonMatching, goodVertexProba, passHLTFilterMuons, trackerAndGlobalMuons, hybridSoftMuons;

	Long64_t totEntries = OniaTree->GetEntries();

	for (Long64_t iEvent = 0; iEvent < totEntries; iEvent++) {
		if (iEvent % 10000 == 0) {
			cout << Form("\rProcessing event %lld / %lld (%.0f%%)", iEvent, totEntries, 100. * iEvent / totEntries) << flush;
		}

		OniaTree->GetEntry(iEvent);

		if (Centrality >= 2 * gCentralityBinMax) continue; // discard events with centrality >= 90% in 2018 data

		firesTrigger = ((HLTriggers & (ULong64_t)(1 << (gUpsilonHLTBit - 1))) == (ULong64_t)(1 << (gUpsilonHLTBit - 1)));
		
        eventWeight = Gen_weight * FindNcoll(Centrality); // * Get_zPV_weight(zVtx);

        // loop over all gen upsilons
        for (int iGen = 0; iGen < Gen_QQ_size; iGen++) {
            genLorentzVector = (TLorentzVector*)Gen_QQ_4mom->At(iGen);

            // fiducial region
            //if (genLorentzVector->Pt() < ptMin || genLorentzVector->Pt() > ptMax) continue; // pt bin of interest

            if (fabs(genLorentzVector->Rapidity()) < gRapidityMin || fabs(genLorentzVector->Rapidity()) > gRapidityMax) continue;

            // single-muon acceptance

            // positive muon first
            genLorentzVector = (TLorentzVector*)Gen_mu_4mom->At(Gen_QQ_mupl_idx[iGen]);

            if (accName == TString("MuonUpsilonTriggerAcc")) {if (!MuonUpsilonTriggerAcc(*genLorentzVector)) continue; MuonAccName = "_TriggerAcc";}
            else if (accName == TString("MuonSimpleAcc")) {if (!MuonSimpleAcc(*genLorentzVector)) continue; MuonAccName = "_SimpleAcc";}
            else if (accName == TString("MuonWithin2018PbPbAcc")) {if (!MuonWithin2018PbPbAcc(*genLorentzVector)) continue; MuonAccName = "_2018Acc";}
            else {
                cout << "Invalid acceptance name. Please choose from 'MuonUpsilonTriggerAcc', 'MuonWithin2018PbPbAcc', or 'MuonSimpleAcc'." << endl;
                return;
            }

            // then negative muon
            genLorentzVector = (TLorentzVector*)Gen_mu_4mom->At(Gen_QQ_mumi_idx[iGen]);

            if (accName == TString("MuonUpsilonTriggerAcc")) {if (!MuonUpsilonTriggerAcc(*genLorentzVector)) continue;}
            else if (accName == TString("MuonSimpleAcc")) {if (!MuonSimpleAcc(*genLorentzVector)) continue;}
            else if (accName == TString("MuonWithin2018PbPbAcc")) {if (!MuonWithin2018PbPbAcc(*genLorentzVector)) continue;}

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

				// get muon coordinates

				TLorentzVector* Reco_mupl_LV = (TLorentzVector*)CloneArr_mu->At(iMuPlus);
				double Reco_mupl_eta = Reco_mupl_LV->Eta();
				double Reco_mupl_pt = Reco_mupl_LV->Pt();

				TLorentzVector* Reco_mumi_LV = (TLorentzVector*)CloneArr_mu->At(iMuMinus);
				double Reco_mumi_eta = Reco_mumi_LV->Eta();
				double Reco_mumi_pt = Reco_mumi_LV->Pt();

				TVector3 muPlus_CS = MuPlusVector_CollinsSoper(*recoLorentzVector, *Reco_mupl_LV);
				double cosThetaCS = muPlus_CS.CosTheta();
				double phiCS = 0;

				if (isPhiFolded == kTRUE) phiCS = fabs(muPlus_CS.Phi() * 180 / TMath::Pi());
				else phiCS = muPlus_CS.Phi() * 180 / TMath::Pi();

				TVector3 muPlus_HX = MuPlusVector_Helicity(*recoLorentzVector, *Reco_mupl_LV);
				double cosThetaHX = muPlus_HX.CosTheta();
				double phiHX = 0;
                
                if (isPhiFolded == kTRUE) phiHX = fabs(muPlus_HX.Phi() * 180 / TMath::Pi());
				else phiHX = muPlus_HX.Phi() * 180 / TMath::Pi();
 
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
                
				dimuTrkWeight_nominal = tnp_weight_trk_pbpb(Reco_mupl_eta, indexNominal) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexNominal);
				dimuMuIdWeight_nominal = tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexNominal);

                int binX = hSF_sum[0][0]->GetXaxis()->FindBin(reco_QQ_pt);
                int binY_CS = hSF_sum[0][0]->GetYaxis()->FindBin(cosThetaCS);
                int binZ_CS = hSF_sum[0][0]->GetZaxis()->FindBin(phiCS);
            
                int binY_HX = hSF_sum[1][0]->GetYaxis()->FindBin(cosThetaHX);
                int binZ_HX = hSF_sum[1][0]->GetZaxis()->FindBin(phiHX);
				
				/// fill total SF histograms
                hSF_sum[0][0]->SetBinContent(binX, binY_CS, binZ_CS, hSF_sum[0][0]->GetBinContent(binX, binY_CS, binZ_CS) + dimuWeight_nominal);
                hSF_count[0][0]->SetBinContent(binX, binY_CS, binZ_CS, (hSF_count[0][0]->GetBinContent(binX, binY_CS, binZ_CS) + 1));

                hSF_sum[1][0]->SetBinContent(binX, binY_HX, binZ_HX, hSF_sum[1][0]->GetBinContent(binX, binY_HX, binZ_HX) + dimuWeight_nominal);
                hSF_count[1][0]->SetBinContent(binX, binY_HX, binZ_HX, (hSF_count[1][0]->GetBinContent(binX, binY_HX, binZ_HX) + 1));
            
				/// fill tracker SF histograms
				hSF_sum[0][1]->SetBinContent(binX, binY_CS, binZ_CS, hSF_sum[0][1]->GetBinContent(binX, binY_CS, binZ_CS) + dimuTrkWeight_nominal);
				hSF_count[0][1]->SetBinContent(binX, binY_CS, binZ_CS, (hSF_count[0][1]->GetBinContent(binX, binY_CS, binZ_CS) + 1));

				hSF_sum[1][1]->SetBinContent(binX, binY_HX, binZ_HX, hSF_sum[1][1]->GetBinContent(binX, binY_HX, binZ_HX) + dimuTrkWeight_nominal);
				hSF_count[1][1]->SetBinContent(binX, binY_HX, binZ_HX, (hSF_count[1][1]->GetBinContent(binX, binY_HX, binZ_HX) + 1));

				/// fill trigger SF histograms
				hSF_sum[0][2]->SetBinContent(binX, binY_CS, binZ_CS, hSF_sum[0][2]->GetBinContent(binX, binY_CS, binZ_CS) + dimuTrigWeight_nominal);
				hSF_count[0][2]->SetBinContent(binX, binY_CS, binZ_CS, (hSF_count[0][2]->GetBinContent(binX, binY_CS, binZ_CS) + 1));

				hSF_sum[1][2]->SetBinContent(binX, binY_HX, binZ_HX, hSF_sum[1][2]->GetBinContent(binX, binY_HX, binZ_HX) + dimuTrigWeight_nominal);
				hSF_count[1][2]->SetBinContent(binX, binY_HX, binZ_HX, (hSF_count[1][2]->GetBinContent(binX, binY_HX, binZ_HX) + 1));

				/// fill muon ID SF histograms
				hSF_sum[0][3]->SetBinContent(binX, binY_CS, binZ_CS, hSF_sum[0][3]->GetBinContent(binX, binY_CS, binZ_CS) + dimuMuIdWeight_nominal);
				hSF_count[0][3]->SetBinContent(binX, binY_CS, binZ_CS, (hSF_count[0][3]->GetBinContent(binX, binY_CS, binZ_CS) + 1));

				hSF_sum[1][3]->SetBinContent(binX, binY_HX, binZ_HX, hSF_sum[1][3]->GetBinContent(binX, binY_HX, binZ_HX) + dimuMuIdWeight_nominal);
				hSF_count[1][3]->SetBinContent(binX, binY_HX, binZ_HX, (hSF_count[1][3]->GetBinContent(binX, binY_HX, binZ_HX) + 1));			
			}
        }
    }

    /// save the histograms
    TString outputFileName = Form("averageSF_MC_%s.root", accName.Data());
    TFile* outputFile = new TFile(outputFileName, "RECREATE");
    
	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < 4; ++j) {
			hSF_sum[i][j]->Write();
			hSF_count[i][j]->Write();
		}
	}

    outputFile->Close();

    cout << "Histograms saved in " << outputFileName << endl;

    file->Close();

    // delete genLorentzVector;
    // delete recoLorentzVector;

    return;
}

void drawAverageSF_MC(TString accName = "MuonUpsilonTriggerAcc") { // accName = "MuonSimpleAcc", "MuonWithin2018PbPbAcc", or "MuonUpsilonTriggerAcc"
    TFile* file = TFile::Open(Form("averageSF_MC_%s.root", accName.Data()), "READ");
    if (!file) {
        cout << "File not found. Check the directory of the file." << endl;
        return;
    }

    TH3D* hSF_sum[2][4]; 
	hSF_sum[0][0] = (TH3D*)file->Get("hSF_sumCS_total");
	hSF_sum[1][0] = (TH3D*)file->Get("hSF_sumHX_total");

	hSF_sum[0][1] = (TH3D*)file->Get("hSF_sumCS_trk");
	hSF_sum[1][1] = (TH3D*)file->Get("hSF_sumHX_trk");

	hSF_sum[0][2] = (TH3D*)file->Get("hSF_sumCS_trig");
	hSF_sum[1][2] = (TH3D*)file->Get("hSF_sumHX_trig");

	hSF_sum[0][3] = (TH3D*)file->Get("hSF_sumCS_muID");
	hSF_sum[1][3] = (TH3D*)file->Get("hSF_sumHX_muID");

    TH3D* hSF_count[2][4]; 
	hSF_count[0][0] = (TH3D*)file->Get("hSF_countCS_total");
	hSF_count[1][0] = (TH3D*)file->Get("hSF_countHX_total");

	hSF_count[0][1] = (TH3D*)file->Get("hSF_countCS_trk");
	hSF_count[1][1] = (TH3D*)file->Get("hSF_countHX_trk");

	hSF_count[0][2] = (TH3D*)file->Get("hSF_countCS_trig");
	hSF_count[1][2] = (TH3D*)file->Get("hSF_countHX_trig");

	hSF_count[0][3] = (TH3D*)file->Get("hSF_countCS_muID");
	hSF_count[1][3] = (TH3D*)file->Get("hSF_countHX_muID");

    TH3D* hSF_average[2][4];
	hSF_average[0][0] = (TH3D*)hSF_sum[0][0]->Clone("hSF_averageCS_total");
	hSF_average[1][0] = (TH3D*)hSF_sum[1][0]->Clone("hSF_averageHX_total");

	hSF_average[0][1] = (TH3D*)hSF_sum[0][1]->Clone("hSF_averageCS_trk");
	hSF_average[1][1] = (TH3D*)hSF_sum[1][1]->Clone("hSF_averageHX_trk");

	hSF_average[0][2] = (TH3D*)hSF_sum[0][2]->Clone("hSF_averageCS_trig");
	hSF_average[1][2] = (TH3D*)hSF_sum[1][2]->Clone("hSF_averageHX_trig");

	hSF_average[0][3] = (TH3D*)hSF_sum[0][3]->Clone("hSF_averageCS_muID");
	hSF_average[1][3] = (TH3D*)hSF_sum[1][3]->Clone("hSF_averageHX_muID");

	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < 4; ++j) {
			hSF_average[i][j]->Divide(hSF_count[i][j]); // divide sum by count to get average   
		}
	}

    // Create a new TH2D to store the result for a specific x-bin
    TH2D* hSF_Avg2D[NPtBins][2][4];
    TCanvas* canvas[NPtBins][2][4];
    
    // Loop over the y and z bins and sum the content for the selected x-bin
	for (int refFrameBin = 0; refFrameBin < 2; refFrameBin++) {
		const char* refFrameName = (refFrameBin == 0) ? "CS" : "HX";

		for (int ptBin = 0; ptBin < NPtBins; ptBin++) {

			for (int SFType = 0; SFType < 4; SFType++) {
				const char* SFTypeName = (SFType == 0) ? "total" : (SFType == 1) ? "trk" : (SFType == 2) ? "trig" : "muID";

				// Create a new TH2D to store the result for a specific x-bin
				hSF_Avg2D[ptBin][refFrameBin][SFType] = new TH2D(Form("hSF_Avg2D_pt%dto%d_%s_%s", (int)gPtBinning[ptBin], (int)gPtBinning[ptBin +1], refFrameName, SFTypeName), "", 
				hSF_average[refFrameBin][SFType]->GetNbinsY(), hSF_average[refFrameBin][SFType]->GetYaxis()->GetXmin(), hSF_average[refFrameBin][SFType]->GetYaxis()->GetXmax(),
				hSF_average[refFrameBin][SFType]->GetNbinsZ(), hSF_average[refFrameBin][SFType]->GetZaxis()->GetXmin(), hSF_average[refFrameBin][SFType]->GetZaxis()->GetXmax());
			
				for (int cosThetaBin = 1; cosThetaBin <= hSF_average[refFrameBin][SFType]->GetNbinsY(); cosThetaBin++) {
					for (int phiBin = 1; phiBin <= hSF_average[refFrameBin][SFType]->GetNbinsZ(); phiBin++) {
						
						// Get the content of the selected x-bin, y-bin, and z-bin
						double content = hSF_average[refFrameBin][SFType]->GetBinContent(ptBin + 1, cosThetaBin, phiBin);

						// Set the content in the new 2D histogram for the selected x-bin
						hSF_Avg2D[ptBin][refFrameBin][SFType]->SetBinContent(cosThetaBin, phiBin, content);
					}
				}
				canvas[ptBin][refFrameBin][SFType] = new TCanvas(Form("canvas_pt%dto%d_%s_%s", (int)gPtBinning[ptBin], (int)gPtBinning[ptBin + 1], refFrameName, SFTypeName), Form("canvas_pt%dto%d_%s", (int)gPtBinning[ptBin], (int)gPtBinning[ptBin + 1], refFrameName), 650, 600);
				canvas[ptBin][refFrameBin][SFType]->SetRightMargin(0.18);
				canvas[ptBin][refFrameBin][SFType]->SetLeftMargin(0.15);
				canvas[ptBin][refFrameBin][SFType]->SetBottomMargin(0.15);

				hSF_Avg2D[ptBin][refFrameBin][SFType]->GetXaxis()->SetTitle(CosThetaVarTitle(refFrameName));
				hSF_Avg2D[ptBin][refFrameBin][SFType]->GetYaxis()->SetTitle(AbsPhiAxisTitle(refFrameName));
				hSF_Avg2D[ptBin][refFrameBin][SFType]->GetZaxis()->SetTitle(Form("average %s SF", SFTypeName));
				
				hSF_Avg2D[ptBin][refFrameBin][SFType]->GetYaxis()->SetTitleOffset(1.1);
				hSF_Avg2D[ptBin][refFrameBin][SFType]->GetZaxis()->SetTitleOffset(1.0);

				hSF_Avg2D[ptBin][refFrameBin][SFType]->GetXaxis()->CenterTitle();
				hSF_Avg2D[ptBin][refFrameBin][SFType]->GetYaxis()->CenterTitle();
				
				hSF_Avg2D[ptBin][refFrameBin][SFType]->GetYaxis()->SetNdivisions(4);
				hSF_Avg2D[ptBin][refFrameBin][SFType]->GetYaxis()->SetRangeUser(0, 240);

				hSF_Avg2D[ptBin][refFrameBin][SFType]->GetZaxis()->SetMaxDigits(1);
				hSF_Avg2D[ptBin][refFrameBin][SFType]->GetZaxis()->SetRangeUser(0.8, 1.2);
				hSF_Avg2D[ptBin][refFrameBin][SFType]->GetZaxis()->SetNdivisions(4);
				hSF_Avg2D[ptBin][refFrameBin][SFType]->Draw("colz");    

				display2DMapContents(hSF_Avg2D[ptBin][refFrameBin][SFType], NCosThetaBins, NPhiBins, kFALSE);

				/// draw legend
				TLatex legend;
				legend.SetTextAlign(22);
				legend.SetTextSize(0.04);
				legend.DrawLatexNDC(.5, .87, Form("%s < 2.4, %s", gDimuonRapidityVarTitle, DimuonPtRangeText(gPtBinning[ptBin], gPtBinning[ptBin + 1])));

				if (strcmp(accName, "MuonUpsilonTriggerAcc") == 0) legend.DrawLatexNDC(.49, .81, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, %s", gUpsilonState, gMuonPtCutText));
				else if (strcmp(accName, "MuonSimpleAcc") == 0) legend.DrawLatexNDC(.49, .81, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, #it{p}_{T}^{ #mu} > 3.5 GeV/#it{c}", gUpsilonState));
				else if (strcmp(accName, "MuonWithin2018PbPbAcc") == 0) legend.DrawLatexNDC(.49, .81, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, #it{p}_{T}^{ #mu} > 2018PbPbAcc", gUpsilonState));

				CMS_lumi(canvas[ptBin][refFrameBin][SFType], gCMSLumiText);

				gSystem->mkdir("averageSF", kTRUE);
				canvas[ptBin][refFrameBin][SFType]->SaveAs(Form("averageSF/averageTotalSF_MC_pt%dto%d_%s_%s.png", (int)gPtBinning[ptBin], (int)gPtBinning[ptBin + 1], refFrameName, SFTypeName));
			}
		}
	}

    return;
}