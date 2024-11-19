#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/Parameters/PhaseSpace.h"

#include "../Tools/Parameters/CentralityValues.h"
#include "../Tools/Parameters/EfficiencyWeights.h"
#include "../Tools/Parameters/MuonScaleFactors.h"

#include "../ReferenceFrameTransform/Transformations.h"

// skim and keep all MC events, but add a category "selected" to know if a given event passes the selection criteria or no

// the reconstructed dimuons are fully weighted, including potential reweighting for polarization

void skimReconstructedMCWeighted(Int_t iState = 1, Double_t lambdaTheta = 0, Double_t lambdaPhi = 0, Double_t lambdaThetaPhi = 0) {
	const char* inputFileName = Form("OniaTree_Y%dS_pThat2_HydjetDrumMB_miniAOD.root", iState);

	const char* outputFileName = Form("Y%dSReconstructedMCWeightedDataset%s_Lambda_Theta%.2f_Phi%.2f_ThetaPhi%.2f.root", iState, gMuonAccName, lambdaTheta, lambdaPhi, lambdaThetaPhi);

	TFile* infile = TFile::Open(inputFileName, "READ");
	TTree* OniaTree = (TTree*)infile->Get("hionia/myTree");

	/// OniaTree variables

	Float_t Gen_weight;

	TClonesArray* Gen_QQ_4mom = nullptr;
	TClonesArray* Gen_mu_4mom = nullptr;
	Short_t Gen_QQ_size;
	Short_t Gen_QQ_whichRec[1000];

	Short_t Gen_QQ_mupl_idx[1000];
	Short_t Gen_QQ_mumi_idx[1000];

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

	Int_t Reco_mu_SelectionType[1000];
	//(parameters for quality cuts)
	Float_t Reco_QQ_VtxProb[1000];
	Int_t Reco_mu_nPixWMea[1000];
	Int_t Reco_mu_nTrkWMea[1000];
	Float_t Reco_mu_dxy[1000];
	Float_t Reco_mu_dz[1000];
	ULong64_t Reco_mu_trig[1000];

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

	OniaTree->SetBranchAddress("Reco_QQ_4mom", &CloneArr_QQ);
	OniaTree->SetBranchAddress("Reco_mu_4mom", &CloneArr_mu);
	OniaTree->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size);
	OniaTree->SetBranchAddress("Reco_QQ_sign", &Reco_QQ_sign);
	OniaTree->SetBranchAddress("Reco_QQ_mupl_idx", Reco_QQ_mupl_idx);
	OniaTree->SetBranchAddress("Reco_QQ_mumi_idx", Reco_QQ_mumi_idx);
	OniaTree->SetBranchAddress("Reco_QQ_whichGen", Reco_QQ_whichGen);

	OniaTree->SetBranchAddress("Reco_mu_trig", Reco_mu_trig);
	OniaTree->SetBranchAddress("Reco_mu_SelectionType", Reco_mu_SelectionType);

	OniaTree->SetBranchAddress("Reco_QQ_VtxProb", &Reco_QQ_VtxProb);
	OniaTree->SetBranchAddress("Reco_mu_nPixWMea", &Reco_mu_nPixWMea);
	OniaTree->SetBranchAddress("Reco_mu_nTrkWMea", &Reco_mu_nTrkWMea);
	OniaTree->SetBranchAddress("Reco_mu_dxy", &Reco_mu_dxy);
	OniaTree->SetBranchAddress("Reco_mu_dz", &Reco_mu_dz);

	/// RooDataSet output: one entry = one dimuon candidate!
	TFile file(outputFileName, "RECREATE");

	// weighting by event directly on the fly
	RooRealVar centVar("centrality", "event centrality", 0, 200);
	RooRealVar eventWeightCSVar("eventWeightCS", "event-by-event weight (Ncoll x MC gen weight x reco pT reweighting x muon scale factors x polarization in CS)", 0, 100000);
	RooRealVar eventWeightHXVar("eventWeightHX", "event-by-event weight (Ncoll x MC gen weight x reco pT reweighting x muon scale factors x polarization in HX)", 0, 100000);

	Float_t lowMassCut = 8, highMassCut = 11;
	RooRealVar massVar("mass", gMassVarTitle, lowMassCut, highMassCut, gMassUnit);
	RooRealVar yVar("rapidity", gDimuonRapidityVarTitle, 0, 2.4);

	RooRealVar ptVar("pt", gDimuonPtVarTitle, 0, 100, gPtUnit);

	char* refFrameName = "Lab";
	RooRealVar cosThetaLabVar(CosThetaVarName(refFrameName), CosThetaVarTitle(refFrameName), -1, 1);
	RooRealVar phiLabVar(PhiVarName(refFrameName), PhiVarTitle(refFrameName), -180, 180, gPhiUnit);
	RooRealVar etaLabMuplVar("etaLabMupl", "eta of positive muon in the lab frame", -2.4, 2.4);
	RooRealVar etaLabMumiVar("etaLabMumi", "eta of negative muon in the lab frame", -2.4, 2.4);

	refFrameName = "CS";
	RooRealVar cosThetaCSVar(CosThetaVarName(refFrameName), CosThetaVarTitle(refFrameName), -1, 1);
	RooRealVar phiCSVar(PhiVarName(refFrameName), PhiVarTitle(refFrameName), -180, 180, gPhiUnit);
	RooRealVar phiTildeCSVar(PhiTildeVarName(refFrameName), PhiVarTitle(refFrameName), -180, 180, gPhiUnit);

	refFrameName = "HX";
	RooRealVar cosThetaHXVar(CosThetaVarName(refFrameName), CosThetaVarTitle(refFrameName), -1, 1);
	RooRealVar phiHXVar(PhiVarName(refFrameName), PhiVarTitle(refFrameName), -180, 180, gPhiUnit);
	RooRealVar phiTildeHXVar(PhiTildeVarName(refFrameName), PhiTildeVarTitle(refFrameName), -180, 180, gPhiUnit);

	RooRealVar lambdaThetaVar("lambdaTheta", "", -1, 1);
	RooRealVar lambdaPhiVar("lambdaPhi", "", -1, 1);
	RooRealVar lambdaThetaPhiVar("lambdaThetaPhi", "", -1, 1);

	// fix polarization extraction parameters
	lambdaThetaVar.setVal(lambdaTheta);
	lambdaThetaVar.setConstant(kTRUE);

	lambdaPhiVar.setVal(lambdaPhi);
	lambdaPhiVar.setConstant(kTRUE);

	lambdaThetaPhiVar.setVal(lambdaThetaPhi);
	lambdaThetaPhiVar.setConstant(kTRUE);

	RooCategory recoCat("recoCategory", "is the dimuon selected?"); // for efficiency PDF
	recoCat.defineType("selected", 1);
	recoCat.defineType("rejected", 0);

	RooDataSet datasetCS("MCdatasetCS", "skimmed MC dataset in CS", RooArgSet(recoCat, centVar, eventWeightCSVar, massVar, yVar, ptVar, cosThetaLabVar, phiLabVar, etaLabMuplVar, etaLabMumiVar, cosThetaCSVar, phiCSVar, phiTildeCSVar, cosThetaHXVar, phiHXVar, phiTildeHXVar), RooFit::WeightVar("eventWeightCS"), RooFit::StoreAsymError(RooArgSet(eventWeightCSVar)));
	RooDataSet datasetHX("MCdatasetHX", "skimmed MC dataset in HX", RooArgSet(recoCat, centVar, eventWeightHXVar, massVar, yVar, ptVar, cosThetaLabVar, phiLabVar, etaLabMuplVar, etaLabMumiVar, cosThetaCSVar, phiCSVar, phiTildeCSVar, cosThetaHXVar, phiHXVar, phiTildeHXVar), RooFit::WeightVar("eventWeightHX"), RooFit::StoreAsymError(RooArgSet(eventWeightHXVar)));

	// loop variables
	TLorentzVector* genLorentzVector = new TLorentzVector();

	Float_t nColl, weight = 0, totalWeightCS = 0, totalWeightHX = 0, polarWeightCS = 0, polarWeightHX = 0, dimuonPtWeight = 0, errorWeightDown = 0, errorWeightUp = 0;

	// for muon scale factors
	int indexNominal = 0;
	int indexSystUp = -1, indexSystDown = -2;
	int indexStatUp = +1, indexStatDown = +2;

	double dimuWeight_nominal = 0;
	double dimuWeight_trk_systUp, dimuWeight_trk_systDown, dimuWeight_trk_statUp, dimuWeight_trk_statDown;
	double dimuWeight_muId_systUp, dimuWeight_muId_systDown, dimuWeight_muId_statUp, dimuWeight_muId_statDown;
	double dimuWeight_trig_systUp, dimuWeight_trig_systDown, dimuWeight_trig_statUp, dimuWeight_trig_statDown;

	double dimuTrigWeight_nominal = -1, dimuTrigWeight_systUp = -1, dimuTrigWeight_systDown = -1, dimuTrigWeight_statUp = -1, dimuTrigWeight_statDown = -1;

	Long64_t totEntries = OniaTree->GetEntries();

	// booleans for the categorisation
	Bool_t allGood, firesTrigger, isRecoMatched, dimuonMatching, goodVertexProba, passHLTFilterMuons, trackerAndGlobalMuons, hybridSoftMuons;

	for (Long64_t iEvent = 0; iEvent < (totEntries); iEvent++) {
		if (iEvent % 10000 == 0) {
			cout << Form("\rProcessing event %lld / %lld (%.0f%%)", iEvent, totEntries, 100. * iEvent / totEntries) << flush;
		}

		OniaTree->GetEntry(iEvent);

		// event selection

		if (Centrality >= 2 * gCentralityBinMax) continue; // discard events with centrality > 90% in 2018 data

		// firesTrigger = ((HLTriggers & (ULong64_t)(1 << (gUpsilonHLTBit - 1))) == (ULong64_t)(1 << (gUpsilonHLTBit - 1)));

		if (!((HLTriggers & (ULong64_t)(1 << (gUpsilonHLTBit - 1))) == (ULong64_t)(1 << (gUpsilonHLTBit - 1)))) continue; // must fire the upsilon HLT path

		nColl = FindNcoll(Centrality);

		// loop over all gen upsilons
		for (int iGen = 0; iGen < Gen_QQ_size; iGen++) {
			genLorentzVector = (TLorentzVector*)Gen_QQ_4mom->At(iGen);

			// fiducial region

			if (fabs(genLorentzVector->Rapidity()) < gRapidityMin || fabs(genLorentzVector->Rapidity()) > gRapidityMax) continue;

			// single-muon acceptance

			// positive muon first
			genLorentzVector = (TLorentzVector*)Gen_mu_4mom->At(Gen_QQ_mupl_idx[iGen]);

			if (!MuonUpsilonTriggerAcc(*genLorentzVector)) continue;

			// then negative muon
			genLorentzVector = (TLorentzVector*)Gen_mu_4mom->At(Gen_QQ_mumi_idx[iGen]);

			if (!MuonUpsilonTriggerAcc(*genLorentzVector)) continue;

			// go to reco level
			Int_t iReco = Gen_QQ_whichRec[iGen];

			if (Reco_QQ_sign[iReco] != 0) continue; // only opposite-sign muon pairs

			isRecoMatched = iReco > -1;

			if (!isRecoMatched) continue;

			/// all the reconstructed upsilons must pass the conditions below!

			// dimuonMatching = (Reco_QQ_trig[iReco] & (ULong64_t)(1 << (gUpsilonHLTBit - 1))) == (ULong64_t)(1 << (gUpsilonHLTBit - 1));

			if (!((Reco_QQ_trig[iReco] & (ULong64_t)(1 << (gUpsilonHLTBit - 1))) == (ULong64_t)(1 << (gUpsilonHLTBit - 1)))) continue; // dimuon matching

			// goodVertexProba = Reco_QQ_VtxProb[iReco] > 0.01;
			if (Reco_QQ_VtxProb[iReco] < 0.01) continue;
			
			TLorentzVector* Reco_QQ_4mom = (TLorentzVector*)CloneArr_QQ->At(iReco);

			/// single-muon selection criteria
			int iMuPlus = Reco_QQ_mupl_idx[iReco];
			int iMuMinus = Reco_QQ_mumi_idx[iReco];

			// global AND tracker muons
			// trackerAndGlobalMuons = (Reco_mu_SelectionType[iMuPlus] & 2) && (Reco_mu_SelectionType[iMuPlus] & 8) && (Reco_mu_SelectionType[iMuMinus] & 2) && (Reco_mu_SelectionType[iMuMinus] & 8);
			// global AND tracker muons
			if (!((Reco_mu_SelectionType[iMuPlus] & 2) && (Reco_mu_SelectionType[iMuPlus] & 8))) continue;
			if (!((Reco_mu_SelectionType[iMuMinus] & 2) && (Reco_mu_SelectionType[iMuMinus] & 8))) continue;

			// passing hybrid-soft Id
			// hybridSoftMuons = (Reco_mu_nTrkWMea[iMuPlus] > 5) && (Reco_mu_nPixWMea[iMuPlus] > 0) && (fabs(Reco_mu_dxy[iMuPlus]) < 0.3) && (fabs(Reco_mu_dz[iMuPlus]) < 20.) && (Reco_mu_nTrkWMea[iMuMinus] > 5) && (Reco_mu_nPixWMea[iMuMinus] > 0) && (fabs(Reco_mu_dxy[iMuMinus]) < 0.3) && (fabs(Reco_mu_dz[iMuMinus]) < 20.);
			if (!((Reco_mu_nTrkWMea[iMuPlus] > 5) && (Reco_mu_nPixWMea[iMuPlus] > 0) && (fabs(Reco_mu_dxy[iMuPlus]) < 0.3) && (fabs(Reco_mu_dz[iMuPlus]) < 20.))) continue;
			if (!((Reco_mu_nTrkWMea[iMuMinus] > 5) && (Reco_mu_nPixWMea[iMuMinus] > 0) && (fabs(Reco_mu_dxy[iMuMinus]) < 0.3) && (fabs(Reco_mu_dz[iMuMinus]) < 20.))) continue;

			TLorentzVector* Reco_mupl_4mom = (TLorentzVector*)CloneArr_mu->At(iMuPlus);
			double Reco_mupl_eta = Reco_mupl_4mom->Eta();
			double Reco_mupl_pt = Reco_mupl_4mom->Pt();

			TLorentzVector* Reco_mumi_4mom = (TLorentzVector*)CloneArr_mu->At(iMuMinus);
			double Reco_mumi_eta = Reco_mumi_4mom->Eta();
			double Reco_mumi_pt = Reco_mumi_4mom->Pt();

			/// reweight for the data/MC reco pT spectrum discrepancies
			dimuonPtWeight = Get_RecoPtWeight(Reco_QQ_4mom->Rapidity(), Reco_QQ_4mom->Pt());

			// get positive muon's coordinates in the studied reference frames
			TVector3 muPlus_CS = MuPlusVector_CollinsSoper(*Reco_QQ_4mom, *Reco_mupl_4mom);

			TVector3 muPlus_HX = MuPlusVector_Helicity(*Reco_QQ_4mom, *Reco_mupl_4mom);

			// HLT filters
			bool mupl_L2Filter = ((Reco_mu_trig[iMuPlus] & ((ULong64_t)pow(2, gL2FilterBit))) == ((ULong64_t)pow(2, gL2FilterBit)));
			bool mupl_L3Filter = ((Reco_mu_trig[iMuPlus] & ((ULong64_t)pow(2, gL3FilterBit))) == ((ULong64_t)pow(2, gL3FilterBit)));
			bool mumi_L2Filter = ((Reco_mu_trig[iMuMinus] & ((ULong64_t)pow(2, gL2FilterBit))) == ((ULong64_t)pow(2, gL2FilterBit)));
			bool mumi_L3Filter = ((Reco_mu_trig[iMuMinus] & ((ULong64_t)pow(2, gL3FilterBit))) == ((ULong64_t)pow(2, gL3FilterBit)));

			bool mupl_isL2 = (mupl_L2Filter && !mupl_L3Filter);
			bool mupl_isL3 = (mupl_L2Filter && mupl_L3Filter);
			bool mumi_isL2 = (mumi_L2Filter && !mumi_L3Filter);
			bool mumi_isL3 = (mumi_L2Filter && mumi_L3Filter);

			passHLTFilterMuons = (mupl_L2Filter && mumi_L3Filter) || (mupl_L3Filter && mumi_L2Filter) || (mupl_L3Filter && mumi_L3Filter);

			// allGood = firesTrigger && dimuonMatching && goodVertexProba && passHLTFilterMuons && trackerAndGlobalMuons && hybridSoftMuons;

			// global AND tracker muons
			if (!((Reco_mu_SelectionType[iMuPlus] & 2) && (Reco_mu_SelectionType[iMuPlus] & 8))) continue;
			if (!((Reco_mu_SelectionType[iMuMinus] & 2) && (Reco_mu_SelectionType[iMuMinus] & 8))) continue;

			// passing hybrid-soft Id
			if (!((Reco_mu_nTrkWMea[iMuPlus] > 5) && (Reco_mu_nPixWMea[iMuPlus] > 0) && (fabs(Reco_mu_dxy[iMuPlus]) < 0.3) && (fabs(Reco_mu_dz[iMuPlus]) < 20.))) continue;
			if (!((Reco_mu_nTrkWMea[iMuMinus] > 5) && (Reco_mu_nPixWMea[iMuMinus] > 0) && (fabs(Reco_mu_dxy[iMuMinus]) < 0.3) && (fabs(Reco_mu_dz[iMuMinus]) < 20.))) continue;


			/// muon scale factors

			// muon trigger SF is tricky, need to know which muon passed which trigger filter

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

			// dimuon weight = product of the total scale factors
			dimuWeight_nominal = tnp_weight_trk_pbpb(Reco_mupl_eta, indexNominal) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexNominal) * dimuTrigWeight_nominal;

			/// variations for muon tracking SF (keeping the nominal values for muon Id and trigger)

			// tracking, syst up
			dimuWeight_trk_systUp = tnp_weight_trk_pbpb(Reco_mupl_eta, indexSystUp) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexSystUp) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexNominal) * dimuTrigWeight_nominal;

			// tracking, syst down
			dimuWeight_trk_systDown = tnp_weight_trk_pbpb(Reco_mupl_eta, indexSystDown) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexSystDown) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexNominal) * dimuTrigWeight_nominal;

			// tracking, stat up
			dimuWeight_trk_statUp = tnp_weight_trk_pbpb(Reco_mupl_eta, indexStatUp) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexStatUp) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexNominal) * dimuTrigWeight_nominal;

			// tracking, stat down
			dimuWeight_trk_statDown = tnp_weight_trk_pbpb(Reco_mupl_eta, indexStatDown) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexStatDown) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexNominal) * dimuTrigWeight_nominal;

			/// variations for muon Id SF (keeping the nominal values for tracking and trigger)

			// Id, syst up
			dimuWeight_muId_systUp = tnp_weight_trk_pbpb(Reco_mupl_eta, indexNominal) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexSystUp) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexSystUp) * dimuTrigWeight_nominal;

			// Id, syst down
			dimuWeight_muId_systDown = tnp_weight_trk_pbpb(Reco_mupl_eta, indexNominal) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexSystDown) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexSystDown) * dimuTrigWeight_nominal;

			// Id, stat up
			dimuWeight_muId_statUp = tnp_weight_trk_pbpb(Reco_mupl_eta, indexNominal) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexStatUp) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexStatUp) * dimuTrigWeight_nominal;

			// Id, stat down
			dimuWeight_muId_statDown = tnp_weight_trk_pbpb(Reco_mupl_eta, indexNominal) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexStatDown) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexStatDown) * dimuTrigWeight_nominal;

			/// variations for trigger SF (keeping the nominal values for tracking and muon Id)

			// trigger, syst up
			dimuWeight_trig_systUp = tnp_weight_trk_pbpb(Reco_mupl_eta, indexNominal) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexNominal) * dimuTrigWeight_systUp;

			// trigger, syst down
			dimuWeight_trig_systDown = tnp_weight_trk_pbpb(Reco_mupl_eta, indexNominal) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexNominal) * dimuTrigWeight_systDown;

			// trigger, stat up
			dimuWeight_trig_statUp = tnp_weight_trk_pbpb(Reco_mupl_eta, indexNominal) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexNominal) * dimuTrigWeight_statUp;

			// trigger, stat down
			dimuWeight_trig_statDown = tnp_weight_trk_pbpb(Reco_mupl_eta, indexNominal) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexNominal) * dimuTrigWeight_statDown;

			/// now the overall event weight
			weight = nColl * Gen_weight * dimuonPtWeight * dimuWeight_nominal;

			// propagate the scale factor uncertainties to the weight

			// uncertainty on the muon scale factors = quadratic sum of the uncertainties from the variation of all the sources
			double dimuWeightDownError_squared = pow(dimuWeight_trk_statDown - dimuWeight_nominal, 2.0) + pow(dimuWeight_muId_statDown - dimuWeight_nominal, 2.0) + pow(dimuWeight_trig_statDown - dimuWeight_nominal, 2.0) + pow(dimuWeight_trk_systDown - dimuWeight_nominal, 2.0) + pow(dimuWeight_muId_systDown - dimuWeight_nominal, 2.0) + pow(dimuWeight_trig_systDown - dimuWeight_nominal, 2.0);

			errorWeightDown = weight * sqrt(dimuWeightDownError_squared);
			// cout << "errorWeightDown: " << errorWeightDown << endl;

			//cout << "Weight from muon scale factors = " << dimuWeight_nominal << " - " << sqrt(dimuWeightDownError_squared) << " (" << 100. * sqrt(dimuWeightDownError_squared) / dimuWeight_nominal << "% relative)" << endl;

			double dimuWeightUpError_squared = pow(dimuWeight_trk_statUp - dimuWeight_nominal, 2.0) + pow(dimuWeight_muId_statUp - dimuWeight_nominal, 2.0) + pow(dimuWeight_trig_statUp - dimuWeight_nominal, 2.0) + pow(dimuWeight_trk_systUp - dimuWeight_nominal, 2.0) + pow(dimuWeight_muId_systUp - dimuWeight_nominal, 2.0) + pow(dimuWeight_trig_systUp - dimuWeight_nominal, 2.0);

			errorWeightUp = weight * sqrt(dimuWeightUpError_squared);

			// cout << "errorWeightUp: " << errorWeightUp << endl;
			// fill the dataset

			recoCat.setLabel((allGood) ? "selected" : "rejected");

			centVar = Centrality;

			polarWeightCS = 1 + lambdaTheta * TMath::Power(muPlus_CS.CosTheta(), 2) + lambdaPhi * TMath::Power(std::sin(muPlus_CS.Theta()), 2) * std::cos(2 * muPlus_CS.Phi()) + lambdaThetaPhi * std::sin(2 * muPlus_CS.Theta()) * std::cos(muPlus_CS.Phi());
			polarWeightHX = 1 + lambdaTheta * TMath::Power(muPlus_HX.CosTheta(), 2) + lambdaPhi * TMath::Power(std::sin(muPlus_HX.Theta()), 2) * std::cos(2 * muPlus_HX.Phi()) + lambdaThetaPhi * std::sin(2 * muPlus_HX.Theta()) * std::cos(muPlus_HX.Phi());

			totalWeightCS = weight * polarWeightCS;
			totalWeightHX = weight * polarWeightHX;

			eventWeightCSVar = weight * polarWeightCS;
			eventWeightHXVar = weight * polarWeightHX;

			eventWeightCSVar.setAsymError(errorWeightDown, errorWeightUp);
			eventWeightHXVar.setAsymError(errorWeightDown, errorWeightUp);

			massVar = Reco_QQ_4mom->M();
			yVar = fabs(Reco_QQ_4mom->Rapidity());
			ptVar = Reco_QQ_4mom->Pt();

			cosThetaLabVar = Reco_mupl_4mom->CosTheta();
			phiLabVar = Reco_mupl_4mom->Phi() * 180 / TMath::Pi();
			etaLabMuplVar = Reco_mupl_4mom->PseudoRapidity();
			etaLabMumiVar = Reco_mumi_4mom->PseudoRapidity();

			cosThetaCSVar = muPlus_CS.CosTheta();
			phiCSVar = muPlus_CS.Phi() * 180 / TMath::Pi();

			if (cosThetaCSVar.getVal() < 0) {
				// if phi value is smaller than -pi, add 2pi
				if ((phiCSVar.getVal() - 135) < -180)
					phiTildeCSVar.setVal(phiCSVar.getVal() + 225);
				else
					phiTildeCSVar.setVal(phiCSVar.getVal() - 135);
			}

			else if (cosThetaCSVar.getVal() > 0) {
				// if phi value is smaller than -pi, add 2pi
				if ((phiCSVar.getVal() - 45) < -180)
					phiTildeCSVar.setVal(phiCSVar.getVal() + 315);
				else
					phiTildeCSVar.setVal(phiCSVar.getVal() - 45);
			}

			cosThetaHXVar = muPlus_HX.CosTheta();
			phiHXVar = muPlus_HX.Phi() * 180 / TMath::Pi();

			if (cosThetaHXVar.getVal() < 0) {
				// if phi value is smaller than -pi, add 2pi
				if ((phiHXVar.getVal() - 135) < -180)
					phiTildeHXVar.setVal(phiHXVar.getVal() + 225);
				else
					phiTildeHXVar.setVal(phiHXVar.getVal() - 135);
			}

			else if (cosThetaHXVar.getVal() > 0) {
				// if phi value is smaller than -pi, add 2pi
				if ((phiHXVar.getVal() - 45) < -180)
					phiTildeHXVar.setVal(phiHXVar.getVal() + 315);
				else
					phiTildeHXVar.setVal(phiHXVar.getVal() - 45);
			}

			// datasetCS.add(RooArgSet(recoCat, centVar, eventWeightCSVar, massVar, yVar, ptVar, cosThetaLabVar, phiLabVar, etaLabMuplVar, etaLabMumiVar, cosThetaCSVar, phiCSVar, phiTildeCSVar, cosThetaHXVar, phiHXVar, phiTildeHXVar), totalWeightCS, errorWeightDown, errorWeightUp);
			// datasetHX.add(RooArgSet(recoCat, centVar, eventWeightHXVar, massVar, yVar, ptVar, cosThetaLabVar, phiLabVar, etaLabMuplVar, etaLabMumiVar, cosThetaCSVar, phiCSVar, phiTildeCSVar, cosThetaHXVar, phiHXVar, phiTildeHXVar), totalWeightHX, errorWeightDown, errorWeightUp);
			datasetCS.add(RooArgSet(centVar, eventWeightCSVar, massVar, yVar, ptVar, cosThetaLabVar, phiLabVar, etaLabMuplVar, etaLabMumiVar, cosThetaCSVar, phiCSVar, phiTildeCSVar, cosThetaHXVar, phiHXVar, phiTildeHXVar), totalWeightCS, errorWeightDown, errorWeightUp);
			datasetHX.add(RooArgSet(centVar, eventWeightHXVar, massVar, yVar, ptVar, cosThetaLabVar, phiLabVar, etaLabMuplVar, etaLabMumiVar, cosThetaCSVar, phiCSVar, phiTildeCSVar, cosThetaHXVar, phiHXVar, phiTildeHXVar), totalWeightHX, errorWeightDown, errorWeightUp);

		}
	}

	datasetCS.Write();
	datasetHX.Write();

	file.Close();

	infile->Close();

	cout << endl
	     << datasetCS.GetName() << " written in " << outputFileName << endl
	     << datasetHX.GetName() << " written in " << outputFileName << endl;

	return;
}

// check the dataset distributions

void draw2DHist(const char* refFrameName = "CS", Double_t lambdaTheta = 0, Double_t lambdaPhi = 0, Double_t lambdaThetaPhi = 0) {
	const char* FileName = Form("Y1SReconstructedMCWeightedDataset%s_Lambda_Theta%.2f_Phi%.2f_ThetaPhi%.2f.root", gMuonAccName, lambdaTheta, lambdaPhi, lambdaThetaPhi);

	TFile* file = openFile(FileName);

	RooDataSet* allDataset = (RooDataSet*)file->Get(Form("MCdataset%s", refFrameName));

	// import the dataset to a workspace
	RooWorkspace wspace("workspace");
	wspace.import(*allDataset);

	RooRealVar cosTheta = *wspace.var(CosThetaVarName(refFrameName));
	RooRealVar phi = *wspace.var(PhiVarName(refFrameName));
	RooRealVar phiTilde = *wspace.var(PhiTildeVarName(refFrameName));

	RooPlot* cosThetaframe = cosTheta.frame();
	RooPlot* phiframe = phi.frame();
	RooPlot* phiTildeframe = phiTilde.frame();

	allDataset->plotOn(cosThetaframe);
	allDataset->plotOn(phiframe);
	allDataset->plotOn(phiTildeframe);

	TCanvas* c = new TCanvas("c", "Canvas", 1200, 600);
	c->Divide(3);

	c->cd(1);
	cosThetaframe->Draw();

	c->cd(2);
	phiframe->Draw();

	c->cd(3);
	phiTildeframe->Draw();
}
