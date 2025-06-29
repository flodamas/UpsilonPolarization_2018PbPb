#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/Parameters/PhaseSpace.h"

#include "../Tools/Parameters/CentralityValues.h"
#include "../Tools/Parameters/EfficiencyWeights.h"
#include "../Tools/Parameters/MuonScaleFactors.h"

#include "../ReferenceFrameTransform/Transformations.h"

// using namespace RooFit;

// the reconstructed dimuons are fully weighted, including potential reweighting for polarization

void skimReconstructedMCWeighted(TString muonAccName = "UpsilonTriggerThresholds", Double_t lambdaTheta = 0, Double_t lambdaPhi = 0, Double_t lambdaThetaPhi = 0, Int_t iState = 1) {
	const char* inputFileName = Form("OniaTree_Y%dS_pThat2_HydjetDrumMB_miniAOD.root", iState);

	TFile* infile = TFile::Open(inputFileName, "READ");
	TTree* OniaTree = (TTree*)infile->Get("hionia/myTree");

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
	Short_t Reco_QQ_whichGen[1000];

	Short_t Gen_QQ_size;
	Short_t Gen_QQ_whichRec[1000];

	Short_t Gen_QQ_mupl_idx[1000];
	Short_t Gen_QQ_mumi_idx[1000];

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
	OniaTree->SetBranchAddress("Reco_QQ_whichGen", Reco_QQ_whichGen);

	OniaTree->SetBranchAddress("Reco_mu_trig", Reco_mu_trig);
	OniaTree->SetBranchAddress("Reco_mu_SelectionType", Reco_mu_SelectionType);

	OniaTree->SetBranchAddress("Reco_QQ_VtxProb", &Reco_QQ_VtxProb);
	OniaTree->SetBranchAddress("Reco_mu_nPixWMea", &Reco_mu_nPixWMea);
	OniaTree->SetBranchAddress("Reco_mu_nTrkWMea", &Reco_mu_nTrkWMea);
	OniaTree->SetBranchAddress("Reco_mu_dxy", &Reco_mu_dxy);
	OniaTree->SetBranchAddress("Reco_mu_dz", &Reco_mu_dz);

	/// RooDataSet output: one entry = one dimuon candidate!

	// weighting by event directly on the fly
	RooRealVar centVar("centrality", "event centrality", 0, 200);
	RooRealVar eventWeightCSVar("eventWeightCS", "event-by-event weight (Ncoll x MC gen weight x reco pT reweighting x muon scale factors x polarization in CS)", 0, 100000);
	RooRealVar eventWeightHXVar("eventWeightHX", "event-by-event weight (Ncoll x MC gen weight x reco pT reweighting x muon scale factors x polarization in HX)", 0, 100000);

	// Add these new RooRealVars above dataset definition
	RooRealVar errorWeightDownCSVar("errorWeightDownCS", "Downward error on event weight in CS", 0, 100000);
	RooRealVar errorWeightUpCSVar("errorWeightUpCS", "Upward error on event weight in CS", 0, 100000);

	RooRealVar errorWeightDownHXVar("errorWeightDownHX", "Downward error on event weight in HX", 0, 100000);
	RooRealVar errorWeightUpHXVar("errorWeightUpHX", "Upward error on event weight in HX", 0, 100000);

	Float_t lowMassCut = 8.8, highMassCut = 10.2;
	RooRealVar massVar("mass", gMassVarTitle, lowMassCut, highMassCut, gMassUnit);
	RooRealVar yVar("rapidity", gDimuonRapidityVarTitle, 0, 2.4);

	RooRealVar ptVar("pt", gDimuonPtVarTitle, 0, 100, gPtUnit);

	const char* refFrameName = "Lab";
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

	RooDataSet datasetCS("MCdatasetCS", "skimmed MC dataset in CS", RooArgSet(centVar, massVar, yVar, ptVar, cosThetaLabVar, phiLabVar, etaLabMuplVar, etaLabMumiVar, cosThetaCSVar, phiCSVar, phiTildeCSVar, cosThetaHXVar, phiHXVar, phiTildeHXVar, errorWeightDownCSVar, errorWeightUpCSVar), RooFit::WeightVar(eventWeightCSVar), RooFit::StoreAsymError(RooArgSet(eventWeightCSVar)));
	RooDataSet datasetHX("MCdatasetHX", "skimmed MC dataset in HX", RooArgSet(centVar, massVar, yVar, ptVar, cosThetaLabVar, phiLabVar, etaLabMuplVar, etaLabMumiVar, cosThetaCSVar, phiCSVar, phiTildeCSVar, cosThetaHXVar, phiHXVar, phiTildeHXVar, errorWeightDownHXVar, errorWeightUpHXVar), RooFit::WeightVar(eventWeightHXVar), RooFit::StoreAsymError(RooArgSet(eventWeightHXVar)));
	
	RooDataSet rawDatasetCS("MCrawDatasetCS", "skimmed unweighted MC dataset in CS", RooArgSet(centVar, massVar, yVar, ptVar, cosThetaLabVar, phiLabVar, etaLabMuplVar, etaLabMumiVar, cosThetaCSVar, phiCSVar, phiTildeCSVar, cosThetaHXVar, phiHXVar, phiTildeHXVar, errorWeightDownCSVar, errorWeightUpCSVar));
	RooDataSet rawDatasetHX("MCrawDatasetHX", "skimmed unweighted MC dataset in HX", RooArgSet(centVar, massVar, yVar, ptVar, cosThetaLabVar, phiLabVar, etaLabMuplVar, etaLabMumiVar, cosThetaCSVar, phiCSVar, phiTildeCSVar, cosThetaHXVar, phiHXVar, phiTildeHXVar, errorWeightDownHXVar, errorWeightUpHXVar));

	// RooDataSet datasetCS("MCdatasetCS", "skimmed MC dataset in CS", RooArgSet(centVar, eventWeightCSVar, massVar, yVar, ptVar, cosThetaLabVar, phiLabVar, etaLabMuplVar, etaLabMumiVar, cosThetaCSVar, phiCSVar, phiTildeCSVar, cosThetaHXVar, phiHXVar, phiTildeHXVar), RooFit::WeightVar("eventWeightCS"), RooFit::StoreAsymError(RooArgSet(eventWeightCSVar)));
	// RooDataSet datasetHX("MCdatasetHX", "skimmed MC dataset in HX", RooArgSet(centVar, eventWeightHXVar, massVar, yVar, ptVar, cosThetaLabVar, phiLabVar, etaLabMuplVar, etaLabMumiVar, cosThetaCSVar, phiCSVar, phiTildeCSVar, cosThetaHXVar, phiHXVar, phiTildeHXVar), RooFit::WeightVar("eventWeightHX"), RooFit::StoreAsymError(RooArgSet(eventWeightHXVar)));
	//datasetCS.Print();

	// TFitResultPtr HydjetGenFitParamCS[NPtBins], HydjetGenFitParamHX[NPtBins];
	// double HydjetFitLambdaThetaCS[NPtBins], HydjetFitLambdaPhiCS[NPtBins], HydjetFitLambdaThetaPhiCS[NPtBins];
	// double HydjetFitLambdaThetaHX[NPtBins], HydjetFitLambdaPhiHX[NPtBins], HydjetFitLambdaThetaPhiHX[NPtBins];
	
	// for (int iptBin = 0; iptBin < NPtBins; iptBin++) {
	// 	HydjetGenFitParamCS[iptBin] = getHydjetGenWeight("CS", lambdaTheta, lambdaPhi, lambdaThetaPhi, gPtBinning[iptBin], gPtBinning[iptBin + 1], 20, -1, 1, 18, -180, 180, kFALSE);
				
	// 	HydjetFitLambdaThetaCS[iptBin] = HydjetGenFitParamCS[iptBin]->Parameter(1);
	// 	HydjetFitLambdaPhiCS[iptBin] = HydjetGenFitParamCS[iptBin]->Parameter(2);
	// 	HydjetFitLambdaThetaPhiCS[iptBin] = HydjetGenFitParamCS[iptBin]->Parameter(3);

	// 	HydjetGenFitParamHX[iptBin] = getHydjetGenWeight("HX", lambdaTheta, lambdaPhi, lambdaThetaPhi, gPtBinning[iptBin], gPtBinning[iptBin + 1], 20, -1, 1, 18, -180, 180, kFALSE);

	// 	HydjetFitLambdaThetaHX[iptBin] = HydjetGenFitParamHX[iptBin]->Parameter(1);
	// 	HydjetFitLambdaPhiHX[iptBin] = HydjetGenFitParamHX[iptBin]->Parameter(2);
	// 	HydjetFitLambdaThetaPhiHX[iptBin] = HydjetGenFitParamHX[iptBin]->Parameter(3);
	// }

	// loop variables
	TLorentzVector* gen_QQ_LV = new TLorentzVector();
	TLorentzVector* gen_mupl_LV = new TLorentzVector();
	TLorentzVector* gen_mumi_LV = new TLorentzVector();

	Float_t nColl, weightCS = 0, weightHX = 0, HydjetGenWeightCS = 0, HydjetGenWeightHX = 0, polarWeightCS = 0, polarWeightHX = 0, dimuonPtWeight = 0, errorWeightDownCS = 0, errorWeightUpCS = 0, errorWeightDownHX = 0, errorWeightUpHX = 0;
	
	Bool_t isRecoMatched, passHLTFilterMuons;

	Int_t hiBin;

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

	for (Long64_t iEvent = 0; iEvent < (totEntries); iEvent++) {
		if (iEvent % 10000 == 0) {
			cout << Form("\rProcessing event %lld / %lld (%.0f%%)", iEvent, totEntries, 100. * iEvent / totEntries) << flush;
		}

		OniaTree->GetEntry(iEvent);

		// if (iEvent > 40) break; // for testing

		if (Centrality >= 2 * gCentralityBinMax) continue; // discard events with centrality > 90% in 2018 data

		if (!((HLTriggers & (ULong64_t)(1 << (gUpsilonHLTBit - 1))) == (ULong64_t)(1 << (gUpsilonHLTBit - 1)))) continue; // must fire the upsilon HLT path

		hiBin = GetHiBinFromhiHF(HFmean);

		nColl = FindNcoll(hiBin);

		for (int iGen = 0; iGen < Gen_QQ_size; iGen++) {
			
			// cut on the muon kinematics at gen level
			gen_QQ_LV = (TLorentzVector*)Gen_QQ_4mom->At(iGen);

			gen_mupl_LV = (TLorentzVector*)Gen_mu_4mom->At(Gen_QQ_mupl_idx[iGen]);
			gen_mumi_LV = (TLorentzVector*)Gen_mu_4mom->At(Gen_QQ_mumi_idx[iGen]);

			if (fabs(gen_QQ_LV->Rapidity()) < gRapidityMin || fabs(gen_QQ_LV->Rapidity()) > gRapidityMax) continue;

			if (gen_QQ_LV->Pt() > gPtMax) continue;

			// pt Weight at gen level
			double gen_QQ_pt = gen_QQ_LV->Pt();
			dimuonPtWeight = Get_GenPtWeight(gen_QQ_LV->Rapidity(), gen_QQ_pt);

			// reference frame transformation at gen level
			double cosThetaLab_gen = gen_mupl_LV->CosTheta();
			
			double phiLab_gen = 0;

			// if (isPhiFolded == kTRUE) {
			// 	phiLab_gen = fabs(gen_mupl_LV->Phi() * 180 / TMath::Pi());
			// } else {
				phiLab_gen = gen_mupl_LV->Phi() * 180 / TMath::Pi();
			// }

			TVector3 muPlus_CS_gen = MuPlusVector_CollinsSoper(*gen_QQ_LV, *gen_mupl_LV);

			double cosThetaCS_gen = muPlus_CS_gen.CosTheta();
			
			double phiCS_gen = 0;

			// if (isPhiFolded == kTRUE) {
			// 	phiCS_gen = fabs(muPlus_CS_gen.Phi() * 180 / TMath::Pi());
			// } else {
				phiCS_gen = muPlus_CS_gen.Phi() * 180 / TMath::Pi();
			// }

			TVector3 muPlus_HX_gen = MuPlusVector_Helicity(*gen_QQ_LV, *gen_mupl_LV);

			double cosThetaHX_gen = muPlus_HX_gen.CosTheta();
			
			double phiHX_gen = 0;

			// if (isPhiFolded == kTRUE) {
			// 	phiHX_gen = fabs(muPlus_HX_gen.Phi() * 180 / TMath::Pi());
			// } else {
				phiHX_gen = muPlus_HX_gen.Phi() * 180 / TMath::Pi();
			// }	

			// go to reco level
			Int_t iReco = Gen_QQ_whichRec[iGen];	

			if (Reco_QQ_sign[iReco] != 0) continue; // only opposite-sign muon pairs

			isRecoMatched = iReco > -1;

			if (!isRecoMatched) continue; // dimuon matching

			TLorentzVector* Reco_QQ_4mom = (TLorentzVector*)CloneArr_QQ->At(iReco);
			
			if (!((Reco_QQ_trig[iReco] & (ULong64_t)(1 << (gUpsilonHLTBit - 1))) == (ULong64_t)(1 << (gUpsilonHLTBit - 1)))) continue; // dimuon matching

			if (Reco_QQ_VtxProb[iReco] < 0.01) continue; // good common vertex proba

			// if (Reco_QQ_4mom->M() < lowMassCut || Reco_QQ_4mom->M() > highMassCut) continue; // speedup!

			/// single-muon selection criteria
			int iMuPlus = Reco_QQ_mupl_idx[iReco];
			int iMuMinus = Reco_QQ_mumi_idx[iReco];

			// global AND tracker muons
			if (!((Reco_mu_SelectionType[iMuPlus] & 2) && (Reco_mu_SelectionType[iMuPlus] & 8))) continue;
			if (!((Reco_mu_SelectionType[iMuMinus] & 2) && (Reco_mu_SelectionType[iMuMinus] & 8))) continue;

			// passing hybrid-soft Id
			if (!((Reco_mu_nTrkWMea[iMuPlus] > 5) && (Reco_mu_nPixWMea[iMuPlus] > 0) && (fabs(Reco_mu_dxy[iMuPlus]) < 0.3) && (fabs(Reco_mu_dz[iMuPlus]) < 20.))) continue;
			if (!((Reco_mu_nTrkWMea[iMuMinus] > 5) && (Reco_mu_nPixWMea[iMuMinus] > 0) && (fabs(Reco_mu_dxy[iMuMinus]) < 0.3) && (fabs(Reco_mu_dz[iMuMinus]) < 20.))) continue;

			// single-muon acceptance
			if (!MuonKinematicsWithinLimits(*gen_mupl_LV, muonAccName)) continue;

			if (!MuonKinematicsWithinLimits(*gen_mumi_LV, muonAccName)) continue;

			/// all selections applied, moving to weight computations

			TLorentzVector* Reco_mupl_4mom = (TLorentzVector*)CloneArr_mu->At(iMuPlus);
			double Reco_mupl_eta = Reco_mupl_4mom->Eta();
			double Reco_mupl_pt = Reco_mupl_4mom->Pt();

			TLorentzVector* Reco_mumi_4mom = (TLorentzVector*)CloneArr_mu->At(iMuMinus);
			double Reco_mumi_eta = Reco_mumi_4mom->Eta();
			double Reco_mumi_pt = Reco_mumi_4mom->Pt();

			/// reweight for the data/MC reco pT spectrum discrepancies
			// dimuonPtWeight = Get_RecoPtWeight(Reco_QQ_4mom->Rapidity(), Reco_QQ_4mom->Pt());

			// get positive muon's coordinates in the studied reference frames
			// TVector3 muPlus_CS_gen = MuPlusVector_CollinsSoper(*gen_QQ_LV, *gen_mupl_LV);
			TVector3 muPlus_CS = MuPlusVector_CollinsSoper(*Reco_QQ_4mom, *Reco_mupl_4mom);

			// TVector3 muPlus_HX_gen = MuPlusVector_Helicity(*gen_QQ_LV, *gen_mupl_LV);
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

			if (!passHLTFilterMuons) continue; 

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

			else {
				dimuTrigWeight_nominal = 1;
				dimuTrigWeight_systUp = 1;
				dimuTrigWeight_systDown = 1;
				dimuTrigWeight_statUp = 1;
				dimuTrigWeight_statDown = 1;
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

			// propagate the scale factor uncertainties to the weight

			// uncertainty on the muon scale factors = quadratic sum of the uncertainties from the variation of all the sources
			double dimuWeightDownError_squared = pow(dimuWeight_trk_statDown - dimuWeight_nominal, 2.0) + pow(dimuWeight_muId_statDown - dimuWeight_nominal, 2.0) + pow(dimuWeight_trig_statDown - dimuWeight_nominal, 2.0) + pow(dimuWeight_trk_systDown - dimuWeight_nominal, 2.0) + pow(dimuWeight_muId_systDown - dimuWeight_nominal, 2.0) + pow(dimuWeight_trig_systDown - dimuWeight_nominal, 2.0);

			double dimuWeightUpError_squared = pow(dimuWeight_trk_statUp - dimuWeight_nominal, 2.0) + pow(dimuWeight_muId_statUp - dimuWeight_nominal, 2.0) + pow(dimuWeight_trig_statUp - dimuWeight_nominal, 2.0) + pow(dimuWeight_trk_systUp - dimuWeight_nominal, 2.0) + pow(dimuWeight_muId_systUp - dimuWeight_nominal, 2.0) + pow(dimuWeight_trig_systUp - dimuWeight_nominal, 2.0);

			// errorWeightUp = weight * sqrt(dimuWeightUpError_squared);

			// cout << "errorWeightUp: " << errorWeightUp << endl;
			
			// fill the dataset
			centVar = Centrality;

			// int idx = 0;

			// if (gen_QQ_LV->Pt() > 0 && gen_QQ_LV->Pt() <= 2) idx = 0;
			// else if (gen_QQ_LV->Pt() > 2 && gen_QQ_LV->Pt() <= 6) idx = 1;
			// else if (gen_QQ_LV->Pt() > 6 && gen_QQ_LV->Pt() <= 12) idx = 2;
			// else if (gen_QQ_LV->Pt() > 12) idx = 3;

			// HydjetGenWeightCS = 1. / (1 + HydjetFitLambdaThetaCS[idx] * TMath::Power(muPlus_CS_gen.CosTheta(), 2) + HydjetFitLambdaPhiCS[idx] * TMath::Power(std::sin(muPlus_CS_gen.Theta()), 2) * std::cos(2 * muPlus_CS_gen.Phi()) + HydjetFitLambdaThetaPhiCS[idx] * std::sin(2 * muPlus_CS_gen.Theta()) * std::cos(muPlus_CS_gen.Phi()));
			// HydjetGenWeightHX = 1. / (1 + HydjetFitLambdaThetaHX[idx] * TMath::Power(muPlus_HX_gen.CosTheta(), 2) + HydjetFitLambdaPhiHX[idx] * TMath::Power(std::sin(muPlus_HX_gen.Theta()), 2) * std::cos(2 * muPlus_HX_gen.Phi()) + HydjetFitLambdaThetaPhiHX[idx] * std::sin(2 * muPlus_HX_gen.Theta()) * std::cos(muPlus_HX_gen.Phi())); 
			
			// HydjetGenWeightCS = 1.;
			// HydjetGenWeightHX = 1.;

			polarWeightCS = 1 + lambdaTheta * TMath::Power(muPlus_CS_gen.CosTheta(), 2) + lambdaPhi * TMath::Power(std::sin(muPlus_CS_gen.Theta()), 2) * std::cos(2 * muPlus_CS_gen.Phi()) + lambdaThetaPhi * std::sin(2 * muPlus_CS_gen.Theta()) * std::cos(muPlus_CS_gen.Phi());
			polarWeightHX = 1 + lambdaTheta * TMath::Power(muPlus_HX_gen.CosTheta(), 2) + lambdaPhi * TMath::Power(std::sin(muPlus_HX_gen.Theta()), 2) * std::cos(2 * muPlus_HX_gen.Phi()) + lambdaThetaPhi * std::sin(2 * muPlus_HX_gen.Theta()) * std::cos(muPlus_HX_gen.Phi());

			// totalWeightCS = weight * polarWeightCS;
			// totalWeightHX = weight * polarWeightHX;

			// totalWeightCSVar = weight * polarWeightCS * HydjetGenWeightCS;
			// totalWeightHXVar = weight * polarWeightHX * HydjetGenWeightHX;

			// totalWeightCSVar = weight * polarWeightCS;
			// totalWeightHXVar = weight * polarWeightHX;

			// totalWeightCSVar.setAsymError(errorWeightDown, errorWeightUp);
			// totalWeightHXVar.setAsymError(errorWeightDown, errorWeightUp);

			// eventWeightCSVar = weight * polarWeightCS * HydjetGenWeightCS;
			// eventWeightHXVar = weight * polarWeightHX * HydjetGenWeightHX;

			// /// now the overall event weight
			// weightCS = nColl * Gen_weight * dimuonPtWeight * dimuWeight_nominal * polarWeightCS;
			// weightHX = nColl * Gen_weight * dimuonPtWeight * dimuWeight_nominal * polarWeightHX;
			
			/// now the overall event weight
			weightCS = nColl * Gen_weight * dimuonPtWeight * polarWeightCS;
			weightHX = nColl * Gen_weight * dimuonPtWeight * polarWeightHX;

			// weightCS = 1;
			// weightHX = 1;

			// errorWeightDownCS = weightCS * sqrt(dimuWeightDownError_squared);
			// errorWeightUpCS = weightCS * sqrt(dimuWeightUpError_squared);
			// errorWeightDownHX = weightHX * sqrt(dimuWeightDownError_squared);
			// errorWeightUpHX = weightHX * sqrt(dimuWeightUpError_squared);

			errorWeightDownCS = 1;
			errorWeightUpCS = 1;
			errorWeightDownHX = 1;
			errorWeightUpHX = 1;

			eventWeightCSVar.setVal(weightCS);
			eventWeightHXVar.setVal(weightHX);

			eventWeightCSVar.setAsymError(errorWeightDownCS, errorWeightUpCS);
			eventWeightHXVar.setAsymError(errorWeightDownHX, errorWeightUpHX);

			errorWeightDownCSVar.setVal(errorWeightDownCS);
			errorWeightUpCSVar.setVal(errorWeightUpCS);
			errorWeightDownHXVar.setVal(errorWeightDownHX);
			errorWeightUpHXVar.setVal(errorWeightUpHX);
			
			massVar = Reco_QQ_4mom->M();
			yVar = fabs(Reco_QQ_4mom->Rapidity());
			ptVar = gen_QQ_LV->Pt();

			// cosThetaLabVar = Reco_mupl_4mom->CosTheta();
			// phiLabVar = Reco_mupl_4mom->Phi() * 180 / TMath::Pi();
			// etaLabMuplVar = Reco_mupl_4mom->PseudoRapidity();
			// etaLabMumiVar = Reco_mumi_4mom->PseudoRapidity();

			cosThetaLabVar = cosThetaLab_gen;
			phiLabVar = phiLab_gen;
			etaLabMuplVar = Reco_mupl_4mom->PseudoRapidity();
			etaLabMumiVar = Reco_mumi_4mom->PseudoRapidity();

			// cosThetaCSVar = muPlus_CS.CosTheta();
			// phiCSVar = muPlus_CS.Phi() * 180 / TMath::Pi();

			cosThetaCSVar = cosThetaCS_gen;
			phiCSVar = phiCS_gen;

			/// calculate phi tilde
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

			// cosThetaHXVar = muPlus_HX.CosTheta();
			// phiHXVar = muPlus_HX.Phi() * 180 / TMath::Pi();
			cosThetaHXVar = cosThetaHX_gen;
			phiHXVar = phiHX_gen;

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

			datasetCS.add(RooArgSet(centVar, massVar, yVar, ptVar, cosThetaLabVar, phiLabVar, etaLabMuplVar, etaLabMumiVar, cosThetaCSVar, phiCSVar, phiTildeCSVar, cosThetaHXVar, phiHXVar, phiTildeHXVar, errorWeightDownCSVar, errorWeightUpCSVar), eventWeightCSVar.getVal(), errorWeightDownCS, errorWeightUpCS);
			datasetHX.add(RooArgSet(centVar, massVar, yVar, ptVar, cosThetaLabVar, phiLabVar, etaLabMuplVar, etaLabMumiVar, cosThetaCSVar, phiCSVar, phiTildeCSVar, cosThetaHXVar, phiHXVar, phiTildeHXVar, errorWeightDownHXVar, errorWeightUpHXVar), eventWeightHXVar.getVal(), errorWeightDownHX, errorWeightUpHX);
			
			rawDatasetCS.add(RooArgSet(centVar, eventWeightCSVar, massVar, yVar, ptVar, cosThetaLabVar, phiLabVar, etaLabMuplVar, etaLabMumiVar, cosThetaCSVar, phiCSVar, phiTildeCSVar, cosThetaHXVar, phiHXVar, phiTildeHXVar));
			rawDatasetHX.add(RooArgSet(centVar, eventWeightHXVar, massVar, yVar, ptVar, cosThetaLabVar, phiLabVar, etaLabMuplVar, etaLabMumiVar, cosThetaCSVar, phiCSVar, phiTildeCSVar, cosThetaHXVar, phiHXVar, phiTildeHXVar));
			
			//	datasetCS.Print();

			// datasetCS.add(RooArgSet(centVar, eventWeightCSVar, massVar, yVar, ptVar, cosThetaLabVar, phiLabVar, etaLabMuplVar, etaLabMumiVar, cosThetaCSVar, phiCSVar, phiTildeCSVar, cosThetaHXVar, phiHXVar, phiTildeHXVar), totalWeightCS, errorWeightDown, errorWeightUp);
			// datasetHX.add(RooArgSet(centVar, eventWeightHXVar, massVar, yVar, ptVar, cosThetaLabVar, phiLabVar, etaLabMuplVar, etaLabMumiVar, cosThetaCSVar, phiCSVar, phiTildeCSVar, cosThetaHXVar, phiHXVar, phiTildeHXVar), totalWeightHX, errorWeightDown, errorWeightUp);
			// datasetCS.Print();
			// if (iEvent % 100 == 0) {
			// cout << "Event " << iEvent << std::endl;
			// 	datasetHX.Print("v");
			// 	datasetCS.Print("v");	

			// 	std::cout << "Passing event with mass = " << Reco_QQ_4mom->M()
			// 	<< ", pt = " << Reco_QQ_4mom->Pt()
			// 	<< ", y = " << Reco_QQ_4mom->Rapidity()
			// 	<< ", weightCS = " << weightCS
			// 	<< ", weightHX = " << weightHX
			// 	<< std::endl;
			// 	std::cout << endl;
			// }

		}
	}

	// const char* outputFileName = Form("Y%dSReconstructedMCWeightedDataset_%s_Lambda_Theta%.2f_Phi%.2f_ThetaPhi%.2f_noMassCut.root", iState, muonAccName.Data(), lambdaTheta, lambdaPhi, lambdaThetaPhi);
	// const char* outputFileName = Form("Y%dSReconstructedMCWeightedDataset_%s_Lambda_Theta%.2f_Phi%.2f_ThetaPhi%.2f.root", iState, muonAccName.Data(), lambdaTheta, lambdaPhi, lambdaThetaPhi);
	// const char* outputFileName = Form("Y%dSReconstructedMCWeightedDataset_%s_Lambda_Theta%.2f_Phi%.2f_ThetaPhi%.2f_gen.root", iState, muonAccName.Data(), lambdaTheta, lambdaPhi, lambdaThetaPhi);
	// const char* outputFileName = Form("Y%dSReconstructedMCWeightedDataset_%s_Lambda_Theta%.2f_Phi%.2f_ThetaPhi%.2f_gen.root", iState, muonAccName.Data(), lambdaTheta, lambdaPhi, lambdaThetaPhi);
	const char* outputFileName = Form("Y%dSReconstructedMCWeightedDataset_%s_Lambda_Theta%.2f_Phi%.2f_ThetaPhi%.2f_gen_noSF.root", iState, muonAccName.Data(), lambdaTheta, lambdaPhi, lambdaThetaPhi);
	// const char* outputFileName = Form("Y%dSReconstructedMCWeightedDataset_%s_Lambda_Theta%.2f_Phi%.2f_ThetaPhi%.2f_gen_noWeights.root", iState, muonAccName.Data(), lambdaTheta, lambdaPhi, lambdaThetaPhi);
	// const char* outputFileName = Form("Y%dSReconstructedMCWeightedDataset_%s_Lambda_Theta%.2f_Phi%.2f_ThetaPhi%.2f_HydjetWeight.root", iState, muonAccName.Data(), lambdaTheta, lambdaPhi, lambdaThetaPhi);

	TFile file(outputFileName, "RECREATE");

	datasetCS.Write();
	datasetHX.Write();

	rawDatasetCS.Write();
	rawDatasetHX.Write();

	file.Close();

	infile->Close();

	cout << endl
	     << datasetCS.GetName() << " written in " << outputFileName << endl
	     << datasetHX.GetName() << " written in " << outputFileName << endl
		 << rawDatasetCS.GetName() << " written in " << outputFileName << endl
		 << rawDatasetHX.GetName() << " written in " << outputFileName << endl;

	return;
}

// check the dataset distributions

void draw2DHist(const char* refFrameName = "CS", Double_t lambdaTheta = 0, Double_t lambdaPhi = 0, Double_t lambdaThetaPhi = 0) {
	const char* FileName = Form("Y1SReconstructedMCWeightedDataset_%s_Lambda_Theta%.2f_Phi%.2f_ThetaPhi%.2f.root", gMuonAccName.Data(), lambdaTheta, lambdaPhi, lambdaThetaPhi);

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
