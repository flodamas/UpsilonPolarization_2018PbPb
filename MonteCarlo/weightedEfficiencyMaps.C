#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "AccEffHelpers.h"

//#include "../Tools/Style/FitDistributions.h"
#include "../Tools/Style/Legends.h"

#include "../Tools/Parameters/CentralityValues.h"
#include "../Tools/Parameters/EfficiencyWeights.h"
#include "../Tools/Parameters/MuonScaleFactors.h"

#include "../ReferenceFrameTransform/Transformations.h"

// return the 2D map of the relative systematic efficiency uncertainty
TH2D* RelSystEffHist(TEfficiency hNominal, TEfficiency hTrk_systUp, TEfficiency hTrk_systDown, TEfficiency hMuId_systUp, TEfficiency hMuId_systDown, TEfficiency hTrig_systUp, TEfficiency hTrig_systDown, TEfficiency hTrk_statUp, TEfficiency hTrk_statDown, TEfficiency hMuId_statUp, TEfficiency hMuId_statDown, TEfficiency hTrig_statUp, TEfficiency hTrig_statDown) {
	// clone one of the input efficiency map, to make sure we get the correct binning
	TH2D* hTotalSyst = (TH2D*)hNominal.GetCopyPassedHisto(); // will rename it outside this function

	// useful loop variables
	int globalBin;
	double nominalEff;
	double trk_systEff, muId_systEff, trig_systEff, systEffSquared;
	double trk_statEff, muId_statEff, trig_statEff, statEffSquared;
	double finalSystEff;

	for (int iCosThetaBin = 1; iCosThetaBin <= hTotalSyst->GetNbinsX(); iCosThetaBin++) {
		for (int iPhiBin = 1; iPhiBin <= hTotalSyst->GetNbinsY(); iPhiBin++) {
			globalBin = hNominal.GetGlobalBin(iCosThetaBin, iPhiBin); // to run though the efficiency maps

			nominalEff = hNominal.GetEfficiency(globalBin);

			// compute the systematic uncertainties associated to the variation of the muon SF uncertainties

			// instructions can be found here: https://twiki.cern.ch/twiki/pub/CMS/HIMuonTagProbe/TnpHeaderFile.pdf#page=5

			// from muon SF systematics
			trk_systEff = max(abs(hTrk_systUp.GetEfficiency(globalBin) - nominalEff), abs(hTrk_systDown.GetEfficiency(globalBin) - nominalEff));

			muId_systEff = max(abs(hMuId_systUp.GetEfficiency(globalBin) - nominalEff), abs(hMuId_systDown.GetEfficiency(globalBin) - nominalEff));

			trig_systEff = max(abs(hTrig_systUp.GetEfficiency(globalBin) - nominalEff), abs(hTrig_systDown.GetEfficiency(globalBin) - nominalEff));

			systEffSquared = trk_systEff * trk_systEff + muId_systEff * muId_systEff + trig_systEff * trig_systEff;

			// from muon SF statistical uncertainties
			trk_statEff = max(abs(hTrk_statUp.GetEfficiency(globalBin) - nominalEff), abs(hTrk_statDown.GetEfficiency(globalBin) - nominalEff));

			muId_statEff = max(abs(hMuId_statUp.GetEfficiency(globalBin) - nominalEff), abs(hMuId_statDown.GetEfficiency(globalBin) - nominalEff));

			trig_statEff = max(abs(hTrig_statUp.GetEfficiency(globalBin) - nominalEff), abs(hTrig_statDown.GetEfficiency(globalBin) - nominalEff));

			statEffSquared = trk_statEff * trk_statEff + muId_statEff * muId_statEff + trig_statEff * trig_statEff;

			finalSystEff = sqrt(systEffSquared + statEffSquared);

			// store the RELATIVE SYSTEMATIC UNCERTAINTIES (not in % though) with respect to the nominal efficiency, more useful at the end
			hTotalSyst->SetBinContent(iCosThetaBin, iPhiBin, (nominalEff == 0) ? 0 : finalSystEff / nominalEff);

			//cout << "Efficiency in global bin " << globalBin << " = " << nominalEff << " +/- " << finalSystEff << endl;
		}
	}

	return hTotalSyst;
}

// return the 3D map of the relative systematic efficiency uncertainty
TH3D* RelSystEffHist3D(TEfficiency* hNominal, TEfficiency* hTrk_systUp, TEfficiency* hTrk_systDown, TEfficiency* hMuId_systUp, TEfficiency* hMuId_systDown, TEfficiency* hTrig_systUp, TEfficiency* hTrig_systDown, TEfficiency* hTrk_statUp, TEfficiency* hTrk_statDown, TEfficiency* hMuId_statUp, TEfficiency* hMuId_statDown, TEfficiency* hTrig_statUp, TEfficiency* hTrig_statDown, const char* refFrameName = "CS") {
	// clone one of the input efficiency map, to make sure we get the correct binning
	TH3D* hTotalSyst = (TH3D*)hNominal->GetCopyPassedHisto(); // will rename it outside this function

	// useful loop variables
	int globalBin;
	double nominalEff;
	double trk_systEff, muId_systEff, trig_systEff, systEffSquared;
	double trk_statEff, muId_statEff, trig_statEff, statEffSquared;
	double finalSystEff;

	for (int iCosThetaBin = 1; iCosThetaBin <= hTotalSyst->GetNbinsX(); iCosThetaBin++) {
		for (int iPhiBin = 1; iPhiBin <= hTotalSyst->GetNbinsY(); iPhiBin++) {
			for (int iPtBin = 1; iPtBin <= hTotalSyst->GetNbinsZ(); iPtBin++) {
				globalBin = hNominal->GetGlobalBin(iCosThetaBin, iPhiBin, iPtBin); // to run though the efficiency maps

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
				hTotalSyst->SetBinContent(iCosThetaBin, iPhiBin, iPtBin, (nominalEff == 0) ? 0 : finalSystEff / nominalEff);

				//cout << "Efficiency in global bin " << globalBin << " = " << nominalEff << " +/- " << finalSystEff << endl;
			}
		}
	}

	hTotalSyst->SetName(RelativeSystTEfficiency3DName(refFrameName));
	hTotalSyst->SetTitle(Form("Relative efficiency uncertainty;%s", TEfficiency3DAxisTitle(refFrameName)));

	return hTotalSyst;
}

const char* EfficiencyLegendText(int ptMin, int ptMax) {
	return Form("cent. %d-%d%%, %s < 2.4, %s", gCentralityBinMin, gCentralityBinMax, gDimuonRapidityVarTitle, DimuonPtRangeText(ptMin, ptMax));
}

void DrawEfficiencyMap(TEfficiency* effMap, Int_t ptMin, Int_t ptMax, int iState = gUpsilonState) {
	TCanvas* canvas = new TCanvas(effMap->GetName(), "", 700, 600);
	effMap->Draw("COLZ");

	CMS_lumi(canvas, Form("#varUpsilon(%dS) Hydjet-embedded MC", iState));

	TLatex legend;
	legend.SetTextAlign(22);
	legend.SetTextSize(0.042);
	legend.DrawLatexNDC(.5, .88, EfficiencyLegendText(ptMin, ptMax));

	gPad->Update();

	effMap->GetPaintedHistogram()->GetYaxis()->SetRangeUser(-190, 240);
	effMap->GetPaintedHistogram()->GetZaxis()->SetRangeUser(0, 0.9);

	effMap->GetPaintedHistogram()->GetXaxis()->CenterTitle();
	effMap->GetPaintedHistogram()->GetYaxis()->CenterTitle();

	gSystem->mkdir(Form("EfficiencyMaps/%dS", iState), kTRUE);
	canvas->SaveAs(Form("EfficiencyMaps/%dS/%s.png", iState, effMap->GetName()), "RECREATE");
}

// create the (cos theta, phi) map of the total efficiency (fully reweighted) for a given pT range

// this macro is based on mapUpsilonEfficiency.C with the difference that only "granular" distributions are made, following the new strategy of correcting the data BEFORE signal extraction from invariant mass fit

void weightedEfficiencyMaps(Int_t ptMin = 0, Int_t ptMax = 2, Int_t iState = gUpsilonState) {
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

	/// (cos theta, phi, pT) 3D distribution maps for CS and HX frames

	// Collins-Soper

	TEfficiency* hNominalEffCS = TEfficiency3D(NominalTEfficiency3DName("CS"), "CS", iState);

	TEfficiency* hCS_trk_systUp = TEfficiency3D("hCS_trk_systUp", "CS", iState);
	TEfficiency* hCS_trk_systDown = TEfficiency3D("hCS_trk_systDown", "CS", iState);
	TEfficiency* hCS_trk_statUp = TEfficiency3D("hCS_trk_statUp", "CS", iState);
	TEfficiency* hCS_trk_statDown = TEfficiency3D("hCS_trk_statDown", "CS", iState);

	TEfficiency* hCS_muId_systUp = TEfficiency3D("hCS_muId_systUp", "CS", iState);
	TEfficiency* hCS_muId_systDown = TEfficiency3D("hCS_muId_systDown", "CS", iState);
	TEfficiency* hCS_muId_statUp = TEfficiency3D("hCS_muId_statUp", "CS", iState);
	TEfficiency* hCS_muId_statDown = TEfficiency3D("hCS_muId_statDown", "CS", iState);

	TEfficiency* hCS_trig_systUp = TEfficiency3D("hCS_trig_systUp", "CS", iState);
	TEfficiency* hCS_trig_systDown = TEfficiency3D("hCS_trig_systDown", "CS", iState);
	TEfficiency* hCS_trig_statUp = TEfficiency3D("hCS_trig_statUp", "CS", iState);
	TEfficiency* hCS_trig_statDown = TEfficiency3D("hCS_trig_statDown", "CS", iState);

	// Helicity
	TEfficiency* hNominalEffHX = TEfficiency3D(NominalTEfficiency3DName("HX"), "HX", iState);

	TEfficiency* hHX_trk_systUp = TEfficiency3D("hHX_trk_systUp", "HX", iState);
	TEfficiency* hHX_trk_systDown = TEfficiency3D("hHX_trk_systDown", "HX", iState);
	TEfficiency* hHX_trk_statUp = TEfficiency3D("hHX_trk_statUp", "HX", iState);
	TEfficiency* hHX_trk_statDown = TEfficiency3D("hHX_trk_statDown", "HX", iState);

	TEfficiency* hHX_muId_systUp = TEfficiency3D("hHX_muId_systUp", "HX", iState);
	TEfficiency* hHX_muId_systDown = TEfficiency3D("hHX_muId_systDown", "HX", iState);
	TEfficiency* hHX_muId_statUp = TEfficiency3D("hHX_muId_statUp", "HX", iState);
	TEfficiency* hHX_muId_statDown = TEfficiency3D("hHX_muId_statDown", "HX", iState);

	TEfficiency* hHX_trig_systUp = TEfficiency3D("hHX_trig_systUp", "HX", iState);
	TEfficiency* hHX_trig_systDown = TEfficiency3D("hHX_trig_systDown", "HX", iState);
	TEfficiency* hHX_trig_statUp = TEfficiency3D("hHX_trig_statUp", "HX", iState);
	TEfficiency* hHX_trig_statDown = TEfficiency3D("hHX_trig_statDown", "HX", iState);

	// (cos theta, phi) for the given pt range
	TEfficiency* hEffCS2D = CosThetaPhiTEfficiency2D(ptMin, ptMax, "CS", iState);

	TEfficiency* hEffHX2D = CosThetaPhiTEfficiency2D(ptMin, ptMax, "HX", iState);

	// vs cos theta, to investigate

	TEfficiency* hEffCS1D = CosThetaTEfficiency1D(ptMin, ptMax, "CS", iState);

	TEfficiency* hEffHX1D = CosThetaTEfficiency1D(ptMin, ptMax, "HX", iState);

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

	// loop variables
	TLorentzVector* genLorentzVector = new TLorentzVector();
	TLorentzVector* recoLorentzVector = new TLorentzVector();

	double eventWeight, dimuonPtWeight, totalWeight;
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

			if (genLorentzVector->Pt() < gMuonPtCut) continue;

			if (fabs(genLorentzVector->Eta()) > 2.4) continue;

			// then negative muon
			genLorentzVector = (TLorentzVector*)Gen_mu_4mom->At(Gen_QQ_mumi_idx[iGen]);

			if (genLorentzVector->Pt() < gMuonPtCut) continue;

			if (fabs(genLorentzVector->Eta()) > 2.4) continue;

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

				allGood = firesTrigger && isRecoMatched && dimuonMatching && goodVertexProba && passHLTFilterMuons && trackerAndGlobalMuons && hybridSoftMuons;

				// get muon coordinates

				TLorentzVector* Reco_mupl_LV = (TLorentzVector*)CloneArr_mu->At(iMuPlus);
				double Reco_mupl_eta = Reco_mupl_LV->Eta();
				double Reco_mupl_pt = Reco_mupl_LV->Pt();

				TLorentzVector* Reco_mumi_LV = (TLorentzVector*)CloneArr_mu->At(iMuMinus);
				double Reco_mumi_eta = Reco_mumi_LV->Eta();
				double Reco_mumi_pt = Reco_mumi_LV->Pt();

				TVector3 muPlus_CS = MuPlusVector_CollinsSoper(*recoLorentzVector, *Reco_mupl_LV);
				double cosThetaCS = muPlus_CS.CosTheta();
				double phiCS = muPlus_CS.Phi() * 180 / TMath::Pi();

				TVector3 muPlus_HX = MuPlusVector_Helicity(*recoLorentzVector, *Reco_mupl_LV);
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

				totalWeight = eventWeight * dimuonPtWeight * dimuWeight_nominal;
				hNominalEffCS->FillWeighted(allGood, totalWeight, cosThetaCS, phiCS, reco_QQ_pt);
				hNominalEffHX->FillWeighted(allGood, totalWeight, cosThetaHX, phiHX, reco_QQ_pt);

				if (reco_QQ_pt > ptMin && reco_QQ_pt < ptMax) { // pT range of interest
					hEffCS2D->FillWeighted(allGood, totalWeight, cosThetaCS, phiCS);
					hEffHX2D->FillWeighted(allGood, totalWeight, cosThetaHX, phiHX);

					hEffCS1D->FillWeighted(allGood, totalWeight, cosThetaCS);
					hEffHX1D->FillWeighted(allGood, totalWeight, cosThetaHX);
				}

				/// variations for muon tracking SF (keeping the nominal efficiency for muon Id and trigger)

				// tracking, syst up
				dimuWeight_trk_systUp = tnp_weight_trk_pbpb(Reco_mupl_eta, indexSystUp) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexSystUp) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexNominal) * dimuTrigWeight_nominal;

				totalWeight = eventWeight * dimuonPtWeight * dimuWeight_trk_systUp;
				hCS_trk_systUp->FillWeighted(allGood, totalWeight, cosThetaCS, phiCS, reco_QQ_pt);
				hHX_trk_systUp->FillWeighted(allGood, totalWeight, cosThetaHX, phiHX, reco_QQ_pt);

				// tracking, syst down
				dimuWeight_trk_systDown = tnp_weight_trk_pbpb(Reco_mupl_eta, indexSystDown) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexSystDown) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexNominal) * dimuTrigWeight_nominal;

				totalWeight = eventWeight * dimuonPtWeight * dimuWeight_trk_systDown;
				hCS_trk_systDown->FillWeighted(allGood, totalWeight, cosThetaCS, phiCS, reco_QQ_pt);
				hHX_trk_systDown->FillWeighted(allGood, totalWeight, cosThetaHX, phiHX, reco_QQ_pt);

				// tracking, stat up
				dimuWeight_trk_statUp = tnp_weight_trk_pbpb(Reco_mupl_eta, indexStatUp) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexStatUp) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexNominal) * dimuTrigWeight_nominal;

				totalWeight = eventWeight * dimuonPtWeight * dimuWeight_trk_statUp;
				hCS_trk_statUp->FillWeighted(allGood, totalWeight, cosThetaCS, phiCS, reco_QQ_pt);
				hHX_trk_statUp->FillWeighted(allGood, totalWeight, cosThetaHX, phiHX, reco_QQ_pt);

				// tracking, stat down
				dimuWeight_trk_statDown = tnp_weight_trk_pbpb(Reco_mupl_eta, indexStatDown) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexStatDown) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexNominal) * dimuTrigWeight_nominal;

				totalWeight = eventWeight * dimuonPtWeight * dimuWeight_trk_statDown;
				hCS_trk_statDown->FillWeighted(allGood, totalWeight, cosThetaCS, phiCS, reco_QQ_pt);
				hHX_trk_statDown->FillWeighted(allGood, totalWeight, cosThetaHX, phiHX, reco_QQ_pt);

				/// variations for muon Id SF (keeping the nominal efficiency for tracking and trigger)

				// Id, syst up
				dimuWeight_muId_systUp = tnp_weight_trk_pbpb(Reco_mupl_eta, indexNominal) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexSystUp) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexSystUp) * dimuTrigWeight_nominal;

				totalWeight = eventWeight * dimuonPtWeight * dimuWeight_muId_systUp;
				hCS_muId_systUp->FillWeighted(allGood, totalWeight, cosThetaCS, phiCS, reco_QQ_pt);
				hHX_muId_systUp->FillWeighted(allGood, totalWeight, cosThetaHX, phiHX, reco_QQ_pt);

				// Id, syst down
				dimuWeight_muId_systDown = tnp_weight_trk_pbpb(Reco_mupl_eta, indexNominal) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexSystDown) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexSystDown) * dimuTrigWeight_nominal;

				totalWeight = eventWeight * dimuonPtWeight * dimuWeight_muId_systDown;
				hCS_muId_systDown->FillWeighted(allGood, totalWeight, cosThetaCS, phiCS, reco_QQ_pt);
				hHX_muId_systDown->FillWeighted(allGood, totalWeight, cosThetaHX, phiHX, reco_QQ_pt);

				// Id, stat up
				dimuWeight_muId_statUp = tnp_weight_trk_pbpb(Reco_mupl_eta, indexNominal) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexStatUp) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexStatUp) * dimuTrigWeight_nominal;

				totalWeight = eventWeight * dimuonPtWeight * dimuWeight_muId_statUp;
				hCS_muId_statUp->FillWeighted(allGood, totalWeight, cosThetaCS, phiCS, reco_QQ_pt);
				hHX_muId_statUp->FillWeighted(allGood, totalWeight, cosThetaHX, phiHX, reco_QQ_pt);

				// Id, stat down
				dimuWeight_muId_statDown = tnp_weight_trk_pbpb(Reco_mupl_eta, indexNominal) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexStatDown) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexStatDown) * dimuTrigWeight_nominal;

				totalWeight = eventWeight * dimuonPtWeight * dimuWeight_muId_statDown;
				hCS_muId_statDown->FillWeighted(allGood, totalWeight, cosThetaCS, phiCS, reco_QQ_pt);
				hHX_muId_statDown->FillWeighted(allGood, totalWeight, cosThetaHX, phiHX, reco_QQ_pt);

				/// variations for trigger SF (keeping the nominal efficiency for tracking and muon Id)

				// trigger, syst up
				dimuWeight_trig_systUp = tnp_weight_trk_pbpb(Reco_mupl_eta, indexNominal) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexNominal) * dimuTrigWeight_systUp;

				totalWeight = eventWeight * dimuonPtWeight * dimuWeight_trig_systUp;
				hCS_trig_systUp->FillWeighted(allGood, totalWeight, cosThetaCS, phiCS, reco_QQ_pt);
				hHX_trig_systUp->FillWeighted(allGood, totalWeight, cosThetaHX, phiHX, reco_QQ_pt);

				// trigger, syst down
				dimuWeight_trig_systDown = tnp_weight_trk_pbpb(Reco_mupl_eta, indexNominal) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexNominal) * dimuTrigWeight_systDown;

				totalWeight = eventWeight * dimuonPtWeight * dimuWeight_trig_systDown;
				hCS_trig_systDown->FillWeighted(allGood, totalWeight, cosThetaCS, phiCS, reco_QQ_pt);
				hHX_trig_systDown->FillWeighted(allGood, totalWeight, cosThetaHX, phiHX, reco_QQ_pt);

				// trigger, stat up
				dimuWeight_trig_statUp = tnp_weight_trk_pbpb(Reco_mupl_eta, indexNominal) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexNominal) * dimuTrigWeight_statUp;

				totalWeight = eventWeight * dimuonPtWeight * dimuWeight_trig_statUp;
				hCS_trig_statUp->FillWeighted(allGood, totalWeight, cosThetaCS, phiCS, reco_QQ_pt);
				hHX_trig_statUp->FillWeighted(allGood, totalWeight, cosThetaHX, phiHX, reco_QQ_pt);

				// trigger, stat down
				dimuWeight_trig_statDown = tnp_weight_trk_pbpb(Reco_mupl_eta, indexNominal) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexNominal) * dimuTrigWeight_statDown;

				totalWeight = eventWeight * dimuonPtWeight * dimuWeight_trig_statDown;
				hCS_trig_statDown->FillWeighted(allGood, totalWeight, cosThetaCS, phiCS, reco_QQ_pt);
				hHX_trig_statDown->FillWeighted(allGood, totalWeight, cosThetaHX, phiHX, reco_QQ_pt);
			}
		} // end of gen upsilon loop
	}

	cout << endl;

	gStyle->SetPadLeftMargin(.15);
	gStyle->SetPadRightMargin(0.18);
	gStyle->SetPalette(kRainBow);
	gStyle->SetNumberContours(256);

	/// display the nominal results

	DrawEfficiencyMap(hEffCS2D, ptMin, ptMax, iState);

	DrawEfficiencyMap(hEffHX2D, ptMin, ptMax, iState);

	/// compute the systematics in this macro since we have all the ingredients for that
	// instructions can be found here: https://twiki.cern.ch/twiki/pub/CMS/HIMuonTagProbe/TnpHeaderFile.pdf#page=5

	// store the RELATIVE SYSTEMATIC UNCERTAINTIES (not in % though) with respect to the nominal efficiency, more useful at the end

	gStyle->SetPadRightMargin(0.19);
	gStyle->SetTitleOffset(1.3, "Z");

	const char* legendText = EfficiencyLegendText(ptMin, ptMax);

	TLatex* legend = new TLatex(.5, .88, legendText);
	legend->SetTextAlign(22);
	legend->SetTextSize(0.042);

	// Collins-Soper
	auto* hSystCS = RelSystEffHist3D(hNominalEffCS, hCS_trk_systUp, hCS_trk_systDown, hCS_muId_systUp, hCS_muId_systDown, hCS_trig_systUp, hCS_trig_systDown, hCS_trk_statUp, hCS_trk_statDown, hCS_muId_statUp, hCS_muId_statDown, hCS_trig_statUp, hCS_trig_statDown, "CS");

	auto* canvasCSsyst = new TCanvas("canvasCSsyst", "", 700, 600);
	hSystCS->Draw("COLZ");

	CMS_lumi(canvasCSsyst, Form("#varUpsilon(%dS) Hydjet-embedded MC", iState));

	legend->DrawLatexNDC(.5, .88, legendText);

	gPad->Update();

	//hSystCS->GetYaxis()->SetRangeUser(-190, 240);

	//canvasCSsyst->SaveAs(Form("EfficiencyMaps/%dS/RelatSystEff_CS.png", gUpsilonState), "RECREATE");

	// Helicity
	auto* hSystHX = RelSystEffHist3D(hNominalEffHX, hHX_trk_systUp, hHX_trk_systDown, hHX_muId_systUp, hHX_muId_systDown, hHX_trig_systUp, hHX_trig_systDown, hHX_trk_statUp, hHX_trk_statDown, hHX_muId_statUp, hHX_muId_statDown, hHX_trig_statUp, hHX_trig_statDown, "HX");

	auto* canvasHXsyst = new TCanvas("canvasHXsyst", "", 700, 600);
	hSystHX->Draw("COLZ");

	CMS_lumi(canvasHXsyst, Form("#varUpsilon(%dS) Hydjet-embedded MC", iState));

	legend->DrawLatexNDC(.5, .88, legendText);

	gPad->Update();

	//hSystHX->GetYaxis()->SetRangeUser(-190, 240);

	//canvasHXsyst->SaveAs(Form("EfficiencyMaps/%dS/RelatSystEff_HX.png", gUpsilonState), "RECREATE");

	/// save the nominal efficiency results and the corresponding systematics in a file for later usage
	gSystem->mkdir(Form("EfficiencyMaps/%dS", gUpsilonState), kTRUE);
	const char* outputFileName = Form("EfficiencyMaps/%dS/EfficiencyResults.root", iState);
	TFile outputFile(outputFileName, "UPDATE");

	hNominalEffCS->Write();
	hNominalEffHX->Write();
	hSystCS->Write();
	hSystHX->Write();

	hEffCS2D->Write();
	hEffHX2D->Write();

	hEffCS1D->Write();
	hEffHX1D->Write();

	outputFile.Close();

	file->Close();

	if (BeVerbose) cout << "\nNominal efficiency and corresponding systematic uncertainty maps saved in " << outputFileName << endl;
}
