#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/Parameters/PhaseSpace.h"

#include "AccEffHelpers.h"

//#include "../Tools/Style/FitDistributions.h"
#include "../Tools/Style/Legends.h"

#include "../Tools/Parameters/CentralityValues.h"
#include "../Tools/Parameters/EfficiencyWeights.h"
#include "../Tools/Parameters/MuonScaleFactors.h"

#include "../ReferenceFrameTransform/Transformations.h"

// return the 2D map of the relative systematic efficiency uncertainty
TH2D* SystEffHist(TEfficiency* hNominal, TEfficiency* hTrk_systUp, TEfficiency* hTrk_systDown, TEfficiency* hMuId_systUp, TEfficiency* hMuId_systDown, TEfficiency* hTrig_systUp, TEfficiency* hTrig_systDown, TEfficiency* hTrk_statUp, TEfficiency* hTrk_statDown, TEfficiency* hMuId_statUp, TEfficiency* hMuId_statDown, TEfficiency* hTrig_statUp, TEfficiency* hTrig_statDown) {
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
			// cout << "globalBin: " << globalBin << endl;

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

			// store the nominal efficiency value + its uncertainty
			// hTotalSyst->SetBinContent(iCosThetaBin, iPhiBin, nominalEff + finalSystEff);
			// hTotalSyst->SetBinContent(iCosThetaBin, iPhiBin, finalSystEff / nominalEff);
			hTotalSyst->SetBinContent(iCosThetaBin, iPhiBin, finalSystEff);

			//cout << "Efficiency in global bin " << globalBin << " = " << nominalEff << " +/- " << finalSystEff << endl;
		}
	}

	return hTotalSyst;
}

// return the 3D map of the relative systematic efficiency uncertainty
TH3D* SystEffHist3D(TEfficiency* hNominal, TEfficiency* hTrk_systUp, TEfficiency* hTrk_systDown, TEfficiency* hMuId_systUp, TEfficiency* hMuId_systDown, TEfficiency* hTrig_systUp, TEfficiency* hTrig_systDown, TEfficiency* hTrk_statUp, TEfficiency* hTrk_statDown, TEfficiency* hMuId_statUp, TEfficiency* hMuId_statDown, TEfficiency* hTrig_statUp, TEfficiency* hTrig_statDown, const char* refFrameName = "CS", Double_t lambdaTheta = 0, Double_t lambdaPhi = 0, Double_t lambdaThetaPhi = 0) {
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

				// store the nominal efficiency value + its uncertainty
				// hTotalSyst->SetBinContent(iCosThetaBin, iPhiBin, iPtBin, nominalEff + finalSystEff);
				// hTotalSyst->SetBinContent(iCosThetaBin, iPhiBin, iPtBin, finalSystEff / nominalEff);
				hTotalSyst->SetBinContent(iCosThetaBin, iPhiBin, iPtBin, finalSystEff);

				// cout << "iCosTheta " << iCosThetaBin << ", iPhi " << iPhiBin << ", iPt " << iPtBin << endl;
				// cout << "systEffSquared = " << systEffSquared << ", statEffSquared = " << statEffSquared << endl;
				// cout << "Efficiency in global bin " << globalBin << " = " << nominalEff << " +/- " << finalSystEff << endl;

				// cout << hNominal->GetPassedHistogram()->GetBinContent(iCosThetaBin, iPhiBin, iPtBin) << endl;
				// cout << hNominal->GetTotalHistogram()->GetBinContent(iCosThetaBin, iPhiBin, iPtBin) << endl;
			}
		}
	}

	hTotalSyst->SetName(SystTEfficiency3DName(refFrameName, lambdaTheta, lambdaPhi, lambdaThetaPhi));
	hTotalSyst->SetTitle(Form("Nominal efficiency + muon SF variation;%s", TEfficiency3DAxisTitle(refFrameName)));

	return hTotalSyst;
}

const char* EfficiencyLegendText(int ptMin, int ptMax) {
	return Form("cent. %d-%d%%, %s < 2.4, %s", gCentralityBinMin, gCentralityBinMax, gDimuonRapidityVarTitle, DimuonPtRangeText(ptMin, ptMax));
}

void DrawEfficiencyMap(TEfficiency* effMap, Int_t ptMin, Int_t ptMax, TString muonAccName, int iState = gUpsilonState, Bool_t isPhiFolded = kFALSE) {
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

	const char* path = EfficiencyResultsPath(muonAccName.Data());

	gSystem->mkdir(path, kTRUE);

	canvas->SaveAs(Form("%s/%s%s.png", path, effMap->GetName(), isPhiFolded ? "" : "_fullPhi"), "RECREATE");
}

// create the (cos theta, phi) map of the total efficiency (fully reweighted) for a given pT range

void weightedEfficiencyMaps(Int_t ptMin = 0, Int_t ptMax = 2, TString muonAccName = "UpsilonTriggerThresholds", Float_t lambdaTheta = 0, Float_t lambdaPhi = 0, Float_t lambdaThetaPhi = 0, Bool_t isPhiFolded = kFALSE, Int_t iState = gUpsilonState) { // please always Bool_t isPhiFolded = kFALSE for the efficiency matrix. Phifolding is applied later when rebinning the efficiency matrix
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

	/// (cos theta, phi, pT) 3D distribution maps for Lab, CS and HX frames
	TEfficiency* hNominalEffLab = TEfficiency3D(NominalTEfficiency3DName("Lab", lambdaTheta, lambdaPhi, lambdaThetaPhi), "Lab", iState);

	// Collins-Soper
	TEfficiency* hNominalEffCS = TEfficiency3D(NominalTEfficiency3DName("CS", lambdaTheta, lambdaPhi, lambdaThetaPhi), "CS", iState);

	TEfficiency* hCS_trk_systUp = TEfficiency3D(Form("hCS_trk_systUp_%s", PolaWeightName(lambdaTheta, lambdaPhi, lambdaThetaPhi)), "CS", iState);
	TEfficiency* hCS_trk_systDown = TEfficiency3D(Form("hCS_trk_systDown_%s", PolaWeightName(lambdaTheta, lambdaPhi, lambdaThetaPhi)), "CS", iState);
	TEfficiency* hCS_trk_statUp = TEfficiency3D(Form("hCS_trk_statUp_%s", PolaWeightName(lambdaTheta, lambdaPhi, lambdaThetaPhi)), "CS", iState);
	TEfficiency* hCS_trk_statDown = TEfficiency3D(Form("hCS_trk_statDown_%s", PolaWeightName(lambdaTheta, lambdaPhi, lambdaThetaPhi)), "CS", iState);

	TEfficiency* hCS_muId_systUp = TEfficiency3D(Form("hCS_muId_systUp_%s", PolaWeightName(lambdaTheta, lambdaPhi, lambdaThetaPhi)), "CS", iState);
	TEfficiency* hCS_muId_systDown = TEfficiency3D(Form("hCS_muId_systDown_%s", PolaWeightName(lambdaTheta, lambdaPhi, lambdaThetaPhi)), "CS", iState);
	TEfficiency* hCS_muId_statUp = TEfficiency3D(Form("hCS_muId_statUp_%s", PolaWeightName(lambdaTheta, lambdaPhi, lambdaThetaPhi)), "CS", iState);
	TEfficiency* hCS_muId_statDown = TEfficiency3D(Form("hCS_muId_statDown_%s", PolaWeightName(lambdaTheta, lambdaPhi, lambdaThetaPhi)), "CS", iState);

	TEfficiency* hCS_trig_systUp = TEfficiency3D(Form("hCS_trig_systUp_%s", PolaWeightName(lambdaTheta, lambdaPhi, lambdaThetaPhi)), "CS", iState);
	TEfficiency* hCS_trig_systDown = TEfficiency3D(Form("hCS_trig_systDown_%s", PolaWeightName(lambdaTheta, lambdaPhi, lambdaThetaPhi)), "CS", iState);
	TEfficiency* hCS_trig_statUp = TEfficiency3D(Form("hCS_trig_statUp_%s", PolaWeightName(lambdaTheta, lambdaPhi, lambdaThetaPhi)), "CS", iState);
	TEfficiency* hCS_trig_statDown = TEfficiency3D(Form("hCS_trig_statDown_%s", PolaWeightName(lambdaTheta, lambdaPhi, lambdaThetaPhi)), "CS", iState);

	TEfficiency* hCS_total_systUp = TEfficiency3D(Form("hCS_total_systUp_%s", PolaWeightName(lambdaTheta, lambdaPhi, lambdaThetaPhi)), "CS", iState);
	TEfficiency* hCS_total_systDown = TEfficiency3D(Form("hCS_total_systDown_%s", PolaWeightName(lambdaTheta, lambdaPhi, lambdaThetaPhi)), "CS", iState);
	TEfficiency* hCS_total_statUp = TEfficiency3D(Form("hCS_total_statUp_%s", PolaWeightName(lambdaTheta, lambdaPhi, lambdaThetaPhi)), "CS", iState);
	TEfficiency* hCS_total_statDown = TEfficiency3D(Form("hCS_total_statDown_%s", PolaWeightName(lambdaTheta, lambdaPhi, lambdaThetaPhi)), "CS", iState);

	// Helicity
	TEfficiency* hNominalEffHX = TEfficiency3D(NominalTEfficiency3DName("HX", lambdaTheta, lambdaPhi, lambdaThetaPhi), "HX", iState);

	TEfficiency* hHX_trk_systUp = TEfficiency3D(Form("hHX_trk_systUp_%s", PolaWeightName(lambdaTheta, lambdaPhi, lambdaThetaPhi)), "HX", iState);
	TEfficiency* hHX_trk_systDown = TEfficiency3D(Form("hHX_trk_systDown_%s", PolaWeightName(lambdaTheta, lambdaPhi, lambdaThetaPhi)), "HX", iState);
	TEfficiency* hHX_trk_statUp = TEfficiency3D(Form("hHX_trk_statUp_%s", PolaWeightName(lambdaTheta, lambdaPhi, lambdaThetaPhi)), "HX", iState);
	TEfficiency* hHX_trk_statDown = TEfficiency3D(Form("hHX_trk_statDown_%s", PolaWeightName(lambdaTheta, lambdaPhi, lambdaThetaPhi)), "HX", iState);

	TEfficiency* hHX_muId_systUp = TEfficiency3D(Form("hHX_muId_systUp_%s", PolaWeightName(lambdaTheta, lambdaPhi, lambdaThetaPhi)), "HX", iState);
	TEfficiency* hHX_muId_systDown = TEfficiency3D(Form("hHX_muId_systDown_%s", PolaWeightName(lambdaTheta, lambdaPhi, lambdaThetaPhi)), "HX", iState);
	TEfficiency* hHX_muId_statUp = TEfficiency3D(Form("hHX_muId_statUp_%s", PolaWeightName(lambdaTheta, lambdaPhi, lambdaThetaPhi)), "HX", iState);
	TEfficiency* hHX_muId_statDown = TEfficiency3D(Form("hHX_muId_statDown_%s", PolaWeightName(lambdaTheta, lambdaPhi, lambdaThetaPhi)), "HX", iState);

	TEfficiency* hHX_trig_systUp = TEfficiency3D(Form("hHX_trig_systUp_%s", PolaWeightName(lambdaTheta, lambdaPhi, lambdaThetaPhi)), "HX", iState);
	TEfficiency* hHX_trig_systDown = TEfficiency3D(Form("hHX_trig_systDown_%s", PolaWeightName(lambdaTheta, lambdaPhi, lambdaThetaPhi)), "HX", iState);
	TEfficiency* hHX_trig_statUp = TEfficiency3D(Form("hHX_trig_statUp_%s", PolaWeightName(lambdaTheta, lambdaPhi, lambdaThetaPhi)), "HX", iState);
	TEfficiency* hHX_trig_statDown = TEfficiency3D(Form("hHX_trig_statDown_%s", PolaWeightName(lambdaTheta, lambdaPhi, lambdaThetaPhi)), "HX", iState);

	TEfficiency* hHX_total_systUp = TEfficiency3D(Form("hHX_total_systUp_%s", PolaWeightName(lambdaTheta, lambdaPhi, lambdaThetaPhi)), "HX", iState);
	TEfficiency* hHX_total_systDown = TEfficiency3D(Form("hHX_total_systDown_%s", PolaWeightName(lambdaTheta, lambdaPhi, lambdaThetaPhi)), "HX", iState);
	TEfficiency* hHX_total_statUp = TEfficiency3D(Form("hHX_total_statUp_%s", PolaWeightName(lambdaTheta, lambdaPhi, lambdaThetaPhi)), "HX", iState);
	TEfficiency* hHX_total_statDown = TEfficiency3D(Form("hHX_total_statDown_%s", PolaWeightName(lambdaTheta, lambdaPhi, lambdaThetaPhi)), "HX", iState);

	// (cos theta, phi) for the given pt range
	TEfficiency* hEffCS2D = CosThetaPhiTEfficiency2D("CS", ptMin, ptMax, iState, lambdaTheta, lambdaPhi, lambdaThetaPhi);

	TEfficiency* hEffHX2D = CosThetaPhiTEfficiency2D("HX", ptMin, ptMax, iState, lambdaTheta, lambdaPhi, lambdaThetaPhi);

	// vs cos theta, to investigate

	TEfficiency* hEffCS1D = CosThetaTEfficiency1D("CS", ptMin, ptMax, iState, false, lambdaTheta, lambdaPhi, lambdaThetaPhi);

	TEfficiency* hEffHX1D = CosThetaTEfficiency1D("HX", ptMin, ptMax, iState, false, lambdaTheta, lambdaPhi, lambdaThetaPhi);

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

	double eventWeight, dimuonPtWeight, totalWeightCS, totalWeightHX, totalWeightLab;
	double dimuTrigWeight_nominal = -1,
	       dimuTrigWeight_systUp = -1, dimuTrigWeight_systDown = -1,
	       dimuTrigWeight_statUp = -1, dimuTrigWeight_statDown = -1;

	Float_t weightCS = 0, weightHX = 0;

	Bool_t allGood, firesTrigger, isRecoMatched, dimuonMatching, goodVertexProba, passHLTFilterMuons, trackerAndGlobalMuons, hybridSoftMuons;

	Int_t hiBin;

	Long64_t totEntries = OniaTree->GetEntries();

	for (Long64_t iEvent = 0; iEvent < totEntries; iEvent++) {
		if (iEvent % 10000 == 0) {
			cout << Form("\rProcessing event %lld / %lld (%.0f%%)", iEvent, totEntries, 100. * iEvent / totEntries) << flush;
		}

		// if (iEvent > 100) break; // for testing

		OniaTree->GetEntry(iEvent);
		//genLorentzVector->Clear();

		// event selection

		if (Centrality >= 2 * gCentralityBinMax) continue; // discard events with centrality >= 90% in 2018 data

		firesTrigger = ((HLTriggers & (ULong64_t)(1 << (gUpsilonHLTBit - 1))) == (ULong64_t)(1 << (gUpsilonHLTBit - 1)));

		hiBin = GetHiBinFromhiHF(HFmean);

		eventWeight = Gen_weight * FindNcoll(hiBin); // * Get_zPV_weight(zVtx); */

		// loop over all gen upsilons
		for (int iGen = 0; iGen < Gen_QQ_size; iGen++) {
			gen_QQ_LV = (TLorentzVector*)Gen_QQ_4mom->At(iGen);

			// fiducial region
			// if (genLorentzVector->Pt() < ptMin || genLorentzVector->Pt() > ptMax) continue; // pt bin of interest

			if (fabs(gen_QQ_LV->Rapidity()) < gRapidityMin || fabs(gen_QQ_LV->Rapidity()) > gRapidityMax) continue;

			if (gen_QQ_LV->Pt() > gPtMax) continue;

			// single-muon acceptance
			// positive muon first
			gen_mupl_LV = (TLorentzVector*)Gen_mu_4mom->At(Gen_QQ_mupl_idx[iGen]);

			if (!MuonKinematicsWithinLimits(*gen_mupl_LV, muonAccName)) continue;

			// then negative muon
			gen_mumi_LV = (TLorentzVector*)Gen_mu_4mom->At(Gen_QQ_mumi_idx[iGen]);

			if (!MuonKinematicsWithinLimits(*gen_mumi_LV, muonAccName)) continue;

			// pt Weight at gen level
			double gen_QQ_pt = gen_QQ_LV->Pt();
			dimuonPtWeight = Get_GenPtWeight(gen_QQ_LV->Rapidity(), gen_QQ_pt);

			// reference frame transformation at gen level
			double cosThetaLab_gen = gen_mupl_LV->CosTheta();
			
			double phiLab_gen = 0;

			if (isPhiFolded == kTRUE) {
				phiLab_gen = fabs(gen_mupl_LV->Phi() * 180 / TMath::Pi());
			} else {
				phiLab_gen = gen_mupl_LV->Phi() * 180 / TMath::Pi();
			}

			TVector3 muPlus_CS_gen = MuPlusVector_CollinsSoper(*gen_QQ_LV, *gen_mupl_LV);

			double cosThetaCS_gen = muPlus_CS_gen.CosTheta();
			
			double phiCS_gen = 0;

			if (isPhiFolded == kTRUE) {
				phiCS_gen = fabs(muPlus_CS_gen.Phi() * 180 / TMath::Pi());
			} else {
				phiCS_gen = muPlus_CS_gen.Phi() * 180 / TMath::Pi();
			}

			TVector3 muPlus_HX_gen = MuPlusVector_Helicity(*gen_QQ_LV, *gen_mupl_LV);

			double cosThetaHX_gen = muPlus_HX_gen.CosTheta();
			
			double phiHX_gen = 0;

			if (isPhiFolded == kTRUE) {
				phiHX_gen = fabs(muPlus_HX_gen.Phi() * 180 / TMath::Pi());
			} else {
				phiHX_gen = muPlus_HX_gen.Phi() * 180 / TMath::Pi();
			}			

			// polarization weights at gen level
			if (isPhiFolded == kTRUE) {
				weightCS = 1 + lambdaTheta * TMath::Power(muPlus_CS_gen.CosTheta(), 2) + lambdaPhi * TMath::Power(std::sin(muPlus_CS_gen.Theta()), 2) * std::cos(2 * fabs(muPlus_CS_gen.Phi())) + lambdaThetaPhi * std::sin(2 * muPlus_CS_gen.Theta()) * std::cos(fabs(muPlus_CS_gen.Phi()));
				weightHX = 1 + lambdaTheta * TMath::Power(muPlus_HX_gen.CosTheta(), 2) + lambdaPhi * TMath::Power(std::sin(muPlus_HX_gen.Theta()), 2) * std::cos(2 * fabs(muPlus_HX_gen.Phi())) + lambdaThetaPhi * std::sin(2 * muPlus_HX_gen.Theta()) * std::cos(fabs(muPlus_HX_gen.Phi()));
			}

			else {
				weightCS = 1 + lambdaTheta * TMath::Power(muPlus_CS_gen.CosTheta(), 2) + lambdaPhi * TMath::Power(std::sin(muPlus_CS_gen.Theta()), 2) * std::cos(2 * muPlus_CS_gen.Phi()) + lambdaThetaPhi * std::sin(2 * muPlus_CS_gen.Theta()) * std::cos(muPlus_CS_gen.Phi());
				weightHX = 1 + lambdaTheta * TMath::Power(muPlus_HX_gen.CosTheta(), 2) + lambdaPhi * TMath::Power(std::sin(muPlus_HX_gen.Theta()), 2) * std::cos(2 * muPlus_HX_gen.Phi()) + lambdaThetaPhi * std::sin(2 * muPlus_HX_gen.Theta()) * std::cos(muPlus_HX_gen.Phi());
			}
			
			/// go to reco level (numerator)
			Int_t iReco = Gen_QQ_whichRec[iGen];

			if (Reco_QQ_sign[iReco] != 0) continue; // only opposite-sign muon pairs			

			/// all the reconstructed upsilons must pass the conditions below!

			isRecoMatched = iReco > -1;

			if (isRecoMatched) {
				
				recoLorentzVector = (TLorentzVector*)CloneArr_QQ->At(iReco);
				double reco_QQ_pt = recoLorentzVector->Pt();

				// dimuonPtWeight = Get_RecoPtWeight(recoLorentzVector->Rapidity(), reco_QQ_pt); // pt Weight at reco level

				dimuonMatching = (Reco_QQ_trig[iReco] & (ULong64_t)(1 << (gUpsilonHLTBit - 1))) == (ULong64_t)(1 << (gUpsilonHLTBit - 1));

				goodVertexProba = Reco_QQ_VtxProb[iReco] > 0.01;

				/// single-muon selection criteria
				int iMuPlus = Reco_QQ_mupl_idx[iReco];
				int iMuMinus = Reco_QQ_mumi_idx[iReco];

				/// HLT filters
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

				// get muon coordinates at reco level
				TLorentzVector* Reco_mupl_LV = (TLorentzVector*)CloneArr_mu->At(iMuPlus);
				double Reco_mupl_eta = Reco_mupl_LV->Eta();
				double Reco_mupl_pt = Reco_mupl_LV->Pt();

				TLorentzVector* Reco_mumi_LV = (TLorentzVector*)CloneArr_mu->At(iMuMinus);
				double Reco_mumi_eta = Reco_mumi_LV->Eta();
				double Reco_mumi_pt = Reco_mumi_LV->Pt();

				// double cosThetaLab = Reco_mupl_LV->CosTheta();

				// double phiLab = 0;

				// if (isPhiFolded == kTRUE) {
				// 	// phiLab_gen = fabs(gen_mupl_LV->Phi() * 180 / TMath::Pi());
				// 	phiLab = fabs(Reco_mupl_LV->Phi() * 180 / TMath::Pi());
				// } else {
				// 	// phiLab_gen = gen_mupl_LV->Phi() * 180 / TMath::Pi();
				// 	phiLab = Reco_mupl_LV->Phi() * 180 / TMath::Pi();
				// }

				// // TVector3 muPlus_CS_gen = MuPlusVector_CollinsSoper(*gen_QQ_LV, *gen_mupl_LV);
				// TVector3 muPlus_CS = MuPlusVector_CollinsSoper(*recoLorentzVector, *Reco_mupl_LV);

				// // double cosThetaCS_gen = muPlus_CS_gen.CosTheta();
				// double cosThetaCS = muPlus_CS.CosTheta();
				// // double phiCS_gen = 0;
				// double phiCS = 0;

				// if (isPhiFolded == kTRUE) {
				// 	// phiCS_gen = fabs(muPlus_CS_gen.Phi() * 180 / TMath::Pi());
				// 	phiCS = fabs(muPlus_CS.Phi() * 180 / TMath::Pi());
				// } else {
				// 	// phiCS_gen = muPlus_CS_gen.Phi() * 180 / TMath::Pi();
				// 	phiCS = muPlus_CS.Phi() * 180 / TMath::Pi();
				// }

				// // TVector3 muPlus_HX_gen = MuPlusVector_Helicity(*gen_QQ_LV, *gen_mupl_LV);
				// TVector3 muPlus_HX = MuPlusVector_Helicity(*recoLorentzVector, *Reco_mupl_LV);

				// // double cosThetaHX_gen = muPlus_HX_gen.CosTheta();
				// double cosThetaHX = muPlus_HX.CosTheta();
				// // double phiHX_gen = 0;
				// double phiHX = 0;

				// if (isPhiFolded == kTRUE) {
				// 	// phiHX_gen = fabs(muPlus_HX_gen.Phi() * 180 / TMath::Pi());
				// 	phiHX = fabs(muPlus_HX.Phi() * 180 / TMath::Pi());
				// } else {
				// 	phiHX = muPlus_HX.Phi() * 180 / TMath::Pi();
				// 	// phiHX_gen = muPlus_HX_gen.Phi() * 180 / TMath::Pi();
				// }

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
				// dimuWeight_nominal = tnp_weight_trk_pbpb(Reco_mupl_eta, indexNominal) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexNominal) * dimuTrigWeight_nominal;

				dimuWeight_nominal = allGood ? tnp_weight_trk_pbpb(Reco_mupl_eta, indexNominal) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexNominal) * dimuTrigWeight_nominal : 1; // if the event is not selected, we do not apply the dimuon weight to the denominator
				
				// total weight
				totalWeightLab = eventWeight * dimuonPtWeight * dimuWeight_nominal;
				totalWeightCS = eventWeight * dimuonPtWeight * dimuWeight_nominal * weightCS;
				totalWeightHX = eventWeight * dimuonPtWeight * dimuWeight_nominal * weightHX;

				// cout << "eventWeight: " << eventWeight << endl;
				// cout << "dimuonPtWeight: " << dimuonPtWeight << endl;
				// cout << "dimuWeight_nominal: " << dimuWeight_nominal << endl;
				// cout << "weightCS: " << weightCS << endl;
				// cout << "totalWeightCS: " << totalWeightCS << endl;
				hNominalEffLab->FillWeighted(allGood, totalWeightLab, cosThetaLab_gen, phiLab_gen, gen_QQ_pt);
				hNominalEffCS->FillWeighted(allGood, totalWeightCS, cosThetaCS_gen, phiCS_gen, gen_QQ_pt);
				hNominalEffHX->FillWeighted(allGood, totalWeightHX, cosThetaHX_gen, phiHX_gen, gen_QQ_pt);

				// cout << "Event " << iEvent << endl;
				// cout << "all good: " << allGood << endl;
				// cout << "total weight HX: " << totalWeightHX << endl;
				// cout << "cosThetaHX: " << cosThetaHX_gen << endl;
				// cout << "phiHX: " << phiHX_gen << endl;
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

				if (gen_QQ_pt > ptMin && gen_QQ_pt < ptMax) { // pT range of interest
					hEffCS2D->FillWeighted(allGood, totalWeightCS, cosThetaCS_gen, phiCS_gen);
					hEffHX2D->FillWeighted(allGood, totalWeightHX, cosThetaHX_gen, phiHX_gen);

					hEffCS1D->FillWeighted(allGood, totalWeightCS, cosThetaCS_gen);
					hEffHX1D->FillWeighted(allGood, totalWeightHX, cosThetaHX_gen);
				}

				/// variations for muon tracking SF (keeping the nominal efficiency for muon Id and trigger)

				// tracking, syst up
				dimuWeight_trk_systUp = tnp_weight_trk_pbpb(Reco_mupl_eta, indexSystUp) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexSystUp) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexNominal) * dimuTrigWeight_nominal;

				totalWeightCS = eventWeight * dimuonPtWeight * dimuWeight_trk_systUp * weightCS;
				totalWeightHX = eventWeight * dimuonPtWeight * dimuWeight_trk_systUp * weightHX;
				hCS_trk_systUp->FillWeighted(allGood, totalWeightCS, cosThetaCS_gen, phiCS_gen, gen_QQ_pt);
				hHX_trk_systUp->FillWeighted(allGood, totalWeightHX, cosThetaHX_gen, phiHX_gen, gen_QQ_pt);

				// cout << "Event " << iEvent << endl;
				// cout << "systUp trk weight CS: " << totalWeightCS << endl;
				// cout << "systUp trk weight HX: " << totalWeightHX << endl;

				// tracking, syst down
				dimuWeight_trk_systDown = tnp_weight_trk_pbpb(Reco_mupl_eta, indexSystDown) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexSystDown) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexNominal) * dimuTrigWeight_nominal;

				totalWeightCS = eventWeight * dimuonPtWeight * dimuWeight_trk_systDown * weightCS;
				totalWeightHX = eventWeight * dimuonPtWeight * dimuWeight_trk_systDown * weightHX;
				hCS_trk_systDown->FillWeighted(allGood, totalWeightCS, cosThetaCS_gen, phiCS_gen, gen_QQ_pt);
				hHX_trk_systDown->FillWeighted(allGood, totalWeightHX, cosThetaHX_gen, phiHX_gen, gen_QQ_pt);

				// cout << "Event " << iEvent << endl;
				// cout << "systDown trk weight CS: " << totalWeightCS << endl;
				// cout << "systDown trk weight HX: " << totalWeightHX << endl;

				// tracking, stat up
				dimuWeight_trk_statUp = tnp_weight_trk_pbpb(Reco_mupl_eta, indexStatUp) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexStatUp) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexNominal) * dimuTrigWeight_nominal;

				totalWeightCS = eventWeight * dimuonPtWeight * dimuWeight_trk_statUp * weightCS;
				totalWeightHX = eventWeight * dimuonPtWeight * dimuWeight_trk_statUp * weightHX;
				hCS_trk_statUp->FillWeighted(allGood, totalWeightCS, cosThetaCS_gen, phiCS_gen, gen_QQ_pt);
				hHX_trk_statUp->FillWeighted(allGood, totalWeightHX, cosThetaHX_gen, phiHX_gen, gen_QQ_pt);

				// tracking, stat down
				dimuWeight_trk_statDown = tnp_weight_trk_pbpb(Reco_mupl_eta, indexStatDown) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexStatDown) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexNominal) * dimuTrigWeight_nominal;

				totalWeightCS = eventWeight * dimuonPtWeight * dimuWeight_trk_statDown * weightCS;
				totalWeightHX = eventWeight * dimuonPtWeight * dimuWeight_trk_statDown * weightHX;
				hCS_trk_statDown->FillWeighted(allGood, totalWeightCS, cosThetaCS_gen, phiCS_gen, gen_QQ_pt);
				hHX_trk_statDown->FillWeighted(allGood, totalWeightHX, cosThetaHX_gen, phiHX_gen, gen_QQ_pt);

				/// variations for muon Id SF (keeping the nominal efficiency for tracking and trigger)

				// Id, syst up
				dimuWeight_muId_systUp = tnp_weight_trk_pbpb(Reco_mupl_eta, indexNominal) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexSystUp) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexSystUp) * dimuTrigWeight_nominal;

				totalWeightCS = eventWeight * dimuonPtWeight * dimuWeight_muId_systUp * weightCS;
				totalWeightHX = eventWeight * dimuonPtWeight * dimuWeight_muId_systUp * weightHX;
				hCS_muId_systUp->FillWeighted(allGood, totalWeightCS, cosThetaCS_gen, phiCS_gen, gen_QQ_pt);
				hHX_muId_systUp->FillWeighted(allGood, totalWeightHX, cosThetaHX_gen, phiHX_gen, gen_QQ_pt);

				// cout << "Event " << iEvent << endl;
				// cout << "systUp muID weight CS: " << totalWeightCS << endl;
				// cout << "systUp muID weight HX: " << totalWeightHX << endl;

				// Id, syst down
				dimuWeight_muId_systDown = tnp_weight_trk_pbpb(Reco_mupl_eta, indexNominal) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexSystDown) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexSystDown) * dimuTrigWeight_nominal;

				totalWeightCS = eventWeight * dimuonPtWeight * dimuWeight_muId_systDown * weightCS;
				totalWeightHX = eventWeight * dimuonPtWeight * dimuWeight_muId_systDown * weightHX;
				hCS_muId_systDown->FillWeighted(allGood, totalWeightCS, cosThetaCS_gen, phiCS_gen, gen_QQ_pt);
				hHX_muId_systDown->FillWeighted(allGood, totalWeightHX, cosThetaHX_gen, phiHX_gen, gen_QQ_pt);

				// cout << "Event " << iEvent << endl;
				// cout << "systDown muID weight CS: " << totalWeightCS << endl;
				// cout << "systDown muID weight HX: " << totalWeightHX << endl;

				// Id, stat up
				dimuWeight_muId_statUp = tnp_weight_trk_pbpb(Reco_mupl_eta, indexNominal) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexStatUp) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexStatUp) * dimuTrigWeight_nominal;

				totalWeightCS = eventWeight * dimuonPtWeight * dimuWeight_muId_statUp * weightCS;
				totalWeightHX = eventWeight * dimuonPtWeight * dimuWeight_muId_statUp * weightHX;
				hCS_muId_statUp->FillWeighted(allGood, totalWeightCS, cosThetaCS_gen, phiCS_gen, gen_QQ_pt);
				hHX_muId_statUp->FillWeighted(allGood, totalWeightHX, cosThetaHX_gen, phiHX_gen, gen_QQ_pt);

				// Id, stat down
				dimuWeight_muId_statDown = tnp_weight_trk_pbpb(Reco_mupl_eta, indexNominal) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexStatDown) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexStatDown) * dimuTrigWeight_nominal;

				totalWeightCS = eventWeight * dimuonPtWeight * dimuWeight_muId_statDown * weightCS;
				totalWeightHX = eventWeight * dimuonPtWeight * dimuWeight_muId_statDown * weightHX;
				hCS_muId_statDown->FillWeighted(allGood, totalWeightCS, cosThetaCS_gen, phiCS_gen, gen_QQ_pt);
				hHX_muId_statDown->FillWeighted(allGood, totalWeightHX, cosThetaHX_gen, phiHX_gen, gen_QQ_pt);

				/// variations for trigger SF (keeping the nominal efficiency for tracking and muon Id)

				// trigger, syst up
				dimuWeight_trig_systUp = tnp_weight_trk_pbpb(Reco_mupl_eta, indexNominal) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexNominal) * dimuTrigWeight_systUp;

				totalWeightCS = eventWeight * dimuonPtWeight * dimuWeight_trig_systUp * weightCS;
				totalWeightHX = eventWeight * dimuonPtWeight * dimuWeight_trig_systUp * weightHX;
				hCS_trig_systUp->FillWeighted(allGood, totalWeightCS, cosThetaCS_gen, phiCS_gen, gen_QQ_pt);
				hHX_trig_systUp->FillWeighted(allGood, totalWeightHX, cosThetaHX_gen, phiHX_gen, gen_QQ_pt);

				// cout << "Event " << iEvent << endl;
				// cout << "systUp trigger weight CS: " << totalWeightCS << endl;
				// cout << "systUp trigger weight HX: " << totalWeightHX << endl;

				// trigger, syst down
				dimuWeight_trig_systDown = tnp_weight_trk_pbpb(Reco_mupl_eta, indexNominal) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexNominal) * dimuTrigWeight_systDown;

				totalWeightCS = eventWeight * dimuonPtWeight * dimuWeight_trig_systDown * weightCS;
				totalWeightHX = eventWeight * dimuonPtWeight * dimuWeight_trig_systDown * weightHX;
				hCS_trig_systDown->FillWeighted(allGood, totalWeightCS, cosThetaCS_gen, phiCS_gen, gen_QQ_pt);
				hHX_trig_systDown->FillWeighted(allGood, totalWeightHX, cosThetaHX_gen, phiHX_gen, gen_QQ_pt);

				// cout << "Event " << iEvent << endl;
				// cout << "systDown trigger weight CS: " << totalWeightCS << endl;
				// cout << "systDown trigger weight HX: " << totalWeightHX << endl;

				// trigger, stat up
				dimuWeight_trig_statUp = tnp_weight_trk_pbpb(Reco_mupl_eta, indexNominal) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexNominal) * dimuTrigWeight_statUp;

				totalWeightCS = eventWeight * dimuonPtWeight * dimuWeight_trig_statUp * weightCS;
				totalWeightHX = eventWeight * dimuonPtWeight * dimuWeight_trig_statUp * weightHX;
				hCS_trig_statUp->FillWeighted(allGood, totalWeightCS, cosThetaCS_gen, phiCS_gen, gen_QQ_pt);
				hHX_trig_statUp->FillWeighted(allGood, totalWeightHX, cosThetaHX_gen, phiHX_gen, gen_QQ_pt);

				// cout << "Event " << iEvent << endl;
				// cout << "systUp trigger weight CS: " << totalWeightCS << endl;
				// cout << "systUp trigger weight HX: " << totalWeightHX << endl;

				// trigger, stat down
				dimuWeight_trig_statDown = tnp_weight_trk_pbpb(Reco_mupl_eta, indexNominal) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexNominal) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexNominal) * dimuTrigWeight_statDown;

				totalWeightCS = eventWeight * dimuonPtWeight * dimuWeight_trig_statDown * weightCS;
				totalWeightHX = eventWeight * dimuonPtWeight * dimuWeight_trig_statDown * weightHX;
				hCS_trig_statDown->FillWeighted(allGood, totalWeightCS, cosThetaCS_gen, phiCS_gen, gen_QQ_pt);
				hHX_trig_statDown->FillWeighted(allGood, totalWeightHX, cosThetaHX_gen, phiHX_gen, gen_QQ_pt);

				// cout << "Event " << iEvent << endl;
				// cout << "systDown trigger weight CS: " << totalWeightCS << endl;
				// cout << "systDown trigger weight HX: " << totalWeightHX << endl;

				// track + muID + trigger, syst up
				dimuWeight_total_systUp = tnp_weight_trk_pbpb(Reco_mupl_eta, indexSystUp) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexSystUp) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexSystUp) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexSystUp) * dimuTrigWeight_systUp;

				totalWeightCS = eventWeight * dimuonPtWeight * dimuWeight_total_systUp * weightCS;
				totalWeightHX = eventWeight * dimuonPtWeight * dimuWeight_total_systUp * weightHX;
				hCS_total_systUp->FillWeighted(allGood, totalWeightCS, cosThetaCS_gen, phiCS_gen, gen_QQ_pt);
				hHX_total_systUp->FillWeighted(allGood, totalWeightHX, cosThetaHX_gen, phiHX_gen, gen_QQ_pt);

				// cout << "Event " << iEvent << endl;
				// cout << "systUp total weight CS: " << totalWeightCS << endl;
				// cout << "systUp total weight HX: " << totalWeightHX << endl;

				// track + muID + trigger, syst down
				dimuWeight_total_systDown = tnp_weight_trk_pbpb(Reco_mupl_eta, indexSystDown) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexSystDown) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexSystDown) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexSystDown) * dimuTrigWeight_systDown;

				totalWeightCS = eventWeight * dimuonPtWeight * dimuWeight_total_systDown * weightCS;
				totalWeightHX = eventWeight * dimuonPtWeight * dimuWeight_total_systDown * weightHX;
				hCS_total_systDown->FillWeighted(allGood, totalWeightCS, cosThetaCS_gen, phiCS_gen, gen_QQ_pt);
				hHX_total_systDown->FillWeighted(allGood, totalWeightHX, cosThetaHX_gen, phiHX_gen, gen_QQ_pt);

				// cout << "Event " << iEvent << endl;
				// cout << "systDown total weight CS: " << totalWeightCS << endl;
				// cout << "systDown total weight HX: " << totalWeightHX << endl;

				// track + muID + trigger, stat up
				dimuWeight_total_statUp = tnp_weight_trk_pbpb(Reco_mupl_eta, indexStatUp) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexStatUp) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexStatUp) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexStatUp) * dimuTrigWeight_statUp;

				totalWeightCS = eventWeight * dimuonPtWeight * dimuWeight_total_statUp * weightCS;
				totalWeightHX = eventWeight * dimuonPtWeight * dimuWeight_total_statUp * weightHX;
				hCS_total_statUp->FillWeighted(allGood, totalWeightCS, cosThetaCS_gen, phiCS_gen, gen_QQ_pt);
				hHX_total_statUp->FillWeighted(allGood, totalWeightHX, cosThetaHX_gen, phiHX_gen, gen_QQ_pt);

				// track + muID + trigger, stat down
				dimuWeight_total_statDown = tnp_weight_trk_pbpb(Reco_mupl_eta, indexStatDown) * tnp_weight_trk_pbpb(Reco_mumi_eta, indexStatDown) * tnp_weight_muid_pbpb(Reco_mupl_pt, Reco_mupl_eta, indexStatDown) * tnp_weight_muid_pbpb(Reco_mumi_pt, Reco_mumi_eta, indexStatDown) * dimuTrigWeight_statDown;

				totalWeightCS = eventWeight * dimuonPtWeight * dimuWeight_total_statDown * weightCS;
				totalWeightHX = eventWeight * dimuonPtWeight * dimuWeight_total_statDown * weightHX;
				hCS_total_statDown->FillWeighted(allGood, totalWeightCS, cosThetaCS_gen, phiCS_gen, gen_QQ_pt);
				hHX_total_statDown->FillWeighted(allGood, totalWeightHX, cosThetaHX_gen, phiHX_gen, gen_QQ_pt);

				TH3D* hPassed = (TH3D*)hNominalEffCS->GetPassedHistogram();

				Int_t iCosThetaBin = hPassed->GetXaxis()->FindBin(cosThetaCS_gen);
				Int_t iPhiBin = hPassed->GetYaxis()->FindBin(phiCS_gen);
				Int_t iPtBin = hPassed->GetZaxis()->FindBin(gen_QQ_pt);

				int globalBin = hNominalEffCS->GetGlobalBin(iCosThetaBin, iPhiBin, iPtBin);

				// cout << "Event " << iEvent << endl;
				// cout << iCosThetaBin << " " << iPhiBin << " " << iPtBin << " " << globalBin << endl;
				// cout << "CosTheta: " << cosThetaCS << ", " << "Phi: " << phiCS << ", " << "pT: " << reco_QQ_pt << endl;
				// cout << "weightCS: " << totalWeightCS << endl;

				// cout << "Nominal efficiency: " << hNominalEffCS->GetEfficiency(globalBin) << endl;

				// cout << "Trk Syst up efficiency: " << hCS_trk_systUp->GetEfficiency(globalBin) << endl;
				// cout << "Trk Syst down efficiency: " << hCS_trk_systDown->GetEfficiency(globalBin) << endl;
				// cout << "Trk Stat up efficiency: " << hCS_trk_statUp->GetEfficiency(globalBin) << endl;
				// cout << "Trk Stat down efficiency: " << hCS_trk_statDown->GetEfficiency(globalBin) << endl;

				// cout << "Trigger Syst up efficiency: " << hCS_trig_statDown->GetEfficiency(globalBin) << endl;
				// cout << "Trigger Syst down efficiency: " << hCS_trig_statDown->GetEfficiency(globalBin) << endl;
				// cout << "Trigger Stat up efficiency: " << hCS_trig_statUp->GetEfficiency(globalBin) << endl;
				// cout << "Trigger Stat down efficiency: " << hCS_trig_statDown->GetEfficiency(globalBin) << endl;

				// cout << "MuId Syst up efficiency: " << hCS_muId_systUp->GetEfficiency(globalBin) << endl;
				// cout << "MuId Syst down efficiency: " << hCS_muId_systDown->GetEfficiency(globalBin) << endl;
				// cout << "MuId Stat up efficiency: " << hCS_muId_statUp->GetEfficiency(globalBin) << endl;
				// cout << "MuId Stat down efficiency: " << hCS_muId_statDown->GetEfficiency(globalBin) << endl;

				// cout << endl;
			}
			// efficiency denominator
			else {
				dimuWeight_nominal = 1.;

				totalWeightLab = eventWeight * dimuonPtWeight * dimuWeight_nominal;
				totalWeightCS = eventWeight * dimuonPtWeight * dimuWeight_nominal * weightCS;
				totalWeightHX = eventWeight * dimuonPtWeight * dimuWeight_nominal * weightHX;

				allGood = 0;

				hNominalEffLab->FillWeighted(allGood, totalWeightLab, cosThetaLab_gen, phiLab_gen, gen_QQ_pt);
				hNominalEffCS->FillWeighted(allGood, totalWeightCS, cosThetaCS_gen, phiCS_gen, gen_QQ_pt);
				hNominalEffHX->FillWeighted(allGood, totalWeightHX, cosThetaHX_gen, phiHX_gen, gen_QQ_pt);
			}
		} // end of gen upsilon loop
	}

	cout << endl;

	gStyle->SetPadLeftMargin(.15);
	gStyle->SetPadRightMargin(0.18);
	SetColorPalette(gEfficiencyColorPaletteName);

	/// display the nominal results

	DrawEfficiencyMap(hEffCS2D, ptMin, ptMax, muonAccName, iState, isPhiFolded);
	DrawEfficiency1DHist(hEffCS1D, ptMin, ptMax, muonAccName, iState, false, true, "CS", isPhiFolded);

	DrawEfficiencyMap(hEffHX2D, ptMin, ptMax, muonAccName, iState, isPhiFolded);
	DrawEfficiency1DHist(hEffHX1D, ptMin, ptMax, muonAccName, iState, false, true, "HX", isPhiFolded);

	/// compute the systematics in this macro since we have all the ingredients for that
	// instructions can be found here: https://twiki.cern.ch/twiki/pub/CMS/HIMuonTagProbe/TnpHeaderFile.pdf#page=5

	// store the RELATIVE SYSTEMATIC UNCERTAINTIES (not in % though) with respect to the nominal efficiency, more useful at the end
	TString path = EfficiencyResultsPath(muonAccName.Data());

	gSystem->mkdir(path, kTRUE);

	gStyle->SetPadRightMargin(0.19);
	gStyle->SetTitleOffset(1.3, "Z");

	const char* legendText = EfficiencyLegendText(ptMin, ptMax);

	TLatex* legend = new TLatex(.5, .88, legendText);
	legend->SetTextAlign(22);
	legend->SetTextSize(0.042);

	// Collins-Soper
	auto* hSystCS3D = SystEffHist3D(hNominalEffCS, hCS_trk_systUp, hCS_trk_systDown, hCS_muId_systUp, hCS_muId_systDown, hCS_trig_systUp, hCS_trig_systDown, hCS_trk_statUp, hCS_trk_statDown, hCS_muId_statUp, hCS_muId_statDown, hCS_trig_statUp, hCS_trig_statDown, "CS", lambdaTheta, lambdaPhi, lambdaThetaPhi);

	auto* canvasCSsyst = new TCanvas("canvasCSsyst", "", 700, 600);
	hSystCS3D->Draw("COLZ");

	CMS_lumi(canvasCSsyst, Form("#varUpsilon(%dS) Hydjet-embedded MC", iState));

	legend->DrawLatexNDC(.5, .88, legendText);

	gPad->Update();

	hSystCS3D->GetYaxis()->SetRangeUser(-190, 240);

	canvasCSsyst->SaveAs(Form("%s/RelatSystEff_CS%s.png", path.Data(), isPhiFolded ? "" : "_fullPhi"), "RECREATE");

	// HelicityTEfficiency
	auto* hSystHX2D = SystEffHist(hNominalEffHX, hHX_trk_systUp, hHX_trk_systDown, hHX_muId_systUp, hHX_muId_systDown, hHX_trig_systUp, hHX_trig_systDown, hHX_trk_statUp, hHX_trk_statDown, hHX_muId_statUp, hHX_muId_statDown, hHX_trig_statUp, hHX_trig_statDown);

	auto* canvasHXsyst2D = new TCanvas("canvasHXsyst2D", "", 700, 600);
	hSystHX2D->Draw();

	auto* hSystHX3D = SystEffHist3D(hNominalEffHX, hHX_trk_systUp, hHX_trk_systDown, hHX_muId_systUp, hHX_muId_systDown, hHX_trig_systUp, hHX_trig_systDown, hHX_trk_statUp, hHX_trk_statDown, hHX_muId_statUp, hHX_muId_statDown, hHX_trig_statUp, hHX_trig_statDown, "HX", lambdaTheta, lambdaPhi, lambdaThetaPhi);

	auto* canvasHXsyst = new TCanvas("canvasHXsyst", "", 700, 600);
	hSystHX3D->Draw("COLZ");

	CMS_lumi(canvasHXsyst, Form("#varUpsilon(%dS) Hydjet-embedded MC", iState));

	legend->DrawLatexNDC(.5, .88, legendText);

	gPad->Update();

	hSystHX3D->GetYaxis()->SetRangeUser(-190, 240);

	canvasHXsyst->SaveAs(Form("%s/RelatSystEff_HX%s.png", path.Data(), isPhiFolded ? "" : "_fullPhi"), "RECREATE");

	/// save the nominal efficiency results and the corresponding systematics in a file for later usage
	TString outputFileName = Form("%s/EfficiencyResults%s.root", path.Data(), isPhiFolded ? "" : "_fullPhi");

	TFile outputFile(outputFileName.Data(), "UPDATE");

	hNominalEffLab->Write();
	hNominalEffCS->Write();
	hNominalEffHX->Write();

	hSystCS3D->Write();
	hSystHX3D->Write();

	hCS_trk_systUp->Write();
	hCS_trk_systDown->Write();
	hCS_muId_systUp->Write();
	hCS_muId_systDown->Write();
	hCS_trig_systUp->Write();
	hCS_trig_systDown->Write();
	hCS_total_systUp->Write();
	hCS_total_systDown->Write();

	hCS_trk_statUp->Write();
	hCS_trk_statDown->Write();
	hCS_muId_statUp->Write();
	hCS_muId_statDown->Write();
	hCS_trig_statUp->Write();
	hCS_trig_statDown->Write();
	hCS_total_statUp->Write();
	hCS_total_statDown->Write();

	hHX_trk_systUp->Write();
	hHX_trk_systDown->Write();
	hHX_muId_systUp->Write();
	hHX_muId_systDown->Write();
	hHX_trig_systUp->Write();
	hHX_trig_systDown->Write();
	hHX_total_systUp->Write();
	hHX_total_systDown->Write();

	hHX_trk_statUp->Write();
	hHX_trk_statDown->Write();
	hHX_muId_statUp->Write();
	hHX_muId_statDown->Write();
	hHX_trig_statUp->Write();
	hHX_trig_statDown->Write();
	hHX_total_statUp->Write();
	hHX_total_statDown->Write();

	hSystHX2D->Write();

	// hEffCS2D->Write();
	// hEffHX2D->Write();

	// hEffCS1D->Write();
	// hEffHX1D->Write();

	outputFile.Close();

	file->Close();

	if (BeVerbose) cout << "\nNominal efficiency and corresponding systematic uncertainty maps saved in " << outputFileName << endl;
}
