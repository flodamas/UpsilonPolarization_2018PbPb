#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/Parameters/PhaseSpace.h"

#include "AccEffHelpers.h"

#include "../Polarization/PolarFitHelpers.h"

//#include "../Tools/Style/FitDistributions.h"
#include "../Tools/Style/Legends.h"

#include "../ReferenceFrameTransform/Transformations.h"

void acceptanceMap_HydjetGen(Int_t ptMin = 0, Int_t ptMax = 30, Bool_t isPhiFolded = kTRUE, TString muonAccName = gMuonAccName, Float_t lambdaTheta = 0, Float_t lambdaPhi = 0, Float_t lambdaThetaPhi = 0, Int_t iState = 1) {
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

	// loop variables
	TLorentzVector* genLorentzVector = new TLorentzVector();
	TLorentzVector* gen_QQ_LV = new TLorentzVector();
	TLorentzVector* gen_mupl_LV = new TLorentzVector();
	TLorentzVector* gen_mumi_LV = new TLorentzVector();
	TLorentzVector* recoLorentzVector = new TLorentzVector();

	// (cos theta, phi, pT) 3D maps for final acceptance correction, variable size binning for the stats
	TEfficiency* accMatrixLab = TEfficiency3D(NominalTEfficiency3DName("Lab", lambdaTheta, lambdaPhi, lambdaThetaPhi), "Lab", iState, isPhiFolded);
	TEfficiency* accMatrixCS = TEfficiency3D(NominalTEfficiency3DName("CS", lambdaTheta, lambdaPhi, lambdaThetaPhi), "CS", iState, isPhiFolded);
	TEfficiency* accMatrixHX = TEfficiency3D(NominalTEfficiency3DName("HX", lambdaTheta, lambdaPhi, lambdaThetaPhi), "HX", iState, isPhiFolded);

    Bool_t withinAcceptance;

	Double_t cosThetaCS, phiCS, cosThetaHX, phiHX, cosThetaLab, phiLab;

	Float_t weightCS = 0, weightHX = 0;

	Long64_t totEntries = OniaTree->GetEntries();
	double counter = 0;

	// Loop over the events
	for (Long64_t iEvent = 0; iEvent < (totEntries); iEvent++) {
		if (iEvent % 10000 == 0) {
			cout << Form("\rProcessing event %lld / %lld (%.0f%%)", iEvent, totEntries, 100. * iEvent / totEntries) << flush;
		}

        OniaTree->GetEntry(iEvent);

		if (Centrality >= 2 * gCentralityBinMax) continue; // discard events with centrality >= 90% in 2018 data

		// firesTrigger = ((HLTriggers & (ULong64_t)(1 << (gUpsilonHLTBit - 1))) == (ULong64_t)(1 << (gUpsilonHLTBit - 1)));

		// hiBin = GetHiBinFromhiHF(HFmean);

        // loop over all gen upsilons
		for (int iGen = 0; iGen < Gen_QQ_size; iGen++) {

            gen_QQ_LV = (TLorentzVector*)Gen_QQ_4mom->At(iGen);
   
            if (fabs(gen_QQ_LV->Rapidity()) < gRapidityMin || fabs(gen_QQ_LV->Rapidity()) > gRapidityMax) continue; // upsilon within fiducial region

           	// single-muon acceptance cuts
            gen_mupl_LV = (TLorentzVector*)Gen_mu_4mom->At(Gen_QQ_mupl_idx[iGen]);
            gen_mumi_LV = (TLorentzVector*)Gen_mu_4mom->At(Gen_QQ_mumi_idx[iGen]);

			withinAcceptance = MuonKinematicsWithinLimits(*gen_mupl_LV, muonAccName) && MuonKinematicsWithinLimits(*gen_mumi_LV, muonAccName);

            TVector3 muPlus_Lab = gen_mupl_LV->Vect();

			cosThetaLab = muPlus_Lab.CosTheta();
			phiLab = muPlus_Lab.Phi() * 180 / TMath::Pi();

			if (isPhiFolded == kTRUE) phiLab = fabs(phiLab);

			accMatrixLab->FillWeighted(withinAcceptance, 1, cosThetaLab, phiLab, gen_QQ_LV->Pt());

			// Reference frame transformations
			TVector3 muPlus_CS = MuPlusVector_CollinsSoper(*gen_QQ_LV, *gen_mupl_LV);

			cosThetaCS = muPlus_CS.CosTheta();

			phiCS = muPlus_CS.Phi() * 180 / TMath::Pi();
			if (isPhiFolded == kTRUE) phiCS = fabs(phiCS);
			// cout << "cosThetaCS: " << cosThetaCS << endl;
			// cout << "phiCS: " << phiCS << endl;

			if (isPhiFolded == kTRUE)
				weightCS = 1 + lambdaTheta * TMath::Power(muPlus_CS.CosTheta(), 2) + lambdaPhi * TMath::Power(std::sin(muPlus_CS.Theta()), 2) * std::cos(2 * fabs(muPlus_CS.Phi())) + lambdaThetaPhi * std::sin(2 * muPlus_CS.Theta()) * std::cos(fabs(muPlus_CS.Phi()));
			else
				weightCS = 1 + lambdaTheta * TMath::Power(muPlus_CS.CosTheta(), 2) + lambdaPhi * TMath::Power(std::sin(muPlus_CS.Theta()), 2) * std::cos(2 * muPlus_CS.Phi()) + lambdaThetaPhi * std::sin(2 * muPlus_CS.Theta()) * std::cos(muPlus_CS.Phi());

			// cout << "weightCS: " << weightCS << endl;
			// cout << "cosThetaCS: " << cosThetaCS << endl;
			// cout << "phiCS: " << phiCS << endl;

			accMatrixCS->FillWeighted(withinAcceptance, weightCS, cosThetaCS, phiCS, gen_QQ_LV->Pt());

			TVector3 muPlus_HX = MuPlusVector_Helicity(*gen_QQ_LV, *gen_mupl_LV);

			cosThetaHX = muPlus_HX.CosTheta();
			phiHX = muPlus_HX.Phi() * 180 / TMath::Pi();
			if (isPhiFolded == kTRUE)
				phiHX = fabs(phiHX);

			if (isPhiFolded == kTRUE)
				weightHX = 1 + lambdaTheta * TMath::Power(muPlus_HX.CosTheta(), 2) + lambdaPhi * TMath::Power(std::sin(muPlus_HX.Theta()), 2) * std::cos(2 * fabs(muPlus_HX.Phi())) + lambdaThetaPhi * std::sin(2 * muPlus_HX.Theta()) * std::cos(fabs(muPlus_HX.Phi()));
			else
				weightHX = 1 + lambdaTheta * TMath::Power(muPlus_HX.CosTheta(), 2) + lambdaPhi * TMath::Power(std::sin(muPlus_HX.Theta()), 2) * std::cos(2 * muPlus_HX.Phi()) + lambdaThetaPhi * std::sin(2 * muPlus_HX.Theta()) * std::cos(muPlus_HX.Phi());

			// cout << "weightHX: " << weightHX << endl;
			// cout << "cosThetaHX: " << cosThetaHX << endl;
			// cout << "phiHX: " << phiHX << endl;

			accMatrixHX->FillWeighted(withinAcceptance, weightHX, cosThetaHX, phiHX, gen_QQ_LV->Pt());
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
    const char* outputFileName = Form("%s/AcceptanceResults_Hydjet%s.root", path, isPhiFolded ? "" : "_fullPhi");
 
    TFile outputFile(outputFileName, "UPDATE");
 
    //accMatrixLab->Write();
    accMatrixCS->Write();
    accMatrixHX->Write();
     
    outputFile.Write();
    outputFile.Close();
 
    if (BeVerbose) cout << "\nAcceptance maps saved in " << outputFileName << endl;     
}

TEfficiency* getAcceptance3DMap(const char* refFrameName, Double_t lambdaTheta = 0, Double_t lambdaPhi = 0, Double_t lambdaThetaPhi = 0, Bool_t isPhiFolded = kFALSE) {
	/// get acceptance and efficiency in 1D
	TString nominalMapName = NominalTEfficiency3DName(refFrameName, lambdaTheta, lambdaPhi, lambdaThetaPhi);

	TString fileName = Form("%s/AcceptanceResults_Hydjet%s.root", AcceptanceResultsPath(gMuonAccName), isPhiFolded ? "" : "_fullPhi");

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

		acceptanceMap_HydjetGen(0, 30, isPhiFolded, "UpsilonTriggerThresholds", lambdaTheta, lambdaPhi, lambdaThetaPhi, gUpsilonState);

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

    accMap = getAcceptance3DMap(refFrameName.Data(), lambdaTheta, lambdaPhi, lambdaThetaPhi, isPhiFolded);

    TEfficiency* accMapCosThetaPhi = rebinTEff3DMap(accMap, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges);

    TH2D* hTotalCosThetaPhi = (TH2D*)accMapCosThetaPhi->GetTotalHistogram();
    // hTotalCosThetaPhi->Scale(1. / hTotalCosThetaPhi->Integral(0, nCosThetaBins + 1, 0, nPhiBins + 1));
    
	TH2D* hPassedCosThetaPhi = (TH2D*)accMapCosThetaPhi->GetPassedHistogram();
    
    /// draw acceptance graph for check
	TCanvas* accCanvas = nullptr;
	TCanvas* accTotalCanvas = nullptr;
    TCanvas* accPassedCanvas = nullptr;

    accCanvas = DrawEfficiency2DHist(accMapCosThetaPhi, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, gUpsilonState, kTRUE, kTRUE, kFALSE, "_TriggerAcc", isPhiFolded, kTRUE, lambdaTheta, lambdaPhi, lambdaThetaPhi);
   
    accTotalCanvas = draw2DMap(hTotalCosThetaPhi, refFrameName.Data(), nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kFALSE, 1, isPhiFolded);
    display2DMapContents(hTotalCosThetaPhi, nCosThetaBins, nPhiBins, kFALSE);

    accPassedCanvas = draw2DMap(hPassedCosThetaPhi, refFrameName.Data(), nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kFALSE, 1, isPhiFolded);
    display2DMapContents(hPassedCosThetaPhi, nCosThetaBins, nPhiBins, kFALSE);
}