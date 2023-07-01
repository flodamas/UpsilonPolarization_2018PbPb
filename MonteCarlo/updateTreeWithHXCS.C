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

#include "../Tools/Parameters/AnalysisParameters.h"

#include "../ReferenceFrameTransform/Transformations.h"

void updateTreeWithHXCS(Int_t iState = 1, Int_t ptMin = 0, Int_t ptMax = 30) {

	// Start measuring time
	clock_t start, end, cpu_time;
	start = clock();

	// Read the MC file
	const char* filename = Form("../Files/OniaTree_Y%dS_GENONLY_NoFilter.root", iState);
	TFile* file = TFile::Open(filename, "Update");
	if (!file) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
		return;
	}

	cout << "File " << filename << " opened" << endl;
	TTree* OniaTree = (TTree*)file->Get("hionia/myTree");


	Int_t Gen_QQ_size;
	TClonesArray* Gen_QQ_4mom = nullptr;
	TClonesArray* Gen_QQ_mumi_4mom = nullptr;
	TClonesArray* Gen_QQ_mupl_4mom = nullptr;
	Float_t muplCosThetaHX;
	Float_t muplPhiHX;
	Float_t muplCosThetaCS;
	Float_t muplPhiCS;


	OniaTree->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size);
	OniaTree->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom);
	OniaTree->SetBranchAddress("Gen_QQ_mumi_4mom", &Gen_QQ_mumi_4mom);
	OniaTree->SetBranchAddress("Gen_QQ_mupl_4mom", &Gen_QQ_mupl_4mom);

	TLorentzVector* gen_QQ_LV = new TLorentzVector();
	TLorentzVector* gen_mumi_LV = new TLorentzVector();
	TLorentzVector* gen_mupl_LV = new TLorentzVector();

	OniaTree->Branch("muplCosThetaHX", &muplCosThetaHX, "muplCosThetaHX/F");
	OniaTree->Branch("muplPhiHX", &muplPhiHX, "muplPhiHX/F");
	OniaTree->Branch("muplCosThetaCS", &muplCosThetaCS, "muplCosThetaCS/F");


	OniaTree -> SetBranchStatus("*", 0);
	OniaTree -> SetBranchStatus("Gen_QQ_size", 1);
	OniaTree -> SetBranchStatus("Gen_QQ_4mom", 1);
	OniaTree -> SetBranchStatus("Gen_QQ_mumi_4mom", 1);
	OniaTree -> SetBranchStatus("Gen_QQ_mupl_4mom", 1);

	Long64_t totEntries = OniaTree->GetEntries();

	// Start the event loop for the reference frame transformations and fill a Ntuple 
	for (Long64_t iEvent = 0; iEvent < (totEntries); iEvent++) {
		if (iEvent % 10000 == 0) {
			cout << Form("\rProcessing event %lld / %lld (%.0f%%)", iEvent, totEntries, 100. * iEvent / totEntries) << flush;
		}

		OniaTree->GetEntry(iEvent);

		// loop over all gen upsilons
		for (int iGen = 0; iGen < Gen_QQ_size; iGen++) {
			gen_QQ_LV = (TLorentzVector*)Gen_QQ_4mom->At(iGen);

			// kinematic region of interest
			if (gen_QQ_LV->Pt() < ptMin || gen_QQ_LV->Pt() > ptMax) continue;

			if (fabs(gen_QQ_LV->Rapidity()) > 2.4) continue;

			gen_mupl_LV = (TLorentzVector*)Gen_QQ_mupl_4mom->At(iGen);
			gen_mumi_LV = (TLorentzVector*)Gen_QQ_mumi_4mom->At(iGen);

			// get positive muon's coordinates in the studied reference frames
			TVector3 muPlus_HX = MuPlusVector_Helicity(*gen_QQ_LV, *gen_mupl_LV);
			TVector3 muPlus_CS = MuPlusVector_CollinsSoper(*gen_QQ_LV, *gen_mupl_LV);
			
			muplCosThetaHX = muPlus_HX.CosTheta();
			muplPhiHX = (muPlus_HX.Phi())*180/TMath::Pi();
			muplCosThetaCS = muPlus_CS.CosTheta();
			muplPhiCS = (muPlus_CS.Phi())*180/TMath::Pi();
		}

		OniaTree -> Fill();
	}

	// Create a file and store the ntuple	
  	OniaTree -> Write();
  	file -> Close();

	// End measuring time
	end = clock();
	cpu_time = (double)(end - start) / CLOCKS_PER_SEC;
	cout << "Time elapsed: " << cpu_time / 60. << "minutes" << endl;


}