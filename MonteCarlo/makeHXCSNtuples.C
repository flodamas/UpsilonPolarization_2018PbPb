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

void makeHXCSNtuples(Int_t iState = 1, Int_t ptMin = 0, Int_t ptMax = 30) {

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

	OniaTree->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size);
	OniaTree->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom);
	OniaTree->SetBranchAddress("Gen_QQ_mumi_4mom", &Gen_QQ_mumi_4mom);
	OniaTree->SetBranchAddress("Gen_QQ_mupl_4mom", &Gen_QQ_mupl_4mom);

	TLorentzVector* gen_QQ_LV = new TLorentzVector();
	TLorentzVector* gen_mumi_LV = new TLorentzVector();
	TLorentzVector* gen_mupl_LV = new TLorentzVector();

	Long64_t totEntries = OniaTree->GetEntries();

	OniaTree -> SetBranchAddress("Gen_QQ_size", &Gen_QQ_size);

	// Create a Ntuple to store kinematics of Upsilon and daughter muons
	gROOT->cd();
	TString varlist = "upsM:upsRapLab:upsPtLab:muplPtLab:muplEtaLab:mumiPtLab:mumiEtaLab:muplCosThetaHX:muplPhiHX:muplCosThetaCS:muplPhiCS";
	TNtuple* HXCSNTuple = new TNtuple("HXCSNTuple", "Upsilon and muons in the HX and CS ntuple", varlist);

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

			// ******** Fill Ntuple with kinematics of upsilon and muons in the Lab, HX, and CS frames******** //
			float tuple[] = {
				static_cast<float>(gen_QQ_LV->M()),
			    static_cast<float>(gen_QQ_LV->Rapidity()),
			    static_cast<float>(gen_QQ_LV->Pt()),

			    static_cast<float>(gen_mupl_LV->Pt()),
			    static_cast<float>(gen_mupl_LV->Eta()),
			    static_cast<float>(gen_mumi_LV->Pt()),
			    static_cast<float>(gen_mumi_LV->Eta()),

			    static_cast<float>(muPlus_HX.CosTheta()),
			    static_cast<float>((muPlus_HX.Phi())*180/TMath::Pi()),
			    static_cast<float>(muPlus_CS.CosTheta()),
			    static_cast<float>((muPlus_CS.Phi())*180/TMath::Pi())
			};
			HXCSNTuple -> Fill(tuple);
		}
	}

	// Create a file and store the ntuple
	// const char* outfilename = Form("../Files/OniaTree_Y%dS_GENONLY_NoFilter_HXCS.root", iState);
	// TFile* outfile = new TFile(outfilename, "RECREATE", "Y 1S in HX and CS");  	
  	file->cd();
  	HXCSNTuple -> Write();
  	file -> Close();

	// End measuring time
	end = clock();
	cpu_time = (double)(end - start) / CLOCKS_PER_SEC;
	cout << "Time elapsed: " << cpu_time / 60. << "minutes" << endl;


}