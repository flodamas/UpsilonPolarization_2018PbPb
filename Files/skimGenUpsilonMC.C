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

#include "../AnalysisParameters.h"

#include "../Tools/Parameters/CentralityValues.h"

#include "../ReferenceFrameTransform/Transformations.h"

void skimGenUpsilonMC(const char* inputFileName = "OniaTree_Y1S_GENONLY_NoFilter.root", const char* outputFileName = "TransverseCSPolarizationGen.root") {
	TFile* infile = TFile::Open(inputFileName, "READ");
	TTree* OniaTree = (TTree*)infile->Get("hionia/myTree");

	/// OniaTree variables, quite old version
	Int_t Gen_QQ_size;

	TClonesArray* Gen_QQ_4mom = nullptr;
	TClonesArray* Gen_QQ_mumi_4mom = nullptr;
	TClonesArray* Gen_QQ_mupl_4mom = nullptr;

	OniaTree->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size);
	OniaTree->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom);

	OniaTree->SetBranchAddress("Gen_QQ_mumi_4mom", &Gen_QQ_mumi_4mom);
	OniaTree->SetBranchAddress("Gen_QQ_mupl_4mom", &Gen_QQ_mupl_4mom);

	/// RooDataSet output: one entry = one dimuon candidate!

	// weighting by event directly on the fly
	RooRealVar eventWeightVar("eventWeight", "polarization weight", 0, 10000);

	RooRealVar yVar("rapidity", "dimuon absolute rapidity", gRapidityMin, gRapidityMax);
	RooRealVar ptVar("pt", "dimuon pT", 0, 100, "GeV/c");

	RooRealVar etaMuPlusVar("etaMuPlus", "", 0, 100);
	RooRealVar ptMuPlusVar("ptMuPlus", "positive muon pT", 0, 100, "GeV/c");

	RooRealVar etaMuMinusVar("etaMuMinus", "", 0, 100);
	RooRealVar ptMuMinusVar("ptMuMinus", "negative muon pT", 0, 100, "GeV/c");

	RooRealVar cosThetaLabVar("cosThetaLab", "cos #theta_{Lab}", -1, 1);
	RooRealVar phiLabVar("phiLab", "#varphi_{Lab}", -180, 180, "#circ");

	RooRealVar cosThetaCSVar("cosThetaCS", "cos #theta_{CS}", -1, 1);
	RooRealVar phiCSVar("phiCS", "#varphi_{CS}", -180, 180, "#circ");

	RooRealVar cosThetaHXVar("cosThetaHX", "cos #theta_{HX}", -1, 1);
	RooRealVar phiHXVar("phiHX", "#varphi_{HX}", -180, 180, "#circ");

	RooDataSet dataset("MCdataset", "skimmed MC dataset", RooArgSet(eventWeightVar, yVar, ptVar, etaMuPlusVar, ptMuPlusVar, etaMuMinusVar, ptMuMinusVar, cosThetaLabVar, phiLabVar, cosThetaCSVar, phiCSVar, cosThetaHXVar, phiHXVar), RooFit::WeightVar("eventWeight"));

	// loop variables
	Float_t weight = 0;

	TLorentzVector* gen_QQ_LV = new TLorentzVector();
	TLorentzVector* gen_mumi_LV = new TLorentzVector();
	TLorentzVector* gen_mupl_LV = new TLorentzVector();

	Long64_t totEntries = OniaTree->GetEntries();

	for (Long64_t iEvent = 0; iEvent < (totEntries); iEvent++) {
		if (iEvent % 10000 == 0) {
			cout << Form("\rProcessing event %lld / %lld (%.0f%%)", iEvent, totEntries, 100. * iEvent / totEntries) << flush;
		}

		OniaTree->GetEntry(iEvent);

		// loop over gen
		for (int iGen = 0; iGen < Gen_QQ_size; iGen++) {
			gen_QQ_LV = (TLorentzVector*)Gen_QQ_4mom->At(iGen);

			gen_mupl_LV = (TLorentzVector*)Gen_QQ_mupl_4mom->At(iGen);

			gen_mumi_LV = (TLorentzVector*)Gen_QQ_mumi_4mom->At(iGen);

			// get positive muon's coordinates in the studied reference frames
			TVector3 muPlus_CS = MuPlusVector_CollinsSoper(*gen_QQ_LV, *gen_mupl_LV);

			TVector3 muPlus_HX = MuPlusVector_Helicity(*gen_QQ_LV, *gen_mupl_LV);

			weight = 1 + TMath::Power(muPlus_CS.CosTheta(), 2);

			// fill the dataset

			eventWeightVar = weight;
			yVar = fabs(gen_QQ_LV->Rapidity());
			ptVar = gen_QQ_LV->Pt();

			ptMuPlusVar = gen_mupl_LV->Pt();
			etaMuPlusVar = fabs(gen_mupl_LV->Eta());

			ptMuMinusVar = gen_mumi_LV->Pt();
			etaMuMinusVar = fabs(gen_mumi_LV->Eta());

			cosThetaLabVar = gen_mupl_LV->CosTheta();
			phiLabVar = gen_mupl_LV->Phi() * 180 / TMath::Pi();

			cosThetaCSVar = muPlus_CS.CosTheta();
			phiCSVar = muPlus_CS.Phi() * 180 / TMath::Pi();

			cosThetaHXVar = muPlus_HX.CosTheta();
			phiHXVar = muPlus_HX.Phi() * 180 / TMath::Pi();

			dataset.add(RooArgSet(eventWeightVar, yVar, ptVar, ptMuPlusVar, etaMuPlusVar, ptMuMinusVar, etaMuMinusVar, cosThetaLabVar, phiLabVar, cosThetaCSVar, phiCSVar, cosThetaHXVar, phiHXVar), weight);
		}
	}

	TFile file(outputFileName, "RECREATE");

	dataset.Write();

	file.Close();

	// check the dataset distributions
}
