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

// skim and keep all MC events, but add a category "accepted" to know if a given dimuon is within our acceptance

void skimAcceptedMC(Int_t iState = 1) {
	const char* inputFileName = Form("OniaTree_Y%dS_GENONLY_NoFilter.root", iState);

	const char* outputFileName = Form("Y%dSGeneratedMCDataset.root", iState);

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
	TFile file(outputFileName, "RECREATE");

	RooCategory eventCat("eventCat", "is the dimuon within acceptance?"); // for efficiency PDF
	eventCat.defineType("accepted", 1);
	eventCat.defineType("rejected", 0);

	RooRealVar yVar("rapidity", "dimuon absolute rapidity", gRapidityMin, gRapidityMax);
	RooRealVar ptVar("pt", "dimuon pT", 0, 100, "GeV/c");

	RooRealVar etaMuPlusVar("etaMuPlus", "", 0, 100);
	RooRealVar ptMuPlusVar("ptMuPlus", "positive muon pT", 0, 100, "GeV/c");

	RooRealVar etaMuMinusVar("etaMuMinus", "", 0, 100);
	RooRealVar ptMuMinusVar("ptMuMinus", "negative muon pT", 0, 100, "GeV/c");

	RooRealVar cosThetaCSVar("cosThetaCS", "cos theta in the Collins-Soper frame", -1, 1);
	RooRealVar phiCSVar("phiCS", "phi angle in the Collins-Soper frame", -180, 180, "#circ");

	RooRealVar cosThetaHXVar("cosThetaHX", "cos theta in the helicity frame", -1, 1);
	RooRealVar phiHXVar("phiHX", "phi angle in the helicity frame", -180, 180, "#circ");

	RooDataSet dataset("MCdataset", "skimmed MC dataset", RooArgSet(eventCat, yVar, ptVar, etaMuPlusVar, ptMuPlusVar, etaMuMinusVar, ptMuMinusVar, cosThetaCSVar, phiCSVar, cosThetaHXVar, phiHXVar));

	// loop variables

	TLorentzVector* gen_QQ_LV = new TLorentzVector();
	TLorentzVector* gen_mumi_LV = new TLorentzVector();
	TLorentzVector* gen_mupl_LV = new TLorentzVector();

	Long64_t totEntries = OniaTree->GetEntries();

	bool isAccepted = false;

	for (Long64_t iEvent = 0; iEvent < (totEntries); iEvent++) {
		if (iEvent % 10000 == 0) {
			cout << Form("\rProcessing event %lld / %lld (%.0f%%)", iEvent, totEntries, 100. * iEvent / totEntries) << flush;
		}

		OniaTree->GetEntry(iEvent);

		// loop over gen
		for (Int_t iGen = 0; iGen < Gen_QQ_size; iGen++) {
			gen_QQ_LV = (TLorentzVector*)Gen_QQ_4mom->At(iGen);

			// fiducial region

			if (fabs(gen_QQ_LV->Rapidity()) < gRapidityMin || fabs(gen_QQ_LV->Rapidity()) > gRapidityMax) continue;

			gen_mupl_LV = (TLorentzVector*)Gen_QQ_mupl_4mom->At(iGen);

			gen_mumi_LV = (TLorentzVector*)Gen_QQ_mumi_4mom->At(iGen);

			// get positive muon's coordinates in the studied reference frames
			TVector3 muPlus_CS = MuPlusVector_CollinsSoper(*gen_QQ_LV, *gen_mupl_LV);

			TVector3 muPlus_HX = MuPlusVector_Helicity(*gen_QQ_LV, *gen_mupl_LV);

			// fill the dataset

			yVar = fabs(gen_QQ_LV->Rapidity());
			ptVar = gen_QQ_LV->Pt();

			ptMuPlusVar = gen_mupl_LV->Pt();
			etaMuPlusVar = fabs(gen_mupl_LV->Eta());

			ptMuMinusVar = gen_mumi_LV->Pt();
			etaMuMinusVar = fabs(gen_mumi_LV->Eta());

			cosThetaCSVar = muPlus_CS.CosTheta();
			phiCSVar = muPlus_CS.Phi() * 180 / TMath::Pi();

			cosThetaHXVar = muPlus_HX.CosTheta();
			phiHXVar = muPlus_HX.Phi() * 180 / TMath::Pi();

			isAccepted = (etaMuPlusVar.getVal() < 2.4 && ptMuPlusVar.getVal() > 3.5) && (etaMuMinusVar.getVal() < 2.4 && ptMuMinusVar.getVal() > 3.5);

			eventCat.setLabel((isAccepted) ? "accepted" : "rejected");

			dataset.add(RooArgSet(eventCat, yVar, ptVar, ptMuPlusVar, etaMuPlusVar, ptMuMinusVar, etaMuMinusVar, cosThetaCSVar, phiCSVar, cosThetaHXVar, phiHXVar));
		}
	}

	dataset.Write();

	file.Close();

	// check the dataset distributions
}
