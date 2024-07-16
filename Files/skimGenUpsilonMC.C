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

#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/Parameters/CentralityValues.h"

#include "../ReferenceFrameTransform/Transformations.h"

void skimGenUpsilonMC(const char* inputFileName = "OniaTree_Y1S_GENONLY_NoFilter.root", Double_t lambdaTheta = 0, Double_t lambdaPhi = 0, Double_t lambdaThetaPhi = 0) {
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
	RooRealVar eventWeightCSVar("eventWeightCS", "polarization weight CS", 0, 10000);
	RooRealVar eventWeightHXVar("eventWeightHX", "polarization weight HX", 0, 10000);

	RooRealVar yVar("rapidity", gDimuonRapidityVarTitle, 0, 2.4);
	RooRealVar ptVar("pt", gDimuonPtVarTitle, 0, 100, gPtUnit);

	RooRealVar etaMuPlusVar("etaMuPlus", "", 0, 100);
	RooRealVar ptMuPlusVar("ptMuPlus", "positive muon pT", 0, 100, gPtUnit);

	RooRealVar etaMuMinusVar("etaMuMinus", "", 0, 100);
	RooRealVar ptMuMinusVar("ptMuMinus", "negative muon pT", 0, 100, gPtUnit);

	const char* refFrameNameLab = "Lab";
	RooRealVar cosThetaLabVar(CosThetaVarName(refFrameNameLab), CosThetaVarTitle(refFrameNameLab), -1, 1);
	RooRealVar phiLabVar(PhiVarName(refFrameNameLab), PhiVarTitle(refFrameNameLab), -180, 180, gPhiUnit);

	const char* refFrameNameCS = "CS";
	RooRealVar cosThetaCSVar(CosThetaVarName(refFrameNameCS), CosThetaVarTitle(refFrameNameCS), -1, 1);
	RooRealVar phiCSVar(PhiVarName(refFrameNameCS), PhiVarTitle(refFrameNameCS), -180, 180, gPhiUnit);
	RooRealVar phiTildeCSVar(PhiTildeVarName(refFrameNameCS), PhiVarTitle(refFrameNameCS), -180, 180, gPhiUnit);

	const char* refFrameNameHX = "HX";
	RooRealVar cosThetaHXVar(CosThetaVarName(refFrameNameHX), CosThetaVarTitle(refFrameNameHX), -1, 1);
	RooRealVar phiHXVar(PhiVarName(refFrameNameHX), PhiVarTitle(refFrameNameHX), -180, 180, gPhiUnit);
	RooRealVar phiTildeHXVar(PhiTildeVarName(refFrameNameHX), PhiTildeVarTitle(refFrameNameHX), -180, 180, gPhiUnit);

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

	RooDataSet datasetCS("MCdatasetCS", "skimmed MC dataset in CS", RooArgSet(eventWeightCSVar, yVar, ptVar, etaMuPlusVar, ptMuPlusVar, etaMuMinusVar, ptMuMinusVar, cosThetaLabVar, phiLabVar, cosThetaCSVar, phiCSVar, phiTildeCSVar, cosThetaHXVar, phiHXVar, phiTildeHXVar), RooFit::WeightVar("eventWeightCS"));
	RooDataSet datasetHX("MCdatasetHX", "skimmed MC dataset in HX", RooArgSet(eventWeightHXVar, yVar, ptVar, etaMuPlusVar, ptMuPlusVar, etaMuMinusVar, ptMuMinusVar, cosThetaLabVar, phiLabVar, cosThetaCSVar, phiCSVar, phiTildeCSVar, cosThetaHXVar, phiHXVar, phiTildeHXVar), RooFit::WeightVar("eventWeightHX"));

	// loop variables
	Float_t weightCS = 0;
	Float_t weightHX = 0;

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

			weightCS = 1 + lambdaTheta * TMath::Power(muPlus_CS.CosTheta(), 2) + lambdaPhi * TMath::Power(std::sin(muPlus_CS.Theta()), 2) * std::cos(2 * muPlus_CS.Phi()) + lambdaThetaPhi * std::sin(2 * muPlus_CS.Theta()) * std::cos(muPlus_CS.Phi());
			weightHX = 1 + lambdaTheta * TMath::Power(muPlus_HX.CosTheta(), 2) + lambdaPhi * TMath::Power(std::sin(muPlus_HX.Theta()), 2) * std::cos(2 * muPlus_HX.Phi()) + lambdaThetaPhi * std::sin(2 * muPlus_HX.Theta()) * std::cos(muPlus_HX.Phi());;

			// fill the dataset

			eventWeightCSVar = weightCS;
			eventWeightHXVar = weightHX;

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
			
			if (cosThetaCSVar.getVal() < 0) {
				// if phi value is smaller than -pi, add 2pi
				if ((phiCSVar.getVal() - 135) < -180) phiTildeCSVar.setVal(phiCSVar.getVal() + 225);
				else phiTildeCSVar.setVal(phiCSVar.getVal() - 135);
			}

			else if (cosThetaCSVar.getVal() > 0) {
				// if phi value is smaller than -pi, add 2pi
				if ((phiCSVar.getVal() - 45) < -180) phiTildeCSVar.setVal(phiCSVar.getVal() + 315);
				else phiTildeCSVar.setVal(phiCSVar.getVal() - 45);
			}

			cosThetaHXVar = muPlus_HX.CosTheta();
			phiHXVar = muPlus_HX.Phi() * 180 / TMath::Pi();
			
			if (cosThetaHXVar.getVal() < 0) {
				// if phi value is smaller than -pi, add 2pi
				if ((phiHXVar.getVal() - 135) < -180) phiTildeHXVar.setVal(phiHXVar.getVal() + 225);
				else phiTildeHXVar.setVal(phiHXVar.getVal() - 135);
			}

			else if (cosThetaHXVar.getVal() > 0) {
				// if phi value is smaller than -pi, add 2pi
				if ((phiHXVar.getVal() - 45) < -180) phiTildeHXVar.setVal(phiHXVar.getVal() + 315);
				else phiTildeHXVar.setVal(phiHXVar.getVal() - 45);
			}

			datasetCS.add(RooArgSet(eventWeightCSVar, yVar, ptVar, ptMuPlusVar, etaMuPlusVar, ptMuMinusVar, etaMuMinusVar, cosThetaLabVar, phiLabVar, cosThetaCSVar, phiCSVar, phiTildeCSVar, cosThetaHXVar, phiHXVar, phiTildeHXVar), weightCS);
			datasetHX.add(RooArgSet(eventWeightHXVar, yVar, ptVar, ptMuPlusVar, etaMuPlusVar, ptMuMinusVar, etaMuMinusVar, cosThetaLabVar, phiLabVar, cosThetaCSVar, phiCSVar, phiTildeCSVar, cosThetaHXVar, phiHXVar, phiTildeHXVar), weightHX);
		}
	}

	const char* outputFileName = Form("Y1SGenNoFilterMCDataset_Lambda_Theta%.2f_Phi%.2f_ThetaPhi%.2f.root", lambdaTheta, lambdaPhi, lambdaThetaPhi);

	TFile file(outputFileName, "RECREATE");

	datasetCS.Write();
	datasetHX.Write();

	file.Close();
}

// check the dataset distributions

void draw2DHist(const char* refFrameName = "CS" , Double_t lambdaTheta = 0, Double_t lambdaPhi = 0, Double_t lambdaThetaPhi = 0){

	const char* FileName = Form("Y1SGenNoFilterMCDataset_Lambda_Theta%.2f_Phi%.2f_ThetaPhi%.2f.root", lambdaTheta, lambdaPhi, lambdaThetaPhi);

	TFile* file = openFile(FileName);

	RooDataSet* allDataset = (RooDataSet*)file->Get(Form("MCdataset%s", refFrameName));

	// import the dataset to a workspace
	RooWorkspace wspace("workspace");
	wspace.import(*allDataset);

	RooRealVar cosTheta = *wspace.var(CosThetaVarName(refFrameName));
	RooRealVar phi = *wspace.var(PhiVarName(refFrameName));
	RooRealVar phiTilde = *wspace.var(PhiTildeVarName(refFrameName));

	RooPlot* xframe = phiTilde.frame();

	allDataset->plotOn(xframe);

	TCanvas* c = new TCanvas("c", "Canvas", 800, 600);
    xframe->Draw();
}
