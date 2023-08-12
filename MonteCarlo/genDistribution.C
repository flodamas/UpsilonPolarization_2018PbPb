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

#include "../AnalysisParameters.h"

#include "../ReferenceFrameTransform/Transformations.h"

void genDistribution(Int_t iState = 1, Int_t ptMin = 0, Int_t ptMax = 30) {
	const char* filename = Form("../Files/OniaTree_Y%dS_GENONLY_NoFilter.root", iState);
	TFile* file = TFile::Open(filename, "READ");
	if (!file) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
		return;
	}

	cout << "File " << filename << " opened" << endl;

	writeExtraText = true; // if extra text
	extraText = "       Internal";

	TTree* OniaTree = (TTree*)file->Get("hionia/myTree");

	/// OniaTree variables, quite old version
	Int_t Gen_QQ_size;

	TClonesArray* Gen_QQ_4mom = nullptr;
	TClonesArray* Gen_QQ_mumi_4mom = nullptr;
	TClonesArray* Gen_QQ_mupl_4mom = nullptr;

	OniaTree->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size);
	OniaTree->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom);

	OniaTree->SetBranchAddress("Gen_QQ_mumi_4mom", &Gen_QQ_mumi_4mom);
	OniaTree->SetBranchAddress("Gen_QQ_mupl_4mom", &Gen_QQ_mupl_4mom);

	// (cos theta, phi) 2D distribution maps for CS and HX frames

	Int_t nCosThetaBins = 20;
	Float_t cosThetaMin = -1, cosThetaMax = 1;

	//	Double_t cosThetaBinning[] = {-1, -0.8, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1};
	//	nCosThetaBins = sizeof(cosThetaBinning) / sizeof(Double_t) - 1;

	Int_t nPhiBins = 21;
	Float_t phiMin = -180, phiMax = 240;

	//Double_t phiBinning[] = {-TMath::Pi(), -2.5, -2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2, 2.5, 180};
	//	nPhiBins = sizeof(cosThetaBinning) / sizeof(Double_t) - 1;

	TH2F* hCS = new TH2F("hCS", ";cos #theta_{CS}; #varphi_{CS} (#circ);number of generated #varUpsilon's", nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);
	TH2F* hHX = new TH2F("hHX", ";cos #theta_{HX}; #varphi_{HX} (#circ);number of generated #varUpsilon's", nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);
	TH2F* hLab = new TH2F("hLab", ";cos #theta_{lab}; #varphi_{lab} (#circ);number of generated #varUpsilon's", nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);

	TLorentzVector* gen_QQ_LV = new TLorentzVector();
	TLorentzVector* gen_mumi_LV = new TLorentzVector();
	TLorentzVector* gen_mupl_LV = new TLorentzVector();

	Long64_t totEntries = OniaTree->GetEntries();

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

			gen_mumi_LV = (TLorentzVector*)Gen_QQ_mumi_4mom->At(iGen);

			gen_mupl_LV = (TLorentzVector*)Gen_QQ_mupl_4mom->At(iGen);

			// get positive muon's coordinates in the studied reference frames
			TVector3 muPlus_CS = MuPlusVector_CollinsSoper(*gen_QQ_LV, *gen_mupl_LV);

			TVector3 muPlus_HX = MuPlusVector_Helicity(*gen_QQ_LV, *gen_mupl_LV);

			hCS->Fill(muPlus_CS.CosTheta(), muPlus_CS.Phi() * 180 / TMath::Pi());
			hHX->Fill(muPlus_HX.CosTheta(), muPlus_HX.Phi() * 180 / TMath::Pi());
			hLab->Fill(gen_mupl_LV->CosTheta(), gen_mupl_LV->Phi() * 180 / TMath::Pi());
		}
	}

	gStyle->SetPadLeftMargin(.15);
	//gStyle->SetTitleYOffset(.9);
	gStyle->SetPadRightMargin(0.18);
	gStyle->SetPalette(kRainBow);
	gStyle->SetNumberContours(256);

	TLatex* legend = new TLatex();
	legend->SetTextAlign(22);
	legend->SetTextSize(0.05);

	auto* canvasCS = new TCanvas("canvasCS", "", 700, 600);
	hCS->Draw("COLZ");

	CMS_lumi(canvasCS, Form("#varUpsilon(%dS) Pythia 8 + EvtGen MC, no filter", iState));

	legend->DrawLatexNDC(.48, .88, Form("|y^{#mu#mu}| < 2.4, %d < p_{T}^{#mu#mu} < %d GeV", ptMin, ptMax));

	gPad->Update();

	hCS->GetYaxis()->SetRangeUser(-190, 240);

	canvasCS->SaveAs(Form("GenMaps/CS_pt%dto%dGeV.png", ptMin, ptMax), "RECREATE");

	auto* canvasHX = new TCanvas("canvasHX", "", 700, 600);
	hHX->Draw("COLZ");

	CMS_lumi(canvasHX, Form("#varUpsilon(%dS) Pythia 8 + EvtGen MC, no filter", iState));

	legend->DrawLatexNDC(.48, .88, Form("|y^{#mu#mu}| < 2.4, %d < p_{T}^{#mu#mu} < %d GeV", ptMin, ptMax));

	gPad->Update();

	hHX->GetYaxis()->SetRangeUser(-190, 240);

	canvasHX->SaveAs(Form("GenMaps/HX_pt%dto%dGeV.png", ptMin, ptMax), "RECREATE");

	auto* canvasLab = new TCanvas("canvasLab", "", 700, 600);
	hLab->Draw("COLZ");

	CMS_lumi(canvasLab, Form("#varUpsilon(%dS) Pythia 8 + EvtGen MC, no filter", iState));

	legend->DrawLatexNDC(.48, .88, Form("|y^{#mu#mu}| < 2.4, %d < p_{T}^{#mu#mu} < %d GeV", ptMin, ptMax));

	gPad->Update();

	hLab->GetYaxis()->SetRangeUser(-190, 240);

	canvasLab->SaveAs(Form("GenMaps/Lab_pt%dto%dGeV.png", ptMin, ptMax), "RECREATE");
}
