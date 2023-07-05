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

#include "../ReferenceFrameTransform/Transformations.h"

// (cos theta, phi) acceptance maps based on Y events generated without any cuts
// MC file available here: /eos/cms/store/group/phys_heavyions/dileptons/MC2015/pp502TeV/TTrees/OniaTree_Ups1SMM_5p02TeV_TuneCUETP8M1_nofilter_pp502Fall15-MCRUN2_71_V1-v1_GENONLY.root

void acceptanceMap_noGenFilter(Int_t iState = 1, Int_t ptMin = 0, Int_t ptMax = 30) {

	// Read GenOnly Nofilter file
	const char* filename = Form("../Files/OniaTree_Y%dS_GENONLY_NoFilter.root", iState);
	TFile* file = TFile::Open(filename, "READ");
	if (!file) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
		return;
	}

	cout << "File " << filename << " opened" << endl;

	// Put the text, "CMS Internal", on the right top of the plot 
	// writeExtraText = true; // if extra text
	// extraText = "       Internal";
	writeExtraText = false;

	TTree* OniaTree = (TTree*)file->Get("hionia/myTree");

	/// OniaTree variables, quite old version (since this genonly file from 2015)
	Int_t Gen_QQ_size;

	TClonesArray* Gen_QQ_4mom = nullptr;
	TClonesArray* Gen_QQ_mumi_4mom = nullptr;
	TClonesArray* Gen_QQ_mupl_4mom = nullptr;

	OniaTree->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size);
	OniaTree->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom);

	OniaTree->SetBranchAddress("Gen_QQ_mumi_4mom", &Gen_QQ_mumi_4mom);
	OniaTree->SetBranchAddress("Gen_QQ_mupl_4mom", &Gen_QQ_mupl_4mom);

	// (cos theta, phi) 2D distribution maps for Lab, CS and HX frames
	Int_t nCosThetaBins = 20;
	Float_t cosThetaMin = -1, cosThetaMax = 1;

	//	Double_t cosThetaBinning[] = {-1, -0.8, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1};
	//	nCosThetaBins = sizeof(cosThetaBinning) / sizeof(Double_t) - 1;

	Int_t nPhiBins = 23;
	Float_t phiMin = -180, phiMax = 280;

	//Double_t phiBinning[] = {-TMath::Pi(), -2.5, -2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2, 2.5, 180};
	//	nPhiBins = sizeof(cosThetaBinning) / sizeof(Double_t) - 1;


	TEfficiency* hLab = new TEfficiency(Form("Lab_pt%dto%dGeV", ptMin, ptMax), ";cos #theta_{Lab}; #varphi_{Lab} (#circ);acceptance", nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);
	TEfficiency* hCS = new TEfficiency(Form("CS_pt%dto%dGeV", ptMin, ptMax), ";cos #theta_{CS}; #varphi_{CS} (#circ);acceptance", nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);
	TEfficiency* hHX = new TEfficiency(Form("HX_pt%dto%dGeV", ptMin, ptMax), ";cos #theta_{HX}; #varphi_{HX} (#circ);acceptance", nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);

	TLorentzVector* gen_QQ_LV = new TLorentzVector();
	TLorentzVector* gen_mumi_LV = new TLorentzVector();
	TLorentzVector* gen_mupl_LV = new TLorentzVector();

	Bool_t withinAcceptance;

	Long64_t totEntries = OniaTree->GetEntries();

	// Loop over the events
	for (Long64_t iEvent = 0; iEvent < (totEntries); iEvent++) {
		if (iEvent % 10000 == 0) {
			cout << Form("\rProcessing event %lld / %lld (%.0f%%)", iEvent, totEntries, 100. * iEvent / totEntries) << flush;
		}

		OniaTree->GetEntry(iEvent);

		// loop over all gen upsilons
		for (int iGen = 0; iGen < Gen_QQ_size; iGen++) {
			gen_QQ_LV = (TLorentzVector*)Gen_QQ_4mom->At(iGen);

			if (gen_QQ_LV->Pt() < ptMin || gen_QQ_LV->Pt() > ptMax) continue; // pt bin of interest

			if (fabs(gen_QQ_LV->Rapidity()) > 2.4) continue; // upsilon within acceptance

			// single-muon acceptance cuts
			gen_mumi_LV = (TLorentzVector*)Gen_QQ_mumi_4mom->At(iGen);
			gen_mupl_LV = (TLorentzVector*)Gen_QQ_mupl_4mom->At(iGen);

			withinAcceptance = (fabs(gen_mupl_LV->Eta()) < 2.4) && (gen_mupl_LV->Pt() > 3.5) && (fabs(gen_mumi_LV->Eta()) < 2.4) && (gen_mumi_LV->Pt() > 3.5);

			// Reference frame transformations
			TVector3 muPlus_CS = MuPlusVector_CollinsSoper(*gen_QQ_LV, *gen_mupl_LV);
			TVector3 muPlus_HX = MuPlusVector_Helicity(*gen_QQ_LV, *gen_mupl_LV);
			
			hLab->Fill(withinAcceptance, gen_mupl_LV->CosTheta(), gen_mupl_LV->Phi()*180/TMath::Pi());	
			hCS->Fill(withinAcceptance, muPlus_CS.CosTheta(), muPlus_CS.Phi() * 180 / TMath::Pi());
			hHX->Fill(withinAcceptance, muPlus_HX.CosTheta(), muPlus_HX.Phi() * 180 / TMath::Pi());
		}
	}

	// Set the plot styles
	gStyle->SetPadLeftMargin(.15);
	//gStyle->SetTitleYOffset(.9);
	gStyle->SetPadRightMargin(0.18);
	gStyle->SetPalette(kCividis);
	gStyle->SetNumberContours(256);

	TLatex* legend = new TLatex();
	legend->SetTextAlign(22);
	legend->SetTextSize(0.05);


	// Draw and save the acceptance map for Lab frame
	auto* canvasLab = new TCanvas("canvasLab", "", 700, 600);
	hLab->Draw("COLZ");

	CMS_lumi(canvasLab, Form("#varUpsilon(%dS) Pythia 8 + EvtGen MC, no filter", iState));

	legend->DrawLatexNDC(.48, .88, Form("|y^{#mu#mu}| < 2.4, %d < p_{T}^{#mu#mu} < %d GeV", ptMin, ptMax));
	legend->DrawLatexNDC(.48, .8, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, p_{T}^{#mu} > 3.5 GeV", iState));

	gPad->Update();

	hLab->GetPaintedHistogram()->GetYaxis()->SetRangeUser(-190, 300);
	hLab->GetPaintedHistogram()->GetZaxis()->SetRangeUser(0, 1);

	canvasLab->SaveAs(Form("AcceptanceMaps/%dS/Lab_pt%dto%dGeV.png", iState, ptMin, ptMax), "RECREATE");


	// Draw and save the acceptance map for CS frame
	auto* canvasCS = new TCanvas("canvasCS", "", 700, 600);
	hCS->Draw("COLZ");

	CMS_lumi(canvasCS, Form("#varUpsilon(%dS) Pythia 8 + EvtGen MC, no filter", iState));

	legend->DrawLatexNDC(.48, .88, Form("|y^{#mu#mu}| < 2.4, %d < p_{T}^{#mu#mu} < %d GeV", ptMin, ptMax));
	legend->DrawLatexNDC(.48, .8, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, p_{T}^{#mu} > 3.5 GeV", iState));

	gPad->Update();

	hCS->GetPaintedHistogram()->GetYaxis()->SetRangeUser(-190, 300);
	hCS->GetPaintedHistogram()->GetZaxis()->SetRangeUser(0, 1);

	canvasCS->SaveAs(Form("AcceptanceMaps/%dS/CS_pt%dto%dGeV.png", iState, ptMin, ptMax), "RECREATE");


	// Draw and save the acceptance map for HX frame
	auto* canvasHX = new TCanvas("canvasHX", "", 700, 600);
	hHX->Draw("COLZ");

	CMS_lumi(canvasHX, Form("#varUpsilon(%dS) Pythia 8 + EvtGen MC, no filter", iState));

	legend->DrawLatexNDC(.48, .88, Form("|y^{#mu#mu}| < 2.4, %d < p_{T}^{#mu#mu} < %d GeV", ptMin, ptMax));
	legend->DrawLatexNDC(.48, .8, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, p_{T}^{#mu} > 3.5 GeV", iState));

	gPad->Update();

	hHX->GetPaintedHistogram()->GetYaxis()->SetRangeUser(-190, 300);
	hHX->GetPaintedHistogram()->GetZaxis()->SetRangeUser(0, 1);

	canvasHX->SaveAs(Form("AcceptanceMaps/%dS/HX_pt%dto%dGeV.png", iState, ptMin, ptMax), "RECREATE");


	/// save the results in a file for later usage (in particular for the polarization extraction)
	const char* outputFileName = Form("AcceptanceMaps/%dS/AcceptanceResults.root", iState);
	TFile outputFile(outputFileName, "RECREATE");

	hLab->Write();
	hCS->Write();
	hHX->Write();

	outputFile.Close();

	cout << endl
	     << "Acceptance maps saved in " << outputFileName << endl;
}