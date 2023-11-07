#include "TEfficiency.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"

// #include "../AnalysisParameters.h"

#include "../Tools/Style/FitDistributions.h"
#include "../Tools/Style/Legends.h"

#include "../ReferenceFrameTransform/Transformations.h"

void DrawAcceptanceMap(TEfficiency* accMap, Int_t ptMin, Int_t ptMax) {
	TCanvas* canvas = new TCanvas(accMap->GetName(), "", 700, 600);
	accMap->Draw("COLZ");

	CMS_lumi(canvas, Form("#varUpsilon(%dS) Pythia 8 MC, no gen filter", gUpsilonState));

	TLatex legend;
	legend.SetTextAlign(22);
	legend.SetTextSize(0.05);
	legend.DrawLatexNDC(.48, .88, Form("|y^{#mu#mu}| < 2.4, %d < p_{T}^{#mu#mu} < %d GeV", ptMin, ptMax));
	legend.DrawLatexNDC(.48, .8, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, p_{T}^{#mu} > 3.5 GeV", gUpsilonState));

	gPad->Update();

	accMap->GetPaintedHistogram()->GetYaxis()->SetRangeUser(-190, 300);
	accMap->GetPaintedHistogram()->GetZaxis()->SetRangeUser(0, 1);

	canvas->SaveAs(Form("AcceptanceMaps/%dS/%s.png", gUpsilonState, accMap->GetName()), "RECREATE");
}

// (cos theta, phi) acceptance maps based on Y events generated without any decay kinematic cut
// MC files available here: /eos/cms/store/group/phys_heavyions/dileptons/MC2015/pp502TeV/TTrees/

void acceptanceMap_noGenFilter(Int_t ptMin = 0, Int_t ptMax = 30) {
	// Read GenOnly Nofilter file
	const char* filename = Form("../Files/OniaTree_Y%dS_GENONLY_NoFilter.root", gUpsilonState);
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

	// can we afford (cos theta, phi, pT) 3D maps for acceptance?
	TEfficiency* accMatrixCS = new TEfficiency("AccMatrixCS", ";cos #theta_{CS}; #varphi_{CS} (#circ);p_{T} (GeV/c);acceptance", NCosThetaBins, gCosThetaMin, gCosThetaMax, NPhiBins, gPhiMin, gPhiMax, gPtBinning[NPtBins] / 2, gPtBinning[0], gPtBinning[NPtBins]); // pT bins 2 GeV wide

	TEfficiency* accMatrixHX = new TEfficiency("AccMatrixHX", ";cos #theta_{HX}; #varphi_{HX} (#circ);p_{T} (GeV/c);acceptance", NCosThetaBins, gCosThetaMin, gCosThetaMax, NPhiBins, gPhiMin, gPhiMax, gPtBinning[NPtBins] / 2, gPtBinning[0], gPtBinning[NPtBins]);

	// (cos theta, phi) 2D distribution maps for Lab, CS and HX frames

	TEfficiency* hGranularLab = new TEfficiency(Form("GranularLab_pt%dto%d", ptMin, ptMax), ";cos #theta_{Lab}; #varphi_{Lab} (#circ);acceptance", NCosThetaBins, gCosThetaMin, gCosThetaMax, NPhiBins, gPhiMin, gPhiMax);
	TEfficiency* hGranularCS = new TEfficiency(Form("GranularCS_pt%dto%d", ptMin, ptMax), ";cos #theta_{CS}; #varphi_{CS} (#circ);acceptance", NCosThetaBins, gCosThetaMin, gCosThetaMax, NPhiBins, gPhiMin, gPhiMax);
	TEfficiency* hGranularHX = new TEfficiency(Form("GranularHX_pt%dto%d", ptMin, ptMax), ";cos #theta_{HX}; #varphi_{HX} (#circ);acceptance", NCosThetaBins, gCosThetaMin, gCosThetaMax, NPhiBins, gPhiMin, gPhiMax);

	// actual analysis binning (defined in AnalysisParameters.h)
	TEfficiency* hAnalysisLab = new TEfficiency(Form("AnalysisLab_pt%dto%d", ptMin, ptMax), ";cos #theta_{Lab}; #varphi_{Lab} (#circ);acceptance", NCosThetaBinsLab, CosThetaBinningLab, NPhiBinsLab, PhiBinningLab);
	TEfficiency* hAnalysisCS = new TEfficiency(Form("AnalysisCS_pt%dto%d", ptMin, ptMax), ";cos #theta_{CS}; #varphi_{CS} (#circ);acceptance", NCosThetaBinsCS, CosThetaBinningCS, NPhiBinsCS, PhiBinningCS);
	TEfficiency* hAnalysisHX = new TEfficiency(Form("AnalysisHX_pt%dto%d", ptMin, ptMax), ";cos #theta_{HX}; #varphi_{HX} (#circ);acceptance", NCosThetaBinsHX, CosThetaBinningHX, NPhiBinsHX, PhiBinningHX);

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

			if (fabs(gen_QQ_LV->Rapidity()) < gRapidityMin || fabs(gen_QQ_LV->Rapidity()) > gRapidityMax) continue; // upsilon within acceptance

			// single-muon acceptance cuts
			gen_mumi_LV = (TLorentzVector*)Gen_QQ_mumi_4mom->At(iGen);
			gen_mupl_LV = (TLorentzVector*)Gen_QQ_mupl_4mom->At(iGen);

			withinAcceptance = (fabs(gen_mupl_LV->Eta()) < 2.4) && (gen_mupl_LV->Pt() > 3.5) && (fabs(gen_mumi_LV->Eta()) < 2.4) && (gen_mumi_LV->Pt() > 3.5);

			// Reference frame transformations
			TVector3 muPlus_CS = MuPlusVector_CollinsSoper(*gen_QQ_LV, *gen_mupl_LV);

			accMatrixCS->Fill(withinAcceptance, muPlus_CS.CosTheta(), muPlus_CS.Phi() * 180 / TMath::Pi(), gen_QQ_LV->Pt());

			TVector3 muPlus_HX = MuPlusVector_Helicity(*gen_QQ_LV, *gen_mupl_LV);

			accMatrixHX->Fill(withinAcceptance, muPlus_HX.CosTheta(), muPlus_HX.Phi() * 180 / TMath::Pi(), gen_QQ_LV->Pt());

			if (gen_QQ_LV->Pt() < ptMin || gen_QQ_LV->Pt() > ptMax) continue; // pt bin of interest

			hGranularCS->Fill(withinAcceptance, muPlus_CS.CosTheta(), muPlus_CS.Phi() * 180 / TMath::Pi());
			hAnalysisCS->Fill(withinAcceptance, muPlus_CS.CosTheta(), muPlus_CS.Phi() * 180 / TMath::Pi());

			hGranularHX->Fill(withinAcceptance, muPlus_HX.CosTheta(), muPlus_HX.Phi() * 180 / TMath::Pi());
			hAnalysisHX->Fill(withinAcceptance, muPlus_HX.CosTheta(), muPlus_HX.Phi() * 180 / TMath::Pi());

			hGranularLab->Fill(withinAcceptance, gen_mupl_LV->CosTheta(), gen_mupl_LV->Phi() * 180 / TMath::Pi());
			hAnalysisLab->Fill(withinAcceptance, gen_mupl_LV->CosTheta(), gen_mupl_LV->Phi() * 180 / TMath::Pi());
		}
	}

	// Set the plot styles
	gStyle->SetPadLeftMargin(.15);
	//gStyle->SetTitleYOffset(.9);
	gStyle->SetPadRightMargin(0.18);
	gStyle->SetPalette(kCividis);
	gStyle->SetNumberContours(256);

	// Draw and save the acceptance map for Lab frame
	DrawAcceptanceMap(hGranularLab, ptMin, ptMax);
	DrawAcceptanceMap(hAnalysisLab, ptMin, ptMax);

	// Draw and save the acceptance map for CS frame
	DrawAcceptanceMap(hGranularCS, ptMin, ptMax);
	DrawAcceptanceMap(hAnalysisCS, ptMin, ptMax);

	// Draw and save the acceptance map for HX frame
	DrawAcceptanceMap(hGranularHX, ptMin, ptMax);
	DrawAcceptanceMap(hAnalysisHX, ptMin, ptMax);

	/// save the results in a file for later usage (in particular for the polarization extraction)
	const char* outputFileName = Form("AcceptanceMaps/%dS/AcceptanceResults.root", gUpsilonState);
	TFile outputFile(outputFileName, "RECREATE");

	accMatrixCS->Write();
	accMatrixHX->Write();

	hGranularLab->Write();
	hAnalysisLab->Write();
	hGranularCS->Write();
	hAnalysisCS->Write();
	hGranularHX->Write();
	hAnalysisHX->Write();

	outputFile.Close();

	cout << endl
	     << "Acceptance maps saved in " << outputFileName << endl;
}
