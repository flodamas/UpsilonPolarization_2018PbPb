#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

//#include "../Tools/Style/FitDistributions.h"
#include "../Tools/Style/Legends.h"

#include "../ReferenceFrameTransform/Transformations.h"

void DrawAcceptanceMap(TEfficiency* accMap, Int_t ptMin, Int_t ptMax, Int_t iState = 1) {
	TCanvas* canvas = new TCanvas(accMap->GetName(), "", 700, 600);
	accMap->Draw("COLZ");

	CMS_lumi(canvas, Form("Unpolarized #varUpsilon(%dS) Pythia 8 MC", iState));

	TLatex legend;
	legend.SetTextAlign(22);
	legend.SetTextSize(0.05);
	legend.DrawLatexNDC(.48, .88, Form("%s < 2.4, %s", gDimuonRapidityVarTitle, DimuonPtRangeText(ptMin, ptMax)));
	legend.DrawLatexNDC(.48, .8, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, %s", iState, gMuonPtCutText));

	gPad->Update();

	accMap->GetPaintedHistogram()->GetXaxis()->CenterTitle();
	accMap->GetPaintedHistogram()->GetYaxis()->CenterTitle();

	accMap->GetPaintedHistogram()->GetYaxis()->SetRangeUser(-190, 300);
	accMap->GetPaintedHistogram()->GetZaxis()->SetRangeUser(0, 1);

	gSystem->mkdir(Form("AcceptanceMaps/%dS", iState), kTRUE);
	canvas->SaveAs(Form("AcceptanceMaps/%dS/%s.png", iState, accMap->GetName()), "RECREATE");
}

// (cos theta, phi) acceptance maps based on Y events generated without any decay kinematic cut
// MC files available here: /eos/cms/store/group/phys_heavyions/dileptons/MC2015/pp502TeV/TTrees/

void acceptanceMap_noGenFilter(Int_t ptMin = 0, Int_t ptMax = 30, Int_t iState = 1) {
	// Read GenOnly Nofilter file
	const char* filename = Form("../Files/OniaTree_Y%dS_GENONLY_NoFilter.root", iState);
	TFile* file = TFile::Open(filename, "READ");
	if (!file) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
		return;
	}

	cout << "File " << filename << " opened" << endl;

	// Put the text, "CMS Internal", on the right top of the plot
	writeExtraText = true; // if extra text
	extraText = "       Internal";

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

	// (cos theta, phi, pT) 3D maps for final acceptance correction, variable size binning for the stats
	TEfficiency* accMatrixCS = new TEfficiency("AccMatrixCS", "acceptance;cos #theta_{CS}; #varphi_{CS} (#circ);p_{T} (GeV/c)", NCosThetaFineBins, gCosThetaFineBinning, NPhiFineBins, gPhiFineBinning, NPtFineBins, gPtFineBinning);

	TEfficiency* accMatrixHX = new TEfficiency("AccMatrixHX", "acceptance;cos #theta_{HX}; #varphi_{HX} (#circ);p_{T} (GeV/c)", NCosThetaFineBins, gCosThetaFineBinning, NPhiFineBins, gPhiFineBinning, NPtFineBins, gPtFineBinning);

	// (cos theta, phi) 2D distribution maps for Lab, CS and HX frames

	TEfficiency* hGranularLab = new TEfficiency(Form("GranularLab_pt%dto%d", ptMin, ptMax), ";cos #theta_{Lab}; #varphi_{Lab} (#circ);acceptance", NCosThetaBins, gCosThetaMin, gCosThetaMax, NPhiBins, gPhiMin, gPhiMax);
	TEfficiency* hGranularCS = new TEfficiency(Form("GranularCS_pt%dto%d", ptMin, ptMax), ";cos #theta_{CS}; #varphi_{CS} (#circ);acceptance", NCosThetaBins, gCosThetaMin, gCosThetaMax, NPhiBins, gPhiMin, gPhiMax);
	TEfficiency* hGranularHX = new TEfficiency(Form("GranularHX_pt%dto%d", ptMin, ptMax), ";cos #theta_{HX}; #varphi_{HX} (#circ);acceptance", NCosThetaBins, gCosThetaMin, gCosThetaMax, NPhiBins, gPhiMin, gPhiMax);

	// actual analysis binning (defined in AnalysisParameters.h)
	TEfficiency* hAnalysisLab = new TEfficiency(Form("AnalysisLab_pt%dto%d", ptMin, ptMax), ";cos #theta_{Lab}; #varphi_{Lab} (#circ);acceptance", NCosThetaBinsLab, CosThetaBinningLab, NPhiBinsLab, PhiBinningLab);
	TEfficiency* hAnalysisCS = new TEfficiency(Form("AnalysisCS_pt%dto%d", ptMin, ptMax), ";cos #theta_{CS}; #varphi_{CS} (#circ);acceptance", NCosThetaBinsCS, CosThetaBinningCS, NPhiBinsCS, PhiBinningCS);
	TEfficiency* hAnalysisHX = new TEfficiency(Form("AnalysisHX_pt%dto%d", ptMin, ptMax), ";cos #theta_{HX}; #varphi_{HX} (#circ);acceptance", NCosThetaBinsHX, CosThetaBinningHX, NPhiBinsHX, PhiBinningHX);

	// vs cos theta, for investigation
	TEfficiency* hAccCS1D = new TEfficiency(Form("AccCosThetaCS_pt%dto%d", ptMin, ptMax), ";cos #theta_{CS}; acceptance", NCosThetaBins, gCosThetaMin, gCosThetaMax);

	TEfficiency* hAccHX1D = new TEfficiency(Form("AccCosThetaHX_pt%dto%d", ptMin, ptMax), ";cos #theta_{HX}; acceptance", NCosThetaBins, gCosThetaMin, gCosThetaMax);

	TLorentzVector* gen_QQ_LV = new TLorentzVector();
	TLorentzVector* gen_mumi_LV = new TLorentzVector();
	TLorentzVector* gen_mupl_LV = new TLorentzVector();

	// set error calculation options (default: kFCP (Frequentist Clopper-Pearson))
	// possible options: kFWilson (Frequentist Wilson interval), kBUniform (Beyesian Uniform prior)...
	// documentation about the options: https://root.cern.ch/doc/master/classTEfficiency.html#a415d9689c9e50f1d6982d043d5fd11ec
	accMatrixCS->SetStatisticOption(TEfficiency::kFCP);
	accMatrixHX->SetStatisticOption(TEfficiency::kFCP);

	hGranularLab->SetStatisticOption(TEfficiency::kFCP);
	hGranularCS->SetStatisticOption(TEfficiency::kFCP);
	hGranularHX->SetStatisticOption(TEfficiency::kFCP);

	hAnalysisLab->SetStatisticOption(TEfficiency::kFCP);
	hAnalysisCS->SetStatisticOption(TEfficiency::kFCP);
	hAnalysisHX->SetStatisticOption(TEfficiency::kFCP);

	hAccCS1D->SetStatisticOption(TEfficiency::kFCP);
	hAccHX1D->SetStatisticOption(TEfficiency::kFCP);

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

			withinAcceptance = (fabs(gen_mupl_LV->Eta()) < 2.4) && (gen_mupl_LV->Pt() > gMuonPtCut) && (fabs(gen_mumi_LV->Eta()) < 2.4) && (gen_mumi_LV->Pt() > gMuonPtCut);

			// Reference frame transformations
			TVector3 muPlus_CS = MuPlusVector_CollinsSoper(*gen_QQ_LV, *gen_mupl_LV);

			accMatrixCS->Fill(withinAcceptance, muPlus_CS.CosTheta(), muPlus_CS.Phi() * 180 / TMath::Pi(), gen_QQ_LV->Pt());

			TVector3 muPlus_HX = MuPlusVector_Helicity(*gen_QQ_LV, *gen_mupl_LV);

			accMatrixHX->Fill(withinAcceptance, muPlus_HX.CosTheta(), muPlus_HX.Phi() * 180 / TMath::Pi(), gen_QQ_LV->Pt());

			if (gen_QQ_LV->Pt() > ptMin && gen_QQ_LV->Pt() < ptMax) { // pt bin of interest for the other distributions

				hGranularCS->Fill(withinAcceptance, muPlus_CS.CosTheta(), muPlus_CS.Phi() * 180 / TMath::Pi());
				hAnalysisCS->Fill(withinAcceptance, muPlus_CS.CosTheta(), muPlus_CS.Phi() * 180 / TMath::Pi());

				hGranularHX->Fill(withinAcceptance, muPlus_HX.CosTheta(), muPlus_HX.Phi() * 180 / TMath::Pi());
				hAnalysisHX->Fill(withinAcceptance, muPlus_HX.CosTheta(), muPlus_HX.Phi() * 180 / TMath::Pi());

				hGranularLab->Fill(withinAcceptance, gen_mupl_LV->CosTheta(), gen_mupl_LV->Phi() * 180 / TMath::Pi());
				hAnalysisLab->Fill(withinAcceptance, gen_mupl_LV->CosTheta(), gen_mupl_LV->Phi() * 180 / TMath::Pi());

				hAccCS1D->Fill(withinAcceptance, muPlus_CS.CosTheta());
				hAccHX1D->Fill(withinAcceptance, muPlus_HX.CosTheta());
			}
		}
	}

	// Set the plot styles
	gStyle->SetPadLeftMargin(.15);
	//gStyle->SetTitleYOffset(.9);
	gStyle->SetPadRightMargin(0.18);
	gStyle->SetPalette(kCividis);
	gStyle->SetNumberContours(256);

	// Draw and save the acceptance map for Lab frame
	DrawAcceptanceMap(hGranularLab, ptMin, ptMax, iState);
	DrawAcceptanceMap(hAnalysisLab, ptMin, ptMax, iState);

	// Draw and save the acceptance map for CS frame
	DrawAcceptanceMap(hGranularCS, ptMin, ptMax, iState);
	DrawAcceptanceMap(hAnalysisCS, ptMin, ptMax, iState);

	// Draw and save the acceptance map for HX frame
	DrawAcceptanceMap(hGranularHX, ptMin, ptMax, iState);
	DrawAcceptanceMap(hAnalysisHX, ptMin, ptMax, iState);

	/// save the results in a file for later usage
	gSystem->mkdir(Form("AcceptanceMaps/%dS", iState), kTRUE);
	const char* outputFileName = Form("AcceptanceMaps/%dS/AcceptanceResults_FCP.root", iState);
	TFile outputFile(outputFileName, "UPDATE");

	accMatrixCS->Write();
	accMatrixHX->Write();

	hGranularLab->Write();
	hAnalysisLab->Write();
	hGranularCS->Write();
	hAnalysisCS->Write();
	hGranularHX->Write();
	hAnalysisHX->Write();

	hAccCS1D->Write();
	hAccHX1D->Write();

	outputFile.Close();

	if (BeVerbose) cout << "\nAcceptance maps saved in " << outputFileName << endl;
}
