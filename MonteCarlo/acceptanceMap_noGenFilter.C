#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/Parameters/PhaseSpace.h"

#include "AccEffHelpers.h"

//#include "../Tools/Style/FitDistributions.h"
#include "../Tools/Style/Legends.h"

#include "../ReferenceFrameTransform/Transformations.h"

void DrawAcceptanceMap(TEfficiency* accMap, Int_t ptMin, Int_t ptMax, Int_t iState = 1, Double_t lambdaTheta = 0, Double_t lambdaPhi = 0, Double_t lambdaThetaPhi = 0) {
	TCanvas* canvas = new TCanvas(accMap->GetName(), "", 700, 600);
	accMap->Draw("COLZ");

	CMS_lumi(canvas, Form("Unpolarized #varUpsilon(%dS) Pythia 8 MC", iState));

	TLatex legend;
	legend.SetTextAlign(22);
	legend.SetTextSize(0.05);
	legend.DrawLatexNDC(.48, .88, Form("%s < 2.4, %s", gDimuonRapidityVarTitle, DimuonPtRangeText(ptMin, ptMax)));
	legend.DrawLatexNDC(.48, .8, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, %s", iState, gMuonPtCutText));
	legend.DrawLatexNDC(.48, .72, Form("#lambda_{#theta} = %.2f, #lambda_{#varphi} = %.2f, #lambda_{#theta#varphi} = %.2f", lambdaTheta, lambdaPhi, lambdaThetaPhi));

	gPad->Update();

	accMap->GetPaintedHistogram()->GetXaxis()->CenterTitle();
	accMap->GetPaintedHistogram()->GetYaxis()->CenterTitle();

	accMap->GetPaintedHistogram()->GetYaxis()->SetRangeUser(-190, 300);
	accMap->GetPaintedHistogram()->GetZaxis()->SetRangeUser(0, 1);

	gSystem->mkdir(Form("AcceptanceMaps/%dS", iState), kTRUE);
	canvas->SaveAs(Form("AcceptanceMaps/%dS/%s_LambdaTheta%.2fPhi%.2fThetaPhi%.2f_2018Acc.png", iState, accMap->GetName(), lambdaTheta, lambdaPhi, lambdaThetaPhi), "RECREATE");
}

void DrawAcceptance1DHist(TEfficiency* accHist, Int_t ptMin, Int_t ptMax, Int_t iState = 1, Double_t lambdaTheta = 0, Double_t lambdaPhi = 0, Double_t lambdaThetaPhi = 0) {
	TCanvas* canvas = new TCanvas(accHist->GetName(), "", 600, 600);
	canvas->SetRightMargin(0.05);

	// empty frame for the axes
	TH1D* frameHist = new TH1D("frameHist", "", NCosThetaBinsHX, CosThetaBinningHX);

	frameHist->Draw();

	accHist->SetLineWidth(3);
	accHist->Draw("PL E0 SAME");

	CMS_lumi(canvas, Form("Unpolarized #varUpsilon(%dS) Pythia 8 MC", iState));

	TLatex legend;
	legend.SetTextAlign(22);
	legend.SetTextSize(0.05);
	legend.DrawLatexNDC(.55, .88, Form("%s < 2.4, %s", gDimuonRapidityVarTitle, DimuonPtRangeText(ptMin, ptMax)));
	legend.DrawLatexNDC(.55, .8, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, %s", iState, gMuonPtCutText));
	legend.DrawLatexNDC(.55, .72, Form("#lambda_{#theta} = %.2f, #lambda_{#varphi} = %.2f, #lambda_{#theta#varphi} = %.2f", lambdaTheta, lambdaPhi, lambdaThetaPhi));

	if (strstr(accHist->GetName(), "CS"))
		frameHist->SetXTitle(CosThetaVarTitle("CS"));
	else
		frameHist->SetXTitle(CosThetaVarTitle("HX"));

	frameHist->SetYTitle(TEfficiencyAccMainTitle(iState));

	frameHist->GetXaxis()->CenterTitle();
	frameHist->GetYaxis()->CenterTitle();

	frameHist->GetXaxis()->SetRangeUser(-1, 1);
	frameHist->GetYaxis()->SetRangeUser(0, 1);

	frameHist->GetXaxis()->SetNdivisions(510, kTRUE);

	gSystem->mkdir(Form("AcceptanceMaps/%dS", iState), kTRUE);
	canvas->SaveAs(Form("AcceptanceMaps/%dS/%s_LambdaTheta%.2fPhi%.2fThetaPhi%.2f_2018Acc.png", iState, accHist->GetName(), lambdaTheta, lambdaPhi, lambdaThetaPhi), "RECREATE");
}

const char* Acceptance2DAxisTitle(const char* refFrameName = "CS") {
	return Form(";%s;%s;acceptance", CosThetaVarTitle(refFrameName), PhiAxisTitle(refFrameName));
}

// (cos theta, phi) acceptance maps based on Y events generated without any decay kinematic cut
// MC files available here: /eos/cms/store/group/phys_heavyions/dileptons/MC2015/pp502TeV/TTrees/ (This file was deleted:/)

void acceptanceMap_noGenFilter(Int_t ptMin = 0, Int_t ptMax = 30, Int_t iState = 1, Double_t lambdaTheta = 0, Double_t lambdaPhi = 0, Double_t lambdaThetaPhi = 0) {
	// // Read GenOnly Nofilter file
	// const char* filename = Form("../Files/OniaTree_Y%dS_GENONLY_NoFilter.root", iState, lambdaTheta, lambdaPhi, lambdaThetaPhi);

	// Read GenOnly Nofilter file with polarization weights
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
	TEfficiency* accMatrixCS = TEfficiency3D(Form("%s_LambdaTheta%.2fPhi%.2fThetaPhi%.2f", NominalTEfficiency3DName("CS"), lambdaTheta, lambdaPhi, lambdaThetaPhi), "CS", iState);

	TEfficiency* accMatrixHX = TEfficiency3D(Form("%s_LambdaTheta%.2fPhi%.2fThetaPhi%.2f", NominalTEfficiency3DName("HX"), lambdaTheta, lambdaPhi, lambdaThetaPhi), "HX", iState);

	// (cos theta, phi) 2D distribution maps for Lab, CS and HX frames

	TEfficiency* hGranularLab = CosThetaPhiAcceptance2D(ptMin, ptMax, "Lab");
	TEfficiency* hGranularCS = CosThetaPhiAcceptance2D(ptMin, ptMax, "CS");
	TEfficiency* hGranularHX = CosThetaPhiAcceptance2D(ptMin, ptMax, "HX");

	// actual analysis binning (defined in AnalysisParameters.h)
	TEfficiency* hAnalysisLab = new TEfficiency(Form("AnalysisLab_pt%dto%d_LambdaTheta%.2fPhi%.2fThetaPhi%.2f", ptMin, ptMax, lambdaTheta, lambdaPhi, lambdaThetaPhi), Acceptance2DAxisTitle("Lab"), NCosThetaBinsLab, CosThetaBinningLab, NPhiBinsLab, PhiBinningLab);
	TEfficiency* hAnalysisCS = new TEfficiency(Form("AnalysisCS_pt%dto%d_LambdaTheta%.2fPhi%.2fThetaPhi%.2f", ptMin, ptMax, lambdaTheta, lambdaPhi, lambdaThetaPhi), Acceptance2DAxisTitle("CS"), NCosThetaBinsCS, CosThetaBinningCS, NPhiBinsCS, PhiBinningCS);
	TEfficiency* hAnalysisHX = new TEfficiency(Form("AnalysisHX_pt%dto%d_LambdaTheta%.2fPhi%.2fThetaPhi%.2f", ptMin, ptMax, lambdaTheta, lambdaPhi, lambdaThetaPhi), Acceptance2DAxisTitle("HX"), NCosThetaBinsHX, CosThetaBinningHX, NPhiBinsHX, PhiBinningHX);

	// vs cos theta, for investigation

	TEfficiency* hAccCS1D = CosThetaTEfficiency1D(ptMin, ptMax, "CS", iState, kTRUE);

	TEfficiency* hAccHX1D = CosThetaTEfficiency1D(ptMin, ptMax, "HX", iState, kTRUE);

	TLorentzVector* gen_QQ_LV = new TLorentzVector();
	TLorentzVector* gen_mumi_LV = new TLorentzVector();
	TLorentzVector* gen_mupl_LV = new TLorentzVector();

	Bool_t withinAcceptance;

	Double_t cosThetaCS, phiCS, cosThetaHX, phiHX;

	Float_t weightCS = 0, weightHX = 0;

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

			if (fabs(gen_QQ_LV->Rapidity()) < gRapidityMin || fabs(gen_QQ_LV->Rapidity()) > gRapidityMax) continue; // upsilon within fiducial region

			// single-muon acceptance cuts
			gen_mumi_LV = (TLorentzVector*)Gen_QQ_mumi_4mom->At(iGen);
			gen_mupl_LV = (TLorentzVector*)Gen_QQ_mupl_4mom->At(iGen);

			// withinAcceptance = MuonSimpleAcc(*gen_mupl_LV) && MuonSimpleAcc(*gen_mumi_LV);
			withinAcceptance = MuonWithin2018PbPbAcc(*gen_mupl_LV) && MuonWithin2018PbPbAcc(*gen_mumi_LV);
			
			// Reference frame transformations
			TVector3 muPlus_CS = MuPlusVector_CollinsSoper(*gen_QQ_LV, *gen_mupl_LV);

			cosThetaCS = muPlus_CS.CosTheta();
			phiCS = muPlus_CS.Phi() * 180 / TMath::Pi();

			weightCS = 1 + lambdaTheta * TMath::Power(muPlus_CS.CosTheta(), 2) + lambdaPhi * TMath::Power(std::sin(muPlus_CS.Theta()), 2) * std::cos(2 * muPlus_CS.Phi()) + lambdaThetaPhi * std::sin(2 * muPlus_CS.Theta()) * std::cos(muPlus_CS.Phi());

			accMatrixCS->FillWeighted(withinAcceptance, weightCS, cosThetaCS, phiCS, gen_QQ_LV->Pt());

			TVector3 muPlus_HX = MuPlusVector_Helicity(*gen_QQ_LV, *gen_mupl_LV);

			cosThetaHX = muPlus_HX.CosTheta();
			phiHX = muPlus_HX.Phi() * 180 / TMath::Pi();

			weightHX = 1 + lambdaTheta * TMath::Power(muPlus_HX.CosTheta(), 2) + lambdaPhi * TMath::Power(std::sin(muPlus_HX.Theta()), 2) * std::cos(2 * muPlus_HX.Phi()) + lambdaThetaPhi * std::sin(2 * muPlus_HX.Theta()) * std::cos(muPlus_HX.Phi());

			accMatrixHX->FillWeighted(withinAcceptance, weightHX, cosThetaHX, phiHX, gen_QQ_LV->Pt());

			if (gen_QQ_LV->Pt() > ptMin && gen_QQ_LV->Pt() < ptMax) { // pt bin of interest for the other distributions

				hGranularCS->FillWeighted(withinAcceptance, weightCS, cosThetaCS, phiCS);
				hAnalysisCS->FillWeighted(withinAcceptance, weightHX, cosThetaCS, phiCS);

				hGranularHX->FillWeighted(withinAcceptance, weightCS, cosThetaHX, phiHX);
				hAnalysisHX->FillWeighted(withinAcceptance, weightHX, cosThetaHX, phiHX);

				hGranularLab->Fill(withinAcceptance, gen_mupl_LV->CosTheta(), gen_mupl_LV->Phi() * 180 / TMath::Pi());
				hAnalysisLab->Fill(withinAcceptance, gen_mupl_LV->CosTheta(), gen_mupl_LV->Phi() * 180 / TMath::Pi());

				hAccCS1D->FillWeighted(withinAcceptance, weightCS, cosThetaCS);
				hAccHX1D->FillWeighted(withinAcceptance, weightHX, cosThetaHX);

				// cout << "withinAcceptance: " << withinAcceptance << endl;
				// cout << "weightHX: " << weightHX << endl;
				// cout << "cosThetaHX: " << cosThetaHX << endl;

				// TH1* passedHist = (TH1*)hAccHX1D->GetPassedHistogram();
				// cout << passedHist->FindBin(cosThetaHX) << endl;
				// cout << passedHist->GetBinContent(passedHist->FindBin(cosThetaHX)) << endl;
			}
		}
	}

	// Set the plot styles
	gStyle->SetPadLeftMargin(.15);
	//gStyle->SetTitleYOffset(.9);
	gStyle->SetPadRightMargin(0.18);
	SetColorPalette(gAcceptanceColorPaletteName);

	// Draw and save the acceptance map for Lab frame
	DrawAcceptanceMap(hGranularLab, ptMin, ptMax, iState);
	DrawAcceptanceMap(hAnalysisLab, ptMin, ptMax, iState);

	// Draw and save the acceptance map for CS frame
	DrawAcceptanceMap(hGranularCS, ptMin, ptMax, iState);
	DrawAcceptanceMap(hAnalysisCS, ptMin, ptMax, iState);
	DrawAcceptance1DHist(hAccCS1D, ptMin, ptMax, iState);

	// Draw and save the acceptance map for HX frame
	DrawAcceptanceMap(hGranularHX, ptMin, ptMax, iState);
	DrawAcceptanceMap(hAnalysisHX, ptMin, ptMax, iState);
	DrawAcceptance1DHist(hAccHX1D, ptMin, ptMax, iState);

	/// save the results in a file for later usage
	gSystem->mkdir(Form("AcceptanceMaps/%dS", iState), kTRUE);
	const char* outputFileName = Form("AcceptanceMaps/%dS/AcceptanceResults_2018Acc.root", iState);
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
