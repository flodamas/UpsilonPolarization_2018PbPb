#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/Parameters/PhaseSpace.h"

#include "AccEffHelpers.h"

//#include "../Tools/Style/FitDistributions.h"
#include "../Tools/Style/Legends.h"
#include "../Tools/Parameters/EfficiencyWeights.h"

#include "../ReferenceFrameTransform/Transformations.h"

void DrawAcceptanceMap(TEfficiency* accMap, Int_t ptMin, Int_t ptMax, TString muonAccName, Double_t lambdaTheta = 0, Double_t lambdaPhi = 0, Double_t lambdaThetaPhi = 0, Int_t iState = 1) {
	TCanvas* canvas = new TCanvas(accMap->GetName(), "", 750, 600);
	accMap->Draw("COLZ");

	CMS_lumi(canvas, Form("#varUpsilon(%dS) Pythia 8 (5.02 TeV)", iState));

	TLatex legend;
	legend.SetTextAlign(22);
	legend.SetTextSize(0.05);
	legend.DrawLatexNDC(.5, .88, Form("%s, %s", DimuonRapidityRangeText(), DimuonPtRangeText(ptMin, ptMax)));
	legend.DrawLatexNDC(.5, .81, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, %s", iState, gMuonPtCutText));
	legend.DrawLatexNDC(.5, .74, Form("inputs: #lambda_{#theta} = %.2f, #lambda_{#varphi} = %.2f, #lambda_{#theta#varphi} = %.2f", lambdaTheta, lambdaPhi, lambdaThetaPhi));

	gPad->Update();

	accMap->GetPaintedHistogram()->GetXaxis()->CenterTitle();
	accMap->GetPaintedHistogram()->GetYaxis()->CenterTitle();

	//accMap->GetPaintedHistogram()->GetYaxis()->SetRangeUser(0, 220);
	accMap->GetPaintedHistogram()->GetZaxis()->SetRangeUser(0, 1);

	gPad->RedrawAxis();

	const char* path = AcceptanceResultsPath(muonAccName.Data());

	gSystem->mkdir(path, kTRUE);
	canvas->SaveAs(Form("%s/%s.png", path, accMap->GetName()), "RECREATE");
	// canvas->SaveAs(Form("AcceptanceMaps/%dS/%s_fullPhi.png", iState, accMap->GetName()), "RECREATE");
}

void DrawAcceptance1DHist(TEfficiency* accHist, Int_t ptMin, Int_t ptMax, TString muonAccName, Double_t lambdaTheta = 0, Double_t lambdaPhi = 0, Double_t lambdaThetaPhi = 0, Int_t iState = 1) {
	TCanvas* canvas = new TCanvas(accHist->GetName(), "", 600, 600);
	canvas->SetRightMargin(0.05);

	// empty frame for the axes
	TH1D* frameHist = new TH1D("frameHist", "", NCosThetaBinsHX, CosThetaBinningHX);

	frameHist->Draw();

	accHist->SetLineWidth(3);
	accHist->Draw("PL E0 SAME");

	CMS_lumi(canvas, Form("#varUpsilon(%dS) Pythia 8 (5.02 TeV)", iState));

	TLatex legend;
	legend.SetTextAlign(22);
	legend.SetTextSize(0.05);
	legend.DrawLatexNDC(.55, .9, Form("%s < 2.4, %s", gDimuonRapidityVarTitle, DimuonPtRangeText(ptMin, ptMax)));
	legend.DrawLatexNDC(.55, .83, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, %s", iState, gMuonPtCutText));
	legend.DrawLatexNDC(.55, .75, Form("#lambda_{#theta} = %.2f, #lambda_{#varphi} = %.2f, #lambda_{#theta#varphi} = %.2f", lambdaTheta, lambdaPhi, lambdaThetaPhi));

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

	const char* path = AcceptanceResultsPath(muonAccName.Data());

	gSystem->mkdir(path, kTRUE);
	canvas->SaveAs(Form("%s/%s.png", path, accHist->GetName()), "RECREATE");
	// canvas->SaveAs(Form("AcceptanceMaps/%dS/%s_fullPhi.png", iState, accHist->GetName()), "RECREATE");
}

const char* Acceptance2DAxisTitle(const char* refFrameName = "CS") {
	return Form(";%s;%s;acceptance", CosThetaVarTitle(refFrameName), PhiAxisTitle(refFrameName));
}

// (cos theta, phi) acceptance maps based on Y events generated without any decay kinematic cut
// MC files available here: /eos/cms/store/group/phys_heavyions/dileptons/MC2015/pp502TeV/TTrees/ (This file was deleted:/)

void acceptanceMap_noGenFilter(Int_t ptMin = 0, Int_t ptMax = 30, Bool_t isPhiFolded = kFALSE, TString muonAccName = gMuonAccName, Float_t lambdaTheta = 0, Float_t lambdaPhi = 0, Float_t lambdaThetaPhi = 0, Int_t iState = 1, Long64_t startEvent = 0, Long64_t endEvent = 10e6, int fileIndex = 0) {
	// Read GenOnly Nofilter file with polarization weights
	// const char* filename = Form("../Files/OniaTree_Y%dS_GENONLY_NoFilter.root", iState);
	// const char* filename = Form("../Files/Oniatree_Upsilon%dS_5p02TeV_TuneCP5_141X_GenOnly_combined.root", iState);
	const char* filename = Form("../Files/Oniatree_Upsilon%dS_5p02TeV_TuneCP5_141X_GenOnly_50M.root", iState);

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
	// Int_t Gen_QQ_size;
	Short_t Gen_QQ_size;

	TClonesArray* Gen_QQ_4mom = nullptr;
    
    TClonesArray* Gen_mu_4mom = nullptr;

    TClonesArray* CloneArr_QQ = nullptr;
    TClonesArray* CloneArr_mu = nullptr;

    Short_t Gen_QQ_mupl_idx[1000];
    Short_t Gen_QQ_mumi_idx[1000];

	// TClonesArray* Gen_QQ_4mom = nullptr;
	// TClonesArray* Gen_QQ_mumi_4mom = nullptr;
	// TClonesArray* Gen_QQ_mupl_4mom = nullptr;

	OniaTree->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size);
	OniaTree->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom);

	OniaTree->SetBranchAddress("Gen_mu_4mom", &Gen_mu_4mom);

    OniaTree->SetBranchAddress("Gen_QQ_mumi_idx", &Gen_QQ_mumi_idx);
    OniaTree->SetBranchAddress("Gen_QQ_mupl_idx", &Gen_QQ_mupl_idx);

	// OniaTree->SetBranchAddress("Gen_QQ_mumi_4mom", &Gen_QQ_mumi_4mom);
	// OniaTree->SetBranchAddress("Gen_QQ_mupl_4mom", &Gen_QQ_mupl_4mom);

	// TLorentzVector* gen_QQ_LV = new TLorentzVector();
	// TLorentzVector* gen_mumi_LV = new TLorentzVector();
	// TLorentzVector* gen_mupl_LV = new TLorentzVector();

	TLorentzVector* gen_QQ_LV = nullptr;
	TLorentzVector* gen_mumi_LV = nullptr;
	TLorentzVector* gen_mupl_LV = nullptr;

    // Long64_t totEntries = OniaTree->GetEntries(); 

	// (cos theta, phi, pT) 3D maps for final acceptance correction, variable size binning for the stats
	TEfficiency* accMatrixLab = TEfficiency3D(NominalTEfficiency3DName("Lab", lambdaTheta, lambdaPhi, lambdaThetaPhi), "Lab", iState, isPhiFolded);
	TEfficiency* accMatrixCS = TEfficiency3D(NominalTEfficiency3DName("CS", lambdaTheta, lambdaPhi, lambdaThetaPhi), "CS", iState, isPhiFolded);
	TEfficiency* accMatrixHX = TEfficiency3D(NominalTEfficiency3DName("HX", lambdaTheta, lambdaPhi, lambdaThetaPhi), "HX", iState, isPhiFolded);

	// (cos theta, phi) 2D distribution maps for Lab, CS and HX frames

	TEfficiency* hGranularLab = CosThetaPhiAcceptance2D("Lab", ptMin, ptMax, lambdaTheta, lambdaPhi, lambdaThetaPhi);
	TEfficiency* hGranularCS = CosThetaPhiAcceptance2D("CS", ptMin, ptMax, lambdaTheta, lambdaPhi, lambdaThetaPhi);
	TEfficiency* hGranularHX = CosThetaPhiAcceptance2D("HX", ptMin, ptMax, lambdaTheta, lambdaPhi, lambdaThetaPhi);

	// actual analysis binning (defined in AnalysisParameters.h)
	TEfficiency* hAnalysisLab = new TEfficiency(Form("Analysis%s", TEfficiencyEndName("Lab", ptMin, ptMax, lambdaTheta, lambdaPhi, lambdaThetaPhi)), Acceptance2DAxisTitle("Lab"), NCosThetaBinsLab, CosThetaBinningLab, NPhiBinsLab, PhiBinningLab);
	TEfficiency* hAnalysisCS = new TEfficiency(Form("Analysis%s", TEfficiencyEndName("CS", ptMin, ptMax, lambdaTheta, lambdaPhi, lambdaThetaPhi)), Acceptance2DAxisTitle("CS"), NCosThetaBinsCS, CosThetaBinningCS, NPhiBinsCS, PhiBinningCS);
	TEfficiency* hAnalysisHX = new TEfficiency(Form("Analysis%s", TEfficiencyEndName("HX", ptMin, ptMax, lambdaTheta, lambdaPhi, lambdaThetaPhi)), Acceptance2DAxisTitle("HX"), NCosThetaBinsHX, CosThetaBinningHX, NPhiBinsHX, PhiBinningHX);

	// vs cos theta, for investigation

	TEfficiency* hAccCS1D = CosThetaTEfficiency1D("CS", ptMin, ptMax, iState, kTRUE, lambdaTheta, lambdaPhi, lambdaThetaPhi);

	TEfficiency* hAccHX1D = CosThetaTEfficiency1D("HX", ptMin, ptMax, iState, kTRUE, lambdaTheta, lambdaPhi, lambdaThetaPhi);

	// TLorentzVector* gen_QQ_LV = new TLorentzVector();
	// TLorentzVector* gen_mumi_LV = new TLorentzVector();
	// TLorentzVector* gen_mupl_LV = new TLorentzVector();

	Bool_t withinAcceptance;

	Double_t cosThetaCS, phiCS, cosThetaHX, phiHX, cosThetaLab, phiLab;

	Double_t dimuonPtWeight = 1.0;

	Float_t weightCS = 0, weightHX = 0;

	Long64_t totEntries = OniaTree->GetEntries();
	double counter = 0;
	// Loop over the events
	// for (Long64_t iEvent = 0; iEvent < (totEntries); iEvent++) {
	for (Long64_t iEvent = startEvent; iEvent < (totEntries); iEvent++) {
		if (iEvent % 10000 == 0) {
			cout << Form("\rProcessing event %lld / %lld (%.0f%%)", iEvent, totEntries, 100. * iEvent / totEntries) << flush;
		}

		OniaTree->GetEntry(iEvent);
		// if (iEvent < 50) continue;
		// if (iEvent == 10000000) break;
		// if (iEvent == 20000000) break;
		// if (iEvent == 30000000) break;

		if (endEvent > 0 && iEvent >= endEvent) break; // limit the number of events to process

		// loop over all gen upsilons
		for (int iGen = 0; iGen < Gen_QQ_size; iGen++) {
			gen_QQ_LV = (TLorentzVector*)Gen_QQ_4mom->At(iGen);

			if (fabs(gen_QQ_LV->Rapidity()) < gRapidityMin || fabs(gen_QQ_LV->Rapidity()) > gRapidityMax) continue; // upsilon within fiducial region

			// single-muon acceptance cuts
			// gen_mumi_LV = (TLorentzVector*)Gen_QQ_mumi_4mom->At(iGen);
			// gen_mupl_LV = (TLorentzVector*)Gen_QQ_mupl_4mom->At(iGen);

			gen_mupl_LV = (TLorentzVector*)Gen_mu_4mom->At(Gen_QQ_mupl_idx[iGen]);
			gen_mumi_LV = (TLorentzVector*)Gen_mu_4mom->At(Gen_QQ_mumi_idx[iGen]);

			withinAcceptance = MuonKinematicsWithinLimits(*gen_mupl_LV, muonAccName) && MuonKinematicsWithinLimits(*gen_mumi_LV, muonAccName);

			// pt Weight at gen level
			double gen_QQ_pt = gen_QQ_LV->Pt();
			dimuonPtWeight = Get_GenPtWeight(gen_QQ_LV->Rapidity(), gen_QQ_pt);

			// cout << "iEvent: " << iEvent << endl;
			// cout << "iGen: " << iGen << endl;
			// cout << "withinAcceptance: " << withinAcceptance << endl;
			// cout << "gen_QQ_LV->Pt(): " << gen_QQ_LV->Pt() << endl;
			// cout << "accName: " << accName << endl;

			TVector3 muPlus_Lab = gen_mupl_LV->Vect();

			cosThetaLab = muPlus_Lab.CosTheta();
			phiLab = muPlus_Lab.Phi() * 180 / TMath::Pi();

			if (isPhiFolded == kTRUE)
				phiLab = fabs(phiLab);

			accMatrixLab->FillWeighted(withinAcceptance, dimuonPtWeight, cosThetaLab, phiLab, gen_QQ_LV->Pt());

			// Reference frame transformations
			TVector3 muPlus_CS = MuPlusVector_CollinsSoper(*gen_QQ_LV, *gen_mupl_LV);

			cosThetaCS = muPlus_CS.CosTheta();

			phiCS = muPlus_CS.Phi() * 180 / TMath::Pi();
			if (isPhiFolded == kTRUE) phiCS = fabs(phiCS);
			// cout << "cosThetaCS: " << cosThetaCS << endl;
			// cout << "phiCS: " << phiCS << endl;

			if (isPhiFolded == kTRUE)
				weightCS = 1 + lambdaTheta * TMath::Power(muPlus_CS.CosTheta(), 2) + lambdaPhi * TMath::Power(std::sin(muPlus_CS.Theta()), 2) * std::cos(2 * fabs(muPlus_CS.Phi())) + lambdaThetaPhi * std::sin(2 * muPlus_CS.Theta()) * std::cos(fabs(muPlus_CS.Phi()));
			else
				weightCS = 1 + lambdaTheta * TMath::Power(muPlus_CS.CosTheta(), 2) + lambdaPhi * TMath::Power(std::sin(muPlus_CS.Theta()), 2) * std::cos(2 * muPlus_CS.Phi()) + lambdaThetaPhi * std::sin(2 * muPlus_CS.Theta()) * std::cos(muPlus_CS.Phi());

			// cout << "weightCS: " << weightCS << endl;
			// cout << "cosThetaCS: " << cosThetaCS << endl;
			// cout << "phiCS: " << phiCS << endl;

			accMatrixCS->FillWeighted(withinAcceptance, weightCS * dimuonPtWeight, cosThetaCS, phiCS, gen_QQ_LV->Pt());

			TVector3 muPlus_HX = MuPlusVector_Helicity(*gen_QQ_LV, *gen_mupl_LV);

			cosThetaHX = muPlus_HX.CosTheta();
			phiHX = muPlus_HX.Phi() * 180 / TMath::Pi();
			if (isPhiFolded == kTRUE)
				phiHX = fabs(phiHX);

			if (isPhiFolded == kTRUE)
				weightHX = 1 + lambdaTheta * TMath::Power(muPlus_HX.CosTheta(), 2) + lambdaPhi * TMath::Power(std::sin(muPlus_HX.Theta()), 2) * std::cos(2 * fabs(muPlus_HX.Phi())) + lambdaThetaPhi * std::sin(2 * muPlus_HX.Theta()) * std::cos(fabs(muPlus_HX.Phi()));
			else
				weightHX = 1 + lambdaTheta * TMath::Power(muPlus_HX.CosTheta(), 2) + lambdaPhi * TMath::Power(std::sin(muPlus_HX.Theta()), 2) * std::cos(2 * muPlus_HX.Phi()) + lambdaThetaPhi * std::sin(2 * muPlus_HX.Theta()) * std::cos(muPlus_HX.Phi());

				// cout << "weightHX: " << weightHX << endl;
			// cout << "cosThetaHX: " << cosThetaHX << endl;
			// cout << "phiHX: " << phiHX << endl;

			accMatrixHX->FillWeighted(withinAcceptance, weightHX * dimuonPtWeight, cosThetaHX, phiHX, gen_QQ_LV->Pt());

			if (gen_QQ_LV->Pt() > ptMin && gen_QQ_LV->Pt() < ptMax) { // pt bin of interest for the other distributions

				hGranularCS->FillWeighted(withinAcceptance, weightCS, cosThetaCS, phiCS);
				hAnalysisCS->FillWeighted(withinAcceptance, weightHX, cosThetaCS, phiCS);

				hGranularHX->FillWeighted(withinAcceptance, weightCS, cosThetaHX, phiHX);
				hAnalysisHX->FillWeighted(withinAcceptance, weightHX, cosThetaHX, phiHX);

				//	hGranularLab->Fill(withinAcceptance, gen_mupl_LV->CosTheta(), gen_mupl_LV->Phi() * 180 / TMath::Pi());
				//	hAnalysisLab->Fill(withinAcceptance, gen_mupl_LV->CosTheta(), gen_mupl_LV->Phi() * 180 / TMath::Pi());

				hAccCS1D->FillWeighted(withinAcceptance, weightCS, cosThetaCS);
				hAccHX1D->FillWeighted(withinAcceptance, weightHX, cosThetaHX);

				// cout << "withinAcceptance: " << withinAcceptance << endl;

				// cout << "cosThetaHX: " << cosThetaHX << endl;

				// TH1* passedHist = (TH1*)hAccHX1D->GetPassedHistogram();
				// cout << passedHist->FindBin(cosThetaHX) << endl;
				// cout << passedHist->GetBinContent(passedHist->FindBin(cosThetaHX)) << endl;
			}
			// cout << "weightHX: " << weightHX << endl;
			if (withinAcceptance && (gen_QQ_LV->Pt() > 2) && (gen_QQ_LV->Pt() < 6) && (cosThetaHX > -0.7) && (cosThetaHX < -0.42) && (phiHX > 60) && (phiHX < 120)) counter++;
			// cout << "" << endl;
		}
	}
	cout << "\ncounter: " << counter << endl;
	// Set the plot styles
	gStyle->SetPadLeftMargin(.15);
	//gStyle->SetTitleYOffset(.9);
	gStyle->SetPadRightMargin(0.18);
	SetColorPalette(gAcceptanceColorPaletteName);

	// Draw and save the acceptance map for Lab frame

	//DrawAcceptanceMap(hGranularLab, ptMin, ptMax, path, lambdaTheta, lambdaPhi, lambdaThetaPhi, iState);
	//DrawAcceptanceMap(hAnalysisLab, ptMin, ptMax, path, lambdaTheta, lambdaPhi, lambdaThetaPhi, iState);

	// Draw and save the acceptance map for CS frame
	DrawAcceptanceMap(hGranularCS, ptMin, ptMax, muonAccName, lambdaTheta, lambdaPhi, lambdaThetaPhi, iState);

	DrawAcceptanceMap(hAnalysisCS, ptMin, ptMax, muonAccName, lambdaTheta, lambdaPhi, lambdaThetaPhi, iState);
	DrawAcceptance1DHist(hAccCS1D, ptMin, ptMax, muonAccName, lambdaTheta, lambdaPhi, lambdaThetaPhi, iState);

	// Draw and save the acceptance map for HX frame
	DrawAcceptanceMap(hGranularHX, ptMin, ptMax, muonAccName, lambdaTheta, lambdaPhi, lambdaThetaPhi, iState);
	DrawAcceptanceMap(hAnalysisHX, ptMin, ptMax, muonAccName, lambdaTheta, lambdaPhi, lambdaThetaPhi, iState);
	DrawAcceptance1DHist(hAccHX1D, ptMin, ptMax, muonAccName, lambdaTheta, lambdaPhi, lambdaThetaPhi, iState);

	// cout << "Phi axis range: " << accMatrixCS->GetTotalHistogram()->GetYaxis()->GetXmin()
	//  << " to " << accMatrixCS->GetTotalHistogram()->GetYaxis()->GetXmax() << endl;

	/// save the results in a file for later usage
	const char* path = AcceptanceResultsPath(muonAccName.Data());

	gSystem->mkdir(path, kTRUE);
	const char* outputFileName = Form("%s/AcceptanceResults_dimuonPtWeight%s_file%d.root", path, isPhiFolded ? "" : "_fullPhi", fileIndex);
	// const char* outputFileName = Form("%s/AcceptanceResults_dimuonPtWeight%s_50M.root", path, isPhiFolded ? "" : "_fullPhi");

	TFile outputFile(outputFileName, "UPDATE");

	//accMatrixLab->Write();
	accMatrixCS->Write();
	accMatrixHX->Write();

	//hGranularLab->Write();
	//hAnalysisLab->Write();
	hGranularCS->Write();
	hAnalysisCS->Write();
	hGranularHX->Write();
	hAnalysisHX->Write();

	hAccCS1D->Write();
	hAccHX1D->Write();

	outputFile.Write();
	outputFile.Close();

	if (BeVerbose) cout << "\nAcceptance maps saved in " << outputFileName << endl;
}

void mergeAcceptanceOutputs(TString outputName = "AcceptanceResults_dimuonPtWeight_fullPhi_50M.root",
    						TString inputPattern = "AcceptanceResults_dimuonPtWeight_fullPhi_file",  // prefix pattern
    						TString muonAccName = "UpsilonTriggerThresholds") {

    const char* path = AcceptanceResultsPath(muonAccName.Data());

	TFileMerger merger(kFALSE); // 'kFALSE' = don't compact memory
    merger.OutputFile(Form("%s/%s", path, outputName.Data()), "RECREATE");

	void* dirp = gSystem->OpenDirectory(path);	

	const char* file;
		
    std::cout << "Merging files in directory: " << path << std::endl;

    while ((file = gSystem->GetDirEntry(dirp))) {
        TString fname = file;
        if (!fname.EndsWith(".root")) continue;

        if (fname.BeginsWith(inputPattern)) {
            TString fullPath = TString::Format("%s/%s", path, fname.Data());
            std::cout << "Adding file: " << fullPath << std::endl;
            merger.AddFile(fullPath);
        }
    }

    gSystem->FreeDirectory(dirp);

    Bool_t success = merger.Merge();
    if (success)
        std::cout << "Merged files saved to " << outputName << std::endl;
    else
        std::cerr << "Merging failed." << std::endl;
}

void run_acceptanceMap_noGenFilter(Long64_t totNEvent = 50000000, Long64_t step = 5000000) {
	
	for (Long64_t startEvent = 0; startEvent < totNEvent; startEvent += step) {
		Long64_t endEvent = startEvent + step;
		if (endEvent > totNEvent) endEvent = totNEvent;

		cout << "Processing events from " << startEvent << " to " << endEvent << endl;

		gSystem->Exec(Form("root -l -b -q 'acceptanceMap_noGenFilter.C+(0,30,kFALSE,\"%s\",0,0,0,1,%lld,%lld,%lld)'", gMuonAccName.Data(), startEvent, endEvent, startEvent / step));
	}

	mergeAcceptanceOutputs("AcceptanceResults_dimuonPtWeight_fullPhi_50M.root", "AcceptanceResults_dimuonPtWeight_fullPhi_file",  gMuonAccName);

	return;
}