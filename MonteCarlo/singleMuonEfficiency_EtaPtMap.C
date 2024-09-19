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

#include "../Tools/Style/Legends.h"

#include "../Tools/Style/FitDistributions.h"
#include "../Tools/Style/Legends.h"

#include "../AnalysisParameters.h"

#include "../Tools/Parameters/CentralityValues.h"

#include "../Tools/Parameters/MuonScaleFactors.h"

void setGuillaumePalette() {
	const Int_t NCont = 999;
	const Int_t NRGBs = 7;
	Double_t stops[NRGBs] = {0.00, 0.0999, 0.1, 0.34, 0.61, 0.84, 1.00};
	Double_t red[NRGBs] = {0.99, 0.99, 0.0, 0.00, 0.87, 1.00, 0.51};
	Double_t green[NRGBs] = {0.86, 0.9, 0.0, 0.81, 1.00, 0.20, 0.00};
	Double_t blue[NRGBs] = {0.99, 0.99, 0.0, 1.00, 0.12, 0.00, 0.00};
	// const Int_t NRGBs = 7;
	// Double_t stops[NRGBs] = { 0.00, 0.025, 0.1, 0.34, 0.61, 0.84, 1.00 };
	// Double_t red[NRGBs]   = { 0.96, 0.99, 0.0, 0.00, 0.87, 1.00, 0.51 };
	// Double_t green[NRGBs] = { 0.85, 0.0, 0.0, 0.81, 1.00, 0.20, 0.00 };
	// Double_t blue[NRGBs]  = { 0.96, 0.99, 0.0, 1.00, 0.12, 0.00, 0.00 };
	TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
	gStyle->SetNumberContours(NCont);
}

TLine* MyLine(float x1, float y1, float x2, float y2, Color_t lineColor = kBlack, int lineStyle = 1, int lineWidth = 2) {
	auto* line = new TLine(x1, y1, x2, y2);
	line->SetLineColor(lineColor);
	line->SetLineStyle(lineStyle);
	line->SetLineWidth(lineWidth);
	return line;
}

TLine* drawTrigger2018acc(float etaMin = 1., float ptMax = 3., Color_t lineColor = kMagenta, int lineStyle = 1, int lineWidth = 3) {
	float eta1 = 1.2, eta2 = 2.1, eta3 = 2.4;

	float pt1 = 3.5, pt1prime = 5.47 - 1.89 * eta1, pt2 = 1.5;

	auto* line1 = MyLine(etaMin, pt1, eta1, pt1, lineColor, lineStyle, lineWidth);
	if (pt1 < ptMax && etaMin < eta1) line1->Draw("l");

	auto* line2prime = MyLine(eta1, pt1, eta1, pt1prime, lineColor, lineStyle, lineWidth);
	if (pt1prime < ptMax) line2prime->Draw("l");	

	auto* line2 = MyLine(eta1, pt1prime, eta2, pt2, lineColor, lineStyle, lineWidth);
	if (pt2 < ptMax) line2->Draw("l");

	auto* line3 = MyLine(eta2, pt2, eta3, pt2, lineColor, lineStyle, lineWidth);
	line3->Draw("l");

	auto* lineVert = MyLine(eta3, pt2, eta3, ptMax, lineColor, lineStyle, lineWidth);
	lineVert->Draw("l");

	return lineVert;
}

void singleMuonEfficiency_EtaPtMap() {
	const char* filename = "../Files/OniaTree_Y1S_pThat2_HydjetDrumMB_miniAOD.root";
	TFile* file = TFile::Open(filename, "READ");
	if (!file) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
		return;
	}

	cout << "File " << filename << " opened" << endl;

	TTree* OniaTree = (TTree*)file->Get("hionia/myTree");

	writeExtraText = true;
	extraText = "       Internal";

	/// OniaTree variables

	Float_t Gen_weight;

	TClonesArray* Gen_mu_4mom = nullptr;

	ULong64_t HLTriggers;
	Int_t Centrality;
	TClonesArray* CloneArr_mu = nullptr;
	Short_t Reco_mu_size;
	Short_t Gen_mu_whichRec[1000];
	ULong64_t Reco_mu_trig[1000];
	Int_t Reco_mu_SelectionType[1000];
	//(parameters for quality cuts)
	Int_t Reco_mu_nPixWMea[1000];
	Int_t Reco_mu_nTrkWMea[1000];
	Float_t Reco_mu_dxy[1000];
	Float_t Reco_mu_dz[1000];

	Short_t Gen_mu_size;

	// event variables
	OniaTree->SetBranchAddress("Gen_weight", &Gen_weight);
	OniaTree->SetBranchAddress("Centrality", &Centrality);
	OniaTree->SetBranchAddress("HLTriggers", &HLTriggers);

	// gen-level variables
	OniaTree->SetBranchAddress("Gen_mu_size", &Gen_mu_size);
	OniaTree->SetBranchAddress("Gen_mu_4mom", &Gen_mu_4mom);

	OniaTree->SetBranchAddress("Gen_mu_whichRec", Gen_mu_whichRec);

	OniaTree->SetBranchAddress("Reco_mu_4mom", &CloneArr_mu);

	OniaTree->SetBranchAddress("Reco_mu_trig", Reco_mu_trig);
	OniaTree->SetBranchAddress("Reco_mu_SelectionType", Reco_mu_SelectionType);

	OniaTree->SetBranchAddress("Reco_mu_nPixWMea", &Reco_mu_nPixWMea);
	OniaTree->SetBranchAddress("Reco_mu_nTrkWMea", &Reco_mu_nTrkWMea);
	OniaTree->SetBranchAddress("Reco_mu_dxy", &Reco_mu_dxy);
	OniaTree->SetBranchAddress("Reco_mu_dz", &Reco_mu_dz);

	// granular binning for acceptance studies

	TEfficiency* hEffMap = new TEfficiency("hEffMap", ";muon |#eta|;muon p_{T} (GeV);(reco + ID + trigger) / generated", 26, 0, 2.6, 60, 0, 6);

	// loop variables
	TLorentzVector* genLorentzVector = new TLorentzVector();

	double eventWeight, genEta, genPt;

	Bool_t allGood, isRecoMatched, passHLTFilter, isTrackerAndGlobal, isHybridSoft;

	Long64_t totEntries = OniaTree->GetEntries();

	for (Long64_t iEvent = 0; iEvent < totEntries; iEvent++) {
		if (iEvent % 10000 == 0) {
			cout << Form("\rProcessing event %lld / %lld (%.0f%%)", iEvent, totEntries, 100. * iEvent / totEntries) << flush;
		}

		OniaTree->GetEntry(iEvent);
		//genLorentzVector->Clear();

		// event selection

		if (Centrality >= 2 * 90) continue; // discard events with centrality >= 90% in 2018 data

		if (!((HLTriggers & (ULong64_t)(1 << (gUpsilonHLTBit - 1))) == (ULong64_t)(1 << (gUpsilonHLTBit - 1)))) continue;

		eventWeight = Gen_weight * FindNcoll(Centrality);

		// loop over all gen muons
		for (int iGen = 0; iGen < Gen_mu_size; iGen++) {
			genLorentzVector = (TLorentzVector*)Gen_mu_4mom->At(iGen);

			genEta = fabs(genLorentzVector->Eta());

			genPt = genLorentzVector->Pt();

			// go to reco level
			Int_t iReco = Gen_mu_whichRec[iGen];

			/// all the reconstructed upsilons must pass the conditions below!

			isRecoMatched = iReco > -1;

			// pass L2 or L3 HLT filters

			passHLTFilter = (Reco_mu_trig[iReco] & ((ULong64_t)pow(2, gL2FilterBit))) == ((ULong64_t)pow(2, gL2FilterBit)) || (Reco_mu_trig[iReco] & ((ULong64_t)pow(2, gL3FilterBit))) == ((ULong64_t)pow(2, gL3FilterBit));

			// global AND tracker muons
			isTrackerAndGlobal = (Reco_mu_SelectionType[iReco] & 2) && (Reco_mu_SelectionType[iReco] & 8);

			// passing hybrid-soft Id
			isHybridSoft = (Reco_mu_nTrkWMea[iReco] > 5) && (Reco_mu_nPixWMea[iReco] > 0) && (fabs(Reco_mu_dxy[iReco]) < 0.3) && (fabs(Reco_mu_dz[iReco]) < 20.);

			allGood = isRecoMatched && passHLTFilter && isTrackerAndGlobal && isHybridSoft;

			eventWeight *= tnp_weight_trk_pbpb(genEta, 0);
			eventWeight *= tnp_weight_muid_pbpb(genPt, genEta, 0);

			hEffMap->FillWeighted(allGood, eventWeight, genEta, genPt);

		} // end of gen upsilon loop
	}

	cout << endl;

	//gStyle->SetPadTopMargin(.02);
	gStyle->SetTitleYOffset(.9);
	gStyle->SetPadLeftMargin(.13);
	gStyle->SetPadRightMargin(0.18);
	//SetColorPalette(gEfficiencyColorPaletteName);

	setGuillaumePalette();

	TCanvas* canvas = new TCanvas("canvas", "", 700, 600);
	hEffMap->Draw("COLZ");

	auto* triggerAcc = drawTrigger2018acc(0, 6);

	auto* legend = new TLegend(.14, .25, .38, .2);
	legend->AddEntry(triggerAcc, "acceptance for 2018 PbPb data", "l");
	legend->Draw();

	CMS_lumi(canvas, "#varUpsilon(1S) Hydjet-embedded MC");

	gPad->Update();

	hEffMap->GetPaintedHistogram()->GetZaxis()->SetRangeUser(0, 1);

	canvas->SaveAs("EfficiencyMaps/SingleMuonTotalEfficiency.png", "RECREATE");
}
