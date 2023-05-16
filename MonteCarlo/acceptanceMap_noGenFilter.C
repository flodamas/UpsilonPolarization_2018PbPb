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

// (cos theta, phi) acceptance maps based on Y events generated without any cuts
// MC file available here: /eos/cms/store/group/phys_heavyions/dileptons/MC2015/pp502TeV/TTrees/OniaTree_Ups1SMM_5p02TeV_TuneCUETP8M1_nofilter_pp502Fall15-MCRUN2_71_V1-v1_GENONLY.root

void acceptanceMap_noGenFilter(Int_t iState = 1, Int_t ptMin = 0, Int_t ptMax = 30) {
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

	Float_t Gen_QQ_phi, Gen_QQ_costheta, Gen_QQ_pt, Gen_QQ_y, Gen_QQ_px, Gen_QQ_py, Gen_QQ_pz, Gen_QQ_E, Gen_QQ_m;
	Float_t Gen_mupl_phi, Gen_mupl_costheta, Gen_mupl_pt, Gen_mupl_y, Gen_mupl_eta, Gen_mupl_px, Gen_mupl_py, Gen_mupl_pz, Gen_mupl_E;
	Float_t Gen_mumi_phi, Gen_mumi_costheta, Gen_mumi_pt, Gen_mumi_y, Gen_mumi_eta, Gen_mumi_px, Gen_mumi_py, Gen_mumi_pz, Gen_mumi_E;

	// ******** Set beam energy for the Collins-Soper reference frame ******** //
	double sqrt_S_NN = 5.02;                 //(Center of mass Energy per nucleon pair in TeV)
	double beam1_p = sqrt_S_NN * 1000. / 2.; //(in GeV) (Note. sqrt_S_NN = sqrt(2*E1*E2+2*p1*p2) = 2E1 when two beams have the same E)
	double beam1_E = beam1_p;
	double beam2_p = -beam1_p;
	double beam2_E = beam1_E;
	//double delta = 0; //(Angle between ZHX(Z-axis in the Helicity frame) and ZCS(Z-axis in the Collins-Soper frame))

	// (cos theta, phi) 2D distribution maps for CS and HX frames

	Int_t nCosThetaBins = 20;
	Float_t cosThetaMin = -1, cosThetaMax = 1;

	//	Double_t cosThetaBinning[] = {-1, -0.8, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1};
	//	nCosThetaBins = sizeof(cosThetaBinning) / sizeof(Double_t) - 1;

	Int_t nPhiBins = 23;
	Float_t phiMin = -180, phiMax = 280;

	//Double_t phiBinning[] = {-TMath::Pi(), -2.5, -2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2, 2.5, 180};
	//	nPhiBins = sizeof(cosThetaBinning) / sizeof(Double_t) - 1;

	TEfficiency* hCS = new TEfficiency(Form("CS_pt%dto%dGeV", ptMin, ptMax), ";cos #theta_{CS}; #varphi_{CS} (#circ);acceptance", nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);
	TEfficiency* hHX = new TEfficiency(Form("HX_pt%dto%dGeV", ptMin, ptMax), ";cos #theta_{HX}; #varphi_{HX} (#circ);acceptance", nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);

	TLorentzVector* gen_QQ_LV = new TLorentzVector();
	TLorentzVector* gen_mumi_LV = new TLorentzVector();
	TLorentzVector* gen_mupl_LV = new TLorentzVector();

	Bool_t withinAcceptance;

	Long64_t totEntries = OniaTree->GetEntries();

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

			// Transformations and rotations

			Gen_QQ_phi = gen_QQ_LV->Phi();
			Gen_QQ_costheta = gen_QQ_LV->CosTheta();
			Gen_QQ_pt = gen_QQ_LV->Pt();
			Gen_QQ_y = gen_QQ_LV->Rapidity();
			Gen_QQ_px = gen_QQ_LV->Px();
			Gen_QQ_py = gen_QQ_LV->Py();
			Gen_QQ_pz = gen_QQ_LV->Pz();
			Gen_QQ_E = gen_QQ_LV->Energy();
			Gen_QQ_m = gen_QQ_LV->M();

			Gen_mupl_phi = gen_mupl_LV->Phi();
			Gen_mupl_costheta = gen_mupl_LV->CosTheta();
			Gen_mupl_pt = gen_mupl_LV->Pt();
			Gen_mupl_y = gen_mupl_LV->Rapidity();
			Gen_mupl_eta = gen_mupl_LV->Eta();
			Gen_mupl_px = gen_mupl_LV->Px();
			Gen_mupl_py = gen_mupl_LV->Py();
			Gen_mupl_pz = gen_mupl_LV->Pz();
			Gen_mupl_E = gen_mupl_LV->Energy();

			Gen_mumi_phi = gen_mumi_LV->Phi();
			Gen_mumi_costheta = gen_mumi_LV->CosTheta();
			Gen_mumi_pt = gen_mumi_LV->Pt();
			Gen_mumi_y = gen_mumi_LV->Rapidity();
			Gen_mumi_eta = gen_mumi_LV->Eta();
			Gen_mumi_px = gen_mumi_LV->Px();
			Gen_mumi_py = gen_mumi_LV->Py();
			Gen_mumi_pz = gen_mumi_LV->Pz();
			Gen_mumi_E = gen_mumi_LV->Energy();

			// ******** Construct 4-momentum vector of upsilon and muons (Lab Frame) ******** //
			// (documetation of TVector3 and TLorentzVector: https://root.cern.ch/root/html534/guides/users-guide/PhysicsVectors.html#lorentz-boost)
			TVector3 upsPvecLab(Gen_QQ_px, Gen_QQ_py, Gen_QQ_pz);
			TLorentzVector ups4MomLab(upsPvecLab, Gen_QQ_E);

			TVector3 muplPvecLab(Gen_mupl_px, Gen_mupl_py, Gen_mupl_pz);
			TLorentzVector mupl4MomLab(muplPvecLab, Gen_mupl_E);

			TVector3 mumiPvecLab(Gen_mumi_px, Gen_mumi_py, Gen_mumi_pz);
			TLorentzVector mumi4MomLab(mumiPvecLab, Gen_mumi_E);

			TVector3 beam1PvecLab(0, 0, beam1_p);
			TLorentzVector beam14MomLab(beam1PvecLab, beam1_E);

			TVector3 beam2PvecLab(0, 0, beam2_p);
			TLorentzVector beam24MomLab(beam2PvecLab, beam2_E);

			// ******** Transform variables of muons from the lab frame to the upsilon's rest frame ******** //
			TLorentzVector ups4MomBoosted(upsPvecLab, Gen_QQ_E);
			TLorentzVector mupl4MomBoosted(muplPvecLab, Gen_mupl_E);
			TLorentzVector mumi4MomBoosted(mumiPvecLab, Gen_mumi_E);

			//(Note. TLorentzVector.BoostVector() gives beta(=px/E,py/E,pz/E) of the parents)
			//(TLorentzVector.Boost() boosts from the rod frame to the lab frame, so plug in -beta to get lab to rod)
			ups4MomBoosted.Boost(-ups4MomLab.BoostVector());
			mupl4MomBoosted.Boost(-ups4MomLab.BoostVector());
			mumi4MomBoosted.Boost(-ups4MomLab.BoostVector());

			// ******** Rotate the coordinate ******** //
			TVector3 muplPvecBoosted(mupl4MomBoosted.Px(), mupl4MomBoosted.Py(), mupl4MomBoosted.Pz());
			TVector3 mumiPvecBoosted(mumi4MomBoosted.Px(), mumi4MomBoosted.Py(), mumi4MomBoosted.Pz());

			//(Note. TVector3.Rotate() rotates the vectors, not the coordinates, so should rotate -phi and -theta)
			muplPvecBoosted.RotateZ(-upsPvecLab.Phi());
			muplPvecBoosted.RotateY(-upsPvecLab.Theta());
			mumiPvecBoosted.RotateZ(-upsPvecLab.Phi());
			mumiPvecBoosted.RotateY(-upsPvecLab.Theta());

			TLorentzVector mupl4MomBoostedRot(muplPvecBoosted, mupl4MomBoosted.E());

			// ******** HX to CS (rotation from HX frame to CS frame) ******** //
			// (1. Boost two beams to upsilon's rest frame)
			// (2. Rotate the coordinate)
			// (3. Get angle between two beams(b1 and -b2), and between b1 and ZHX in the upsilon's rest frame)
			// (4. Calculate delta (angle btw ZHX and ZCS))

			// ******** Transform variables of beams from the lab frame to the upsilon's rest frame ******** //
			TLorentzVector beam14MomBoosted(beam1PvecLab, beam1_E);
			TLorentzVector beam24MomBoosted(beam2PvecLab, beam2_E);

			beam14MomBoosted.Boost(-ups4MomLab.BoostVector());
			beam24MomBoosted.Boost(-ups4MomLab.BoostVector());

			// ******** Rotate the coordinate ******** //
			TVector3 beam1PvecBoosted(beam14MomBoosted.Px(), beam14MomBoosted.Py(), beam14MomBoosted.Pz());
			TVector3 beam2PvecBoosted(beam24MomBoosted.Px(), beam24MomBoosted.Py(), beam24MomBoosted.Pz());

			beam1PvecBoosted.RotateZ(-upsPvecLab.Phi());
			beam1PvecBoosted.RotateY(-upsPvecLab.Theta());
			beam2PvecBoosted.RotateZ(-upsPvecLab.Phi());
			beam2PvecBoosted.RotateY(-upsPvecLab.Theta());

			// ******** Calculate the angle between z_HX and z_CS ******** //
			TVector3 ZHXunitVec(0, 0, 1);                                    //(define z_HX unit vector)
			double Angle_B1ZHX = beam1PvecBoosted.Angle(ZHXunitVec);         //(angle between beam1 and z_HX)
			double Angle_B2ZHX = beam2PvecBoosted.Angle(-ZHXunitVec);        //(angle between beam2 and -z_HX =(-beam2 and z_HX) )
			double Angle_B1miB2 = beam1PvecBoosted.Angle(-beam2PvecBoosted); //(angle between beam1 and -beam2)

			double delta = 0; //(define and initialize the angle between z_HX and z_CS)

			// Maths for calculating the angle between z_HX and z_CS is different depending on the sign of the beam1's z-coordinate)
			if (Angle_B1ZHX > Angle_B2ZHX)
				delta = Angle_B2ZHX + Angle_B1miB2 / 2.;
			else if (Angle_B1ZHX < Angle_B2ZHX)
				delta = Angle_B1ZHX + Angle_B1miB2 / 2.;
			else
				cout << "beam1PvecBoosted.Pz() = 0?" << endl;

			// ******** Rotate the coordinate along the y-axis by the angle between z_HX and z_CS ******** //
			TVector3 muplPvecBoostedCS(muplPvecBoosted.Px(), muplPvecBoosted.Py(), muplPvecBoosted.Pz());
			TVector3 mumiPvecBoostedCS(mumiPvecBoosted.Px(), mumiPvecBoosted.Py(), mumiPvecBoosted.Pz());

			ZHXunitVec.RotateY(delta);
			muplPvecBoostedCS.RotateY(delta);
			mumiPvecBoostedCS.RotateY(delta);

			hCS->Fill(withinAcceptance, muplPvecBoostedCS.CosTheta(), muplPvecBoostedCS.Phi() * 180 / TMath::Pi());
			hHX->Fill(withinAcceptance, muplPvecBoosted.CosTheta(), muplPvecBoosted.Phi() * 180 / TMath::Pi());
		}
	}

	gStyle->SetPadLeftMargin(.15);
	//gStyle->SetTitleYOffset(.9);
	gStyle->SetPadRightMargin(0.18);
	gStyle->SetPalette(kCividis);
	gStyle->SetNumberContours(256);

	TLatex* legend = new TLatex();
	legend->SetTextAlign(22);
	legend->SetTextSize(0.05);

	auto* canvasCS = new TCanvas("canvasCS", "", 700, 600);
	hCS->Draw("COLZ");

	CMS_lumi(canvasCS, Form("#varUpsilon(%dS) Pythia 8 + EvtGen MC, no filter", iState));

	legend->DrawLatexNDC(.48, .88, Form("|y^{#mu#mu}| < 2.4, %d < p_{T}^{#mu#mu} < %d GeV", ptMin, ptMax));
	legend->DrawLatexNDC(.48, .8, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, p_{T}^{#mu} > 3.5 GeV", iState));

	gPad->Update();

	hCS->GetPaintedHistogram()->GetYaxis()->SetRangeUser(-190, 300);
	hCS->GetPaintedHistogram()->GetZaxis()->SetRangeUser(0, 1);

	canvasCS->SaveAs(Form("AcceptanceMaps/%dS/CS_pt%dto%dGeV.png", iState, ptMin, ptMax), "RECREATE");

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

	hCS->Write();
	hHX->Write();

	outputFile.Close();

	cout << endl
	     << "Acceptance maps saved in " << outputFileName << endl;
}
