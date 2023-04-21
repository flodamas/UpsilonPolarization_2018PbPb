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

void acceptanceMap(Int_t iState = 1, Int_t ptMin = 0, Int_t ptMax = 30) {
	const char* filename = Form("../Files/OniaTree_Y%dS_miniAOD_HydjetEmbeddedMC.root", iState);
	TFile* file = TFile::Open(filename, "READ");
	if (!file) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
		return;
	}

	cout << "File " << filename << " opened" << endl;

	writeExtraText = true; // if extra text
	extraText = "       Internal";

	Int_t centMin = 0, centMax = 90;

	TTree* OniaTree = (TTree*)file->Get("hionia/myTree");

	/// OniaTree variables

	Float_t Gen_weight;
	Float_t SumET_HF;

	TClonesArray* Gen_QQ_4mom = nullptr;

	Int_t Centrality;
	TClonesArray* CloneArr_mu = nullptr;
	Short_t Gen_QQ_size;
	Short_t Gen_QQ_sign[1000];
	Short_t Gen_QQ_mupl_idx[1000];
	Short_t Gen_QQ_mumi_idx[1000];
	Short_t Gen_QQ_whichRec[1000];

	// event variables
	OniaTree->SetBranchAddress("Gen_weight", &Gen_weight);
	OniaTree->SetBranchAddress("Centrality", &Centrality);

	// gen-level variables
	OniaTree->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size);
	OniaTree->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom);

	OniaTree->SetBranchAddress("Gen_mu_4mom", &CloneArr_mu);
	OniaTree->SetBranchAddress("Gen_QQ_mupl_idx", Gen_QQ_mupl_idx);
	OniaTree->SetBranchAddress("Gen_QQ_mumi_idx", Gen_QQ_mumi_idx);

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

	Int_t nPhiBins = 25;
	Float_t phiMin = -180, phiMax = 280;

	//Double_t phiBinning[] = {-TMath::Pi(), -2.5, -2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2, 2.5, 180};
	//	nPhiBins = sizeof(cosThetaBinning) / sizeof(Double_t) - 1;

	TEfficiency* hCS = new TEfficiency("hCS", ";cos #theta_{CS}; #varphi_{CS} (#circ);acceptance", nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);
	TEfficiency* hHX = new TEfficiency("hHX", ";cos #theta_{HX}; #varphi_{HX} (#circ);acceptance", nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);

	TLorentzVector* genLorentzVector = new TLorentzVector();

	double weight;
	Bool_t withinAcceptance;

	Long64_t totEntries = OniaTree->GetEntries();

	for (Long64_t iEvent = 0; iEvent < (totEntries); iEvent++) {
		if (iEvent % 10000 == 0) {
			cout << Form("\rProcessing event %lld / %lld (%.0f%%)", iEvent, totEntries, 100. * iEvent / totEntries) << flush;
		}

		OniaTree->GetEntry(iEvent);
		//genLorentzVector->Clear();

		// event selection

		//if (Centrality >= 2 * 90) continue; // discard events with centrality > 90% in 2018 data

		// get N_coll from HF
		//Int_t hiBin = GetHiBinFromhiHF(SumET_HF);

		weight = Gen_weight;

		// loop over all gen upsilons
		for (int iGen = 0; iGen < Gen_QQ_size; iGen++) {
			genLorentzVector = (TLorentzVector*)Gen_QQ_4mom->At(iGen);

			if (genLorentzVector->Pt() < ptMin || genLorentzVector->Pt() > ptMax) continue; // pt bin of interest

			if (fabs(genLorentzVector->Rapidity()) > 2.4) continue; // upsilon within acceptance

			/// single-muon selection criteria
			int iMuPlus = Gen_QQ_mupl_idx[iGen];
			int iMuMinus = Gen_QQ_mumi_idx[iGen];

			// acceptance
			TLorentzVector* Gen_mupl_4mom = (TLorentzVector*)CloneArr_mu->At(iMuPlus);

			TLorentzVector* Gen_mumi_4mom = (TLorentzVector*)CloneArr_mu->At(iMuMinus);

			withinAcceptance = (fabs(Gen_mupl_4mom->Eta()) < 2.4) && (Gen_mupl_4mom->Pt() > 3.5) && (fabs(Gen_mumi_4mom->Eta()) < 2.4) && (Gen_mumi_4mom->Pt() > 3.5);

			// Transformations and rotations

			Gen_QQ_phi = genLorentzVector->Phi();
			Gen_QQ_costheta = genLorentzVector->CosTheta();
			Gen_QQ_pt = genLorentzVector->Pt();
			Gen_QQ_y = genLorentzVector->Rapidity();
			Gen_QQ_px = genLorentzVector->Px();
			Gen_QQ_py = genLorentzVector->Py();
			Gen_QQ_pz = genLorentzVector->Pz();
			Gen_QQ_E = genLorentzVector->Energy();
			Gen_QQ_m = genLorentzVector->M();

			Gen_mupl_phi = Gen_mupl_4mom->Phi();
			Gen_mupl_costheta = Gen_mupl_4mom->CosTheta();
			Gen_mupl_pt = Gen_mupl_4mom->Pt();
			Gen_mupl_y = Gen_mupl_4mom->Rapidity();
			Gen_mupl_eta = Gen_mupl_4mom->Eta();
			Gen_mupl_px = Gen_mupl_4mom->Px();
			Gen_mupl_py = Gen_mupl_4mom->Py();
			Gen_mupl_pz = Gen_mupl_4mom->Pz();
			Gen_mupl_E = Gen_mupl_4mom->Energy();

			Gen_mumi_phi = Gen_mumi_4mom->Phi();
			Gen_mumi_costheta = Gen_mumi_4mom->CosTheta();
			Gen_mumi_pt = Gen_mumi_4mom->Pt();
			Gen_mumi_y = Gen_mumi_4mom->Rapidity();
			Gen_mumi_eta = Gen_mumi_4mom->Eta();
			Gen_mumi_px = Gen_mumi_4mom->Px();
			Gen_mumi_py = Gen_mumi_4mom->Py();
			Gen_mumi_pz = Gen_mumi_4mom->Pz();
			Gen_mumi_E = Gen_mumi_4mom->Energy();

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

			// cout << "<<In the lab frame>>" << endl;
			// cout << "ups: p = (" << upsPvecLab.Px() << ", " << upsPvecLab.Py()  << ", " << upsPvecLab.Pz() << ")" << endl;
			// cout << "mu+: p = (" << muplPvecLab.Px() << ", " << muplPvecLab.Py()  << ", " << muplPvecLab.Pz() << ")" << endl;
			// cout << "mu-: p = (" << mumiPvecLab.Px() << ", " << mumiPvecLab.Py()  << ", " << mumiPvecLab.Pz() << ")" << endl;

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

			// ******** Print out momentums of daughter muons in the upsilon's rest frame after coordinate rotation ******** //
			// cout << endl;
			// cout << "<<Rotated the quarkonium rest frame>>" << endl;
			// cout << "mu+: p = (" << muplPvecBoosted.Px() << ", " << muplPvecBoosted.Py()  << ", " << muplPvecBoosted.Pz() << ")" << endl;
			// cout << "mu-: p = (" << mumiPvecBoosted.Px() << ", " << mumiPvecBoosted.Py()  << ", " << mumiPvecBoosted.Pz() << ")" << endl;

			TLorentzVector mupl4MomBoostedRot(muplPvecBoosted, mupl4MomBoosted.E());

			// ******** HX to CS (rotation from HX frame to CS frame) ******** //
			// (1. Boost two beams to upsilon's rest frame)
			// (2. Rotate the coordinate)
			// (3. Get angle between two beams(b1 and -b2), and between b1 and ZHX in the upsilon's rest frame)
			// (4. Calculate delta (angle btw ZHX and ZCS))

			// ******** Transform variables of beams from the lab frame to the upsilon's rest frame ******** //
			TLorentzVector beam14MomBoosted(beam1PvecLab, beam1_E);
			TLorentzVector beam24MomBoosted(beam2PvecLab, beam2_E);

			// ups4MomLab.SetX(-1);
			// ups4MomLab.SetY(0);
			// ups4MomLab.SetZ(0);

			beam14MomBoosted.Boost(-ups4MomLab.BoostVector());
			beam24MomBoosted.Boost(-ups4MomLab.BoostVector());

			// ******** Rotate the coordinate ******** //
			TVector3 beam1PvecBoosted(beam14MomBoosted.Px(), beam14MomBoosted.Py(), beam14MomBoosted.Pz());
			TVector3 beam2PvecBoosted(beam24MomBoosted.Px(), beam24MomBoosted.Py(), beam24MomBoosted.Pz());

			// upsPvecLab.SetX(-1);
			// upsPvecLab.SetY(0);
			// upsPvecLab.SetZ(0);

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

			// // (The math for caculating the angle between z_HX and z_CS is different depending on the sign of the beam1's z-coordinate)
			// if(beam1PvecBoosted.Pz()>0) delta = Angle_B1ZHX + Angle_B1miB2/2.;
			// else if(beam1PvecBoosted.Pz()<0) delta = Angle_B1ZHX - Angle_B1miB2/2.;
			// else cout <<  "beam1PvecBoosted.Pz() = 0?" << endl;
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

			hCS->FillWeighted(withinAcceptance, weight, muplPvecBoostedCS.CosTheta(), muplPvecBoostedCS.Phi() * 180 / TMath::Pi());
			hHX->FillWeighted(withinAcceptance, weight, muplPvecBoosted.CosTheta(), muplPvecBoosted.Phi() * 180 / TMath::Pi());
		}
	}

	gStyle->SetPadLeftMargin(.15);
	//gStyle->SetTitleYOffset(.9);
	gStyle->SetPadRightMargin(0.18);
	gStyle->SetPalette(kRainBow);

	TLatex* legend = new TLatex();
	legend->SetTextAlign(22);
	legend->SetTextSize(0.05);

	auto* canvasCS = new TCanvas("canvasCS", "", 700, 600);
	hCS->Draw("COLZ");

	CMS_lumi(canvasCS, Form("#varUpsilon(%dS) Hydjet-embedded PbPb MC", iState));

	legend->DrawLatexNDC(.48, .88, Form("|y^{#mu#mu}| < 2.4, %d < p_{T}^{#mu#mu} < %d GeV", ptMin, ptMax));
	legend->DrawLatexNDC(.48, .8, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, p_{T}^{#mu} > 3.5 GeV", iState));

	gPad->Update();

	hCS->GetPaintedHistogram()->GetYaxis()->SetRangeUser(-190, 300);
	hCS->GetPaintedHistogram()->GetZaxis()->SetRangeUser(0, 1);

	canvasCS->SaveAs(Form("AcceptanceMaps/%dS/CS_cent%dto%d_pt%dto%dGeV.png", iState, centMin, centMax, ptMin, ptMax), "RECREATE");

	auto* canvasHX = new TCanvas("canvasHX", "", 700, 600);
	hHX->Draw("COLZ");

	CMS_lumi(canvasHX, Form("#varUpsilon(%dS) Hydjet-embedded PbPb MC", iState));

	legend->DrawLatexNDC(.48, .88, Form("|y^{#mu#mu}| < 2.4, %d < p_{T}^{#mu#mu} < %d GeV", ptMin, ptMax));
	legend->DrawLatexNDC(.48, .8, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, p_{T}^{#mu} > 3.5 GeV", iState));

	gPad->Update();

	hHX->GetPaintedHistogram()->GetYaxis()->SetRangeUser(-190, 300);
	hHX->GetPaintedHistogram()->GetZaxis()->SetRangeUser(0, 1);

	canvasHX->SaveAs(Form("AcceptanceMaps/%dS/HX_cent%dto%d_pt%dto%dGeV.png", iState, centMin, centMax, ptMin, ptMax), "RECREATE");
}
