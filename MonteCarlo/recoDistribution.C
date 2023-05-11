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

#include "../Tools/Parameters/AnalysisParameters.h"

#include "../Tools/Parameters/CentralityValues.h"

#include "../Tools/Parameters/MuonScaleFactors.h"

void recoDistribution(Int_t ptMin = 0, Int_t ptMax = 30) {
	const char* filename = "../Files/OniaTree_Y1S_miniAOD_HydjetEmbeddedMC.root";
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

	Float_t SumET_HF;

	TClonesArray* Reco_QQ_4mom = nullptr;

	ULong64_t HLTriggers;
	ULong64_t Reco_QQ_trig[1000];
	Int_t Centrality;
	TClonesArray* CloneArr_QQ = nullptr;
	TClonesArray* CloneArr_mu = nullptr;
	Short_t Reco_QQ_size;
	Short_t Reco_QQ_sign[1000];
	Short_t Reco_QQ_mupl_idx[1000];
	Short_t Reco_QQ_mumi_idx[1000];
	Short_t Reco_QQ_whichGen[1000];

	ULong64_t Reco_mu_trig[1000];
	Int_t Reco_mu_SelectionType[1000];
	//(parameters for quality cuts)
	Float_t Reco_QQ_VtxProb[1000];
	Int_t Reco_mu_nPixWMea[1000];
	Int_t Reco_mu_nTrkWMea[1000];
	Float_t Reco_mu_dxy[1000];
	Float_t Reco_mu_dz[1000];

	// event variables
	OniaTree->SetBranchAddress("SumET_HF", &SumET_HF);
	OniaTree->SetBranchAddress("HLTriggers", &HLTriggers);
	OniaTree->SetBranchAddress("Centrality", &Centrality);

	OniaTree->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig);
	OniaTree->SetBranchAddress("Reco_QQ_whichGen", Reco_QQ_whichGen);
	OniaTree->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom);
	OniaTree->SetBranchAddress("Reco_mu_4mom", &CloneArr_mu);
	OniaTree->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size);
	OniaTree->SetBranchAddress("Reco_QQ_sign", &Reco_QQ_sign);
	OniaTree->SetBranchAddress("Reco_QQ_mupl_idx", Reco_QQ_mupl_idx);
	OniaTree->SetBranchAddress("Reco_QQ_mumi_idx", Reco_QQ_mumi_idx);

	OniaTree->SetBranchAddress("Reco_mu_trig", Reco_mu_trig);
	OniaTree->SetBranchAddress("Reco_mu_SelectionType", Reco_mu_SelectionType);

	OniaTree->SetBranchAddress("Reco_QQ_VtxProb", &Reco_QQ_VtxProb);
	OniaTree->SetBranchAddress("Reco_mu_nPixWMea", &Reco_mu_nPixWMea);
	OniaTree->SetBranchAddress("Reco_mu_nTrkWMea", &Reco_mu_nTrkWMea);
	OniaTree->SetBranchAddress("Reco_mu_dxy", &Reco_mu_dxy);
	OniaTree->SetBranchAddress("Reco_mu_dz", &Reco_mu_dz);

	Float_t Reco_QQ_phi, Reco_QQ_costheta, Reco_QQ_pt, Reco_QQ_y, Reco_QQ_px, Reco_QQ_py, Reco_QQ_pz, Reco_QQ_E, Reco_QQ_m;
	Float_t Reco_mupl_phi, Reco_mupl_costheta, Reco_mupl_pt, Reco_mupl_y, Reco_mupl_eta, Reco_mupl_px, Reco_mupl_py, Reco_mupl_pz, Reco_mupl_E;
	Float_t Reco_mumi_phi, Reco_mumi_costheta, Reco_mumi_pt, Reco_mumi_y, Reco_mumi_eta, Reco_mumi_px, Reco_mumi_py, Reco_mumi_pz, Reco_mumi_E;

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

	Int_t nPhiBins = 21;
	Float_t phiMin = -180, phiMax = 240;

	//Double_t phiBinning[] = {-TMath::Pi(), -2.5, -2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2, 2.5, 180};
	//	nPhiBins = sizeof(cosThetaBinning) / sizeof(Double_t) - 1;

	TH2F* hCS = new TH2F("hCS", ";cos #theta_{CS}; #varphi_{CS} (#circ);# reco QQ", nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);
	TH2F* hHX = new TH2F("hHX", ";cos #theta_{HX}; #varphi_{HX} (#circ);# reco QQ", nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);

	TLorentzVector* recoLV = new TLorentzVector();

	double weight, muTrkSF, muIdSF, muTrigSF;
	Bool_t allGood, firesTrigger, isRecoMatched, goodVertexProba, withinAcceptance, trackerAndGlobalMuons, hybridSoftMuons;

	Long64_t totEntries = OniaTree->GetEntries();

	for (Long64_t iEvent = 0; iEvent < (totEntries); iEvent++) {
		if (iEvent % 10000 == 0) {
			cout << Form("\rProcessing event %lld / %lld (%.0f%%)", iEvent, totEntries, 100. * iEvent / totEntries) << flush;
		}

		OniaTree->GetEntry(iEvent);
		//recoLV->Clear();

		// event selection

		if (Centrality >= 2 * 90) continue; // discard events with centrality >= 90% in 2018 data

		firesTrigger = ((HLTriggers & (ULong64_t)(1 << (UpsilonHLTBit - 1))) == (ULong64_t)(1 << (UpsilonHLTBit - 1)));

		if (!firesTrigger) continue;

		// get N_coll from HF
		Int_t hiBin = GetHiBinFromhiHF(SumET_HF);

		weight = FindNcoll(hiBin);

		// loop over all reco QQ candidates
		for (int iReco = 0; iReco < Reco_QQ_size; iReco++) {
			recoLV = (TLorentzVector*)Reco_QQ_4mom->At(iReco);

			if (recoLV->Pt() < ptMin || recoLV->Pt() > ptMax) continue; // pt bin of interest

			if (fabs(recoLV->Rapidity()) > 2.4 || fabs(recoLV->Rapidity()) < 1.2) continue; // upsilon within acceptance

			if (Reco_QQ_whichGen[iReco] < 0) continue; //matching gen

			if (Reco_QQ_sign[iReco] != 0) continue; // only opposite-sign muon pairs

			if (Reco_QQ_VtxProb[iReco] < 0.01) continue;

			/// single-muon selection criteria
			int iMuPlus = Reco_QQ_mupl_idx[iReco];
			int iMuMinus = Reco_QQ_mumi_idx[iReco];

			// global AND tracker muons
			trackerAndGlobalMuons = (Reco_mu_SelectionType[iMuPlus] & 2) && (Reco_mu_SelectionType[iMuPlus] & 8) && (Reco_mu_SelectionType[iMuMinus] & 2) && (Reco_mu_SelectionType[iMuMinus] & 8);

			if (!trackerAndGlobalMuons) continue;

			// passing hybrid-soft Id
			hybridSoftMuons = (Reco_mu_nTrkWMea[iMuPlus] > 5) && (Reco_mu_nPixWMea[iMuPlus] > 0) && (fabs(Reco_mu_dxy[iMuPlus]) < 0.3) && (fabs(Reco_mu_dz[iMuPlus]) < 20.) && (Reco_mu_nTrkWMea[iMuMinus] > 5) && (Reco_mu_nPixWMea[iMuMinus] > 0) && (fabs(Reco_mu_dxy[iMuMinus]) < 0.3) && (fabs(Reco_mu_dz[iMuMinus]) < 20.);

			if (!hybridSoftMuons) continue;

			// acceptance
			TLorentzVector* Reco_mupl_4mom = (TLorentzVector*)CloneArr_mu->At(iMuPlus);

			TLorentzVector* Reco_mumi_4mom = (TLorentzVector*)CloneArr_mu->At(iMuMinus);

			withinAcceptance = (1.2 < fabs(Reco_mupl_4mom->Eta()) < 2.4) && (Reco_mupl_4mom->Pt() > 3.) && (1.2 < fabs(Reco_mumi_4mom->Eta()) < 2.4) && (Reco_mumi_4mom->Pt() > 3.);

			if (!withinAcceptance) continue;

			// Transformations and rotations

			Reco_QQ_phi = recoLV->Phi();
			Reco_QQ_costheta = recoLV->CosTheta();
			Reco_QQ_pt = recoLV->Pt();
			Reco_QQ_y = recoLV->Rapidity();
			Reco_QQ_px = recoLV->Px();
			Reco_QQ_py = recoLV->Py();
			Reco_QQ_pz = recoLV->Pz();
			Reco_QQ_E = recoLV->Energy();
			Reco_QQ_m = recoLV->M();

			Reco_mupl_phi = Reco_mupl_4mom->Phi();
			Reco_mupl_costheta = Reco_mupl_4mom->CosTheta();
			Reco_mupl_pt = Reco_mupl_4mom->Pt();
			Reco_mupl_y = Reco_mupl_4mom->Rapidity();
			Reco_mupl_eta = Reco_mupl_4mom->Eta();
			Reco_mupl_px = Reco_mupl_4mom->Px();
			Reco_mupl_py = Reco_mupl_4mom->Py();
			Reco_mupl_pz = Reco_mupl_4mom->Pz();
			Reco_mupl_E = Reco_mupl_4mom->Energy();

			Reco_mumi_phi = Reco_mumi_4mom->Phi();
			Reco_mumi_costheta = Reco_mumi_4mom->CosTheta();
			Reco_mumi_pt = Reco_mumi_4mom->Pt();
			Reco_mumi_y = Reco_mumi_4mom->Rapidity();
			Reco_mumi_eta = Reco_mumi_4mom->Eta();
			Reco_mumi_px = Reco_mumi_4mom->Px();
			Reco_mumi_py = Reco_mumi_4mom->Py();
			Reco_mumi_pz = Reco_mumi_4mom->Pz();
			Reco_mumi_E = Reco_mumi_4mom->Energy();

			// ******** Construct 4-momentum vector of upsilon and muons (Lab Frame) ******** //
			// (documetation of TVector3 and TLorentzVector: https://root.cern.ch/root/html534/guides/users-guide/PhysicsVectors.html#lorentz-boost)
			TVector3 upsPvecLab(Reco_QQ_px, Reco_QQ_py, Reco_QQ_pz);
			TLorentzVector ups4MomLab(upsPvecLab, Reco_QQ_E);

			TVector3 muplPvecLab(Reco_mupl_px, Reco_mupl_py, Reco_mupl_pz);
			TLorentzVector mupl4MomLab(muplPvecLab, Reco_mupl_E);

			TVector3 mumiPvecLab(Reco_mumi_px, Reco_mumi_py, Reco_mumi_pz);
			TLorentzVector mumi4MomLab(mumiPvecLab, Reco_mumi_E);

			TVector3 beam1PvecLab(0, 0, beam1_p);
			TLorentzVector beam14MomLab(beam1PvecLab, beam1_E);

			TVector3 beam2PvecLab(0, 0, beam2_p);
			TLorentzVector beam24MomLab(beam2PvecLab, beam2_E);

			// cout << "<<In the lab frame>>" << endl;
			// cout << "ups: p = (" << upsPvecLab.Px() << ", " << upsPvecLab.Py()  << ", " << upsPvecLab.Pz() << ")" << endl;
			// cout << "mu+: p = (" << muplPvecLab.Px() << ", " << muplPvecLab.Py()  << ", " << muplPvecLab.Pz() << ")" << endl;
			// cout << "mu-: p = (" << mumiPvecLab.Px() << ", " << mumiPvecLab.Py()  << ", " << mumiPvecLab.Pz() << ")" << endl;

			// ******** Transform variables of muons from the lab frame to the upsilon's rest frame ******** //
			TLorentzVector ups4MomBoosted(upsPvecLab, Reco_QQ_E);
			TLorentzVector mupl4MomBoosted(muplPvecLab, Reco_mupl_E);
			TLorentzVector mumi4MomBoosted(mumiPvecLab, Reco_mumi_E);

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

			hCS->Fill(muplPvecBoostedCS.CosTheta(), muplPvecBoostedCS.Phi() * 180 / TMath::Pi(), weight);
			hHX->Fill(muplPvecBoosted.CosTheta(), muplPvecBoosted.Phi() * 180 / TMath::Pi(), weight);
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

	CMS_lumi(canvasCS, "#varUpsilon(1S) Hydjet-embedded PbPb MC");

	legend->DrawLatexNDC(.5, .88, Form("centrality %d-%d%%, %d < p_{T}^{#mu#mu} < %d GeV", centMin, centMax, ptMin, ptMax));

	gPad->Update();

	hCS->GetYaxis()->SetRangeUser(-190, 240);
	//	hCS->GetPaintedHistogram()->GetZaxis()->SetRangeUser(0, 0.9);

	canvasCS->SaveAs(Form("RecoMaps/CS_cent%dto%d_pt%dto%dGeV.png", centMin, centMax, ptMin, ptMax), "RECREATE");

	auto* canvasHX = new TCanvas("canvasHX", "", 700, 600);
	hHX->Draw("COLZ");

	CMS_lumi(canvasHX, "#varUpsilon(1S) Hydjet-embedded PbPb MC");

	legend->DrawLatexNDC(.5, .88, Form("centrality %d-%d%%, %d < p_{T}^{#mu#mu} < %d GeV", centMin, centMax, ptMin, ptMax));

	gPad->Update();

	hHX->GetYaxis()->SetRangeUser(-190, 240);
	//	hHX->GetPaintedHistogram()->GetZaxis()->SetRangeUser(0, 0.9);

	canvasHX->SaveAs(Form("RecoMaps/HX_cent%dto%d_pt%dto%dGeV.png", centMin, centMax, ptMin, ptMax), "RECREATE");
}
