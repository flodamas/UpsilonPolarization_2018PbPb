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

#include "../AnalysisParameters.h"

// (https://twiki.cern.ch/twiki/bin/viewauth/CMS/UpsilonPolarizationInPbPb5TeV)

void skimUpsilonCandidates(const char* inputFileName = "OniaTree_miniAOD_PbPbPrompt_112X_DATA_ep.root", const char* outputFileName = "upsilonSkimmedDataset.root") {
	TFile* infile = TFile::Open(inputFileName, "READ");
	TTree* OniaTree = (TTree*)infile->Get("hionia/myTree");

	/// OniaTree variables
	Float_t zVtx;
	ULong64_t HLTriggers;
	ULong64_t Reco_QQ_trig[1000];
	Int_t Centrality;
	TClonesArray* CloneArr_QQ = nullptr;
	TClonesArray* CloneArr_mu = nullptr;
	Short_t Reco_QQ_size;
	Short_t Reco_QQ_sign[1000];
	Short_t Reco_QQ_mupl_idx[1000];
	Short_t Reco_QQ_mumi_idx[1000];

	Int_t Reco_mu_SelectionType[1000];
	//(parameters for quality cuts)
	Float_t Reco_QQ_VtxProb[1000];
	Int_t Reco_mu_nPixWMea[1000];
	Int_t Reco_mu_nTrkWMea[1000];
	Float_t Reco_mu_dxy[1000];
	Float_t Reco_mu_dz[1000];

	OniaTree->SetBranchAddress("zVtx", &zVtx);

	OniaTree->SetBranchAddress("HLTriggers", &HLTriggers);
	OniaTree->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig);
	OniaTree->SetBranchAddress("Centrality", &Centrality);
	OniaTree->SetBranchAddress("Reco_QQ_4mom", &CloneArr_QQ);
	OniaTree->SetBranchAddress("Reco_mu_4mom", &CloneArr_mu);
	OniaTree->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size);
	OniaTree->SetBranchAddress("Reco_QQ_sign", &Reco_QQ_sign);
	OniaTree->SetBranchAddress("Reco_QQ_mupl_idx", Reco_QQ_mupl_idx);
	OniaTree->SetBranchAddress("Reco_QQ_mumi_idx", Reco_QQ_mumi_idx);
	OniaTree->SetBranchAddress("Reco_mu_SelectionType", Reco_mu_SelectionType);

	OniaTree->SetBranchAddress("Reco_QQ_VtxProb", &Reco_QQ_VtxProb);
	OniaTree->SetBranchAddress("Reco_mu_nPixWMea", &Reco_mu_nPixWMea);
	OniaTree->SetBranchAddress("Reco_mu_nTrkWMea", &Reco_mu_nTrkWMea);
	OniaTree->SetBranchAddress("Reco_mu_dxy", &Reco_mu_dxy);
	OniaTree->SetBranchAddress("Reco_mu_dz", &Reco_mu_dz);

	double Reco_QQ_E, Reco_mupl_E, Reco_mumi_E;

	/// RooDataSet output: one entry = one dimuon candidate!
	RooRealVar centVar("centrality", "event centrality", 0, 200);

	Float_t lowMassCut = 8, highMassCut = 14;
	RooRealVar massVar("mass", "m_{#mu^{#plus}#mu^{#minus}}", lowMassCut, highMassCut, "GeV/c^{2}");
	RooRealVar yVar("rapidity", "dimuon absolute rapidity", gRapidityMin, gRapidityMax);
	RooRealVar ptVar("pt", "dimuon pT", 0, 100, "GeV/c");

	RooRealVar cosThetaLabVar("cosThetaLab", "cos theta in the lab frame", -1, 1);
	RooRealVar phiLabVar("phiLab", "phi angle in the lab frame", -180, 180, "#circ");

	RooRealVar cosThetaCSVar("cosThetaCS", "cos theta in the Collins-Soper frame", -1, 1);
	RooRealVar phiCSVar("phiCS", "phi angle in the Collins-Soper frame", -180, 180, "#circ");

	RooRealVar cosThetaHXVar("cosThetaHX", "cos theta in the helicity frame", -1, 1);
	RooRealVar phiHXVar("phiHX", "phi angle in the helicity frame", -180, 180, "#circ");

	RooDataSet dataset("dataset", "skimmed dataset", RooArgSet(centVar, massVar, yVar, ptVar, cosThetaLabVar, phiLabVar, cosThetaCSVar, phiCSVar, cosThetaHXVar, phiHXVar));

	// ******** Set beam energy for the Collins-Soper reference frame ******** //
	double sqrt_S_NN = 5.02;                 //(Center of mass Energy per nucleon pair in TeV)
	double beam1_p = sqrt_S_NN * 1000. / 2.; //(in GeV) (Note. sqrt_S_NN = sqrt(2*E1*E2+2*p1*p2) = 2E1 when two beams have the same E)
	double beam1_E = beam1_p;
	double beam2_p = -beam1_p;
	double beam2_E = beam1_E;
	double delta = 0; //(Angle between ZHX(Z-axis in the Helicity frame) and ZCS(Z-axis in the Collins-Soper frame))

	Long64_t totEntries = OniaTree->GetEntries();

	for (Long64_t iEvent = 0; iEvent < (totEntries); iEvent++) {
		if (iEvent % 10000 == 0) {
			cout << Form("\rProcessing event %lld / %lld (%.0f%%)", iEvent, totEntries, 100. * iEvent / totEntries) << flush;
		}

		OniaTree->GetEntry(iEvent);

		// event selection
		//if (fabs(zVtx) > 15) continue;

		if (Centrality >= 2 * gCentralityBinMax) continue; // discard events with centrality >= 90% in 2018 data

		if (!((HLTriggers & (ULong64_t)(1 << (gUpsilonHLTBit - 1))) == (ULong64_t)(1 << (gUpsilonHLTBit - 1)))) continue; // must fire the upsilon HLT path

		// loop over reconstructed dimuon candidates
		for (int iQQ = 0; iQQ < Reco_QQ_size; iQQ++) {
			if (!((Reco_QQ_trig[iQQ] & (ULong64_t)(1 << (gUpsilonHLTBit - 1))) == (ULong64_t)(1 << (gUpsilonHLTBit - 1)))) continue; // dimuon matching

			if (Reco_QQ_sign[iQQ] != 0) continue; // only opposite-sign muon pairs

			if (Reco_QQ_VtxProb[iQQ] < 0.01) continue; // good common vertex proba

			TLorentzVector* Reco_QQ_4mom = (TLorentzVector*)CloneArr_QQ->At(iQQ);

			if (Reco_QQ_4mom->M() < lowMassCut || Reco_QQ_4mom->M() > highMassCut) continue; // speedup!

			if (fabs(Reco_QQ_4mom->Rapidity()) < gRapidityMin || fabs(Reco_QQ_4mom->Rapidity()) > gRapidityMax) continue;

			//if (Reco_QQ_4mom->Pt() > 50) continue;

			/// single-muon selection criteria
			int iMuPlus = Reco_QQ_mupl_idx[iQQ];
			int iMuMinus = Reco_QQ_mumi_idx[iQQ];

			// global AND tracker muons
			if (!((Reco_mu_SelectionType[iMuPlus] & 2) && (Reco_mu_SelectionType[iMuPlus] & 8))) continue;
			if (!((Reco_mu_SelectionType[iMuMinus] & 2) && (Reco_mu_SelectionType[iMuMinus] & 8))) continue;

			// passing hybrid-soft Id
			if (!((Reco_mu_nTrkWMea[iMuPlus] > 5) && (Reco_mu_nPixWMea[iMuPlus] > 0) && (fabs(Reco_mu_dxy[iMuPlus]) < 0.3) && (fabs(Reco_mu_dz[iMuPlus]) < 20.))) continue;
			if (!((Reco_mu_nTrkWMea[iMuMinus] > 5) && (Reco_mu_nPixWMea[iMuMinus] > 0) && (fabs(Reco_mu_dxy[iMuMinus]) < 0.3) && (fabs(Reco_mu_dz[iMuMinus]) < 20.))) continue;

			// acceptance

			TLorentzVector* Reco_mupl_4mom = (TLorentzVector*)CloneArr_mu->At(iMuPlus);

			if (fabs(Reco_mupl_4mom->Eta()) > 2.4) continue;
			if (Reco_mupl_4mom->Pt() < 3.5) continue;

			TLorentzVector* Reco_mumi_4mom = (TLorentzVector*)CloneArr_mu->At(iMuMinus);

			if (fabs(Reco_mumi_4mom->Eta()) > 2.4) continue;
			if (Reco_mumi_4mom->Pt() < 3.5) continue;

			// ******** Store kinematics of upsilon and muons (Lab Frame) into variables ******** //

			Reco_QQ_E = Reco_QQ_4mom->Energy();

			Reco_mupl_E = Reco_mupl_4mom->Energy();

			Reco_mumi_E = Reco_mumi_4mom->Energy();

			// ******** Construct 4-momentum vector of upsilon and muons (Lab Frame) ******** //
			// (documetation of TVector3 and TLorentzVector: https://root.cern.ch/root/html534/guides/users-guide/PhysicsVectors.html#lorentz-boost)
			TVector3 upsPvecLab(Reco_QQ_4mom->Px(), Reco_QQ_4mom->Py(), Reco_QQ_4mom->Pz());
			TLorentzVector ups4MomLab(upsPvecLab, Reco_QQ_E);

			TVector3 muplPvecLab(Reco_mupl_4mom->Px(), Reco_mupl_4mom->Py(), Reco_mupl_4mom->Pz());
			TLorentzVector mupl4MomLab(muplPvecLab, Reco_mupl_E);

			TVector3 mumiPvecLab(Reco_mumi_4mom->Px(), Reco_mumi_4mom->Py(), Reco_mumi_4mom->Pz());
			TLorentzVector mumi4MomLab(mumiPvecLab, Reco_mumi_E);

			TVector3 beam1PvecLab(0, 0, beam1_p);
			TLorentzVector beam14MomLab(beam1PvecLab, beam1_E);

			TVector3 beam2PvecLab(0, 0, beam2_p);
			TLorentzVector beam24MomLab(beam2PvecLab, beam2_E);

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

			delta = 0; // reset the angle between z_HX and z_CS

			// // (The math for caculating the angle between z_HX and z_CS is different depending on the sign of the beam1's z-coordinate)
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

			// fill the dataset

			centVar = Centrality;

			massVar = Reco_QQ_4mom->M();
			yVar = fabs(Reco_QQ_4mom->Rapidity());
			ptVar = Reco_QQ_4mom->Pt();

			cosThetaLabVar = Reco_mupl_4mom->CosTheta();
			phiLabVar = Reco_mupl_4mom->Phi() * 180 / TMath::Pi();

			cosThetaCSVar = muplPvecBoostedCS.CosTheta();
			phiCSVar = muplPvecBoostedCS.Phi() * 180 / TMath::Pi();

			cosThetaHXVar = muplPvecBoosted.CosTheta();
			phiHXVar = muplPvecBoosted.Phi() * 180 / TMath::Pi();

			dataset.add(RooArgSet(centVar, massVar, yVar, ptVar, cosThetaLabVar, phiLabVar, cosThetaCSVar, phiCSVar, cosThetaHXVar, phiHXVar));

		} // end of reco QQ loop
	}   // enf of event loop

	TFile file(outputFileName, "RECREATE");

	dataset.Write();

	file.Close();

	return;
}
