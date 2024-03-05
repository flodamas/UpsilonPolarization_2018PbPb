#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../ReferenceFrameTransform/Transformations.h"

// (https://twiki.cern.ch/twiki/bin/viewauth/CMS/UpsilonPolarizationInPbPb5TeV)

void skimUpsilonCandidates(const char* inputFileName = "OniaTree_miniAOD_PbPbPrompt_112X_DATA_ep.root", const char* outputFileName = "UpsilonSkimmedDataset.root") {
	/// Open OniaTree
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

	/// RooDataSet output: one entry = one dimuon candidate!
	RooRealVar centVar("centrality", "event centrality", 0, 200);

	Float_t lowMassCut = 6.5, highMassCut = 14.5; // large invariant mass window to ease the integration for RooFit
	RooRealVar massVar("mass", "m_{#mu^{#plus}#mu^{#minus}}", lowMassCut, highMassCut, "GeV/c^{2}");
	RooRealVar yVar("rapidity", "dimuon absolute rapidity", 0, 2.4);
	RooRealVar ptVar("pt", "dimuon pT", 0, gPtBinning[NPtBins], "GeV/c");

	RooRealVar cosThetaLabVar("cosThetaLab", "cos theta in the lab frame", -1, 1);
	RooRealVar phiLabVar("phiLab", "phi angle in the lab frame", -180, 180, "#circ");
	RooRealVar etaLabMuplVar("etaLabMupl", "eta of positive muon in the lab frame", -2.4, 2.4);
	RooRealVar etaLabMumiVar("etaLabMumi", "eta of negative muon in the lab frame", -2.4, 2.4);

	RooDataSet datasetLab("datasetLab", "skimmed dataset for the Lab frame", RooArgSet(centVar, massVar, yVar, ptVar, cosThetaLabVar, phiLabVar, etaLabMuplVar, etaLabMumiVar));

	RooRealVar cosThetaCSVar("cosThetaCS", "cos theta in the Collins-Soper frame", -1, 1);
	RooRealVar phiCSVar("phiCS", "phi angle in the Collins-Soper frame", 0, 180, "#circ");

	RooDataSet datasetCS("datasetCS", "skimmed dataset for the CS frame", RooArgSet(centVar, massVar, yVar, ptVar, cosThetaCSVar, phiCSVar));

	RooRealVar cosThetaHXVar("cosThetaHX", "cos theta in the helicity frame", -1, 1);
	RooRealVar phiHXVar("phiHX", "phi angle in the helicity frame", -180, 180, "#circ");

	RooDataSet datasetHX("datasetHX", "skimmed dataset for the HX frame", RooArgSet(centVar, massVar, yVar, ptVar, cosThetaHXVar, phiHXVar));

	// loop variables
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

			if (fabs(Reco_QQ_4mom->Rapidity()) > 2.4) continue; // cut on the dimuon rapidity when reducing the dataset later

			if (Reco_QQ_4mom->Pt() > gPtBinning[NPtBins]) continue;

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

			// get positive muon's coordinates in the studied reference frames
			TVector3 muPlus_CS = MuPlusVector_CollinsSoper(*Reco_QQ_4mom, *Reco_mupl_4mom);

			TVector3 muPlus_HX = MuPlusVector_Helicity(*Reco_QQ_4mom, *Reco_mupl_4mom);

			// fill the datasets

			centVar = Centrality;

			massVar = Reco_QQ_4mom->M();
			yVar = fabs(Reco_QQ_4mom->Rapidity());
			ptVar = Reco_QQ_4mom->Pt();

			// Lab
			cosThetaLabVar = Reco_mupl_4mom->CosTheta();
			phiLabVar = Reco_mupl_4mom->Phi() * 180 / TMath::Pi();
			etaLabMuplVar = Reco_mupl_4mom->PseudoRapidity();
			etaLabMumiVar = Reco_mumi_4mom->PseudoRapidity();

			datasetLab.add(RooArgSet(centVar, massVar, yVar, ptVar, cosThetaLabVar, phiLabVar, etaLabMuplVar, etaLabMumiVar));

			// Collins-Soper
			cosThetaCSVar = muPlus_CS.CosTheta();
			phiCSVar = muPlus_CS.Phi() * 180 / TMath::Pi();

			datasetCS.add(RooArgSet(centVar, massVar, yVar, ptVar, cosThetaCSVar, phiCSVar));

			// Helicity
			cosThetaHXVar = muPlus_HX.CosTheta();
			phiHXVar = muPlus_HX.Phi() * 180 / TMath::Pi();

			datasetHX.add(RooArgSet(centVar, massVar, yVar, ptVar, cosThetaHXVar, phiHXVar));

		} // end of reco QQ loop
	}   // enf of event loop

	TFile file(outputFileName, "RECREATE");

	datasetCS.Write();

	datasetHX.Write();

	file.Close();

	return;
}
