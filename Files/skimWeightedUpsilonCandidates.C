#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../ReferenceFrameTransform/Transformations.h"

// (https://twiki.cern.ch/twiki/bin/viewauth/CMS/UpsilonPolarizationInPbPb5TeV)

void skimWeightedUpsilonCandidates(const char* inputFileName = "OniaTree_miniAOD_PbPbPrompt_112X_DATA_ep.root", const char* outputFileName = "WeightedUpsilonSkimmedDataset.root") {
	/// Files

	// acceptance maps
	TFile* acceptanceFile = TFile::Open("../MonteCarlo/AcceptanceMaps/1S/AcceptanceResults.root", "READ");
	if (!acceptanceFile) {
		cout << "Acceptance file not found. Check the directory of the file." << endl;
		return;
	}

	auto* accMapCS = (TEfficiency*)acceptanceFile->Get("AccMatrixCS");
	auto* accMapHX = (TEfficiency*)acceptanceFile->Get("AccMatrixHX");

	// efficiency maps
	TFile* efficiencyFile = TFile::Open("../MonteCarlo/EfficiencyMaps/1S/EfficiencyResults.root", "READ");
	if (!efficiencyFile) {
		cout << "Efficiency file not found. Check the directory of the file." << endl;
		return;
	}

	auto* effMapCS = (TEfficiency*)efficiencyFile->Get("NominalEff_CS");
	auto* systEffCS = (TH3D*)efficiencyFile->Get("RelatSystEff_CS");
	auto* effMapHX = (TEfficiency*)efficiencyFile->Get("NominalEff_HX");
	auto* systEffHX = (TH3D*)efficiencyFile->Get("RelatSystEff_HX");

	// Oniatree
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
	RooRealVar massVar("mass", gMassVarTitle, lowMassCut, highMassCut, gMassUnit);
	RooRealVar yVar("rapidity", gDimuonRapidityVarTitle, 0, 2.4);

	RooRealVar ptVar("pt", gDimuonPtVarTitle, 0, gPtBinning[NPtBins], gPtUnit);

	// separate into two datasets (easier for entry weighting) and forget about lab variables

	char* refFrameName = "CS";
	RooRealVar cosThetaCSVar(CosThetaVarName(refFrameName), CosThetaVarTitle(refFrameName), -1, 1);
	RooRealVar phiCSVar(PhiVarName(refFrameName), PhiVarTitle(refFrameName), -180, 180, gPhiUnit);

	RooRealVar dimuonWeightCSVar("dimuonWeightCS", "1 / (acc x eff) weight in the CS frame", 0, 1000000);

	RooDataSet datasetCS("datasetCS", "skimmed weighted dataset for the CS frame", RooArgSet(centVar, massVar, yVar, ptVar, cosThetaCSVar, phiCSVar, dimuonWeightCSVar), RooFit::WeightVar("dimuonWeightCS"), RooFit::StoreAsymError(RooArgSet(dimuonWeightCSVar)));

	refFrameName = "HX";
	RooRealVar cosThetaHXVar(CosThetaVarName(refFrameName), CosThetaVarTitle(refFrameName), -1, 1);
	RooRealVar phiHXVar(PhiVarName(refFrameName), PhiVarTitle(refFrameName), -180, 180, gPhiUnit);

	RooRealVar dimuonWeightHXVar("dimuonWeightHX", "1 / (acc x eff) weight in the HX frame", 0, 1000000);

	RooDataSet datasetHX("datasetHX", "skimmed weighted dataset for the HX frame", RooArgSet(centVar, massVar, yVar, ptVar, cosThetaHXVar, phiHXVar, dimuonWeightHXVar), RooFit::WeightVar("dimuonWeightHX"), RooFit::StoreAsymError(RooArgSet(dimuonWeightHXVar)));

	// loop variables
	Long64_t totEntries = OniaTree->GetEntries();

	Double_t weightCS = 0, weightHX = 0;
	Double_t errorWeightLowCS = 0, errorWeightHighCS = 0, errorWeightLowHX = 0, errorWeightHighHX = 0;

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
			if (Reco_mupl_4mom->Pt() < gMuonPtCut) continue;

			TLorentzVector* Reco_mumi_4mom = (TLorentzVector*)CloneArr_mu->At(iMuMinus);

			if (fabs(Reco_mumi_4mom->Eta()) > 2.4) continue;
			if (Reco_mumi_4mom->Pt() < gMuonPtCut) continue;

			// get positive muon's coordinates in the studied reference frames
			TVector3 muPlus_CS = MuPlusVector_CollinsSoper(*Reco_QQ_4mom, *Reco_mupl_4mom);

			TVector3 muPlus_HX = MuPlusVector_Helicity(*Reco_QQ_4mom, *Reco_mupl_4mom);

			// get the corresponding weights
			int accBinCS = accMapCS->FindFixBin(muPlus_CS.CosTheta(), muPlus_CS.Phi() * 180 / TMath::Pi(), Reco_QQ_4mom->Pt());
			double acceptanceCS = accMapCS->GetEfficiency(accBinCS);

			int effBinCS = effMapCS->FindFixBin(muPlus_CS.CosTheta(), muPlus_CS.Phi() * 180 / TMath::Pi(), Reco_QQ_4mom->Pt());
			double efficiencyCS = effMapCS->GetEfficiency(effBinCS);

			if ((acceptanceCS == 0) || (efficiencyCS == 0)) { // IMPORTANT!
				weightCS = 0;
				errorWeightLowCS = 1;
				errorWeightHighCS = 1;
			} else {
				weightCS = 1 / (acceptanceCS * efficiencyCS);

				// propagate both scale factor uncertainties and efficiency stat errors to the weight
				errorWeightLowCS = weightCS * TMath::Hypot(systEffCS->GetBinContent(effBinCS), effMapCS->GetEfficiencyErrorUp(effBinCS) / efficiencyCS);

				errorWeightHighCS = weightCS * TMath::Hypot(systEffCS->GetBinContent(effBinCS), effMapCS->GetEfficiencyErrorLow(effBinCS) / efficiencyCS);
			}

			int accBinHX = accMapHX->FindFixBin(muPlus_HX.CosTheta(), muPlus_HX.Phi() * 180 / TMath::Pi(), Reco_QQ_4mom->Pt());
			double acceptanceHX = accMapHX->GetEfficiency(accBinHX);

			int effBinHX = effMapHX->FindFixBin(muPlus_HX.CosTheta(), muPlus_HX.Phi() * 180 / TMath::Pi(), Reco_QQ_4mom->Pt());
			double efficiencyHX = effMapHX->GetEfficiency(effBinHX);

			if ((acceptanceHX == 0) || (efficiencyHX == 0)) { // IMPORTANT!
				weightHX = 0;
				errorWeightLowHX = 1;
				errorWeightHighHX = 1;
			} else {
				weightHX = 1 / (acceptanceHX * efficiencyHX);

				// propagate both scale factor uncertainties and efficiency stat errors to the weight
				errorWeightLowHX = weightHX * TMath::Hypot(systEffHX->GetBinContent(effBinHX), effMapHX->GetEfficiencyErrorUp(effBinHX) / efficiencyHX);

				errorWeightHighHX = weightHX * TMath::Hypot(systEffHX->GetBinContent(effBinHX), effMapHX->GetEfficiencyErrorLow(effBinHX) / efficiencyHX);
			}

			// fill the datasets

			centVar = Centrality;

			massVar = Reco_QQ_4mom->M();
			yVar = fabs(Reco_QQ_4mom->Rapidity());
			ptVar = Reco_QQ_4mom->Pt();

			// Collins-Soper
			cosThetaCSVar = muPlus_CS.CosTheta();
			phiCSVar = fabs(muPlus_CS.Phi() * 180 / TMath::Pi());

			dimuonWeightCSVar = weightCS;

			dimuonWeightCSVar.setAsymError(errorWeightLowCS, errorWeightHighCS);

			datasetCS.add(RooArgSet(centVar, massVar, yVar, ptVar, cosThetaCSVar, phiCSVar, dimuonWeightCSVar), weightCS, errorWeightLowCS, errorWeightHighCS);

			// Helicity
			cosThetaHXVar = muPlus_HX.CosTheta();
			phiHXVar = fabs(muPlus_HX.Phi() * 180 / TMath::Pi());

			dimuonWeightHXVar = weightHX;
			dimuonWeightHXVar.setAsymError(errorWeightLowHX, errorWeightHighHX);

			datasetHX.add(RooArgSet(centVar, massVar, yVar, ptVar, cosThetaHXVar, phiHXVar, dimuonWeightHXVar), weightHX, errorWeightLowHX, errorWeightHighHX);

		} // end of reco QQ loop
	}   // enf of event loop

	TFile file(outputFileName, "RECREATE");

	datasetCS.Write();

	datasetHX.Write();

	file.Close();

	infile->Close();

	return;
}
