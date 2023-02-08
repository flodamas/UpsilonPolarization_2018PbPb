#include "../Tools/Style/tdrStyle.C"
#include "../Tools/Style/CMS_lumi.C"

TPaveText* createCustomLegend(RooRealVar mean, RooRealVar sigma, RooRealVar alpha, RooRealVar order, RooRealVar alphaR, RooRealVar orderR) {
	auto* text = new TPaveText(0.58, 0.9, 0.93, 0.6, "NDCNB");
	text->SetFillStyle(4000);
	text->SetBorderSize(0);

	text->AddText(Form("m = %.3f #pm %.3f GeV", mean.getVal(), mean.getError()));
	text->AddText(Form("#sigma = %.2f #pm %.2f MeV", 1000 * sigma.getVal(), 1000 * sigma.getError()));
	text->AddText(Form("#alpha = %.3f #pm %.3f", alpha.getVal(), alpha.getError()));
	text->AddText(Form("n = %.3f #pm %.3f", order.getVal(), order.getError()));
	text->AddText(Form("#alpha' = %.3f #pm %.3f", alphaR.getVal(), alphaR.getError()));
	text->AddText(Form("n' = %.3f #pm %.3f", orderR.getVal(), orderR.getError()));

	text->SetAllWith(" ", "align", 12);
	return text;
}

void fitMCSignalShape_asymDoubleCBplusGauss(Int_t minPt = 0, Int_t maxPt = 30, Int_t centMin = 0, Int_t centMax = 90, Long64_t nEvents = -1) {
	const char* filename = "../../Run3_muons/files/UpsilonEmbedded_2018pbpb_Oniatree_10_3_2.root";

	TFile* file = TFile::Open(filename, "READ");
	if (!file) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
		return;
	}

	cout << "File " << filename << " opened" << endl;

	TTree* OniaTree = (TTree*)file->Get("hionia/myTree");

	writeExtraText = true; // if extra text
	extraText = "      Simulation Internal";

	//tdrStyle->SetTitleYOffset(1.2);
	Float_t binMin = 8, binMax = 11;
	Int_t nBins = 60;

	// ******** Select Upsilon mass region bits ******** //
	// 2018
	// Bit1: HLT_HIL1DoubleMuOpen_v1       (Double muon inclusive)
	// Bit13: HLT_HIL3MuONHitQ10_L2MuO_MAXdR3p5_M1to5_v1  (J/psi region)
	// Bit14: HLT_HIL3Mu2p5NHitQ10_L2Mu2_M7toinf_v1 (Upsilon + high masses)
	const Int_t NTriggers = 3;
	const Int_t Bits[NTriggers] = {1, 13, 14};
	Int_t SelectedBit = 2; //(This will be used in the loop for HLTrigger and Reco_QQ_Trig)

	/// OniaTree variables

	// gen level

	Float_t Gen_weight;

	// reco level
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

	OniaTree->SetBranchAddress("Gen_weight", &Gen_weight);

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

	Float_t Reco_QQ_phi, Reco_QQ_costheta, Reco_QQ_pt, Reco_QQ_y, Reco_QQ_px, Reco_QQ_py, Reco_QQ_pz, Reco_QQ_E, Reco_QQ_m;
	Float_t Reco_mupl_phi, Reco_mupl_costheta, Reco_mupl_pt, Reco_mupl_y, Reco_mupl_eta, Reco_mupl_px, Reco_mupl_py, Reco_mupl_pz, Reco_mupl_E;
	Float_t Reco_mumi_phi, Reco_mumi_costheta, Reco_mumi_pt, Reco_mumi_y, Reco_mumi_eta, Reco_mumi_px, Reco_mumi_py, Reco_mumi_pz, Reco_mumi_E;

	// ******** Create a Ntuple to store kinematics of Upsilon and daughter muons ******** //
	gROOT->cd();
	TString varlist = "centrality:upsM:upsRap:upsPt:upsPz:upsPhi:upsCosTheta:muplPt:muplPz:muplEta:muplPhi:muplCosTheta:mumiPt:mumiPz:mumiEta:mumiPhi:muplM:muplCosThetaPrimeHX:muplPhiPrimeHX:muplCosThetaPrimeCS:muplPhiPrimeCS";
	TNtuple* UpsMuNTuple = new TNtuple("UpsMuKinematics", "Upsilon in the lab frame ntuple", varlist);

	// ******** Set beam energy for the Collins-Soper reference frame ******** //
	double sqrt_S_NN = 5.02;                 //(Center of mass Energy per nucleon pair in TeV)
	double beam1_p = sqrt_S_NN * 1000. / 2.; //(in GeV) (Note. sqrt_S_NN = sqrt(2*E1*E2+2*p1*p2) = 2E1 when two beams have the same E)
	double beam1_E = beam1_p;
	double beam2_p = -beam1_p;
	double beam2_E = beam1_E;
	//double delta = 0; //(Angle between ZHX(Z-axis in the Helicity frame) and ZCS(Z-axis in the Collins-Soper frame))

	if (nEvents == -1) nEvents = OniaTree->GetEntries();

	TTree* massTree = new TTree("massTree", "");
	massTree->SetDirectory(0);
	Float_t invMass = 0.;
	massTree->Branch("mass", &invMass, "invMass/F");

	for (Long64_t iEvent = 0; iEvent < nEvents; iEvent++) {
		if (iEvent % 10000 == 0) {
			cout << Form("\rProcessing event %lld / %lld (%.0f%%)", iEvent, nEvents, 100. * iEvent / nEvents) << flush;
		}

		// ******** Load the values ******** //
		OniaTree->GetEntry(iEvent);

		// event selection
		if (Centrality < 2 * centMin) continue;

		if (Centrality > 2 * centMax) continue;

		if (!((HLTriggers & (ULong64_t)(1 << (Bits[SelectedBit] - 1))) == (ULong64_t)(1 << (Bits[SelectedBit] - 1)))) continue; // must fire the upsilon HLT path

		// loop over reconstructed dimuon candidates
		for (int iQQ = 0; iQQ < Reco_QQ_size; iQQ++) {
			if (!((Reco_QQ_trig[iQQ] & (ULong64_t)(1 << (Bits[SelectedBit] - 1))) == (ULong64_t)(1 << (Bits[SelectedBit] - 1)))) continue; // dimuon matching

			if (Reco_QQ_sign[iQQ] != 0) continue; // only opposite-sign muon pairs

			if (Reco_QQ_VtxProb[iQQ] < 0.01) continue; // good common vertex proba

			TLorentzVector* Reco_QQ_4mom = (TLorentzVector*)CloneArr_QQ->At(iQQ);

			if (Reco_QQ_4mom->M() < binMin || Reco_QQ_4mom->M() > binMax) continue; // speedup!

			if (fabs(Reco_QQ_4mom->Rapidity()) > 2.4) continue;

			if (Reco_QQ_4mom->Pt() < minPt) continue;

			if (Reco_QQ_4mom->Pt() > maxPt) continue;

			/// single-muon selection criteria
			int iMuPlus = Reco_QQ_mupl_idx[iQQ];
			int iMuMinus = Reco_QQ_mumi_idx[iQQ];

			// global AND tracker muons
			//if (!((Reco_mu_SelectionType[iMuPlus] & 2) && (Reco_mu_SelectionType[iMuPlus] & 8))) continue;
			//if (!((Reco_mu_SelectionType[iMuMinus] & 2) && (Reco_mu_SelectionType[iMuMinus] & 8))) continue;

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
			Reco_QQ_phi = Reco_QQ_4mom->Phi();
			Reco_QQ_costheta = Reco_QQ_4mom->CosTheta();
			Reco_QQ_pt = Reco_QQ_4mom->Pt();
			Reco_QQ_y = Reco_QQ_4mom->Rapidity();
			Reco_QQ_px = Reco_QQ_4mom->Px();
			Reco_QQ_py = Reco_QQ_4mom->Py();
			Reco_QQ_pz = Reco_QQ_4mom->Pz();
			Reco_QQ_E = Reco_QQ_4mom->Energy();
			Reco_QQ_m = Reco_QQ_4mom->M();

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

			// ******** Print out momentums of two beams in the upsilon's rest frame ******** //
			// cout << endl;
			// cout << "<<Boosted to the quarkonium rest frame>>" << endl;
			// cout << "ups: p = (" << ups4MomBoosted.Px() << ", " << ups4MomBoosted.Py()  << ", " << ups4MomBoosted.Pz() << ")" << endl;
			// cout << "beam1: p = (" << beam14MomBoosted.Px() << ", " << beam14MomBoosted.Py()  << ", " << beam14MomBoosted.Pz() << ")" << endl;
			// cout << "beam2: p = (" << beam24MomBoosted.Px() << ", " << beam24MomBoosted.Py()  << ", " << beam24MomBoosted.Pz() << ")" << endl;

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

			// ******** Print out momentums of daughter muons in the upsilon's rest frame after coordinate rotation ******** //
			// cout << endl;
			// cout << "<<Rotated the quarkonium rest frame>>" << endl;
			// cout << "beam1: p = (" << beam1PvecBoosted.Px() << ", " << beam1PvecBoosted.Py()  << ", " << beam1PvecBoosted.Pz() << ")" << endl;
			// cout << "beam2: p = (" << beam2PvecBoosted.Px() << ", " << beam2PvecBoosted.Py()  << ", " << beam2PvecBoosted.Pz() << ")" << endl;

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

			invMass = Reco_QQ_m;

			massTree->Fill();
		}

	} // end of the loop on events

	cout << endl;

	Long64_t nEntries = massTree->GetEntries();

	using namespace RooFit;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	// The RooRealVar name must be the same as the branch name!!!
	RooRealVar mass("mass", "m_{#mu^{#plus}#mu^{#minus}}", binMin, binMax, "GeV");
	RooDataSet data("data", "data", massTree, RooArgSet(mass));

	auto* canvas = new TCanvas("canvas", "", 600, 650);
	//canvas->Divide(2);
	TPad* pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1.0);
	pad1->SetBottomMargin(0.01);
	pad1->SetLogy();
	pad1->Draw();
	pad1->cd();

	RooPlot* frame = mass.frame(Title(" "), Range(binMin, binMax));
	//frame->SetYTitle(Form("Candidates / (%d MeV)", (int)1000 * (binMax - binMin) / nBins));
	data.plotOn(frame, Name("data"), Binning(nBins), DrawOption("P0Z"));

	// Asymetric double CB
	RooRealVar mean("mean", "", 9., 10.);
	RooRealVar sigmaInf("sigmaInf", "", .05, .12);
	RooRealVar alphaInf("alphaInf", "", 0.1, 10);
	RooRealVar orderInf("orderInf", "", 0.1, 10);
	RooRealVar sigmaSup("sigmaSup", "", .05, .12);
	RooRealVar alphaSup("alphaSup", "", 0.1, 10);
	RooRealVar orderSup("orderSup", "", 0.1, 20);

	RooCrystalBall asymCB("asymCB", "", mass, mean, sigmaInf, sigmaSup, alphaInf, orderInf, alphaSup, orderSup);

	// Gaussian
	//RooRealVar resFraction("resFraction", "", 0.01, 1);
	RooRealVar sigma_gauss("sigma_gauss", "", .02, .2);

	RooGaussian gauss("gauss", "", mass, mean, sigma_gauss);

	// add the two
	RooRealVar normFraction("normFraction", "", 0.01, 1);
	RooAddPdf signal("signal", "sum", RooArgList(asymCB, gauss), RooArgList(normFraction));

	auto* fitResult = signal.fitTo(data, Save(), Extended(kTRUE), PrintLevel(-1), Minos(kTRUE), NumCPU(3), Range(binMin, binMax));

	fitResult->Print("v");

	signal.plotOn(frame, LineColor(kRed));

	frame->Draw();

	frame->SetMaximum(nEntries);
	frame->SetMinimum(.8);

	//auto* text = createCustomLegend(mean, sigma, alphaInf, orderInf, alphaSup, orderSup);
	//text->Draw("SAME");

	TPaveText* pt = new TPaveText(0.18, 0.9, 0.5, 0.65, "NDCNB");
	pt->SetFillColor(4000);
	pt->SetBorderSize(0);
	pt->AddText(Form("Centrality %d-%d%%", centMin, centMax));
	pt->AddText("|#eta^{#mu}| < 2.4, p_{T}^{#mu} > 3.5 GeV");
	pt->AddText("|y^{#mu#mu}| < 2.4");
	pt->AddText(Form("%d < p_{T}^{#mu#mu} < %d GeV", minPt, maxPt));

	pt->SetAllWith("", "align", 12);
	pt->Draw("SAME");

	CMS_lumi(pad1, "Hydjet-embedded PbPb MC");

	// pull distribution
	canvas->cd();
	TPad* pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, .25);
	pad2->SetTopMargin(0.015);
	pad2->SetBottomMargin(0.34);
	pad2->SetTicks(1, 1);
	pad2->Draw();
	pad2->cd();

	RooHist* hpull = frame->pullHist();
	//hpull->SetMarkerSize(0.8);
	RooPlot* pullFrame = mass.frame(Title(" "));
	pullFrame->addPlotable(hpull, "PZ");
	//pullFrame->SetTitleSize(0);
	pullFrame->GetYaxis()->SetTitleOffset(0.4);
	pullFrame->GetYaxis()->SetTitle("pull");
	pullFrame->GetYaxis()->SetTitleSize(0.17);
	pullFrame->GetYaxis()->SetLabelSize(0.15);
	//pullFrame->GetYaxis()->SetRangeUser(-4.5, 4.5);
	//  pullFrame->GetYaxis()->SetLimits(-6,6) ;
	//	pullFrame->GetYaxis()->CenterTitle();

	pullFrame->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{#minus}} (GeV)");
	//pullFrame->GetXaxis()->SetTitleOffset(1.20);
	//pullFrame->GetXaxis()->SetLabelOffset(0.1);
	pullFrame->GetXaxis()->SetLabelSize(0.15);
	pullFrame->GetXaxis()->SetTitleSize(0.17);
	//	pullFrame->GetXaxis()->CenterTitle();
	// pullFrame->GetXaxis()->SetTitleFont(43);
	// pullFrame->GetYaxis()->SetTitleFont(43);

	pullFrame->GetYaxis()->SetTickSize(0.03);
	//pullFrame->GetYaxis()->SetNdivisions(505);
	pullFrame->GetXaxis()->SetTickSize(0.03);
	pullFrame->Draw();

	TLatex textChi2;
	textChi2.SetTextAlign(12);
	textChi2.SetTextSize(0.15);
	textChi2.DrawLatexNDC(0.7, 0.85, Form("#chi^{2} / n_{d.o.f.} = %.1f", frame->chiSquare(9)));

	//canvas->Modified();
	//canvas->Update();
	canvas->cd();
	pad1->Draw();
	pad2->Draw();

	canvas->SaveAs(Form("signal/asymDoubleCBplusGauss_cent%dto%d_pt%dto%dGeV.png", centMin, centMax, minPt, maxPt), "RECREATE");
}
