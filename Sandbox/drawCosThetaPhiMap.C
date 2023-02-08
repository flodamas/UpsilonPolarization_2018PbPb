#include "../Tools/Style/tdrStyle.C"
#include "../Tools/Style/CMS_lumi.C"

void drawCosThetaPhiMap(Int_t minPt = 0, Int_t maxPt = 30, Int_t centMin = 0, Int_t centMax = 90,
                        const char* filename = "../Files/miniAOD_upsilonTree.root") {
	TFile* f = TFile::Open(filename, "READ");
	if (!f) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
		return;
	}

	cout << "File " << filename << " opened" << endl;

	writeExtraText = true; // if extra text
	extraText = "       Internal";

	TTreeReader reader("UpsMuKinematics", f);
	TTreeReaderValue<Float_t> centrality(reader, "centrality");
	TTreeReaderValue<Float_t> upsM(reader, "upsM"), upsPt(reader, "upsPt");

	TTreeReaderValue<Float_t> muplCosThetaPrimeHX(reader, "muplCosThetaPrimeHX"), muplPhiPrimeHX(reader, "muplPhiPrimeHX"), muplCosThetaPrimeCS(reader, "muplCosThetaPrimeCS"), muplPhiPrimeCS(reader, "muplPhiPrimeCS");

	// (cos theta, phi) 2D distribution maps for CS and HX frames

	Int_t nCosThetaBins = 20;
	Float_t cosThetaMin = -1, cosThetaMax = 1;

	//	Double_t cosThetaBinning[] = {-1, -0.8, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1};
	//	nCosThetaBins = sizeof(cosThetaBinning) / sizeof(Double_t) - 1;

	Int_t nPhiBins = 22;
	Float_t phiMin = -220, phiMax = 220;

	//Double_t phiBinning[] = {-TMath::Pi(), -2.5, -2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2, 2.5, 180};
	//	nPhiBins = sizeof(cosThetaBinning) / sizeof(Double_t) - 1;

	TH2F* hCS = new TH2F("hCS", ";cos #theta_{CS}; #varphi_{CS} (#circ)", nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);
	TH2F* hHX = new TH2F("hHX", ";cos #theta_{HX}; #varphi_{HX} (#circ)", nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);

	// loop over the number of events in the TTree
	while (reader.Next()) {
		if (*centrality < 2 * centMin) continue;

		if (*centrality > 2 * centMax) continue;

		if (*upsPt < minPt) continue;

		if (*upsPt > maxPt) continue;

		// select dimuon candidates within the Y(1S) mass window
		if (*upsM < 9) continue;
		if (*upsM > 10) continue;

		hCS->Fill(*muplCosThetaPrimeCS, *muplPhiPrimeCS * 180 / TMath::Pi());
		hHX->Fill(*muplCosThetaPrimeHX, *muplPhiPrimeHX * 180 / TMath::Pi());

	} // end of the loop on events

	gStyle->SetPadLeftMargin(.15);
	gStyle->SetTitleYOffset(.9);
	gStyle->SetPadRightMargin(0.15);
	gStyle->SetPalette(kRainBow);

	TLatex* legend = new TLatex();
	legend->SetTextAlign(22);
	legend->SetTextSize(0.05);

	auto* canvasCS = new TCanvas("canvasCS", "", 700, 600);
	hCS->Draw("COLZ");
	legend->DrawLatexNDC(.5, .87, Form("centrality %d-%d%%, %d < p_{T}^{#mu#mu} < %d GeV", centMin, centMax, minPt, maxPt));

	hCS->GetYaxis()->SetRangeUser(-180, 240);

	gPad->Update();

	CMS_lumi(canvasCS, "2018 PbPb data (1.6 nb^{-1}, 5.02 TeV)");

	canvasCS->SaveAs(Form("frame_distrib/CS_cent%dto%d_pt%dto%dGeV.png", centMin, centMax, minPt, maxPt), "RECREATE");

	auto* canvasHX = new TCanvas("canvasHX", "", 700, 600);
	hHX->Draw("COLZ");
	legend->DrawLatexNDC(.5, .87, Form("centrality %d-%d%%, %d < p_{T}^{#mu#mu} < %d GeV", centMin, centMax, minPt, maxPt));

	gPad->Update();

	hHX->GetYaxis()->SetRangeUser(-180, 240);

	CMS_lumi(canvasHX, "2018 PbPb data (1.6 nb^{-1}, 5.02 TeV)");

	canvasHX->SaveAs(Form("frame_distrib/HX_cent%dto%d_pt%dto%dGeV.png", centMin, centMax, minPt, maxPt), "RECREATE");
}
