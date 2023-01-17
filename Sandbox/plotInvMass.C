#include "../Tools/Style/tdrStyle.C"
#include "../Tools/Style/CMS_lumi.C"

void plotInvMass(Int_t minPt = 0, Int_t maxPt = 30,
                 const char* filename = "../Files/upsilonSkimmedTree.root") {
	TFile* f = TFile::Open(filename, "READ");
	if (!f) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
		return;
	}

	cout << "File " << filename << " opened" << endl;

	writeExtraText = true; // if extra text
	extraText = "       Internal";

	//tdrStyle->SetTitleYOffset(1.2);
	Float_t binMin = 8, binMax = 14;
	Int_t nBins = 80;

	Int_t minCent = 0, maxCent = 90;

	TTreeReader reader("UpsMuKinematics", f);
	TTreeReaderValue<Float_t> centrality(reader, "centrality");
	TTreeReaderValue<Float_t> upsM(reader, "upsM"), upsPt(reader, "upsPt");

	//	TTreeReaderValue<Float_t> muonPlusPt(reader, "muonPlusPt"), muonMinusPt(reader, "muonMinusPt"), muonPlusEta(reader, "muonPlusEta"), muonMinusEta(reader, "muonMinusEta");

	TTree* massTree = new TTree("massTree", "");
	massTree->SetDirectory(0);
	Float_t invMass = 0.;
	massTree->Branch("mass", &invMass, "invMass/F");

	// loop over the number of events in the TTree
	while (reader.Next()) {
		if (*upsPt < minPt) continue;

		if (*upsPt > maxPt) continue;

		invMass = *upsM;

		massTree->Fill();

	} // end of the loop on events

	Long64_t nEntries = massTree->GetEntries();

	using namespace RooFit;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	// The RooRealVar name must be the same as the branch name!!!
	RooRealVar mass("mass", "m_{#mu^{#plus}#mu^{#minus}}", binMin, binMax, "GeV");
	RooDataSet data("data", "data", massTree, RooArgSet(mass));

	auto* canvas = new TCanvas("canvas", "", 600, 600);

	RooPlot* frame = mass.frame(Title(" "), Range(binMin, binMax));
	//frame->SetYTitle(Form("Candidates / (%d MeV)", (int)1000 * (binMax - binMin) / nBins));
	data.plotOn(frame, Name("data"), Binning(nBins), DrawOption("P0Z"));

	frame->Draw();

	//frame->SetMaximum(nEntries / 150);
	//frame->SetMinimum(0.8);

	TPaveText* pt = new TPaveText(0.5, 0.9, 0.9, 0.7, "NDCNB");
	pt->SetFillColor(4000);
	pt->SetBorderSize(0);
	pt->AddText(Form("Centrality %d-%d%%", minCent, maxCent));
	pt->AddText("|#eta^{#mu}| < 2.4, p_{T}^{#mu} > 3.5 GeV");
	pt->AddText(Form("|y^{#mu#mu}| < 2.4, %d < p_{T}^{#mu#mu} < %d GeV", minPt, maxPt));

	pt->SetAllWith("", "align", 12);
	pt->Draw();

	CMS_lumi(canvas, "2018 PbPb data (1.6 nb^{-1}, 5.02 TeV)");

	canvas->Modified();
	canvas->Update();

	canvas->SaveAs("mass_distrib/test.png", "RECREATE");
}
