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

void fitMCSignalShape_tripleCB(Int_t minPt = 0, Int_t maxPt = 30, Int_t centMin = 0, Int_t centMax = 90) {
	const char* filename = "../Files/MCUpsilonSkimmedTree.root";

	TFile* file = TFile::Open(filename, "READ");
	if (!file) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
		return;
	}

	cout << "File " << filename << " opened" << endl;

	writeExtraText = true; // if extra text
	extraText = "      Simulation Internal";

	//tdrStyle->SetTitleYOffset(1.2);
	Float_t massMin = 8.5, massMax = 10.5;
	Int_t nBins = (massMax - massMin) * 20; // to get 50 MeV / bin

	// ******** Select Upsilon mass region bits ******** //
	// 2018
	// Bit1: HLT_HIL1DoubleMuOpen_v1       (Double muon inclusive)
	// Bit13: HLT_HIL3MuONHitQ10_L2MuO_MAXdR3p5_M1to5_v1  (J/psi region)
	// Bit14: HLT_HIL3Mu2p5NHitQ10_L2Mu2_M7toinf_v1 (Upsilon + high masses)
	const Int_t NTriggers = 3;
	const Int_t Bits[NTriggers] = {1, 13, 14};
	Int_t SelectedBit = 2; //(This will be used in the loop for HLTrigger and Reco_QQ_Trig)

	TTreeReader reader("UpsMuKinematics", file);
	TTreeReaderValue<Float_t> centrality(reader, "centrality"), nColl(reader, "nColl"), genWeight(reader, "genWeight");
	TTreeReaderValue<Float_t> upsM(reader, "upsM"), upsPt(reader, "upsPt");

	using namespace RooFit;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	RooRealVar mass("mass", "m_{#mu^{#plus}#mu^{#minus}}", massMin, massMax, "GeV");
	RooRealVar weight("weight", "event weight", 0, 10000);

	RooDataSet data("data", "data", RooArgSet(mass, weight), WeightVar("weight"));
	//RooDataSet data("data", "data", RooArgSet(mass));

	Float_t eventWeight = 0;

	Long64_t iEvent = 0, nEvents = reader.GetEntries(false);
	// loop over the number of events in the TTree
	while (reader.Next()) {
		iEvent++;
		if (iEvent % 10000 == 0) {
			cout << Form("\rProcessing event %lld / %lld (%.0f%%)", iEvent, nEvents, 100. * iEvent / nEvents) << flush;
		}
		if (*centrality < 2 * centMin) continue;

		if (*centrality > 2 * centMax) continue;

		if (*upsPt < minPt) continue;

		if (*upsPt > maxPt) continue;

		if (*upsM < massMin) continue;

		if (*upsM > massMax) continue;

		mass = *upsM;

		eventWeight = *genWeight * *nColl;

		weight = eventWeight;

		data.add(RooArgSet(mass, weight), eventWeight);

		//cout << "dataset filled with weight " << eventWeight << endl;

	} // end of the loop on events

	Long64_t nEntries = data.sumEntries();

	auto* canvas = new TCanvas("canvas", "", 600, 650);
	//canvas->Divide(2);
	TPad* pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1.0);
	pad1->SetBottomMargin(0.01);
	pad1->SetLogy();
	pad1->Draw();
	pad1->cd();

	RooPlot* frame = mass.frame(Title(" "), Range(massMin, massMax));
	//frame->SetYTitle(Form("Candidates / (%d MeV)", (int)1000 * (massMax - massMin) / nBins));
	data.plotOn(frame, Name("data"), Binning(nBins), DrawOption("P0Z"));

	// first CB
	RooRealVar mean("mean", "", 9., 10.);
	RooRealVar sigma_1("sigma_1", "", .05, .5);

	RooRealVar alphaPar("alphaPar", "", 0.1, 10);
	RooRealVar orderPar("orderPar", "", 0.1, 10);

	RooCBShape CB_1("CB_1", "", mass, mean, sigma_1, alphaPar, orderPar);

	// second CB
	RooRealVar resFraction_2("resFraction_2", "", 0.01, 1);
	RooFormulaVar sigma_2("sigma_2", "resFraction_2*sigma_1", RooArgSet(resFraction_2, sigma_1));

	RooCBShape CB_2("CB_2", "", mass, mean, sigma_2, alphaPar, orderPar);
	RooRealVar normFraction_2("normFraction_2", "", 0.0001, 1);

	// third CB
	RooRealVar resFraction_3("resFraction_3", "", 0.0001, 1);
	RooFormulaVar sigma_3("sigma_3", "resFraction_3*sigma_1", RooArgSet(resFraction_3, sigma_1));

	RooCBShape CB_3("CB_3", "", mass, mean, sigma_3, alphaPar, orderPar);
	RooRealVar normFraction_3("normFraction_3", "", 0.0001, 1);

	// add the two

	RooAddPdf signal("signal", "sum of two CBs PDF", RooArgList(CB_1, CB_2, CB_3), RooArgList(normFraction_2, normFraction_3), kTRUE);

	// fit
	bool doWeightedError = false;
	auto* fitResult = signal.fitTo(data, Save(), Extended(true), PrintLevel(-1), Minos(!doWeightedError), NumCPU(3), Range(massMin, massMax), SumW2Error(doWeightedError)); // quoting RooFit: "sum-of-weights and asymptotic error correction do not work with MINOS errors"

	fitResult->Print("v");

	signal.plotOn(frame, Components("CB_1"), LineColor(kRed));
	signal.plotOn(frame, Components("CB_2"), LineColor(kOrange + 1));
	signal.plotOn(frame, Components("CB_3"), LineColor(kGreen + 1));
	signal.plotOn(frame, LineColor(kBlue));

	frame->Draw();

	frame->SetMaximum(nEntries);
	frame->SetMinimum(8);

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
	pullFrame->GetYaxis()->SetTitle("Pull");
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
	textChi2.DrawLatexNDC(0.7, 0.85, Form("#chi^{2} / n_{d.o.f.} = %.1f", frame->chiSquare(fitResult->floatParsFinal().getSize())));

	//canvas->Modified();
	//canvas->Update();
	canvas->cd();
	pad1->Draw();
	pad2->Draw();

	canvas->SaveAs(Form("Figures/tripleCB_cent%dto%d_pt%dto%dGeV.png", centMin, centMax, minPt, maxPt), "RECREATE");
}
