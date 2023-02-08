#include "../Tools/Style/tdrStyle.C"
#include "../Tools/Style/CMS_lumi.C"

#include "../Tools/Constants.h"

#include "../Tools/CustomRoofitPDFs/ErrorFuncTimesExp.h"

void fitAsymDoubleCB(Int_t minPt = 0, Int_t maxPt = 30, Int_t centMin = 0, Int_t centMax = 90,
                     const char* filename = "../Files/miniAOD_upsilonTree.root") {
	TFile* f = TFile::Open(filename, "READ");
	if (!f) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
		return;
	}

	cout << "File " << filename << " opened" << endl;

	writeExtraText = true; // if extra text
	extraText = "      Internal";

	//tdrStyle->SetTitleYOffset(1.2);
	Float_t binMin = 8, binMax = 14;
	Int_t nBins = 80;

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
		if (*centrality < 2 * centMin) continue;

		if (*centrality > 2 * centMax) continue;

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
	TPad* pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1.0);
	pad1->SetBottomMargin(0.015);
	pad1->Draw();
	pad1->cd();

	RooPlot* frame = mass.frame(Title(" "), Range(binMin, binMax));
	//frame->SetYTitle(Form("Candidates / (%d MeV)", (int)1000 * (binMax - binMin) / nBins));
	data.plotOn(frame, Name("data"), Binning(nBins), DrawOption("P0Z"));

	/// fitting model

	// signal: one asymmetrical double-sided Crystal Ball PDF per resonance
	// tail parameters fixed to MC extracted values, and identical for the three resonances
	RooConstVar alphaInf("alphaInf", "", 1.49);
	RooConstVar orderInf("orderInf", "", 2.41);
	RooConstVar alphaSup("alphaSup", "", 1.16);
	RooConstVar orderSup("orderSup", "", 36.4);

	// Y(1S) signal shape
	RooRealVar mean_1S("mean_1S", "mean 1S", PDGmass_1S, 9.3, 9.6);
	RooRealVar sigmaInf_1S("sigmaInf_1S", "", .05, .15);
	RooRealVar sigmaSup_1S("sigmaSup_1S", "", .05, .15);

	RooCrystalBall signal_1S("signal_1S", "", mass, mean_1S, sigmaInf_1S, sigmaSup_1S, alphaInf, orderInf, alphaSup, orderSup);
	RooRealVar nSignal_1S("nSignal_1S", "N 1S", 10000, 0, nEntries);

	// Y(2S) signal shape, mass scaling for mean and widths
	RooConstVar massScaling_2S("massScaling_2S", "", PDGmass_2S / PDGmass_1S);

	RooFormulaVar mean_2S("mean_2S", "massScaling_2S*mean_1S", RooArgSet(massScaling_2S, mean_1S));
	RooFormulaVar sigmaInf_2S("sigmaInf_2S", "massScaling_2S*sigmaInf_1S", RooArgSet(massScaling_2S, sigmaInf_1S));
	RooFormulaVar sigmaSup_2S("sigmaSup_2S", "massScaling_2S*sigmaSup_1S", RooArgSet(massScaling_2S, sigmaSup_1S));

	RooCrystalBall signal_2S("signal_2S", "", mass, mean_2S, sigmaInf_2S, sigmaSup_2S, alphaInf, orderInf, alphaSup, orderSup);
	RooRealVar nSignal_2S("nSignal_2S", "N 2S", 1000, 0, nEntries);

	// Y(3S) signal shape, mass scaling for mean and widths
	RooConstVar massScaling_3S("massScaling_3S", "", PDGmass_3S / PDGmass_1S);

	RooFormulaVar mean_3S("mean_3S", "massScaling_3S*mean_1S", RooArgSet(massScaling_3S, mean_1S));
	RooFormulaVar sigmaInf_3S("sigmaInf_3S", "massScaling_3S*sigmaInf_1S", RooArgSet(massScaling_3S, sigmaInf_1S));
	RooFormulaVar sigmaSup_3S("sigmaSup_3S", "massScaling_3S*sigmaSup_1S", RooArgSet(massScaling_3S, sigmaSup_1S));

	RooCrystalBall signal_3S("signal_3S", "", mass, mean_3S, sigmaInf_3S, sigmaSup_3S, alphaInf, orderInf, alphaSup, orderSup);
	RooRealVar nSignal_3S("nSignal_3S", "N 3S", 500, 0, nEntries);

	// background: error function x exponential
	RooRealVar err_mu("err_mu", "err_mu", 0, 15);
	RooRealVar err_sigma("err_sigma", "err_sigma", 0, 10);
	RooRealVar exp_lambda("exp_lambda", "m_lambda", 0, 15);

	ErrorFuncTimesExp bkgPDF("bkgPDF", "", mass, err_mu, err_sigma, exp_lambda);
	RooRealVar nBkg("nBkg", "N background events", 0, nEntries);

	RooAddPdf fitModel("fitModel", "", RooArgList(signal_1S, signal_2S, signal_3S, bkgPDF), RooArgList(nSignal_1S, nSignal_2S, nSignal_3S, nBkg));

	auto* fitResult = fitModel.fitTo(data, Save(), Extended(kTRUE), PrintLevel(-1), Minos(kTRUE), NumCPU(3), Range(binMin, binMax));

	fitResult->Print("v");

	fitModel.plotOn(frame, Components(bkgPDF), LineColor(kGray + 2), LineStyle(kDashed));
	fitModel.plotOn(frame, Components(signal_1S), LineColor(kRed));
	fitModel.plotOn(frame, Components(signal_2S), LineColor(kRed));
	fitModel.plotOn(frame, Components(signal_3S), LineColor(kRed));
	fitModel.plotOn(frame, LineColor(kBlue));

	frame->Draw();
	gPad->RedrawAxis();

	//frame->SetMaximum(nEntries / 150);
	//frame->SetMinimum(0.8);

	TPaveText* pt = new TPaveText(0.5, 0.9, 0.9, 0.7, "NDCNB");
	pt->SetFillColor(4000);
	pt->SetBorderSize(0);
	pt->AddText(Form("Centrality %d-%d%%", centMin, centMax));
	pt->AddText("|#eta^{#mu}| < 2.4, p_{T}^{#mu} > 3.5 GeV");
	pt->AddText(Form("|y^{#mu#mu}| < 2.4, %d < p_{T}^{#mu#mu} < %d GeV", minPt, maxPt));

	pt->SetAllWith("", "align", 12);
	pt->Draw();

	CMS_lumi(pad1, "2018 PbPb data (1.6 nb^{-1}, 5.02 TeV)");

	// pull distribution
	canvas->cd();
	TPad* pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, .25);
	pad2->SetTopMargin(0.02);
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

	canvas->SaveAs(Form("mass_distrib/asymDoubleCB_cent%dto%d_pt%dto%dGeV.png", centMin, centMax, minPt, maxPt), "RECREATE");
}
