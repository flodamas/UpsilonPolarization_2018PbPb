#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/FitShortcuts.h"
#include "../Tools/Style/Legends.h"

#include "../Tools/Style/FitDistributions.h"

#include "../Tools/RooFitPDFs/GeneralPolarizationPDF.h"

RooDataSet* CosThetaPhiDataset(RooDataSet* allDataset, RooWorkspace& wspace, Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS", Int_t phiMin = 0, Int_t phiMax = 180) {
	if (allDataset == nullptr) {
		cerr << "Null RooDataSet provided to the reducer method!!" << endl;
		return nullptr;
	}

	const char* kinematicCut = Form("(centrality >= %d && centrality < %d) && (rapidity > %f && rapidity < %f) && (pt > %d && pt < %d) && (phi%s > %d && phi%s < %d)", 2 * gCentralityBinMin, 2 * gCentralityBinMax, gRapidityMin, gRapidityMax, ptMin, ptMax, refFrameName, phiMin, refFrameName, phiMax);

	RooDataSet* reducedDataset = (RooDataSet*)allDataset->reduce(RooArgSet(*(wspace.var(Form("cosTheta%s", refFrameName))), *(wspace.var(Form("phi%s", refFrameName)))), kinematicCut);

	wspace.import(*reducedDataset, RooFit::Rename("(cos theta, phi) dataset"));

	return reducedDataset;
}

void selectedMC_2D(Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS", Int_t phiMin = -180, Int_t phiMax = 180, Int_t iState = 1) {
	writeExtraText = true; // if extra text
	extraText = "       Internal";

	Int_t nCosThetaBins = 40;
	Float_t cosThetaMin = -1, cosThetaMax = 1.;

	/// Set up the data
	using namespace RooFit;
	using namespace RooStats;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	const char* filename = Form("../Files/Y%dSSelectedMCWeightedDataset.root", iState);

	TFile* f = TFile::Open(filename, "READ");
	if (!f) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
		return;
	}

	cout << "File " << filename << " opened" << endl;

	RooDataSet* allDataset = (RooDataSet*)f->Get("MCdataset");

	// import the dataset to a workspace
	RooWorkspace wspace("workspace");
	wspace.import(*allDataset);

	auto* data = CosThetaPhiDataset(allDataset, wspace, ptMin, ptMax, refFrameName, phiMin, phiMax);

	RooRealVar cosTheta = *wspace.var(Form("cosTheta%s", refFrameName));
	RooRealVar phi = *wspace.var(Form("phi%s", refFrameName));

	// polarization fit

	RooRealVar lambdaTheta("lambdaTheta", "lambdaTheta", -1.0, 1.0);
	RooRealVar lambdaPhi("lambdaPhi", "lambdaPhi", -1.0, 1.0);
	RooRealVar lambdaThetaPhi("lambdaThetaPhi", "lambdaThetaPhi", -1.0, 1.0);

	auto polarizationPDF = GeneralPolarizationPDF("polarizationPDF", " ", cosTheta, phi, lambdaTheta, lambdaPhi, lambdaThetaPhi);

	auto* polarizationFitResult = polarizationPDF.fitTo(*data, Save(), Extended(kTRUE), PrintLevel(+1), NumCPU(NCPUs), Range(cosThetaMin, cosThetaMax), AsymptoticError(DoAsymptoticError));

	polarizationFitResult->Print("v");

	// draw
	//gStyle->SetPadLeftMargin(.15);
	//gStyle->SetTitleYOffset(1.2);
	gStyle->SetPadRightMargin(0.18);
	gStyle->SetTitleOffset(1.2, "z");
	gStyle->SetPalette(kRainBow);
	gStyle->SetNumberContours(256);

	TCanvas* canvas = new TCanvas("canvas", "canvas", 650, 600);

	Int_t nPhiBins = 20; // just for the vizualization

	TH1* histo = data->createHistogram("histo", cosTheta, Binning(nCosThetaBins, cosThetaMin, cosThetaMax), YVar(phi, Binning(nPhiBins, phiMin, phiMax)));

	histo->SetTitle(Form(";cos #theta_{%s}; #varphi_{%s} (#circ)", refFrameName, refFrameName));
	histo->GetXaxis()->CenterTitle();
	histo->GetYaxis()->SetRangeUser(-180, 300);
	histo->GetYaxis()->CenterTitle();
	histo->GetZaxis()->SetMaxDigits(3);
	histo->Draw("COLZ");

	TH1* histoPDF = polarizationPDF.createHistogram("histoPDF", cosTheta, Binning(nCosThetaBins), Binning(nPhiBins));
	histoPDF->SetLineColor(kPink);
	histoPDF->Draw("SAME cont3");

	// cosmetics

	TLegend legend(.18, .9, .45, .7);
	legend.SetTextSize(.05);
	legend.SetHeader(Form("centrality %d-%d%%, %d < p_{T}^{#mu#mu} < %d GeV/c", gCentralityBinMin, gCentralityBinMax, ptMin, ptMax));
	legend.AddEntry(histoPDF, Form("distribution fit: #lambda_{#theta} = %.2f #pm %.2f", lambdaTheta.getVal(), lambdaTheta.getError()), "l");

	legend.DrawClone();

	gPad->Update();

	CMS_lumi(canvas, Form("#varUpsilon(%dS) Hydjet-embedded MC", gUpsilonState));

	gSystem->mkdir("DistributionFits/2D", kTRUE);
	canvas->SaveAs(Form("DistributionFits/2D/SelectedY%dSMC_%s_cent%dto%d_pt%dto%dGeV_phi%dto%d.png", iState, refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, phiMin, phiMax), "RECREATE");
}
