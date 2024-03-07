#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/FitShortcuts.h"
#include "../Tools/Style/Legends.h"

#include "../Tools/RooFitPDFs/InvariantMassModels.h"
#include "../Tools/Style/FitDistributions.h"

#include "../Tools/RooFitPDFs/CosThetaPolarizationPDF.h"
#include "../Tools/RooFitPDFs/cosThetaPolarFunc.h"

RooDataSet* InvMassDataset(RooDataSet* allDataset, RooWorkspace& wspace, Int_t ptMin = 0, Int_t ptMax = 30, Float_t cosThetaMin = -0.1, Float_t cosThetaMax = 0.1, const char* refFrameName = "CS", Int_t phiMin = 0, Int_t phiMax = 180) {
	if (allDataset == nullptr) {
		cerr << "Null RooDataSet provided to the reducer method!!" << endl;
		return nullptr;
	}

	const char* kinematicCut = Form("(centrality >= %d && centrality < %d) && (rapidity > %f && rapidity < %f) && (pt > %d && pt < %d) && (cosTheta%s > %f && cosTheta%s < %f) && (phi%s > %d && phi%s < %d)", 2 * gCentralityBinMin, 2 * gCentralityBinMax, gRapidityMin, gRapidityMax, ptMin, ptMax, refFrameName, cosThetaMin, refFrameName, cosThetaMax, refFrameName, phiMin, refFrameName, phiMax);

	RooDataSet* reducedDataset = (RooDataSet*)allDataset->reduce(RooArgSet(*(wspace.var("mass"))), kinematicCut);

	wspace.import(*reducedDataset, Rename(Form("dataset_cosTheta_%.1fto%.1f", cosThetaMin, cosThetaMax)));

	return reducedDataset;
}

void correctedYield_1D(Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS", Int_t nCosThetaBins = 10, Float_t cosThetaMin = -1, Float_t cosThetaMax = 1., Int_t phiMin = 0, Int_t phiMax = 180, const char* filename = "../Files/WeightedUpsilonSkimmedDataset.root") {
	writeExtraText = true; // if extra text
	extraText = "      Internal";

	/// Set up the data
	using namespace RooFit;
	using namespace RooStats;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	TFile* f = TFile::Open(filename, "READ");
	if (!f) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
		return;
	}

	cout << "File " << filename << " opened" << endl;

	const char* datasetName = Form("dataset%s", refFrameName);
	RooDataSet* allDataset = (RooDataSet*)f->Get(datasetName);

	// import the dataset to a workspace
	RooWorkspace wspace(Form("workspace_%s", datasetName));
	wspace.import(*allDataset);

	RooRealVar invMass = *wspace.var("mass");

	RooRealVar cosTheta = *wspace.var(Form("cosTheta%s", refFrameName));

	Long64_t nEntries = allDataset->sumEntries() / 1000;

	/// Invariant mass model

	// signal: one double-sided Crystal Ball PDF (symmetric Gaussian core) per Y resonance
	// tail parameters fixed to MC extracted values, and identical for the three resonances

	const char* signalShapeName = "SymDSCB";

	// background: Chebychev polynomial

	int order = 2;
	const char* bkgShapeName = Form("ChebychevOrder%d", order);

	auto* invMassModel = MassFitModel(wspace, signalShapeName, bkgShapeName, ptMin, ptMax, nEntries);

	RooRealVar* yield1S = wspace.var("yield1S");
	RooRealVar* yield2S = wspace.var("yield2S");
	RooRealVar* yield3S = wspace.var("yield3S");

	/// "Standard" procedure: extract the yields per bin

	TH1D* standardCorrectedHist = new TH1D("standardCorrectedHist", " ", nCosThetaBins, cosThetaMin, cosThetaMax);

	Float_t cosThetaStep = ((cosThetaMax - cosThetaMin) / nCosThetaBins);

	TCanvas* massCanvas = 0;
	Bool_t isCSframe = (strcmp(refFrameName, "CS") == 0) ? kTRUE : kFALSE;

	Float_t maxYield = 0;

	for (Int_t iCosTheta = 0; iCosTheta < nCosThetaBins; iCosTheta++) {
		Float_t cosThetaVal = cosThetaMin + iCosTheta * cosThetaStep;

		cout << "Invariant mass fit for cos theta = [" << cosThetaVal << ", " << cosThetaVal + cosThetaStep << "]" << endl;

		RooDataSet* reducedDataset = InvMassDataset(allDataset, wspace, ptMin, ptMax, cosThetaVal, cosThetaVal + cosThetaStep, refFrameName, phiMin, phiMax);

		auto* fitResult = invMassModel->fitTo(*reducedDataset, Save(), Extended(kTRUE), PrintLevel(-1), NumCPU(NCPUs), Range(MassBinMin, MassBinMax), AsymptoticError(DoAsymptoticError), SumW2Error(!DoAsymptoticError));

		fitResult->Print("v");

		// save the invariant mass distribution fit for further checks
		// one pad for the invariant mass data distribution with fit components, one for the pull distribution
		// TCanvas* massCanvas = new TCanvas("massCanvas", "", 600, 600);
		if (massCanvas) delete massCanvas;
		massCanvas = new TCanvas("massCanvas", "", 600, 600);
		TPad* pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1.0);
		pad1->SetBottomMargin(0.03);
		pad1->Draw();
		pad1->cd();

		// RooPlot* frame = InvariantMassRooPlot(wspace, reducedDataset, invMassModel);

		RooPlot* frame = (*wspace.var("mass")).frame(Title(" "), Range(MassBinMin, MassBinMax));
		frame->GetXaxis()->SetLabelOffset(1); // to make it disappear under the pull distribution pad
		reducedDataset->plotOn(frame, Name("data"), Binning(NMassBins), DrawOption("P0Z"), DataError(RooAbsData::SumW2));

		// auto* fitModel = wspace.pdf(invMassModel);
		// fitModel->plotOn(frame, Components(*wspace.pdf("bkgPDF")), LineColor(kGray + 2), LineStyle(kDashed));
		// fitModel->plotOn(frame, Components(*wspace.pdf("signalPDF_1S")), LineColor(kRed));
		// fitModel->plotOn(frame, Components(*wspace.pdf("signalPDF_2S")), LineColor(kRed));
		// fitModel->plotOn(frame, Components(*wspace.pdf("signalPDF_3S")), LineColor(kRed));
		// fitModel->plotOn(frame, LineColor(kBlue));

		invMassModel->plotOn(frame, Components(*wspace.pdf("bkgPDF")), LineColor(kGray + 2), LineStyle(kDashed));
		invMassModel->plotOn(frame, Components(*wspace.pdf("signalPDF_1S")), LineColor(kRed));
		invMassModel->plotOn(frame, Components(*wspace.pdf("signalPDF_2S")), LineColor(kRed));
		invMassModel->plotOn(frame, Components(*wspace.pdf("signalPDF_3S")), LineColor(kRed));
		invMassModel->plotOn(frame, LineColor(kBlue));

		frame->GetYaxis()->SetMaxDigits(3);

		frame->addObject(KinematicsText(gCentralityBinMin, gCentralityBinMax, ptMin, ptMax));

		frame->addObject(RefFrameText(isCSframe, cosThetaVal, cosThetaVal + cosThetaStep, phiMin, phiMax));

		frame->addObject(FitResultText(*wspace.var("yield1S"), ComputeSignalSignificance(wspace, 1), *wspace.var("yield2S"), ComputeSignalSignificance(wspace, 2)));

		frame->Draw();

		gPad->RedrawAxis();

		// pull distribution
		massCanvas->cd();

		TPad* pad2 = GetPadPullDistribution(frame, fitResult->floatParsFinal().getSize());

		//canvas->Modified();
		//canvas->Update();
		massCanvas->cd();
		pad1->Draw();
		pad2->Draw();

		const char* fitModelName = GetFitModelName(signalShapeName, ptMin, ptMax, isCSframe, cosThetaVal, cosThetaVal + cosThetaStep, phiMin, phiMax);

		gSystem->mkdir("InvMassFits", kTRUE);
		massCanvas->SaveAs(Form("InvMassFits/CorrectedData_ChebychevOrder%d_%s.png", order, fitModelName), "RECREATE");

		// frame->Clear();
		standardCorrectedHist->SetBinContent(iCosTheta + 1, yield1S->getVal());
		standardCorrectedHist->SetBinError(iCosTheta + 1, yield1S->getError());

		if (yield1S->getVal() > maxYield) maxYield = yield1S->getVal();
	}

	RooDataHist correctedHist("correctedHist", " ", cosTheta, standardCorrectedHist);

	/// Draw the cos theta distributions

	TCanvas* canvas = new TCanvas("canvas", "canvas", 650, 600);

	RooPlot* frame = cosTheta.frame(Title(" "), Range(cosThetaMin, cosThetaMax));

	correctedHist.plotOn(frame, DrawOption("P0Z"), MarkerColor(kAzure + 2), Name("dataPoints"));

	//frame->GetYaxis()->SetRangeUser(0, 1000);
	frame->GetYaxis()->SetMaxDigits(3);

	frame->Draw();

	gPad->RedrawAxis();

	frame->SetMaximum(2 * maxYield);

	/// Polarization fit

	cout << endl
	     << "Distribution fit for polarization paramaters extraction" << endl;

	RooRealVar lambdaTheta("lambdaTheta", "lambdaTheta", -2, 2);

	auto cosThetaPDF_1S = CosThetaPolarizationPDF("cosThetaPDF_1S", " ", cosTheta, lambdaTheta);

	enableBinIntegrator(cosThetaPDF_1S, nCosThetaBins);
	auto* polarizationFitResult = cosThetaPDF_1S.fitTo(correctedHist, Save(), Extended(kTRUE), PrintLevel(+1), NumCPU(NCPUs), Range(cosThetaMin, cosThetaMax), SumW2Error(false));

	polarizationFitResult->Print("v");

	cosThetaPDF_1S.plotOn(frame, LineColor(kRed + 1), Name("polaResult"));

	frame->Draw();

	// cosmetics

	TLegend legend(.22, .88, .5, .68);
	legend.SetTextSize(.05);
	legend.SetHeader(Form("centrality %d-%d%%, %d < p_{T}^{#mu#mu} < %d GeV/c", gCentralityBinMin, gCentralityBinMax, ptMin, ptMax));
	legend.AddEntry(frame->findObject("dataPoints"), "#varUpsilon(1S) corrected yield", "lp");
	legend.AddEntry(frame->findObject("polaResult"), Form("distribution fit: #lambda_{#theta} = %.2f #pm %.2f", lambdaTheta.getVal(), lambdaTheta.getError()), "l");

	legend.DrawClone();

	gPad->Update();

	//CMS_lumi(canvas, gCMSLumiText);

	gSystem->mkdir("DistributionFits/1D", kTRUE);
	canvas->SaveAs(Form("DistributionFits/1D/%s_cent%dto%d_pt%dto%dGeV_phi%dto%d.png", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, phiMin, phiMax), "RECREATE");
}
