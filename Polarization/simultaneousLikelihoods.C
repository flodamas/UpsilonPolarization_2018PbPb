#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/FitShortcuts.h"
#include "../Tools/Style/Legends.h"

#include "../Tools/RooFitPDFs/InvariantMassModels.h"
#include "../Tools/Style/FitDistributions.h"

#include "../Tools/RooFitPDFs/CosThetaPolarizationPDF.h"

#include "RooStats/SPlot.h"

RooDataSet* CosThetaDataset(RooDataSet* allDataset, RooWorkspace& wspace, Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS", Int_t phiMin = -180, Int_t phiMax = 180) {
	if (allDataset == nullptr) {
		cerr << "Null RooDataSet provided to the reducer method!!" << endl;
		return nullptr;
	}

	const char* kinematicCut = Form("(centrality >= %d && centrality < %d) && (rapidity > %f && rapidity < %f) && (pt > %d && pt < %d) && (phi%s > %d && phi%s < %d)", 2 * gCentralityBinMin, 2 * gCentralityBinMax, gRapidityMin, gRapidityMax, ptMin, ptMax, refFrameName, phiMin, refFrameName, phiMax);

	RooDataSet* reducedDataset = (RooDataSet*)allDataset->reduce(RooArgSet(*(wspace.var(Form("cosTheta%s", refFrameName)))), kinematicCut);

	wspace.import(*reducedDataset, RooFit::Rename("cos theta dataset"));

	return reducedDataset;
}

void simultaneousLikelihoods(Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS", Int_t nCosThetaBins = 10, Float_t cosThetaMin = -1, Float_t cosThetaMax = 1., Int_t phiMin = -180, Int_t phiMax = 180, Int_t iState = 1) {
	using namespace RooFit;
	using namespace RooStats;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	writeExtraText = true; // if extra text
	extraText = "      Internal";

	/// 1. create the signal-sWeighted dataset and prepare the NLL for the data

	// set up the data
	const char* filename = "../Files/UpsilonSkimmedDataset.root";

	TFile* f = TFile::Open(filename, "READ");
	if (!f) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
		return;
	}

	cout << "File " << filename << " opened" << endl;

	const char* datasetName = Form("dataset%s", refFrameName);
	RooDataSet* allDataset = (RooDataSet*)f->Get(datasetName);

	// import the dataset to a workspace
	RooWorkspace wspace(Form("workspace_%s", refFrameName));
	wspace.import(*allDataset);

	auto* data = InvMassCosThetaPhiDataset(allDataset, wspace, ptMin, ptMax, refFrameName, phiMin, phiMax);

	std::cout << "\n------------------------------------------\nThe dataset before creating sWeights:\n";
	data->Print();

	// read variables in the reduced dataset in the workspace
	RooRealVar* invMass = wspace.var("mass");

	RooRealVar* cosTheta = wspace.var(Form("cosTheta%s", refFrameName));

	Long64_t nEntries = data->sumEntries();

	/// Invariant mass model

	// signal: one double-sided Crystal Ball PDF (symmetric Gaussian core) per Y resonance
	// tail parameters fixed to MC extracted values, and identical for the three resonances

	const char* signalShapeName = "SymDSCB";

	// get the tail parameters of the signal shape first in case the MC fit is needed
	RooRealVar* alphaInf = new RooRealVar("alphaInf", "", 1);
	RooRealVar* orderInf = new RooRealVar("orderInf", "", 1);
	RooRealVar* alphaSup = new RooRealVar("alphaSup", "", 1);
	RooRealVar* orderSup = new RooRealVar("orderSup", "", 1);

	RooArgSet tailParams = GetMCSignalTailParameters(alphaInf, orderInf, alphaSup, orderSup, signalShapeName, ptMin, ptMax);

	auto signalModel = NominalSignalModel(wspace, alphaInf, orderInf, alphaSup, orderSup, nEntries);

	RooAbsPdf* signalPDF_1S = wspace.pdf("signalPDF_1S");
	RooAbsPdf* signalPDF_2S = wspace.pdf("signalPDF_2S");
	RooAbsPdf* signalPDF_3S = wspace.pdf("signalPDF_3S");

	RooRealVar* yield1S = wspace.var("yield1S");
	RooRealVar* yield2S = wspace.var("yield2S");
	RooRealVar* yield3S = wspace.var("yield3S");

	// background: Chebychev polynomial
	int order = 2;
	auto bkgModel = NominalBkgModel(wspace, Form("ChebychevOrder%d", order), nEntries);

	RooAbsPdf* bkgPDF = wspace.pdf("bkgPDF");

	RooRealVar* yieldBkg = wspace.var("yieldBkg");

	RooAddPdf* invMassModel = new RooAddPdf("fitModel", "", RooArgList(*signalPDF_1S, *signalPDF_2S, *signalPDF_3S, *bkgPDF), {*yield1S, *yield2S, *yield3S, *yieldBkg});

	auto* invMassFitResult = invMassModel->fitTo(*data, Save(), Extended(kTRUE), PrintLevel(-1), NumCPU(NCPUs), Range(MassBinMin, MassBinMax), Minos(kTRUE));

	invMassFitResult->Print("v");

	wspace.import(*invMassModel, RecycleConflictNodes());

	// draw the invariant mass distribution, to check the fit
	TCanvas* massCanvas = DrawMassFitDistributions(wspace, data, invMassFitResult->floatParsFinal().getSize(), ptMin, ptMax);

	massCanvas->SaveAs(Form("InvMassFits/rawInvMassFit_ChebychevBkg_cent%dto%d_pt%dto%dGeV.png", gCentralityBinMin, gCentralityBinMax, ptMin, ptMax), "RECREATE");

	// sPlot time!

	// yields from invariant mass distribution fit as sWeights
	// see https://root.cern/doc/master/classRooStats_1_1SPlot.html#a5b30f5b1b2a3723bbebef17ffb6507b2 constructor for the arguments
	RooStats::SPlot sData{"sData", "", *data, invMassModel, RooArgList(*yield1S, *yield2S, *yield3S, *yieldBkg), RooArgSet(), true, false, "dataWithSWeights", Range(MassBinMin, MassBinMax), NumCPU(NCPUs)};

	std::cout << "\n\nThe dataset after creating sWeights:\n";
	data->Print();

	// check that our weights have the desired properties

	std::cout << "\n------------------------------------------\n\nCheck SWeights:" << std::endl;

	std::cout << std::endl
	          << "Yield of Y(1S) is\t" << yield1S->getVal() << ".  From sWeights it is "
	          << sData.GetYieldFromSWeight("yield1S") << std::endl;

	std::cout << "Yield of background is\t" << yieldBkg->getVal() << ".  From sWeights it is "
	          << sData.GetYieldFromSWeight("yieldBkg") << std::endl
	          << std::endl;

	// import this new dataset with sWeights
	//std::cout << "import new dataset with sWeights" << std::endl;
	wspace.import(*data, Rename("dataWithSWeights"));

	RooDataSet data_weight1S{data->GetName(), data->GetTitle(), data, *data->get(), nullptr, "yield1S_sw"};

	// set up the NLL

	RooRealVar lambdaTheta("lambdaTheta", "lambdaTheta", -1.2, 1.2);

	auto cosThetaPDF = CosThetaPolarizationPDF("cosThetaPDF", " ", *cosTheta, lambdaTheta);

	RooNLLVar dataNLL("dataNLL", "", cosThetaPDF, data_weight1S);
	wspace.import(dataNLL);

	/// 2. create the selected MC dataset and prepare the NLL
	const char* MCfilename = Form("../Files/Y%dSSelectedMCWeightedDataset.root", iState);

	TFile* MCfile = TFile::Open(MCfilename, "READ");
	if (!MCfile) {
		cout << "File " << MCfilename << " not found. Check the directory of the file." << endl;
		return;
	}

	cout << "File " << MCfilename << " opened" << endl;

	RooDataSet* MCdataset = (RooDataSet*)MCfile->Get("MCdataset");

	// import the dataset to the workspace
	wspace.import(*MCdataset);

	auto* selectedMCDataset = CosThetaDataset(MCdataset, wspace, ptMin, ptMax, refFrameName, phiMin, phiMax);

	RooNLLVar selectedMCNLL("selectedMCNLL", "", cosThetaPDF, *selectedMCDataset);
	wspace.import(selectedMCNLL);

	/// 3. combine the NLL functions
	//	wspace.factory(sum::nllJoint(dataNLL, selectedMCNLL));

	//	RooFormulaVar totalNLL("totalNLL", "", "@0 - @1", {dataNLL, selectedMCNLL}); // ??!!

	RooAddition totalNLL("totalNLL", "totalNLL", RooArgSet(dataNLL, selectedMCNLL));

	/// 4. polarization fit!

	cout << endl
	     << "Simultaneous fit!" << endl
	     << endl;

	RooMinimizer min(totalNLL);
	//auto polarizationFitResult = min.fit(Save(), Extended(kTRUE), PrintLevel(+1), NumCPU(NCPUs), Range(cosThetaMin, cosThetaMax), AsymptoticError(DoAsymptoticError));
	min.setMinimizerType("Minuit2");
	min.migrad();
	//min.hesse(); //

	auto* polarizationFitResult = min.save();

	//auto* polarizationFitResult = cosThetaPDF.fitTo(data_weight1S, Save(), Extended(kTRUE), PrintLevel(+1), NumCPU(NCPUs), Range(cosThetaMin, cosThetaMax), AsymptoticError(DoAsymptoticError));

	polarizationFitResult->Print("v");

	/// Draw everything!

	TCanvas* canvas = new TCanvas("canvas", "canvas", 650, 600);

	RooPlot* frame = cosTheta->frame(Title(" "), Range(cosThetaMin, cosThetaMax));

	data_weight1S.plotOn(frame, Binning(nCosThetaBins), DrawOption("P0Z"), DataError(RooAbsData::SumW2), Name("data1S"));
	cosThetaPDF.plotOn(frame, LineColor(kRed + 1), Name("polaResult"));

	frame->GetYaxis()->SetMaxDigits(3);

	frame->Draw();

	gPad->RedrawAxis();

	// cosmetics

	TLegend legend(.22, .88, .5, .65);
	legend.SetTextSize(.05);
	legend.SetHeader(Form("centrality %d-%d%%, %d < p_{T}^{#mu#mu} < %d GeV/c", gCentralityBinMin, gCentralityBinMax, ptMin, ptMax));
	legend.AddEntry(frame->findObject("data1S"), Form("#varUpsilon(%dS) signal sWeights", iState), "lep");
	legend.AddEntry(frame->findObject("polaResult"), Form("distribution fit: #lambda_{#theta} = %.2f #pm %.2f", lambdaTheta.getVal(), lambdaTheta.getError()), "l");

	legend.DrawClone();

	gPad->Update();

	//CMS_lumi(canvas, gCMSLumiText);
	canvas->SaveAs(Form("DistributionFits/1D/SimultaneousFit_Y%dSCosTheta%s_cent%dto%d_pt%dto%dGeV_phi%dto%d.png", iState, refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, phiMin, phiMax), "RECREATE");
}
