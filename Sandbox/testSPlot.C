#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/FitShortcuts.h"
#include "../Tools/Parameters/PhysicsConstants.h"
#include "../Tools/CustomRoofitPDFs/ErrorFuncTimesExp.h"

#include "../Tools/Style/Legends.h"

#include "RooStats/SPlot.h"

// reduce the whole dataset (N dimensions)
RooDataSet* InvMassCosThetaWeightedDataset(RooDataSet* allDataset, RooWorkspace* wspace, Int_t ptMin = 0, Int_t ptMax = 30) {
	if (allDataset == nullptr) {
		cerr << "Null RooDataSet provided to the reducer method!!" << endl;
		return nullptr;
	}

	const char* kinematicCut = Form("(centrality >= %d && centrality < %d) && (rapidity > %f && rapidity < %f) && (pt > %d && pt < %d)", 2 * gCentralityBinMin, 2 * gCentralityBinMax, gRapidityMin, gRapidityMax, ptMin, ptMax);

	RooDataSet* reducedDataset = (RooDataSet*)allDataset->reduce(RooArgSet(*(wspace->var("mass")), *(wspace->var("cosThetaCS"))), kinematicCut);

	wspace->import(*reducedDataset, RooFit::Rename("(inv mass, cos theta CS) dataset"));

	return reducedDataset;
}

void testSPlot(Int_t ptMin = 0, Int_t ptMax = 30, const char* filename = "../Files/WeightedUpsilonSkimmedDataset.root") {
	TFile* f = TFile::Open(filename, "READ");
	if (!f) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
		return;
	}

	cout << "File " << filename << " opened" << endl;

	writeExtraText = true; // if extra text
	extraText = "      Internal";

	Int_t nCosThetaBins = 20;
	Float_t cosThetaMin = -1, cosThetaMax = 1;

	/// Set up the data
	using namespace RooFit;
	using namespace RooStats;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	// Read skimmed dataset (contains angular distributions in CS and HX after kinematic cuts due to acceptance)
	RooDataSet* allDatasetCS = (RooDataSet*)f->Get("datasetCS");

	// import the dataset to a workspace
	RooWorkspace* wspace = new RooWorkspace("workspace");
	wspace->import(*allDatasetCS);

	auto* data = InvMassCosThetaWeightedDataset(allDatasetCS, wspace, ptMin, ptMax);

	std::cout << "\n\n------------------------------------------\nThe dataset before creating sWeights:\n";
	data->Print();

	// read variables in the reduced dataset in the workspace
	RooRealVar* invMass = wspace->var("mass");

	RooRealVar* cosThetaCS = wspace->var("cosThetaCS");

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

	// Y(1S) signal shape
	RooRealVar mean_1S("mean_1S", "mean 1S", PDGmass_1S, 9.3, 9.6);
	RooRealVar sigma_1S("sigma_1S", "", .01, .15);

	RooCrystalBall signal_1S("signal_1S", "", *invMass, mean_1S, sigma_1S, *alphaInf, *orderInf, *alphaSup, *orderSup);
	RooRealVar yield1S("yield1S", "N 1S", nEntries / 5, 0, nEntries);

	// Y(2S) signal shape, mass scaling for mean and widths
	RooConstVar massScaling_2S("massScaling_2S", "", PDGmass_2S / PDGmass_1S);

	RooFormulaVar mean_2S("mean_2S", "massScaling_2S*mean_1S", RooArgSet(massScaling_2S, mean_1S));
	RooFormulaVar sigma_2S("sigma_2S", "massScaling_2S*sigma_1S", RooArgSet(massScaling_2S, sigma_1S));

	RooCrystalBall signal_2S("signal_2S", "", *invMass, mean_2S, sigma_2S, *alphaInf, *orderInf, *alphaSup, *orderSup);
	RooRealVar yield2S("yield2S", "N 2S", nEntries / 10, 0, nEntries / 2);

	// Y(3S) signal shape, mass scaling for mean and widths
	RooConstVar massScaling_3S("massScaling_3S", "", PDGmass_3S / PDGmass_1S);

	RooFormulaVar mean_3S("mean_3S", "massScaling_3S*mean_1S", RooArgSet(massScaling_3S, mean_1S));
	RooFormulaVar sigma_3S("sigma_3S", "massScaling_3S*sigma_1S", RooArgSet(massScaling_3S, sigma_1S));

	RooCrystalBall signal_3S("signal_3S", "", *invMass, mean_3S, sigma_3S, *alphaInf, *orderInf, *alphaSup, *orderSup);
	RooRealVar yield3S("yield3S", "N 3S", nEntries / 20, 0, nEntries / 2);

	// background: error function x exponential
	RooRealVar err_mu("err_mu", "err_mu", 0, 10);
	RooRealVar err_sigma("err_sigma", "err_sigma", 0, 10);
	RooRealVar exp_lambda("exp_lambda", "m_lambda", 0, 10);

	ErrorFuncTimesExp bkgPDF("bkgPDF", "", *invMass, err_mu, err_sigma, exp_lambda);
	RooRealVar yieldBkg("yieldBkg", "N background events", 0, nEntries);

	RooAddPdf* invMassModel = new RooAddPdf("invMassModel", "", RooArgList(signal_1S, signal_2S, signal_3S, bkgPDF), RooArgList(yield1S, yield2S, yield3S, yieldBkg));

	/// SPlot time!

	// yields from invariant mass distribution fit as sWeights
	// see https://root.cern/doc/master/classRooStats_1_1SPlot.html#a5b30f5b1b2a3723bbebef17ffb6507b2 constructor for the arguments
	RooStats::SPlot sData{"sData", "", *data, invMassModel, RooArgList(yield1S, yield2S, yield3S, yieldBkg), RooArgSet(), true, false, "dataWithSWeights", Range(MassBinMin, MassBinMax), AsymptoticError(DoAsymptoticError), SumW2Error(!DoAsymptoticError), NumCPU(NCPUs)};

	std::cout << "\n\nThe dataset after creating sWeights:\n";
	data->Print();

	// check that our weights have the desired properties

	std::cout << "\n\n------------------------------------------\n\nCheck SWeights:" << std::endl;

	std::cout << std::endl
	          << "Yield of Y(1S) is\t" << yield1S.getVal() << ".  From sWeights it is "
	          << sData.GetYieldFromSWeight("yield1S") << std::endl;

	std::cout << std::endl
	          << "Yield of Y(2S) is\t" << yield2S.getVal() << ".  From sWeights it is "
	          << sData.GetYieldFromSWeight("yield2S") << std::endl;

	std::cout << "Yield of background is\t" << yieldBkg.getVal() << ".  From sWeights it is "
	          << sData.GetYieldFromSWeight("yieldBkg") << std::endl
	          << std::endl;

	// import this new dataset with sWeights
	std::cout << "import new dataset with sWeights" << std::endl;
	wspace->import(*data, Rename("dataWithSWeights"));

	/// Draw the cos theta distribution with and without sWeights

	TCanvas* canvas = new TCanvas("canvas", "canvas", 650, 600);

	RooPlot* frame = cosThetaCS->frame(Title(" "), Range(cosThetaMin, cosThetaMax));
	frame->SetXTitle("cos #theta_{CS}");
	data->plotOn(frame, Binning(nCosThetaBins), DrawOption("P0Z"), Name("data"));

	// create weighted data sets
	RooDataSet data_weight1S{data->GetName(), data->GetTitle(), data, *data->get(), nullptr, "yield1S_sw"};

	data_weight1S.plotOn(frame, Binning(nCosThetaBins), DrawOption("P0Z"), MarkerColor(kRed), DataError(RooAbsData::SumW2), Name("dataWithSWeights"));
	frame->Draw();

	TLatex text;
	text.SetTextAlign(22);
	text.SetTextSize(0.05);
	text.DrawLatexNDC(.55, .85, Form("centrality %d-%d%%, %d < p_{T}^{#mu#mu} < %d GeV/c", gCentralityBinMin, gCentralityBinMax, ptMin, ptMax));
	//text.DrawLatexNDC(.48, .8, Form("%.0f < m_{#mu#mu} < %.0f GeV/c^{2}", massMin, massMax));

	TLegend legend(.25, .8, .5, .65);
	legend.SetTextSize(.05);

	legend.AddEntry(frame->findObject("data"), "efficiency-corrected dimuon events", "lp");
	legend.AddEntry(frame->findObject("dataWithSWeights"), "with #varUpsilon(1S) yield sWeights", "lp");
	legend.DrawClone();

	gPad->Update();

	CMS_lumi(canvas, gCMSLumiText);
	canvas->SaveAs(Form("sPlotTests/cosThetaCS_cent%dto%d_pt%dto%dGeV.png", gCentralityBinMin, gCentralityBinMax, ptMin, ptMax), "RECREATE");

	/// Draw the invariant mass distribution, to check the fit
	auto* massCanvas = new TCanvas("massCanvas", "", 650, 600);

	RooPlot* massFrame = invMass->frame(Title(" "), Range(MassBinMin, MassBinMax));
	//frame->SetYTitle(Form("Candidates / (%d MeV)", (int)1000 * (binMax - binMin) / nBins));
	data->plotOn(massFrame, Binning(NMassBins), DrawOption("P0Z"));

	invMassModel->plotOn(massFrame, Components(bkgPDF), LineColor(kGray + 2), LineStyle(kDashed));
	invMassModel->plotOn(massFrame, Components(signal_1S), LineColor(kRed));
	invMassModel->plotOn(massFrame, Components(signal_2S), LineColor(kRed));
	invMassModel->plotOn(massFrame, Components(signal_3S), LineColor(kRed));
	invMassModel->plotOn(massFrame, LineColor(kBlue));

	massFrame->addObject(KinematicsText(gCentralityBinMin, gCentralityBinMax, ptMin, ptMax));

	massFrame->addObject(FitResultText(yield1S, 0, yield2S, 0));

	massFrame->Draw();
	massFrame->GetYaxis()->SetMaxDigits(3);

	gPad->RedrawAxis();

	massCanvas->SaveAs(Form("sPlotTests/invMassFit_cent%dto%d_pt%dto%dGeV.png", gCentralityBinMin, gCentralityBinMax, ptMin, ptMax), "RECREATE");
}
