#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/Datasets/RooDataSetHelpers.h"
#include "../Tools/Datasets/SPlotHelpers.h"

#include "../Tools/FitShortcuts.h"
#include "../Tools/Style/Legends.h"

#include "../Tools/RooFitPDFs/InvariantMassModels.h"
#include "../Tools/Style/FitDistributions.h"

#include "../Tools/RooFitPDFs/CosThetaPolarizationPDF.h"

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
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	writeExtraText = true; // if extra text
	extraText = "      Internal";

	/// 1. create the signal-sWeighted dataset and prepare the NLL for the data

	// set up the data
	const char* filename = Form("../Files/UpsilonSkimmedDataset%s.root", gMuonAccName);

	RooWorkspace wspace = SetUpWorkspace(filename, refFrameName);

	auto* data = InvMassCosThetaPhiDataset(allDataset, wspace, ptMin, ptMax, refFrameName, phiMin, phiMax);

	// read variables in the reduced dataset in the workspace
	RooRealVar* invMass = wspace.var("mass");

	RooRealVar* cosTheta = wspace.var(Form("cosTheta%s", refFrameName));

	Long64_t nEntries = data->sumEntries();

	/// Invariant mass model

	// signal: one double-sided Crystal Ball PDF (symmetric Gaussian core) per Y resonance

	const char* signalShapeName = "SymDSCB";

	// background
	//int order = 2;
	//const char* bkgShapeName = Form("ChebychevOrder%d", order);
	const char* bkgShapeName = "ExpTimesErr";

	auto* invMassModel = MassFitModel(wspace, signalShapeName, bkgShapeName, ptMin, ptMax, nEntries);

	auto* invMassFitResult = invMassModel->fitTo(*data, Save(), Extended(kTRUE), PrintLevel(-1), NumCPU(NCPUs), Range(MassBinMin, MassBinMax), Minos(kTRUE));

	invMassFitResult->Print("v");

	wspace.import(*invMassModel, RecycleConflictNodes());

	// draw the invariant mass distribution, to check the fit
	TCanvas* massCanvas = DrawMassFitDistributions(wspace, data, invMassFitResult->floatParsFinal().getSize(), ptMin, ptMax);

	massCanvas->SaveAs(Form("InvMassFits/rawInvMassFit_%s_cent%dto%d_pt%dto%dGeV.png", bkgShapeName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax), "RECREATE");

	// sPlot time!

	SPlot sData = CreateSPlot(wspace, data, invMassModel);

	RooDataSet data_weight1S = GetSWeightedDataset(data, "1S");

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
