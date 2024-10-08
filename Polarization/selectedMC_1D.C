#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/FitShortcuts.h"
#include "../Tools/Style/Legends.h"

#include "../Tools/Style/FitDistributions.h"

#include "../Tools/RooFitPDFs/CosThetaPolarizationPDF.h"

void selectedMC_1D(Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS", Int_t phiMin = -180, Int_t phiMax = 180, Int_t iState = 1) {
	writeExtraText = true; // if extra text
	extraText = "      Internal";

	Int_t nCosThetaBins = 50;
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

	auto* data = InvMassCosThetaPhiDataset(allDataset, wspace, ptMin, ptMax, refFrameName, phiMin, phiMax);

	RooRealVar cosTheta = *wspace.var(Form("cosTheta%s", refFrameName));

	// polarization fit

	RooRealVar lambdaTheta("lambdaTheta", "lambdaTheta", -1.0, 1.0);

	auto cosThetaPDF = CosThetaPolarizationPDF("cosThetaPDF", " ", cosTheta, lambdaTheta);

	auto* polarizationFitResult = cosThetaPDF.fitTo(*data, Save(), Extended(kTRUE), PrintLevel(+1), NumCPU(NCPUs), Range(cosThetaMin, cosThetaMax), AsymptoticError(DoAsymptoticError));

	polarizationFitResult->Print("v");

	// draw

	TCanvas* canvas = new TCanvas("canvas", "canvas", 650, 600);

	RooPlot* frame = cosTheta.frame(Title(" "), Range(cosThetaMin, cosThetaMax));

	data->plotOn(frame, Name("data"), Binning(nCosThetaBins), DrawOption("P0Z"), DataError(RooAbsData::SumW2));

	cosThetaPDF.plotOn(frame, LineColor(gColor1S), Name("polaResult"));

	//frame->GetYaxis()->SetRangeUser(0, 1000);
	frame->GetYaxis()->SetMaxDigits(3);
	frame->SetMaximum(2 * data->sumEntries() / nCosThetaBins);

	frame->Draw();

	// cosmetics

	TLegend legend(.22, .88, .5, .65);
	legend.SetTextSize(.05);
	legend.SetHeader(Form("centrality %d-%d%%, %d < p_{T}^{#mu#mu} < %d GeV/c", gCentralityBinMin, gCentralityBinMax, ptMin, ptMax));
	legend.AddEntry(frame->findObject("data"), Form("selected #varUpsilon(%dS) MC candidates", iState), "ep");
	legend.AddEntry(frame->findObject("polaResult"), Form("distribution fit: #lambda_{#theta} = %.2f #pm %.2f", lambdaTheta.getVal(), lambdaTheta.getError()), "l");

	legend.DrawClone();

	gPad->Update();

	//CMS_lumi(canvas, gCMSLumiText);

	gSystem->mkdir("DistributionFits/1D", kTRUE);
	canvas->SaveAs(Form("DistributionFits/1D/SelectedY%dSMC_%s_cent%dto%d_pt%dto%dGeV_phi%dto%d.png", iState, refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, phiMin, phiMax), "RECREATE");
}
