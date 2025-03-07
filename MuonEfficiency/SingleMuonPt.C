#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/Style/Legends.h"

// void drawPtDistribution(RooDataSet* dataset, RooWorkspace& wspace, const char* refFrameName = "CS", int ptMin = 0, int ptMax = 30, const char* extraString = "") {
TH1* drawPtDistribution(RooDataSet* dataset, RooWorkspace& wspace, const char* refFrameName = "CS", int ptMin = 0, int ptMax = 30, const char* extraString = "", double cosThetaMin = -1, double cosThetaMax = 1, int phiMin = -180, int phiMax = 180) {

    double massMin = 9., massMax = 10.;

    const char* histoName = Form("pt%dto%d_cent%dto%d_%s%s", ptMin, ptMax, gCentralityBinMin, gCentralityBinMax, refFrameName, extraString);

    // if (wspace.data("dataset")) {
    //     wspace.remove(wspace.data("dataset"));
    // }

    wspace.import(*dataset);

    RooRealVar ptVar = *wspace.var("pt");
    RooRealVar massVar = *wspace.var("mass");
    RooRealVar ptLabMuplVar = *wspace.var("ptLabMupl");
    RooRealVar ptLabMumiVar = *wspace.var("ptLabMumi");

	// reduce the whole dataset (N dimensions) to 2D (cos theta, phi)
	const char* kinematicCut = Form("(centrality >= %d && centrality < %d) && (mass > %f && mass < %f) && (rapidity > %f && rapidity < %f) && (pt > %d && pt < %d) && (%s > %f && %s < %f) && (%s > fabs(%d) && %s < fabs(%d))", 2 * gCentralityBinMin, 2 * gCentralityBinMax, massMin, massMax, gRapidityMin, gRapidityMax, ptMin, ptMax, CosThetaVarName(refFrameName), cosThetaMin, CosThetaVarName(refFrameName), cosThetaMax, PhiVarName(refFrameName), phiMin, PhiVarName(refFrameName), phiMax);

    RooDataSet* reducedDataset = (RooDataSet*)dataset->reduce(RooArgSet(ptVar, massVar, ptLabMuplVar, ptLabMumiVar), kinematicCut);

    TCanvas* canvas = new TCanvas(Form("canvas%s", refFrameName), "canvas", 600, 600);

    canvas->SetRightMargin(0.05);
    canvas->SetLeftMargin(0.18);
    canvas->SetTopMargin(0.06);
    canvas->SetBottomMargin(0.15);

    TH1* histoMupl = dynamic_cast<TH1*>(reducedDataset->createHistogram(histoName, ptLabMuplVar, RooFit::Binning(20, 0, 20)));
    TH1* histoMumi = dynamic_cast<TH1*>(reducedDataset->createHistogram(histoName, ptLabMumiVar, RooFit::Binning(20, 0, 20)));

    // Set the y-axis maximum value based on the maximum value of histoMupl or histoMumi
    double maxVal = std::max(histoMupl->GetMaximum(), histoMumi->GetMaximum());
    histoMupl->GetYaxis()->SetRangeUser(0, 1.1 * maxVal);  // Add a little extra space above the maximum value

    histoMupl->SetTitle(Form(";p_{T}^{#mu} (GeV/c); Events / (1 GeV/c)"));
    histoMupl->Draw("PE");

    histoMupl->SetMarkerStyle(21);  // Set marker style (e.g., 20 = solid circle)
    histoMupl->SetMarkerSize(1.2);  // Set marker size
    histoMupl->SetMarkerColor(TColor::GetColor("#FF1F5B"));  // Set marker color (e.g., blue)
    histoMupl->SetLineColor(TColor::GetColor("#FF1F5B"));  // Set line color (e.g., blue)

    histoMupl->GetXaxis()->CenterTitle();
    histoMupl->GetXaxis()->SetNdivisions(-20);
    histoMupl->GetXaxis()->SetTitleSize(0.05);

    histoMupl->GetYaxis()->SetTitleSize(0.05);

    histoMumi->Draw("PE same");

    histoMumi->SetMarkerStyle(22);  // Set marker style (e.g., 20 = solid circle)
    histoMumi->SetMarkerSize(1.2);  // Set marker scD6Cize
    histoMumi->SetMarkerColor(TColor::GetColor("#00CD6C"));  // Set marker color (e.g., blue)
    histoMumi->SetLineColor(TColor::GetColor("#00CD6C"));  // Set line color (e.g., blue)

    gPad->Update();

    TPaveStats* ptstats = (TPaveStats*)histoMupl->FindObject("stats");
    if (ptstats) {
        ptstats->SetX1NDC(0.73);  // Set new x1 position
        ptstats->SetY1NDC(0.76);  // Set new y1 position
        ptstats->SetX2NDC(0.92);  // Set new x2 position
        ptstats->SetY2NDC(0.92);  // Set new y2 position
        ptstats->SetBorderSize(1);
        ptstats->SetFillColor(0);
        ptstats->SetTextAlign(12);
        ptstats->SetTextFont(42);
        ptstats->SetOptStat(1111);
        ptstats->SetOptFit(0);
    }
    else {
        cout << "ptstats not found" << endl;
    }
	// frame->addObject(KinematicsText(gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, 0.67, 0.70, 0.94, 0.95));  // (HX, 2to6, -0.42to-0.14, 60to120): 0.67, 0.70, 0.94, 0.95 // (HX, 12to20, 0.42to0.70, 60to120): 0.63, 0.695, 0.94, 0.95

	// if (!isPhiFolded)
	// 	frame->addObject(RefFrameText(isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax));
	// else  
	// 	frame->addObject(RefFrameTextPhiFolded(isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax, 0.595, 0.47, 0.945, 0.66, 32)); // (HX, 2to6, -0.42to-0.14, 60to120): 0.595, 0.47, 0.945, 0.66, 32 // (HX, 12to20, 0.42to0.70, 60to120): 0.605, 0.47, 0.945, 0.66

    TPaveText* kinematicsText = KinematicsText(gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, 0.66, 0.53, 0.93, 0.73);  // (HX, 2to6, -0.42to-0.14,
    kinematicsText->Draw();

    TPaveText* refFrameText;

    if (strcmp(refFrameName, "CS") == 0)
    refFrameText = RefFrameTextPhiFolded(true, cosThetaMin, cosThetaMax, phiMin, phiMax, 0.58, 0.27, 0.935, 0.41, 32);
    else
    refFrameText = RefFrameTextPhiFolded(false, cosThetaMin, cosThetaMax, phiMin, phiMax, 0.58, 0.27, 0.935, 0.41, 32);
    
    refFrameText->Draw();

	TLatex latex;
	latex.SetTextAlign(31);
	latex.SetTextSize(0.035);
	// latex.DrawLatexNDC(.92, .70, Form("centrality %d-%d%%", gCentralityBinMin, gCentralityBinMax));
	// latex.DrawLatexNDC(.92, .64, Form("%d < p_{T}^{#mu#mu} < %d GeV/c", ptMin, ptMax));
    latex.DrawLatexNDC(.93, .49, Form("%.0f < m^{#mu#mu} < %.0f GeV/#it{c^{2}}", massMin, massMax));
    // latex.DrawLatexNDC(.92, .52, Form("%.2f < cos#theta_{%s} < %.2f", cosThetaMin, refFrameName, cosThetaMax));
    // latex.DrawLatexNDC(.92, .46, Form("%d < #varphi_{%s} < %d #circ", phiMin, refFrameName, phiMax));

	gPad->Update();

    TLegend* legend = new TLegend(0.61, 0.80, 0.76, 0.90);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(histoMupl, "#mu^{+}", "lep");
    legend->AddEntry(histoMumi, "#mu^{-}", "lep");
    legend->Draw();

	gPad->Update();

    CMS_lumi(canvas, gCMSLumiText);
    canvas->SaveAs(Form("singleMuonPt_%s_cosTheta%.2fto%.2f_absphi%dto%d.png", histoName, cosThetaMin, cosThetaMax, phiMin, phiMax), "RECREATE");
	// canvas->SaveAs(Form("frame_distrib/%s.png", histoName), "RECREATE");

    delete canvas;

    return histoMupl;
}


void SingleMuonPt(Int_t ptMin = 0, Int_t ptMax = 30, const char* extraString = "", const char* refFrameName = "HX") {

    const char* filename = Form("../Files/UpsilonSkimmedDataset%s_pTmu.root", gMuonAccName);

    TFile* f = TFile::Open(filename, "READ");
	if (!f) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
		return;
	}
	
    cout << "File " << filename << " opened" << endl;

	writeExtraText = true; // if extra text
	extraText = "       Internal";

	using namespace RooFit;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	RooWorkspace wspace("workspace");

	// // for (int ptBin = 0; ptBin < NPtBins; ++ptBin) {
    //     for (int ptBin = 1; ptBin < NPtBins; ++ptBin) {
    //         for (int cosThetaBin = 3; cosThetaBin < NCosThetaBins; ++cosThetaBin) {
    //             for (int refFrame = 0; refFrame < 2; ++refFrame) {
    //                 accEffplots_3Dto1D_comparison(gPtBinning[ptBin], gPtBinning[ptBin + 1], refFrame == 0 ? "CS" : "HX", 1, gCosThetaBinning[cosThetaBin], gCosThetaBinning[cosThetaBin + 1], 6, -180, 180, gUpsilonState, kFALSE, kFALSE, "MuonUpsilonTriggerAcc");
    //             }
    //         }
    //     }
    // drawPtDistribution((RooDataSet*)f->Get("dataset"), wspace, "HX", ptMin, ptMax, extraString);
    for (int cosThetaBin = 0; cosThetaBin < NCosThetaBins; ++cosThetaBin) {
        for (int phiBin = 0; phiBin < NPhiBins; ++phiBin) {
            drawPtDistribution((RooDataSet*)f->Get("dataset"), wspace, refFrameName, ptMin, ptMax, extraString, gCosThetaBinning[cosThetaBin], gCosThetaBinning[cosThetaBin + 1], (int)gPhiBinning[phiBin], (int)gPhiBinning[phiBin + 1]);
        }
    }
    return;
}