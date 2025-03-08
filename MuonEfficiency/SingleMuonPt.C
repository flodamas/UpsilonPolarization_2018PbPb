#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/Style/Legends.h"

// #include "../Tools/Parameters/MuonScaleFactors.h"

#include "../Tools/Parameters/tnp_weight_lowptPbPb.h"

TLine* drawLine(double x1, double y1, double x2, double y2) {
    TLine* line = new TLine(x1, y1, x2, y2);
    line->SetLineStyle(2);
    line->SetLineColor(kBlack);

    line->Draw("SAME");
    return line;
}

TH1D* drawSF(const double* ptBinEdges, int nBins = 10, double etaMin = 0, double etaMax = 1.2, Bool_t isL3 = kTRUE) {
    double pt;

    double eta = (etaMin + etaMax) / 2;

    int filterID;

    if (isL3) filterID = 3;
    else filterID = 1;

    int idx = 0;

    gStyle->SetOptStat(0);

    TH1D* frameHist = new TH1D("frameHist", "", 60, 0, 30);

    frameHist->GetXaxis()->SetTitle("p_{T}^{#mu} [GeV/c]");
    frameHist->GetXaxis()->SetTitleSize(0.1);
    frameHist->GetXaxis()->SetLabelSize(0.1);
    frameHist->GetXaxis()->CenterTitle();

    frameHist->GetYaxis()->SetTitle("SF");
    frameHist->GetYaxis()->SetTitleOffset(0.5);
    frameHist->GetYaxis()->SetTitleSize(0.1);
    frameHist->GetYaxis()->SetLabelSize(0.1);
    frameHist->GetYaxis()->CenterTitle();

    frameHist->GetXaxis()->SetNdivisions(506);

    frameHist->GetYaxis()->SetRangeUser(0.7, 1.5);
    frameHist->GetYaxis()->SetNdivisions(404);

    TGraphAsymmErrors* sfGraph = new TGraphAsymmErrors(nBins);

    // Styling the graph
    sfGraph->SetMarkerStyle(20);
    sfGraph->SetMarkerSize(1.2);
    sfGraph->SetMarkerColor(TColor::GetColor("#045275"));
    sfGraph->SetLineColor(TColor::GetColor("#045275"));

    for (int ptBin = 0; ptBin < nBins; ptBin++) {
        double pt = (ptBinEdges[ptBin] + ptBinEdges[ptBin + 1]) / 2;  // Bin center
        double sf = tnp_weight_trg_pbpb(pt, eta, filterID, idx);      // Scale factor
    
        double statUp = tnp_weight_trg_pbpb(pt, eta, filterID, 1) - sf;
        double statDown = sf - tnp_weight_trg_pbpb(pt, eta, filterID, 2);
    
        double systUp = tnp_weight_trg_pbpb(pt, eta, filterID, -1) - sf;
        double systDown = sf - tnp_weight_trg_pbpb(pt, eta, filterID, -2);
    
        double errUp = sqrt(statUp * statUp + systUp * systUp);
        double errDown = sqrt(statDown * statDown + systDown * systDown);

        // X-axis errors (bin width) - assuming equal bin widths
        double xLow = (pt - ptBinEdges[ptBin]); 
        double xHigh = (ptBinEdges[ptBin + 1] - pt);
    
        // Set the point and errors in TGraphAsymmErrors
        sfGraph->SetPoint(ptBin, pt, sf);
        sfGraph->SetPointError(ptBin, xLow, xHigh, statUp, statDown); // Using statistical errors
    
        cout << "pt: " << pt << ", SF: " << sf << ", statUp: " << statUp << ", statDown: " << statDown << endl;
    }

    TCanvas* canvasSF = new TCanvas("canvasSF", "", 600, 250);

    canvasSF->SetRightMargin(0.05);
    canvasSF->SetLeftMargin(0.12);
    canvasSF->SetBottomMargin(0.22);

    CMS_lumi(canvasSF, gCMSLumiText);

    frameHist->Draw();
    sfGraph->Draw("P SAME");  // Draw asymmetric error bars on top

    drawLine(0, 0.8, 30, 0.8);
    drawLine(0, 1, 30, 1);    
    drawLine(0, 1.2, 30, 1.2);
    drawLine(0, 1.4, 30, 1.4);

    canvasSF->Update();

    return nullptr;
}

// void drawPtDistribution(RooDataSet* dataset, RooWorkspace& wspace, const char* refFrameName = "CS", int ptMin = 0, int ptMax = 30, const char* extraString = "") {
TH1* drawPtDistribution(const char* refFrameName = "CS", double etaMin = 0, double etaMax = 1.2, int ptMin = 0, int ptMax = 30, const char* extraString = "", double cosThetaMin = -1, double cosThetaMax = 1, int phiMin = -180, int phiMax = 180, Bool_t isL3 = kTRUE) {

    const char* filename = Form("../Files/UpsilonSkimmedDataset%s_pTmu.root", gMuonAccName);

    TFile* f = TFile::Open(filename, "READ");
	if (!f) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
		return nullptr;
	}
	
    cout << "File " << filename << " opened" << endl;

	writeExtraText = true; // if extra text
	extraText = "       Internal";

	using namespace RooFit;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

    RooDataSet* dataset = (RooDataSet*)f->Get("dataset");

	RooWorkspace wspace("workspace");

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
	const char* kinematicCutMupl = Form("(centrality >= %d && centrality < %d) && (mass > %f && mass < %f) && (rapidity > %f && rapidity < %f) && (fabs(etaLabMupl) > %f && fabs(etaLabMupl) < %f) && (pt > %d && pt < %d) && (%s > %f && %s < %f) && (fabs(%s) > %d && fabs(%s) < %d)", 2 * gCentralityBinMin, 2 * gCentralityBinMax, massMin, massMax, gRapidityMin, gRapidityMax, etaMin, etaMax, ptMin, ptMax, CosThetaVarName(refFrameName), cosThetaMin, CosThetaVarName(refFrameName), cosThetaMax, PhiVarName(refFrameName), phiMin, PhiVarName(refFrameName), phiMax);
	const char* kinematicCutMumi = Form("(centrality >= %d && centrality < %d) && (mass > %f && mass < %f) && (rapidity > %f && rapidity < %f) && (fabs(etaLabMupl) > %f && fabs(etaLabMupl) < %f) && (pt > %d && pt < %d) && (%s > %f && %s < %f) && (fabs(%s) > %d && fabs(%s) < %d)", 2 * gCentralityBinMin, 2 * gCentralityBinMax, massMin, massMax, gRapidityMin, gRapidityMax, etaMin, etaMax, ptMin, ptMax, CosThetaVarName(refFrameName), cosThetaMin, CosThetaVarName(refFrameName), cosThetaMax, PhiVarName(refFrameName), phiMin, PhiVarName(refFrameName), phiMax);

    RooDataSet* reducedDatasetMupl = (RooDataSet*)dataset->reduce(RooArgSet(ptVar, massVar, ptLabMuplVar, ptLabMumiVar), kinematicCutMupl);
    RooDataSet* reducedDatasetMumi = (RooDataSet*)dataset->reduce(RooArgSet(ptVar, massVar, ptLabMuplVar, ptLabMumiVar), kinematicCutMumi);

    const int nBins = 10;
    double ptBinEdges[nBins + 1] = {3.5, 4, 4.5, 5, 5.5, 6.5, 8, 10.5, 14, 18, 30};
    RooBinning ptBinning(nBins, ptBinEdges);

    drawSF(ptBinEdges, nBins, etaMin, etaMax, isL3);

    TCanvas* canvas = new TCanvas(Form("canvas%s", refFrameName), "canvas", 600, 600);

    canvas->SetRightMargin(0.05);
    canvas->SetLeftMargin(0.18);
    canvas->SetTopMargin(0.06);
    canvas->SetBottomMargin(0.15);

    // double framehisto_ptBinEdges[12] = {0, 3.5, 4, 4.5, 5, 5.5, 6.5, 8, 10.5, 14, 18, 30};
    TH1D* framehisto = new TH1D("framehisto", "framehisto", 60, 0, 30);

    TH1* histoMupl = dynamic_cast<TH1*>(reducedDatasetMupl->createHistogram("histoMupl", ptLabMuplVar, RooFit::Binning(ptBinning)));
    TH1* histoMumi = dynamic_cast<TH1*>(reducedDatasetMumi->createHistogram("histoMumi", ptLabMumiVar, RooFit::Binning(ptBinning)));

    framehisto->GetXaxis()->SetNdivisions(-506);
    framehisto->GetXaxis()->SetRangeUser(0, 30);

    // Set the y-axis maximum value based on the maximum value of histoMupl or histoMumi
    double maxVal = std::max(histoMupl->GetMaximum(), histoMumi->GetMaximum());
    
    framehisto->GetYaxis()->SetRangeUser(0, 1.1 * maxVal);  // Add a little extra space above the maximum value

    framehisto->SetTitle(Form(";p_{T}^{#mu} (GeV/c); Events"));
    framehisto->GetXaxis()->CenterTitle();
    framehisto->GetXaxis()->SetTitleSize(0.05);
    framehisto->GetYaxis()->SetTitleSize(0.05);
    
    framehisto->SetBins(60, 0, 30);  // Explicitly set bins from 0 to 30
    framehisto->GetXaxis()->SetLimits(0, 30);  // Force X-axis limits

    gPad->SetTicks(1, 1);
    gStyle->SetOptStat(0);

    framehisto->Draw();

    histoMupl->Draw("HIST PL same");
    histoMupl->Draw("E same");

    histoMupl->SetMarkerStyle(21);  // Set marker style: solid square
    histoMupl->SetMarkerSize(1.2);  // Set marker size
    histoMupl->SetMarkerColor(TColor::GetColor("#FF1F5B"));  // red
    histoMupl->SetLineColor(TColor::GetColor("#FF1F5B"));  // 
    histoMupl->SetLineWidth(3);

    histoMumi->Draw("HIST PL same");
    histoMumi->Draw("E same");

    histoMumi->SetMarkerStyle(22);  // Set marker style: solid triangle
    histoMumi->SetMarkerSize(1.2);  // Set marker scD6Cize
    histoMumi->SetMarkerColor(TColor::GetColor("#009ADE"));  // blue
    histoMumi->SetLineColor(TColor::GetColor("#009ADE"));  // 
    histoMumi->SetLineWidth(3);

    gPad->Update();

    // TPaveStats* ptstats = (TPaveStats*)histoMupl->FindObject("stats");
    // if (ptstats) {
    //     ptstats->SetX1NDC(0.73);  // Set new x1 position
    //     ptstats->SetY1NDC(0.76);  // Set new y1 position
    //     ptstats->SetX2NDC(0.92);  // Set new x2 position
    //     ptstats->SetY2NDC(0.92);  // Set new y2 position
    //     ptstats->SetBorderSize(1);
    //     ptstats->SetFillColor(0);
    //     ptstats->SetTextAlign(12);
    //     ptstats->SetTextFont(42);
    //     ptstats->SetOptStat(1111);
    //     ptstats->SetOptFit(0);
    // }
    // else {
    //     cout << "ptstats not found" << endl;
    // }
	// frame->addObject(KinematicsText(gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, 0.67, 0.70, 0.94, 0.95));  // (HX, 2to6, -0.42to-0.14, 60to120): 0.67, 0.70, 0.94, 0.95 // (HX, 12to20, 0.42to0.70, 60to120): 0.63, 0.695, 0.94, 0.95

	// if (!isPhiFolded)
	// 	frame->addObject(RefFrameText(isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax));
	// else  
	// 	frame->addObject(RefFrameTextPhiFolded(isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax, 0.595, 0.47, 0.945, 0.66, 32)); // (HX, 2to6, -0.42to-0.14, 60to120): 0.595, 0.47, 0.945, 0.66, 32 // (HX, 12to20, 0.42to0.70, 60to120): 0.605, 0.47, 0.945, 0.66

    TPaveText* kinematicsText = KinematicsText(gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, 0.67, 0.73, 0.93, 0.92);  // (HX, 2to6, -0.42to-0.14,
    kinematicsText->Draw();

    TPaveText* refFrameText;

    if (strcmp(refFrameName, "CS") == 0)
    refFrameText = RefFrameTextPhiFolded(true, cosThetaMin, cosThetaMax, phiMin, phiMax, 0.59, 0.40, 0.94, 0.56, 32);
    else
    refFrameText = RefFrameTextPhiFolded(false, cosThetaMin, cosThetaMax, phiMin, phiMax,0.59, 0.40, 0.94, 0.56, 32);

    refFrameText->Draw();

	TLatex latex;
	latex.SetTextAlign(31);
	latex.SetTextSize(0.035);
	// latex.DrawLatexNDC(.92, .70, Form("centrality %d-%d%%", gCentralityBinMin, gCentralityBinMax));
	// latex.DrawLatexNDC(.92, .64, Form("%d < p_{T}^{#mu#mu} < %d GeV/c", ptMin, ptMax));
    latex.DrawLatexNDC(.92, .68, Form("%.0f < m^{#mu#mu} < %.0f GeV/#it{c^{2}}", massMin, massMax));
    // latex.DrawLatexNDC(.92, .52, Form("%.2f < cos#theta_{%s} < %.2f", cosThetaMin, refFrameName, cosThetaMax));
    // latex.DrawLatexNDC(.92, .46, Form("%d < #varphi_{%s} < %d #circ", phiMin, refFrameName, phiMax));
    latex.DrawLatexNDC(.92, .63, Form("%.1f < |#eta^{#mu}| < %.1f", etaMin, etaMax));

	gPad->Update();

    TLegend* legend = new TLegend(0.54, 0.82, 0.69, 0.92);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(histoMupl, "#mu^{+}", "lep");
    legend->AddEntry(histoMumi, "#mu^{-}", "lep");
    legend->Draw();

	gPad->Update();

    CMS_lumi(canvas, gCMSLumiText);
    canvas->SaveAs(Form("singleMuonPt_%s_cosTheta%.2fto%.2f_absphi%dto%d.png", histoName, cosThetaMin, cosThetaMax, phiMin, phiMax), "RECREATE");
	// canvas->SaveAs(Form("frame_distrib/%s.png", histoName), "RECREATE");

    // TFile outputFile(Form("singleMuonPt_%s_cosTheta%.2fto%.2f_absphi%dto%d.root", histoName, cosThetaMin, cosThetaMax, phiMin, phiMax), "UPDATE");

	// histoMupl->Write();
	// histoMumi->Write();

	// outputFile.Close();

    // cout << endl
    // << "CosTheta-Phi maps saved in " << outputFile.GetName() << endl;
    
    // delete canvas;

    return histoMupl;
}


void SingleMuonPt(Int_t ptMin = 0, Int_t ptMax = 30, double etaMin = 0, double etaMax = 1.2, const char* extraString = "", const char* refFrameName = "HX", Bool_t isL3 = kTRUE) {

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
            drawPtDistribution(refFrameName, etaMin, etaMax, ptMin, ptMax, extraString, gCosThetaBinning[cosThetaBin], gCosThetaBinning[cosThetaBin + 1], (int)gPhiBinning[phiBin], (int)gPhiBinning[phiBin + 1], isL3);
        }
    }
    return;
}