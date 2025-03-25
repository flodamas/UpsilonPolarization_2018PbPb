/// This macro is to draw single muon efficiency using the TnP method
/// Basically it repoduces the plots for Upsilon SFs in AN-2018-316 using the header file obtained in 2020 (though not exactly the same as the one used in the analysis note due to updates)
/// The header file was obtained from the official repository of the muon analysis:
/// https://github.com/CMS-HIN-dilepton/MuonAnalysis-TagAndProbe/blob/62e934924e9a9cc8e47db0d1657c5b57650fc852/macros/tnp_weight_lowptPbPb.h

#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/Style/Legends.h"

// #include "../Tools/Parameters/MuonScaleFactors_extractEff.h"

#include "../Tools/Parameters/tnp_weight_lowptPbPb_officialRepo_extractEff.h" // plots in the analysis note (AN-2018-316)

TLine* drawLine(double x1, double y1, double x2, double y2) {
    TLine* line = new TLine(x1, y1, x2, y2);
    line->SetLineStyle(2);
    line->SetLineColor(kBlack);

    line->Draw("SAME");
    return line;
}

double getSF(TString effType = "trg", double pt = 1, double eta = 1, int filterID = 3, int idx = 0, double* numEff = nullptr, double* denEff = nullptr) {
    if (effType == TString("trg")) return tnp_weight_trg_pbpb(pt, eta, filterID, idx, numEff, denEff);
    else if (effType == TString("muid")) return tnp_weight_muid_pbpb(pt, eta, idx, numEff, denEff);
    else if (effType == TString("trk")) return tnp_weight_trk_pbpb(eta, idx, numEff, denEff);
    else {
        cout << "Invalid efficiency type" << endl;
        return -1;
    }
}

/// Draw the single muon TnP efficiency (muID and Trigger) for a given eta range (trk efficiency is in the function below due to different binning)
TH1D* drawSFEff(double etaMin = 0, double etaMax = 1.2, TString effType = "trg", Bool_t isL3 = kTRUE, Bool_t isANGraph = kFALSE) {

    /// Define the bin edges for the pT bins
    int nBins = 0;
    double* ptBinEdges = nullptr;
  
    if (effType == TString("trg")) {
        if (isANGraph == kFALSE) {
            /// Jul 2020
            if (etaMin == 0 && etaMax == 1.2) {
                nBins = 10;
                ptBinEdges = new double[nBins + 1]{3.5, 4, 4.5, 5, 5.5, 6.5, 8, 10.5, 14, 18, 30};
            }
            else if (etaMin == 1.2 && etaMax == 1.8) {
                nBins = 10;
                ptBinEdges = new double[nBins + 1]{2.07, 3, 3.5, 4, 4.5, 5, 6, 7.5, 10, 15, 30};
            }
            else if (etaMin == 1.8 && etaMax == 2.1) {
                nBins = 11;
                ptBinEdges = new double[nBins + 1]{1.5, 2.5, 3, 3.5, 4, 4.5, 5.5, 6.5, 8, 9.5, 13, 30};
            }
            else if (etaMin == 2.1 && etaMax == 2.4) {
                nBins = 7;
                ptBinEdges = new double[nBins + 1]{1.5, 2.5, 3, 4.5, 6.5, 8.5, 11, 30};
            }
            else {
                cout << "Invalid eta range" << endl;
                return nullptr;
            }
        }

        else if (isANGraph == kTRUE) {
            /// Oct 2019 (Upsilon) (plots in the analysis note, AN-2018-316)
            if (etaMin == 0 && etaMax == 1.2) {
                nBins = 10;
                ptBinEdges = new double[nBins + 1]{3.5, 4, 4.5, 5, 5.5, 6.5, 8, 10.5, 14, 18, 30};
            }
            else if (etaMin == 1.2 && etaMax == 1.8) {
                nBins = 10;
                ptBinEdges = new double[nBins + 1]{2.37, 3, 3.5, 4, 4.5, 5, 6, 7.5, 10, 15, 30};
            }
            else if (etaMin == 1.8 && etaMax == 2.1) {
                nBins = 12;
                ptBinEdges = new double[nBins + 1]{1.8, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7.5, 10, 15, 30};
            }
            else if (etaMin == 2.1 && etaMax == 2.4) {
                nBins = 9;
                ptBinEdges = new double[nBins + 1]{1.8, 2.2, 2.7, 3.2, 3.7, 4.7, 6.5, 8.5, 11, 30};
            }
            else {
                cout << "Invalid eta range" << endl;
                return nullptr;
            }
        }
    }

    else if (effType == TString("muid")) {
        if (isANGraph == kFALSE) {
            /// Jul 2020
            if (etaMin == 0 && etaMax == 1.2) {
                nBins = 10;
                ptBinEdges = new double[nBins + 1]{3.5, 4, 4.5, 5, 5.5, 6.5, 8, 10.5, 14, 18, 30};
            }
            else if (etaMin == 1.2 && etaMax == 1.8) {
                nBins = 10;
                ptBinEdges = new double[nBins + 1]{2.07, 3, 3.5, 4, 4.5, 5, 6, 7.5, 10, 15, 30};
            }
            else if (etaMin == 1.8 && etaMax == 2.1) {
                nBins = 10;
                ptBinEdges = new double[nBins + 1]{1.5, 2.5, 3, 3.5, 4, 4.5, 5.5, 7, 9, 12, 30};
            }
            else if (etaMin == 2.1 && etaMax == 2.4) {
                nBins = 9;
                ptBinEdges = new double[nBins + 1]{1.5, 2.2, 2.7, 3.2, 3.7, 4.7, 8, 11, 14, 30};
            }
            else {
                cout << "Invalid eta range" << endl;
                return nullptr;
            }
        }

        else if (isANGraph == kTRUE) {
            /// Oct 2019 (Upsilon) (plots in the analysis note, AN-2018-316)
            if (etaMin == 0 && etaMax == 1.2) {
                nBins = 10;
                ptBinEdges = new double[nBins + 1]{3.5, 4, 4.5, 5, 5.5, 6.5, 8, 10.5, 14, 18, 30};
            }
            else if (etaMin == 1.2 && etaMax == 1.8) {
                nBins = 10;
                ptBinEdges = new double[nBins + 1]{2.37, 3, 3.5, 4, 4.5, 5, 6, 7.5, 10, 15, 30};
            }
            else if (etaMin == 1.8 && etaMax == 2.1) {
                nBins = 12;
                ptBinEdges = new double[nBins + 1]{1.8, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7.5, 10, 15, 30};
            }
            else if (etaMin == 2.1 && etaMax == 2.4) {
                nBins = 9;
                ptBinEdges = new double[nBins + 1]{1.8, 2.2, 2.7, 3.2, 3.7, 4.7, 6.5, 8.5, 11, 30};
            }
            else {
                cout << "Invalid eta range" << endl;
                return nullptr;
            }
        }
    }

    double eta = (etaMin + etaMax) / 2;

    int filterID;

    if (isL3) filterID = 3; // Upsilon L3
    else filterID = 2; // Upsilon L2

    int idx = 0;

    /// set frame histogram for Efficiency plot
    gStyle->SetOptStat(0);
    gStyle->SetEndErrorSize(0);  // Remove cross-line (caps) at the end of error bars

    TH1D* frameHistEff = new TH1D("frameHistEff", "", 60, 0, 30);

    frameHistEff->GetXaxis()->SetTitleSize(0.0);
    frameHistEff->GetXaxis()->SetLabelSize(0.0);

    frameHistEff->GetYaxis()->SetTitle("Single #mu Efficiency");
    frameHistEff->GetYaxis()->SetTitleOffset(1.2);
    frameHistEff->GetYaxis()->SetTitleSize(0.05);
    frameHistEff->GetYaxis()->SetLabelSize(0.05);
    // frameHistEff->GetYaxis()->CenterTitle();

    frameHistEff->GetYaxis()->SetRangeUser(0, 1.05);
    frameHistEff->GetYaxis()->SetNdivisions(409);

    /// set frame histogram for SF plot
    TH1D* frameHistSF = new TH1D("frameHistSF", "", 60, 0, 30);

    frameHistSF->GetXaxis()->SetTitle("p_{T}^{#mu} [GeV/c]");
    frameHistSF->GetXaxis()->SetTitleSize(0.15);
    frameHistSF->GetXaxis()->SetLabelSize(0.15);
    frameHistSF->GetXaxis()->CenterTitle();

    frameHistSF->GetYaxis()->SetTitle("SF");
    frameHistSF->GetYaxis()->SetTitleOffset(0.4);
    frameHistSF->GetYaxis()->SetTitleSize(0.15);
    frameHistSF->GetYaxis()->SetLabelSize(0.15);
    frameHistSF->GetYaxis()->CenterTitle();

    frameHistSF->GetXaxis()->SetNdivisions(506);

    frameHistSF->GetYaxis()->SetRangeUser(0.7, 1.5);
    frameHistSF->GetYaxis()->SetNdivisions(404);

    /// set points and errors for the Efficiency and SF plots
    TGraphAsymmErrors* sfNumEffGraph = new TGraphAsymmErrors(nBins);
    TGraphAsymmErrors* sfDenEffGraph = new TGraphAsymmErrors(nBins);
    TGraphAsymmErrors* sfGraph = new TGraphAsymmErrors(nBins);
    
    TGraphAsymmErrors* systEffGraph = new TGraphAsymmErrors(nBins);
    TGraphAsymmErrors* systRatioGraph = new TGraphAsymmErrors(nBins);

    /// Styling the graph
    sfNumEffGraph->SetMarkerStyle(25);
    sfNumEffGraph->SetMarkerSize(1.2);
    sfNumEffGraph->SetMarkerColor(4);
    sfNumEffGraph->SetLineColor(4);

    sfDenEffGraph->SetMarkerStyle(20);
    sfDenEffGraph->SetMarkerSize(1.2);
    sfDenEffGraph->SetMarkerColor(2);
    sfDenEffGraph->SetLineColor(2);

    sfGraph->SetMarkerStyle(20);
    sfGraph->SetMarkerSize(1.2);
    sfGraph->SetMarkerColor(TColor::GetColor("#333333"));
    sfGraph->SetLineColor(TColor::GetColor("#333333"));

    systEffGraph->SetMarkerStyle(22);
    systEffGraph->SetMarkerSize(1.2);
    systEffGraph->SetMarkerColor(TColor::GetColor("#46AEA0"));
    systEffGraph->SetLineColor(TColor::GetColor("#46AEA0"));

    systRatioGraph->SetMarkerStyle(20);
    systRatioGraph->SetMarkerSize(1.2);
    systRatioGraph->SetMarkerColor(TColor::GetColor("#333333"));
    systRatioGraph->SetLineColor(TColor::GetColor("#333333"));

    /// Read the efficiencies and scale factors from the header file
    for (int ptBin = 0; ptBin < nBins; ptBin++) {
        double numEff = 1;
        double denEff = 1;

        double pt = (ptBinEdges[ptBin] + ptBinEdges[ptBin + 1]) / 2;  // Bin center
        
        double sf = getSF(effType, pt, eta, filterID, idx, &numEff, &denEff); // nominal scale factor

        double statUp = getSF(effType, pt, eta, filterID, 1, &numEff, &denEff) - sf;
        double statDown = sf - getSF(effType, pt, eta, filterID, 2, &numEff, &denEff);
    
        double systUp = getSF(effType, pt, eta, filterID, -1, &numEff, &denEff) * denEff;
        double systDown = getSF(effType, pt, eta, filterID, -2, &numEff, &denEff) * denEff;
    
        double errUp = sqrt(statUp * statUp + systUp * systUp);
        double errDown = sqrt(statDown * statDown + systDown * systDown);

        /// X-axis errors (bin width) - assuming equal bin widths
        double xLow = (pt - ptBinEdges[ptBin]); 
        double xHigh = (ptBinEdges[ptBin + 1] - pt);
    
        /// Set the point and errors in TGraphAsymmErrors
        sfNumEffGraph->SetPoint(ptBin, pt, numEff); // data TnP efficiency
        sfNumEffGraph->SetPointError(ptBin, xLow, xHigh, statDown * numEff, statUp * numEff); // No errors for data TnP efficiency

        sfDenEffGraph->SetPoint(ptBin, pt, denEff); // MC TnP efficiency
        sfDenEffGraph->SetPointError(ptBin, xLow, xHigh, 0, 0); // No errors for MC TnP efficiency

        sfGraph->SetPoint(ptBin, pt, sf);
        sfGraph->SetPointError(ptBin, xLow, xHigh, statDown, statUp); // Using statistical errors

        systEffGraph->SetPoint(ptBin, pt, systUp);
        systEffGraph->SetPointError(ptBin, xLow, xHigh, 0, 0); // No errors info for systematic uncertainties

        systRatioGraph->SetPoint(ptBin, pt, systUp / numEff);
        systRatioGraph->SetPointError(ptBin, xLow, xHigh, 0, 0); // No errors info for systematic uncertainties
    
        // cout << "pt: " << pt << ", SF: " << sf << ", statUp: " << statUp << ", statDown: " << statDown << endl;
    }

    TCanvas* canvasSFEff = new TCanvas("canvasSFEff", "", 600, 600);

    TPad* padEff = new TPad("padEff", "pad for efficiency", 0.0, 0.25, 1.0, 1.0);

    padEff->Draw();
    padEff->cd();

    padEff->SetRightMargin(0.05);
    padEff->SetLeftMargin(0.12);
    padEff->SetTopMargin(0.07);
    padEff->SetBottomMargin(0.02);

    frameHistEff->Draw();
    sfNumEffGraph->Draw("P SAME"); 
    sfDenEffGraph->Draw("P SAME");

	TLegend* legendEff = new TLegend(0.5, 0.2, 0.7, 0.5);
    legendEff->SetTextSize(0.045);
    legendEff->SetBorderSize(0);
	legendEff->SetFillStyle(0);
    if (effType == TString("trg")) legendEff->SetHeader(Form("#splitline{%s Upsilon Trigger Efficiency}{(p_{T}^{#mu} > %.1f GeV, |#eta| #in [%.1f, %.1f])} ", isL3 ? "L3" : "L2", ptBinEdges[0], etaMin, etaMax));
    else if (effType == TString("muid")) legendEff->SetHeader(Form("#splitline{Hybrid Soft ID 2018 Efficiency}{(p_{T}^{#mu} > %.1f GeV, |#eta| #in [%.1f, %.1f])} ", ptBinEdges[0], etaMin, etaMax));
    else if (effType == TString("trk")) legendEff->SetHeader(Form("#splitline{Inner tracking Efficiency}{(p_{T}^{#mu} > 0.0 GeV, |#eta| #in [%.1f, %.1f])} ", etaMin, etaMax));
    legendEff->AddEntry((TObject*)0, "", "");
    legendEff->AddEntry(sfDenEffGraph, "MC PYTHIA+HYDJET", "lep");
    legendEff->AddEntry(sfNumEffGraph, "DATA", "lep");
	legendEff->Draw("SAME");

    padEff->Update();
    // gPad->Update();
    canvasSFEff->Update();

    canvasSFEff->cd();

    TPad* padSF = new TPad("padSF", "pad for SF", 0.0, 0.0, 1.0, 0.25);

    padSF->Draw();
    padSF->cd();

    padSF->SetTopMargin(0.01);
    padSF->SetRightMargin(0.05);
    padSF->SetLeftMargin(0.12);
    padSF->SetBottomMargin(0.40);

    frameHistSF->Draw();
    sfGraph->Draw("P SAME");  // Draw asymmetric error bars on top

    drawLine(0, 0.8, 30, 0.8);
    drawLine(0, 1, 30, 1);    
    drawLine(0, 1.2, 30, 1.2);
    drawLine(0, 1.4, 30, 1.4);

    CMS_lumi(padEff, gCMSLumiText);
    
    gPad->Update();
    canvasSFEff->Update();

    gSystem->mkdir("SingleMuonEfficiency/", kTRUE);
    
    if (isANGraph == kFALSE) canvasSFEff->SaveAs(Form("SingleMuonEfficiency/SingleMuonEfficiency_SF_%s%s_%.1f_%.1f.png", effType.Data(), effType == TString("trg") ? (isL3 ? "_UpsL3" : "_UpsL2") : "", etaMin, etaMax));
    else canvasSFEff->SaveAs(Form("SingleMuonEfficiency/SingleMuonEfficiency_SF_%s%s_%.1f_%.1f_old.png", effType.Data(), effType == TString("trg") ? (isL3 ? "_UpsL3" : "_UpsL2") : "", etaMin, etaMax));

    /// draw the efficiency with the systematic variation
    TCanvas* canvasSystEff = new TCanvas("canvasSystEff", "", 600, 600);

    TPad* padSystEff = new TPad("padSystEff", "pad for systematic efficiency", 0.0, 0.25, 1.0, 1.0);

    padSystEff->Draw();
    padSystEff->cd();

    padSystEff->SetRightMargin(0.05);
    padSystEff->SetLeftMargin(0.12);
    padSystEff->SetTopMargin(0.07);
    padSystEff->SetBottomMargin(0.02);

    frameHistEff->Draw();
    sfNumEffGraph->Draw("P SAME"); 
    systEffGraph->Draw("P SAME");

	TLegend* legendSystEff = new TLegend(0.5, 0.2, 0.7, 0.5);
    legendSystEff->SetTextSize(0.045);
    legendSystEff->SetBorderSize(0);
	legendSystEff->SetFillStyle(0);
    if (effType == TString("trg")) legendSystEff->SetHeader(Form("#splitline{%s Upsilon Trigger Efficiency}{(p_{T}^{#mu} > %.1f GeV, |#eta| #in [%.1f, %.1f])} ", isL3 ? "L3" : "L2", ptBinEdges[0], etaMin, etaMax));
    else if (effType == TString("muid")) legendSystEff->SetHeader(Form("#splitline{Hybrid Soft ID 2018 Efficiency}{(p_{T}^{#mu} > %.1f GeV, |#eta| #in [%.1f, %.1f])} ", ptBinEdges[0], etaMin, etaMax));
    else if (effType == TString("trk")) legendSystEff->SetHeader(Form("#splitline{Inner tracking Efficiency}{(p_{T}^{#mu} > 0.0 GeV, |#eta| #in [%.1f, %.1f])} ", etaMin, etaMax));
    legendSystEff->AddEntry((TObject*)0, "", "");
    legendSystEff->AddEntry(sfNumEffGraph, "nominal", "lep");
    legendSystEff->AddEntry(systEffGraph, "max syst", "lep");
	legendSystEff->Draw("SAME");

    padSystEff->Update();
    // gPad->Update();
    canvasSystEff->Update();

    canvasSystEff->cd();

    TPad* padRatio = new TPad("padRatio", "pad for Ratio", 0.0, 0.0, 1.0, 0.25);

    padRatio->Draw();
    padRatio->cd();

    padRatio->SetTopMargin(0.01);
    padRatio->SetRightMargin(0.05);
    padRatio->SetLeftMargin(0.12);
    padRatio->SetBottomMargin(0.40);

    frameHistSF->GetYaxis()->SetTitle("Ratio");
    frameHistSF->Draw();
    systRatioGraph->Draw("P SAME");  // Draw asymmetric error bars on top

    drawLine(0, 0.8, 30, 0.8);
    drawLine(0, 1, 30, 1);    
    drawLine(0, 1.2, 30, 1.2);
    drawLine(0, 1.4, 30, 1.4);

    CMS_lumi(padSystEff, gCMSLumiText);
    
    gPad->Update();
    canvasSystEff->Update();

    if (isANGraph == kFALSE) canvasSystEff->SaveAs(Form("SingleMuonEfficiency/SingleMuonEfficiency_syst_%s%s_%.1f_%.1f.png", effType.Data(), effType == TString("trg") ? (isL3 ? "_UpsL3" : "_UpsL2") : "", etaMin, etaMax));
    else canvasSystEff->SaveAs(Form("SingleMuonEfficiency/SingleMuonEfficiency_syst_%s%s_%.1f_%.1f_old.png", effType.Data(), effType == TString("trg") ? (isL3 ? "_UpsL3" : "_UpsL2") : "", etaMin, etaMax));

    return nullptr;
}

/// Draw trcking TnP efficiency and SF (It's separate from the trigger and muon ID efficiencies due to the different binning)
void drawSFEff_trk(Bool_t isANGraph = kFALSE) {

    const int nEtaBins_trk = 11;
    double etaBinning_trk[nEtaBins_trk + 1] = {-2.4, -1.6, -1.2, -0.9, -0.6, -0.3, 0.3, 0.6, 0.9, 1.2, 1.6, 2.4};

    TString effType = "trk";
    
    int idx = 0;

    int dummy = 1;

    /// set frame histogram for Efficiency plot
    gStyle->SetOptStat(0);
    gStyle->SetEndErrorSize(0);  // Remove cross-line (caps) at the end of error bars
    
    TH1D* frameHistEff = new TH1D("frameHistEff", "", nEtaBins_trk, etaBinning_trk[0], etaBinning_trk[nEtaBins_trk]);

    frameHistEff->GetXaxis()->SetTitleSize(0.0);
    frameHistEff->GetXaxis()->SetLabelSize(0.0);

    frameHistEff->GetXaxis()->SetNdivisions(510);

    frameHistEff->GetYaxis()->SetTitle("Single #mu Efficiency");
    frameHistEff->GetYaxis()->SetTitleOffset(1.2);
    frameHistEff->GetYaxis()->SetTitleSize(0.05);
    frameHistEff->GetYaxis()->SetLabelSize(0.05);
    // frameHistEff->GetYaxis()->CenterTitle();

    frameHistEff->GetYaxis()->SetRangeUser(0.6, 1.05);
    frameHistEff->GetYaxis()->SetNdivisions(509);

    /// set frame histogram for SF plot
    TH1D* frameHistSF = new TH1D("frameHistSF", "", nEtaBins_trk, etaBinning_trk[0], etaBinning_trk[nEtaBins_trk]);

    frameHistSF->GetXaxis()->SetTitle("#eta^{#mu}");
    frameHistSF->GetXaxis()->SetTitleSize(0.15);
    frameHistSF->GetXaxis()->SetLabelSize(0.15);
    frameHistSF->GetXaxis()->CenterTitle();

    frameHistSF->GetYaxis()->SetTitle("SF");
    frameHistSF->GetYaxis()->SetTitleOffset(0.4);
    frameHistSF->GetYaxis()->SetTitleSize(0.15);
    frameHistSF->GetYaxis()->SetLabelSize(0.15);
    frameHistSF->GetYaxis()->CenterTitle();

    frameHistSF->GetXaxis()->SetNdivisions(510);

    frameHistSF->GetYaxis()->SetRangeUser(0.96, 1.019);
    frameHistSF->GetYaxis()->SetNdivisions(503);

    /// set points and errors for the Efficiency and SF plots
    TGraphAsymmErrors* sfNumEffGraph = new TGraphAsymmErrors(nEtaBins_trk);
    TGraphAsymmErrors* sfDenEffGraph = new TGraphAsymmErrors(nEtaBins_trk);
    TGraphAsymmErrors* sfGraph = new TGraphAsymmErrors(nEtaBins_trk);

    TGraphAsymmErrors* systEffGraph = new TGraphAsymmErrors(nEtaBins_trk);
    TGraphAsymmErrors* systRatioGraph = new TGraphAsymmErrors(nEtaBins_trk);

    /// Styling the graph
    sfNumEffGraph->SetMarkerStyle(25);
    sfNumEffGraph->SetMarkerSize(1.2);
    sfNumEffGraph->SetMarkerColor(4);
    sfNumEffGraph->SetLineColor(4);

    sfDenEffGraph->SetMarkerStyle(20);
    sfDenEffGraph->SetMarkerSize(1.2);
    sfDenEffGraph->SetMarkerColor(2);
    sfDenEffGraph->SetLineColor(2);

    sfGraph->SetMarkerStyle(20);
    sfGraph->SetMarkerSize(1.2);
    sfGraph->SetMarkerColor(TColor::GetColor("#333333"));
    sfGraph->SetLineColor(TColor::GetColor("#333333"));

    systEffGraph->SetMarkerStyle(22);
    systEffGraph->SetMarkerSize(1.2);
    systEffGraph->SetMarkerColor(TColor::GetColor("#46AEA0"));
    systEffGraph->SetLineColor(TColor::GetColor("#46AEA0"));

    systRatioGraph->SetMarkerStyle(20);
    systRatioGraph->SetMarkerSize(1.2);
    systRatioGraph->SetMarkerColor(TColor::GetColor("#333333"));
    systRatioGraph->SetLineColor(TColor::GetColor("#333333"));
    
    /// Read the efficiencies and scale factors from the header file
    for (int etaBin = 0; etaBin < nEtaBins_trk; etaBin++) {

        double numEff = 1;
        double denEff = 1;

        double eta = (etaBinning_trk[etaBin] + etaBinning_trk[etaBin + 1]) / 2;  // Bin center

        double sf = getSF(effType, dummy, eta, dummy, idx, &numEff, &denEff); // Scale factor

        double statUp = getSF(effType, dummy, eta, dummy, 1, &numEff, &denEff) - sf;
        double statDown = sf - getSF(effType, dummy, eta, dummy, 2, &numEff, &denEff);

        double systUp = getSF(effType, dummy, eta, dummy, -1, &numEff, &denEff) * denEff;
        double systDown = getSF(effType, dummy, eta, dummy, -2, &numEff, &denEff) * denEff;

        double errUp = sqrt(statUp * statUp + systUp * systUp);
        double errDown = sqrt(statDown * statDown + systDown * systDown);

        /// X-axis errors (bin width) - assuming equal bin widths
        double xLow = (eta - etaBinning_trk[etaBin]);
        double xHigh = (etaBinning_trk[etaBin + 1] - eta);

        /// Set the point and errors in TGraphAsymmErrors
        sfNumEffGraph->SetPoint(etaBin, eta, numEff); // data TnP efficiency
        sfNumEffGraph->SetPointError(etaBin, xLow, xHigh, statDown * numEff, statUp * numEff); // No errors for data TnP efficiency

        sfDenEffGraph->SetPoint(etaBin, eta, denEff); // MC TnP efficiency
        sfDenEffGraph->SetPointError(etaBin, xLow, xHigh, 0, 0); // No errors for MC TnP efficiency

        sfGraph->SetPoint(etaBin, eta, sf);
        sfGraph->SetPointError(etaBin, xLow, xHigh, statDown, statUp); // Using statistical errors

        systEffGraph->SetPoint(etaBin, eta, systUp);
        systEffGraph->SetPointError(etaBin, xLow, xHigh, 0, 0); // No errors info for systematic uncertainties

        systRatioGraph->SetPoint(etaBin, eta, systUp / numEff);
        systRatioGraph->SetPointError(etaBin, xLow, xHigh, 0, 0); // No errors info for systematic uncertainties
    }

    TCanvas* canvasSFEff = new TCanvas("canvasSFEff", "", 600, 600);

    TPad* padEff = new TPad("padEff", "pad for efficiency", 0.0, 0.25, 1.0, 1.0);

    padEff->Draw();
    padEff->cd();

    padEff->SetRightMargin(0.05);
    padEff->SetLeftMargin(0.12);
    padEff->SetTopMargin(0.07);
    padEff->SetBottomMargin(0.02);

    frameHistEff->Draw();
    sfNumEffGraph->Draw("P SAME");
    sfDenEffGraph->Draw("P SAME");

    TLegend* legendEff = new TLegend(0.5, 0.2, 0.7, 0.5);
    legendEff->SetTextSize(0.045);
    legendEff->SetBorderSize(0);
    legendEff->SetFillStyle(0);
    legendEff->SetHeader("#splitline{Inner tracking Efficiency}{(p_{T}^{#mu} > 0.0 GeV)}");
    legendEff->AddEntry((TObject*)0, "", "");
    legendEff->AddEntry(sfDenEffGraph, "MC PYTHIA+HYDJET", "lep");
    legendEff->AddEntry(sfNumEffGraph, "DATA", "lep");
    legendEff->Draw("SAME");

    padEff->Update();
    // gPad->Update();
    canvasSFEff->Update();

    canvasSFEff->cd();

    TPad* padSF = new TPad("padSF", "pad for SF", 0.0, 0.0, 1.0, 0.25);

    padSF->Draw();
    padSF->cd();

    padSF->SetTopMargin(0.01);
    padSF->SetRightMargin(0.05);
    padSF->SetLeftMargin(0.12);
    padSF->SetBottomMargin(0.40);

    frameHistSF->Draw();
    sfGraph->Draw("P SAME");  // Draw asymmetric error bars on top

    drawLine(etaBinning_trk[0], 0.97, etaBinning_trk[nEtaBins_trk], 0.97);
    drawLine(etaBinning_trk[0], 0.98, etaBinning_trk[nEtaBins_trk], 0.98);
    drawLine(etaBinning_trk[0], 0.99, etaBinning_trk[nEtaBins_trk], 0.99);
    drawLine(etaBinning_trk[0], 1, etaBinning_trk[nEtaBins_trk], 1);
    drawLine(etaBinning_trk[0], 1.01, etaBinning_trk[nEtaBins_trk], 1.01);

    CMS_lumi(padEff, gCMSLumiText);

    gPad->Update();
    canvasSFEff->Update();

    gSystem->mkdir("SingleMuonEfficiency/", kTRUE);

    if (isANGraph == kFALSE) canvasSFEff->SaveAs(Form("SingleMuonEfficiency/SingleMuonEfficiency_SF_%s.png", effType.Data()));
    else canvasSFEff->SaveAs(Form("SingleMuonEfficiency/SingleMuonEfficiency_SF_%s_old.png", effType.Data()));

    /// draw the efficiency with the systematic variation
    TCanvas* canvasSystEff = new TCanvas("canvasSystEff", "", 600, 600);

    TPad* padSystEff = new TPad("padSystEff", "pad for systematic efficiency", 0.0, 0.25, 1.0, 1.0);

    padSystEff->Draw();
    padSystEff->cd();

    padSystEff->SetRightMargin(0.05);
    padSystEff->SetLeftMargin(0.12);
    padSystEff->SetTopMargin(0.07);
    padSystEff->SetBottomMargin(0.02);

    frameHistEff->Draw();
    sfNumEffGraph->Draw("P SAME"); 
    systEffGraph->Draw("P SAME");

	TLegend* legendSystEff = new TLegend(0.5, 0.2, 0.7, 0.5);
    legendSystEff->SetTextSize(0.045);
    legendSystEff->SetBorderSize(0);
	legendSystEff->SetFillStyle(0);
    legendSystEff->SetHeader("#splitline{Inner tracking Efficiency}{(p_{T}^{#mu} > 0.0 GeV)}");
    legendSystEff->AddEntry((TObject*)0, "", "");
    legendSystEff->AddEntry(sfNumEffGraph, "nominal", "lep");
    legendSystEff->AddEntry(systEffGraph, "max syst", "lep");
	legendSystEff->Draw("SAME");

    padSystEff->Update();
    // gPad->Update();
    canvasSystEff->Update();

    canvasSystEff->cd();

    TPad* padRatio = new TPad("padRatio", "pad for Ratio", 0.0, 0.0, 1.0, 0.25);

    padRatio->Draw();
    padRatio->cd();

    padRatio->SetTopMargin(0.01);
    padRatio->SetRightMargin(0.05);
    padRatio->SetLeftMargin(0.12);
    padRatio->SetBottomMargin(0.40);

    frameHistSF->GetYaxis()->SetTitle("Ratio");
    frameHistSF->Draw();
    systRatioGraph->Draw("P SAME");  // Draw asymmetric error bars on top

    drawLine(etaBinning_trk[0], 0.97, etaBinning_trk[nEtaBins_trk], 0.97);
    drawLine(etaBinning_trk[0], 0.98, etaBinning_trk[nEtaBins_trk], 0.98);
    drawLine(etaBinning_trk[0], 0.99, etaBinning_trk[nEtaBins_trk], 0.99);
    drawLine(etaBinning_trk[0], 1, etaBinning_trk[nEtaBins_trk], 1);
    drawLine(etaBinning_trk[0], 1.01, etaBinning_trk[nEtaBins_trk], 1.01);

    CMS_lumi(padSystEff, gCMSLumiText);
    
    gPad->Update();
    canvasSystEff->Update();

    if (isANGraph == kFALSE) canvasSystEff->SaveAs(Form("SingleMuonEfficiency/SingleMuonEfficiency_syst_%s.png", effType.Data()));
    else canvasSystEff->SaveAs(Form("SingleMuonEfficiency/SingleMuonEfficiency_syst_%s_old.png", effType.Data()));


    return;
}

/// scan drawSFEff
void SingleMuonTnPEff(Bool_t isANGraph = kFALSE) { // plots in the analysis note (AN-2018-316) if isANGraph = kTRUE (need to change the header file)

    const int nEtaBins = 4;
    double etaBinning[nEtaBins + 1] = {0, 1.2, 1.8, 2.1, 2.4};

    Bool_t isL3 = kTRUE;

    // const int nEffTypes = 3;
    // TString effType[3] = {"trg", "muid", "trk"};
    const int nEffTypes = 2;
    TString effType[nEffTypes] = {"trg", "muid"};

    for (int etaBin = 0; etaBin < nEtaBins; etaBin++) {
        for (int effTypeBin = 0; effTypeBin < nEffTypes; effTypeBin++) {
            if (effTypeBin == 0) {
                for (int filterBin = 0; filterBin < 2; filterBin++) {
                    if (filterBin == 0) isL3 = kTRUE; // Upsilon L3
                    else isL3 = kFALSE; // Upsilon L2  

                    drawSFEff(etaBinning[etaBin], etaBinning[etaBin + 1], effType[effTypeBin], isL3, isANGraph);
                }
            }

            else drawSFEff(etaBinning[etaBin], etaBinning[etaBin + 1], effType[effTypeBin], isL3, isANGraph);
        }    
    }

    drawSFEff_trk(isANGraph);

    return;
}