// This macro tests the 2D polarization fit, skipping the acceptance forbiden region. 

#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/FitShortcuts.h"
#include "../Tools/Style/Legends.h"

#include "../Tools/Style/FitDistributions.h"

#include "../Tools/RooFitPDFs/PolarFunc.h"

#include "../Tools/RooFitPDFs/CosThetaPolarizationPDF.h"
#include "../Tools/RooFitPDFs/CosThetaPolarizationPDF.cxx"

#include "../Tools/RooFitPDFs/PhiPolarizationPDF.h"
#include "../Tools/RooFitPDFs/PhiPolarizationPDF.cxx"

#include "../Tools/RooFitPDFs/GeneralPolarizationPDF.h"
#include "../Tools/RooFitPDFs/GeneralPolarizationPDF.cxx"

using namespace RooFit;

TPaveStats* positionStatBox(TH1* hist) {

	// Access the stat box
    auto *statBox = dynamic_cast<TPaveStats*>(hist->FindObject("stats"));

    if (statBox) {
        statBox->SetX1NDC(0.26); // New X1 coordinate in NDC
        statBox->SetX2NDC(0.87); // New X2 coordinate in NDC
        statBox->SetY1NDC(0.09); // New Y1 coordinate in NDC
        statBox->SetY2NDC(0.32); // New Y2 coordinate in NDC
        statBox->SetLineColor(0);
    }
    else cout << "not found stat box!!" << endl;

    return statBox;
}

TPad* GetPadPullDistribution(TH1* dataHist, TF1* fitFunction, TFitResultPtr fitResult) {

    Int_t nBins = dataHist->GetNbinsX();
    
    TH1D* pullHist = new TH1D("pullHist", "Pull Distribution", nBins, dataHist->GetXaxis()->GetXmin(), dataHist->GetXaxis()->GetXmax());

    for (int iBin=1; iBin<=nBins; iBin++) {
    	double xValue = dataHist->GetXaxis()->GetBinCenter(iBin);
        double observed = dataHist->GetBinContent(iBin);
        double fitted = fitFunction->Eval(xValue);
        double error = dataHist->GetBinError(iBin);

        if (error != 0) {
            double pull = (observed - fitted) / error;
            pullHist->SetBinContent(iBin, pull);
        } else {
            pullHist->SetBinContent(iBin, 0);
        }
    }

    TPad* bottomPad = new TPad("bottomPad", "bottomPad", 0, 0.0, 1, .25);
    bottomPad->SetTopMargin(0.015);
	bottomPad->SetBottomMargin(0.4);
	bottomPad->SetRightMargin(0.01);
	bottomPad->SetLeftMargin(0.13);
	bottomPad->SetTicks(1, 1);
	bottomPad->Draw();
	bottomPad->cd();

	pullHist->SetTitle(" ");
	
	pullHist->GetYaxis()->SetTitleOffset(0.3);
	pullHist->GetYaxis()->SetTitle("Pull");
	
	pullHist->GetYaxis()->SetTitleSize(0.17);
	pullHist->GetYaxis()->SetLabelSize(0.15);
	pullHist->GetYaxis()->CenterTitle();

	pullHist->GetXaxis()->SetTitle(dataHist->GetXaxis()->GetTitle());
	pullHist->GetXaxis()->SetTitleSize(0.17);	
	pullHist->GetXaxis()->SetLabelSize(0.15);
	pullHist->GetXaxis()->CenterTitle();
	
	pullHist->GetYaxis()->SetTickSize(0.03);
	pullHist->GetXaxis()->SetTickSize(0.1);

	pullHist->SetMaximum(9.5);
	pullHist->SetMinimum(-9.5);
    
	pullHist->SetMarkerStyle(8);
	pullHist->SetMarkerSize(1);
	pullHist->SetMarkerColor(kBlack);

    pullHist->Draw("PE");

    TLine zeroLine(bottomPad->GetUxmin(), 0, bottomPad->GetUxmax(), 0);
	zeroLine.SetLineStyle(kDashed);
	zeroLine.Draw("SAME");

	bottomPad->Update();

	TLatex textChi2;
	textChi2.SetTextAlign(12);
	textChi2.SetTextSize(0.15);
	textChi2.DrawLatexNDC(0.75, 0.12, Form("#chi^{2} / n_{dof} = %.1f", fitResult->Chi2() / fitResult->Ndf()));

    return bottomPad;
}


TH2D* generateGeneralPolarizationHist(Int_t nCosThetaBins, Float_t cosThetaMin, Float_t cosThetaMax, Int_t nPhiBins, Float_t phiMin, Float_t phiMax, Double_t n, Double_t lambdaTheta, Double_t lambdaPhi, Double_t lambdaThetaPhi){
	
	/// dummy histogram to adjust the plot range
	TH2F* hdummy = new TH2F("hdummy", ";cos #theta; #varphi (#circ);Number of generated #varUpsilons", nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins + 5, phiMin - 20, phiMax + 100); 

	/// 2D Angular distribution function
	TF2* generalPolarFunc = getGeneralPolarFunc(n); // no need to be n for the argument

	generalPolarFunc->FixParameter(0, 1); // (parameter of normalization)

	generalPolarFunc->FixParameter(1, lambdaTheta); // (input parameter of lambda theta)
	generalPolarFunc->FixParameter(2, lambdaPhi); // (input parameter of lambda phi)
	generalPolarFunc->FixParameter(3, lambdaThetaPhi); // (input parameter of lambda theta phi)
	
	generalPolarFunc->SetTitle(";cos #theta; #varphi (#circ);Number of generated #varUpsilons");

	/// histogram for Random sampling from the angular distribution function
	TH2D* generalPolarHist = new TH2D("generalPolarHist", ";cos #theta; #varphi (#circ);Number of generated #varUpsilons", nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax); 

	generalPolarHist->FillRandom("generalPolarFunc", n); // convert the function to histogram for fit procedure

	// TH2D* discontinuousHist = new TH2D("discontinuousHist", ";cos #theta; #varphi (#circ);Number of generated #varUpsilons", nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax); 

	TH2D* discontinuousHist = (TH2D*) generalPolarHist->Clone();

	Int_t iEmptyBin = 3;
	
	for (Int_t iCosTheta = 1; iCosTheta <= nCosThetaBins; iCosTheta++) {
		discontinuousHist->SetBinContent(iCosTheta, iEmptyBin, 0);
	}

	// fit the histogram to see the change in the input values
	// TFitResultPtr fitResults = fitGeneralPolarizationHist(generalPolarHist);

	Float_t maxYield = discontinuousHist->GetEntries();

	// TF2* generalPolarFuncFit = getGeneralPolarFunc(maxYield);

	TF2* discontinuousFuncFit = getDiscontPolarFunc(maxYield);

	TFitResultPtr fitResults = discontinuousHist->Fit("generalPolarFunc", "ESVIM0"); // chi2 fit to the integrated bin 

	// draw plots
	TCanvas *polarCanvas2D = new TCanvas("polarCanvas2D", "polarCanvas2D", 1250, 600);

	polarCanvas2D->Divide(2);

	polarCanvas2D->cd(1);

    gPad->SetRightMargin(0.18);

	hdummy->Draw("COLZ");

	discontinuousHist->Draw("COLZ SAME");

	// Styles of the texts in the plot
	TLatex* legend1 = new TLatex();
	legend1->SetTextAlign(22);
	legend1->SetTextSize(0.05);

	// Put texts inside the plot
	legend1->DrawLatexNDC(.50, .88, "Idea 2D polarization distribution");
	legend1->DrawLatexNDC(.50, .80, Form("#lambda_{#theta} = %.2f, #lambda_{#varphi} = %.2f, #lambda_{#theta#varphi} = %.2f", generalPolarFunc->GetParameter(1), generalPolarFunc->GetParameter(2), generalPolarFunc->GetParameter(3)));

	polarCanvas2D->cd(2);

	gPad->SetTopMargin(0.05);

	hdummy->Draw("LEGO");

	discontinuousHist->Draw("LEGO SAME");

	discontinuousFuncFit->Draw("SURFACE SAME");

	/// cosmetics	
	polarCanvas2D->SetTopMargin(.15);
	polarCanvas2D->SetLeftMargin(.1);

	Double_t maxYaxis = discontinuousHist->GetMaximum() * 1.4;

	hdummy->GetZaxis()->SetRangeUser(0, maxYaxis);

	discontinuousHist->GetZaxis()->SetMaxDigits(3);
	discontinuousHist->GetZaxis()->SetRangeUser(0, maxYaxis);

	discontinuousFuncFit->SetRange(-1, -180, 0, 1, 180, maxYaxis);
	discontinuousFuncFit->SetMaximum(maxYaxis);
	discontinuousFuncFit->SetMinimum(0);

	// Styles of the texts in the plot
	TLatex* legend2 = new TLatex();

	// legend2->SetTextAlign(22);
	legend2->SetTextSize(0.05);

	// Put texts inside the plot
	legend2->DrawLatexNDC(.5, .90, Form("#lambda_{#theta, fit} = %.4f #pm %.4f", fitResults->Parameter(1), fitResults->ParError(1)));
	legend2->DrawLatexNDC(.5, .84, Form("#lambda_{#varphi, fit} = %.4f #pm %.4f", fitResults->Parameter(2), fitResults->ParError(2)));
	legend2->DrawLatexNDC(.5, .78, Form("#lambda_{#theta#varphi, fit} = %.4f #pm %.4f", fitResults->Parameter(3), fitResults->ParError(3)));

	double chi2 = fitResults->Chi2();
	double nDOF = nCosThetaBins * (nPhiBins - 1) - discontinuousFuncFit->GetNpar();

	TLatex textChi2;
	textChi2.SetTextAlign(12);
	textChi2.SetTextSize(0.045);
	textChi2.DrawLatexNDC(0.75, 0.08, Form("#chi^{2} / n_{dof} = %.1f", fitResults->Chi2() / fitResults->Ndf()));

	cout << "ndf from fitresult: " << fitResults->Ndf() << endl;
	cout << "my calculation: " << nDOF << endl;

	// Set the plot styles
	hdummy->GetZaxis()->SetTitleOffset(1.);
	hdummy->GetZaxis()->SetMaxDigits(3);

	hdummy->GetXaxis()->SetTitleOffset(1.);
	hdummy->GetXaxis()->CenterTitle();

	hdummy->GetYaxis()->SetTitleOffset(1.5);
	hdummy->GetYaxis()->CenterTitle();

	SetColorPalette(gPreferredColorPaletteName);

	gPad->Update();

	polarCanvas2D->Update();

	// save the plot
	gSystem->mkdir("DistributionFitsMC", kTRUE);
	polarCanvas2D->SaveAs(Form("DistributionFitsMC/IdealDiscontinuousDistributionTheta%.2f_Phi%.2f.png", generalPolarFunc->GetParameter(1), generalPolarFunc->GetParameter(2)), "RECREATE");

	return generalPolarHist;
}

void testPolarizationFit() {
	
	Int_t nCosThetaBins = 5;
	Float_t cosThetaMin = -0.7, cosThetaMax = 0.7;

	Int_t nPhiBins = 5;
	Float_t phiMin = -180, phiMax = 180;

	Double_t n = 1e7, lambdaTheta = 0., lambdaPhi = 0.1, lambdaThetaPhi = 0;

	generateGeneralPolarizationHist(nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax, n, lambdaTheta, lambdaPhi, lambdaThetaPhi);

}