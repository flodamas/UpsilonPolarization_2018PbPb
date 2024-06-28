// This code extracts the polarization parameters 
// from the 2D angular distribution function (ideal MC distribution for the test) in the costheta phi space 
// using 1D fit

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
        double fitted = fitFunction->Eval(xValue); // Adjust index if necessary
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
	bottomPad->SetTicks(1, 1);
	bottomPad->Draw();
	bottomPad->cd();

	pullHist->SetTitle(" ");
	
	pullHist->GetYaxis()->SetTitleOffset(0.35);
	pullHist->GetYaxis()->SetTitle("Pull");
	
	pullHist->GetYaxis()->SetTitleSize(0.17);
	pullHist->GetYaxis()->SetLabelSize(0.15);
	pullHist->GetYaxis()->CenterTitle();

	pullHist->GetXaxis()->SetTitleSize(0.17);	
	pullHist->GetXaxis()->SetLabelSize(0.15);
	pullHist->GetXaxis()->CenterTitle();
	
	pullHist->GetYaxis()->SetTickSize(0.03);
	pullHist->GetXaxis()->SetTickSize(0.1);

	pullHist->SetMaximum(10.5);
	pullHist->SetMinimum(-10.5);
    
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
	TH2D *generalPolarHist = new TH2D("generalPolarHist", ";cos #theta; #varphi (#circ);Number of generated #varUpsilons", nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax); 

	generalPolarHist->FillRandom("generalPolarFunc", n); // convert the function to histogram for fit procedure

	// fit the histogram to see the change in the input values
	// TFitResultPtr fitResults = fitGeneralPolarizationHist(generalPolarHist);

	Float_t maxYield = generalPolarHist->GetEntries();

	TF2* generalPolarFuncFit = getGeneralPolarFunc(maxYield);

	TFitResultPtr fitResults = generalPolarHist->Fit("generalPolarFunc", "ESVIM0"); // chi2 fit to the integrated bin 

	// draw plots
	TCanvas *polarCanvas2D = new TCanvas("polarCanvas2D", "polarCanvas2D", 1250, 600);

	polarCanvas2D->Divide(2);

	polarCanvas2D->cd(1);

    gPad->SetRightMargin(0.18);

	hdummy->Draw("COLZ");

	generalPolarHist->Draw("COLZ SAME");

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

	generalPolarHist->Draw("LEGO SAME");

	generalPolarFuncFit->Draw("SURFACE SAME");

	/// cosmetics	
	polarCanvas2D->SetTopMargin(.15);
	polarCanvas2D->SetLeftMargin(.1);

	hdummy->GetZaxis()->SetRangeUser(0, generalPolarHist->GetMaximum());

	generalPolarHist->GetZaxis()->SetMaxDigits(3);

	// Styles of the texts in the plot
	TLatex* legend2 = new TLatex();

	// legend2->SetTextAlign(22);
	legend2->SetTextSize(0.05);

	// Put texts inside the plot
	legend2->DrawLatexNDC(.5, .90, Form("#lambda_{#theta, fit} = %.4f #pm %.4f", fitResults->Parameter(1), fitResults->ParError(1)));
	legend2->DrawLatexNDC(.5, .84, Form("#lambda_{#varphi, fit} = %.4f #pm %.4f", fitResults->Parameter(2), fitResults->ParError(2)));
	legend2->DrawLatexNDC(.5, .78, Form("#lambda_{#theta#varphi, fit} = %.4f #pm %.4f", fitResults->Parameter(3), fitResults->ParError(3)));

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
	polarCanvas2D->SaveAs(Form("DistributionFitsMC/IdealDistributionTheta%.2f_Phi%.2f.png", generalPolarFunc->GetParameter(1), generalPolarFunc->GetParameter(2)), "RECREATE");

	return generalPolarHist;
}

void extractPolarizationParameters1D(Double_t lambdaTheta0 = 0.88, Double_t lambdaPhi0 = -0.8) {  

	/// Generate a Toy Data (This part can be replaced by data)
	
	// set binning and min, max of cosTheta and phi
	Int_t nCosThetaBins = 20;
	Float_t cosThetaMin = -1, cosThetaMax = 1;

	Int_t nPhiBins = 18;
	Float_t phiMin = -180, phiMax = 180;

	/// set input values
	Double_t n0 = 1e7; // normalization
	Double_t lambdaThetaPhi0 = 0;

	// generate the data
	TH2D* generalPolarHist = generateGeneralPolarizationHist(nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax, n0, lambdaTheta0, lambdaPhi0, lambdaThetaPhi0);

	Double_t nEntries = generalPolarHist->GetEntries();

	cout << "--------------------------------------" << endl;
	cout << "number of entries: " <<  nEntries << endl;
	cout << "--------------------------------------" << endl;

	/// Make 2D histograms to 1D 
	// (integrate over phi, that is, costheta graph)
	TH1D *polarHistCosTheta = generalPolarHist->ProjectionX("cos #theta", 0, nPhiBins); // arguments: (name, firstybin, lastybin)

	// (integrate over cosTheta, that is, phi graph)
	TH1D *polarHistPhi = generalPolarHist->ProjectionY("#varphi", 0, nCosThetaBins);

	// /// Define variables
	// // x and y axis
	// RooRealVar cosThetaVar("cosTheta", "cos#theta", cosThetaMin, cosThetaMax); 
	// RooRealVar phiVar("phi", "#varphi (#circ)", phiMin, phiMax); 

	// // polarization parameters
	// RooRealVar lambdaThetaVar("lambdaTheta", "lambda Theta", 0, -2, 2);
	// RooRealVar lambdaPhiVar("lambdaPhi", "lambda Phi", 0, -2, 2);
	// RooRealVar lambdaThetaPhiVar("lambdaThetaPhi", "lambda Theta Phi", 0);

	// // normalization factor
	// RooRealVar normCosThetaVar("normCosTheta", "normalization factor for costheta graph",  n0, n0 * 0.1, n0 * 2.5); 	
	// RooRealVar normPhiVar("normPhi", "normalization factor for phi graph", 5e5, n0 * 0.001, n0 * 2.5);

	// /// Extract Polarization Parameters with Fit
	// // Import histogram into RooFit
	// RooDataHist rooHistCosTheta("rooHistCosTheta","angular distribution with cosTheta", cosThetaVar, polarHistCosTheta);
	
	// RooDataHist rooHistPhi("rooHistPhi","angular distribution with phi", phiVar, polarHistPhi);
	
    // RooArgList varList(cosThetaVar, phiVar);

	// // Define model function and apply normalization factor
	// CosThetaPolarizationPDF rooPdfCosTheta("rooPdfCosTheta", "rooPdfCosTheta", cosThetaVar, lambdaThetaVar);
	// RooExtendPdf extendedPdfCosTheta("extendedPdfCosTheta", "extended CosTheta PDF", rooPdfCosTheta, normCosThetaVar);

	// PhiPolarizationPDF rooPdfPhi("rooPdfPhi", "rooPdfPhi", phiVar, lambdaThetaVar, lambdaPhiVar);
	// RooExtendPdf extendedPdfPhi("extendedPdfPhi", "extended Phi PDF", rooPdfPhi, normPhiVar);

	// // Make frames
	// RooPlot* frameCosTheta = cosThetaVar.frame(cosThetaMin, cosThetaMax, nCosThetaBins);

	// RooPlot* framePhi = phiVar.frame(phiMin, phiMax, nPhiBins);
	
	// // locate imported histogram on the frame
	// rooHistCosTheta.plotOn(frameCosTheta, Name("cosTheta Hist"), MarkerSize(1.5), DrawOption("P0Z"));

	// rooHistPhi.plotOn(framePhi, Name("rooHistPhi"), MarkerSize(1.5), DrawOption("P0Z"));

	// // Fit the model to the histogram

	// enableBinIntegrator(extendedPdfCosTheta, nCosThetaBins);
	
	// auto* fitResultCosTheta = extendedPdfCosTheta.chi2FitTo(rooHistCosTheta, Save(), Extended(kTRUE), Minos(kTRUE), NumCPU(3), Range(cosThetaMin, cosThetaMax), SumW2Error(kFALSE), Integrate(kTRUE));
	// fitResultCosTheta->Print("v");

	// // Put the fit model on the frame
	// extendedPdfCosTheta.plotOn(frameCosTheta, Name("rooPdfCosTheta"));
	
	// // Draw the histogram and the fit
	// TCanvas* polarCanvas = new TCanvas("polarCanvas", "", 1200, 600);
	
	// polarCanvas->Divide(2, 1);

	// polarCanvas->cd(1);
	
	// TPad* padCosTheta = new TPad("padCosTheta", "padCosTheta", 0, 0.25, 1, 1.0);
	// padCosTheta->SetBottomMargin(0.03);
	// padCosTheta->SetTicks(1, 1);
	// padCosTheta->Draw();
	// padCosTheta->cd();
	
	// // add legends
	// frameCosTheta->addObject(PolarParamsText(lambdaTheta0, lambdaPhi0, normCosThetaVar, lambdaThetaVar, normPhiVar, lambdaPhiVar, false));
	// frameCosTheta->GetXaxis()->SetLabelOffset(1); // to make it disappear under the pull distribution pad
	// frameCosTheta->SetTitle("");
	// frameCosTheta->Draw();
	// gPad->RedrawAxis();

	// // pull Distribution
	// TPad* padCosThetaPull = GetPadPullDistribution(frameCosTheta, fitResultCosTheta->floatParsFinal().getSize());
  	
	// polarCanvas->cd(1);
  	
  	// padCosTheta->Draw();
  	// padCosThetaPull->Draw();

	// cout << "--------------------------------------" << endl;
	// cout << "Done costheta fit!" << endl;
	// cout << "normalization: " << normCosThetaVar << ", lambdaTheta: " << lambdaThetaVar << endl;
	// cout << "--------------------------------------" << endl;

	// // Fix lambda Theta and use this value for the phi graph fit
	// // lambdaThetaVar.setVal(0.88);
	// lambdaThetaVar.setConstant(kTRUE) ;

	// /// Extract Polarization Parameters with 1D Fit (phi)
	// enableBinIntegrator(extendedPdfPhi, nPhiBins);
	// // Fit the model to the histogram
	// auto* fitResultPhi = extendedPdfPhi.chi2FitTo(rooHistPhi, {Save(), Extended(kTRUE), Minos(kTRUE), NumCPU(3), Range(phiMin, phiMax), IntegrateBins(1000)/*, SumW2Error(kFALSE), IntegrateBins(10000)*/});
	// fitResultPhi->Print("v");
	
	// // Put the histogram on the frame
	// extendedPdfPhi.plotOn(framePhi, Name("rooPdfPhi"));

	// // Draw the histogram and the fit
	// polarCanvas->cd(2);

	// TPad *padPhi = new TPad("padPhi", "padPhi", 0, 0.25, 1, 1.0);
	// padPhi->SetBottomMargin(0.03);
	// padPhi->SetTicks(1, 1);
	// padPhi->Draw();
	// padPhi->cd();
	
	// // add legends
	// framePhi->addObject(PolarParamsText(lambdaTheta0, lambdaPhi0, normCosThetaVar, lambdaThetaVar, normPhiVar, lambdaPhiVar, true));
	// framePhi->GetXaxis()->SetLabelOffset(1); // to make it disappear under the pull distribution pad
	// framePhi->SetTitle("");
	// framePhi->Draw();
	// gPad->RedrawAxis();

	// // Pull Distribution
	// TPad* padPhiPull = GetPadPullDistribution(framePhi, fitResultPhi->floatParsFinal().getSize());

	// polarCanvas->cd(2);

  	// padPhi->Draw();
  	// padPhiPull->Draw();

  	// // save the canvas
  	// gSystem->mkdir("DistributionFitsMC", kTRUE);
  	// polarCanvas-> SaveAs(Form("DistributionFitsMC/RooFit_fit1D_lambdaCosTheta%.2fPhi%.2f.png", lambdaTheta0, lambdaPhi0), "RECREATE");

	// cout << "--------------------------------------" << endl;
	// cout << "Done phi fit!" << endl;
	// cout << "normalization: " << normPhiVar << ", lambdaPhi: " << lambdaPhiVar << endl;
	// cout << "--------------------------------------" << endl;
	
	// /// Print out the results
	// cout << "input      : lambdaTheta = " << lambdaTheta0 << ", lambdaPhi = " << lambdaPhi0 << endl;
	// cout << "fit results: lambdaTheta = " << lambdaThetaVar << " ± " << lambdaThetaVar.getError() << ", lambdaPhi = " << lambdaPhiVar << " ± " << lambdaPhiVar.getError() << endl;

	/// using ROOT::Fit()

	// perform fit of cosTheta graph 
	Float_t maxYield = polarHistCosTheta->GetEntries();

	TF1* cosThetaPolarFuncFit = getCosThetaPolarFunc(maxYield);

	TFitResultPtr cosThetafitResults = polarHistCosTheta->Fit("cosThetaPolarFunc", "ESVIMR0"); // chi2 fit to the integrated bin 

	cout << "--------------------------------------" << endl;
	cout << "Done costheta fit!" << endl;
	cout << "normalization: " << cosThetafitResults->Parameter(0) << ", lambdaTheta: " << cosThetafitResults->Parameter(1) << endl;
	cout << "--------------------------------------" << endl;

	cosThetafitResults->Print("v");

	TF1* phiPolarFuncFit = getPhiPolarFunc(maxYield);

	// fix the extracted lambdaTheta value
	// phiPolarFuncFit->SetParameter(1, cosThetafitResults->Parameter(1));
	phiPolarFuncFit->FixParameter(1, cosThetafitResults->Parameter(1));

	// perform fit of phi graph 
	TFitResultPtr phifitResults = polarHistPhi->Fit("phiPolarFunc", "ESVIMR0"); // chi2 fit to the integrated bin 

	cout << "--------------------------------------" << endl;
	cout << "Done phi fit!" << endl;
	cout << "normalization: " << phifitResults->Parameter(0) << ", lambdaPhi: " << phifitResults->Parameter(1) << endl;
	cout << "--------------------------------------" << endl;
	
	phifitResults->Print("v");

	// draw the histogram and the fit
	TCanvas* polarCanvas1D = new TCanvas("polarCanvas1D", "", 1200, 600);
	
	polarCanvas1D->Divide(2, 1);

	polarCanvas1D->cd(1);

	TPad* padCosTheta = new TPad("padCosTheta", "padCosTheta", 0, 0.25, 1, 1.0);
	padCosTheta->SetBottomMargin(0.03);
	padCosTheta->SetTicks(1, 1);
	padCosTheta->Draw();
	padCosTheta->cd();

	polarHistCosTheta->SetMinimum(0);

	polarHistCosTheta->SetMarkerStyle(8);
	polarHistCosTheta->SetMarkerSize(1.5);
	polarHistCosTheta->SetMarkerColor(kBlack);
	
	polarHistCosTheta->GetXaxis()->SetLabelSize(0);

	polarHistCosTheta->Draw("PE");

	cosThetaPolarFuncFit->Draw("SAME");

  	padCosTheta->Draw();

  	TPaveStats* cosThetaStatbox = positionStatBox(polarHistCosTheta);
  	cosThetaStatbox->Draw();

	polarCanvas1D->cd(1);

	// pull Distribution
	TPad* padCosThetaPull = GetPadPullDistribution(polarHistCosTheta, cosThetaPolarFuncFit, cosThetafitResults);

	// phi graph
	polarCanvas1D->cd(2);

	TPad* padPhi = new TPad("padPhi", "padPhi", 0, 0.25, 1, 1.0);
	padPhi->SetBottomMargin(0.03);
	padPhi->SetTicks(1, 1);
	padPhi->Draw();
	padPhi->cd();

	polarHistPhi->SetMinimum(0);

	polarHistPhi->SetMarkerStyle(8);
	polarHistPhi->SetMarkerSize(1.5);
	polarHistPhi->SetMarkerColor(kBlack);
	
	polarHistPhi->GetXaxis()->SetLabelSize(0);

	polarHistPhi->Draw("PE");

	phiPolarFuncFit->Draw("SAME");
  	
  	padPhi->Draw();

  	TPaveStats* phiStatbox = positionStatBox(polarHistPhi);
  	phiStatbox->Draw();

	polarCanvas1D->cd(2);

	// pull Distribution
	TPad* padPhiPull = GetPadPullDistribution(polarHistPhi, phiPolarFuncFit, phifitResults);

	gSystem->mkdir("DistributionFitsMC", kTRUE);
  	polarCanvas1D->SaveAs(Form("DistributionFitsMC/fit1D_lambdaCosTheta%.2fPhi%.2f.png", lambdaTheta0, lambdaPhi0), "RECREATE");

	return;
}