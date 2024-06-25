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
	TH1D *polarHistCosTheta = generalPolarHist->ProjectionX("cos #theta", 1, nPhiBins); // arguments: (name, firstybin, lastybin)

	// (integrate over cosTheta, that is, phi graph)
	TH1D *polarHistPhi = generalPolarHist->ProjectionY("#varphi", 1, nCosThetaBins);

	/// Define variables
	// x and y axis
	RooRealVar cosThetaVar("cosTheta", "cos#theta", cosThetaMin, cosThetaMax); 
	RooRealVar phiVar("phi", "#varphi (#circ)", phiMin, phiMax); 

	// polarization parameters
	RooRealVar lambdaThetaVar("lambdaTheta", "lambda Theta", 0, -2, 2);
	RooRealVar lambdaPhiVar("lambdaPhi", "lambda Phi", 0, -2, 2);
	RooRealVar lambdaThetaPhiVar("lambdaThetaPhi", "lambda Theta Phi", 0);

	// normalization factor
	RooRealVar normCosThetaVar("normCosTheta", "normalization factor for costheta graph",  n0, n0 * 0.1, n0 * 2.5); 	
	RooRealVar normPhiVar("normPhi", "normalization factor for phi graph", n0, n0 * 0.1, n0 * 2.5);

	/// Extract Polarization Parameters with Fit
	// Import histogram into RooFit
	RooDataHist rooHistCosTheta("rooHistCosTheta","angular distribution with cosTheta", cosThetaVar, polarHistCosTheta);
	
	RooDataHist rooHistPhi("rooHistPhi","angular distribution with phi", phiVar, polarHistPhi);
	
    RooArgList varList(cosThetaVar, phiVar);

	// Define model function and apply normalization factor
	CosThetaPolarizationPDF rooPdfCosTheta("rooPdfCosTheta", "rooPdfCosTheta", cosThetaVar, lambdaThetaVar);
	RooExtendPdf extendedPdfCosTheta("extendedPdfCosTheta", "extended CosTheta PDF", rooPdfCosTheta, normCosThetaVar);

	PhiPolarizationPDF rooPdfPhi("rooPdfPhi", "rooPdfPhi", phiVar, lambdaThetaVar, lambdaPhiVar);
	RooExtendPdf extendedPdfPhi("extendedPdfPhi", "extended Phi PDF", rooPdfPhi, normPhiVar);

	// Make frames
	RooPlot* frameCosTheta = cosThetaVar.frame(cosThetaMin, cosThetaMax, nCosThetaBins);

	RooPlot* framePhi = phiVar.frame(phiMin, phiMax, nPhiBins);
	
	// locate imported histogram on the frame
	rooHistCosTheta.plotOn(frameCosTheta, Name("cosTheta Hist"), MarkerSize(1.5), DrawOption("P0Z"));

	rooHistPhi.plotOn(framePhi, Name("rooHistPhi"), MarkerSize(1.5), DrawOption("P0Z"));

	// Fit the model to the histogram

	auto* fitResultCosTheta = extendedPdfCosTheta.chi2FitTo(rooHistCosTheta, Save(), Extended(kTRUE), Minos(kTRUE), NumCPU(3), Range(cosThetaMin, cosThetaMax), SumW2Error(kFALSE), IntegrateBins(10));
	fitResultCosTheta->Print("v");

	// Put the fit model on the frame
	extendedPdfCosTheta.plotOn(frameCosTheta, Name("rooPdfCosTheta"));
	
	// Draw the histogram and the fit
	TCanvas* polarCanvas = new TCanvas(generalPolarHist->GetName(), "", 1200, 600);
	
	polarCanvas->Divide(2, 1);

	polarCanvas->cd(1);
	
	TPad* padCosTheta = new TPad("padCosTheta", "padCosTheta", 0, 0.25, 1, 1.0);
	padCosTheta->SetBottomMargin(0.03);
	padCosTheta->SetTicks(1, 1);
	padCosTheta->Draw();
	padCosTheta->cd();
	
	// add legends
	frameCosTheta->addObject(PolarParamsText(lambdaTheta0, lambdaPhi0, normCosThetaVar, lambdaThetaVar, normPhiVar, lambdaPhiVar, false));
	frameCosTheta->GetXaxis()->SetLabelOffset(1); // to make it disappear under the pull distribution pad
	frameCosTheta->SetTitle("");
	frameCosTheta->Draw();
	gPad->RedrawAxis();

	// pull Distribution
	TPad* padCosThetaPull = GetPadPullDistribution(frameCosTheta, fitResultCosTheta->floatParsFinal().getSize());
  	
	polarCanvas->cd(1);
  	
  	padCosTheta->Draw();
  	padCosThetaPull->Draw();

	cout << "--------------------------------------" << endl;
	cout << "Done costheta fit!" << endl;
	cout << "normalization: " << normCosThetaVar << ", lambdaTheta: " << lambdaThetaVar << endl;
	cout << "--------------------------------------" << endl;

	// Fix lambda Theta and use this value for the phi graph fit
	// lambdaThetaVar.setConstant(kTRUE) ;

	/// Extract Polarization Parameters with 1D Fit (phi)

	// Fit the model to the histogram
	auto* fitResultPhi = extendedPdfPhi.chi2FitTo(rooHistPhi, Save(), Extended(kTRUE), Minos(kTRUE), NumCPU(3), Range(phiMin, phiMax), SumW2Error(kFALSE), IntegrateBins(10));
	fitResultPhi->Print("v");
	
	// Put the histogram on the frame
	extendedPdfPhi.plotOn(framePhi, Name("rooPdfPhi"));

	// Draw the histogram and the fit
	polarCanvas->cd(2);

	TPad *padPhi = new TPad("padPhi", "padPhi", 0, 0.25, 1, 1.0);
	padPhi->SetBottomMargin(0.03);
	padPhi->SetTicks(1, 1);
	padPhi->Draw();
	padPhi->cd();
	
	// add legends
	framePhi->addObject(PolarParamsText(lambdaTheta0, lambdaPhi0, normCosThetaVar, lambdaThetaVar, normPhiVar, lambdaPhiVar, true));
	framePhi->GetXaxis()->SetLabelOffset(1); // to make it disappear under the pull distribution pad
	framePhi->SetTitle("");
	framePhi->Draw();
	gPad->RedrawAxis();

	// Pull Distribution
	TPad* padPhiPull = GetPadPullDistribution(framePhi, fitResultPhi->floatParsFinal().getSize());

	polarCanvas->cd(2);

  	padPhi->Draw();
  	padPhiPull->Draw();

  	// save the canvas
  	gSystem->mkdir("DistributionFitsMC", kTRUE);
  	polarCanvas-> SaveAs(Form("DistributionFitsMC/fit1D_lambdaCosTheta%.2fPhi%.2f.png", lambdaTheta0, lambdaPhi0), "RECREATE");

	cout << "--------------------------------------" << endl;
	cout << "Done phi fit!" << endl;
	cout << "normalization: " << normPhiVar << ", lambdaPhi: " << lambdaPhiVar << endl;
	cout << "--------------------------------------" << endl;
	
	/// Print out the results
	cout << "input      : lambdaTheta = " << lambdaTheta0 << ", lambdaPhi = " << lambdaPhi0 << endl;
	cout << "fit results: lambdaTheta = " << lambdaThetaVar << " ± " << lambdaThetaVar.getError() << ", lambdaPhi = " << lambdaPhiVar << " ± " << lambdaPhiVar.getError() << endl;

	// using ROOT::Fit()

	// Float_t maxYield = polarHistCosTheta->GetEntries();

	// TF1* cosThetaPolarFuncFit = getCosThetaPolarFunc(maxYield);

	// TFitResultPtr cosThetafitResults = polarHistCosTheta->Fit("cosThetaPolarFunc", "ESVIM0"); // chi2 fit to the integrated bin 

	// cout << "--------------------------------------" << endl;
	// cout << "Done costheta fit!" << endl;
	// cout << "normalization: " << cosThetafitResults->Parameter(0) << ", lambdaTheta: " << cosThetafitResults->Parameter(1) << endl;
	// cout << "--------------------------------------" << endl;

	// cosThetafitResults->Print("v");

	// TF1* phiPolarFuncFit = getPhiPolarFunc(maxYield);

	// phiPolarFuncFit->SetParameter(1, cosThetafitResults->Parameter(1));

	// TFitResultPtr phifitResults = polarHistPhi->Fit("phiPolarFunc", "ESVIM0"); // chi2 fit to the integrated bin 

	// cout << "--------------------------------------" << endl;
	// cout << "Done phi fit!" << endl;
	// cout << "normalization: " << phifitResults->Parameter(0) << ", lambdaPhi: " << phifitResults->Parameter(1) << endl;
	// cout << "--------------------------------------" << endl;
	
	// phifitResults->Print("v");


	return;
	}