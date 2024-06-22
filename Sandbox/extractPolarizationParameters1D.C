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

/// Define and draw ideal 2D polarization distribution in (costheta, phi) phase space

TH2D* generate2DPolarization(Int_t nCosThetaBins, Float_t cosThetaMin, Float_t cosThetaMax, Int_t nPhiBins, Float_t phiMin, Float_t phiMax, Double_t nIn, Double_t lambdaThetaIn, Double_t lambdaPhiIn, Double_t lambdaThetaPhiIn){
	
	/// (dummy histogram to adjust the plot range)
	TH2F* hdummy = new TH2F("hdummy", ";cos #theta; #varphi (#circ);Number of generated #varUpsilons", nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins+5, phiMin-20, phiMax+100); 

	/// (2D Angular distribution function)
	TF2* polarFunc2D = generalPolarFunc(nIn);

	polarFunc2D->FixParameter(0, 1); // (parameter of normalization)

	polarFunc2D->FixParameter(1, lambdaThetaIn); // (parameter of lambda theta)
	polarFunc2D->FixParameter(2, lambdaPhiIn); // (parameter of lambda phi)
	polarFunc2D->FixParameter(3, lambdaThetaPhiIn); // (additional factor since y(#phi) is in the unit of degree)
	
	polarFunc2D->SetParNames("NormFactor", "lambdaTheta", "lambdaPhi",  "lambdaThetaPhi");
	polarFunc2D->SetTitle(";cos #theta; #varphi (#circ);Number of generated #varUpsilons");

	/// (histogram for Random sampling from the angular distribution function)
	TH2D *polarHist2D = new TH2D("polarHist2D", ";cos #theta; #varphi (#circ);Number of generated #varUpsilons", nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax); 

	polarHist2D->FillRandom("polarFunc2D", nIn); // convert the fuction to histogram for fit procedure

	// draw plots
	TCanvas *polarCanvas2D = new TCanvas("polarCanvas2D", "polarCanvas2D", 650, 600);

	hdummy->Draw("COLZ");

	polarHist2D->Draw("COLZ same");

	/// cosmetics	
	polarCanvas2D->SetLeftMargin(.15);
	polarCanvas2D->SetRightMargin(.20);

	hdummy->GetZaxis()->SetRangeUser(0, polarHist2D->GetMaximum());

	polarHist2D->GetZaxis()->SetMaxDigits(3);

	// Styles of the texts in the plot
	TLatex* legend = new TLatex();
	legend->SetTextAlign(22);
	legend->SetTextSize(0.05);

	// Put texts inside the plot
	legend->DrawLatexNDC(.48, .88, "Idea 2D polarization distribution");
	legend->DrawLatexNDC(.48, .80, Form("#lambda_{#theta} = %.2f, #lambda_{#varphi} = %.2f, #lambda_{#theta#varphi} = %.2f", polarFunc2D->GetParameter(1), polarFunc2D->GetParameter(2), polarFunc2D->GetParameter(3)));
	
	// Set the plot styles
	//gStyle->SetTitleYOffset(.9);
	gStyle->SetTitleOffset(1.5, "z"); // gap between color bar and z title
	SetColorPalette(gPreferredColorPaletteName);
	gStyle->SetNumberContours(256);

	gPad->Update();

	// save the plot
	gSystem->mkdir("DistributionFitsMC", kTRUE);
	polarCanvas2D->SaveAs(Form("DistributionFitsMC/IdealDistributionTheta%.2f_Phi%.2f.png", polarFunc2D->GetParameter(1), polarFunc2D->GetParameter(2)), "RECREATE");

	return polarHist2D;
}

void extractPolarizationParameters1D(Double_t lambdaThetaIn = 0.88, Double_t lambdaPhiIn = -0.8) {  

	/// Generate a Toy Data (This part can be replaced by data)
	
	// set binning and min, max of cosTheta and phi
	Int_t nCosThetaBins = 20;
	Float_t cosThetaMin = -1, cosThetaMax = 1;

	Int_t nPhiBins = 18;
	Float_t phiMin = -180, phiMax = 180;

	/// set input values
	Double_t nIn = 1e8; // normalization
	Double_t lambdaThetaPhiIn = 0;

	TH2D* polarHist2D = generate2DPolarization(nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax, nIn, lambdaThetaIn, lambdaPhiIn, lambdaThetaPhiIn);

	Double_t nEntries = polarHist2D->GetEntries();

	cout << "--------------------------------------" << endl;
	cout << "number of entries: " <<  nEntries << endl;
	cout << "--------------------------------------" << endl;

	/// Make 2D histograms to 1D 
	// (integrate over phi, that is, costheta graph)
	TH1D *polarHistCosTheta = polarHist2D->ProjectionX("cos #theta", 1, nPhiBins); // arguments: (name, firstybin, lastybin)

	// (integrate over cosTheta, that is, phi graph)
	TH1D *polarHistPhi = polarHist2D->ProjectionY("#varphi", 1, nCosThetaBins);

	/// Define variables
	// x and y axis
	RooRealVar cosThetaVar("cosTheta", "cos#theta", cosThetaMin, cosThetaMax); 
	RooRealVar phiVar("phi", "#varphi (#circ)", phiMin, phiMax); 

	// polarization parameters
	RooRealVar lambdaThetaVar("lambdaTheta", "lambda Theta", 0, -2, 2);
	RooRealVar lambdaPhiVar("lambdaPhi", "lambda Phi", 0, -2, 2);
	RooRealVar lambdaThetaPhiVar("lambdaThetaPhi", "lambda Theta Phi", 0);

	// normalization factor
	RooRealVar normCosThetaVar("normCosTheta", "normalization factor for costheta graph",  nIn, nIn * 0.1, nIn * 2.5); 	
	RooRealVar normPhiVar("normPhi", "normalization factor for phi graph", nIn, nIn * 0.1, nIn * 2.5);
	RooRealVar norm2DVar("norm2D", "normalization factor for 2D graph", nIn, nIn * 0.1, nIn * 2.5);

	/// Extract Polarization Parameters with Fit
	// Import histogram into RooFit
	RooDataHist rooHistCosTheta("rooHistCosTheta","angular distribution with cosTheta", cosThetaVar, polarHistCosTheta);
	
	RooDataHist rooHistPhi("rooHistPhi","angular distribution with phi", phiVar, polarHistPhi);
	
    RooArgList varList(cosThetaVar, phiVar);

	RooDataHist rooHist2D("rooHist2D", "rooHist2D", varList, polarHist2D);

	// Define model function and apply normalization factor
	CosThetaPolarizationPDF rooPdfCosTheta("rooPdfCosTheta", "rooPdfCosTheta", cosThetaVar, lambdaThetaVar);
	RooExtendPdf extendedPdfCosTheta("extendedPdfCosTheta", "extended CosTheta PDF", rooPdfCosTheta, normCosThetaVar);

	PhiPolarizationPDF rooPdfPhi("rooPdfPhi", "rooPdfPhi", phiVar, lambdaThetaVar, lambdaPhiVar);
	RooExtendPdf extendedPdfPhi("extendedPdfPhi", "extended Phi PDF", rooPdfPhi, normPhiVar);

	GeneralPolarizationPDF rooPdf2D("rooPdf2D", "rooPdf2D", cosThetaVar, phiVar, lambdaThetaVar, lambdaPhiVar, lambdaThetaPhiVar);
	RooExtendPdf extendedPdf2D("extendedPdf2D", "extendedPdf2D", rooPdf2D, norm2DVar);

	// Make frames
	RooPlot* frameCosTheta = cosThetaVar.frame(cosThetaMin, cosThetaMax, nCosThetaBins);

	RooPlot* framePhi = phiVar.frame(phiMin, phiMax, nPhiBins);

	// RooPlot* frame2D = 
	
	// locate imported histogram on the frame
	rooHistCosTheta.plotOn(frameCosTheta, Name("cosTheta Hist"), MarkerSize(1.5), DrawOption("P0Z"));

	rooHistPhi.plotOn(framePhi, Name("rooHistPhi"), MarkerSize(1.5), DrawOption("P0Z"));

	enableBinIntegrator(extendedPdfCosTheta, nCosThetaBins);

	// Fit the model to the histogram
	auto* fitResult2D = extendedPdf2D.fitTo(rooHist2D, Save(), Extended(kTRUE), Minos(kTRUE), NumCPU(3));
	fitResult2D->Print("v");


	auto* fitResultCosTheta = extendedPdfCosTheta.fitTo(rooHistCosTheta, Save(), Extended(kTRUE), Minos(kTRUE), NumCPU(3), Range(cosThetaMin, cosThetaMax));
	fitResultCosTheta->Print("v");

	// Put the fit model on the frame
	extendedPdfCosTheta.plotOn(frameCosTheta, Name("rooPdfCosTheta"));
	
	// Draw the histogram and the fit
	TCanvas* polarCanvas = new TCanvas(polarHist2D->GetName(), "", 1200, 600);
	
	polarCanvas->Divide(2, 1);

	polarCanvas->cd(1);
	
	TPad* padCosTheta = new TPad("padCosTheta", "padCosTheta", 0, 0.25, 1, 1.0);
	padCosTheta->SetBottomMargin(0.03);
	padCosTheta->SetTicks(1, 1);
	padCosTheta->Draw();
	padCosTheta->cd();
	
	// add legends
	frameCosTheta->addObject(PolarParamsText(lambdaThetaIn, lambdaPhiIn, normCosThetaVar, lambdaThetaVar, normPhiVar, lambdaPhiVar, false));
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
	lambdaThetaVar.setConstant(kTRUE) ;

	/// Extract Polarization Parameters with 1D Fit (phi)	

	// Import histogram into RooFit
	
	// Define model function

	// Make a frame

	enableBinIntegrator(extendedPdfPhi, nPhiBins);

	// Fit the model to the histogram
	auto* fitResultPhi = extendedPdfPhi.fitTo(rooHistPhi, Save(), Extended(kTRUE), Minos(kTRUE), NumCPU(3), Range(phiMin, phiMax));
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
	framePhi->addObject(PolarParamsText(lambdaThetaIn, lambdaPhiIn, normCosThetaVar, lambdaThetaVar, normPhiVar, lambdaPhiVar, true));
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
  	polarCanvas-> SaveAs(Form("DistributionFitsMC/fit1D_lambdaCosTheta%.2fPhi%.2f.png", lambdaThetaIn, lambdaPhiIn), "RECREATE");

	cout << "--------------------------------------" << endl;
	cout << "Done phi fit!" << endl;
	cout << "normalization: " << normPhiVar << ", lambdaPhi: " << lambdaPhiVar << endl;
	cout << "--------------------------------------" << endl;
	
	/// Print out the results
	cout << "input      : lambdaTheta = " << lambdaThetaIn << ", lambdaPhi = " << lambdaPhiIn << endl;
	cout << "fit results: lambdaTheta = " << lambdaThetaVar << ", lambdaPhi = " << lambdaPhiVar << endl;

	return;
	}