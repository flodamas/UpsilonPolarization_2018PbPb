// This code extracts the polarization parameters 
// from the 2D angular distribution function (ideal distribution for the test) in the costheta phi space 
// using 2D fit and RooFit


/// C++ 
#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>

/// ROOT
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TF1.h"
#include "TH2.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TRotation.h"
#include "TLorentzVector.h"
#include "TGraph.h"
#include "TString.h"
#include "TLatex.h"

///RooFit
#include "RooDataSet.h"
#include "RooWorkspace.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooPlot.h"
#include "RooCBShape.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooChebychev.h"
#include "RooGenericPdf.h"
#include "RooFitResult.h"
#include "RooHist.h"
#include "RooCFunction1Binding.h"
#include "RooCFunction3Binding.h"

///Reference files
#include "../Tools/Style/tdrStyle.C"
#include "../Tools/Style/CMS_lumi.C"
#include "../Tools/Style/FitDistributions.h"
#include "../Tools/Style/Legends.h"

using namespace RooFit;


void extractPolarizationParameters2D() {  
// void extractPolarizationParameters2D(Int_t iState = 1, Int_t ptMin = 0, Int_t ptMax = 30) {  

	/// Start measuring time 
	// clock_t start, end, cpu_time;
	// start = clock();
	auto start = std::chrono::high_resolution_clock::now();

	/// Put the text, "CMS Internal", on the right top of the plot 
	// writeExtraText = true; // if extra text
	// extraText = "       Internal";

	writeExtraText = false; // if don't want to put a text

	/// Generate a Toy Data (This part can be replaced by data)
	// (cos theta, phi) 2D distribution maps for CS and HX frames
	Int_t nCosThetaBins = 20;
	Float_t cosThetaMin = -1, cosThetaMax = 1;

	//	Double_t cosThetaBinning[] = {-1, -0.8, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1};
	//	nCosThetaBins = sizeof(cosThetaBinning) / sizeof(Double_t) - 1;

	Int_t nPhiBins = 18;
	Float_t phiMin = -180, phiMax = 180;

	//Double_t phiBinning[] = {-TMath::Pi(), -2.5, -2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2, 2.5, 180};
	//	nPhiBins = sizeof(cosThetaBinning) / sizeof(Double_t) - 1;

	/// Set the plot styles
	gStyle -> SetPadLeftMargin(.15);
	gStyle -> SetPadRightMargin(.2);
	//gStyle->SetTitleYOffset(.9);
	gStyle -> SetTitleOffset(1.2,"z"); // gap between color bar and z title
	gStyle -> SetPalette(kCividis);
	gStyle -> SetNumberContours(256);

	/// Styles of the texts in the plot
	TLatex* legend = new TLatex();
	legend -> SetTextAlign(22);
	legend -> SetTextSize(0.05);

	/// Define variables 
	RooRealVar costheta("costheta", "cos #theta", cosThetaMin, cosThetaMax);
	RooRealVar phi("phi", "#phi", phiMin, phiMax);
	RooRealVar norm("norm", "normalization factor", 0, 5);
	RooRealVar lambTheta("lambTheta", "#lambda_{#theta}", 0, -1, 1);
	RooRealVar lambPhi("lambPhi", "#lambda_{#phi}", 0, -1, 1);
	RooRealVar pi("pi", "#pi", M_PI);


	/// Define and draw 2D Angular distribution function in (costheta, phi) phase space
	/// (I'm testing 2D fit with an ideal case. This part could be replaced with MC or Data later)	
    ///----------------------------------------------------------------------------------------
	const int NumIdx = 6;
	double weights[NumIdx][2][2] = {{{0.88, -0.99}, {0.00, -0.50}}, //test1 --idx:0
							   {{-1.00, 0.00}, {-1.00, 0.00}}, //test2 --idx:1
							   {{1.00, 1.00}, {1.00, 1.00}}, //start1 --idx:2
							   {{0.00, 0.00}, {0.00, 0.00}}, //start2 --idx:3
							   {{1.00, 0.00}, {1.00, 0.00}}, //start3 --idx:4
							   {{0.00, -0.50}, {0.00, -0.50}}}; //start4 --idx:5

	int chooseIdx = 0; // choose the index (0-5) above
	int HXorCS = 0; // HX is 0, CS is 1

	RooRealVar lambThetaTest("lambThetaTest", "lambda theta for generating data", weights[chooseIdx][HXorCS][0]);
	RooRealVar lambPhiTest("lambPhiTest", "lambda phi for generating data", weights[chooseIdx][HXorCS][1]);

	/// Write the function with the polarization parameters chosen above
	RooGenericPdf AngDisFuncPdfData("AngDisFuncPdfData", "norm*(1.+lambThetaTest*costheta*costheta+lambPhiTest*(1.-costheta*costheta)*(cos(2.*phi*pi/180.)))/(3+lambThetaTest)", RooArgSet(costheta, phi, norm, lambThetaTest, lambPhiTest, pi)); 
	/// Generate the data from the function
	RooDataSet* IdealData = AngDisFuncPdfData.generate(RooArgSet(costheta, phi), 1e7);

	// /// Read GenOnly Nofilter file with HX and CS frames (Run "makeHXCSNtuple.C" to get this file)
	// const char* filename = Form("./AcceptanceMaps/%dS/withWeights/AcceptanceResults_withWeights_pt%dto%dGeV.root", iState, ptMin, ptMax);
	// TFile* file = TFile::Open(filename, "READ");
	// if (!file) {
	// 	cout << "File " << filename << " not found. Check the directory of the file." << endl;
	// 	return;
	// }

	// cout << "File " << filename << " opened" << endl;


	// const char* filename = "../Files/AngDisHistTheta0.00Phi-0.50.root";
	// TFile* file = TFile::Open(filename, "READ");
	// if (!file) {
	// 	cout << "File " << filename << " not found. Check the directory of the file." << endl;
	// 	return;
	// }

	// cout << "File " << filename << " opened" << endl;

    ///----------------------------------------------------------------------------------------


	/// Define model function
	RooGenericPdf AngDisFuncPdf("AngDisFuncPdf", "norm*(1.+lambTheta*costheta*costheta+lambPhi*(1.-costheta*costheta)*(cos(2.*phi*pi/180.)))/(3+lambTheta)", RooArgSet(costheta, phi, norm, lambTheta, lambPhi, pi)); 
	AngDisFuncPdf.getVal(RooArgSet(costheta,phi)); // this line is needed for identifying axis variables (x, y)

	/// Fit the model to the data
	AngDisFuncPdf.fitTo(*IdealData) ;

	// Plot the x distribution of data(x,y) and f(x,y)
	RooPlot* framex = costheta.frame();
	IdealData->plotOn(framex);
	AngDisFuncPdf.plotOn(framex);
	TCanvas* cx = new TCanvas("cx", "cx", 700, 600);
	framex->Draw();

	// Plot the y distribution of data(x,y) and f(x,y)
	RooPlot* framey = phi.frame();
	IdealData->plotOn(framey);
	AngDisFuncPdf.plotOn(framey);
	TCanvas* cy = new TCanvas("cy", "cy", 700, 600);
	framey -> Draw();

	// Create histogram for the data
	TCanvas* c2Ddata = new TCanvas("c2Ddata", "c2Ddata", 700, 600);
	TH2* h2Ddata = dynamic_cast<TH2*>(IdealData->createHistogram("h2Ddata", costheta, Binning(nCosThetaBins,cosThetaMin,cosThetaMax), YVar(phi, Binning(nPhiBins,phiMin,phiMax))));

	// Draw histogram
	h2Ddata -> Draw("COLZ");
	// h2Ddata -> Draw("Surf"); //mesh 3D plot


	// Create histogram for the fit
	TCanvas* c2Dfit = new TCanvas("c2Dfit", "c2Dfit", 700, 600);
	TH2* h2Dfit = dynamic_cast<TH2*>(AngDisFuncPdf.createHistogram("h2Dfit", costheta, Binning(nCosThetaBins,cosThetaMin,cosThetaMax), YVar(phi, Binning(nPhiBins,phiMin,phiMax))));

	// Draw histogram
	h2Dfit -> Draw("COLZ");
	// h2Dfit -> Draw("Surf"); //mesh 3D plot

	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> duration = end - start;
	cout << "\nEllapsed Time: " << duration.count() / 60.0 << " mins (" << duration.count() << " secs)\n";

	return;



	}