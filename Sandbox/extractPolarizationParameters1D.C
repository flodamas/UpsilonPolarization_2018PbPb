// This code extracts the polarization parameters 
// from the 2D angular distribution function (ideal distribution for the test) in the costheta phi space 
// using 1D fit

/// C++ 
#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>

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
#include "TAxis.h"

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

void extractPolarizationParameters1D() {  
// void extractPolarizationParameters1D(Int_t iState = 1, Int_t ptMin = 0, Int_t ptMax = 30) {  

	/// Start measuring time 
	clock_t start, end, cpu_time;
	start = clock();


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
	//gStyle -> SetTitleYOffset(.9);
	gStyle -> SetTitleOffset(1.5,"z"); // gap between color bar and z title
	gStyle -> SetPalette(kCividis);
	gStyle -> SetNumberContours(256);

	/// Define and draw 2D Angular distribution function in (costheta, phi) phase space
	/// (I'm testing 2D fit with an ideal case. This part could be replaced with MC or Data later)	
    ///----------------------------------------------------------------------------------------
	/// (dummy histogram to adjust the plot range)
	TH2F *hdummy = new TH2F("hdummy", ";cos #theta; #varphi (#circ);Number of generated #varUpsilons", nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins+5, phiMin-20, phiMax+100); 

	/// Set the weights
	// Test1 values: λ_θ(HX)=0.88, λ_φ(HX)=-0.99, λ_θ(CS)=0.00, λ_φ(CS)=-0.50
	// Test2 values: λ_θ(HX)=-1.00, λ_φ(HX)=0.00, λ_θ(CS)=-1.00, λ_φ(CS)=0.00	

	// Start1 values: λ_θ(HX)=1.00, λ_φ(HX)=1.00, λ_θ(CS)=1.00, λ_φ(CS)=1.00
	// Start2 values: λ_θ(HX)=0.00, λ_φ(HX)=0.00, λ_θ(CS)=0.00, λ_φ(CS)=0.00
	// Start3 values: λ_θ(HX)=1.00, λ_φ(HX)=0.00, λ_θ(CS)=1.00, λ_φ(CS)=0.00
	// Start4 values: λ_θ(HX)=0.00, λ_φ(HX)=-0.50, λ_θ(CS)=0.00, λ_φ(CS)=-0.50	

	const int NumIdx = 6;
	double weights[NumIdx][2][2] = {{{0.88, -0.99}, {0.00, -0.50}}, //test1 --idx:0
							   {{-1.00, 0.00}, {-1.00, 0.00}}, //test2 --idx:1
							   {{1.00, 1.00}, {1.00, 1.00}}, //start1 --idx:2
							   {{0.00, 0.00}, {0.00, 0.00}}, //start2 --idx:3
							   {{1.00, 0.00}, {1.00, 0.00}}, //start3 --idx:4
							   {{0.00, -0.50}, {0.00, -0.50}}}; //start4 --idx:5

	int chooseIdx = 0; // choose the index (0-5) above
	int HXorCS = 1; // HX is 0, CS is 1

	/// (2D Angular distribution function)
	TF2  *AngDisFunc = new TF2("AngDisFunc","[2]/(3+[0])*(1+[0]*x*x+[1]*(1-x*x)*(cos(2*y*[3]/180.)))", cosThetaMin, cosThetaMax, phiMin, phiMax); 
	AngDisFunc -> FixParameter(0, weights[chooseIdx][HXorCS][0]); // (parameter of lambda theta)
	AngDisFunc -> FixParameter(1, weights[chooseIdx][HXorCS][1]); // (parameter of lambda phi)
	AngDisFunc -> FixParameter(2, 1e2); // (parameter of normalization)
	AngDisFunc -> FixParameter(3, M_PI); // (additional factor since y(#phi) is in the unit of degree)
	AngDisFunc -> SetParNames("lambda_theta", "lambda_phi", "NormFactor", "pi");
	AngDisFunc -> SetTitle(";cos #theta; #varphi (#circ);Number of generated #varUpsilons");

	/// (histogram for Random sampling from the angular distribution function)
	TH2D *hAngDis = new TH2D("hAngDis", ";cos #theta; #varphi (#circ);Number of generated #varUpsilons", nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax); 

	hAngDis -> FillRandom("AngDisFunc", 1e7); // convert the fuction to histogram for fit procedure


	/// Draw graphs
	TCanvas *cAngDisFunc = new TCanvas("cAngDisFunc", "cAngDisFunc", 700, 600);
	
	/// Set canvas margin
	cAngDisFunc -> SetLeftMargin(.15);
	cAngDisFunc -> SetRightMargin(.22);

	hdummy -> GetZaxis() ->SetRangeUser(0, hAngDis->GetMaximum());
	hdummy -> Draw("COLZ");

	hAngDis -> Draw("COLZ same");

	/// Put text upper right conner of the plot
	// CMS_lumi(FitFunc, "Angular distribution function");

	/// Styles of the texts in the plot
	TLatex* legend = new TLatex();
	legend -> SetTextAlign(22);
	legend -> SetTextSize(0.05);

	/// Put texts inside the plot
	legend -> DrawLatexNDC(.48, .88, "Angular distribution function (ideal)");
	legend -> DrawLatexNDC(.48, .80, Form("#lambda_{#theta} = %.2f, #lambda_{#phi} = %.2f", AngDisFunc->GetParameter(0), AngDisFunc->GetParameter(1)));
	
	gPad -> Update();

	// cAngDisFunc -> SaveAs(Form("Polarization2D/AngularDisributionFuctionTheta%.2f_Phi%.2f.png", AngDisFunc->GetParameter(0), AngDisFunc->GetParameter(1)), "RECREATE");

    ///----------------------------------------------------------------------------------------


	/// Extract Polarization Parameters with 1D Fit (costheta)
	/// Make 2D histogram to 1D (integrate over phi, that is, costheta graph)
	TH1D *hAngDis1DCos = hAngDis -> ProjectionX("cos #theta", 1, nPhiBins);

	/// Define variables
	RooRealVar costheta("costheta", "cos#theta", cosThetaMin, cosThetaMax); 
	RooRealVar phi("phi", "#phi (#circ)", phiMin, phiMax); //y

	RooRealVar normCos("normCos", "normalization factor for costheta graph", hAngDis->GetEntries(), 9.5e6, 1.5e7); //normalization factor
	RooRealVar normPhi("normPhi", "normalization factor for phi graph", hAngDis->GetEntries(), 9.5e6, 1.5e7); //normalization factor
	cout << "--------------------------------------" << endl;
	cout << "entries: " << hAngDis -> GetEntries() << endl;
	cout << "--------------------------------------" << endl;
	RooRealVar lambTheta("lambTheta", "lambda Theta", 1, -1, 1);
	RooRealVar lambPhi("lambPhi", "lambda Phi", 1, -1, 1);
	RooRealVar pi("pi", "pi", M_PI);

	/// Import histogram into RooFit
	RooDataHist cosHist("cosHist","angular distribution with costheta", costheta, hAngDis1DCos);

	/// Define model function
	RooGenericPdf AngDisFuncCosPdf("AngDisFuncCosPdf", "normCos*(1+lambTheta*costheta*costheta)/(3+lambTheta)", RooArgSet(costheta, normCos, lambTheta)); 

	/// Make a frame
	RooPlot* cosframe = costheta.frame(cosThetaMin, cosThetaMax, nCosThetaBins);
	cosHist.plotOn(cosframe, Name("cosHist"), MarkerSize(1.5), DrawOption("P0Z"));

	/// Fit the model to the histogram
	auto* fitResultCos = AngDisFuncCosPdf.fitTo(cosHist, Save(), Extended(kTRUE), Minos(kTRUE), NumCPU(3), Range(cosThetaMin, cosThetaMax));
	fitResultCos -> Print("v");
	/// Put the histogram on the frame
	AngDisFuncCosPdf.plotOn(cosframe, Name("AngDisFuncCosPdf"));
	

	/// Draw the histogram and the fit
	TCanvas* cCosHist = new TCanvas("cCosHist", "cCosHist", 600, 600);
	cCosHist -> cd();
	
	TPad *padCos1 = new TPad("padCos1", "padCos1", 0, 0.25, 1, 1.0);
	padCos1 -> SetBottomMargin(0.03);
	padCos1 -> SetTicks(1, 1);
	padCos1 -> Draw();
	padCos1 -> cd();
	
	// add legends
	cosframe -> addObject(PolarParamsText(weights[chooseIdx][HXorCS][0], weights[chooseIdx][HXorCS][1], lambTheta, lambPhi, false));
	cosframe -> GetXaxis() -> SetLabelOffset(1); // to make it disappear under the pull distribution pad
	cosframe -> SetTitle("");
	cosframe -> Draw();
	gPad -> RedrawAxis();

	// Pull Distribution
	cCosHist -> cd();
	
	TPad* padCos2 = GetPadPullDistribution(cosframe, fitResultCos->floatParsFinal().getSize());

	cCosHist -> cd();
  	padCos1 -> Draw();
  	padCos2 -> Draw();

	cCosHist -> SaveAs(Form("CosTheta1D_lambCos%flambPhi%f.png", weights[chooseIdx][HXorCS][0], weights[chooseIdx][HXorCS][1]), "RECREATE");

	cout << "--------------------------------------" << endl;
	cout << "Done costheta fit!" << endl;
	cout << "normalization: " << normCos << ", lambdaTheta: " << lambTheta << endl;
	cout << "--------------------------------------" << endl;




	/// Fix lambda Theta and use this value for the phi graph fit
	lambTheta.setConstant(kTRUE) ;

	/// Extract Polarization Parameters with 1D Fit (phi)	
	/// Make 2D histogram to 1D (integrate over costheta, that is, phi graph)
	TH1D *hAngDis1DPhi = hAngDis -> ProjectionY("#phi", 1, nCosThetaBins);

	/// Import histogram into RooFit
	RooDataHist phiHist("phiHist","angular distribution with phi", phi, hAngDis1DPhi);

	/// Define model function
	RooGenericPdf AngDisFuncPhiPdf("AngDisFuncPhiPdf", "normPhi*(1+lambPhi*2.*cos(2*phi* pi/180.)/(3+lambTheta))", RooArgSet(phi, normPhi, lambPhi, pi, lambTheta)); 

	/// Make a frame
	RooPlot* phiframe = phi.frame();
	phiHist.plotOn(phiframe, Name("phiHist"), MarkerSize(1.5), DrawOption("P0Z"));

	/// Fit the model to the histogram
	auto* fitResultPhi = AngDisFuncPhiPdf.fitTo(phiHist, Save(), Extended(kTRUE), Minos(kTRUE), NumCPU(3), Range(phiMin, phiMax));
	fitResultPhi -> Print("v");
	/// Put the histogram on the frame
	AngDisFuncPhiPdf.plotOn(phiframe, Name("AngDisFuncPhiPdf"));

	


	/// Draw the histogram and the fit
	TCanvas* cPhiHist = new TCanvas("cPhiHist", "cPhiHist", 600, 600);

	TPad *padPhi1 = new TPad("padPhi1", "padPhi1", 0, 0.25, 1, 1.0);
	padPhi1 -> SetBottomMargin(0.03);
	padPhi1 -> SetTicks(1, 1);
	padPhi1 -> Draw();
	padPhi1 -> cd();
	
	// add legends
	phiframe -> addObject(PolarParamsText(weights[chooseIdx][HXorCS][0], weights[chooseIdx][HXorCS][1], lambTheta, lambPhi, true));
	phiframe -> GetXaxis() -> SetLabelOffset(1); // to make it disappear under the pull distribution pad
	phiframe -> SetTitle("");
	phiframe -> Draw();
	gPad -> RedrawAxis();

	// Pull Distribution
	cCosHist -> cd();
	
	TPad* padPhi2 = GetPadPullDistribution(phiframe, fitResultPhi->floatParsFinal().getSize());

	cPhiHist -> cd();
  	padPhi1 -> Draw();
  	padPhi2 -> Draw();

  	cPhiHist-> SaveAs(Form("Phi1D_lambCos%flambPhi%f.png", weights[chooseIdx][HXorCS][0], weights[chooseIdx][HXorCS][1]), "RECREATE");

	cout << "--------------------------------------" << endl;
	cout << "Done phi fit!" << endl;
	cout << "normalization: " << normPhi << ", lambdaPhi: " << lambPhi << endl;
	cout << "--------------------------------------" << endl;
	
	/// Print out the results
	cout << "input      : lambdaTheta = " << weights[chooseIdx][HXorCS][0] << ", lambdaPhi = " << weights[chooseIdx][HXorCS][1] << endl;
	cout << "fit resulta: lambdaTheta = " << lambTheta << ", lambdaPhi = " << lambPhi << endl;


	/// End measuring time
	end = clock();
	cpu_time = (double)(end - start)/CLOCKS_PER_SEC;
	cout << "Time elapsed: " << cpu_time/60. << "minutes" << endl;


	return;



	}