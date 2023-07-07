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

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

#include "../Tools/Style/tdrStyle.C"
#include "../Tools/Style/CMS_lumi.C"

#include "../Tools/Style/FitDistributions.h"
#include "../Tools/Style/Legends.h"

#include "../Tools/Parameters/AnalysisParameters.h"

#include "../ReferenceFrameTransform/Transformations.h"


void acceptanceMap_GenWeights(Int_t iState = 1, Int_t ptMin = 0, Int_t ptMax = 30) {  

	// Start measuring time
	clock_t start, end, cpu_time;
	start = clock();

	// Read GenOnly Nofilter file with HX and CS frames (Run "makeHXCSNtuple.C" to get this file)
	const char* filename = Form("../Files/OniaTree_Y%dS_GENONLY_NoFilter.root", iState);
	TFile* file = TFile::Open(filename, "READ");
	if (!file) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
		return;
	}

	cout << "File " << filename << " opened" << endl;

	TNtuple *HXCSNTuple = (TNtuple*)gDirectory -> Get("HXCSNTuple");

	// Put the text, "CMS Internal", on the right top of the plot 
	// writeExtraText = true; // if extra text
	// extraText = "       Internal";

	writeExtraText = false; // if don't want to put a text


	// Define variables in the Ntuple
	Float_t upsPtLab;
	Float_t upsRapLab;
	Float_t muplPtLab;
	Float_t muplEtaLab;
	Float_t mumiPtLab;
	Float_t mumiEtaLab;
	Float_t muplCosThetaHX;
	Float_t muplPhiHX;
	Float_t muplCosThetaCS;
	Float_t muplPhiCS;

	HXCSNTuple -> SetBranchAddress("upsPtLab", &upsPtLab);
	HXCSNTuple -> SetBranchAddress("upsRapLab", &upsRapLab);
	HXCSNTuple -> SetBranchAddress("muplPtLab", &muplPtLab);
	HXCSNTuple -> SetBranchAddress("muplEtaLab", &muplEtaLab);
	HXCSNTuple -> SetBranchAddress("mumiPtLab", &mumiPtLab);
	HXCSNTuple -> SetBranchAddress("mumiEtaLab", &mumiEtaLab);
	HXCSNTuple -> SetBranchAddress("muplCosThetaHX", &muplCosThetaHX);
	HXCSNTuple -> SetBranchAddress("muplPhiHX", &muplPhiHX);
	HXCSNTuple -> SetBranchAddress("muplCosThetaCS", &muplCosThetaCS);
	HXCSNTuple -> SetBranchAddress("muplPhiCS", &muplPhiCS);

	
	// Generate a Toy Data (This part can be replaced by data)
	// (cos theta, phi) 2D distribution maps for CS and HX frames
	Int_t nCosThetaBins = 20;
	Float_t cosThetaMin = -1, cosThetaMax = 1;

	//	Double_t cosThetaBinning[] = {-1, -0.8, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1};
	//	nCosThetaBins = sizeof(cosThetaBinning) / sizeof(Double_t) - 1;

	Int_t nPhiBins = 23;
	Float_t phiMin = -180, phiMax = 280;

	//Double_t phiBinning[] = {-TMath::Pi(), -2.5, -2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2, 2.5, 180};
	//	nPhiBins = sizeof(cosThetaBinning) / sizeof(Double_t) - 1;

	// Define Muon acceptance cuts
	const char* SingleMuonAccCuts = "(muplPtLab>3.5)&&(mumiPtLab>3.5)&&(abs(muplEtaLab)<2.4)&&(abs(mumiEtaLab)<2.4)";

	char DoubleMuonCuts[50] = {};
	sprintf(DoubleMuonCuts, "(upsPtLab>%d)&&(upsPtLab<%d)&&(abs(upsRapLab)<2.4)", ptMin, ptMax);

	// Set the plot styles
	gStyle -> SetPadLeftMargin(.15);
	gStyle -> SetPadRightMargin(.2);
	//gStyle->SetTitleYOffset(.9);
	gStyle -> SetTitleOffset(1.2,"z"); // gap between color bar and z title
	gStyle -> SetPalette(kCividis);
	gStyle -> SetNumberContours(256);

	// Styles of the texts in the plot
	TLatex* legend = new TLatex();
	legend -> SetTextAlign(22);
	legend -> SetTextSize(0.04);

	// Set the weights
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

	char* weightsLabel[NumIdx] = {"Test1", "Test2", "Start1", "Start2", "Start3", "Start4"};
	cout << weightsLabel[1] << endl;
	// Define histograms for numerator and denominator in acceptance//
	TH2F *hHXNum[NumIdx] = {0};
	TH2F *hHXDen[NumIdx] = {0};
	TH2F *hHXAccp[NumIdx] = {0};
	TCanvas *cHXNum[NumIdx] = {0};
	TCanvas *cHXDen[NumIdx] = {0};
	TCanvas *cHXAccp[NumIdx] = {0};

	TH2F *hCSNum[NumIdx] = {0};
	TH2F *hCSDen[NumIdx] = {0};
	TH2F *hCSAccp[NumIdx] = {0};
	TCanvas *cCSNum[NumIdx] = {0};
	TCanvas *cCSDen[NumIdx] = {0};
	TCanvas *cCSAccp[NumIdx] = {0};

	// Fill histograms
	for(int idx=0; idx<NumIdx; idx++){

		// Define histogram and canvas names

		// (For HX with the first, second order weight (lambda_theta, lambda_phi) (1+ lambda_theta*cos^2(theta)+ lambda_phi*sin^2(theta) * cos(2*phi)))
		char hHXNumName[50]={},hHXDenName[50]={},hHXAccpName[50]={}, cHXNumName[50]={}, cHXDenName[50]={}, cHXAccpName[50]={};
		sprintf(hHXNumName, "hHXNum%s", weightsLabel[idx]);
		sprintf(hHXDenName, "hHXDen%s", weightsLabel[idx]);
		sprintf(hHXAccpName, "hHXAccp%s", weightsLabel[idx]);
		sprintf(cHXNumName, "cHXNum%s", weightsLabel[idx]);
		sprintf(cHXDenName, "cHXDen%s", weightsLabel[idx]);
		sprintf(cHXAccpName, "cHXAccp%s", weightsLabel[idx]);

		// (For CS with the first, second order weight (lambda_theta, lambda_phi) (1+ lambda_theta*cos^2(theta)+ lambda_phi*sin^2(theta) * cos(2*phi)))
		char hCSNumName[50]={},hCSDenName[50]={},hCSAccpName[50]={}, cCSNumName[50]={}, cCSDenName[50]={}, cCSAccpName[50]={};
		sprintf(hCSNumName, "hCSNum%s", weightsLabel[idx]);
		sprintf(hCSDenName, "hCSDen%s", weightsLabel[idx]);
		sprintf(hCSAccpName, "hCSAccp%s", weightsLabel[idx]);
		sprintf(cCSNumName, "cCSNum%s", weightsLabel[idx]);
		sprintf(cCSDenName, "cCSDen%s", weightsLabel[idx]);
		sprintf(cCSAccpName, "cCSAccp%s", weightsLabel[idx]);


		// (For HX with the first, second order weight (lambda_theta, lambda_phi) (1+ lambda_theta*cos^2(theta)+ lambda_phi*sin^2(theta) * cos(2*phi)))
		// (numerator)
		cHXNum[idx] = new TCanvas(cHXNumName, cHXNumName, 700, 600);
		hHXNum[idx] = new TH2F(hHXNumName, ";cos #theta_{HX}; #varphi_{HX} (#circ);Number of generated #varUpsilons", nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);
		
		// (apply weights, cuts and draw 2D graph on the phi-costheta phase space)
		HXCSNTuple -> Draw(Form("muplPhiHX:muplCosThetaHX>>hHXNum%s", weightsLabel[idx]), Form("(1 + %f*pow(muplCosThetaHX, 2.0) + %f*(1- pow(muplCosThetaHX, 2.0))*cos(2*muplPhiHX) )*(%s&&%s)", weights[idx][0][0], weights[idx][0][1], SingleMuonAccCuts, DoubleMuonCuts), "COLZ");
		// HXCSNTuple -> Draw(Form("muplPhiHX:muplCosThetaHX>>hHXNum%s", weightsLabel[idx]), "(1 + 0.88*pow(muplCosThetaHX, 2.0) -0.99*(1- pow(muplCosThetaHX, 2.0))*cos(2*muplPhiHX) )*((muplPtLab>3.5)&&(mumiPtLab>3.5)&&(abs(muplEtaLab)<2.4)&&(abs(mumiEtaLab)<2.4))", "COLZ");
		// Put text upper right conner of the plot
		CMS_lumi(cHXNum[idx], Form("#varUpsilon(%dS) Pythia 8 + EvtGen MC, no filter", iState));

		// Put texts inside the plot
		legend -> DrawLatexNDC(.48, .90, Form("|y^{#mu#mu}| < 2.4, %d < p_{T}^{#mu#mu} < %d GeV", ptMin, ptMax));
		legend -> DrawLatexNDC(.48, .84, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, p_{T}^{#mu} > 3.5 GeV", iState));
		legend -> DrawLatexNDC(.48, .78, Form("#lambda_{#theta, HX} = %.2f, #lambda_{#phi, HX} = %.2f", weights[idx][0][0], weights[idx][0][1]));

		gPad -> Update();

		hHXNum[idx] -> GetYaxis() -> SetRangeUser(-190, 300);
	
		// (Save the canvas)
		cHXNum[idx] -> SaveAs(Form("AcceptanceMaps/%dS/withWeights/HXNum_Theta%.2f_Phi%.2f_pt%dto%dGeV.png", iState, weights[idx][0][0], weights[idx][0][1], ptMin, ptMax), "RECREATE");


		// (dinominator)
		cHXDen[idx] = new TCanvas(cHXDenName, cHXDenName, 700, 600);
		hHXDen[idx] = new TH2F(hHXDenName, ";cos #theta_{HX}; #varphi_{HX} (#circ);Number of generated #varUpsilons", nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);

		HXCSNTuple -> Draw(Form("muplPhiHX:muplCosThetaHX>>hHXDen%s", weightsLabel[idx]), Form("(1+ %f*pow(muplCosThetaHX, 2.0)+ %f*(1- pow(muplCosThetaHX, 2.0)) * cos(2*muplPhiHX))*(%s)", weights[idx][0][0], weights[idx][0][1], DoubleMuonCuts), "COLZ");
		// Put text upper right conner of the plot
		CMS_lumi(cHXDen[idx], Form("#varUpsilon(%dS) Pythia 8 + EvtGen MC, no filter", iState));

		// Put texts inside the plot
		legend -> DrawLatexNDC(.48, .90, Form("|y^{#mu#mu}| < 2.4, %d < p_{T}^{#mu#mu} < %d GeV", ptMin, ptMax));
		// legend -> DrawLatexNDC(.48, .84, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, p_{T}^{#mu} > 3.5 GeV", iState));
		legend -> DrawLatexNDC(.48, .84, Form("#lambda_{#theta, HX} = %.2f, #lambda_{#phi, HX} = %.2f", weights[idx][0][0], weights[idx][0][1]));

		gPad -> Update();

		hHXDen[idx] -> GetYaxis() -> SetRangeUser(-190, 300);

		cHXDen[idx] -> SaveAs(Form("AcceptanceMaps/%dS/withWeights/HXDen_Theta%.2f_Phi%.2f_pt%dto%dGeV.png", iState, weights[idx][0][0], weights[idx][0][1], ptMin, ptMax), "RECREATE");



		// (For CS with the first, second order weight (lambda_theta, lambda_phi) (1+ lambda_theta*cos^2(theta)+ lambda_phi*sin^2(theta) * cos(2*phi)))
		// (numerator)
		cCSNum[idx] = new TCanvas(cCSNumName, cCSNumName, 700, 600);
		hCSNum[idx] = new TH2F(hCSNumName, ";cos #theta_{CS}; #varphi_{CS} (#circ);Number of generated #varUpsilons", nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);
		
		// (apply weights, cuts and draw 2D graph on the phi-costheta phase space)
		HXCSNTuple -> Draw(Form("muplPhiCS:muplCosThetaCS>>hCSNum%s", weightsLabel[idx]), Form("(1 + %f*pow(muplCosThetaCS, 2.0) + %f*(1- pow(muplCosThetaCS, 2.0)) * cos(2*muplPhiCS))*(%s&&%s)", weights[idx][1][0], weights[idx][1][1], SingleMuonAccCuts, DoubleMuonCuts), "COLZ");
		// Put text upper right conner of the plot
		CMS_lumi(cCSNum[idx], Form("#varUpsilon(%dS) Pythia 8 + EvtGen MC, no filter", iState));

		// Put texts inside the plot
		legend -> DrawLatexNDC(.48, .90, Form("|y^{#mu#mu}| < 2.4, %d < p_{T}^{#mu#mu} < %d GeV", ptMin, ptMax));
		legend -> DrawLatexNDC(.48, .84, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, p_{T}^{#mu} > 3.5 GeV", iState));
		legend -> DrawLatexNDC(.48, .78, Form("#lambda_{#theta, CS} = %.2f, #lambda_{#phi, CS} = %.2f", weights[idx][1][0], weights[idx][1][1]));

		gPad -> Update();

		hCSNum[idx] -> GetYaxis() -> SetRangeUser(-190, 300);
	
		// (Save the canvas)
		cCSNum[idx] -> SaveAs(Form("AcceptanceMaps/%dS/withWeights/CSNum_Theta%.2f_Phi%.2f_pt%dto%dGeV.png", iState, weights[idx][1][0], weights[idx][1][1], ptMin, ptMax), "RECREATE");

		// (dinominator)
		cCSDen[idx] = new TCanvas(cCSDenName, cCSDenName, 700, 600);
		hCSDen[idx] = new TH2F(hCSDenName, ";cos #theta_{CS}; #varphi_{CS} (#circ);Number of generated #varUpsilons", nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);

		HXCSNTuple -> Draw(Form("muplPhiCS:muplCosThetaCS>>hCSDen%s", weightsLabel[idx]), Form("(1+ %f*pow(muplCosThetaCS, 2.0)+ %f*(1- pow(muplCosThetaCS, 2.0)) * cos(2*muplPhiCS))*(%s)", weights[idx][1][0], weights[idx][1][1], DoubleMuonCuts), "COLZ");
		// Put text upper right conner of the plot
		CMS_lumi(cCSDen[idx], Form("#varUpsilon(%dS) Pythia 8 + EvtGen MC, no filter", iState));

		// Put texts inside the plot
		legend -> DrawLatexNDC(.48, .90, Form("|y^{#mu#mu}| < 2.4, %d < p_{T}^{#mu#mu} < %d GeV", ptMin, ptMax));
		// legend -> DrawLatexNDC(.48, .84, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, p_{T}^{#mu} > 3.5 GeV", iState));
		legend -> DrawLatexNDC(.48, .84, Form("#lambda_{#theta, CS} = %.2f, #lambda_{#phi, CS} = %.2f", weights[idx][1][0], weights[idx][1][1]));

		gPad -> Update();

		hCSDen[idx] -> GetYaxis() -> SetRangeUser(-190, 300);
	
		cCSDen[idx] -> SaveAs(Form("AcceptanceMaps/%dS/withWeights/CSDen_Theta%.2f_Phi%.2f_pt%dto%dGeV.png", iState, weights[idx][1][0], weights[idx][1][1], ptMin, ptMax), "RECREATE");



		// Calculate the acceptance
		hHXAccp[idx] = new TH2F(hHXAccpName, ";cos #theta_{HX}; #varphi_{HX} (#circ);Acceptance", nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);
		hCSAccp[idx] = new TH2F(hCSAccpName, ";cos #theta_{CS}; #varphi_{CS} (#circ);Acceptance", nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);


		for(int iThetaBin=1; iThetaBin<=nCosThetaBins; iThetaBin++){
			for(int iPhiBin=1; iPhiBin<=nPhiBins; iPhiBin++){
				// (For HX with the first, second order weight (lambda_theta, lambda_phi) (1+ lambda_theta*cos^2(theta)+ lambda_phi*sin^2(theta) * cos(2*phi)))
				if((hHXDen[idx]->GetBinContent(iThetaBin, iPhiBin))==0) continue;
				hHXAccp[idx] -> SetBinContent(iThetaBin, iPhiBin, hHXNum[idx]->GetBinContent(iThetaBin,iPhiBin)/hHXDen[idx]->GetBinContent(iThetaBin,iPhiBin));

				// (For CS with the first, second order weight (lambda_theta, lambda_phi) (1+ lambda_theta*cos^2(theta)+ lambda_phi*sin^2(theta) * cos(2*phi)))
				if((hCSDen[idx]->GetBinContent(iThetaBin,iPhiBin))==0) continue;
				hCSAccp[idx] -> SetBinContent(iThetaBin,iPhiBin, hCSNum[idx]->GetBinContent(iThetaBin,iPhiBin)/hCSDen[idx]->GetBinContent(iThetaBin,iPhiBin));
			}
		}

		cHXAccp[idx] = new TCanvas(cHXAccpName, cHXAccpName, 700, 600);
		hHXAccp[idx] -> Draw("COLZ");
		// Put text upper right conner of the plot
		CMS_lumi(cHXAccp[idx], Form("#varUpsilon(%dS) Pythia 8 + EvtGen MC, no filter", iState));

		// Put texts inside the plot
		legend -> DrawLatexNDC(.48, .90, Form("|y^{#mu#mu}| < 2.4, %d < p_{T}^{#mu#mu} < %d GeV", ptMin, ptMax));
		legend -> DrawLatexNDC(.48, .84, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, p_{T}^{#mu} > 3.5 GeV", iState));
		legend -> DrawLatexNDC(.48, .78, Form("#lambda_{#theta, HX} = %.2f, #lambda_{#phi, HX} = %.2f", weights[idx][0][0], weights[idx][0][1]));

		gPad -> Update();

		hHXAccp[idx] -> GetYaxis() -> SetRangeUser(-190, 300);
		hHXAccp[idx] -> GetZaxis() -> SetRangeUser(0, 1);
		cHXAccp[idx] -> SaveAs(Form("AcceptanceMaps/%dS/withWeights/HXAccp_Theta%.2f_Phi%.2f_pt%dto%dGeV.png", iState, weights[idx][0][0], weights[idx][0][1], ptMin, ptMax), "RECREATE");


		cCSAccp[idx] = new TCanvas(cCSAccpName, cCSAccpName, 700, 600);
		hCSAccp[idx] -> Draw("COLZ");
		// Put text upper right conner of the plot
		CMS_lumi(cCSAccp[idx], Form("#varUpsilon(%dS) Pythia 8 + EvtGen MC, no filter", iState));

		// Put texts inside the plot
		legend -> DrawLatexNDC(.48, .90, Form("|y^{#mu#mu}| < 2.4, %d < p_{T}^{#mu#mu} < %d GeV", ptMin, ptMax));
		legend -> DrawLatexNDC(.48, .84, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, p_{T}^{#mu} > 3.5 GeV", iState));
		legend -> DrawLatexNDC(.48, .78, Form("#lambda_{#theta, CS} = %.2f, #lambda_{#phi, CS} = %.2f", weights[idx][1][0], weights[idx][1][1]));

		gPad -> Update();

		hCSAccp[idx] -> GetYaxis() -> SetRangeUser(-190, 300);
		hCSAccp[idx] -> GetZaxis() -> SetRangeUser(0, 1);
		cCSAccp[idx] -> SaveAs(Form("AcceptanceMaps/%dS/withWeights/CSAccp_Theta%.2f_Phi%.2f_pt%dto%dGeV.png", iState, weights[idx][1][0], weights[idx][1][1], ptMin, ptMax), "RECREATE");


	}

	/// save the results in a file for later usage (in particular for the polarization extraction)
	const char* outputFileName = Form("AcceptanceMaps/%dS/withWeights/AcceptanceResults_withWeights.root", iState);
	TFile outputFile(outputFileName, "RECREATE");

	for(int idx=0; idx<NumIdx; idx++){
		hCSNum[idx] -> Write();
		hHXNum[idx] -> Write();
		hCSAccp[idx] -> Write();
		hHXAccp[idx] -> Write();
	}

	outputFile.Close();

	cout << endl
	     << "Acceptance maps saved in " << outputFileName << endl;

	// End measuring time
	end = clock();
	cpu_time = (double)(end - start)/CLOCKS_PER_SEC;
	cout << "Time elapsed: " << cpu_time/60. << "minutes" << endl;


	return;

}