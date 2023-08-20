#include "../Tools/Style/tdrStyle.C"
#include "../Tools/Style/CMS_lumi.C"

#include "../Tools/Style/FitDistributions.h"
#include "../Tools/Style/Legends.h"

#include "../Tools/FitShortcuts.h"

using namespace RooFit;

Double_t calculateAIC(RooFitResult* fitResult){

	/// calculate the AIC(Akaike Information Criterion) value
	/// (ref: https://en.wikipedia.org/wiki/Akaike_information_criterion)
	Double_t k = fitResult -> floatParsFinal().getSize();/// number of parameters
	Double_t minlogL = fitResult -> minNll(); ///minimum negative log liklihood 
	Double_t AIC = 2.*k + 2.*minlogL;
	
	cout << "k: " << k << endl;
	cout << "-log(L) at minimum: " << minlogL << endl; 
	cout << "AIC :" << AIC << endl;

	return AIC;
}

void AICtest(Int_t centMin = 0, Int_t centMax = 90, Int_t ptMin = 0, Int_t ptMax = 30){

	/// Start measuring time
	clock_t start, end, cpu_time;
	start = clock();


	/// open the skimmed MC file
	const char* filename = "../Files/MCUpsilonSkimmedWeightedDataset.root";

	TFile* file = TFile::Open(filename, "READ");
	if (!file) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
	}
	cout << "File " << filename << " opened" << endl;

	writeExtraText = true; // if extra text
	extraText = "      Simulation Internal";

	// we only extract the tail parameters for a specific pt bin, they do not vary significantly with cos theta or phi within this pt bin
	// Bool_t isCSframe = 1; //or kTRUE
	Bool_t isCSframe = 0; //or kFALSE //need integer expression for Legend

	Float_t cosThetaMin = -1, cosThetaMax = 1;
	Float_t phiMin = -180, phiMax = 180;

	Float_t massMin = 8.5, massMax = 10.5;
	Int_t nBins = 80;

	/// create RooRealVar to contain cuts to save values in a workspace
	RooRealVar centMinVar("centMinVar", "centMinVar", centMin); 
	RooRealVar centMaxVar("centMaxVar", "centMaxVar", centMax); 
	RooRealVar ptMinVar("ptMinVar", "ptMinVar", ptMin); 
	RooRealVar ptMaxVar("ptMaxVar", "ptMaxVar", ptMax); 
	RooRealVar isCSframeVar("isCSframeVar", "isCSframeVar", isCSframe); 
	RooRealVar cosThetaMinVar("cosThetaMinVar", "cosThetaMinVar", cosThetaMin); 
	RooRealVar cosThetaMaxVar("cosThetaMaxVar", "cosThetaMaxVar", cosThetaMax); 
	RooRealVar phiMinVar("phiMinVar", "phiMinVar", phiMin); 
	RooRealVar phiMaxVar("phiMaxVar", "phiMaxVar", phiMax); 

	/// reduce RooDataset to only mass data in the range we are interested in
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	RooDataSet* allDataset = (RooDataSet*)file->Get("MCdataset");

	RooWorkspace* wspace = new RooWorkspace("workspace");
	wspace->import(*allDataset);

	RooDataSet* massDataset = ReducedMassDataset(allDataset, wspace, centMin, centMax, ptMin, ptMax, isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax);
	RooRealVar* massVar = wspace->var("mass");

	wspace -> import(RooArgSet(centMinVar, centMaxVar, ptMinVar, ptMaxVar, isCSframeVar, cosThetaMinVar, cosThetaMaxVar, phiMinVar, phiMaxVar));
	// wspace -> Print();

	/// choose function 1 and function 2 among 
	/// symDSCB, AsymDSCB, symDSCBGauss

	/// The order of pdfs could be swapped
	TString pdf1 = "AsymDSCB";
	TString pdf2 = "SymDSCBGauss";
	TString pdfName[2];
	
	/// perform the fits
	RooFitResult* fitResult[2];
	int idx = 0;
	
	if(pdf1=="SymDSCB" || pdf2 =="SymDSCB"){
		pdfName[idx] = "SymDSCB";
		fitResult[idx] = SymDSCBfit(massVar, wspace, massDataset, massMin, massMax);	
		drawFitGraph(massVar, wspace, massDataset, pdfName[idx], fitResult[idx], massMin, massMax, nBins);
		idx++;
	}
	
	if(pdf1=="AsymDSCB" || pdf2 =="AsymDSCB"){
		pdfName[idx] = "AsymDSCB";
		fitResult[idx] = AsymDSCBfit(massVar, wspace, massDataset, massMin, massMax);	
		drawFitGraph(massVar, wspace, massDataset, pdfName[idx], fitResult[idx], massMin, massMax, nBins);
		idx++;
	}

	if(pdf1=="SymDSCBGauss" || pdf2 =="SymDSCBGauss"){
		pdfName[idx] = "SymDSCBGauss";
		fitResult[idx] = SymDSCBGaussfit(massVar, wspace, massDataset, massMin, massMax);	
		drawFitGraph(massVar, wspace, massDataset, pdfName[idx], fitResult[idx], massMin, massMax, nBins);
		idx++;
	}

	// fitResult[0] -> Print("v");
	// fitResult[1] -> Print("v");

	wspace -> Print();
	// // /// calculate AIC values
	// // Double_t AIC_SymDSCB = calculateAIC(fitResult_SymDSCB);
	// // // Double_t AIC_AsymDSCB = calculateAIC(fitResult_AsymDSCB);
 	// // Double_t AIC_SymDSCBGauss = calculateAIC(fitResult_SymDSCBGauss);

	/// get minNLL values
	Double_t minNLL_f1 = fitResult[0] -> minNll();
	Double_t minNLL_f2 = fitResult[1] -> minNll();

 	// // cout << AIC_SymDSCB << endl;
 	// // cout << AIC_AsymDSCB << endl;

	/// check check :) 
 	cout << "minNll of " << pdfName[0] << ": " << minNLL_f1 << endl;
 	cout << "minNll of " << pdfName[1] << ": " << minNLL_f2 << endl;

 	Double_t Chi2Score = 2.*(minNLL_f1 - minNLL_f2);
 	Int_t ndf = (fitResult[0]->floatParsFinal().getSize()) - (fitResult[1]->floatParsFinal().getSize());

 	cout << "chi2score: " << Chi2Score << endl;
 	cout << "ndf: " << abs(ndf) << endl;

	/// End measuring time
	end = clock();
	cpu_time = (double)(end - start) / CLOCKS_PER_SEC;
	cout << "Time elapsed: " << cpu_time / 60. << "minutes" << endl;

}