#include "../Tools/Style/FitDistributions.h"
#include "../Tools/Style/Legends.h"

using namespace RooFit;

void DrawFitGraph(RooRealVar* massVar, RooWorkspace* wspace, RooDataSet* massDataset, TString pdfName, RooFitResult* fitResult, Float_t massMin, Float_t massMax, Int_t nBins){
	
	/// draw the fit to see if the fit is reasonable
	string spdfName = pdfName.Data(); // convert TString to string
	auto* canvas = new TCanvas(Form("canvas%s",spdfName.c_str()), Form("canvas%s",spdfName.c_str()), 600, 600);
	TPad* pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1.0);
	pad1->SetBottomMargin(0.03);
	pad1->Draw();
	pad1->cd();

	// writeExtraText = true; // if extra text
	writeExtraText = false;
	// extraText = "      Simulation Internal";

	/// define frame for the invariant mass plot with a fit
	RooPlot* frame = massVar -> frame(Title(" "), Range(massMin, massMax));
	frame->GetXaxis()->SetLabelOffset(1); // to make it disappear under the pull distribution pad
	// frame->SetYTitle(Form("Candidates / (%d MeV)", (int)1000 * (binMax - binMin) / nBins));
	massDataset -> plotOn(frame, Name("data"), Binning(nBins), DrawOption("P0Z"));

	RooAbsPdf* signal = nullptr;
	/// get the signal pdf depending on the pdf model
	if(pdfName=="SymDSCB") signal = wspace -> pdf("SymDSCB"); // symmetric core Double-sided Crystal Ball
	else if(pdfName=="SymDSCBGauss") signal = wspace -> pdf("DSCBGauss"); // Double-sided Crystal Ball + Gaussian
	else if(pdfName=="AsymDSCB") signal = wspace -> pdf("AsymDSCB"); // asymmetric core Double-sided Crystal Ball
	else if(pdfName=="Hypatia") signal = wspace -> pdf("Hypatia"); // Hypatia 
	else( "Error! can't find the pdf model..");

	wspace -> Print();	

	/// put the signal pdf on the frame
	signal -> plotOn(frame, LineColor(kBlue));

	/// convert RooRealVals in wspace to Int_t
	RooRealVar* centMinVar = wspace->var("centMinVar");
	Int_t centMin = centMinVar -> getVal();
	RooRealVar* centMaxVar = wspace->var("centMaxVar");
	Int_t centMax = centMaxVar -> getVal();
	
	RooRealVar* ptMinVar = wspace->var("ptMinVar");
	Int_t ptMin = ptMinVar -> getVal();
	RooRealVar* ptMaxVar = wspace->var("ptMaxVar");
	Int_t ptMax = ptMaxVar -> getVal();
	
	RooRealVar* cosThetaMinVar = wspace->var("cosThetaMinVar");
	Int_t cosThetaMin = cosThetaMinVar -> getVal();
	RooRealVar* cosThetaMaxVar = wspace->var("cosThetaMaxVar");
	Int_t cosThetaMax = cosThetaMaxVar -> getVal();

	RooRealVar* isCSframeVar = wspace->var("isCSframeVar");
	Int_t isCSframe = isCSframeVar -> getVal();

	RooRealVar* phiMinVar = wspace->var("phiMinVar");
	Int_t phiMin = phiMinVar -> getVal();
	RooRealVar* phiMaxVar = wspace->var("phiMaxVar");
	Int_t phiMax = phiMaxVar -> getVal();

	/// Get variables in the signal
	RooArgSet* params = nullptr;
	params = signal -> getVariables() ;
	cout << "params of " << pdfName << ": " << *params << endl;

	/// texts on the plot
	frame->addObject(KinematicsText(centMin, centMax, ptMin, ptMax));
	frame->addObject(RefFrameText(isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax));
	if(pdfName=="SymDSCB"){
		RooRealVar* meanVar = (RooRealVar*) params->find("meanSymDSCB");
		RooRealVar* sigmaVar = (RooRealVar*) params->find("sigmaSymDSCB");
		RooRealVar* alphaInfVar = (RooRealVar*) params->find("alphaInfSymDSCB");
		RooRealVar* orderInfVar = (RooRealVar*) params->find("orderInfSymDSCB");
		RooRealVar* alphaSupVar = (RooRealVar*) params->find("alphaSupSymDSCB");
		RooRealVar* orderSupVar = (RooRealVar*) params->find("orderSupSymDSCB");
		frame->addObject(SymCoreDoubleCBParamsText(*meanVar, *sigmaVar, *alphaInfVar, *orderInfVar, *alphaSupVar, *orderSupVar));
	}
	else if(pdfName=="SymDSCBGauss"){
		RooRealVar* meanVar = (RooRealVar*) params->find("meanDSCBGauss");
		RooRealVar* sigmaVar = (RooRealVar*) params->find("sigmaDSCBGauss");
		RooRealVar* alphaInfVar = (RooRealVar*) params->find("alphaInfDSCBGauss");
		RooRealVar* orderInfVar = (RooRealVar*) params->find("orderInfDSCBGauss");
		RooRealVar* alphaSupVar = (RooRealVar*) params->find("alphaSupDSCBGauss");
		RooRealVar* orderSupVar = (RooRealVar*) params->find("orderSupDSCBGauss");
		RooRealVar* sigma_gaussVar = (RooRealVar*) params->find("sigma_gauss");
		RooRealVar* normFractionVar = (RooRealVar*) params->find("normFraction");
		frame->addObject(SymCoreDoubleCBGaussParamsText(*meanVar, *sigmaVar, *alphaInfVar, *orderInfVar, *alphaSupVar, *orderSupVar, *sigma_gaussVar, *normFractionVar));
	}
	else if(pdfName=="AsymDSCB"){
		RooRealVar* meanVar = (RooRealVar*) params->find("meanAsymDSCB");
		RooRealVar* sigmaInfVar = (RooRealVar*) params->find("sigmaInfAsymDSCB");
		RooRealVar* alphaInfVar = (RooRealVar*) params->find("alphaInfAsymDSCB");
		RooRealVar* orderInfVar = (RooRealVar*) params->find("orderInfAsymDSCB");
		RooRealVar* sigmaSupVar = (RooRealVar*) params->find("sigmaSupAsymDSCB");
		RooRealVar* alphaSupVar = (RooRealVar*) params->find("alphaSupAsymDSCB");
		RooRealVar* orderSupVar = (RooRealVar*) params->find("orderSupAsymDSCB");
		frame->addObject(AsymDoubleCBParamsText(*meanVar, *sigmaInfVar, *alphaInfVar, *orderInfVar, *sigmaSupVar, *alphaSupVar, *orderSupVar));	
	}
	else if(pdfName=="Hypatia"){
		RooRealVar* meanVar = (RooRealVar*) params->find("meanHypatia");
		RooRealVar* lambdaVar = (RooRealVar*) params->find("lambdaHypatia");
		RooRealVar* zetaVar = (RooRealVar*) params->find("zetaHypatia");
		RooRealVar* betaVar = (RooRealVar*) params->find("betaHypatia");
		RooRealVar* sigmaVar = (RooRealVar*) params->find("sigmaHypatia");
		RooRealVar* alphaInfVar = (RooRealVar*) params->find("alphaInfHypatia");
		RooRealVar* orderInfVar = (RooRealVar*) params->find("orderInfHypatia");
		RooRealVar* alphaSupVar = (RooRealVar*) params->find("alphaSupHypatia");
		RooRealVar* orderSupVar = (RooRealVar*) params->find("orderSupHypatia");
		frame->addObject(HypatiaParamsText(*meanVar, *lambdaVar, *zetaVar, *betaVar, *sigmaVar, *alphaInfVar, *orderInfVar, *alphaSupVar, *orderSupVar));	
	}
	else("Error! can't find the pdf model..");
	frame->Draw();

	canvas->cd();

	TPad* pad2 = GetPadPullDistribution(frame, fitResult->floatParsFinal().getSize());

	canvas->cd();
	pad1->Draw();
	pad2->Draw();
	canvas->SaveAs(Form("SignalParameters/plots/MCfit_%s_%dto%d.png", spdfName.c_str(), ptMin, ptMax), "RECREATE");
}