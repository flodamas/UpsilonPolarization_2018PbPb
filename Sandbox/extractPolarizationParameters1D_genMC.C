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

void getPolarizedMCHist(TH2D* generalPolarHist, TH2D* generalPolarTildeHist, Int_t nCosThetaBins, Float_t cosThetaMin, Float_t cosThetaMax, Int_t nPhiBins, Float_t phiMin, Float_t phiMax, Double_t lambdaTheta = 0.88, Double_t lambdaPhi = -0.8, Double_t lambdaThetaPhi = 0){

	// Gen without any filter
	const char* mcFileName = Form("../Files/Y1SGenNoFilterMCDataset_Lambda_Theta%.2f_Phi%.2f_ThetaPhi%.2f.root", lambdaTheta, lambdaPhi, lambdaThetaPhi);

	TFile* mcFile = openFile(mcFileName);

	const char* refFrameName = "HX";

	RooDataSet* allDataset = (RooDataSet*)mcFile->Get(Form("MCdataset%s", refFrameName));

	// import the dataset to a workspace
	RooWorkspace wspace("workspace");
	wspace.import(*allDataset);

	RooRealVar cosTheta = *wspace.var(CosThetaVarName(refFrameName));
	RooRealVar phi = *wspace.var(PhiVarName(refFrameName));
	RooRealVar phiTilde = *wspace.var(PhiTildeVarName(refFrameName));

	/// dummy histogram to adjust the plot range
	TH2F* hdummy = new TH2F("hdummy", ";cos #theta; #varphi (#circ);Number of generated #varUpsilons", nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins + 5, phiMin - 20, phiMax + 100); 

	// TH2D* generalPolarHist = new TH2D("generalPolarHist", ";cos #theta; #varphi (#circ);Number of generated #varUpsilons", nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax); 
	// TH2D* generalPolarTildeHist = new TH2D("generalPolarTildeHist", ";cos #theta; #tilde{#varphi} (#circ);Number of generated #varUpsilons", nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax); 

	for (Int_t iEvent = 0; iEvent < allDataset->numEntries(); iEvent++) {
		const RooArgSet* iRooArgSet = allDataset->get(iEvent);
		
		Double_t cosThetaVal = ((RooRealVar*)iRooArgSet->find(CosThetaVarName(refFrameName)))->getVal();
		Double_t phiVal = ((RooRealVar*)iRooArgSet->find(PhiVarName(refFrameName)))->getVal();
		Double_t phiTildeVal = ((RooRealVar*)iRooArgSet->find(PhiTildeVarName(refFrameName)))->getVal();
		Double_t weight = allDataset->weight();

		generalPolarHist->Fill(cosThetaVal, phiVal, weight);
		generalPolarTildeHist->Fill(cosThetaVal, phiTildeVal, weight);
	}	

	Float_t maxYield = generalPolarHist->GetEntries();

	TF2* generalPolarFuncFit = getGeneralPolarFunc(maxYield);

	TFitResultPtr fitResults = generalPolarHist->Fit("generalPolarFunc", "ESVIM0"); // chi2 fit to the integrated bin 

	// draw plots

	TCanvas* mc2DCanvas = new TCanvas("mc2DCanvas", "mc2DCanvas", 1250, 600);

	mc2DCanvas->Divide(2);

	mc2DCanvas->cd(1);

    gPad->SetRightMargin(0.18);

	hdummy->Draw("COLZ");

	generalPolarHist->Draw("COLZ SAME");

	// Styles of the texts in the plot
	TLatex* legend1 = new TLatex();
	legend1->SetTextAlign(22);
	legend1->SetTextSize(0.05);

	// Put texts inside the plot
	legend1->DrawLatexNDC(.50, .88, "No filter Gen MC");
	legend1->DrawLatexNDC(.50, .80, Form("#lambda_{#theta} = %.2f, #lambda_{#varphi} = %.2f, #lambda_{#theta#varphi} = %.2f", lambdaTheta, lambdaPhi, lambdaThetaPhi));

	mc2DCanvas->cd(2);

	gPad->SetTopMargin(0.05);

	hdummy->Draw("LEGO");

	generalPolarHist->Draw("LEGO SAME");

	generalPolarFuncFit->Draw("SURFACE SAME");

	/// cosmetics	
	mc2DCanvas->SetTopMargin(.15);
	mc2DCanvas->SetLeftMargin(.1);

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

	mc2DCanvas->Update();

	// save the plot
	gSystem->mkdir("DistributionFitsMC", kTRUE);
	mc2DCanvas->SaveAs(Form("DistributionFitsMC/noFilterGenTheta%.2f_Phi%.2f_ThetaPhi%.2f.png", lambdaTheta, lambdaPhi, lambdaThetaPhi), "RECREATE");

	// return generalPolarHist;
}

void extractPolarizationParameters1D(Double_t lambdaTheta0 = 0.88, Double_t lambdaPhi0 = -0.8, Double_t lambdaThetaPhi0 = 0, Double_t n0 = 1e7) {  

	/// Generate a Toy Data (This part can be replaced by data)
	
	// set binning and min, max of cosTheta and phi
	Int_t nCosThetaBins = 20;
	Float_t cosThetaMin = -1, cosThetaMax = 1;

	Int_t nPhiBins = 18;
	Float_t phiMin = -180, phiMax = 180;

	// generate the data
	// (sampled from the general polarization function)
	// TH2D* generalPolarHist = generateGeneralPolarizationHist(nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax, n0, lambdaTheta0, lambdaPhi0, lambdaThetaPhi0);
	
	// (the sample from Gen only MC without any filters)
	TH2D* generalPolarHist = new TH2D("generalPolarHist", ";cos #theta; #varphi (#circ);Number of generated #varUpsilons", nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax); 
	TH2D* generalPolarTildeHist = new TH2D("generalPolarTildeHist", ";cos #theta; #tilde{#varphi} (#circ);Number of generated #varUpsilons", nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax); 

	getPolarizedMCHist(generalPolarHist, generalPolarTildeHist, nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax, lambdaTheta0, lambdaPhi0, lambdaThetaPhi0);

	Double_t nEntries = generalPolarHist->GetEntries();

	cout << "--------------------------------------" << endl;
	cout << "number of entries: " <<  nEntries << endl;
	cout << "--------------------------------------" << endl;

	/// Make 2D histograms to 1D 
	// (integrate over phi, that is, costheta graph)
	TH1D* polarHistCosTheta = generalPolarHist->ProjectionX("cos #theta", 0, nPhiBins); // arguments: (name, firstybin, lastybin)

	// (integrate over cosTheta, that is, phi graph)
	TH1D* polarHistPhi = generalPolarHist->ProjectionY("#varphi", 0, nCosThetaBins);

	TH1D* polarHistPhiTilde = generalPolarTildeHist->ProjectionY("#tilde{#varphi}", 0, nCosThetaBins);

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
	cout << "normalization: " << phifitResults->Parameter(0) << ", lambdaPhi: " << phifitResults->Parameter(2) << endl;
	cout << "--------------------------------------" << endl;
	
	phifitResults->Print("v");

	TF1* phiTildePolarFuncFit = getPhiTildePolarFunc(maxYield);

	// fix the extracted lambdaTheta value
	phiTildePolarFuncFit->FixParameter(1, cosThetafitResults->Parameter(1));

	// perform fit of phi graph 
	TFitResultPtr phiTildefitResults = polarHistPhiTilde->Fit("phiTildePolarFunc", "ESVIMR0"); // chi2 fit to the integrated bin 

	cout << "--------------------------------------" << endl;
	cout << "Done phiTilde fit!" << endl;
	cout << "normalization: " << phiTildefitResults->Parameter(0) << ", lambdaPhiTilde: " << phiTildefitResults->Parameter(2) << endl;
	cout << "--------------------------------------" << endl;
	
	phiTildefitResults->Print("v");

	// draw the histogram and the fit
	TCanvas* polarCanvas1D = new TCanvas("polarCanvas1D", "", 1350, 500);

	polarCanvas1D->Divide(3, 1);

	polarCanvas1D->cd(1);

	TPad* padCosTheta = new TPad("padCosTheta", "padCosTheta", 0, 0.25, 1, 1.0);
	padCosTheta->SetBottomMargin(0.03);
	padCosTheta->SetRightMargin(0.01);
	padCosTheta->SetLeftMargin(0.13);
	padCosTheta->SetTicks(1, 1);
	padCosTheta->Draw();
	padCosTheta->cd();

	polarHistCosTheta->SetMinimum(0);

	polarHistCosTheta->SetMarkerStyle(8);
	polarHistCosTheta->SetMarkerSize(1.5);
	polarHistCosTheta->SetMarkerColor(kBlack);
	
	polarHistCosTheta->GetXaxis()->SetLabelSize(0);

	polarHistCosTheta->GetYaxis()->SetTitle(Form("Events / (%.1f)", (cosThetaMax - cosThetaMin) / nCosThetaBins));
	polarHistCosTheta->GetYaxis()->SetTitleOffset(1.1);

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
	padPhi->SetRightMargin(0.01);
	padPhi->SetLeftMargin(0.13);
	padPhi->SetTicks(1, 1);
	padPhi->Draw();
	padPhi->cd();

	polarHistPhi->SetMinimum(0);

	polarHistPhi->SetMarkerStyle(8);
	polarHistPhi->SetMarkerSize(1.5);
	polarHistPhi->SetMarkerColor(kBlack);
	
	polarHistPhi->GetXaxis()->SetLabelSize(0);

	polarHistPhi->GetYaxis()->SetTitle(Form("Events / (%.1f)", (phiMax - phiMin) / nPhiBins));
	polarHistPhi->GetYaxis()->SetTitleOffset(1.1);

	polarHistPhi->Draw("PE");

	phiPolarFuncFit->Draw("SAME");
  	
  	padPhi->Draw();

  	TPaveStats* phiStatbox = positionStatBox(polarHistPhi);
  	phiStatbox->Draw();

	polarCanvas1D->cd(2);

	// pull Distribution
	TPad* padPhiPull = GetPadPullDistribution(polarHistPhi, phiPolarFuncFit, phifitResults);

	// phi tilde graph
	polarCanvas1D->cd(3);

	TPad* padPhiTilde = new TPad("padPhiTilde", "padPhiTilde", 0, 0.25, 1, 1.0);
	padPhiTilde->SetBottomMargin(0.03);
	padPhiTilde->SetRightMargin(0.01);
	padPhiTilde->SetLeftMargin(0.13);
	padPhiTilde->SetTicks(1, 1);
	padPhiTilde->Draw();
	padPhiTilde->cd();

	polarHistPhiTilde->SetMinimum(0);

	polarHistPhiTilde->SetMarkerStyle(8);
	polarHistPhiTilde->SetMarkerSize(1.5);
	polarHistPhiTilde->SetMarkerColor(kBlack);
	
	polarHistPhiTilde->GetXaxis()->SetLabelSize(0);

	polarHistPhiTilde->GetYaxis()->SetTitle(Form("Events / (%.1f)", (phiMax - phiMin) / nPhiBins));
	polarHistPhiTilde->GetYaxis()->SetTitleOffset(1.1);

	polarHistPhiTilde->Draw("PE");

	phiTildePolarFuncFit->Draw("SAME");
  	
  	padPhiTilde->Draw();

  	TPaveStats* phiTildeStatbox = positionStatBox(polarHistPhiTilde);
  	phiTildeStatbox->Draw();

	polarCanvas1D->cd(3);

	// pull Distribution
	TPad* padPhiTildePull = GetPadPullDistribution(polarHistPhiTilde, phiTildePolarFuncFit, phiTildefitResults);

	gSystem->mkdir("DistributionFitsMC", kTRUE);
  	polarCanvas1D->SaveAs(Form("DistributionFitsMC/fit1D_lambdaCosTheta%.2fPhi%.2fThetaPhi%.2f.png", lambdaTheta0, lambdaPhi0, lambdaThetaPhi0), "RECREATE");

	return;
}