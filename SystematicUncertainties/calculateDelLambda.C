#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/FitShortcuts.h"
#include "../Tools/Style/Legends.h"

void calculateDelLambda(Bool_t isCSframe = true, const char* bkgShapeName = "ExpTimesErr", const Int_t nCosThetaBins = 5, Double_t cosThetaMin = -0.7, Double_t cosThetaMax = 0.7, const Int_t nPhiBins = 3, Int_t phiMin = 0, Int_t phiMax = 180) {
	writeExtraText = true; // if extra text
	extraText = "       Internal";

	RooRealVar* lambdaTheta = new RooRealVar("lambdaTheta", "lambdaTheta", 0);
	RooRealVar* lambdaPhi = new RooRealVar("lambdaPhi", "lambdaPhi", 0);
	RooRealVar* lambdaThetaPhi = new RooRealVar("lambdaThetaPhi", "lambdaThetaPhi", 0);
	RooRealVar* lambdaTilde = new RooRealVar("lambdaTilde", "lambdaTilde", 0);

	RooRealVar* lambdaThetaSys = new RooRealVar("lambdaTheta", "lambdaTheta", 0);
	RooRealVar* lambdaPhiSys = new RooRealVar("lambdaPhi", "lambdaPhi", 0);
	RooRealVar* lambdaThetaPhiSys = new RooRealVar("lambdaThetaPhi", "lambdaThetaPhi", 0);
	RooRealVar* lambdaTildeSys = new RooRealVar("lambdaTilde", "lambdaTilde", 0);

	// TH1D* lambdaThetaHist = new TH1D("lambdaThetaHist", "; p_{T} (GeV/c); #lambda_{#theta}", NPtBins, gPtBinning);
	// TH1D* lambdaPhiHist = new TH1D("lambdaPhiHist", "; p_{T} (GeV/c); #lambda_{#varphi}", NPtBins, gPtBinning);
	// TH1D* lambdaThetaPhiHist = new TH1D("lambdaThetaPhiHist", "; p_{T} (GeV/c); #lambda_{#theta#varphi}", NPtBins, gPtBinning);
	// TH1D* lambdaTildeHist = new TH1D("lambdaTildeHist", "; p_{T} (GeV/c); #tilde{#lambda}", NPtBins, gPtBinning);

	const char* refFrameName = (isCSframe ? "CS" : "HX");
	const char* signalShapeName = "SymDSCB";
	const char* methodName = "rootFit";

    // const char* methodNameSys = "rootFit_SysFitModels_altSig";
    // // const char* methodNameSys = "rootFit_SysFitModels_altBkg";
    // // const char* methodNameSys = "rootFit_SysFitModels_nominal";

    const char* methodNameSys[] = {
        "rootFit_SysFitModels_altSig_nominal",
        "rootFit_SysFitModels_altBkg_nominal",
        "rootFit_SysFitModels_nominal_nominal",
        // "rootFit_SysFitModels__muIDSysDown",
        // "rootFit_SysFitModels__muIDSysUp",
        // "rootFit_SysFitModels__trigSysDown",
        // "rootFit_SysFitModels__trigSysUp",
        // "rootFit_SysFitModels__trkSysDown",
        // "rootFit_SysFitModels__trkSysUp",
        "rootFit_SysFitModels__totalSysUp",
        "rootFit_SysFitModels__totalSysDown",
    };

    double totalUncLambdaTheta[NPtBins] = {0};
    double totalUncLambdaPhi[NPtBins] = {0};
    double totalUncLambdaThetaPhi[NPtBins] = {0};

    for (Int_t index = 0 ; index < 5; index++) {

        cout << "------------------------------------------------------------------------" << endl;
        cout << "** " << refFrameName << " frame" << endl;
        cout << "** " << methodNameSys[index] << endl;
        cout << "| pT (GeV/c) | deltaLambdaTheta | deltaLambdaPhi | deltaLambdaThetaPhi |" << endl;

        for (Int_t ibin = 1; ibin <= NPtBins; ibin++) {
            // get polarization parameters
            const char* fitModelName = GetFitModelName(signalShapeName, gPtBinning[ibin - 1], gPtBinning[ibin], refFrameName, cosThetaMin, cosThetaMax, phiMin, phiMax);

            RooArgSet polarParams = GetPolarParams(lambdaTheta, lambdaPhi, lambdaThetaPhi, lambdaTilde, methodName, fitModelName, bkgShapeName, false);

            double lambdaThetaVal = lambdaTheta->getVal();
            double lambdaThetaUnc = lambdaTheta->getError();

            double lambdaPhiVal = lambdaPhi->getVal();
            double lambdaPhiUnc = lambdaPhi->getError();

            double lambdaThetaPhiVal = lambdaThetaPhi->getVal();
            double lambdaThetaPhiUnc = lambdaThetaPhi->getError();

            double lambdaTildeVal = lambdaTilde->getVal();
            double lambdaTildeUnc = lambdaTilde->getError();


            // get polarization parameters
            // const char* fitModelNameSys = GetFitModelNameSys(signalShapeName, gPtBinning[ibin - 1], gPtBinning[ibin], refFrameName, cosThetaMin, cosThetaMax, phiMin, phiMax);

            RooArgSet polarParamsSys = GetPolarParams(lambdaThetaSys, lambdaPhiSys, lambdaThetaPhiSys, lambdaTildeSys, methodNameSys[index], fitModelName, bkgShapeName, false);

            double lambdaThetaValSys = lambdaThetaSys->getVal();
            double lambdaThetaUncSys = lambdaThetaSys->getError();

            double lambdaPhiValSys = lambdaPhiSys->getVal();
            double lambdaPhiUncSys = lambdaPhiSys->getError();

            double lambdaThetaPhiValSys = lambdaThetaPhiSys->getVal();
            double lambdaThetaPhiUncSys = lambdaThetaPhiSys->getError();

            double lambdaTildeValSys = lambdaTildeSys->getVal();
            double lambdaTildeUncSys = lambdaTildeSys->getError();

            // cout << "pt bin: " << gPtBinning[ibin - 1] << "to" << gPtBinning[ibin] << endl;

            // cout << "lambdaThetaVal: " << lambdaThetaVal << endl;
            // cout << "lambdaThetaUnc: " << lambdaThetaVal << endl;

            // cout << "lambdaPhiVal: " << lambdaPhiVal << endl;
            // cout << "lambdaPhiUnc: " << lambdaPhiUnc << endl;

            // cout << "lambdaThetaPhiVal: " << lambdaThetaPhiVal << endl;
            // cout << "lambdaThetaPhiUnc: " << lambdaThetaPhiUnc << endl;

            // cout << "lambdaTildeVal: " << lambdaTildeVal << endl;
            // cout << "lambdaTildeUnc: " << lambdaTildeUnc << endl;

            // cout << "lambdaThetaValSys: " << lambdaThetaValSys << endl;
            // cout << "lambdaThetaUncSys: " << lambdaThetaValSys << endl;

            // cout << "lambdaPhiValSys: " << lambdaPhiValSys << endl;
            // cout << "lambdaPhiUncSys: " << lambdaPhiUncSys << endl;

            // cout << "lambdaThetaPhiValSys: " << lambdaThetaPhiValSys << endl;
            // cout << "lambdaThetaPhiUncSys: " << lambdaThetaPhiUncSys << endl;

            // cout << "lambdaTildeValSys: " << lambdaTildeValSys << endl;
            // cout << "lambdaTildeUncSys: " << lambdaTildeUncSys << endl;

            // cout << "deltaLambdaTheta: " << lambdaThetaVal - lambdaThetaValSys << endl;
            // cout << "deltaLambdaPhi: " << lambdaPhiVal - lambdaPhiValSys << endl;
            // cout << "deltaLambdaThetaPhi: " << lambdaThetaPhiVal - lambdaThetaPhiValSys << endl;
            // cout << "deltaLambdaTilde: " << lambdaTildeVal - lambdaTildeValSys << endl;

            cout << "|   " << gPtBinning[ibin - 1] << "-" << gPtBinning[ibin] << " & \\num{" << lambdaThetaVal - lambdaThetaValSys << "} & \\num{" << lambdaPhiVal - lambdaPhiValSys << "} & \\num{" << lambdaThetaPhiVal - lambdaThetaPhiValSys << "}      | " << endl;

            // cout << "------------------------------------------------------------------------" << endl;

            if (index != 2) { 
            totalUncLambdaTheta[ibin] += TMath::Power(lambdaThetaVal - lambdaThetaValSys, 2);
            totalUncLambdaPhi[ibin] += TMath::Power(lambdaPhiVal - lambdaPhiValSys, 2);
            totalUncLambdaThetaPhi[ibin] += TMath::Power(lambdaThetaPhiVal - lambdaThetaPhiValSys, 2);
            }
        }

        if (index == 4){
            for (Int_t ibin = 1; ibin <= NPtBins; ibin++) {
                cout << "------------------------------------------------------------------------" << endl;
                cout << "** " << refFrameName << " frame" << endl;
                cout << "| pT (GeV/c) | totalUncLambdaTheta | totalUncLambdaPhi | totalUncLambdaThetaPhi |" << endl;
                cout << "|   " << gPtBinning[ibin - 1] << "-" << gPtBinning[ibin] << " & \\num{" << TMath::Sqrt(totalUncLambdaTheta[ibin]) << "} & \\num{" << TMath::Sqrt(totalUncLambdaPhi[ibin]) << "} & \\num{" << TMath::Sqrt(totalUncLambdaThetaPhi[ibin]) << "}      | " << endl;
                cout << "------------------------------------------------------------------------" << endl;
            }
        }    
    }
}
