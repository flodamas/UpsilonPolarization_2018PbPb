#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/FitShortcuts.h"
#include "../Tools/Style/Legends.h"

void calculateDelLambda(const char* bkgShapeName = "ExpTimesErr", const Int_t nCosThetaBins = 5, Double_t cosThetaMin = -0.7, Double_t cosThetaMax = 0.7, const Int_t nPhiBins = 3, Int_t phiMin = 0, Int_t phiMax = 180) {
	writeExtraText = true; // if extra text
	extraText = "       Internal";

    /// define the polarization parameters for nominal and systematic variations
	RooRealVar* lambdaTheta = new RooRealVar("lambdaTheta", "lambdaTheta", 0);
	RooRealVar* lambdaPhi = new RooRealVar("lambdaPhi", "lambdaPhi", 0);
	RooRealVar* lambdaThetaPhi = new RooRealVar("lambdaThetaPhi", "lambdaThetaPhi", 0);
	RooRealVar* lambdaTilde = new RooRealVar("lambdaTilde", "lambdaTilde", 0);

	RooRealVar* lambdaThetaSys = new RooRealVar("lambdaTheta", "lambdaTheta", 0);
	RooRealVar* lambdaPhiSys = new RooRealVar("lambdaPhi", "lambdaPhi", 0);
	RooRealVar* lambdaThetaPhiSys = new RooRealVar("lambdaThetaPhi", "lambdaThetaPhi", 0);
	RooRealVar* lambdaTildeSys = new RooRealVar("lambdaTilde", "lambdaTilde", 0);

    /// some variables for the for loops
	const char* signalShapeName = "SymDSCB";
	const char* methodName = "rootFit";

    const char* refFrameName =  "";

    const int NRefFrame = 2;

    const char* methodNameSys[] = {
        "rootFit_SysFitModels_altSig_nominal",
        "rootFit_SysFitModels_altBkg_nominal",
        "rootFit_SysFitModels_nominal_nominal",
        // "rootFit_SysFitModels_altSig_nominalEff_nominalAcc",
        // "rootFit_SysFitModels_altBkg_nominalEff_nominalAcc",
        // "rootFit_SysFitModels_nominal_nominalEff_nominalAcc",
        
        // "rootFit_SysFitModels__muIDSysDown",
        // "rootFit_SysFitModels__muIDSysUp",
        // "rootFit_SysFitModels__trigSysDown",
        // "rootFit_SysFitModels__trigSysUp",
        // "rootFit_SysFitModels__trkSysDown",
        // "rootFit_SysFitModels__trkSysUp",
        
        // "rootFit_SysFitModels__totalSysUp",
        // "rootFit_SysFitModels__totalSysDown",
        "rootFit_SysFitModels_totalSysUpEff_nominalAcc",
        "rootFit_SysFitModels_totalSysDownEff_nominalAcc",

        "rootFit_SysFitModels_nominalEff_statUpAcc",
        "rootFit_SysFitModels_nominalEff_statDownAcc",
        "rootFit_SysFitModels_statUpEff_nominalAcc",
        "rootFit_SysFitModels_statDownEff_nominalAcc",
    };

    /// arrays to store the systematic uncertainties and calculate the total uncertainty
    double totalUncLambdaTheta[NPtBins][NRefFrame] = {0};
    double totalUncLambdaPhi[NPtBins][NRefFrame] = {0};
    double totalUncLambdaThetaPhi[NPtBins][NRefFrame] = {0};
    double totalUncLambdaTilde[NPtBins][NRefFrame] = {0};

    /// output file for the systematic uncertainty table
    const char* outFileName = "systematic_errors_latex_table.txt";
    std::ofstream outFile(outFileName);

    /// write table headers for the systematic uncertainty table
    outFile << "\\begin{table}[htpb]" << std::endl;
    outFile << "    \\centering" << std::endl;
    outFile << "    \\resizebox{\\textwidth}{!}{" << std::endl;
    outFile << "        \\begin{tabular}{|c|c|c|c|c|c|c|c|}" << std::endl;
    outFile << "        \\hline" << std::endl;
    outFile << "        \\multirow{2}{*}{\\begin{tabular}{@{}c@{}}Uncertainty \\\\ source\\end{tabular}} & " << std::endl;
    outFile << "        \\multirow{2}{*}{\\pt (\\GeVc)} & " << std::endl;
    outFile << "        \\multicolumn{3}{c|}{CS frame} & \\multicolumn{3}{c|}{HX frame} \\\\ " << std::endl;
    outFile << "        \\cline{3-8}" << std::endl;
    outFile << "        && $\\Delta\\lambda_{\\theta}$ & $\\Delta\\lambda_{\\varphi}$ & $\\Delta\\lambda_{\\theta\\varphi}$ & $\\Delta\\lambda_{\\theta}$ & $\\Delta\\lambda_{\\varphi}$ & $\\Delta\\lambda_{\\theta\\varphi}$ \\\\ " << std::endl;
    outFile << "        \\hline" << std::endl;
    outFile << "        \\hline" << std::endl;

    /// loop over the systematic uncertainty types
    for (Int_t index = 0 ; index < 9; index++) {

        /// set the first column of the table for the systematic uncertainty source
        if (index == 0) outFile << "        \\multirow{4}*{\\begin{tabular}[l]{@{}c@{}}Alternative \\\\ signal shape\\end{tabular}}" << std::endl;
        if (index == 1) outFile << "        \\multirow{4}*{\\begin{tabular}[l]{@{}c@{}}Alternative \\\\ background shape\\end{tabular}}" << std::endl;
        if (index == 2) outFile << "        \\multirow{4}*{\\begin{tabular}[l]{@{}c@{}}Background peak \\\\ under signal peak\\end{tabular}}" << std::endl;
        if (index == 3) outFile << "        \\multirow{4}*{\\begin{tabular}[l]{@{}c@{}}Muon scale factor \\\\ systematic up\\end{tabular}}" << std::endl;
        if (index == 4) outFile << "        \\multirow{4}*{\\begin{tabular}[l]{@{}c@{}}Muon scale factor \\\\ systematic down\\end{tabular}}" << std::endl;
        if (index == 5) outFile << "        \\multirow{4}*{\\begin{tabular}[l]{@{}c@{}}Acceptance \\\\ statistical up\\end{tabular}}" << std::endl;
        if (index == 6) outFile << "        \\multirow{4}*{\\begin{tabular}[l]{@{}c@{}}Acceptance \\\\ statistical down\\end{tabular}}" << std::endl;
        if (index == 7) outFile << "        \\multirow{4}*{\\begin{tabular}[l]{@{}c@{}}Efficiency \\\\ statistical up\\end{tabular}}" << std::endl;
        if (index == 8) outFile << "        \\multirow{4}*{\\begin{tabular}[l]{@{}c@{}}Efficiency \\\\ statistical down\\end{tabular}}" << std::endl;

        /// loop over the pt bins
        for (Int_t ibin = 1; ibin < NPtBins; ibin++) {

            /// set the pt bin range (second column of the table)
            outFile << "            & " << gPtBinning[ibin] << "-" << gPtBinning[ibin + 1];

            /// loop over the reference frames
            for (Int_t iRefFrame = 0; iRefFrame < NRefFrame; iRefFrame++) {

                if (iRefFrame == 0) {
                    refFrameName = "CS";
                } else {
                    refFrameName = "HX";
                }

                /// get polarization parameters for the nominal
                const char* fitModelName = GetFitModelName(signalShapeName, gPtBinning[ibin], gPtBinning[ibin + 1], refFrameName, cosThetaMin, cosThetaMax, phiMin, phiMax);

                RooArgSet polarParams = GetPolarParams(lambdaTheta, lambdaPhi, lambdaThetaPhi, lambdaTilde, methodName, fitModelName, bkgShapeName, false);

                double lambdaThetaVal = lambdaTheta->getVal();
                double lambdaThetaUnc = lambdaTheta->getError();

                double lambdaPhiVal = lambdaPhi->getVal();
                double lambdaPhiUnc = lambdaPhi->getError();

                double lambdaThetaPhiVal = lambdaThetaPhi->getVal();
                double lambdaThetaPhiUnc = lambdaThetaPhi->getError();

                double lambdaTildeVal = lambdaTilde->getVal();
                double lambdaTildeUnc = lambdaTilde->getError();

                /// get polarization parameters for the systematic variations
                RooArgSet polarParamsSys = GetPolarParams(lambdaThetaSys, lambdaPhiSys, lambdaThetaPhiSys, lambdaTildeSys, methodNameSys[index], fitModelName, bkgShapeName, false);

                double lambdaThetaValSys = lambdaThetaSys->getVal();
                double lambdaThetaUncSys = lambdaThetaSys->getError();

                double lambdaPhiValSys = lambdaPhiSys->getVal();
                double lambdaPhiUncSys = lambdaPhiSys->getError();

                double lambdaThetaPhiValSys = lambdaThetaPhiSys->getVal();
                double lambdaThetaPhiUncSys = lambdaThetaPhiSys->getError();

                double lambdaTildeValSys = lambdaTildeSys->getVal();
                double lambdaTildeUncSys = lambdaTildeSys->getError();

                /// print the values for cross-check
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

                /// write the values to the output file
                outFile << std::fixed << std::setprecision(8); // Set fixed-point notation with 8 digits after the decimal point 
                outFile << " & \\num{" << lambdaThetaVal - lambdaThetaValSys << "} & \\num{" << lambdaPhiVal - lambdaPhiValSys << "} & \\num{" << lambdaThetaPhiVal - lambdaThetaPhiValSys << "}" ;
                
                /// calculate the total uncertainty for the systematic variations (except for the systematic variations of background peak under signal peak)
                if (index != 2) { 
                totalUncLambdaTheta[ibin][iRefFrame] += TMath::Power(lambdaThetaVal - lambdaThetaValSys, 2);
                totalUncLambdaPhi[ibin][iRefFrame] += TMath::Power(lambdaPhiVal - lambdaPhiValSys, 2);
                totalUncLambdaThetaPhi[ibin][iRefFrame] += TMath::Power(lambdaThetaPhiVal - lambdaThetaPhiValSys, 2);
                totalUncLambdaTilde[ibin][iRefFrame] += TMath::Power(lambdaTildeVal - lambdaTildeValSys, 2);
                }
            }
            outFile << " \\\\ " << std::endl;
            outFile.unsetf(std::ios::fixed); // Unset fixed-point notation
        }

        /// write the total uncertainty for the systematic variations
        if (index == 8){
            
            outFile << "        \\hline" << std::endl;
            outFile << "        \\multirow{4}*{\\begin{tabular}[l]{@{}c@{}}Total systematic \\\\ uncertainty \\end{tabular}}" << std::endl;

            for (Int_t ibin = 1; ibin < NPtBins; ibin++) {
                 outFile << "            & " << gPtBinning[ibin] << "-" << gPtBinning[ibin + 1] ;

                for (Int_t iRefFrame = 0; iRefFrame < NRefFrame; iRefFrame++) {
                    outFile << " & \\num{" << TMath::Sqrt(totalUncLambdaTheta[ibin][iRefFrame]) << "} & \\num{" << TMath::Sqrt(totalUncLambdaPhi[ibin][iRefFrame]) << "} & \\num{" << TMath::Sqrt(totalUncLambdaThetaPhi[ibin][iRefFrame]) << "}";  
                }
                outFile << " \\\\" << std::endl;
            }
        } 
        outFile << "        \\hline" << std::endl; 
    }

    /// write the table footer
    outFile << "        \\end{tabular}" << std::endl;
    outFile << "    }" << std::endl;
    outFile << "    \\caption{Summary of the systematic uncertainties.}" << std::endl;
    outFile << "    \\label{tab:sysUnc}" << std::endl;
    outFile << "\\end{table}" << std::endl;
    
    /// close the file
    outFile.close();

    /// print out the work done
    cout << Form("... Systematic uncertainties in an Overleaf table form written to \"%s\".", outFileName) << endl;
    cout << "... Copy and paste the whole portion of the txt to your Overleaf project!" << endl;  

    /// save the total systematic uncertainties to a txt file
    const char* outFileName_total = "systematic_errors_total.txt";
    std::ofstream outFile_total("systematic_errors_total.txt");

    outFile_total << "lambdaTheta lambdaPhi lambdaThetaPhi lambdaTilde" << endl;

    for (int iRF = 0; iRF < NRefFrame; iRF++) {
        iRF == 0 ? refFrameName = "CS" : refFrameName = "HX";
        for (int ipt = 0; ipt < NPtBins; ipt++) {
            outFile_total << refFrameName << " pT" << gPtBinning[ipt] << "-" << gPtBinning[ipt + 1] <<  " " << TMath::Sqrt(totalUncLambdaTheta[ipt][iRF]) << " " << TMath::Sqrt(totalUncLambdaPhi[ipt][iRF]) << " " << TMath::Sqrt(totalUncLambdaThetaPhi[ipt][iRF]) << " " << TMath::Sqrt(totalUncLambdaTilde[ipt][iRF]) << endl;
        }
    }

    outFile_total.close();

    cout << "" << endl;
    cout << Form("... Total systematic uncertainties also saved in \"%s\".", outFileName_total) << endl; 

    return;
}   
