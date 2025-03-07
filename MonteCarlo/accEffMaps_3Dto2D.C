#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "AccEffHelpers.h"

#include "../Tools/Datasets/RooDataSetHelpers.h"
// #include "../sPlot/SPlotHelpers.h"

#include "../Tools/FitShortcuts.h"
#include "../Tools/Style/Legends.h"

#include "../Polarization/PolarFitHelpers.h"

#include "../Tools/RooFitPDFs/InvariantMassModels.h"
#include "../Tools/Style/FitDistributions.h"

#include "../Tools/RooFitPDFs/CosThetaPolarizationPDF.h"
#include "../Tools/RooFitPDFs/PolarFunc.h"

#include "../ReferenceFrameTransform/Transformations.h"


void accEffMaps_3Dto2D(Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS", const Int_t nCosThetaBins = 5, Double_t cosThetaMin = -0.7, Double_t cosThetaMax = 0.7, const Int_t nPhiBins = 5, Int_t phiMin = -180, Int_t phiMax = 180, Int_t iState = gUpsilonState, Bool_t isPhiFolded = kFALSE, TString accName = "MuonUpsilonTriggerAcc") { // accName = "MuonSimpleAcc", "MuonWithin2018PbPbAcc", or "MuonUpsilonTriggerAcc"
	writeExtraText = true; // if extra text
	extraText = "      Internal";

	using namespace std;
	using namespace RooFit;
	using namespace RooStats;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	/// Bin edges and width
	// Set the bin edges along the cosTheta/phi axis depending on the number of bins, min and max values
	// (If want to use non-uniform bin width, the bin edges should be pre-defined in PolarFitHelpers.h)
	std::vector<Double_t> cosThetaBinEdges = setCosThetaBinEdges(nCosThetaBins, cosThetaMin, cosThetaMax);

	std::vector<Double_t> phiBinEdges = setPhiBinEdges(nPhiBins, phiMin, phiMax);

	// bin width
	Double_t cosThetaStep = (cosThetaBinEdges[nCosThetaBins] - cosThetaBinEdges[0]) / nCosThetaBins;

	Double_t phiStep = (phiBinEdges[nPhiBins] - phiBinEdges[0]) / nPhiBins;

	const char* nominalMapName = NominalTEfficiency3DName(refFrameName);

	// get acceptance maps
	Double_t lambdaTheta = 0., lambdaPhi = 0., lambdaThetaPhi = 0.;
	
	TString MuonAccName = "";
	TString accFileName = "";

	if (accName == TString("MuonUpsilonTriggerAcc")) MuonAccName = "_TriggerAcc";
	else if (accName == TString("MuonWithin2018PbPbAcc")) MuonAccName = "_2018PbPbAcc";
	else if (accName == TString("MuonSimpleAcc")) MuonAccName = "_SimpleAcc";
	else {
		cout << "Invalid acceptance name. Please choose from 'MuonUpsilonTriggerAcc', 'MuonWithin2018PbPbAcc', or 'MuonSimpleAcc'." << endl;
		return;
	}

	if (isPhiFolded == kTRUE) accFileName = Form("./AcceptanceMaps/1S/AcceptanceResults%s.root", MuonAccName.Data());
	else accFileName = Form("./AcceptanceMaps/1S/AcceptanceResults%s_fullPhi.root", MuonAccName.Data());

	TFile* acceptanceFile = openFile(accFileName);
	cout << accFileName << " opened" << endl;
	
	auto* accMap = (TEfficiency*)acceptanceFile->Get(nominalMapName);
	// cout << "Nominal map name: " << nominalMapName << endl;

	if (!accMap) {
    	std::cerr << "Error: accMap is null." << std::endl;
    	return;
	}

	// rebin acceptance maps based on costheta, phi, and pT selection
	TEfficiency* accMapCosThetaPhi = rebinTEff3DMap(accMap, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges);

	// get efficiency maps
	// TFile* efficiencyFile = openFile(Form("./EfficiencyMaps/1S/EfficiencyResults%s.root", gMuonAccName));
	TString effFileName = "";

	if (isPhiFolded == kTRUE) effFileName = Form("./EfficiencyMaps/1S/EfficiencyResults%s.root", MuonAccName.Data());
	else effFileName = Form("./EfficiencyMaps/1S/EfficiencyResults%s_fullPhi.root", MuonAccName.Data());

	TFile* efficiencyFile = openFile(effFileName);
	auto* effMap = (TEfficiency*)efficiencyFile->Get(nominalMapName);
	// cout << "Nominal map name: " << nominalMapName << endl;

	if (!effMap) {
    	std::cerr << "Error: effMap is null." << std::endl;
    	return;
	}

	// rebin efficiency maps based on costheta, phi, and pT selection
	TEfficiency* effMapCosThetaPhi = rebinTEff3DMap(effMap, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges);

	// get relative systematic uncertainty of efficiency
	auto* systEff = (TH3D*)efficiencyFile->Get(SystTEfficiency3DName(refFrameName));

	// rebin uncertainty map based on costheta, phi, and pT selection
	TH2D* systEffCosThetaPhi = rebinRel3DUncMap(effMap, systEff, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges);

	Bool_t isCSframe = (strcmp(refFrameName, "CS") == 0) ? kTRUE : kFALSE;

	// initialize variables used in the loop below
	Double_t totalRelUncHigh = 0, totalRelUncLow = 0;
	Double_t totalUncHigh = 0, totalUncLow = 0;

	Float_t maxYield = 0;

	TCanvas* massCanvas = 0;

	// draw acc and eff histograms to check if the rebinning works well
 	// (last five variables: isAcc, displayValues, displayErrors, extraString, isPhiFolded)
	DrawEfficiency2DHist(accMapCosThetaPhi, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, iState, kTRUE, kTRUE, kTRUE, MuonAccName.Data(), isPhiFolded);

	DrawEfficiency2DHist(effMapCosThetaPhi, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, iState, kFALSE, kTRUE, kTRUE, MuonAccName.Data(), isPhiFolded);
 	
	// (last five variables: isAcc, displayValues, displayErrors, displayYieldValues, extraString, isPhiFolded)
	DrawEffxAcc2DHist(accMapCosThetaPhi, effMapCosThetaPhi, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, iState, kTRUE, kTRUE, kFALSE, MuonAccName.Data(), isPhiFolded);
}

void accEffMaps_3Dto2D_scan(const char* refFrameName = "CS", Bool_t isPhiFolded = kFALSE, TString accName = "MuonUpsilonTriggerAcc") { // accName = "MuonSimpleAcc", "MuonWithin2018PbPbAcc", or "MuonUpsilonTriggerAcc") {

	const int NptBins = 4;
	int ptBinEdges[NptBins + 1] = {0, 2, 6, 12, 20};

	for (int ipt = 0; ipt < NptBins; ipt++) {
		// accEffMaps_3Dto2D(ptBinEdges[ipt], ptBinEdges[ipt + 1], refFrameName, 20, -1, 1, 23, -180, 280, 1, isPhiFolded, accName);
		accEffMaps_3Dto2D(ptBinEdges[ipt], ptBinEdges[ipt + 1], refFrameName, 5, -0.7, 0.7, 7, -180, 240, 1, isPhiFolded, accName);
	}

}