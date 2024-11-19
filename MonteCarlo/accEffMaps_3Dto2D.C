#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "AccEffHelpers.h"

#include "../Tools/Datasets/RooDataSetHelpers.h"
#include "../sPlot/SPlotHelpers.h"

#include "../Tools/FitShortcuts.h"
#include "../Tools/Style/Legends.h"

#include "../Polarization/PolarFitHelpers.h"

#include "../Tools/RooFitPDFs/InvariantMassModels.h"
#include "../Tools/Style/FitDistributions.h"

#include "../Tools/RooFitPDFs/CosThetaPolarizationPDF.h"
#include "../Tools/RooFitPDFs/PolarFunc.h"

#include "../ReferenceFrameTransform/Transformations.h"


void accEffMaps_3Dto2D(Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS", const Int_t nCosThetaBins = 5, Double_t cosThetaMin = -0.7, Double_t cosThetaMax = 0.7, const Int_t nPhiBins = 5, Int_t phiMin = -180, Int_t phiMax = 180, Int_t iState = gUpsilonState) {
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

	TFile* acceptanceFile = openFile(Form("./AcceptanceMaps/1S/AcceptanceResults%s.root", gMuonAccName));
	cout << Form("./AcceptanceMaps/1S/AcceptanceResults%s.root", gMuonAccName) << endl;
	
	auto* accMap = (TEfficiency*)acceptanceFile->Get(nominalMapName);
	cout << nominalMapName << endl;

	if (!accMap) {
    	std::cerr << "Error: accMap is null." << std::endl;
    	return;
	}

	// rebin acceptance maps based on costheta, phi, and pT selection
	TEfficiency* accMapCosThetaPhi = rebinTEff3DMap(accMap, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges);

	// get efficiency maps
	TFile* efficiencyFile = openFile(Form("./EfficiencyMaps/1S/EfficiencyResults%s.root", gMuonAccName));
	auto* effMap = (TEfficiency*)efficiencyFile->Get(nominalMapName);
	cout << nominalMapName << endl;

	if (!effMap) {
    	std::cerr << "Error: effMap is null." << std::endl;
    	return;
	}

	// rebin efficiency maps based on costheta, phi, and pT selection
	TEfficiency* effMapCosThetaPhi = rebinTEff3DMap(effMap, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges);

	// get relative systematic uncertainty of efficiency
	auto* systEff = (TH3D*)efficiencyFile->Get(RelativeSystTEfficiency3DName(refFrameName));

	// rebin uncertainty map based on costheta, phi, and pT selection
	TH2D* systEffCosThetaPhi = rebinRel3DUncMap(effMap, systEff, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges);

	Bool_t isCSframe = (strcmp(refFrameName, "CS") == 0) ? kTRUE : kFALSE;

	// initialize variables used in the loop below
	Double_t totalRelUncHigh = 0, totalRelUncLow = 0;
	Double_t totalUncHigh = 0, totalUncLow = 0;

	Float_t maxYield = 0;

	TCanvas* massCanvas = 0;

	// draw acc and eff histograms to check if the rebinning works well
 	// (last three variables: isAcc, displayValues, extraString)
	DrawEfficiency2DHist(accMapCosThetaPhi, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, iState, kTRUE, kFALSE, gMuonAccName);

	DrawEfficiency2DHist(effMapCosThetaPhi, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, iState, kFALSE, kFALSE, gMuonAccName);
}

void accEffMaps_3Dto2D_scan(const char* refFrameName = "CS") {

	const int NptBins = 4;
	int ptBinEdges[NptBins + 1] = {0, 2, 6, 12, 20};

	for (int ipt = 0; ipt < NptBins; ipt++) {

		accEffMaps_3Dto2D(ptBinEdges[ipt], ptBinEdges[ipt + 1], refFrameName, 20, -1, 1, 24, -200, 280, 1);

	}

}