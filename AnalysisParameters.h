#include "Tools/Style/tdrStyle.C"
#include "Tools/Style/CMS_lumi.C"

/// Definition of all global parameters to be used within this analysis framework

// CANNOT BE CHANGED!
const int gUpsilonHLTBit = 14;
const int gL2FilterBit = 19;
const int gL3FilterBit = 20;

/// Kinematic intervals
const int gCentralityBinMin = 0;
const int gCentralityBinMax = 90;

const float gRapidityMin = 0;
const float gRapidityMax = 2.4;

/// cos theta and phi binning, specific to each pt interval

// select one and comment out the others!!
/*
const int gPtMin = 0;
const int gPtMax = 2;

// 0 < pt < 2 GeV, Collins-Soper
const int NCosThetaBinsCS = 6;
//const double CosThetaBinningCS[NCosThetaBinsCS + 1] = {-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6};
const double CosThetaBinningCS[NCosThetaBinsCS + 1] = {-0.5, -0.3, -0.15, 0, 0.15, 0.3, 0.5};

const int NPhiBinsCS = 8;
const double PhiBinningCS[NPhiBinsCS + 1] = {-180, -120, -80, -40, 0, 40, 80, 120, 180};

// 0 < pt < 2 GeV, Helicity
const int NCosThetaBinsHX = 8;
const double CosThetaBinningHX[NCosThetaBinsHX + 1] = {-1, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 1};

const int NPhiBinsHX = 6;
const double PhiBinningHX[NPhiBinsHX + 1] = {-180, -120, -60, 0, 60, 120, 180};
*/

const int gPtMin = 16;
const int gPtMax = 30;

// 16 < pt < 30 GeV, Collins-Soper
const int NCosThetaBinsCS = 10;
//const double CosThetaBinningCS[NCosThetaBinsCS + 1] = {-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6};
const double CosThetaBinningCS[NCosThetaBinsCS + 1] = {-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1};

const int NPhiBinsCS = 1; // 1D analysis
const double PhiBinningCS[NPhiBinsCS + 1] = {-180, 180};

/// Settings for invariant mass fits

const int NCPUs = 2;

const float MassBinMin = 8;
const float MassBinMax = 13;
const int NMassBins = 80;

// MC signal shape
bool DoMCWeightedError = true;
