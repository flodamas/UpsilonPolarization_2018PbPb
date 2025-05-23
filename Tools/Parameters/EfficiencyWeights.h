// Data / MC factors for the z position of the primary vertex (to account for the distribution shift)
// when the header is included several times, to avoid the redefinition error
#ifndef EFFICIENCY_WEIGHTS_H
#define EFFICIENCY_WEIGHTS_H

#include "TF1.h"

const Int_t nBins_zPV = 2 * 40;

const Float_t zPV_weight[nBins_zPV] = {1.28619, 0.687722, 1.03597, 0.880682, 1.11021, 1.50242, 0.995876, 1.18442, 0.971516, 0.939498, 1.19896, 1.02778, 1.18657, 1.13459, 1.2129, 1.12699, 1.16711, 1.18659, 1.18977, 1.24434, 1.15642, 1.17536, 1.20528, 1.16062, 1.16346, 1.18841, 1.14715, 1.14755, 1.13624, 1.133, 1.10685, 1.10358, 1.08736, 1.0798, 1.07459, 1.04775, 1.05677, 1.03215, 1.02656, 1.03072, 1.01536, 0.99726, 1.00044, 0.980783, 0.972381, 0.970434, 0.949301, 0.944801, 0.936345, 0.929792, 0.92977, 0.918105, 0.897718, 0.9107, 0.885819, 0.888816, 0.858479, 0.86341, 0.837018, 0.814051, 0.815967, 0.767701, 0.797429, 0.772583, 0.775842, 0.724849, 0.69916, 0.693069, 0.642584, 0.67468, 0.616945, 0.629109, 0.539328, 0.538292, 0.478972, 0.589674, 0.511554, 0.399289, 0.408257, 0.370139};

Float_t Get_zPV_weight(Float_t zpos) {
	if (fabs(zpos) >= 20) return 1.; // weights not estimated beyond |z| = 20 cm, but it should not be an issue if one discards such MC events

	return zPV_weight[(Int_t)((nBins_zPV / 2) + (2 * zpos))];
}

// Data / MC factors for the reconstructed pT spectra, in two rapidity regions (affordable)

// in form of a fit function (how to deal with parameters uncertainties TBD)

Double_t PtWeight_absy0to1p2(Float_t pT) {
	TF1* PtWeightFunc_absy0to1p2 = new TF1("PtWeightFunc_absy0to1p2", "[0]/([1] + x)", gPtBinning[0], gPtBinning[NPtBins]);
	PtWeightFunc_absy0to1p2->SetParameter(0, 14.1);
	PtWeightFunc_absy0to1p2->SetParError(0, 0.8);
	PtWeightFunc_absy0to1p2->SetParameter(1, 7.6);
	PtWeightFunc_absy0to1p2->SetParError(1, 0.7);

	return PtWeightFunc_absy0to1p2->Eval(pT);
}

Double_t PtWeight_absy1p2to2p4(Float_t pT) {
	TF1* PtWeightFunc_absy1p2to2p4 = new TF1("PtWeightFunc_absy1p2to2p4", "[0]/([1] + x)", gPtBinning[0], gPtBinning[NPtBins]);
	PtWeightFunc_absy1p2to2p4->SetParameter(0, 11.2);
	PtWeightFunc_absy1p2to2p4->SetParError(0, 0.6);
	PtWeightFunc_absy1p2to2p4->SetParameter(1, 5.8);
	PtWeightFunc_absy1p2to2p4->SetParError(1, 0.6);

	return PtWeightFunc_absy1p2to2p4->Eval(pT);
}
/*
Double_t PtWeight_absy0to1p2(Float_t pT) {
	TF1* PtWeightFunc_absy0to1p2 = new TF1("PtWeightFunc_absy0to1p2", "([0] + [1]*x*x) / ( x - [2])^3", gPtBinning[0], gPtBinning[NPtBins]);
	PtWeightFunc_absy0to1p2->SetParameter(0, 23.2);
	PtWeightFunc_absy0to1p2->SetParError(0, 7.9);
	PtWeightFunc_absy0to1p2->SetParameter(1, 13.5);
	PtWeightFunc_absy0to1p2->SetParError(1, 0.7);
	PtWeightFunc_absy0to1p2->SetParameter(2, -1.8);
	PtWeightFunc_absy0to1p2->SetParError(2, 0.2);

	return PtWeightFunc_absy0to1p2->Eval(pT);
}

Double_t PtWeight_absy1p2to2p4(Float_t pT) {
	TF1* PtWeightFunc_absy1p2to2p4 = new TF1("PtWeightFunc_absy1p2to2p4", "([0] + [1]*x*x) / ( x - [2])^3", gPtBinning[0], gPtBinning[NPtBins]);
	PtWeightFunc_absy1p2to2p4->SetParameter(0, 8.4);
	PtWeightFunc_absy1p2to2p4->SetParError(0, 2.8);
	PtWeightFunc_absy1p2to2p4->SetParameter(1, 10.4);
	PtWeightFunc_absy1p2to2p4->SetParError(1, 0.6);
	PtWeightFunc_absy1p2to2p4->SetParameter(2, -1.3);
	PtWeightFunc_absy1p2to2p4->SetParError(2, 0.14);

	return PtWeightFunc_absy1p2to2p4->Eval(pT);
}
*/

Double_t PtGenWeight_absy0to2p4(Float_t pT) {
	TF1* PtWeightFunc_absy0to2p4 = new TF1("PtGenWeightFunc_absy0to2p4", "[0]/([1] + x)", gPtBinning[0], gPtBinning[NPtBins]);
	PtWeightFunc_absy0to2p4->SetParameter(0, 3.7);
	PtWeightFunc_absy0to2p4->SetParError(0, 0.1);
	PtWeightFunc_absy0to2p4->SetParameter(1, 0.4);
	PtWeightFunc_absy0to2p4->SetParError(1, 0.2);

	return PtWeightFunc_absy0to2p4->Eval(pT);
}

Float_t Get_RecoPtWeight(Float_t absY, Float_t pT) {
	if (fabs(absY) < 1.2)
		return PtWeight_absy0to1p2(pT);

	else if (fabs(absY) > 1.2 && fabs(absY) < 2.4)
		return PtWeight_absy1p2to2p4(pT);

	else
		return 0.; // for safety
}

Float_t Get_GenPtWeight(Float_t absY, Float_t pT) {
	if (fabs(absY) < 2.4)
		return PtGenWeight_absy0to2p4(pT);

	else
		return 0.; // for safety
}

#endif // EFFICIENCY_WEIGHTS_H