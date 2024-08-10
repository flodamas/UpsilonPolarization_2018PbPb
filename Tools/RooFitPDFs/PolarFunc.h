#include "../BasicHeaders.h"

TF1* getCosThetaPolarFunc(Float_t maxYield) {
	TF1* cosThetaPolarFunc = new TF1("cosThetaPolarFunc", "[0] / (3 + [1]) * (1 + [1] * x * x)", -1, 1);

	cosThetaPolarFunc->SetParameter(0, 0.5 * maxYield);
	cosThetaPolarFunc->SetParameter(1, 0);

	cosThetaPolarFunc->SetParLimits(1, -2, 2);

	cosThetaPolarFunc->SetParName(0, "normFactor");
	cosThetaPolarFunc->SetParName(1, "lambdaTheta");

	cosThetaPolarFunc->SetLineWidth(3);
	cosThetaPolarFunc->SetLineColor(kRed + 1);

	return cosThetaPolarFunc;
}

TF1* getPhiPolarFunc(Float_t maxYield) {
	TF1* phiPolarFunc = new TF1("phiPolarFunc", "[0] * (1. + 2. * [2] / (3. + [1]) * TMath::Cos(2. * x * pi / 180.))", -180, 180);

	phiPolarFunc->SetParameter(0, 0.5 * maxYield);
	phiPolarFunc->SetParameter(1, 0);
	phiPolarFunc->SetParameter(2, 0);

	phiPolarFunc->SetParLimits(1, -2, 2);
	phiPolarFunc->SetParLimits(2, -2, 2);

	phiPolarFunc->SetParName(0, "normFactor");
	phiPolarFunc->SetParName(1, "lambdaTheta2");
	phiPolarFunc->SetParName(2, "lambdaPhi");

	phiPolarFunc->SetLineWidth(3);
	phiPolarFunc->SetLineColor(kRed + 1);

	return phiPolarFunc;
}

TF1* getPhiTildePolarFunc(Float_t maxYield) {
	TF1* phiTildePolarFunc = new TF1("phiTildePolarFunc", "[0] * (1. + sqrt(2.) * [2] / (3. + [1]) * TMath::Cos(x * pi / 180.))", -180, 180);

	phiTildePolarFunc->SetParameter(0, 0.5 * maxYield);
	phiTildePolarFunc->SetParameter(1, 0);
	phiTildePolarFunc->SetParameter(2, 0);

	phiTildePolarFunc->SetParLimits(1, -2, 2);
	phiTildePolarFunc->SetParLimits(2, -2, 2);

	phiTildePolarFunc->SetParName(0, "normFactor");
	phiTildePolarFunc->SetParName(1, "lambdaTheta3");
	phiTildePolarFunc->SetParName(2, "lambdaPhiTilde");

	phiTildePolarFunc->SetLineWidth(3);
	phiTildePolarFunc->SetLineColor(kRed + 1);

	return phiTildePolarFunc;
}

TF2* getGeneralPolarFunc(Float_t maxYield) {
	// TF2* generalPolarFunc = new TF2("generalPolarFunc", "[0] / (3 + [1]) * (1 + [1] * x * x + [2] * (1 - x * x) * TMath::Cos(2. * y * pi / 180.) + [3] * TMath::Sqrt(1 - x * x) * x * TMath::Cos(y * pi / 180.))", -1, 1, -180, 180);
	TF2* generalPolarFunc = new TF2("generalPolarFunc", "[0] / (3 + [1]) * (1 + [1] * x * x + [2] * TMath::Sin(TMath::ACos(x)) * TMath::Sin(TMath::ACos(x)) * TMath::Cos(2. * y * pi / 180.) + [3] * TMath::Sin(2 * TMath::ACos(x)) * TMath::Cos(y * pi / 180.))", -1, 1, -180, 180);

	generalPolarFunc->SetParameter(0, 0.5 * maxYield);
	generalPolarFunc->SetParameter(1, 0);
	generalPolarFunc->SetParameter(2, 0);
	generalPolarFunc->SetParameter(3, 0);
 
	generalPolarFunc->SetParLimits(1, -2, 2);
	generalPolarFunc->SetParLimits(2, -2, 2);
	generalPolarFunc->SetParLimits(3, -2, 2);

	generalPolarFunc->SetParName(0, "normFactor");
	generalPolarFunc->SetParName(1, "lambdaTheta");
	generalPolarFunc->SetParName(2, "lambdaPhi");
	generalPolarFunc->SetParName(3, "lambdaThetaPhi");

	generalPolarFunc->SetLineWidth(3);
	generalPolarFunc->SetLineColor(kRed + 1);

	return generalPolarFunc;
}

double discontPolarFunc(double *xy, double *par) {
	double x = xy[0];
	double y = xy[1];

	double normFactor = par[0];
	double lambdaTheta = par[1];
	double lambdaPhi = par[2];
	double lambdaThetaPhi = par[3];

	if (fabs(y) < 36) {
		return 0;
	}
	else {
		return normFactor / (3 + lambdaTheta) * (1 + lambdaTheta * x * x + lambdaPhi * (1 - x * x) * TMath::Cos(2. * y * TMath::Pi() / 180.) + lambdaThetaPhi * TMath::Sqrt(1 - x * x) * x * TMath::Cos(y * TMath::Pi() / 180.));
	}
}

TF2* getDiscontPolarFunc(Float_t maxYield) {
	TF2* generalPolarFunc = new TF2("generalPolarFunc", &discontPolarFunc,  -1, 1, -180, 180, 4);

	generalPolarFunc->SetParameter(0, 0.5 * maxYield);
	generalPolarFunc->SetParameter(1, 0);
	generalPolarFunc->SetParameter(2, 0);
	generalPolarFunc->SetParameter(3, 0);
 
	generalPolarFunc->SetParLimits(1, -2, 2);
	generalPolarFunc->SetParLimits(2, -2, 2);
	generalPolarFunc->SetParLimits(3, -2, 2);

	generalPolarFunc->SetParName(0, "normFactor");
	generalPolarFunc->SetParName(1, "lambdaTheta");
	generalPolarFunc->SetParName(2, "lambdaPhi");
	generalPolarFunc->SetParName(3, "lambdaThetaPhi");

	TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 500);

	generalPolarFunc->Draw("COLZ");

	return generalPolarFunc;
}