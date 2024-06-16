#include "../BasicHeaders.h"

TF1* cosThetaPolarFunc(Float_t maxYield) {
	TF1* cosThetaPolarFunc = new TF1("PolarFunc", "[0] / (3 + [1]) * (1 + [1] * x * x)", -1, 1);

	cosThetaPolarFunc->SetParameter(0, 0.5 * maxYield);
	cosThetaPolarFunc->SetParameter(1, 0);

	cosThetaPolarFunc->SetParLimits(1, -2, 2);

	cosThetaPolarFunc->SetParName(0, "normFactor");
	cosThetaPolarFunc->SetParName(1, "lambdaTheta");

	cosThetaPolarFunc->SetLineWidth(3);
	cosThetaPolarFunc->SetLineColor(kRed + 1);

	return cosThetaPolarFunc;
}

TF1* phiPolarFunc(Float_t maxYield) {
	TF1* phiPolarFunc = new TF1("phiPolarFunc", "[0] * (1. + 2. * [2] / (3. + [1]) * TMath::Cos(2. * x * pi / 180.))", -180, 180);

	phiPolarFunc->SetParameter(0, 0.5 * maxYield);
	phiPolarFunc->SetParameter(1, 0);
	phiPolarFunc->SetParameter(2, 0);

	phiPolarFunc->SetParLimits(1, -2, 2);
	phiPolarFunc->SetParLimits(2, -2, 2);

	phiPolarFunc->SetParName(0, "normFactor");
	phiPolarFunc->SetParName(1, "lambdaTheta");
	phiPolarFunc->SetParName(2, "lambdaPhi");

	phiPolarFunc->SetLineWidth(3);
	phiPolarFunc->SetLineColor(kRed + 1);

	return phiPolarFunc;
}

TF2* generalPolarFunc(Float_t maxYield) {
	TF2* polarFunc2D = new TF2("polarFunc2D", "[0] / (3 + [1]) * (1 + [1] * x * x + [2] * (1 - x * x) * TMath::Cos(2. * y * pi / 180.) + [3] * TMath::Sqrt(1 - x * x) * x * TMath::Cos(y * pi / 180.))", -1, 1, -180, 180);

	polarFunc2D->SetParameter(0, 0.5 * maxYield);
	polarFunc2D->SetParameter(1, 0);
	polarFunc2D->SetParameter(2, 0);
	polarFunc2D->SetParameter(3, 0);
 
	polarFunc2D->SetParLimits(1, -2, 2);
	polarFunc2D->SetParLimits(2, -2, 2);
	polarFunc2D->SetParLimits(3, -2, 2);

	polarFunc2D->SetParName(0, "normFactor");
	polarFunc2D->SetParName(1, "lambdaTheta");
	polarFunc2D->SetParName(2, "lambdaPhi");
	polarFunc2D->SetParName(3, "lambdaThetaPhi");

	polarFunc2D->SetLineWidth(3);
	polarFunc2D->SetLineColor(kRed + 1);

	return polarFunc2D;
}