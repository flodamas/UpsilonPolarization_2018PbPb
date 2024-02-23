#include "../BasicHeaders.h"

TF1* cosThetaPolarFunc(Float_t maxYield){
	TF1* cosThetaPolarFunc = new TF1("PolarFunc", "[1]*(1+[0]*x*x)", -1, 1);
	
	cosThetaPolarFunc->SetParameter(0, 0);
	cosThetaPolarFunc->SetParameter(1, 0.5 * maxYield);
	cosThetaPolarFunc->SetParLimits(0, -1, 1);

	cosThetaPolarFunc->SetParName(0, "lambdaTheta");
	cosThetaPolarFunc->SetParName(1, "normFactor");

	cosThetaPolarFunc->SetLineWidth(3);
	cosThetaPolarFunc->SetLineColor(kRed + 1);

	return cosThetaPolarFunc;
}