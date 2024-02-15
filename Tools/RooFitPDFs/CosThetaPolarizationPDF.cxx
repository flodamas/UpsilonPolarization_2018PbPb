#include "RooFit.h"

#include "Riostream.h"

#include "CosThetaPolarizationPDF.h"
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"

CosThetaPolarizationPDF::CosThetaPolarizationPDF(const char* name, const char* title,
                                                 RooAbsReal& _cosTheta,
                                                 RooAbsReal& _lambdaTheta) :
  RooAbsPdf(name, title),
  cosTheta("cosTheta", "cos theta variable", this, _cosTheta),
  lambdaTheta("lambdaTheta", "lambda theta POI", this, _lambdaTheta) {}

CosThetaPolarizationPDF::CosThetaPolarizationPDF(const CosThetaPolarizationPDF& other, const char* name) :
  RooAbsPdf(other, name),
  cosTheta("cosTheta", this, other.cosTheta),
  lambdaTheta("lambdaTheta", this, other.lambdaTheta) {}

Double_t CosThetaPolarizationPDF::evaluate() const {
	return 1 + lambdaTheta * cosTheta * cosTheta;
}
