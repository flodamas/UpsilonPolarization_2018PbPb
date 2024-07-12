#include "../BasicHeaders.h"

#include "RooFit.h"

#include "Riostream.h"

#include "PhiPolarizationPDF.h"
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"

PhiPolarizationPDF::PhiPolarizationPDF(const char* name, const char* title,
                                                 RooAbsReal& _phi,
                                                 RooAbsReal& _lambdaTheta,
                                                 RooAbsReal& _lambdaPhi) :
  RooAbsPdf(name, title),
  phi("phi", "phi variable", this, _phi),
  lambdaTheta("lambdaTheta", "lambda theta POI", this, _lambdaTheta),
  lambdaPhi("lambdaPhi", "lambda phi POI", this, _lambdaPhi) {}

PhiPolarizationPDF::PhiPolarizationPDF(const PhiPolarizationPDF& other, const char* name) :
  RooAbsPdf(other, name),
  phi("phi", this, other.phi),
  lambdaTheta("lambdaTheta", this, other.lambdaTheta),
  lambdaPhi("lambdaPhi", this, other.lambdaPhi) {}

Double_t PhiPolarizationPDF::evaluate() const {
	return 1. + 2. * lambdaPhi / (3. + lambdaTheta) * std::cos(2. * phi * M_PI / 180);
}
