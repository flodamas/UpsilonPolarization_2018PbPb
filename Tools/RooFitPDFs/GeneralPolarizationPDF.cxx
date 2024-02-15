#include "RooFit.h"

#include "Riostream.h"

#include "GeneralPolarizationPDF.h"
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"

#include "math.h"

GeneralPolarizationPDF::GeneralPolarizationPDF(const char* name, const char* title,
                                               RooAbsReal& _cosTheta, RooAbsReal& _phi,
                                               RooAbsReal& _lambdaTheta, RooAbsReal& _lambdaPhi, RooAbsReal& _lambdaThetaPhi) :
  RooAbsPdf(name, title),
  cosTheta("cosTheta", "cos theta variable", this, _cosTheta),
  phi("phi", "phi variable", this, _phi),
  lambdaTheta("lambdaTheta", "lambda theta POI", this, _lambdaTheta),
  lambdaPhi("lambdaPhi", "lambda phi POI", this, _lambdaPhi),
  lambdaThetaPhi("lambdaThetaPhi", "lambda theta-phi POI", this, _lambdaThetaPhi) {}

GeneralPolarizationPDF::GeneralPolarizationPDF(const GeneralPolarizationPDF& other, const char* name) :
  RooAbsPdf(other, name),
  cosTheta("cosTheta", this, other.cosTheta),
  phi("phi", this, other.phi),
  lambdaTheta("lambdaTheta", this, other.lambdaTheta),
  lambdaPhi("lambdaPhi", this, other.lambdaPhi),
  lambdaThetaPhi("lambdaThetaPhi", this, other.lambdaThetaPhi) {}

Double_t GeneralPolarizationPDF::evaluate() const {
	double theta = acos(cosTheta);
	double phiInRadian = phi * M_PI / 180;

	return 1 + lambdaTheta * pow(cosTheta, 2.0) + lambdaPhi * pow(sin(theta), 2.0) * cos(2 * phiInRadian) + lambdaThetaPhi * sin(2 * theta) * cos(phiInRadian);
}
