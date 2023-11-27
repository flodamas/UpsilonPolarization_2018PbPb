#include "RooFit.h"

#include "Riostream.h"
#include <math.h>

#include "ErrorFuncTimesExp.h"
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "TMath.h"

ErrorFuncTimesExp::ErrorFuncTimesExp(const char* name, const char* title,
                                     RooAbsReal& _m, RooAbsReal& _mu,
                                     RooAbsReal& _sigma,
                                     RooAbsReal& _lambda) :
  RooAbsPdf(name, title),
  m("m", "mass", this, _m),
  mu("mu", "pole", this, _mu),
  sigma("sigma", "Sigma", this, _sigma),
  lambda("lambda", "exponential rate", this, _lambda) {}

ErrorFuncTimesExp::ErrorFuncTimesExp(const ErrorFuncTimesExp& other, const char* name) :
  RooAbsPdf(other, name),
  m("m", this, other.m),
  mu("mu", this, other.mu),
  sigma("sigma", this, other.sigma),
  lambda("lambda", this, other.lambda) {}

Double_t ErrorFuncTimesExp::evaluate() const {
	return TMath::Exp(-m / lambda) * (TMath::Erf((m - mu) / (TMath::Sqrt(2) * sigma)) + 1) * 0.5;
}
