#include "ErrorFuncTimesExp.h"

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
	return (RooMath::erf((m - mu) * sigma) + 1) * exp(-m / lambda);
}
