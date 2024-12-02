
#ifndef ERRORFUNCTIMESEXP
#define ERRORFUNCTIMESEXP

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "TMath.h"
#include "RooMath.h"

class ErrorFuncTimesExp : public RooAbsPdf {
 public:
	ErrorFuncTimesExp(){};
	ErrorFuncTimesExp(const char* name, const char* title, RooAbsReal& _m,
	                  RooAbsReal& _mu, RooAbsReal& _sigma, RooAbsReal& _lambda);

	ErrorFuncTimesExp(const ErrorFuncTimesExp& other, const char* name);
	inline virtual TObject* clone(const char* newname) const {
		return new ErrorFuncTimesExp(*this, newname);
	}

	inline ~ErrorFuncTimesExp() {}

	Double_t evaluate() const;

	ClassDef(ErrorFuncTimesExp, 2)

	  protected : RooRealProxy m;
	RooRealProxy mu;
	RooRealProxy sigma;
	RooRealProxy lambda;
};

#endif
