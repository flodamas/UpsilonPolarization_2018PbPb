
#ifndef ERRORFUNCTIMESEXP
#define ERRORFUNCTIMESEXP

using namespace RooFit;

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"

class ErrorFuncTimesExp : public RooAbsPdf {
 public:
	ErrorFuncTimesExp(){};
	ErrorFuncTimesExp(const char* name, const char* title, RooAbsReal& _m,
	                  RooAbsReal& _mu, RooAbsReal& _sigma, RooAbsReal& _lambda);

	ErrorFuncTimesExp(const ErrorFuncTimesExp& other, const char* name = 0);
	virtual TObject* clone(const char* newname) const {
		return new ErrorFuncTimesExp(*this, newname);
	}

	inline virtual ~ErrorFuncTimesExp() {}

 protected:
	RooRealProxy m;
	RooRealProxy mu;
	RooRealProxy sigma;
	RooRealProxy lambda;

	Double_t evaluate() const;

 private:
	ClassDef(ErrorFuncTimesExp, 1)
};

#endif
