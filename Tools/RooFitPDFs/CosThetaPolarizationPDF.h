
#ifndef COSTHETAPOLARIZATIONPDF
#define COSTHETAPOLARIZATIONPDF

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"

class CosThetaPolarizationPDF : public RooAbsPdf {
 public:
	CosThetaPolarizationPDF(){};
	CosThetaPolarizationPDF(const char* name, const char* title, RooAbsReal& _cosTheta, RooAbsReal& _lambdaTheta);

	CosThetaPolarizationPDF(const CosThetaPolarizationPDF& other, const char* name = 0);
	virtual TObject* clone(const char* newname) const {
		return new CosThetaPolarizationPDF(*this, newname);
	}

	inline virtual ~CosThetaPolarizationPDF() {}

 protected:
	RooRealProxy cosTheta;
	RooRealProxy lambdaTheta;

	Double_t evaluate() const;

 private:
	ClassDef(CosThetaPolarizationPDF, 1)
};

#endif
