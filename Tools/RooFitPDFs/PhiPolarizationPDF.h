
#ifndef PHIPOLARIZATIONPDF
#define PHIPOLARIZATIONPDF

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"

class PhiPolarizationPDF : public RooAbsPdf {
 public:
	PhiPolarizationPDF(){};
	PhiPolarizationPDF(const char* name, const char* title, RooAbsReal& _phi, RooAbsReal& _lambdaTheta, RooAbsReal& _lambdaPhi);

	PhiPolarizationPDF(const PhiPolarizationPDF& other, const char* name = 0);
	virtual TObject* clone(const char* newname) const {
		return new PhiPolarizationPDF(*this, newname);
	}

	inline virtual ~PhiPolarizationPDF() {}

 protected:
	RooRealProxy phi;
	RooRealProxy lambdaTheta;
	RooRealProxy lambdaPhi;

	Double_t evaluate() const;

 private:
	ClassDef(PhiPolarizationPDF, 1)
};

#endif
