
#ifndef GENERALPOLARIZATIONPDF
#define GENERALPOLARIZATIONPDF

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"

class GeneralPolarizationPDF : public RooAbsPdf {
 public:
	GeneralPolarizationPDF(){};
	GeneralPolarizationPDF(const char* name, const char* title, RooAbsReal& _cosTheta, RooAbsReal& _phi, RooAbsReal& _lambdaTheta, RooAbsReal& _lambdaPhi, RooAbsReal& _lambdaThetaPhi);

	GeneralPolarizationPDF(const GeneralPolarizationPDF& other, const char* name = 0);
	virtual TObject* clone(const char* newname) const {
		return new GeneralPolarizationPDF(*this, newname);
	}

	inline virtual ~GeneralPolarizationPDF() {}

 protected:
	RooRealProxy cosTheta;
	RooRealProxy phi;
	RooRealProxy lambdaTheta;
	RooRealProxy lambdaPhi;
	RooRealProxy lambdaThetaPhi;

	Double_t evaluate() const;

 private:
	ClassDef(GeneralPolarizationPDF, 1)
};

#endif
