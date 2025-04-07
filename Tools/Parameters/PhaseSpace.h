// when the header is included several times, to avoid the redefinition error
#ifndef PHASESPACE_H
#define PHASESPACE_H

Bool_t MuonWithin2018PbPbAcc(const TLorentzVector& muonLV) {
	float eta = fabs(muonLV.Eta());
	float pt = muonLV.Pt();

	return (fabs(eta) < 2.4 &&
	        ((fabs(eta) < 1.2 && pt >= 3.5) ||
	         (1.2 <= fabs(eta) && fabs(eta) < 2.1 && pt >= 5.47 - 1.89 * fabs(eta)) ||
	         (2.1 <= fabs(eta) && pt >= 1.5)));
}

Bool_t MuonSimpleAcc(const TLorentzVector& muonLV) {
	float eta = fabs(muonLV.Eta());
	float pt = muonLV.Pt();
	return (fabs(eta) < 2.4 && pt > 3.5);
}

Bool_t MuonUpsilonTriggerAcc(const TLorentzVector& muonLV) {
	if (!MuonWithin2018PbPbAcc(muonLV)) return false; // minimal condition in any case

	return muonLV.Pt() > 2.5; // threshold applied at L3
}

Bool_t MuonEffStepAcc(const TLorentzVector& muonLV) {
	float eta = fabs(muonLV.Eta());
	float pt = muonLV.Pt();

	return (fabs(eta) < 2.4 &&
	        ((fabs(eta) < 1.3 && pt >= 4.0) ||
	         (1.3 <= fabs(eta) && fabs(eta) < 1.5 && pt >= 3.5) ||
	         (1.5 <= fabs(eta) && pt >= 3.0)));
}

#endif