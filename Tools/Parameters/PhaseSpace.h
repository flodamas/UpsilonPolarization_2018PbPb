// when the header is included several times, to avoid the redefinition error
#ifndef PHASESPACE_H
#define PHASESPACE_H

Bool_t PbPb2018MuonLimits(float absEta, float pt) {
	return (absEta < 2.4 &&
	        ((absEta < 1.2 && pt >= 3.5) ||
	         (1.2 <= absEta && absEta < 2.1 && pt >= 5.47 - 1.89 * absEta) ||
	         (2.1 <= absEta && pt >= 1.5)));
}

Bool_t UpsilonTriggerThresholds(float absEta, float pt) {
	if (!PbPb2018MuonLimits(absEta, pt)) return false; // minimal condition in any case

	return pt > 2.5; // threshold applied at L3
}

Bool_t FlatMuonPt3p5Cut(float absEta, float pt) {
	return (absEta < 2.4 && pt > 3.5);
}

Bool_t GoodMuonEffSteps(float absEta, float pt) {
	return (absEta < 2.4 &&
	        ((absEta < 1.4 && pt >= 4.0) ||
	         (1.4 <= absEta && pt >= 3.0)));
}

// single place caller for the muon regions define above
Bool_t MuonKinematicsWithinLimits(const TLorentzVector& muonLV, TString name = gMuonAccName) {
	float absEta = fabs(muonLV.Eta());
	if (absEta > 2.4) return false;

	float pt = muonLV.Pt();

	if (name == "PbPb2018MuonLimits") {
		return PbPb2018MuonLimits(absEta, pt);
	}

	else if (name == "UpsilonTriggerThresholds") {
		return UpsilonTriggerThresholds(absEta, pt);
	}

	else if (name == "FlatMuonPt3p5Cut") {
		return FlatMuonPt3p5Cut(absEta, pt);
	}

	else if (name == "GoodMuonEffSteps") {
		return GoodMuonEffSteps(absEta, pt);
	}

	else {
		cout << "\n[Phase Space] Invalid muon acceptance name. Please choose from 'UpsilonTriggerThresholds', 'PbPb2018MuonLimits', 'FlatMuonPt3p5Cut', or 'GoodMuonEffSteps'." << endl;
		return false;
	}
}

#endif
