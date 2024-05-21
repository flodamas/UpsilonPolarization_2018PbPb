/// print out error values of 1D acceptance and efficiency obtained from the skimming codes

void print1DEffErrors(int ptMin = 2, int ptMax = 4, const char* refFrameName = "HX"){

	/// open the acceptance and efficiency result files
	TFile* acceptanceFile = TFile::Open("../MonteCarlo/AcceptanceMaps/1S/AcceptanceResults.root");

	TFile* efficiencyFile = TFile::Open("../MonteCarlo/EfficiencyMaps/1S/EfficiencyResults.root");

	/// set the 1D efficiency plot name
	const char* nominalMapName = Form("CosTheta%s_pt%dto%d", refFrameName, ptMin, ptMax);

	/// load the corresponding plot
	auto* accEffHist = (TEfficiency*)acceptanceFile->Get(nominalMapName);
	auto* effEffHist = (TEfficiency*)efficiencyFile->Get(nominalMapName);

	/// obtain the number of bins of the plots
	int NBins = accEffHist->GetTotalHistogram()->GetNbinsX();

	/// loop over the costheta bins and print error values
	for (int iCosTheta = 1; iCosTheta <= NBins; iCosTheta++){

		cout << "eff error up: " << effEffHist->GetEfficiencyErrorUp(iCosTheta) << endl;
		cout << "eff error low: " << effEffHist->GetEfficiencyErrorLow(iCosTheta) << endl;

		cout << "acc error up: " << accEffHist->GetEfficiencyErrorUp(iCosTheta) << endl;
		cout << "acc error low: " << accEffHist->GetEfficiencyErrorLow(iCosTheta) << endl;

		cout << endl;
	}
	
	return;

}