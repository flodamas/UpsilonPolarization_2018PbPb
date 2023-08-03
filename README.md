# UpsilonPolarization_2018PbPb (writing in progress)
	
	Data 

	1.	Skim data (apply cuts) and transform reference frame from OniaTree 
		Files > skimUpsilonCandidates.C

		headers: 1) AnalysisParameters.h
				 2) ReferenceFrameTransformation > Transformations.h 
		input: OniaTree_miniAOD_PbPbPrompt_112X_DATA_ep.root
		output: upsilonSkimmedDataset.root

	2. 	Draw Costheta vs Phi map  
		Sandbox > drawCosThetaPhiMap.C

	3. 	Extract signal in a given (pT, costheta, pi) bin 
	   	SignalExtraction > nominalFit_hightPt.C
	   	SignalExtraction > nominalFit_lowPt.C

	   	input: upsilonSkimmedDataset.root
	   	output: fit results


	4.  Extract polarization parameters


	MC

	1. Skim MC (apply cuts) and transform reference frame from OniaTree
	   Files > skimMCUpsilon.C


	   input: OniaTree_Y1S_pThat2_HydjetDrumMB_miniAOD.root
	   		  OniaTree_Y2S_pThat2_HydjetDrumMB_miniAOD.root
	   output: MCUpsilonSkimmedWeightedDataset.root



