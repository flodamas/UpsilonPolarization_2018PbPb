# UpsilonPolarization_2018PbPb (writing in progress)
	
## Data 

**1.**	Skim data (apply cuts) and transform reference frame from OniaTree <br>
  - Files > skimUpsilonCandidates.C
    - headers: 1) AnalysisParameters.h <br>
    $\qquad$ $\qquad$ 2) ReferenceFrameTransformation > Transformations.h <br>
    - input: OniaTree_miniAOD_PbPbPrompt_112X_DATA_ep.root <br>
    - output: upsilonSkimmedDataset.root <br>

**2.** 	Draw Costheta vs Phi map   <br>
  - Sandbox > drawCosThetaPhiMap.C

**3.** 	Extract signal in a given (pT, costheta, pi) bin  <br>
  - SignalExtraction > nominalFit_hightPt.C <br>
  - SignalExtraction > nominalFit_lowPt.C <br>
    - input: upsilonSkimmedDataset.root
    - output: fit results


**4.**  Extract polarization parameters


## MC

**1.** Skim MC (apply cuts) and transform reference frame from OniaTree  <br>
   - Files > skimMCUpsilon.C
     - input: OniaTree_Y1S_pThat2_HydjetDrumMB_miniAOD.root <br>
     $\qquad$ OniaTree_Y2S_pThat2_HydjetDrumMB_miniAOD.root <br>
     - output: MCUpsilonSkimmedWeightedDataset.root



