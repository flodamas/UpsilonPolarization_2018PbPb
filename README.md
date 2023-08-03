# UpsilonPolarization_2018PbPb (writing in progress)
	
## Data 

**1.**	Skim data (apply cuts) and transform reference frame from OniaTree <br>
  - Files / skimUpsilonCandidates.C
    - headers: 1) AnalysisParameters.h <br>
    $\qquad$ $\qquad$ 2) ReferenceFrameTransformation / Transformations.h <br>
    - input: OniaTree_miniAOD_PbPbPrompt_112X_DATA_ep.root <br>
    - output: upsilonSkimmedDataset.root <br>

**2.** 	Draw Cos&theta; vs &phi; map   <br>
  - Sandbox / rawCosThetaPhiMap.C

**3.** Acceptance & Efficiency correction

**4.** 	Extract signal in a given (pT, costheta, pi) bin  <br>
1) Get the tail parameters before fit from MC 
  - MonteCarlo / extractMCSignalTails_symCoreDSCB.C
    - input: MCUpsilonSkimmedWeightedDataset.root
    - output: MonteCarlo / SignalParameters.txt
2) Fit data and extract the signal
  - SignalExtraction / nominalFit_hightPt.C <br>
  - SignalExtraction / nominalFit_lowPt.C <br>
    - input: upsilonSkimmedDataset.root <br>
    $\qquad$ MonteCarlo / SignalParameters.txt <br>
    $\qquad$ (read by the function in /Tools/Shortcuts.h)
    - output: fit results


**4.**  Extract polarization parameters


## MC

**1.** Skim MC (apply cuts) and transform reference frame from OniaTree  <br>
   - Files / skimMCUpsilon.C
     - input: OniaTree_Y1S_pThat2_HydjetDrumMB_miniAOD.root <br>
     $\qquad$ OniaTree_Y2S_pThat2_HydjetDrumMB_miniAOD.root <br>
     - output: MCUpsilonSkimmedWeightedDataset.root



