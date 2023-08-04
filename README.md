# UpsilonPolarization_2018PbPb
	
## Data 

**1.**	Skim data (apply cuts) and transform reference frame from OniaTree <br>
  - Files / skimUpsilonCandidates.C
    ```
    - headers: 1) AnalysisParameters.h
               2) ReferenceFrameTransformation/Transformations.h
    - input: OniaTree_miniAOD_PbPbPrompt_112X_DATA_ep.root 
    - output: upsilonSkimmedDataset.root
    ```
    
**2.** 	Draw Cos&theta; vs &phi; map   <br>
  - Sandbox / rawCosThetaPhiMap.C
    ```
    - input: Files/upsilonSkimmedDataset.root
    - output: Sandbox/frame_distrib/CosThetaPhiMap_pt%dto%GeV.root
    ```
    
**3.** Acceptance & Efficiency correction 
1) Get acceptance correction map
  - MonteCarlo / acceptanceMap_noGenFilter.C
    ```
    - input: Files/OniaTree_Y%dS_GENONLY_NoFilter.root
    - output: MonteCarlo/acceptanceMaps/%S/AcceptanceResults_pt%dto%dGeV.root 
    ```
2) Get efficiency map
3) Do Acceptance & Efficiency correction
   
**4.** 	Extract signal in a given (pT, costheta, pi) bin  <br>
1) Get the tail parameters from MC before fit  
  - MonteCarlo / extractMCSignalTails_symCoreDSCB.C
    ```
    - input: Files/MCUpsilonSkimmedWeightedDataset.root
    - output: MonteCarlo/SignalParameters.txt
    ```
  - MonteCarlo / drawTailParaPtDepend.C
    ```
    - input: MonteCarlo/SignalParameters/symCoreDSCB_cent%dto%d_pt%dto%d.txt
    - output: MonteCarlo/SignalParameters/Alpha(n)_pT_dependence.png
    ``` 
2) Fit data and extract the signal
  - SignalExtraction / nominalFit_hightPt.C <br>
  - SignalExtraction / nominalFit_lowPt.C <br>
    ```
    - input: upsilonSkimmedDataset.root
             MonteCarlo/SignalParameters.txt 
             (read by the function in /Tools/Shortcuts.h)
    - output: fit results
    ```
    
**5.**  Extract polarization parameters


## MC

**1.** Skim MC (apply cuts) and transform reference frame from OniaTree  <br>
   - Files / skimMCUpsilon.C
     ```
     - input: OniaTree_Y1S_pThat2_HydjetDrumMB_miniAOD.root
              OniaTree_Y2S_pThat2_HydjetDrumMB_miniAOD.root
     - output: MCUpsilonSkimmedWeightedDataset.root
     ```



