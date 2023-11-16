# UpsilonPolarization_2018PbPb
  
## Data 

**1.**  Skim data (apply cuts) and transform reference frame from OniaTree <br>
  - Files / skimUpsilonCandidates.C
    ```
    - headers: 1) AnalysisParameters.h
               2) ReferenceFrameTransformation/Transformations.h
    - input: OniaTree_miniAOD_PbPbPrompt_112X_DATA_ep.root 
    - output: upsilonSkimmedDataset.root
    ```
    
**2.**  Draw Cos&theta; vs &phi; map   <br>
  - Sandbox / rawCosThetaPhiMap.C
    ```
    - input: Files/upsilonSkimmedDataset.root
    - output: Sandbox/frame_distrib/CosThetaPhiMap_pt%dto%GeV.root
    ```
    
**3.** Acceptance & Efficiency correction <br>
  i) Get acceptance correction map
  - MonteCarlo / acceptanceMap_noGenFilter.C
    ```
    - input: Files/OniaTree_Y%dS_GENONLY_NoFilter.root
    - output: MonteCarlo/acceptanceMaps/%S/AcceptanceResults.root 
    ```
  ii) Get efficiency map
  - MonteCarlo / mapUpsilonEfficiency.C
    ```
    - input: Files/OniaTree_Y%dS_pThat2_HydjetDrumMB_miniAOD.root
    - output: MonteCarlo/EfficiencyMaps/%S/EfficiencyResults_pt%dto%dGeV.root 
    ```
  iii) Do Acceptance & Efficiency correction
  - Files / skimWeightedUpsilonCandidates.C
    ```
    - inputs: 1) Files/OniaTree_miniAOD_PbPbPrompt_112X_DATA_ep.root
          2) MonteCarlo/AcceptanceMaps/1S/AcceptanceResults.root
          3) MonteCarlo/EfficiencyMaps/1S/EfficiencyResults_pt0to30.root
    - output: Files/WeightedUpsilonSkimmedDataset.root 
    ```
  iv) Get Upsilon yield 2D map after Acceptance & Efficiency correction
  - Sandbox / drawCosThetaPhiMapWeightedDataset.C
    ```
    - input: Files/WeightedUpsilonSkimmedDataset.root
    - outputs: 1) Sandbox/frame_distrib/WeightedCosThetaPhiMap_pt%dto%dGeV.root
        2) Sandbox/frame_distrib/WeightedCS_cent%dto%d_pt%dto%dGeV.png
        3) Sandbox/frame_distrib/WeightedHX_cent%dto%d_pt%dto%dGeV.png
    ```
       
**4.**  Extract signal in a given ($p_{T}$, cos&theta;, &phi;) bin  <br>
  i) Get the tail parameters from MC before fit  
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
  ii) Fit data and extract the signal
  - SignalExtraction / nominalFit_hightPt.C <br>
  - SignalExtraction / nominalFit_lowPt.C <br>
    ```
    - inputs: 1) upsilonSkimmedDataset.root
              2) MonteCarlo/SignalParameters.txt 
              (read by the function in Tools/Shortcuts.h)
    - output: fit results
    ```
    
**5.**  Extract polarization parameters


## MC

**1.** Skim MC (apply cuts, weights) and transform reference frame from OniaTree  <br>
       weight = (NColl) x (MC gen weight) x (Data / MC vertex z position) <br>
   - Files / skimRecoUpsilonMC.C
     ```
     - inputs: 1) OniaTree_Y1S_pThat2_HydjetDrumMB_miniAOD.root
               2) OniaTree_Y2S_pThat2_HydjetDrumMB_miniAOD.root
     - output: MCUpsilonSkimmedWeightedDataset.root
     ```



