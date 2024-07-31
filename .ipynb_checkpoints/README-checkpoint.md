<!-- #region -->
# Belle II Summer Workshop 2023

## 1. Basf2 hands-on session

### Samples

- Some mdst files were downloaded from the grid MC15ri_b for local tests, stored in the `mdst` folder. The local test steering script is `1_Reconstruction_local.py`.

- The final version of the steering script is `1_Reconstruction_grid.py` for generic and misID MC, `1_Reconstruction_grid_signal.py` for signal MC. 2M signal (`B± -> K± pi0`) and 2M misID (`B± -> pi± pi0`), 200/fb generic mixed, charged, ssbar and ccbar samples. The output Ntuples are stored in the `Ntuples` folder.

- Only events in the signal region (`abs(B_deltaE)<0.5`) will eventually be used for the continuum suppression and template fitting. Those events are stored in the `Ntuples/signal_region/`.

### Components

|Components|Selection|
|:---|:---|
|Signal|B_isSignal==1|
|misID|B_mcErrors==128|
|BBbar|B_mcErrors>0 and B_mcErrors!=128|
|qqbar|B_isContinuumEvent==True|

- misID could be selected also by: K_isSignal!=1 and pi_isSignal==1 and K_genMotherPDG==pi_genMotherPDG and abs(K_genMotherPDG)==521

- BBbar could be selected also by: (K_isSignal==1 and B_isSignal!=1) or (K_isSignal!=1 and (K_genMotherPDG!=pi_genMotherPDG or abs(K_genMotherPDG)!=521))


### Selections in the steering script

||Steps|Cuts|
|:---|:---|:---|
|Signal B|1. Tracks|[abs(dz)<2] and [dr<0.5] and thetaInCDCAcceptance and [nCDCHits>20]|
||2. Clusters|[clusterE>0.05] and [abs(clusterTiming)<formula(2*clusterErrorTiming)] and [abs(clusterTiming)<200] <br>and inECLAcceptance and [beamBackgroundSuppression>0.5] and [fakePhotonSuppression>0.8] |
||3. K±| good track and [kaonID_noSVD>0.9]|
||4. Pi0| good cluster daughters and [0.115<M<0.15]|
||5. B±|[Mbc>5.26] and [abs(deltaE)<0.5]|
|-|-|-|
|ROE| 1. Tracks|[abs(dz)<20] and [dr<10] and thetaInCDCAcceptance and [nCDCHits>0] and [pt>0.075]|
||2. Clusters|[clusterE>0.05] and [clusterNHits>1.5] and [abs(clusterTiming)<formula(2*clusterErrorTiming)] and <br>[abs(clusterTiming)<200] and [beamBackgroundSuppression>0.2] and [fakePhotonSuppression>0.2]|
<!-- #endregion -->

```python

```
