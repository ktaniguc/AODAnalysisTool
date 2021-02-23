# AODAnalysisTool
## setup and compile
------------------------------------
- default athena version : Athena,22.0.15
- From my experience, AOD files derived from Athena,21.X.XX as input may be available without any change Athena version in setup.sh
```sh
$ git clone https://gitlab.cern.ch/ktaniguc/AODAnalysisTool.git
$ cd AODAnalysisTool
#export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase 
#alias setupATLAS='source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh'
$ source setup.sh
$ ./compile.sh clean #(or "$./compile cmake" if your first compile)
```

## run
**If you want to set up the Athena version of more than 22.0.16(maybe..), then the script for submitting batch job (e.g. "/run/run.sh") does not work with following error**
```txt
coral.Exception: Connection on "ATLASDD" cannot be established ( CORAL : "ConnectionPool::getSessionFromNewConnection" from "CORAL/Services/ConnectionService" )
```
- I prepared new script for walkaround. Please try "$cd run_bsub and $./bsub.sh" instead of "$cd run/ and $./run.sh"
------------------------------------
```sh
$ cd run/ #($cd run_bsub/  in case of more than 22.0.16)
```
- please change the CalcEffAlg_optionsMT.py in run/ (or run_bsub/) directory before you run
- if you want to get the informations of containers, then please make "makeSimpleNtupe" option true. 
- "makeSimpleNtuple": create another ntuple root file of which the name is automaticaly added "..._Ntuple.root"
- "applyMuonVeto": data reduction study with MC ttbar alhad, then this flag should be true to get closer condition with real data
```sh
$ LocalInputFileList="the file list name with PATH as input"
$ ......
$ job += CalcEffAlg( ....
$                    OutputFile = "the file name of output with existing PATH",
$                    ....
$                    TapMethod="NoTagJpsi",     #described at "TapMethod"
$                    doNSWMon=True,             #described at "NSW monitoring"
$                    monTrigName="HLT_mu6_L1MU6"#described at "NSW monitoring"
$                    isAsymNSW=True,            #described at "NSW monitoring"
$                    makeSimpleNtuple=True,     #described at "Branch names (if makeSimpleNtuple=True)"
$                    applyMuonVeto=False,       #skip the events containing muons associated with c, b, tau
$                   )
```
- if you want to analyse the AOD file made by Athena 21.X.XX, please change the following line in "run.sh"
```sh
$ # bsub -q $BSUB_QUE -o ./LOGS/log_${DATE}.txt athena.py CalcEffAlg_optionsMT.py
$ bsub -q $BSUB_QUE -o ./LOGS/log_${DATE}.txt athena.py CalcEffAlg_options.py 
```
```sh
$ ./run.sh # (./bsub.sh    if you setted up more than 22.0.16)
$   # finished running, the file for simple ntuple is made with the name "(named by you)_Ntuple.root"
```

### TapMethod
- the parent particles(Z or J/psi) are selected by this item
- If the case when you want to do only offline-trigger matching without considering trigger acceptance for tag-muon, then please choose "NoTag"
----------------------------------
| method name | description |
|:------------:|:------------:|
| Jpsi | |
| Ztap | | 
| JpsiMC | |
| NoMass | same requirements with the ones of Ztap except for mass requirement. |
| NoTag | no requirements for tan and probe. so only trigger-offline matching are done|
| NoTagJpsi | same setup with NoTag, except for mass requirement. mass requirement is same with the one of Jpsi |  

### NSW monitoring
- if you run AODAnalysisTool with "doNSWMon=True" in jobOption, then "NSW-monitoring.root" file is created automatically.
- NSW-monitoring.root file include some histograms for NSW variables in L2MuonSA
- there are two types for NSW sample which are "Asymmetric:A-side only" and "Symmetric:A and C side". "isAsymNSW" option is prepared for the case if you use Asymmetric sample.
- if "isAsymNSW" option is true, then &eta;<sub>offline probe muon</sub> > 0 and sAddress=-1(Endcap L2SA) is required when NSW-monitoring histograms are filled
- trigger chain associated with trigger objects can be selected by "monTrigName" option.(default is "HLT_mu6_L1MU6")
- trigger efficiencies by p<sub>T</sub> are prepared. the checked chains are automatically setted in initialize of CalcEffAlg.cxx. but &eta;<sub>offline probe muon</sub> are not required even if "isAsymNSW" is setted as true.
----------------------------------
| histogram's name | description |
|:------------:|:------------:|
| m_h_stgcClusterZR | the 2D distribution of stgcClusters|
| m_h_mmClusterZR | the 2D distribution of mmClusters| 
| m_h_residual_offvsinnSP | the shortest distance between inner offline trajectory and position of inner superpoint|
| m_h_residual_roadvsinnSP | the shortest distance between endcap inner road and position of inner superpoint|
| m_h_dtheta_offvsinnSP | &delta;&theta; distribution between inner offline trajectory and inner superpoint|
| m_h_XXpass_(chain) | offline p<sub>T</sub> distribution only when the events passed the trigger step|  
| m_eff_XXpass_(chain) | the efficiencies by offline p<sub>T</sub> distribution |  
| m_h_trigPassEvents_(chain) | the trigger acceptance by each step |

## Branch names (tag-and-probe)
------------------------------------
| branch name | description |
|:------------:|:------------:|
| n_trig | number of trigger chains to measure |
| trigname | names of trigger chains |
| tpsumReqdRl1 | sum of dR threshold between Level1 and offline (tag + probe) |
| tpextdR | offline dR between tag and probe muons( extrapolated to Muon Spectrometer ) |
| invMass | reconstructed mass of tag and probe muons |
| tag_ReqdRL1 | dR threshold between Level1 and offline of tag muon |
| tag_pt | offline p<sub>T</sub> of tag muon |
| tag_eta | offline &eta; of tag muon |
| tag_extEta | offline &eta; of tag muon ( extrapolated to Muon Spectrometer ) |
| tag_phi | offline &phi; of tag muon |
| tag_extPhi | offline &phi; of tag muon ( extrapolated to Muon Spectrometer ) |
| tag_d0 | offline d<sub>0</sub> of tag muon |
| tag_z0 | offline z<sub>0</sub> of tag muon |
| probe_pt | offline p<sub>T</sub> of probe muon |
| probe_eta | offline &eta; of probe muon |
| probe_extEta | offline &eta; of probe muon ( extrapolated to Muon Spectrometer ) |
| probe_phi | offline &phi; of probe muon |
| probe_extPhi | offline &phi; of probe muon ( extrapolated to Muon Spectrometer ) |
| probe_d0 | offline d<sub>0</sub> of probe muon |
| probe_z0 | offline z<sub>0</sub> of probe muon |
| probe_charge | offline charge of probe muon |
| probe_segment... | offline segment. description is [here](https://acode-browser.usatlas.bnl.gov/lxr/source/athena/Event/xAOD/xAODMuon/xAODMuon/versions/MuonSegment_v1.h) |

Branch names beginning with probe_L1, probe_SA, probe_CB, probe_EF are variables of probe muon measured in Level1, muonSA, muComb, EventFilter and respectively.  
For example,

### L1
| branch name | description |
|:------------:|:------------:|
| probe_L1_pass | whether pass L1 or not |
| probe_L1_eta | &eta; L1RoI |
| probe_L1_phi | &phi; L1RoI |
| probe_L1_dR | dR between probe muon and L1RoI |
| probe_L1_thrValue | The highest threshold value (in MeV) passed by the muon candidate |
| probe_L1_roiNum | RoI number |
| probe_L1_thrNumber | the logic number of the highest threshold this RoI passed |
| probe_L1_isMoreCandInRoI | the flag is true when RPC Pad found more than two candidates(if endcap RoI, this value is not correct) |

### SA
| branch name | description |
|:------------:|:------------:|
| probe_SA_pass | whether pass muonSA or not |
| probe_SA_dR | dR between muonSA and offline muons |
| probe_SA_pt | p<sub>T</sub> measured in muonSA |
| probe_SA_eta | &eta; measured in muonSA ( in case of muonSA, etaIP ) |
| probe_SA_phi | &phi; measured in muonSA ( in case of muonSA, etaIP ) |
| probe_SA_etams | &eta; measured in muonSA ( etaMS ) |
| probe_SA_phims | &phi; measured in muonSA ( phiMS ) |
| probe_SA_etabe | &eta; measured in muonSA ( &eta; back-extrapolated to IP with same method of muComb ) |
| probe_SA_phibe | &phi; measured in muonSA ( &phi; back-extrapolated to IP with same method of muComb ) |
| probe_SA_sAddress | the station address of the muon (-1:Endcap, other:Barrel) |
| probe_SA_isRpcFailure | flag to see if RPC is properly read (if false, then falied to fit RPCs) |
| probe_SA_isTgcFailure | flag to see if TGC is properly read |
| probe_SA_superPointR_XX | the measured radious of the muon in one particular super point |
| probe_SA_superPointZ_XX | the measured Z position of the muon in one particular super point |
| probe_SA_superPointSlope_XX | the measured slope of the muon in one particular super point |
| probe_SA_superPointIntercept_XX | the measured intercept of the muon in one particular super point |
| probe_SA_superPointChi2_XX | the chi2 of the fit in one particular super point |
| probe_SA_rpcHitMeasPhi | true = &phi; strip, false = z strip |
| probe_SA_rpcHitLayer | (0,1 = low-pT plane), (2,3 = pivot plane), (4,5 = high-pT plane), (6,7 = feet's outermost plane) |
| probe_SA_mdtHitIsOutlier | 0 = used to make SuperPoint, 1 = not used |
| probe_SA_roadAw | road slope (Z/R), index (0, 1, 2) = (inner, middle, outer) |
| probe_SA_roadBw | road offset (R) |
| probe_SA_zMin | mdt region minimum z, maybe not used in Run3 (Roi based instead) |
| probe_SA_stgcClusterXX | information of stgc clusters|
| probe_SA_mmClusterXX | information of mm clusters|

## Branch names (if makeSimpleNtuple=True)
------------------------------------
| branch name | description |
|:------------:|:------------:|
|  eventNumber | |
|  runNumber||
|  trigname| trigger name for TDT input|
|  n_trig| size of trigger name|
|  isPassedTrig | isPassed the trigger chain. please access like isPassedTrig->at(i_trig)|
|  isPassedL1_evt| isPassed L1 by event. |
|  isPassedSA_evt||
|  isPassedCB_evt||
|  isPassedCBIO_evt|if you analyze l2io chain, then look this info|
|  isPassedEF_evt||
|  isPassedL1|isPassed L1 by L1 objects|
|  isPassedSA||
|  isPassedCB||
|  isPassedCBIO||
|  isPassedEF||
|  tdt_L1RoINumber|roiNumber retrieved by TDT. same structure of index with "isPassedL1".| 
|  tdt_SARoINumber||
|  tdt_CBRoINumber||
|  tdt_SAIORoINumber||
|  tdt_CBIORoINumber||
|  tdt_L1RoISector|roiSector retrieved by TDT. same structure of index with "isPassedL1".| 
|  tdt_SARoISector||
|  tdt_CBRoISector||
|  tdt_SAIORoISector||   
|  tdt_CBIORoISector||
|  tdt_L1_isMoreCandInRoI | the flag is true when RPC Pad found more than two candidates |
|  tdt_SApt|SA pt retrieved by TDT. same structure of index with "isPassedSA"|
|  tdt_SAeta||
|  tdt_SAphi||
|  tdt_SAetaMS||
|  tdt_SAphiMS||
|  tdt_SAsAddress||
|  tdt_SAroiEta||
|  tdt_SAroiPhi||
|  tdt_SAsuperPointAXIS_SECTOR||
|  tdt_SAIOpt|SA IOmode pt retrieved by TDT. same structure of index with "isPassedSAIO"|
|  tdt_SAIOeta||
|  tdt_SAIOphi||
|  tdt_SAIOetaMS||
|  tdt_SAIOphiMS||
|  tdt_SAIOsAddress||
|  tdt_SAIOroiEta||
|  tdt_SAIOroiPhi||
|  tdt_SAIOsuperPointAXIS_SECTOR|superpoint AXIS=R  or Z,  SECTOR=BI etc...|
|  tdt_CBpt||
|  tdt_CBeta||
|  tdt_CBphi||
|  tdt_CBIOpt||
|  tdt_CBIOeta||
|  tdt_CBIOphi||
|  muon_pt|pt calculated in offline|
|  muon_eta||
|  muon_phi||
|  muon_extEta||
|  muon_extPhi||
|  L1_XX|just dump MuonRoIContainer's informations without TDT.  XX="roiEta"  etc...|
|  SA_XX|just dump informations in L2StandAloneMuonContainer "HLT_MuonL2SAInfo"|
|  SAIO_XX|just dump informations in L2StandAloneMuonContainer "HLT_MuonL2SAInfoIOmode". XX="pt" etc....|
|  CB_XX||
|  CBIO_XX||

## How to add the Ntuple contents
------------------------------------
- update soon
- add valiables in CalcEfficiency/??Objects.h, src/TrigMatchingTool.cxx, EventTree
### in SAObjects.h
```sh
$  public :
$    SAObject() :
$    isPassed(false),
$    pt(0), 
$    eta(0),
$    phi(0),
$    roiNum(-1),
$    MyNewValue(initial value)   //add!
$    {};
$    ~SAObject() {};
$  public :
$    bool isPassed;
$    double pt; 
$    double eta;
$    double phi;
$    int roiNum;
$    int MyNewValue;      //add!
```
### in TrigMatchingToolMT ( MT analysis )
```sh
$ bool TrigMatchingToolMT::matchSA( std::string& trig,
$                                   const L1Object l1obj,
$                                   SAObject& saobj )
$ {  .......
$       saobj.pt = (*l2sa)->pt();
$       saobj.eta = (*l2sa)->eta();
$       saobj.phi = (*l2sa)->phi();
$       saobj.MyNewValue = (*l2sa)->MyNewValue();    //add!
```
### in TrigMatchingTool ( Run2 analysis )
- if you want to analyse only MT file, you do not need to change here(then the value of "probe_SA_MyNewValue" in the Ntuple of output is always initial one)
- for setProbes function, matchXX get the TrigCombination variable. On the other hand, for doProbeMatching, matchXX get the FeatureContainer. Therefore, it's ok to add the new value only in matchXX for doProbeMatching
```sh
$ bool TrigMatchingTool::matchSA( const Trig::FeatureContainer& fc,
$                                 std::string& mesSATEName,
$                                 const L1Object l1obj,
$                                 SAObject& saobj )
$ {  .......
$       saobj.pt = l2sa->pt();
$       saobj.eta = l2sa->eta();
$       saobj.phi = l2sa->phi();
$       saobj.MyNewValue = l2sa->MyNewValue();   //add!
```
### in EventTreeMT.h
```sh
$    //
$    std::vector < bool >* probe_SA_pass;
$    std::vector < double >* probe_SA_pt;
$    std::vector < int >* probe_SA_MyNewValue;  //add!
```
### in EventTreeMT.cxx
```sh
$ int EventTreeMT::initialize( TString outfile = "test.root" ) {
$   probe_SA_MyNewValue   = new std::vector < int > ();    //add!
$   //
$   m_tree->Branch( "probe_SA_pass",   &probe_SA_pass );
$   m_tree->Branch( "probe_SA_MyNewValue",   &probe_SA_MyNewValue );   //add!

$ void EventTreeMT::clear() {
$   //
$   probe_SA_pass->clear();
$   probe_SA_pt->clear();
$   probe_SA_MyNewValue->clear();    //add!

$ template<typename TAP> void EventTreeMT::filltree( TAP& tap )
$ {
$ // fill the variable vectors
$   int probe_n = tap.m_vL1objects_probe.size();
$       //
$       probe_SA_pass->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].isPassed );
$       probe_SA_pt->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].pt );
$       probe_SA_MyNewValue->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].MyNewValue );  //add!
```
## Notification
### Tag-and-probe requirements
- I implement the requirements following in [ATLAS trigger support note](https://cds.cern.ch/record/2638544/files/ATL-COM-DAQ-2018-140.pdf)
- Traditional CalcEffTool was not following this requirements but [L2MuonSA's TWiki page](https://twiki.cern.ch/twiki/bin/view/Atlas/L2MuonSA#Cut_conditions_for_efficiency_ca)
- if you want to change traditional requirements, please change m_run3evtSelection's value to false in TagAndProbe(MT).h
- the bit-calculation of L1 RoiSector's reference : [L1TGCNtuple](https://twiki.cern.ch/twiki/bin/view/Main/L1TGCNtuple#sectorAddress_8_bit_information)

## To do
------------------------------------
- enrich the Ntuple contents
- use athena message service // make this available by inheriting AgsTool to MyClass?
