#ifndef _TrigMuonAnalysis_EventTreeMT
#define _TrigMuonAnalysis_EventTreeMT

#include <iostream>
#include <fstream>
#include <vector>
#include <stdint.h>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TLorentzVector.h"

#include "CalcEfficiency/TagAndProbeMT.h"
#include "CalcEfficiency/TagAndProbe.h"
#include "CalcEfficiency/Utils.h"

class EventTreeMT {
	
  public:
    EventTreeMT();
    ~EventTreeMT();

    TFile* m_file; //!
    TTree* m_tree; //!
  
  public:
    int initialize( TString outfile );
    void clear();
    template<typename TAP> void filltree( TAP& tap, 
                                          int evtNum, 
                                          int runNum, 
                                          int lumiBlk, 
                                          float avgIntPerXing);
    int finalize();
    
    //event info
    int eventNumber;
    int runNumber;
    int lumiBlock;
    float averageInteractionsPerCrossing;
    //
    double invMass;
    std::vector<std::string>* trigname;
    std::vector<int>* passtrig;
    std::vector<std::string>* L1trigname;
    std::vector < bool >* L1_passbyEvent;
    int n_trig;
    double tpsumReqdRl1;
    double tpextdR;
    //
    double tag_pt;
    double tag_eta;
    double tag_phi;
    double tag_extEta;
    double tag_extPhi;
    double tag_d0;
    double tag_z0;
    std::vector < int >* tag_L1_pass;
    std::vector < int >* tag_L1_roiNum;
    std::vector < int >* tag_SA_pass;
    std::vector < int >* tag_CB_pass;
    std::vector < int >* tag_EF_pass;
    //
    double probe_pt;
    double probe_eta;
    double probe_phi;
    double probe_extEta;
    double probe_extPhi;
    double probe_d0;
    double probe_z0;
    double probe_segment_n;
    double probe_segment_x[10];
    double probe_segment_y[10];
    double probe_segment_z[10];
    double probe_segment_px[10];
    double probe_segment_py[10];
    double probe_segment_pz[10];
    double probe_segment_chiSquared[10];
    double probe_segment_numberDoF[10];
    double probe_segment_sector[10];
    double probe_segment_chamberIndex[10];
    double probe_segment_etaIndex[10];
    double probe_segment_nPrecisionHits[10];
    double probe_segment_nPhiLayers[10];
    double probe_segment_nTrigEtaLayers[10];
    //
    std::vector < int >* probe_L1_pass;
    std::vector < double >* probe_L1_eta;
    std::vector < double >* probe_L1_phi;
    std::vector < double >* probe_L1_dR;
    std::vector < double >* probe_L1_thrValue;
    std::vector < int >* probe_L1_roiNum;
    std::vector < int >* probe_L1_thrNumber;
    std::vector < bool >* probe_L1_isMoreCandInRoI;
    //
    std::vector < int >* probe_SA_pass;
    std::vector < double >* probe_SA_pt;
    std::vector < double >* probe_SA_eta;
    std::vector < double >* probe_SA_phi;
    std::vector < double >* probe_SA_etaBE;
    std::vector < double >* probe_SA_phiBE;
    std::vector < double >* probe_SA_etaMS;
    std::vector < double >* probe_SA_phiMS;
    std::vector < double >* probe_SA_tgcPt;
    std::vector < double >* probe_SA_ptBarrelRadius;
    std::vector < double >* probe_SA_ptBarrelSagitta;
    std::vector < double >* probe_SA_ptEndcapAlpha;
    std::vector < double >* probe_SA_ptEndcapBeta;
    std::vector < double >* probe_SA_ptEndcapRadius;
    std::vector < double >* probe_SA_ptCSC;
    std::vector < int >* probe_SA_sAddress;
    std::vector < float >* probe_SA_roiEta;
    std::vector < float >* probe_SA_roiPhi;
    std::vector < int >* probe_SA_isRpcFailure;
    std::vector < int >* probe_SA_isTgcFailure;
    //the measured radious of the muon in one particular super point
    std::vector < double >* probe_SA_superPointR_BI;
    std::vector < double >* probe_SA_superPointR_BM;
    std::vector < double >* probe_SA_superPointR_BO;
    std::vector < double >* probe_SA_superPointR_EI;
    std::vector < double >* probe_SA_superPointR_EM;
    std::vector < double >* probe_SA_superPointR_EO;
    std::vector < double >* probe_SA_superPointR_EE;
    std::vector < double >* probe_SA_superPointR_CSC;
    std::vector < double >* probe_SA_superPointR_BEE;
    std::vector < double >* probe_SA_superPointR_BME;
    //the measured Z position of the muon in one particular super point
    std::vector < double >* probe_SA_superPointZ_BI;
    std::vector < double >* probe_SA_superPointZ_BM;
    std::vector < double >* probe_SA_superPointZ_BO;
    std::vector < double >* probe_SA_superPointZ_EI;
    std::vector < double >* probe_SA_superPointZ_EM;
    std::vector < double >* probe_SA_superPointZ_EO;
    std::vector < double >* probe_SA_superPointZ_EE;
    std::vector < double >* probe_SA_superPointZ_CSC;
    std::vector < double >* probe_SA_superPointZ_BEE;
    std::vector < double >* probe_SA_superPointZ_BME;
    //the measured slope of the muon in one particular super point
    std::vector < double >* probe_SA_superPointSlope_BI;
    std::vector < double >* probe_SA_superPointSlope_BM;
    std::vector < double >* probe_SA_superPointSlope_BO;
    std::vector < double >* probe_SA_superPointSlope_EI;
    std::vector < double >* probe_SA_superPointSlope_EM;
    std::vector < double >* probe_SA_superPointSlope_EO;
    std::vector < double >* probe_SA_superPointSlope_EE;
    std::vector < double >* probe_SA_superPointSlope_CSC;
    std::vector < double >* probe_SA_superPointSlope_BEE;
    std::vector < double >* probe_SA_superPointSlope_BME;
    //the measured intercept of the muon in one particular super point
    std::vector < double >* probe_SA_superPointIntercept_BI;
    std::vector < double >* probe_SA_superPointIntercept_BM;
    std::vector < double >* probe_SA_superPointIntercept_BO;
    std::vector < double >* probe_SA_superPointIntercept_EI;
    std::vector < double >* probe_SA_superPointIntercept_EM;
    std::vector < double >* probe_SA_superPointIntercept_EO;
    std::vector < double >* probe_SA_superPointIntercept_EE;
    std::vector < double >* probe_SA_superPointIntercept_CSC;
    std::vector < double >* probe_SA_superPointIntercept_BEE;
    std::vector < double >* probe_SA_superPointIntercept_BME;
    //the chi2 of the fit in one particular super point
    std::vector < double >* probe_SA_superPointChi2_BI;
    std::vector < double >* probe_SA_superPointChi2_BM;
    std::vector < double >* probe_SA_superPointChi2_BO;
    std::vector < double >* probe_SA_superPointChi2_EI;
    std::vector < double >* probe_SA_superPointChi2_EM;
    std::vector < double >* probe_SA_superPointChi2_EO;
    std::vector < double >* probe_SA_superPointChi2_EE;
    std::vector < double >* probe_SA_superPointChi2_CSC;
    std::vector < double >* probe_SA_superPointChi2_BEE;
    std::vector < double >* probe_SA_superPointChi2_BME;
    //RPC hits
    std::vector < std::vector < float > >* probe_SA_rpcHitX;
    std::vector < std::vector < float > >* probe_SA_rpcHitY;
    std::vector < std::vector < float > >* probe_SA_rpcHitZ;
    std::vector < std::vector < float > >* probe_SA_rpcHitR;
    std::vector < std::vector < float > >* probe_SA_rpcHitEta;
    std::vector < std::vector < float > >* probe_SA_rpcHitPhi;
    std::vector < std::vector < uint32_t > >* probe_SA_rpcHitLayer;
    std::vector < std::vector < std::string > >* probe_SA_rpcHitStationName;
    std::vector < std::vector < uint32_t > >* probe_SA_rpcHitMeasPhi;
    //TGC hits
    std::vector < std::vector < float > >* probe_SA_tgcHitZ;
    std::vector < std::vector < float > >* probe_SA_tgcHitR;
    std::vector < std::vector < float > >* probe_SA_tgcHitEta;
    std::vector < std::vector < float > >* probe_SA_tgcHitPhi;
    std::vector < std::vector < float > >* probe_SA_tgcHitWidth;
    std::vector < std::vector < int > >* probe_SA_tgcHitStationNum;
    std::vector < std::vector < bool > >* probe_SA_tgcHitIsStrip;
    std::vector < std::vector < int > >* probe_SA_tgcHitBCTag;
    std::vector < std::vector < bool > >* probe_SA_tgcHitInRoad;
    //MDT hits
    std::vector < std::vector < int > >* probe_SA_mdtHitIsOutlier;
    std::vector < std::vector < int > >* probe_SA_mdtHitChamber;
    std::vector < std::vector < float > >* probe_SA_mdtHitR;
    std::vector < std::vector < float > >* probe_SA_mdtHitZ;
    std::vector < std::vector < float > >* probe_SA_mdtHitPhi;
    std::vector < std::vector < float > >* probe_SA_mdtHitResidual;

    std::vector < std::vector < float > >* probe_SA_roadAw;
    std::vector < std::vector < float > >* probe_SA_roadBw;
    std::vector < std::vector < float > >* probe_SA_zMin;
    std::vector < std::vector < float > >* probe_SA_zMax;
    std::vector < std::vector < float > >* probe_SA_rMin;
    std::vector < std::vector < float > >* probe_SA_rMax;
    std::vector < std::vector < float > >* probe_SA_etaMin;
    std::vector < std::vector < float > >* probe_SA_etaMax;
    //NSW hits
    std::vector < std::vector < float > >* probe_SA_stgcClusterR;
    std::vector < std::vector < float > >* probe_SA_stgcClusterZ;
    std::vector < std::vector < float > >* probe_SA_stgcClusterEta;
    std::vector < std::vector < float > >* probe_SA_stgcClusterPhi;
    std::vector < std::vector < float > >* probe_SA_stgcClusterResidualR;
    std::vector < std::vector < float > >* probe_SA_stgcClusterResidualPhi;
    std::vector < std::vector < int > >* probe_SA_stgcClusterStationEta;
    std::vector < std::vector < int > >* probe_SA_stgcClusterStationPhi;
    std::vector < std::vector < int > >* probe_SA_stgcClusterStationName;
    std::vector < std::vector < int > >* probe_SA_stgcClusterType;
    std::vector < std::vector < int > >* probe_SA_stgcClusterIsOutlier;
    std::vector < std::vector < unsigned int > >* probe_SA_stgcClusterLayer;
    std::vector < std::vector < float > >* probe_SA_mmClusterR;
    std::vector < std::vector < float > >* probe_SA_mmClusterZ;
    std::vector < std::vector < float > >* probe_SA_mmClusterEta;
    std::vector < std::vector < float > >* probe_SA_mmClusterPhi;
    std::vector < std::vector < float > >* probe_SA_mmClusterResidualR;
    std::vector < std::vector < float > >* probe_SA_mmClusterResidualPhi;
    std::vector < std::vector < int > >* probe_SA_mmClusterStationEta;
    std::vector < std::vector < int > >* probe_SA_mmClusterStationPhi;
    std::vector < std::vector < int > >* probe_SA_mmClusterStationName;
    std::vector < std::vector < int > >* probe_SA_mmClusterIsOutlier;
    std::vector < std::vector < unsigned int > >* probe_SA_mmClusterLayer;
    //
    //
    std::vector < int >* probe_CB_pass;
    std::vector < double >* probe_CB_pt;
    std::vector < double >* probe_CB_eta;
    std::vector < double >* probe_CB_phi;
    //
    std::vector < int >* probe_EF_pass;
    std::vector < double >* probe_EF_pt;
    std::vector < double >* probe_EF_eta;
    std::vector < double >* probe_EF_phi;

};

#endif

