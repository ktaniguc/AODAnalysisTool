#include <iostream>
#include <fstream>
#include <vector>

#include "CalcEfficiency/EventTreeMT.h"
#include "CalcEfficiency/Utils.h"

#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TLorentzVector.h"

EventTreeMT::EventTreeMT() {
}

EventTreeMT::~EventTreeMT() {
}

int EventTreeMT::initialize( TString outfile = "test.root" ) {

  //initialize the event tree
  m_file	= new TFile( outfile, "recreate" );
  //m_file	= outfile;
  m_tree 	= new TTree( "t_tap", "TrigMuonTagAndProbe" );
  //m_tree    = outtree;
  //m_tree->SetDirectory( m_file );

  //--------------------------------------------------
  // VARIABLE SET UP
  //--------------------------------------------------
	
  // Event info
  eventNumber = -1;
  runNumber   = -1;
  lumiBlock   = -1;
  averageInteractionsPerCrossing = -1;
  // Tag and probe variables
  tpsumReqdRl1      = 99999;
  tpextdR           = 99999;
  invMass           = 99999;
  n_trig            = 0;
  trigname          = new std::vector<std::string>();   
  passtrig          = new std::vector<int>();
  L1trigname        = new std::vector<std::string>();   
  L1_passbyEvent    = new std::vector < bool > ();
  //
  tag_pt           = -99999;
  tag_eta          = -99999;
  tag_phi          = -99999;
  tag_extEta       = -99999;
  tag_extPhi       = -99999;
  tag_d0           = -99999;
  tag_z0           = -99999;
  tag_L1_pass  = new std::vector < int > ();
  tag_L1_roiNum = new std::vector < int > ();
  tag_SA_pass  = new std::vector < int > ();
  tag_CB_pass  = new std::vector < int > ();
  tag_EF_pass  = new std::vector < int > ();
  //
  probe_pt           = -99999;
  probe_eta          = -99999;
  probe_phi          = -99999;
  probe_extEta       = -99999;
  probe_extPhi       = -99999;
  probe_d0           = -99999;
  probe_z0           = -99999;
  probe_segment_n   = -99999.;
  for(int i = 0 ; i < 10 ; i++ ){
    probe_segment_x[i]              = -99999.;
    probe_segment_y[i]              = -99999.;
    probe_segment_z[i]              = -99999.;
    probe_segment_px[i]             = -99999.;
    probe_segment_py[i]             = -99999.;
    probe_segment_pz[i]             = -99999.;
    probe_segment_chiSquared[i]     = -99999.;
    probe_segment_numberDoF[i]      = -99999.;
    probe_segment_sector[i]         = -99999.;
    probe_segment_chamberIndex[i]   = -99999.;
    probe_segment_etaIndex[i]       = -99999.;
    probe_segment_nPrecisionHits[i] = -99999.;
    probe_segment_nPhiLayers[i]     = -99999.;
    probe_segment_nTrigEtaLayers[i] = -99999.;
  }
  //
  probe_L1_pass  = new std::vector < int > ();
  probe_L1_eta   = new std::vector < double > ();
  probe_L1_phi   = new std::vector < double > ();
  probe_L1_dR   = new std::vector < double > ();
  probe_L1_roiNum   = new std::vector < int > ();
  probe_L1_thrValue   = new std::vector < double > ();
  probe_L1_thrNumber   = new std::vector < int > ();
  probe_L1_isMoreCandInRoI   = new std::vector < bool > ();
  //
  probe_SA_pass  = new std::vector < int > ();
  probe_SA_pt    = new std::vector < double > ();
  probe_SA_eta   = new std::vector < double > ();
  probe_SA_phi   = new std::vector < double > ();
  probe_SA_etaBE   = new std::vector < double > ();
  probe_SA_phiBE   = new std::vector < double > ();
  probe_SA_etaMS   = new std::vector < double > ();
  probe_SA_phiMS   = new std::vector < double > ();
  //pT from different algorithms
  probe_SA_tgcPt = new std::vector < double > ();
  probe_SA_ptBarrelRadius = new std::vector < double > ();
  probe_SA_ptBarrelSagitta = new std::vector < double > ();
  probe_SA_ptEndcapAlpha = new std::vector < double > ();
  probe_SA_ptEndcapBeta = new std::vector < double > ();
  probe_SA_ptEndcapRadius = new std::vector < double > ();
  probe_SA_ptCSC = new std::vector < double > ();
  probe_SA_sAddress = new std::vector < int > ();

  probe_SA_roiEta = new std::vector < float > ();
  probe_SA_roiPhi = new std::vector < float > ();
  probe_SA_isRpcFailure = new std::vector < int > ();
  probe_SA_isTgcFailure = new std::vector < int > ();

  //the measured radious of the muon in one particular super point
  probe_SA_superPointR_BI  = new std::vector < double > ();
  probe_SA_superPointR_BM  = new std::vector < double > ();
  probe_SA_superPointR_BO  = new std::vector < double > ();
  probe_SA_superPointR_EI  = new std::vector < double > ();
  probe_SA_superPointR_EM  = new std::vector < double > ();
  probe_SA_superPointR_EO  = new std::vector < double > ();
  probe_SA_superPointR_EE  = new std::vector < double > ();
  probe_SA_superPointR_CSC = new std::vector < double > ();
  probe_SA_superPointR_BEE = new std::vector < double > ();
  probe_SA_superPointR_BME = new std::vector < double > ();
  //the measured Z position of the muon in one particular super point
  probe_SA_superPointZ_BI  = new std::vector < double > ();
  probe_SA_superPointZ_BM  = new std::vector < double > ();
  probe_SA_superPointZ_BO  = new std::vector < double > ();
  probe_SA_superPointZ_EI  = new std::vector < double > ();
  probe_SA_superPointZ_EM  = new std::vector < double > ();
  probe_SA_superPointZ_EO  = new std::vector < double > ();
  probe_SA_superPointZ_EE  = new std::vector < double > ();
  probe_SA_superPointZ_CSC = new std::vector < double > ();
  probe_SA_superPointZ_BEE = new std::vector < double > ();
  probe_SA_superPointZ_BME = new std::vector < double > ();
  //the measured slope of the muon in one particular super point
  probe_SA_superPointSlope_BI  = new std::vector < double > ();
  probe_SA_superPointSlope_BM  = new std::vector < double > ();
  probe_SA_superPointSlope_BO  = new std::vector < double > ();
  probe_SA_superPointSlope_EI  = new std::vector < double > ();
  probe_SA_superPointSlope_EM  = new std::vector < double > ();
  probe_SA_superPointSlope_EO  = new std::vector < double > ();
  probe_SA_superPointSlope_EE  = new std::vector < double > ();
  probe_SA_superPointSlope_CSC = new std::vector < double > ();
  probe_SA_superPointSlope_BEE = new std::vector < double > ();
  probe_SA_superPointSlope_BME = new std::vector < double > ();
  //the measured intercept of the muon in one particular super point
  probe_SA_superPointIntercept_BI  = new std::vector < double > ();
  probe_SA_superPointIntercept_BM  = new std::vector < double > ();
  probe_SA_superPointIntercept_BO  = new std::vector < double > ();
  probe_SA_superPointIntercept_EI  = new std::vector < double > ();
  probe_SA_superPointIntercept_EM  = new std::vector < double > ();
  probe_SA_superPointIntercept_EO  = new std::vector < double > ();
  probe_SA_superPointIntercept_EE  = new std::vector < double > ();
  probe_SA_superPointIntercept_CSC = new std::vector < double > ();
  probe_SA_superPointIntercept_BEE = new std::vector < double > ();
  probe_SA_superPointIntercept_BME = new std::vector < double > ();
  //the chi2 of the muon in one particular super point
  probe_SA_superPointChi2_BI  = new std::vector < double > ();
  probe_SA_superPointChi2_BM  = new std::vector < double > ();
  probe_SA_superPointChi2_BO  = new std::vector < double > ();
  probe_SA_superPointChi2_EI  = new std::vector < double > ();
  probe_SA_superPointChi2_EM  = new std::vector < double > ();
  probe_SA_superPointChi2_EO  = new std::vector < double > ();
  probe_SA_superPointChi2_EE  = new std::vector < double > ();
  probe_SA_superPointChi2_CSC = new std::vector < double > ();
  probe_SA_superPointChi2_BEE = new std::vector < double > ();
  probe_SA_superPointChi2_BME = new std::vector < double > ();
  //
  probe_SA_rpcHitX   = new std::vector < std::vector < float > > ();
  probe_SA_rpcHitY   = new std::vector < std::vector < float > > ();
  probe_SA_rpcHitZ   = new std::vector < std::vector < float > > ();
  probe_SA_rpcHitR   = new std::vector < std::vector < float > > ();
  probe_SA_rpcHitEta = new std::vector < std::vector < float > > ();
  probe_SA_rpcHitPhi = new std::vector < std::vector < float > > ();
  probe_SA_rpcHitMeasPhi = new std::vector < std::vector < uint32_t > > ();
  probe_SA_rpcHitLayer = new std::vector < std::vector < uint32_t > > ();
  probe_SA_rpcHitStationName = new std::vector < std::vector < string > > ();
  //
  probe_SA_tgcHitZ = new std::vector < std::vector < float > > ();
  probe_SA_tgcHitR = new std::vector < std::vector < float > > ();
  probe_SA_tgcHitEta = new std::vector < std::vector < float > > ();
  probe_SA_tgcHitPhi = new std::vector < std::vector < float > > ();
  probe_SA_tgcHitWidth = new std::vector < std::vector < float > > ();
  probe_SA_tgcHitStationNum = new std::vector < std::vector < int > > ();
  probe_SA_tgcHitIsStrip = new std::vector < std::vector < bool > > ();
  probe_SA_tgcHitBCTag = new std::vector < std::vector < int > > ();
  probe_SA_tgcHitInRoad = new std::vector < std::vector < bool > > ();
  //
  probe_SA_mdtHitIsOutlier = new std::vector < std::vector < int > > ();
  probe_SA_mdtHitChamber   = new std::vector < std::vector < int > > ();
  probe_SA_mdtHitR         = new std::vector < std::vector < float > > ();
  probe_SA_mdtHitZ         = new std::vector < std::vector < float > > ();
  probe_SA_mdtHitPhi       = new std::vector < std::vector < float > > ();
  probe_SA_mdtHitResidual  = new std::vector < std::vector < float > > ();

  probe_SA_roadAw = new std::vector < std::vector < float > > ();
  probe_SA_roadBw = new std::vector < std::vector < float > > ();
  probe_SA_zMin   = new std::vector < std::vector < float > > ();
  probe_SA_zMax   = new std::vector < std::vector < float > > ();
  probe_SA_rMin   = new std::vector < std::vector < float > > ();
  probe_SA_rMax   = new std::vector < std::vector < float > > ();
  probe_SA_etaMin = new std::vector < std::vector < float > > ();
  probe_SA_etaMax = new std::vector < std::vector < float > > ();
  //NSW
  probe_SA_stgcClusterR = new std::vector < std::vector < float > > ();
  probe_SA_stgcClusterZ = new std::vector < std::vector < float > > ();
  probe_SA_stgcClusterEta = new std::vector < std::vector < float > > ();
  probe_SA_stgcClusterPhi = new std::vector < std::vector < float > > ();
  probe_SA_stgcClusterResidualR = new std::vector < std::vector < float > > ();
  probe_SA_stgcClusterResidualPhi = new std::vector < std::vector < float > > ();
  probe_SA_stgcClusterStationEta = new std::vector < std::vector < int > > ();
  probe_SA_stgcClusterStationPhi = new std::vector < std::vector < int > > ();
  probe_SA_stgcClusterStationName = new std::vector < std::vector < int > > ();
  probe_SA_stgcClusterType = new std::vector < std::vector < int > > ();
  probe_SA_stgcClusterIsOutlier = new std::vector < std::vector < int > > ();
  probe_SA_stgcClusterLayer = new std::vector < std::vector < unsigned int > > ();
  probe_SA_mmClusterR = new std::vector < std::vector < float > > ();
  probe_SA_mmClusterZ = new std::vector < std::vector < float > > ();
  probe_SA_mmClusterEta = new std::vector < std::vector < float > > ();
  probe_SA_mmClusterPhi = new std::vector < std::vector < float > > ();
  probe_SA_mmClusterResidualR = new std::vector < std::vector < float > > ();
  probe_SA_mmClusterResidualPhi = new std::vector < std::vector < float > > ();
  probe_SA_mmClusterStationEta = new std::vector < std::vector < int > > ();
  probe_SA_mmClusterStationPhi = new std::vector < std::vector < int > > ();
  probe_SA_mmClusterStationName = new std::vector < std::vector < int > > ();
  probe_SA_mmClusterIsOutlier = new std::vector < std::vector < int > > ();
  probe_SA_mmClusterLayer = new std::vector < std::vector < unsigned int > > ();
  //
  probe_CB_pass  = new std::vector < int > ();
  probe_CB_pt    = new std::vector < double > ();
  probe_CB_eta   = new std::vector < double > ();
  probe_CB_phi   = new std::vector < double > ();
  //
  probe_EF_pass  = new std::vector < int > ();
  probe_EF_pt    = new std::vector < double > ();
  probe_EF_eta   = new std::vector < double > ();
  probe_EF_phi   = new std::vector < double > ();
    
  //--------------------------------------------------
  // BRANCH SET UP
  //--------------------------------------------------
 
  //Event info 
  m_tree->Branch( "eventNumber",       &eventNumber );
  m_tree->Branch( "runNumber",       &runNumber );
  m_tree->Branch( "lumiBlock",       &lumiBlock );
  m_tree->Branch( "averageInteractionsPerCrossing",       &averageInteractionsPerCrossing );

  //tag and probe
  m_tree->Branch( "tpsumReqdRl1",       &tpsumReqdRl1 );
  m_tree->Branch( "tpextdR",       &tpextdR );
  m_tree->Branch( "invMass",       &invMass );
  m_tree->Branch( "trigname",       &trigname );
  m_tree->Branch( "passtrig",       &passtrig );
  m_tree->Branch( "L1trigname",       &L1trigname );
  m_tree->Branch( "L1_passbyEvent",   &L1_passbyEvent );
  m_tree->Branch( "n_trig",       &n_trig );
  //tag
  m_tree->Branch( "tag_pt",     &tag_pt );
  m_tree->Branch( "tag_eta",    &tag_eta );
  m_tree->Branch( "tag_phi",    &tag_phi );
  m_tree->Branch( "tag_extEta",    &tag_extEta );
  m_tree->Branch( "tag_extPhi",    &tag_extPhi );
  m_tree->Branch( "tag_d0",     &tag_d0 );
  m_tree->Branch( "tag_z0",     &tag_z0 );
  m_tree->Branch( "tag_L1_pass",   &tag_L1_pass );
  m_tree->Branch( "tag_L1_roiNum",   &tag_L1_roiNum );
  m_tree->Branch( "tag_SA_pass",   &tag_SA_pass );
  m_tree->Branch( "tag_CB_pass",   &tag_CB_pass );
  m_tree->Branch( "tag_EF_pass",   &tag_EF_pass );
  //probe
  m_tree->Branch( "probe_pt",     &probe_pt );
  m_tree->Branch( "probe_eta",    &probe_eta );
  m_tree->Branch( "probe_phi",    &probe_phi );
  m_tree->Branch( "probe_extEta",    &probe_extEta );
  m_tree->Branch( "probe_extPhi",    &probe_extPhi );
  m_tree->Branch( "probe_d0",     &probe_d0 );
  m_tree->Branch( "probe_z0",     &probe_z0 );
  m_tree->Branch( "probe_segment_n",    &probe_segment_n );
  m_tree->Branch( "probe_segment_x",              probe_segment_x,              Form("probe_segment_x[%d]/D",10) );
  m_tree->Branch( "probe_segment_y",              probe_segment_y,              Form("probe_segment_y[%d]/D",10) );
  m_tree->Branch( "probe_segment_z",              probe_segment_z,              Form("probe_segment_z[%d]/D",10) );
  m_tree->Branch( "probe_segment_px",             probe_segment_px,             Form("probe_segment_px[%d]/D",10) );
  m_tree->Branch( "probe_segment_py",             probe_segment_py,             Form("probe_segment_py[%d]/D",10) );
  m_tree->Branch( "probe_segment_pz",             probe_segment_pz,             Form("probe_segment_pz[%d]/D",10) );
  m_tree->Branch( "probe_segment_chiSquared",     probe_segment_chiSquared,     Form("probe_segment_chiSquared[%d]/D",10) );
  m_tree->Branch( "probe_segment_numberDoF",      probe_segment_numberDoF,      Form("probe_segment_numberDoF[%d]/D",10) );
  m_tree->Branch( "probe_segment_sector",         probe_segment_sector,         Form("probe_segment_sector[%d]/D",10) );
  m_tree->Branch( "probe_segment_chamberIndex",   probe_segment_chamberIndex,   Form("probe_segment_chamberIndex[%d]/D",10) );
  m_tree->Branch( "probe_segment_etaIndex",       probe_segment_etaIndex,       Form("probe_segment_etaIndex[%d]/D",10) );
  m_tree->Branch( "probe_segment_nPrecisionHits", probe_segment_nPrecisionHits, Form("probe_segment_nPrecisionHits[%d]/D",10) );
  m_tree->Branch( "probe_segment_nPhiLayers",     probe_segment_nPhiLayers,     Form("probe_segment_nPhiLayers[%d]/D",10) );
  m_tree->Branch( "probe_segment_nTrigEtaLayers", probe_segment_nTrigEtaLayers, Form("probe_segment_nTrigEtaLayers[%d]/D",10) );
  //L1
  m_tree->Branch( "probe_L1_pass",   &probe_L1_pass );
  m_tree->Branch( "probe_L1_eta",    &probe_L1_eta );
  m_tree->Branch( "probe_L1_phi",    &probe_L1_phi );
  m_tree->Branch( "probe_L1_dR",    &probe_L1_dR );
  m_tree->Branch( "probe_L1_thrValue",    &probe_L1_thrValue );
  m_tree->Branch( "probe_L1_roiNum",    &probe_L1_roiNum );
  m_tree->Branch( "probe_L1_thrNumber",    &probe_L1_thrNumber );
  m_tree->Branch( "probe_L1_isMoreCandInRoI",    &probe_L1_isMoreCandInRoI );
  //SA
  m_tree->Branch( "probe_SA_pass",   &probe_SA_pass );
  m_tree->Branch( "probe_SA_pt",     &probe_SA_pt );
  m_tree->Branch( "probe_SA_eta",    &probe_SA_eta );
  m_tree->Branch( "probe_SA_phi",    &probe_SA_phi );
  m_tree->Branch( "probe_SA_etaBE",    &probe_SA_etaBE );
  m_tree->Branch( "probe_SA_phiBE",    &probe_SA_phiBE );
  m_tree->Branch( "probe_SA_etaMS",    &probe_SA_etaMS );
  m_tree->Branch( "probe_SA_phiMS",    &probe_SA_phiMS );
  //
  m_tree->Branch( "probe_SA_tgcPt",  &probe_SA_tgcPt );
  m_tree->Branch( "probe_SA_ptBarrelRadius",  &probe_SA_ptBarrelRadius );
  m_tree->Branch( "probe_SA_ptBarrelSagitta",  &probe_SA_ptBarrelSagitta );
  m_tree->Branch( "probe_SA_ptEndcapAlpha",  &probe_SA_ptEndcapAlpha );
  m_tree->Branch( "probe_SA_ptEndcapBeta",  &probe_SA_ptEndcapBeta );
  m_tree->Branch( "probe_SA_ptEndcapRadius",  &probe_SA_ptEndcapRadius );
  m_tree->Branch( "probe_SA_ptCSC",  &probe_SA_ptCSC );
  m_tree->Branch( "probe_SA_sAddress",  &probe_SA_sAddress );

  m_tree->Branch( "probe_SA_roiEta",  &probe_SA_roiEta );
  m_tree->Branch( "probe_SA_roiPhi",  &probe_SA_roiPhi );
  m_tree->Branch( "probe_SA_isRpcFailure",  &probe_SA_isRpcFailure );
  m_tree->Branch( "probe_SA_isTgcFailure",  &probe_SA_isTgcFailure );
  //
  m_tree->Branch( "probe_SA_superPointR_BI",   &probe_SA_superPointR_BI );
  m_tree->Branch( "probe_SA_superPointR_BM",   &probe_SA_superPointR_BM );
  m_tree->Branch( "probe_SA_superPointR_BO",   &probe_SA_superPointR_BO );
  m_tree->Branch( "probe_SA_superPointR_EI",   &probe_SA_superPointR_EI );
  m_tree->Branch( "probe_SA_superPointR_EM",   &probe_SA_superPointR_EM );
  m_tree->Branch( "probe_SA_superPointR_EO",   &probe_SA_superPointR_EO );
  m_tree->Branch( "probe_SA_superPointR_EE",   &probe_SA_superPointR_EE );
  m_tree->Branch( "probe_SA_superPointR_CSC",  &probe_SA_superPointR_CSC );
  m_tree->Branch( "probe_SA_superPointR_BEE",  &probe_SA_superPointR_BEE );
  m_tree->Branch( "probe_SA_superPointR_BME",  &probe_SA_superPointR_BME );
  //
  m_tree->Branch( "probe_SA_superPointZ_BI",   &probe_SA_superPointZ_BI );
  m_tree->Branch( "probe_SA_superPointZ_BM",   &probe_SA_superPointZ_BM );
  m_tree->Branch( "probe_SA_superPointZ_BO",   &probe_SA_superPointZ_BO );
  m_tree->Branch( "probe_SA_superPointZ_EI",   &probe_SA_superPointZ_EI );
  m_tree->Branch( "probe_SA_superPointZ_EM",   &probe_SA_superPointZ_EM );
  m_tree->Branch( "probe_SA_superPointZ_EO",   &probe_SA_superPointZ_EO );
  m_tree->Branch( "probe_SA_superPointZ_EE",   &probe_SA_superPointZ_EE );
  m_tree->Branch( "probe_SA_superPointZ_CSC",  &probe_SA_superPointZ_CSC );
  m_tree->Branch( "probe_SA_superPointZ_BEE",  &probe_SA_superPointZ_BEE );
  m_tree->Branch( "probe_SA_superPointZ_BME",  &probe_SA_superPointZ_BME );
  //
  m_tree->Branch( "probe_SA_superPointSlope_BI",   &probe_SA_superPointSlope_BI );
  m_tree->Branch( "probe_SA_superPointSlope_BM",   &probe_SA_superPointSlope_BM );
  m_tree->Branch( "probe_SA_superPointSlope_BO",   &probe_SA_superPointSlope_BO );
  m_tree->Branch( "probe_SA_superPointSlope_EI",   &probe_SA_superPointSlope_EI );
  m_tree->Branch( "probe_SA_superPointSlope_EM",   &probe_SA_superPointSlope_EM );
  m_tree->Branch( "probe_SA_superPointSlope_EO",   &probe_SA_superPointSlope_EO );
  m_tree->Branch( "probe_SA_superPointSlope_EE",   &probe_SA_superPointSlope_EE );
  m_tree->Branch( "probe_SA_superPointSlope_CSC",  &probe_SA_superPointSlope_CSC );
  m_tree->Branch( "probe_SA_superPointSlope_BEE",  &probe_SA_superPointSlope_BEE );
  m_tree->Branch( "probe_SA_superPointSlope_BME",  &probe_SA_superPointSlope_BME );
  //
  m_tree->Branch( "probe_SA_superPointIntercept_BI",   &probe_SA_superPointIntercept_BI );
  m_tree->Branch( "probe_SA_superPointIntercept_BM",   &probe_SA_superPointIntercept_BM );
  m_tree->Branch( "probe_SA_superPointIntercept_BO",   &probe_SA_superPointIntercept_BO );
  m_tree->Branch( "probe_SA_superPointIntercept_EI",   &probe_SA_superPointIntercept_EI );
  m_tree->Branch( "probe_SA_superPointIntercept_EM",   &probe_SA_superPointIntercept_EM );
  m_tree->Branch( "probe_SA_superPointIntercept_EO",   &probe_SA_superPointIntercept_EO );
  m_tree->Branch( "probe_SA_superPointIntercept_EE",   &probe_SA_superPointIntercept_EE );
  m_tree->Branch( "probe_SA_superPointIntercept_CSC",  &probe_SA_superPointIntercept_CSC );
  m_tree->Branch( "probe_SA_superPointIntercept_BEE",  &probe_SA_superPointIntercept_BEE );
  m_tree->Branch( "probe_SA_superPointIntercept_BME",  &probe_SA_superPointIntercept_BME );
  //
  m_tree->Branch( "probe_SA_superPointChi2_BI",   &probe_SA_superPointChi2_BI );
  m_tree->Branch( "probe_SA_superPointChi2_BM",   &probe_SA_superPointChi2_BM );
  m_tree->Branch( "probe_SA_superPointChi2_BO",   &probe_SA_superPointChi2_BO );
  m_tree->Branch( "probe_SA_superPointChi2_EI",   &probe_SA_superPointChi2_EI );
  m_tree->Branch( "probe_SA_superPointChi2_EM",   &probe_SA_superPointChi2_EM );
  m_tree->Branch( "probe_SA_superPointChi2_EO",   &probe_SA_superPointChi2_EO );
  m_tree->Branch( "probe_SA_superPointChi2_EE",   &probe_SA_superPointChi2_EE );
  m_tree->Branch( "probe_SA_superPointChi2_CSC",  &probe_SA_superPointChi2_CSC );
  m_tree->Branch( "probe_SA_superPointChi2_BEE",  &probe_SA_superPointChi2_BEE );
  m_tree->Branch( "probe_SA_superPointChi2_BME",  &probe_SA_superPointChi2_BME );
  // 
  m_tree->Branch( "probe_SA_rpcHitX",  &probe_SA_rpcHitX );
  m_tree->Branch( "probe_SA_rpcHitY",  &probe_SA_rpcHitY );
  m_tree->Branch( "probe_SA_rpcHitZ",  &probe_SA_rpcHitZ );
  m_tree->Branch( "probe_SA_rpcHitR",  &probe_SA_rpcHitR );
  m_tree->Branch( "probe_SA_rpcHitEta",  &probe_SA_rpcHitEta );
  m_tree->Branch( "probe_SA_rpcHitPhi",  &probe_SA_rpcHitPhi );
  m_tree->Branch( "probe_SA_rpcHitMeasPhi",  &probe_SA_rpcHitMeasPhi );
  m_tree->Branch( "probe_SA_rpcHitLayer",  &probe_SA_rpcHitLayer );
  m_tree->Branch( "probe_SA_rpcHitStationName",  &probe_SA_rpcHitStationName );
  //
  m_tree->Branch( "probe_SA_tgcHitZ",  &probe_SA_tgcHitZ );
  m_tree->Branch( "probe_SA_tgcHitR",  &probe_SA_tgcHitR );
  m_tree->Branch( "probe_SA_tgcHitEta",  &probe_SA_tgcHitEta );
  m_tree->Branch( "probe_SA_tgcHitPhi",  &probe_SA_tgcHitPhi );
  m_tree->Branch( "probe_SA_tgcHitWidth",  &probe_SA_tgcHitWidth );
  m_tree->Branch( "probe_SA_tgcHitStationNum",  &probe_SA_tgcHitStationNum );
  m_tree->Branch( "probe_SA_tgcHitIsStrip",  &probe_SA_tgcHitIsStrip );
  m_tree->Branch( "probe_SA_tgcHitBCTag",  &probe_SA_tgcHitBCTag );
  m_tree->Branch( "probe_SA_tgcHitInRoad",  &probe_SA_tgcHitInRoad );
  //
  m_tree->Branch( "probe_SA_mdtHitIsOutlier", &probe_SA_mdtHitIsOutlier );
  m_tree->Branch( "probe_SA_mdtHitChamber",   &probe_SA_mdtHitChamber   );
  m_tree->Branch( "probe_SA_mdtHitR",         &probe_SA_mdtHitR         );
  m_tree->Branch( "probe_SA_mdtHitZ",         &probe_SA_mdtHitZ         );
  m_tree->Branch( "probe_SA_mdtHitPhi",       &probe_SA_mdtHitPhi       );
  m_tree->Branch( "probe_SA_mdtHitResidual",  &probe_SA_mdtHitResidual  );
  //
  m_tree->Branch( "probe_SA_roadAw",  &probe_SA_roadAw );
  m_tree->Branch( "probe_SA_roadBw",  &probe_SA_roadBw );
  m_tree->Branch( "probe_SA_zMin",  &probe_SA_zMin );
  m_tree->Branch( "probe_SA_zMax",  &probe_SA_zMax );
  m_tree->Branch( "probe_SA_rMin",  &probe_SA_rMin );
  m_tree->Branch( "probe_SA_rMax",  &probe_SA_rMax );
  m_tree->Branch( "probe_SA_etaMin",  &probe_SA_etaMin );
  m_tree->Branch( "probe_SA_etaMax",  &probe_SA_etaMax );
  //NSW
  m_tree->Branch( "probe_SA_stgcClusterZ",  &probe_SA_stgcClusterZ );
  m_tree->Branch( "probe_SA_stgcClusterR",  &probe_SA_stgcClusterR );
  m_tree->Branch( "probe_SA_stgcClusterEta",  &probe_SA_stgcClusterEta );
  m_tree->Branch( "probe_SA_stgcClusterPhi",  &probe_SA_stgcClusterPhi );
  m_tree->Branch( "probe_SA_stgcClusterResidualR",  &probe_SA_stgcClusterResidualR );
  m_tree->Branch( "probe_SA_stgcClusterResidualPhi",  &probe_SA_stgcClusterResidualPhi );
  m_tree->Branch( "probe_SA_stgcClusterStationEta",  &probe_SA_stgcClusterStationEta );
  m_tree->Branch( "probe_SA_stgcClusterStationPhi",  &probe_SA_stgcClusterStationPhi );
  m_tree->Branch( "probe_SA_stgcClusterStationName",  &probe_SA_stgcClusterStationName );
  m_tree->Branch( "probe_SA_stgcClusterType",  &probe_SA_stgcClusterType );
  m_tree->Branch( "probe_SA_stgcClusterIsOutlier",  &probe_SA_stgcClusterIsOutlier );
  m_tree->Branch( "probe_SA_stgcClusterLayer",  &probe_SA_stgcClusterLayer );
  m_tree->Branch( "probe_SA_mmClusterZ",  &probe_SA_mmClusterZ );
  m_tree->Branch( "probe_SA_mmClusterR",  &probe_SA_mmClusterR );
  m_tree->Branch( "probe_SA_mmClusterEta",  &probe_SA_mmClusterEta );
  m_tree->Branch( "probe_SA_mmClusterPhi",  &probe_SA_mmClusterPhi );
  m_tree->Branch( "probe_SA_mmClusterResidualR",  &probe_SA_mmClusterResidualR );
  m_tree->Branch( "probe_SA_mmClusterResidualPhi",  &probe_SA_mmClusterResidualPhi );
  m_tree->Branch( "probe_SA_mmClusterStationEta",  &probe_SA_mmClusterStationEta );
  m_tree->Branch( "probe_SA_mmClusterStationPhi",  &probe_SA_mmClusterStationPhi );
  m_tree->Branch( "probe_SA_mmClusterStationName",  &probe_SA_mmClusterStationName );
  m_tree->Branch( "probe_SA_mmClusterIsOutlier",  &probe_SA_mmClusterIsOutlier );
  m_tree->Branch( "probe_SA_mmClusterLayer",  &probe_SA_mmClusterLayer );
  //CB
  m_tree->Branch( "probe_CB_pass",   &probe_CB_pass );
  m_tree->Branch( "probe_CB_pt",     &probe_CB_pt );
  m_tree->Branch( "probe_CB_eta",    &probe_CB_eta );
  m_tree->Branch( "probe_CB_phi",    &probe_CB_phi );
  //EF
  m_tree->Branch( "probe_EF_pass",   &probe_EF_pass );
  m_tree->Branch( "probe_EF_pt",     &probe_EF_pt );
  m_tree->Branch( "probe_EF_eta",    &probe_EF_eta );
  m_tree->Branch( "probe_EF_phi",    &probe_EF_phi );
  return 1;
}

void EventTreeMT::clear() {
  // clear the vector for branch
  trigname->clear();
  passtrig->clear();
  L1trigname->clear();
  L1_passbyEvent->clear();
  //
  tag_L1_pass->clear();
  tag_L1_roiNum->clear();
  tag_SA_pass->clear();
  tag_CB_pass->clear();
  tag_EF_pass->clear();
  //
  probe_L1_pass->clear();
  probe_L1_eta->clear();
  probe_L1_phi->clear();
  probe_L1_dR->clear();
  probe_L1_thrValue->clear();
  probe_L1_roiNum->clear();
  probe_L1_thrNumber->clear();
  probe_L1_isMoreCandInRoI->clear();
  //
  probe_SA_pass->clear();
  probe_SA_pt->clear();
  probe_SA_eta->clear();
  probe_SA_phi->clear();
  probe_SA_etaBE->clear();
  probe_SA_phiBE->clear();
  probe_SA_etaMS->clear();
  probe_SA_phiMS->clear();
  probe_SA_tgcPt->clear();
  probe_SA_ptBarrelRadius->clear();
  probe_SA_ptBarrelSagitta->clear();
  probe_SA_ptEndcapAlpha->clear();
  probe_SA_ptEndcapBeta->clear();
  probe_SA_ptEndcapRadius->clear();
  probe_SA_ptCSC->clear();
  probe_SA_sAddress->clear();
  probe_SA_roiEta->clear();
  probe_SA_roiPhi->clear();
  probe_SA_isRpcFailure->clear();
  probe_SA_isTgcFailure->clear();

  //the measured radious of the muon in one particular super point
  probe_SA_superPointR_BI->clear();
  probe_SA_superPointR_BM->clear();
  probe_SA_superPointR_BO->clear();
  probe_SA_superPointR_EI->clear();
  probe_SA_superPointR_EM->clear();
  probe_SA_superPointR_EO->clear();
  probe_SA_superPointR_EE->clear();
  probe_SA_superPointR_CSC->clear();
  probe_SA_superPointR_BEE->clear();
  probe_SA_superPointR_BME->clear();
  //the measured Z position of the muon in one particular super point
  probe_SA_superPointZ_BI->clear();
  probe_SA_superPointZ_BM->clear();
  probe_SA_superPointZ_BO->clear();
  probe_SA_superPointZ_EI->clear();
  probe_SA_superPointZ_EM->clear();
  probe_SA_superPointZ_EO->clear();
  probe_SA_superPointZ_EE->clear();
  probe_SA_superPointZ_CSC->clear();
  probe_SA_superPointZ_BEE->clear();
  probe_SA_superPointZ_BME->clear();
  //the measured slope of the muon in one particular super point
  probe_SA_superPointSlope_BI->clear();
  probe_SA_superPointSlope_BM->clear();
  probe_SA_superPointSlope_BO->clear();
  probe_SA_superPointSlope_EI->clear();
  probe_SA_superPointSlope_EM->clear();
  probe_SA_superPointSlope_EO->clear();
  probe_SA_superPointSlope_EE->clear();
  probe_SA_superPointSlope_CSC->clear();
  probe_SA_superPointSlope_BEE->clear();
  probe_SA_superPointSlope_BME->clear();
  //the measured intercept of the muon in one particular super point
  probe_SA_superPointIntercept_BI->clear();
  probe_SA_superPointIntercept_BM->clear();
  probe_SA_superPointIntercept_BO->clear();
  probe_SA_superPointIntercept_EI->clear();
  probe_SA_superPointIntercept_EM->clear();
  probe_SA_superPointIntercept_EO->clear();
  probe_SA_superPointIntercept_EE->clear();
  probe_SA_superPointIntercept_CSC->clear();
  probe_SA_superPointIntercept_BEE->clear();
  probe_SA_superPointIntercept_BME->clear();
  //the chi2 of the muon in one particular super point
  probe_SA_superPointChi2_BI->clear();
  probe_SA_superPointChi2_BM->clear();
  probe_SA_superPointChi2_BO->clear();
  probe_SA_superPointChi2_EI->clear();
  probe_SA_superPointChi2_EM->clear();
  probe_SA_superPointChi2_EO->clear();
  probe_SA_superPointChi2_EE->clear();
  probe_SA_superPointChi2_CSC->clear();
  probe_SA_superPointChi2_BEE->clear();
  probe_SA_superPointChi2_BME->clear();
  //
  probe_SA_rpcHitX->clear();
  probe_SA_rpcHitY->clear();
  probe_SA_rpcHitZ->clear();
  probe_SA_rpcHitR->clear();
  probe_SA_rpcHitEta->clear();
  probe_SA_rpcHitPhi->clear();
  probe_SA_rpcHitMeasPhi->clear();
  probe_SA_rpcHitLayer->clear();
  probe_SA_rpcHitStationName->clear();
  //
  probe_SA_tgcHitZ->clear();
  probe_SA_tgcHitR->clear();
  probe_SA_tgcHitEta->clear();
  probe_SA_tgcHitPhi->clear();
  probe_SA_tgcHitWidth->clear();
  probe_SA_tgcHitStationNum->clear();
  probe_SA_tgcHitIsStrip->clear();
  probe_SA_tgcHitBCTag->clear();
  probe_SA_tgcHitInRoad->clear();
  //
  probe_SA_mdtHitIsOutlier->clear();
  probe_SA_mdtHitChamber->clear();
  probe_SA_mdtHitR->clear();
  probe_SA_mdtHitZ->clear();
  probe_SA_mdtHitPhi->clear();
  probe_SA_mdtHitResidual->clear();
  //
  probe_SA_roadAw->clear();
  probe_SA_roadBw->clear();
  probe_SA_zMin->clear();
  probe_SA_zMax->clear();
  probe_SA_rMin->clear();
  probe_SA_rMax->clear();
  probe_SA_etaMin->clear();
  probe_SA_etaMax->clear();
  //
  probe_SA_stgcClusterR->clear();
  probe_SA_stgcClusterZ->clear();
  probe_SA_stgcClusterEta->clear();
  probe_SA_stgcClusterPhi->clear();
  probe_SA_stgcClusterResidualR->clear();
  probe_SA_stgcClusterResidualPhi->clear();
  probe_SA_stgcClusterStationEta->clear();
  probe_SA_stgcClusterStationPhi->clear();
  probe_SA_stgcClusterStationName->clear();
  probe_SA_stgcClusterType->clear();
  probe_SA_stgcClusterIsOutlier->clear();
  probe_SA_stgcClusterLayer->clear();
  probe_SA_mmClusterR->clear();
  probe_SA_mmClusterZ->clear();
  probe_SA_mmClusterEta->clear();
  probe_SA_mmClusterPhi->clear();
  probe_SA_mmClusterResidualR->clear();
  probe_SA_mmClusterResidualPhi->clear();
  probe_SA_mmClusterStationEta->clear();
  probe_SA_mmClusterStationPhi->clear();
  probe_SA_mmClusterStationName->clear();
  probe_SA_mmClusterIsOutlier->clear();
  probe_SA_mmClusterLayer->clear();
  //
  probe_CB_pass->clear();
  probe_CB_pt->clear();
  probe_CB_eta->clear();
  probe_CB_phi->clear();
  //
  probe_EF_pass->clear();
  probe_EF_pt->clear();
  probe_EF_eta->clear();
  probe_EF_phi->clear();
}

template void EventTreeMT::filltree<TagAndProbeMT>(TagAndProbeMT& tap, 
                                                   int evtNum,
                                                   int runNum,
                                                   int lumiBlk,
                                                   float avgIntPerXing);
template void EventTreeMT::filltree<TagAndProbe>(TagAndProbe& tap,
                                                 int evtNum,
                                                 int runNum,
                                                 int lumiBlk,
                                                 float avgIntPerXing);
template<typename TAP> void EventTreeMT::filltree( TAP& tap,
                                                   int evtNum,
                                                   int runNum,
                                                   int lumiBlk,
                                                   float avgIntPerXing)
{
// fill the variable vectors
  int probe_n = tap.m_vL1objects_probe.size();
  for( int i_probe = 0; i_probe < probe_n; i_probe++ ) {
    this->clear();
    //event info
    eventNumber = evtNum;
    runNumber = runNum;
    lumiBlock = lumiBlk;
    averageInteractionsPerCrossing = avgIntPerXing;
    //tag and probe
    n_trig = tap.m_nmesChain; 
    tpsumReqdRl1 = tap.m_requirements_tag[i_probe].reqdRl1 + tap.m_requirements_probe[i_probe].reqdRl1;
    tpextdR = tap.m_probe[i_probe].tpextdR;
    tag_pt = tap.m_tag[i_probe].pt;
    tag_eta = tap.m_tag[i_probe].eta;
    tag_phi = tap.m_tag[i_probe].phi;
    tag_extEta = tap.m_tag[i_probe].extEta;
    tag_extPhi = tap.m_tag[i_probe].extPhi;
    tag_d0 = tap.m_tag[i_probe].d0;
    tag_z0 = tap.m_tag[i_probe].z0;

    if((int)tap.m_tag_L1_pass.size() == n_trig) {
      for( int i_trig = 0; i_trig < n_trig; i_trig++ ){
	tag_L1_pass->push_back( tap.m_tag_L1_pass[i_trig] );
	tag_L1_roiNum->push_back( tap.m_tag_L1_roiNum[i_trig] );
	tag_SA_pass->push_back( tap.m_tag_SA_pass[i_trig] );
	tag_CB_pass->push_back( tap.m_tag_CB_pass[i_trig] );
	tag_EF_pass->push_back( tap.m_tag_EF_pass[i_trig] );
      }
    }

    probe_pt = tap.m_probe[i_probe].pt;
    probe_eta = tap.m_probe[i_probe].eta;
    probe_phi = tap.m_probe[i_probe].phi;
    probe_extEta = tap.m_probe[i_probe].extEta;
    probe_extPhi = tap.m_probe[i_probe].extPhi;
    probe_d0 = tap.m_probe[i_probe].d0;
    probe_z0 = tap.m_probe[i_probe].z0;
    probe_segment_n   = tap.m_probe[i_probe].segmentN;
    for(int i_seg = 0 ; i_seg < 10 ; i_seg++ ){
      probe_segment_x[i_seg]              = tap.m_probe[i_probe].segmentX[i_seg];
      probe_segment_y[i_seg]              = tap.m_probe[i_probe].segmentY[i_seg];
      probe_segment_z[i_seg]              = tap.m_probe[i_probe].segmentZ[i_seg];
      probe_segment_px[i_seg]             = tap.m_probe[i_probe].segmentPx[i_seg];
      probe_segment_py[i_seg]             = tap.m_probe[i_probe].segmentPy[i_seg];
      probe_segment_pz[i_seg]             = tap.m_probe[i_probe].segmentPz[i_seg];
      probe_segment_chiSquared[i_seg]     = tap.m_probe[i_probe].segmentChiSquared[i_seg];
      probe_segment_numberDoF[i_seg]      = tap.m_probe[i_probe].segmentNumberDoF[i_seg];
      probe_segment_sector[i_seg]         = tap.m_probe[i_probe].segmentSector[i_seg];
      probe_segment_chamberIndex[i_seg]   = tap.m_probe[i_probe].segmentChamberIndex[i_seg];
      probe_segment_etaIndex[i_seg]       = tap.m_probe[i_probe].segmentEtaIndex[i_seg];
      probe_segment_nPrecisionHits[i_seg] = tap.m_probe[i_probe].segmentNPrecisionHits[i_seg];
      probe_segment_nPhiLayers[i_seg]     = tap.m_probe[i_probe].segmentNPhiLayers[i_seg];
      probe_segment_nTrigEtaLayers[i_seg] = tap.m_probe[i_probe].segmentNTrigEtaLayers[i_seg];
    }
    for( int i_trig = 0; i_trig < n_trig; i_trig++ ){
      trigname->push_back(tap.m_HLTtrigmesName[i_trig]);
      L1trigname->push_back(tap.m_L1trigmesName[i_trig]);
      L1_passbyEvent->push_back( tap.m_vL1objects_probe[i_probe][i_trig].isPassedByEvent_TDT );
      if((int)tap.m_passTrigmes.size() == n_trig)
	passtrig->push_back(tap.m_passTrigmes[i_trig]);
      //
      probe_L1_pass->push_back( tap.m_vL1objects_probe[i_probe][i_trig].isPassed );
      probe_L1_eta->push_back( tap.m_vL1objects_probe[i_probe][i_trig].eta );
      probe_L1_phi->push_back( tap.m_vL1objects_probe[i_probe][i_trig].phi );
      probe_L1_dR->push_back( tap.m_vL1objects_probe[i_probe][i_trig].dRl1 );
      probe_L1_thrValue->push_back( tap.m_vL1objects_probe[i_probe][i_trig].thrValue );
      probe_L1_roiNum->push_back( tap.m_vL1objects_probe[i_probe][i_trig].roiNum );
      probe_L1_thrNumber->push_back( tap.m_vL1objects_probe[i_probe][i_trig].thrNumber );
      probe_L1_isMoreCandInRoI->push_back( tap.m_vL1objects_probe[i_probe][i_trig].isMoreCandInRoI );
      //
      probe_SA_pass->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].isPassed );
      probe_SA_pt->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].pt );
      probe_SA_eta->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].eta );
      probe_SA_phi->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].phi );
      probe_SA_etaBE->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].etaBE );
      probe_SA_phiBE->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].phiBE );
      probe_SA_etaMS->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].etaMS );
      probe_SA_phiMS->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].phiMS );
      probe_SA_tgcPt->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].tgcPt );
      probe_SA_ptBarrelRadius->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].ptBarrelRadius );
      probe_SA_ptBarrelSagitta->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].ptBarrelSagitta );
      probe_SA_ptEndcapAlpha->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].ptEndcapAlpha );
      probe_SA_ptEndcapBeta->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].ptEndcapBeta );
      probe_SA_ptEndcapRadius->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].ptEndcapRadius );
      probe_SA_ptCSC->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].ptCSC );
      probe_SA_sAddress->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].sAddress );
      probe_SA_roiEta->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].roiEta );
      probe_SA_roiPhi->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].roiPhi );
      probe_SA_isRpcFailure->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].isRpcFailure );
      probe_SA_isTgcFailure->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].isTgcFailure );

      //the measured radious of the muon in one particular super point
      probe_SA_superPointR_BI->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointR_BI );
      probe_SA_superPointR_BM->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointR_BM );
      probe_SA_superPointR_BO->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointR_BO );
      probe_SA_superPointR_EI->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointR_EI );
      probe_SA_superPointR_EM->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointR_EM );
      probe_SA_superPointR_EO->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointR_EO );
      probe_SA_superPointR_EE->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointR_EE );
      probe_SA_superPointR_CSC->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointR_CSC );
      probe_SA_superPointR_BEE->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointR_BEE );
      probe_SA_superPointR_BME->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointR_BME );
      //the measured Z position of the muon in one particular super point
      probe_SA_superPointZ_BI->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointZ_BI );
      probe_SA_superPointZ_BM->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointZ_BM );
      probe_SA_superPointZ_BO->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointZ_BO );
      probe_SA_superPointZ_EI->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointZ_EI );
      probe_SA_superPointZ_EM->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointZ_EM );
      probe_SA_superPointZ_EO->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointZ_EO );
      probe_SA_superPointZ_EE->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointZ_EE );
      probe_SA_superPointZ_CSC->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointZ_CSC );
      probe_SA_superPointZ_BEE->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointZ_BEE );
      probe_SA_superPointZ_BME->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointZ_BME );
      //the measured slope of the muon in one particular super point
      probe_SA_superPointSlope_BI->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointSlope_BI );
      probe_SA_superPointSlope_BM->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointSlope_BM );
      probe_SA_superPointSlope_BO->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointSlope_BO );
      probe_SA_superPointSlope_EI->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointSlope_EI );
      probe_SA_superPointSlope_EM->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointSlope_EM );
      probe_SA_superPointSlope_EO->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointSlope_EO );
      probe_SA_superPointSlope_EE->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointSlope_EE );
      probe_SA_superPointSlope_CSC->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointSlope_CSC );
      probe_SA_superPointSlope_BEE->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointSlope_BEE );
      probe_SA_superPointSlope_BME->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointSlope_BME );
      //the measured intercept of the muon in one particular super point
      probe_SA_superPointIntercept_BI->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointIntercept_BI );
      probe_SA_superPointIntercept_BM->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointIntercept_BM );
      probe_SA_superPointIntercept_BO->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointIntercept_BO );
      probe_SA_superPointIntercept_EI->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointIntercept_EI );
      probe_SA_superPointIntercept_EM->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointIntercept_EM );
      probe_SA_superPointIntercept_EO->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointIntercept_EO );
      probe_SA_superPointIntercept_EE->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointIntercept_EE );
      probe_SA_superPointIntercept_CSC->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointIntercept_CSC );
      probe_SA_superPointIntercept_BEE->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointIntercept_BEE );
      probe_SA_superPointIntercept_BME->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointIntercept_BME );
      //the chi2 of the muon in one particular super point
      probe_SA_superPointChi2_BI->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointChi2_BI );
      probe_SA_superPointChi2_BM->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointChi2_BM );
      probe_SA_superPointChi2_BO->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointChi2_BO );
      probe_SA_superPointChi2_EI->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointChi2_EI );
      probe_SA_superPointChi2_EM->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointChi2_EM );
      probe_SA_superPointChi2_EO->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointChi2_EO );
      probe_SA_superPointChi2_EE->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointChi2_EE );
      probe_SA_superPointChi2_CSC->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointChi2_CSC );
      probe_SA_superPointChi2_BEE->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointChi2_BEE );
      probe_SA_superPointChi2_BME->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].superPointChi2_BME );
      //
      probe_SA_rpcHitX->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].rpcHitX );
      probe_SA_rpcHitY->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].rpcHitY );
      probe_SA_rpcHitZ->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].rpcHitZ );
      probe_SA_rpcHitR->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].rpcHitR );
      probe_SA_rpcHitEta->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].rpcHitEta );
      probe_SA_rpcHitPhi->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].rpcHitPhi );
      probe_SA_rpcHitMeasPhi->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].rpcHitMeasPhi );
      probe_SA_rpcHitLayer->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].rpcHitLayer );
      probe_SA_rpcHitStationName->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].rpcHitStationName );
      //
      probe_SA_mdtHitIsOutlier->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].mdtHitIsOutlier );
      probe_SA_mdtHitChamber->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].mdtHitChamber );
      probe_SA_mdtHitR->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].mdtHitR );
      probe_SA_mdtHitZ->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].mdtHitZ );
      probe_SA_mdtHitPhi->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].mdtHitPhi );
      probe_SA_mdtHitResidual->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].mdtHitResidual );

      probe_SA_roadAw->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].roadAw );
      probe_SA_roadBw->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].roadBw );
      probe_SA_zMin->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].zMin );
      probe_SA_zMax->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].zMax );
      probe_SA_rMin->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].rMin );
      probe_SA_rMax->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].rMax );
      probe_SA_etaMin->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].etaMin );
      probe_SA_etaMax->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].etaMax );
      //
      probe_SA_stgcClusterR->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].stgcClusterR );
      probe_SA_stgcClusterZ->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].stgcClusterZ );
      probe_SA_stgcClusterEta->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].stgcClusterEta );
      probe_SA_stgcClusterPhi->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].stgcClusterPhi );
      probe_SA_stgcClusterResidualR->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].stgcClusterResidualR );
      probe_SA_stgcClusterResidualPhi->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].stgcClusterResidualPhi );
      probe_SA_stgcClusterStationEta->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].stgcClusterStationEta );
      probe_SA_stgcClusterStationPhi->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].stgcClusterStationPhi );
      probe_SA_stgcClusterStationName->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].stgcClusterStationName );
      probe_SA_stgcClusterType->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].stgcClusterType );
      probe_SA_stgcClusterIsOutlier->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].stgcClusterIsOutlier );
      probe_SA_stgcClusterLayer->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].stgcClusterLayer );
      probe_SA_mmClusterR->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].mmClusterR );
      probe_SA_mmClusterZ->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].mmClusterZ );
      probe_SA_mmClusterEta->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].mmClusterEta );
      probe_SA_mmClusterPhi->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].mmClusterPhi );
      probe_SA_mmClusterResidualR->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].mmClusterResidualR );
      probe_SA_mmClusterResidualPhi->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].mmClusterResidualPhi );
      probe_SA_mmClusterStationEta->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].mmClusterStationEta );
      probe_SA_mmClusterStationPhi->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].mmClusterStationPhi );
      probe_SA_mmClusterStationName->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].mmClusterStationName );
      probe_SA_mmClusterIsOutlier->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].mmClusterIsOutlier );
      probe_SA_mmClusterLayer->push_back( tap.m_vSAobjects_probe[i_probe][i_trig].mmClusterLayer );
      //
      probe_CB_pass->push_back( tap.m_vCBobjects_probe[i_probe][i_trig].isPassed );
      probe_CB_pt->push_back( tap.m_vCBobjects_probe[i_probe][i_trig].pt );
      probe_CB_eta->push_back( tap.m_vCBobjects_probe[i_probe][i_trig].eta );
      probe_CB_phi->push_back( tap.m_vCBobjects_probe[i_probe][i_trig].phi );
      //
      probe_EF_pass->push_back( tap.m_vEFobjects_probe[i_probe][i_trig].isPassed );
      probe_EF_pt->push_back( tap.m_vEFobjects_probe[i_probe][i_trig].pt );
      probe_EF_eta->push_back( tap.m_vEFobjects_probe[i_probe][i_trig].eta );
      probe_EF_phi->push_back( tap.m_vEFobjects_probe[i_probe][i_trig].phi );
    }
    m_tree->Fill();
  }
}

int EventTreeMT::finalize() {
  // write into file and close	
  m_file->Write();	
  delete  trigname;
  delete  L1trigname;
  delete  L1_passbyEvent;
  delete  passtrig;
  delete  tag_L1_pass;
  delete  tag_L1_roiNum;
  delete  tag_SA_pass;
  delete  tag_CB_pass;
  delete  tag_EF_pass;
  delete  probe_L1_pass;  
  delete  probe_L1_eta; 
  delete  probe_L1_phi; 
  delete  probe_L1_dR; 
  delete  probe_L1_thrValue; 
  delete  probe_L1_roiNum; 
  delete  probe_L1_thrNumber; 
  delete  probe_L1_isMoreCandInRoI; 
  delete  probe_SA_pass; 
  delete  probe_SA_pt; 
  delete  probe_SA_eta; 
  delete  probe_SA_phi; 
  delete  probe_SA_etaMS; 
  delete  probe_SA_phiMS; 
  delete  probe_SA_etaBE; 
  delete  probe_SA_phiBE; 
  delete  probe_SA_tgcPt;
  delete  probe_SA_ptBarrelRadius;
  delete  probe_SA_ptBarrelSagitta;
  delete  probe_SA_ptEndcapAlpha;
  delete  probe_SA_ptEndcapBeta;
  delete  probe_SA_ptEndcapRadius;
  delete  probe_SA_ptCSC;
  delete  probe_SA_sAddress;
  delete  probe_SA_roiEta;
  delete  probe_SA_roiPhi;
  delete  probe_SA_isRpcFailure;
  delete  probe_SA_isTgcFailure;

  //the measured radious of the muon in one particular super point
  delete  probe_SA_superPointR_BI;
  delete  probe_SA_superPointR_BM;
  delete  probe_SA_superPointR_BO;
  delete  probe_SA_superPointR_EI;
  delete  probe_SA_superPointR_EM;
  delete  probe_SA_superPointR_EO;
  delete  probe_SA_superPointR_EE;
  delete  probe_SA_superPointR_CSC;
  delete  probe_SA_superPointR_BEE;
  delete  probe_SA_superPointR_BME;
  //the measured Z position of the muon in one particular super point
  delete  probe_SA_superPointZ_BI;
  delete  probe_SA_superPointZ_BM;
  delete  probe_SA_superPointZ_BO;
  delete  probe_SA_superPointZ_EI;
  delete  probe_SA_superPointZ_EM;
  delete  probe_SA_superPointZ_EO;
  delete  probe_SA_superPointZ_EE;
  delete  probe_SA_superPointZ_CSC;
  delete  probe_SA_superPointZ_BEE;
  delete  probe_SA_superPointZ_BME;
  //the measured slope of the muon in one particular super point
  delete  probe_SA_superPointSlope_BI;
  delete  probe_SA_superPointSlope_BM;
  delete  probe_SA_superPointSlope_BO;
  delete  probe_SA_superPointSlope_EI;
  delete  probe_SA_superPointSlope_EM;
  delete  probe_SA_superPointSlope_EO;
  delete  probe_SA_superPointSlope_EE;
  delete  probe_SA_superPointSlope_CSC;
  delete  probe_SA_superPointSlope_BEE;
  delete  probe_SA_superPointSlope_BME;
  //the measured intercept of the muon in one particular super point
  delete  probe_SA_superPointIntercept_BI;
  delete  probe_SA_superPointIntercept_BM;
  delete  probe_SA_superPointIntercept_BO;
  delete  probe_SA_superPointIntercept_EI;
  delete  probe_SA_superPointIntercept_EM;
  delete  probe_SA_superPointIntercept_EO;
  delete  probe_SA_superPointIntercept_EE;
  delete  probe_SA_superPointIntercept_CSC;
  delete  probe_SA_superPointIntercept_BEE;
  delete  probe_SA_superPointIntercept_BME;
  //the chi2 of the muon in one particular super point
  delete  probe_SA_superPointChi2_BI;
  delete  probe_SA_superPointChi2_BM;
  delete  probe_SA_superPointChi2_BO;
  delete  probe_SA_superPointChi2_EI;
  delete  probe_SA_superPointChi2_EM;
  delete  probe_SA_superPointChi2_EO;
  delete  probe_SA_superPointChi2_EE;
  delete  probe_SA_superPointChi2_CSC;
  delete  probe_SA_superPointChi2_BEE;
  delete  probe_SA_superPointChi2_BME;
  //
  delete  probe_SA_rpcHitX;
  delete  probe_SA_rpcHitY;
  delete  probe_SA_rpcHitZ;
  delete  probe_SA_rpcHitR;
  delete  probe_SA_rpcHitEta;
  delete  probe_SA_rpcHitPhi;
  delete  probe_SA_rpcHitMeasPhi;
  delete  probe_SA_rpcHitLayer;
  delete  probe_SA_rpcHitStationName;
  //
  delete  probe_SA_tgcHitZ;
  delete  probe_SA_tgcHitR;
  delete  probe_SA_tgcHitEta;
  delete  probe_SA_tgcHitPhi;
  delete  probe_SA_tgcHitWidth;
  delete  probe_SA_tgcHitStationNum;
  delete  probe_SA_tgcHitIsStrip;
  delete  probe_SA_tgcHitBCTag;
  delete  probe_SA_tgcHitInRoad;
  //
  delete  probe_SA_mdtHitIsOutlier;
  delete  probe_SA_mdtHitChamber;
  delete  probe_SA_mdtHitR;
  delete  probe_SA_mdtHitZ;
  delete  probe_SA_mdtHitPhi;
  delete  probe_SA_mdtHitResidual;

  delete  probe_SA_roadAw;
  delete  probe_SA_roadBw;
  delete  probe_SA_zMin;
  delete  probe_SA_zMax;
  delete  probe_SA_rMin;
  delete  probe_SA_rMax;
  delete  probe_SA_etaMin;
  delete  probe_SA_etaMax;
  delete  probe_SA_stgcClusterR;
  delete  probe_SA_stgcClusterZ;
  delete  probe_SA_stgcClusterEta;
  delete  probe_SA_stgcClusterPhi;
  delete  probe_SA_stgcClusterResidualR;
  delete  probe_SA_stgcClusterResidualPhi;
  delete  probe_SA_stgcClusterStationEta;
  delete  probe_SA_stgcClusterStationPhi;
  delete  probe_SA_stgcClusterStationName;
  delete  probe_SA_stgcClusterType;
  delete  probe_SA_stgcClusterIsOutlier;
  delete  probe_SA_stgcClusterLayer;
  delete  probe_SA_mmClusterR;
  delete  probe_SA_mmClusterZ;
  delete  probe_SA_mmClusterEta;
  delete  probe_SA_mmClusterPhi;
  delete  probe_SA_mmClusterResidualR;
  delete  probe_SA_mmClusterResidualPhi;
  delete  probe_SA_mmClusterStationEta;
  delete  probe_SA_mmClusterStationPhi;
  delete  probe_SA_mmClusterStationName;
  delete  probe_SA_mmClusterIsOutlier;
  delete  probe_SA_mmClusterLayer;
  delete  probe_CB_pass; 
  delete  probe_CB_pt; 
  delete  probe_CB_eta; 
  delete  probe_CB_phi; 
  delete  probe_EF_pass; 
  delete  probe_EF_pt; 
  delete  probe_EF_eta; 
  delete  probe_EF_phi;
  
  return 1;
}            
