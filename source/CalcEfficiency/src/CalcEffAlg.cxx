// Gaudi
#include "GaudiKernel/IAlgTool.h"

#include "GaudiKernel/ITHistSvc.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>

// xAOD 
#include "EventInfo/EventInfo.h"
#include <EventInfo/EventID.h>
#include "xAODEventInfo/EventInfo.h"
#include "xAODMuon/MuonContainer.h"
#include "xAODMuon/MuonSegmentContainer.h"
#include "xAODTrigMuon/L2StandAloneMuon.h"
#include "xAODTrigMuon/L2StandAloneMuonContainer.h"
#include "xAODTrigMuon/L2CombinedMuon.h"
#include "xAODTrigMuon/L2CombinedMuonContainer.h"
#include "xAODTrigger/MuonRoIContainer.h"
#include "xAODTracking/TrackParticleContainer.h"
#include "xAODTracking/VertexContainer.h"
#include "GoodRunsLists/IGoodRunsListSelectorTool.h"
#include "TrigConfxAOD/xAODConfigTool.h"
#include "TrigConfInterfaces/ITrigConfigTool.h"
#include "TrkParameters/TrackParameters.h"
#include "TrkExInterfaces/IExtrapolator.h"
#include "TrkVKalVrtFitter/TrkVKalVrtFitter.h"
#include "TrkVertexFitterInterfaces/IVertexFitter.h"

#include "TrigConfHLTData/HLTTriggerElement.h"

#include "TrigConfHLTData/HLTUtils.h"
#include "TrigConfHLTData/HLTTriggerElement.h"


#include "TrigConfL1Data/TriggerThreshold.h"
#include "TrigConfL1Data/CTPConfig.h"

//TriggerDecisionTool
#include "TrigDecisionTool/ChainGroup.h"
#include "TrigDecisionTool/Feature.h"
#include "TrigDecisionTool/FeatureContainer.h"

//MyTool
#include "CalcEfficiency/CalcEffAlg.h"

using namespace SG;
using namespace std;
using namespace TrigCompositeUtils;


//  m_trigConfTool("TrigConf::xAODConfigTool/xAODConfigTool"),
///THis function is constructor.
CalcEffAlg::CalcEffAlg(const std::string& name, ISvcLocator* pSvcLocator)
  : AthAlgorithm(name, pSvcLocator),
  m_configSvc( "TrigConf::TrigConfigSvc/TrigConfigSvc", name )
{
  ///(https://twiki.cern.ch/twiki/bin/view/Sandbox/WritingYourOwnAthenaAlgorithms)Notice the calls of the declareProperty("electron_Et_min_cut", m_electron_Et_min_cut = 20*GeV) in the constructor. This makes the C++ m_electron_Et_min_cut variable configurable from the Python job options. The first argument is a string containing the Python name of the variable. The second is the C++ variable with an optional equals-sign followed by a default value. Configuration from Python will be explained more later when we get to the job options.
  declareProperty( "message", m_message);
  declareProperty( "OutputFile", m_etname );
  declareProperty( "TapMethod", m_tapmethod );
  declareProperty( "doNSWMon", m_doNSWMon );
  declareProperty( "monTrigName", m_monTrigName );
  declareProperty( "isAsymNSW", m_isAsymNSW );
  declareProperty( "makeSimpleNtuple", m_makeNtuple );
  declareProperty( "applyMuonVeto", m_applyMuonVeto );
  declareProperty( "Extrapolate", m_useExt );
  declareProperty( "GRL", m_useGRL );
  declareProperty( "DataType", m_dataType );
}

StatusCode CalcEffAlg::initialize() {
  ATH_MSG_INFO("initialize()");
  ATH_MSG_INFO("Tag-and-probe method: " << m_tapmethod);
  ATH_MSG_INFO("do NSW montoring: " << m_doNSWMon << " with " << m_monTrigName);
  ATH_MSG_INFO("create simple ntuple file?: " << m_makeNtuple);
  ATH_MSG_INFO("apply muon associated with c, b, tau veto : " << m_applyMuonVeto);
  //==============================================================
  //==  GRL Tools
  //==============================================================
  ATH_CHECK( m_grlTool.retrieve() );
  //==============================================================
  //==  Trigger Tools
  //==============================================================
  ATH_CHECK(m_configSvc.retrieve());
  ATH_CHECK(m_trigDecTool.retrieve());
  if(m_trigDecTool->getNavigationFormat() == "TrigComposite") m_run3 = true;
  ATH_MSG_INFO("TrigDecTool->getNavigationFormat() == " << m_trigDecTool->getNavigationFormat());
  m_trigDecTool->ExperimentalAndExpertMethods()->enable();

  ////==============================================================
  ////==  Tools
  ////==============================================================
  ATH_CHECK(m_extrapolator.retrieve());
  //ATH_CHECK(m_vrtfitter.retrieve());

  ////==============================================================
  ////==  MuonExtrapolatorUtils Class
  ////==============================================================
  m_ext.initialize( m_extrapolator );
  m_isFirstEvent = true; 
  //==============================================================
  //==  TagAndProbe Class
  if(m_run3){  // Run3
    m_tapMT.initialize(m_message, m_useExt, m_tapmethod, m_ext, m_trigDecTool, m_dataType);
    m_tapMT.addMesChain( "L1_MU4", "HLT_mu4_l2io_L1MU4" );
    m_tapMT.addMesChain( "L1_MU4", "HLT_mu4_L1MU4" );
    m_tapMT.addMesChain( "L1_2MU10", "HLT_2mu14_l2io_L12MU10" );
    m_tapMT.addMesChain( "L1_2MU10", "HLT_2mu14_L12MU10" );
    m_tapMT.addMesChain( "L1_2MU6", "HLT_2mu6_l2io_L12MU6" );
    m_tapMT.addMesChain( "L1_2MU6", "HLT_2mu6_L12MU6" );
    m_tapMT.addMesChain( "L1_MU6", "HLT_mu6_L1MU6" );
    m_tapMT.addMesChain( "L1_MU10", "HLT_mu14_L1MU10" );
    m_tapMT.addMesChain( "L1_MU20", "HLT_mu26_ivarmedium_L1MU20" );
    if(m_doNSWMon){
      m_histsMT.initialize( "NSW-monitoring.root", m_tapMT, m_isAsymNSW, m_monTrigName );
    }
    if(m_makeNtuple){
      m_ntupMakerTool.initialize(m_trigDecTool, m_ext);
      m_ntupMakerTool.setTriggerName("L1_MU4", "HLT_mu4_l2io_L1MU4");
      m_ntupMakerTool.setTriggerName("L1_2MU6", "HLT_2mu6_l2io_L12MU6");
      m_ntupMakerTool.setTriggerName("L1_2MU10", "HLT_2mu14_l2io_L12MU10");
      m_ntupMakerTool.setTriggerName("L1_MU4", "HLT_mu4_L1MU4");
      m_ntupMakerTool.setTriggerName("L1_MU20", "HLT_mu24_L1MU20");
      m_ntupMakerTool.setTriggerName("L1_MU10", "HLT_mu14_L1MU10");
      m_ntupMakerTool.setTriggerName("L1_2MU6", "HLT_2mu6_L12MU6");
      std::cout << "--> make simple ntuple. TDT input chains : ";
      for(int i=0; i< (int)m_ntupMakerTool.m_hltchainName.size(); i++){
        std::cout << "[ " <<  m_ntupMakerTool.m_hltchainName.at(i) << " ], ";
      }
      std::cout << std::endl;
    }
  } else {     // Run2
    m_tap.initialize(m_message, m_useExt, m_tapmethod, m_ext, m_trigDecTool, m_dataType);
    m_tap.addMesChain( "L1_MU4", "HLT_mu4" );
    m_tap.addMesChain( "L1_MU6", "HLT_mu6" );
    m_tap.addMesChain( "L1_MU10", "HLT_mu10" );
    m_tap.addMesChain( "L1_MU10", "HLT_mu14" );
    m_tap.addMesChain( "L1_MU20", "HLT_mu26_ivarmedium" );
    // please initalize HistNtuple after addMesChain finished
//    m_histsMT.initialize( "efficiency-monitoring.root", m_tap );
  }

  ////==============================================================
  ////==  VrtFitterUtils Class
  ////==============================================================
  //m_vft.initialize( m_vrtfitter );

  //==============================================================
  //==  EventTree Class
  //==============================================================
  m_etMT.initialize( m_etname );
  if(m_makeNtuple){
    m_ntuple.initialize( m_etname ); 
  }

  ATH_CHECK( m_EventInfoKey.initialize() );
  ATH_CHECK( m_truthParticleContainerKey.initialize() );
  ATH_CHECK( m_MuonContainerKey.initialize() );
  ATH_CHECK( m_MuonRoIContainerKey.initialize() );
  ATH_CHECK( m_L2SAKey.initialize() );
  ATH_CHECK( m_L2CBKey.initialize() );
  ATH_CHECK( m_L2SAIOKey.initialize() );
  ATH_CHECK( m_L2CBIOKey.initialize() );
  
  ATH_MSG_INFO("CalcEffAlg::initialize() end");

  return StatusCode::SUCCESS;
}

StatusCode CalcEffAlg::finalize() {
  ATH_MSG_INFO("finalize()");
  if(m_run3 && m_doNSWMon){
    m_histsMT.finalize( m_tapMT );
  } //else {
//    m_histsMT.finalize( m_tap );
//  }
  ATH_MSG_INFO("HistNtuple successfully finished");
  m_etMT.finalize();
  ATH_MSG_INFO("EventTree successfully finished");
  if(m_run3 && m_makeNtuple){
    m_ntuple.finalize();
    ATH_MSG_INFO("simple ntuple file was successfully created");
  }
  std::cout << "CalcEffAlg successfully finished" << std::endl;
  return StatusCode::SUCCESS;
}

StatusCode CalcEffAlg::execute() {
  ATH_MSG_INFO("execute()");
  if(m_run3) m_tapMT.clear();
  else       m_tap.clear();

  const EventContext& ctx = getContext();
  ATH_MSG_INFO("Get event context << " << ctx );
  ///main function to do "Tag and Probe" .
  //==============================================================
  //=  Event information
  //==============================================================
  SG::ReadHandle<xAOD::EventInfo> eventInfo(m_EventInfoKey, ctx);
  uint32_t runNumber = eventInfo->runNumber();
  unsigned long long eventNumber = eventInfo->eventNumber();
  int lumiBlock = eventInfo-> lumiBlock();
  double averageInteractionsPerCrossing =  eventInfo -> averageInteractionsPerCrossing();
  ATH_MSG_INFO("===================================================");
  ATH_MSG_INFO("eventNumber==========#" << eventNumber << "========" );
  ATH_MSG_INFO("===================================================");
  ATH_MSG_INFO("Run = " << runNumber << " : Event = " << eventNumber << " : mu = " << averageInteractionsPerCrossing );
  bool isMC = true;
  if(!eventInfo->eventType(xAOD::EventInfo::IS_SIMULATION ) ){
    isMC = false;
  }
  // GRL
  if( !isMC && m_useGRL ){
    ATH_MSG_INFO("Skip this event via GRL");
    if(!m_grlTool->passRunLB(*eventInfo)) return StatusCode::SUCCESS; //checks the GRL and skips to next event if not passing
  } // end if not MC
  
  // do event cleaning if not MC
  if( !isMC ){
    bool isClean = true;
    if( eventInfo->errorState(xAOD::EventInfo::LAr) == xAOD::EventInfo::Error ) isClean = false;
    if( eventInfo->errorState(xAOD::EventInfo::Tile) == xAOD::EventInfo::Error ) isClean = false;
    if( eventInfo->errorState(xAOD::EventInfo::SCT) == xAOD::EventInfo::Error ) isClean = false;
    if( eventInfo->isEventFlagBitSet(xAOD::EventInfo::Core, 18) ) isClean = false;

    if(!isClean) return StatusCode::SUCCESS;
  }
  // ==== truth muon =====
  //pdgId : +-15 = tau, XX4XX = c meson, X4XXX = c baryon, XX5XX = b meson, X5XXX = b baryon
  if(m_applyMuonVeto){
    ATH_MSG_INFO("# check existance of muons associated with b, c, tau");
    SG::ReadHandle<xAOD::TruthParticleContainer> mctruths( m_truthParticleContainerKey, ctx);
    bool fromTau = false;
    bool fromC = false;
    bool fromB = false;
    for (auto mctruth : *mctruths) {
      if(!mctruth) continue;
      float mctruth_pt = mctruth->pt();
      if(mctruth_pt<0.00001) continue;

      if(std::abs(mctruth->pdgId()) == 13) {
        ATH_MSG_INFO(" MC truth muon barcode" << mctruth->barcode() << " pdgId/status/pt = " << mctruth->pdgId() << "/" << mctruth->status() << "/" << mctruth_pt);

        bool hasProdVtx = mctruth->hasProdVtx();
        if(hasProdVtx) {
          const xAOD::TruthVertex* vtx = mctruth->prodVtx();
          for(unsigned int i=0; i<vtx->nIncomingParticles(); i++) {
            const xAOD::TruthParticle* parent = vtx->incomingParticle(i);
            if(!parent) continue;
            float parent_pt = parent->pt();
            int parentId = std::abs(parent->pdgId());
            if(parent_pt<0.00001) continue;
            ATH_MSG_INFO("... parent of " << mctruth->barcode() << ", barcode" << parent->barcode() << " pdgId/status/pt = " << parent->pdgId() << "/" << parent->status() << "/" << parent_pt);
            if(parentId == 15) fromTau = true;
            int parId[2] = {0, 0};
            parId[0] = (parentId % 1000);
            parId[0] = (parId[0]/100);
            parId[1] = (parentId % 10000);
            parId[1] = (parId[1]/1000);
            if(parId[0] == 4 || parId[1] == 4) fromC = true;
            if(parId[0] == 5 || parId[1] == 5) fromB = true;
          }
        }
      }
    }
    if(fromTau || fromB || fromC){
      ATH_MSG_DEBUG("---> this event contains muons assiciated with b/c/tau = " << fromB << "/" << fromC << "/" << fromTau << ", skip");
      return StatusCode::SUCCESS;
    }
  }

  // ===== retrieve offline muons and lvl1 rois
  SG::ReadHandle<xAOD::MuonContainer> muons( m_MuonContainerKey, ctx );
  if( !muons.isValid() ){
    ATH_MSG_DEBUG("No valid MuonContainer with tag : " << m_MuonContainerKey);
    return StatusCode::SUCCESS;
  }
  SG::ReadHandle<xAOD::MuonRoIContainer> rois( m_MuonRoIContainerKey, ctx );
  if( !rois.isValid() ){
    ATH_MSG_DEBUG("No valid MuonRoIContainer with tag : " << m_MuonRoIContainerKey);
    return StatusCode::SUCCESS;
  }
  SG::ReadHandle<xAOD::L2StandAloneMuonContainer> l2sa( m_L2SAKey, ctx );
  if( !l2sa.isValid() ){
    ATH_MSG_DEBUG("No valid L2SA with tag : " << m_L2SAKey);
  }
  SG::ReadHandle<xAOD::L2CombinedMuonContainer> l2cb( m_L2CBKey, ctx );
  if( !l2cb.isValid() ){
    ATH_MSG_DEBUG("No valid L2CB with tag : " << m_L2CBKey);
  }
  SG::ReadHandle<xAOD::L2StandAloneMuonContainer> l2sa_io( m_L2SAIOKey, ctx );
  if( !l2sa_io.isValid() ){
    ATH_MSG_DEBUG("No valid L2SA with tag : " << m_L2SAIOKey);
    if(m_makeNtuple)m_ntupMakerTool.m_isvalidIOmode = false;
  }
  SG::ReadHandle<xAOD::L2CombinedMuonContainer> l2cb_io( m_L2CBIOKey, ctx );
  if( !l2cb_io.isValid() ){
    ATH_MSG_DEBUG("No valid L2CB with tag : " << m_L2CBIOKey);
    if(m_makeNtuple)m_ntupMakerTool.m_isvalidIOmode = false;
  }
  // ==== check run2API or run3API via trigDecTol
  if(m_trigDecTool->getNavigationFormat() == "TriggerElement"){ //Run2 TriggerElement
    if(m_isFirstEvent){
      auto cgs = m_trigDecTool->getChainGroup("HLT.*mu.*|L1_.*MU.*|HLT_noalg_L1.*MU.*");
      for( auto &trig : cgs->getListOfTriggers() ){
        auto cg = m_trigDecTool->getChainGroup(trig);
        bool isPassedCurrent = cg->isPassed();
        ATH_MSG_INFO(trig << " is passed ==>" << isPassedCurrent);
      }
    }
    m_isFirstEvent = false;
    if(!m_tap.isPassedTrigger()) {
      ATH_MSG_DEBUG("trigger not passeed");
      return StatusCode::SUCCESS;
    }
    if(!m_tap.setProbes( *muons )){
      ATH_MSG_DEBUG("could not set tag and probes");
      return StatusCode::SUCCESS;
    }
    m_tap.doProbeMatching( rois );
    m_etMT.filltree( m_tap, eventNumber, runNumber, lumiBlock, averageInteractionsPerCrossing );
//    m_histsMT.FillHist( m_tap );
  
  } else { // Run3 TrigComposite

    auto cgs = m_trigDecTool->getChainGroup("HLT.*mu.*|L1_.*MU.*|HLT_noalg_L1.*MU.*");

    if(m_isFirstEvent){
      for( auto &trig : cgs->getListOfTriggers() ){
        auto cg = m_trigDecTool->getChainGroup(trig);
        bool isPassedCurrent = cg->isPassed( TrigDefs::eventAccepted );
        ATH_MSG_INFO(trig << " is passed ===> " << isPassedCurrent);
      }
      m_isFirstEvent = false;
    }

    if(m_makeNtuple){
      m_ntupMakerTool.Clear();
      m_ntupMakerTool.withdrawInfo(muons, rois, l2sa, l2sa_io, l2cb, l2cb_io);
      m_ntuple.filltree(eventNumber, runNumber, m_ntupMakerTool);
    }
    //start tag-and-probe
    m_tapMT.checkMesChain();
    if(!m_tapMT.isPassedTrigger()) {
      ATH_MSG_DEBUG("trigger not passeed");
      return StatusCode::SUCCESS;
    }
    if(!m_tapMT.setProbes( muons, rois )){
      ATH_MSG_DEBUG("could not set tag and probes");
      return StatusCode::SUCCESS;
    }
    m_tapMT.doProbeMatching( rois );
    m_tapMT.doTagMatching( rois );
    m_etMT.filltree( m_tapMT, eventNumber, runNumber, lumiBlock, averageInteractionsPerCrossing );
    if(m_doNSWMon) m_histsMT.FillHist( m_tapMT );

  } // TrigComposite end

  return StatusCode::SUCCESS;
}

