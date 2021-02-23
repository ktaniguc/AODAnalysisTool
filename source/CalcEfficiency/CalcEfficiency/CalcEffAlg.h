#ifndef CALCEFFALG_H // Include guard
#define CALCEFFALG_H

// Base
#include "AthenaBaseComps/AthAlgorithm.h"

class ITHistSvc;
class TH1D;
class TH2D;

// ASG Tools
#include "AsgTools/AsgMetadataTool.h"
#include "AsgTools/AsgTool.h"
#include "AsgTools/ToolHandle.h"

#include "xAODTruth/TruthParticleContainer.h"         
#include "xAODTruth/TruthVertexContainer.h"         
#include "xAODMuon/MuonContainer.h"
#include "xAODMuon/MuonSegmentContainer.h"
#include "xAODTrigMuon/L2StandAloneMuon.h"
#include "xAODTrigMuon/L2StandAloneMuonContainer.h"
#include "xAODTrigMuon/L2CombinedMuon.h"
#include "xAODTrigMuon/L2CombinedMuonContainer.h"
#include "xAODTrigger/MuonRoIContainer.h"
#include "xAODTracking/TrackParticleContainer.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TFile.h"

// Official Tools
#include "TrigConfxAOD/xAODConfigTool.h"
#include "TrigConfInterfaces/ITrigConfigSvc.h"
//#include "GoodRunsLists/IGoodRunsListSelectionTool.h"
#include "AsgAnalysisInterfaces/IGoodRunsListSelectionTool.h"
#include <TrigConfInterfaces/ITrigConfigTool.h>
#include "TrigDecisionTool/TrigDecisionTool.h"

#include "xAODTrigger/TrigCompositeContainer.h"
#include "xAODTrigger/TrigCompositeAuxContainer.h"
#include "xAODTrigger/MuonRoIContainer.h"

#include "StoreGate/ReadHandleKey.h"
#include "GaudiKernel/ITHistSvc.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODTrigger/versions/TrigComposite_v1.h"
#include "TrkExInterfaces/IExtrapolator.h"
            
// My class
#include "CalcEfficiency/MuonExtUtils.h"
#include "CalcEfficiency/HistNtupleMT.h"
#include "CalcEfficiency/TagAndProbe.h"
#include "CalcEfficiency/TagAndProbeMT.h"
#include "CalcEfficiency/EventTreeMT.h"
#include "CalcEfficiency/NtupleManager.h"
#include "CalcEfficiency/NtupleMakerTool.h"

using namespace TrigCompositeUtils;
using namespace SG;

// foaward
class GoodRunsListSelectionTool;
namespace CP {
  class MuonCalibrationAndSmearingTool; 
}
namespace TrigConf {
  class xAODConfigTool;
  class ITrigConfigTool;
}
namespace Trk { 
  class IExtrapolator; 
  class IVertexFitter; 
}


//class CalcEffAlg : public AthReentrantAlgorithm
class CalcEffAlg : public AthAlgorithm
{

  public:
    CalcEffAlg(const std::string& name, ISvcLocator* pSvcLocator); // Constructor
    enum L1Items{ NOTHING, L1_MU4, L1_MU6, L1_MU10, L1_MU11, L1_MU15, L1_MU20, L1_MU21, ERROR };
    int L1trigThr(std::string l1item) const {
      if( "L1_MU4"==l1item) return L1Items::L1_MU4;
      if( "L1_MU6"==l1item) return L1Items::L1_MU6;
      if( "L1_MU10"==l1item) return L1Items::L1_MU10;
      if( "L1_MU11"==l1item) return L1Items::L1_MU11;
      if( "L1_MU15"==l1item) return L1Items::L1_MU15;
      if( "L1_MU20"==l1item) return L1Items::L1_MU20;
      if( "L1_MU21"==l1item) return L1Items::L1_MU21;
      return L1Items::ERROR;   
    }
    int dimuL1trigThr(std::string l1item) const {
      if( "L1_2MU4"==l1item) return L1Items::L1_MU4;
      if( "L1_2MU6"==l1item) return L1Items::L1_MU6;
      if( "L1_2MU10"==l1item) return L1Items::L1_MU10;
      if( "L1_2MU11"==l1item) return L1Items::L1_MU11;
      if( "L1_2MU15"==l1item) return L1Items::L1_MU15;
      if( "L1_2MU20"==l1item) return L1Items::L1_MU20;
      if( "L1_2MU21"==l1item) return L1Items::L1_MU21;
      return L1Items::ERROR;   
    }
    enum TapType{ L1, L2, EF, ALL };

    StatusCode initialize();
    StatusCode finalize();
    StatusCode execute();

  private:
    int m_message;
    std::string m_etname;
    std::string m_dataType;
    bool m_isFirstEvent; //!
    HistNtupleMT m_histsMT;
    EventTreeMT m_etMT;
    NtupleManager m_ntuple;
    NtupleMakerTool m_ntupMakerTool;
    
    TagAndProbe m_tap;
    TagAndProbeMT m_tapMT;
    // GRL tool
    ToolHandle<IGoodRunsListSelectionTool> m_grlTool{"GoodRunsListSelectionTool/MyGRLTool"};
    bool m_useGRL; //!
    bool m_run3{false}; //if getNavigationFormat is TrigComposite, then this flag turn true
    bool m_afterIOmode{true}; // check the version if IOmode is implemented ?

    // The decision tool
    ToolHandle<Trig::TrigDecisionTool> m_trigDecTool{"Trig::TrigDecisionTool/TrigDecisionTool"}; //!
    ServiceHandle<TrigConf::ITrigConfigSvc> m_configSvc;
    SG::ReadHandleKey<xAOD::EventInfo> m_EventInfoKey{this, "xAODEventInfo", "EventInfo", "Name of EventInfo object"};
    SG::ReadHandleKey<xAOD::TruthParticleContainer> m_truthParticleContainerKey{this, "xAODTruthParticleContainer", "TruthParticles", "Name of TruthParticleContainer object"};
    SG::ReadHandleKey<xAOD::MuonContainer> m_MuonContainerKey{this, "xAODMuonContainer", "Muons", "Name of MuonContainer object"};
    SG::ReadHandleKey<xAOD::MuonRoIContainer> m_MuonRoIContainerKey{this, "xAODMuonRoIContainer", "LVL1MuonRoIs", "Name of MuonRoIContainer object"};
    SG::ReadHandleKey<xAOD::L2StandAloneMuonContainer> m_L2SAKey{this, "xAODL2SAContainer", "HLT_MuonL2SAInfo", "Name of L2SAContainer object"};
    SG::ReadHandleKey<xAOD::L2CombinedMuonContainer> m_L2CBKey{this, "xAODL2CBContainer", "HLT_MuonL2CBInfo", "Name of L2CBContainer object"};
    SG::ReadHandleKey<xAOD::L2StandAloneMuonContainer> m_L2SAIOKey{this, "xAODL2SAContainer_IO", "HLT_MuonL2SAInfoIOmode", "Name of L2SAIOContainer object"};
    SG::ReadHandleKey<xAOD::L2CombinedMuonContainer> m_L2CBIOKey{this, "xAODL2CBContainer_IO", "HLT_MuonL2CBInfoIOmode", "Name of L2CBIOContainer object"};
    std::string m_trigDecisionKey; //!< SG key of the trigger data (TrigDecision object) add 0602

    // Extrapolator Class
    ToolHandle<Trk::IExtrapolator> m_extrapolator{"Trk::Extrapolator/AtlasExtrapolator"};
    MuonExtUtils m_ext; //!
    std::string m_tapmethod; //!
    bool m_doNSWMon{false};
    bool m_makeNtuple{false};
    bool m_applyMuonVeto{false};
    std::string m_monTrigName{"HLT_mu6_L1MU6"};
    bool m_isAsymNSW{false};
    bool m_useExt; //!
};

#endif // CALCALG_H
