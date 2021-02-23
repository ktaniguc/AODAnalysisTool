#ifndef _INC_TRIGMATCHINGTOOLMT
#define _INC_TRIGMATCHINGTOOLMT

#include "CalcEfficiency/OfflineObjects.h"
#include "CalcEfficiency/Requirements.h"
#include "CalcEfficiency/L1Objects.h"
#include "CalcEfficiency/SAObjects.h"
#include "CalcEfficiency/CBObjects.h"
#include "CalcEfficiency/EFObjects.h"
#include "CalcEfficiency/MuonExtUtils.h"
#include "CalcEfficiency/Utils.h"
// stdlib
#include <iostream>
#include <fstream>
#include <vector>

#include "TrigDecisionTool/TrigDecisionTool.h"
#include "TrkExInterfaces/IExtrapolator.h"
#include "xAODMuon/MuonContainer.h"
#include "xAODMuon/MuonSegmentContainer.h"
#include "xAODTracking/TrackParticleContainer.h"
#include "xAODTrigMuon/L2StandAloneMuon.h"
#include "xAODTrigMuon/L2StandAloneMuonContainer.h"
#include "xAODTrigMuon/L2CombinedMuon.h"
#include "xAODTrigMuon/L2CombinedMuonContainer.h"
#include "xAODTrigger/MuonRoIContainer.h"
#include "StoreGate/ReadHandle.h"
#include "TrigCompositeUtils/TrigCompositeUtils.h"

class TrigMatchingToolMT {
  public:
    TrigMatchingToolMT();

    void initialize(ToolHandle<Trig::TrigDecisionTool> tdt, MuonExtUtils ext)
    { 
      m_trigDecTool = tdt;
      m_ext = ext;
    };

    enum L1Items{ L1_MU4=1, L1_MU6, L1_MU10, L1_MU11, L1_MU20, L1_MU21, ERROR};
    int L1trigThr(std::string l1item) const {
      if( "L1_MU4"==l1item) return L1Items::L1_MU4;
      if( "L1_MU6"==l1item) return L1Items::L1_MU6;
      if( "L1_MU10"==l1item) return L1Items::L1_MU10;
      if( "L1_MU11"==l1item) return L1Items::L1_MU11;
      if( "L1_MU20"==l1item) return L1Items::L1_MU20;
      if( "L1_MU21"==l1item) return L1Items::L1_MU21;
      return L1Items::ERROR;   
    }
    int dimuL1trigThr(std::string l1item) const {
      if( "L1_2MU4"==l1item) return L1Items::L1_MU4;
      if( "L1_2MU6"==l1item) return L1Items::L1_MU6;
      if( "L1_2MU10"==l1item) return L1Items::L1_MU10;
      if( "L1_2MU11"==l1item) return L1Items::L1_MU11;
      if( "L1_2MU20"==l1item) return L1Items::L1_MU20;
      if( "L1_2MU21"==l1item) return L1Items::L1_MU21;
      return L1Items::ERROR;   
    }

    bool matchL1( SG::ReadHandle<xAOD::MuonRoIContainer> &murois,
                  const xAOD::Muon* muon, 
                  L1Object& l1obj,
                  std::string trig,
                  std::string HLTtrig );
    bool matchSA( std::string& trig,
                  const L1Object l1obj,
                  SAObject& saobj );
    bool matchCB( std::string& trig,
                  const SAObject saobj,
                  CBObject& cbobj );
    bool matchEF( std::string& trig, 
                  const xAOD::Muon* muon,
                  EFObject& efobj );
//    bool matchEFFS( const Trig::FeatureContainer& fc,
//                    const xAOD::Muon* muon,
//                    EFObject& efobj );
//    //following match** function is for probe's matching
//    bool matchL1( const Trig::FeatureContainer& fc, 
//                  const xAOD::Muon* muon,
//                  L1Object& l1obj ); 
//    bool matchSA( const Trig::FeatureContainer& fc,
//                  std::string& mesSATEName,
//                  const L1Object l1obj,
//                  SAObject& saobj);
//    bool matchCB( const Trig::FeatureContainer& fc,
//                  std::string& mesCBTEName,
//                  const L1Object l1obj,
//                  CBObject& cbobj);
//    bool matchEF( const Trig::FeatureContainer& fc,
//                  const xAOD::Muon* muon,
//                  EFObject& efobj );

  private:
    MuonExtUtils m_ext;
    dk::Utils m_utils;
    ToolHandle<Trig::TrigDecisionTool> m_trigDecTool{"Trig::TrigDecisionTool/TrigDecisionTool"}; //!


};

#endif
