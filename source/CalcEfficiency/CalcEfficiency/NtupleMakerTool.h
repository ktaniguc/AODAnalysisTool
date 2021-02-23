#ifndef _INC_NTUPLEMAKERTOOL
#define _INC_NTUPLEMAKERTOOL

#include "CalcEfficiency/OfflineObjects.h"
#include "CalcEfficiency/Requirements.h"
#include "CalcEfficiency/L1Objects.h"
#include "CalcEfficiency/SAObjects.h"
#include "CalcEfficiency/CBObjects.h"
#include "CalcEfficiency/EFObjects.h"
#include "CalcEfficiency/TDTObjects.h"
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

class NtupleMakerTool {
  public:
    NtupleMakerTool();

    void initialize(ToolHandle<Trig::TrigDecisionTool> tdt, MuonExtUtils ext)
    { 
      m_trigDecTool = tdt;
      m_ext = ext;
    };
    void setTriggerName(std::string l1, std::string hlt);

    void withdrawInfo(SG::ReadHandle<xAOD::MuonContainer>       &muons,
                      SG::ReadHandle<xAOD::MuonRoIContainer>    &rois,
                      SG::ReadHandle<xAOD::L2StandAloneMuonContainer> &l2sa,
                      SG::ReadHandle<xAOD::L2StandAloneMuonContainer> &l2saio,
                      SG::ReadHandle<xAOD::L2CombinedMuonContainer> &l2cb,
                      SG::ReadHandle<xAOD::L2CombinedMuonContainer> &l2cbio);
    void Clear(){
      m_muons.clear();
      m_TDTobjects.clear();
      m_L1objects.clear();
      m_SAobjects.clear();
      m_SAIOobjects.clear();
      m_CBobjects.clear();
      m_CBIOobjects.clear();
      m_EFobjects.clear();
    }
    std::vector<std::string> m_hltchainName;
    std::vector<std::string> m_l1chainName;

    OfflineObjects m_muons;
    TDTObjects m_TDTobjects;
    L1Objects m_L1objects;
    SAObjects m_SAobjects;
    SAObjects m_SAIOobjects;
    CBObjects m_CBobjects;
    CBObjects m_CBIOobjects;
    SAObjects m_EFobjects;
    bool m_isvalidIOmode{true};
  private:
    MuonExtUtils m_ext;
    dk::Utils m_utils;
    ToolHandle<Trig::TrigDecisionTool> m_trigDecTool{"Trig::TrigDecisionTool/TrigDecisionTool"}; //!

};

#endif
