#ifndef _INC_TAGANDPROBEMT
#define _INC_TAGANDPROBEMT

// stdlib
#include <iostream>
#include <fstream>
#include <vector>

// ASG Tools
#include "AsgTools/AsgMetadataTool.h"
#include "AsgTools/AsgTool.h"
#include "AsgTools/ToolHandle.h"
#include "AthenaBaseComps/AthMsgStreamMacros.h"
#include "AthenaBaseComps/AthMessaging.h"

//StoreGate
#include "StoreGate/ReadHandleKey.h"

// xAOD
#include "TrkExInterfaces/IExtrapolator.h"
#include "TrkVertexFitterInterfaces/IVertexFitter.h"
#include "TrigDecisionTool/TrigDecisionTool.h"
#include "xAODMuon/MuonContainer.h"
#include "xAODMuon/MuonSegmentContainer.h"
#include "xAODTracking/TrackParticleContainer.h"
#include "xAODTrigMuon/L2StandAloneMuon.h"
#include "xAODTrigMuon/L2StandAloneMuonContainer.h"
#include "xAODTrigMuon/L2CombinedMuon.h"
#include "xAODTrigMuon/L2CombinedMuonContainer.h"
#include "xAODTrigger/MuonRoIContainer.h"

// my classes
#include "CalcEfficiency/ParticleSelecterTool.h"
#include "CalcEfficiency/TrigMatchingToolMT.h"
#include "CalcEfficiency/MuonExtUtils.h"
#include "CalcEfficiency/Requirements.h"
#include "CalcEfficiency/OfflineObjects.h"
#include "CalcEfficiency/L1Objects.h"
#include "CalcEfficiency/SAObjects.h"
#include "CalcEfficiency/CBObjects.h"
#include "CalcEfficiency/EFObjects.h"
#include "CalcEfficiency/Utils.h"

class TagAndProbeMT {
	
  public:
    TagAndProbeMT();

    virtual int initialize( const int& message, 
                            const bool& useExt,
                            const TString method,
                            MuonExtUtils ext,
                            ToolHandle<Trig::TrigDecisionTool> tdt,
                            const std::string dataType ); //!
    virtual void clear(); //!
    virtual void addMesChain( const string& L1trig, const string& HLTtrig );
    virtual void checkMesChain();
    virtual bool isPassedTrigger();
    virtual bool setProbes( SG::ReadHandle<xAOD::MuonContainer> &muons, SG::ReadHandle<xAOD::MuonRoIContainer> &rois );
    virtual void doProbeMatching( SG::ReadHandle<xAOD::MuonRoIContainer> &rois );
    virtual void doTagMatching( SG::ReadHandle<xAOD::MuonRoIContainer> &rois );
    double dRl1bypt( double mupt ) const;
    bool passDimuSelection( const xAOD::Muon* tag, 
                            const xAOD::Muon* probe, 
                            int chain,
                            Requirement& req );
//    Bool_t getSATEName( const std::string& mesHLT, std::string& teName );
//    Bool_t getCBTEName( const std::string& mesHLT, std::string& teName );

  private:
    ToolHandle<Trig::TrigDecisionTool> m_trigDecTool{"Trig::TrigDecisionTool/TrigDecisionTool"}; //!
    //My tools
    MuonExtUtils m_ext;
    ParticleSelecterTool m_pst;
    dk::Utils m_utils;
    TrigMatchingToolMT m_matchToolMT;
  public:    
    std::vector< std::string > m_trigchainList; //required chain for tag muon
    std::vector< std::string > m_hltchainList; //required chain for tag muon in L1~CB, when trigchainList = JpsimumuFS
    std::vector< std::string > m_l1chainList; 
    std::vector< string > m_trigTagSATEName; //!
    std::vector< string > m_trigTagCBTEName; //!
    
    std::vector< std::string > m_L1trigmesName;
    std::vector< std::string > m_HLTtrigmesName; //the chain for probe muon
    std::vector< int > m_passTrigmes;
    std::vector< string > m_trigmesSATEname; //!
    std::vector< string > m_trigmesCBTEname; //!
    
    std::vector< double > m_massMin;
    std::vector< double > m_massMax;
    std::vector< bool > m_isPassed;
    int m_nChain;
    int m_nmesChain;
    float m_chi2cut{20};
    bool m_isnoL1{false};
    bool m_reqdR{true}; //still not implement --> always true
    bool m_ignoreMassreq{false};
    bool m_useExt{true};
    int m_message; //!
    bool m_run3evtSelection{true};
    static constexpr int SegmentMaxNumber{10};
    TString m_method; //!
    std::string dataType; //!
 
  public:   
    std::vector< std::pair< const xAOD::Muon*, const xAOD::Muon* > > m_tappairs; //!
    Requirements m_requirements_tag; //requirements tag and probe
    L1Objects m_L1objects_tag;
    SAObjects m_SAobjects_tag;
    CBObjects m_CBobjects_tag;
    EFObjects m_EFobjects_tag;
    OfflineObjects m_tag;
    Requirements m_requirements_probe; //requirements tag and probe
    std::vector<L1Objects> m_vL1objects_probe;
    std::vector<SAObjects> m_vSAobjects_probe;
    std::vector<CBObjects> m_vCBobjects_probe;
    std::vector<EFObjects> m_vEFobjects_probe;
    OfflineObjects m_probe;
    std::vector< int > m_tag_L1_pass;
    std::vector< int > m_tag_L1_roiNum;
    std::vector< int > m_tag_SA_pass;
    std::vector< int > m_tag_CB_pass;
    std::vector< int > m_tag_EF_pass;

};
#endif
