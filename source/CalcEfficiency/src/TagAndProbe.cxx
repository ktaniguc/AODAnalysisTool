#include <iostream>
#include <fstream>
#include <vector>

#include "TrigConfxAOD/xAODConfigTool.h"
#include "xAODTrigMuon/L2StandAloneMuon.h"
#include "xAODTrigMuon/L2StandAloneMuonContainer.h"
#include "xAODTrigMuon/L2CombinedMuon.h"
#include "xAODTrigMuon/L2CombinedMuonContainer.h"
#include "xAODTrigger/MuonRoIContainer.h"
#include "xAODTracking/TrackParticleContainer.h"
#include "AsgTools/AsgTool.h"

#include "CalcEfficiency/TagAndProbe.h"

#include "TMath.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TVector3.h"

TagAndProbe::TagAndProbe(){ }

int TagAndProbe::initialize( const int& message, 
                             const bool& useExt,
                             const TString method,
                             MuonExtUtils ext,
                             ToolHandle<Trig::TrigDecisionTool> tdt,
                             const std::string dataType
                           )
{
/// This function depends on message, useExt, method, ext vft, tdt, dataType.
/// vecTrigEvent, m_trigEvent, ... is push_backed for each method.
/// m_tapType is ALL, L2 or EFF. 
  std::cout << "== tag and probe execute" << std::endl;
  // m_message : 0 -> INFO, 1 -> DEBUG, 2 -> WARNING
  m_message     = message;
  m_method      = method;
  m_trigDecTool = tdt;
  m_ext         = ext;
  m_useExt      = useExt; //still not implement --> always use
  this -> dataType = dataType;
  m_chi2cut = 20.; //still not implement
  m_tmt.initialize(tdt, ext);

  if(m_method=="Jpsi"){
//    m_l1chainList.push_back("L1_MU6");
    m_hltchainList.push_back("HLT_mu20");
//    m_trigchainList.push_back("HLT_mu6");
    m_trigchainList.push_back("HLT_mu20_2mu2noL1_JpsimumuFS");
    m_isPassed.push_back(false);
    m_reqdR = true;
    m_isnoL1 = true;
    m_massMin.push_back(3099.-300.);
    m_massMax.push_back(3099.+300.);
    std::string tsSATE = "";
    getSATEName( "HLT_mu20", tsSATE );
    m_trigTagSATEName.push_back( tsSATE );
    std::string tsCBTE = "";
    getCBTEName( "HLT_mu20", tsCBTE );
    m_trigTagCBTEName.push_back( tsCBTE );
  }
  else if(m_method=="JpsiMC"){
//    m_l1chainList.push_back("L1_MU6");
    m_hltchainList.push_back("HLT_mu6");
//    m_trigchainList.push_back("HLT_mu6");
    m_trigchainList.push_back("HLT_mu6");
    m_isPassed.push_back(false);
    m_reqdR = true;
    m_isnoL1 = true;
    m_massMin.push_back(3099.-300.);
    m_massMax.push_back(3099.+300.);
    std::string tsSATE = "";
    getSATEName( "HLT_mu6", tsSATE );
    m_trigTagSATEName.push_back( tsSATE );
    std::string tsCBTE = "";
    getCBTEName( "HLT_mu6", tsCBTE );
    m_trigTagCBTEName.push_back( tsCBTE );
  } 
  else if(m_method=="Ztap"){
    m_hltchainList.push_back("HLT_mu26_ivarmedium");
    m_trigchainList.push_back("HLT_mu26_ivarmedium");
    m_isPassed.push_back(false);
    m_reqdR = true;
    m_massMin.push_back(90000.-10000.);
    m_massMax.push_back(90000.+10000.);
    std::string tsSATE = "";
    getSATEName( "HLT_mu26_ivarmedium", tsSATE );
    m_trigTagSATEName.push_back( tsSATE );
    std::string tsCBTE = "";
    getCBTEName( "HLT_mu26_ivarmedium", tsCBTE );
    m_trigTagCBTEName.push_back( tsCBTE );
  }
  else if(m_method=="NoMass"){
//    m_l1chainList.push_back("L1_MU6");
    m_hltchainList.push_back("HLT_mu6");
    m_trigchainList.push_back("HLT_mu6");
    m_isPassed.push_back(false);
    m_reqdR = false;
    m_ignoreMassreq = true;
    m_massMin.push_back(0);
    m_massMax.push_back(0);
    std::string tsSATE = "";
    getSATEName( "HLT_mu6", tsSATE );
    m_trigTagSATEName.push_back( tsSATE );
    std::string tsCBTE = "";
    getCBTEName( "HLT_mu6", tsCBTE );
    m_trigTagCBTEName.push_back( tsCBTE );
  }
  else if(m_method=="NoTag"){
//    m_l1chainList.push_back("L1_MU6");
    m_hltchainList.push_back("dummy");
    m_trigchainList.push_back("dummy");
    m_isPassed.push_back(true);
    m_reqdR = false;
    m_ignoreMassreq = true;
    m_massMin.push_back(0);
    m_massMax.push_back(0);
    m_trigTagSATEName.push_back( "dummy" );
    m_trigTagCBTEName.push_back( "dummy" );
  }
  m_nChain = m_trigchainList.size();
  return 0;
}

//---------------------------------------------------------//
//---------------------------------------------------------//

void TagAndProbe::clear() {
  m_tappairs.clear();
  m_tag.clear();
  m_probe.clear();
  m_requirements_tag.clear();
  m_L1objects_tag.clear();
  m_SAobjects_tag.clear();
  m_CBobjects_tag.clear();
  m_EFobjects_tag.clear();
  m_requirements_probe.clear();
  m_vL1objects_probe.clear();
  m_vSAobjects_probe.clear();
  m_vCBobjects_probe.clear();
  m_vEFobjects_probe.clear();
  m_passTrigmes.clear();
}

//---------------------------------------------------------//
//---------------------------------------------------------//

void TagAndProbe::addMesChain( const string& L1trig, const string& HLTtrig )
{
  m_L1trigmesName.push_back(L1trig);
  m_HLTtrigmesName.push_back(HLTtrig);
  std::string trigSATE = "";
  getSATEName( HLTtrig, trigSATE );
  m_trigmesSATEname.push_back( trigSATE );
  std::string trigCBTE = "";
  getCBTEName( HLTtrig, trigCBTE );
  m_trigmesCBTEname.push_back( trigCBTE );

  m_nmesChain = m_HLTtrigmesName.size();
}

//---------------------------------------------------------//
//---------------------------------------------------------//

void TagAndProbe::checkMesChain()
{

  m_passTrigmes.clear();
  for(size_t i=0; i<m_HLTtrigmesName.size(); i++) {
    m_passTrigmes.push_back( m_trigDecTool->isPassed( m_HLTtrigmesName.at(i), TrigDefs::eventAccepted ) );
  }

}

//---------------------------------------------------------//
//---------------------------------------------------------//

bool TagAndProbe::isPassedTrigger()
{
  if(m_message <= DEBUG)
    std::cout << "TagAndProbe::isPassedTrigger =====" << std::endl;
  bool isTriggered = false;
  if(m_method=="NoTag"){
    return true;
  }
  for(int i = 0; i < m_nChain; i++){
    if(m_trigDecTool->isPassed(m_trigchainList.at(i))){
      std::cout << m_trigchainList.at(i) << " is passed" << std::endl; 
      m_isPassed.at(i) = true;
      isTriggered = true;
    }
  }
  return isTriggered;
}

//---------------------------------------------------------//
//---------------------------------------------------------//
//when dRl1tag-probe is 0.12 ? --> m_method is "NoMass"
bool TagAndProbe::setProbes( const xAOD::MuonContainer& muons )
{
  for(const auto& tag : muons){
    if(m_pst.getQualityOfTrack( tag ) < 0) continue;

    for(const auto& probe : muons){
      if(tag==probe) continue;
      if(m_pst.getQualityOfTrack( probe ) < 0) continue;
      
      // tag offline -- trigger matching start
      for(int ichain = 0; ichain < m_nChain; ichain++){
        if(m_isPassed.at(ichain) == false) continue;
        if( m_message <= DEBUG ) std::cout << ">>================================<<" << std::endl;
        Requirement req;
        if(!passDimuSelection(tag, probe, ichain, req)) continue;
        Trig::FeatureContainer fc = m_trigDecTool->features( m_hltchainList.at(ichain), TrigDefs::alsoDeactivateTEs );
        Trig::FeatureContainer fc_fullScan = m_trigDecTool->features( m_trigchainList.at(ichain), TrigDefs::alsoDeactivateTEs );
        req.reqdRl1 = m_reqdR ? dRl1bypt( tag->pt() ) : 0.12;
        req.EFmatchingdR = m_reqdR ? 0.01 : 0.05;
        L1Object l1obj;
        l1obj.dRl1 = m_reqdR ? dRl1bypt( tag->pt() ) : 0.12;
        SAObject saobj;
        CBObject cbobj;
        EFObject efobj;
        efobj.dRef = m_reqdR ? 0.01 : 0.05;
        if(m_method=="NoTag"){
          std::cout << ">> ACCEPT AS TAG MUON <<" << std::endl;
          std::pair<const xAOD::Muon*, const xAOD::Muon*> tappair = make_pair( tag, probe );
          //push_back infos of tag
          m_tappairs.push_back(tappair);
          m_requirements_tag.push_back(req);
          m_L1objects_tag.push_back(l1obj);
          m_SAobjects_tag.push_back(saobj);
          m_CBobjects_tag.push_back(cbobj);
          m_EFobjects_tag.push_back(efobj);
        } else {
          if(m_isnoL1){
            for( const Trig::Combination comb : fc.getCombinations()){
              if(!m_tmt.matchL1(comb, tag, l1obj)) continue;
              if(!m_tmt.matchSA(comb, m_trigTagSATEName.at(ichain), l1obj, saobj)) continue;
              if(!m_tmt.matchCB(comb, m_trigTagCBTEName.at(ichain), saobj, cbobj)) continue;
              if(!m_tmt.matchEF(fc_fullScan, tag, efobj)) continue;
            }
          } else {
            for( const Trig::Combination comb : fc.getCombinations()){
              // start matching
              if(!m_tmt.matchL1(comb, tag, l1obj)) continue;
              if(!m_tmt.matchSA(comb, m_trigTagSATEName.at(ichain), l1obj, saobj)) continue;
              if(!m_tmt.matchCB(comb, m_trigTagCBTEName.at(ichain), saobj, cbobj)) continue;
              if(!m_tmt.matchEF(comb, tag, efobj)) continue;
            }
            std::cout << "--> matched L1 roiNum/eta/phi/dR(mu-roi) = " << l1obj.roiNum << "/" << l1obj.eta << "/" << l1obj.phi << "/" << l1obj.dRl1 << std::endl;
            std::cout << "--> matched SA pt/eta/phi = " << saobj.pt << "/" << saobj.eta << "/" << saobj.phi << std::endl;
            std::cout << "--> matched CB pt/eta/phi = " << cbobj.pt << "/" << cbobj.eta << "/" << cbobj.phi << std::endl;
            std::cout << "--> matched EF pt/eta/phi = " << efobj.pt << "/" << efobj.eta << "/" << efobj.phi << std::endl;
          }
          // check if trigger passed
          if(l1obj.isPassed == 1 && saobj.isPassed == 1 && cbobj.isPassed == 1 && efobj.isPassed == 1){
            std::cout << ">> ACCEPT AS TAG MUON <<" << std::endl;
            std::pair<const xAOD::Muon*, const xAOD::Muon*> tappair = make_pair( tag, probe );
            //push_back infos of tag
            m_tappairs.push_back(tappair);
            m_requirements_tag.push_back(req);
            m_L1objects_tag.push_back(l1obj);
            m_SAobjects_tag.push_back(saobj);
            m_CBobjects_tag.push_back(cbobj);
            m_EFobjects_tag.push_back(efobj);
          }
        }
      }
    }
  }
  if( !m_tappairs.empty() ) return true;
  return false;
}

//---------------------------------------------------------//
//---------------------------------------------------------//

void TagAndProbe::doProbeMatching( SG::ReadHandle<xAOD::MuonRoIContainer> &rois )
{
  std::cout << "=====> doProbeMatching <=====" << std::endl;
  for(int ipair = 0; ipair < (int)m_tappairs.size(); ipair++){
    const xAOD::Muon* tag = m_tappairs.at(ipair).first;
    const xAOD::Muon* probe = m_tappairs.at(ipair).second;    
    const xAOD::TrackParticle* tagtrk = tag->trackParticle( xAOD::Muon::InnerDetectorTrackParticle );
    pair< double, double > tagextEtaAndPhi = m_ext.extTrack( tagtrk );
    const double tagExtEta = tagextEtaAndPhi.first;
    const double tagExtPhi = tagextEtaAndPhi.second;
    const xAOD::TrackParticle* probetrk = probe->trackParticle( xAOD::Muon::InnerDetectorTrackParticle );
    pair< double, double > probeextEtaAndPhi = m_ext.extTrack( probetrk );
    const double probeExtEta = probeextEtaAndPhi.first;
    const double probeExtPhi = probeextEtaAndPhi.second;
    const double dRext = m_utils.deltaR( tagExtEta, tagExtPhi, probeExtEta, probeExtPhi );
    OfflineObject probeobj, tagobj;
    probeobj.pt = probe->pt();
    probeobj.eta = probe->eta();
    probeobj.phi = probe->phi();
    probeobj.extEta = probeExtEta;
    probeobj.extPhi = probeExtPhi;
    probeobj.tpextdR = dRext;
    probeobj.d0 = (probetrk)? probetrk->d0() : -99999;
    probeobj.z0 = (probetrk)? probetrk->z0() : -99999;
    
    probeobj.segmentN = probe->nMuonSegments();
    for(int i_seg = 0; i_seg < probeobj.segmentN; i_seg++){
      const xAOD::MuonSegment* segment = probe->muonSegment(i_seg);
      if(!segment) continue;
      probeobj.segmentX[i_seg]              = segment->x();
      probeobj.segmentY[i_seg]              = segment->y();
      probeobj.segmentZ[i_seg]              = segment->z();
      probeobj.segmentPx[i_seg]             = segment->px();
      probeobj.segmentPy[i_seg]             = segment->py();
      probeobj.segmentPz[i_seg]             = segment->pz();
      probeobj.segmentChiSquared[i_seg]     = segment->chiSquared();
      probeobj.segmentNumberDoF[i_seg]      = segment->numberDoF();
      probeobj.segmentSector[i_seg]         = segment->sector();
      probeobj.segmentChamberIndex[i_seg]   = segment->chamberIndex();
      probeobj.segmentEtaIndex[i_seg]       = segment->etaIndex();
      probeobj.segmentNPrecisionHits[i_seg] = segment->nPrecisionHits();
      probeobj.segmentNPhiLayers[i_seg]     = segment->nPhiLayers();
      probeobj.segmentNTrigEtaLayers[i_seg] = segment->nTrigEtaLayers();
    }
    m_probe.push_back(probeobj);
    tagobj.pt = tag->pt();
    tagobj.eta = tag->eta();
    tagobj.phi = tag->phi();
    tagobj.extEta = tagExtEta;
    tagobj.extPhi = tagExtPhi;
    tagobj.tpextdR = dRext;
    m_tag.push_back(tagobj);
    std::cout << "#tag-probe extdR : " << dRext << "<<=======" << std::endl;
    Requirement req;
    req.reqdRl1 = m_reqdR ? dRl1bypt( probe->pt() ) : 0.12;
    req.EFmatchingdR = m_reqdR ? 0.01 : 0.05;
    m_requirements_probe.push_back(req);
    
    L1Objects L1objects;
    SAObjects SAobjects;
    CBObjects CBobjects;
    EFObjects EFobjects;
    for(int ichain = 0; ichain < m_nmesChain; ichain++){
      L1Object l1obj;
      l1obj.dRl1 = m_reqdR ? dRl1bypt( probe->pt() ) : 0.12;
      SAObject saobj;
      CBObject cbobj;
      EFObject efobj;
      efobj.dRef = m_reqdR ? 0.01 : 0.05;
      std::cout << "##Trigger chain --> " << "L1: " << m_L1trigmesName.at(ichain) << "/ HLT: " << m_HLTtrigmesName.at(ichain) << "<<=======" << std::endl;
      Trig::FeatureContainer fc = m_trigDecTool->features( m_HLTtrigmesName.at(ichain), TrigDefs::alsoDeactivateTEs );
            
      m_tmt.matchL1(rois, probe, l1obj, m_L1trigmesName.at(ichain));
      m_tmt.matchSA(fc, m_trigmesSATEname.at(ichain), l1obj, saobj);
      m_tmt.matchCB(fc, m_trigmesCBTEname.at(ichain), saobj, cbobj);
      m_tmt.matchEF(fc, probe, efobj);
    
      std::cout << "--> matched L1 roiNum/eta/phi/dR(mu-roi) = " << l1obj.roiNum << "/" << l1obj.eta << "/" << l1obj.phi << "/" << l1obj.dRl1 << std::endl;
      std::cout << "--> matched SA pt/eta/phi = " << saobj.pt << "/" << saobj.eta << "/" << saobj.phi << std::endl;
      std::cout << "--> matched CB pt/eta/phi = " << cbobj.pt << "/" << cbobj.eta << "/" << cbobj.phi << std::endl;
      std::cout << "--> matched EF pt/eta/phi = " << efobj.pt << "/" << efobj.eta << "/" << efobj.phi << std::endl;
      std::cout << "Trigger status : " << "L1 --> " << l1obj.isPassed << std::endl; 
      std::cout << "                 " << "SA --> " << saobj.isPassed << std::endl;
      std::cout << "                 " << "CB --> " << cbobj.isPassed << std::endl;
      std::cout << "                 " << "EF --> " << efobj.isPassed << std::endl;
      L1objects.push_back(l1obj);
      SAobjects.push_back(saobj);
      CBobjects.push_back(cbobj);
      EFobjects.push_back(efobj);
    } // ichain loop end
    m_vL1objects_probe.push_back(L1objects);
    m_vSAobjects_probe.push_back(SAobjects);
    m_vCBobjects_probe.push_back(CBobjects);
    m_vEFobjects_probe.push_back(EFobjects);
  } // muon pair loop end

}

//---------------------------------------------------------//
//---------------------------------------------------------//

double TagAndProbe::dRl1bypt( double mupt ) const
{
  double dR = 0.08;
  if( mupt < 10000. ) {
    dR = -0.00001*mupt + 0.18;
  }
  return dR;
}

//---------------------------------------------------------//
//---------------------------------------------------------//

bool TagAndProbe::passDimuSelection( const xAOD::Muon* tag, 
                                       const xAOD::Muon* probe, 
                                       int chain, 
                                       Requirement& req )
{
  if(m_message <= DEBUG) std::cout << "start dimuon selection" << std::endl;
  bool passMassreq = false;
  bool passdphireq = false;
  bool passdRreq = false;
  bool passOSreq = false;
  bool passMuqual = false;
  bool passchi2req = false;
  float chi2tag = 20.;
  float chi2probe = 20.;
  tag->parameter(chi2tag, xAOD::Muon::msInnerMatchChi2);
  probe->parameter(chi2probe, xAOD::Muon::msInnerMatchChi2);
  if( (chi2tag < m_chi2cut) && (chi2probe < m_chi2cut) ) passchi2req = true;
  if(m_run3evtSelection){
    if(tag->quality() == xAOD::Muon::Tight) passMuqual = true;
  } else {
    if(m_pst.getQualityOfTrack( tag ) == 0) passMuqual = true;
  }
  if(m_run3evtSelection){
    if(probe->quality() == xAOD::Muon::Medium || 
       probe->quality() == xAOD::Muon::Tight) passMuqual = true;
  } else {
    if(m_pst.getQualityOfTrack( probe ) == 0) passMuqual = true;
  }
      
  // extrapolate the muon's (eta, phi) to MS
  const xAOD::TrackParticle* tagtrk = tag->trackParticle( xAOD::Muon::InnerDetectorTrackParticle );
  if( !tagtrk ) return false;
  pair< double, double > tagextEtaAndPhi = m_ext.extTrack( tagtrk );
  const double tagExtEta = tagextEtaAndPhi.first;
  const double tagExtPhi = tagextEtaAndPhi.second;
  const xAOD::TrackParticle* probetrk = probe->trackParticle( xAOD::Muon::InnerDetectorTrackParticle );
  if( !probetrk ) return false;
  pair< double, double > probeextEtaAndPhi = m_ext.extTrack( probetrk );
  const double probeExtEta = probeextEtaAndPhi.first;
  const double probeExtPhi = probeextEtaAndPhi.second;
  const double dRext = m_utils.deltaR( tagExtEta, tagExtPhi, probeExtEta, probeExtPhi );
  if( m_message <= DEBUG ) std::cout << "... tag-probe extdR : " << dRext << std::endl;
  if(dRext > 0.2 || m_method=="NoTag") passdRreq = true;
  if(tag->charge()*probe->charge() < 0) passOSreq = true;
  TLorentzVector lvtag, lvprobe;
  lvtag.SetPhi( tag->phi() ); 
  lvprobe.SetPhi( probe->phi() );

  double dphi = TVector2::Phi_mpi_pi( tag->phi() - probe->phi() );
  double mass = (tag->p4() + probe->p4()).M();
  if( m_message <= DEBUG ) std::cout << "... dimuon mass : " << mass << std::endl;
  if( m_message <= DEBUG ) std::cout << "... tag-probe dphi : " << dphi << std::endl;
  req.tpdPhi = dphi;
  req.invMass = mass;
  if( (mass > m_massMin.at(chain) && mass < m_massMax.at(chain)) || m_ignoreMassreq ) passMassreq = true;
  if(std::fabs(dphi) < 3.0 || m_method=="NoTag") passdphireq = true;
  if(passMassreq && passdphireq && passOSreq && passdRreq && passMuqual && passchi2req)
    return true;
  return false;
}
//---------------------------------------------------------//
//---------------------------------------------------------//

Bool_t TagAndProbe::getSATEName( const std::string& mesHLT, std::string& teName )
{
  ///This function is used in "int TagAndProbe::addMesChain ".
  ///mesHLT and teName is parameters.

  TString hlt = mesHLT.c_str();
  TObjArray* toa = hlt . Tokenize( "_" );
  Int_t thrValue = 0;
  Bool_t isBarrelOnly = kFALSE;
  for( Int_t i = 0; i < toa -> GetEntries(); i++ ){
    TString tsToken = (static_cast<TObjString*>( toa -> At(i) ) ) -> String();
    if ( tsToken . Contains( "0eta105" ) ){
      isBarrelOnly = kTRUE;
    }
    if ( !tsToken . Contains( "mu" ) )
      continue;
    tsToken = tsToken.ReplaceAll( "noL1", "" );
    tsToken = tsToken.ReplaceAll( "3mu", "" );
    tsToken = tsToken.ReplaceAll( "2mu", "" );
    tsToken = tsToken.ReplaceAll( "mu", "" );
    if ( !tsToken . IsDec() )
      continue;
    thrValue = tsToken . Atoi();
  }
  if ( thrValue == 0 )
    return 1;
  if ( thrValue == 4 || thrValue == 2 || isBarrelOnly ){
  } else {
    thrValue = 6;
  }
  teName = Form( "%dGeV%s_v15a", thrValue, (isBarrelOnly)?("_barrelOnly"):("") );
  return 0;
}

//---------------------------------------------------------//
//---------------------------------------------------------//

Bool_t TagAndProbe::getCBTEName( const std::string& mesHLT, std::string& teName )
{
  TString hlt = mesHLT.c_str();
  TObjArray* toa = hlt . Tokenize( "_" );
  Int_t thrValue = 0;
  for( Int_t i = 0; i < toa -> GetEntries(); i++ ){
    TString tsToken = (static_cast<TObjString*>( toa -> At(i) ) ) -> String();
    if ( !tsToken . Contains( "mu" ) )
      continue;
    tsToken = tsToken.ReplaceAll( "noL1", "" );
    tsToken = tsToken.ReplaceAll( "3mu", "" );
    tsToken = tsToken.ReplaceAll( "2mu", "" );
    tsToken = tsToken.ReplaceAll( "mu", "" );
    if ( !tsToken . IsDec() )
      continue;
    thrValue = tsToken . Atoi();
  }
  if ( thrValue == 0 )
    return 1;
  if ( thrValue >= 24 )
    thrValue = 22;
  teName = Form( "_mu%d_", thrValue );
  return 0;
}

//---------------------------------------------------------//
//---------------------------------------------------------//

