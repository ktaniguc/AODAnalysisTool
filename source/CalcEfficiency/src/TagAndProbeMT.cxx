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

#include "CalcEfficiency/TagAndProbeMT.h"

#include "TMath.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TVector3.h"

TagAndProbeMT::TagAndProbeMT(){ }

int TagAndProbeMT::initialize( const int& message, 
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
  m_matchToolMT.initialize(tdt, ext);

  if(m_method=="Jpsi"){
    m_l1chainList.push_back("L1_MU20");
    m_hltchainList.push_back("HLT_mu20_L1MU20");
//    m_trigchainList.push_back("HLT_mu6");
    m_trigchainList.push_back("HLT_mu20_2mu2noL1_JpsimumuFS");
    m_isPassed.push_back(false);
    m_reqdR = true;
    m_isnoL1 = true;
    m_massMin.push_back(2700.);
    m_massMax.push_back(3500.);
    std::string tsSATE = "";
//    getSATEName( "HLT_mu20", tsSATE );
    m_trigTagSATEName.push_back( tsSATE );
    std::string tsCBTE = "";
//    getCBTEName( "HLT_mu20", tsCBTE );
    m_trigTagCBTEName.push_back( tsCBTE );
  }
  else if(m_method=="JpsiMC"){
    m_l1chainList.push_back("L1_MU6");
    m_hltchainList.push_back("HLT_mu6_L1MU6");
//    m_trigchainList.push_back("HLT_mu6");
    m_trigchainList.push_back("HLT_mu6_L1MU6");
    m_isPassed.push_back(false);
    m_reqdR = true;
    m_isnoL1 = true;
    m_massMin.push_back(2700.);
    m_massMax.push_back(3500.);
    std::string tsSATE = "";
//    getSATEName( "HLT_mu20", tsSATE );
    m_trigTagSATEName.push_back( tsSATE );
    std::string tsCBTE = "";
//    getCBTEName( "HLT_mu20", tsCBTE );
    m_trigTagCBTEName.push_back( tsCBTE );
  }
  else if(m_method=="Ztap"){
    m_l1chainList.push_back("L1_MU20");
    m_hltchainList.push_back("HLT_mu26_ivarmedium_L1MU20");
    m_trigchainList.push_back("HLT_mu26_ivarmedium_L1MU20");
    m_isPassed.push_back(false);
    m_reqdR = true;
    m_massMin.push_back(90000.-10000.);
    m_massMax.push_back(90000.+10000.);
    std::string tsSATE = "";
//    getSATEName( "HLT_mu26_ivarmedium", tsSATE );
    m_trigTagSATEName.push_back( tsSATE );
    std::string tsCBTE = "";
//    getCBTEName( "HLT_mu26_ivarmedium", tsCBTE );
    m_trigTagCBTEName.push_back( tsCBTE );
  }
  else if(m_method=="NoMass"){
    m_l1chainList.push_back("L1_MU20");
    m_hltchainList.push_back("HLT_mu26_ivarmedium_L1MU20");
    m_trigchainList.push_back("HLT_mu26_ivarmedium_L1MU20");
//    m_l1chainList.push_back("L1_MU6");
//    m_hltchainList.push_back("HLT_mu6_L1MU6");
//    m_trigchainList.push_back("HLT_mu6_L1MU6");
    m_isPassed.push_back(false);
    m_reqdR = false;
    m_ignoreMassreq = true;
    m_massMin.push_back(0);
    m_massMax.push_back(0);
    std::string tsSATE = "";
//    getSATEName( "HLT_mu6", tsSATE );
    m_trigTagSATEName.push_back( tsSATE );
    std::string tsCBTE = "";
//    getCBTEName( "HLT_mu6", tsCBTE );
    m_trigTagCBTEName.push_back( tsCBTE );
  }
  else if(m_method=="NoTag"){
    m_l1chainList.push_back("dummy");
    m_hltchainList.push_back("dummy");
    m_trigchainList.push_back("dummy");
    m_isPassed.push_back(true);
    m_reqdR = false;
    m_ignoreMassreq = true;
    m_massMin.push_back(0);
    m_massMax.push_back(0);
    std::string tsSATE = "";
//    getSATEName( "HLT_mu6", tsSATE );
    m_trigTagSATEName.push_back( tsSATE );
    std::string tsCBTE = "";
//    getCBTEName( "HLT_mu6", tsCBTE );
    m_trigTagCBTEName.push_back( tsCBTE );
  }
  else if(m_method=="NoTagJpsi"){
    m_l1chainList.push_back("dummy");
    m_hltchainList.push_back("dummy");
    m_trigchainList.push_back("dummy");
    m_isPassed.push_back(true);
    m_reqdR = false;
    m_ignoreMassreq = false;
    m_massMin.push_back(2700.);
    m_massMax.push_back(3500.);
    std::string tsSATE = "";
    m_trigTagSATEName.push_back( tsSATE );
    std::string tsCBTE = "";
    m_trigTagCBTEName.push_back( tsCBTE );
  }
  m_nChain = m_trigchainList.size();
  return 0;
}

//---------------------------------------------------------//
//---------------------------------------------------------//

void TagAndProbeMT::addMesChain( const string& L1trig, const string& HLTtrig )
{
  if(HLTtrig.find("b") != std::string::npos) 
    std::cout << "WARNING         this chain maybe 'nomucomb' trigger. if so, isPassedCB are not correct" << std::endl;
  m_L1trigmesName.push_back(L1trig);
  m_HLTtrigmesName.push_back(HLTtrig);

  m_nmesChain = m_HLTtrigmesName.size();
}

//---------------------------------------------------------//
//---------------------------------------------------------//

void TagAndProbeMT::checkMesChain()
{

  m_passTrigmes.clear();
  for(size_t i=0; i<m_HLTtrigmesName.size(); i++) {
    m_passTrigmes.push_back( m_trigDecTool->isPassed( m_HLTtrigmesName.at(i), TrigDefs::eventAccepted ) );
  }

}

//---------------------------------------------------------//
//---------------------------------------------------------//

void TagAndProbeMT::clear() {
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
  m_tag_L1_pass.clear();
  m_tag_L1_roiNum.clear();
  m_tag_SA_pass.clear();
  m_tag_CB_pass.clear();
  m_tag_EF_pass.clear();
}

//---------------------------------------------------------//
//---------------------------------------------------------//

bool TagAndProbeMT::isPassedTrigger()
{
  std::cout << "TagAndProbeMT::isPassedTrigger =====" << std::endl;
  bool isTriggered = false;
  if(m_method.Contains("NoTag")){
    return true;
  }
  for(int i = 0; i < m_nChain; i++){
    if(m_trigDecTool->isPassed( m_trigchainList.at(i), TrigDefs::eventAccepted )){
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
bool TagAndProbeMT::setProbes( SG::ReadHandle<xAOD::MuonContainer> &muons,
                               SG::ReadHandle<xAOD::MuonRoIContainer> &rois 
                              )
{
  for(const xAOD::Muon* tag : *muons){
    for(const xAOD::Muon* probe : *muons){
      if(tag==probe) continue;
      // tag offline -- trigger matching start
      for(int ichain = 0; ichain < m_nChain; ichain++){
        if(m_isPassed.at(ichain) == false) continue;
        std::cout << ">>================================<<" << std::endl;
        Requirement req;
        if(!passDimuSelection(tag, probe, ichain, req)) continue;
        req.reqdRl1 = m_reqdR ? dRl1bypt( tag->pt() ) : 0.12;
        req.EFmatchingdR = m_reqdR ? 0.01 : 0.05;
        L1Object l1obj;
        l1obj.dRl1 = m_reqdR ? dRl1bypt( tag->pt() ) : 0.12;
        SAObject saobj;
        CBObject cbobj;
        EFObject efobj;
        efobj.dRef = m_reqdR ? 0.01 : 0.05;
        // start matching for tag
	if( m_matchToolMT.matchL1( rois, tag, l1obj, m_l1chainList.at(ichain), m_hltchainList.at(ichain) ) )
	  std::cout << "--> matched L1 roiNum/eta/phi/dR(mu-roi) = " << l1obj.roiNum << "/" << l1obj.eta << "/" << l1obj.phi << "/" << l1obj.dRl1 << std::endl;
	if( m_matchToolMT.matchSA( m_hltchainList.at(ichain), l1obj, saobj ) )
	  std::cout << "--> matched SA pt/eta/phi = " << saobj.pt << "/" << saobj.eta << "/" << saobj.phi << std::endl;
	if( m_matchToolMT.matchCB( m_hltchainList.at(ichain), saobj, cbobj ) )
	  std::cout << "--> matched CB pt/eta/phi = " << cbobj.pt << "/" << cbobj.eta << "/" << cbobj.phi << std::endl;
	if( m_matchToolMT.matchEF( m_hltchainList.at(ichain), tag, efobj ) )
	  std::cout << "--> matched EF pt/eta/phi = " << efobj.pt << "/" << efobj.eta << "/" << efobj.phi << std::endl;
	// check if tag passes trigger
        if(m_method.Contains("NoTag")
	   ||
	   (l1obj.isPassed  == 1 && saobj.isPassed == 1 && cbobj.isPassed == 1 && efobj.isPassed == 1))
	  std::cout << ">> ACCEPT AS TAG MUON <<" << std::endl;
	else
	  continue;
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
  return true;
}

//---------------------------------------------------------//
//---------------------------------------------------------//

void TagAndProbeMT::doProbeMatching( SG::ReadHandle<xAOD::MuonRoIContainer> &rois )
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
    std::cout << "# of MuonSegment = " << probeobj.segmentN << std::endl;
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
    tagobj.d0 = (tagtrk)? tagtrk->d0() : -99999;
    tagobj.z0 = (tagtrk)? tagtrk->z0() : -99999;
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
            
      m_matchToolMT.matchL1( rois, probe, l1obj, m_L1trigmesName.at(ichain), m_HLTtrigmesName.at(ichain) );
      std::cout << "--> matched L1 roiNum/eta/phi/dR(mu-roi) = " << l1obj.roiNum << "/" << l1obj.eta << "/" << l1obj.phi << "/" << l1obj.dRl1 << std::endl;
      m_matchToolMT.matchSA( m_HLTtrigmesName.at(ichain), l1obj, saobj );
      std::cout << "--> matched SA pt/eta/phi = " << saobj.pt << "/" << saobj.eta << "/" << saobj.phi << std::endl;
      m_matchToolMT.matchCB( m_HLTtrigmesName.at(ichain), saobj, cbobj );
      std::cout << "--> matched CB pt/eta/phi = " << cbobj.pt << "/" << cbobj.eta << "/" << cbobj.phi << std::endl;
      m_matchToolMT.matchEF( m_HLTtrigmesName.at(ichain), probe, efobj );
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
void TagAndProbeMT::doTagMatching( SG::ReadHandle<xAOD::MuonRoIContainer> &rois )
{
  std::cout << "=====> doTagMatching <=====" << std::endl;
  if(m_tappairs.size() == 0) return;

  const xAOD::Muon* tag = m_tappairs.at(0).first;

  for(int ichain = 0; ichain < m_nmesChain; ichain++){
    L1Object l1obj;
    l1obj.dRl1 = m_reqdR ? dRl1bypt( tag->pt() ) : 0.12;
    SAObject saobj;
    CBObject cbobj;
    EFObject efobj;
    efobj.dRef = m_reqdR ? 0.01 : 0.05;
    std::cout << "TagMatching ##Trigger chain --> " << "L1: " << m_L1trigmesName.at(ichain) << "/ HLT: " << m_HLTtrigmesName.at(ichain) << "<<=======" << std::endl;

    m_matchToolMT.matchL1( rois, tag, l1obj, m_L1trigmesName.at(ichain), m_HLTtrigmesName.at(ichain) );
    std::cout << "TagMatching --> matched L1 roiNum/eta/phi/dR(mu-roi) = " << l1obj.roiNum << "/" << l1obj.eta << "/" << l1obj.phi << "/" << l1obj.dRl1 << std::endl;
    m_matchToolMT.matchSA( m_HLTtrigmesName.at(ichain), l1obj, saobj );
    std::cout << "TagMatching --> matched SA pt/eta/phi = " << saobj.pt << "/" << saobj.eta << "/" << saobj.phi << std::endl;
    m_matchToolMT.matchCB( m_HLTtrigmesName.at(ichain), saobj, cbobj );
    std::cout << "TagMatching --> matched CB pt/eta/phi = " << cbobj.pt << "/" << cbobj.eta << "/" << cbobj.phi << std::endl;
    m_matchToolMT.matchEF( m_HLTtrigmesName.at(ichain), tag, efobj );
    std::cout << "TagMatching --> matched EF pt/eta/phi = " << efobj.pt << "/" << efobj.eta << "/" << efobj.phi << std::endl;
    std::cout << "                             " << "L1 --> " << saobj.isPassed << std::endl;
    std::cout << "                             " << "SA --> " << saobj.isPassed << std::endl;
    std::cout << "                             " << "CB --> " << cbobj.isPassed << std::endl;
    std::cout << "                             " << "EF --> " << efobj.isPassed << std::endl;
    m_tag_L1_pass.push_back(l1obj.isPassed);
    m_tag_L1_roiNum.push_back(l1obj.roiNum);
    m_tag_SA_pass.push_back(saobj.isPassed);
    m_tag_CB_pass.push_back(cbobj.isPassed);
    m_tag_EF_pass.push_back(efobj.isPassed);
  } // ichain loop end

}

//---------------------------------------------------------//
//---------------------------------------------------------//

double TagAndProbeMT::dRl1bypt( double mupt ) const
{
  double dR = 0.08;
  if( mupt < 10000. ) {
    dR = -0.00001*mupt + 0.18;
  }
  return dR;
}

//---------------------------------------------------------//
//---------------------------------------------------------//

bool TagAndProbeMT::passDimuSelection( const xAOD::Muon* tag, 
                                       const xAOD::Muon* probe, 
                                       int chain, 
                                       Requirement& req )
{
  std::cout << "start dimuon selection" << std::endl;
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
  std::cout << "... tag-probe extdR : " << dRext << std::endl;
  if(dRext > 0.2 || m_method.Contains("NoTag")) passdRreq = true;
  if(tag->charge()*probe->charge() < 0) passOSreq = true;
  TLorentzVector lvtag, lvprobe;
  lvtag.SetPhi( tag->phi() ); 
  lvprobe.SetPhi( probe->phi() );

  double dphi = TVector2::Phi_mpi_pi( tag->phi() - probe->phi() );
  double mass = (tag->p4() + probe->p4()).M();
  std::cout << "... dimuon mass : " << mass << std::endl;
  std::cout << "... tag-probe dphi : " << dphi << std::endl;
  req.tpdPhi = dphi;
  req.invMass = mass;
  if( (mass > m_massMin.at(chain) && mass < m_massMax.at(chain)) || m_ignoreMassreq ) passMassreq = true;
  if(std::fabs(dphi) < 3.0 || m_method.Contains("NoTag")) passdphireq = true;
  if(passchi2req && passMuqual && passMassreq && passdphireq && passOSreq && passdRreq)
    return true;
  return false;
}
//---------------------------------------------------------//
//---------------------------------------------------------//

