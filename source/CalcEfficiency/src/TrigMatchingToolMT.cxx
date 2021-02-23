#include "CalcEfficiency/TrigMatchingToolMT.h"
#include "TVector3.h"

TrigMatchingToolMT::TrigMatchingToolMT(){ }

//---------------------------------------------------------//
//---------------------------------------------------------//
//in the end, L1 requirement dRl1tag = minimum dR between l1RoI and offline
//if there are no RoIs fulfill the initial req, dRl1tag val isn't changed
bool TrigMatchingToolMT::matchL1( SG::ReadHandle<xAOD::MuonRoIContainer> &murois,
                                  const xAOD::Muon* muon, 
                                  L1Object& l1obj,
                                  std::string trig,
                                  std::string HLTtrig )
{
  std::cout << "start L1 matching, trigger ----> " << trig << std::endl;
  l1obj.isPassedByEvent_TDT = m_trigDecTool->isPassed(trig, TrigDefs::eventAccepted);
  l1obj.isPassedByEvent = m_trigDecTool->isPassed(trig, TrigDefs::eventAccepted);
  
  L1Objects L1objs_forTrigDec;
  std::vector< TrigCompositeUtils::LinkInfo<xAOD::L2StandAloneMuonContainer> > l2saLinks = m_trigDecTool->features<xAOD::L2StandAloneMuonContainer>( HLTtrig, TrigDefs::includeFailedDecisions );
  for(const TrigCompositeUtils::LinkInfo<xAOD::L2StandAloneMuonContainer>& l2saLinkInfo : l2saLinks){
    if( !l2saLinkInfo.isValid() ) continue;
    const ElementLink<xAOD::L2StandAloneMuonContainer> l2sa = l2saLinkInfo.link;
    if( !l2sa.isValid() ) continue;
    // start L1 -- SA matching
    const TrigCompositeUtils::Decision* muDecision = l2saLinkInfo.source; 
    const TrigCompositeUtils::LinkInfo<TrigRoiDescriptorCollection> roiLinkInfo = TrigCompositeUtils::findLink<TrigRoiDescriptorCollection>(muDecision, "initialRoI");
    if( !roiLinkInfo.isValid() ) continue;
    L1Object l1obj;
    const ElementLink<TrigRoiDescriptorCollection> roiEL = roiLinkInfo.link;
    l1obj.eta = (*roiEL)->eta();
    l1obj.phi = (*roiEL)->phi();
    l1obj.roiWord = (*roiEL)->roiWord();
    //std::cout << "... [TrigRoiDescriptor] L1RoIs seeded by muFast: roiEta/roiPhi/roiWord = " << (*roiEL)->eta() << "/" << (*roiEL)->phi() << "/" << (*roiEL)->roiWord() << std::endl;
    L1objs_forTrigDec.push_back(l1obj);
  }
  // extrapolate the muon's (eta, phi) to MS
  const xAOD::TrackParticle* mutrk = muon->trackParticle( xAOD::Muon::InnerDetectorTrackParticle );
  std::pair< double, double > extEtaAndPhi = m_ext.extTrack( mutrk );
  const double muExtEta = extEtaAndPhi.first;
  const double muExtPhi = extEtaAndPhi.second;

  for(const auto &roi : *murois){
    std::cout << "... L1 RoI eta/phi/thre/roiNum/roiWord = " << roi->eta() << "/" << roi->phi() << "/" << roi->thrValue() << "/" << roi->getRoI() << "/" << roi->roiWord() << std::endl;
    const double l1eta = roi->eta();
    const double l1phi = roi->phi();
    const int    l1thr = roi->getThrNumber();
    const double roidRext = m_utils.deltaR( muExtEta, muExtPhi, l1eta, l1phi );
    std::cout << "... L1RoI vs offline mu dR = " << roidRext << std::endl;
    bool passdR = (l1obj.dRl1 > roidRext);
    bool isSeededSA = false;
    for(int i_link = 0; i_link < (int)L1objs_forTrigDec.size(); i_link++){
      if(roi->roiWord() == L1objs_forTrigDec.at(i_link).roiWord) isSeededSA = true;
    }
    if( passdR ){
      l1obj.dRl1 = roidRext;
      (isSeededSA) ? l1obj.isPassed = 1 : l1obj.isPassed = 0;
      l1obj.eta = roi->eta();
      l1obj.phi = roi->phi();
      l1obj.thrValue = roi->thrValue();
      l1obj.roiNum = roi->getRoI();
      l1obj.isMoreCandInRoI = roi->isMoreCandInRoI();
      l1obj.thrNumber = l1thr;
      //how to calculate the roiSector : https://twiki.cern.ch/twiki/bin/view/Main/L1TGCNtuple#sectorAddress_8_bit_information
      uint32_t roiSec = roi->getSectorAddress() >> 1;
      (roiSec >= 64 ) ? (l1obj.roiSector = roiSec & 0x3f) : (l1obj.roiSector = roiSec & 0x1f);
      std::cout << "... matched" << std::endl;
    }
    else { std::cout << "... not matched" << std::endl; }
  }
  if( l1obj.isPassed > -1 ) return true;
  return false;
}

//---------------------------------------------------------//
//---------------------------------------------------------//

bool TrigMatchingToolMT::matchSA( std::string& trig,
                                  const L1Object l1obj,
                                  SAObject& saobj )
{
  std::cout << "start SA matching , L1 roiNum = " << l1obj.roiNum << std::endl;
  // get SA cont by using TrigDecTool
  std::vector< TrigCompositeUtils::LinkInfo<xAOD::L2StandAloneMuonContainer> > l2saLinks = m_trigDecTool->features<xAOD::L2StandAloneMuonContainer>( trig, TrigDefs::includeFailedDecisions );
  if(trig.find("l2mp") != std::string::npos){
    l2saLinks = m_trigDecTool->features<xAOD::L2StandAloneMuonContainer>( trig, TrigDefs::includeFailedDecisions, "HLT_MuonL2SAInfoMPmode" );
  }
  for(const TrigCompositeUtils::LinkInfo<xAOD::L2StandAloneMuonContainer>& l2saLinkInfo : l2saLinks){
    if( !l2saLinkInfo.isValid() ) continue;
    const ElementLink<xAOD::L2StandAloneMuonContainer> l2sa = l2saLinkInfo.link;
    if( !l2sa.isValid() ) continue;
    // start L1 -- SA matching
    std::cout << "... L2SA roiNum/pt/eta/phi/roiThr: " << (*l2sa)->roiNumber() << "/" << (*l2sa)->pt() << "/" << (*l2sa)->eta() << "/" << (*l2sa)->phi() << "/" << (*l2sa)->roiThreshold() << std::endl;
    std::cout << "... RoI's sector address: L1/SA = " << l1obj.roiSector << "/" << (*l2sa)->roiSector() << std::endl;
    if( (int)(*l2sa)->roiNumber() == l1obj.roiNum && (int)(*l2sa)->roiSector() == l1obj.roiSector ){
      ( l2saLinkInfo.state == TrigCompositeUtils::ActiveState::ACTIVE ) ? saobj.isPassed = 1 : saobj.isPassed = 0;
      saobj.roiNum = (*l2sa)->roiNumber();
      saobj.roiSector = (*l2sa)->roiSector();
      saobj.pt = (*l2sa)->pt();
      saobj.eta = (*l2sa)->eta();
      saobj.phi = (*l2sa)->phi();
      pair< double, double > extsa = m_ext.extTrackMuComb( *l2sa );
      saobj.etaBE = extsa.first;
      saobj.phiBE = extsa.second;
      saobj.etaMS = (*l2sa)->etaMS();
      saobj.phiMS = (*l2sa)->phiMS();
      saobj.sAddress = (*l2sa)->sAddress();
      //
      saobj.tgcPt  = (*l2sa)->tgcPt();
      saobj.ptBarrelRadius  = (*l2sa)->ptBarrelRadius();
      saobj.ptBarrelSagitta  = (*l2sa)->ptBarrelSagitta();
      saobj.ptEndcapAlpha  = (*l2sa)->ptEndcapAlpha();
      saobj.ptEndcapBeta  = (*l2sa)->ptEndcapBeta();
      saobj.ptEndcapRadius  = (*l2sa)->ptEndcapRadius();
      saobj.ptCSC = (*l2sa)->ptCSC();
      saobj.sAddress  = (*l2sa)->sAddress();

      saobj.roiEta  = (*l2sa)->roiEta();
      saobj.roiPhi  = (*l2sa)->roiPhi();
      saobj.isRpcFailure  = (*l2sa)->isRpcFailure();
      saobj.isTgcFailure  = (*l2sa)->isTgcFailure();
      //
      //float superPointR( int chamber ) const; ( defined in L2StandAloneMuon.h )
      //enum Chamber {BarrelInner = 0, BarrelMiddle = 1, BarrelOuter =2, EndcapInner =3, EndcapMiddle = 4, EndcapOuter = 5,EndcapExtra =6, CSC = 7, BEE = 8, BME = 9, Backup = 10, MaxChamber = 11} ( defined in TrigMuonDefs.h )
      saobj.superPointR_BI  = (*l2sa)->superPointR(0);
      saobj.superPointR_BM  = (*l2sa)->superPointR(1);
      saobj.superPointR_BO  = (*l2sa)->superPointR(2);
      saobj.superPointR_EI  = (*l2sa)->superPointR(3);
      saobj.superPointR_EM  = (*l2sa)->superPointR(4);
      saobj.superPointR_EO  = (*l2sa)->superPointR(5);
      saobj.superPointR_EE  = (*l2sa)->superPointR(6);
      saobj.superPointR_CSC = (*l2sa)->superPointR(7);
      saobj.superPointR_BEE = (*l2sa)->superPointR(8);
      saobj.superPointR_BME = (*l2sa)->superPointR(9);
      //float superPointZ( int chamber ) const; ( defined in L2StandAloneMuon.h )
      saobj.superPointZ_BI  = (*l2sa)->superPointZ(0);
      saobj.superPointZ_BM  = (*l2sa)->superPointZ(1);
      saobj.superPointZ_BO  = (*l2sa)->superPointZ(2);
      saobj.superPointZ_EI  = (*l2sa)->superPointZ(3);
      saobj.superPointZ_EM  = (*l2sa)->superPointZ(4);
      saobj.superPointZ_EO  = (*l2sa)->superPointZ(5);
      saobj.superPointZ_EE  = (*l2sa)->superPointZ(6);
      saobj.superPointZ_CSC = (*l2sa)->superPointZ(7);
      saobj.superPointZ_BEE = (*l2sa)->superPointZ(8);
      saobj.superPointZ_BME = (*l2sa)->superPointZ(9);
      //float superPointSlope( int chamber ) const; ( defined in L2StandAloneMuon.h )
      saobj.superPointSlope_BI  = (*l2sa)->superPointSlope(0);
      saobj.superPointSlope_BM  = (*l2sa)->superPointSlope(1);
      saobj.superPointSlope_BO  = (*l2sa)->superPointSlope(2);
      saobj.superPointSlope_EI  = (*l2sa)->superPointSlope(3);
      saobj.superPointSlope_EM  = (*l2sa)->superPointSlope(4);
      saobj.superPointSlope_EO  = (*l2sa)->superPointSlope(5);
      saobj.superPointSlope_EE  = (*l2sa)->superPointSlope(6);
      saobj.superPointSlope_CSC = (*l2sa)->superPointSlope(7);
      saobj.superPointSlope_BEE = (*l2sa)->superPointSlope(8);
      saobj.superPointSlope_BME = (*l2sa)->superPointSlope(9);
      //float superPointIntercept( int chamber ) const; ( defined in L2StandAloneMuon.h )
      saobj.superPointIntercept_BI  = (*l2sa)->superPointIntercept(0);
      saobj.superPointIntercept_BM  = (*l2sa)->superPointIntercept(1);
      saobj.superPointIntercept_BO  = (*l2sa)->superPointIntercept(2);
      saobj.superPointIntercept_EI  = (*l2sa)->superPointIntercept(3);
      saobj.superPointIntercept_EM  = (*l2sa)->superPointIntercept(4);
      saobj.superPointIntercept_EO  = (*l2sa)->superPointIntercept(5);
      saobj.superPointIntercept_EE  = (*l2sa)->superPointIntercept(6);
      saobj.superPointIntercept_CSC = (*l2sa)->superPointIntercept(7);
      saobj.superPointIntercept_BEE = (*l2sa)->superPointIntercept(8);
      saobj.superPointIntercept_BME = (*l2sa)->superPointIntercept(9);
      //float superPointChi2( int chamber ) const; ( defined in L2StandAloneMuon.h )
      saobj.superPointChi2_BI  = (*l2sa)->superPointChi2(0);
      saobj.superPointChi2_BM  = (*l2sa)->superPointChi2(1);
      saobj.superPointChi2_BO  = (*l2sa)->superPointChi2(2);
      saobj.superPointChi2_EI  = (*l2sa)->superPointChi2(3);
      saobj.superPointChi2_EM  = (*l2sa)->superPointChi2(4);
      saobj.superPointChi2_EO  = (*l2sa)->superPointChi2(5);
      saobj.superPointChi2_EE  = (*l2sa)->superPointChi2(6);
      saobj.superPointChi2_CSC = (*l2sa)->superPointChi2(7);
      saobj.superPointChi2_BEE = (*l2sa)->superPointChi2(8);
      saobj.superPointChi2_BME = (*l2sa)->superPointChi2(9);
      //
      saobj.rpcHitX = (*l2sa)->rpcHitX();
      saobj.rpcHitY = (*l2sa)->rpcHitY();
      saobj.rpcHitZ = (*l2sa)->rpcHitZ();
      saobj.rpcHitMeasPhi = (*l2sa)->rpcHitMeasuresPhi();
      saobj.rpcHitStationName = (*l2sa)->rpcHitStationName();
      saobj.rpcHitLayer = (*l2sa)->rpcHitLayer();
      for(int i = 0; i < (int)saobj.rpcHitX.size(); i++){
        TVector3 rpcHit(saobj.rpcHitX.at(i), saobj.rpcHitY.at(i), saobj.rpcHitZ.at(i));
        saobj.rpcHitEta.push_back(rpcHit.Eta());
        saobj.rpcHitPhi.push_back(rpcHit.Phi());
        saobj.rpcHitR.push_back(rpcHit.Perp());
      }
      saobj.tgcHitEta = (*l2sa)->tgcHitEta();
      saobj.tgcHitPhi = (*l2sa)->tgcHitPhi();
      saobj.tgcHitR = (*l2sa)->tgcHitEta();
      saobj.tgcHitZ = (*l2sa)->tgcHitZ();
      saobj.tgcHitWidth = (*l2sa)->tgcHitWidth();
      saobj.tgcHitStationNum = (*l2sa)->tgcHitStationNum();
      saobj.tgcHitIsStrip = (*l2sa)->tgcHitIsStrip();
      saobj.tgcHitBCTag = (*l2sa)->tgcHitBCTag();
      saobj.tgcHitInRoad = (*l2sa)->tgcHitInRoad();
      //
      for(int i = 0; i < (int)(*l2sa)->nMdtHits(); i++){
        saobj.mdtHitIsOutlier.push_back((*l2sa)->mdtHitIsOutlier(i));
        saobj.mdtHitChamber.push_back((*l2sa)->mdtHitChamber(i));
        saobj.mdtHitR.push_back((*l2sa)->mdtHitR(i));
        saobj.mdtHitZ.push_back((*l2sa)->mdtHitZ(i));
        saobj.mdtHitPhi.push_back((*l2sa)->mdtHitPhi(i));
        saobj.mdtHitResidual.push_back((*l2sa)->mdtHitResidual(i));
      }
      for(int i = 0; i < 9; /*<-chamber*/ i++){
        saobj.roadAw.push_back((*l2sa)->roadAw(i,0));
        saobj.roadBw.push_back((*l2sa)->roadBw(i,0));
        saobj.zMin.push_back((*l2sa)->zMin(i,0));
        saobj.zMax.push_back((*l2sa)->zMax(i,0));
        saobj.rMin.push_back((*l2sa)->rMin(i,0));
        saobj.rMax.push_back((*l2sa)->rMax(i,0));
        saobj.etaMin.push_back((*l2sa)->etaMin(i,0));
        saobj.etaMax.push_back((*l2sa)->etaMax(i,0));
      }
      saobj.stgcClusterR = (*l2sa)->stgcClusterR();
      saobj.stgcClusterZ = (*l2sa)->stgcClusterZ();
      saobj.stgcClusterEta = (*l2sa)->stgcClusterEta();
      saobj.stgcClusterPhi = (*l2sa)->stgcClusterPhi();
      saobj.stgcClusterResidualR = (*l2sa)->stgcClusterResidualR();
      saobj.stgcClusterResidualPhi = (*l2sa)->stgcClusterResidualPhi();
      saobj.stgcClusterStationEta = (*l2sa)->stgcClusterStationEta();
      saobj.stgcClusterStationPhi = (*l2sa)->stgcClusterStationPhi();
      saobj.stgcClusterStationName = (*l2sa)->stgcClusterStationName();
      saobj.stgcClusterType = (*l2sa)->stgcClusterType();
      saobj.stgcClusterIsOutlier = (*l2sa)->stgcClusterIsOutlier();
      saobj.stgcClusterLayer = (*l2sa)->stgcClusterLayer();
      saobj.mmClusterR = (*l2sa)->mmClusterR();
      saobj.mmClusterZ = (*l2sa)->mmClusterZ();
      saobj.mmClusterEta = (*l2sa)->mmClusterEta();
      saobj.mmClusterPhi = (*l2sa)->mmClusterPhi();
      saobj.mmClusterResidualR = (*l2sa)->mmClusterResidualR();
      saobj.mmClusterResidualPhi = (*l2sa)->mmClusterResidualPhi();
      saobj.mmClusterStationEta = (*l2sa)->mmClusterStationEta();
      saobj.mmClusterStationPhi = (*l2sa)->mmClusterStationPhi();
      saobj.mmClusterStationName = (*l2sa)->mmClusterStationName();
      saobj.mmClusterIsOutlier = (*l2sa)->mmClusterIsOutlier();
      saobj.mmClusterLayer = (*l2sa)->mmClusterLayer();
      std::cout << "... matched" << std::endl;

      break;
    }
    else {std::cout << "... not matched" << std::endl;}
  }
  if(saobj.isPassed > -1) return true;
  return false;
}

//---------------------------------------------------------//
//---------------------------------------------------------//

bool TrigMatchingToolMT::matchCB( std::string& trig,
                                  const SAObject saobj,
                                  CBObject& cbobj )
{
  std::cout << "start CB matching , SA roiNum/Sector = " << saobj.roiNum << "/" << saobj.roiSector << std::endl;
  // get CB cont by using TrigDecTool
  std::vector< TrigCompositeUtils::LinkInfo<xAOD::L2CombinedMuonContainer> > l2cbLinks;
  if(trig.find("l2io") != std::string::npos){
    l2cbLinks = m_trigDecTool->features<xAOD::L2CombinedMuonContainer>( trig, TrigDefs::includeFailedDecisions, "HLT_MuonL2CBInfoIOmode" );
  } 
  else if(trig.find("l2mp") != std::string::npos){
    l2cbLinks = m_trigDecTool->features<xAOD::L2CombinedMuonContainer>( trig, TrigDefs::includeFailedDecisions, "HLT_MuonL2CBInfoMPmode" );
  } else {
    l2cbLinks = m_trigDecTool->features<xAOD::L2CombinedMuonContainer>( trig, TrigDefs::includeFailedDecisions, "HLT_MuonL2CBInfo" );
  }
  //
  for(const TrigCompositeUtils::LinkInfo<xAOD::L2CombinedMuonContainer>& l2cbLinkInfo : l2cbLinks){
    if( !l2cbLinkInfo.isValid() ) continue;
    const ElementLink<xAOD::L2CombinedMuonContainer> l2cb = l2cbLinkInfo.link;
    if( !l2cb.isValid() ) continue;
    //if errorFlag = 1, then muFast pT is exact 0GeV and link to muFast can't be available.
    if( (*l2cb)->errorFlag() == 1){ 
      std::cout << "... L2CB: link to muFast is not available." << std::endl;
      std::cout << "... L2CB pt/eta/phi: " << (*l2cb)->pt() << "/" << (*l2cb)->eta() << "/" << (*l2cb)->phi() << std::endl;
    } else {
      std::cout << "... L2CB SAroiNum/pt/eta/phi: " << (*l2cb)->muSATrack()->roiNumber() << "/" << (*l2cb)->pt() << "/" << (*l2cb)->eta() << "/" << (*l2cb)->phi() << std::endl;
      if( (int)(*l2cb)->muSATrack()->roiNumber() == saobj.roiNum && (int)(*l2cb)->muSATrack()->roiSector() == saobj.roiSector ){ 
        cbobj.roiNum = (int)(*l2cb)->muSATrack()->roiNumber();
        cbobj.roiSector = (int)(*l2cb)->muSATrack()->roiSector();
        cbobj.pt = (*l2cb)->pt();
        cbobj.eta = (*l2cb)->eta();
        cbobj.phi = (*l2cb)->phi();
        ( l2cbLinkInfo.state == TrigCompositeUtils::ActiveState::ACTIVE ) ? cbobj.isPassed = 1 : cbobj.isPassed = 0;
        std::cout << "... matched" << std::endl;
        break;
      }
      else { std::cout << "... not matched" << std::endl;}
    }
  }
  if(cbobj.isPassed > -1) return true;
  return false;
}

//---------------------------------------------------------//
//---------------------------------------------------------//

bool TrigMatchingToolMT::matchEF( std::string& trig,
                                  const xAOD::Muon* muon,
                                  EFObject& efobj )
{
  std::cout << "start EF matching" << std::endl;
  bool match = false;
  // get EF cont by using TrigDecTool
  std::vector< TrigCompositeUtils::LinkInfo<xAOD::MuonContainer> > efLinks = m_trigDecTool->features<xAOD::MuonContainer>( trig, TrigDefs::includeFailedDecisions );
  //
  for(const TrigCompositeUtils::LinkInfo<xAOD::MuonContainer>& efLinkInfo : efLinks){
    if( !efLinkInfo.isValid() ) continue;
    const ElementLink<xAOD::MuonContainer> ef = efLinkInfo.link;
    if( !ef.isValid() ) continue;
    std::cout << "... ef pt/eta/phi: " << (*ef)->pt() << "/" << (*ef)->eta() << "/" << (*ef)->phi() << std::endl;
    if( efobj.dRef > m_utils.deltaR( muon->eta(), muon->phi(), (*ef)->eta(), (*ef)->phi() ) ){
      match = true;
      ( efLinkInfo.state == TrigCompositeUtils::ActiveState::ACTIVE ) ? efobj.isPassed = 1 : efobj.isPassed = 0;
      efobj.pt      = (*ef)->pt();
      efobj.eta     = (*ef)->eta();
      efobj.phi     = (*ef)->phi();
      efobj.dRef = m_utils.deltaR( muon->eta(), muon->phi(), (*ef)->eta(), (*ef)->phi() );
      std::cout << "... matched, dR = " << efobj.dRef << std::endl;
    }
    else{ std::cout << "... not matched" << std::endl; }
  }
  return match;
}

//---------------------------------------------------------//
//---------------------------------------------------------//
