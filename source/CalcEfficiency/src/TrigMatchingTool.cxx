#include "CalcEfficiency/TrigMatchingTool.h"
#include "TVector3.h"

TrigMatchingTool::TrigMatchingTool(){ }

//---------------------------------------------------------//
//---------------------------------------------------------//
//in the end, L1 requirement dRl1tag = minimum dR between l1RoI and offline
//if there are no RoIs fulfill the initial req, dRl1tag val isn't changed
bool TrigMatchingTool::matchL1( const Trig::Combination& comb,
                                const xAOD::Muon* muon, 
                                L1Object& l1obj) 
{
  std::cout << "start L1 matching" << std::endl;
  // extrapolate the muon's (eta, phi) to MS
  const xAOD::TrackParticle* mutrk = muon->trackParticle( xAOD::Muon::InnerDetectorTrackParticle );
  std::pair< double, double > extEtaAndPhi = m_ext.extTrack( mutrk );
  const double muExtEta = extEtaAndPhi.first;
  const double muExtPhi = extEtaAndPhi.second;
  // get L1 cont by using TrigRoiDescriptor 
  std::vector< Trig::Feature<TrigRoiDescriptor> > initRois = comb.elementFeature<TrigRoiDescriptorCollection>();
  if(initRois.size() < 1){
    std::cout << "L1Matching : No rois!" << std::endl;
    return false;
  }
  for(const auto &initRoi : initRois){
    auto itMu = m_trigDecTool->ancestor<xAOD::MuonRoI>(initRoi);
    const auto *l1 = itMu.cptr();
    std::cout << "... L1 RoI eta/phi/thre/roiNum = " << l1->eta() << "/" << l1->phi() << "/" << l1->thrValue() << "/" << l1->getRoI() << std::endl;
    const double l1eta = l1->eta();
    const double l1phi = l1->phi();
    const double roidRext = m_utils.deltaR( muExtEta, muExtPhi, l1eta, l1phi );
    std::cout << "... L1RoI vs offline mu dR = " << roidRext << std::endl;
    bool pass = (l1obj.dRl1 > roidRext);
    if( pass ){
      l1obj.dRl1 = roidRext;
      initRoi.te()->getActiveState() ? l1obj.isPassed = 1 : l1obj.isPassed = 0;
      l1obj.eta = l1->eta();
      l1obj.phi = l1->phi();
      l1obj.thrValue = l1->thrValue();
      l1obj.roiNum = l1->getRoI();
      uint32_t roiSec = l1->getSectorAddress() >> 1;
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
// in L1-SA matching, RoI number and RoISectorAddress matching is done.
// ===description of RoISectorAddress===
/// The sector address is an 8-bit identifier of the sector. For its detailed
/// description, see page 38 of https://edms.cern.ch/file/248757/1/mirod.pdf
///
/// The least significant bit shown which hemisphere the sector is in
/// (0: positive side, 1: negative side), the upper (1 or 2) bits show
/// what kind of sector it is, the rest of the address is the number of the
/// sector.
///
/// @return An 8-bit identifier
///=====================================

bool TrigMatchingTool::matchSA( const Trig::Combination& comb,
                                std::string& mesSATEName,
                                const L1Object l1obj,
                                SAObject& saobj )
{
  std::cout << "start SA matching , L1 roiNum = " << l1obj.roiNum << std::endl;
  std::cout << " #SATEname = " << mesSATEName << std::endl;
  // get SA cont by using TrigDecTool
  std::vector< Trig::Feature <xAOD::L2StandAloneMuonContainer> > l2samufeats = comb.containerFeature<xAOD::L2StandAloneMuonContainer>();
  //
  Trig::ExpertMethods* expert = m_trigDecTool -> ExperimentalAndExpertMethods();
  expert->enable();
  //
  for(const Trig::Feature<xAOD::L2StandAloneMuonContainer> mus:l2samufeats){
    //getTE and check active state
    bool isActiveTE = false;
    const HLT::TriggerElement *trigElement1 = mus.te();
    std::vector<HLT::TriggerElement*> TEsuccessors = expert->getNavigation()->getDirectSuccessors(trigElement1);
    for(auto te2 : TEsuccessors){
      std::cout << "  #TrigElementName=" << Trig::getTEName( *te2 ) << " -> " << te2->getActiveState() << std::endl;
      if ( te2 -> getActiveState() ){
        TString teName = Trig::getTEName( *te2 );
        if ( (teName.Contains( "L2_mu_SAhyp") && teName.Contains( mesSATEName.c_str() ) )|| ( teName.Contains("L2_mu_hypo1" ) ) )
          isActiveTE = true;
      }
    }
    // start L1 -- SA matching
    const xAOD::L2StandAloneMuonContainer *fL2SAs = mus.cptr();
    for(auto l2sa : *fL2SAs){
      std::cout << "... L2SA roiNum/pt/eta/phi/roiThr: " << l2sa->roiNumber() << "/" << l2sa->pt() << "/" << l2sa->eta() << "/" << l2sa->phi() << "/" << l2sa->roiThreshold() << std::endl;
      std::cout << "... RoI's sector address: L1/SA = " << l1obj.roiSector << "/" << l2sa->roiSector() << std::endl;
      int saroiNum = l2sa->roiNumber();
      int saroiSec = l2sa->roiSector();
      if( saroiNum == l1obj.roiNum && saroiSec == l1obj.roiSector ){
        saobj.roiNum = saroiNum;
        saobj.roiSector = saroiSec;
        saobj.pt = l2sa->pt();
        saobj.eta = l2sa->eta();
        saobj.phi = l2sa->phi();
        isActiveTE ? saobj.isPassed = 1 : saobj.isPassed = 0;
        std::cout << "... matched" << std::endl;

        break;
      }
      else {std::cout << "... not matched" << std::endl;}
    }
  }
  if(saobj.isPassed > -1) return true;
  return false;
}

//---------------------------------------------------------//
//---------------------------------------------------------//

bool TrigMatchingTool::matchCB( const Trig::Combination& comb,
                                std::string& mesCBTEName,
                                const SAObject saobj,
                                CBObject& cbobj)
{
  std::cout << "start CB matching , SA roiNum = " << saobj.roiNum << std::endl;
  std::cout << " #CBTEname = " << mesCBTEName << std::endl;
  // get SA cont by using TrigDecTool
  std::vector< Trig::Feature <xAOD::L2CombinedMuonContainer> > l2cbmufeats = comb.containerFeature<xAOD::L2CombinedMuonContainer>();
  //
  Trig::ExpertMethods* expert = m_trigDecTool -> ExperimentalAndExpertMethods();
  expert->enable();
  //
  for(const Trig::Feature<xAOD::L2CombinedMuonContainer> mus:l2cbmufeats){
    //getTE and check active state
    bool isActiveTE = false;
    const HLT::TriggerElement *trigElement1 = mus.te();
    std::vector<HLT::TriggerElement*> TEsuccessors = expert->getNavigation()->getDirectSuccessors(trigElement1);
    for(auto te2 : TEsuccessors){
      std::cout << "  #TrigElementName=" << Trig::getTEName( *te2 ) << " -> " << te2->getActiveState() << std::endl;
      if ( te2 -> getActiveState() ){
        TString teName = Trig::getTEName( *te2 );
        if ( (teName.Contains( "L2_mucombhyp") && teName.Contains( mesCBTEName.c_str() ) )|| ( teName.Contains("L2_mu_hypo2" ) ) )
          isActiveTE = true;
      }
    }
//    if(!isActiveTE) continue;
    const xAOD::L2CombinedMuonContainer *fL2CBs = mus.cptr();
    for( const auto& l2cb : *fL2CBs ) {
      std::cout << "... L2CB SAroiNum/pt/eta/phi: " << l2cb->muSATrack()->roiNumber() << "/" << l2cb->pt() << "/" << l2cb->eta() << "/" << l2cb->phi() << std::endl;
      if( (int)l2cb->muSATrack()->roiNumber() == saobj.roiNum && (int)l2cb->muSATrack()->roiSector() == saobj.roiSector ){ 
        cbobj.roiNum = l2cb->muSATrack()->roiNumber();
        cbobj.roiSector = l2cb->muSATrack()->roiSector();
        cbobj.pt = l2cb->pt();
        cbobj.eta = l2cb->eta();
        cbobj.phi = l2cb->phi();
        isActiveTE ? cbobj.isPassed = 1 : cbobj.isPassed = 0;
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

bool TrigMatchingTool::matchEFFS( const Trig::FeatureContainer& fc,
                                  const xAOD::Muon* muon,
                                  EFObject& efobj )
{
  std::cout << "start EF matching" << std::endl;
  bool match = false;
  const std::vector< Trig::Feature<xAOD::MuonContainer> > fEFs = fc.get<xAOD::MuonContainer>( "" );
  for ( auto& fEF : fEFs ){
    const HLT::TriggerElement* efTE = ( fEF.te() );
    const xAOD::MuonContainer* cont = fEF.cptr();
    for( const auto& ef : *cont ) {
      std::cout << "... EF pt/eta/phi: " << ef->pt() << "/" << ef->eta() << "/" << ef->phi() << std::endl;
      if( efobj.dRef > m_utils.deltaR( muon->eta(), muon->phi(), ef->eta(), ef->phi() ) ){
        match = true;
        efobj.pt      = ef->pt();
        efobj.eta     = ef->eta();
        efobj.phi     = ef->phi();
        efobj.dRef = m_utils.deltaR( muon->eta(), muon->phi(), ef->eta(), ef->phi() );
        efTE->getActiveState() ? efobj.isPassed = 1 : efobj.isPassed = 0;
        std::cout << "... matched, dR = " << efobj.dRef << std::endl;
      }
      else{ std::cout << "... not matched" << std::endl; }
    }
  }
  return match;
}

//---------------------------------------------------------//
//---------------------------------------------------------//

bool TrigMatchingTool::matchEF( const Trig::Combination& comb,
                                const xAOD::Muon* muon,
                                EFObject& efobj )
{
  std::cout << "start EF matching" << std::endl;
  bool match = false;
  const std::vector< Trig::Feature<xAOD::MuonContainer> > fEFs = comb.containerFeature<xAOD::MuonContainer>();
  for ( auto& fEF : fEFs ){
    const HLT::TriggerElement* efTE = ( fEF.te() );
    const xAOD::MuonContainer* cont = fEF.cptr();
    for( const auto& ef : *cont ) {
      std::cout << "... EF pt/eta/phi: " << ef->pt() << "/" << ef->eta() << "/" << ef->phi() << std::endl;
      if( efobj.dRef > m_utils.deltaR( muon->eta(), muon->phi(), ef->eta(), ef->phi() ) ){
        match = true;
        efTE->getActiveState() ? efobj.isPassed = 1 : efobj.isPassed = 0;
        efobj.pt      = ef->pt();
        efobj.eta     = ef->eta();
        efobj.phi     = ef->phi();
        efobj.dRef = m_utils.deltaR( muon->eta(), muon->phi(), ef->eta(), ef->phi() );
        std::cout << "... matched, dR = " << efobj.dRef << std::endl;
      }
      else{ std::cout << "... not matched" << std::endl; }
    }
  }
  return match;
}

//---------------------------------------------------------//
//---------------------------------------------------------//
//following match** function is for checking probe trigger status

bool TrigMatchingTool::matchL1( SG::ReadHandle<xAOD::MuonRoIContainer> &murois, 
                                const xAOD::Muon* muon, 
                                L1Object& l1obj,
                                std::string& trig ) 
{
  std::cout << "start L1 matching" << std::endl;
  // extrapolate the muon's (eta, phi) to MS
  const xAOD::TrackParticle* mutrk = muon->trackParticle( xAOD::Muon::InnerDetectorTrackParticle );
  std::pair< double, double > extEtaAndPhi = m_ext.extTrack( mutrk );
  const double muExtEta = extEtaAndPhi.first;
  const double muExtPhi = extEtaAndPhi.second;
//  std::vector< Trig::Feature<TrigRoiDescriptor> > initRois = fc.elementFeature<TrigRoiDescriptorCollection>("initialRoI");
//  if(initRois.size() < 1){
//    std::cout << "L1Matching : No rois!" << std::endl;
//    return false;
//  }
  for(const auto &l1 : *murois){
    std::cout << "... L1 RoI eta/phi/thre/roiNum = " << l1->eta() << "/" << l1->phi() << "/" << l1->thrValue() << "/" << l1->getRoI() << std::endl;
    const double l1eta = l1->eta();
    const double l1phi = l1->phi();
    const int l1thr = l1->getThrNumber();
    const double roidRext = m_utils.deltaR( muExtEta, muExtPhi, l1eta, l1phi );
    std::cout << "... L1RoI vs offline mu dR = " << roidRext << std::endl;
    bool passdR = (l1obj.dRl1 > roidRext);
    bool passThr = ( l1thr >= L1trigThr(trig) );
    if( passdR ){
      l1obj.dRl1 = roidRext;
      passThr ? l1obj.isPassed = 1 : l1obj.isPassed = 0;
      l1obj.eta = l1->eta();
      l1obj.phi = l1->phi();
      l1obj.thrValue = l1->thrValue();
      l1obj.roiNum = l1->getRoI();
      l1obj.thrNumber = l1thr;
      uint32_t roiSec = l1->getSectorAddress() >> 1;
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

bool TrigMatchingTool::matchSA( const Trig::FeatureContainer& fc,
                                std::string& mesSATEName,
                                const L1Object l1obj,
                                SAObject& saobj )
{
  std::cout << "start SA matching , L1 roiNum = " << l1obj.roiNum << std::endl;
  std::cout << " #SATEname = " << mesSATEName << std::endl;
  // get SA cont by using TrigDecTool
  std::vector< Trig::Feature <xAOD::L2StandAloneMuonContainer> > l2samufeats = fc.get<xAOD::L2StandAloneMuonContainer>("", TrigDefs::alsoDeactivateTEs);
  //
  Trig::ExpertMethods* expert = m_trigDecTool -> ExperimentalAndExpertMethods();
  expert->enable();
  //
  for(const Trig::Feature<xAOD::L2StandAloneMuonContainer> mus:l2samufeats){
    //getTE and check active state
    bool isActiveTE = false;
    const HLT::TriggerElement *trigElement1 = mus.te();
    std::vector<HLT::TriggerElement*> TEsuccessors = expert->getNavigation()->getDirectSuccessors(trigElement1);
    for(auto te2 : TEsuccessors){
      std::cout << "  #TrigElementName=" << Trig::getTEName( *te2 ) << " -> " << te2->getActiveState() << std::endl;
      if ( te2 -> getActiveState() ){
        TString teName = Trig::getTEName( *te2 );
        if ( (teName.Contains( "L2_mu_SAhyp") && teName.Contains( mesSATEName.c_str() ) )|| ( teName.Contains("L2_mu_hypo1" ) ) )
          isActiveTE = true;
      }
    }
    // start L1 -- SA matching
    const xAOD::L2StandAloneMuonContainer *fL2SAs = mus.cptr();
    for(auto l2sa : *fL2SAs){
      std::cout << "... L2SA roiNum/pt/eta/phi/roiThr: " << l2sa->roiNumber() << "/" << l2sa->pt() << "/" << l2sa->eta() << "/" << l2sa->phi() << "/" << l2sa->roiThreshold() << std::endl;
      std::cout << "... RoI's sector address: L1/SA = " << l1obj.roiSector << "/" << l2sa->roiSector() << std::endl;
      if( (int)l2sa->roiNumber() == l1obj.roiNum && (int)l2sa->roiSector() == l1obj.roiSector ){
        isActiveTE ? saobj.isPassed = 1 : saobj.isPassed = 0;
        saobj.roiNum = l2sa->roiNumber();
        saobj.roiSector = l2sa->roiSector();
        saobj.pt = l2sa->pt();
        saobj.eta = l2sa->eta();
        saobj.phi = l2sa->phi();
        pair< double, double > extsa = m_ext.extTrackMuComb( l2sa );
        saobj.etaBE = extsa.first;
        saobj.phiBE = extsa.second;
        saobj.etaMS = l2sa->etaMS();
        saobj.phiMS = l2sa->phiMS();
        //
        saobj.tgcPt  = l2sa->tgcPt();
        saobj.ptBarrelRadius  = l2sa->ptBarrelRadius();
        saobj.ptBarrelSagitta  = l2sa->ptBarrelSagitta();
        saobj.ptEndcapAlpha  = l2sa->ptEndcapAlpha();
        saobj.ptEndcapBeta  = l2sa->ptEndcapBeta();
        saobj.ptEndcapRadius  = l2sa->ptEndcapRadius();
        saobj.ptCSC = l2sa->ptCSC();
        saobj.sAddress  = l2sa->sAddress();

        saobj.roiEta  = l2sa->roiEta();
        saobj.roiPhi  = l2sa->roiPhi();
        saobj.isRpcFailure  = l2sa->isRpcFailure();
        saobj.isTgcFailure  = l2sa->isTgcFailure();
        //
        //float superPointR( int chamber ) const; ( defined in L2StandAloneMuon.h )
        //enum Chamber {BarrelInner = 0, BarrelMiddle = 1, BarrelOuter =2, EndcapInner =3, EndcapMiddle = 4, EndcapOuter = 5,EndcapExtra =6, CSC = 7, BEE = 8, BME = 9, Backup = 10, MaxChamber = 11} ( defined in TrigMuonDefs.h )
        saobj.superPointR_BI  = l2sa->superPointR(0);
        saobj.superPointR_BM  = l2sa->superPointR(1);
        saobj.superPointR_BO  = l2sa->superPointR(2);
        saobj.superPointR_EI  = l2sa->superPointR(3);
        saobj.superPointR_EM  = l2sa->superPointR(4);
        saobj.superPointR_EO  = l2sa->superPointR(5);
        saobj.superPointR_EE  = l2sa->superPointR(6);
        saobj.superPointR_CSC = l2sa->superPointR(7);
        saobj.superPointR_BEE = l2sa->superPointR(8);
        saobj.superPointR_BME = l2sa->superPointR(9);
        //float superPointZ( int chamber ) const; ( defined in L2StandAloneMuon.h )
        saobj.superPointZ_BI  = l2sa->superPointZ(0);
        saobj.superPointZ_BM  = l2sa->superPointZ(1);
        saobj.superPointZ_BO  = l2sa->superPointZ(2);
        saobj.superPointZ_EI  = l2sa->superPointZ(3);
        saobj.superPointZ_EM  = l2sa->superPointZ(4);
        saobj.superPointZ_EO  = l2sa->superPointZ(5);
        saobj.superPointZ_EE  = l2sa->superPointZ(6);
        saobj.superPointZ_CSC = l2sa->superPointZ(7);
        saobj.superPointZ_BEE = l2sa->superPointZ(8);
        saobj.superPointZ_BME = l2sa->superPointZ(9);
        //float superPointSlope( int chamber ) const; ( defined in L2StandAloneMuon.h )
        saobj.superPointSlope_BI  = l2sa->superPointSlope(0);
        saobj.superPointSlope_BM  = l2sa->superPointSlope(1);
        saobj.superPointSlope_BO  = l2sa->superPointSlope(2);
        saobj.superPointSlope_EI  = l2sa->superPointSlope(3);
        saobj.superPointSlope_EM  = l2sa->superPointSlope(4);
        saobj.superPointSlope_EO  = l2sa->superPointSlope(5);
        saobj.superPointSlope_EE  = l2sa->superPointSlope(6);
        saobj.superPointSlope_CSC = l2sa->superPointSlope(7);
        saobj.superPointSlope_BEE = l2sa->superPointSlope(8);
        saobj.superPointSlope_BME = l2sa->superPointSlope(9);
        //float superPointIntercept( int chamber ) const; ( defined in L2StandAloneMuon.h )
        saobj.superPointIntercept_BI  = l2sa->superPointIntercept(0);
        saobj.superPointIntercept_BM  = l2sa->superPointIntercept(1);
        saobj.superPointIntercept_BO  = l2sa->superPointIntercept(2);
        saobj.superPointIntercept_EI  = l2sa->superPointIntercept(3);
        saobj.superPointIntercept_EM  = l2sa->superPointIntercept(4);
        saobj.superPointIntercept_EO  = l2sa->superPointIntercept(5);
        saobj.superPointIntercept_EE  = l2sa->superPointIntercept(6);
        saobj.superPointIntercept_CSC = l2sa->superPointIntercept(7);
        saobj.superPointIntercept_BEE = l2sa->superPointIntercept(8);
        saobj.superPointIntercept_BME = l2sa->superPointIntercept(9);
        //float superPointChi2( int chamber ) const; ( defined in L2StandAloneMuon.h )
        saobj.superPointChi2_BI  = l2sa->superPointChi2(0);
        saobj.superPointChi2_BM  = l2sa->superPointChi2(1);
        saobj.superPointChi2_BO  = l2sa->superPointChi2(2);
        saobj.superPointChi2_EI  = l2sa->superPointChi2(3);
        saobj.superPointChi2_EM  = l2sa->superPointChi2(4);
        saobj.superPointChi2_EO  = l2sa->superPointChi2(5);
        saobj.superPointChi2_EE  = l2sa->superPointChi2(6);
        saobj.superPointChi2_CSC = l2sa->superPointChi2(7);
        saobj.superPointChi2_BEE = l2sa->superPointChi2(8);
        saobj.superPointChi2_BME = l2sa->superPointChi2(9);
        //
        saobj.rpcHitX = l2sa->rpcHitX();
        saobj.rpcHitY = l2sa->rpcHitY();
        saobj.rpcHitZ = l2sa->rpcHitZ();
        saobj.rpcHitMeasPhi = l2sa->rpcHitMeasuresPhi();
        saobj.rpcHitStationName = l2sa->rpcHitStationName();
        saobj.rpcHitLayer = l2sa->rpcHitLayer();
        for(int i = 0; i < (int)saobj.rpcHitX.size(); i++){
          TVector3 rpcHit(saobj.rpcHitX.at(i), saobj.rpcHitY.at(i), saobj.rpcHitZ.at(i));
          saobj.rpcHitEta.push_back(rpcHit.Eta());
          saobj.rpcHitPhi.push_back(rpcHit.Phi());
          saobj.rpcHitR.push_back(rpcHit.Perp());
        }
        saobj.tgcHitEta = l2sa->tgcHitEta();
        saobj.tgcHitPhi = l2sa->tgcHitPhi();
        saobj.tgcHitR = l2sa->tgcHitEta();
        saobj.tgcHitZ = l2sa->tgcHitZ();
        saobj.tgcHitWidth = l2sa->tgcHitWidth();
        saobj.tgcHitStationNum = l2sa->tgcHitStationNum();
        saobj.tgcHitIsStrip = l2sa->tgcHitIsStrip();
        saobj.tgcHitBCTag = l2sa->tgcHitBCTag();
        saobj.tgcHitInRoad = l2sa->tgcHitInRoad();
        //
        for(int i = 0; i < (int)l2sa->nMdtHits(); i++){
          saobj.mdtHitIsOutlier.push_back(l2sa->mdtHitIsOutlier(i));
          saobj.mdtHitChamber.push_back(l2sa->mdtHitChamber(i));
          saobj.mdtHitR.push_back(l2sa->mdtHitR(i));
          saobj.mdtHitZ.push_back(l2sa->mdtHitZ(i));
          saobj.mdtHitPhi.push_back(l2sa->mdtHitPhi(i));
          saobj.mdtHitResidual.push_back(l2sa->mdtHitResidual(i));
        }
        for(int i = 0; i < 9; /*<-chamber*/ i++){
          saobj.roadAw.push_back(l2sa->roadAw(i,0));
          saobj.roadBw.push_back(l2sa->roadBw(i,0));
          saobj.zMin.push_back(l2sa->zMin(i,0));
          saobj.zMax.push_back(l2sa->zMax(i,0));
          saobj.rMin.push_back(l2sa->rMin(i,0));
          saobj.rMax.push_back(l2sa->rMax(i,0));
          saobj.etaMin.push_back(l2sa->etaMin(i,0));
          saobj.etaMax.push_back(l2sa->etaMax(i,0));
        }
        std::cout << "... matched" << std::endl;

        break;
      }
      else {std::cout << "... not matched" << std::endl;}
    }
  }
  if(saobj.isPassed > -1) return true;
  return false;
}

//---------------------------------------------------------//
//---------------------------------------------------------//

bool TrigMatchingTool::matchCB( const Trig::FeatureContainer& fc,
                                std::string& mesCBTEName,
                                const SAObject saobj,
                                CBObject& cbobj)
{
  std::cout << "start CB matching , SA roiNum/Sector = " << saobj.roiNum << "/" << saobj.roiSector << std::endl;
  std::cout << " #CBTEname = " << mesCBTEName << std::endl;
  // get SA cont by using TrigDecTool
  std::vector< Trig::Feature <xAOD::L2CombinedMuonContainer> > l2cbmufeats = fc.get<xAOD::L2CombinedMuonContainer>("", TrigDefs::alsoDeactivateTEs);
  //
  Trig::ExpertMethods* expert = m_trigDecTool -> ExperimentalAndExpertMethods();
  expert->enable();
  //
  for(const Trig::Feature<xAOD::L2CombinedMuonContainer> mus:l2cbmufeats){
    //getTE and check active state
    bool isActiveTE = false;
    const HLT::TriggerElement *trigElement1 = mus.te();
    std::vector<HLT::TriggerElement*> TEsuccessors = expert->getNavigation()->getDirectSuccessors(trigElement1);
    for(auto te2 : TEsuccessors){
      std::cout << "  #TrigElementName=" << Trig::getTEName( *te2 ) << " -> " << te2->getActiveState() << std::endl;
      if ( te2 -> getActiveState() ){
        TString teName = Trig::getTEName( *te2 );
        if ( (teName.Contains( "L2_mucombhyp") && teName.Contains( mesCBTEName.c_str() ) )|| ( teName.Contains("L2_mu_hypo2" ) ) )
          isActiveTE = true;
      }
    }
//    if(!isActiveTE) continue;
    const xAOD::L2CombinedMuonContainer *fL2CBs = mus.cptr();
    for( const auto& l2cb : *fL2CBs ) {
      std::cout << "... L2CB SAroiNum/pt/eta/phi: " << l2cb->muSATrack()->roiNumber() << "/" << l2cb->pt() << "/" << l2cb->eta() << "/" << l2cb->phi() << std::endl;
      if( (int)l2cb->muSATrack()->roiNumber() == saobj.roiNum && (int)l2cb->muSATrack()->roiSector() == saobj.roiSector ){ 
        cbobj.roiNum = l2cb->muSATrack()->roiNumber();
        cbobj.roiSector = l2cb->muSATrack()->roiSector();
        cbobj.pt = l2cb->pt();
        cbobj.eta = l2cb->eta();
        cbobj.phi = l2cb->phi();
        isActiveTE ? cbobj.isPassed = 1 : cbobj.isPassed = 0;
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

bool TrigMatchingTool::matchEF( const Trig::FeatureContainer& fc,
                                const xAOD::Muon* muon,
                                EFObject& efobj )
{
  std::cout << "start EF matching" << std::endl;
  bool match = false;
  const std::vector< Trig::Feature<xAOD::MuonContainer> > fEFs = fc.get<xAOD::MuonContainer>("");
  for ( auto& fEF : fEFs ){
    const HLT::TriggerElement* efTE = ( fEF.te() );
    const xAOD::MuonContainer* cont = fEF.cptr();
    for( const auto& ef : *cont ) {
      std::cout << "... EF pt/eta/phi: " << ef->pt() << "/" << ef->eta() << "/" << ef->phi() << std::endl;
      if( efobj.dRef > m_utils.deltaR( muon->eta(), muon->phi(), ef->eta(), ef->phi() ) ){
        match = true;
        efTE->getActiveState() ? efobj.isPassed = 1 : efobj.isPassed = 0;
        efobj.pt      = ef->pt();
        efobj.eta     = ef->eta();
        efobj.phi     = ef->phi();
        efobj.dRef = m_utils.deltaR( muon->eta(), muon->phi(), ef->eta(), ef->phi() );
        std::cout << "... matched, dR = " << efobj.dRef << std::endl;
      }
      else{ std::cout << "... not matched" << std::endl; }
    }
  }
  return match;
}

