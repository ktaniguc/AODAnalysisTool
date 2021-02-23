#include "CalcEfficiency/NtupleMakerTool.h"
#include "TVector3.h"

NtupleMakerTool::NtupleMakerTool(){ }

void NtupleMakerTool::setTriggerName(std::string l1, std::string hlt)
{
  if(hlt.find("b") != std::string::npos) 
    std::cout << "WARNING         this chain maybe 'nomucomb' trigger. if so, isPassedCB are not correct" << std::endl;
  m_hltchainName.push_back(hlt);
  m_l1chainName.push_back(l1);
}

void NtupleMakerTool::withdrawInfo(SG::ReadHandle<xAOD::MuonContainer>       &muons,
                                   SG::ReadHandle<xAOD::MuonRoIContainer>    &rois,
                                   SG::ReadHandle<xAOD::L2StandAloneMuonContainer> &l2sa,
                                   SG::ReadHandle<xAOD::L2StandAloneMuonContainer> &l2saio,
                                   SG::ReadHandle<xAOD::L2CombinedMuonContainer> &l2cb,
                                   SG::ReadHandle<xAOD::L2CombinedMuonContainer> &l2cbio){

  //======== offline =======
  for( const auto& muon : *muons){
    const xAOD::TrackParticle* muontrk = muon->trackParticle( xAOD::Muon::InnerDetectorTrackParticle );
    if(!muontrk) continue;
    std::pair< double, double > muonextEtaAndPhi = m_ext.extTrack( muontrk );
    OfflineObject muonobj;
    muonobj.pt = muon->pt();
    muonobj.eta = muon->eta();
    muonobj.phi = muon->phi();
    muonobj.extEta = muonextEtaAndPhi.first;
    muonobj.extPhi = muonextEtaAndPhi.second;
    m_muons.push_back(muonobj);
  }
  //======== L1 =======
  for(const auto& roi : *rois){
    L1Object l1obj;
    std::cout << "... L1RoI roiNum/thrNum/eta/phi/roiWord = " << roi->getRoI() << "/" << roi->getThrNumber() << "/" << roi->eta() << "/" << roi->phi() << "/" << roi->roiWord() << std::endl;
    l1obj.eta = roi->eta();
    l1obj.phi = roi->phi();
    l1obj.thrValue = roi->thrValue();
    l1obj.roiNum = roi->getRoI();
    l1obj.isMoreCandInRoI = roi->isMoreCandInRoI();
    l1obj.thrNumber = roi->getThrNumber();
    //how to calculate the roiSector : https://twiki.cern.ch/twiki/bin/view/Main/L1TGCNtuple#sectorAddress_8_bit_information
    uint32_t roiSec = roi->getSectorAddress() >> 1;
    (roiSec >= 64 ) ? (l1obj.roiSector = roiSec & 0x3f) : (l1obj.roiSector = roiSec & 0x1f);
    m_L1objects.push_back(l1obj);
  }
  //======== default L2SA =======
  for(const auto& sa : *l2sa){
    SAObject saobj;
    std::cout << "... HLT_MuonL2SAInfo roiNum/pt/etaMS/phiMS/roiEta/roiPhi = " << sa->roiNumber() << "/" << sa->pt() << "/" << sa->etaMS() << "/" << sa->phiMS() << "/" << sa->roiEta() << "/" << sa->roiPhi() << std::endl;
    saobj.roiNum = sa->roiNumber();
    saobj.roiSector = sa->roiSector();
    saobj.pt = sa->pt();
    saobj.eta = sa->eta();
    saobj.phi = sa->phi();
    saobj.etaMS = sa->etaMS();
    saobj.phiMS = sa->phiMS();
    saobj.sAddress = sa->sAddress();
    saobj.roiEta  = sa->roiEta();
    saobj.roiPhi  = sa->roiPhi();
    saobj.superPointR_BI  = sa->superPointR(0);
    saobj.superPointR_BM  = sa->superPointR(1);
    saobj.superPointR_BO  = sa->superPointR(2);
    saobj.superPointR_EI  = sa->superPointR(3);
    saobj.superPointR_EM  = sa->superPointR(4);
    saobj.superPointR_EO  = sa->superPointR(5);
    saobj.superPointR_EE  = sa->superPointR(6);
    saobj.superPointR_CSC = sa->superPointR(7);
    saobj.superPointR_BEE = sa->superPointR(8);
    saobj.superPointR_BME = sa->superPointR(9);
    //float superPointZ int chamber ) const;  defined in L2StandAloneMuon.h )
    saobj.superPointZ_BI  = sa->superPointZ(0);
    saobj.superPointZ_BM  = sa->superPointZ(1);
    saobj.superPointZ_BO  = sa->superPointZ(2);
    saobj.superPointZ_EI  = sa->superPointZ(3);
    saobj.superPointZ_EM  = sa->superPointZ(4);
    saobj.superPointZ_EO  = sa->superPointZ(5);
    saobj.superPointZ_EE  = sa->superPointZ(6);
    saobj.superPointZ_CSC = sa->superPointZ(7);
    saobj.superPointZ_BEE = sa->superPointZ(8);
    saobj.superPointZ_BME = sa->superPointZ(9);
    m_SAobjects.push_back(saobj);
  }
  //======== L2SA IOmode =======
  for(const auto& saio : *l2saio){
    SAObject saobj;
    std::cout << "... HLT_MuonL2SAInfoIOmode roiNum/pt/etaMS/phiMS/roiEta/roiPhi = " << saio->roiNumber() << "/" << saio->pt() << "/" << saio->etaMS() << "/" << saio->phiMS() << "/" << saio->roiEta() << "/" << saio->roiPhi() << std::endl;
    saobj.roiNum = saio->roiNumber();
    saobj.roiSector = saio->roiSector();
    saobj.pt = saio->pt();
    saobj.eta = saio->eta();
    saobj.phi = saio->phi();
    saobj.etaMS = saio->etaMS();
    saobj.phiMS = saio->phiMS();
    saobj.sAddress = saio->sAddress();
    saobj.roiEta  = saio->roiEta();
    saobj.roiPhi  = saio->roiPhi();
    saobj.superPointR_BI  = saio->superPointR(0);
    saobj.superPointR_BM  = saio->superPointR(1);
    saobj.superPointR_BO  = saio->superPointR(2);
    saobj.superPointR_EI  = saio->superPointR(3);
    saobj.superPointR_EM  = saio->superPointR(4);
    saobj.superPointR_EO  = saio->superPointR(5);
    saobj.superPointR_EE  = saio->superPointR(6);
    saobj.superPointR_CSC = saio->superPointR(7);
    saobj.superPointR_BEE = saio->superPointR(8);
    saobj.superPointR_BME = saio->superPointR(9);
    //float superPointZ int chamber ) const;  defined in L2StandAloneMuon.h )
    saobj.superPointZ_BI  = saio->superPointZ(0);
    saobj.superPointZ_BM  = saio->superPointZ(1);
    saobj.superPointZ_BO  = saio->superPointZ(2);
    saobj.superPointZ_EI  = saio->superPointZ(3);
    saobj.superPointZ_EM  = saio->superPointZ(4);
    saobj.superPointZ_EO  = saio->superPointZ(5);
    saobj.superPointZ_EE  = saio->superPointZ(6);
    saobj.superPointZ_CSC = saio->superPointZ(7);
    saobj.superPointZ_BEE = saio->superPointZ(8);
    saobj.superPointZ_BME = saio->superPointZ(9);
    m_SAIOobjects.push_back(saobj);
  }
  //======== default CB =======
  for(const auto& cb : *l2cb){
    CBObject cbobj;
    if((cb)->errorFlag() != 1){
      std::cout << "... HLT_MuonL2CBInfo roiNum/pt/eta/phi = " << cb->muSATrack()->roiNumber() << "/" << cb->pt() << "/" << cb->eta() << "/" << cb->phi() << std::endl;
      cbobj.roiNum = (int)cb->muSATrack()->roiNumber();
      cbobj.roiSector = (int)cb->muSATrack()->roiSector();
    } else {
      std::cout << "... HLT_MuonL2CBInfo pt/eta/phi = " << cb->pt() << "/" << cb->eta() << "/" << cb->phi() << std::endl;
    }
    if((cb)->errorFlag() == 0){
      cbobj.idpt = (cb)->idTrack()->pt();
      cbobj.ideta = (cb)->idTrack()->eta();
      cbobj.idphi = (cb)->idTrack()->phi();
    }
    cbobj.pt = cb->pt();
    cbobj.eta = cb->eta();
    cbobj.phi = cb->phi();
    m_CBobjects.push_back(cbobj);
  }
  //======== CB IOmode =======
  for(const auto& cbio : *l2cbio){
    CBObject cbobj;
    if((cbio)->errorFlag() != 1){
      std::cout << "... HLT_MuonL2CBInfoMPmode roiNum/pt/eta/phi = " << cbio->muSATrack()->roiNumber() << "/" << cbio->pt() << "/" << cbio->eta() << "/" << cbio->phi() << std::endl;
      cbobj.roiNum = (int)cbio->muSATrack()->roiNumber();
      cbobj.roiSector = (int)cbio->muSATrack()->roiSector();
    } else {
      std::cout << "... HLT_MuonL2CBInfoMPmode pt/eta/phi = " << cbio->pt() << "/" << cbio->eta() << "/" << cbio->phi() << std::endl;
    }
    if((cbio)->errorFlag() == 0){
      cbobj.idpt = (cbio)->idTrack()->pt();
      cbobj.ideta = (cbio)->idTrack()->eta();
      cbobj.idphi = (cbio)->idTrack()->phi();
    }
    cbobj.pt = cbio->pt();
    cbobj.eta = cbio->eta();
    cbobj.phi = cbio->phi();
    m_CBIOobjects.push_back(cbobj);
  }
  //========= TDT info =========
  std::cout << "----------- start retrieving trigger objects from tdt -----------" << std::endl;
  for(int i_trig=0; i_trig< (int)m_hltchainName.size(); i_trig++){
    TDTObject tdtobj;
    tdtobj.trigChainName = m_hltchainName.at(i_trig);
    auto cg = m_trigDecTool->getChainGroup(m_hltchainName.at(i_trig));
    tdtobj.isPassedChain = cg->isPassed( TrigDefs::eventAccepted );
    std::cout << m_hltchainName.at(i_trig) << " is Passed? ===> " << tdtobj.isPassedChain << std::endl;
    //L1
    // by event
    tdtobj.isPassedL1_evt = m_trigDecTool->isPassed(m_l1chainName.at(i_trig), TrigDefs::eventAccepted);
    std::cout << m_l1chainName.at(i_trig) << " is Passed ===> " << tdtobj.isPassedL1_evt << std::endl;

    //get initialRoIs seeded by L2MuonSA for judging if passing L1 trigger or not
    std::vector<uint32_t> link_roiWords; link_roiWords.clear();
    std::vector< TrigCompositeUtils::LinkInfo<xAOD::L2StandAloneMuonContainer> > l2saLinks_forL1 = m_trigDecTool->features<xAOD::L2StandAloneMuonContainer>( m_hltchainName.at(i_trig), TrigDefs::includeFailedDecisions );
    for(const TrigCompositeUtils::LinkInfo<xAOD::L2StandAloneMuonContainer>& l2saLinkInfo : l2saLinks_forL1){
      if( !l2saLinkInfo.isValid() ) continue;
      const ElementLink<xAOD::L2StandAloneMuonContainer> l2sa = l2saLinkInfo.link;
      if( !l2sa.isValid() ) continue;
      const TrigCompositeUtils::Decision* muDecision = l2saLinkInfo.source; 
      const TrigCompositeUtils::LinkInfo<TrigRoiDescriptorCollection> roiLinkInfo = TrigCompositeUtils::findLink<TrigRoiDescriptorCollection>(muDecision, "initialRoI");
      if( !roiLinkInfo.isValid() ) continue;
      const ElementLink<TrigRoiDescriptorCollection> roiEL = roiLinkInfo.link;
      link_roiWords.push_back((*roiEL)->roiWord());
      //std::cout << "[TrigRoiDescriptor] L1RoIs seeded by muFast: roiEta/roiPhi/roiWord = " << (*roiEL)->eta() << "/" << (*roiEL)->phi() << "/" << (*roiEL)->roiWord() << std::endl;
    }
    // by object
    for(const auto& roi : *rois){
      L1Object l1obj;
      bool isSeededSA = false;
      for(int i_link = 0; i_link < (int)link_roiWords.size(); i_link++){
        if(roi->roiWord() == link_roiWords.at(i_link)) isSeededSA = true;
      }
      uint32_t roiSec = roi->getSectorAddress() >> 1;
      int roiSector = -999;
      (roiSec >= 64 ) ? (roiSector = roiSec & 0x3f) : (roiSector = roiSec & 0x1f);
      l1obj.isPassed = isSeededSA;
      l1obj.roiNum = roi->getRoI();
      l1obj.roiSector = roiSector;
      l1obj.isMoreCandInRoI = roi->isMoreCandInRoI();
      tdtobj.m_l1.push_back(l1obj);
      std::cout << "L1 w/ TDT: roiNumber/roiSector/isPassed = " << l1obj.roiNum << "/" << l1obj.roiSector << "/" << isSeededSA << std::endl;
    }
    //default SA
    std::vector< TrigCompositeUtils::LinkInfo<xAOD::L2StandAloneMuonContainer> > l2saLinks = m_trigDecTool->features<xAOD::L2StandAloneMuonContainer>( m_hltchainName.at(i_trig), TrigDefs::includeFailedDecisions, "HLT_MuonL2SAInfo" );
    for(const TrigCompositeUtils::LinkInfo<xAOD::L2StandAloneMuonContainer>& l2saLinkInfo : l2saLinks){
      if( !l2saLinkInfo.isValid() ) continue;
      const ElementLink<xAOD::L2StandAloneMuonContainer> l2sa = l2saLinkInfo.link;
      if( !l2sa.isValid() ) continue;
      if( l2saLinkInfo.state == TrigCompositeUtils::ActiveState::ACTIVE ) tdtobj.isPassedSA_evt = true;

      SAObject saobj;
      bool ispassed_sa = false;
      if(l2saLinkInfo.state == TrigCompositeUtils::ActiveState::ACTIVE) ispassed_sa = true;
      saobj.isPassed = ispassed_sa;
      saobj.roiNum = (*l2sa)->roiNumber();
      saobj.roiSector = (*l2sa)->roiSector();
      saobj.pt = (*l2sa)->pt();
      saobj.eta = (*l2sa)->eta();
      saobj.phi = (*l2sa)->phi();
      saobj.etaMS = (*l2sa)->etaMS();
      saobj.phiMS = (*l2sa)->phiMS();
      saobj.sAddress = (*l2sa)->sAddress();
      saobj.roiEta  = (*l2sa)->roiEta();
      saobj.roiPhi  = (*l2sa)->roiPhi();
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
      std::cout << "SA default w/ TDT: pt/eta/phi/isPassed = " << saobj.pt << "/" << saobj.eta << "/"
        << saobj.phi << "/" << ispassed_sa << std::endl;
      tdtobj.m_sa.push_back(saobj);
    }
    //SAIO
    if(m_isvalidIOmode){
      std::vector< TrigCompositeUtils::LinkInfo<xAOD::L2StandAloneMuonContainer> > l2saLinksIO = m_trigDecTool->features<xAOD::L2StandAloneMuonContainer>( m_hltchainName.at(i_trig), TrigDefs::includeFailedDecisions, "HLT_MuonL2SAInfoIOmode" );
      for(const TrigCompositeUtils::LinkInfo<xAOD::L2StandAloneMuonContainer>& l2saLinkInfo : l2saLinksIO){
        if( !l2saLinkInfo.isValid() ) continue;
        const ElementLink<xAOD::L2StandAloneMuonContainer> l2sa = l2saLinkInfo.link;
        if( !l2sa.isValid() ) continue;
        if( l2saLinkInfo.state == TrigCompositeUtils::ActiveState::ACTIVE ) tdtobj.isPassedSAIO_evt = true;

        SAObject saioobj;
        bool ispassed_saio = false;
        if(l2saLinkInfo.state == TrigCompositeUtils::ActiveState::ACTIVE) ispassed_saio = true;
        saioobj.isPassed = ispassed_saio;
        saioobj.roiNum = (*l2sa)->roiNumber();
        saioobj.roiSector = (*l2sa)->roiSector();
        saioobj.pt = (*l2sa)->pt();
        saioobj.eta = (*l2sa)->eta();
        saioobj.phi = (*l2sa)->phi();
        saioobj.etaMS = (*l2sa)->etaMS();
        saioobj.phiMS = (*l2sa)->phiMS();
        saioobj.sAddress = (*l2sa)->sAddress();
        saioobj.roiEta  = (*l2sa)->roiEta();
        saioobj.roiPhi  = (*l2sa)->roiPhi();
        saioobj.superPointR_BI  = (*l2sa)->superPointR(0);
        saioobj.superPointR_BM  = (*l2sa)->superPointR(1);
        saioobj.superPointR_BO  = (*l2sa)->superPointR(2);
        saioobj.superPointR_EI  = (*l2sa)->superPointR(3);
        saioobj.superPointR_EM  = (*l2sa)->superPointR(4);
        saioobj.superPointR_EO  = (*l2sa)->superPointR(5);
        saioobj.superPointR_EE  = (*l2sa)->superPointR(6);
        saioobj.superPointR_CSC = (*l2sa)->superPointR(7);
        saioobj.superPointR_BEE = (*l2sa)->superPointR(8);
        saioobj.superPointR_BME = (*l2sa)->superPointR(9);
        saioobj.superPointZ_BI  = (*l2sa)->superPointZ(0);
        saioobj.superPointZ_BM  = (*l2sa)->superPointZ(1);
        saioobj.superPointZ_BO  = (*l2sa)->superPointZ(2);
        saioobj.superPointZ_EI  = (*l2sa)->superPointZ(3);
        saioobj.superPointZ_EM  = (*l2sa)->superPointZ(4);
        saioobj.superPointZ_EO  = (*l2sa)->superPointZ(5);
        saioobj.superPointZ_EE  = (*l2sa)->superPointZ(6);
        saioobj.superPointZ_CSC = (*l2sa)->superPointZ(7);
        saioobj.superPointZ_BEE = (*l2sa)->superPointZ(8);
        saioobj.superPointZ_BME = (*l2sa)->superPointZ(9);
        std::cout << "SA IOmode w/ TDT: pt/eta/phi/isPassed = " << saioobj.pt << "/" << saioobj.eta << "/"
            << saioobj.phi << "/" << ispassed_saio << std::endl;
        tdtobj.m_saio.push_back(saioobj);
      }
    }
    //CB
    std::vector< TrigCompositeUtils::LinkInfo<xAOD::L2CombinedMuonContainer> > l2cbLinks = m_trigDecTool->features<xAOD::L2CombinedMuonContainer>( m_hltchainName.at(i_trig), TrigDefs::includeFailedDecisions, "HLT_MuonL2CBInfo");
    //
    for(const TrigCompositeUtils::LinkInfo<xAOD::L2CombinedMuonContainer>& l2cbLinkInfo : l2cbLinks){
      CBObject cbobj;
      bool isPassedCBobj = false;
      if( !l2cbLinkInfo.isValid() ) continue;
      const ElementLink<xAOD::L2CombinedMuonContainer> l2cb = l2cbLinkInfo.link;
      if( !l2cb.isValid() ) continue;
      if( l2cbLinkInfo.state == TrigCompositeUtils::ActiveState::ACTIVE ){ 
        isPassedCBobj = true;
        tdtobj.isPassedCB_evt = true;
      }
      cbobj.isPassed = isPassedCBobj;
      if((*l2cb)->errorFlag() == 1){
        cbobj.roiNum = -99999;
        cbobj.roiSector = -99999;
      } else {
        cbobj.roiNum = (*l2cb)->muSATrack()->roiNumber();
        cbobj.roiSector = (*l2cb)->muSATrack()->roiSector();
      }
      if((*l2cb)->errorFlag() == 0){
        cbobj.idpt = (*l2cb)->idTrack()->pt();
        cbobj.ideta = (*l2cb)->idTrack()->eta();
        cbobj.idphi = (*l2cb)->idTrack()->phi();
      }
      cbobj.pt = (*l2cb)->pt();
      cbobj.eta = (*l2cb)->eta();
      cbobj.phi = (*l2cb)->phi();
      std::cout << "CB default w/ TDT: pt/eta/phi/isPassed = " << (*l2cb)->pt() << "/" << (*l2cb)->eta() << "/"
          << (*l2cb)->phi() << "/" << isPassedCBobj << std::endl;
      cbobj.isPassed = isPassedCBobj;
      tdtobj.m_cb.push_back(cbobj);
    }
    //CBIO
    //
    if(m_isvalidIOmode){
      std::vector< TrigCompositeUtils::LinkInfo<xAOD::L2CombinedMuonContainer> > l2cbLinksIO = m_trigDecTool->features<xAOD::L2CombinedMuonContainer>( m_hltchainName.at(i_trig), TrigDefs::includeFailedDecisions, "HLT_MuonL2CBInfoIOmode");
      //
      for(const TrigCompositeUtils::LinkInfo<xAOD::L2CombinedMuonContainer>& l2cbLinkInfo : l2cbLinksIO){
        CBObject cbobj;
        bool isPassedCBobj = false;
        if( !l2cbLinkInfo.isValid() ) continue;
        const ElementLink<xAOD::L2CombinedMuonContainer> l2cb = l2cbLinkInfo.link;
        if( !l2cb.isValid() ) continue;
        if( l2cbLinkInfo.state == TrigCompositeUtils::ActiveState::ACTIVE ){ 
          isPassedCBobj = true;
          tdtobj.isPassedCBIO_evt = true;
        }
        cbobj.isPassed = isPassedCBobj;
        if((*l2cb)->errorFlag() == 1){
          cbobj.roiNum = -99999;
          cbobj.roiSector = -99999;
        } else {
          cbobj.roiNum = (*l2cb)->muSATrack()->roiNumber();
          cbobj.roiSector = (*l2cb)->muSATrack()->roiSector();
        }
        if((*l2cb)->errorFlag() == 0){
          cbobj.idpt = (*l2cb)->idTrack()->pt();
          cbobj.ideta = (*l2cb)->idTrack()->eta();
          cbobj.idphi = (*l2cb)->idTrack()->phi();
        }
        cbobj.pt = (*l2cb)->pt();
        cbobj.eta = (*l2cb)->eta();
        cbobj.phi = (*l2cb)->phi();
        std::cout << "CB IOmode w/ TDT: pt/eta/phi/isPassed = " << (*l2cb)->pt() << "/" << (*l2cb)->eta() << "/"
          << (*l2cb)->phi() << "/" << isPassedCBobj << std::endl;
        cbobj.isPassed = isPassedCBobj;
        tdtobj.m_cbio.push_back(cbobj);
      }
    }
    std::vector< TrigCompositeUtils::LinkInfo<xAOD::MuonContainer> > efLinks = m_trigDecTool->features<xAOD::MuonContainer>( m_hltchainName.at(i_trig), TrigDefs::includeFailedDecisions );
    //
    for(const TrigCompositeUtils::LinkInfo<xAOD::MuonContainer>& efLinkInfo : efLinks){
      if( !efLinkInfo.isValid() ) continue;
      const ElementLink<xAOD::MuonContainer> ef = efLinkInfo.link;
      if( !ef.isValid() ) continue;
      EFObject efobj;
      if( efLinkInfo.state == TrigCompositeUtils::ActiveState::ACTIVE ){
        tdtobj.isPassedEF_evt = true;
        efobj.isPassed = true;
      }
      efobj.pt = (*ef)->pt();
      efobj.eta = (*ef)->eta();
      efobj.phi = (*ef)->phi();
      tdtobj.m_ef.push_back(efobj);
      std::cout << "EF w/ TDT: pt/eta/phi/isPassed = " << (*ef)->pt() << "/" << (*ef)->eta() << "/" << (*ef)->phi() << "/" << efobj.isPassed << std::endl;
    }
    
    std::cout << "-----> " << m_hltchainName.at(i_trig) << " ispassed L1/SA/SAIO/CB/CBIO/EF = " << 
    tdtobj.isPassedL1_evt << "/" << tdtobj.isPassedSA_evt << "/" << tdtobj.isPassedSAIO_evt << "/" << 
    tdtobj.isPassedCB_evt << "/" << tdtobj.isPassedCBIO_evt << "/" << tdtobj.isPassedEF_evt << std::endl;
    m_TDTobjects.push_back(tdtobj);
  }

}
//---------------------------------------------------------//
//---------------------------------------------------------//
