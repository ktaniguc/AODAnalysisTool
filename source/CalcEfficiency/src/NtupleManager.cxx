#include <iostream>
#include <fstream>
#include <vector>

#include "CalcEfficiency/NtupleManager.h"
#include "CalcEfficiency/Utils.h"

#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "CalcEfficiency/SAObjects.h"
#include "CalcEfficiency/CBObjects.h"

NtupleManager::NtupleManager() {
}

NtupleManager::~NtupleManager() {
}

int NtupleManager::initialize( std::string outfile = "test.root" ) {
  std::string name4ntuple = outfile;
  int i_name = name4ntuple.find(".");
  name4ntuple.erase(i_name);
  std::string ntupname = name4ntuple+"_Ntuple.root";

  //initialize the event tree
  m_file	= new TFile( ntupname.c_str(), "recreate" );
  //m_file	= outfile;
  m_tree 	= new TTree( "ntuple", "AODVariables" );
  //m_tree    = outtree;
  //m_tree->SetDirectory( m_file );

  //--------------------------------------------------
  // VARIABLE SET UP
  //--------------------------------------------------
	
  // Event info
  eventNumber = -1;
  runNumber   = -1;
  // trigger variables
  n_trig            = 0;
  trigname          = new std::vector<std::string>();   
  isPassedTrig      = new std::vector<bool> ();
  isPassedL1_evt      = new std::vector<bool> ();
  isPassedSA_evt      = new std::vector<bool> ();
  isPassedCB_evt      = new std::vector<bool> ();
  isPassedSAIO_evt      = new std::vector<bool> ();
  isPassedCBIO_evt      = new std::vector<bool> ();
  isPassedEF_evt      = new std::vector<bool> ();
  isPassedL1      = new std::vector<std::vector<bool>> ();
  isPassedSA      = new std::vector<std::vector<bool>> ();
  isPassedCB      = new std::vector<std::vector<bool>> ();
  isPassedSAIO      = new std::vector<std::vector<bool>> ();
  isPassedCBIO      = new std::vector<std::vector<bool>> ();
  isPassedEF      = new std::vector<std::vector<bool>> ();
  tdt_L1RoINumber      = new std::vector<std::vector<int>> ();
  tdt_SARoINumber      = new std::vector<std::vector<int>> ();
  tdt_CBRoINumber      = new std::vector<std::vector<int>> ();
  tdt_SAIORoINumber      = new std::vector<std::vector<int>> ();
  tdt_CBIORoINumber      = new std::vector<std::vector<int>> ();
  tdt_L1RoISector      = new std::vector<std::vector<int>> ();
  tdt_SARoISector      = new std::vector<std::vector<int>> ();
  tdt_CBRoISector      = new std::vector<std::vector<int>> ();
  tdt_SAIORoISector      = new std::vector<std::vector<int>> ();
  tdt_CBIORoISector      = new std::vector<std::vector<int>> ();
  tdt_L1isMoreCandInRoI      = new std::vector<std::vector<bool>> ();
  tdt_SApt    = new std::vector < std::vector<double> > ();
  tdt_SAeta   = new std::vector < std::vector<double> > ();
  tdt_SAphi   = new std::vector < std::vector<double> > ();
  tdt_SAetaMS   = new std::vector < std::vector<double> > ();
  tdt_SAphiMS   = new std::vector < std::vector<double> > ();
  tdt_SAsAddress = new std::vector < std::vector<int> > ();
  tdt_SAroiEta = new std::vector < std::vector<float> > ();
  tdt_SAroiPhi = new std::vector < std::vector<float> > ();
  tdt_SAsuperPointR_BI  = new std::vector < std::vector<double> > ();
  tdt_SAsuperPointR_BM  = new std::vector < std::vector<double> > ();
  tdt_SAsuperPointR_BO  = new std::vector < std::vector<double> > ();
  tdt_SAsuperPointR_EI  = new std::vector < std::vector<double> > ();
  tdt_SAsuperPointR_EM  = new std::vector < std::vector<double> > ();
  tdt_SAsuperPointR_EO  = new std::vector < std::vector<double> > ();
  tdt_SAsuperPointR_EE  = new std::vector < std::vector<double> > ();
  tdt_SAsuperPointR_CSC = new std::vector < std::vector<double> > ();
  tdt_SAsuperPointR_BEE = new std::vector < std::vector<double> > ();
  tdt_SAsuperPointR_BME = new std::vector < std::vector<double> > ();
  tdt_SAsuperPointZ_BI  = new std::vector < std::vector<double> > ();
  tdt_SAsuperPointZ_BM  = new std::vector < std::vector<double> > ();
  tdt_SAsuperPointZ_BO  = new std::vector < std::vector<double> > ();
  tdt_SAsuperPointZ_EI  = new std::vector < std::vector<double> > ();
  tdt_SAsuperPointZ_EM  = new std::vector < std::vector<double> > ();
  tdt_SAsuperPointZ_EO  = new std::vector < std::vector<double> > ();
  tdt_SAsuperPointZ_EE  = new std::vector < std::vector<double> > ();
  tdt_SAsuperPointZ_CSC = new std::vector < std::vector<double> > ();
  tdt_SAsuperPointZ_BEE = new std::vector < std::vector<double> > ();
  tdt_SAsuperPointZ_BME = new std::vector < std::vector<double> > ();
  tdt_SAIOpt    = new std::vector < std::vector<double> > ();
  tdt_SAIOeta   = new std::vector < std::vector<double> > ();
  tdt_SAIOphi   = new std::vector < std::vector<double> > ();
  tdt_SAIOetaMS   = new std::vector < std::vector<double> > ();
  tdt_SAIOphiMS   = new std::vector < std::vector<double> > ();
  tdt_SAIOsAddress = new std::vector < std::vector<int> > ();
  tdt_SAIOroiEta = new std::vector < std::vector<float> > ();
  tdt_SAIOroiPhi = new std::vector < std::vector<float> > ();
  tdt_SAIOsuperPointR_BI  = new std::vector < std::vector<double> > ();
  tdt_SAIOsuperPointR_BM  = new std::vector < std::vector<double> > ();
  tdt_SAIOsuperPointR_BO  = new std::vector < std::vector<double> > ();
  tdt_SAIOsuperPointR_EI  = new std::vector < std::vector<double> > ();
  tdt_SAIOsuperPointR_EM  = new std::vector < std::vector<double> > ();
  tdt_SAIOsuperPointR_EO  = new std::vector < std::vector<double> > ();
  tdt_SAIOsuperPointR_EE  = new std::vector < std::vector<double> > ();
  tdt_SAIOsuperPointR_CSC = new std::vector < std::vector<double> > ();
  tdt_SAIOsuperPointR_BEE = new std::vector < std::vector<double> > ();
  tdt_SAIOsuperPointR_BME = new std::vector < std::vector<double> > ();
  tdt_SAIOsuperPointZ_BI  = new std::vector < std::vector<double> > ();
  tdt_SAIOsuperPointZ_BM  = new std::vector < std::vector<double> > ();
  tdt_SAIOsuperPointZ_BO  = new std::vector < std::vector<double> > ();
  tdt_SAIOsuperPointZ_EI  = new std::vector < std::vector<double> > ();
  tdt_SAIOsuperPointZ_EM  = new std::vector < std::vector<double> > ();
  tdt_SAIOsuperPointZ_EO  = new std::vector < std::vector<double> > ();
  tdt_SAIOsuperPointZ_EE  = new std::vector < std::vector<double> > ();
  tdt_SAIOsuperPointZ_CSC = new std::vector < std::vector<double> > ();
  tdt_SAIOsuperPointZ_BEE = new std::vector < std::vector<double> > ();
  tdt_SAIOsuperPointZ_BME = new std::vector < std::vector<double> > ();
  tdt_CBpt    = new std::vector < std::vector<double> > ();
  tdt_CBeta   = new std::vector < std::vector<double> > ();
  tdt_CBphi   = new std::vector < std::vector<double> > ();
  tdt_CBIOpt    = new std::vector < std::vector<double> > ();
  tdt_CBIOeta   = new std::vector < std::vector<double> > ();
  tdt_CBIOphi   = new std::vector < std::vector<double> > ();
  tdt_CBidpt    = new std::vector < std::vector<double> > ();
  tdt_CBideta   = new std::vector < std::vector<double> > ();
  tdt_CBidphi   = new std::vector < std::vector<double> > ();
  tdt_CBIOidpt    = new std::vector < std::vector<double> > ();
  tdt_CBIOideta   = new std::vector < std::vector<double> > ();
  tdt_CBIOidphi   = new std::vector < std::vector<double> > ();
  //
  muon_pt           = new std::vector < double > ();
  muon_eta          = new std::vector < double > ();
  muon_phi          = new std::vector < double > ();
  muon_extEta       = new std::vector < double > ();
  muon_extPhi       = new std::vector < double > ();
  //
  L1_eta   = new std::vector < double > ();
  L1_phi   = new std::vector < double > ();
  L1_roiNum   = new std::vector < int > ();
  L1_roiSector   = new std::vector < int > ();
  L1_thrValue   = new std::vector < double > ();
  L1_thrNumber   = new std::vector < int > ();
  L1_isMoreCandInRoI   = new std::vector < bool > ();
  //
  SA_pt    = new std::vector < double > ();
  SA_eta   = new std::vector < double > ();
  SA_phi   = new std::vector < double > ();
  SA_etaMS   = new std::vector < double > ();
  SA_phiMS   = new std::vector < double > ();
  SA_sAddress = new std::vector < int > ();
  SA_roiEta = new std::vector < float > ();
  SA_roiPhi = new std::vector < float > ();
  SA_roiNum = new std::vector < int > ();
  SA_roiSector = new std::vector < int > ();

  //the measured radious of the muon in one particular super point
  SA_superPointR_BI  = new std::vector < double > ();
  SA_superPointR_BM  = new std::vector < double > ();
  SA_superPointR_BO  = new std::vector < double > ();
  SA_superPointR_EI  = new std::vector < double > ();
  SA_superPointR_EM  = new std::vector < double > ();
  SA_superPointR_EO  = new std::vector < double > ();
  SA_superPointR_EE  = new std::vector < double > ();
  SA_superPointR_CSC = new std::vector < double > ();
  SA_superPointR_BEE = new std::vector < double > ();
  SA_superPointR_BME = new std::vector < double > ();
  //the measured Z position of the muon in one particular super point
  SA_superPointZ_BI  = new std::vector < double > ();
  SA_superPointZ_BM  = new std::vector < double > ();
  SA_superPointZ_BO  = new std::vector < double > ();
  SA_superPointZ_EI  = new std::vector < double > ();
  SA_superPointZ_EM  = new std::vector < double > ();
  SA_superPointZ_EO  = new std::vector < double > ();
  SA_superPointZ_EE  = new std::vector < double > ();
  SA_superPointZ_CSC = new std::vector < double > ();
  SA_superPointZ_BEE = new std::vector < double > ();
  SA_superPointZ_BME = new std::vector < double > ();
  SAIO_pt    = new std::vector < double > ();
  SAIO_eta   = new std::vector < double > ();
  SAIO_phi   = new std::vector < double > ();
  SAIO_etaMS   = new std::vector < double > ();
  SAIO_phiMS   = new std::vector < double > ();
  SAIO_sAddress = new std::vector < int > ();
  SAIO_roiEta = new std::vector < float > ();
  SAIO_roiPhi = new std::vector < float > ();
  SAIO_roiNum = new std::vector < int > ();
  SAIO_roiSector = new std::vector < int > ();

  //the measured radious of the muon in one particular super point
  SAIO_superPointR_BI  = new std::vector < double > ();
  SAIO_superPointR_BM  = new std::vector < double > ();
  SAIO_superPointR_BO  = new std::vector < double > ();
  SAIO_superPointR_EI  = new std::vector < double > ();
  SAIO_superPointR_EM  = new std::vector < double > ();
  SAIO_superPointR_EO  = new std::vector < double > ();
  SAIO_superPointR_EE  = new std::vector < double > ();
  SAIO_superPointR_CSC = new std::vector < double > ();
  SAIO_superPointR_BEE = new std::vector < double > ();
  SAIO_superPointR_BME = new std::vector < double > ();
  //the measured Z position of the muon in one particular super point
  SAIO_superPointZ_BI  = new std::vector < double > ();
  SAIO_superPointZ_BM  = new std::vector < double > ();
  SAIO_superPointZ_BO  = new std::vector < double > ();
  SAIO_superPointZ_EI  = new std::vector < double > ();
  SAIO_superPointZ_EM  = new std::vector < double > ();
  SAIO_superPointZ_EO  = new std::vector < double > ();
  SAIO_superPointZ_EE  = new std::vector < double > ();
  SAIO_superPointZ_CSC = new std::vector < double > ();
  SAIO_superPointZ_BEE = new std::vector < double > ();
  SAIO_superPointZ_BME = new std::vector < double > ();
  //
  CB_pt    = new std::vector < double > ();
  CB_eta   = new std::vector < double > ();
  CB_phi   = new std::vector < double > ();
  CB_idpt    = new std::vector < double > ();
  CB_ideta   = new std::vector < double > ();
  CB_idphi   = new std::vector < double > ();
  CB_roiNumber   = new std::vector < int > ();
  CB_roiSector   = new std::vector < int > ();
  //
  CBIO_pt    = new std::vector < double > ();
  CBIO_eta   = new std::vector < double > ();
  CBIO_phi   = new std::vector < double > ();
  CBIO_idpt    = new std::vector < double > ();
  CBIO_ideta   = new std::vector < double > ();
  CBIO_idphi   = new std::vector < double > ();
  CBIO_roiNumber   = new std::vector < int > ();
  CBIO_roiSector   = new std::vector < int > ();
  //
  EF_pt    = new std::vector < double > ();
  EF_eta   = new std::vector < double > ();
  EF_phi   = new std::vector < double > ();
    
  //--------------------------------------------------
  // BRANCH SET UP
  //--------------------------------------------------
 
  //Event info 
  m_tree->Branch( "eventNumber",       &eventNumber );
  m_tree->Branch( "runNumber",       &runNumber );
  m_tree->Branch( "trigname",       &trigname );
  m_tree->Branch( "isPassedTrig",       &isPassedTrig );
  m_tree->Branch( "isPassedL1_evt",       &isPassedL1_evt );
  m_tree->Branch( "isPassedSA_evt",       &isPassedSA_evt );
  m_tree->Branch( "isPassedCB_evt",       &isPassedCB_evt );
  m_tree->Branch( "isPassedSAIO_evt",       &isPassedSAIO_evt );
  m_tree->Branch( "isPassedCBIO_evt",       &isPassedCBIO_evt );
  m_tree->Branch( "isPassedEF_evt",       &isPassedEF_evt );
  m_tree->Branch( "isPassedL1",       &isPassedL1 );
  m_tree->Branch( "isPassedSA",       &isPassedSA );
  m_tree->Branch( "isPassedCB",       &isPassedCB );
  m_tree->Branch( "isPassedSAIO",       &isPassedSAIO );
  m_tree->Branch( "isPassedCBIO",       &isPassedCBIO );
  m_tree->Branch( "isPassedEF",       &isPassedEF );
  m_tree->Branch( "tdt_L1RoINumber",       &tdt_L1RoINumber );
  m_tree->Branch( "tdt_SARoINumber",       &tdt_SARoINumber );
  m_tree->Branch( "tdt_CBRoINumber",       &tdt_CBRoINumber );
  m_tree->Branch( "tdt_SAIORoINumber",       &tdt_SAIORoINumber );
  m_tree->Branch( "tdt_CBIORoINumber",       &tdt_CBIORoINumber );
  m_tree->Branch( "tdt_L1RoISector",       &tdt_L1RoISector );
  m_tree->Branch( "tdt_SARoISector",       &tdt_SARoISector );
  m_tree->Branch( "tdt_CBRoISector",       &tdt_CBRoISector );
  m_tree->Branch( "tdt_SAIORoISector",       &tdt_SAIORoISector );
  m_tree->Branch( "tdt_CBIORoISector",       &tdt_CBIORoISector );
  m_tree->Branch( "tdt_L1isMoreCandInRoI",       &tdt_L1isMoreCandInRoI );
  m_tree->Branch( "tdt_SApt",     &tdt_SApt );
  m_tree->Branch( "tdt_SAeta",    &tdt_SAeta );
  m_tree->Branch( "tdt_SAphi",    &tdt_SAphi );
  m_tree->Branch( "tdt_SAetaMS",    &tdt_SAetaMS );
  m_tree->Branch( "tdt_SAphiMS",    &tdt_SAphiMS );
  //
  m_tree->Branch( "tdt_SAsAddress",  &tdt_SAsAddress );
  m_tree->Branch( "tdt_SAroiEta",  &tdt_SAroiEta );
  m_tree->Branch( "tdt_SAroiPhi",  &tdt_SAroiPhi );
  //
  m_tree->Branch( "tdt_SAsuperPointR_BI",   &tdt_SAsuperPointR_BI );
  m_tree->Branch( "tdt_SAsuperPointR_BM",   &tdt_SAsuperPointR_BM );
  m_tree->Branch( "tdt_SAsuperPointR_BO",   &tdt_SAsuperPointR_BO );
  m_tree->Branch( "tdt_SAsuperPointR_EI",   &tdt_SAsuperPointR_EI );
  m_tree->Branch( "tdt_SAsuperPointR_EM",   &tdt_SAsuperPointR_EM );
  m_tree->Branch( "tdt_SAsuperPointR_EO",   &tdt_SAsuperPointR_EO );
  m_tree->Branch( "tdt_SAsuperPointR_EE",   &tdt_SAsuperPointR_EE );
  m_tree->Branch( "tdt_SAsuperPointR_CSC",  &tdt_SAsuperPointR_CSC );
  m_tree->Branch( "tdt_SAsuperPointR_BEE",  &tdt_SAsuperPointR_BEE );
  m_tree->Branch( "tdt_SAsuperPointR_BME",  &tdt_SAsuperPointR_BME );
  //
  m_tree->Branch( "tdt_SAsuperPointZ_BI",   &tdt_SAsuperPointZ_BI );
  m_tree->Branch( "tdt_SAsuperPointZ_BM",   &tdt_SAsuperPointZ_BM );
  m_tree->Branch( "tdt_SAsuperPointZ_BO",   &tdt_SAsuperPointZ_BO );
  m_tree->Branch( "tdt_SAsuperPointZ_EI",   &tdt_SAsuperPointZ_EI );
  m_tree->Branch( "tdt_SAsuperPointZ_EM",   &tdt_SAsuperPointZ_EM );
  m_tree->Branch( "tdt_SAsuperPointZ_EO",   &tdt_SAsuperPointZ_EO );
  m_tree->Branch( "tdt_SAsuperPointZ_EE",   &tdt_SAsuperPointZ_EE );
  m_tree->Branch( "tdt_SAsuperPointZ_CSC",  &tdt_SAsuperPointZ_CSC );
  m_tree->Branch( "tdt_SAsuperPointZ_BEE",  &tdt_SAsuperPointZ_BEE );
  m_tree->Branch( "tdt_SAsuperPointZ_BME",  &tdt_SAsuperPointZ_BME );
  m_tree->Branch( "tdt_SAIOpt",     &tdt_SAIOpt );
  m_tree->Branch( "tdt_SAIOeta",    &tdt_SAIOeta );
  m_tree->Branch( "tdt_SAIOphi",    &tdt_SAIOphi );
  m_tree->Branch( "tdt_SAIOetaMS",    &tdt_SAIOetaMS );
  m_tree->Branch( "tdt_SAIOphiMS",    &tdt_SAIOphiMS );
  //
  m_tree->Branch( "tdt_SAIOsAddress",  &tdt_SAIOsAddress );
  m_tree->Branch( "tdt_SAIOroiEta",  &tdt_SAIOroiEta );
  m_tree->Branch( "tdt_SAIOroiPhi",  &tdt_SAIOroiPhi );
  //
  m_tree->Branch( "tdt_SAIOsuperPointR_BI",   &tdt_SAIOsuperPointR_BI );
  m_tree->Branch( "tdt_SAIOsuperPointR_BM",   &tdt_SAIOsuperPointR_BM );
  m_tree->Branch( "tdt_SAIOsuperPointR_BO",   &tdt_SAIOsuperPointR_BO );
  m_tree->Branch( "tdt_SAIOsuperPointR_EI",   &tdt_SAIOsuperPointR_EI );
  m_tree->Branch( "tdt_SAIOsuperPointR_EM",   &tdt_SAIOsuperPointR_EM );
  m_tree->Branch( "tdt_SAIOsuperPointR_EO",   &tdt_SAIOsuperPointR_EO );
  m_tree->Branch( "tdt_SAIOsuperPointR_EE",   &tdt_SAIOsuperPointR_EE );
  m_tree->Branch( "tdt_SAIOsuperPointR_CSC",  &tdt_SAIOsuperPointR_CSC );
  m_tree->Branch( "tdt_SAIOsuperPointR_BEE",  &tdt_SAIOsuperPointR_BEE );
  m_tree->Branch( "tdt_SAIOsuperPointR_BME",  &tdt_SAIOsuperPointR_BME );
  //
  m_tree->Branch( "tdt_SAIOsuperPointZ_BI",   &tdt_SAIOsuperPointZ_BI );
  m_tree->Branch( "tdt_SAIOsuperPointZ_BM",   &tdt_SAIOsuperPointZ_BM );
  m_tree->Branch( "tdt_SAIOsuperPointZ_BO",   &tdt_SAIOsuperPointZ_BO );
  m_tree->Branch( "tdt_SAIOsuperPointZ_EI",   &tdt_SAIOsuperPointZ_EI );
  m_tree->Branch( "tdt_SAIOsuperPointZ_EM",   &tdt_SAIOsuperPointZ_EM );
  m_tree->Branch( "tdt_SAIOsuperPointZ_EO",   &tdt_SAIOsuperPointZ_EO );
  m_tree->Branch( "tdt_SAIOsuperPointZ_EE",   &tdt_SAIOsuperPointZ_EE );
  m_tree->Branch( "tdt_SAIOsuperPointZ_CSC",  &tdt_SAIOsuperPointZ_CSC );
  m_tree->Branch( "tdt_SAIOsuperPointZ_BEE",  &tdt_SAIOsuperPointZ_BEE );
  m_tree->Branch( "tdt_SAIOsuperPointZ_BME",  &tdt_SAIOsuperPointZ_BME );
  m_tree->Branch( "tdt_CBpt",     &tdt_CBpt );
  m_tree->Branch( "tdt_CBeta",    &tdt_CBeta );
  m_tree->Branch( "tdt_CBphi",    &tdt_CBphi );
  m_tree->Branch( "tdt_CBIOpt",     &tdt_CBIOpt );
  m_tree->Branch( "tdt_CBIOeta",    &tdt_CBIOeta );
  m_tree->Branch( "tdt_CBIOphi",    &tdt_CBIOphi );
  m_tree->Branch( "tdt_CBidpt",     &tdt_CBidpt );
  m_tree->Branch( "tdt_CBideta",    &tdt_CBideta );
  m_tree->Branch( "tdt_CBidphi",    &tdt_CBidphi );
  m_tree->Branch( "tdt_CBIOidpt",     &tdt_CBIOidpt );
  m_tree->Branch( "tdt_CBIOideta",    &tdt_CBIOideta );
  m_tree->Branch( "tdt_CBIOidphi",    &tdt_CBIOidphi );
  m_tree->Branch( "n_trig",       &n_trig );
  //offline muon
  m_tree->Branch( "muon_pt",     &muon_pt );
  m_tree->Branch( "muon_eta",    &muon_eta );
  m_tree->Branch( "muon_phi",    &muon_phi );
  m_tree->Branch( "muon_extEta",    &muon_extEta );
  m_tree->Branch( "muon_extPhi",    &muon_extPhi );
  //L1
  m_tree->Branch( "L1_eta",    &L1_eta );
  m_tree->Branch( "L1_phi",    &L1_phi );
  m_tree->Branch( "L1_thrValue",    &L1_thrValue );
  m_tree->Branch( "L1_roiNum",    &L1_roiNum );
  m_tree->Branch( "L1_roiSector",    &L1_roiSector );
  m_tree->Branch( "L1_thrNumber",    &L1_thrNumber );
  m_tree->Branch( "L1_isMoreCandInRoI",    &L1_isMoreCandInRoI );
  //SA
  m_tree->Branch( "SAIO_roiNumber",  &SAIO_roiNum );
  m_tree->Branch( "SAIO_roiSector",  &SAIO_roiSector );
  m_tree->Branch( "SAIO_pt",     &SAIO_pt );
  m_tree->Branch( "SAIO_eta",    &SAIO_eta );
  m_tree->Branch( "SAIO_phi",    &SAIO_phi );
  m_tree->Branch( "SAIO_etaMS",    &SAIO_etaMS );
  m_tree->Branch( "SAIO_phiMS",    &SAIO_phiMS );
  //
  m_tree->Branch( "SAIO_sAddress",  &SAIO_sAddress );
  m_tree->Branch( "SAIO_roiEta",  &SAIO_roiEta );
  m_tree->Branch( "SAIO_roiPhi",  &SAIO_roiPhi );
  //
  m_tree->Branch( "SAIO_superPointR_BI",   &SAIO_superPointR_BI );
  m_tree->Branch( "SAIO_superPointR_BM",   &SAIO_superPointR_BM );
  m_tree->Branch( "SAIO_superPointR_BO",   &SAIO_superPointR_BO );
  m_tree->Branch( "SAIO_superPointR_EI",   &SAIO_superPointR_EI );
  m_tree->Branch( "SAIO_superPointR_EM",   &SAIO_superPointR_EM );
  m_tree->Branch( "SAIO_superPointR_EO",   &SAIO_superPointR_EO );
  m_tree->Branch( "SAIO_superPointR_EE",   &SAIO_superPointR_EE );
  m_tree->Branch( "SAIO_superPointR_CSC",  &SAIO_superPointR_CSC );
  m_tree->Branch( "SAIO_superPointR_BEE",  &SAIO_superPointR_BEE );
  m_tree->Branch( "SAIO_superPointR_BME",  &SAIO_superPointR_BME );
  //
  m_tree->Branch( "SAIO_superPointZ_BI",   &SAIO_superPointZ_BI );
  m_tree->Branch( "SAIO_superPointZ_BM",   &SAIO_superPointZ_BM );
  m_tree->Branch( "SAIO_superPointZ_BO",   &SAIO_superPointZ_BO );
  m_tree->Branch( "SAIO_superPointZ_EI",   &SAIO_superPointZ_EI );
  m_tree->Branch( "SAIO_superPointZ_EM",   &SAIO_superPointZ_EM );
  m_tree->Branch( "SAIO_superPointZ_EO",   &SAIO_superPointZ_EO );
  m_tree->Branch( "SAIO_superPointZ_EE",   &SAIO_superPointZ_EE );
  m_tree->Branch( "SAIO_superPointZ_CSC",  &SAIO_superPointZ_CSC );
  m_tree->Branch( "SAIO_superPointZ_BEE",  &SAIO_superPointZ_BEE );
  m_tree->Branch( "SAIO_superPointZ_BME",  &SAIO_superPointZ_BME );
  // 
  m_tree->Branch( "SA_roiNumber",  &SA_roiNum );
  m_tree->Branch( "SA_roiSector",  &SA_roiSector );
  m_tree->Branch( "SA_pt",     &SA_pt );
  m_tree->Branch( "SA_eta",    &SA_eta );
  m_tree->Branch( "SA_phi",    &SA_phi );
  m_tree->Branch( "SA_etaMS",    &SA_etaMS );
  m_tree->Branch( "SA_phiMS",    &SA_phiMS );
  //
  m_tree->Branch( "SA_sAddress",  &SA_sAddress );
  m_tree->Branch( "SA_roiEta",  &SA_roiEta );
  m_tree->Branch( "SA_roiPhi",  &SA_roiPhi );
  //
  m_tree->Branch( "SA_superPointR_BI",   &SA_superPointR_BI );
  m_tree->Branch( "SA_superPointR_BM",   &SA_superPointR_BM );
  m_tree->Branch( "SA_superPointR_BO",   &SA_superPointR_BO );
  m_tree->Branch( "SA_superPointR_EI",   &SA_superPointR_EI );
  m_tree->Branch( "SA_superPointR_EM",   &SA_superPointR_EM );
  m_tree->Branch( "SA_superPointR_EO",   &SA_superPointR_EO );
  m_tree->Branch( "SA_superPointR_EE",   &SA_superPointR_EE );
  m_tree->Branch( "SA_superPointR_CSC",  &SA_superPointR_CSC );
  m_tree->Branch( "SA_superPointR_BEE",  &SA_superPointR_BEE );
  m_tree->Branch( "SA_superPointR_BME",  &SA_superPointR_BME );
  //
  m_tree->Branch( "SA_superPointZ_BI",   &SA_superPointZ_BI );
  m_tree->Branch( "SA_superPointZ_BM",   &SA_superPointZ_BM );
  m_tree->Branch( "SA_superPointZ_BO",   &SA_superPointZ_BO );
  m_tree->Branch( "SA_superPointZ_EI",   &SA_superPointZ_EI );
  m_tree->Branch( "SA_superPointZ_EM",   &SA_superPointZ_EM );
  m_tree->Branch( "SA_superPointZ_EO",   &SA_superPointZ_EO );
  m_tree->Branch( "SA_superPointZ_EE",   &SA_superPointZ_EE );
  m_tree->Branch( "SA_superPointZ_CSC",  &SA_superPointZ_CSC );
  m_tree->Branch( "SA_superPointZ_BEE",  &SA_superPointZ_BEE );
  m_tree->Branch( "SA_superPointZ_BME",  &SA_superPointZ_BME );
  // 
  //CB
  m_tree->Branch( "CB_pt",     &CB_pt );
  m_tree->Branch( "CB_eta",    &CB_eta );
  m_tree->Branch( "CB_phi",    &CB_phi );
  m_tree->Branch( "CB_idpt",     &CB_idpt );
  m_tree->Branch( "CB_ideta",    &CB_ideta );
  m_tree->Branch( "CB_idphi",    &CB_idphi );
  m_tree->Branch( "CB_roiNumber",    &CB_roiNumber );
  m_tree->Branch( "CB_roiSector",    &CB_roiSector );
  //CB
  m_tree->Branch( "CBIO_pt",     &CBIO_pt );
  m_tree->Branch( "CBIO_eta",    &CBIO_eta );
  m_tree->Branch( "CBIO_phi",    &CBIO_phi );
  m_tree->Branch( "CBIO_idpt",     &CBIO_idpt );
  m_tree->Branch( "CBIO_ideta",    &CBIO_ideta );
  m_tree->Branch( "CBIO_idphi",    &CBIO_idphi );
  m_tree->Branch( "CBIO_roiNumber",    &CBIO_roiNumber );
  m_tree->Branch( "CBIO_roiSector",    &CBIO_roiSector );
  //EF
  m_tree->Branch( "EF_pt",     &EF_pt );
  m_tree->Branch( "EF_eta",    &EF_eta );
  m_tree->Branch( "EF_phi",    &EF_phi );
  return 1;
}

void NtupleManager::clear() {
  // clear the vector for branch
  std::cout << "NtupleManager::clear() start" <<  std::endl;
  muon_pt->clear();
  muon_eta->clear();
  muon_phi->clear();
  muon_extEta->clear();
  muon_extPhi->clear();
  trigname->clear();
  isPassedTrig->clear();
  isPassedL1_evt->clear();
  isPassedSA_evt->clear();
  isPassedCB_evt->clear();
  isPassedSAIO_evt->clear();
  isPassedCBIO_evt->clear();
  isPassedEF_evt->clear();
  isPassedL1->clear();
  isPassedSA->clear();
  isPassedCB->clear();
  isPassedSAIO->clear();
  isPassedCBIO->clear();
  isPassedEF->clear();
  tdt_L1RoINumber->clear();
  tdt_SARoINumber->clear();
  tdt_CBRoINumber->clear();
  tdt_SAIORoINumber->clear();
  tdt_CBIORoINumber->clear();
  tdt_L1RoISector->clear();
  tdt_SARoISector->clear();
  tdt_CBRoISector->clear();
  tdt_SAIORoISector->clear();
  tdt_CBIORoISector->clear();
  tdt_L1isMoreCandInRoI->clear();
  tdt_SApt->clear();
  tdt_SAeta->clear();
  tdt_SAphi->clear();
  tdt_SAetaMS->clear();
  tdt_SAphiMS->clear();
  tdt_SAsAddress->clear();
  tdt_SAroiEta->clear();
  tdt_SAroiPhi->clear();

  tdt_SAsuperPointR_BI->clear();
  tdt_SAsuperPointR_BM->clear();
  tdt_SAsuperPointR_BO->clear();
  tdt_SAsuperPointR_EI->clear();
  tdt_SAsuperPointR_EM->clear();
  tdt_SAsuperPointR_EO->clear();
  tdt_SAsuperPointR_EE->clear();
  tdt_SAsuperPointR_CSC->clear();
  tdt_SAsuperPointR_BEE->clear();
  tdt_SAsuperPointR_BME->clear();
  tdt_SAsuperPointZ_BI->clear();
  tdt_SAsuperPointZ_BM->clear();
  tdt_SAsuperPointZ_BO->clear();
  tdt_SAsuperPointZ_EI->clear();
  tdt_SAsuperPointZ_EM->clear();
  tdt_SAsuperPointZ_EO->clear();
  tdt_SAsuperPointZ_EE->clear();
  tdt_SAsuperPointZ_CSC->clear();
  tdt_SAsuperPointZ_BEE->clear();
  tdt_SAsuperPointZ_BME->clear();
  tdt_SAIOpt->clear();
  tdt_SAIOeta->clear();
  tdt_SAIOphi->clear();
  tdt_SAIOetaMS->clear();
  tdt_SAIOphiMS->clear();
  tdt_SAIOsAddress->clear();
  tdt_SAIOroiEta->clear();
  tdt_SAIOroiPhi->clear();

  tdt_SAIOsuperPointR_BI->clear();
  tdt_SAIOsuperPointR_BM->clear();
  tdt_SAIOsuperPointR_BO->clear();
  tdt_SAIOsuperPointR_EI->clear();
  tdt_SAIOsuperPointR_EM->clear();
  tdt_SAIOsuperPointR_EO->clear();
  tdt_SAIOsuperPointR_EE->clear();
  tdt_SAIOsuperPointR_CSC->clear();
  tdt_SAIOsuperPointR_BEE->clear();
  tdt_SAIOsuperPointR_BME->clear();
  tdt_SAIOsuperPointZ_BI->clear();
  tdt_SAIOsuperPointZ_BM->clear();
  tdt_SAIOsuperPointZ_BO->clear();
  tdt_SAIOsuperPointZ_EI->clear();
  tdt_SAIOsuperPointZ_EM->clear();
  tdt_SAIOsuperPointZ_EO->clear();
  tdt_SAIOsuperPointZ_EE->clear();
  tdt_SAIOsuperPointZ_CSC->clear();
  tdt_SAIOsuperPointZ_BEE->clear();
  tdt_SAIOsuperPointZ_BME->clear();
  tdt_CBpt->clear();
  tdt_CBeta->clear();
  tdt_CBphi->clear();
  tdt_CBIOpt->clear();
  tdt_CBIOeta->clear();
  tdt_CBIOphi->clear();
  tdt_CBidpt->clear();
  tdt_CBideta->clear();
  tdt_CBidphi->clear();
  tdt_CBIOidpt->clear();
  tdt_CBIOideta->clear();
  tdt_CBIOidphi->clear();
  //
  L1_eta->clear();
  L1_phi->clear();
  L1_thrValue->clear();
  L1_roiNum->clear();
  L1_roiSector->clear();
  L1_thrNumber->clear();
  L1_isMoreCandInRoI->clear();
  //
  SA_roiNum->clear();
  SA_roiSector->clear();
  SA_pt->clear();
  SA_eta->clear();
  SA_phi->clear();
  SA_etaMS->clear();
  SA_phiMS->clear();
  SA_sAddress->clear();
  SA_roiEta->clear();
  SA_roiPhi->clear();

  //the measured radious of the muon in one particular super point
  SA_superPointR_BI->clear();
  SA_superPointR_BM->clear();
  SA_superPointR_BO->clear();
  SA_superPointR_EI->clear();
  SA_superPointR_EM->clear();
  SA_superPointR_EO->clear();
  SA_superPointR_EE->clear();
  SA_superPointR_CSC->clear();
  SA_superPointR_BEE->clear();
  SA_superPointR_BME->clear();
  //the measured Z position of the muon in one particular super point
  SA_superPointZ_BI->clear();
  SA_superPointZ_BM->clear();
  SA_superPointZ_BO->clear();
  SA_superPointZ_EI->clear();
  SA_superPointZ_EM->clear();
  SA_superPointZ_EO->clear();
  SA_superPointZ_EE->clear();
  SA_superPointZ_CSC->clear();
  SA_superPointZ_BEE->clear();
  SA_superPointZ_BME->clear();
  //
  SAIO_roiNum->clear();
  SAIO_roiSector->clear();
  SAIO_pt->clear();
  SAIO_eta->clear();
  SAIO_phi->clear();
  SAIO_etaMS->clear();
  SAIO_phiMS->clear();
  SAIO_sAddress->clear();
  SAIO_roiEta->clear();
  SAIO_roiPhi->clear();

  //the measured radious of the muon in one particular super point
  SAIO_superPointR_BI->clear();
  SAIO_superPointR_BM->clear();
  SAIO_superPointR_BO->clear();
  SAIO_superPointR_EI->clear();
  SAIO_superPointR_EM->clear();
  SAIO_superPointR_EO->clear();
  SAIO_superPointR_EE->clear();
  SAIO_superPointR_CSC->clear();
  SAIO_superPointR_BEE->clear();
  SAIO_superPointR_BME->clear();
  //the measured Z position of the muon in one particular super point
  SAIO_superPointZ_BI->clear();
  SAIO_superPointZ_BM->clear();
  SAIO_superPointZ_BO->clear();
  SAIO_superPointZ_EI->clear();
  SAIO_superPointZ_EM->clear();
  SAIO_superPointZ_EO->clear();
  SAIO_superPointZ_EE->clear();
  SAIO_superPointZ_CSC->clear();
  SAIO_superPointZ_BEE->clear();
  SAIO_superPointZ_BME->clear();
  //
  CB_pt->clear();
  CB_eta->clear();
  CB_phi->clear();
  CB_idpt->clear();
  CB_ideta->clear();
  CB_idphi->clear();
  CB_roiSector->clear();
  CB_roiNumber->clear();
  //
  CBIO_pt->clear();
  CBIO_eta->clear();
  CBIO_phi->clear();
  CBIO_idpt->clear();
  CBIO_ideta->clear();
  CBIO_idphi->clear();
  CBIO_roiSector->clear();
  CBIO_roiNumber->clear();
  //
  EF_pt->clear();
  EF_eta->clear();
  EF_phi->clear();
  std::cout << "NtupleManager::clear() end" <<  std::endl;
}

void NtupleManager::filltree( int evtNum,
                              int runNum,
                              NtupleMakerTool ntupTool )
{
// fill the variable vectors
  this->clear();
  eventNumber = evtNum;
  runNumber = runNum;
  n_trig = ntupTool.m_TDTobjects.size();
  for( int i_tdt = 0; i_tdt < (int)ntupTool.m_TDTobjects.size(); i_tdt++ ) {
    trigname->push_back(ntupTool.m_TDTobjects[i_tdt].trigChainName);
    isPassedTrig->push_back(ntupTool.m_TDTobjects[i_tdt].isPassedChain);
    isPassedL1_evt->push_back(ntupTool.m_TDTobjects[i_tdt].isPassedL1_evt);
    isPassedSA_evt->push_back(ntupTool.m_TDTobjects[i_tdt].isPassedSA_evt);
    isPassedCB_evt->push_back(ntupTool.m_TDTobjects[i_tdt].isPassedCB_evt);
    isPassedSAIO_evt->push_back(ntupTool.m_TDTobjects[i_tdt].isPassedSAIO_evt);
    isPassedCBIO_evt->push_back(ntupTool.m_TDTobjects[i_tdt].isPassedCBIO_evt);
    isPassedEF_evt->push_back(ntupTool.m_TDTobjects[i_tdt].isPassedEF_evt);
    std::vector<bool> L1, isMoreCand, SA, CB, SAIO, CBIO, EF;
    std::vector<int> L1Num, SANum, CBNum, SAIONum, CBIONum;
    std::vector<int> L1Sec, SASec, CBSec, SAIOSec, CBIOSec;
    std::vector < double > tmp_sa_pt;
    std::vector < double > tmp_sa_eta;
    std::vector < double > tmp_sa_phi;
    std::vector < double > tmp_sa_etaMS;
    std::vector < double > tmp_sa_phiMS;
    std::vector < int > tmp_sa_sAddress;
    std::vector < float > tmp_sa_roiEta;
    std::vector < float > tmp_sa_roiPhi;
    std::vector < int > tmp_sa_roiNum;
    std::vector < int > tmp_sa_roiSector;
    std::vector < double > tmp_sa_superPointR_BI;
    std::vector < double > tmp_sa_superPointR_BM;
    std::vector < double > tmp_sa_superPointR_BO;
    std::vector < double > tmp_sa_superPointR_EI;
    std::vector < double > tmp_sa_superPointR_EM;
    std::vector < double > tmp_sa_superPointR_EO;
    std::vector < double > tmp_sa_superPointR_EE;
    std::vector < double > tmp_sa_superPointR_CSC;
    std::vector < double > tmp_sa_superPointR_BEE;
    std::vector < double > tmp_sa_superPointR_BME;
    std::vector < double > tmp_sa_superPointZ_BI;
    std::vector < double > tmp_sa_superPointZ_BM;
    std::vector < double > tmp_sa_superPointZ_BO;
    std::vector < double > tmp_sa_superPointZ_EI;
    std::vector < double > tmp_sa_superPointZ_EM;
    std::vector < double > tmp_sa_superPointZ_EO;
    std::vector < double > tmp_sa_superPointZ_EE;
    std::vector < double > tmp_sa_superPointZ_CSC;
    std::vector < double > tmp_sa_superPointZ_BEE;
    std::vector < double > tmp_sa_superPointZ_BME;
    std::vector < double > tmp_saio_pt;
    std::vector < double > tmp_saio_eta;
    std::vector < double > tmp_saio_phi;
    std::vector < double > tmp_saio_etaMS;
    std::vector < double > tmp_saio_phiMS;
    std::vector < int > tmp_saio_sAddress;
    std::vector < float > tmp_saio_roiEta;
    std::vector < float > tmp_saio_roiPhi;
    std::vector < int > tmp_saio_roiNum;
    std::vector < int > tmp_saio_roiSector;
    std::vector < double > tmp_saio_superPointR_BI;
    std::vector < double > tmp_saio_superPointR_BM;
    std::vector < double > tmp_saio_superPointR_BO;
    std::vector < double > tmp_saio_superPointR_EI;
    std::vector < double > tmp_saio_superPointR_EM;
    std::vector < double > tmp_saio_superPointR_EO;
    std::vector < double > tmp_saio_superPointR_EE;
    std::vector < double > tmp_saio_superPointR_CSC;
    std::vector < double > tmp_saio_superPointR_BEE;
    std::vector < double > tmp_saio_superPointR_BME;
    std::vector < double > tmp_saio_superPointZ_BI;
    std::vector < double > tmp_saio_superPointZ_BM;
    std::vector < double > tmp_saio_superPointZ_BO;
    std::vector < double > tmp_saio_superPointZ_EI;
    std::vector < double > tmp_saio_superPointZ_EM;
    std::vector < double > tmp_saio_superPointZ_EO;
    std::vector < double > tmp_saio_superPointZ_EE;
    std::vector < double > tmp_saio_superPointZ_CSC;
    std::vector < double > tmp_saio_superPointZ_BEE;
    std::vector < double > tmp_saio_superPointZ_BME;
    std::vector < double > tmp_cb_pt;
    std::vector < double > tmp_cb_eta;
    std::vector < double > tmp_cb_phi;
    std::vector < double > tmp_cb_idpt;
    std::vector < double > tmp_cb_ideta;
    std::vector < double > tmp_cb_idphi;
    std::vector < double > tmp_cbio_pt;
    std::vector < double > tmp_cbio_eta;
    std::vector < double > tmp_cbio_phi;
    std::vector < double > tmp_cbio_idpt;
    std::vector < double > tmp_cbio_ideta;
    std::vector < double > tmp_cbio_idphi;
    for(int i_obj=0; i_obj<(int)ntupTool.m_TDTobjects[i_tdt].m_l1.size(); i_obj++){
      L1Object* l1obj = &ntupTool.m_TDTobjects.at(i_tdt).m_l1.at(i_obj);
      L1.push_back(l1obj->isPassed);
      isMoreCand.push_back(l1obj->isMoreCandInRoI);
      L1Num.push_back(l1obj->roiNum);
      L1Sec.push_back(l1obj->roiSector);
    }
    for(int i_obj=0; i_obj<(int)ntupTool.m_TDTobjects[i_tdt].m_sa.size(); i_obj++){
      SAObject* sa = &ntupTool.m_TDTobjects[i_tdt].m_sa.at(i_obj);
      SA.push_back(sa->isPassed);
      SANum.push_back(sa->roiNum);
      SASec.push_back(sa->roiSector);
      tmp_sa_pt.push_back(sa->pt);
      tmp_sa_eta.push_back(sa->eta);
      tmp_sa_phi.push_back(sa->phi);
      tmp_sa_etaMS.push_back(sa->etaMS);
      tmp_sa_phiMS.push_back(sa->phiMS);
      tmp_sa_sAddress.push_back(sa->sAddress);
      tmp_sa_roiEta.push_back(sa->roiEta);
      tmp_sa_roiPhi.push_back(sa->roiPhi);
      tmp_sa_superPointR_BI.push_back(sa->superPointR_BI);
      tmp_sa_superPointR_BM.push_back(sa->superPointR_BM);
      tmp_sa_superPointR_BO.push_back(sa->superPointR_BO);
      tmp_sa_superPointR_EI.push_back(sa->superPointR_EI);
      tmp_sa_superPointR_EM.push_back(sa->superPointR_EM);
      tmp_sa_superPointR_EO.push_back(sa->superPointR_EO);
      tmp_sa_superPointR_EE.push_back(sa->superPointR_EE);
      tmp_sa_superPointR_CSC.push_back(sa->superPointR_CSC);
      tmp_sa_superPointR_BEE.push_back(sa->superPointR_BEE);
      tmp_sa_superPointR_BME.push_back(sa->superPointR_BME);
      tmp_sa_superPointZ_BI.push_back(sa->superPointZ_BI);
      tmp_sa_superPointZ_BM.push_back(sa->superPointZ_BM);
      tmp_sa_superPointZ_BO.push_back(sa->superPointZ_BO);
      tmp_sa_superPointZ_EI.push_back(sa->superPointZ_EI);
      tmp_sa_superPointZ_EM.push_back(sa->superPointZ_EM);
      tmp_sa_superPointZ_EO.push_back(sa->superPointZ_EO);
      tmp_sa_superPointZ_EE.push_back(sa->superPointZ_EE);
      tmp_sa_superPointZ_CSC.push_back(sa->superPointZ_CSC);
      tmp_sa_superPointZ_BEE.push_back(sa->superPointZ_BEE);
      tmp_sa_superPointZ_BME.push_back(sa->superPointZ_BME);
    }
    for(int i_obj=0; i_obj<(int)ntupTool.m_TDTobjects[i_tdt].m_cb.size(); i_obj++){
      CBObject* cb = &ntupTool.m_TDTobjects[i_tdt].m_cb.at(i_obj);
      CB.push_back(cb->isPassed);
      CBNum.push_back(cb->roiNum);
      CBSec.push_back(cb->roiSector);
      tmp_cb_pt.push_back(cb->pt);
      tmp_cb_eta.push_back(cb->eta);
      tmp_cb_phi.push_back(cb->phi);
      tmp_cb_idpt.push_back(cb->idpt);
      tmp_cb_ideta.push_back(cb->ideta);
      tmp_cb_idphi.push_back(cb->idphi);
      CB.push_back(ntupTool.m_TDTobjects.at(i_tdt).m_cb.at(i_obj).isPassed);
      CBNum.push_back(ntupTool.m_TDTobjects.at(i_tdt).m_cb.at(i_obj).roiNum);
      CBSec.push_back(ntupTool.m_TDTobjects.at(i_tdt).m_cb.at(i_obj).roiSector);
    }
    for(int i_obj=0; i_obj<(int)ntupTool.m_TDTobjects[i_tdt].m_saio.size(); i_obj++){
      SAObject* saio = &ntupTool.m_TDTobjects[i_tdt].m_saio.at(i_obj);
      SAIO.push_back(saio->isPassed);
      SAIONum.push_back(saio->roiNum);
      SAIOSec.push_back(saio->roiSector);
      tmp_saio_pt.push_back(saio->pt);
      tmp_saio_eta.push_back(saio->eta);
      tmp_saio_phi.push_back(saio->phi);
      tmp_saio_etaMS.push_back(saio->etaMS);
      tmp_saio_phiMS.push_back(saio->phiMS);
      tmp_saio_sAddress.push_back(saio->sAddress);
      tmp_saio_roiEta.push_back(saio->roiEta);
      tmp_saio_roiPhi.push_back(saio->roiPhi);
      tmp_saio_superPointR_BI.push_back(saio->superPointR_BI);
      tmp_saio_superPointR_BM.push_back(saio->superPointR_BM);
      tmp_saio_superPointR_BO.push_back(saio->superPointR_BO);
      tmp_saio_superPointR_EI.push_back(saio->superPointR_EI);
      tmp_saio_superPointR_EM.push_back(saio->superPointR_EM);
      tmp_saio_superPointR_EO.push_back(saio->superPointR_EO);
      tmp_saio_superPointR_EE.push_back(saio->superPointR_EE);
      tmp_saio_superPointR_CSC.push_back(saio->superPointR_CSC);
      tmp_saio_superPointR_BEE.push_back(saio->superPointR_BEE);
      tmp_saio_superPointR_BME.push_back(saio->superPointR_BME);
      tmp_saio_superPointZ_BI.push_back(saio->superPointZ_BI);
      tmp_saio_superPointZ_BM.push_back(saio->superPointZ_BM);
      tmp_saio_superPointZ_BO.push_back(saio->superPointZ_BO);
      tmp_saio_superPointZ_EI.push_back(saio->superPointZ_EI);
      tmp_saio_superPointZ_EM.push_back(saio->superPointZ_EM);
      tmp_saio_superPointZ_EO.push_back(saio->superPointZ_EO);
      tmp_saio_superPointZ_EE.push_back(saio->superPointZ_EE);
      tmp_saio_superPointZ_CSC.push_back(saio->superPointZ_CSC);
      tmp_saio_superPointZ_BEE.push_back(saio->superPointZ_BEE);
      tmp_saio_superPointZ_BME.push_back(saio->superPointZ_BME);
    }
    for(int i_obj=0; i_obj<(int)ntupTool.m_TDTobjects[i_tdt].m_cbio.size(); i_obj++){
      CBObject* cbio = &ntupTool.m_TDTobjects[i_tdt].m_cbio.at(i_obj);
      CBIO.push_back(cbio->isPassed);
      CBIONum.push_back(cbio->roiNum);
      CBIOSec.push_back(cbio->roiSector);
      tmp_cbio_pt.push_back(cbio->pt);
      tmp_cbio_eta.push_back(cbio->eta);
      tmp_cbio_phi.push_back(cbio->phi);
      tmp_cbio_idpt.push_back(cbio->idpt);
      tmp_cbio_ideta.push_back(cbio->ideta);
      tmp_cbio_idphi.push_back(cbio->idphi);
    }
    for(int i_obj=0; i_obj<(int)ntupTool.m_TDTobjects[i_tdt].m_ef.size(); i_obj++){
      EFObject* ef = &ntupTool.m_TDTobjects.at(i_tdt).m_ef.at(i_obj);
      EF.push_back(ef->isPassed);
    }
    isPassedL1->push_back(L1);
    isPassedSA->push_back(SA);
    isPassedCB->push_back(CB);
    isPassedSAIO->push_back(SAIO);
    isPassedCBIO->push_back(CBIO);
    isPassedEF->push_back(EF);
    tdt_L1RoINumber->push_back(L1Num);
    tdt_SARoINumber->push_back(SANum);
    tdt_CBRoINumber->push_back(CBNum);
    tdt_SAIORoINumber->push_back(SAIONum);
    tdt_CBIORoINumber->push_back(CBIONum);
    tdt_L1RoISector->push_back(L1Sec);
    tdt_SARoISector->push_back(SASec);
    tdt_CBRoISector->push_back(CBSec);
    tdt_SAIORoISector->push_back(SAIOSec);
    tdt_CBIORoISector->push_back(CBIOSec);
    tdt_L1isMoreCandInRoI->push_back(isMoreCand);
    tdt_SApt->push_back(tmp_sa_pt);
    tdt_SAeta->push_back(tmp_sa_eta);
    tdt_SAphi->push_back(tmp_sa_phi);
    tdt_SAetaMS->push_back(tmp_sa_etaMS);
    tdt_SAphiMS->push_back(tmp_sa_phiMS);
    tdt_SAsAddress->push_back(tmp_sa_sAddress);
    tdt_SAroiEta->push_back(tmp_sa_roiEta);
    tdt_SAroiPhi->push_back(tmp_sa_roiPhi);
    tdt_SAsuperPointR_BI->push_back(tmp_sa_superPointR_BI);
    tdt_SAsuperPointR_BM->push_back(tmp_sa_superPointR_BM);
    tdt_SAsuperPointR_BO->push_back(tmp_sa_superPointR_BO);
    tdt_SAsuperPointR_EI->push_back(tmp_sa_superPointR_EI);
    tdt_SAsuperPointR_EM->push_back(tmp_sa_superPointR_EM);
    tdt_SAsuperPointR_EO->push_back(tmp_sa_superPointR_EO);
    tdt_SAsuperPointR_EE->push_back(tmp_sa_superPointR_EE);
    tdt_SAsuperPointR_CSC->push_back(tmp_sa_superPointR_CSC);
    tdt_SAsuperPointR_BEE->push_back(tmp_sa_superPointR_BEE);
    tdt_SAsuperPointR_BME->push_back(tmp_sa_superPointR_BME);
    tdt_SAsuperPointZ_BI->push_back(tmp_sa_superPointZ_BI);
    tdt_SAsuperPointZ_BM->push_back(tmp_sa_superPointZ_BM);
    tdt_SAsuperPointZ_BO->push_back(tmp_sa_superPointZ_BO);
    tdt_SAsuperPointZ_EI->push_back(tmp_sa_superPointZ_EI);
    tdt_SAsuperPointZ_EM->push_back(tmp_sa_superPointZ_EM);
    tdt_SAsuperPointZ_EO->push_back(tmp_sa_superPointZ_EO);
    tdt_SAsuperPointZ_EE->push_back(tmp_sa_superPointZ_EE);
    tdt_SAsuperPointZ_CSC->push_back(tmp_sa_superPointZ_CSC);
    tdt_SAsuperPointZ_BEE->push_back(tmp_sa_superPointZ_BEE);
    tdt_SAsuperPointZ_BME->push_back(tmp_sa_superPointZ_BME);
    tdt_SAIOpt->push_back(tmp_saio_pt);
    tdt_SAIOeta->push_back(tmp_saio_eta);
    tdt_SAIOphi->push_back(tmp_saio_phi);
    tdt_SAIOetaMS->push_back(tmp_saio_etaMS);
    tdt_SAIOphiMS->push_back(tmp_saio_phiMS);
    tdt_SAIOsAddress->push_back(tmp_saio_sAddress);
    tdt_SAIOroiEta->push_back(tmp_saio_roiEta);
    tdt_SAIOroiPhi->push_back(tmp_saio_roiPhi);
    tdt_SAIOsuperPointR_BI->push_back(tmp_saio_superPointR_BI);
    tdt_SAIOsuperPointR_BM->push_back(tmp_saio_superPointR_BM);
    tdt_SAIOsuperPointR_BO->push_back(tmp_saio_superPointR_BO);
    tdt_SAIOsuperPointR_EI->push_back(tmp_saio_superPointR_EI);
    tdt_SAIOsuperPointR_EM->push_back(tmp_saio_superPointR_EM);
    tdt_SAIOsuperPointR_EO->push_back(tmp_saio_superPointR_EO);
    tdt_SAIOsuperPointR_EE->push_back(tmp_saio_superPointR_EE);
    tdt_SAIOsuperPointR_CSC->push_back(tmp_saio_superPointR_CSC);
    tdt_SAIOsuperPointR_BEE->push_back(tmp_saio_superPointR_BEE);
    tdt_SAIOsuperPointR_BME->push_back(tmp_saio_superPointR_BME);
    tdt_SAIOsuperPointZ_BI->push_back(tmp_saio_superPointZ_BI);
    tdt_SAIOsuperPointZ_BM->push_back(tmp_saio_superPointZ_BM);
    tdt_SAIOsuperPointZ_BO->push_back(tmp_saio_superPointZ_BO);
    tdt_SAIOsuperPointZ_EI->push_back(tmp_saio_superPointZ_EI);
    tdt_SAIOsuperPointZ_EM->push_back(tmp_saio_superPointZ_EM);
    tdt_SAIOsuperPointZ_EO->push_back(tmp_saio_superPointZ_EO);
    tdt_SAIOsuperPointZ_EE->push_back(tmp_saio_superPointZ_EE);
    tdt_SAIOsuperPointZ_CSC->push_back(tmp_saio_superPointZ_CSC);
    tdt_SAIOsuperPointZ_BEE->push_back(tmp_saio_superPointZ_BEE);
    tdt_SAIOsuperPointZ_BME->push_back(tmp_saio_superPointZ_BME);
    tdt_CBpt->push_back(tmp_cb_pt);
    tdt_CBeta->push_back(tmp_cb_eta);
    tdt_CBphi->push_back(tmp_cb_phi);
    tdt_CBidpt->push_back(tmp_cb_idpt);
    tdt_CBideta->push_back(tmp_cb_ideta);
    tdt_CBidphi->push_back(tmp_cb_idphi);
    tdt_CBIOpt->push_back(tmp_cbio_pt);
    tdt_CBIOeta->push_back(tmp_cbio_eta);
    tdt_CBIOphi->push_back(tmp_cbio_phi);
    tdt_CBIOidpt->push_back(tmp_cbio_idpt);
    tdt_CBIOideta->push_back(tmp_cbio_ideta);
    tdt_CBIOidphi->push_back(tmp_cbio_idphi);
  }
  for( int i_muon = 0; i_muon < (int)ntupTool.m_muons.size(); i_muon++ ) {
    muon_pt->push_back(ntupTool.m_muons[i_muon].pt);
    muon_eta->push_back(ntupTool.m_muons[i_muon].eta);
    muon_phi->push_back(ntupTool.m_muons[i_muon].phi);
    muon_extEta->push_back(ntupTool.m_muons[i_muon].extEta);
    muon_extPhi->push_back(ntupTool.m_muons[i_muon].extPhi);
  }
  for( int i_l1 = 0; i_l1 < (int)ntupTool.m_L1objects.size(); i_l1++ ) {
    L1_thrValue->push_back(ntupTool.m_L1objects[i_l1].thrValue);
    L1_thrNumber->push_back(ntupTool.m_L1objects[i_l1].thrNumber);
    L1_eta->push_back(ntupTool.m_L1objects[i_l1].eta);
    L1_phi->push_back(ntupTool.m_L1objects[i_l1].phi);
    L1_roiNum->push_back(ntupTool.m_L1objects[i_l1].roiNum);
    L1_roiSector->push_back(ntupTool.m_L1objects[i_l1].roiSector);
    L1_isMoreCandInRoI->push_back(ntupTool.m_L1objects[i_l1].isMoreCandInRoI);
  }
  for( int i_sa = 0; i_sa < (int)ntupTool.m_SAobjects.size(); i_sa++ ) {
    SA_pt->push_back(ntupTool.m_SAobjects[i_sa].pt);
    SA_eta->push_back(ntupTool.m_SAobjects[i_sa].eta);
    SA_phi->push_back(ntupTool.m_SAobjects[i_sa].phi);
    SA_roiEta->push_back(ntupTool.m_SAobjects[i_sa].roiEta);
    SA_roiPhi->push_back(ntupTool.m_SAobjects[i_sa].roiPhi);
    SA_etaMS->push_back(ntupTool.m_SAobjects[i_sa].etaMS);
    SA_phiMS->push_back(ntupTool.m_SAobjects[i_sa].phiMS);
    SA_roiNum->push_back(ntupTool.m_SAobjects[i_sa].roiNum);
    SA_roiSector->push_back(ntupTool.m_SAobjects[i_sa].roiSector);
    //the measured radious of the muon in one particular super point
    SA_superPointR_BI->push_back( ntupTool.m_SAobjects[i_sa].superPointR_BI );
    SA_superPointR_BM->push_back( ntupTool.m_SAobjects[i_sa].superPointR_BM );
    SA_superPointR_BO->push_back( ntupTool.m_SAobjects[i_sa].superPointR_BO );
    SA_superPointR_EI->push_back( ntupTool.m_SAobjects[i_sa].superPointR_EI );
    SA_superPointR_EM->push_back( ntupTool.m_SAobjects[i_sa].superPointR_EM );
    SA_superPointR_EO->push_back( ntupTool.m_SAobjects[i_sa].superPointR_EO );
    SA_superPointR_EE->push_back( ntupTool.m_SAobjects[i_sa].superPointR_EE );
    SA_superPointR_CSC->push_back( ntupTool.m_SAobjects[i_sa].superPointR_CSC );
    SA_superPointR_BEE->push_back( ntupTool.m_SAobjects[i_sa].superPointR_BEE );
    SA_superPointR_BME->push_back( ntupTool.m_SAobjects[i_sa].superPointR_BME );
    //the measured Z position of the muon in one particular super point
    SA_superPointZ_BI->push_back( ntupTool.m_SAobjects[i_sa].superPointZ_BI );
    SA_superPointZ_BM->push_back( ntupTool.m_SAobjects[i_sa].superPointZ_BM );
    SA_superPointZ_BO->push_back( ntupTool.m_SAobjects[i_sa].superPointZ_BO );
    SA_superPointZ_EI->push_back( ntupTool.m_SAobjects[i_sa].superPointZ_EI );
    SA_superPointZ_EM->push_back( ntupTool.m_SAobjects[i_sa].superPointZ_EM );
    SA_superPointZ_EO->push_back( ntupTool.m_SAobjects[i_sa].superPointZ_EO );
    SA_superPointZ_EE->push_back( ntupTool.m_SAobjects[i_sa].superPointZ_EE );
    SA_superPointZ_CSC->push_back( ntupTool.m_SAobjects[i_sa].superPointZ_CSC );
    SA_superPointZ_BEE->push_back( ntupTool.m_SAobjects[i_sa].superPointZ_BEE );
    SA_superPointZ_BME->push_back( ntupTool.m_SAobjects[i_sa].superPointZ_BME );
  }
  for( int i_sa = 0; i_sa < (int)ntupTool.m_SAIOobjects.size(); i_sa++ ) {
    SAIO_pt->push_back(ntupTool.m_SAIOobjects[i_sa].pt);
    SAIO_eta->push_back(ntupTool.m_SAIOobjects[i_sa].eta);
    SAIO_phi->push_back(ntupTool.m_SAIOobjects[i_sa].phi);
    SAIO_roiEta->push_back(ntupTool.m_SAIOobjects[i_sa].roiEta);
    SAIO_roiPhi->push_back(ntupTool.m_SAIOobjects[i_sa].roiPhi);
    SAIO_etaMS->push_back(ntupTool.m_SAIOobjects[i_sa].etaMS);
    SAIO_phiMS->push_back(ntupTool.m_SAIOobjects[i_sa].phiMS);
    SAIO_roiNum->push_back(ntupTool.m_SAIOobjects[i_sa].roiNum);
    SAIO_roiSector->push_back(ntupTool.m_SAIOobjects[i_sa].roiSector);
    SAIO_superPointR_BI->push_back( ntupTool.m_SAIOobjects[i_sa].superPointR_BI );
    SAIO_superPointR_BM->push_back( ntupTool.m_SAIOobjects[i_sa].superPointR_BM );
    SAIO_superPointR_BO->push_back( ntupTool.m_SAIOobjects[i_sa].superPointR_BO );
    SAIO_superPointR_EI->push_back( ntupTool.m_SAIOobjects[i_sa].superPointR_EI );
    SAIO_superPointR_EM->push_back( ntupTool.m_SAIOobjects[i_sa].superPointR_EM );
    SAIO_superPointR_EO->push_back( ntupTool.m_SAIOobjects[i_sa].superPointR_EO );
    SAIO_superPointR_EE->push_back( ntupTool.m_SAIOobjects[i_sa].superPointR_EE );
    SAIO_superPointR_CSC->push_back( ntupTool.m_SAIOobjects[i_sa].superPointR_CSC );
    SAIO_superPointR_BEE->push_back( ntupTool.m_SAIOobjects[i_sa].superPointR_BEE );
    SAIO_superPointR_BME->push_back( ntupTool.m_SAIOobjects[i_sa].superPointR_BME );
    SAIO_superPointZ_BI->push_back( ntupTool.m_SAIOobjects[i_sa].superPointZ_BI );
    SAIO_superPointZ_BM->push_back( ntupTool.m_SAIOobjects[i_sa].superPointZ_BM );
    SAIO_superPointZ_BO->push_back( ntupTool.m_SAIOobjects[i_sa].superPointZ_BO );
    SAIO_superPointZ_EI->push_back( ntupTool.m_SAIOobjects[i_sa].superPointZ_EI );
    SAIO_superPointZ_EM->push_back( ntupTool.m_SAIOobjects[i_sa].superPointZ_EM );
    SAIO_superPointZ_EO->push_back( ntupTool.m_SAIOobjects[i_sa].superPointZ_EO );
    SAIO_superPointZ_EE->push_back( ntupTool.m_SAIOobjects[i_sa].superPointZ_EE );
    SAIO_superPointZ_CSC->push_back( ntupTool.m_SAIOobjects[i_sa].superPointZ_CSC );
    SAIO_superPointZ_BEE->push_back( ntupTool.m_SAIOobjects[i_sa].superPointZ_BEE );
    SAIO_superPointZ_BME->push_back( ntupTool.m_SAIOobjects[i_sa].superPointZ_BME );
  }
  for( int i_CB = 0; i_CB < (int)ntupTool.m_CBobjects.size(); i_CB++ ) {
    CB_pt->push_back(ntupTool.m_CBobjects[i_CB].pt);
    CB_eta->push_back(ntupTool.m_CBobjects[i_CB].eta);
    CB_phi->push_back(ntupTool.m_CBobjects[i_CB].phi);
    CB_idpt->push_back(ntupTool.m_CBobjects[i_CB].idpt);
    CB_ideta->push_back(ntupTool.m_CBobjects[i_CB].ideta);
    CB_idphi->push_back(ntupTool.m_CBobjects[i_CB].idphi);
    CB_roiNumber->push_back(ntupTool.m_CBobjects[i_CB].roiNum);
    CB_roiSector->push_back(ntupTool.m_CBobjects[i_CB].roiSector);
  }
  for( int i_CBIO = 0; i_CBIO < (int)ntupTool.m_CBIOobjects.size(); i_CBIO++ ) {
    CBIO_pt->push_back(ntupTool.m_CBIOobjects[i_CBIO].pt);
    CBIO_eta->push_back(ntupTool.m_CBIOobjects[i_CBIO].eta);
    CBIO_phi->push_back(ntupTool.m_CBIOobjects[i_CBIO].phi);
    CBIO_idpt->push_back(ntupTool.m_CBIOobjects[i_CBIO].idpt);
    CBIO_ideta->push_back(ntupTool.m_CBIOobjects[i_CBIO].ideta);
    CBIO_idphi->push_back(ntupTool.m_CBIOobjects[i_CBIO].idphi);
    CBIO_roiNumber->push_back(ntupTool.m_CBIOobjects[i_CBIO].roiNum);
    CBIO_roiSector->push_back(ntupTool.m_CBIOobjects[i_CBIO].roiSector);
  }
  m_tree->Fill();
  std::cout << "NtupleManager::filltree  end" << std::endl;
}

int NtupleManager::finalize() {
  // write into file and close	
  m_file->Write();	
  return 1;
}            
