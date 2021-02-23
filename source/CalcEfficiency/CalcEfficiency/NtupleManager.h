#ifndef _TrigMuonAnalysis_NtupleManager
#define _TrigMuonAnalysis_NtupleManager

#include <iostream>
#include <fstream>
#include <vector>
#include <stdint.h>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TLorentzVector.h"

#include "CalcEfficiency/NtupleMakerTool.h"
#include "CalcEfficiency/OfflineObjects.h"
#include "CalcEfficiency/L1Objects.h"
#include "CalcEfficiency/SAObjects.h"
#include "CalcEfficiency/CBObjects.h"
#include "CalcEfficiency/EFObjects.h"
#include "CalcEfficiency/TDTObjects.h"
#include "CalcEfficiency/Utils.h"

class NtupleManager {
	
  public:
    NtupleManager();
    ~NtupleManager();

    TFile* m_file; //!
    TTree* m_tree; //!
    OfflineObjects m_muon;
    TDTObjects m_TDTobjects;
    L1Objects m_L1objects;
    SAObjects m_SAobjects;
    CBObjects m_CBobjects;
    EFObjects m_EFobjects;
  
  public:
    int initialize( std::string outfile );
    void clear();
    void filltree( int evtNum, 
                   int runNum,
                   NtupleMakerTool ntupTool );
    int finalize();
    
    //event info
    int eventNumber;
    int runNumber;
    //
    std::vector<std::string>* trigname;
    std::vector<bool>* isPassedTrig;
    std::vector<bool>* isPassedL1_evt;
    std::vector<bool>* isPassedSA_evt;
    std::vector<bool>* isPassedCB_evt;
    std::vector<bool>* isPassedSAIO_evt;
    std::vector<bool>* isPassedCBIO_evt;
    std::vector<bool>* isPassedEF_evt;
    std::vector<std::vector<bool>>* isPassedL1;
    std::vector<std::vector<bool>>* isPassedSA;
    std::vector<std::vector<bool>>* isPassedCB;
    std::vector<std::vector<bool>>* isPassedSAIO;
    std::vector<std::vector<bool>>* isPassedCBIO;
    std::vector<std::vector<bool>>* isPassedEF;
    std::vector<std::vector<int>>* tdt_L1RoINumber;
    std::vector<std::vector<int>>* tdt_SARoINumber;
    std::vector<std::vector<int>>* tdt_CBRoINumber;
    std::vector<std::vector<int>>* tdt_SAIORoINumber;
    std::vector<std::vector<int>>* tdt_CBIORoINumber;
    std::vector<std::vector<int>>* tdt_L1RoISector;
    std::vector<std::vector<int>>* tdt_SARoISector;
    std::vector<std::vector<int>>* tdt_CBRoISector;
    std::vector<std::vector<int>>* tdt_SAIORoISector;
    std::vector<std::vector<int>>* tdt_CBIORoISector;
    std::vector<std::vector<bool>>* tdt_L1isMoreCandInRoI;
    std::vector < std::vector<double> >* tdt_SApt;
    std::vector < std::vector<double> >* tdt_SAeta;
    std::vector < std::vector<double> >* tdt_SAphi;
    std::vector < std::vector<double> >* tdt_SAetaMS;
    std::vector < std::vector<double> >* tdt_SAphiMS;
    std::vector < std::vector<int> >* tdt_SAsAddress;
    std::vector < std::vector<float> >* tdt_SAroiEta;
    std::vector < std::vector<float> >* tdt_SAroiPhi;
    std::vector < std::vector<double> >* tdt_SAsuperPointR_BI;
    std::vector < std::vector<double> >* tdt_SAsuperPointR_BM;
    std::vector < std::vector<double> >* tdt_SAsuperPointR_BO;
    std::vector < std::vector<double> >* tdt_SAsuperPointR_EI;
    std::vector < std::vector<double> >* tdt_SAsuperPointR_EM;
    std::vector < std::vector<double> >* tdt_SAsuperPointR_EO;
    std::vector < std::vector<double> >* tdt_SAsuperPointR_EE;
    std::vector < std::vector<double> >* tdt_SAsuperPointR_CSC;
    std::vector < std::vector<double> >* tdt_SAsuperPointR_BEE;
    std::vector < std::vector<double> >* tdt_SAsuperPointR_BME;
    //the measured Z position of the muon in one particular super point
    std::vector < std::vector<double> >* tdt_SAsuperPointZ_BI;
    std::vector < std::vector<double> >* tdt_SAsuperPointZ_BM;
    std::vector < std::vector<double> >* tdt_SAsuperPointZ_BO;
    std::vector < std::vector<double> >* tdt_SAsuperPointZ_EI;
    std::vector < std::vector<double> >* tdt_SAsuperPointZ_EM;
    std::vector < std::vector<double> >* tdt_SAsuperPointZ_EO;
    std::vector < std::vector<double> >* tdt_SAsuperPointZ_EE;
    std::vector < std::vector<double> >* tdt_SAsuperPointZ_CSC;
    std::vector < std::vector<double> >* tdt_SAsuperPointZ_BEE;
    std::vector < std::vector<double> >* tdt_SAsuperPointZ_BME;
    std::vector < std::vector<double> >* tdt_SAIOpt;
    std::vector < std::vector<double> >* tdt_SAIOeta;
    std::vector < std::vector<double> >* tdt_SAIOphi;
    std::vector < std::vector<double> >* tdt_SAIOetaMS;
    std::vector < std::vector<double> >* tdt_SAIOphiMS;
    std::vector < std::vector<int> >* tdt_SAIOsAddress;
    std::vector < std::vector<float> >* tdt_SAIOroiEta;
    std::vector < std::vector<float> >* tdt_SAIOroiPhi;
    std::vector < std::vector<double> >* tdt_SAIOsuperPointR_BI;
    std::vector < std::vector<double> >* tdt_SAIOsuperPointR_BM;
    std::vector < std::vector<double> >* tdt_SAIOsuperPointR_BO;
    std::vector < std::vector<double> >* tdt_SAIOsuperPointR_EI;
    std::vector < std::vector<double> >* tdt_SAIOsuperPointR_EM;
    std::vector < std::vector<double> >* tdt_SAIOsuperPointR_EO;
    std::vector < std::vector<double> >* tdt_SAIOsuperPointR_EE;
    std::vector < std::vector<double> >* tdt_SAIOsuperPointR_CSC;
    std::vector < std::vector<double> >* tdt_SAIOsuperPointR_BEE;
    std::vector < std::vector<double> >* tdt_SAIOsuperPointR_BME;
    //the measured Z position of the muon in one particular super point
    std::vector < std::vector<double> >* tdt_SAIOsuperPointZ_BI;
    std::vector < std::vector<double> >* tdt_SAIOsuperPointZ_BM;
    std::vector < std::vector<double> >* tdt_SAIOsuperPointZ_BO;
    std::vector < std::vector<double> >* tdt_SAIOsuperPointZ_EI;
    std::vector < std::vector<double> >* tdt_SAIOsuperPointZ_EM;
    std::vector < std::vector<double> >* tdt_SAIOsuperPointZ_EO;
    std::vector < std::vector<double> >* tdt_SAIOsuperPointZ_EE;
    std::vector < std::vector<double> >* tdt_SAIOsuperPointZ_CSC;
    std::vector < std::vector<double> >* tdt_SAIOsuperPointZ_BEE;
    std::vector < std::vector<double> >* tdt_SAIOsuperPointZ_BME;
    std::vector < std::vector<double> >* tdt_CBIOpt;
    std::vector < std::vector<double> >* tdt_CBIOeta;
    std::vector < std::vector<double> >* tdt_CBIOphi;
    std::vector < std::vector<double> >* tdt_CBIOidpt;
    std::vector < std::vector<double> >* tdt_CBIOideta;
    std::vector < std::vector<double> >* tdt_CBIOidphi;
    std::vector < std::vector<double> >* tdt_CBpt;
    std::vector < std::vector<double> >* tdt_CBeta;
    std::vector < std::vector<double> >* tdt_CBphi;
    std::vector < std::vector<double> >* tdt_CBidpt;
    std::vector < std::vector<double> >* tdt_CBideta;
    std::vector < std::vector<double> >* tdt_CBidphi;
    int n_trig;
    //
    std::vector< double >* muon_pt;
    std::vector< double >* muon_eta;
    std::vector< double >* muon_phi;
    std::vector< double >* muon_extEta;
    std::vector< double >* muon_extPhi;
    //
    std::vector < double >* L1_eta;
    std::vector < double >* L1_phi;
    std::vector < double >* L1_thrValue;
    std::vector < int >* L1_roiNum;
    std::vector < int >* L1_roiSector;
    std::vector < int >* L1_thrNumber;
    std::vector < bool >* L1_isMoreCandInRoI;
    //
    std::vector < double >* SA_pt;
    std::vector < double >* SA_eta;
    std::vector < double >* SA_phi;
    std::vector < double >* SA_etaMS;
    std::vector < double >* SA_phiMS;
    std::vector < int >* SA_sAddress;
    std::vector < float >* SA_roiEta;
    std::vector < float >* SA_roiPhi;
    std::vector < int >* SA_roiNum;
    std::vector < int >* SA_roiSector;
    //the measured radious of the muon in one particular super point
    std::vector < double >* SA_superPointR_BI;
    std::vector < double >* SA_superPointR_BM;
    std::vector < double >* SA_superPointR_BO;
    std::vector < double >* SA_superPointR_EI;
    std::vector < double >* SA_superPointR_EM;
    std::vector < double >* SA_superPointR_EO;
    std::vector < double >* SA_superPointR_EE;
    std::vector < double >* SA_superPointR_CSC;
    std::vector < double >* SA_superPointR_BEE;
    std::vector < double >* SA_superPointR_BME;
    //the measured Z position of the muon in one particular super point
    std::vector < double >* SA_superPointZ_BI;
    std::vector < double >* SA_superPointZ_BM;
    std::vector < double >* SA_superPointZ_BO;
    std::vector < double >* SA_superPointZ_EI;
    std::vector < double >* SA_superPointZ_EM;
    std::vector < double >* SA_superPointZ_EO;
    std::vector < double >* SA_superPointZ_EE;
    std::vector < double >* SA_superPointZ_CSC;
    std::vector < double >* SA_superPointZ_BEE;
    std::vector < double >* SA_superPointZ_BME;
    //
    std::vector < double >* SAIO_pt;
    std::vector < double >* SAIO_eta;
    std::vector < double >* SAIO_phi;
    std::vector < double >* SAIO_etaMS;
    std::vector < double >* SAIO_phiMS;
    std::vector < int >* SAIO_sAddress;
    std::vector < float >* SAIO_roiEta;
    std::vector < float >* SAIO_roiPhi;
    std::vector < int >* SAIO_roiNum;
    std::vector < int >* SAIO_roiSector;
    //the measured radious of the muon in one particular super point
    std::vector < double >* SAIO_superPointR_BI;
    std::vector < double >* SAIO_superPointR_BM;
    std::vector < double >* SAIO_superPointR_BO;
    std::vector < double >* SAIO_superPointR_EI;
    std::vector < double >* SAIO_superPointR_EM;
    std::vector < double >* SAIO_superPointR_EO;
    std::vector < double >* SAIO_superPointR_EE;
    std::vector < double >* SAIO_superPointR_CSC;
    std::vector < double >* SAIO_superPointR_BEE;
    std::vector < double >* SAIO_superPointR_BME;
    //the measured Z position of the muon in one particular super point
    std::vector < double >* SAIO_superPointZ_BI;
    std::vector < double >* SAIO_superPointZ_BM;
    std::vector < double >* SAIO_superPointZ_BO;
    std::vector < double >* SAIO_superPointZ_EI;
    std::vector < double >* SAIO_superPointZ_EM;
    std::vector < double >* SAIO_superPointZ_EO;
    std::vector < double >* SAIO_superPointZ_EE;
    std::vector < double >* SAIO_superPointZ_CSC;
    std::vector < double >* SAIO_superPointZ_BEE;
    std::vector < double >* SAIO_superPointZ_BME;
    //
    std::vector < double >* CB_pt;
    std::vector < double >* CB_eta;
    std::vector < double >* CB_phi;
    std::vector < double >* CB_idpt;
    std::vector < double >* CB_ideta;
    std::vector < double >* CB_idphi;
    std::vector < int >* CB_roiNumber;
    std::vector < int >* CB_roiSector;
    //
    std::vector < double >* CBIO_pt;
    std::vector < double >* CBIO_eta;
    std::vector < double >* CBIO_phi;
    std::vector < double >* CBIO_idpt;
    std::vector < double >* CBIO_ideta;
    std::vector < double >* CBIO_idphi;
    std::vector < int >* CBIO_roiNumber;
    std::vector < int >* CBIO_roiSector;
    //
    std::vector < double >* EF_pt;
    std::vector < double >* EF_eta;
    std::vector < double >* EF_phi;

};

#endif

