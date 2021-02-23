#ifndef HISTNTUPLEMT_H
#define HISTNTUPLEMT_H

// ASG Tools
#include "AthenaBaseComps/AthAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"
#include "GaudiKernel/ITHistSvc.h"
#include "TFile.h"
#include "TEfficiency.h"
#include "CalcEfficiency/TagAndProbeMT.h"
#include "CalcEfficiency/TagAndProbe.h"

class TTree;
class TH1D;
class TH2D;

class HistNtupleMT {
  public:
    HistNtupleMT();
    virtual ~HistNtupleMT();

    template <typename TAP>
    int initialize( std::string histname, TAP& tap, const bool isAsymNSW, std::string montrigName="HLT_mu6_L1MU6" );
    template <typename TAP> int finalize( TAP& tap );
    template <typename TAP> void FillHist( TAP& tap );
    void CalcEfficiency( TH1D* h_num, TH1D* h_den, TH1D* h_set );
    inline float calc_residual(float x, float y, float aw, float bw);

    bool m_isAsymNSW{false};
    int m_monChain{0};
    std::string m_monTrigName{"HLT_mu6_L1MU6"};

  public:
    TFile* m_FILE;
    //check valiables
    TH1D* m_h_probept;
    TH1D* m_h_probeexteta;
    TH1D* m_h_probeextphi;
    TH2D* m_h_stgcClusterZR; //should be "m"
    TH2D* m_h_mmClusterZR;
    TH1D* m_h_residual_offsegvsinnSP;
    TH1D* m_h_dtheta_offvsinnSP;
    TH1D* m_h_residual_roadvsinnSP;
    // for isMoreCand study 
    TH1D* m_h_passedisMoreCand;
    //efficiency
    TH1D* m_h_probeEvents_offpt;
    std::vector< TH1D* > m_h_L1pass_offpt;
    std::vector< TH1D* > m_h_SApass_offpt;
    std::vector< TH1D* > m_h_CBpass_offpt;
    std::vector< TH1D* > m_h_EFpass_offpt;
    TEfficiency* m_peffL1_offpt;
    std::vector< TH1D* > m_eff_L1pass_offpt;
    std::vector< TH1D* > m_eff_SApass_offpt;
    std::vector< TH1D* > m_eff_CBpass_offpt;
    std::vector< TH1D* > m_eff_EFpass_offpt;
    std::vector< TH1D* > m_h_trigPassEvents;
};
#endif  //HISTNTUPLEMT_H
