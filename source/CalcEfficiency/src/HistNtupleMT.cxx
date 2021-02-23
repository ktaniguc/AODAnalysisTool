//Gaudi
#include "GaudiKernel/ITHistSvc.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ToolHandle.h"
#include "StoreGate/StoreGateSvc.h"

#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TVector3.h"
#include "CalcEfficiency/HistNtupleMT.h"

HistNtupleMT::HistNtupleMT() {
}

HistNtupleMT::~HistNtupleMT() {
}

template void HistNtupleMT::FillHist<TagAndProbeMT>(TagAndProbeMT& tap);
template void HistNtupleMT::FillHist<TagAndProbe>(TagAndProbe& tap);
template <typename TAP> void HistNtupleMT::FillHist( TAP& tap )
{
  std::cout << "---> start filling histogram for NSW monitoring <---" << std::endl;
  int probe_n = tap.m_vL1objects_probe.size();
  if(tap.m_nmesChain < 1){
     std::cout << "No chains required to probe muon! Skip to fill NSW histogram." << std::endl;
     return;
  }
  if(tap.m_HLTtrigmesName.at(m_monChain) != m_monTrigName){
    std::cout << "choosed chain is inconsistent monitored chain! someting wrong, Skip..." << std::endl;
    return;
  }
  std::cout << "using chain: " << tap.m_HLTtrigmesName.at(m_monChain) << ", tapMethod: " << tap.m_method << std::endl;
  for( int i_probe = 0; i_probe < probe_n; i_probe++ ) {
    // NSW validation =============================
    // using HLT_mu6_L1MU6
    if( !m_isAsymNSW || 
        (tap.m_probe[i_probe].extEta > 0 && 
        tap.m_vSAobjects_probe[i_probe][m_monChain].sAddress == -1 /*endcap L2SA*/))
      {
      m_h_probept->Fill(tap.m_probe[i_probe].pt/1000.);
      m_h_probeexteta->Fill(tap.m_probe[i_probe].extEta);
      m_h_probeextphi->Fill(tap.m_probe[i_probe].extPhi);
      if(tap.m_vL1objects_probe[i_probe][m_monChain].isPassed < 1){
        std::cout << "L1 failed, continue" << std::endl;
        continue;
      }
      std::vector<float>* stgcClusterR = &tap.m_vSAobjects_probe[i_probe][m_monChain].stgcClusterR;
      std::vector<float>* stgcClusterZ = &tap.m_vSAobjects_probe[i_probe][m_monChain].stgcClusterZ;
      std::vector<float>* mmClusterR = &tap.m_vSAobjects_probe[i_probe][m_monChain].mmClusterR;
      std::vector<float>* mmClusterZ = &tap.m_vSAobjects_probe[i_probe][m_monChain].mmClusterZ;
      float innSPR = (tap.m_vSAobjects_probe[i_probe][m_monChain].superPointR_EI); //[mm]
      float innSPZ = (tap.m_vSAobjects_probe[i_probe][m_monChain].superPointZ_EI); //[mm]
      float innSPSlope = (tap.m_vSAobjects_probe[i_probe][m_monChain].superPointSlope_EI);
      //float innSPOffset = (tap.m_vSAobjects_probe[i_probe][m_monChain].superPointIntercept_EI); currently unused
      //road index [0, 1, 2] = [inner, middle, outer]
      float innroadAw = tap.m_vSAobjects_probe[i_probe][m_monChain].roadAw.at(0);
      float innroadBw = tap.m_vSAobjects_probe[i_probe][m_monChain].roadBw.at(0);
      std::vector<TVector3> offsegInner;
      std::vector<TVector3> offsegInnerDir;
      for(int i_seg = 0; i_seg < tap.m_probe[i_probe].segmentN; i_seg++){
        TVector3 segment;
        segment.SetXYZ(tap.m_probe[i_probe].segmentX[i_seg], tap.m_probe[i_probe].segmentY[i_seg], tap.m_probe[i_probe].segmentZ[i_seg] );
        TVector3 segmentDir;
        segmentDir.SetXYZ(tap.m_probe[i_probe].segmentPx[i_seg], tap.m_probe[i_probe].segmentPy[i_seg], tap.m_probe[i_probe].segmentPz[i_seg] );
        if(fabs(segment.Z()) > 5000. && fabs(segment.Z()) < 8500. && fabs(tap.m_probe[i_probe].extEta) > 1.05){//choose offline segment in Endcap inner
          offsegInner.push_back(segment);
          offsegInnerDir.push_back(segmentDir);
        }
      }
      for(int i_stgc = 0; i_stgc < (int)stgcClusterR->size(); i_stgc++){
        m_h_stgcClusterZR->Fill(stgcClusterZ->at(i_stgc)/1000. , stgcClusterR->at(i_stgc)/1000.);
      }
      for(int i_mm = 0; i_mm < (int)mmClusterR->size(); i_mm++){
        m_h_mmClusterZR->Fill(mmClusterZ->at(i_mm)/1000. , mmClusterR->at(i_mm)/1000.);
      }
      if(offsegInner.size() == 1){
        float aw_offseg = (fabs(offsegInnerDir.at(0).Z()) > 1e-5) ? offsegInnerDir.at(0).Perp()/offsegInnerDir.at(0).Z() : 0;
        float bw_offseg = offsegInner.at(0).Perp() - aw_offseg*offsegInner.at(0).Z(); 
        float SPvsOffresidual = calc_residual(innSPR, innSPZ, aw_offseg, bw_offseg);
        float SPvsRoadresidual = calc_residual(innSPR, innSPZ, innroadAw, innroadBw);
        float theta_SP = atan2(innSPSlope, 1);
        float theta_OFF = atan2(offsegInnerDir.at(0).Perp(), offsegInnerDir.at(0).Z());
        float SPvsOffdtheta = theta_OFF - theta_SP;
        m_h_residual_offsegvsinnSP->Fill(SPvsOffresidual);
        m_h_residual_roadvsinnSP->Fill(SPvsRoadresidual);
        m_h_dtheta_offvsinnSP->Fill(SPvsOffdtheta);
      }
    }
    // NSW validation end =========================
    double tpsumReqdRl1 = tap.m_requirements_tag[i_probe].reqdRl1 + tap.m_requirements_probe[i_probe].reqdRl1;
    if( tpsumReqdRl1 < tap.m_probe[i_probe].tpextdR ){
      m_h_probeEvents_offpt->Fill(tap.m_probe[i_probe].pt/1000.);
      
      // efficiency by pT check =========================
      for(int i_trig = 0; i_trig < tap.m_nmesChain; i_trig++){
        m_h_trigPassEvents[i_trig]->Fill("offline", 1); //offline
        if(tap.m_vL1objects_probe[i_probe][i_trig].isPassed == 1){
          m_h_trigPassEvents[i_trig]->Fill("L1", 1); //L1
          m_h_L1pass_offpt[i_trig]->Fill(tap.m_probe[i_probe].pt/1000);
          if(tap.m_vSAobjects_probe[i_probe][i_trig].isPassed == 1){
            m_h_trigPassEvents[i_trig]->Fill("SA", 1); //SA
            m_h_SApass_offpt[i_trig]->Fill(tap.m_probe[i_probe].pt/1000);
            if(tap.m_vCBobjects_probe[i_probe][i_trig].isPassed == 1){
              m_h_trigPassEvents[i_trig]->Fill("CB", 1); //CB
              m_h_CBpass_offpt[i_trig]->Fill(tap.m_probe[i_probe].pt/1000);
              if(tap.m_vEFobjects_probe[i_probe][i_trig].isPassed == 1){
                m_h_trigPassEvents[i_trig]->Fill("EF", 1); //EF
                m_h_EFpass_offpt[i_trig]->Fill(tap.m_probe[i_probe].pt/1000);
              }
            }
          }
        }
      } //m_nmesChain loop end
      // efficiency by pT check end =========================
    }
  }
}

template int HistNtupleMT::initialize<TagAndProbeMT>( std::string histname, TagAndProbeMT& tap, const bool isAsymNSWi, std::string montrigName );
template int HistNtupleMT::initialize<TagAndProbe>( std::string histname, TagAndProbe& tap, const bool isAsymNSW, std::string montrigName );
template <typename TAP>
int HistNtupleMT::initialize( std::string histname, TAP& tap, const bool isAsymNSW, std::string montrigName ) {
  m_isAsymNSW = isAsymNSW;
  m_monTrigName = montrigName;
  for(int i_trig = 0; i_trig < tap.m_nmesChain; i_trig++){ 
    if(tap.m_HLTtrigmesName.at(i_trig).find(montrigName) != std::string::npos){
      m_monChain = i_trig;
    }
  }
  TString etname = histname.c_str();
  m_FILE 	= new TFile( etname, "recreate" );
  // check valiables
  m_h_probept = new TH1D("m_h_probept", "probe p_{T};p_{T}^{offline}[GeV];events", 55, 0, 110);
  m_h_probeexteta = new TH1D("m_h_probeexteta", "probe #eta at MuonSpectrometer;#eta_{#mu} at MuonSpectrometer;events", 50, -2.5, 2.5);
  m_h_probeextphi = new TH1D("m_h_probeextphi", "probe #phi at MuonSpectrometer;#phi_{#mu} at MuonSpectrometer;events", 64, -3.2, 3.2);
  //
  //NSW validation
  m_h_stgcClusterZR = new TH2D("m_h_stgcClusterZR", "stgcCluster (Z, R) [m];Z[m];R[m];events", 80, -10, 10, 48, 0, 12);
  m_h_mmClusterZR = new TH2D("m_h_mmClusterZR", "mmCluster (Z, R) [m];Z[m];R[m];events", 80, -10, 10, 48, 0, 12);
  m_h_residual_offsegvsinnSP = new TH1D("m_h_residual_offvsinnSP", "residual (inner offline track vs innerSP);residual[mm];events", 100, -500, 500);
  m_h_residual_roadvsinnSP = new TH1D("m_h_residual_roadvsinnSP", "residual (inner Road vs innerSP);residual[mm];events", 100, -500, 500);
  m_h_dtheta_offvsinnSP = new TH1D("m_h_dtheta_offvsinnSP", "dtheta (inner offline track - innerSPslope);dtheta (rad);events", 40, -0.1, 0.1);
  //
  m_h_passedisMoreCand = new TH1D("m_h_passedisMoreCand", "N_{isMoreCandflag = true}, barrel events;;events", 2, 0, 2);
  m_h_passedisMoreCand->GetXaxis()->SetBinLabel(1, "L1 roi");
  m_h_passedisMoreCand->GetXaxis()->SetBinLabel(2, "isMoreCand true");
  m_h_passedisMoreCand->SetMinimum(0);
  //efficiency
  m_h_probeEvents_offpt = new TH1D("m_h_probeEvents_offpt", "probe events;p_{T}^{offline}[GeV];events", 100, 0, 50);
  for(int i_trig = 0; i_trig < tap.m_nmesChain; i_trig++){
    //L1 pass
    m_h_L1pass_offpt.push_back(new TH1D(Form("m_h_L1pass_%s", tap.m_L1trigmesName.at(i_trig).data()), Form("%s passed events;p_{T}^{offline}[GeV];events", tap.m_L1trigmesName.at(i_trig).data()), 100, 0, 50));
    m_eff_L1pass_offpt.push_back(new TH1D(Form("m_eff_L1pass_%s", tap.m_L1trigmesName.at(i_trig).data()), Form("%s passed efficiency;p_{T}^{offline}[GeV];#epsilon", tap.m_L1trigmesName.at(i_trig).data()), 100, 0, 50));
    //SA pass
    m_h_SApass_offpt.push_back(new TH1D(Form("m_h_SApass_%s", tap.m_HLTtrigmesName.at(i_trig).data()), Form("%s passed events;p_{T}^{offline}[GeV];events", tap.m_HLTtrigmesName.at(i_trig).data()), 100, 0, 50));
    m_eff_SApass_offpt.push_back(new TH1D(Form("m_eff_SApass_%s", tap.m_HLTtrigmesName.at(i_trig).data()), Form("%s passed efficiency SA/L1;p_{T}^{offline}[GeV];#epsilon", tap.m_HLTtrigmesName.at(i_trig).data()), 100, 0, 50));
    //CB pass
    m_h_CBpass_offpt.push_back(new TH1D(Form("m_h_CBpass_%s", tap.m_HLTtrigmesName.at(i_trig).data()), Form("%s passed events;p_{T}^{offline}[GeV];events", tap.m_HLTtrigmesName.at(i_trig).data()), 100, 0, 50));
    m_eff_CBpass_offpt.push_back(new TH1D(Form("m_eff_CBpass_%s", tap.m_HLTtrigmesName.at(i_trig).data()), Form("%s passed efficiency CB/SA;p_{T}^{offline}[GeV];#epsilon", tap.m_HLTtrigmesName.at(i_trig).data()), 100, 0, 50));
    //EF pass
    m_h_EFpass_offpt.push_back(new TH1D(Form("m_h_EFpass_%s", tap.m_HLTtrigmesName.at(i_trig).data()), Form("%s passed events;p_{T}^{offline}[GeV];events", tap.m_HLTtrigmesName.at(i_trig).data()), 100, 0, 50));
    m_eff_EFpass_offpt.push_back(new TH1D(Form("m_eff_EFpass_%s", tap.m_HLTtrigmesName.at(i_trig).data()), Form("%s passed efficiency EF/CB;p_{T}^{offline}[GeV];#epsilon", tap.m_HLTtrigmesName.at(i_trig).data()), 100, 0, 50));
    // each trigger step passed events
    m_h_trigPassEvents.push_back(new TH1D(Form("m_h_trigPassEvents_%s", tap.m_HLTtrigmesName.at(i_trig).data()), Form("%s passed events;trigger step;events", tap.m_HLTtrigmesName.at(i_trig).data()), 5, 0, 5));
    m_h_trigPassEvents.at(i_trig)->GetXaxis()->SetBinLabel(1, "offline");
    m_h_trigPassEvents.at(i_trig)->GetXaxis()->SetBinLabel(2, "L1");
    m_h_trigPassEvents.at(i_trig)->GetXaxis()->SetBinLabel(3, "SA");
    m_h_trigPassEvents.at(i_trig)->GetXaxis()->SetBinLabel(4, "CB");
    m_h_trigPassEvents.at(i_trig)->GetXaxis()->SetBinLabel(5, "EF");
  }
  return 0;
}

template int HistNtupleMT::finalize<TagAndProbeMT>(TagAndProbeMT& tap);
template int HistNtupleMT::finalize<TagAndProbe>(TagAndProbe& tap);
template <typename TAP> int HistNtupleMT::finalize( TAP& tap )
{
  for(int i_trig = 0; i_trig < tap.m_nmesChain; i_trig++){
    CalcEfficiency(m_h_L1pass_offpt.at(i_trig), m_h_probeEvents_offpt, m_eff_L1pass_offpt.at(i_trig));
    CalcEfficiency(m_h_SApass_offpt.at(i_trig), m_h_L1pass_offpt.at(i_trig), m_eff_SApass_offpt.at(i_trig));
    CalcEfficiency(m_h_CBpass_offpt.at(i_trig), m_h_SApass_offpt.at(i_trig), m_eff_CBpass_offpt.at(i_trig));
    CalcEfficiency(m_h_EFpass_offpt.at(i_trig), m_h_CBpass_offpt.at(i_trig), m_eff_EFpass_offpt.at(i_trig));
  }
  m_FILE->Write();
  delete m_h_passedisMoreCand;
  delete m_h_probeEvents_offpt;
  for(int i_trig = 0; i_trig < tap.m_nmesChain; i_trig++){
    delete m_h_L1pass_offpt[i_trig];
    delete m_eff_L1pass_offpt[i_trig];
    delete m_h_SApass_offpt[i_trig];
    delete m_eff_SApass_offpt[i_trig];
    delete m_h_CBpass_offpt[i_trig];
    delete m_eff_CBpass_offpt[i_trig];
    delete m_h_EFpass_offpt[i_trig];
    delete m_eff_EFpass_offpt[i_trig];
    delete m_h_trigPassEvents[i_trig];
  }

  return 0;
}

void HistNtupleMT::CalcEfficiency(TH1D* h_num, TH1D* h_den, TH1D* h_set){
  double eff_y, eff_yerr;
  unsigned int nbins = h_den->GetXaxis()->GetNbins();
  for(unsigned int ibin = 1; ibin <= nbins; ibin++){
    double denominator = h_den->GetBinContent(ibin);
    double numerator = h_num->GetBinContent(ibin);
    eff_y = (denominator != 0) ? numerator/denominator : 0;
    eff_yerr = (denominator != 0) ? sqrt((numerator*(1.-2*eff_y)+(denominator*pow(eff_y, 2)))/pow(denominator, 2)):0;
    h_set->SetBinContent(ibin, eff_y);
    h_set->SetBinError(ibin, eff_yerr);
  }//ibin loop end
}

inline float HistNtupleMT::calc_residual(float x, float y, float aw, float bw)
{
  const float ZERO_LIMIT = 1e-4;
  if( fabs(aw) < ZERO_LIMIT ) return y-bw;
  float ia  = 1/aw;
  float iaq = ia*ia;
  float dz  = x - (y-bw)*ia;
  return dz/sqrt(1.+iaq);
}

