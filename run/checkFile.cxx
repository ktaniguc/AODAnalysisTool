#include "TH1D.h"
#include <iostream>
#include <TH1.h>
#include <TFile.h>
#include "/home/ktaniguc/RootUtils/RootUtils/TCanvas_opt.h"
#include "/home/ktaniguc/RootUtils/RootUtils/TLegend_addfunc.h"
#include "/home/ktaniguc/RootUtils/src/rootlogon.C"
#include "TEfficiency.h"

using namespace std;

void checkFile(){
  gStyle->SetOptStat(0);
  rootlogon(); 
  TFile *_file0 = TFile::Open("./wofix/efficiency-monitoring.root");
  TCanvas_opt *c1 = new TCanvas_opt();
  string histname[1];
  //histname[0] = "m_h_SARoadtheta";
  //histname[0] = "m_h_SARoadoffset";
  //histname[0] = "m_h_dtheta_offvsroad";
  //histname[0] = "h_Offpt";
//  histname[0] = "m_h_passedisMoreCand;1";
  histname[0] = "m_h_trigPassEvents_HLT_mu26_ivarmedium_L1MU20;1";
  TH1D *h0 = (TH1D*)_file0->Get(histname[0].c_str());
  c1->cd();
  h0->Draw("hist text");
  TLegend_addfunc *legptres = new TLegend_addfunc(2);
  legptres->AddSelection("Z->#mu#mu");
  legptres->AddSelection("roiNum only");
  legptres->Draw("same");
}


