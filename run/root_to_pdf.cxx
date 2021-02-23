
void root_to_pdf(){
  rootlogon();
  TCanvas *c1 = new TCanvas("c1", "c1", 10, 10, 1020, 700);
  c1->SetGrid();
  c1->SetRightMargin(0.20);
  c1->SetLeftMargin(0.23);
  c1->SetBottomMargin(0.20);

  gStyle->SetLegendBorderSize(0);


  //TString label = "./plot/plot_tmp_short.pdf";
  //TFile *file0 = TFile::Open("./outroot/short_bsub_reco12_max_barrel.root");

  //TString label = "./plot/plot_tmp.pdf";
  //TFile *file0 = TFile::Open("short_bsub_reco12_max.root");

  //TString label = "./plot/plot_tmp_barrel.pdf";
  //TFile *file0 = TFile::Open("./outroot/bsub_reco12_max_barrel.root");

  //TString label = "./plot/plot_tmp_barrel1.pdf";
  //TFile *file0 = TFile::Open("./outroot/bsub_reco12_max_barrel1.root");

  //TString label = "./plot/plot_tmp_barrel_pt.pdf";
  //TFile *file0 = TFile::Open("./outroot/bsub_reco12_max_barrel_pt.root");

  //TString label = "./plot/plot_tmp_max.pdf";
  //TFile *file0 = TFile::Open("./outroot/bsub_reco12_max.root");

  //TString label = "short_bsub_reco12_max_barrel1.pdf";
  //TFile *file0 = TFile::Open("./outroot/short_bsub_reco12_max_barrel1.root");

  //TString label = "./plot/plot_bsub_reco12_max_barrel2.pdf";
  //TFile *file0 = TFile::Open("./outroot/bsub_reco12_max_barrel2.root");

  //TString label = "./plot/plot_short_bsub_reco12_max_barrel2.pdf";
  //TFile *file0 = TFile::Open("./outroot/short_bsub_reco12_max_barrel2.root");

  //TString label = "./plot/plot_short_bsub_reco12_max_barrel3.pdf";
  //TFile *file0 = TFile::Open("./outroot/short_bsub_reco12_max_barrel3.root");

  //TString label = "./plot/plot_short_bsub_reco12_max_barrel4.pdf";
  //TFile *file0 = TFile::Open("./outroot/short_bsub_reco12_max_barrel4.root");

  //TString label = "./plot/plot_short_bsub_reco12_max_barrel5.pdf";
  //TFile *file0 = TFile::Open("./outroot/short_bsub_reco12_max_barrel5.root");

  //TString label = "./plot/plot_bsub_reco12_max_barrel3.pdf";
  //TFile *file0 = TFile::Open("./outroot/bsub_reco12_max_barrel3.root");

  //TString label = "./plot/plot_short_bsub_reco12_max_barrel7.pdf";
  //TFile *file0 = TFile::Open("./outroot/short_bsub_reco12_max_barrel7.root");

  //TString label = "./plot/plot_short_bsub_reco12_max_barrel8.pdf";
  //TFile *file0 = TFile::Open("./outroot/short_bsub_reco12_max_barrel8.root");

  //TString label = "./plot/plot_short_bsub_reco12_max_barrel9.pdf";
  //TFile *file0 = TFile::Open("./outroot/short_bsub_reco12_max_barrel9.root");

  //TString label = "./plot/plot_short_bsub_reco12_max_barrel10.pdf";
  //TFile *file0 = TFile::Open("./outroot/short_bsub_reco12_max_barrel10.root");

  //TString label = "./plot/plot_short_bsub_reco12_max_barrel11.pdf";
  //TFile *file0 = TFile::Open("./outroot/short_bsub_reco12_max_barrel11.root");

  //TString label = "./plot/plot_short_bsub_reco12_max_barrel12.pdf";
  //TFile *file0 = TFile::Open("./outroot/short_bsub_reco12_max_barrel12.root");

  //TString label = "./plot/plot_bsub_reco12_max_barrel5.pdf";
  //TFile *file0 = TFile::Open("./outroot/bsub_reco12_max_barrel5.root");

  //TString label = "./plot/plot_bsub_reco14_max_barrel.pdf";
  //TFile *file0 = TFile::Open("./outroot/bsub_reco14_max_barrel.root");

  //TString label = "./plot/plot_short_bsub_reco14_max_barrel.pdf";
  //TFile *file0 = TFile::Open("./outroot/short_bsub_reco14_max_barrel.root");

  //TString label = "./plot/plot_short_bsub_reco14_max_barrel1.pdf";
  //TFile *file0 = TFile::Open("./outroot/short_bsub_reco14_max_barrel1.root");

  //t_349014.DrawHist("../plot/DrawHist_" + PdfLabel + ".pdf");
  //TFile *fout = new TFile(Form("outroot/%s.root", PdfLabel.Data()), "recreate");

  //TString label = "./plot/plot_bsub_reco14_max_barrel1.pdf";
  //TFile *file0 = TFile::Open("./outroot/bsub_reco14_max_barrel1.root");

  //TString label = "./plot/plot_short_bsub_reco12_max_barrel13.pdf";
  //TFile *file0 = TFile::Open("./outroot/short_bsub_reco12_max_barrel13.root");

  //TString label = "./plot/plot_short_bsub_reco12_max_barrel14.pdf";
  //TFile *file0 = TFile::Open("./outroot/short_bsub_reco12_max_barrel14.root");

  //TString label = "./plot/plot_bsub_reco14_max_barrel_new.pdf";
  //TFile *file0 = TFile::Open("./outroot/bsub_reco14_max_barrel_new.root");

  //TString label = "./plot/plot_bsub_reco12_max_barrel_new.pdf";
  //TFile *file0 = TFile::Open("./outroot/bsub_reco12_max_barrel_new.root");

  //TString label = "./plot/plot_bsub_reco14_max_barrel_new1.pdf";
  //TFile *file0 = TFile::Open("./outroot/bsub_reco14_max_barrel_new1.root");

  //TString label = "./plot/plot_short_bsub_reco14_max_barrel2.pdf";
  //TFile *file0 = TFile::Open("./outroot/short_bsub_reco14_max_barrel2.root");

  //TString label = "./plot/plot_bsub_reco14_max_barrel_new2.pdf";
  //TFile *file0 = TFile::Open("./outroot/bsub_reco14_max_barrel_new2.root");

  //TString label = "./plot/plot_bsub_reco16_max_barrel_new.pdf";
  //TFile *file0 = TFile::Open("./outroot/bsub_reco16_max_barrel_new.root");

  //TString label = "./plot/plot_bsub_reco16_max_barrel_new1.pdf";
  //TFile *file0 = TFile::Open("./outroot/bsub_reco16_max_barrel_new1.root");

  //TString label = "./plot/plot_bsub_reco14_max_barrel_new3.pdf";
  //TFile *file0 = TFile::Open("./outroot/bsub_reco14_max_barrel_new3.root");

  //TString label = "./plot/plot_bsub_reco14_max_barrel_new6.pdf";
  //TFile *file0 = TFile::Open("./outroot/bsub_reco14_max_barrel_new6.root");

  //TString label = "./plot/plot_bsub_reco14_max_barrel_new7.pdf";
  //TFile *file0 = TFile::Open("./outroot/bsub_reco14_max_barrel_new7.root");

  //TString label = "./plot/plot_bsub_reco14_max_barrel_new8.pdf";
  //TFile *file0 = TFile::Open("./outroot/bsub_reco14_max_barrel_new8.root");

  //TString label = "./plot/plot_bsub_reco14_max_barrel_new9.pdf";
  //TFile *file0 = TFile::Open("./outroot/bsub_reco14_max_barrel_new9.root");

  //TString label = "./plot/tmp_reco18_short.pdf";
  //TFile *file0 = TFile::Open("./outroot/tmp_reco18_short.root");

  //TString label = "./plot/plot_bsub_reco18_max_barrel_new.pdf";
  //TFile *file0 = TFile::Open("./outroot/bsub_reco18_max_barrel_new.root");

  //TString label = "./plot/plot_bsub_reco18_max_barrel_new1.pdf";
  //TFile *file0 = TFile::Open("./outroot/bsub_reco18_max_barrel_new1.root");

  //TString label = "./plot/plot_bsub_reco18_max_barrel_new2.pdf";
  //TFile *file0 = TFile::Open("./outroot/bsub_reco18_max_barrel_new2.root");

  //TString label = "./plot/plot_bsub_reco18_max_barrel_new3.pdf";
  //TFile *file0 = TFile::Open("./outroot/bsub_reco18_max_barrel_new3.root");

  //TString label = "./plot/plot_bsub_reco18_max_barrel_new4.pdf";
  //TFile *file0 = TFile::Open("./outroot/bsub_reco18_max_barrel_new4.root");

  //TString label = "./plot/plot_bsub_reco19_max_barrel_new.pdf";
  //TFile *file0 = TFile::Open("./outroot/bsub_reco19_max_barrel_new.root");

  //TString label = "./plot/plot_bsub_reco19_max_barrel_new3.pdf";
  //TFile *file0 = TFile::Open("./outroot/bsub_reco19_max_barrel_new3.root");

  //TString label = "./plot/plot_bsub_reco19_max_barrel_new4.pdf";
  //TFile *file0 = TFile::Open("./outroot/bsub_reco19_max_barrel_new.root");

  //TString label = "./plot/plot_bsub_reco19_max_barrel_new5.pdf";
  //TFile *file0 = TFile::Open("./outroot/bsub_reco19_max_barrel_new5.root");

  //TString label = "./plot/plot_bsub_reco19_max_barrel_new6.pdf";
  //TFile *file0 = TFile::Open("./outroot/bsub_reco19_max_barrel_new6.root");

  //TString label = "./plot/plot_bsub_reco19_max_barrel_new7.pdf";
  //TFile *file0 = TFile::Open("./outroot/bsub_reco19_max_barrel_new7.root");

  TString label = "./plot/plot_bsub_reco19_max_barrel_new8.pdf";
  TFile *file0 = TFile::Open("./outroot/bsub_reco19_max_barrel_new8.root");

  TH1F* h_pt_L1 = (TH1F*)file0->Get("m_h_pt_muL1");
  TH1F* h_pt_SA = (TH1F*)file0->Get("m_h_pt_muSA");
  TH1F* h_pt_CB = (TH1F*)file0->Get("m_h_pt_muCB");
  TH1F* h_pt_SAIO = (TH1F*)file0->Get("m_h_pt_muSAIO");
  TH1F* h_pt_CBIO = (TH1F*)file0->Get("m_h_pt_muCBIO");

  TH1F* h_pt_L1_1 = (TH1F*)file0->Get("m_h_pt_muL1_1");
  TH1F* h_pt_SA_1 = (TH1F*)file0->Get("m_h_pt_muSA_1");
  TH1F* h_pt_CB_1 = (TH1F*)file0->Get("m_h_pt_muCB_1");
  TH1F* h_pt_SAIO_1 = (TH1F*)file0->Get("m_h_pt_muSAIO_1");
  TH1F* h_pt_CBIO_1 = (TH1F*)file0->Get("m_h_pt_muCBIO_1");

  TH1F* h_L1_sep = (TH1F*)file0->Get("m_h_dR_2muL1_sep");
  TH1F* h_SA_sep = (TH1F*)file0->Get("m_h_dR_2muSA_sep");
  TH1F* h_CB_sep = (TH1F*)file0->Get("m_h_dR_2muCB_sep");
  TH1F* h_SAIO_sep = (TH1F*)file0->Get("m_h_dR_2muSAIO_sep");
  TH1F* h_CBIO_sep = (TH1F*)file0->Get("m_h_dR_2muCBIO_sep");

  TH1F* h_L1_1_sep = (TH1F*)file0->Get("m_h_dR_2muL1_sep1");
  TH1F* h_SA_1_sep = (TH1F*)file0->Get("m_h_dR_2muSA_sep1");
  TH1F* h_CB_1_sep = (TH1F*)file0->Get("m_h_dR_2muCB_sep1");
  TH1F* h_SAIO_1_sep = (TH1F*)file0->Get("m_h_dR_2muSAIO_sep1");
  TH1F* h_CBIO_1_sep = (TH1F*)file0->Get("m_h_dR_2muCBIO_sep1");

  TH1F* h_L1_2_sep = (TH1F*)file0->Get("m_h_dR_2muL1_sep2");
  TH1F* h_SA_2_sep = (TH1F*)file0->Get("m_h_dR_2muSA_sep2");
  TH1F* h_CB_2_sep = (TH1F*)file0->Get("m_h_dR_2muCB_sep2");
  TH1F* h_SAIO_2_sep = (TH1F*)file0->Get("m_h_dR_2muSAIO_sep2");
  TH1F* h_CBIO_2_sep = (TH1F*)file0->Get("m_h_dR_2muCBIO_sep2");

  TH1F* h_L1_3_sep = (TH1F*)file0->Get("m_h_dR_2muL1_sep3");
  TH1F* h_SA_3_sep = (TH1F*)file0->Get("m_h_dR_2muSA_sep3");
  TH1F* h_CB_3_sep = (TH1F*)file0->Get("m_h_dR_2muCB_sep3");
  TH1F* h_SAIO_3_sep = (TH1F*)file0->Get("m_h_dR_2muSAIO_sep3");
  TH1F* h_CBIO_3_sep = (TH1F*)file0->Get("m_h_dR_2muCBIO_sep3");

  TH1F* h_L1_4_sep = (TH1F*)file0->Get("m_h_dR_2muL1_sep4");
  TH1F* h_SA_4_sep = (TH1F*)file0->Get("m_h_dR_2muSA_sep4");
  TH1F* h_CB_4_sep = (TH1F*)file0->Get("m_h_dR_2muCB_sep4");
  TH1F* h_SAIO_4_sep = (TH1F*)file0->Get("m_h_dR_2muSAIO_sep4");
  TH1F* h_CBIO_4_sep = (TH1F*)file0->Get("m_h_dR_2muCBIO_sep4");

  TH1F* h_L1_5_sep = (TH1F*)file0->Get("m_h_dR_2muL1_sep5");
  TH1F* h_SA_5_sep = (TH1F*)file0->Get("m_h_dR_2muSA_sep5");
  TH1F* h_CB_5_sep = (TH1F*)file0->Get("m_h_dR_2muCB_sep5");
  TH1F* h_SAIO_5_sep = (TH1F*)file0->Get("m_h_dR_2muSAIO_sep5");
  TH1F* h_CBIO_5_sep = (TH1F*)file0->Get("m_h_dR_2muCBIO_sep5");

  TH1F* h_L1_6_sep = (TH1F*)file0->Get("m_h_dR_2muL1_sep6");
  TH1F* h_SA_6_sep = (TH1F*)file0->Get("m_h_dR_2muSA_sep6");
  TH1F* h_CB_6_sep = (TH1F*)file0->Get("m_h_dR_2muCB_sep6");
  TH1F* h_SAIO_6_sep = (TH1F*)file0->Get("m_h_dR_2muSAIO_sep6");
  TH1F* h_CBIO_6_sep = (TH1F*)file0->Get("m_h_dR_2muCBIO_sep6");

  TH2F* h_SA_pt_res = (TH2F*)file0->Get("m_h_pt_res_muSA");
  TH2F* h_CB_pt_res = (TH2F*)file0->Get("m_h_pt_res_muCB");
  TH2F* h_SAIO_pt_res = (TH2F*)file0->Get("m_h_pt_res_muSAIO");
  TH2F* h_CBIO_pt_res = (TH2F*)file0->Get("m_h_pt_res_muCBIO");

  TH2F* h_SA_pt_res_1 = (TH2F*)file0->Get("m_h_pt_res_muSA_1");
  TH2F* h_CB_pt_res_1 = (TH2F*)file0->Get("m_h_pt_res_muCB_1");
  TH2F* h_SAIO_pt_res_1 = (TH2F*)file0->Get("m_h_pt_res_muSAIO_1");
  TH2F* h_CBIO_pt_res_1 = (TH2F*)file0->Get("m_h_pt_res_muCBIO_1");

  TH2F* h_SA_pt_res_2 = (TH2F*)file0->Get("m_h_pt_res_muSA_2");
  TH2F* h_CB_pt_res_2 = (TH2F*)file0->Get("m_h_pt_res_muCB_2");
  TH2F* h_SAIO_pt_res_2 = (TH2F*)file0->Get("m_h_pt_res_muSAIO_2");
  TH2F* h_CBIO_pt_res_2 = (TH2F*)file0->Get("m_h_pt_res_muCBIO_2");

  TH2F* h_SA_pt_res_3 = (TH2F*)file0->Get("m_h_pt_res_muSA_3");
  TH2F* h_CB_pt_res_3 = (TH2F*)file0->Get("m_h_pt_res_muCB_3");
  TH2F* h_SAIO_pt_res_3 = (TH2F*)file0->Get("m_h_pt_res_muSAIO_3");
  TH2F* h_CBIO_pt_res_3 = (TH2F*)file0->Get("m_h_pt_res_muCBIO_3");

  TH2F* m_h_pt_numberOfSPs_muSA = (TH2F*)file0->Get("m_h_pt_numberOfSPs_muSA");
  //TH2F* h_SA_pt_res = (TH2F*)file0->Get("m_h_numberOfSPs_muCB");
  TH2F* m_h_pt_numberOfSPs_muSAIO = (TH2F*)file0->Get("m_h_pt_numberOfSPs_muSAIO");
  TH2F* m_h_pt_numberOfSPs_muSAIO_1 = (TH2F*)file0->Get("m_h_pt_numberOfSPs_muSAIO_1");
  //TH2F* h_SA_pt_res = (TH2F*)file0->Get("m_h_numberOfSPs_muCBIO");

  TH2F* m_h_pt_numberOfFTFs = (TH2F*)file0->Get("m_h_pt_numberOfFTFs");
  TH2F* m_h_pt_numberOfFTFs1 = (TH2F*)file0->Get("m_h_pt_numberOfFTFs1");
  TH2F* m_h_pt_numberOfFTFs2 = (TH2F*)file0->Get("m_h_pt_numberOfFTFs2");

  std::cout << "debug" << std::endl;
  c1 -> Print(label + "[", "pdf");

  //gStyle->SetPalette(kSolar);
  gStyle->SetPalette(kNeon);
  TColor::InvertPalette();
  //TH2F* m_h_pt_alpha1 = (TH2F*)file0->Get("m_h_pt_alpha1");
  //m_h_pt_alpha1->GetYaxis()->SetRangeUser(0,0.1);
  //m_h_pt_alpha1->GetXaxis()->SetTitle("1/p_{T}[1/GeV]");
  //m_h_pt_alpha1->GetYaxis()->SetTitle("#alpha");
  //m_h_pt_alpha1->GetZaxis()->SetTitle("Counts");
  //m_h_pt_alpha1->Draw("colz");
  //c1 -> Print(label, "pdf");


  //TH2F* m_h_pt_beta1 = (TH2F*)file0->Get("m_h_pt_beta1");
  //m_h_pt_beta1->GetYaxis()->SetRangeUser(0,0.1);
  //m_h_pt_beta1->GetXaxis()->SetTitle("1/p_{T}[1/GeV]");
  //m_h_pt_beta1->GetYaxis()->SetTitle("#beta");
  //m_h_pt_beta1->GetZaxis()->SetTitle("Counts");
  //m_h_pt_beta1->Draw("colz");
  //c1 -> Print(label, "pdf");


  TH2F* m_h_di_extIP_deltaR = (TH2F*)file0->Get("m_h_di_extIP_deltaR");
  m_h_di_extIP_deltaR->GetYaxis()->SetRangeUser(0,0.4);
  m_h_di_extIP_deltaR->GetXaxis()->SetRangeUser(0,0.4);
  m_h_di_extIP_deltaR->Draw("colz");
  m_h_di_extIP_deltaR->Draw("colz");
  TColor::InvertPalette();
  c1 -> Print(label, "pdf");


  h_CBIO_pt_res->SetLineColor(kGreen+2);
  h_CBIO_pt_res->ProjectionY()->GetXaxis()->SetRangeUser(-1,1);
  h_CBIO_pt_res->ProjectionY()->Draw();
  h_SA_pt_res->SetLineColor(kBlack);
  h_SA_pt_res->ProjectionY()->Draw("same");
  h_CB_pt_res->SetLineColor(kRed);
  h_CB_pt_res->ProjectionY()->Draw("same");
  h_SAIO_pt_res->SetLineColor(kBlue);
  h_SAIO_pt_res->ProjectionY()->Draw("same");
  c1 -> Print(label, "pdf");

  h_CBIO_pt_res->ProjectionY()->DrawNormalized("",1./(h_CBIO_pt_res->ProjectionY()->GetBinWidth(1)));
  h_SA_pt_res->ProjectionY()->DrawNormalized("same",1./(h_SA_pt_res->ProjectionY()->GetBinWidth(1)));
  h_CB_pt_res->ProjectionY()->DrawNormalized("same",1./(h_CB_pt_res->ProjectionY()->GetBinWidth(1)));
  h_SAIO_pt_res->ProjectionY()->DrawNormalized("same",1./(h_SAIO_pt_res->ProjectionY()->GetBinWidth(1)));
  c1 -> Print(label, "pdf");

  h_CBIO_pt_res_1->SetLineColor(kGreen+2);
  h_CBIO_pt_res_1->ProjectionY()->Draw();
  h_SA_pt_res_1->SetLineColor(kBlack);
  h_SA_pt_res_1->ProjectionY()->Draw("same");
  h_CB_pt_res_1->SetLineColor(kRed);
  h_CB_pt_res_1->ProjectionY()->Draw("same");
  h_SAIO_pt_res_1->SetLineColor(kBlue);
  h_SAIO_pt_res_1->ProjectionY()->Draw("same");
  c1 -> Print(label, "pdf");

  h_CBIO_pt_res_1->ProjectionY()->DrawNormalized("",1./(h_CBIO_pt_res_1->ProjectionY()->GetBinWidth(1)));
  h_SA_pt_res_1->ProjectionY()->DrawNormalized("same",1./(h_SA_pt_res_1->ProjectionY()->GetBinWidth(1)));
  h_CB_pt_res_1->ProjectionY()->DrawNormalized("same",1./(h_CB_pt_res_1->ProjectionY()->GetBinWidth(1)));
  h_SAIO_pt_res_1->ProjectionY()->DrawNormalized("same",1./(h_SAIO_pt_res_1->ProjectionY()->GetBinWidth(1)));
  c1 -> Print(label, "pdf");

  h_CBIO_pt_res_2->SetLineColor(kGreen+2);
  h_CBIO_pt_res_2->ProjectionY()->Draw();
  h_SA_pt_res_2->SetLineColor(kBlack);
  h_SA_pt_res_2->ProjectionY()->Draw("same");
  h_CB_pt_res_2->SetLineColor(kRed);
  h_CB_pt_res_2->ProjectionY()->Draw("same");
  h_SAIO_pt_res_2->SetLineColor(kBlue);
  h_SAIO_pt_res_2->ProjectionY()->Draw("same");
  c1 -> Print(label, "pdf");

  h_CBIO_pt_res_2->ProjectionY()->DrawNormalized("",1./(h_CBIO_pt_res_2->ProjectionY()->GetBinWidth(1)));
  h_SA_pt_res_2->ProjectionY()->DrawNormalized("same",1./(h_SA_pt_res_2->ProjectionY()->GetBinWidth(1)));
  h_CB_pt_res_2->ProjectionY()->DrawNormalized("same",1./(h_CB_pt_res_2->ProjectionY()->GetBinWidth(1)));
  h_SAIO_pt_res_2->ProjectionY()->DrawNormalized("same",1./(h_SAIO_pt_res_2->ProjectionY()->GetBinWidth(1)));
  c1 -> Print(label, "pdf");

  h_CBIO_pt_res_3->GetYaxis()->SetRangeUser(-1,1);
  h_CBIO_pt_res_3->SetLineColor(kBlue);
  h_CBIO_pt_res_3->ProjectionY()->Draw();
  h_SA_pt_res_3->SetLineColor(kRed);
  h_SA_pt_res_3->ProjectionY()->Draw("same");
  //h_CB_pt_res_3->SetLineColor(kRed);
  //h_CB_pt_res_3->ProjectionY()->Draw("same");
  //h_SAIO_pt_res_3->SetLineColor(kBlue);
  //h_SAIO_pt_res_3->ProjectionY()->Draw("same");
  c1 -> Print(label, "pdf");

  h_CBIO_pt_res_3->ProjectionY()->DrawNormalized("",1./(h_CBIO_pt_res_3->ProjectionY()->GetBinWidth(1)));
  h_SA_pt_res_3->ProjectionY()->DrawNormalized("same",1./(h_SA_pt_res_3->ProjectionY()->GetBinWidth(1)));
  //h_CB_pt_res_3->ProjectionY()->DrawNormalized("same",1./(h_CB_pt_res_3->ProjectionY()->GetBinWidth(1)));
  //h_SAIO_pt_res_3->ProjectionY()->DrawNormalized("same",1./(h_SAIO_pt_res_3->ProjectionY()->GetBinWidth(1)));
  c1 -> Print(label, "pdf");

  // ================fraction of #FTFs==============
  TH1D* nFTFs_all = m_h_pt_numberOfFTFs->ProjectionX("nFTFs_all");
  TH1D* nFTFs_0   = m_h_pt_numberOfFTFs->ProjectionX("nFTFs_0", m_h_pt_numberOfFTFs->GetYaxis()->FindBin(-2),  m_h_pt_numberOfFTFs->GetYaxis()->FindBin(0.5));
  TH1D* nFTFs_1   = m_h_pt_numberOfFTFs->ProjectionX("nFTFs_1", m_h_pt_numberOfFTFs->GetYaxis()->FindBin(0.6), m_h_pt_numberOfFTFs->GetYaxis()->FindBin(2.5));
  TH1D* nFTFs_2   = m_h_pt_numberOfFTFs->ProjectionX("nFTFs_2", m_h_pt_numberOfFTFs->GetYaxis()->FindBin(2.6), m_h_pt_numberOfFTFs->GetYaxis()->FindBin(4.5));
  TH1D* nFTFs_3   = m_h_pt_numberOfFTFs->ProjectionX("nFTFs_3", m_h_pt_numberOfFTFs->GetYaxis()->FindBin(4.6), m_h_pt_numberOfFTFs->GetYaxis()->FindBin(6.5));
  TH1D* nFTFs_4   = m_h_pt_numberOfFTFs->ProjectionX("nFTFs_4", m_h_pt_numberOfFTFs->GetYaxis()->FindBin(6.6), m_h_pt_numberOfFTFs->GetYaxis()->FindBin(11));

  //nSPs_all -> Draw();
  //nSPs_0 -> Draw("same");
  //nSPs_1 -> Draw("same");
  //nSPs_2 -> Draw("same");
  //nSPs_3 -> Draw("same");
  //nSPs_4 -> Draw("same");
  //c1 -> Print( pdf, "pdf" );

  // divided by pass-hist and not-pass-hist
  nFTFs_0->Divide(nFTFs_all);
  nFTFs_1->Divide(nFTFs_all);
  nFTFs_2->Divide(nFTFs_all);
  nFTFs_3->Divide(nFTFs_all);
  nFTFs_4->Divide(nFTFs_all);

  THStack *nFTFs_pt = new THStack("hs_pt",Form(";%s;Fraction of #FTFs (p_{T} > 3GeV)",m_h_pt_numberOfFTFs->GetXaxis()->GetTitle()));
  nFTFs_0->SetFillColor(kYellow);
  nFTFs_1->SetFillColor(kGreen+2);
  nFTFs_2->SetFillColor(kBlue);
  nFTFs_3->SetFillColor(kRed);
  nFTFs_4->SetFillColor(kBlack);

  nFTFs_pt->Add(nFTFs_4);
  nFTFs_pt->Add(nFTFs_3);
  nFTFs_pt->Add(nFTFs_2);
  nFTFs_pt->Add(nFTFs_1);
  nFTFs_pt->Add(nFTFs_0);

  gStyle->SetHistTopMargin(0);
  nFTFs_pt->Draw();
  nFTFs_pt->SetMaximum(1.1);
  nFTFs_pt->Draw();

  TLegend *leg_nFTFs_pt = new TLegend(0.02,0.00,0.1,0.4);
  leg_nFTFs_pt->AddEntry(nFTFs_0,"0","f");
  leg_nFTFs_pt->AddEntry(nFTFs_1,"1,2","f");
  leg_nFTFs_pt->AddEntry(nFTFs_2,"3,4","f");
  leg_nFTFs_pt->AddEntry(nFTFs_3,"5,6","f");
  leg_nFTFs_pt->AddEntry(nFTFs_4,"7~","f");
  leg_nFTFs_pt->Draw();

  c1 -> Print( label, "pdf" );


  // ================fraction of high pT #FTFs==============
  TH1D* nFTFs1_all = m_h_pt_numberOfFTFs1->ProjectionX("nFTFs_all");
  TH1D* nFTFs1_0   = m_h_pt_numberOfFTFs1->ProjectionX("nFTFs_0", m_h_pt_numberOfFTFs1->GetYaxis()->FindBin(-1),  m_h_pt_numberOfFTFs1->GetYaxis()->FindBin(0.5));
  TH1D* nFTFs1_1   = m_h_pt_numberOfFTFs1->ProjectionX("nFTFs_1", m_h_pt_numberOfFTFs1->GetYaxis()->FindBin(0.6), m_h_pt_numberOfFTFs1->GetYaxis()->FindBin(2.5));
  TH1D* nFTFs1_2   = m_h_pt_numberOfFTFs1->ProjectionX("nFTFs_2", m_h_pt_numberOfFTFs1->GetYaxis()->FindBin(2.6), m_h_pt_numberOfFTFs1->GetYaxis()->FindBin(4.5));
  TH1D* nFTFs1_3   = m_h_pt_numberOfFTFs1->ProjectionX("nFTFs_3", m_h_pt_numberOfFTFs1->GetYaxis()->FindBin(4.6), m_h_pt_numberOfFTFs1->GetYaxis()->FindBin(6.5));
  TH1D* nFTFs1_4   = m_h_pt_numberOfFTFs1->ProjectionX("nFTFs_4", m_h_pt_numberOfFTFs1->GetYaxis()->FindBin(6.6), m_h_pt_numberOfFTFs1->GetYaxis()->FindBin(11));

  //nSPs_all -> Draw();
  //nSPs_0 -> Draw("same");
  //nSPs_1 -> Draw("same");
  //nSPs_2 -> Draw("same");
  //nSPs_3 -> Draw("same");
  //nSPs_4 -> Draw("same");
  //c1 -> Print( pdf, "pdf" );

  // divided by pass-hist and not-pass-hist
  nFTFs1_0->Divide(nFTFs1_all);
  nFTFs1_1->Divide(nFTFs1_all);
  nFTFs1_2->Divide(nFTFs1_all);
  nFTFs1_3->Divide(nFTFs1_all);
  nFTFs1_4->Divide(nFTFs1_all);

  THStack *nFTFs1_pt = new THStack("hs_pt",Form(";%s;Fraction of #FTFs (p_{T} > 10GeV)",m_h_pt_numberOfFTFs1->GetXaxis()->GetTitle()));
  nFTFs1_0->SetFillColor(kYellow);
  nFTFs1_1->SetFillColor(kGreen+2);
  nFTFs1_2->SetFillColor(kBlue);
  nFTFs1_3->SetFillColor(kRed);
  nFTFs1_4->SetFillColor(kBlack);

  nFTFs1_pt->Add(nFTFs1_4);
  nFTFs1_pt->Add(nFTFs1_3);
  nFTFs1_pt->Add(nFTFs1_2);
  nFTFs1_pt->Add(nFTFs1_1);
  nFTFs1_pt->Add(nFTFs1_0);

  gStyle->SetHistTopMargin(0);
  nFTFs1_pt->Draw();
  nFTFs1_pt->SetMaximum(1.1);
  nFTFs1_pt->Draw();

  TLegend *leg_nFTFs1_pt = new TLegend(0.02,0.00,0.1,0.4);
  leg_nFTFs1_pt->AddEntry(nFTFs1_0,"0","f");
  leg_nFTFs1_pt->AddEntry(nFTFs1_1,"1,2","f");
  leg_nFTFs1_pt->AddEntry(nFTFs1_2,"3,4","f");
  leg_nFTFs1_pt->AddEntry(nFTFs1_3,"5,6","f");
  leg_nFTFs1_pt->AddEntry(nFTFs1_4,"7~","f");
  leg_nFTFs1_pt->Draw();

  c1 -> Print( label, "pdf" );

  // ================fraction of low dR #FTFs==============
  TH1D* nFTFs2_all = m_h_pt_numberOfFTFs2->ProjectionX("nFTFs_all");
  TH1D* nFTFs2_0   = m_h_pt_numberOfFTFs2->ProjectionX("nFTFs_0", m_h_pt_numberOfFTFs2->GetYaxis()->FindBin(-1),  m_h_pt_numberOfFTFs2->GetYaxis()->FindBin(0.5));
  TH1D* nFTFs2_1   = m_h_pt_numberOfFTFs2->ProjectionX("nFTFs_1", m_h_pt_numberOfFTFs2->GetYaxis()->FindBin(0.6), m_h_pt_numberOfFTFs2->GetYaxis()->FindBin(2.5));
  TH1D* nFTFs2_2   = m_h_pt_numberOfFTFs2->ProjectionX("nFTFs_2", m_h_pt_numberOfFTFs2->GetYaxis()->FindBin(2.6), m_h_pt_numberOfFTFs2->GetYaxis()->FindBin(4.5));
  TH1D* nFTFs2_3   = m_h_pt_numberOfFTFs2->ProjectionX("nFTFs_3", m_h_pt_numberOfFTFs2->GetYaxis()->FindBin(4.6), m_h_pt_numberOfFTFs2->GetYaxis()->FindBin(6.5));
  TH1D* nFTFs2_4   = m_h_pt_numberOfFTFs2->ProjectionX("nFTFs_4", m_h_pt_numberOfFTFs2->GetYaxis()->FindBin(6.6), m_h_pt_numberOfFTFs2->GetYaxis()->FindBin(11));

  //nSPs_all -> Draw();
  //nSPs_0 -> Draw("same");
  //nSPs_1 -> Draw("same");
  //nSPs_2 -> Draw("same");
  //nSPs_3 -> Draw("same");
  //nSPs_4 -> Draw("same");
  //c1 -> Print( pdf, "pdf" );

  // divided by pass-hist and not-pass-hist
  nFTFs2_0->Divide(nFTFs2_all);
  nFTFs2_1->Divide(nFTFs2_all);
  nFTFs2_2->Divide(nFTFs2_all);
  nFTFs2_3->Divide(nFTFs2_all);
  nFTFs2_4->Divide(nFTFs2_all);

  THStack *nFTFs2_pt = new THStack("hs_pt",Form(";%s;Fraction of #FTFs (dR < 0.1)",m_h_pt_numberOfFTFs2->GetXaxis()->GetTitle()));
  nFTFs2_0->SetFillColor(kYellow);
  nFTFs2_1->SetFillColor(kGreen+2);
  nFTFs2_2->SetFillColor(kBlue);
  nFTFs2_3->SetFillColor(kRed);
  nFTFs2_4->SetFillColor(kBlack);

  nFTFs2_pt->Add(nFTFs2_4);
  nFTFs2_pt->Add(nFTFs2_3);
  nFTFs2_pt->Add(nFTFs2_2);
  nFTFs2_pt->Add(nFTFs2_1);
  nFTFs2_pt->Add(nFTFs2_0);

  gStyle->SetHistTopMargin(0);
  nFTFs2_pt->Draw();
  nFTFs2_pt->SetMaximum(1.1);
  nFTFs2_pt->Draw();

  TLegend *leg_nFTFs2_pt = new TLegend(0.02,0.00,0.1,0.4);
  leg_nFTFs2_pt->AddEntry(nFTFs2_0,"0","f");
  leg_nFTFs2_pt->AddEntry(nFTFs2_1,"1,2","f");
  leg_nFTFs2_pt->AddEntry(nFTFs2_2,"3,4","f");
  leg_nFTFs2_pt->AddEntry(nFTFs2_3,"5,6","f");
  leg_nFTFs2_pt->AddEntry(nFTFs2_4,"7~","f");
  leg_nFTFs2_pt->Draw();

  c1 -> Print( label, "pdf" );

  // ================fraction of #SPs SA==============
  TH1D* nSPs_all = m_h_pt_numberOfSPs_muSA->ProjectionX("nSPs_all");
  TH1D* nSPs_0   = m_h_pt_numberOfSPs_muSA->ProjectionX("nSPs_0", m_h_pt_numberOfSPs_muSA->GetYaxis()->FindBin(-0.5),  m_h_pt_numberOfSPs_muSA->GetYaxis()->FindBin(0.5));
  TH1D* nSPs_1   = m_h_pt_numberOfSPs_muSA->ProjectionX("nSPs_1", m_h_pt_numberOfSPs_muSA->GetYaxis()->FindBin(0.6), m_h_pt_numberOfSPs_muSA->GetYaxis()->FindBin(1.5));
  TH1D* nSPs_2   = m_h_pt_numberOfSPs_muSA->ProjectionX("nSPs_2", m_h_pt_numberOfSPs_muSA->GetYaxis()->FindBin(1.6), m_h_pt_numberOfSPs_muSA->GetYaxis()->FindBin(2.5));
  TH1D* nSPs_3   = m_h_pt_numberOfSPs_muSA->ProjectionX("nSPs_3", m_h_pt_numberOfSPs_muSA->GetYaxis()->FindBin(2.6), m_h_pt_numberOfSPs_muSA->GetYaxis()->FindBin(3.5));
  TH1D* nSPs_4   = m_h_pt_numberOfSPs_muSA->ProjectionX("nSPs_4", m_h_pt_numberOfSPs_muSA->GetYaxis()->FindBin(3.6), m_h_pt_numberOfSPs_muSA->GetYaxis()->FindBin(10));

  //nSPs_all -> Draw();
  //nSPs_0 -> Draw("same");
  //nSPs_1 -> Draw("same");
  //nSPs_2 -> Draw("same");
  //nSPs_3 -> Draw("same");
  //nSPs_4 -> Draw("same");
  //c1 -> Print( pdf, "pdf" );

  // divided by pass-hist and not-pass-hist
  nSPs_0->Divide(nSPs_all);
  nSPs_1->Divide(nSPs_all);
  nSPs_2->Divide(nSPs_all);
  nSPs_3->Divide(nSPs_all);
  nSPs_4->Divide(nSPs_all);

  THStack *nSPs_pt = new THStack("hs_pt",Form(";%s;Fraction of #SPs for each RoI",m_h_pt_numberOfSPs_muSA->GetXaxis()->GetTitle()));
  nSPs_0->SetFillColor(kYellow);
  nSPs_1->SetFillColor(kGreen+2);
  nSPs_2->SetFillColor(kBlue);
  nSPs_3->SetFillColor(kRed);
  nSPs_4->SetFillColor(kBlack);

  nSPs_pt->Add(nSPs_4);
  nSPs_pt->Add(nSPs_3);
  nSPs_pt->Add(nSPs_2);
  nSPs_pt->Add(nSPs_1);
  nSPs_pt->Add(nSPs_0);

  gStyle->SetHistTopMargin(0);
  nSPs_pt->Draw();
  nSPs_pt->GetXaxis()->SetRangeUser(0,20);
  nSPs_pt->SetMaximum(1.1);
  nSPs_pt->Draw();

  TLegend *leg_nSPs_pt = new TLegend(0.02,0.00,0.1,0.4);
  leg_nSPs_pt->AddEntry(nSPs_0,"0","f");
  leg_nSPs_pt->AddEntry(nSPs_1,"1","f");
  leg_nSPs_pt->AddEntry(nSPs_2,"2","f");
  leg_nSPs_pt->AddEntry(nSPs_3,"3","f");
  leg_nSPs_pt->AddEntry(nSPs_4,"4~","f");
  leg_nSPs_pt->Draw();

  c1 -> Print( label, "pdf" );



  // ================fraction of #SPs SAIO==============
  TH1D* nSPs1_all = m_h_pt_numberOfSPs_muSAIO->ProjectionX("nSPs1_all");
  TH1D* nSPs1_0   = m_h_pt_numberOfSPs_muSAIO->ProjectionX("nSPs1_0", m_h_pt_numberOfSPs_muSAIO->GetYaxis()->FindBin(-10),  m_h_pt_numberOfSPs_muSAIO->GetYaxis()->FindBin(0.5));
  TH1D* nSPs1_1   = m_h_pt_numberOfSPs_muSAIO->ProjectionX("nSPs1_1", m_h_pt_numberOfSPs_muSAIO->GetYaxis()->FindBin(0.6), m_h_pt_numberOfSPs_muSAIO->GetYaxis()->FindBin(1.5));
  TH1D* nSPs1_2   = m_h_pt_numberOfSPs_muSAIO->ProjectionX("nSPs1_2", m_h_pt_numberOfSPs_muSAIO->GetYaxis()->FindBin(1.6), m_h_pt_numberOfSPs_muSAIO->GetYaxis()->FindBin(2.5));
  TH1D* nSPs1_3   = m_h_pt_numberOfSPs_muSAIO->ProjectionX("nSPs1_3", m_h_pt_numberOfSPs_muSAIO->GetYaxis()->FindBin(2.6), m_h_pt_numberOfSPs_muSAIO->GetYaxis()->FindBin(3.5));
  TH1D* nSPs1_4   = m_h_pt_numberOfSPs_muSAIO->ProjectionX("nSPs1_4", m_h_pt_numberOfSPs_muSAIO->GetYaxis()->FindBin(3.6), m_h_pt_numberOfSPs_muSAIO->GetYaxis()->FindBin(12));

  //nSPs_all -> Draw();
  //nSPs_0 -> Draw("same");
  //nSPs_1 -> Draw("same");
  //nSPs_2 -> Draw("same");
  //nSPs_3 -> Draw("same");
  //nSPs_4 -> Draw("same");
  //c1 -> Print( pdf, "pdf" );

  // divided by pass-hist and not-pass-hist
  nSPs1_0->Divide(nSPs1_all);
  nSPs1_1->Divide(nSPs1_all);
  nSPs1_2->Divide(nSPs1_all);
  nSPs1_3->Divide(nSPs1_all);
  nSPs1_4->Divide(nSPs1_all);

  THStack *nSPs1_pt = new THStack("hs1_pt",Form(";%s;Fraction of #SPs for each RoI",m_h_pt_numberOfSPs_muSAIO->GetXaxis()->GetTitle()));
  nSPs1_0->SetFillColor(kYellow);
  nSPs1_1->SetFillColor(kGreen+2);
  nSPs1_2->SetFillColor(kBlue);
  nSPs1_3->SetFillColor(kRed);
  nSPs1_4->SetFillColor(kBlack);

  nSPs1_pt->Add(nSPs1_4);
  nSPs1_pt->Add(nSPs1_3);
  nSPs1_pt->Add(nSPs1_2);
  nSPs1_pt->Add(nSPs1_1);
  nSPs1_pt->Add(nSPs1_0);

  gStyle->SetHistTopMargin(0);
  nSPs1_pt->Draw();
  nSPs1_pt->GetXaxis()->SetRangeUser(0,20);
  nSPs1_pt->SetMaximum(1.1);
  nSPs1_pt->Draw();

  TLegend *leg_nSPs1_pt = new TLegend(0.02,0.00,0.1,0.4);
  leg_nSPs1_pt->AddEntry(nSPs1_0,"0","f");
  leg_nSPs1_pt->AddEntry(nSPs1_1,"1","f");
  leg_nSPs1_pt->AddEntry(nSPs1_2,"2","f");
  leg_nSPs1_pt->AddEntry(nSPs1_3,"3","f");
  leg_nSPs1_pt->AddEntry(nSPs1_4,"4~","f");
  leg_nSPs1_pt->Draw();

  c1 -> Print( label, "pdf" );

  // ================fraction of #SPs SAIO==============
  TH1D* nSPs2_all = m_h_pt_numberOfSPs_muSAIO_1->ProjectionX("nSPs2_all");
  TH1D* nSPs2_0   = m_h_pt_numberOfSPs_muSAIO_1->ProjectionX("nSPs2_0", m_h_pt_numberOfSPs_muSAIO_1->GetYaxis()->FindBin(-10), m_h_pt_numberOfSPs_muSAIO_1->GetYaxis()->FindBin(0.5));
  TH1D* nSPs2_1   = m_h_pt_numberOfSPs_muSAIO_1->ProjectionX("nSPs2_1", m_h_pt_numberOfSPs_muSAIO_1->GetYaxis()->FindBin(0.6), m_h_pt_numberOfSPs_muSAIO_1->GetYaxis()->FindBin(1.5));
  TH1D* nSPs2_2   = m_h_pt_numberOfSPs_muSAIO_1->ProjectionX("nSPs2_2", m_h_pt_numberOfSPs_muSAIO_1->GetYaxis()->FindBin(1.6), m_h_pt_numberOfSPs_muSAIO_1->GetYaxis()->FindBin(2.5));
  TH1D* nSPs2_3   = m_h_pt_numberOfSPs_muSAIO_1->ProjectionX("nSPs2_3", m_h_pt_numberOfSPs_muSAIO_1->GetYaxis()->FindBin(2.6), m_h_pt_numberOfSPs_muSAIO_1->GetYaxis()->FindBin(3.5));
  TH1D* nSPs2_4   = m_h_pt_numberOfSPs_muSAIO_1->ProjectionX("nSPs2_4", m_h_pt_numberOfSPs_muSAIO_1->GetYaxis()->FindBin(3.6), m_h_pt_numberOfSPs_muSAIO_1->GetYaxis()->FindBin(12));

  //nSPs_all -> Draw();
  //nSPs_0 -> Draw("same");
  //nSPs_1 -> Draw("same");
  //nSPs_2 -> Draw("same");
  //nSPs_3 -> Draw("same");
  //nSPs_4 -> Draw("same");
  //c1 -> Print( pdf, "pdf" );

  // divided by pass-hist and not-pass-hist
  nSPs2_0->Divide(nSPs2_all);
  nSPs2_1->Divide(nSPs2_all);
  nSPs2_2->Divide(nSPs2_all);
  nSPs2_3->Divide(nSPs2_all);
  nSPs2_4->Divide(nSPs2_all);

  THStack *nSPs2_pt = new THStack("hs1_pt",Form(";%s;Fraction of #SPs for each RoI",m_h_pt_numberOfSPs_muSAIO->GetXaxis()->GetTitle()));
  nSPs2_0->SetFillColor(kYellow);
  nSPs2_1->SetFillColor(kGreen+2);
  nSPs2_2->SetFillColor(kBlue);
  nSPs2_3->SetFillColor(kRed);
  nSPs2_4->SetFillColor(kBlack);

  nSPs2_pt->Add(nSPs2_4);
  nSPs2_pt->Add(nSPs2_3);
  nSPs2_pt->Add(nSPs2_2);
  nSPs2_pt->Add(nSPs2_1);
  nSPs2_pt->Add(nSPs2_0);

  gStyle->SetHistTopMargin(0);
  nSPs2_pt->Draw();
  nSPs2_pt->GetXaxis()->SetRangeUser(0,20);
  nSPs2_pt->SetMaximum(1.1);
  nSPs2_pt->Draw();

  TLegend *leg_nSPs2_pt = new TLegend(0.02,0.00,0.1,0.4);
  leg_nSPs2_pt->AddEntry(nSPs2_0,"0","f");
  leg_nSPs2_pt->AddEntry(nSPs2_1,"1","f");
  leg_nSPs2_pt->AddEntry(nSPs2_2,"2","f");
  leg_nSPs2_pt->AddEntry(nSPs2_3,"3","f");
  leg_nSPs2_pt->AddEntry(nSPs2_4,"4~","f");
  leg_nSPs2_pt->Draw();

  c1 -> Print( label, "pdf" );







  //=============single muon eff========

  TH1F* frame_pt = c1->DrawFrame(0,0,20,1.2);
  TEfficiency *h_pt_L1SA = new TEfficiency(*h_pt_SA, *h_pt_L1);
  h_pt_L1SA->Draw("same");
  TEfficiency *h_pt_L1SAIO = new TEfficiency(*h_pt_CB, *h_pt_SA);
  h_pt_L1SAIO->SetLineColor(kRed);
  h_pt_L1SAIO->SetMarkerColor(kRed);
  h_pt_L1SAIO->Draw("same");
  std::cout << "debug" << std::endl;
  c1 -> Print(label, "pdf");

  TH1F* frame_pt_1 = c1->DrawFrame(0,0,20,1.2);
  TEfficiency *h_pt_L1SA_1 = new TEfficiency(*h_pt_SA_1, *h_pt_L1_1);
  h_pt_L1SA_1->Draw("same");
  TEfficiency *h_pt_L1SAIO_1 = new TEfficiency(*h_pt_CB_1, *h_pt_SA_1);
  h_pt_L1SAIO_1->SetLineColor(kRed);
  h_pt_L1SAIO_1->SetMarkerColor(kRed);
  h_pt_L1SAIO_1->Draw("same");
  std::cout << "debug" << std::endl;
  c1 -> Print(label, "pdf");



  //========sep==========

  TH1F* frame0 = c1->DrawFrame(0,0,0.4,1000);
  h_L1_sep->SetLineColor(kBlack);
  h_SA_sep->SetLineColor(kBlue);
  h_SAIO_sep->SetLineColor(kRed);
  h_L1_sep->SetMarkerColor(kBlack);
  h_SA_sep->SetMarkerColor(kBlue);
  h_SAIO_sep->SetMarkerColor(kRed);
  h_L1_sep->Draw("same");
  h_SA_sep->Draw("same");
  h_SAIO_sep->Draw("same");
  h_CBIO_sep->Draw("same");
  std::cout << "debug" << std::endl;
  c1 -> Print(label, "pdf");


  TH1F* frame1 = c1->DrawFrame(0,0,0.4,1.2);
  frame1->GetXaxis()->SetTitle("p_{T}[GeV]");
  frame1->GetYaxis()->SetTitle("/separated-L1");
  TEfficiency *hL1SA_sep = new TEfficiency(*h_SA_sep, *h_L1_sep); //SApass+SAsep/L1pass+L1sep
  hL1SA_sep->Draw("same");
  //hL1SA_sep->GetYaxis()->SetTitle("/separated-L1");
  //hL1SA_sep->SetTitle("p_{T} [GeV];/separated-L1;");
  //TEfficiency *hL1CB_sep = new TEfficiency(*h_CB_sep, *h_L1_sep); //CBpass+CBsep/L1pass+L1sep
  //hL1CB_sep->SetLineColor(kBlue);
  //hL1CB_sep->SetMarkerColor(kBlue);
  //hL1CB_sep->Draw("same");
  TEfficiency *hL1CBIO_sep = new TEfficiency(*h_CBIO_sep, *h_L1_sep); //CBIOpass+CBIOsep/L1pass+L1sep
  hL1CBIO_sep->SetLineColor(kRed);
  hL1CBIO_sep->SetMarkerColor(kRed);
  //hL1CBIO_sep->SetTitle("p_{T} [GeV];/separated-L1;");
  hL1CBIO_sep->Draw("same");
  //TEfficiency *hL1SAIO = new TEfficiency(*h_SAIO, *h_L1);
  //hL1SAIO->SetLineColor(kRed);
  //hL1SAIO->SetMarkerColor(kRed);
  //hL1SAIO->Draw("same");
  std::cout << "debug" << std::endl;
  c1 -> Print(label, "pdf");

  TH1F* frame2 = c1->DrawFrame(0,0,0.4,1.2);
  frame2->GetXaxis()->SetTitle("p_{T}[GeV]");
  frame2->GetYaxis()->SetTitle("/L1");
  TEfficiency *hL1SA_1_sep = new TEfficiency(*h_SA_1_sep, *h_L1_1_sep); //SAsep+SAsep/L1pass
  hL1SA_1_sep->Draw("same");
  TEfficiency *hL1CB_1_sep = new TEfficiency(*h_CB_1_sep, *h_L1_1_sep); //CBsep+SAsep/L1pass
  hL1CB_1_sep->SetLineColor(kBlack);
  hL1CB_1_sep->SetMarkerColor(kBlack);
  hL1CB_1_sep->Draw("same");
  TEfficiency *hL1CBIO_1_sep = new TEfficiency(*h_CBIO_1_sep, *h_L1_1_sep); //CBIOpass+SAsep/L1pass
  hL1CBIO_1_sep->SetLineColor(kRed);
  hL1CBIO_1_sep->SetMarkerColor(kRed);
  hL1CBIO_1_sep->Draw("same");
  //TEfficiency *hL1SAIO_1 = new TEfficiency(*h_SAIO_1, *h_L1_1);
  //hL1SAIO_1->SetLineColor(kBlue);
  //hL1SAIO_1->SetMarkerColor(kBlue);
  //hL1SAIO_1->Draw("same");
  std::cout << "debug" << std::endl;
  c1 -> Print(label, "pdf");

  TH1F* frame3 = c1->DrawFrame(0,0,0.4,1.2);
  TEfficiency *hL1SA_2_sep = new TEfficiency(*h_SA_2_sep, *h_L1_2_sep); //SAsep+SAsep/L1pass+!L1sep
  hL1SA_2_sep->Draw("same");
  TEfficiency *hL1CBIO_2_sep = new TEfficiency(*h_CBIO_2_sep, *h_L1_2_sep); //CBIOsep+CBIOsep/L1pass+!L1sep
  hL1CBIO_2_sep->SetLineColor(kRed);
  hL1CBIO_2_sep->SetMarkerColor(kRed);
  hL1CBIO_2_sep->Draw("same");
  //TEfficiency *hL1SAIO_2 = new TEfficiency(*h_SAIO_2, *h_L1_2);
  //hL1SAIO_2->SetLineColor(kBlue);
  //hL1SAIO_2->SetMarkerColor(kBlue);
  //hL1SAIO_2->Draw("same");
  std::cout << "debug" << std::endl;
  c1 -> Print(label, "pdf");


  TH1F* frame4 = c1->DrawFrame(0,0,0.4,1.2);
  //TEfficiency *hL1SA_3 = new TEfficiency(*h_SA_3, *h_L1_3);
  //hL1SA_3->Draw("same");
  TEfficiency *hL1CBIO_3_sep = new TEfficiency(*h_CB_3_sep, *h_SA_3_sep);
  hL1CBIO_3_sep->SetLineColor(kRed);
  hL1CBIO_3_sep->SetMarkerColor(kRed);
  hL1CBIO_3_sep->Draw("same");
  //TEfficiency *hSACBIO_3_sep = new TEfficiency(*h_CBIO_3_sep, *h_SA_3_sep);
  //hSACBIO_3_sep->SetLineColor(kBlue);
  //hSACBIO_3_sep->SetMarkerColor(kBlue);
  //hSACBIO_3_sep->Draw("same");
  //TEfficiency *hL1SAIO_2 = new TEfficiency(*h_SAIO_2, *h_L1_2);
  //hL1SAIO_2->SetLineColor(kBlue);
  //hL1SAIO_2->SetMarkerColor(kBlue);
  //hL1SAIO_2->Draw("same");
  std::cout << "debug" << std::endl;
  c1 -> Print(label, "pdf");


  TH1F* frame5 = c1->DrawFrame(0,0,0.4,1.2);
  //TEfficiency *hL1SA_4 = new TEfficiency(*h_SA_4, *h_L1_4);
  //hL1SA_4->Draw("same");
  TEfficiency *hL1CBIO_4_sep = new TEfficiency(*h_CBIO_4_sep, *h_SA_4_sep); //CBIOsep+CBIOsep/SApass+!SAsep+L1pass+L1sep
  hL1CBIO_4_sep->SetLineColor(kRed);
  hL1CBIO_4_sep->SetMarkerColor(kRed);
  hL1CBIO_4_sep->Draw();
  //TEfficiency *hSACBIO_4_sep = new TEfficiency(*h_CBIO_4_sep, *h_SA_4_sep); //SAsep+SAsep/L1pass+L1sep
  //hSACBIO_4_sep->SetLineColor(kBlue);
  //hSACBIO_4_sep->SetMarkerColor(kBlue);
  //hSACBIO_4_sep->Draw("same");
  //TEfficiency *hL1SAIO_2 = new TEfficiency(*h_SAIO_2, *h_L1_2);
  //hL1SAIO_2->SetLineColor(kBlue);
  //hL1SAIO_2->SetMarkerColor(kBlue);
  //hL1SAIO_2->Draw("same");
  std::cout << "debug" << std::endl;
  c1 -> Print(label, "pdf");


  TH1F* frame6 = c1->DrawFrame(0,0,0.4,1.2);
  //TEfficiency *hL1SA_4 = new TEfficiency(*h_SA_4, *h_L1_4);
  //hL1SA_4->Draw("same");
  TEfficiency *hL1CBIO_5_sep = new TEfficiency(*h_CBIO_5_sep, *h_L1_5_sep); //CBIOsep+CBIOsep/SApass+!SAsep+L1pass
  hL1CBIO_5_sep->SetLineColor(kRed);
  hL1CBIO_5_sep->SetMarkerColor(kRed);
  hL1CBIO_5_sep->Draw();
  //TEfficiency *hSACBIO_5_sep = new TEfficiency(*h_CBIO_5_sep, *h_SA_5_sep); //CBIOsep+CBIOsep/SApass+!SAsep+L1pass+L1sep
  //hSACBIO_5_sep->SetLineColor(kBlue);
  //hSACBIO_5_sep->SetMarkerColor(kBlue);
  //hSACBIO_5_sep->Draw("same");
  //TEfficiency *hL1SAIO_2 = new TEfficiency(*h_SAIO_2, *h_L1_2);
  //hL1SAIO_2->SetLineColor(kBlue);
  //hL1SAIO_2->SetMarkerColor(kBlue);
  //hL1SAIO_2->Draw("same");
  std::cout << "debug" << std::endl;
  c1 -> Print(label, "pdf");


  TH1F* frame7 = c1->DrawFrame(0,0,0.4,1.2);
  //TEfficiency *hL1SA_4 = new TEfficiency(*h_SA_4, *h_L1_4);
  //hL1SA_4->Draw("same");
  TEfficiency *hL1CB_6_sep = new TEfficiency(*h_CB_6_sep, *h_L1_6_sep); //SApass+CBpass+CBsep/L1pass
  hL1CB_6_sep->SetLineColor(kRed);
  hL1CB_6_sep->SetMarkerColor(kRed);
  hL1CB_6_sep->Draw("same");
  TEfficiency *hL1SA_6_sep = new TEfficiency(*h_SA_6_sep, *h_L1_6_sep); //SApass+CBpass+CBsep/L1pass
  hL1SA_6_sep->SetLineColor(kBlue);
  hL1SA_6_sep->SetMarkerColor(kBlue);
  hL1SA_6_sep->Draw("same");
  //TEfficiency *hSACBIO_5_sep = new TEfficiency(*h_CBIO_5_sep, *h_SA_5_sep); //CBIOsep+CBIOsep/SApass+!SAsep+L1pass+L1sep
  //hSACBIO_5_sep->SetLineColor(kBlue);
  //hSACBIO_5_sep->SetMarkerColor(kBlue);
  //hSACBIO_5_sep->Draw("same");
  //TEfficiency *hL1SAIO_2 = new TEfficiency(*h_SAIO_2, *h_L1_2);
  //hL1SAIO_2->SetLineColor(kBlue);
  //hL1SAIO_2->SetMarkerColor(kBlue);
  //hL1SAIO_2->Draw("same");
  std::cout << "debug" << std::endl;
  c1 -> Print(label, "pdf");









  std::cout << "debug" << std::endl;
  c1 -> Print(label + "]", "pdf");
}
