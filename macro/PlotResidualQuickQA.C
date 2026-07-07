#include "TCanvas.h"
#include "TChain.h"
#include "TColor.h"
#include "TError.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TObject.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"

#include <iostream>
#include <string>
#include <vector>

// -------------------------------
// Cut and plotting configuration
// -------------------------------
const double kVzMin = -15.0;
const double kVzMax = 15.0;
const double kVzMinDCA = -10.0;
const double kVzMaxDCA = 10.0;
const double kPtMinResidualVsR = 2.0;
const double kPtMinResidualEtaSplit = 0.5;
const double kEtaAbsSplit = 0.5;
const double kEtaAbsMaxResidualVsR = 1;
const double kQualityPtMin = 0.9;
const int kMinTpcClusters = 30;

const int kCanvasWidth = 900;
const int kCanvasHeight = 760;
const bool kBatchMode = true;

std::string AndCuts(const std::vector<std::string>& cuts)
{
  std::string out;
  for (const auto& cut : cuts)
  {
    if (cut.empty())
    {
      continue;
    }
    if (!out.empty())
    {
      out += "&&";
    }
    out += "(" + cut + ")";
  }
  return out;
}

void DrawLabel(const std::string& label, const std::string& cut)
{
  TLatex latex;
  latex.SetNDC(true);
  latex.SetTextFont(42);
  latex.SetTextSize(0.038);
  latex.DrawLatex(0.13, 0.935, label.c_str());

  latex.SetTextSize(0.027);
  latex.DrawLatex(0.13, 0.895, cut.c_str());
}

void StyleHist(TH1* hist, const char* xTitle, const char* yTitle)
{
  hist->SetTitle("");
  hist->SetLineColor(kAzure + 2);
  hist->SetMarkerColor(kAzure + 2);
  hist->SetLineWidth(2);
  hist->GetXaxis()->SetTitle(xTitle);
  hist->GetYaxis()->SetTitle(yTitle);
  hist->GetXaxis()->CenterTitle(true);
  hist->GetYaxis()->CenterTitle(true);
  hist->GetYaxis()->SetTitleOffset(1.35);
}

void Style2D(TH2* hist, const char* xTitle, const char* yTitle)
{
  StyleHist(hist, xTitle, yTitle);
  hist->GetZaxis()->SetTitle("Entries");
  hist->GetZaxis()->CenterTitle(true);
}

void SaveCanvasAndObject(TFile* fout, TCanvas* canvas, TH1* hist,
                         const std::string& outdir, const std::string& name)
{
  fout->cd();
  hist->Write("", TObject::kOverwrite);
  canvas->Write("", TObject::kOverwrite);

  const std::string png = outdir + "/" + name + ".png";
  canvas->SaveAs(png.c_str());
}

void DrawAndSave2D(TChain& tree, TFile* fout,
                   const std::string& histName,
                   const std::string& drawExpr,
                   const std::string& cut,
                   const std::string& plotLabel,
                   const std::string& outdir,
                   const char* xTitle,
                   const char* yTitle,
                   int nBinsX,
                   double xMin,
                   double xMax,
                   int nBinsY,
                   double yMin,
                   double yMax)
{
  TH2F* hist = new TH2F(histName.c_str(), "", nBinsX, xMin, xMax, nBinsY, yMin, yMax);
  Style2D(hist, xTitle, yTitle);

  const std::string command = drawExpr + ">>" + histName;
  tree.Draw(command.c_str(), cut.c_str(), "goff");

  const std::string canvasName = "c_" + histName;
  TCanvas* canvas = new TCanvas(canvasName.c_str(), "", kCanvasWidth, kCanvasHeight);
  canvas->SetRightMargin(0.16);
  canvas->SetTopMargin(0.10);
  hist->Draw("COLZ");
  DrawLabel(plotLabel, cut);
  SaveCanvasAndObject(fout, canvas, hist, outdir, histName);
}
// .x PlotResidualQuickQA.C("/sphenix/tg/tg01/hf/mitrankova/PatternReco/79513_keff_m1/outout_resid/tpc_poly_cluster_residualsHITS_ppFieldOn_79513*.root","79513_keff_m1")
// .x PlotResidualQuickQA.C("/sphenix/tg/tg01/hf/mitrankova/PatternReco/79513_keff_1/outout_resid/tpc_poly_cluster_residualsHITS_ppFieldOn_79513*.root","79513_keff_1")
// .x PlotResidualQuickQA.C("/sphenix/tg/tg01/hf/mitrankova/PatternReco/79513_keff_2/outout_resid/tpc_poly_cluster_residualsHITS_ppFieldOn_79513*.root","79513_keff_2")
// .x PlotResidualQuickQA.C("/sphenix/tg/tg01/hf/mitrankova/PatternReco/79513_keff_4/outout_resid/tpc_poly_cluster_residualsHITS_ppFieldOn_79513*.root","79513_keff_4")
// .x PlotResidualQuickQA.C("/sphenix/tg/tg01/hf/mitrankova/PatternReco/79513/outout_resid/tpc_poly_cluster_residualsHITS_ppFieldOn_79513*.root","79513")
void PlotResidualQuickQA(const std::string& inputfile = "/sphenix/tg/tg01/hf/mitrankova/PatternReco/79513/outout_resid/tpc_poly_cluster_residualsHITS_ppFieldOn_79513*.root",
                         const std::string& runnumber = "79513",
                         const std::string& collisionSystem = "p+p",
                         const std::string& outputBaseDir = "ResidualQuickQA")
{
  if (kBatchMode)
  {
    gROOT->SetBatch(kTRUE);
  }

  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBird);
  gStyle->SetNumberContours(255);
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetTitleSize(0.045, "XYZ");
  gStyle->SetLabelSize(0.038, "XYZ");

  const std::string label = collisionSystem + "+" + runnumber;
  const std::string outdir = outputBaseDir + "/" + collisionSystem + "_" + runnumber;
  gSystem->mkdir(outputBaseDir.c_str(), kTRUE);
  gSystem->mkdir(outdir.c_str(), kTRUE);

  TChain residuals("residuals");
  const int nfiles = residuals.Add(inputfile.c_str());
  if (nfiles <= 0 || residuals.GetEntries() <= 0)
  {
    std::cerr << "No entries found in tree 'residuals' from input: " << inputfile << std::endl;
    return;
  }

  const std::string outfile = outdir + "/" + collisionSystem + "_" + runnumber + "_residual_quickQA.root";
  TFile* fout = TFile::Open(outfile.c_str(), "RECREATE");
  if (!fout || fout->IsZombie())
  {
    std::cerr << "Could not create output ROOT file: " << outfile << std::endl;
    return;
  }

  const std::string vzCut = Form("vertex_z>%g&&vertex_z<%g", kVzMin, kVzMax);
  const std::string vzCutDCA = Form("vertex_z>%g&&vertex_z<%g", kVzMinDCA, kVzMaxDCA);
  const std::string tpcClusterCut = Form("ntpc_clusters>%d", kMinTpcClusters);
  const std::string ptResidualVsRCut = Form("pt>%g", kPtMinResidualVsR);
  const std::string ptResidualEtaSplitCut = Form("pt>%g", kPtMinResidualEtaSplit);
  const std::string etaLtSplitCut = Form("abs(eta)<%g", kEtaAbsSplit);
  const std::string etaGtSplitCut = Form("abs(eta)>%g", kEtaAbsSplit);
  const std::string etaResidualVsRCut = Form("abs(eta)<%g", kEtaAbsMaxResidualVsR);
  const std::string qualityPtCut = Form("pt>%g", kQualityPtMin);
  const std::string positiveChargeCut = "charge>0";
  const std::string negativeChargeCut = "charge<0";
  const std::string side0Cut = "side==0";
  const std::string side1Cut = "side==1";

  const std::string cutResidualR = AndCuts({vzCut, ptResidualVsRCut, etaResidualVsRCut});
  //  const std::string cutResidualR = AndCuts({vzCut, ptResidualVsRCut});//GenFit
  const std::string cutResidualRPositive = AndCuts({cutResidualR, positiveChargeCut});
  const std::string cutResidualRNegative = AndCuts({cutResidualR, negativeChargeCut});
  const std::string cutResidualRSide0 = AndCuts({cutResidualR, side0Cut});
  const std::string cutResidualRSide1 = AndCuts({cutResidualR, side1Cut});
  const std::string cutResidualRPositiveSide0 = AndCuts({cutResidualRPositive, side0Cut});
  const std::string cutResidualRPositiveSide1 = AndCuts({cutResidualRPositive, side1Cut});
  const std::string cutResidualRNegativeSide0 = AndCuts({cutResidualRNegative, side0Cut});
  const std::string cutResidualRNegativeSide1 = AndCuts({cutResidualRNegative, side1Cut});
  const std::string cutResidualREtaLtBase = AndCuts({vzCut, ptResidualEtaSplitCut, etaLtSplitCut});
  const std::string cutResidualREtaGtBase = AndCuts({vzCut, ptResidualEtaSplitCut, etaGtSplitCut});
  const std::string cutResidualRPositiveSide0EtaLt = AndCuts({cutResidualREtaLtBase, positiveChargeCut, side0Cut});
  const std::string cutResidualRPositiveSide1EtaLt = AndCuts({cutResidualREtaLtBase, positiveChargeCut, side1Cut});
  const std::string cutResidualRNegativeSide0EtaLt = AndCuts({cutResidualREtaLtBase, negativeChargeCut, side0Cut});
  const std::string cutResidualRNegativeSide1EtaLt = AndCuts({cutResidualREtaLtBase, negativeChargeCut, side1Cut});
  const std::string cutResidualRPositiveSide0EtaGt = AndCuts({cutResidualREtaGtBase, positiveChargeCut, side0Cut});
  const std::string cutResidualRPositiveSide1EtaGt = AndCuts({cutResidualREtaGtBase, positiveChargeCut, side1Cut});
  const std::string cutResidualRNegativeSide0EtaGt = AndCuts({cutResidualREtaGtBase, negativeChargeCut, side0Cut});
  const std::string cutResidualRNegativeSide1EtaGt = AndCuts({cutResidualREtaGtBase, negativeChargeCut, side1Cut});
  const std::string cutResidualPt = vzCut;
  const std::string cutVertexAll = tpcClusterCut;
  const std::string cutVertexPt = AndCuts({tpcClusterCut, ptResidualVsRCut});
  const std::string cutQuality = AndCuts({tpcClusterCut, qualityPtCut});
  const std::string cutDca = AndCuts({tpcClusterCut, vzCutDCA});
  const std::string cutDcaSide0 = AndCuts({cutDca, side0Cut});
  const std::string cutDcaSide1 = AndCuts({cutDca, side1Cut});

  std::cout << "Input: " << inputfile << std::endl;
  std::cout << "Entries: " << residuals.GetEntries() << std::endl;
  std::cout << "Output directory: " << outdir << std::endl;

  // residual_rphi vs cluster_r
  TH2F* h_residual_rphi_vs_cluster_r = new TH2F("h_residual_rphi_vs_cluster_r", "",
                                                55, 20, 75, 200, -0.8, 0.8);
  Style2D(h_residual_rphi_vs_cluster_r, "cluster r [cm]", "r#Delta#phi residual [cm]");
  residuals.Draw("residual_rphi:cluster_r>>h_residual_rphi_vs_cluster_r", cutResidualR.c_str(), "goff");
  TCanvas* c_residual_rphi_vs_cluster_r = new TCanvas("c_residual_rphi_vs_cluster_r", "",
                                                      kCanvasWidth, kCanvasHeight);
  c_residual_rphi_vs_cluster_r->SetRightMargin(0.16);
  c_residual_rphi_vs_cluster_r->SetTopMargin(0.10);
  h_residual_rphi_vs_cluster_r->Draw("COLZ");
  DrawLabel(label, cutResidualR);
  SaveCanvasAndObject(fout, c_residual_rphi_vs_cluster_r, h_residual_rphi_vs_cluster_r,
                      outdir, "residual_rphi_vs_cluster_r");

  // residual_rphi vs cluster_r, positive charge
  TH2F* h_residual_rphi_vs_cluster_r_positive =
      new TH2F("h_residual_rphi_vs_cluster_r_positive", "", 55, 20, 75, 200, -0.8, 0.8);
  Style2D(h_residual_rphi_vs_cluster_r_positive, "cluster r [cm]", "r#Delta#phi residual [cm]");
  residuals.Draw("residual_rphi:cluster_r>>h_residual_rphi_vs_cluster_r_positive",
                 cutResidualRPositive.c_str(), "goff");
  TCanvas* c_residual_rphi_vs_cluster_r_positive =
      new TCanvas("c_residual_rphi_vs_cluster_r_positive", "", kCanvasWidth, kCanvasHeight);
  c_residual_rphi_vs_cluster_r_positive->SetRightMargin(0.16);
  c_residual_rphi_vs_cluster_r_positive->SetTopMargin(0.10);
  h_residual_rphi_vs_cluster_r_positive->Draw("COLZ");
  DrawLabel(label + ", charge > 0", cutResidualRPositive);
  SaveCanvasAndObject(fout, c_residual_rphi_vs_cluster_r_positive,
                      h_residual_rphi_vs_cluster_r_positive, outdir,
                      "residual_rphi_vs_cluster_r_positive");

  // residual_rphi vs cluster_r, negative charge
  TH2F* h_residual_rphi_vs_cluster_r_negative =
      new TH2F("h_residual_rphi_vs_cluster_r_negative", "", 55, 20, 75, 200, -0.8, 0.8);
  Style2D(h_residual_rphi_vs_cluster_r_negative, "cluster r [cm]", "r#Delta#phi residual [cm]");
  residuals.Draw("residual_rphi:cluster_r>>h_residual_rphi_vs_cluster_r_negative",
                 cutResidualRNegative.c_str(), "goff");
  TCanvas* c_residual_rphi_vs_cluster_r_negative =
      new TCanvas("c_residual_rphi_vs_cluster_r_negative", "", kCanvasWidth, kCanvasHeight);
  c_residual_rphi_vs_cluster_r_negative->SetRightMargin(0.16);
  c_residual_rphi_vs_cluster_r_negative->SetTopMargin(0.10);
  h_residual_rphi_vs_cluster_r_negative->Draw("COLZ");
  DrawLabel(label + ", charge < 0", cutResidualRNegative);
  SaveCanvasAndObject(fout, c_residual_rphi_vs_cluster_r_negative,
                      h_residual_rphi_vs_cluster_r_negative, outdir,
                      "residual_rphi_vs_cluster_r_negative");

  DrawAndSave2D(residuals, fout, "residual_rphi_vs_cluster_r_side0",
                "residual_rphi:cluster_r", cutResidualRSide0, label + ", side == 0", outdir,
                "cluster r [cm]", "r#Delta#phi residual [cm]", 55, 20, 75, 200, -0.8, 0.8);
  DrawAndSave2D(residuals, fout, "residual_rphi_vs_cluster_r_side1",
                "residual_rphi:cluster_r", cutResidualRSide1, label + ", side == 1", outdir,
                "cluster r [cm]", "r#Delta#phi residual [cm]", 55, 20, 75, 200, -0.8, 0.8);
  DrawAndSave2D(residuals, fout, "residual_rphi_vs_cluster_r_positive_side0",
                "residual_rphi:cluster_r", cutResidualRPositiveSide0,
                label + ", charge > 0, side == 0", outdir,
                "cluster r [cm]", "r#Delta#phi residual [cm]", 55, 20, 75, 200, -0.8, 0.8);
  DrawAndSave2D(residuals, fout, "residual_rphi_vs_cluster_r_positive_side1",
                "residual_rphi:cluster_r", cutResidualRPositiveSide1,
                label + ", charge > 0, side == 1", outdir,
                "cluster r [cm]", "r#Delta#phi residual [cm]", 55, 20, 75, 200, -0.8, 0.8);
  DrawAndSave2D(residuals, fout, "residual_rphi_vs_cluster_r_negative_side0",
                "residual_rphi:cluster_r", cutResidualRNegativeSide0,
                label + ", charge < 0, side == 0", outdir,
                "cluster r [cm]", "r#Delta#phi residual [cm]", 55, 20, 75, 200, -0.8, 0.8);
  DrawAndSave2D(residuals, fout, "residual_rphi_vs_cluster_r_negative_side1",
                "residual_rphi:cluster_r", cutResidualRNegativeSide1,
                label + ", charge < 0, side == 1", outdir,
                "cluster r [cm]", "r#Delta#phi residual [cm]", 55, 20, 75, 200, -0.8, 0.8);

  DrawAndSave2D(residuals, fout, "residual_rphi_vs_cluster_r_positive_side0_eta_lt_0p5",
                "residual_rphi:cluster_r", cutResidualRPositiveSide0EtaLt,
                label + ", charge > 0, side == 0, |eta| < 0.5", outdir,
                "cluster r [cm]", "r#Delta#phi residual [cm]", 55, 20, 75, 200, -0.8, 0.8);
  DrawAndSave2D(residuals, fout, "residual_rphi_vs_cluster_r_positive_side1_eta_lt_0p5",
                "residual_rphi:cluster_r", cutResidualRPositiveSide1EtaLt,
                label + ", charge > 0, side == 1, |eta| < 0.5", outdir,
                "cluster r [cm]", "r#Delta#phi residual [cm]", 55, 20, 75, 200, -0.8, 0.8);
  DrawAndSave2D(residuals, fout, "residual_rphi_vs_cluster_r_negative_side0_eta_lt_0p5",
                "residual_rphi:cluster_r", cutResidualRNegativeSide0EtaLt,
                label + ", charge < 0, side == 0, |eta| < 0.5", outdir,
                "cluster r [cm]", "r#Delta#phi residual [cm]", 55, 20, 75, 200, -0.8, 0.8);
  DrawAndSave2D(residuals, fout, "residual_rphi_vs_cluster_r_negative_side1_eta_lt_0p5",
                "residual_rphi:cluster_r", cutResidualRNegativeSide1EtaLt,
                label + ", charge < 0, side == 1, |eta| < 0.5", outdir,
                "cluster r [cm]", "r#Delta#phi residual [cm]", 55, 20, 75, 200, -0.8, 0.8);
  DrawAndSave2D(residuals, fout, "residual_rphi_vs_cluster_r_positive_side0_eta_gt_0p5",
                "residual_rphi:cluster_r", cutResidualRPositiveSide0EtaGt,
                label + ", charge > 0, side == 0, |eta| > 0.5", outdir,
                "cluster r [cm]", "r#Delta#phi residual [cm]", 55, 20, 75, 200, -0.8, 0.8);
  DrawAndSave2D(residuals, fout, "residual_rphi_vs_cluster_r_positive_side1_eta_gt_0p5",
                "residual_rphi:cluster_r", cutResidualRPositiveSide1EtaGt,
                label + ", charge > 0, side == 1, |eta| > 0.5", outdir,
                "cluster r [cm]", "r#Delta#phi residual [cm]", 55, 20, 75, 200, -0.8, 0.8);
  DrawAndSave2D(residuals, fout, "residual_rphi_vs_cluster_r_negative_side0_eta_gt_0p5",
                "residual_rphi:cluster_r", cutResidualRNegativeSide0EtaGt,
                label + ", charge < 0, side == 0, |eta| > 0.5", outdir,
                "cluster r [cm]", "r#Delta#phi residual [cm]", 55, 20, 75, 200, -0.8, 0.8);
  DrawAndSave2D(residuals, fout, "residual_rphi_vs_cluster_r_negative_side1_eta_gt_0p5",
                "residual_rphi:cluster_r", cutResidualRNegativeSide1EtaGt,
                label + ", charge < 0, side == 1, |eta| > 0.5", outdir,
                "cluster r [cm]", "r#Delta#phi residual [cm]", 55, 20, 75, 200, -0.8, 0.8);

  // residual_z vs cluster_z, both charges
  TH2F* h_residual_z_vs_cluster_z = new TH2F("h_residual_z_vs_cluster_z", "",
                                             204, -102, 102, 200, -1, 1);
  Style2D(h_residual_z_vs_cluster_z, "cluster z [cm]", "#Deltaz residual [cm]");
  residuals.Draw("residual_z:cluster_z>>h_residual_z_vs_cluster_z", cutResidualR.c_str(), "goff");
  TCanvas* c_residual_z_vs_cluster_z = new TCanvas("c_residual_z_vs_cluster_z", "",
                                                   kCanvasWidth, kCanvasHeight);
  c_residual_z_vs_cluster_z->SetRightMargin(0.16);
  c_residual_z_vs_cluster_z->SetTopMargin(0.10);
  h_residual_z_vs_cluster_z->Draw("COLZ");
  DrawLabel(label, cutResidualR);
  SaveCanvasAndObject(fout, c_residual_z_vs_cluster_z, h_residual_z_vs_cluster_z,
                      outdir, "residual_z_vs_cluster_z");

  // residual_z vs cluster_z, positive charge
  TH2F* h_residual_z_vs_cluster_z_positive =
      new TH2F("h_residual_z_vs_cluster_z_positive", "", 204, -102, 102, 200, -10, 10);
  Style2D(h_residual_z_vs_cluster_z_positive, "cluster z [cm]", "#Deltaz residual [cm]");
  residuals.Draw("residual_z:cluster_z>>h_residual_z_vs_cluster_z_positive",
                 cutResidualRPositive.c_str(), "goff");
  TCanvas* c_residual_z_vs_cluster_z_positive =
      new TCanvas("c_residual_z_vs_cluster_z_positive", "", kCanvasWidth, kCanvasHeight);
  c_residual_z_vs_cluster_z_positive->SetRightMargin(0.16);
  c_residual_z_vs_cluster_z_positive->SetTopMargin(0.10);
  h_residual_z_vs_cluster_z_positive->Draw("COLZ");
  DrawLabel(label + ", charge > 0", cutResidualRPositive);
  SaveCanvasAndObject(fout, c_residual_z_vs_cluster_z_positive,
                      h_residual_z_vs_cluster_z_positive, outdir,
                      "residual_z_vs_cluster_z_positive");

  // residual_z vs cluster_z, negative charge
  TH2F* h_residual_z_vs_cluster_z_negative =
      new TH2F("h_residual_z_vs_cluster_z_negative", "", 204, -102, 102, 200, -10, 10);
  Style2D(h_residual_z_vs_cluster_z_negative, "cluster z [cm]", "#Deltaz residual [cm]");
  residuals.Draw("residual_z:cluster_z>>h_residual_z_vs_cluster_z_negative",
                 cutResidualRNegative.c_str(), "goff");
  TCanvas* c_residual_z_vs_cluster_z_negative =
      new TCanvas("c_residual_z_vs_cluster_z_negative", "", kCanvasWidth, kCanvasHeight);
  c_residual_z_vs_cluster_z_negative->SetRightMargin(0.16);
  c_residual_z_vs_cluster_z_negative->SetTopMargin(0.10);
  h_residual_z_vs_cluster_z_negative->Draw("COLZ");
  DrawLabel(label + ", charge < 0", cutResidualRNegative);
  SaveCanvasAndObject(fout, c_residual_z_vs_cluster_z_negative,
                      h_residual_z_vs_cluster_z_negative, outdir,
                      "residual_z_vs_cluster_z_negative");

  // residual_rphi vs pt
  TH2F* h_residual_rphi_vs_pt = new TH2F("h_residual_rphi_vs_pt", "",
                                         55, 0, 5, 200, -0.8, 0.8);
  Style2D(h_residual_rphi_vs_pt, "p_{T} [GeV]", "r#Delta#phi residual [cm]");
  residuals.Draw("residual_rphi:pt>>h_residual_rphi_vs_pt", cutResidualPt.c_str(), "goff");
  TCanvas* c_residual_rphi_vs_pt = new TCanvas("c_residual_rphi_vs_pt", "",
                                               kCanvasWidth, kCanvasHeight);
  c_residual_rphi_vs_pt->SetRightMargin(0.16);
  c_residual_rphi_vs_pt->SetTopMargin(0.10);
  h_residual_rphi_vs_pt->Draw("COLZ");
  DrawLabel(label, cutResidualPt);
  SaveCanvasAndObject(fout, c_residual_rphi_vs_pt, h_residual_rphi_vs_pt,
                      outdir, "residual_rphi_vs_pt");

  // vertex xy, all accepted tracks
  TH2F* h_vertex_xy = new TH2F("h_vertex_xy", "", 100, -10, 10, 100, -10, 10);
  Style2D(h_vertex_xy, "vertex y [cm]", "vertex x [cm]");
  residuals.Draw("vertex_x:vertex_y>>h_vertex_xy", cutVertexAll.c_str(), "goff");
  TCanvas* c_vertex_xy = new TCanvas("c_vertex_xy", "", kCanvasWidth, kCanvasHeight);
  c_vertex_xy->SetRightMargin(0.16);
  c_vertex_xy->SetTopMargin(0.10);
  h_vertex_xy->Draw("COLZ");
  DrawLabel(label, cutVertexAll);
  SaveCanvasAndObject(fout, c_vertex_xy, h_vertex_xy, outdir, "vertex_xy");

  // vertex xy, pt selected
  TH2F* h_vertex_xy_pt = new TH2F("h_vertex_xy_pt", "", 100, -10, 10, 100, -10, 10);
  Style2D(h_vertex_xy_pt, "vertex y [cm]", "vertex x [cm]");
  residuals.Draw("vertex_x:vertex_y>>h_vertex_xy_pt", cutVertexPt.c_str(), "goff");
  TCanvas* c_vertex_xy_pt = new TCanvas("c_vertex_xy_pt", "", kCanvasWidth, kCanvasHeight);
  c_vertex_xy_pt->SetRightMargin(0.16);
  c_vertex_xy_pt->SetTopMargin(0.10);
  h_vertex_xy_pt->Draw("COLZ");
  DrawLabel(label, cutVertexPt);
  SaveCanvasAndObject(fout, c_vertex_xy_pt, h_vertex_xy_pt, outdir, "vertex_xy_pt_gt_2");

  // quality
  TH1F* h_quality = new TH1F("h_quality", "", 100, 0, 0.5);
  StyleHist(h_quality, "quality", "Tracks");
  h_quality->SetFillColorAlpha(kAzure + 2, 0.30);
  residuals.Draw("quality>>h_quality", cutQuality.c_str(), "goff");
  TCanvas* c_quality = new TCanvas("c_quality", "", kCanvasWidth, kCanvasHeight);
  c_quality->SetTopMargin(0.10);
  h_quality->Draw("HIST");
  DrawLabel(label, cutQuality);
  SaveCanvasAndObject(fout, c_quality, h_quality, outdir, "quality");

  // rDCA vs charge/pt
  TH2F* h_rDCA_vs_chargeOverPt = new TH2F("h_rDCA_vs_chargeOverPt", "",
                                          100, -5, 5, 200, -10, 10);
  Style2D(h_rDCA_vs_chargeOverPt, "charge / p_{T} [1/GeV]", "rDCA [cm]");
  residuals.Draw("rDCA:charge/pt>>h_rDCA_vs_chargeOverPt", cutDca.c_str(), "goff");
  TCanvas* c_rDCA_vs_chargeOverPt = new TCanvas("c_rDCA_vs_chargeOverPt", "",
                                                kCanvasWidth, kCanvasHeight);
  c_rDCA_vs_chargeOverPt->SetRightMargin(0.16);
  c_rDCA_vs_chargeOverPt->SetTopMargin(0.10);
  h_rDCA_vs_chargeOverPt->Draw("COLZ");
  DrawLabel(label, cutDca);
  SaveCanvasAndObject(fout, c_rDCA_vs_chargeOverPt, h_rDCA_vs_chargeOverPt,
                      outdir, "rDCA_vs_chargeOverPt");

  DrawAndSave2D(residuals, fout, "rDCA_vs_chargeOverPt_side0",
                "rDCA:charge/pt", cutDcaSide0, label + ", side == 0", outdir,
                "charge / p_{T} [1/GeV]", "rDCA [cm]", 100, -5, 5, 200, -10, 10);
  DrawAndSave2D(residuals, fout, "rDCA_vs_chargeOverPt_side1",
                "rDCA:charge/pt", cutDcaSide1, label + ", side == 1", outdir,
                "charge / p_{T} [1/GeV]", "rDCA [cm]", 100, -5, 5, 200, -10, 10);

  // rDCA_zero vs charge/pt
  TH2F* h_rDCA_zero_vs_chargeOverPt = new TH2F("h_rDCA_zero_vs_chargeOverPt", "",
                                               100, -5, 5, 200, -10, 10);
  Style2D(h_rDCA_zero_vs_chargeOverPt, "charge / p_{T} [1/GeV]", "rDCA_zero [cm]");
  residuals.Draw("rDCA_zero:charge/pt>>h_rDCA_zero_vs_chargeOverPt", cutDca.c_str(), "goff");
  TCanvas* c_rDCA_zero_vs_chargeOverPt = new TCanvas("c_rDCA_zero_vs_chargeOverPt", "",
                                                     kCanvasWidth, kCanvasHeight);
  c_rDCA_zero_vs_chargeOverPt->SetRightMargin(0.16);
  c_rDCA_zero_vs_chargeOverPt->SetTopMargin(0.10);
  h_rDCA_zero_vs_chargeOverPt->Draw("COLZ");
  DrawLabel(label, cutDca);
  SaveCanvasAndObject(fout, c_rDCA_zero_vs_chargeOverPt, h_rDCA_zero_vs_chargeOverPt,
                      outdir, "rDCA_zero_vs_chargeOverPt");

  DrawAndSave2D(residuals, fout, "rDCA_zero_vs_chargeOverPt_side0",
                "rDCA_zero:charge/pt", cutDcaSide0, label + ", side == 0", outdir,
                "charge / p_{T} [1/GeV]", "rDCA_zero [cm]", 100, -5, 5, 200, -10, 10);
  DrawAndSave2D(residuals, fout, "rDCA_zero_vs_chargeOverPt_side1",
                "rDCA_zero:charge/pt", cutDcaSide1, label + ", side == 1", outdir,
                "charge / p_{T} [1/GeV]", "rDCA_zero [cm]", 100, -5, 5, 200, -10, 10);

  fout->Close();
  delete fout;

  std::cout << "Wrote ROOT file: " << outfile << std::endl;
  std::cout << "Wrote PNG files in: " << outdir << std::endl;
}
