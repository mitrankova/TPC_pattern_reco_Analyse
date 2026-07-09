#include "TCanvas.h"
#include "TChain.h"
#include "TColor.h"
#include "TError.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TLeaf.h"
#include "TObject.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

enum class BranchType
{
  kFloat,
  kDouble,
  kInt,
  kUInt,
  kShort,
  kUShort,
  kLong64,
  kULong64
};

struct ScalarBranch
{
  std::string name;
  BranchType type = BranchType::kFloat;
  Float_t f = 0;
  Double_t d = 0;
  Int_t i = 0;
  UInt_t ui = 0;
  Short_t s = 0;
  UShort_t us = 0;
  Long64_t l = 0;
  ULong64_t ul = 0;

  bool Bind(TChain& tree, const char* branchName)
  {
    name = branchName;
    TLeaf* leaf = tree.GetLeaf(branchName);
    if (!leaf)
    {
      std::cerr << "Missing branch: " << branchName << std::endl;
      return false;
    }

    const std::string typeName = leaf->GetTypeName();
    void* address = &f;
    if (typeName == "Float_t" || typeName == "float")
    {
      type = BranchType::kFloat;
      address = &f;
    }
    else if (typeName == "Double_t" || typeName == "double")
    {
      type = BranchType::kDouble;
      address = &d;
    }
    else if (typeName == "Int_t" || typeName == "int")
    {
      type = BranchType::kInt;
      address = &i;
    }
    else if (typeName == "UInt_t" || typeName == "unsigned int")
    {
      type = BranchType::kUInt;
      address = &ui;
    }
    else if (typeName == "Short_t" || typeName == "short")
    {
      type = BranchType::kShort;
      address = &s;
    }
    else if (typeName == "UShort_t" || typeName == "unsigned short")
    {
      type = BranchType::kUShort;
      address = &us;
    }
    else if (typeName == "Long64_t" || typeName == "long long")
    {
      type = BranchType::kLong64;
      address = &l;
    }
    else if (typeName == "ULong64_t" || typeName == "unsigned long long")
    {
      type = BranchType::kULong64;
      address = &ul;
    }
    else
    {
      std::cerr << "Unsupported branch type for " << branchName << ": "
                << typeName << std::endl;
      return false;
    }

    return tree.SetBranchAddress(branchName, address) >= 0;
  }

  double Value() const
  {
    switch (type)
    {
    case BranchType::kFloat:
      return f;
    case BranchType::kDouble:
      return d;
    case BranchType::kInt:
      return i;
    case BranchType::kUInt:
      return ui;
    case BranchType::kShort:
      return s;
    case BranchType::kUShort:
      return us;
    case BranchType::kLong64:
      return l;
    case BranchType::kULong64:
      return ul;
    }
    return 0;
  }
};

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

void PrintEntryProgress(Long64_t current, Long64_t total)
{
  if (total <= 0)
  {
    return;
  }

  const int barWidth = 50;
  const double fraction = static_cast<double>(current) / static_cast<double>(total);
  const int filled = static_cast<int>(barWidth * fraction);
  const int percent = static_cast<int>(100.0 * fraction);

  std::cout << "\rReading entries [";
  for (int i = 0; i < barWidth; ++i)
  {
    if (i < filled)
    {
      std::cout << "=";
    }
    else if (i == filled && current < total)
    {
      std::cout << ">";
    }
    else
    {
      std::cout << " ";
    }
  }
  std::cout << "] " << percent << "% (" << current << "/" << total << ")" << std::flush;

  if (current >= total)
  {
    std::cout << std::endl;
  }
}

void DrawAndSaveHist(TFile* fout, TH1* hist,
                     const std::string& plotLabel,
                     const std::string& cut,
                     const std::string& outdir,
                     const std::string& name,
                     const char* drawOption,
                     bool setRightMargin)
{
  const std::string canvasName = "c_" + name;
  TCanvas* canvas = new TCanvas(canvasName.c_str(), "", kCanvasWidth, kCanvasHeight);
  if (setRightMargin)
  {
    canvas->SetRightMargin(0.16);
  }
  canvas->SetTopMargin(0.10);
  hist->Draw(drawOption);
  DrawLabel(plotLabel, cut);
  SaveCanvasAndObject(fout, canvas, hist, outdir, name);
}
// .x PlotResidualQuickQA.C("/sphenix/tg/tg01/hf/mitrankova/PatternReco/79513_keff_m1/outout_resid/tpc_poly_cluster_residualsHITS_ppFieldOn_79513*.root","79513_keff_m1")
// .x PlotResidualQuickQA.C("/sphenix/tg/tg01/hf/mitrankova/PatternReco/79513_keff_1/outout_resid/tpc_poly_cluster_residualsHITS_ppFieldOn_79513*.root","79513_keff_1")
// .x PlotResidualQuickQA.C("/sphenix/tg/tg01/hf/mitrankova/PatternReco/79513_keff_2/outout_resid/tpc_poly_cluster_residualsHITS_ppFieldOn_79513*.root","79513_keff_2")
// .x PlotResidualQuickQA.C("/sphenix/tg/tg01/hf/mitrankova/PatternReco/79513_keff_4/outout_resid/tpc_poly_cluster_residualsHITS_ppFieldOn_79513*.root","79513_keff_4")
// .x PlotResidualQuickQA.C("/sphenix/tg/tg01/hf/mitrankova/PatternReco/79513/outout_resid/tpc_poly_cluster_residualsHITS_ppFieldOn_79513*.root","79513")
// .x PlotResidualQuickQA.C("/sphenix/tg/tg01/hf/mitrankova/PatternReco/76905_6t0/outout_resid/tpc_poly_cluster_residualsHITS_ppFieldOn_769051*.root","76905_6t0","auau_lowInt")
// .x PlotResidualQuickQA.C("/sphenix/tg/tg01/hf/mitrankova/PatternReco/76905_85t0/outout_resid/tpc_poly_cluster_residualsHITS_ppFieldOn_769051*.root","76905_85t0","auau_lowInt")
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
  const std::string dedxCut = "ntpc_clusters>20";
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

  // Book histograms once, then fill all of them in a single pass over the tree.
  TH2F* h_residual_rphi_vs_cluster_r = new TH2F("h_residual_rphi_vs_cluster_r", "",
                                                55, 20, 75, 200, -0.8, 0.8);
  Style2D(h_residual_rphi_vs_cluster_r, "cluster r [cm]", "r#Delta#phi residual [cm]");
  TH2F* h_residual_rphi_vs_cluster_r_positive =
      new TH2F("h_residual_rphi_vs_cluster_r_positive", "", 55, 20, 75, 200, -0.8, 0.8);
  Style2D(h_residual_rphi_vs_cluster_r_positive, "cluster r [cm]", "r#Delta#phi residual [cm]");
  TH2F* h_residual_rphi_vs_cluster_r_negative =
      new TH2F("h_residual_rphi_vs_cluster_r_negative", "", 55, 20, 75, 200, -0.8, 0.8);
  Style2D(h_residual_rphi_vs_cluster_r_negative, "cluster r [cm]", "r#Delta#phi residual [cm]");

  TH2F* h_residual_rphi_vs_cluster_r_side0 =
      new TH2F("residual_rphi_vs_cluster_r_side0", "", 55, 20, 75, 200, -0.8, 0.8);
  Style2D(h_residual_rphi_vs_cluster_r_side0, "cluster r [cm]", "r#Delta#phi residual [cm]");
  TH2F* h_residual_rphi_vs_cluster_r_side1 =
      new TH2F("residual_rphi_vs_cluster_r_side1", "", 55, 20, 75, 200, -0.8, 0.8);
  Style2D(h_residual_rphi_vs_cluster_r_side1, "cluster r [cm]", "r#Delta#phi residual [cm]");
  TH2F* h_residual_rphi_vs_cluster_r_positive_side0 =
      new TH2F("residual_rphi_vs_cluster_r_positive_side0", "", 55, 20, 75, 200, -0.8, 0.8);
  Style2D(h_residual_rphi_vs_cluster_r_positive_side0, "cluster r [cm]", "r#Delta#phi residual [cm]");
  TH2F* h_residual_rphi_vs_cluster_r_positive_side1 =
      new TH2F("residual_rphi_vs_cluster_r_positive_side1", "", 55, 20, 75, 200, -0.8, 0.8);
  Style2D(h_residual_rphi_vs_cluster_r_positive_side1, "cluster r [cm]", "r#Delta#phi residual [cm]");
  TH2F* h_residual_rphi_vs_cluster_r_negative_side0 =
      new TH2F("residual_rphi_vs_cluster_r_negative_side0", "", 55, 20, 75, 200, -0.8, 0.8);
  Style2D(h_residual_rphi_vs_cluster_r_negative_side0, "cluster r [cm]", "r#Delta#phi residual [cm]");
  TH2F* h_residual_rphi_vs_cluster_r_negative_side1 =
      new TH2F("residual_rphi_vs_cluster_r_negative_side1", "", 55, 20, 75, 200, -0.8, 0.8);
  Style2D(h_residual_rphi_vs_cluster_r_negative_side1, "cluster r [cm]", "r#Delta#phi residual [cm]");

  TH2F* h_residual_rphi_vs_cluster_r_positive_side0_eta_lt_0p5 =
      new TH2F("residual_rphi_vs_cluster_r_positive_side0_eta_lt_0p5", "", 55, 20, 75, 200, -0.8, 0.8);
  Style2D(h_residual_rphi_vs_cluster_r_positive_side0_eta_lt_0p5, "cluster r [cm]", "r#Delta#phi residual [cm]");
  TH2F* h_residual_rphi_vs_cluster_r_positive_side1_eta_lt_0p5 =
      new TH2F("residual_rphi_vs_cluster_r_positive_side1_eta_lt_0p5", "", 55, 20, 75, 200, -0.8, 0.8);
  Style2D(h_residual_rphi_vs_cluster_r_positive_side1_eta_lt_0p5, "cluster r [cm]", "r#Delta#phi residual [cm]");
  TH2F* h_residual_rphi_vs_cluster_r_negative_side0_eta_lt_0p5 =
      new TH2F("residual_rphi_vs_cluster_r_negative_side0_eta_lt_0p5", "", 55, 20, 75, 200, -0.8, 0.8);
  Style2D(h_residual_rphi_vs_cluster_r_negative_side0_eta_lt_0p5, "cluster r [cm]", "r#Delta#phi residual [cm]");
  TH2F* h_residual_rphi_vs_cluster_r_negative_side1_eta_lt_0p5 =
      new TH2F("residual_rphi_vs_cluster_r_negative_side1_eta_lt_0p5", "", 55, 20, 75, 200, -0.8, 0.8);
  Style2D(h_residual_rphi_vs_cluster_r_negative_side1_eta_lt_0p5, "cluster r [cm]", "r#Delta#phi residual [cm]");
  TH2F* h_residual_rphi_vs_cluster_r_positive_side0_eta_gt_0p5 =
      new TH2F("residual_rphi_vs_cluster_r_positive_side0_eta_gt_0p5", "", 55, 20, 75, 200, -0.8, 0.8);
  Style2D(h_residual_rphi_vs_cluster_r_positive_side0_eta_gt_0p5, "cluster r [cm]", "r#Delta#phi residual [cm]");
  TH2F* h_residual_rphi_vs_cluster_r_positive_side1_eta_gt_0p5 =
      new TH2F("residual_rphi_vs_cluster_r_positive_side1_eta_gt_0p5", "", 55, 20, 75, 200, -0.8, 0.8);
  Style2D(h_residual_rphi_vs_cluster_r_positive_side1_eta_gt_0p5, "cluster r [cm]", "r#Delta#phi residual [cm]");
  TH2F* h_residual_rphi_vs_cluster_r_negative_side0_eta_gt_0p5 =
      new TH2F("residual_rphi_vs_cluster_r_negative_side0_eta_gt_0p5", "", 55, 20, 75, 200, -0.8, 0.8);
  Style2D(h_residual_rphi_vs_cluster_r_negative_side0_eta_gt_0p5, "cluster r [cm]", "r#Delta#phi residual [cm]");
  TH2F* h_residual_rphi_vs_cluster_r_negative_side1_eta_gt_0p5 =
      new TH2F("residual_rphi_vs_cluster_r_negative_side1_eta_gt_0p5", "", 55, 20, 75, 200, -0.8, 0.8);
  Style2D(h_residual_rphi_vs_cluster_r_negative_side1_eta_gt_0p5, "cluster r [cm]", "r#Delta#phi residual [cm]");

  TH2F* h_residual_z_vs_cluster_z = new TH2F("h_residual_z_vs_cluster_z", "",
                                             204, -102, 102, 200, -1, 1);
  Style2D(h_residual_z_vs_cluster_z, "cluster z [cm]", "#Deltaz residual [cm]");
  TH2F* h_residual_z_vs_cluster_z_positive =
      new TH2F("h_residual_z_vs_cluster_z_positive", "", 204, -102, 102, 200, -10, 10);
  Style2D(h_residual_z_vs_cluster_z_positive, "cluster z [cm]", "#Deltaz residual [cm]");
  TH2F* h_residual_z_vs_cluster_z_negative =
      new TH2F("h_residual_z_vs_cluster_z_negative", "", 204, -102, 102, 200, -10, 10);
  Style2D(h_residual_z_vs_cluster_z_negative, "cluster z [cm]", "#Deltaz residual [cm]");

  TH2F* h_residual_rphi_vs_pt = new TH2F("h_residual_rphi_vs_pt", "",
                                         55, 0, 5, 200, -0.8, 0.8);
  Style2D(h_residual_rphi_vs_pt, "p_{T} [GeV]", "r#Delta#phi residual [cm]");

  TH2F* h_dedx_vs_charge_momentum = new TH2F("h_dedx_vs_charge_momentum", "",
                                                200, -3, 3, 100, 0, 1000);
  Style2D(h_dedx_vs_charge_momentum, "charge #times p [GeV]", "dE/dx");

  TH2F* h_vertex_xy = new TH2F("h_vertex_xy", "", 100, -10, 10, 100, -10, 10);
  Style2D(h_vertex_xy, "vertex y [cm]", "vertex x [cm]");
  TH2F* h_vertex_xy_pt = new TH2F("h_vertex_xy_pt", "", 100, -10, 10, 100, -10, 10);
  Style2D(h_vertex_xy_pt, "vertex y [cm]", "vertex x [cm]");

  TH1F* h_quality = new TH1F("h_quality", "", 100, 0, 0.5);
  StyleHist(h_quality, "quality", "Tracks");
  h_quality->SetFillColorAlpha(kAzure + 2, 0.30);

  TH2F* h_rDCA_vs_chargeOverPt = new TH2F("h_rDCA_vs_chargeOverPt", "",
                                          100, -5, 5, 200, -10, 10);
  Style2D(h_rDCA_vs_chargeOverPt, "charge / p_{T} [1/GeV]", "rDCA [cm]");
  TH2F* h_rDCA_vs_chargeOverPt_side0 = new TH2F("rDCA_vs_chargeOverPt_side0", "",
                                                100, -5, 5, 200, -10, 10);
  Style2D(h_rDCA_vs_chargeOverPt_side0, "charge / p_{T} [1/GeV]", "rDCA [cm]");
  TH2F* h_rDCA_vs_chargeOverPt_side1 = new TH2F("rDCA_vs_chargeOverPt_side1", "",
                                                100, -5, 5, 200, -10, 10);
  Style2D(h_rDCA_vs_chargeOverPt_side1, "charge / p_{T} [1/GeV]", "rDCA [cm]");

  TH2F* h_rDCA_zero_vs_chargeOverPt = new TH2F("h_rDCA_zero_vs_chargeOverPt", "",
                                               100, -5, 5, 200, -10, 10);
  Style2D(h_rDCA_zero_vs_chargeOverPt, "charge / p_{T} [1/GeV]", "rDCA_zero [cm]");
  TH2F* h_rDCA_zero_vs_chargeOverPt_side0 = new TH2F("rDCA_zero_vs_chargeOverPt_side0", "",
                                                     100, -5, 5, 200, -10, 10);
  Style2D(h_rDCA_zero_vs_chargeOverPt_side0, "charge / p_{T} [1/GeV]", "rDCA_zero [cm]");
  TH2F* h_rDCA_zero_vs_chargeOverPt_side1 = new TH2F("rDCA_zero_vs_chargeOverPt_side1", "",
                                                     100, -5, 5, 200, -10, 10);
  Style2D(h_rDCA_zero_vs_chargeOverPt_side1, "charge / p_{T} [1/GeV]", "rDCA_zero [cm]");

  residuals.SetBranchStatus("*", 0);
  residuals.SetBranchStatus("residual_rphi", 1);
  residuals.SetBranchStatus("residual_z", 1);
  residuals.SetBranchStatus("cluster_r", 1);
  residuals.SetBranchStatus("cluster_z", 1);
  residuals.SetBranchStatus("pt", 1);
  residuals.SetBranchStatus("pz", 1);
  residuals.SetBranchStatus("py", 1);
  residuals.SetBranchStatus("px", 1);
  residuals.SetBranchStatus("dedx", 1);
  residuals.SetBranchStatus("eta", 1);
  residuals.SetBranchStatus("vertex_z", 1);
  residuals.SetBranchStatus("vertex_x", 1);
  residuals.SetBranchStatus("vertex_y", 1);
  residuals.SetBranchStatus("quality", 1);
  residuals.SetBranchStatus("rDCA", 1);
  residuals.SetBranchStatus("rDCA_zero", 1);
  residuals.SetBranchStatus("ntpc_clusters", 1);
  residuals.SetBranchStatus("charge", 1);
  residuals.SetBranchStatus("side", 1);

  ScalarBranch b_residual_rphi;
  ScalarBranch b_residual_z;
  ScalarBranch b_cluster_r;
  ScalarBranch b_cluster_z;
  ScalarBranch b_pt;
  ScalarBranch b_dedx;
  ScalarBranch b_px;
  ScalarBranch b_py;
  ScalarBranch b_pz;
  ScalarBranch b_eta;
  ScalarBranch b_vertex_z;
  ScalarBranch b_vertex_x;
  ScalarBranch b_vertex_y;
  ScalarBranch b_quality;
  ScalarBranch b_rDCA;
  ScalarBranch b_rDCA_zero;
  ScalarBranch b_ntpc_clusters;
  ScalarBranch b_charge;
  ScalarBranch b_side;

  const bool branchesOk =
      b_residual_rphi.Bind(residuals, "residual_rphi") &&
      b_residual_z.Bind(residuals, "residual_z") &&
      b_cluster_r.Bind(residuals, "cluster_r") &&
      b_cluster_z.Bind(residuals, "cluster_z") &&
      b_pt.Bind(residuals, "pt") &&
      b_dedx.Bind(residuals, "dedx") &&
      b_px.Bind(residuals, "px") &&
      b_py.Bind(residuals, "py") &&
      b_pz.Bind(residuals, "pz") &&
      b_eta.Bind(residuals, "eta") &&
      b_vertex_z.Bind(residuals, "vertex_z") &&
      b_vertex_x.Bind(residuals, "vertex_x") &&
      b_vertex_y.Bind(residuals, "vertex_y") &&
      b_quality.Bind(residuals, "quality") &&
      b_rDCA.Bind(residuals, "rDCA") &&
      b_rDCA_zero.Bind(residuals, "rDCA_zero") &&
      b_ntpc_clusters.Bind(residuals, "ntpc_clusters") &&
      b_charge.Bind(residuals, "charge") &&
      b_side.Bind(residuals, "side");

  if (!branchesOk)
  {
    std::cerr << "Input tree is missing at least one required branch or has an unsupported branch type." << std::endl;
    fout->Close();
    delete fout;
    return;
  }

  const Long64_t nEntries = residuals.GetEntries();
  const Long64_t progressEvery = nEntries > 1000 ? nEntries / 1000 : 1;
  PrintEntryProgress(0, nEntries);

  for (Long64_t ientry = 0; ientry < nEntries; ++ientry)
  {
    residuals.GetEntry(ientry);

    if ((ientry + 1) % progressEvery == 0 || ientry + 1 == nEntries)
    {
      PrintEntryProgress(ientry + 1, nEntries);
    }

    const double residual_rphi = b_residual_rphi.Value();
    const double residual_z = b_residual_z.Value();
    const double cluster_r = b_cluster_r.Value();
    const double cluster_z = b_cluster_z.Value();
    const double pt = b_pt.Value();
    const double dedx = b_dedx.Value();
    const double px = b_px.Value();
    const double py = b_py.Value();
    const double pz = b_pz.Value();
    const double eta = b_eta.Value();
    const double vertex_z = b_vertex_z.Value();
    const double vertex_x = b_vertex_x.Value();
    const double vertex_y = b_vertex_y.Value();
    const double quality = b_quality.Value();
    const double rDCA = b_rDCA.Value();
    const double rDCA_zero = b_rDCA_zero.Value();
    const int ntpc_clusters = b_ntpc_clusters.Value();
    const int charge = b_charge.Value();
    const int side = b_side.Value();

    const double absEta = std::abs(eta);
    const bool passVz = vertex_z > kVzMin && vertex_z < kVzMax;
    const bool passVzDca = vertex_z > kVzMinDCA && vertex_z < kVzMaxDCA;
    const bool passTpcClusters = ntpc_clusters > kMinTpcClusters;
    const bool passResidualR = passVz && pt > kPtMinResidualVsR && absEta < kEtaAbsMaxResidualVsR;
    const bool passResidualEtaLt = passVz && pt > kPtMinResidualEtaSplit && absEta < kEtaAbsSplit;
    const bool passResidualEtaGt = passVz && pt > kPtMinResidualEtaSplit && absEta > kEtaAbsSplit;
    const bool passResidualPt = passVz;
    const bool passDedx = ntpc_clusters > 20 && pt > 0.5;
    const bool passVertexAll = passTpcClusters;
    const bool passVertexPt = passTpcClusters && pt > kPtMinResidualVsR;
    const bool passQuality = passTpcClusters && pt > kQualityPtMin;
    const bool passDca = passTpcClusters && passVzDca;
    const bool positiveCharge = charge > 0;
    const bool negativeCharge = charge < 0;
    const bool side0 = side == 0;
    const bool side1 = side == 1;

    if (passResidualR)
    {
      h_residual_rphi_vs_cluster_r->Fill(cluster_r, residual_rphi);
      h_residual_z_vs_cluster_z->Fill(cluster_z, residual_z);
      if (positiveCharge)
      {
        h_residual_rphi_vs_cluster_r_positive->Fill(cluster_r, residual_rphi);
        h_residual_z_vs_cluster_z_positive->Fill(cluster_z, residual_z);
      }
      if (negativeCharge)
      {
        h_residual_rphi_vs_cluster_r_negative->Fill(cluster_r, residual_rphi);
        h_residual_z_vs_cluster_z_negative->Fill(cluster_z, residual_z);
      }
      if (side0)
      {
        h_residual_rphi_vs_cluster_r_side0->Fill(cluster_r, residual_rphi);
      }
      if (side1)
      {
        h_residual_rphi_vs_cluster_r_side1->Fill(cluster_r, residual_rphi);
      }
      if (positiveCharge && side0)
      {
        h_residual_rphi_vs_cluster_r_positive_side0->Fill(cluster_r, residual_rphi);
      }
      if (positiveCharge && side1)
      {
        h_residual_rphi_vs_cluster_r_positive_side1->Fill(cluster_r, residual_rphi);
      }
      if (negativeCharge && side0)
      {
        h_residual_rphi_vs_cluster_r_negative_side0->Fill(cluster_r, residual_rphi);
      }
      if (negativeCharge && side1)
      {
        h_residual_rphi_vs_cluster_r_negative_side1->Fill(cluster_r, residual_rphi);
      }
    }

    if (passResidualEtaLt)
    {
      if (positiveCharge && side0)
      {
        h_residual_rphi_vs_cluster_r_positive_side0_eta_lt_0p5->Fill(cluster_r, residual_rphi);
      }
      if (positiveCharge && side1)
      {
        h_residual_rphi_vs_cluster_r_positive_side1_eta_lt_0p5->Fill(cluster_r, residual_rphi);
      }
      if (negativeCharge && side0)
      {
        h_residual_rphi_vs_cluster_r_negative_side0_eta_lt_0p5->Fill(cluster_r, residual_rphi);
      }
      if (negativeCharge && side1)
      {
        h_residual_rphi_vs_cluster_r_negative_side1_eta_lt_0p5->Fill(cluster_r, residual_rphi);
      }
    }

    if (passResidualEtaGt)
    {
      if (positiveCharge && side0)
      {
        h_residual_rphi_vs_cluster_r_positive_side0_eta_gt_0p5->Fill(cluster_r, residual_rphi);
      }
      if (positiveCharge && side1)
      {
        h_residual_rphi_vs_cluster_r_positive_side1_eta_gt_0p5->Fill(cluster_r, residual_rphi);
      }
      if (negativeCharge && side0)
      {
        h_residual_rphi_vs_cluster_r_negative_side0_eta_gt_0p5->Fill(cluster_r, residual_rphi);
      }
      if (negativeCharge && side1)
      {
        h_residual_rphi_vs_cluster_r_negative_side1_eta_gt_0p5->Fill(cluster_r, residual_rphi);
      }
    }

    if (passResidualPt)
    {
      h_residual_rphi_vs_pt->Fill(pt, residual_rphi);
    }
    if (passDedx)
    {
      const double chargeMomentum = charge * std::sqrt(px * px + py * py + pz * pz);
      h_dedx_vs_charge_momentum->Fill(chargeMomentum, dedx);
    }
    if (passVertexAll)
    {
      h_vertex_xy->Fill(vertex_y, vertex_x);
    }
    if (passVertexPt)
    {
      h_vertex_xy_pt->Fill(vertex_y, vertex_x);
    }
    if (passQuality)
    {
      h_quality->Fill(quality);
    }
    if (passDca && pt != 0)
    {
      const double chargeOverPt = charge / pt;
      h_rDCA_vs_chargeOverPt->Fill(chargeOverPt, rDCA);
      h_rDCA_zero_vs_chargeOverPt->Fill(chargeOverPt, rDCA_zero);
      if (side0)
      {
        h_rDCA_vs_chargeOverPt_side0->Fill(chargeOverPt, rDCA);
        h_rDCA_zero_vs_chargeOverPt_side0->Fill(chargeOverPt, rDCA_zero);
      }
      if (side1)
      {
        h_rDCA_vs_chargeOverPt_side1->Fill(chargeOverPt, rDCA);
        h_rDCA_zero_vs_chargeOverPt_side1->Fill(chargeOverPt, rDCA_zero);
      }
    }
  }

  DrawAndSaveHist(fout, h_residual_rphi_vs_cluster_r, label, cutResidualR,
                  outdir, "residual_rphi_vs_cluster_r", "COLZ", true);
  DrawAndSaveHist(fout, h_residual_rphi_vs_cluster_r_positive, label + ", charge > 0", cutResidualRPositive,
                  outdir, "residual_rphi_vs_cluster_r_positive", "COLZ", true);
  DrawAndSaveHist(fout, h_residual_rphi_vs_cluster_r_negative, label + ", charge < 0", cutResidualRNegative,
                  outdir, "residual_rphi_vs_cluster_r_negative", "COLZ", true);
  DrawAndSaveHist(fout, h_residual_rphi_vs_cluster_r_side0, label + ", side == 0", cutResidualRSide0,
                  outdir, "residual_rphi_vs_cluster_r_side0", "COLZ", true);
  DrawAndSaveHist(fout, h_residual_rphi_vs_cluster_r_side1, label + ", side == 1", cutResidualRSide1,
                  outdir, "residual_rphi_vs_cluster_r_side1", "COLZ", true);
  DrawAndSaveHist(fout, h_residual_rphi_vs_cluster_r_positive_side0, label + ", charge > 0, side == 0", cutResidualRPositiveSide0,
                  outdir, "residual_rphi_vs_cluster_r_positive_side0", "COLZ", true);
  DrawAndSaveHist(fout, h_residual_rphi_vs_cluster_r_positive_side1, label + ", charge > 0, side == 1", cutResidualRPositiveSide1,
                  outdir, "residual_rphi_vs_cluster_r_positive_side1", "COLZ", true);
  DrawAndSaveHist(fout, h_residual_rphi_vs_cluster_r_negative_side0, label + ", charge < 0, side == 0", cutResidualRNegativeSide0,
                  outdir, "residual_rphi_vs_cluster_r_negative_side0", "COLZ", true);
  DrawAndSaveHist(fout, h_residual_rphi_vs_cluster_r_negative_side1, label + ", charge < 0, side == 1", cutResidualRNegativeSide1,
                  outdir, "residual_rphi_vs_cluster_r_negative_side1", "COLZ", true);

  DrawAndSaveHist(fout, h_residual_rphi_vs_cluster_r_positive_side0_eta_lt_0p5, label + ", charge > 0, side == 0, |eta| < 0.5", cutResidualRPositiveSide0EtaLt,
                  outdir, "residual_rphi_vs_cluster_r_positive_side0_eta_lt_0p5", "COLZ", true);
  DrawAndSaveHist(fout, h_residual_rphi_vs_cluster_r_positive_side1_eta_lt_0p5, label + ", charge > 0, side == 1, |eta| < 0.5", cutResidualRPositiveSide1EtaLt,
                  outdir, "residual_rphi_vs_cluster_r_positive_side1_eta_lt_0p5", "COLZ", true);
  DrawAndSaveHist(fout, h_residual_rphi_vs_cluster_r_negative_side0_eta_lt_0p5, label + ", charge < 0, side == 0, |eta| < 0.5", cutResidualRNegativeSide0EtaLt,
                  outdir, "residual_rphi_vs_cluster_r_negative_side0_eta_lt_0p5", "COLZ", true);
  DrawAndSaveHist(fout, h_residual_rphi_vs_cluster_r_negative_side1_eta_lt_0p5, label + ", charge < 0, side == 1, |eta| < 0.5", cutResidualRNegativeSide1EtaLt,
                  outdir, "residual_rphi_vs_cluster_r_negative_side1_eta_lt_0p5", "COLZ", true);
  DrawAndSaveHist(fout, h_residual_rphi_vs_cluster_r_positive_side0_eta_gt_0p5, label + ", charge > 0, side == 0, |eta| > 0.5", cutResidualRPositiveSide0EtaGt,
                  outdir, "residual_rphi_vs_cluster_r_positive_side0_eta_gt_0p5", "COLZ", true);
  DrawAndSaveHist(fout, h_residual_rphi_vs_cluster_r_positive_side1_eta_gt_0p5, label + ", charge > 0, side == 1, |eta| > 0.5", cutResidualRPositiveSide1EtaGt,
                  outdir, "residual_rphi_vs_cluster_r_positive_side1_eta_gt_0p5", "COLZ", true);
  DrawAndSaveHist(fout, h_residual_rphi_vs_cluster_r_negative_side0_eta_gt_0p5, label + ", charge < 0, side == 0, |eta| > 0.5", cutResidualRNegativeSide0EtaGt,
                  outdir, "residual_rphi_vs_cluster_r_negative_side0_eta_gt_0p5", "COLZ", true);
  DrawAndSaveHist(fout, h_residual_rphi_vs_cluster_r_negative_side1_eta_gt_0p5, label + ", charge < 0, side == 1, |eta| > 0.5", cutResidualRNegativeSide1EtaGt,
                  outdir, "residual_rphi_vs_cluster_r_negative_side1_eta_gt_0p5", "COLZ", true);

  DrawAndSaveHist(fout, h_residual_z_vs_cluster_z, label, cutResidualR,
                  outdir, "residual_z_vs_cluster_z", "COLZ", true);
  DrawAndSaveHist(fout, h_residual_z_vs_cluster_z_positive, label + ", charge > 0", cutResidualRPositive,
                  outdir, "residual_z_vs_cluster_z_positive", "COLZ", true);
  DrawAndSaveHist(fout, h_residual_z_vs_cluster_z_negative, label + ", charge < 0", cutResidualRNegative,
                  outdir, "residual_z_vs_cluster_z_negative", "COLZ", true);
  DrawAndSaveHist(fout, h_residual_rphi_vs_pt, label, cutResidualPt,
                  outdir, "residual_rphi_vs_pt", "COLZ", true);
  DrawAndSaveHist(fout, h_dedx_vs_charge_momentum, label, dedxCut,
                  outdir, "dedx_vs_charge_momentum", "COLZ", true);
  DrawAndSaveHist(fout, h_vertex_xy, label, cutVertexAll,
                  outdir, "vertex_xy", "COLZ", true);
  DrawAndSaveHist(fout, h_vertex_xy_pt, label, cutVertexPt,
                  outdir, "vertex_xy_pt_gt_2", "COLZ", true);
  DrawAndSaveHist(fout, h_quality, label, cutQuality,
                  outdir, "quality", "HIST", false);
  DrawAndSaveHist(fout, h_rDCA_vs_chargeOverPt, label, cutDca,
                  outdir, "rDCA_vs_chargeOverPt", "COLZ", true);
  DrawAndSaveHist(fout, h_rDCA_vs_chargeOverPt_side0, label + ", side == 0", cutDcaSide0,
                  outdir, "rDCA_vs_chargeOverPt_side0", "COLZ", true);
  DrawAndSaveHist(fout, h_rDCA_vs_chargeOverPt_side1, label + ", side == 1", cutDcaSide1,
                  outdir, "rDCA_vs_chargeOverPt_side1", "COLZ", true);
  DrawAndSaveHist(fout, h_rDCA_zero_vs_chargeOverPt, label, cutDca,
                  outdir, "rDCA_zero_vs_chargeOverPt", "COLZ", true);
  DrawAndSaveHist(fout, h_rDCA_zero_vs_chargeOverPt_side0, label + ", side == 0", cutDcaSide0,
                  outdir, "rDCA_zero_vs_chargeOverPt_side0", "COLZ", true);
  DrawAndSaveHist(fout, h_rDCA_zero_vs_chargeOverPt_side1, label + ", side == 1", cutDcaSide1,
                  outdir, "rDCA_zero_vs_chargeOverPt_side1", "COLZ", true);

  fout->Close();
  delete fout;

  std::cout << "Wrote ROOT file: " << outfile << std::endl;
  std::cout << "Wrote PNG files in: " << outdir << std::endl;
}
