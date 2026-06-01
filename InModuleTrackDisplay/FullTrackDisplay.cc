#include "FullTrackDisplay.h"

#include "FullTrack.h"
#include "FullTrackContainer.h"
#include "TpcPadMap.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>

#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <TCanvas.h>
#include <TColor.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TH3D.h>
#include <TMath.h>
#include <TPolyLine3D.h>
#include <TString.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace
{
  const double TWO_PI = 2.0 * TMath::Pi();

  double wrap_to_pi(double phi)
  {
    while (phi > TMath::Pi())  phi -= TWO_PI;
    while (phi <= -TMath::Pi()) phi += TWO_PI;
    return phi;
  }


  bool is_good_number(const double x)
  {
    return (x == x && std::fabs(x) < 1.0e30);
  }

  struct ClippedFitSegment
  {
    double tbin0;
    double phi0;
    double radius0;
    double tbin1;
    double phi1;
    double radius1;
  };

  bool clip_phi_line_to_window(double& t0, double& p0, double& r0,
                               double& t1, double& p1, double& r1)
  {
    const double pmin = -TMath::Pi();
    const double pmax =  TMath::Pi();

    if ((p0 < pmin && p1 < pmin) || (p0 > pmax && p1 > pmax)) return false;

    if (std::fabs(p1 - p0) < 1.0e-12)
    {
      return (p0 >= pmin && p0 <= pmax);
    }

    if (p0 < pmin)
    {
      const double f = (pmin - p0) / (p1 - p0);
      t0 += f * (t1 - t0);
      r0 += f * (r1 - r0);
      p0 = pmin;
    }
    else if (p0 > pmax)
    {
      const double f = (pmax - p0) / (p1 - p0);
      t0 += f * (t1 - t0);
      r0 += f * (r1 - r0);
      p0 = pmax;
    }

    if (p1 < pmin)
    {
      const double f = (pmin - p0) / (p1 - p0);
      t1 = t0 + f * (t1 - t0);
      r1 = r0 + f * (r1 - r0);
      p1 = pmin;
    }
    else if (p1 > pmax)
    {
      const double f = (pmax - p0) / (p1 - p0);
      t1 = t0 + f * (t1 - t0);
      r1 = r0 + f * (r1 - r0);
      p1 = pmax;
    }

    return (p0 >= pmin && p0 <= pmax && p1 >= pmin && p1 <= pmax);
  }

  void make_clipped_fit_segments(const std::vector<double>& tbin,
                                 const std::vector<double>& phi_unwrapped,
                                 const std::vector<double>& radius,
                                 std::vector<ClippedFitSegment>& segments)
  {
    segments.clear();
    if (tbin.size() < 2) return;

    for (unsigned int i = 1; i < tbin.size(); ++i)
    {
      const double p0raw = phi_unwrapped[i - 1];
      const double p1raw = phi_unwrapped[i];

      const int kmin = static_cast<int>(std::floor((-TMath::Pi() - std::max(p0raw, p1raw)) / TWO_PI)) - 1;
      const int kmax = static_cast<int>(std::ceil( ( TMath::Pi() - std::min(p0raw, p1raw)) / TWO_PI)) + 1;

      for (int k = kmin; k <= kmax; ++k)
      {
        double t0 = tbin[i - 1];
        double t1 = tbin[i];
        double r0 = radius[i - 1];
        double r1 = radius[i];
        double p0 = p0raw + static_cast<double>(k) * TWO_PI;
        double p1 = p1raw + static_cast<double>(k) * TWO_PI;

        if (!clip_phi_line_to_window(t0, p0, r0, t1, p1, r1)) continue;
        if (!is_good_number(t0) || !is_good_number(t1) ||
            !is_good_number(p0) || !is_good_number(p1) ||
            !is_good_number(r0) || !is_good_number(r1)) continue;

        ClippedFitSegment s;
        s.tbin0 = t0;
        s.phi0 = p0;
        s.radius0 = r0;
        s.tbin1 = t1;
        s.phi1 = p1;
        s.radius1 = r1;
        segments.push_back(s);
      }
    }
  }

  int track_color(const unsigned int itrk)
  {
    static const int colors[] = {
      kRed + 1,
      kBlue + 1,
      kGreen + 2,
      kMagenta + 1,
      kCyan + 2,
      kOrange + 7,
      kViolet + 1,
      kAzure + 1,
      kPink + 7,
      kTeal + 3
    };

    const unsigned int ncolors = sizeof(colors) / sizeof(colors[0]);
    return colors[itrk % ncolors];
  }

  int layer_to_region(const unsigned int layer)
  {
    if (layer < 7 || layer > 54) return -1;
    const int region = static_cast<int>((layer - 7) / 16);
    if (region < 0 || region >= 3) return -1;
    return region;
  }

  unsigned long long make_unique_hit_id(const TrkrDefs::hitsetkey hsk,
                                        const TrkrDefs::hitkey hk)
  {
    return (static_cast<unsigned long long>(hsk) << 32) |
           static_cast<unsigned long long>(hk);
  }

  void style_hit_graph(TGraph* g, const int color)
  {
    if (!g) return;
    g->SetMarkerStyle(20);
    g->SetMarkerSize(0.65);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetLineWidth(1);
  }

  void style_fit_graph(TGraph* g, const int color)
  {
    if (!g) return;
    g->SetMarkerStyle(1);
    g->SetMarkerSize(0.0);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetLineWidth(3);
    g->SetLineStyle(2);
  }

  void style_fit_line_3d(TPolyLine3D* line, const int color)
  {
    if (!line) return;
    line->SetLineColor(color);
    line->SetLineWidth(4);
    line->SetLineStyle(1);
  }

  TGraph* make_graph(const std::vector<double>& x,
                     const std::vector<double>& y,
                     const std::string& name,
                     const int color,
                     const bool fit_style)
  {
    TGraph* g = new TGraph(static_cast<int>(x.size()));
    g->SetName(name.c_str());

    for (unsigned int i = 0; i < x.size(); ++i)
    {
      g->SetPoint(static_cast<int>(i), x[i], y[i]);
    }

    if (fit_style) style_fit_graph(g, color);
    else style_hit_graph(g, color);

    return g;
  }

  struct SideBundle
  {
    SideBundle()
      : h3_adc_hits(0)
      , c3_adc_hits_fits(0)
      , c_tbin_phi(0)
      , c_tbin_layer(0)
      , mg_tbin_phi_hits(0)
      , mg_tbin_phi_fits(0)
      , mg_tbin_layer_hits(0)
      , mg_tbin_layer_fits(0)
    {}

    TH3D* h3_adc_hits;
    TCanvas* c3_adc_hits_fits;
    TCanvas* c_tbin_phi;
    TCanvas* c_tbin_layer;

    TMultiGraph* mg_tbin_phi_hits;
    TMultiGraph* mg_tbin_phi_fits;
    TMultiGraph* mg_tbin_layer_hits;
    TMultiGraph* mg_tbin_layer_fits;

    std::vector<TPolyLine3D*> fit_lines_3d;
    std::vector<std::string> fit_line_names_3d;
    std::set<unsigned long long> filled_hit_ids;
  };

  void make_side_bundle(SideBundle& b,
                        const int side,
                        const unsigned int evt)
  {
    const TString tag = Form("side%d", side);

    b.h3_adc_hits =
      new TH3D(Form("h3_%s_fulltrack_adc_hits", tag.Data()),
               Form("event %u side %d: FULLTRACK ADC-weighted hits;timebin;global #phi;radius",
                    evt, side),
               512, -0.5, 511.5,
               720, -M_PI, M_PI,
               108, 20.0, 80.0);
    b.h3_adc_hits->SetStats(0);

    b.mg_tbin_phi_hits = new TMultiGraph();
    b.mg_tbin_phi_hits->SetName(Form("mg_%s_fulltrack_tbin_phi_hits", tag.Data()));
    b.mg_tbin_phi_hits->SetTitle(Form("event %u side %d: FULLTRACK hits;timebin;global #phi", evt, side));

    b.mg_tbin_phi_fits = new TMultiGraph();
    b.mg_tbin_phi_fits->SetName(Form("mg_%s_fulltrack_tbin_phi_fits", tag.Data()));
    b.mg_tbin_phi_fits->SetTitle(Form("event %u side %d: FULLTRACK fits;timebin;global #phi", evt, side));

    b.mg_tbin_layer_hits = new TMultiGraph();
    b.mg_tbin_layer_hits->SetName(Form("mg_%s_fulltrack_tbin_layer_hits", tag.Data()));
    b.mg_tbin_layer_hits->SetTitle(Form("event %u side %d: FULLTRACK hits;timebin;radius", evt, side));

    b.mg_tbin_layer_fits = new TMultiGraph();
    b.mg_tbin_layer_fits->SetName(Form("mg_%s_fulltrack_tbin_layer_fits", tag.Data()));
    b.mg_tbin_layer_fits->SetTitle(Form("event %u side %d: FULLTRACK fits;timebin;radius", evt, side));
  
  }

  SideBundle& get_side_bundle(std::map<int, SideBundle>& bundles,
                              const int side,
                              const unsigned int evt)
  {
    std::map<int, SideBundle>::iterator it = bundles.find(side);
    if (it == bundles.end())
    {
      const std::pair<std::map<int, SideBundle>::iterator, bool> inserted =
        bundles.insert(std::make_pair(side, SideBundle()));
      it = inserted.first;
      make_side_bundle(it->second, side, evt);
    }
    return it->second;
  }

  void write_side_bundle(SideBundle& b, const int side)
  {
    const TString tag = Form("side%d", side);

    if (b.h3_adc_hits)
    {
      b.c3_adc_hits_fits =
        new TCanvas(Form("c3_%s_fulltrack_adc_hits_fits", tag.Data()),
                    Form("side %d FULLTRACK ADC hits and full-track fits", side),
                    1200, 900);

      b.h3_adc_hits->Draw("BOX2Z");
      for (unsigned int iline = 0; iline < b.fit_lines_3d.size(); ++iline)
      {
        if (!b.fit_lines_3d[iline]) continue;
        b.fit_lines_3d[iline]->Draw("same");
      }

      b.c3_adc_hits_fits->Modified();
      b.c3_adc_hits_fits->Update();
      b.c3_adc_hits_fits->Write();
    }
    /*
    if (b.mg_tbin_phi_hits && b.mg_tbin_phi_fits)
    {
      b.c_tbin_phi =
        new TCanvas(Form("c_%s_fulltrack_tbin_phi", tag.Data()),
                    Form("side %d FULLTRACK timebin vs global phi", side),
                    1200, 900);
      b.mg_tbin_phi_hits->Draw("AP");
      b.mg_tbin_phi_fits->Draw("L");
      b.c_tbin_phi->Modified();
      b.c_tbin_phi->Update();
      b.c_tbin_phi->Write();
    }*/

    if (b.mg_tbin_layer_hits && b.mg_tbin_layer_fits)
    {
      b.c_tbin_layer =
        new TCanvas(Form("c_%s_fulltrack_tbin_layer", tag.Data()),
                    Form("side %d FULLTRACK timebin vs radius", side),
                    1200, 900);
      b.mg_tbin_layer_hits->Draw("AP");
      b.mg_tbin_layer_fits->Draw("L");
      b.c_tbin_layer->Modified();
      b.c_tbin_layer->Update();
      b.c_tbin_layer->Write();
    }

    /*
    if (b.mg_tbin_phi_hits) b.mg_tbin_phi_hits->Write();
    if (b.mg_tbin_phi_fits) b.mg_tbin_phi_fits->Write();
    if (b.mg_tbin_layer_hits) b.mg_tbin_layer_hits->Write();
    if (b.mg_tbin_layer_fits) b.mg_tbin_layer_fits->Write();
    */
  }
}

FullTrackDisplay::HitPoint::HitPoint()
  : ok(false)
  , hitsetkey(0)
  , hitkey(0)
  , layer(0)
  , region(0)
  , sector(0)
  , side(0)
  , pad(0)
  , tbin(0)
  , adc(0)
  , global_phi(0.0)
  , radius(0.0)
{
}

FullTrackDisplay::FullTrackDisplay(const std::string& name,
                                   const std::string& outfilename,
                                   const std::string& trackNodeName,
                                   const unsigned int maxEventDisplays)
  : SubsysReco(name)
  , m_outfilename(outfilename)
  , m_trackNodeName(trackNodeName)
  , m_maxEventDisplays(maxEventDisplays)
  , m_evt(0)
  , m_eventsSaved(0)
  , m_outfile(0)
  , m_tracks(0)
  , m_hits(0)
  , m_padMap(0)
{
}

FullTrackDisplay::~FullTrackDisplay()
{
  // TPC_PADMAP is owned by the node tree, not by this display module.
  m_padMap = 0;
}

int FullTrackDisplay::Init(PHCompositeNode*)
{
  m_outfile = new TFile(m_outfilename.c_str(), "RECREATE");
  if (!m_outfile || m_outfile->IsZombie())
  {
    std::cerr << "FullTrackDisplay::Init - cannot open output file "
              << m_outfilename << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  std::cout << "FullTrackDisplay::Init - writing up to "
            << m_maxEventDisplays << " events to "
            << m_outfilename << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

int FullTrackDisplay::process_event(PHCompositeNode* topNode)
{
  if (!get_nodes(topNode))
  {
    ++m_evt;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  if (!m_outfile || m_eventsSaved >= m_maxEventDisplays)
  {
    ++m_evt;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  m_outfile->cd();
  TDirectory* event_dir =
    m_outfile->mkdir(Form("event_%u_fulltracks", m_evt));

  if (!event_dir)
  {
    std::cerr << "FullTrackDisplay::process_event - failed to create event directory"
              << std::endl;
    ++m_evt;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  event_dir->cd();

  std::map<int, SideBundle> bundles;

  for (unsigned int itrk = 0; itrk < m_tracks->size(); ++itrk)
  {
    const FullTrack* trk = m_tracks->get_track(itrk);
    if (!trk) continue;

    const int side = trk->get_side();
    SideBundle& b = get_side_bundle(bundles, side, m_evt);
    const int color = track_color(itrk);
    const unsigned int tid = trk->get_track_id();

    std::vector<HitPoint> pts;
    pts.reserve(trk->size_hit_indices());

    std::vector<double> hit_tbin;
    std::vector<double> hit_phi;
    std::vector<double> hit_layer;

    for (unsigned int ihit = 0; ihit < trk->size_hit_indices(); ++ihit)
    {
      const FullTrack::HitIndex idx = trk->get_hit_index(ihit);
      const HitPoint p = make_hit_point(idx.first, idx.second);
      if (!p.ok) continue;
      if (p.side != side) continue;

      const unsigned long long uid = make_unique_hit_id(idx.first, idx.second);
      if (b.filled_hit_ids.find(uid) == b.filled_hit_ids.end())
      {
        if (b.h3_adc_hits)
        {
          b.h3_adc_hits->Fill(static_cast<double>(p.tbin),
                              p.global_phi,
                              p.radius,
                              static_cast<double>(p.adc));
        }
        b.filled_hit_ids.insert(uid);
      }

      pts.push_back(p);
      hit_tbin.push_back(static_cast<double>(p.tbin));
      hit_phi.push_back(p.global_phi);
      hit_layer.push_back(p.radius);
    }

    if (pts.empty()) continue;

    const std::string hit_phi_name =
      Form("g_side%d_fulltrk%u_tbin_phi_hits", side, tid);
    const std::string hit_layer_name =
      Form("g_side%d_fulltrk%u_tbin_layer_hits", side, tid);

    b.mg_tbin_phi_hits->Add(make_graph(hit_tbin, hit_phi,
                                       hit_phi_name, color, false), "P");
    b.mg_tbin_layer_hits->Add(make_graph(hit_tbin, hit_layer,
                                         hit_layer_name, color, false), "P");

    const unsigned int first_layer = trk->get_first_layer();
    const unsigned int last_layer = trk->get_last_layer();

    if (last_layer >= first_layer &&
        is_good_number(trk->get_phi_slope()) &&
        is_good_number(trk->get_phi_intercept()) &&
        is_good_number(trk->get_tbin_slope()) &&
        is_good_number(trk->get_tbin_intercept()))
    {
      std::vector<double> fit_tbin;
      std::vector<double> fit_phi_unwrapped;
      std::vector<double> fit_layer;

      for (unsigned int il = first_layer; il <= last_layer; ++il)
      {
        const int reg = layer_to_region(il);
        if (reg < 0) continue;

        const double radius = m_padMap->get_local_radius(static_cast<unsigned int>(reg), il);
        const double tbin = trk->get_tbin_slope() * radius +
                            trk->get_tbin_intercept();
        const double phi = trk->get_phi_slope() * radius +
                           trk->get_phi_intercept();

        if (!is_good_number(radius) || !is_good_number(tbin) || !is_good_number(phi)) continue;

        fit_tbin.push_back(tbin);
        fit_phi_unwrapped.push_back(phi);
        fit_layer.push_back(radius);
      }

      std::vector<ClippedFitSegment> fit_segments;
      make_clipped_fit_segments(fit_tbin, fit_phi_unwrapped, fit_layer, fit_segments);

      for (unsigned int iseg = 0; iseg < fit_segments.size(); ++iseg)
      {
        const ClippedFitSegment& seg = fit_segments[iseg];

        TPolyLine3D* line = new TPolyLine3D(2);
        const std::string line_name =
          Form("line3_side%d_fulltrk%u_fit_seg%u", side, tid, iseg);

        line->SetPoint(0, seg.tbin0, seg.phi0, seg.radius0);
        line->SetPoint(1, seg.tbin1, seg.phi1, seg.radius1);

        style_fit_line_3d(line, color);
        b.fit_lines_3d.push_back(line);
        b.fit_line_names_3d.push_back(line_name);

        std::vector<double> seg_tbin;
        std::vector<double> seg_phi;
        std::vector<double> seg_layer;
        seg_tbin.push_back(seg.tbin0);
        seg_tbin.push_back(seg.tbin1);
        seg_phi.push_back(seg.phi0);
        seg_phi.push_back(seg.phi1);
        seg_layer.push_back(seg.radius0);
        seg_layer.push_back(seg.radius1);

        const std::string fit_phi_name =
          Form("g_side%d_fulltrk%u_tbin_phi_fit_seg%u", side, tid, iseg);
        const std::string fit_layer_name =
          Form("g_side%d_fulltrk%u_tbin_layer_fit_seg%u", side, tid, iseg);

        b.mg_tbin_phi_fits->Add(make_graph(seg_tbin, seg_phi,
                                           fit_phi_name, color, true), "L");
        b.mg_tbin_layer_fits->Add(make_graph(seg_tbin, seg_layer,
                                             fit_layer_name, color, true), "L");
      }
    }
  }

  for (std::map<int, SideBundle>::iterator it = bundles.begin();
       it != bundles.end(); ++it)
  {
    write_side_bundle(it->second, it->first);
  }

  ++m_eventsSaved;
  ++m_evt;
  return Fun4AllReturnCodes::EVENT_OK;
}

int FullTrackDisplay::End(PHCompositeNode*)
{
  if (!m_outfile) return Fun4AllReturnCodes::EVENT_OK;

  m_outfile->cd();
  m_outfile->Write();
  m_outfile->Close();

  delete m_outfile;
  m_outfile = 0;

  std::cout << "FullTrackDisplay::End - events seen: " << m_evt
            << ", events written: " << m_eventsSaved
            << ", output file: " << m_outfilename << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

bool FullTrackDisplay::get_nodes(PHCompositeNode* topNode)
{
  m_hits = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  m_tracks = findNode::getClass<FullTrackContainer>(topNode,
                                                   m_trackNodeName.c_str());

  if (!m_tracks)
  {
    const char* candidate_names[] = {
      "FULLTRACKS",
      "FullTracks",
      "FullTrackContainer",
      "FULLTRACKCONTAINER",
      "FULLTRACKS_CONTAINER"
    };

    for (unsigned int i = 0;
         i < sizeof(candidate_names) / sizeof(candidate_names[0]) && !m_tracks;
         ++i)
    {
      m_tracks = findNode::getClass<FullTrackContainer>(topNode,
                                                       candidate_names[i]);
      if (m_tracks) m_trackNodeName = candidate_names[i];
    }
  }

  if (!m_tracks)
  {
    std::cerr << "FullTrackDisplay - could not find FullTrackContainer node"
              << std::endl;
    return false;
  }

  if (!m_hits)
  {
    std::cerr << "FullTrackDisplay - missing TRKR_HITSET" << std::endl;
    return false;
  }

  m_padMap = findNode::getClass<TpcPadMap>(topNode, "TPC_PADMAP");
  if (!m_padMap)
  {
    std::cerr << "FullTrackDisplay - missing TPC_PADMAP" << std::endl;
    return false;
  }

  if (!m_padMap->isValid())
  {
    std::cerr << "FullTrackDisplay - TPC_PADMAP is invalid" << std::endl;
    return false;
  }

  return true;
}

FullTrackDisplay::HitPoint
FullTrackDisplay::make_hit_point(const TrkrDefs::hitsetkey hsk,
                                 const TrkrDefs::hitkey hk) const
{
  HitPoint p;

  TrkrHitSet* hitset = m_hits->findHitSet(hsk);
  if (!hitset) return p;

  TrkrHit* hit = hitset->getHit(hk);
  if (!hit) return p;

  p.hitsetkey = hsk;
  p.hitkey = hk;

  p.layer = TrkrDefs::getLayer(hsk);
  const int region = layer_to_region(p.layer);
  if (region < 0) return p;

  p.region = static_cast<unsigned int>(region);
  p.side = static_cast<int>(TpcDefs::getSide(hsk));
  p.sector = static_cast<unsigned int>(TpcDefs::getSectorId(hsk));
  p.pad = TpcDefs::getPad(hk);
  p.tbin = TpcDefs::getTBin(hk);
  p.adc = hit->getAdc();
  if (!m_padMap) return p;
  p.global_phi =
    wrap_to_pi(m_padMap->get_global_phi(static_cast<unsigned int>(p.side),
                                        p.region,
                                        static_cast<double>(p.pad)));
  p.radius = m_padMap->get_local_radius(p.region, p.layer);

  if (!is_good_number(p.global_phi) || !is_good_number(p.radius)) return p;

  p.ok = true;
  return p;
}
