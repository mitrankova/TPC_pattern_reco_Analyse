#include "InModuleTrackDisplay.h"

#include "TpcPadMap.h"

#include "InModuleTrack.h"
#include "InModuleTrackContainer.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>

#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <g4detectors/PHG4TpcGeom.h>
#include <g4detectors/PHG4TpcGeomContainer.h>

#include <TDirectory.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TH3D.h>
#include <TPolyLine3D.h>
#include <TMath.h>
#include <TString.h>
#include <TColor.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <vector>

namespace
{
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

  void style_hit_graph(TGraph* g, const int color)
  {
    if (!g) return;

    g->SetMarkerStyle(20);
    g->SetMarkerSize(0.8);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetLineWidth(2);
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

  int layer_to_module(const unsigned int layer)
  {
    // TPC layer groups requested by the user:
    // module 0: layers  7--22
    // module 1: layers 23--38
    // module 2: layers 39--54
    if (layer < 7 || layer > 54) return -1;
    const int module = static_cast<int>((layer - 7) / 16);
    if (module < 0 || module >= 3) return -1;
    return module;
  }

  unsigned long long make_unique_hit_id(const TrkrDefs::hitsetkey hsk,
                                        const TrkrDefs::hitkey hk)
  {
    return (static_cast<unsigned long long>(hsk) << 32) |
           static_cast<unsigned long long>(hk);
  }


  // Convert the in-module fit coordinate (sector-local phi) to global TPC phi.
  // InModuleTracks still fits in local module coordinates; only the display is global.
  double local_phi_to_global_phi(const TpcPadMap* padMap,
                                 const int side,
                                 const int sector,
                                 const unsigned int region,
                                 const double local_phi)
  {
    if (!padMap) return 0.0;

    const double phi_width = padMap->get_phi_bin_width(region);
    const unsigned int pads_per_sector = padMap->get_pads_per_sector(region);

    if (phi_width <= 0.0 || pads_per_sector == 0U) return 0.0;

    const double local_pad = local_phi / phi_width;
    const double full_pad =
      static_cast<double>(sector) * static_cast<double>(pads_per_sector) +
      local_pad;

    return padMap->get_global_phi(static_cast<unsigned int>(side),
                                  region,
                                  full_pad);
  }

  double local_phi_to_full_pad(const TpcPadMap* padMap,
                               const int sector,
                               const unsigned int region,
                               const double local_phi)
  {
    if (!padMap) return 0.0;

    const double phi_width = padMap->get_phi_bin_width(region);
    const unsigned int pads_per_sector = padMap->get_pads_per_sector(region);

    if (phi_width <= 0.0 || pads_per_sector == 0U) return 0.0;

    return static_cast<double>(sector) * static_cast<double>(pads_per_sector) +
           local_phi / phi_width;
  }


  bool is_good_number(const double x)
  {
    return (x == x && std::fabs(x) < 1.0e30);
  }

  double wrap_delta_phi_display(double dphi)
  {
    while (dphi > TMath::Pi()) dphi -= 2.0 * TMath::Pi();
    while (dphi < -TMath::Pi()) dphi += 2.0 * TMath::Pi();
    return dphi;
  }

  double sagitta_model_display(const double xrot,
                               const double S,
                               const double x0,
                               const double invR)
  {
    const double dx = xrot - x0;
    const double invR2 = invR * invR;
    const double dx2 = dx * dx;

    // Same truncated expansion as the Python notebook:
    // _sagitta = S - 1/2 invR dx^2 - 1/8 invR^3 dx^4 - 1/16 invR^5 dx^6
    return S
      - 0.5 * invR * dx2
      - 0.125 * invR * invR2 * dx2 * dx2
      - 0.0625 * invR * invR2 * invR2 * dx2 * dx2 * dx2;
  }

  bool sagitta_phi_at_radius_display(const double rloc,
                                     const double mean_phi,
                                     const double S,
                                     const double x0,
                                     const double invR,
                                     const double theta,
                                     const double bline,
                                     double& phi_out)
  {
    if (rloc <= 0.0) return false;
    if (!is_good_number(S) || !is_good_number(x0) || !is_good_number(invR) ||
        !is_good_number(theta) || !is_good_number(bline) ||
        !is_good_number(mean_phi)) return false;

    const double c = std::cos(theta);
    const double s = std::sin(theta);

    // Sagitta parameters are fitted in the 2D plane:
    //   x = radius, y = local phi.
    // weighted_sagitta_fit first rotates that plane:
    //   xrot =  cos(theta)*x + sin(theta)*(y - bline)
    //   yrot = -sin(theta)*x + cos(theta)*(y - bline)
    // and fits yrot = sagitta_model(xrot).
    // For display we know x = rloc and solve only for y = local phi.
    const double phimin = mean_phi - 0.8;
    const double phimax = mean_phi + 0.8;
    const int nscan = 240;

    bool have_best = false;
    double best_phi = 0.0;
    double best_dist = 1.0e99;

    double pprev = phimin;
    double xprev_rot = c * rloc + s * (pprev - bline);
    double yprev_rot = -s * rloc + c * (pprev - bline);
    double gprev = yprev_rot - sagitta_model_display(xprev_rot, S, x0, invR);

    double best_abs_g = std::fabs(gprev);
    double best_abs_phi = pprev;

    for (int iscan = 1; iscan <= nscan; ++iscan)
    {
      const double pcur = phimin + (phimax - phimin) * static_cast<double>(iscan) / static_cast<double>(nscan);
      const double xcur_rot = c * rloc + s * (pcur - bline);
      const double ycur_rot = -s * rloc + c * (pcur - bline);
      const double gcur = ycur_rot - sagitta_model_display(xcur_rot, S, x0, invR);

      if (std::fabs(gcur) < best_abs_g)
      {
        best_abs_g = std::fabs(gcur);
        best_abs_phi = pcur;
      }

      bool bracket = false;
      if (gprev == 0.0 || gcur == 0.0) bracket = true;
      if ((gprev < 0.0 && gcur > 0.0) || (gprev > 0.0 && gcur < 0.0)) bracket = true;

      if (bracket)
      {
        double pa = pprev;
        double pb = pcur;
        double ga = gprev;

        for (int ib = 0; ib < 50; ++ib)
        {
          const double pm = 0.5 * (pa + pb);
          const double xm_rot = c * rloc + s * (pm - bline);
          const double ym_rot = -s * rloc + c * (pm - bline);
          const double gm = ym_rot - sagitta_model_display(xm_rot, S, x0, invR);

          if ((ga < 0.0 && gm <= 0.0) || (ga > 0.0 && gm >= 0.0))
          {
            pa = pm;
            ga = gm;
          }
          else
          {
            pb = pm;
          }
        }

        const double phi = 0.5 * (pa + pb);
        const double dphi = std::fabs(phi - mean_phi);
        if (!have_best || dphi < best_dist)
        {
          have_best = true;
          best_dist = dphi;
          best_phi = phi;
        }
      }

      pprev = pcur;
      gprev = gcur;
    }

    if (have_best)
    {
      phi_out = best_phi;
      return true;
    }

    // Last-resort near-tangent case: accept the closest point only if the
    // fitted curve is numerically very close in local-phi units.
    if (best_abs_g < 2.0e-2)
    {
      phi_out = best_abs_phi;
      return true;
    }

    return false;
  }

  struct GroupKey
  {
    GroupKey()
      : side(-1)
      , sector(-1)
      , module(-1)
    {}

    GroupKey(const int side_in, const int sector_in, const int module_in)
      : side(side_in)
      , sector(sector_in)
      , module(module_in)
    {}

    int side;
    int sector;
    int module;

    bool operator<(const GroupKey& rhs) const
    {
      if (side != rhs.side) return side < rhs.side;
      if (sector != rhs.sector) return sector < rhs.sector;
      return module < rhs.module;
    }
  };

  struct GraphBundle
  {
    GraphBundle()
      : mg_tbin_pad_hits(0)
      , mg_tbin_layer_hits(0)
      , mg_pad_layer_hits(0)
      , mg_tbin_pad_fits(0)
      , mg_tbin_layer_fits(0)
      , mg_pad_layer_fits(0)
      , h3_adc_hits(0)
      , h3_adc_unassociated_hits(0)
      , h3_adc_globalphi_hits(0)
      , c3_adc_hits_fits(0)
      , c3_adc_unassociated_hits(0)
      , c3_adc_globalphi_hits_fits(0)
    {}

    TMultiGraph* mg_tbin_pad_hits;
    TMultiGraph* mg_tbin_layer_hits;
    TMultiGraph* mg_pad_layer_hits;

    TMultiGraph* mg_tbin_pad_fits;
    TMultiGraph* mg_tbin_layer_fits;
    TMultiGraph* mg_pad_layer_fits;

    // 3D hit display:
    //   x = timebin, y = pad, z = layer
    //   bin content = ADC-weighted hit entries
    TH3D* h3_adc_hits;                  // hits associated with tracks
    TH3D* h3_adc_unassociated_hits;     // TRKR_HITSET hits not associated with any track

    // Global-phi 3D hit display:
    //   x = timebin, y = global phi, z = radius
    //   bin content = ADC-weighted hit entries
    TH3D* h3_adc_globalphi_hits;

    TCanvas* c3_adc_hits_fits;
    TCanvas* c3_adc_unassociated_hits;
    TCanvas* c3_adc_globalphi_hits_fits;

    std::vector<TPolyLine3D*> fit_lines_3d;
    std::vector<std::string> fit_line_names_3d;

    std::vector<TPolyLine3D*> globalphi_fit_lines_3d;
    std::vector<std::string> globalphi_fit_line_names_3d;

    // avoid filling the same hit more than once if tracks share hits
    std::set<unsigned long long> filled_hit_ids;
  };

  void make_bundle_graphs(GraphBundle& b,
                          const GroupKey& key,
                          const unsigned int evt,
                          const TpcPadMap* padMap)
  {
    const TString tag = Form("s%d_sec%02d_mod%d",
                             key.side, key.sector, key.module);

    b.mg_tbin_pad_hits = new TMultiGraph();
    b.mg_tbin_pad_hits->SetName(Form("mg_%s_tbin_pad_hits", tag.Data()));
    b.mg_tbin_pad_hits->SetTitle(Form("event %u side %d sector %d module %d: hits;timebin;pad",
                                      evt, key.side, key.sector, key.module));

    b.mg_tbin_layer_hits = new TMultiGraph();
    b.mg_tbin_layer_hits->SetName(Form("mg_%s_tbin_layer_hits", tag.Data()));
    b.mg_tbin_layer_hits->SetTitle(Form("event %u side %d sector %d module %d: hits;timebin;layer",
                                        evt, key.side, key.sector, key.module));

    b.mg_pad_layer_hits = new TMultiGraph();
    b.mg_pad_layer_hits->SetName(Form("mg_%s_pad_layer_hits", tag.Data()));
    b.mg_pad_layer_hits->SetTitle(Form("event %u side %d sector %d module %d: hits;pad;layer",
                                       evt, key.side, key.sector, key.module));

    b.mg_tbin_pad_fits = new TMultiGraph();
    b.mg_tbin_pad_fits->SetName(Form("mg_%s_tbin_pad_fits", tag.Data()));
    b.mg_tbin_pad_fits->SetTitle(Form("event %u side %d sector %d module %d: fits;timebin;pad",
                                      evt, key.side, key.sector, key.module));

    b.mg_tbin_layer_fits = new TMultiGraph();
    b.mg_tbin_layer_fits->SetName(Form("mg_%s_tbin_layer_fits", tag.Data()));
    b.mg_tbin_layer_fits->SetTitle(Form("event %u side %d sector %d module %d: fits;timebin;layer",
                                        evt, key.side, key.sector, key.module));

    b.mg_pad_layer_fits = new TMultiGraph();
    b.mg_pad_layer_fits->SetName(Form("mg_%s_pad_layer_fits", tag.Data()));
    b.mg_pad_layer_fits->SetTitle(Form("event %u side %d sector %d module %d: fits;pad;layer",
                                       evt, key.side, key.sector, key.module));

    const int layer_min = 7 + 16 * key.module;
    const int layer_max = layer_min + 15;
    const unsigned int nPads = padMap ? padMap->get_pads_per_sector(static_cast<unsigned int>(key.module)) : 0U;
    const int pad_min = static_cast<int>(nPads) * key.sector;
    const int pad_max = pad_min + static_cast<int>(nPads) - 1;

    b.h3_adc_hits = new TH3D(Form("h3_%s_adc_hits", tag.Data()),
                             Form("event %u side %d sector %d module %d: ADC-weighted hits;timebin;pad;layer",
                                  evt, key.side, key.sector, key.module),
                             512, -0.5, 511.5,
                             static_cast<int>(nPads), static_cast<double>(pad_min) - 0.5, static_cast<double>(pad_max) + 0.5,
                             16, static_cast<double>(layer_min) - 0.5,
                             static_cast<double>(layer_max) + 0.5);

    b.h3_adc_hits->SetStats(0);

    const double radius_min =
      padMap ? padMap->get_local_radius(static_cast<unsigned int>(key.module), static_cast<unsigned int>(layer_min)) - 0.5 : 0.0;
    const double radius_max =
      padMap ? padMap->get_local_radius(static_cast<unsigned int>(key.module), static_cast<unsigned int>(layer_max)) + 0.5 : 1.0;

    const double global_phi_min = -TMath::Pi();
    const double global_phi_max =  TMath::Pi();

    b.h3_adc_globalphi_hits =
      new TH3D(Form("h3_%s_adc_global_hits", tag.Data()),
               Form("event %u side %d sector %d module %d: ADC-weighted hits in global coordinates;timebin;global #phi;radius",
                    evt, key.side, key.sector, key.module),
               512, -0.5, 511.5,
               720, global_phi_min, global_phi_max,
               16, radius_min, radius_max);

    b.h3_adc_globalphi_hits->SetStats(0);


    b.h3_adc_unassociated_hits =
    new TH3D(Form("h3_%s_adc_unassociated_hits", tag.Data()),
            Form("event %u side %d sector %d module %d: ADC-weighted unassociated TRKR_HITSET hits;timebin;pad;layer",
                  evt, key.side, key.sector, key.module),
            512, -0.5, 511.5,
             static_cast<int>(nPads), static_cast<double>(pad_min) - 0.5, static_cast<double>(pad_max) + 0.5,
            16, static_cast<double>(layer_min) - 0.5,
            static_cast<double>(layer_max) + 0.5);

    b.h3_adc_unassociated_hits->SetStats(0);


  }


  GraphBundle& get_bundle(std::map<GroupKey, GraphBundle>& bundles,
                          const GroupKey& key,
                          const unsigned int evt,
                          const TpcPadMap* padMap)
  {
    std::map<GroupKey, GraphBundle>::iterator it = bundles.find(key);

    if (it == bundles.end())
    {
      const std::pair<std::map<GroupKey, GraphBundle>::iterator, bool> inserted =
        bundles.insert(std::make_pair(key, GraphBundle()));

      it = inserted.first;
      make_bundle_graphs(it->second, key, evt, padMap);
    }

    return it->second;
  }

  void write_bundle(GraphBundle& b, const GroupKey& key)
  {
    //if (b.h3_adc_hits) b.h3_adc_hits->Write();
    //if (b.h3_adc_unassociated_hits) b.h3_adc_unassociated_hits->Write();

    /*for (unsigned int iline = 0; iline < b.fit_lines_3d.size(); ++iline)
    {
      if (b.fit_lines_3d[iline]) b.fit_lines_3d[iline]->Write(b.fit_line_names_3d[iline].c_str());
    }*/

      if (b.h3_adc_hits)
      {
        const TString tag = Form("s%d_sec%02d_mod%d",
                                key.side, key.sector, key.module);

        b.c3_adc_hits_fits = new TCanvas(Form("c3_%s_adc_hits_fits", tag.Data()),
                                        Form("side %d sector %d module %d ADC hits and fits",
                                              key.side, key.sector, key.module),
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

      if (b.h3_adc_globalphi_hits)
      {
        const TString tag = Form("s%d_sec%02d_mod%d",
                                key.side, key.sector, key.module);

        b.c3_adc_globalphi_hits_fits =
          new TCanvas(Form("c3_%s_adc_globalphi_hits_fits", tag.Data()),
                      Form("side %d sector %d module %d ADC hits and global-phi fits",
                            key.side, key.sector, key.module),
                      1200, 900);

        b.h3_adc_globalphi_hits->Draw("BOX2Z");

        for (unsigned int iline = 0; iline < b.globalphi_fit_lines_3d.size(); ++iline)
        {
          if (!b.globalphi_fit_lines_3d[iline]) continue;
          b.globalphi_fit_lines_3d[iline]->Draw("same");
          //b.globalphi_fit_lines_3d[iline]->Write(b.globalphi_fit_line_names_3d[iline].c_str());
        }

        b.c3_adc_globalphi_hits_fits->Modified();
        b.c3_adc_globalphi_hits_fits->Update();
        b.c3_adc_globalphi_hits_fits->Write();
      }

    if (b.h3_adc_unassociated_hits)
    {
      const TString tag = Form("s%d_sec%02d_mod%d",
                              key.side, key.sector, key.module);

      b.c3_adc_unassociated_hits =
        new TCanvas(Form("c3_%s_adc_unassociated_hits", tag.Data()),
                    Form("side %d sector %d module %d unassociated TRKR_HITSET hits",
                        key.side, key.sector, key.module),
                    1200, 900);

      b.h3_adc_unassociated_hits->Draw("BOX2Z");

      b.c3_adc_unassociated_hits->Modified();
      b.c3_adc_unassociated_hits->Update();
      b.c3_adc_unassociated_hits->Write();
    }

/*
    if (b.mg_tbin_pad_hits) b.mg_tbin_pad_hits->Write();
    //if (b.mg_tbin_layer_hits) b.mg_tbin_layer_hits->Write();
    //if (b.mg_pad_layer_hits) b.mg_pad_layer_hits->Write();

    if (b.mg_tbin_pad_fits) b.mg_tbin_pad_fits->Write();
    if (b.mg_tbin_layer_fits) b.mg_tbin_layer_fits->Write();
    if (b.mg_pad_layer_fits) b.mg_pad_layer_fits->Write();
    */
  }
}

InModuleTrackDisplay::HitPoint::HitPoint()
  : ok(false)
  , hitsetkey(0)
  , hitkey(0)
  , layer(0)
  , pad(0)
  , tbin(0)
  , adc(0)
  , x(0)
  , y(0)
  , r(0)
  , phi(0)
  , local_phi(0)
{}

InModuleTrackDisplay::InModuleTrackDisplay(const std::string& name,
                                           const std::string& outfilename,
                                           const std::string& trackNodeName,
                                           unsigned int maxEventDisplays)
  : SubsysReco(name)
  , m_outfilename(outfilename)
  , m_trackNodeName(trackNodeName)
  , m_maxEventDisplays(maxEventDisplays)
  , m_evt(0)
  , m_eventsSaved(0)
  , m_outfile(0)
  , m_tracks(0)
  , m_hits(0)
  , m_tpcGeom(0)
  , m_padMap(0)
{}

int InModuleTrackDisplay::Init(PHCompositeNode*)
{
  m_outfile = new TFile(m_outfilename.c_str(), "RECREATE");
  if (!m_outfile || m_outfile->IsZombie())
  {
    std::cerr << "InModuleTrackDisplay::Init - cannot open output file "
              << m_outfilename << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_outfile->mkdir("events");

  std::cout << "InModuleTrackDisplay::Init - writing only event multigraphs to "
            << m_outfilename << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

int InModuleTrackDisplay::process_event(PHCompositeNode* topNode)
{
  ++m_evt;

  if (!get_nodes(topNode))
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  const unsigned int ntracks = m_tracks ? m_tracks->size() : 0;
  const bool saveThisEvent = (m_eventsSaved < m_maxEventDisplays);

  if (!saveThisEvent)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  TDirectory* eventsTop = m_outfile->GetDirectory("events");
  if (!eventsTop)
  {
    eventsTop = m_outfile->mkdir("events");
  }

  eventsTop->cd();
  TDirectory* eventDir = eventsTop->mkdir(Form("event_%06u", m_evt));
  if (!eventDir)
  {
    std::cerr << "InModuleTrackDisplay::process_event - failed to create event directory"
              << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  eventDir->cd();

  std::map<GroupKey, GraphBundle> bundles;
  std::set<unsigned long long> associated_hit_ids;

  for (unsigned int itrk = 0; itrk < ntracks; ++itrk)
  {
    const InModuleTrack* trk = m_tracks->get_track(itrk);
    if (!trk) continue;

    const unsigned int tid = trk->get_track_id();
    const int side = static_cast<int>(trk->get_side());
    const int sector = static_cast<int>(trk->get_sector());
    const int col = track_color(itrk);

    std::map<GroupKey, std::vector<HitPoint> > points_by_group;

    for (unsigned int ih = 0; ih < trk->size_hit_indices(); ++ih)
    {
      const InModuleTrack::HitIndex idx = trk->get_hit_index(ih);
      const HitPoint p = make_hit_point(idx.first, idx.second);

      if (!p.ok) continue;
      associated_hit_ids.insert(make_unique_hit_id(p.hitsetkey, p.hitkey));

      const int module = layer_to_module(p.layer);
      if (module < 0) continue;

      const GroupKey key(side, sector, module);
      points_by_group[key].push_back(p);
    }

    std::map<GroupKey, std::vector<HitPoint> >::iterator pg;
    for (pg = points_by_group.begin(); pg != points_by_group.end(); ++pg)
    {
      const GroupKey& key = pg->first;
      const std::vector<HitPoint>& pts = pg->second;

      if (pts.empty()) continue;

      GraphBundle& b = get_bundle(bundles, key, m_evt, m_padMap);

      TGraph* g_tbin_pad_hits = new TGraph();
      g_tbin_pad_hits->SetName(Form("g_s%d_sec%02d_mod%d_trk%u_tbin_pad_hits",
                                    key.side, key.sector, key.module, tid));
      g_tbin_pad_hits->SetTitle(Form("event %u track %u hits;timebin;pad",
                                     m_evt, tid));
      style_hit_graph(g_tbin_pad_hits, col);

      TGraph* g_tbin_layer_hits = new TGraph();
      g_tbin_layer_hits->SetName(Form("g_s%d_sec%02d_mod%d_trk%u_tbin_layer_hits",
                                      key.side, key.sector, key.module, tid));
      g_tbin_layer_hits->SetTitle(Form("event %u track %u hits;timebin;layer",
                                       m_evt, tid));
      style_hit_graph(g_tbin_layer_hits, col);

      TGraph* g_pad_layer_hits = new TGraph();
      g_pad_layer_hits->SetName(Form("g_s%d_sec%02d_mod%d_trk%u_pad_layer_hits",
                                     key.side, key.sector, key.module, tid));
      g_pad_layer_hits->SetTitle(Form("event %u track %u hits;pad;layer",
                                      m_evt, tid));
      style_hit_graph(g_pad_layer_hits, col);

      unsigned int min_layer = pts[0].layer;
      unsigned int max_layer = pts[0].layer;
      double min_r = pts[0].r;
      double max_r = pts[0].r;

      for (unsigned int ip = 0; ip < pts.size(); ++ip)
      {
        const HitPoint& p = pts[ip];

        min_layer = std::min(min_layer, p.layer);
        max_layer = std::max(max_layer, p.layer);
        min_r = std::min(min_r, p.r);
        max_r = std::max(max_r, p.r);

        const int n = g_tbin_pad_hits->GetN();
        g_tbin_pad_hits->SetPoint(n,
                                  static_cast<double>(p.tbin),
                                  static_cast<double>(p.pad));

        g_tbin_layer_hits->SetPoint(n,
                                    static_cast<double>(p.tbin),
                                    static_cast<double>(p.layer));

        g_pad_layer_hits->SetPoint(n,
                                   static_cast<double>(p.pad),
                                   static_cast<double>(p.layer));

        const unsigned long long uid = make_unique_hit_id(p.hitsetkey, p.hitkey);
        if (b.filled_hit_ids.insert(uid).second && b.h3_adc_hits)
        {
          b.h3_adc_hits->Fill(static_cast<double>(p.tbin),
                              static_cast<double>(p.pad),
                              static_cast<double>(p.layer),
                              static_cast<double>(p.adc));

          if (b.h3_adc_globalphi_hits)
          {
            b.h3_adc_globalphi_hits->Fill(static_cast<double>(p.tbin),
                                      p.phi,
                                      p.r,
                                      static_cast<double>(p.adc));
          }
        }
      }

      b.mg_tbin_pad_hits->Add(g_tbin_pad_hits, "LP");
      b.mg_tbin_layer_hits->Add(g_tbin_layer_hits, "LP");
      b.mg_pad_layer_hits->Add(g_pad_layer_hits, "LP");

      const double l1 = static_cast<double>(min_layer);
      const double l2 = static_cast<double>(max_layer);

      const double pad1 = trk->get_pad_slope() * l1 + trk->get_pad_intercept();
      const double pad2 = trk->get_pad_slope() * l2 + trk->get_pad_intercept();

      const double tbin1 = trk->get_tbin_slope() * l1 + trk->get_tbin_intercept();
      const double tbin2 = trk->get_tbin_slope() * l2 + trk->get_tbin_intercept();

      TGraph* g_tbin_pad_fit = new TGraph();
      g_tbin_pad_fit->SetName(Form("g_s%d_sec%02d_mod%d_trk%u_tbin_pad_fit",
                                   key.side, key.sector, key.module, tid));
      g_tbin_pad_fit->SetTitle(Form("event %u track %u fit;timebin;pad",
                                    m_evt, tid));
      style_fit_graph(g_tbin_pad_fit, col);
      g_tbin_pad_fit->SetPoint(0, tbin1, pad1);
      g_tbin_pad_fit->SetPoint(1, tbin2, pad2);

      TGraph* g_tbin_layer_fit = new TGraph();
      g_tbin_layer_fit->SetName(Form("g_s%d_sec%02d_mod%d_trk%u_tbin_layer_fit",
                                     key.side, key.sector, key.module, tid));
      g_tbin_layer_fit->SetTitle(Form("event %u track %u fit;timebin;layer",
                                      m_evt, tid));
      style_fit_graph(g_tbin_layer_fit, col);
      g_tbin_layer_fit->SetPoint(0, tbin1, l1);
      g_tbin_layer_fit->SetPoint(1, tbin2, l2);

      TGraph* g_pad_layer_fit = new TGraph();
      g_pad_layer_fit->SetName(Form("g_s%d_sec%02d_mod%d_trk%u_pad_layer_fit",
                                    key.side, key.sector, key.module, tid));
      g_pad_layer_fit->SetTitle(Form("event %u track %u fit;pad;layer",
                                     m_evt, tid));
      style_fit_graph(g_pad_layer_fit, col);
      g_pad_layer_fit->SetPoint(0, pad1, l1);
      g_pad_layer_fit->SetPoint(1, pad2, l2);

      // Build 3D fit lines.
      // For field-on tracks, prefer the stored sagitta fit from the Python
      // notebook convention.  This avoids the unstable explicit circle
      // conversion for almost-straight tracks.  Fall back to the old circle
      // drawing only if no sagitta is stored, then fall back to straight lines.
      const double phi_width = m_padMap ? m_padMap->get_phi_bin_width(static_cast<unsigned int>(key.module)) : 0.0;

      double mean_phi = 0.0;
      double mean_pad = 0.0;
      double mean_tbin = 0.0;
      for (unsigned int ip = 0; ip < pts.size(); ++ip)
      {
        mean_phi  += pts[ip].local_phi;
        mean_pad  += static_cast<double>(pts[ip].pad);
        mean_tbin += static_cast<double>(pts[ip].tbin);
      }
      mean_phi  /= static_cast<double>(pts.size());
      mean_pad  /= static_cast<double>(pts.size());
      mean_tbin /= static_cast<double>(pts.size());

      const bool has_sagitta_fit =
        (trk->get_has_sagitta_fit() != 0 &&
         is_good_number(trk->get_sagitta_S()) &&
         is_good_number(trk->get_sagitta_x0()) &&
         is_good_number(trk->get_sagitta_invR()) &&
         is_good_number(trk->get_sagitta_theta()) &&
         is_good_number(trk->get_sagitta_b()) &&
         phi_width > 0.0);

      bool drew_curved_fit = false;

      if (has_sagitta_fit)
      {
        std::vector<double> v_tbin;
        std::vector<double> v_pad;
        std::vector<double> v_layer;
        std::vector<double> v_phi;
        std::vector<double> v_r;

        for (int il = static_cast<int>(min_layer); il <= static_cast<int>(max_layer); ++il)
        {
          const unsigned int layer = static_cast<unsigned int>(il);
          const double rloc = m_padMap->get_local_radius(static_cast<unsigned int>(key.module), layer);
          if (rloc <= 0.0) continue;

          double philoc = 0.0;
          if (!sagitta_phi_at_radius_display(rloc,
                                             mean_phi,
                                             trk->get_sagitta_S(),
                                             trk->get_sagitta_x0(),
                                             trk->get_sagitta_invR(),
                                             trk->get_sagitta_theta(),
                                             trk->get_sagitta_b(),
                                             philoc))
          {
            continue;
          }

          const double padloc = philoc / phi_width;
          const double tloc =
            trk->get_radius_tbin_slope() * rloc +
            trk->get_radius_tbin_intercept();

          if (!is_good_number(tloc) || !is_good_number(padloc) || !is_good_number(philoc)) continue;

          v_tbin.push_back(tloc);
          v_pad.push_back(local_phi_to_full_pad(m_padMap, key.sector, static_cast<unsigned int>(key.module), philoc));
          v_layer.push_back(static_cast<double>(layer));
          v_phi.push_back(local_phi_to_global_phi(m_padMap, key.side, key.sector, static_cast<unsigned int>(key.module), philoc));
          v_r.push_back(rloc);
        }

        if (v_tbin.size() >= 2)
        {
          TPolyLine3D* fit_line_3d = new TPolyLine3D(v_tbin.size());
          const std::string line_name =
            Form("line3_s%d_sec%02d_mod%d_trk%u_sagitta_fit",
                 key.side, key.sector, key.module, tid);

          for (unsigned int ip = 0; ip < v_tbin.size(); ++ip)
            fit_line_3d->SetPoint(ip, v_tbin[ip], v_pad[ip], v_layer[ip]);

          style_fit_line_3d(fit_line_3d, col);
          b.fit_lines_3d.push_back(fit_line_3d);
          b.fit_line_names_3d.push_back(line_name);

          TPolyLine3D* globalphi_fit_line_3d = new TPolyLine3D(v_tbin.size());
          const std::string globalphi_line_name =
            Form("line3_s%d_sec%02d_mod%d_trk%u_globalphi_sagitta_fit",
                 key.side, key.sector, key.module, tid);

          for (unsigned int ip = 0; ip < v_tbin.size(); ++ip)
            globalphi_fit_line_3d->SetPoint(ip, v_tbin[ip], v_phi[ip], v_r[ip]);

          style_fit_line_3d(globalphi_fit_line_3d, col);
          b.globalphi_fit_lines_3d.push_back(globalphi_fit_line_3d);
          b.globalphi_fit_line_names_3d.push_back(globalphi_line_name);

          drew_curved_fit = true;
        }
        else
        {
          std::cout << "InModuleTrackDisplay: no drawable sagitta segment for track "
                    << tid << " side=" << key.side
                    << " sector=" << key.sector
                    << " module=" << key.module
                    << " S=" << trk->get_sagitta_S()
                    << " x0=" << trk->get_sagitta_x0()
                    << " invR=" << trk->get_sagitta_invR()
                    << " theta=" << trk->get_sagitta_theta()
                    << " b=" << trk->get_sagitta_b()
                    << std::endl;
        }
      }

      const bool has_circle_fit =
        (!drew_curved_fit &&
         trk->get_circle_radius() > 0.0 &&
         is_good_number(trk->get_circle_radius()) &&
         is_good_number(trk->get_circle_x0()) &&
         is_good_number(trk->get_circle_y0()) &&
         phi_width > 0.0);

      if (has_circle_fit)
      {
        const double x0 = trk->get_circle_x0();
        const double y0 = trk->get_circle_y0();
        const double rc = trk->get_circle_radius();
        const double d = std::sqrt(x0 * x0 + y0 * y0);

        std::vector<double> v_tbin;
        std::vector<double> v_pad;
        std::vector<double> v_layer;
        std::vector<double> v_phi;
        std::vector<double> v_r;

        for (int il = static_cast<int>(min_layer); il <= static_cast<int>(max_layer); ++il)
        {
          const unsigned int layer = static_cast<unsigned int>(il);
          const double rloc = m_padMap->get_local_radius(static_cast<unsigned int>(key.module), layer);
          if (rloc <= 0.0) continue;
          if (d <= 1.0e-12) continue;
          if (d > rloc + rc) continue;
          if (d < std::fabs(rloc - rc)) continue;

          const double a = (rloc * rloc - rc * rc + d * d) / (2.0 * d);
          double h2 = rloc * rloc - a * a;
          if (h2 < 0.0 && h2 > -1.0e-8) h2 = 0.0;
          if (h2 < 0.0) continue;

          const double h = std::sqrt(h2);
          const double ux = x0 / d;
          const double uy = y0 / d;

          const double xb = a * ux;
          const double yb = a * uy;

          const double xint1 = xb - h * uy;
          const double yint1 = yb + h * ux;
          const double xint2 = xb + h * uy;
          const double yint2 = yb - h * ux;

          double phi1 = std::atan2(yint1, xint1);
          double phi2 = std::atan2(yint2, xint2);
          if (phi1 < 0.0) phi1 += 2.0 * TMath::Pi();
          if (phi2 < 0.0) phi2 += 2.0 * TMath::Pi();

          const double pad_seed =
            trk->get_pad_slope() * static_cast<double>(layer) +
            trk->get_pad_intercept();
          double seed = pad_seed * phi_width;
          if (!is_good_number(seed)) seed = mean_phi;

          double dphi1 = std::fabs(wrap_delta_phi_display(phi1 - seed));
          double dphi2 = std::fabs(wrap_delta_phi_display(phi2 - seed));

          const double philoc = (dphi1 <= dphi2) ? phi1 : phi2;
          const double padloc = philoc / phi_width;
          const double tloc =
            trk->get_radius_tbin_slope() * rloc +
            trk->get_radius_tbin_intercept();

          if (!is_good_number(tloc) || !is_good_number(padloc) || !is_good_number(philoc)) continue;

          v_tbin.push_back(tloc);
          v_pad.push_back(local_phi_to_full_pad(m_padMap, key.sector, static_cast<unsigned int>(key.module), philoc));
          v_layer.push_back(static_cast<double>(layer));
          v_phi.push_back(local_phi_to_global_phi(m_padMap, key.side, key.sector, static_cast<unsigned int>(key.module), philoc));
          v_r.push_back(rloc);
        }

        if (v_tbin.size() >= 2)
        {
          TPolyLine3D* fit_line_3d = new TPolyLine3D(v_tbin.size());
          const std::string line_name =
            Form("line3_s%d_sec%02d_mod%d_trk%u_circle_fit",
                 key.side, key.sector, key.module, tid);

          for (unsigned int ip = 0; ip < v_tbin.size(); ++ip)
            fit_line_3d->SetPoint(ip, v_tbin[ip], v_pad[ip], v_layer[ip]);

          style_fit_line_3d(fit_line_3d, col);
          b.fit_lines_3d.push_back(fit_line_3d);
          b.fit_line_names_3d.push_back(line_name);

          TPolyLine3D* globalphi_fit_line_3d = new TPolyLine3D(v_tbin.size());
          const std::string globalphi_line_name =
            Form("line3_s%d_sec%02d_mod%d_trk%u_globalphi_circle_fit",
                 key.side, key.sector, key.module, tid);

          for (unsigned int ip = 0; ip < v_tbin.size(); ++ip)
            globalphi_fit_line_3d->SetPoint(ip, v_tbin[ip], v_phi[ip], v_r[ip]);

          style_fit_line_3d(globalphi_fit_line_3d, col);
          b.globalphi_fit_lines_3d.push_back(globalphi_fit_line_3d);
          b.globalphi_fit_line_names_3d.push_back(globalphi_line_name);

          drew_curved_fit = true;
        }
        else
        {
          std::cout << "InModuleTrackDisplay: no drawable circle segment for track "
                    << tid << " side=" << key.side
                    << " sector=" << key.sector
                    << " module=" << key.module
                    << " cx=" << x0 << " cy=" << y0 << " R=" << rc
                    << " d=" << d
                    << " min_layer=" << min_layer
                    << " max_layer=" << max_layer
                    << std::endl;
        }
      }

      if (!drew_curved_fit)
      {
        TPolyLine3D* fit_line_3d = new TPolyLine3D(2);
        const std::string line_name =
          Form("line3_s%d_sec%02d_mod%d_trk%u_line_fit",
               key.side, key.sector, key.module, tid);

        fit_line_3d->SetPoint(0, tbin1, pad1, l1);
        fit_line_3d->SetPoint(1, tbin2, pad2, l2);
        style_fit_line_3d(fit_line_3d, col);
        b.fit_lines_3d.push_back(fit_line_3d);
        b.fit_line_names_3d.push_back(line_name);

        TPolyLine3D* globalphi_fit_line_3d = new TPolyLine3D(51);
        const std::string globalphi_line_name =
          Form("line3_s%d_sec%02d_mod%d_trk%u_globalphi_line_fit",
               key.side, key.sector, key.module, tid);

        for (int iplocal = 0; iplocal <= 50; ++iplocal)
        {
          const double frac = static_cast<double>(iplocal) / 50.0;
          const double rloc = min_r + frac * (max_r - min_r);
          const double tloc =
            trk->get_radius_tbin_slope() * rloc +
            trk->get_radius_tbin_intercept();
          const double philoc =
            trk->get_radius_phi_slope() * rloc +
            trk->get_radius_phi_intercept();

          globalphi_fit_line_3d->SetPoint(iplocal,
                                      tloc,
                                      local_phi_to_global_phi(m_padMap, key.side, key.sector, static_cast<unsigned int>(key.module), philoc),
                                      rloc);
        }

        style_fit_line_3d(globalphi_fit_line_3d, col);
        b.globalphi_fit_lines_3d.push_back(globalphi_fit_line_3d);
        b.globalphi_fit_line_names_3d.push_back(globalphi_line_name);
      }

      b.mg_tbin_pad_fits->Add(g_tbin_pad_fit, "L");
      b.mg_tbin_layer_fits->Add(g_tbin_layer_fit, "L");
      b.mg_pad_layer_fits->Add(g_pad_layer_fit, "L");
    }
  }

      // Fill all TRKR_HITSET hits that are not associated with any reconstructed in-module track
    TrkrHitSetContainer::ConstRange hitset_range = m_hits->getHitSets();

    for (TrkrHitSetContainer::ConstIterator hsiter = hitset_range.first;
        hsiter != hitset_range.second;
        ++hsiter)
    {
      TrkrHitSet* hitset = hsiter->second;
      if (!hitset) continue;

      const TrkrDefs::hitsetkey hsk = hsiter->first;

      const unsigned int layer = TrkrDefs::getLayer(hsk);
      const int module = layer_to_module(layer);
      if (module < 0) continue;

      const int side = static_cast<int>(TpcDefs::getSide(hsk));
      const int sector = static_cast<int>(TpcDefs::getSectorId(hsk));

      const GroupKey key(side, sector, module);
      GraphBundle& b = get_bundle(bundles, key, m_evt, m_padMap);

      TrkrHitSet::ConstRange hit_range = hitset->getHits();

      for (TrkrHitSet::ConstIterator hiter = hit_range.first;
          hiter != hit_range.second;
          ++hiter)
      {
        const TrkrDefs::hitkey hk = hiter->first;
        TrkrHit* hit = hiter->second;
        if (!hit) continue;

        const unsigned long long uid = make_unique_hit_id(hsk, hk);

        // Skip hits already used by any track
        if (associated_hit_ids.find(uid) != associated_hit_ids.end()) continue;

        const unsigned int pad = TpcDefs::getPad(hk);
        const unsigned int tbin = TpcDefs::getTBin(hk);
        const unsigned short adc = hit->getAdc();

        if (b.h3_adc_unassociated_hits)
        {
          b.h3_adc_unassociated_hits->Fill(static_cast<double>(tbin),
                                          static_cast<double>(pad),
                                          static_cast<double>(layer),
                                          static_cast<double>(adc));
        }
      }
    }

  std::map<GroupKey, GraphBundle>::iterator ib;
  for (ib = bundles.begin(); ib != bundles.end(); ++ib)
  {
    eventDir->cd();
    write_bundle(ib->second, ib->first);
  }

  std::cout << "InModuleTrackDisplay - saved event "
            << m_evt << " with " << ntracks
            << " tracks and " << bundles.size()
            << " side/sector/module graph groups"
            << std::endl;

  ++m_eventsSaved;

  return Fun4AllReturnCodes::EVENT_OK;
}

int InModuleTrackDisplay::End(PHCompositeNode*)
{
  if (m_outfile)
  {
    std::cout << "InModuleTrackDisplay::End - events seen: "
              << m_evt << std::endl;

    m_outfile->cd();
    m_outfile->Write();
    m_outfile->Close();

    delete m_outfile;
    m_outfile = 0;
  }

  std::cout << "InModuleTrackDisplay::End - wrote "
            << m_outfilename << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

bool InModuleTrackDisplay::get_nodes(PHCompositeNode* topNode)
{
  m_hits = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  m_tpcGeom = findNode::getClass<PHG4TpcGeomContainer>(topNode, "TPCGEOMCONTAINER");
  m_padMap = findNode::getClass<TpcPadMap>(topNode, "TPC_PADMAP");

  m_tracks = findNode::getClass<InModuleTrackContainer>(topNode,
                                                        m_trackNodeName.c_str());

  if (!m_tracks)
  {
    const char* candidate_names[] = {
      "INMODULETRACKS",
      "InModuleTracks",
      "InModuleTrackContainer",
      "INMODULETRACKCONTAINER",
      "INMODULETRACKS_CONTAINER"
    };

    for (unsigned int i = 0;
         i < sizeof(candidate_names) / sizeof(candidate_names[0]) && !m_tracks;
         ++i)
    {
      m_tracks = findNode::getClass<InModuleTrackContainer>(topNode,
                                                            candidate_names[i]);

      if (m_tracks)
      {
        m_trackNodeName = candidate_names[i];
      }
    }
  }

  if (!m_tracks)
  {
    std::cerr << "InModuleTrackDisplay - could not find InModuleTrackContainer node"
              << std::endl;
    return false;
  }

  if (!m_hits)
  {
    std::cerr << "InModuleTrackDisplay - missing TRKR_HITSET"
              << std::endl;
    return false;
  }

  if (!m_padMap)
  {
    std::cerr << "InModuleTrackDisplay - missing TPC_PADMAP"
              << std::endl;
    return false;
  }

  if (!m_padMap->isValid())
  {
    std::cerr << "InModuleTrackDisplay - TPC_PADMAP is invalid"
              << std::endl;
    return false;
  }

  return true;
}

InModuleTrackDisplay::HitPoint
InModuleTrackDisplay::make_hit_point(const TrkrDefs::hitsetkey hsk,
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
  p.pad = TpcDefs::getPad(hk);
  p.tbin = TpcDefs::getTBin(hk);
  p.adc = hit->getAdc();

  const int module = layer_to_module(p.layer);
  if (module < 0) return p;

  if (!m_padMap) return p;

  const unsigned int region = static_cast<unsigned int>(module);
  const unsigned int side = static_cast<unsigned int>(TpcDefs::getSide(hsk));

  p.r = m_padMap->get_local_radius(region, p.layer);
  p.local_phi = m_padMap->get_local_phi(region, static_cast<double>(p.pad));
  p.phi = m_padMap->get_global_phi(side, region, static_cast<double>(p.pad));

  p.x = p.r * std::cos(p.phi);
  p.y = p.r * std::sin(p.phi);

  p.ok = true;
  return p;
}
