#include "InModuleTrackDisplay.h"

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
    if (layer < 7) return -1;
    return static_cast<int>((layer - 7) / 16);
  }

  unsigned long long make_unique_hit_id(const TrkrDefs::hitsetkey hsk,
                                        const TrkrDefs::hitkey hk)
  {
    return (static_cast<unsigned long long>(hsk) << 32) |
           static_cast<unsigned long long>(hk);
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
      , c3_adc_hits_fits(0)
      , c3_adc_unassociated_hits(0)
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

    TCanvas* c3_adc_hits_fits;
    TCanvas* c3_adc_unassociated_hits;
    std::vector<TPolyLine3D*> fit_lines_3d;
    std::vector<std::string> fit_line_names_3d;

    // avoid filling the same hit more than once if tracks share hits
    std::set<unsigned long long> filled_hit_ids;
  };

  void make_bundle_graphs(GraphBundle& b,
                          const GroupKey& key,
                          const unsigned int evt)
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
    unsigned int Npads[3] = {94, 128, 192};
    const int pad_min =  Npads[key.module]*key.sector;
    const int pad_max =  pad_min + (Npads[key.module]-1);

    b.h3_adc_hits = new TH3D(Form("h3_%s_adc_hits", tag.Data()),
                             Form("event %u side %d sector %d module %d: ADC-weighted hits;timebin;pad;layer",
                                  evt, key.side, key.sector, key.module),
                             512, -0.5, 511.5,
                             Npads[key.module], static_cast<double>(pad_min) - 0.5, static_cast<double>(pad_max) + 0.5,
                             16, static_cast<double>(layer_min) - 0.5,
                             static_cast<double>(layer_max) + 0.5);

    b.h3_adc_hits->SetStats(0);


    b.h3_adc_unassociated_hits =
    new TH3D(Form("h3_%s_adc_unassociated_hits", tag.Data()),
            Form("event %u side %d sector %d module %d: ADC-weighted unassociated TRKR_HITSET hits;timebin;pad;layer",
                  evt, key.side, key.sector, key.module),
            512, -0.5, 511.5,
             Npads[key.module], static_cast<double>(pad_min) - 0.5, static_cast<double>(pad_max) + 0.5,
            16, static_cast<double>(layer_min) - 0.5,
            static_cast<double>(layer_max) + 0.5);

    b.h3_adc_unassociated_hits->SetStats(0);


  }


  GraphBundle& get_bundle(std::map<GroupKey, GraphBundle>& bundles,
                          const GroupKey& key,
                          const unsigned int evt)
  {
    std::map<GroupKey, GraphBundle>::iterator it = bundles.find(key);

    if (it == bundles.end())
    {
      const std::pair<std::map<GroupKey, GraphBundle>::iterator, bool> inserted =
        bundles.insert(std::make_pair(key, GraphBundle()));

      it = inserted.first;
      make_bundle_graphs(it->second, key, evt);
    }

    return it->second;
  }

  void write_bundle(GraphBundle& b, const GroupKey& key)
  {
    if (b.h3_adc_hits) b.h3_adc_hits->Write();
    if (b.h3_adc_unassociated_hits) b.h3_adc_unassociated_hits->Write();

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

    if (b.mg_tbin_pad_hits) b.mg_tbin_pad_hits->Write();
    //if (b.mg_tbin_layer_hits) b.mg_tbin_layer_hits->Write();
    //if (b.mg_pad_layer_hits) b.mg_pad_layer_hits->Write();

    if (b.mg_tbin_pad_fits) b.mg_tbin_pad_fits->Write();
    if (b.mg_tbin_layer_fits) b.mg_tbin_layer_fits->Write();
    if (b.mg_pad_layer_fits) b.mg_pad_layer_fits->Write();
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

      GraphBundle& b = get_bundle(bundles, key, m_evt);

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

      for (unsigned int ip = 0; ip < pts.size(); ++ip)
      {
        const HitPoint& p = pts[ip];

        min_layer = std::min(min_layer, p.layer);
        max_layer = std::max(max_layer, p.layer);

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

      TPolyLine3D* fit_line_3d = new TPolyLine3D(2);
      //fit_line_3d->SetName(Form("line3_s%d_sec%02d_mod%d_trk%u_fit",
      //                          key.side, key.sector, key.module, tid));
      const std::string line_name =  Form("line3_s%d_sec%02d_mod%d_trk%u_fit", key.side, key.sector, key.module, tid);
      fit_line_3d->SetPoint(0, tbin1, pad1, l1);
      fit_line_3d->SetPoint(1, tbin2, pad2, l2);
      style_fit_line_3d(fit_line_3d, col);
      b.fit_lines_3d.push_back(fit_line_3d);
      b.fit_line_names_3d.push_back(line_name);

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
      GraphBundle& b = get_bundle(bundles, key, m_evt);

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

  if (!m_tpcGeom)
  {
    std::cerr << "InModuleTrackDisplay - missing TPCGEOMCONTAINER"
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

  PHG4TpcGeom* layergeom = m_tpcGeom->GetLayerCellGeom(p.layer);
  if (!layergeom) return p;

  const unsigned int side = TpcDefs::getSide(hsk);

  p.r = layergeom->get_radius();
  p.phi = layergeom->get_phicenter(p.pad, side);

  p.x = p.r * std::cos(p.phi);
  p.y = p.r * std::sin(p.phi);

  p.ok = true;
  return p;
}
