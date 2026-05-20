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

#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4detectors/PHG4TpcGeom.h>
#include <g4detectors/PHG4TpcGeomContainer.h>

#include <TDirectory.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TMultiGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TMath.h>
#include <TColor.h>
#include <TString.h>

#include <cmath>
#include <iostream>

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

  void style_track_graph(TGraph* g, const int color)
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
}

InModuleTrackDisplay::HitPoint::HitPoint()
  : ok(false)
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
  , h_ntracks(0)
  , h_nhits_per_track(0)
  , h_nrawhits_field(0)
  , h_nblobs(0)
  , h_chi2pad_ndf(0)
  , h_chi2tbin_ndf(0)
  , h_layer_first_last(0)
  , h_pad_slope_vs_tbin_slope(0)
  , h_sector_side(0)
  , h_xy_all(0)
  , h_pad_tbin_all(0)
  , h_layer_tbin_all(0)
  , h_layer_pad_all(0)
  , h_xyz_all(0)
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

  m_outfile->mkdir("summary");
  m_outfile->mkdir("events");
  m_outfile->cd("summary");

  h_ntracks = new TH1D("h_ntracks",
                       "In-module tracks per event;N tracks;Events",
                       300, 0, 300);

  h_nhits_per_track = new TH1D("h_nhits_per_track",
                               "Referenced raw hits per track;N hit indices;Tracks",
                               300, 0, 300);

  h_nrawhits_field = new TH1D("h_nrawhits_field",
                              "track->get_nrawhits();N raw hits;Tracks",
                              300, 0, 300);

  h_nblobs = new TH1D("h_nblobs",
                      "Blobs per track;N blobs;Tracks",
                      100, 0, 100);

  h_layer_first_last = new TH2D("h_layer_first_last",
                                "First vs last layer;first layer;last layer;Tracks",
                                60, 0, 60, 60, 0, 60);

  h_chi2pad_ndf = new TH1D("h_chi2pad_ndf",
                           "Pad fit #chi^{2}/ndf;#chi^{2}/ndf;Tracks",
                           200, 0, 50);

  h_chi2tbin_ndf = new TH1D("h_chi2tbin_ndf",
                            "Tbin fit #chi^{2}/ndf;#chi^{2}/ndf;Tracks",
                            200, 0, 50);

  h_pad_slope_vs_tbin_slope = new TH2D("h_pad_slope_vs_tbin_slope",
                                       "Fit slopes;tbin slope;pad slope;Tracks",
                                       200, -20, 20, 200, -20, 20);

  h_sector_side = new TH2D("h_sector_side",
                           "Track module side/sector;sector;side;Tracks",
                           24, 0, 24, 5, -2, 3);

  h_xy_all = new TH2D("h_xy_all",
                      "All referenced track hits, x-y;x [cm];y [cm];ADC-weighted hits",
                      800, -90, 90, 800, -90, 90);

  h_pad_tbin_all = new TH2D("h_pad_tbin_all",
                            "All referenced track hits, pad-tbin;pad;tbin;ADC-weighted hits",
                            3000, 0, 3000, 3000, 0, 3000);

  h_layer_tbin_all = new TH2D("h_layer_tbin_all",
                              "All referenced track hits, layer-tbin;layer;tbin;ADC-weighted hits",
                              60, 0, 60, 3000, 0, 3000);

  h_layer_pad_all = new TH2D("h_layer_pad_all",
                             "All referenced track hits, layer-pad;layer;pad;ADC-weighted hits",
                             60, 0, 60, 3000, 0, 3000);

  h_xyz_all = new TH3D("h_xyz_all",
                       "All referenced track hits, x-y-tbin;x [cm];y [cm];tbin",
                       240, -90, 90, 240, -90, 90, 300, 0, 3000);

  std::cout << "InModuleTrackDisplay::Init - writing to "
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
  h_ntracks->Fill(ntracks);

  const bool saveThisEvent = (m_eventsSaved < m_maxEventDisplays);

  TDirectory* eventDir = 0;

  TH2D* h_xy_evt = 0;
  TH2D* h_xy_trackid_evt = 0;
  TH2D* h_pad_tbin_evt = 0;
  TH2D* h_layer_tbin_evt = 0;
  TH2D* h_layer_pad_evt = 0;
  TH3D* h_xyz_evt = 0;

  TMultiGraph* mg_xy_tracks = 0;
  TMultiGraph* mg_pad_tbin_tracks = 0;
  TMultiGraph* mg_layer_pad_tracks = 0;
  TMultiGraph* mg_layer_tbin_tracks = 0;

  TMultiGraph* mg_layer_pad_fits = 0;
  TMultiGraph* mg_layer_tbin_fits = 0;

  if (saveThisEvent)
  {
    TDirectory* eventsTop = m_outfile->GetDirectory("events");
    eventsTop->cd();

    eventDir = eventsTop->mkdir(Form("event_%06u", m_evt));
    eventDir->cd();

    h_xy_evt = new TH2D("h_xy",
                        Form("event %u: x-y;x [cm];y [cm]", m_evt),
                        800, -90, 90, 800, -90, 90);

    h_xy_trackid_evt = new TH2D("h_xy_trackid",
                                Form("event %u: x-y weighted by track id;x [cm];y [cm]", m_evt),
                                800, -90, 90, 800, -90, 90);

    h_pad_tbin_evt = new TH2D("h_pad_tbin",
                              Form("event %u: pad-tbin;pad;tbin", m_evt),
                              3000, 0, 3000, 3000, 0, 3000);

    h_layer_tbin_evt = new TH2D("h_layer_tbin",
                                Form("event %u: layer-tbin;layer;tbin", m_evt),
                                60, 0, 60, 3000, 0, 3000);

    h_layer_pad_evt = new TH2D("h_layer_pad",
                               Form("event %u: layer-pad;layer;pad", m_evt),
                               60, 0, 60, 3000, 0, 3000);

    h_xyz_evt = new TH3D("h_xyz",
                         Form("event %u: x-y-tbin;x [cm];y [cm];tbin", m_evt),
                         240, -90, 90, 240, -90, 90, 300, 0, 3000);

    mg_xy_tracks = new TMultiGraph();
    mg_xy_tracks->SetName("mg_xy_tracks");
    mg_xy_tracks->SetTitle(Form("event %u: colored tracks, x-y;x [cm];y [cm]", m_evt));

    mg_pad_tbin_tracks = new TMultiGraph();
    mg_pad_tbin_tracks->SetName("mg_pad_tbin_tracks");
    mg_pad_tbin_tracks->SetTitle(Form("event %u: colored tracks, pad-tbin;pad;tbin", m_evt));

    mg_layer_pad_tracks = new TMultiGraph();
    mg_layer_pad_tracks->SetName("mg_layer_pad_tracks");
    mg_layer_pad_tracks->SetTitle(Form("event %u: colored hits, layer-pad;layer;pad", m_evt));

    mg_layer_tbin_tracks = new TMultiGraph();
    mg_layer_tbin_tracks->SetName("mg_layer_tbin_tracks");
    mg_layer_tbin_tracks->SetTitle(Form("event %u: colored hits, layer-tbin;layer;tbin", m_evt));

    mg_layer_pad_fits = new TMultiGraph();
    mg_layer_pad_fits->SetName("mg_layer_pad_fits");
    mg_layer_pad_fits->SetTitle(Form("event %u: colored pad fit lines;layer;pad", m_evt));

    mg_layer_tbin_fits = new TMultiGraph();
    mg_layer_tbin_fits->SetName("mg_layer_tbin_fits");
    mg_layer_tbin_fits->SetTitle(Form("event %u: colored tbin fit lines;layer;tbin", m_evt));
  }

  for (unsigned int itrk = 0; itrk < ntracks; ++itrk)
  {
    const InModuleTrack* trk = m_tracks->get_track(itrk);
    if (!trk) continue;

    const unsigned int tid = trk->get_track_id();

    h_nhits_per_track->Fill(trk->size_hit_indices());
    h_nrawhits_field->Fill(trk->get_nrawhits());
    h_nblobs->Fill(trk->get_nblobs());
    h_layer_first_last->Fill(trk->get_first_layer(), trk->get_last_layer());
    h_sector_side->Fill(trk->get_sector(), trk->get_side());

    if (trk->get_ndof_pad() > 0)
    {
      h_chi2pad_ndf->Fill(trk->get_chi2_pad() / static_cast<double>(trk->get_ndof_pad()));
    }

    if (trk->get_ndof_tbin() > 0)
    {
      h_chi2tbin_ndf->Fill(trk->get_chi2_tbin() / static_cast<double>(trk->get_ndof_tbin()));
    }

    h_pad_slope_vs_tbin_slope->Fill(trk->get_tbin_slope(), trk->get_pad_slope());

    TGraph* gxy = 0;
    TGraph* gpad_tbin = 0;
    TGraph* glayer_pad = 0;
    TGraph* glayer_tbin = 0;
    TGraph* gfit_layer_pad = 0;
    TGraph* gfit_layer_tbin = 0;


    if (saveThisEvent)
    {
      eventDir->cd();

      const int col = track_color(itrk);

      gxy = new TGraph();
      gxy->SetName(Form("gxy_track_%u", tid));
      gxy->SetTitle(Form("event %u track %u;x [cm];y [cm]", m_evt, tid));
      style_track_graph(gxy, col);

      gpad_tbin = new TGraph();
      gpad_tbin->SetName(Form("gpad_tbin_track_%u", tid));
      gpad_tbin->SetTitle(Form("event %u track %u;pad;tbin", m_evt, tid));
      style_track_graph(gpad_tbin, col);

      glayer_pad = new TGraph();
      glayer_pad->SetName(Form("glayer_pad_track_%u", tid));
      glayer_pad->SetTitle(Form("event %u track %u;layer;pad", m_evt, tid));
      style_track_graph(glayer_pad, col);

      glayer_tbin = new TGraph();
      glayer_tbin->SetName(Form("glayer_tbin_track_%u", tid));
      glayer_tbin->SetTitle(Form("event %u track %u;layer;tbin", m_evt, tid));
      style_track_graph(glayer_tbin, col);



      const double firstLayer = static_cast<double>(trk->get_first_layer());
      const double lastLayer = static_cast<double>(trk->get_last_layer());

      const double pad1 = trk->get_pad_slope() * firstLayer + trk->get_pad_intercept();
      const double pad2 = trk->get_pad_slope() * lastLayer + trk->get_pad_intercept();

      const double tbin1 = trk->get_tbin_slope() * firstLayer + trk->get_tbin_intercept();
      const double tbin2 = trk->get_tbin_slope() * lastLayer + trk->get_tbin_intercept();

      gfit_layer_pad = new TGraph();
      gfit_layer_pad->SetName(Form("gfit_layer_pad_track_%u", tid));
      gfit_layer_pad->SetTitle(Form("event %u track %u pad fit;layer;pad", m_evt, tid));
      style_fit_graph(gfit_layer_pad, col);
      gfit_layer_pad->SetPoint(0, firstLayer, pad1);
      gfit_layer_pad->SetPoint(1, lastLayer, pad2);

      gfit_layer_tbin = new TGraph();
      gfit_layer_tbin->SetName(Form("gfit_layer_tbin_track_%u", tid));
      gfit_layer_tbin->SetTitle(Form("event %u track %u tbin fit;layer;tbin", m_evt, tid));
      style_fit_graph(gfit_layer_tbin, col);
      gfit_layer_tbin->SetPoint(0, firstLayer, tbin1);
      gfit_layer_tbin->SetPoint(1, lastLayer, tbin2);
    }

    for (unsigned int ih = 0; ih < trk->size_hit_indices(); ++ih)
    {
      const InModuleTrack::HitIndex idx = trk->get_hit_index(ih);
      const HitPoint p = make_hit_point(idx.first, idx.second);

      if (!p.ok) continue;

      h_xy_all->Fill(p.x, p.y, p.adc);
      h_pad_tbin_all->Fill(p.pad, p.tbin, p.adc);
      h_layer_tbin_all->Fill(p.layer, p.tbin, p.adc);
      h_layer_pad_all->Fill(p.layer, p.pad, p.adc);
      h_xyz_all->Fill(p.x, p.y, p.tbin, p.adc);

      if (saveThisEvent)
      {
        h_xy_evt->Fill(p.x, p.y, p.adc);
        h_xy_trackid_evt->Fill(p.x, p.y, tid + 1.0);
        h_pad_tbin_evt->Fill(p.pad, p.tbin, p.adc);
        h_layer_tbin_evt->Fill(p.layer, p.tbin, p.adc);
        h_layer_pad_evt->Fill(p.layer, p.pad, p.adc);
        h_xyz_evt->Fill(p.x, p.y, p.tbin, p.adc);

        const int n = gxy->GetN();

        gxy->SetPoint(n, p.x, p.y);
        gpad_tbin->SetPoint(n, p.pad, p.tbin);
        glayer_pad->SetPoint(n, p.layer, p.pad);
        glayer_tbin->SetPoint(n, p.layer, p.tbin);

      }
    }

    if (saveThisEvent)
    {
      if (gxy && gxy->GetN() > 0)
      {
        mg_xy_tracks->Add(gxy, "LP");
      }

      if (gpad_tbin && gpad_tbin->GetN() > 0)
      {
        mg_pad_tbin_tracks->Add(gpad_tbin, "LP");
      }

      if (glayer_pad && glayer_pad->GetN() > 0)
      {
        mg_layer_pad_tracks->Add(glayer_pad, "LP");
      }

      if (glayer_tbin && glayer_tbin->GetN() > 0)
      {
        mg_layer_tbin_tracks->Add(glayer_tbin, "LP");
      }

      if (gfit_layer_pad)
      {
        mg_layer_pad_fits->Add(gfit_layer_pad, "L");
      }

      if (gfit_layer_tbin)
      {
        mg_layer_tbin_fits->Add(gfit_layer_tbin, "L");
      }
    }
  }

  if (saveThisEvent)
  {
    eventDir->cd();

    if (mg_xy_tracks) mg_xy_tracks->Write();
    if (mg_pad_tbin_tracks) mg_pad_tbin_tracks->Write();
    if (mg_layer_pad_tracks) mg_layer_pad_tracks->Write();
    if (mg_layer_tbin_tracks) mg_layer_tbin_tracks->Write();

    if (mg_layer_pad_fits) mg_layer_pad_fits->Write();
    if (mg_layer_tbin_fits) mg_layer_tbin_fits->Write();

    std::cout << "InModuleTrackDisplay - saving event "
              << m_evt << " with " << ntracks << " tracks"
              << std::endl;

    eventDir->Write();
    ++m_eventsSaved;
  }

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

  m_tracks = findNode::getClass<InModuleTrackContainer>(topNode, m_trackNodeName.c_str());

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
      m_tracks = findNode::getClass<InModuleTrackContainer>(topNode, candidate_names[i]);

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