#include "TpcPolyClusterDisplay.h"

#include "TpcPolyClusterTrack.h"
#include "TpcPolyClusterTrackContainer.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH3D.h>
#include <TPolyMarker3D.h>
#include <TString.h>

#include <cmath>
#include <iostream>
#include <vector>

namespace
{
  int cluster_color(const unsigned int icluster)
  {
    static const int colors[] = {
      kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kCyan + 2,
      kOrange + 7, kViolet + 1, kAzure + 1, kPink + 7, kTeal + 3
    };
    return colors[icluster % (sizeof(colors) / sizeof(colors[0]))];
  }

  TPolyMarker3D* make_cluster_marker(const double z,
                                     const double x,
                                     const double y,
                                     const int color)
  {
    TPolyMarker3D* marker = new TPolyMarker3D(1);
    marker->SetPoint(0, z, x, y);
    marker->SetMarkerColor(color);
    marker->SetMarkerStyle(20);
    marker->SetMarkerSize(1.4);
    return marker;
  }
}

TpcPolyClusterDisplay::TpcPolyClusterDisplay(const std::string& name,
                                             const std::string& outfilename,
                                             const std::string& clusterNodeName,
                                             const unsigned int maxEventDisplays)
  : SubsysReco(name)
  , m_outfilename(outfilename)
  , m_clusterNodeName(clusterNodeName)
  , m_maxEventDisplays(maxEventDisplays)
  , m_evt(0)
  , m_eventsSaved(0)
  , m_zmin(-102.0)
  , m_zmax(102.0)
  , m_xymax(85.0)
  , m_outfile(nullptr)
  , m_clusterTracks(nullptr)
{
}

TpcPolyClusterDisplay::~TpcPolyClusterDisplay()
{
  if (m_outfile)
  {
    delete m_outfile;
    m_outfile = nullptr;
  }
}

int TpcPolyClusterDisplay::Init(PHCompositeNode*)
{
  m_outfile = new TFile(m_outfilename.c_str(), "RECREATE");
  if (!m_outfile || m_outfile->IsZombie())
  {
    std::cerr << Name() << "::Init - cannot open output file " << m_outfilename << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  m_outfile->mkdir("events");
  return Fun4AllReturnCodes::EVENT_OK;
}

bool TpcPolyClusterDisplay::get_nodes(PHCompositeNode* topNode)
{
  m_clusterTracks = findNode::getClass<TpcPolyClusterTrackContainer>(topNode, m_clusterNodeName.c_str());
  if (!m_clusterTracks)
  {
    std::cerr << Name() << " - missing " << m_clusterNodeName << std::endl;
    return false;
  }
  return true;
}

int TpcPolyClusterDisplay::process_event(PHCompositeNode* topNode)
{
  ++m_evt;
  if (!get_nodes(topNode)) return Fun4AllReturnCodes::EVENT_OK;
  if (!m_outfile || m_eventsSaved >= m_maxEventDisplays) return Fun4AllReturnCodes::EVENT_OK;

  TDirectory* eventsTop = m_outfile->GetDirectory("events");
  if (!eventsTop) eventsTop = m_outfile->mkdir("events");
  eventsTop->cd();

  TDirectory* eventDir = eventsTop->mkdir(Form("event_%06u", m_evt));
  if (!eventDir) return Fun4AllReturnCodes::EVENT_OK;
  eventDir->cd();

  TH3D* h3 = new TH3D(Form("h3_evt%06u_tpcpolycluster_centroids_bothsides", m_evt),
                      Form("event %u both sides TpcPolyCluster centroids;z [cm];x [cm];y [cm]", m_evt),
                      204, m_zmin, m_zmax,
                      170, -m_xymax, m_xymax,
                      170, -m_xymax, m_xymax);
  h3->SetStats(0);
  h3->SetDirectory(nullptr);

  std::vector<TPolyMarker3D*> markers;
  const unsigned int ncluster_tracks = m_clusterTracks->size();
  unsigned int nclusters = 0;
  unsigned int nclustersPlotted = 0;
  for (unsigned int itrack = 0; itrack < ncluster_tracks; ++itrack)
  {
    const TpcPolyClusterTrack* trk = m_clusterTracks->get_track(itrack);
    if (!trk || !trk->isValid()) continue;

    for (unsigned int icluster = 0; icluster < trk->size_clusters(); ++icluster)
    {
      ++nclusters;
      const double x = trk->get_cluster_x(icluster);
      const double y = trk->get_cluster_y(icluster);
      const double z = trk->get_cluster_z(icluster);
      if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) continue;
      if (z < m_zmin || z > m_zmax) continue;

      h3->Fill(z, x, y);
      markers.push_back(make_cluster_marker(z, x, y, cluster_color(itrack)));
      ++nclustersPlotted;
    }
  }

  TCanvas* c3 = new TCanvas(Form("c3_evt%06u_tpcpolycluster_z_x_y_centroids_bothsides", m_evt),
                            Form("event %u both sides TpcPolyCluster centroids", m_evt),
                            1200, 900);
  h3->Draw("BOX2Z");
  for (TPolyMarker3D* marker : markers)
  {
    if (marker) marker->Draw("same");
  }
  c3->Modified();
  c3->Update();
  c3->Write();

  std::cout << Name() << " - saved event " << m_evt
            << " cluster_tracks=" << ncluster_tracks
            << " layer_clusters=" << nclusters
            << " plotted=" << nclustersPlotted << std::endl;

  ++m_eventsSaved;
  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcPolyClusterDisplay::End(PHCompositeNode*)
{
  if (m_outfile)
  {
    m_outfile->Close();
    delete m_outfile;
    m_outfile = nullptr;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
