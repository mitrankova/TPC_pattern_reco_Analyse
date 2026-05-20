#ifndef INMODULETRACKDISPLAY_H
#define INMODULETRACKDISPLAY_H

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

#include <string>

class PHCompositeNode;
class TFile;
class TH1D;
class TH2D;
class TH3D;
class InModuleTrackContainer;
class TrkrHitSetContainer;
class PHG4TpcGeom;
class PHG4TpcGeomContainer;
class PHG4CylinderGeomContainer;

class InModuleTrackDisplay : public SubsysReco
{
 public:
  InModuleTrackDisplay(const std::string& name = "InModuleTrackDisplay",
                       const std::string& outfilename = "inmodule_display.root",
                       const std::string& trackNodeName = "INMODULETRACKS",
                       unsigned int maxEventDisplays = 5);
  ~InModuleTrackDisplay() override {}

  int Init(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int End(PHCompositeNode* topNode) override;

 private:
  struct HitPoint
  {
    HitPoint();
    bool ok;
    unsigned int layer;
    unsigned int pad;
    unsigned int tbin;
    unsigned short adc;
    double x;
    double y;
    double r;
    double phi;
  };

  bool get_nodes(PHCompositeNode* topNode);
  HitPoint make_hit_point(const TrkrDefs::hitsetkey hsk, const TrkrDefs::hitkey hk) const;

  std::string m_outfilename;
  std::string m_trackNodeName;
  unsigned int m_maxEventDisplays;
  unsigned int m_evt;
  unsigned int m_eventsSaved;

  TFile* m_outfile;
  InModuleTrackContainer* m_tracks;
  TrkrHitSetContainer* m_hits;
  PHG4TpcGeomContainer* m_tpcGeom;

  TH1D* h_ntracks;
  TH1D* h_nhits_per_track;
  TH1D* h_nrawhits_field;
  TH1D* h_nblobs;
  TH1D* h_chi2pad_ndf;
  TH1D* h_chi2tbin_ndf;
  TH2D* h_layer_first_last;
  TH2D* h_pad_slope_vs_tbin_slope;
  TH2D* h_sector_side;
  TH2D* h_xy_all;
  TH2D* h_pad_tbin_all;
  TH2D* h_layer_tbin_all;
  TH2D* h_layer_pad_all;
  TH3D* h_xyz_all;
};

#endif
