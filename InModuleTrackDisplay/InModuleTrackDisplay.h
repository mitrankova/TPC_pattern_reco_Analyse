#ifndef INMODULETRACKDISPLAY_H
#define INMODULETRACKDISPLAY_H

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

#include <string>

class PHCompositeNode;
class TFile;
class InModuleTrackContainer;
class TrkrHitSetContainer;
class PHG4TpcGeomContainer;
class TpcPadMap;

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
    TrkrDefs::hitsetkey hitsetkey;
    TrkrDefs::hitkey hitkey;
    unsigned int layer;
    unsigned int pad;
    unsigned int tbin;
    unsigned short adc;
    double x;
    double y;
    double r;
    double phi;       // global phi
    double local_phi; // sector-local phi used by in-module fit formulas
  };

  bool get_nodes(PHCompositeNode* topNode);
  HitPoint make_hit_point(const TrkrDefs::hitsetkey hsk,
                          const TrkrDefs::hitkey hk) const;

  std::string m_outfilename;
  std::string m_trackNodeName;
  unsigned int m_maxEventDisplays;
  unsigned int m_evt;
  unsigned int m_eventsSaved;

  TFile* m_outfile;
  InModuleTrackContainer* m_tracks;
  TrkrHitSetContainer* m_hits;
  PHG4TpcGeomContainer* m_tpcGeom;
  TpcPadMap* m_padMap;
};

#endif
