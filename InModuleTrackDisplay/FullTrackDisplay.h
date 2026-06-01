#ifndef FULLTRACKDISPLAY_H
#define FULLTRACKDISPLAY_H

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

#include <string>

class PHCompositeNode;
class TFile;
class FullTrackContainer;
class TrkrHitSetContainer;
class TpcPadMap;

class FullTrackDisplay : public SubsysReco
{
 public:
  FullTrackDisplay(const std::string& name = "FullTrackDisplay",
                   const std::string& outfilename = "fulltrack_display.root",
                   const std::string& trackNodeName = "FULLTRACKS",
                   unsigned int maxEventDisplays = 5);
  ~FullTrackDisplay() override;

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
    unsigned int region;
    unsigned int sector;
    int side;

    unsigned int pad;
    unsigned int tbin;
    unsigned short adc;

    double global_phi;
    double radius;
  };

  bool get_nodes(PHCompositeNode* topNode);
  HitPoint make_hit_point(TrkrDefs::hitsetkey hsk,
                          TrkrDefs::hitkey hk) const;

  std::string m_outfilename;
  std::string m_trackNodeName;
  unsigned int m_maxEventDisplays;
  unsigned int m_evt;
  unsigned int m_eventsSaved;

  TFile* m_outfile;
  FullTrackContainer* m_tracks;
  TrkrHitSetContainer* m_hits;
  TpcPadMap* m_padMap;
};

#endif
