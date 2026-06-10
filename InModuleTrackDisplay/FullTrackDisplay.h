#ifndef FULLTRACKDISPLAY_H
#define FULLTRACKDISPLAY_H

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

#include <string>

class PHCompositeNode;
class TFile;
class FullTrackContainer;
class TrkrHitSetContainer;
class IdealPadMap;

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

  // Display-only fit mode. Use Fitter::FIT_SAGITTA for field-on curvature
  // or Fitter::FIT_LINEAR for straight-line fits. Default is sagitta.
  void setFitMode(int mode) { m_fitMode = mode; }
  void setFitWithSagitta(bool v) { m_fitMode = v ? 1 : 0; }

  void setFitWeightPower(double v) { m_fitWeightPower = v; }
  void setFitWeightFloorFrac(double v) { m_fitWeightFloorFrac = v; }
  void setFitWeights(double power, double floor_frac)
  {
    m_fitWeightPower = power;
    m_fitWeightFloorFrac = floor_frac;
  }

 private:
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
  IdealPadMap* m_idealPadMap;
  int m_fitMode;
  double m_fitWeightPower;
  double m_fitWeightFloorFrac;
};

#endif
