/*
 * This macro shows a minimum working example of running the tracking
 * hit unpackers with some basic seeding algorithms to try to put together
 * tracks. There are some analysis modules run at the end which package
 * hits, clusters, and clusters on tracks into trees for analysis.
 */

// leave the GlobalVariables.C at the beginning, an empty line afterwards
// protects its position against reshuffling by clang-format
#include <GlobalVariables.C>

#include <G4_ActsGeom.C>
#include <G4_Global.C>
#include <G4_Magnet.C>
#include <G4_Mbd.C>
#include <QA.C>
#include <Trkr_Clustering.C>
#include <Trkr_LaserClustering.C>
#include <Trkr_Reco.C>
#include <Trkr_RecoInit.C>
#include <Trkr_TpcReadoutInit.C>

#include <cdbobjects/CDBTTree.h>

#include <ffamodules/CDBInterface.h>
#include <ffamodules/FlagHandler.h>

#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllUtils.h>

#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/recoConsts.h>

//#include </sphenix/user/mitrankova/F4A/TPC_pattern_reco/install/include/inmoduletracks/TpcPadMapBuilder.h>
#include </sphenix/user/mitrankova/F4A/TPC_pattern_reco/install/include/inmoduletracks/InModuleTracks.h>
#include </sphenix/user/mitrankova/F4A/TPC_pattern_reco/install/include/inmoduletracks/FullTrackConnector.h>
//#include </sphenix/user/mitrankova/F4A/TPC_pattern_reco/install/include/inmoduletracks/TpcPolyTrackReco.h>
//#include </sphenix/user/mitrankova/F4A/TPC_pattern_reco/install/include/inmoduletracks/FullTrackVertexer.h>
#include </sphenix/user/mitrankova/F4A/TPC_pattern_reco/install/include/inmoduletracks/TpcPolyClusterizer.h>
#include </sphenix/user/mitrankova/F4A/TPC_pattern_reco/install/include/inmoduletracks/FinalTrackVertexer.h>
/*
#include </sphenix/user/mitrankova/F4A/InModuleTrackDisplay/install/include/inmoduletrackdisplay/InModuleTrackDisplay.h>
#include </sphenix/user/mitrankova/F4A/InModuleTrackDisplay/install/include/inmoduletrackdisplay/FullTrackDisplay.h>
#include </sphenix/user/mitrankova/F4A/InModuleTrackDisplay/install/include/inmoduletrackdisplay/TpcPolyTrackDisplay.h>
*/
#include </sphenix/user/mitrankova/F4A/InModuleTrackDisplay/install/include/inmoduletrackdisplay/TpcPolyClusterDisplay.h>

#include </sphenix/user/mitrankova/F4A/InModuleTrackDisplay/install/include/inmoduletrackdisplay/TpcPolyClusterResiduals.h>


R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libphool.so)
R__LOAD_LIBRARY(libcdbobjects.so)
R__LOAD_LIBRARY(libmvtx.so)
R__LOAD_LIBRARY(libintt.so)
R__LOAD_LIBRARY(libtpc.so)
R__LOAD_LIBRARY(libmicromegas.so)
R__LOAD_LIBRARY(/sphenix/user/mitrankova/F4A/TPC_pattern_reco/install/lib/libInModuleTracks.so)
R__LOAD_LIBRARY(/sphenix/user/mitrankova/F4A/InModuleTrackDisplay/install/lib/libInModuleTrackDisplay.so)


class SkipFirstN : public SubsysReco {
 public:
  explicit SkipFirstN(int n) : SubsysReco("SkipFirstN"), target_(n) {}
  int process_event(PHCompositeNode*) override {
    if (count_ < target_) { ++count_; return Fun4AllReturnCodes::ABORTEVENT; }
    return Fun4AllReturnCodes::EVENT_OK;
  }
 private:
  int target_ = 0;
  int count_  = 0;
};
//xxxxxxxxxxxxxxxxxxxxxxxxxx  .x Fun4All_raw_hit_ZeroField_frog.C(10, 81000, 0, ".", 0, "run3pp", "ana537_nocdbtag_v001","HITS_ppZeroField")
//xxxxxxxxxxxxxxxxxxxxxxxxxx  .x Fun4All_raw_hit_ZeroField_frog.C(10, 81000, 0, ".", 0, "run3pp", "ana537_nocdbtag_v001","HITS_ppZeroField")


//xxxxxxxxxxxxxxxxxxxxxxxxxx  .x Fun4All_raw_hit_ZeroField_frog.C(3, 79513, 0, ".", 0, "run3pp", "ana537_nocdbtag_v001","HITS_ppFieldOn")
//.x Fun4All_raw_hit_ZeroField_frog.C(2, 79513, 0, ".", 0, "run3pp", "ana532_nocdbtag_v001","HITS_ppFieldOn")
//.x Fun4All_raw_hit_ZeroField_frog.C(2, 79516, 0, ".", 0, "run3pp", "ana532_nocdbtag_v001","HITS_ppFieldOn")

//==========AuAu Zero Field
// 111x111   75555-75557 - 1mrad; 75557 - 0mrad
//.x Fun4All_raw_hit_ZeroField_frog.C(10, 75555, 0, ".", 0, "run3auau", "ana514_nocdbtag_v001","HITS_AuAu_ZeroField_1mrad")

// 6x6       75570-75573 - 1mrad; 75574 - 0mrad
//.x Fun4All_raw_hit_ZeroField_frog.C(10, 75570, 0, ".", 0, "run3auau", "ana514_nocdbtag_v001","HITS_AuAu_ZeroField_1mrad")
//.x Fun4All_raw_hit_ZeroField_frog.C(10, 75574, 0, ".", 0, "run3auau", "ana514_nocdbtag_v001","HITS_AuAu_ZeroField_0mrad")

//.x Fun4All_raw_hit_ZeroField_frog.C(10, 75396, 0, ".", 0, "run3auau", "ana514_nocdbtag_v001","HITS_AuAu")
//.x Fun4All_raw_hit_ZeroField_frog.C(2, 75405, 0, ".", 0, "run3auau", "ana514_nocdbtag_v001","HITS_AuAu")

//.x Fun4All_raw_hit_ZeroField_frog.C(2, 82626, 0, ".", 0, "run3oo", "ana537_nocdbtag_v001","HITS_OO")
//.x Fun4All_raw_hit_ZeroField_frog.C(10, 79767, 0, ".", 0, "run3cosmics", "ana523_nocdbtag_v001","HITS_cosmics")
//.x Fun4All_raw_hit_ZeroField_frog.C(10, 83174, 0, ".", 0, "run3line_laser", "ana540_nocdbtag_v001","HITS_line_laser")


//.x Fun4All_TPC_Tracking.C(2,75405, 0, "/sphenix/tg/tg01/hf/mitrankova/PatternReco/75405/helixfit", 0, "run3auau", "ana514_nocdbtag_v001","HITS_clusters_seeds","HITS_HELIX_AuAu")
void Fun4All_TPC_Tracking(
    const int nEvents = 10,
    const int runnumber = 79513,
    const int segment = 0,
    const std::string outdir = ".",
    const std::string indir = ".",
    const int nSkip = 0,
    const std::string collision = "run3pp",
    const std::string production = "ana532_nocdbtag_v001",
    const std::string& infilename = "HITS_clusters_seeds",
    const std::string& outfilename = "HITS_ppFieldOn")
{
  const bool convertSeeds = true;
  auto *se = Fun4AllServer::instance();
  se->Verbosity(1);
  auto *rc = recoConsts::instance();
  se->registerSubsystem(new SkipFirstN(nSkip));

  const std::string dsttype = "STREAMING_EVENT";

 // std::string filename = indir+"/output_DST/"+infilename+"_"+to_string(runnumber)+"_"+to_string(segment)+".root";
   std::string filename = indir+"/"+infilename+"_"+to_string(runnumber)+"_"+to_string(segment)+".root";

  auto *hitsinclus = new Fun4AllDstInputManager("ClusterInputManager");
  hitsinclus->fileopen(filename);
  se->registerInputManager(hitsinclus);

  rc->set_IntFlag("RUNNUMBER", runnumber);
  rc->set_IntFlag("RUNSEGMENT", segment);

  Enable::QA = false;
  Enable::CDB = true;
  rc->set_StringFlag("CDB_GLOBALTAG", "newcdbtag");
  rc->set_uint64Flag("TIMESTAMP", runnumber);

 

  G4TRACKING::convert_seeds_to_svtxtracks = convertSeeds;


  std::cout << " run: " << runnumber
            << " samples: " << TRACKING::reco_tpc_maxtime_sample
            << " pre: " << TRACKING::reco_tpc_time_presample
            << " vdrift: " << G4TPC::tpc_drift_velocity_reco
            << std::endl;

  TRACKING::streaming_mode = true;


  FlagHandler *flag = new FlagHandler();
  se->registerSubsystem(flag);

  std::string geofile = CDBInterface::instance()->getUrl("Tracking_Geometry");

  Fun4AllRunNodeInputManager *ingeo = new Fun4AllRunNodeInputManager("GeoIn");
  ingeo->AddFile(geofile);
  se->registerInputManager(ingeo);



  TpcReadoutInit(runnumber);
  G4TPC::REJECT_LASER_EVENTS = true;
  // Flag for running the tpc hit unpacker with zero suppression on
  TRACKING::tpc_zero_supp = true;




 // se->registerSubsystem( new InModuleTracks());
  //se->registerSubsystem( new FullTrackConnector());
  //se->registerSubsystem( new FullTrackVertexer());
  //se->registerSubsystem( new TpcPolyTrackReco());
   //se->registerSubsystem( new TpcPolyClusterizer());
  // se->registerSubsystem( new FinalTrackVertexer());
   //se->registerSubsystem( new TpcPolyClusterDisplay("TpcPolyClusterDisplay", "tpc_poly_cluster_display_"+outfilename+"_" + to_string(runnumber) + ".root" ));
      auto resid = new TpcPolyClusterResiduals("TpcPolyClusterResiduals",
                                         outdir+"/outout_resid/tpc_poly_cluster_residuals"+outfilename+"_" + to_string(runnumber) + to_string(segment) + ".root" );
resid->setMinPt(0.2);
resid->setMinTpcClusters(20);
se->registerSubsystem(resid);

/*
    FullTrackConnector* fullconn = new FullTrackConnector();
    fullconn->Verbosity(0);
    se->registerSubsystem(fullconn);
    */
    
    /*FullTrackVertexer* vtx = new FullTrackVertexer();
    vtx->Verbosity(0);
    se->registerSubsystem(vtx);
    */
  //se->registerSubsystem(new TpcPadMapBuilder("TpcPadMapBuilder", "TPC_PADMAP"));

  //  se->registerSubsystem( new InModuleTrackDisplay( "InModuleTrackDisplay", "inmodule_display_"+outfilename+"_" + to_string(runnumber) + ".root" ));
  //  InModuleTrackDisplay* imt = new InModuleTrackDisplay("InModuleTrackDisplay", "inmodule_display_"+outfilename+"_" + to_string(runnumber) + ".root" );
  //  imt->Verbosity(10);   // or 2, 3, etc.
  //  se->registerSubsystem(imt);

   //se->registerSubsystem( new FullTrackDisplay( "FullTrackDisplay", "full_track_display_"+outfilename+"_" + to_string(runnumber) + ".root" ));

  // auto display = new FullTrackDisplay( "FullTrackDisplay", "34_full_track_display_"+outfilename+"_" + to_string(runnumber) + ".root" );
   //display->setPlotNModulesRange(1, 2);
   // display->setPlotNModulesRange(3, 4);
   // display->setPlotDirectionRange(-0.02, 0.02);  // radius-phi angle
   //display->setPlotThetaRange(-1.05, 1.05);        // radius-timebin angle
   //display->setPlotCurvatureRange(-0.005, 0.005);
   // se->registerSubsystem(display);
    


//se->registerSubsystem( new TpcPolyTrackDisplay( "TpcPolyTrackDisplay", "tpc_poly_track_display_"+outfilename+"_" + to_string(runnumber) + ".root",  "TPCPOLYTRACKS", 2));
  
/*
auto poly = new TpcPolyTrackReco();
poly->Verbosity(10);
se->registerSubsystem(poly);
*/
    /* 
    auto display = new TpcPolyTrackDisplay( "TpcPolyTrackDisplay", "tpc_poly_track_display_"+outfilename+"_" + to_string(runnumber) + ".root",  "TPCPOLYTRACKS", 2);
    display->setHighPtMin(0.5);
    display->setHighPtMaxAbsD0(1.0);
    se->registerSubsystem(display);
    */

/*
    FullTrackConnector* fullconn = new FullTrackConnector("FullTrackConnector", "FullTracks.root");
    fullconn->Verbosity(10);
    fullconn->setConnectMaxLayerGap(17);
    //fullconn->setConnectWindow(0.08, 20.0);
    //fullconn->setConnectSlopeWindow(5.0, 0.03);
    se->registerSubsystem(fullconn);


    FullTrackVertexer* vtx = new FullTrackVertexer();
    vtx->Verbosity(10);   // or 2, 3, etc.
    se->registerSubsystem(vtx);
*/
   
 //std::cout<< "Output DST "<<Form("%s/output_DST/%s_%d_%d.root",outdir.c_str(), outfilename.c_str(), runnumber, segment)  << std::endl;

  
  se->run(nEvents+nSkip);
  se->Print("NODETREE");
  se->End();
  se->PrintTimer();


  CDBInterface::instance()->Print();
  delete se;
  std::cout << "Finished" << std::endl;
  gSystem->Exit(0);
}