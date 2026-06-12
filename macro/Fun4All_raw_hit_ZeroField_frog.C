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

#include </sphenix/user/mitrankova/F4A/TPC_pattern_reco/install/include/inmoduletracks/TpcPadMapBuilder.h>
#include </sphenix/user/mitrankova/F4A/TPC_pattern_reco/install/include/inmoduletracks/InModuleTracks.h>
#include </sphenix/user/mitrankova/F4A/TPC_pattern_reco/install/include/inmoduletracks/FullTrackConnector.h>
#include </sphenix/user/mitrankova/F4A/TPC_pattern_reco/install/include/inmoduletracks/TpcPolyTrackReco.h>
#include </sphenix/user/mitrankova/F4A/TPC_pattern_reco/install/include/inmoduletracks/FullTrackVertexer.h>
#include </sphenix/user/mitrankova/F4A/InModuleTrackDisplay/install/include/inmoduletrackdisplay/InModuleTrackDisplay.h>
#include </sphenix/user/mitrankova/F4A/InModuleTrackDisplay/install/include/inmoduletrackdisplay/FullTrackDisplay.h>
#include </sphenix/user/mitrankova/F4A/InModuleTrackDisplay/install/include/inmoduletrackdisplay/TpcPolyTrackDisplay.h>



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

void Fun4All_raw_hit_ZeroField_frog(
    const int nEvents = 10,
    const int runnumber = 81000,
    const int segment = 0,
    const std::string outdir = ".",
    const int nSkip = 0,
    const std::string collision = "run3pp",
    const std::string production = "ana537_nocdbtag_v001",
    const std::string& outfilename = "HITS_clusters_seeds")
{
  const bool convertSeeds = true;
  auto *se = Fun4AllServer::instance();
  se->Verbosity(1);
  auto *rc = recoConsts::instance();
  se->registerSubsystem(new SkipFirstN(nSkip));

  const std::string dsttype = "STREAMING_EVENT";


  std::ostringstream runstr;
  runstr << std::setw(8) << std::setfill('0') << runnumber;

  std::ostringstream segstr;
  segstr << std::setw(5) << std::setfill('0') << segment;

  const int runbase = (runnumber / 100) * 100;
  const int runnext = runbase + 100;

  std::ostringstream rundir;
  rundir << "run_" << std::setw(8) << std::setfill('0') << runbase << "_" << std::setw(8) << std::setfill('0') << runnext;

  std::vector<std::string> streams;

  // TPC EBDC streams
  for (int ebdc = 0; ebdc < 24; ++ebdc)
  {
    for (int endpoint = 0; endpoint < 2; ++endpoint)
    {
      std::ostringstream s;
      s << "ebdc" << std::setw(2) << std::setfill('0') << ebdc << "_" << endpoint;
      streams.push_back(s.str());
    }
  }
if(collision!="run3line_laser"&&collision!="run3cosmics")
{
  // TPOT
  streams.push_back("ebdc39");

  // INTT streams
  for (int server = 0; server < 8; ++server)
  {
    std::ostringstream s;
    s << "intt" << server;
    streams.push_back(s.str());
  }

  // MVTX streams
  for (int felix = 0; felix < 6; ++felix)
  {
    std::ostringstream s;
    s << "mvtx" << felix;
    streams.push_back(s.str());
  }
}

  int i = 0;
    std::stringstream nice_runnumber;
  nice_runnumber << std::setw(8) << std::setfill('0') << std::to_string(runnumber);

  int rounded_up = 100 * (std::ceil((float) runnumber / 100));
  std::stringstream nice_rounded_up;
  nice_rounded_up << std::setw(8) << std::setfill('0') << std::to_string(rounded_up);

  int rounded_down = 100 * (std::floor((float) runnumber / 100));
  std::stringstream nice_rounded_down;
  nice_rounded_down << std::setw(8) << std::setfill('0') << std::to_string(rounded_down);

  for (unsigned int is = 0; is < streams.size(); ++is)
  {
    const std::string stream = streams[is];

    std::string filename = "DST_" + dsttype + "_" + stream + "_" + collision + "_" + production + "-" +  runstr.str() + "-" + segstr.str() + ".root";
   std::string filepath = "/sphenix/lustre01/sphnxpro/production/" + collision + "/physics/" + production + "/DST_" + dsttype + "_" + stream + "/run_" + nice_rounded_down.str()  + "_" + nice_rounded_up.str()  + "/" + filename;

    std::ifstream testfile(filepath.c_str());
    if (!testfile.good())
    {
      std::cout << "Skipping missing DST: " << filepath << std::endl;
      continue;
    }

    std::cout << "Adding DST: " << filepath << std::endl;
    if (i == 0)
    {
      rc->set_IntFlag("RUNNUMBER", runnumber);
      rc->set_uint64Flag("TIMESTAMP", runnumber);
    }

    std::string inputname = "InputManager" + std::to_string(i);
    auto *hitsin = new Fun4AllDstInputManager(inputname);
    hitsin->fileopen(filepath);
    se->registerInputManager(hitsin);

    ++i;
  }

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

  TRACKING::pp_mode = true;



  FlagHandler *flag = new FlagHandler();
  se->registerSubsystem(flag);

  std::string geofile = CDBInterface::instance()->getUrl("Tracking_Geometry");

  Fun4AllRunNodeInputManager *ingeo = new Fun4AllRunNodeInputManager("GeoIn");
  ingeo->AddFile(geofile);
  se->registerInputManager(ingeo);

  TpcReadoutInit(runnumber);
  G4MAGNET::magfield = "0.01"; 
  G4MAGNET::magfield_tracking = G4MAGNET::magfield; 
  G4MAGNET::magfield_rescale = 1;

  G4TPC::REJECT_LASER_EVENTS = true;
  G4TPC::ENABLE_MODULE_EDGE_CORRECTIONS = false;
  // Flag for running the tpc hit unpacker with zero suppression on
  TRACKING::tpc_zero_supp = true;

  // MVTX
  Enable::MVTX_APPLYMISALIGNMENT = true;
  ACTSGEOM::mvtx_applymisalignment = Enable::MVTX_APPLYMISALIGNMENT;


  TrackingInit();

  for (int felix = 0; felix < 6; felix++)
  {
    Mvtx_HitUnpacking(std::to_string(felix));
  }
  for (int server = 0; server < 8; server++)
  {
    Intt_HitUnpacking(std::to_string(server));
  }
  std::ostringstream ebdcname;
  for (int ebdc = 0; ebdc < 24; ebdc++)
  {
    for (int endpoint = 0; endpoint < 2; endpoint++)
      {
        ebdcname.str("");
        if (ebdc < 10)
	  {
	    ebdcname << "0";
	  }
        ebdcname << ebdc << "_" << endpoint;
        Tpc_HitUnpacking(ebdcname.str());
      }
  }


  Micromegas_HitUnpacking();

  Tpc_LaserEventIdentifying();

  Reject_Laser_Events();



  se->registerSubsystem( new InModuleTracks());
  se->registerSubsystem( new FullTrackConnector());
  se->registerSubsystem( new TpcPolyTrackReco());

  //se->registerSubsystem(new TpcPadMapBuilder("TpcPadMapBuilder", "TPC_PADMAP"));

  //  se->registerSubsystem( new InModuleTrackDisplay( "InModuleTrackDisplay", "inmodule_display_"+outfilename+"_" + to_string(runnumber) + ".root" ));
  //  InModuleTrackDisplay* imt = new InModuleTrackDisplay("InModuleTrackDisplay", "inmodule_display_"+outfilename+"_" + to_string(runnumber) + ".root" );
  //  imt->Verbosity(10);   // or 2, 3, etc.
  //  se->registerSubsystem(imt);


  // se->registerSubsystem( new FullTrackDisplay( "FullTrackDisplay", "full_track_display_"+outfilename+"_" + to_string(runnumber) + ".root" ));

  // auto display = new FullTrackDisplay( "FullTrackDisplay", "34_full_track_display_"+outfilename+"_" + to_string(runnumber) + ".root" );
   //display->setPlotNModulesRange(1, 2);
   // display->setPlotNModulesRange(3, 4);
   // display->setPlotDirectionRange(-0.02, 0.02);  // radius-phi angle
   //display->setPlotThetaRange(-1.05, 1.05);        // radius-timebin angle
   //display->setPlotCurvatureRange(-0.005, 0.005);
   // se->registerSubsystem(display);
    


//se->registerSubsystem( new TpcPolyTrackDisplay( "TpcPolyTrackDisplay", "tpc_poly_track_display_"+outfilename+"_" + to_string(runnumber) + ".root",  "TPCPOLYTRACKS", 2));
    auto display = new TpcPolyTrackDisplay( "TpcPolyTrackDisplay", "tpc_poly_track_display_"+outfilename+"_" + to_string(runnumber) + ".root",  "TPCPOLYTRACKS", 2);
    display->setHighPtMin(0.5);
    display->setHighPtMaxAbsD0(1.0);
    se->registerSubsystem(display);
/*
auto poly = new TpcPolyTrackReco();
poly->setInputNodeName("FULLTRACKS");
poly->setOutputNodeName("TPCPOLYTRACKS");
poly->setT0(360.0);
poly->setTpcAdcClock(50.037280);
poly->setStartZ(-105.5, 105.5); // side 0, side 1
se->registerSubsystem(poly);
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
  Fun4AllOutputManager *out = new Fun4AllDstOutputManager("out", Form("%s/output_DST/%s_%d_%d.root",outdir.c_str(), outfilename.c_str(), runnumber, segment));

  out->AddNode("Sync");
  out->AddNode("EventHeader");
  out->AddNode("TRKR_HITSET");
  out->AddNode("INMODULETRACKS");
  out->AddRunNode("TPCGEOMCONTAINER");

  //se->registerOutputManager(out);
  
  se->run(nEvents+nSkip);
  se->Print("NODETREE");
  se->End();
  se->PrintTimer();


  CDBInterface::instance()->Print();
  delete se;
  std::cout << "Finished" << std::endl;
  gSystem->Exit(0);
}