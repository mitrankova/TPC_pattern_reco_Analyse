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

#include <tpccalib/PHTpcResiduals.h>

#include <mvtxrawhitqa/MvtxRawHitQA.h>

#include <inttrawhitqa/InttRawHitQA.h>

#include <trackingqa/InttClusterQA.h>
#include <trackingqa/MicromegasClusterQA.h>
#include <trackingqa/MvtxClusterQA.h>
#include <trackingqa/SiliconSeedsQA.h>
#include <trackingqa/TpcClusterQA.h>
#include <trackingqa/TpcSeedsQA.h>
#include <trackingqa/TpcSiliconQA.h>

#include <tpcqa/TpcRawHitQA.h>

#include <trackingdiagnostics/KshortReconstruction.h>
#include <trackingdiagnostics/TrackResiduals.h>
#include <trackingdiagnostics/TrkrNtuplizer.h>



#include <ffamodules/CDBInterface.h>
#include <ffamodules/FlagHandler.h>

#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllUtils.h>

#include <phool/recoConsts.h>


R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libphool.so)
R__LOAD_LIBRARY(libcdbobjects.so)
R__LOAD_LIBRARY(libmvtx.so)
R__LOAD_LIBRARY(libintt.so)
R__LOAD_LIBRARY(libtpc.so)
R__LOAD_LIBRARY(libmicromegas.so)
R__LOAD_LIBRARY(libTrackingDiagnostics.so)
R__LOAD_LIBRARY(libtrackingqa.so)
R__LOAD_LIBRARY(libtpcqa.so)


#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllReturnCodes.h>

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

void Fun4All_raw_hit_ZeroField_std(
    const int nEvents = 10,
    const std::string& filelist = "/sphenix/user/mitrankova/pp2025/ZeroField/macro/filelists_81000/filelist_00081000_00000.list",
    const std::string outdir = "./",
    const int nSkip = 0,
    const std::string& outfilename = "STD_clusters_seeds",
    const bool convertSeeds = true)
{
  auto *se = Fun4AllServer::instance();
  se->Verbosity(1);
  auto *rc = recoConsts::instance();
  se->registerSubsystem(new SkipFirstN(nSkip));
  // input manager for QM production raw hit DST file
  std::ifstream ifs(filelist);
  std::string filepath;

  int i = 0;
  int runnumber = std::numeric_limits<int>::quiet_NaN();
  int segment = std::numeric_limits<int>::quiet_NaN();
  

  while (std::getline(ifs, filepath))
  {
    std::cout << "Adding DST with filepath: " << filepath << std::endl;
    if (i == 0)
    {
      std::pair<int, int>
          runseg = Fun4AllUtils::GetRunSegment(filepath);
      runnumber = runseg.first;
      segment = runseg.second;
      rc->set_IntFlag("RUNNUMBER", runnumber);
      rc->set_uint64Flag("TIMESTAMP", runnumber);
    }
  
    std::string inputname = "InputManager" + std::to_string(i);
    auto *hitsin = new Fun4AllDstInputManager(inputname);
    hitsin->fileopen(filepath);
    se->registerInputManager(hitsin);
    i++;
  }

  rc->set_IntFlag("RUNNUMBER", runnumber);
  rc->set_IntFlag("RUNSEGMENT", segment);

  Enable::QA = false;
  Enable::CDB = true;
  rc->set_StringFlag("CDB_GLOBALTAG", "newcdbtag");
  rc->set_uint64Flag("TIMESTAMP", runnumber);

  std::stringstream nice_runnumber;
  nice_runnumber << std::setw(8) << std::setfill('0') << std::to_string(runnumber);

 

  G4TRACKING::convert_seeds_to_svtxtracks = convertSeeds;
  std::cout << "Converting to seeds : " << G4TRACKING::convert_seeds_to_svtxtracks << std::endl;

  std::cout << " run: " << runnumber
            << " samples: " << TRACKING::reco_tpc_maxtime_sample
            << " pre: " << TRACKING::reco_tpc_time_presample
            << " vdrift: " << G4TPC::tpc_drift_velocity_reco
            << std::endl;

  TRACKING::pp_mode = true;

  // distortion calibration mode
  /*
   * set to true to enable residuals in the TPC with
   * TPC clusters not participating to the ACTS track fit
   */
  /*TString outfile = Form(
      "%s%s_%s-%s_%dskip_%devts",
      outdir.Data(),
      outfilename.Data(),
      runnumber.Data(),
      segment.Data(),
      nSkip,
      nEvents
  );*/

  TString outfile = outdir + outfilename + "_" + to_string(nEvents) + "evts_"  + to_string(nSkip) + "nSkip_" + runnumber + "-" + segment;
  std::string theOutfile = outfile.Data();


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

  // to turn on the default static corrections, enable the two lines below
  //G4TPC::ENABLE_STATIC_CORRECTIONS = true;
  //G4TPC::USE_PHI_AS_RAD_STATIC_CORRECTIONS = false;

  // to turn on the average corrections derived from simulation, enable the three lines below
  // note: these are designed to be used only if static corrections are also applied
   //G4TPC::ENABLE_AVERAGE_CORRECTIONS = true;
   //G4TPC::USE_PHI_AS_RAD_AVERAGE_CORRECTIONS = false;
   //G4TPC::average_correction_filename = CDBInterface::instance()->getUrl("TPC_LAMINATION_FIT_CORRECTION");

  

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

  Mvtx_Clustering();

  Intt_Clustering();

  Tpc_LaserEventIdentifying();

  auto *tpcclusterizer = new TpcClusterizer;
  tpcclusterizer->Verbosity(0);
  tpcclusterizer->set_do_hit_association(G4TPC::DO_HIT_ASSOCIATION);
  tpcclusterizer->set_rawdata_reco();
  tpcclusterizer->set_reject_event(G4TPC::REJECT_LASER_EVENTS);
  se->registerSubsystem(tpcclusterizer);

  Micromegas_Clustering();

  Reject_Laser_Events();
  /*
   * Begin Track Seeding
   */
  
  //Tracking_Reco_TrackSeed_run2pp();
  Tracking_Reco_TrackSeed_ZeroField();
 // Tracking_Reco_TrackMatching_run2pp();

  auto silicon_match = new PHSiliconTpcTrackMatching;
  silicon_match->Verbosity(0);
  //silicon_match->set_use_legacy_windowing(false);
  // set search windows matching Silicon to TPC seeds
  // Selected for tracks with ntpc>34,|z_Si-z_TPC|<30,crossing==0
  // see https://indico.bnl.gov/event/26202/attachments/59703/102575/2025_01_31_ZeroField.pdf
  // Will probably be narrowed from improved alignment soon.
  silicon_match->window_dx.set_QoverpT_maxabs({2.6,0,0});
  silicon_match->window_dy.set_QoverpT_maxabs({2.3,0,0});
  silicon_match->window_dz.set_QoverpT_range({-2.9,0,0},{4.2,0,0});
  silicon_match->window_deta.set_QoverpT_maxabs({0.06,0,0});
  silicon_match->window_dphi.set_QoverpT_maxabs({0.11,0,0});
  silicon_match->set_test_windows_printout(false);
  silicon_match->set_pp_mode(TRACKING::pp_mode);
  silicon_match->zeroField(true);
  se->registerSubsystem(silicon_match);

  auto mm_match = new PHMicromegasTpcTrackMatching;
  mm_match->Verbosity(0);
  mm_match->set_pp_mode(TRACKING::pp_mode);
  mm_match->set_rphi_search_window_lyr1(6.);
  mm_match->set_rphi_search_window_lyr2(15.0);
  mm_match->set_z_search_window_lyr1(30.0);
  mm_match->set_z_search_window_lyr2(6.);

  mm_match->set_use_silicon(true);
  mm_match->set_min_tpc_layer(38);            // layer in TPC to start projection fit
  mm_match->set_test_windows_printout(true);  // used for tuning search windows only
  mm_match->zeroField(true);
  se->registerSubsystem(mm_match);

  /*
   * Either converts seeds to tracks with a straight line/helix fit
   * or run the full Acts track kalman filter fit
   */
  if (G4TRACKING::convert_seeds_to_svtxtracks)
  {
    auto *converter = new TrackSeedTrackMapConverter;
    // Default set to full SvtxTrackSeeds. Can be set to
    // SiliconTrackSeedContainer or TpcTrackSeedContainer
    converter->setTrackSeedName("SvtxTrackSeedContainer");
    converter->setFieldMap(G4MAGNET::magfield_tracking);
    converter->Verbosity(0);
    converter->constField(); 
    se->registerSubsystem(converter);
  }
  else
  {
    auto *deltazcorr = new PHTpcDeltaZCorrection;
    deltazcorr->Verbosity(0);
    se->registerSubsystem(deltazcorr);

    // perform final track fit with ACTS
    auto *actsFit = new PHActsTrkFitter;
    actsFit->Verbosity(0);
    actsFit->commissioning(G4TRACKING::use_alignment);
    // in calibration mode, fit only Silicons and Micromegas hits
    actsFit->fitSiliconMMs(G4TRACKING::SC_CALIBMODE);
    actsFit->setUseMicromegas(G4TRACKING::SC_USE_MICROMEGAS);
    actsFit->set_pp_mode(TRACKING::pp_mode);
    actsFit->set_use_clustermover(true);  // default is true for now
    actsFit->useActsEvaluator(false);
    actsFit->useOutlierFinder(false);
    actsFit->setFieldMap(G4MAGNET::magfield_tracking);
    se->registerSubsystem(actsFit);

    auto *cleaner = new PHTrackCleaner();
    cleaner->Verbosity(0);
    cleaner->set_pp_mode(TRACKING::pp_mode);
    se->registerSubsystem(cleaner);

    if (G4TRACKING::SC_CALIBMODE)
    {
      /*
       * in calibration mode, calculate residuals between TPC and fitted tracks,
       * store in dedicated structure for distortion correction
       */
      auto *residuals = new PHTpcResiduals;
      const TString tpc_residoutfile = theOutfile + "_PhTpcResiduals.root";
      residuals->setOutputfile(tpc_residoutfile.Data());
      residuals->setUseMicromegas(G4TRACKING::SC_USE_MICROMEGAS);

      // matches Tony's analysis
      residuals->setMinPt(0.2);

      // reconstructed distortion grid size (phi, r, z)
      residuals->setGridDimensions(36, 48, 80);
      se->registerSubsystem(residuals);
    }
  }

  /*
  auto *finder = new PHSimpleVertexFinder;
  finder->Verbosity(0);
  finder->zeroField(true); 
  finder->setDcaCut(0.5);
  finder->setTrackPtCut(0.3);
  finder->setBeamLineCut(1);
  finder->setTrackQualityCut(1000);
  finder->setNmvtxRequired(3);
  finder->setOutlierPairCut(0.1);
  se->registerSubsystem(finder);
  */ 



  if (!G4TRACKING::convert_seeds_to_svtxtracks)
  {
    // Propagate track positions to the vertex position
    auto *vtxProp = new PHActsVertexPropagator;
    vtxProp->Verbosity(0);
    vtxProp->fieldMap(G4MAGNET::magfield_tracking);
    se->registerSubsystem(vtxProp);
  }

  TString residoutfile = theOutfile + "_resid.root";
  std::string residstring(residoutfile.Data());

  auto *resid = new TrackResiduals("TrackResiduals");
  resid->outfileName(residstring);
  resid->alignment(false);
  resid->vertexTree();

  resid->clusterTree();
  resid->hitTree();
  resid->zeroField();
  resid->linefitAll(); 
  resid->convertSeeds(G4TRACKING::convert_seeds_to_svtxtracks);

  //   resid->Verbosity(0);
  std::cout << "!!!!!!!!!!!!!!!!!!!!! Residual output file: " << residstring <<std::endl;
  se->registerSubsystem(resid);

  // auto ntuplizer = new TrkrNtuplizer("TrkrNtuplizer");
  // se->registerSubsystem(ntuplizer);

  /*
    // To write an output DST
    TString dstfile = theOutfile;
   std::string theDSTFile = dstfile.Data();
   Fun4AllOutputManager *out = new Fun4AllDstOutputManager("DSTOUT", theDSTFile.c_str());
   out->AddNode("Sync");
   out->AddNode("EventHeader");
   out->AddNode("TRKR_CLUSTER");
   out->AddNode("SiliconTrackSeedContainer");
   out->AddNode("TpcTrackSeedContainer");
   out->AddNode("SvtxTrackSeedContainer");
   out->AddNode("SvtxTrackMap");
   out->AddNode("SvtxVertexMap");
   out->AddNode("MbdVertexMap");
   out->AddNode("GL1RAWHIT");
   se->registerOutputManager(out);

  */
  if (Enable::QA)
  {
    //se->registerSubsystem(new MvtxRawHitQA);
    //se->registerSubsystem(new InttRawHitQA);
    se->registerSubsystem(new TpcRawHitQA);
    //se->registerSubsystem(new MvtxClusterQA);
    //se->registerSubsystem(new InttClusterQA);
    se->registerSubsystem(new TpcClusterQA);
    //se->registerSubsystem(new MicromegasClusterQA);

 /*
    auto *converter = new TrackSeedTrackMapConverter("SiliconSeedConverter");
    // Default set to full SvtxTrackSeeds. Can be set to
    // SiliconTrackSeedContainer or TpcTrackSeedContainer
    converter->setTrackSeedName("SiliconTrackSeedContainer");
    converter->setTrackMapName("SiliconSvtxTrackMap");
    converter->setFieldMap(G4MAGNET::magfield_tracking);
    converter->Verbosity(0);
    se->registerSubsystem(converter);

    auto *finder_svx = new PHSimpleVertexFinder("SiliconVertexFinder");
    finder_svx->Verbosity(0);
    finder_svx->setDcaCut(0.1);
    finder_svx->setTrackPtCut(0.2);
    finder_svx->setBeamLineCut(1);
    finder_svx->setTrackQualityCut(500);
    finder_svx->setNmvtxRequired(3);
    finder_svx->setOutlierPairCut(0.1);
    finder_svx->setTrackMapName("SiliconSvtxTrackMap");
    finder_svx->setVertexMapName("SiliconSvtxVertexMap");
    se->registerSubsystem(finder_svx);

    auto *siliconqa = new SiliconSeedsQA;
    siliconqa->setTrackMapName("SiliconSvtxTrackMap");
    siliconqa->setVertexMapName("SiliconSvtxVertexMap");
    se->registerSubsystem(siliconqa);
    */

    auto *convertertpc = new TrackSeedTrackMapConverter("TpcSeedConverter");
    // Default set to full SvtxTrackSeeds. Can be set to
    // SiliconTrackSeedContainer or TpcTrackSeedContainer
    convertertpc->setTrackSeedName("TpcTrackSeedContainer");
    convertertpc->setTrackMapName("TpcSvtxTrackMap");
    convertertpc->setFieldMap(G4MAGNET::magfield_tracking);
    convertertpc->Verbosity(0);
    se->registerSubsystem(convertertpc);

    auto *findertpc = new PHSimpleVertexFinder("TpcSimpleVertexFinder");
    findertpc->Verbosity(0);
    findertpc->setDcaCut(0.5);
    //findertpc->setTrackPtCut(0.2);
    findertpc->zeroField(true);
    findertpc->setBeamLineCut(1);
    findertpc->setBeamSpotCutX(-1, 1);
    findertpc->setBeamSpotCutY(-1, 1);
    findertpc->setTrackQualityCut(1000000000);
    // findertpc->setNmvtxRequired(3);
    findertpc->setRequireMVTX(false);
    findertpc->setOutlierPairCut(0.1);
    findertpc->setTrackMapName("TpcSvtxTrackMap");
    findertpc->setVertexMapName("TpcSvtxVertexMap");
    se->registerSubsystem(findertpc);

    auto *tpcqa = new TpcSeedsQA;
    tpcqa->setTrackMapName("TpcSvtxTrackMap");
    tpcqa->setVertexMapName("TpcSvtxVertexMap");
    tpcqa->setSegment(rc->get_IntFlag("RUNSEGMENT"));
    se->registerSubsystem(tpcqa);

   // se->registerSubsystem(new TpcSiliconQA);
  }

  
  //se->skip(nSkip);
  se->run(nEvents+nSkip);
  se->End();
  se->PrintTimer();


  if (Enable::QA)
  {
    std::string qaOutputFileName = theOutfile + "_qa.root";
    std::cout << "!!!!!!!!!!!!!!!!!!!!! QA output file: " << qaOutputFileName <<std::endl;
    QAHistManagerDef::saveQARootFile(qaOutputFileName);
  }
  CDBInterface::instance()->Print();
  delete se;
  std::cout << "Finished" << std::endl;
  gSystem->Exit(0);
}