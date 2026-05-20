#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllDstInputManager.h>

#include </sphenix/user/mitrankova/F4A/InModuleTrackDisplay/install/include/inmoduletrackdisplay/InModuleTrackDisplay.h>

#include <TSystem.h>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libphool.so)
R__LOAD_LIBRARY(/sphenix/user/mitrankova/F4A/TPC_pattern_reco/install/lib/libInModuleTracks.so)
R__LOAD_LIBRARY(/sphenix/user/mitrankova/F4A/InModuleTrackDisplay/install/lib/libInModuleTrackDisplay.so)

namespace
{
  std::vector<std::string> read_file_list_or_single(const std::string& input)
  {
    std::vector<std::string> files;
    const bool looks_like_list = input.size() > 5 &&
      (input.substr(input.size() - 5) == ".list" || input.substr(input.size() - 4) == ".txt");

    if (looks_like_list)
    {
      std::ifstream fin(input.c_str());
      if (!fin)
      {
        std::cerr << "Cannot open input list: " << input << std::endl;
        return files;
      }
      std::string line;
      while (std::getline(fin, line))
      {
        if (line.empty()) continue;
        if (line[0] == '#') continue;
        files.push_back(line);
      }
    }
    else
    {
      files.push_back(input);
    }
    return files;
  }
}

void Fun4All_InModuleTrackDisplay(const std::string& input = "my_dst_list.txt",
                                  const std::string& output = "inmodule_display.root",
                                  const int nEvents = 5,
                                  const std::string& trackNodeName = "INMODULETRACKS",
                                  const unsigned int maxEventDisplays = 5)
{
  gSystem->Load("libInModuleTracks.so");
  gSystem->Load("libInModuleTrackDisplay.so");

  Fun4AllServer* se = Fun4AllServer::instance();
  se->Verbosity(1);

  Fun4AllDstInputManager* in = new Fun4AllDstInputManager("DSTin");
  std::vector<std::string> files = read_file_list_or_single(input);
  std::cout << "Input files: " << files.size() << std::endl;
  for (std::vector<std::string>::const_iterator it = files.begin(); it != files.end(); ++it)
  {
    std::cout << "  adding " << *it << std::endl;
    in->AddFile(*it);
  }
  se->registerInputManager(in);

  se->registerSubsystem( new InModuleTrackDisplay());



  if (nEvents > 0) se->run(nEvents);
  else se->run();

  se->End();
  delete se;
  gSystem->Exit(0);
}
