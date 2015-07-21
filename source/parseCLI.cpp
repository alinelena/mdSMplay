#include <iostream>
#include <fstream>
#include "colours.h"
#include "customException.hpp"
#include "parseCLI.h"

boost::program_options::variables_map parseCLiandConfig(int ac, char** av, generalInputs *inps)
{
  namespace clio = boost::program_options;
  clio::options_description clionly("Command Line Only Options");
  clionly.add_options()
  ("help,h","Prints help messages")
  ("version,v","Show me the version")
  ("config,c",clio::value<std::string>(&(inps->configFile))->default_value("options.conf"),"Config file from where options can be read");
  clio::options_description clic("Config File Options");
  clic.add_options()
  ("steps",clio::value<int>(&(inps->nSteps))->default_value(100000),"Number of steps")
  ("frequency",clio::value<int>(&(inps->frequency))->default_value(100),"Frequency of printing")
  ("fancy",clio::value<bool>(&(inps->fancy))->default_value(true),"Fancy printing")
  ("box",clio::value<double>(&(inps->box))->default_value(7.55952),"box size")
  ("tbath",clio::value<double>(&(inps->TBath))->default_value(1.0),"Temperature for bath")
  ("timestep",clio::value<double>(&(inps->dt))->default_value(0.001),"Time step")
  ("cutoff",clio::value<double>(&(inps->cut))->default_value(2.5),"Cut off")
  ("filename",clio::value< std::string >(&(inps->filename))->default_value("ar216.xyz"),"Coordinate file name")
  ("xyz",clio::value< std::string >(&(inps->xyz))->default_value("test.xyz"),"Coordinate file name")
  ("name",clio::value< std::string >(&(inps->name))->default_value("LJ"),"Name of the system")
  ("bin",clio::value<double>(&(inps->bin))->default_value(0.001),"Bin size")
  ("tcv",clio::value<double>(&(inps->Tcv))->default_value(10.0),"Temperature for collective variable")
  ("mccycles",clio::value<int>(&(inps->mccycles))->default_value(1000),"MC cycles")
  ("mcsteps",clio::value<int>(&(inps->mcsteps))->default_value(200),"MC steps")
  ("dr",clio::value<double>(&(inps->dr))->default_value(0.2),"dr MC displacement?")
  ("t0",clio::value<double>(&(inps->T0))->default_value(1.0),"Temperature")
  ("gamma",clio::value<double>(&(inps->gamma))->default_value(1.0),"Gamma")
  ("gcv",clio::value<double>(&(inps->gCV))->default_value(10.0),"Gamma CV")
  ("logf",clio::value< std::string >(&(inps->logf))->default_value("mcplay.log"),"Log file name")
  ("debug",clio::value< std::string >(&(inps->dbgf))->default_value("mcplay.debug"),"Debug file name")
  ("shift",clio::value<bool>(&(inps->shift))->default_value(false),"Force shift")
  ("natoms",clio::value<int>(&(inps->N))->default_value(216),"Number of atoms")
  ("nz",clio::value<int>(&(inps->Nz))->default_value(10),"Number of collective variables")
  ("nequil",clio::value<int>(&(inps->nEquil))->default_value(300),"Number of steps for equilibration")
  ("rho",clio::value<double>(&(inps->rho))->default_value(0.5),"density")
  ("mass",clio::value<double>(&(inps->mass))->default_value(1.0),"mass")
  ("mu",clio::value<double>(&(inps->mu))->default_value(1.0),"colllective variable mass")
  ("zk",clio::value<double>(&(inps->zk))->default_value(50.0),"colllective spring constant")
  ("dtz",clio::value<double>(&(inps->dtz))->default_value(0.002),"time step for CV")
  ("el",clio::value< std::string >(&(inps->el))->default_value("Ar"),"Element")
  ("seed",clio::value<int>(&(inps->seed))->default_value(2015),"default seed")
  ("sampler",clio::value< std::string >(&(inps->sampler))->default_value("TAMC"),"Sampler to be employed: VV, VEC, TAMD, MMC, TAMC, TAHMC");
  clio::options_description clics("Command Line and Config File Options");
  clics.add(clionly).add(clic);

  clio::variables_map clioMap;
  try {
    clio::store(clio::parse_command_line(ac,av,clics),clioMap);
    clio::notify(clioMap);
  }
  catch(clio::error &e) {
    std::cerr<<colours::red<<"Error: "<< e.what() << colours::end<<"\n";
    std::cout<<colours::green<<"Command line Options:\n"<<colours::end<<clics<<"\n";
    throw customException(ERROR_CLI);
  }
  try {
    if (inps->configFile == "options.conf") {
      std::ifstream icfile(inps->configFile.c_str());
      if(!icfile) {
        std::cout<<colours::yellow<<"Creating empty config file\n"<<colours::end;
        icfile.close();
        std::fstream icfile(inps->configFile.c_str(),std::ios::out);
        icfile.close();
      }
    }
    clio::store(clio::parse_config_file<char>(inps->configFile.c_str(),clic),clioMap);
    clio::notify(clioMap);
  }
  catch(clio::error &e) {
    std::cerr<<colours::red<<"Error: "<< e.what() << colours::end<<"\n";
    std::cout<<colours::green<<"Config FilesOptions:\n"<<colours::end<<clic<<"\n";
    throw customException(ERROR_CFILE);
  }
  if (clioMap.count("help")) {
    std::cout<<colours::green<<"Command line Options:\n"<<colours::end<<clics<<"\n";
  }
  return clioMap;
}

void printOptions(boost::program_options::variables_map *clioMap, std::string *dbgf)
{
  std::fstream dbg(dbgf->c_str(),std::ios::out);
  for (const auto& it : *clioMap) {
    std::cout << it.first.c_str() << " ";
    dbg << it.first.c_str() << " = ";
    auto& value = it.second.value();
    if (auto v = boost::any_cast<int>(&value)) {
      std::cout << colours::yellow<<*v<<"\n"<<colours::end;
      dbg <<*v<<"\n";
    }
    else if (auto v = boost::any_cast<std::string>(&value)) {
      std::cout << colours::yellow<<*v<<"\n"<<colours::end;
      dbg <<*v<<"\n";
    }
    else if (auto v = boost::any_cast<bool>(&value)) {
      std::cout << colours::yellow<<*v<<"\n"<<colours::end;
      dbg <<*v<<"\n";
    }
    else if (auto v = boost::any_cast<double>(&value)) {
      std::cout << colours::yellow<<*v<<"\n"<<colours::end;
      dbg <<*v<<"\n";
    }
    else {
      std::cout << colours::red<<"No idea about this type\n"<<colours::end;
    }
  }
  dbg.close();
}

