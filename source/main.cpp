#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>

#include "atom.h"
#include "utils.h"
#include "solution.h"
#include "colours.h"
#include "gutenberg.h"
#include "generalinputs.h"
#include "collectivevariable.h"
#include "reading.h"
#include "parseCLI.h"
#include "customException.hpp"


using namespace std;

int main ( int argc, char **argv ) {
    generalInputs inps;

  try {
    boost::program_options::variables_map clioMap=parseCLiandConfig(argc,argv,&inps);
    if (clioMap.count("help")) return SUCCESS;
    if (clioMap.count("version")) {
      std::cout<<colours::green<<"version: "<<inps.version<<colours::end<<"\n";
      return SUCCESS;
    }

    std::cout<<colours::green<<"Echo the options used\n"<<colours::end;
    std::string debugFile=inps.configFile+".dbg";
    printOptions(&clioMap,&debugFile);
  }
  catch(customException &e) {
    std::cerr<<colours::red<<"Error: "<< e.what() << colours::end<<"\n";
    return e.id;
  }

    if ( inps.seeds.size() ==0 ) inps.seeds.push_back ( inps.seed );
    inps.cut=min<double> ( inps.cut,inps.box/2.0 );
    vector<atom> a;
    inps.log.open ( inps.logf.c_str(),ios::out );
    gutenberg::timeStamp ( "start time:",inps.log );
    for ( int i=0; i<argc; i++ ) {
        utils::Print ( argv[i],inps.log,colours::yellow );
        }
    inps.log<<endl;
    if ( inps.sampler== "VV" ) {
        solution::InitializeSimulation ( a,inps.T0,inps.box,inps.seeds[0], inps.filename );
        gutenberg::logInput ( inps );
        solution::MDDriverVV ( a,inps );
        }
    else if ( inps.sampler== "VEC" ) {
        solution::InitializeSimulation ( a,inps.T0,inps.box,inps.seeds[0], inps.filename,inps.gamma );
        gutenberg::logInput ( inps );
        solution::MDDriverVEC ( a,inps );
        }
    else if ( inps.sampler== "TAMD" ) {
        solution::InitializeSimulation ( a,inps.T0,inps.box,inps.seeds[0], inps.filename,inps.gamma );
        vector<colVar> cv ( 2 );
        cv[0].mu=inps.mu;
        cv[0].gamma=inps.gCV;
        cv[0].k=inps.zk;
        cv[0].i=1;
        cv[0].j=2;
        cv[0].z=cv[0].theta ( a,inps.box );
        cv[0].th=cv[0].z;
        cv[0].v=0.0;
        cv[0].f=0.0;

        cv[1].mu=inps.mu;
        cv[1].gamma=inps.gCV;
        cv[1].k=inps.zk;
        cv[1].i=5;
        cv[1].j=6;
        cv[1].z=cv[1].theta ( a,inps.box );
        cv[1].th=cv[1].z;
        cv[1].v=0.0;
        cv[1].f=0.0;
        gutenberg::logInput ( inps );
        gutenberg::logCV ( cv,inps.log );
        solution::MDDriverTAMD ( a,cv,inps );
        }
    else if ( inps.sampler== "MMC" ) {
        solution::InitializeSimulation ( a,inps.T0,inps.box,inps.seeds[0], inps.filename,inps.gamma );
        // inps.Temperature should contain the temperature of the system... so be sure it is computed if not an input
        gutenberg::logInput ( inps );
//        double vir;
//         int i=105;
//         double en=a[i].myPotentialEnergy ( a,inps.cut,inps.box,inps.shift,vir);
//         utils::Print(en,cout,32,16);
//         cout <<endl;
//         vector<atom> b;
//
//         int l=reading::readXYZ("fort.111",b);
//         for(int k=0;k<b.size();k++) b[k].setID(k);
//         double em=b[i].myPotentialEnergy ( a,inps.cut,inps.box,inps.shift,vir);
//         utils::Print(em,cout,32,16);cout<<endl;
        solution::MMCDriver ( a,inps );
        }
    else if ( inps.sampler== "TAMC" ) {
        solution::InitializeSimulation ( a,inps.T0,inps.box,inps.seeds[0], inps.filename,inps.gamma );
        vector<colVar> cv ( 2 );
        cv[0].mu=inps.mu;
        cv[1].mu=inps.mu;
        cv[0].gamma=inps.gCV;
        cv[1].gamma=inps.gCV;
        cv[0].k=50.0;
        cv[1].k=50.0;
        cv[0].i=1;
        cv[0].j=2;
        cv[1].i=5;
        cv[1].j=6;
        cv[0].z=cv[0].theta ( a,inps.box );
        cv[1].z=cv[1].theta ( a,inps.box );
        cv[0].th=cv[0].z;
        cv[1].th=cv[1].z;
        cv[0].v=0.0;
        cv[1].v=0.0;
        cv[0].f=0.0;
        cv[1].f=0.0;
        gutenberg::logInput ( inps );
        solution::MDDriverTAMC ( a,cv,inps );
        }
    else if ( inps.sampler== "TAHMC" ) {
        solution::InitializeSimulation ( a,inps.T0,inps.box,inps.seeds[0], inps.filename,inps.gamma );
        vector<colVar> cv ( 2 );
        cv[0].mu=inps.mu;
        cv[1].mu=inps.mu;
        cv[0].gamma=inps.gCV;
        cv[1].gamma=inps.gCV;
        cv[0].k=50.0;
        cv[1].k=50.0;
        cv[0].i=1;
        cv[0].j=2;
        cv[1].i=5;
        cv[1].j=6;
        cv[0].z=cv[0].theta ( a,inps.box );
        cv[1].z=cv[1].theta ( a,inps.box );
        cv[0].th=cv[0].z;
        cv[1].th=cv[1].z;
        cv[0].v=0.0;
        cv[1].v=0.0;
        cv[0].f=0.0;
        cv[1].f=0.0;
        gutenberg::logInput ( inps );
        solution::MDDriverTAHMC ( a,cv,inps );
        }    
    else {
        solution::GenerateLattice ( a,inps.rho,inps.mass, inps.N,inps.box,inps.el );
        gutenberg::printXYZ ( inps.xyz,utils::stringPlusDouble ( "LJ box ",inps.box ),a,true );

        utils::Print ( "Lattice generated\n", cout,colours::yellow );
        utils::Print ( "use -h for more info on usage\n", cout,colours::red );
        }
    gutenberg::timeStamp ( "end time:",inps.log );
    inps.log.close();
    return 0;
    }
// kate: indent-mode cstyle; space-indent on; indent-width 4;
