#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>

#include "atom.h"
#include "utils.h"
#include "solution.h"
#include "colours.h"
#include "gutenberg.h"
#include "generalinputs.h"
#include "collectivevariable.h"



using namespace std;

int main ( int argc, char **argv ) {
    generalInputs inps;

    inps.fancy=true;
    inps.nSteps=100000;
    inps.box=7.55952;
    inps.TBath=1.0;
    inps.dt=0.001;
    inps.frequency=100;
    inps.cut=2.5;
    inps.filename="ar216.xyz";
    inps.xyz="test.xyz";
    inps.name="LJ";
    inps.bin=0.001;
    inps.Tcv=10.0;
    inps.mccycles=1000;
    inps.mcsteps=200;
    inps.dr=0.2;
    inps.T0=1.0;
    inps.TBath=1.0;
    inps.sampler="VV";
    inps.gamma=1.0;
    inps.gCV=10.0;
    inps.logf="mcplay.log";
    inps.dbgf="mcplay.debug";


    for ( int i=1; i<argc; i++ ) {
        if ( !strcmp ( argv[i],"-steps" ) ) inps.nSteps=atoi ( argv[++i] );
        else if ( !strcmp ( argv[i],"-L" ) ) inps.box=atof ( argv[++i] );
        else if ( !strcmp ( argv[i],"-seed" ) ) inps.seeds.push_back ( atoi ( argv[++i] ) );
        else if ( !strcmp ( argv[i],"-dt" ) ) inps.dt=atof ( argv[++i] );
        else if ( !strcmp ( argv[i],"-freq" ) ) inps.frequency=atoi ( argv[++i] );
        else if ( !strcmp ( argv[i],"-xyz" ) ) inps.xyz=argv[++i];
        else if ( !strcmp ( argv[i],"-inp" ) ) inps.filename=argv[++i];
        else if ( !strcmp ( argv[i],"-cut" ) ) inps.cut=atof ( argv[++i] );
        else if ( !strcmp ( argv[i],"-T0" ) ) inps.T0=atof ( argv[++i] );
        else if ( !strcmp ( argv[i],"-Tb" ) ) inps.TBath=atof ( argv[++i] );
        else if ( !strcmp ( argv[i],"-Tvc" ) ) inps.Tcv=atof ( argv[++i] );
        else if ( !strcmp ( argv[i],"-gamma" ) ) inps.gamma=atof ( argv[++i] );
        else if ( !strcmp ( argv[i],"-gCV" ) ) inps.gCV=atof ( argv[++i] );
        else if ( !strcmp ( argv[i],"-bin" ) ) inps.bin=atof ( argv[++i] );
        else if ( !strcmp ( argv[i],"-fancy" ) ) inps.fancy= atoi (argv[++i])!=0?true:false;
        else if ( !strcmp ( argv[i],"-samp" ) ) inps.sampler=argv[++i];
        else if ( !strcmp ( argv[i],"-log" ) ) inps.logf=argv[++i];
        else if ( !strcmp ( argv[i],"-dbg" ) ) inps.dbgf=argv[++i];
        else if ( !strcmp ( argv[i],"-name" ) ) inps.name=argv[++i];
        else if ( !strcmp ( argv[i],"-mcsteps" ) ) inps.mcsteps=atoi ( argv[++i] );
        else if ( !strcmp ( argv[i],"-mccycles" ) ) inps.mccycles=atoi ( argv[++i] );
        else if ( !strcmp ( argv[i],"-dr" ) ) inps.dr=atof ( argv[++i] );
        else if ( !strcmp ( argv[i],"-h" ) ) {
            utils::Print ( "usage: ", cout,colours::red );
            cout<<endl;
            utils::Print ( argv[0],cout,colours::green );
            cout <<endl;
            utils::Print ( "-steps",cout,colours::green );
            utils::Print ( "<int> number of steps",cout,colours::purple );
            cout<<endl;
            utils::Print ( "-L",cout,colours::green );
            utils::Print ( "<double> box size",cout,colours::purple );
            cout<<endl;
            utils::Print ( "-seed",cout,colours::green );
            utils::Print ( "<int> seed for random generator",cout,colours::purple );
            cout<<endl;
            utils::Print ( "-dt" ,cout,colours::green );
            utils::Print ( "<double> timestep for integrator",cout,colours::purple );
            cout<<endl;
            utils::Print ( "-freq",cout,colours::green );
            utils::Print ( "<int> frequency to print xyz files",cout,colours::purple );
            cout<<endl;
            utils::Print ( "-xyz" ,cout,colours::green );
            utils::Print ( "<string> name of the file for animation",cout,colours::purple );
            cout<<endl;
            utils::Print ( "-inp" ,cout,colours::green );
            utils::Print ( "<string> name of the file with initial postitions(xyz fmt)",cout,colours::purple );
            cout<<endl;
            utils::Print ( "-cut" ,cout,colours::green );
            utils::Print ( "<double> cutoff for LJ",cout,colours::purple );
            cout<<endl;
            utils::Print ( "-bin" ,cout,colours::green );
            utils::Print ( "<double> width of bin in histograms",cout,colours::purple );
            cout<<endl;
            utils::Print ( "-T0" ,cout,colours::green );
            utils::Print ( "<double> initial temperature of the sample",cout,colours::purple );
            cout<<endl;
            utils::Print ( "-Tb" ,cout,colours::green );
            utils::Print ( "<double> bath temperature",cout,colours::purple );
            cout<<endl;
            utils::Print ( "-Tcv" ,cout,colours::green );
            utils::Print ( "<double> collective variables bath temperature",cout,colours::purple );
            cout<<endl;
            utils::Print ( "-gamma" ,cout,colours::green );
            utils::Print ( "<double> friction coefficient ",cout,colours::purple );
            cout<<endl;
            utils::Print ( "-gCV" ,cout,colours::green );
            utils::Print ( "<double> friction coefficient for CV",cout,colours::purple );
            cout<<endl;
            utils::Print ( "-fancy" ,cout,colours::green );
            utils::Print ( "<int> as logical,  enable/discable colours",cout,colours::purple );
            cout<<endl;
            utils::Print ( "-samp" ,cout,colours::green );
            utils::Print ( "<string> sampler available VV, VEC, TAMD, MMC, TAMC",cout,colours::purple );
            cout<<endl;
            utils::Print ( "-log" ,cout,colours::green );
            utils::Print ( "<string> name for log file ",cout,colours::purple );
            cout<<endl;
            utils::Print ( "-dbg" ,cout,colours::green );
            utils::Print ( "<string> name for dbg file ",cout,colours::purple );
            cout<<endl;
            utils::Print ( "-name" ,cout,colours::green );
            utils::Print ( "<string> systems name ",cout,colours::purple );
            cout<<endl;
            utils::Print ( "-mcsteps",cout,colours::green );
            utils::Print ( "<int> number of monte carlo moves",cout,colours::purple );
            cout<<endl;
            utils::Print ( "-mccycles",cout,colours::green );
            utils::Print ( "<int> number of monte carlo cycles",cout,colours::purple );
            cout<<endl;
            utils::Print ( "-dr" ,cout,colours::green );
            utils::Print ( "<double> maximum MC displacement in one direction",cout,colours::purple );
            cout<<endl;
            utils::Print ( "-h" ,cout,colours::green );
            utils::Print ( "this help ",cout,colours::purple );
            cout<<endl;
            return -1;
            }
        }
        if (inps.seeds.size()==0) inps.seeds.push_back ( 12345 );

    vector<atom> a;
//     vector<colVar> cv(2);

    inps.log.open ( inps.logf.c_str(),ios::out );
    inps.debug.open ( inps.dbgf.c_str(),ios::out );

    gutenberg::timeStamp("start time:",inps.log);
    if ( inps.sampler== "VV" ) {
        solution::InitializeSimulation ( a,inps.T0,inps.box,inps.seeds[0], inps.filename );
        gutenberg::logInput(inps);
        solution::MDDriverVV ( a,inps );
        }
    else if ( inps.sampler== "VEC" ) {
        solution::InitializeSimulation ( a,inps.T0,inps.box,inps.seeds[0], inps.filename,inps.gamma );
        gutenberg::logInput(inps);
        solution::MDDriverVEC ( a,inps );
        }
    else if ( inps.sampler== "TAMD" ) {
        solution::InitializeSimulation ( a,inps.T0,inps.box,inps.seeds[0], inps.filename,inps.gamma );
        vector<colVar> cv ( 2 );
        cv[0].mu=10.0;
        cv[1].mu=10.0;
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
        gutenberg::logInput(inps);
        solution::MDDriverTAMD ( a,cv,inps );
        }
    else if ( inps.sampler== "MMC" ) {
        solution::InitializeSimulation ( a,inps.T0,inps.box,inps.seeds[0], inps.filename,inps.gamma );
        // inps.Temperature should contain the temperature of the system... so be sure it is computed if not an input
        gutenberg::logInput(inps);
        solution::MMCDriver ( a,inps );
        }
    else if ( inps.sampler== "TAMC" ) {
        solution::InitializeSimulation ( a,inps.T0,inps.box,inps.seeds[0], inps.filename,inps.gamma );
        vector<colVar> cv ( 2 );
        cv[0].mu=10.0;
        cv[1].mu=10.0;
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
        gutenberg::logInput(inps);
        solution::MDDriverTAMC ( a,cv,inps );
        }
    else {
        utils::Print ( "Unknown sampler\n", cout,colours::red );
        }
    gutenberg::timeStamp("end time:",inps.log);
    inps.log.close();
    inps.debug.close();
    return 0;
    }
// kate: indent-mode cstyle; space-indent on; indent-width 4;
