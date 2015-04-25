/*
    Copyright (c) 2011 Alin Marin Elena <alinm.elena@gmail.com>

    Permission is hereby granted, free of charge, to any person
    obtaining a copy of this software and associated documentation
    files (the "Software"), to deal in the Software without
    restriction, including without limitation the rights to use,
    copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the
    Software is furnished to do so, subject to the following
    conditions:

    The above copyright notice and this permission notice shall be
    included in all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
    OTHER DEALINGS IN THE SOFTWARE.
*/


#include <ctime>
#include "gutenberg.h"
#include "colours.h"
#include "utils.h"
#include "collectivevariable.h"
using namespace std;

void gutenberg::printAtoms ( vector<atom> const& a, ostream& o ) {

    for ( int i=0; i<a.size(); i++ ) {
        utils::Print ( i+1,o );
        a[i].print ( o );
        o<<endl;
        }
    }

void gutenberg::printForces ( vector<atom> const& a, ostream& o ) {
    double fx=0.0,fy=0.0,fz=0.0,f;
    for ( int i=0; i<2* ( utils::sw+1 ) +4* ( utils::w+1 ); i++ ) o<<"=";
    o<<endl;
    o<<colours::green;
    utils::Print ( "|No",o,utils::sw );
    utils::Print ( "|El",o,utils::sw );
    utils::Print ( "|Fx",o );
    utils::Print ( "|Fy",o );
    utils::Print ( "|Fz",o );
    utils::Print ( "|F",o );
    o<<colours::end;
    o<<endl;
    for ( int i=0; i<2* ( utils::sw+1 ) +4* ( utils::w+1 ); i++ ) o<<"-";
    o<<endl;
    for ( int i=0; i<a.size(); i++ ) {
        utils::Print ( i+1,o,utils::sw );
        a[i].printForce ( o );
        fx+=a[i].rfx();
        fy+=a[i].rfy();
        fz+=a[i].rfz();
        utils::Print ( a[i].force(),o );
        o<<endl;
        }
    for ( int i=0; i<2* ( utils::sw+1 ) +4* ( utils::w+1 ); i++ ) o<<"-";
    o<<endl;
    o.width ( utils::sw+utils::sw+2 );
    o<<"Total force: ";
    if ( fx<utils::zerotol ) {
        utils::Print ( fx,o,colours::green );
        }
    else {
        utils::Print ( fx,o,colours::red );
        }
    if ( fy<utils::zerotol ) {
        utils::Print ( fy,o,colours::green );
        }
    else {
        utils::Print ( fy,o,colours::red );
        }
    if ( fz<utils::zerotol ) {
        utils::Print ( fz,o,colours::green );
        }
    else {
        utils::Print ( fz,o,colours::red );
        }
    o<<endl;
    for ( int i=0; i<2* ( utils::sw+1 ) +4* ( utils::w+1 ); i++ ) o<<"=";
    o<<endl;
    }


int gutenberg::printXYZ ( string const& filename , string const& label, std::vector<atom> const& a, bool const& append ) {
    ofstream file;

    if ( append ) {
        file.open ( filename.c_str(),ios::out | ios::app );
        }
    else {
        file.open ( filename.c_str(),ios::out );
        }
    if ( file.is_open() ) {
        file << a.size() <<endl;
        file << label<<endl;
        for ( int  i=0; i<a.size(); i++ ) {
            a[i].printXYZ ( file );
            file<<endl;
            }
        file.close();
        return 0;
        }
    else {
        return -1;
        }
    }

void gutenberg::mdHeader ( ostream& o, bool const& fancy ) {
    if ( fancy ) {
        o << colours::green;
        }
    utils::Print ( "#Time",o );
    utils::Print ( "|E_pot/N",o );
    utils::Print ( "|E_kin/N",o );
    utils::Print ( "|iT",o );
    utils::Print ( "|T",o );
    utils::Print ( "|p",o );
    utils::Print ( "|Q1 (E)/N",o );
    utils::Print ( "|Q2(Vcm)",o );
    if ( fancy ) {
        o<<colours::end;
        }
    cout<<endl;
    }

void gutenberg::mdRepLine ( const double & step , const double & ep, const double & ek,  const double & it, const double & tt,const double & pressure, const double & q1, const double & q2, const bool& fancy, ostream& o ) {
    utils::Print ( step,o );
    utils::Print ( ep,o );
    utils::Print ( ek,o );
    utils::Print ( it,o );
    utils::Print ( tt,o );
    utils::Print ( pressure,o );
    if ( fancy ) {
        o<<colours::red;
        }
    utils::Print ( q1,o );
    utils::Print ( q2,o );
    if ( fancy ) {
        o<<colours::end;
        }
    o<<endl;
    }


void gutenberg::mdTAMDHeader ( const vector<colVar>&z, ostream& o, bool const& fancy ) {
    if ( fancy ) {
        o << colours::green;
        }
    utils::Print ( "#Time",o );
    utils::Print ( "|E_pot/N",o );
    utils::Print ( "|E_kin/N",o );
    utils::Print ( "|T_X",o );
    utils::Print ( "|T_Z",o );
    utils::Print ( "|pressure",o );
    utils::Print ( "|K_Z",o );
    for ( int i=0; i<z.size(); i++ ) utils::Print ( utils::stringPlusInt ( "|z_",i+1 ),o );
    for ( int i=0; i<z.size(); i++ ) utils::Print ( utils::stringPlusInt ( "|theta_",i+1 ),o );
    utils::Print ( "|Q1 (E)/N",o );
    utils::Print ( "|Q2(Vcm)",o );

    if ( fancy ) {
        o<<colours::end;
        }
    cout<<endl;
    }

void gutenberg::mdTAMDRepLine ( const double & step , const double & ep, const double & ek,  const double & tt, const double & tz, const double& pressure,const double& kinz,const double & q1, const double & q2,const vector<colVar>&z, const bool& fancy, ostream& o ) {
    utils::Print ( step,o );
    utils::Print ( ep,o );
    utils::Print ( ek,o );
    utils::Print ( tt,o );
    utils::Print ( tz,o );
    utils::Print ( pressure,o );
    utils::Print ( kinz,o );
    if ( fancy ) {
        o<<colours::purple;
        }
    for ( int i=0; i<z.size(); i++ ) utils::Print ( z[i].z,o );
    for ( int i=0; i<z.size(); i++ ) utils::Print ( z[i].th,o );
    if ( fancy ) {
        o<<colours::end;
        }
    if ( fancy ) {
        o<<colours::red;
        }
    utils::Print ( q1,o );
    utils::Print ( q2,o );
    if ( fancy ) {
        o<<colours::end;
        }
    o<<endl;
    }



int gutenberg::printHistogram ( std::vector<double> const& h, double const& min, double const&step,string const& filename ) {
    ofstream file ( filename.c_str() );
    if ( file.is_open() ) {
        for ( int i=0; i<h.size(); i++ ) {
            utils::Print ( min+i*step,file );
            utils::Print ( h[i],file );
            file<<endl;
            }
        file.close();
        return 0;
        }
    else {
        return -1;
        }


    }


void gutenberg::MMCHeader ( ostream& o, bool const& fancy ) {
    if ( fancy ) {
        o << colours::green;
        }
    utils::Print ( "#MC Steps",o );
    utils::Print ( "|E_pot/N",o );
    utils::Print ( "|p",o );
    if ( fancy ) {
        o<<colours::end;
        }
    cout<<endl;
    }

void gutenberg::MMCRepLine ( const int& i, const double& ep, const double& pressure, const bool& fancy, ostream& o ) {
    utils::Print ( i,o );
    utils::Print ( ep,o );
    if ( fancy ) {
        utils::Print ( pressure,o,colours::red );
        }
    else {
        utils::Print ( pressure,o );
        }
    o<<endl;
    }

void gutenberg::mdTAMCHeader ( const vector<colVar>&z, ostream& o, bool const& fancy ) {
    if ( fancy ) {
        o << colours::green;
        }
    utils::Print ( "#Time",o );
    utils::Print ( "|E_pot/N",o );
    utils::Print ( "|T_Z",o );
    utils::Print ( "|pressure",o );
    utils::Print ( "|KZ",o );
    for ( int i=0; i<z.size(); i++ ) utils::Print ( utils::stringPlusInt ( "|z_",i+1 ),o );
    for ( int i=0; i<z.size(); i++ ) utils::Print ( utils::stringPlusInt ( "|theta_",i+1 ),o );
    if ( fancy ) {
        o<<colours::end;
        }
    cout<<endl;
    }

void gutenberg::mdTAMCRepLine ( const double & step , const double & ep, const double& tz, const double& pressure,const double& kinz,const vector<colVar>&z, const bool& fancy, ostream& o ) {
    utils::Print ( step,o );
    utils::Print ( ep,o );
    utils::Print ( tz,o );
    utils::Print ( pressure,o );
    utils::Print ( kinz,o );
    if ( fancy ) {
        o<<colours::purple;
        }
    for ( int i=0; i<z.size(); i++ ) utils::Print ( z[i].z,o );
    for ( int i=0; i<z.size(); i++ ) utils::Print ( z[i].th,o );
    if ( fancy ) {
        o<<colours::end;
        }
    if ( fancy ) {
        o<<colours::red;
        }
    if ( fancy ) {
        o<<colours::end;
        }
    o<<endl;
    }
void gutenberg::logInput ( generalInputs& in ) {

    gutenberg::logLine ( "system name:",in.name,in.log );
    gutenberg::logLine ( "sampler: ",in.sampler,in.log );
    gutenberg::logLine ( "logfile: ",in.logf,in.log );
    gutenberg::logLine ( "debugfile: ",in.dbgf,in.log );
    gutenberg::logLine ( "timestep: ",in.dt,in.log );
    gutenberg::logLine ( "NSteps: ",in.nSteps,in.log );
    gutenberg::logLine ( "Initial T: ",in.T0,in.log );
    gutenberg::logLine ( "Bath T: ",in.TBath,in.log );
    gutenberg::logLine ( "CV T: ",in.Tcv,in.log );
    gutenberg::logLine ( "Gamma: ",in.gamma,in.log );
    gutenberg::logLine ( "CV Gamma: ",in.gCV,in.log );
    gutenberg::logLine ( "Sampling Freq: ",in.frequency,in.log );
    gutenberg::logLine ( "Hist bin: ",in.bin,in.log );
    gutenberg::logLine ( "LJ cutoff: ",in.cut,in.log );
    gutenberg::logLine ( "sigma: ",utils::sigma,in.log );
    gutenberg::logLine ( "epsilon: ",utils::epsilon,in.log );
    gutenberg::logLine ( "E cutoff: ",utils::LennardJones ( in.cut ),in.log );
    gutenberg::logLine ( "use shift: ",in.shift,in.log );

    gutenberg::logLine ( "MC cycles: ",in.mccycles,in.log );
    gutenberg::logLine ( "MC Steps: ",in.mcsteps,in.log );
    gutenberg::logLine ( "MC displacement: ",in.dr,in.log );
    gutenberg::logLine ( "XYZ input: ",in.filename,in.log );
    gutenberg::logLine ( "XYZ output: ",in.xyz,in.log );
    gutenberg::logLine ( "Coloured output: ",in.fancy,in.log );
    gutenberg::logLine ( "Element name: ",in.el,in.log );
    utils::Print ( "random seeds: ",in.log ,colours::green );
    for ( int i=0; i<in.seeds.size(); i++ ) utils::Print ( in.seeds[i],in.log,colours::red );
    in.log<<endl;
    gutenberg::logLine ( "box size",in.box,in.log );

    }

void gutenberg::logLine ( const string& s, const double& a, ostream& o ) {
    gutenberg::logInfo ( s,a,o );
    o<<endl;
    }

void gutenberg::logLine ( const string& s, const int& a, ostream& o ) {
    gutenberg::logInfo ( s,a,o );
    o<<endl;
    }
void gutenberg::logLine ( const string& s, const string& a, ostream& o ) {
    gutenberg::logInfo ( s,a,o );
    o<<endl;
    }

void gutenberg::logLine ( const string& s, const string& a, ostream& o, const int& w ) {
    gutenberg::logInfo ( s,a,o,w );
    o<<endl;
    }

void gutenberg::logInfo ( const string& s, const double& a, ostream& o ) {
    utils::Print ( s,o,colours::green );
    utils::Print ( a,o,colours::red );
    }

void gutenberg::logInfo ( const string& s, const int& a, ostream& o ) {
    utils::Print ( s,o,colours::green );
    utils::Print ( a,o,colours::red );
    }
void gutenberg::logInfo ( const string& s, const string& a, ostream& o ) {
    utils::Print ( s,o,colours::green );
    utils::Print ( a,o,colours::red );
    }

void gutenberg::logInfo ( const string& s, const string& a, ostream& o, const int& w ) {
    utils::Print ( s,o,colours::green );
    utils::Print ( a,o,w,colours::red );
    }


void gutenberg::timeStamp ( string const& s,ostream& o ) {
    time_t t;
    time ( &t );
    tm *ct= localtime ( &t );
    gutenberg::logLine ( s, asctime ( ct ),o,40 );

    }

void gutenberg::logCV (  vector< colVar >& cv, ostream& o ) {

    utils::Print ( "|no",o,colours::green );
    utils::Print ( "|z",o,colours::green );
    utils::Print ( "|mu",o,colours::green );
    utils::Print ( "|k",o,colours::green );
    utils::Print ( "|gamma",o,colours::green );
    utils::Print ( "|v",o,colours::green );
    utils::Print ( "|f",o,colours::green );
    utils::Print ( "|th",o,colours::green );
    o<<endl;
    int k=1;
    for ( vector<colVar>::iterator it=cv.begin(); it<cv.end(); it++ ) {
        utils::Print ( k,o );
        utils::Print ( ( *it ).z,o );
        utils::Print ( ( *it ).mu,o );
        utils::Print ( ( *it ).k,o );
        utils::Print ( ( *it ).gamma,o );
        utils::Print ( ( *it ).v,o );
        utils::Print ( ( *it ).f,o );
        utils::Print ( ( *it ).th,o );o<<endl;
        k++;
        }

    }

// kate: indent-mode cstyle; space-indent on; indent-width 4;
