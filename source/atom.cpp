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


#include "atom.h"

#include <iostream>
#include <string>
#include <cmath>
#include "utils.h"
#include "atom.h"
#include <stdlib.h>
using namespace std;

void atom::print ( ostream &o ) const {
    utils::Print ( name,o,utils::sw );
    utils::Print ( mass,o );
    utils::Print ( x,o );
    utils::Print ( y,o );
    utils::Print ( z,o );
    utils::Print ( vx,o );
    utils::Print ( vy,o );
    utils::Print ( vz,o );
    utils::Print ( fx,o );
    utils::Print ( fy,o );
    utils::Print ( fz,o );
    }

void atom::printXYZ ( ostream &o ) const {

    utils::Print ( name,o,utils::sw );
    utils::Print ( x,o );
    utils::Print ( y,o );
    utils::Print ( z,o );
    }


void atom::printForce ( ostream&o ) const {
    utils::Print ( name,o,utils::sw );
    utils::Print ( fx,o );
    utils::Print ( fy,o );
    utils::Print ( fz,o );
    }


void atom::setElement ( string const& el ) {
    name=el;
    }

void atom::setPositions ( double const& p, double const& q, double const& r ) {
    x=p;
    y=q;
    z=r;
    }

void atom::setForces ( const double& p, const double& q, const double& r ) {
    fx=p;
    fy=q;
    fz=r;
    }

void atom::addtoForces ( const double& p, const double& q, const double& r ) {
    fx+=p;
    fy+=q;
    fz+=r;
    }

void atom::setVelocities ( const double& p, const double& q, const double& r ) {
    vx=p;
    vy=q;
    vz=r;
    }

void atom::updatePositionsVV ( const double& dt ) {
    x+=vx*dt+fx/ ( mass*2.0 ) *dt*dt;
    y+=vy*dt+fy/ ( mass*2.0 ) *dt*dt;
    z+=vz*dt+fz/ ( mass*2.0 ) *dt*dt;
    }

void atom::updateVelocitiesVV ( const double& dt ) {
    vx+= ( fx ) / ( mass*2.0 ) *dt;
    vy+= ( fy ) / ( mass*2.0 ) *dt;
    vz+= ( fz ) / ( mass*2.0 ) *dt;
    }

void atom::setMass ( const double& m ) {
    mass=m;
    }


void atom::putInBox ( const double& a ) {
    if ( ( x<=0.0 ) || ( x>a ) ) x-=utils::nint ( x/a ) *a;
    if ( ( y<=0.0 ) || ( y>a ) ) y-=utils::nint ( y/a ) *a;
    if ( ( z<=0.0 ) || ( z>a ) ) z-=utils::nint ( z/a ) *a;
    }

void atom::scaleVelocity ( const double& a, const double& b, const double& c ) {
    vx*=a;
    vy*=b;
    vz*=c;
    }

double atom::v2() {
    return ( vx*vx+vy*vy+vz*vz );
    }

double atom::getMass() {
    return mass;
    }

double atom::px() {
    return mass*vx;
    }

double atom::py() {
    return mass*vy;
    }
double atom::pz() {
    return mass*vz;
    }
void atom::addtoVelocities ( const double& p, const double& q, const double& r ) {
    vx+=p;
    vy+=q;
    vz+=r;
    }
double atom::force() const {
    return sqrt ( fx*fx+fy*fy+fz*fz );
    }


double r ( atom const& i, atom const& j ) {
    return sqrt ( ( j.x-i.x ) * ( j.x-i.x ) + ( j.y-i.y ) * ( j.y-i.y ) + ( j.z-i.z ) * ( j.z-i.z ) );
    }

double r ( atom const& i, atom const& j, double const &a ) {

    vector<double> R ( 3 );
    if ( a<0.0 ) {
        return r ( i,j );
        }

    dr ( R,i,j,a );
    return sqrt ( R[0]*R[0]+R[1]*R[1]+R[2]*R[2] );
    }



void dr ( vector<double>& x, atom const &i, atom const &j, double const& a ) {
    x[0]=j.x-i.x;
    x[0]=x[0]-utils::nint ( x[0]/a ) *a;
    x[1]=j.y-i.y;
    x[1]=x[1]-utils::nint ( x[1]/a ) *a;
    x[2]=j.z-i.z;
    x[2]=x[2]-utils::nint ( x[2]/a ) *a;
    }


double rdr ( vector<double>& g, atom const &i, atom const &j, double const& a ) {
    dr ( g,i,j,a );

    return sqrt ( g[0]*g[0]+g[1]*g[1]+g[2]*g[2] );
    }


void atom::dr ( vector<double>& x,  atom const &j, double const& a ) {
    x[0]=j.x-this->x;
    x[0]=x[0]-utils::nint ( x[0]/a ) *a;
    x[1]=j.y-this->y;
    x[1]=x[1]-utils::nint ( x[1]/a ) *a;
    x[2]=j.z-this->z;
    x[2]=x[2]-utils::nint ( x[2]/a ) *a;
    }

double atom::rdr ( vector<double>& g, atom const &j, double const& a ) {
    this->dr ( g,j,a );
    return sqrt ( g[0]*g[0]+g[1]*g[1]+g[2]*g[2] );
    }

double KineticEnergy ( vector< atom >& a ) {
    double kin=0.0;
    for ( vector<atom>::iterator i=a.begin(); i<a.end(); i++ ) {
        kin+= ( *i ).v2() * ( *i ).mass;
        }
    return kin/2.0;
    }

double VCM ( vector< atom >& a ) {
    double mt=0.0,a1=0.0,a2=0.0,a3=0.0;
    for ( vector<atom>::iterator i=a.begin(); i<a.end(); i++ ) {
        a1+= ( *i ).px();
        a2+= ( *i ).py();
        a3+= ( *i ).pz();
        mt+= ( *i ).mass;
        }
    a1=a1/mt;
    a2=a2/mt;
    a3=a3/mt;
    return sqrt ( a1*a1+a2*a2+a3*a3 );
    }

double PotentialEnergy ( vector<atom> const& a, const double& rc, const double& box, const bool& shift,double & vir ) {
    double R,energy;
    vector<double> dR ( 3 );
    vir =0.0;
    double ec=0.0;
    if ( shift ) ec=utils::LennardJones ( rc );
    energy=0.0;
    for ( int i=0; i<a.size()-1; i++ ) {
        for ( int j=i+1; j<a.size(); j++ ) {
            R=rdr ( dR,a[i],a[j],box );
            if ( R<rc ) {
                energy+=utils::LennardJones ( R )-ec;
                vir+=utils::LennardJonesdR ( R ) *R;
                }
            }
        }
    return energy;

    }

double atom::myPotentialEnergy ( vector< atom >& a , const double& rc, const double& box, const bool& shift,double& vir ) {

    double R,energy;
    vector<double> dR ( 3 );
    double ec=0.0;
    if ( shift ) ec=utils::LennardJones ( rc );
    energy=0.0;
    vir =0.0;
    for ( int it=0; it<id; it++ ) {
        R=rdr ( dR,a[it],box );
//         if ( R<utils::zero ) return utils::infinity;
        if ( R<rc ) {
            energy+=utils::LennardJones ( R )-ec;
            vir+=utils::LennardJonesdR ( R ) *R;
            }
        }

    for ( int it=id+1; it<a.size(); it++ ) {
        R=rdr ( dR,a[it],box );
//         if ( R<utils::zero ) return utils::infinity;
        if ( R<rc ) {
            energy+=utils::LennardJones ( R )-ec;
            vir+=utils::LennardJonesdR ( R ) *R;
            }
        }
    return energy;
    }


void atom::setNoise ( const vector< double >& b ) {
    xi.erase ( xi.begin(),xi.end() );
    eta.erase ( eta.begin(),eta.end() );
    xi.push_back ( b[0] );
    xi.push_back ( b[1] );
    xi.push_back ( b[2] );
    eta.push_back ( b[3] );
    eta.push_back ( b[4] );
    eta.push_back ( b[5] );
    }

void atom::updateVelocitiesVEC ( const double& h, const double&s ) {

    double ss=s*sqrt ( gamma/mass );
    double ax=0.5*h* ( 1.0-0.25*h*gamma );
    vx+= ( fx/mass-gamma*vx ) *ax+utils::AH ( h,xi[0],eta[0],gamma,ss );
    vy+= ( fy/mass-gamma*vy ) *ax+utils::AH ( h,xi[1],eta[1],gamma,ss );
    vz+= ( fz/mass-gamma*vz ) *ax+utils::AH ( h,xi[2],eta[2],gamma,ss );
    }

void atom::updatePositionsVEC ( const double& h, const double& s ) {
    double ss=s*sqrt ( gamma/mass );
    x+=h*vx+ss*eta[0];
    y+=h*vy+ss*eta[1];
    z+=h*vz+ss*eta[2];
    }


void atom::randomMove ( const double& dr, const double& a, const double& b, const double& c ) {
    x+= ( a-0.5 ) *dr;
    y+= ( b-0.5 ) *dr;
    z+= ( c-0.5 ) *dr;
    }

void atom::saveOldPositions ( double& xo, double& yo, double& zo ) {
    xo=x;
    yo=y;
    zo=z;
    }

double density ( vector< atom >& a, const double& box ) {
    double tm=0.0;
    for ( vector<atom>::iterator it=a.begin(); it<a.end(); it++ ) {
        tm+= ( *it ).mass;
        }
    return tm/ ( box*box*box );
    }


// kate: indent-mode cstyle; space-indent on; indent-width 4;
