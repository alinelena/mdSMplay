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

#include <cmath>
#include <iostream>

#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>

#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_01.hpp>

#include "solution.h"
#include "utils.h"
#include "reading.h"
#include "colours.h"
#include "gutenberg.h"
#include "collectivevariable.h"

using namespace std;
using namespace colours;

void solution::InitializeSimulation ( std::vector<atom>& a, const double& temp, const double& box, const int& seed, string const& filename , double const& gamma )
{
     InitializeSimulation ( a,temp,box,seed,filename );
     int k=0;
     for ( vector<atom>::iterator i=a.begin(); i<a.end(); i++ ) {
          ( *i ).setGamma ( gamma );
          ( *i ).setID ( k );
          k++;
     }
}

void solution::InitializeSimulation ( std::vector<atom>& a, const double& temp, const double& box, const int& seed, string const& filename )
{
// equal spaced in the box
     boost::mt19937 igen ( seed );
     boost::variate_generator<boost::mt19937, boost::normal_distribution<> >    rgauss ( igen, boost::normal_distribution<> ( 0.0,1.0 ) );
     double x,y,z,fact;

     int k=reading::readXYZ ( filename.c_str(),a );
     for ( int i=0; i<a.size(); i++ ) {
          a[i].setMass ( 1.0 );
          a[i].setID ( i );
          fact=sqrt ( utils::Kb*temp/a[i].getMass() );
          a[i].setVelocities ( fact*rgauss(),fact*rgauss(),fact*rgauss() );
          a[i].setForces ( 0.0,0.0,0.0 );
     }

//remove the movement of the centre of mass
     centerOfMassV ( a );
     scaleVelocities ( a,temp );
}

void solution::RandomVelocities ( std::vector<atom>& a, const double& temp, const double& box, const int& seed )
{
// equal spaced in the box
     boost::mt19937 igen ( seed );
     boost::variate_generator<boost::mt19937, boost::normal_distribution<> >    rgauss ( igen, boost::normal_distribution<> ( 0.0,1.0 ) );
     double x,y,z,fact;

     for ( int i=0; i<a.size(); i++ ) {
          fact=sqrt ( utils::Kb*temp/a[i].getMass() );
          a[i].setVelocities ( fact*rgauss(),fact*rgauss(),fact*rgauss() );
          a[i].setForces ( 0.0,0.0,0.0 );
     }
//remove the movement of the centre of mass
     centerOfMassV ( a );
     scaleVelocities ( a,temp );
}


void solution::InitializeSimulation ( std::vector<atom>& a, string const& filename )
{
     int k=reading::readFull ( filename.c_str(),a );
     if ( k!=0 ) {
          utils::Print ( "Error!!!",cout,colours::red );
     }
}


void solution::GenerateLattice ( vector< atom >&a, const double&rho,
                                 const int&mass, const double&nAtoms, double& box, const string& el )
{

     double volume=nAtoms*mass/rho;
     box=pow ( volume,1.0/3.0 );
     int n=int ( pow ( nAtoms,1.0/3.0 ) );
     if ( n*n*n!=nAtoms ) n++;
     double h=box/ ( n );
     double x;
     double y;
     double z;
     int m=0;
     for ( int i=0; i<n; i++ ) {
          x=i*h;
          for ( int j=0; j<n; j++ ) {
               y=j*h;
               for ( int k=0; k<n; k++ ) {
                    z=k*h;
                    if ( m<nAtoms ) {
                         a.push_back ( atom() );
                         a[m].setElement ( el );
                         a[m].setPositions ( x,y,z );
                         a[m].setMass ( mass );
                         a[m].setID ( m );
                    }
                    m++;
               }
          }
     }

}



void solution::centerOfMassV ( vector< atom >& a )
{
     double vmx=0.0, vmy=0.0, vmz=0.0,mt=0.0;
     for ( vector<atom>::iterator it=a.begin(); it<a.end(); it++ ) {
          vmx+= ( *it ).px();
          vmy+= ( *it ).py();
          vmz+= ( *it ).pz();
          mt+= ( *it ).getMass();
     }
     vmx=-vmx/mt;
     vmy=-vmy/mt;
     vmz=-vmz/mt;
     for ( vector<atom>::iterator it=a.begin(); it<a.end(); it++ ) {
          ( *it ).addtoVelocities ( vmx,vmy,vmz );
     }
}

void solution::scaleVelocities ( vector< atom >& a, const double& T )
{
     double tx=0.0,ty=0.0, tz=0.0;
     if ( T>0.0 ) {
          for ( vector<atom>::iterator it=a.begin(); it<a.end(); it++ ) {
               tx+= ( *it ).getMass() * ( *it ).vx2();
               ty+= ( *it ).getMass() * ( *it ).vy2();
               tz+= ( *it ).getMass() * ( *it ).vz2();
          }
          tx=tx/ ( a.size() *utils::Kb ) ;
          ty=ty/ ( a.size() *utils::Kb ) ;
          tz=tz/ ( a.size() *utils::Kb ) ;
          tx=sqrt ( T/tx );
          ty=sqrt ( T/ty );
          tz=sqrt ( T/tz );
     }
     for ( vector<atom>::iterator it=a.begin(); it<a.end(); it++ ) {
          ( *it ).scaleVelocity ( tx,ty,tz );
     }
}

double solution::computeForces ( vector< atom >& a, const double & rc, const double& box, double& vir )
{
     for ( vector<atom>::iterator it=a.begin(); it<a.end(); it++ ) {
          ( *it ).setForces ( 0.0,0.0,0.0 );
     }

     double R,fR,fx,fy,fz,energy;
     vector<double> dR ( 3 );
     double ec=utils::LennardJones ( rc );
     energy=0.0;
     for ( vector<atom>::iterator i=a.begin(); i<a.end()-1; i++ ) {
          for ( vector<atom>::iterator j=i+1; j<a.end(); j++ ) {
               R=rdr ( dR,*i,*j,box );
               if ( R<rc ) {
                    fR=utils::LennardJonesdR ( R );
                    fx=-fR*dR[0]/R;
                    fy=-fR*dR[1]/R;
                    fz=-fR*dR[2]/R;
                    ( *i ).addtoForces ( fx,fy,fz );
                    ( *j ).addtoForces ( -fx,-fy,-fz );
                    energy+=utils::LennardJones ( R )-ec;
                    vir+=fR*R;
               }
          }
     }
     return energy;
}

double solution::velocityVerlet ( vector< atom >& a, const double& rc, const double& box, const double& dt , double& kin, double &vir )
{
     kin=0.0;
     vir=0.0;
     for ( vector<atom>::iterator i=a.begin(); i<a.end(); i++ ) {
          ( *i ).updatePositionsVV ( dt );
          ( *i ).putInBox ( box );
     }
     for ( vector<atom>::iterator i=a.begin(); i<a.end(); i++ ) {
          ( *i ).updateVelocitiesVV ( dt );
     }
     double energy = solution::computeForces ( a,rc,box,vir );
     for ( vector<atom>::iterator i=a.begin(); i<a.end(); i++ ) {
          ( *i ).updateVelocitiesVV ( dt );
          kin+= ( *i ).v2() * ( *i ).getMass();
     }
     kin=kin/2.0;
     return energy;
}


void solution::MDDriverVV ( vector< atom >& a, generalInputs& in )
{

     double potEnergy, kinEnergy, vir,pressure;
     potEnergy = PotentialEnergy ( a,in.cut,in.box, in.shift,vir );
     kinEnergy = KineticEnergy ( a );
     double iT,tT;
     double consQ;
     double rho=density ( a, in.box );
     double volume = in.box*in.box*in.box;
     pressure=rho* ( utils::Kb*in.TBath ) +vir/volume/3.0;
     consQ=potEnergy+kinEnergy;
     iT=2.0/3.0*kinEnergy/a.size() /utils::Kb;
     tT=iT;
     gutenberg::mdHeader ( cout,in.fancy );
     gutenberg::mdRepLine ( 0.0,potEnergy/a.size(),kinEnergy/a.size(),iT,tT,pressure,consQ/a.size(),VCM ( a ),in.fancy,cout );
     gutenberg::printXYZ ( in.xyz,in.name,a,true );
     gutenberg::printAtoms ( a,in.debug );
     vector<double> avg ( 2 );
     int nsamp=0;
     for ( int k=0; k<avg.size(); k++ ) avg[k]=0.0;
     for ( int i=1; i<=in.nSteps; i++ ) {
          potEnergy = solution::velocityVerlet ( a,in.cut,in.box,in.dt,kinEnergy,vir );

          iT=2.0/3.0*kinEnergy/a.size() /utils::Kb;
          tT+= ( iT-tT ) / ( i+1.0 );

          if ( i%in.frequency==0 ) {
               gutenberg::printForces ( a,in.log );
               gutenberg::printXYZ ( in.xyz,in.name,a,true );
               gutenberg::printAtoms ( a,in.debug );
          }
          consQ=potEnergy+kinEnergy;
          avg[0]+=consQ;
          avg[1]+=consQ*consQ;
          nsamp++;
          pressure=rho* ( utils::Kb*in.TBath ) +vir/volume/3.0;
          gutenberg::mdRepLine ( i*in.dt,potEnergy/a.size(),kinEnergy/a.size(),iT,tT,pressure,consQ/a.size(),VCM ( a ),in.fancy,cout );
     }
     for ( int k=0; k<avg.size(); k++ ) avg[k]=avg[k]/nsamp;
     utils::Print ( "CV:",cout,red );
     utils::Print ( ( avg[1]-avg[0]*avg[0] ) /utils::Kb/tT/tT,cout, green );
     cout<<endl;

}

double solution::VECIntegrator ( vector< atom >& a, const double&rc , const double& box , const double& dt, const double& s1, const double& s2, double& kin, double& vir )
{

     kin=0.0;
     vir = 0.0;
     for ( vector<atom>::iterator i=a.begin(); i<a.end(); i++ ) {
          ( *i ).updateVelocitiesVEC ( dt,s2 );
     }

     for ( vector<atom>::iterator i=a.begin(); i<a.end(); i++ ) {
          ( *i ).updatePositionsVEC ( dt,s1 );
     }
     double energy = solution::computeForces ( a,rc,box, vir );
     for ( vector<atom>::iterator i=a.begin(); i<a.end(); i++ ) {
          ( *i ).updateVelocitiesVEC ( dt,s2 );
          kin+= ( *i ).v2() * ( *i ).getMass();
     }
     kin=kin/2.0;
     return energy;
}

void solution::MDDriverVEC ( vector< atom >& a, generalInputs & in )
{
     //set the random generators...
     boost::mt19937 igen ( in.seeds[0] );
     boost::variate_generator<boost::mt19937, boost::normal_distribution<> >  xi ( igen, boost::normal_distribution<> ( 0.0,1.0 ) );
     vector<double> b ( 6 );
     double potEnergy, kinEnergy, vir;


     potEnergy = PotentialEnergy ( a,in.cut,in.box,in.shift,vir );
     kinEnergy = KineticEnergy ( a );

     double iT,tT;
     double consQ;
     double rho=density ( a, in.box );
     double volume = in.box*in.box*in.box;
     double pressure=rho* ( utils::Kb*in.TBath ) +vir/volume/3.0;
     gutenberg::logLine ( "Volume: ",volume,in.log );
     gutenberg::logLine ( "Density: ",rho,in.log );
     consQ=potEnergy+kinEnergy;
     iT=2.0/3.0*kinEnergy/a.size() /utils::Kb;
     tT=iT;
     gutenberg::mdHeader ( cout,in.fancy );
     gutenberg::mdRepLine ( 0.0,potEnergy/a.size(),kinEnergy/a.size(),iT,tT,pressure,consQ/a.size(),VCM ( a ),in.fancy,in.log );
     gutenberg::printXYZ ( in.xyz,in.name,a,true );
     gutenberg::printAtoms ( a,in.debug );

     double s1=sqrt ( utils::Kb*in.TBath*in.dt/6.0 ) *in.dt;
     double s2=sqrt ( 2.0*utils::Kb*in.TBath );
     double v2=0.0;

     vector<double> K,E,vx,vy,vz,P,v;
//     K.push_back ( kinEnergy/a.size() );
//     E.push_back ( consQ/a.size() );
//     P.push_back ( potEnergy/a.size() );
//     for ( vector<atom>::iterator it=a.begin(); it<a.end(); it++ ) {
//         vx.push_back ( ( *it ).gvx() );
//         vy.push_back ( ( *it ).gvy() );
//         vz.push_back ( ( *it ).gvz() );
//         v.push_back ( sqrt ( ( *it ).v2() ) );
//         }
     double aene=0.0, aene2=0.0,cv;
     int nsamp=0;
     for ( int i=1; i<=in.nSteps; i++ ) {
          for ( vector<atom>::iterator it=a.begin(); it<a.end(); it++ ) {
               b[0]=xi();
               b[1]=xi();
               b[2]=xi();
               b[3]=xi();
               b[4]=xi();
               b[5]=xi();
               ( *it ).setNoise ( b );
          }
          potEnergy = solution::VECIntegrator ( a,in.cut,in.box,in.dt,s1,s2,kinEnergy,vir );



//         pressure=rho* ( utils::Kb*tT ) +vir/volume/3.0;
          gutenberg::mdRepLine ( i*in.dt,potEnergy/a.size(),kinEnergy/a.size(),iT,tT,pressure,consQ/a.size(),VCM ( a ),in.fancy,in.log );

          if ( i%in.frequency==0 )  {
//             gutenberg::printForces ( a,in.log );
//             gutenberg::printXYZ ( in.xyz,in.name,a,true );
//             gutenberg::printAtoms ( a,in.debug );
               iT=2.0*kinEnergy/ ( 3.0*a.size() *utils::Kb );

               consQ=potEnergy+kinEnergy;
               if ( i == in.nEquil ) {
                    nsamp=0;
                    E.erase ( E.begin(),E.end() );
               }
               if ( nsamp==0 ) {

                    aene=consQ;
                    tT=iT;
                    nsamp=0;
               } else {
                    nsamp++;
                    aene+= ( aene-consQ ) / ( nsamp+1.0 );
                    tT+= ( iT-tT ) / ( nsamp+1.0 );
               }
               E.push_back ( consQ/a.size() );

               cv= utils::sigma2 ( E ) / ( utils::Kb*tT*tT );
               if ( i>in.nEquil ) {
                    gutenberg::logInfo ( "step ",i-in.nEquil,cout );

               } else {
                    gutenberg::logInfo ( "step ",i,cout );
               }
               gutenberg::logInfo ( "T ",tT,cout );
               gutenberg::logInfo ( "<E> ",aene,cout );
               gutenberg::logLine ( "Cv: ", cv*a.size(),cout );
               gutenberg::logLine ( "Cv/N: ", cv,cout );

               if ( i>in.nEquil ) {

                    P.push_back ( potEnergy/a.size() );
                    for ( vector<atom>::iterator it=a.begin(); it<a.end(); it++ ) {
                         vx.push_back ( ( *it ).gvx() );
                         vy.push_back ( ( *it ).gvy() );
                         vz.push_back ( ( *it ).gvz() );
                         v2= ( *it ).v2();
                         v.push_back ( sqrt ( v2 ) );
                         K.push_back ( v2* ( *it ).getMass() /2.0 );
                    }
               }
          }
     }

     cv=utils::sigma2 ( E );
     utils::Print ( "CV:",cout,red );
     utils::Print ( cv/ ( utils::Kb*tT*tT ) *a.size(),cout, green );
     cout<<endl;
     utils::Print ( "sigma^2 K",cout, red );
     utils::Print ( utils::sigma2 ( K ),cout,green );
     cout<<endl;
     utils::Print ( "sigma^2 E",cout, red );
     utils::Print ( cv,cout,green );
     cout<<endl;
     solution::histogramWrapper ( K,in.bin,"Kvec.hist" );
     solution::histogramWrapper ( E,in.bin,"Evec.hist" );
     solution::histogramWrapper ( P,in.bin,"Pvec.hist" );
     solution::histogramWrapper ( vx,in.bin,"vxvec.hist" );
     solution::histogramWrapper ( vy,in.bin,"vyvec.hist" );
     solution::histogramWrapper ( vz,in.bin,"vzvec.hist" );
     solution::histogramWrapper ( v,in.bin,"vvec.hist" );
}

double solution::computeForcesTAMD ( vector< atom >& a, vector<colVar>& cv,const double & rc, const double& box, double& vir )
{
     for ( vector<atom>::iterator it=a.begin(); it<a.end(); it++ ) {
          ( *it ).setForces ( 0.0,0.0,0.0 );
     }

     double R,fR,fx,fy,fz,energy;
     vector<double> dR ( 3 );
     double ec=utils::LennardJones ( rc );
     int k=0;
     energy=0.0;
     for ( vector<colVar>::iterator i=cv.begin(); i<cv.end(); i++ ) {
          ( *i ).force ( a,box );
     }
     for ( vector<atom>::iterator i=a.begin(); i<a.end(); i++ ) {
          fx=0.0;
          fy=0.0;
          fz=0.0;
          for ( vector<colVar>::iterator it=cv.begin(); it<cv.end(); it++ ) {
               ( *it ).dtheta ( dR,k,a,box );
               fx=fx+ ( *it ).f*dR[0];
               fy=fy+ ( *it ).f*dR[1];
               fz=fz+ ( *it ).f*dR[2];
          }
          k++;
          ( *i ).addtoForces ( -fx,-fy,-fz );
          for ( vector<atom>::iterator j=i+1; j<a.end(); j++ ) {
               R=rdr ( dR,*i,*j,box );
               if ( R<rc ) {
                    fR=utils::LennardJonesdR ( R );
                    fx=-fR*dR[0]/R;
                    fy=-fR*dR[1]/R;
                    fz=-fR*dR[2]/R;
                    ( *i ).addtoForces ( fx,fy,fz );
                    ( *j ).addtoForces ( -fx,-fy,-fz );
                    energy+=utils::LennardJones ( R )-ec;
                    vir+=fR*R;
               }
          }
     }
     return energy;
}

double solution::TAMDIntegrator ( vector< atom >& a, vector<colVar>&cv, const double&rc , const double& box , const double& dt, const double& s1,const double& c1, const double& s2, double& kin, double&kinz, double&vir )
{

     kin=0.0;
     kinz=0.0;
     vir =0.0;

     for ( vector<atom>::iterator i=a.begin(); i<a.end(); i++ ) {
          ( *i ).updateVelocitiesVEC ( dt,s2 );
     }
     for ( vector<colVar>::iterator i=cv.begin(); i<cv.end(); i++ ) {
          ( *i ).updatezvVEC ( dt );
     }

     for ( vector<atom>::iterator i=a.begin(); i<a.end(); i++ ) {
          ( *i ).updatePositionsVEC ( dt,s1 );
     }
     for ( vector<colVar>::iterator i=cv.begin(); i<cv.end(); i++ ) {
          ( *i ).updatezpVEC ( dt,c1 );
     }
     double energy = solution::computeForcesTAMD ( a,cv,rc,box, vir );
     for ( vector<atom>::iterator i=a.begin(); i<a.end(); i++ ) {
          ( *i ).updateVelocitiesVEC ( dt,s2 );
          kin+= ( *i ).v2() * ( *i ).getMass();
     }

     for ( vector<colVar>::iterator i=cv.begin(); i<cv.end(); i++ ) {
          ( *i ).updatezvVEC ( dt );
          kinz+= ( *i ).v* ( *i ).v* ( *i ).mu;
     }

     kin=kin/2,0;
     kinz=kinz/2.0;
     return energy;
}

///
///  @brief the driver for TAMD
///  @details http://dx.doi.org/10.1016/j.cplett.2006.05.062
///   it uses the VEC integrator for EOM
///  @param a a vector of atoms
///  @param cv the vector of collective variables
///  @param in contains the input parameters
///
void solution::MDDriverTAMD ( vector< atom >& a, vector<colVar>& cv, generalInputs & in )
{
     //set the random generators...

     boost::mt19937 igen ( in.seeds[0] );
     boost::variate_generator<boost::mt19937, boost::normal_distribution<> >  xi ( igen, boost::normal_distribution<> ( 0.0,1.0 ) );
     vector<double> b ( 6 );
     double potEnergy, kinEnergy, vir;

     potEnergy = PotentialEnergy ( a,in.cut,in.box,in.shift,vir );
     kinEnergy = KineticEnergy ( a );

     double iT,tT;
     double iTz,tTz;
     double consQ;
     double rho=density ( a, in.box );
     double volume = in.box*in.box*in.box;
     double pressure=rho* ( utils::Kb*in.TBath ) +vir/volume/3.0;

     consQ=potEnergy+kinEnergy;
     iT=2.0*kinEnergy/ ( 3.0*a.size() *utils::Kb );
     tT=iT;
     iTz=0.0;
     tTz=iTz;
     double kinz=0.0;
     gutenberg::mdTAMDHeader ( cv,cout,in.fancy );
     gutenberg::mdTAMDRepLine ( 0.0,potEnergy/a.size(),kinEnergy/a.size(),tT,tTz,pressure,kinz/cv.size(),consQ/a.size(),VCM ( a ),cv,in.fancy,cout );
     gutenberg::printXYZ ( in.xyz,in.name,a,true );
     gutenberg::printAtoms ( a,in.debug );

     double s1=sqrt ( utils::Kb*in.TBath*in.dt/6.0 ) *in.dt;
     double c1=sqrt ( utils::Kb*in.Tcv*in.dt/6.0 ) *in.dt;
     double s2=sqrt ( 2.0*utils::Kb*in.TBath );

     vector<double> KX,KZ,vxX,vyX,vyZ,vX,vzX,vZ;
     double v2;
     for ( vector<colVar>::iterator it=cv.begin(); it<cv.end(); it++ ) {
          ( *it ).s=sqrt ( 2.0*utils::Kb*in.Tcv* ( *it ).gamma/ ( *it ).mu );
     }
     for ( int i=1; i<=in.nSteps; i++ ) {
          for ( vector<atom>::iterator it=a.begin(); it<a.end(); it++ ) {
               b[0]=xi();
               b[1]=xi();
               b[2]=xi();
               b[3]=xi();
               b[4]=xi();
               b[5]=xi();
               ( *it ).setNoise ( b );
          }
          for ( vector<colVar>::iterator it=cv.begin(); it<cv.end(); it++ ) {
               b[0]=xi();
               b[1]=xi();
               ( *it ).setNoise ( b );
          }
          potEnergy = solution::TAMDIntegrator ( a,cv,in.cut,in.box,in.dt,s1,c1,s2,kinEnergy, kinz,vir );

          iT=2.0*kinEnergy/ ( 3.0*a.size() *utils::Kb );
          tT+= ( iT-tT ) / ( i+1.0 );

          iTz=2.0*kinz/ ( cv.size() *utils::Kb );
          tTz+= ( iTz-tTz ) / ( i+1.0 );
          if ( i%in.frequency==0 && i>in.nEquil ) {
               gutenberg::printForces ( a,in.log );
               gutenberg::printXYZ ( in.xyz,in.name,a,true );
               gutenberg::printAtoms ( a,in.debug );
               for ( vector<atom>::iterator it=a.begin(); it<a.end(); it++ ) {
                    vxX.push_back ( ( *it ).gvx() );
                    vyX.push_back ( ( *it ).gvy() );
                    vzX.push_back ( ( *it ).gvz() );
                    v2= ( *it ).v2();
                    vX.push_back ( sqrt ( v2 ) );
                    KX.push_back ( v2* ( *it ).getMass() /2.0 );
               }
               for ( vector<colVar>::iterator it=cv.begin(); it<cv.end(); it++ ) {
                    vZ.push_back ( ( *it ).v );
                    KZ.push_back ( ( *it ).mu* ( *it ).v* ( *it ).v/2.0 );
               }
          }
          consQ=potEnergy+kinEnergy;
          pressure=rho* ( utils::Kb*in.TBath ) +vir/volume/3.0;
          gutenberg::mdTAMDRepLine ( i*in.dt,potEnergy/a.size(),kinEnergy/a.size(),tT,tTz,pressure,kinz/cv.size(),consQ/a.size(),VCM ( a ),cv,in.fancy,cout );
     }
     solution::histogramWrapper ( vxX,in.bin,"vxX.hist" );
     solution::histogramWrapper ( vyX,in.bin,"vyX.hist" );
     solution::histogramWrapper ( vzX,in.bin,"vzX.hist" );
     solution::histogramWrapper ( vX,in.bin,"vX.hist" );
     solution::histogramWrapper ( KX,in.bin,"KX.hist" );
     solution::histogramWrapper ( vZ,in.bin,"vZ.hist" );
     solution::histogramWrapper ( KZ,in.bin,"KZ.hist" );
}


void solution::optimalDr ( double& d , int& otry, int& ogood, const int& itry, const int& igood, const double& box )
{
     if ( itry!=otry ) {
          double f= ( double ) ( igood-ogood ) / ( double ) ( itry-otry );
          double od=d;
          d=d*f/0.5; //0.5 is the target accteptance ratio;
          if ( d/od>1.5 ) d=od*1.5;
          if ( d/od<0.5 ) d=od*0.5;
          if ( d>box/2 ) d=box/2;
          otry=itry;
          ogood=igood;
     }

}


///
///  @brief the driver for a metropolis montecarlo simulation
///  @param a a vector of atoms
///  @param in contains the input parameters
///
void solution::MMCDriver ( vector< atom >& a, generalInputs& in )
{
     boost::mt19937 igen ( in.seeds[0] );
     boost::uniform_01<boost::mt19937> rn ( igen );

     double ei,ef;
     int rAtom;
     double rx,ry,rz;

     int itry=0;
     int igood=0;
     int otry, ogood;
     otry=itry;
     ogood=igood;
     int nsamp=0;
     double fact;
     double et,vt;
     double aene=0.0, ap=0.0,aene2=0.0;
     double rho=density ( a,in.box );
     double volume=in.box*in.box*in.box;
     double vir, viri,virf;
     double etot=PotentialEnergy ( a,in.cut,in.box,in.shift,vir );
     gutenberg::logLine ( "Volume: ",volume,in.log );
     gutenberg::logLine ( "Density: ",rho,in.log );
     utils::Print ( etot,cout,32,16 );
     cout<<endl;

     vector<double>ene;

     double pressure;
     int cycles;
     int startsamp;
     pressure=rho* ( utils::Kb*in.TBath ) +vir/volume/3.0;
     gutenberg::MMCHeader ( cout, in.fancy );
     gutenberg::MMCRepLine ( 0,etot/a.size(),pressure,in.fancy,cout );
     gutenberg::logLine ( "Initial etot ", etot,cout );
     gutenberg::logLine ( "Initial etot/N ", etot/a.size(),cout );
     gutenberg::logLine ( "Initial virial ", vir,cout );
     solution::optimalDr ( in.dr,otry,ogood,itry,igood,in.box );
     for ( int stage=0; stage<2; stage++ ) {
          if ( stage==0 ) {
               cycles=in.nEquil;
               startsamp=cycles/10;
               utils::Print ( "Equilibration stage:", cout, blue );
               cout<<endl;
          } else {
               startsamp=0;
               cycles=in.mccycles;
               utils::Print ( "Production stage:", cout, blue );
               cout <<endl;
               etot=PotentialEnergy ( a,in.cut,in.box,in.shift,vir );
               ene.push_back ( etot/a.size() );
               gutenberg::MMCHeader ( cout, in.fancy );
               gutenberg::MMCRepLine ( 0,etot/a.size(),pressure,in.fancy,cout );
          }
          gutenberg::logInfo ( "old dr:",in.dr,cout );
          solution::optimalDr ( in.dr,otry,ogood,itry,igood,in.box );
          gutenberg::logLine ( "new dr:",in.dr,cout );
          itry=0;
          igood=0;
          aene=0.0;
          aene2=0.0;
          ap=0.0;
          nsamp=0;
          for ( int i=0; i<cycles; i++ ) {
               // one should compute the optimal dr here...
               for ( int j=0; j<in.mcsteps; j++ ) {
                    itry++;
                    rAtom=int ( rn() *a.size() );
                    ei=a[rAtom].myPotentialEnergy ( a,in.cut,in.box,in.shift,viri );
                    a[rAtom].saveOldPositions ( rx,ry,rz );
                    a[rAtom].randomMove ( in.dr,rn(),rn(),rn() );
                    ef=a[rAtom].myPotentialEnergy ( a,in.cut,in.box,in.shift,virf );
                    fact=- ( ef-ei ) / ( utils::Kb*in.TBath );
                    if ( rn() <exp ( fact ) ) {
                         igood++;
                         etot=etot+ef-ei;
                         vir=vir+virf-viri;
                         a[rAtom].putInBox ( in.box );
                    } else { //if rejected
                         a[rAtom].setPositions ( rx,ry,rz );
                    }
               }
               if ( ( ( i+1 ) %in.frequency==0 ) && ( i+1>startsamp ) ) {
                    //sampling point...
                    pressure=rho* ( utils::Kb*in.TBath ) +vir/volume/3.0;
                    gutenberg::MMCRepLine ( ( i+1 ) *in.mcsteps,etot/a.size(),pressure,in.fancy,in.log );
//                 gutenberg::printXYZ ( in.xyz,in.name,a,true );
                    ene.push_back ( etot/a.size() );
                    aene+=etot;
                    aene2+=etot*etot;
                    ap+=pressure;
                    nsamp++;

                    gutenberg::logInfo ( "|cycle: ",i+1,in.log );
                    gutenberg::logInfo ( "|attempts: ",itry, in.log );
                    gutenberg::logInfo ( "|accepted: ",igood, in.log );
                    gutenberg::logLine ( "|ratio", ( double ) ( igood ) / ( double ) ( itry ),in.log );
               }
               if ( ( i+1 ) % ( cycles/5 ) ==0 ) {
                    gutenberg::logInfo ( "|cycle: ",i+1,cout );
                    gutenberg::logInfo ( "|attempts: ",itry-otry, cout );
                    gutenberg::logInfo ( "|accepted: ",igood-ogood, cout );
                    gutenberg::logLine ( "|ratio", ( double ) ( igood-ogood ) / ( double ) ( itry-otry ),cout );
                    gutenberg::logInfo ( "old dr:",in.dr,cout );
                    solution::optimalDr ( in.dr,otry,ogood,itry,igood,in.box );
                    gutenberg::logInfo ( "new dr:",in.dr,cout );
                    cout<<endl;
               }

          }
          gutenberg::logLine ( "Etot: ",etot,cout );
          gutenberg::logLine ( "Erot/N",etot/a.size(),cout );
          pressure=rho* ( utils::Kb*in.TBath ) +vir/volume/3.0;
          gutenberg::logLine ( "virial: ",vir,cout );
          gutenberg::logLine ( "pressure: ",pressure,cout );
          et=PotentialEnergy ( a,in.cut,in.box,in.shift,vt );
          gutenberg::logLine ( "diff energy: ",et-etot,cout );
          gutenberg::logLine ( "diff virial: ", vt-vir,cout );
          gutenberg::logInfo ( "|cycle: ",cycles,cout );
          gutenberg::logInfo ( "|attempts: ",itry, cout );
          gutenberg::logInfo ( "|accepted: ",igood, cout );
          gutenberg::logLine ( "|ratio", ( double ) ( igood ) / ( double ) ( itry ),cout );
          aene2=aene2/nsamp;
          aene=aene/nsamp;
          gutenberg::logLine ( "Avg energy:", aene,cout );
          gutenberg::logLine ( "Avg energy/N:", aene/a.size(),cout );
          gutenberg::logLine ( "Avg pressure:", ap/nsamp,cout );
          gutenberg::logLine ( "CV: ",1.5*utils::Kb*a.size() + ( aene2-aene*aene ) / ( utils::Kb*in.TBath*in.TBath ),cout );
     }

     solution::histogramWrapper ( ene,in.bin,"Pmmc.hist" );

}



///
///  @brief the driver for TAMC
///  @details TAMC paper
///   it uses the VEC integrator for EOM
///  @param a a vector of atoms
///  @param cv the vector of collective variables
///  @param in contains the input parameters
///
void solution::MDDriverTAMC ( vector< atom >& a, vector<colVar>& cv, generalInputs & in )
{
     //set the random generators...

     boost::mt19937 igen ( in.seeds[0] );
     boost::variate_generator<boost::mt19937, boost::normal_distribution<> >  xi ( igen, boost::normal_distribution<> ( 0.0,1.0 ) );
     boost::uniform_01<boost::mt19937> rn ( igen );
     vector<double> b ( 2 );
     double potEnergy, vir,fact,etot;

     potEnergy = PotentialEnergy ( a,in.cut,in.box,in.shift,vir ) +zpot ( a,cv,in.box );

     double iT=0.0,tT=0.0;
     double consQ,kin;
     double rho=density ( a, in.box );
     double volume = in.box*in.box*in.box;
     double pressure=rho* ( utils::Kb*in.TBath ) +vir/volume/3.0;

     int itry=0, otry=0,igood=0, ogood=0;
     double ei,ef, viri,virf;
     int rAtom;
     double rx,ry,rz;
     double kinz=0.0;

     gutenberg::mdTAMCHeader ( cv,cout,in.fancy );
     gutenberg::mdTAMCRepLine ( 0.0,potEnergy/a.size(),tT,pressure,kinz,cv,in.fancy,cout );
     gutenberg::printXYZ ( in.xyz,in.name,a,true );

     double c1=sqrt ( utils::Kb*in.Tcv*in.dt/6.0 ) *in.dt;


     double aene=0.0, ap=0.0;
     int nsamp=0;
     vector<double> avgfz ( cv.size() ), KZ,vZ;
     for ( vector<colVar>::iterator it=cv.begin(); it<cv.end(); it++ ) {
          ( *it ).s=sqrt ( 2.0*utils::Kb*in.Tcv* ( *it ).gamma/ ( *it ).mu );
     }
     for ( int step=1; step<=in.nSteps; step++ ) {
          for ( vector<colVar>::iterator it=cv.begin(); it<cv.end(); it++ ) {
               b[0]=xi();
               b[1]=xi();
               ( *it ).setNoise ( b );
          }
          vir =0.0;
          aene=0.0;
          ap=0.0;
          for ( int k=0; k<cv.size(); k++ ) avgfz[k]=0.0;
          for ( vector<colVar>::iterator i=cv.begin(); i<cv.end(); i++ ) {
               ( *i ).updatezvVEC ( in.dt );
          }
          for ( vector<colVar>::iterator i=cv.begin(); i<cv.end(); i++ ) {
               ( *i ).updatezpVEC ( in.dt,c1 );
          }
          //steps MC

          for ( int i=0; i<in.mccycles; i++ ) {
               for ( int j=0; j<in.mcsteps; j++ ) {
                    itry++;
                    rAtom=int ( rn() *a.size() );
                    ei=a[rAtom].myPotentialEnergy ( a,in.cut,in.box,in.shift,viri ) +zpot ( a,cv,in.box );
                    a[rAtom].saveOldPositions ( rx,ry,rz );
                    a[rAtom].randomMove ( in.dr,rn(),rn(),rn() );
                    a[rAtom].putInBox ( in.box );
                    ef=a[rAtom].myPotentialEnergy ( a,in.cut,in.box,in.shift,virf ) +zpot ( a,cv,in.box );
                    fact=- ( ef-ei ) / ( utils::Kb*in.TBath );
                    if ( rn() <exp ( fact ) ) {
                         igood++;
                         etot=etot+ef-ei;
                         vir=vir+virf-viri;
                    } else { //if rejected
                         a[rAtom].setPositions ( rx,ry,rz );
                    }
               }
               //sampling point...

               pressure=rho* ( utils::Kb*in.TBath ) +vir/volume/3.0;
               aene+=etot/a.size();
               ap+=pressure;
               nsamp++;
               for ( int k=0; k<cv.size(); k++ ) {
                    cv[k].force ( a,in.box );
                    avgfz[k]+=cv[k].f;
               }
               solution::optimalDr ( in.dr,otry,ogood,itry,igood,in.box );
          }
          for ( int k=0; k<cv.size(); k++ ) {
               cv[k].f=avgfz[k]/nsamp;
          }
          for ( vector<colVar>::iterator i=cv.begin(); i<cv.end(); i++ ) {
               ( *i ).updatezvVEC ( in.dt );
          }

          if ( step%in.frequency==0 ) {
               gutenberg::printXYZ ( in.xyz,in.name,a,true );
          }
          kinz=zKin ( cv );

          iT=2.0*kinz/ ( cv.size() *utils::Kb );
          if ( step>=in.nEquil ) {
               tT+= ( iT-tT ) / ( step-in.nEquil+1.0 );
               if ( step%in.frequency==0 ) {
                    for ( vector<colVar>::iterator it=cv.begin(); it<cv.end(); it++ ) {
                         vZ.push_back ( ( *it ).v );
                         KZ.push_back ( ( *it ).mu* ( *it ).v* ( *it ).v/2.0 );
                    }
               }
          }
          gutenberg::mdTAMCRepLine ( step*in.dt,aene/nsamp,tT,pressure/nsamp,kinz,cv,in.fancy,cout );
          nsamp=0;
     }

     solution::histogramWrapper ( vZ,in.bin,"vZtdmc.hist" );
     solution::histogramWrapper ( KZ,in.bin,"KZtdmc.hist" );

}


///
///  @brief the driver for TAHMC
///  @details TAMC paper and HMC paper
///   it uses the VEC integrator for EOM and velocity verlet for atoms
///  @param a a vector of atoms
///  @param cv the vector of collective variables
///  @param in contains the input parameters
///
void solution::MDDriverTAHMC ( vector< atom >& a, vector<colVar>& cv, generalInputs & in )
{
     //set the random generators...

     boost::mt19937 igen ( in.seeds[0] );
     boost::variate_generator<boost::mt19937, boost::normal_distribution<> >  xi ( igen, boost::normal_distribution<> ( 0.0,1.0 ) );
     boost::uniform_01<boost::mt19937> rn ( igen );
     vector<double> b ( 2 );
     double aPot, vir,fact,etot;
     double etoti,etotf;

     aPot = PotentialEnergy ( a,in.cut,in.box,in.shift,vir ) + zpot ( a,cv,in.box );

     double iT=0.0,tT=0.0;
     double consQ,kin;
     double rho=density ( a, in.box );
     double volume = in.box*in.box*in.box;
     double pressure=rho* ( utils::Kb*in.TBath ) +vir/volume/3.0;

     int itry=0, otry=0,igood=0, ogood=0;
     double zpoti,zpotf, viri,virf;
     int irseed;

     double kinz=0.0;
     double kini,kinf;
     vector<double> pos ( 3* a.size() );
     cout << pos.size() <<" "<<a.size() <<endl;
     for ( int i=0; i<a.size(); i++ ) {
          a[i].saveOldPositions ( pos[3*i],pos[3*i+1],pos[3*i+2] );
     }
     gutenberg::mdTAMCHeader ( cv,cout,in.fancy );
     gutenberg::mdTAMCRepLine ( 0.0,aPot/a.size(),tT,pressure,kinz,cv,in.fancy,cout );
     gutenberg::printXYZ ( in.xyz,in.name,a,true );

     double c1=sqrt ( utils::Kb*in.Tcv*in.dtz/6.0 ) *in.dt;


     double aene=0.0, ap=0.0,acc;
     int nsamp;
     vector<double> avgfz ( cv.size() ), KZ,vZ;
     irseed=rn();
     solution::zHMCSampler ( a, cv, in, irseed,pos, avgfz );
     for ( vector<colVar>::iterator it=cv.begin(); it<cv.end(); it++ ) {
          ( *it ).s=sqrt ( 2.0*utils::Kb*in.Tcv* ( *it ).gamma/ ( *it ).mu );
     }
     for ( int step=1; step<=in.Nz; step++ ) {
          for ( vector<colVar>::iterator it=cv.begin(); it<cv.end(); it++ ) {
               b[0]=xi();
               b[1]=xi();
               ( *it ).setNoise ( b );
          }
          for ( vector<colVar>::iterator i=cv.begin(); i<cv.end(); i++ ) {
               ( *i ).updatezvVEC ( in.dtz );
          }
          for ( vector<colVar>::iterator i=cv.begin(); i<cv.end(); i++ ) {
               ( *i ).updatezpVEC ( in.dtz,c1 );
          }
          //steps HMC
          irseed=rn();
          solution::zHMCSampler ( a, cv, in, irseed,pos, avgfz );
//         solution::zSample


          for ( vector<colVar>::iterator i=cv.begin(); i<cv.end(); i++ ) {
               ( *i ).updatezvVEC ( in.dtz );
          }
          if ( step%in.frequency==0 ) {
               gutenberg::printXYZ ( in.xyz,in.name,a,true );
          }
          kinz=zKin ( cv );

          iT=2.0*kinz/ ( cv.size() *utils::Kb );

          tT+= ( iT-tT ) / ( step+1.0 );
          if ( step%in.frequency==0 ) {
               for ( vector<colVar>::iterator it=cv.begin(); it<cv.end(); it++ ) {
                    vZ.push_back ( ( *it ).v );
                    KZ.push_back ( ( *it ).mu* ( *it ).v* ( *it ).v/2.0 );
               }
          }

          gutenberg::mdTAMCRepLine ( step*in.dt,aene/nsamp,tT,pressure/nsamp,kinz,cv,in.fancy,cout );
          nsamp=0;
     }

//     solution::histogramWrapper ( vZ,in.bin,"vZtdmc.hist" );
//     solution::histogramWrapper ( KZ,in.bin,"KZtdmc.hist" );

}


void solution::zHMCSampler ( vector< atom >& a, vector<colVar>& cv, generalInputs & in, const int& irseed,
                             vector<double>& pos, vector<double>& avgfz )
{

     boost::mt19937 igen ( irseed );
     boost::uniform_01<boost::mt19937> rn ( igen );

     double epoti,epotf,zpoti,zpotf,kini,kinf;
     int igood,iseed,nsamp;
     double aHg,aH,dHg,fact,vir,acc;

     aHg=0.0;
     aH=0.0;
     nsamp=0;
     for ( int k=0; k<cv.size(); k++ ) avgfz[k]=0.0;
     epoti=PotentialEnergy ( a,in.cut,in.box,in.shift,vir );
     zpoti=zpot ( a,cv,in.box );
     igood=0;
     for ( int i=0; i<in.mccycles; i++ ) {
          iseed=int ( rn() *10000 );
          utils::Print ( "Start Mini-Trajectory ",cout,colours::green );
          utils::Print ( i,cout,colours::green );
          cout <<endl;

          solution::RandomVelocities ( a, in.TBath, in.box,iseed );
          kini=KineticEnergy ( a );
          solution::MDDriverVV ( a,in );
          epotf=PotentialEnergy ( a,in.cut,in.box,in.shift,vir );
          zpotf=zpot ( a,cv,in.box );
          kinf=KineticEnergy ( a );
          fact=- ( epotf-epoti+zpotf-zpoti+kinf-kini ) / ( utils::Kb*in.TBath );
          dHg=- ( epotf-epoti+kinf-kini ) / ( utils::Kb*in.TBath );
          if ( rn() <exp ( fact ) ) {
               igood++;
               epoti=epotf;
               zpoti=zpotf;
               for ( int j=0; j<a.size(); j++ ) {
                    a[j].saveOldPositions ( pos[3*j],pos[3*j+1],pos[3*j+2] );
               }
          } else { //if rejected
               for ( int j=0; j<a.size(); j++ ) {
                    a[j].setPositions ( pos[3*j],pos[3*j+1],pos[3*j+2] );
               }
          }
          utils::Print ( "End Mini-Trajectory ",cout,colours::green );
          utils::Print ( i,cout,colours::green );
          utils::Print ( fact,cout,colours::green );
          cout <<endl;
          if ( ( i%in.frequency==0 ) && ( i>=in.nEquil ) ) {
               nsamp++;
               for ( int k=0; k<cv.size(); k++ ) {
                    cv[k].force ( a,in.box );
                    avgfz[k]+=cv[k].f;
               }
               aHg+= ( exp ( dHg )-aHg ) / ( i-in.nEquil+1.0 );
               aH+= ( exp ( fact )-aH ) / ( i-in.nEquil+1.0 );
               utils::Print ( "runing averages ",cout,colours::red );
               utils::Print ( aHg,cout,colours::red );
               utils::Print ( aH,cout,colours::red );
               cout<<endl;
          }
     }
     utils::Print ( "runing averages ",cout,colours::red );
     utils::Print ( aHg,cout,colours::red );
     utils::Print ( aH,cout,colours::red );
     acc= ( double ) igood/ ( double ) in.mccycles;
     utils::Print ( "acceptance ",cout,colours::green );
     utils::Print ( acc,cout,colours::green );
     cout <<endl;

     for ( int k=0; k<cv.size(); k++ ) {
          cv[k].f=avgfz[k]/nsamp;
     }

}



void solution::histogramWrapper ( std::vector< double >& a, const double& bin, const string& filename )
{

     vector<double> aHist;
     double amin;
     utils::histogram ( a,bin,aHist,amin );
     utils::NormaliseHistogram ( aHist,bin );
     gutenberg::printHistogram ( aHist,amin,bin,filename );

}




// kate: indent-mode cstyle; indent-width 5; replace-tabs on; 
