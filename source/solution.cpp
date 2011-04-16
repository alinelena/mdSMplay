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

void solution::InitializeSimulation ( std::vector<atom>& a, const double& temp, const double& box, const int& seed, string const& filename , double const& gamma ) {
    InitializeSimulation ( a,temp,box,seed,filename );
    int k=0;
    for ( vector<atom>::iterator i=a.begin(); i<a.end(); i++ ) {
        ( *i ).setGamma ( gamma );
        ( *i ).setID ( k );
        k++;
        }
    }

void solution::InitializeSimulation ( std::vector<atom>& a, const double& temp, const double& box, const int& seed, string const& filename ) {
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

void solution::InitializeSimulation ( std::vector<atom>& a, string const& filename ) {
    int k=reading::readFull ( filename.c_str(),a );
    if ( k!=0 ) {
        utils::Print ( "Error!!!",cout,colours::red );
        }
    }

void solution::centerOfMassV ( vector< atom >& a ) {
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

void solution::scaleVelocities ( vector< atom >& a, const double& T ) {
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

double solution::computeForces ( vector< atom >& a, const double & rc, const double& box, double& vir ) {
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
                fx=fR*dR[0]/R;
                fy=fR*dR[1]/R;
                fz=fR*dR[2]/R;
                ( *i ).addtoForces ( fx,fy,fz );
                ( *j ).addtoForces ( -fx,-fy,-fz );
                energy+=utils::LennardJones ( R )-ec;
                vir+=fR*R;
                }
            }
        }
    return energy;
    }

double solution::velocityVerlet ( vector< atom >& a, const double& rc, const double& box, const double& dt , double& kin, double &vir ) {
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


void solution::MDDriverVV ( vector< atom >& a, generalInputs& in ) {

    double potEnergy, kinEnergy, vir,pressure;
    potEnergy = PotentialEnergy ( a,in.cut,in.box, vir );
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

double solution::VECIntegrator ( vector< atom >& a, const double&rc , const double& box , const double& dt, const double& s1, const double& s2, double& kin, double& vir ) {

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

void solution::MDDriverVEC ( vector< atom >& a, generalInputs & in ) {
    //set the random generators...
    boost::mt19937 igen ( in.seeds[0] );
    boost::variate_generator<boost::mt19937, boost::normal_distribution<> >  xi ( igen, boost::normal_distribution<> ( 0.0,1.0 ) );
    vector<double> b ( 6 );
    double potEnergy, kinEnergy, vir;


    potEnergy = PotentialEnergy ( a,in.cut,in.box,vir );
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
    gutenberg::mdRepLine ( 0.0,potEnergy/a.size(),kinEnergy/a.size(),iT,tT,pressure,consQ/a.size(),VCM ( a ),in.fancy,cout );
    gutenberg::printXYZ ( in.xyz,in.name,a,true );
    gutenberg::printAtoms ( a,in.debug );

    double s1=sqrt ( utils::Kb*in.TBath*in.dt/6.0 ) *in.dt;
    double s2=sqrt ( 2.0*utils::Kb*in.TBath );

    vector<double> K,E,vx,vy,vz,P;
    K.push_back ( kinEnergy/a.size() );
    E.push_back ( consQ/a.size() );
    P.push_back ( potEnergy/a.size() );
    for ( vector<atom>::iterator it=a.begin(); it<a.end(); it++ ) {
        vx.push_back ( ( *it ).gvx() );
        vy.push_back ( ( *it ).gvy() );
        vz.push_back ( ( *it ).gvz() );
        }
    double aene=0.0, aene2=0.0;
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

        iT=2.0*kinEnergy/ ( 3.0*a.size() *utils::Kb );
        tT+= ( iT-tT ) / ( i+1.0 );


        consQ=potEnergy+kinEnergy;

        pressure=rho* ( utils::Kb*tT ) +vir/volume/3.0;
        gutenberg::mdRepLine ( i*in.dt,potEnergy/a.size(),kinEnergy/a.size(),iT,tT,pressure,consQ/a.size(),VCM ( a ),in.fancy,cout );

        if ( i%in.frequency==0 ) {
            gutenberg::printForces ( a,in.log );
            gutenberg::printXYZ ( in.xyz,in.name,a,true );
            gutenberg::printAtoms ( a,in.debug );
            aene=+consQ;
            aene2+=consQ*consQ;
            nsamp++;
            K.push_back ( kinEnergy/a.size() );
            E.push_back ( consQ/a.size() );
            P.push_back ( potEnergy/a.size() );
            for ( vector<atom>::iterator it=a.begin(); it<a.end(); it++ ) {
                vx.push_back ( ( *it ).gvx() );
                vy.push_back ( ( *it ).gvy() );
                vz.push_back ( ( *it ).gvz() );
                }
            }
        }
    aene2=aene2/nsamp;
    aene=aene/nsamp;
    utils::Print ( "CV:",cout,red );
    utils::Print ( ( aene2-aene*aene ) /utils::Kb/tT/tT,cout, green );
    cout<<endl;
    utils::Print ( "sigma^2 K",cout, red );
    utils::Print ( utils::sigma2 ( K ),cout,green );
    cout<<endl;
    utils::Print ( "sigma^2 E",cout, red );
    utils::Print ( utils::sigma2 ( E ),cout,green );
    cout<<endl;
    solution::histogramWrapper ( K,in.bin,"Kvec.hist" );
    solution::histogramWrapper ( E,in.bin,"Evec.hist" );
    solution::histogramWrapper ( P,in.bin,"Pvec.hist" );
    solution::histogramWrapper ( vx,in.bin,"vxvec.hist" );
    solution::histogramWrapper ( vy,in.bin,"vyvec.hist" );
    solution::histogramWrapper ( vz,in.bin,"vzvec.hist" );
    }

double solution::computeForcesTAMD ( vector< atom >& a, vector<colVar>& cv,const double & rc, const double& box, double& vir ) {
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
                fx=fR*dR[0]/R;
                fy=fR*dR[1]/R;
                fz=fR*dR[2]/R;
                ( *i ).addtoForces ( fx,fy,fz );
                ( *j ).addtoForces ( -fx,-fy,-fz );
                energy+=utils::LennardJones ( R )-ec;
                vir+=fR*R;
                }
            }
        }
    return energy;
    }

double solution::TAMDIntegrator ( vector< atom >& a, vector<colVar>&cv, const double&rc , const double& box , const double& dt, const double& s1,const double& c1, const double& s2, double& kin, double&vir ) {

    kin=0.0;
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
        }
    kin=kin/2.0;
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
void solution::MDDriverTAMD ( vector< atom >& a, vector<colVar>& cv, generalInputs & in ) {
    //set the random generators...

    boost::mt19937 igen ( in.seeds[0] );
    boost::variate_generator<boost::mt19937, boost::normal_distribution<> >  xi ( igen, boost::normal_distribution<> ( 0.0,1.0 ) );
    vector<double> b ( 6 );
    double potEnergy, kinEnergy, vir;

    potEnergy = PotentialEnergy ( a,in.cut,in.box,vir );
    kinEnergy = KineticEnergy ( a );

    double iT,tT;
    double consQ;
    double rho=density ( a, in.box );
    double volume = in.box*in.box*in.box;
    double pressure=rho* ( utils::Kb*in.TBath ) +vir/volume/3.0;

    consQ=potEnergy+kinEnergy;
    iT=2.0/3.0*kinEnergy/a.size() /utils::Kb;
    tT=iT;
    gutenberg::mdTAMDHeader ( cv,cout,in.fancy );
    gutenberg::mdTAMDRepLine ( 0.0,potEnergy/a.size(),kinEnergy/a.size(),iT,tT,pressure,consQ/a.size(),VCM ( a ),cv,in.fancy,cout );
    gutenberg::printXYZ ( in.xyz,in.name,a,true );
    gutenberg::printAtoms ( a,in.debug );

    double s1=sqrt ( utils::Kb*in.TBath*in.dt/6.0 ) *in.dt;
    double c1=sqrt ( utils::Kb*in.Tcv*in.dt/6.0 ) *in.dt;
    double s2=sqrt ( 2.0*utils::Kb*in.TBath );

    vector<double> z;
    z.push_back ( cv[0].z );
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
        potEnergy = solution::TAMDIntegrator ( a,cv,in.cut,in.box,in.dt,s1,c1,s2,kinEnergy, vir );

        iT=2.0/3.0*kinEnergy/a.size() /utils::Kb;
        tT+= ( iT-tT ) / ( i+1.0 );

        if ( i%in.frequency==0 ) {
            gutenberg::printForces ( a,in.log );
            gutenberg::printXYZ ( in.xyz,in.name,a,true );
            gutenberg::printAtoms ( a,in.debug );
            }
        consQ=potEnergy+kinEnergy;
        pressure=rho* ( utils::Kb*in.TBath ) +vir/volume/3.0;
        gutenberg::mdTAMDRepLine ( i*in.dt,potEnergy/a.size(),kinEnergy/a.size(),iT,tT,pressure,consQ/a.size(),VCM ( a ),cv,in.fancy,cout );
        z.push_back ( cv[0].z );
        }
    solution::histogramWrapper ( z,in.bin,"z1.hist" );
    }


void solution::optimalDr ( double& d , int& otry, int& ogood, const int& itry, const int& igood, const double& box ) {

    double f= ( double ) ( igood-ogood ) / ( double ) ( itry-otry );
    double od=d;
    d=d*f/0.5; //0.5 is the target accteptance ratio;
    if ( d/od>1.5 ) d=od*1.5;
    if ( d/od<0.5 ) d=od*0.5;
    if ( d>box/2 ) d=box/2;
    otry=itry;
    ogood=igood;

    }


///
///  @brief the driver for a metropolis montecarlo simulation
///  @param a a vector of atoms
///  @param in contains the input parameters
///
void solution::MMCDriver ( vector< atom >& a, generalInputs& in ) {
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
    double aene=0.0, ap=0.0,aene2=0.0;
    double rho=density ( a,in.box );
    double volume=in.box*in.box*in.box;
    double vir, viri,virf;
    double etot=PotentialEnergy ( a,in.cut,in.box,vir );
//     gutenberg::logLine ( "Volume: ",volume,in.log );
//     gutenberg::logLine ( "Density: ",rho,in.log );


    vector<double>ene;
    ene.push_back ( etot/a.size() );
    double pressure;

    pressure=rho* ( utils::Kb*in.TBath ) +vir/volume/3.0;
    gutenberg::MMCHeader ( cout, in.fancy );
    gutenberg::MMCRepLine ( 0,etot/a.size(),pressure,in.fancy,cout );
    for ( int i=0; i<in.mccycles; i++ ) {
        // one should compute the optimal dr here...
        for ( int j=0; j<in.mcsteps; j++ ) {
            itry++;
            rAtom=int ( rn() *a.size() );
            ei=a[rAtom].myPotentialEnergy ( a,in.cut,in.box,viri );
            a[rAtom].saveOldPositions ( rx,ry,rz );
            a[rAtom].randomMove ( in.dr,rn(),rn(),rn() );
            a[rAtom].putInBox ( in.box );
            ef=a[rAtom].myPotentialEnergy ( a,in.cut,in.box,virf );
            fact=- ( ef-ei ) / ( utils::Kb*in.TBath );
            if ( rn() <exp ( fact ) ) {
                igood++;
                etot=etot+ef-ei;
                vir=vir+virf-viri;
                }
            else {//if rejected
                a[rAtom].setPositions ( rx,ry,rz );
                }
            }
        if ( ( i+1 ) %in.frequency==0 ) {
            //sampling point...
            pressure=rho* ( utils::Kb*in.TBath ) +vir/volume/3.0;
            gutenberg::MMCRepLine ( ( i+1 ) *in.mcsteps,etot/a.size(),pressure,in.fancy,cout );
            gutenberg::printXYZ ( in.xyz,in.name,a,true );
            ene.push_back ( etot/a.size() );
            aene+=etot;
            aene2+=etot*etot;
            ap+=pressure;
            nsamp++;
            solution::optimalDr ( in.dr,otry,ogood,itry,igood,in.box );
            }
        utils::Print ( "|cycle: ",in.log );
        utils::Print ( i+1,in.log );
        utils::Print ( "|attempts: ", in.log );
        utils::Print ( itry,in.log );
        utils::Print ( "|accepted: ", in.log );
        utils::Print ( igood,in.log );
        utils::Print ( "|ratio: ", in.log,red );
        utils::Print ( ( double ) igood/ ( double ) itry,in.log,red );
        in.log<<endl;
        }
    aene2=aene2/nsamp;
    aene=aene/nsamp;
    utils::Print ( "Avg energy:",cout,red );
    utils::Print ( aene/a.size(),cout,green );
    cout<<endl;
    utils::Print ( "Avg pressure:",cout,red );
    utils::Print ( ap/nsamp,cout,green );
    cout<<endl;
    utils::Print ( "CV:",cout,red );
    utils::Print ( 1.5*utils::Kb*a.size() + ( aene2-aene*aene ) /utils::Kb/in.TBath/in.TBath,cout,green );
    cout<<endl;
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
void solution::MDDriverTAMC ( vector< atom >& a, vector<colVar>& cv, generalInputs & in ) {
    //set the random generators...

    boost::mt19937 igen ( in.seeds[0] );
    boost::variate_generator<boost::mt19937, boost::normal_distribution<> >  xi ( igen, boost::normal_distribution<> ( 0.0,1.0 ) );
    boost::uniform_01<boost::mt19937> rn ( igen );
    vector<double> b ( 2 );
    double potEnergy, vir,fact,etot;

    potEnergy = PotentialEnergy ( a,in.cut,in.box,vir ) +zpot ( a,cv,in.box );

    double iT,tT;
    double consQ,kin;
    double rho=density ( a, in.box );
    double volume = in.box*in.box*in.box;
    double pressure=rho* ( utils::Kb*in.TBath ) +vir/volume/3.0;

    int itry=0, otry=0,igood=0, ogood=0;
    double ei,ef, viri,virf;
    int rAtom;
    double rx,ry,rz;

    gutenberg::mdTAMCHeader ( cv,cout,in.fancy );
    gutenberg::mdTAMCRepLine ( 0.0,potEnergy/a.size(),pressure,cv,in.fancy,cout );
    gutenberg::printXYZ ( in.xyz,in.name,a,true );

    double c1=sqrt ( utils::Kb*in.Tcv*in.dt/6.0 ) *in.dt;


    double aene=0.0, ap=0.0;
    int nsamp=0;
    vector<double> avgfz ( cv.size() );
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
                ei=a[rAtom].myPotentialEnergy ( a,in.cut,in.box,viri ) +zpot ( a,cv,in.box );
//                 ei=PotentialEnergy ( a,in.cut,in.box,vir )+zpot ( a,cv,in.box );
                a[rAtom].saveOldPositions ( rx,ry,rz );
                a[rAtom].randomMove ( in.dr,rn(),rn(),rn() );
                a[rAtom].putInBox ( in.box );
                ef=a[rAtom].myPotentialEnergy ( a,in.cut,in.box,virf ) +zpot ( a,cv,in.box );
//                 ef=PotentialEnergy ( a,in.cut,in.box,vir )+zpot ( a,cv,in.box );
                fact=- ( ef-ei ) / ( utils::Kb*in.TBath );
                if ( rn() <exp ( fact ) ) {
                    igood++;
                    etot=etot+ef-ei;
                    vir=vir+virf-viri;
                    }
                else {//if rejected
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

        gutenberg::mdTAMCRepLine ( step*in.dt,etot/nsamp,pressure/nsamp,cv,in.fancy,cout );
        nsamp=0;
        }

    }

void solution::histogramWrapper ( std::vector< double >& a, const double& bin, const string& filename ) {

    vector<double> aHist;
    double amin;
    utils::histogram ( a,bin,aHist,amin );
    utils::NormaliseHistogram ( aHist,bin );
    gutenberg::printHistogram ( aHist,amin,bin,filename );

    }




// kate: indent-mode cstyle; space-indent on; indent-width 4;
