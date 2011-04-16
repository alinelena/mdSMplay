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


#ifndef ATOM_H
#define ATOM_H
#include <string>
#include <vector>
#include <fstream>
class atom {
        double x,y,z;
        double vx,vy,vz;
        double fx,fy,fz;

        double gamma;
        std::vector<double> xi;
        std::vector<double> eta;

        double mass;
        std::string name;
        int type;
        int id;
    public:
        atom() {};
        void print ( std::ostream & ) const;
        void setID ( const int& i ) {
            id=i;
            };
        void printXYZ ( std::ostream & ) const;
        void printForce ( std::ostream & ) const;
        void setElement ( std::string const& );
        void setMass ( double const& );
        void setGamma ( double const&a ) {
            gamma=a;
            };
        void setPositions ( double const&,double const&,double const& );
        void setVelocities ( double const&,double const&,double const& );
        void setForces ( double const&,double const&,double const& );
        void addtoForces ( double const&,double const&,double const& );
        void addtoVelocities ( double const&,double const&,double const& );
        void setpPositions ( );
        void setpVelocities ();
        void setpForces ( );
        void updatePositionsVV ( const double& );
        void updateVelocitiesVV ( const double& );
        void putInBox ( const double& );
        void scaleVelocity ( const double&,const double&,const double& );
        void setNoise ( const std::vector<double>& );
        void updateVelocitiesVEC ( const double &, const double& );
        double v2();
        double gvx() {
            return vx;
            };
        double gvy() {
            return vy;
            };
        double gvz() {
            return vz;
            };
        double vx2() {
            return vx*vx;
            };
        double vy2() {
            return vy*vy;
            };
        double vz2() {
            return vz*vz;
            };
        double getMass();
        double px();
        double py();
        double pz();
        double rfx() const {
            return fx;
            };
        double rfy() const {
            return fy;
            };
        double rfz() const {
            return fz;
            };
        double force() const;
        void randomMove ( double const&, double const&, double const&, double const& );
        friend double r ( atom const&, atom const& );
        friend double r ( atom const&, atom const&, double const& );
        friend void dr ( std::vector<double>&, atom const&, atom const&, double const& );
        friend double rdr ( std::vector<double>&, atom const&, atom const&, double const& );
        friend double KineticEnergy ( std::vector<atom>& );
        friend double PotentialEnergy ( std::vector<atom> const& , const double& , const double&, double& );
        friend double VCM ( std::vector<atom>& );
        friend double density ( std::vector<atom>&, const double& );
        virtual ~atom() {
            xi.erase ( xi.begin(),xi.end() );
            eta.erase ( eta.begin(),eta.end() );
            };
        void updatePositionsVEC ( const double &, const double & );
        void saveOldPositions ( double&, double&, double& );
        void dr ( std::vector<double>&,atom const&, double const& );
        double rdr ( std::vector<double>&, atom const&, double const& );
        double myPotentialEnergy ( std::vector<atom> & , const double& , const double&, double& );
    };

#endif // ATOM_H
// kate: indent-mode cstyle; space-indent on; indent-width 4;
