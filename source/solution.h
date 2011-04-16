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


#ifndef SOLUTION_H
#define SOLUTION_H
#include <vector>
#include <iostream>
#include <string>
#include "atom.h"
#include "generalinputs.h"
#include "collectivevariable.h"

///
/// @brief just a simple namescape containing mainly drivers and init functions
/// @author Alin M Elena
///
namespace solution {


void InitializeSimulation ( std::vector<atom>&, const double& ,  const double &, const int&, std::string const& );
void InitializeSimulation ( std::vector<atom>& , std::string const& );
void InitializeSimulation ( std::vector<atom>&, const double& ,  const double &, const int&, std::string const&, const double& );



/// computes the center of mass velocity
void centerOfMassV ( std::vector<atom>& );
/// scales the velocities to a certain temperature/kinetic energy
void scaleVelocities ( std::vector<atom>&, const double& );

/// comnputes forces on oll atoms
double computeForces ( std::vector<atom> &, const double &, const double &, double& );
///velocity verlet integrator
double velocityVerlet ( std::vector<atom> &, const double &, const double &,const double&, double &, double& );
/// md with velocity verlet
void MDDriverVV ( std::vector<atom>&, generalInputs& );

/// VEC integrator for langevin dynamics
double VECIntegrator ( std::vector<atom> &, const double &, const double &,const double&,
                       const double &, const double&,double &, double& );
/// driver for langevin dynamics with VEC integrator
void MDDriverVEC ( std::vector<atom>&, generalInputs& );

/// driver for TAMD
void MDDriverTAMD ( std::vector<atom>&, std::vector<colVar>&, generalInputs& );
/// TAMD integrator
double TAMDIntegrator ( std::vector<atom> &, std::vector<colVar>&, const double &, const double &,const double&, const double &, const double&,const double&,double &, double& );
/// Computes  forces for TAMD method
double computeForcesTAMD ( std::vector<atom> &, std::vector<colVar>&,const double &, const double &, double& );
///  metropolis monte carlo driver
void MMCDriver ( std::vector<atom>&, generalInputs& );
/// oprtimal displacementin MC so acceptacne rate is around 1/2
void optimalDr(double&, int&,int&, const int&,const int&, double const &);

/// driver for TAMC
void MDDriverTAMC ( std::vector<atom>&, std::vector<colVar>&, generalInputs& );

/// prints to a file the normalised histogram of a set of values
void histogramWrapper(std::vector<double>&, const double&, const std::string& );

}

#endif // SOLUTION_H
// kate: indent-mode cstyle; space-indent on; indent-width 4;
