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

#ifndef UTILS_H
#define UTILS_H
#include <iostream>
#include <vector>
#include <string>

namespace utils {

double const Kb=1.0;
double const zerotol=1.0e-15;
int const w=16;
int const sw=6;
int const pre=8;

//LJ params
double const sigma=1.0;
double const epsilon=1.0;
//small and big
double const zero=1.0e-14;
double const infinity=1.0e69;



int nint ( double const& );
double LennardJones ( const double & );
double LennardJonesdR ( const double& );

void Print ( const double& , std::ostream & ) ;
void Print ( const double& , std::ostream &, const std::string& ) ;
void Print ( const double& , std::ostream &, const int& ,const int& ) ;

void Print ( const int& , std::ostream &, const int& ) ;
void Print ( const int& , std::ostream & ) ;
void Print ( const int& , std::ostream &, const std::string& ) ;
void Print ( const int& , std::ostream &, const std::string&, const int& ) ;


void Print ( const std::string&, std::ostream &, const int& );
void Print ( const std::string& , std::ostream& );
void Print ( const std::string&, std::ostream &, const std::string& );
void Print ( const std::string&, std::ostream &, const int&, const std::string& );

// void Print ( const bool& , std::ostream & ) ;

double AH ( const double&, const double&, const double&, const double&, const double& );

double maximum ( std::vector<double> const& );
double minimum ( std::vector<double> const& );
void histogram ( std::vector<double> const&, double const&, std::vector<double> &, double & );
void NormaliseHistogram ( std::vector<double>&, double const& );
std::string int2string ( const int& );
std::string stringPlusInt ( const std::string &, const int & );
std::string stringPlusDouble ( const std::string &, const double & );

double average ( std::vector<double>& );
double average2 ( std::vector<double>& );
double sigma2 ( std::vector<double>& );
}
#endif // UTILS_H

// kate: indent-mode cstyle; space-indent on; indent-width 4;
