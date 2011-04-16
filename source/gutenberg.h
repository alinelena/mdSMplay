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


#ifndef GUTENBERG_H
#define GUTENBERG_H
#include <string>
#include <vector>
#include <fstream>

#include "atom.h"
#include "collectivevariable.h"
#include "generalinputs.h"

namespace gutenberg {
int printXYZ ( std::string const&, std::string const&, std::vector<atom> const&, bool const& );
void printAtoms ( std::vector<atom> const&, std::ostream & );
void printForces ( std::vector<atom> const&, std::ostream & );
void mdHeader ( std::ostream&, bool const & );
void mdRepLine ( double const&, double const&,double const&,double const&,double const&,double const&,double const&,double const&,bool const&, std::ostream& );
int printHistogram ( std::vector<double> const&, double const&,double const&,std::string const& );

void mdTAMDHeader ( const std::vector<colVar>&, std::ostream&, bool const & );
void mdTAMDRepLine ( double const&, double const&,double const&,double const&,double const&,double const&,double const&,double const&,const std::vector<colVar>&,bool const&, std::ostream& );

void MMCRepLine ( const int&,double const&,double const&,const bool&,std::ostream& );
void MMCHeader ( std::ostream&, bool const & );

void mdTAMCHeader ( const std::vector<colVar>&, std::ostream&, bool const & );
void mdTAMCRepLine ( double const&, double const&,double const&,const std::vector<colVar>&,bool const&, std::ostream& );
void logInput ( generalInputs& );
void logLine ( std::string const&, const double&, std::ostream& );
void logLine ( std::string const&, const int&, std::ostream& );
void logLine ( std::string const&, std::string const&, std::ostream& );
void logLine ( std::string const&, std::string const&, std::ostream&, const int& );
void timeStamp ( std::string const&,std::ostream& );
}

#endif // GUTENBERG_H

// kate: indent-mode cstyle; space-indent on; indent-width 4;
