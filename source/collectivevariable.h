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


#ifndef COLLECTIVEVARIABLE_H
#define COLLECTIVEVARIABLE_H
#include <vector>
#include "atom.h"

class colVar
{
public:
    double z;
    double v;
    double f;
    double mu,gamma,k;
    double s;
    double xi,eta;
    double th;
    int i,j;// the collective variable is the distance r_ij
    int type;
    colVar() {};
    void setNoise ( std::vector<double>const &a )
    {
        xi=a[0];
        eta=a[1];
    };
    double theta ( std::vector<atom> const&, const double& );
    void dtheta ( std::vector<double>, const int&, const std::vector<atom>& , const double& );

    void updatezvVEC ( const double & );
    void updatezpVEC ( const double &, const double& );
    void force ( std::vector<atom> const&, const double& );
    friend double zpot(std::vector<atom> const& ,std::vector<colVar> &,const double&);
    virtual ~colVar() {};
};

#endif //COLLECTIVEVARIABLE_H
// kate: indent-mode cstyle; space-indent on; indent-width 4;
