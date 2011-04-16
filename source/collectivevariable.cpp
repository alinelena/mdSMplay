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

#include "collectivevariable.h"
#include "utils.h"


using namespace std;

void colVar::force ( const std::vector< atom >& a, const double& box )
{
    th=this->theta(a,box);
    f=-k* ( z-th );
}

double colVar::theta ( const std::vector< atom >& a, const double& box )
{
    return r ( a[i],a[j],box );
}

void colVar::dtheta ( vector<double> w, const int&t,const std::vector< atom >&a , const double&box )
{
    double aux=0.0,q;
    q=rdr ( w,a[i],a[j],box );
    if ( t==i )
    {
        w[0]=-w[0]/q;
        w[1]=-w[1]/q;
        w[2]=-w[2]/q;
    }
    else if ( t==j )
    {
        w[0]=-w[0]/q;
        w[1]=-w[1]/q;
        w[2]=-w[2]/q;
    }
    else
    {
        w[0]=0.0;
        w[1]=0.0;
        w[2]=0.0;
    }
}

void colVar::updatezvVEC ( const double& h )
{
    v=v+ ( f/mu-gamma*v ) *0.5*h* ( 1.0-0.25*h*gamma ) +utils::AH ( h,xi,eta,gamma,s );
}

void colVar::updatezpVEC ( const double& h, const double& s )
{
    double ss=s*sqrt ( gamma/mu );
    z=z+h*v+ss*eta;
}

double zpot(const vector< atom >&a , vector< colVar >& cv, const double& box){
    double mp,th;
    mp=0.0;
    for (vector<colVar>::iterator it=cv.begin();it<cv.end();it++){
        th=(*it).z-(*it).theta(a,box);
        mp+=(*it).k*th*th;
    }
 return mp/2.0;
}



// kate: indent-mode cstyle; space-indent on; indent-width 4;
