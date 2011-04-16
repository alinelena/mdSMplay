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
#include <algorithm>
#include <sstream>


#include "utils.h"
#include "colours.h"
using namespace std;
int utils::nint ( const double& a ) {

    int i;

    if ( a>=0.0 ) {
        i = int ( a+0.5 );

        }
    else {
        i = int ( a-0.5 );
        }
    return i;

    }

double utils::LennardJones ( const double& r ) {
    double r3=r*r*r/ ( utils::sigma*utils::sigma*utils::sigma );
    double r6=r3*r3;
    double r12=r6*r6;
    return ( 4.0*utils::epsilon* ( 1.0/r12-1.0/r6 ) );
    }

double utils::LennardJonesdR ( const double& r ) {
    double r3=r*r*r/ ( utils::sigma*utils::sigma*utils::sigma );
    double r6=r3*r3;
    double r12=r6*r6;
    return ( -48.0*utils::epsilon/r* ( 1.0/r12-0.5/r6 ) );
    }

void utils::Print ( const double& a, ostream &o, const int& ww, const int&wpre ) {
    o.unsetf ( ios::left|ios::right );
    o.precision ( wpre );
    o.width ( ww );
    if ( ( abs ( a ) <1.0e-3 ) || ( abs ( a ) >1.0e3 ) ) {
        o.setf ( ios_base::floatfield,ios_base::right );
        o<< scientific<<a<<" ";
        }
    else {
        o.setf ( ios_base::floatfield,ios_base::right );
        o<<fixed<<a<<" ";
        }
    }

void utils::Print ( const double& a, ostream &o ) {
    utils::Print ( a,o,utils::w,utils::pre );
    }

void utils::Print ( const double& i, ostream& o, const string& c ) {
    o<<c;
    utils::Print ( i,o );
    o<<colours::end;
    }



void utils::Print ( const int& a, ostream &o ) {
    utils::Print ( a, o, utils::w );
    }

void utils::Print ( const int& a, ostream &o, const int& ww ) {
    o.unsetf ( ios::left|ios::right );
    o.setf ( ios_base::fixed, ios_base::right );
    o.width ( ww );
    o<<dec<<a<<" ";
    }



void utils::Print ( const int& i, ostream& o, const string& c ) {
    o<<c;
    utils::Print ( i,o );
    o<<colours::end;
    }

void utils::Print ( const int& i, ostream& o, const string& c, const int& ww ) {
    o<<c;
    utils::Print ( i,o,ww );
    o<<colours::end;
    }

void utils::Print ( const string& s, ostream& o,const int& ww ) {
    o.unsetf ( ios_base::left|ios_base::right );
    o.setf ( ios_base::left );
    o.width ( ww );
    o << s<<" ";
    }

void utils::Print(const string& s, ostream& o){
    utils::Print ( s,o,utils::w );
}


void utils::Print ( const string& s, ostream& o, const string& c ) {
    o<<c;
    utils::Print ( s,o,utils::w );
    o<<colours::end;
    }

void utils::Print ( const string& s, ostream& o, const int& ww, const string& c ) {
    o<<c;
    utils::Print ( s,o,ww );
    o<<colours::end;
    }

double utils::AH ( const double& h, const double& xi, const double& eta, const double& g, const double& s ) {
    double aux;
    aux=sqrt ( h ) *s;
    return ( 0.5*aux*xi-0.25*h*aux*g* ( 0.5*xi+eta/sqrt ( 3.0 ) ) );
    }

double utils::maximum ( vector<double> const& a ) {
    return *max_element ( a.begin(),a.end() );
    }

double utils::minimum ( vector<double> const& a ) {
    return *min_element ( a.begin(),a.end() );
    }
void utils::histogram ( const std::vector< double >&a , const double&h , std::vector< double >&hist, double& min ) {
    double max=utils::maximum ( a );
    min=utils::minimum ( a );
    int n;
    if ( hist.size() !=0 ) {
        hist.erase ( hist.begin(),hist.end() );
        }
    n= ( int ) ( ( max-min ) /h );
    for ( int i=0; i<=n; i++ ) {
        hist.push_back ( 0.0 );
        }
    int m;
    for ( int i=0; i<a.size(); i++ ) {
        m = ( int ) ( ( a[i]-min ) /h );
        hist[m]++;
        }
    }

void utils::NormaliseHistogram ( vector< double >& h, const double&  bin ) {
    double s=0.0;
    for ( int i=0; i<h.size(); i++ ) {
        s+=bin*h[i];
        }
    for ( int i=0; i<h.size(); i++ ) {
        h[i]=h[i]/s;
        }
    }


string utils::int2string ( const int& n ) {
    stringstream s;
    s << n;
    return s.str();
    }

string utils::stringPlusInt ( const string& w,const int& n ) {
    stringstream s;
    s <<w<< n;
    return s.str();
    }

double utils::average ( std::vector< double >& a ) {
    double sum=0.0;
    for ( vector<double>::iterator it=a.begin(); it<a.end(); it++ ) {
        sum+= ( *it );
        }
    return sum/a.size();
    }

double utils::average2 ( std::vector< double >& a ) {
    double sum=0.0;
    for ( vector<double>::iterator it=a.begin(); it<a.end(); it++ ) {
        sum+= ( *it ) * ( *it );
        }
    return sum/a.size();
    }
double utils::sigma2 ( std::vector< double >& a ) {

    double av=utils::average ( a );
    return utils::average2 ( a )-av*av;;
    }


// kate: indent-mode cstyle; space-indent on; indent-width 4;
