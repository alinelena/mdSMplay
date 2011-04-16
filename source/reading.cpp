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


#include "reading.h"
#include <iostream>


using namespace std;
int reading::readXYZ ( string const& filename, vector<atom>& a ) {
    ifstream file ( filename.c_str(), ios::in );
    if ( file.is_open() ) {
        string line;
        int m;
        file >> m;
        getline ( file,line );
        getline ( file,line );
        string el;
        double x,y,z;
        for ( int i=0; i<m; i++ ) {
            file >> el >> x>> y>>z;
            a.push_back ( atom() );
            a[i].setElement ( el );
            a[i].setPositions ( x,y,z );
            }
        file.close();
        return 0;
        }
    else {
        return -1;
        }
    }


int reading::readFull ( string const& filename, vector<atom>& a ) {
    ifstream file ( filename.c_str(), ios::in );
    if ( file.is_open() ) {
        string line;
        int m;
        file >> m;
        getline ( file,line );
        getline ( file,line );
        string el;
        double x,y,z,vx,vy,vz,mass;
        for ( int i=0; i<m; i++ ) {
            file >> el >> mass>>x>>y>>z>>vx>>vy>>vz;
            a.push_back ( atom() );
            a[i].setElement ( el );
            a[i].setPositions ( x,y,z );
            a[i].setMass ( mass );
            a[i].setVelocities ( vx,vy,vz );
            a[i].setForces ( 0.0,0.0,0.0 );
            }
        file.close();
        return 0;
        }
    else {
        return -1;
        }
    }
// kate: indent-mode cstyle; space-indent on; indent-width 4;
