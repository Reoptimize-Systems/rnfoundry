/***************************************************************************
 *   Copyright (C) 2009 by Luigi Alberti                                   *
 *   luigi.alberti@unipd.it                                                *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "coil.h"
#include "constant.h"
#include <math.h>

coil::coil() {LENGTH_SET = false;}

coil::coil(int s, int e, wire w, int n)
{
    start = s;
    end   = e;
    mWire = w;
    nc    = n;
    LENGTH_SET = false;
}

coil::coil(int s, int e, wire w, int n, double _L, double _hew)
{
    start = s;
    end   = e;
    mWire = w;
    nc    = n;
    Lstk  = _L;
    hew   = _hew;
    LENGTH_SET = true;
}



double coil::get_resistance(){

return  mWire.get_resistivity() *  get_length()  / mWire.get_section();

}

double coil::get_length(){

    int yq = fabs(end-start);
    if (yq > Q/2) yq = fabs(Q-yq);

//! If the length has been not assigned returns 0

    if (LENGTH_SET) return 2 * nc * (Lstk + pi/Q*(D+hs)*yq + 2*hew);
    else return 0;
}

double coil::get_weight(){

return mWire.get_specific_weight()*get_length();

}

void coil::set_data(double _Lstk, double _D, double _hs, int _Q, double _hew){

    Lstk     = _Lstk;
    D        = _D;
    hs       = _hs;
    Q        = _Q;
    hew      = _hew;
    LENGTH_SET = true;
}

int coil::s(){return start;}
int coil::e(){return end;}
int coil::n(){return nc;}


void coil::set_wire(wire _w){ mWire = _w;}
