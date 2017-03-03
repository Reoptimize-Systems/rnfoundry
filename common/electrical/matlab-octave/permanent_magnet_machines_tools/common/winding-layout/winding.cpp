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

#include "winding.h"
#include <math.h>

using namespace Koil;

winding::winding(int _Q, int _p) {
    Q = _Q;
    p = _p;

    for (int i=0; i<Q; ++i){
      slot_matrix.push_back(0);
    }
}

void winding::add_coil(coil c) {coils.push_back(c);}

double winding::get_resistance(){

    double R,r = 0;

    for(unsigned int i=0; i<coils.size(); ++i){
        r = coils[i].get_resistance();
        if (r==0) return 0;    //!< return 0 if one coil has 0 resistance (length not assigned)
        else R +=r;
    }
    return R;
}

double winding::get_length(){

    double L,l=0;

    for(unsigned int i=0; i<coils.size(); ++i){
        l = coils[i].get_length();
        if (l==0) return 0;    //!< return 0 if one coil length is not assigned
        else L +=l;
    }
    return L;
}

double winding::get_weight(){

    double W,w=0;

    for(unsigned int i=0; i<coils.size(); ++i){
        w = coils[i].get_weight();
        if (w==0) return 0;    //!< return 0 if one coil length is not assigned
        else W +=w;
    }
    return W;
}

double winding::get_kw(){

    if (p<=0) return 0;
    else return get_kw(p);
}

double winding::get_kw(int nu){

  compute_slot_matrix();
  int val = 0;
  double X=0, Y=0, R=0, N=0, angle = 0;
  double ase = 2 * pi / Q * nu;

  for(int i=0; i<Q; ++i){
    val = slot_matrix[i];
    angle = ase * i;
    if(val<0){val = -val; angle = angle+pi;}
    X = X + val * cos(angle);
    Y = Y + val * sin(angle);
    N = N + val;
  }

  R = sqrt(X*X+Y*Y);
  if(R/N<1e-8) return 0;
  return R/N;
}

bool winding::compute_slot_matrix(){

    if (coils.size()<1) return false;
    for(int i=0; i<Q; ++i) slot_matrix[i] = 0;

    for(int i=0; i<coils.size(); ++i){
        slot_matrix[coils[i].s() - 1] += coils[i].n();
        slot_matrix[coils[i].e() - 1] -= coils[i].n();
    }
    return true;
}

std::vector<int> winding::get_slot_matrix(){

    compute_slot_matrix();
    return slot_matrix;
}

std::vector <coil> winding::get_coils(){ return coils; }

void winding::set_wire(wire _w){

    for(unsigned int i=0; i<coils.size(); ++i)         coils[i].set_wire(_w);
}
