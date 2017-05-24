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

#include "wire.h"
#include "wire_size.h"

using namespace Koil;

wire::wire(){ }

void wire::add_wire_size(wire_size size) {wires.push_back(size);}

double wire::get_resistivity(){

    double a,b=0;
    for (unsigned int i=0; i<wires.size(); ++i){
        a += wires[i].get_section();
        b += wires[i].get_section()/wires[i].get_resistivity();
        }
return a/b;
}

double wire::get_section(){

    double S = 0;
    for (unsigned int i=0; i<wires.size(); ++i)  S += wires[i].get_section();
    return S;
}

double wire::get_specific_weight(){

    double W = 0;
    for (unsigned int i=0; i<wires.size(); ++i)  W += wires[i].get_specific_weight();
    return W;
}
