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

#include "wire_size.h"

using namespace Koil;

wire_size::wire_size()
{
}

wire_size::wire_size(double r, double d, double ins_d,  double sw)
{
    resistivity        = r;
    diameter           = d;
    insulated_diameter = ins_d;
    specific_weight    = sw;
}

double wire_size::get_resistivity()     {return resistivity;}
double wire_size::get_specific_weight() {return specific_weight;}
double wire_size::get_section()         {return pi*diameter*diameter/4;}
