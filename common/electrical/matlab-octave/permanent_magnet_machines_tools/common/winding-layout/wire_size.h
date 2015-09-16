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

#ifndef WIRE_SIZE_H
#define WIRE_SIZE_H

#include "constant.h"

/*! \class coil
 *  \brief This class represents a wire size which compose a wire
 *  \author Luigi Alberti
 *
 *  \version 1.1
 *  \date    July 2009
 *
 *  This class represents a wire size of the coil \c wire. A wire can be composed by many \c wire_size
 */

class wire_size
{
public:
    wire_size();
    wire_size(double resistivity, double diameter, double insulated_diameter, double specific_weight);
    double get_resistivity();               //!< Ohm mm^2/m
    double get_section();                   //!< mm^2
    double get_specific_weight();           //!< kg/m


private:
    double resistivity;                     //!< Ohm mm^2/m
    double diameter;                        //!< mm
    double insulated_diameter;              //!< mm
    double specific_weight;                 //!< kg/m


};

#endif // WIRE_SIZE_H
