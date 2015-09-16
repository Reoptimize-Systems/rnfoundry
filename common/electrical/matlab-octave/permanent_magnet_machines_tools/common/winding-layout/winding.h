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

#ifndef WINDING_H
#define WINDING_H

#include <vector>
#include "coil.h"

/*! \class winding
 *  \brief This class represents a set of coils.
 *  \author Luigi Alberti
 *
 *  \version 1.1
 *  \date    July 2009
 *
 *  This class represents a winding of an electrical machine.
 *  A winding is intended as a set of coils all connected in series.
 *  For example it is used to represent one phase of a poliphase machine
 */

class winding
{
public:
    winding(int Q, int p);
    void add_coil(coil c);                     //!< Insert a coil into the winding
    double get_resistance();                   //!< Compute the resistance of the winding
    double get_length();                       //!< Compute the length of the winding
    double get_weight();                       //!< Compute the weight of the winding
    double get_kw(int nu);                     //!< Compute the winding factor of the winding
    double get_kw();                           //!< Compute the winding factor of the winding
    bool compute_slot_matrix();                //!< Compute the slot matrix starting from the coils;
    vector<int> get_slot_matrix();             //!< Get the slot matrix starting
    vector<coil> get_coils();                  //!< Get the vector of the coils
    void set_wire(wire _w);                     //!< Set the wire of all the coils to the giceb wire _w

private:

    int Q;                                     //!< The number of slots
    int p;                                     //!< The pole pairs number
    vector <int> slot_matrix;                  //!< The slot matrix of the winding. The element are the slot conductors
    vector<coil> coils;                        //!< The vector of coils composing the winding


};

#endif // WINDING_H
