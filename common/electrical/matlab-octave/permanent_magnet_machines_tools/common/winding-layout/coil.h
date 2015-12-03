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

#ifndef COIL_H
#define COIL_H

#include "wire.h"

namespace Koil {
    
/*! \class coil
 *  \brief This class represents a coil of the winding.
 *  \author Luigi Alberti
 *
 *  \version 1.1
 *  \date    July 2009
 *
 *  This class represents a coil of the winding. All the data of the coil are stored.
 *  The coil is used to compose a winding.
 *  An m-phase winding is composed by m windings.
 *
 * The Coil resistance is computed as
   \f[ R     = \rho \dfrac{L}{S}  \f]
   where \f$ \rho \f$ and \f$ S \f$ are given from the actual wire of the coil.

   \f$ L \f$ is computed as:
   \f[ L= 2 \cdot n_c \cdot \left(L_{stk}+\dfrac{\pi(D+h_s)}{Q}\cdot y_q +2 h_{ew} \right)\f]
   where
   \f$ n_c\f$ is the number of coil conductors,
   \f$ L_{stk} \f$ is the stack lenght,
   \f$ D \f$ is the stator inner diameter,
   \f$ h_s \f$ is the stator slot height,
   \f$ Q \f$ is the number of slots,
   \f$ y_q \f$ is the coil throw and
   \f$ h_{ew} \f$ is the average end winding lenght.

 */

class coil{
public:
    coil();
    coil(int start, int end, wire mWire, int nc);
    coil(int start, int end, wire mWire, int nc, double Lstk, double hew);

    void set_data(double Lstk, double D, double hs, int Q, double hew);
    double get_resistance();                     //!< Compute the resistance of the coil
    double get_length();                         //!< Compute the length of the coil
    double get_weight();                         //!< Compute the weight of the coil
    void   set_wire(wire _w);                    //!< Set the wire of the coil
    int s();
    int e();
    int n();

private:
    int  start;                                //!< The beginning slot
    int  end;                                  //!< The final slot
    wire mWire;                                //!< The wire description
    int  nc;                                   //!< Number of turns
    double Lstk;                               //!< The stack length
    double hew;                                //!< The average height of the endwindings (m)
    double D;                                  //!< The stator inner diameter (m)
    double hs;                                 //!< The stator slot height (m)
    int Q;                                     //!< The number of slots
    bool LENGTH_SET;                           //!< if false the resistance &co are not computed

};

}; // namespace Koil {

#endif // COIL_H
