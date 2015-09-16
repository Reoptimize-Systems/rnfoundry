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

#ifndef STAROFSLOT_H
#define STAROFSLOT_H

#include <vector>
#include "wire.h"
class mPhaseWinding;

//! A structure to represent the spokes of the star of slot
typedef struct { int       id; //!< The index of the spoke (the related slot number)
                 double angle;  //<! The angle of the spoke
} spoke;


//! A structure to represent a sector (positive or negative) of a phase measured in counter-clock wise direction
typedef struct sector{ double StartAngle;         //!< The start angle of the phasor
                       double EndAngle;           //!< The end angle of the phasor
                       bool AngleInside(double);  //!< Return true if the given angle is inside the sector
                       void NormalizeAngles();    //!< Makes the angles in the range 0-2*pi
                       vector <int> slot;         //!< The list of slots that belong to the sector
};


/*! \class StarOfSlot
 *  \brief This class represents the star of slot of a m-phase winding.
 *  \author Luigi Alberti
 *
 *  \version 1.1
 *  \date    July 2009
 *
 *  A simmetrical m-phase system is considered with 2*m sectors
 */
class StarOfSlot
{
public:
    StarOfSlot(int m, int Q, int p);
    StarOfSlot();
    void CreateStar();        //!< Create the star of slots
    void CreateSectors();     //!< Create the sectors for each phase
    void PopulateWinding(mPhaseWinding * Win, bool SL = false); //!< Create the winding populating the given mPhaseWinding. If SL is true a single layer winding is created
    void PopulateWinding(mPhaseWinding * Win, int yq);          //!< Create the winding populating the given mPhaseWinding with the given coil throw yq. A double layer winding is created
    int  get_t();             //!< Return the machine periodicity
    bool get_SL_feasible();   //!< Return the single layer feasibility of the winding
    bool get_M_zero();        //!< Return true if the mutual inductance between phases is zero
    int  get_yq();            //!< Return the coil throw
    vector<spoke> get_star(); //!< Return the star

private:
    int m,           //!< The number of phase
        Q,           //!< The number of slots
        p,           //!< The number of poles pair
        yq,          //!< The coil throw
        t;           //!< The machine periodicity
    bool single_layer_feasible;                    //!< Single layer feasibility (computed by this class)
    bool single_layer_wanted;                      //!< If true the single layer winding is computed (if feasible)
    bool mutual_inductance_zero;                   //!< if true the mutual inductance between two phases is zero with yq=1 (computed by this class)
    static constexpr double zero = 1e-4;               //!< The zero for the angles in the star

    vector<spoke>    star;      //!< The vector containing the spoke number label sequence
    vector<sector>   p_sec;     //!< The vector containing the positive sectors of the star. The sectors are m
    vector<sector>   n_sec;     //!< The vector containing the negative sectors of the star. The sectors are m

    int gcd(int a, int b); //!< Returns the great common divisor of two numbers
};

#endif // STAROFSLOT_H
