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

#include "starofslot.h"
#include "m_phase_winding.h"
#include <math.h>

using namespace Koil;

const double StarOfSlot::zero = 1e-4;

bool sector::AngleInside(double angle){

//    qDebug()<<StartAngle<<"  "<<angle<<"  "<<EndAngle;
    if(StartAngle<0){ //! The first sector is always centered on the origin axis and has to be trated carefully
        double add = EndAngle;
        angle += add;
        while (angle > 2*pi) angle -= 2*pi;
        if (angle > 0 && angle < EndAngle+add) return true;
    }
    if (angle > StartAngle && angle < EndAngle) return true;
    return false;
}

void sector::NormalizeAngles(){

    while (StartAngle > 2*pi) StartAngle -= 2*pi;
    while (EndAngle   > 2*pi) EndAngle   -= 2*pi;
}

StarOfSlot::StarOfSlot(){

    m                      = -1;
    Q                      = -1;
    p                      = -1;
    yq                     = -1;
    t                      = -1;
    single_layer_feasible  = false;
    single_layer_wanted    = false;
    mutual_inductance_zero = false;
//     zero = 1e-4;
}

StarOfSlot::StarOfSlot(int _m, int _Q, int _p){

//     zero = 1e-4;
    m = _m;
    Q = _Q;
    p = _p;

    t = gcd(Q,p); //! Compute the machine periodicity
    yq = Q/2/p;   //! Compute the theoretical coil throw
    if (yq<1) yq = 1;
    mutual_inductance_zero = false;

  //! Winding feasibility check
  if ( Q/(m*t) != double(Q)/(m*t))    t = -1;

  else { //!<  Q/ (m*t) is integer and  the winding is feasible

    //! Single layer feasibility and mutual inductance check
    if(gcd(Q,2)!=2) single_layer_feasible = false; //! Q is odd. SL winding is feasible only if Q is even

    else{//! Q is even
          if ( gcd (Q/t , 2 ) != 2 ){ //!  Q/t odd
              mutual_inductance_zero = false;
              if   ( (t%2)!=0 )        single_layer_feasible = false; //! Q/t is odd and t is odd
              else                     single_layer_feasible = true;  //! Q/t is odd and t is even
            }
          else{ //!  Q/t even
            if ( gcd(Q/(2*t) , 2) ==2 )  { //!  Q / ( 2t ) even
                single_layer_feasible = true;              //!Q/t is even and Q/(2t) is even";
                if (yq==1) mutual_inductance_zero = true;  //! M = to 0 for both SL and DL winding if yq = 1
                }
            else { //!  Q / ( 2t ) odd
                    single_layer_feasible = true;          //! Q/t is even and Q/(2t) is odd
                    if (yq==1) mutual_inductance_zero = true;  //! M = 0 if yq = 1
                    if (single_layer_wanted) mutual_inductance_zero = false;  //! M = 0 only for double layer winding
                   }
                }
        }
    CreateStar(); //! If the winding is feasible the star is created
    }
}

//! Create the star of slot populating the spoke
void StarOfSlot::CreateStar(){

//    int p = nu;

    //! clear all the quantities that will be computed
    star.clear();

    double alpha_se = 2 * pi / Q * p;         /*!< Slot electrical angle*/
    double alpha_ph = 2 * pi / Q * t;         //!< Angle between two adiacent phasor (or star spoke)

    double epsilon = 0;
    //! Shift of the angle if there is overlapping between spokes and the sector border
    if ( Q/(2*m*t) == double(Q)/(2*m*t)) epsilon = - alpha_ph /4;


    /*! First an array with all the angle and the sequence of label is create.
    *  Then the angle are sorted in order to achieve the correct sequence of spoke labels.
    *  At the same time, also the array of labels is sorted.
    */
    for(int i=0; i<Q; ++i){         //!< Array creation
        spoke s;
        s.angle = alpha_se*i+epsilon;
        s.id    = i+1;
        star.push_back(s);
    }

    for(int i=0; i<Q; ++i){         //!< Makes all the angles in the rang [0 2*pi]
      while (star[i].angle>=2*pi) star[i].angle = star[i].angle - 2 * pi;
      if (fabs(star[i].angle-2*pi) < zero) star[i].angle = 0.0; //!< makes 0 the angle ~ 2pi
    }


    spoke swap;

    for(int i=0; i<Q; ++i){         //!< Sorting of the arrays
      for(int j=i; j<Q; ++j){
        if (star[i].angle < star[i].angle) {
          swap     = star[i];
          star[i]  = star[j];
          star[j]  = swap;
          }
        }
    }
}

void StarOfSlot::CreateSectors(){

    p_sec.clear();
    n_sec.clear();

    for (int i=0; i<m; ++i){ //! Create the 2*m sectors
        sector sec;
        sec.StartAngle = pi/m*(2*i-0.5);
        sec.EndAngle   = pi/m*(2*i+0.5);
        sec.NormalizeAngles();
        p_sec.push_back(sec); //! positive
        sec.StartAngle += pi;
        sec.EndAngle   += pi;
        sec.NormalizeAngles();
        n_sec.push_back(sec); //! negative
    }

    for (unsigned int i=0; i<star.size() ; ++i){ //! Populate the sectors with the slots number
        //! Check for positive sectors
        for (unsigned int j=0; j<p_sec.size() ; ++j){
            if (p_sec[j].AngleInside(star[i].angle))
                p_sec[j].slot.push_back(star[i].id);
        }
        for (unsigned int j=0; j<n_sec.size() ; ++j){
            if (n_sec[j].AngleInside(star[i].angle))
                n_sec[j].slot.push_back(star[i].id);
        }
    }

}


int StarOfSlot::gcd(int a,int b) {

  while (b > 0) {
    a = a % b;
    a ^= b;
    b ^= a;
    a ^= b;
    }
  return a;
}

void StarOfSlot::PopulateWinding(mPhaseWinding * _Win, int _yq){

    yq = _yq;
    PopulateWinding(_Win, false);
}

void StarOfSlot::PopulateWinding(mPhaseWinding * Win, bool SL){

    if (single_layer_feasible == false) SL = false;
    if (SL && (yq %2 == 0 )) yq --; //!< yq must be odd

    Win->clear();//! Erase all previous data
    wire _w;
    _w.add_wire_size(wire_size(1,1,1,1));
    int nc = 1;

    for (int i=0; i<m; ++i){
        winding w(Q, p);
        //! Add the coils of the positive sectors
        for (unsigned int j=0; j<p_sec[i].slot.size(); ++j){
            if (SL && (p_sec[i].slot[j]%2==0) ) {
//                qDebug()<<p_sec[i].slot[j];
                continue;
            }
            int end = p_sec[i].slot[j]+yq;
            while (end  > Q) end -= Q;
            coil c(p_sec[i].slot[j], end, _w, nc);
            w.add_coil(c);
        }
        //! Add the coils of the negative sectors
        for (unsigned int j=0; j<n_sec[i].slot.size(); ++j){
            if (SL && (n_sec[i].slot[j]%2==0) ) {
//            qDebug()<<n_sec[i].slot[j];
            continue;
        }
            int end = n_sec[i].slot[j]+yq;
            while (end  > Q) end -= Q;
            coil c(n_sec[i].slot[j], end, _w, -nc);
            w.add_coil(c);
        }
        Win->AddWinding(w);
    }
}

int StarOfSlot::get_t(){return t;}

int StarOfSlot::get_yq(){return yq;}

bool StarOfSlot::get_M_zero(){return mutual_inductance_zero;}

bool StarOfSlot::get_SL_feasible(){return single_layer_feasible;}

std::vector<spoke> StarOfSlot::get_star(){return star;}
