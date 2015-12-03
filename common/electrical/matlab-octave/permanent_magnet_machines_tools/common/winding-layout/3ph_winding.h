/***************************************************************************
 *   Copyright (C) 2008 by Luigi Alberti and Nicola Bianchi                *
 *  <luigi.alberti@unipd.it> <bianchi@die.unipd.it>                        *
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
#ifndef THREEPH_WINDING
#define THREEPH_WINDING

#include <vector>
#include <string>

#include "3ph_coil.h"

namespace Koil {
    
//using namespace std;


class Threephase_winding{

public:
  Threephase_winding();
  Threephase_winding(int Q, int p); //!< The class constructor; slot_star is called if winding is feasible
  ~Threephase_winding();
  void slot_star();         //!< Compute the star of slot of the main harmonic of order p
  void slot_star(int nu);   //!< Compute the star of slot of the harmonic of order nu. slot_matrix is called.
  void slot_matrix();       //!< Compute the slot matrix and generate the list of phase coils coils_A, coils_B and coils_C
  void set_winding(int Q, int p); //!< Set the winding data and check its feasibility. slot_star is called if winding is feasible 
  void set_yq(int yq);     //!< Change the coil pitch and recompute the slot_matrix
  double kw( );            //!< Compute the winding factor of the main harmonics of order p
  double kw(int nu);       //!< Compute the winding factor of the harmonics of order nu
  int get_Q();
  int get_p();
  int get_t();
  int get_yq();
  void write_slot_matrix(const char * Filename);
  void write_report(char * Filename);

//  private:
  int Q, p, t, yq;
  std::vector<int> star;                              //!< The vector containing the spoke number label sequence
  std::vector<double> angle;                          //!< The vector of the angles (1 2 3 4,...,Q)
  std::vector <double> mat_A, mat_B, mat_C;
  int gcd(int a,int b);
  static const double pi = 3.1415926535897924 ;
  

  std::vector < coil_3ph > coils_A; //!< The list of the coils of phase A
  std::vector < coil_3ph > coils_B; //!< The list of the coils of phase B
  std::vector < coil_3ph > coils_C; //!< The list of the coils of phase C
//  vector < wire > m_wires; //!< The list of the wire present in the winding. This shoul be only one element

};

}; // namespace Koil {

#endif
