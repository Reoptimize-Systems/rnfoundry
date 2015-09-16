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

#ifndef THREEPH_COIL_H
#define THREEPH_COIL_H

#include <iostream>

using namespace std;

/** \brief This struct represents the coil of a winding of an electrical machine.
  *
  * This class represent the coil of a winding of an electrical machine.
  * A coils is defined by the beginning slot, the end slot and the number
  * of turns.
  * 
  */
struct coil_3ph{
  
  coil_3ph ( );                           //!< The constructor for initialization
  coil_3ph (int, int, int, int, int );    //!< The constructor for creation
  int s;                              //!< Beginning slot. Current IN
  int e;                              //!< End slot. Current OUT  
  int nc;                             //!< Number of turns of the coil
  int w;                              //!< The wire that forms the coil
  int sec;                            //!< The slot star sector  of the coil 1: positive, -1: negative

  friend ostream& operator<<(ostream& os, coil_3ph c);
};

#endif // THREEPH_COIL_H
