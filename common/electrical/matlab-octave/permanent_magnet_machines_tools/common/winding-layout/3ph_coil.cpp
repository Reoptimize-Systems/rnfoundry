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


#include "3ph_coil.h"


/*!
 * \brief This constructor makes some initialization for the struct
 * 
 * The coil variable are set equal to -1 and the wire constructor is invoked
 */
coil_3ph::coil_3ph ( ) {

  s   = -1;
  e   = -1;
  nc  = -1;
  w   = -1;
  sec =  0;

}

coil_3ph::coil_3ph (int start, int end, int n, int id_wire, int _sec ) {

  s   = start;
  e   = end;
  nc  = n;
  w   = id_wire;
  sec = _sec;
}

ostream& operator<<(ostream& os, coil_3ph c)
{
    os<<"("<<c.s<<","<<c.e<<","<<c.nc<<","<<c.w<<","<<c.sec<<")"<<endl;
    return os;
}

