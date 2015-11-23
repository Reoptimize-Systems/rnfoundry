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

#ifndef M_PHASE_WINDING_H
#define M_PHASE_WINDING_H

#include "winding.h"
#include "starofslot.h"


/*! \class mPhaseWinding
 *  \brief This class represents an m-phase winding.
 *  \author Luigi Alberti
 *
 *  \version 1.1
 *  \date    July 2009
 *
 *  This class represents a m-phase winding of an electrical machine.
 *  An m-phase winding is composed by m windings.
 *
 *  The three foundamental data for a m-phase winding are the number of slots stored in Q,
 *  the number of poles pair stored in p and the phase number. The last is stored
 *  in windings.size()
 *
 *
 * The \f$ MMF \f$  of the m-phase winding is expressed as:
   \f[ MMF(x) = \dfrac{ D}{2 \nu}
             \sum_{\nu=1}^\infty
               a_{\nu} \sin \dfrac{2 \nu}{D}  x -
               b_{\nu} \cos \dfrac{2 \nu}{D}  x
  \f]
  where \f$ \nu \f$ is the harmonic order, \f$ D \f$ is the airgap diameter.

 */

class mPhaseWinding
{
public:
    mPhaseWinding();
    ~mPhaseWinding();
    void ComputeWinding(int m, int _Q, int _p, bool SL=false); //!< Fill the winding from the star of slot. If SL is true a single layer winding is create when possible
    void ComputeWinding(int m, int _Q, int _p, int yq);      //!< Fill the winding from the star of slot. A double layer winding with the given coil throw yq is created
    void setData(int Q, int p, double ws = 0.01, double D = 1.0/pi, double hs=-1);             //!< Set the input data
    void setData(double ws = 0.01, double D = 1.0/pi);             //!< Set the input data
    void AddWinding(winding win);           //!< Add a winding to the m-phase winding
    vector<double> Get_MMF_a(int Nmax);     //!< Returns the \f$a_{\nu}\f$ coefficients of the MMF Fourier series expansion
    vector<double> Get_MMF_b(int Nmax);     //!< Returns the \f$b_{\nu}\f$ coefficients of the MMF Fourier series expansion
    void SetCurrents(vector<double>);       //!< Set the m value current which are used to compute the \f$ MMF \f$. If the current vector is empty, a symmetrical current set is created.
    void clear();                           //!< Clear the winding
    vector<double> Get_slot_cur_matrix();   //!< Return the slot_cur_matrix
    int getQ();
    int getp();
    int getm();
    int gett();
    double get_RL_index(double k, double a, double b, double g, double mu, double sigma, int Nmax, double f);
    vector<int> get_nu_symmetrical(int Nmax);           //!< Return the harmonic order of a balanced m-phase winding computed as (+-)k * m +1
    vector<double> get_frnu(double f, int Nmax);        //!< Return the harmonic frequency in the rotor reference frame. f is the stator current frequency
    vector<double> get_omega_rnu(double f, int Nmax);   //!< Return the harmonic mechanical angular speed in the rotor reference frame. f is the stator current frequency
    double get_phase_axis(int m);                       //!< Return the angle of the phasor of phase m, for sequence component computation

    vector <winding> windings;

private:
    int Q,                                  //!< The number of slots
        p;                                  //!< The number of poles pair
    double ws,                              //!< The slot opening width (m)
           D;                               //!< The airgap diameter (m)
    vector<double> currents;                //!< The phases currents to compute the \f$ MMF \f$ (istantaneous values) (A)
    vector<double> slot_cur_matrix;         //!< The slot matrix with the total current of all the phases
    void fill_slot_cur_matrix();            //!< Fill the slot_cur_matrix starting from the windings slot matrix and the given phase currents
    double zero = 1.0e-10;                  //!< The tolerance for zero coefficients in an and bn
    StarOfSlot  star;                       //!< The star of slot object
};

#endif // M_PHASE_WINDING_H
