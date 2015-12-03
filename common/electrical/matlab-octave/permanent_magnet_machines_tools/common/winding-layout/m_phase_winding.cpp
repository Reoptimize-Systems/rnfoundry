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

//#include <math>
#include <cmath>
#include "m_phase_winding.h"

using namespace Koil;

const double mPhaseWinding::zero = 1.0e-10;

mPhaseWinding::mPhaseWinding(){

    Q  = -1;
    ws = -1;
    D  = -1;
//     zero = 1.0e-10;
}

void mPhaseWinding::ComputeWinding(int m, int _Q, int _p, bool SL){

    Q = _Q;
    p = _p;
    star =  StarOfSlot(m, Q, p);
    star.CreateSectors();
    star.PopulateWinding(this,SL);
}

void mPhaseWinding::ComputeWinding(int m, int _Q, int _p, int yq){

    Q = _Q;
    p = _p;
    star =  StarOfSlot(m, Q, p);
    star.CreateSectors();
    star.PopulateWinding(this,yq);
}

mPhaseWinding::~mPhaseWinding(){ }

void mPhaseWinding::AddWinding(winding win){     windings.push_back(win);  }

int mPhaseWinding::getm(){ return windings.size(); }
int mPhaseWinding::getQ(){ return Q; }
int mPhaseWinding::getp(){ return p; }
int mPhaseWinding::gett(){ return star.get_t(); }

/*!
  This function returns the \f$ a_{\nu} \f$ coefficients of the MMF
  Fourier series expansion. The maximum harmonic order desidered is given
  as input parameter.
  The coefficients are computed as:
  \f[  a_{\nu} = \dfrac{2}{\pi D} \int_0^C K(x) \cos \dfrac{2 \nu x}{\pi D} dx \f]
*/
std::vector<double> mPhaseWinding::Get_MMF_a(int Nmax){

    std::vector<double> a;
    double sum;

    if (ws <= 0) ws=1e-3; //Use a default value for the slop opening

    if (currents.size()==0){ //! If no current is set, create the default current set to compute the MMF harmonics
        for(int m=0; m < windings.size(); ++m){ currents.push_back(cos(2*pi/windings.size()*m)); }
    }

    fill_slot_cur_matrix();

    for(int nu=1; nu < Nmax+1; ++nu){
        sum = 0;
//        for(int j=0; j < Q; ++j) sum += slot_cur_matrix[j]*cos(nu*2*pi/Q*(j+0.5));
        for(int j=0; j < Q; ++j) sum += slot_cur_matrix[j] / ws * cos(nu*2*pi/Q*(j+0.5));
        //qDebug()<<nu<<sum<<2/pi/nu*sin(nu*ws/D);
        sum = sum * 2/pi/nu*sin(nu*ws/D);
        if (fabs(sum) > zero) a.push_back(sum);
        else                  a.push_back(0);
    }
return a;
}

/*!
  This function returns the \f$ b_{\nu} \f$ coefficients of the MMF
  Fourier series expansion. The maximum harmonic order desidered is given
  as input parameter.
  The coefficients are computed as:
  \f[  b_{\nu} = \dfrac{2}{\pi D} \int_0^C K(x) \sin \dfrac{2 \nu x}{\pi D} dx \f]
*/
std::vector<double> mPhaseWinding::Get_MMF_b(int Nmax){

    std::vector<double> b;
    double sum;

    if (ws <= 0) ws=1e-3;

    if (currents.size()==0){ //! If no current is set, create the default current set to compute the MMF harmonics
        for(int m=0; m < windings.size(); ++m){ currents.push_back(cos(2*pi/windings.size()*m)); }
    }

    fill_slot_cur_matrix();

    for(int nu=1; nu < Nmax+1; ++nu){
        sum = 0;
//        for(int j=0; j < Q; ++j) sum += slot_cur_matrix[j]*sin(nu*2*pi/Q*(j+0.5));
        for(int j=0; j < Q; ++j) sum += slot_cur_matrix[j] / ws * sin(nu*2*pi/Q*(j+0.5));
        sum = sum * 2/pi/nu*sin(nu*ws/D);
        if (fabs(sum) > zero) b.push_back(sum);
        else                  b.push_back(0);
    }
return b;
}

void mPhaseWinding::clear(){

//    Q  = -1;
//    ws = -1;
//    D  = -1;
    currents.clear();
    windings.clear();
    slot_cur_matrix.clear();
}



void mPhaseWinding::SetCurrents(std::vector<double> cur){

    currents.clear(); //!< removes previous values
    for(int i=0; i<cur.size(); ++i) currents.push_back(cur[i]);
}


/*!
This function take the slot matrix of each winding in windings,
compute for each slot the total current considering the current values
in currents and the right number of conductors, and sum all windings in each slot
*/
void mPhaseWinding::fill_slot_cur_matrix(){

    slot_cur_matrix.clear(); //!< removes previous values

    for (int i=0; i<Q; ++i) slot_cur_matrix.push_back(0); //!< build an empty arrat
    int m = windings.size();
    for (int j=0; j<m; ++j){// for each winding/phase
        std::vector<int> slot_matrix = windings[j].get_slot_matrix();
        for (int i=0; i<Q; ++i){// for each slot
            slot_cur_matrix[i] += slot_matrix[i] * currents[j];
        }
    }
}

void mPhaseWinding::setData(int _Q, int _p, double _ws, double _D, double _hs){

    Q = _Q;
    p = _p;
    ws = _ws;
    D = _D;
}

void mPhaseWinding::setData(double _ws, double _D){

    ws = _ws;
    D = _D;
}

std::vector<double> mPhaseWinding::Get_slot_cur_matrix(){ return slot_cur_matrix;  }

double mPhaseWinding::get_RL_index(double k, double a, double b, double g, double mu, double sigma, int Nmax, double f){

    double index = 0;
    std::vector<double> frnu = get_frnu(f,Nmax);
    double kw = windings[0].get_kw();

    for(int i=1; i<Nmax+1; ++i){
        double kgap = k*exp(-((a*g/D)+b )*i);
        double kxi = pi*D*sqrt(pi*mu*mu0*sigma*0.5);
        double xi   = kxi*sqrt(frnu[i-1])/double(i);
        double kwnu = windings[0].get_kw(i);
        index += pow(xi,4)/pow((pow(pow(xi,4)+pow(pi,4),3) ),1./4.)*
                 pow(kwnu/kw,2)*double(i)/p*kgap;
//        qDebug()<< xi<< pow(xi,4)/pow((pow(pow(xi,4)+pow(pi,4),3) ),0.25)<<"   "<<
//                pow(kwnu/kw,2)<<"   " <<double(i)/p<<"   "<<kgap<<"   "<<index;
    }

    return index;
}



std::vector<int> mPhaseWinding::get_nu_symmetrical(int Nmax){
    
    std::vector<int> temp,nu;
    int m       = windings.size();
    int sign = 0;
    int t       = star.get_t();

    temp.push_back(t);     if(p==t) sign = 1.;

    for (int k=1; k<Nmax; ++k){
        temp.push_back((-k * m +1)*t); //! Negative
        temp.push_back(( k * m +1)*t); //! Positive
        if ( (-k * m +1)*t == -p ) sign = -1.;
        if ( ( k * m +1)*t ==  p ) sign =  1.;
        //qDebug()<<k<<p<<(-k * m +1)<<(k * m +1);
        }

    int nu_index = 0;
    for (int k=0; k<Nmax; ++k){
        if (fabs(temp[nu_index])== k+1 ) {
            nu.push_back(temp[nu_index]*sign);
            nu_index++;
            }
        else nu.push_back(0);
    }

    return nu;

}

std::vector<double> mPhaseWinding::get_frnu(double f, int Nmax){

    std::vector<int>    nu   = get_nu_symmetrical(Nmax);
    std::vector<double> wrnu = get_omega_rnu(f,Nmax);
    std::vector<double> frnu;
    for (int k=0; k<Nmax; ++k)  frnu.push_back( fabs(nu[k]*wrnu[k]/2./pi));

    return frnu;
}

std::vector<double> mPhaseWinding::get_omega_rnu(double f, int Nmax){

    std::vector<int>    nu = get_nu_symmetrical(Nmax);
    std::vector<double>  a = Get_MMF_a(Nmax);
    std::vector<double>  b = Get_MMF_b(Nmax);
    std::vector<double>  wrnu;
    for (int k=0; k<Nmax; ++k){
        double c    = sqrt(a[k]*a[k]+b[k]*b[k]);
        double temp = 0;
        if (c>zero) temp = 2.*pi*f*(1./nu[k]-1./p);
        wrnu.push_back(temp);
//        if (nu[k]==0) wrnu.push_back(0);// The harmonic is not present
//        else wrnu.push_back(2*pi*f*(1./nu[k]-1./p));
    }
    return wrnu;


}



double mPhaseWinding::get_phase_axis(int m){

    if (m>windings.size()-1) return -1; //!< check dimensions
    winding Win=windings[m];

    double ase=2*pi/Q*p;                  //!< Electrical slot angle
    std::vector<coil> coils = Win.get_coils(); // Recover all the coils of the winding of index m
    double x=0;
    double y=0;

    for(int i=0; i<coils.size(); ++i){
        coil c =coils[i];
        int n = c.n();
        if (n>0){
            x += int(std::abs(float(n)))*std::cos(c.s()*ase);
            y += int(std::abs(float(n)))*std::sin(c.s()*ase);
            x += int(std::abs(float(n)))*std::cos(c.e()*ase+pi);
            y += int(std::abs(float(n)))*std::sin(c.e()*ase+pi);
            }
        else if(n<0){
            x += int(std::abs(float(n)))*std::cos(c.s()*ase+pi);
            y += int(std::abs(float(n)))*std::sin(c.s()*ase+pi);
            x += int(std::abs(float(n)))*std::cos(c.e()*ase);
            y += int(std::abs(float(n)))*std::sin(c.e()*ase);
        }
    }
return atan2(y,x)*180/pi;



}
