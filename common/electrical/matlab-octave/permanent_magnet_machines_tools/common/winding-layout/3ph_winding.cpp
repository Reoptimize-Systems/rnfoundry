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

#include "3ph_winding.h"
#include <math.h>
#include <fstream>

using namespace Koil;

Threephase_winding::Threephase_winding()
{

    Q  = -1;
    t  = -1;
    p  = -1;
    yq = -1;
}

Threephase_winding::Threephase_winding(int _Q, int _p)
{

    set_winding(_Q , _p);

}

Threephase_winding::~Threephase_winding() {}

void Threephase_winding::set_winding(int _Q, int _p)
{

    Q  = _Q; // the total number of coils Qc
    p  = _p; // the total number of poles
    t = gcd(Q,p); // the greatest common denominator of coils and poles
    yq = Q/2/p; // slot pitch of coils
    if (yq<1) yq = 1;


    //! Winding feasibility check (only for three phase winding)
    if ( Q/(3*t) != (Q+0.0)/(3*t))
    {
        t = -1;
    }
    else slot_star();

}



void Threephase_winding::slot_star()
{
    slot_star(p);
}

void Threephase_winding::slot_star(int nu)
{

    int p = nu;


//! clear all the quantities that will be computed
    star.clear();
    angle.clear();

    double alpha_se = 2 * pi / Q * p;              //!< Slot electrical angle
    double alpha_ph = 2*pi / Q *t;            //!< Angle between two adiacent phasor (or star spoke)

    double epsilon;
// fractpart = modf (param , &intpart);
    epsilon = 0;
    if ( Q/(6*t) == (Q+0.0)/(6*t))                 //!< shift of the angle if there is overlapping with the sector
    {
        epsilon = - alpha_ph /4;
    }



    /*! First an array with all the angle and the sequence of label is create.
     *  Then the angle are sorted in order to achieve the correct sequence of spoke labels.
     *  At the same time, also the array of label are sorted.
     */
    for(int i=0; i<Q; ++i)          //!< Array creation
    {
        angle.push_back(alpha_se*i+epsilon);
        star.push_back(i+1);
    }

    for(int i=0; i<Q; ++i)          //!< Makes all the angles in the rang [0 2*pi]
    {
        while (angle[i]>=2*pi) angle[i] = angle[i] - 2 * pi;
        if (fabs(angle[i]-2*pi)<0.001) angle[i]=0; //!< makes 0 the angle ~ 2pi
    }


    double swap;
    int swap_int;

    for(int i=0; i<Q; ++i)          //!< Sorting of the arrays
    {
        for(int j=i; j<Q; ++j)
        {
            if (angle[j]<angle[i])
            {
                swap     = angle[i];
                angle[i] = angle[j];
                angle[j] = swap;
                swap_int = star[i];
                star[i]  = star[j];
                star[j]  = swap_int;
            }
        }
    }
    slot_matrix();
}// slot_star


void Threephase_winding::slot_matrix()
{

//! Clean up of the previous data
    mat_A.clear();
    mat_B.clear();
    mat_C.clear();
    coils_A.clear();
    coils_B.clear();
    coils_C.clear();

//! slot matrix initialization
    for(int i=0; i<Q; ++i)
    {
        mat_A.push_back(0); // both the layer
        mat_B.push_back(0);
        mat_C.push_back(0);
    }


//! Subdivision of the star in the six sector
//! and coil arrays population
    const double cos_30 = cos(pi/6);
    int index;
    for(int i=0; i<Q; ++i)
    {
        index = star[i]+yq;
        if (index>Q) index=index-Q;

        if (cos(angle[i]) > cos_30)
        {
            coils_A.push_back(/*coil_3ph::*/coil_3ph(star[i],index,1,1,1));   // A+
        }
        if (cos(angle[i]) < -cos_30)
        {
            coils_A.push_back(/*coil_3ph::*/coil_3ph(index,star[i],1,1,-1));   // A-
        }
        if ( (cos(angle[i]) < 0) && (sin(angle[i]) > 0.5) )
        {
            coils_B.push_back(/*coil_3ph::*/coil_3ph(star[i],index,1,1,1));   // B+
        }
        if ( (cos(angle[i]) > 0) && (sin(angle[i]) <-0.5) )
        {
            coils_B.push_back(/*coil_3ph::*/coil_3ph(index,star[i],1,1,-1));   // B-
        }
        if ( (cos(angle[i]) < 0) && (sin(angle[i]) <-0.5) )
        {
            coils_C.push_back(/*coil_3ph::*/coil_3ph(star[i],index,1,1,1));   // C+
        }
        if ( (cos(angle[i]) > 0) && (sin(angle[i]) > 0.5) )
        {
            coils_C.push_back(/*coil_3ph::*/coil_3ph(index,star[i],1,1,-1));   // C-
        }

    }

//! Slot matrix computation
    for(int i=0; i<coils_A.size(); ++i)
    {
        mat_A[coils_A[i].s-1]=mat_A[coils_A[i].s-1] + 0.5;// coils_A[i].nc;
        mat_A[coils_A[i].e-1]=mat_A[coils_A[i].e-1] - 0.5;// coils_A[i].nc;
    }

    for(int i=0; i<coils_B.size(); ++i)
    {
        mat_B[coils_B[i].s-1]=mat_B[coils_B[i].s-1] + 0.5;// coils_B[i].nc;
        mat_B[coils_B[i].e-1]=mat_B[coils_B[i].e-1] - 0.5;// coils_B[i].nc;
    }

    for(int i=0; i<coils_C.size(); ++i)
    {
        mat_C[coils_C[i].s-1]=mat_C[coils_C[i].s-1] + 0.5;// coils_C[i].nc;
        mat_C[coils_C[i].e-1]=mat_C[coils_C[i].e-1] - 0.5;// coils_C[i].nc;
    }

}

void Threephase_winding::write_slot_matrix(const char * Filename)
{

// export the slot matrix in a file

    ofstream out(Filename);
// Setup the output format
// out.setf(ios::scientific, ios::floatfield);
    out.setf(ios::fixed, ios::floatfield);
    out.setf(ios::right, ios::adjustfield);
// out.setf(ios::showpos);
    out.precision(1);
    out.fill('0');

    out<< "-- slot matrix of phase A. Elements: "<<mat_A.size()<<endl;
    out<< "ka = {"<<endl;
    int cont = 1;
    for(int i=0; i<mat_A.size()-1; i++)
    {
        if(cont<10)
        {
            out<<mat_A[i]<<", ";
            cont++;
        }
        else
        {
            out<<mat_A[i]<<","<<endl;
            cont =1;
        }
    }
    out<< mat_A[mat_A.size()-1]<<"} \n\n";

    out<< "-- slot matrix of phase B. Elements: "<<mat_B.size()<<endl;
    out<< "kb = {"<<endl;
    cont = 1;
    for(int i=0; i<mat_B.size()-1; i++)
    {
        if(cont<10)
        {
            out<<mat_B[i]<<", ";
            cont++;
        }
        else
        {
            out<<mat_B[i]<<","<<endl;
            cont =1;
        }
    }
    out<< mat_B[mat_B.size()-1]<<"} \n\n";

    out<< "-- slot matrix of phase C. Elements: "<<mat_C.size()<<endl;
    out<< "kc = {"<<endl;
    cont = 1;
    for(int i=0; i<mat_C.size()-1; i++)
    {
        if(cont<10)
        {
            out<<mat_C[i]<<", ";
            cont++;
        }
        else
        {
            out<<mat_C[i]<<","<<endl;
            cont =1;
        }
    }
    out<< mat_C[mat_C.size()-1]<<"} \n";
}


void Threephase_winding::write_report(char * Filename)
{

    ofstream out(Filename);

    out.setf(ios::fixed, ios::floatfield);
    out.setf(ios::right, ios::adjustfield);
// out.setf(ios::showpos);
    out.precision(3);
    out.fill('0');

    out<<"--------------------------------------------------------------------------------"<<endl;
    out<<"Winding report generated from Koil"<<endl;
    out<<"--------------------------------------------------------------------------------"<<endl<<endl;

    out<<"Phase A coils ("<<coils_A.size()<<" elements):"<<endl;
    for(int i=0; i<coils_A.size(); ++i)
        out<<coils_A[i];//.s<<"  "<<coils_A[i].e<<endl;
    out<<endl;

    out<<"Phase B coils ("<<coils_B.size()<<" elements):"<<endl;
    for(int i=0; i<coils_B.size(); ++i)
        out<<coils_B[i];//.s<<"  "<<coils_B[i].e<<endl;
    out<<endl;

    out<<"Phase C coils ("<<coils_C.size()<<" elements):"<<endl;
    for(int i=0; i<coils_C.size(); ++i)
        out<<coils_C[i];//.s<<"  "<<coils_C[i].e<<endl;
    out<<endl;

}

void Threephase_winding::set_yq(int _yq)
{

    yq = _yq;
    slot_matrix();

}

double Threephase_winding::kw()
{

    return kw(p);
}

double Threephase_winding::kw(int nu)
{

    if(fmod (nu ,t* 3)<0.1) return 0; // return zero for harmonic of order n*3


    double X = 0;
    double Y = 0;
    double R = 0;
    double N = 0;
    double ase = 2 * pi / Q * nu;
    double val = 0;
    double angle;

    for(int i=0; i<Q; ++i)
    {
        val = mat_A[i];
        angle = ase * i;
        if(val<0)
        {
            val=-val;
            angle=angle+pi;
        }
        X = X + val * cos(angle);
        Y = Y + val * sin(angle);
        N = N + val;
    }

    R = sqrt(X*X+Y*Y);

    if(R/N<1e-8) return 0;

    return R/N;
}


int Threephase_winding::get_yq()
{
    return yq;
}


int Threephase_winding::get_t()
{
    return t;
}

int Threephase_winding::gcd(int a,int b)
{

    while (b > 0)
    {
        a = a % b;
        a ^= b;
        b ^= a;
        a ^= b;
    }
    return a;
}









