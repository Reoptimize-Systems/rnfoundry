#include <cmath>
#include "hydrobody.h"

#ifdef DEBUG
#include <iostream>

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#include <string>       // std::string
#include <sstream>      // std::stringstream
std::stringstream mexoutstream;
#define MEXSTREAM(X)                               \
{                                                  \
    mexoutstream << X << std::endl;                \
    std::string mexoutstr = mexoutstream.str();    \
    mexPrintf ("%s", mexoutstr.c_str ());          \
    mexoutstream.str("");                          \
    mexoutstream.clear();                          \
}
#endif // MATLAB_MEX_FILE

#endif // DEBUG

using namespace Eigen;

namespace wsim
{

    hydroBody::hydroBody()
    {
        _isInititialised = false;
        _stepCount = 0;
        _timeStepHist = convArray1Hlend::Zero ();
        _radForceOldF_FM = Vector61d::Zero ();
        _prev_waveElevation = 0;

        _doMorrisonElementViscousDrag = false;
        _doNonLinearFKExcitation = false;

        F_Total = Vector61d::Zero ();
        F_ExcitLin = Vector61d::Zero ();
        F_Excit = Vector61d::Zero ();
        F_ExcitRamp = Vector61d::Zero ();
        F_LinearDamping = Vector61d::Zero ();
        F_ViscousDamping = Vector61d::Zero ();
        F_AddedMass = Vector61d::Zero ();
        F_RadiationDamping = Vector61d::Zero ();
        F_Restoring = Vector61d::Zero ();
            //Eigen::Matrix6Nd BodyHSPressure;
        F_ExcitNonLin = Vector61d::Zero ();
            //Matrix6Nd WaveNonLinearPressure;
            //Eigen::MatrixXd WaveLinearPressure;
        F_MorrisonElement = Vector61d::Zero ();

    }

    hydroBody::~hydroBody()
    {
        //dtor
    }

    void hydroBody::init( const unsigned int &bodyNumber,
                          const bool &doViscousDamping,
                          const bool &doLinearDamping,
                          const bool &doNonLinearFKExcitation,
                          const bool &doMorrisonElementViscousDrag,
                          const ExcitMethod &excitationMethod,
                          const FreeSurfaceMethod &freeSurfaceMethod,
                          const RadiationMethod &radiationMethod,
                          const hydroForceData &hydroForce,
                          const waveData &waves,
                          const Matrix66d &linearDamping,
                          const Vector31d &cg,
                          const Vector31d &cb,
                          const double &CIdt,
                          const std::vector<Eigen::ArrayXXd> &radForce_IRKB_interp,
                          const double &g,
                          const double &rho,
                          const double &dtFeNonlin,
                          const double &dt,
                          const double &startTime,
                          const bool &bodyToBody,
                          const unsigned int &numBodies,
                          const Array1Nd &CTTime,
                          const Matrix1Nd &lenJ,
                          const double &rampT,
                          const double &mass,
                          const double &dispVol )
    {
#if defined(DEBUG) && defined(MATLAB_MEX_FILE)
mexPrintf ("linearDamping:\n");
MEXSTREAM(linearDamping)
#endif // DEBUG

        _bodyNumber = bodyNumber;
        _doViscousDamping = doViscousDamping;
        _doLinearDamping = doLinearDamping;
        _doNonLinearFKExcitation = doNonLinearFKExcitation;
        _doMorrisonElementViscousDrag = doMorrisonElementViscousDrag;
        _excitationMethod = excitationMethod;
        _freeSurfaceMethod = freeSurfaceMethod;
        _radiationMethod = radiationMethod;
        _hydroForce = hydroForce;
        _waves = waves;
        _linearDamping = linearDamping;
        _cg = cg;
        _cb = cb;
        _mass = mass;
        _dispVol = dispVol;

        // radiation related
        _CIdt = CIdt;
        _radForce_IRKB_interp = radForce_IRKB_interp;

        // simu
        _g = g;
        _rho = rho;
        _dtFeNonlin = dtFeNonlin;
        _dt = dt;
        _startTime = startTime;
        _bodyToBody = bodyToBody;
        _numBodies = numBodies;
        _CTTime = CTTime;
        _lenJ = lenJ;
        _rampT = rampT;

        _accelHist = convArrayHlenNd::Zero (CONV_HIST_LEN,6*_numBodies);
        _velHist = convArrayHlenNd::Zero (CONV_HIST_LEN,6*_numBodies);
        _radForceVelocity = ArrayXXd::Zero (_radForce_IRKB_interp[0].rows (), _radForce_IRKB_interp[0].cols ());
        _radForceOldTime = _startTime;
#if defined(DEBUG) && defined(MATLAB_MEX_FILE)
mexPrintf ("init : _accelHist:\n");
MEXSTREAM(_accelHist)
mexPrintf ("init : _velHist:\n");
MEXSTREAM(_velHist)
mexPrintf ("init : _timeStepHist:\n");
MEXSTREAM(_timeStepHist)
#endif // DEBUG

        _isInititialised = true;

    }

    void hydroBody::hydroForces(const double& t, const Vector61d pos, const Matrix6Nd vel, const Matrix6Nd accel)
    {
#if defined(DEBUG) && defined(MATLAB_MEX_FILE)
mexPrintf ("calling linearExcitationForces\n");
#endif // DEBUG
        // always do linear excitation forces
        linearExcitationForces (t);

        if (_doViscousDamping)
        {
#if defined(DEBUG) && defined(MATLAB_MEX_FILE)
mexPrintf ("calling viscousDampingForces\n");
#endif // DEBUG
            viscousDampingForces (vel);
        }
        else
        {
            F_ViscousDamping = Vector61d::Zero ();
        }

        // linear damping
        if (_doLinearDamping)
        {
#if defined(DEBUG) && defined(MATLAB_MEX_FILE)
mexPrintf ("calling linearDampingForces\n");
#endif // DEBUG
            linearDampingForces (vel);
        }
        else
        {
            F_LinearDamping = Vector61d::Zero ();
        }

#if defined(DEBUG) && defined(MATLAB_MEX_FILE)
mexPrintf ("calling radiationForces\n");
#endif // DEBUG
        // always do radiation forces
        radiationForces (t, vel, accel);

#if defined(DEBUG) && defined(MATLAB_MEX_FILE)
mexPrintf ("calling hydrostaticForces\n");
#endif // DEBUG
        // always do hydrostatic restoring forces
        hydrostaticForces (t, pos);

        if (_doNonLinearFKExcitation)
        {
#if defined(DEBUG) && defined(MATLAB_MEX_FILE)
mexPrintf ("calling nonlinearExcitationForces\n");
#endif // DEBUG
//            if (isempty (waveElv))
//            {
//                // calculate the wave elevation if (it has not already
//                // been done by the hydrostaticForces method
//                waveElevation(t, pos);
//            }
//
//            [ F_ExcitNonLin,
//              WaveNonLinearPressure,
//              WaveLinearPressure ] = nonlinearExcitationForces (obj, t, pos, waveElv);
        }
        else
        {
            F_ExcitNonLin = Vector61d::Zero ();
        }

        if (_doMorrisonElementViscousDrag)
        {
#if defined(DEBUG) && defined(MATLAB_MEX_FILE)
mexPrintf ("calling morrisonElementForce\n");
#endif // DEBUG
//            F_MorrisonElement = morrisonElementForce ( obj, t, pos,
//                                                                 vel(:,_bodyNumber),
//                                                                 accel(:,_bodyNumber) );
        }
        else
        {
            F_MorrisonElement = Vector61d::Zero ();

        }

        F_Excit = F_ExcitLin + F_ExcitNonLin;
#if defined(DEBUG) && defined(MATLAB_MEX_FILE)
mexPrintf ("calling applyRamp\n");
#endif // DEBUG
        applyRamp (t, F_Excit, F_ExcitRamp);

        F_Total = F_ExcitRamp
                 - F_ViscousDamping
                 - F_AddedMass
                 - F_Restoring
                 - F_RadiationDamping
                 - F_MorrisonElement
                 - F_LinearDamping;
#if defined(DEBUG) && defined(MATLAB_MEX_FILE)
mexPrintf ("hydroForces : F_Excit:\n");
MEXSTREAM(F_Excit)
mexPrintf ("hydroForces : F_ExcitRamp:\n");
MEXSTREAM(F_ExcitRamp)
mexPrintf ("hydroForces : F_ViscousDamping:\n");
MEXSTREAM(F_ViscousDamping)
mexPrintf ("hydroForces : F_AddedMass:\n");
MEXSTREAM(F_AddedMass)
mexPrintf ("hydroForces : F_Restoring:\n");
MEXSTREAM(F_Restoring)
mexPrintf ("hydroForces : F_RadiationDamping:\n");
MEXSTREAM(F_RadiationDamping)
mexPrintf ("hydroForces : F_MorrisonElement:\n");
MEXSTREAM(F_MorrisonElement)
mexPrintf ("hydroForces : F_LinearDamping:\n");
MEXSTREAM(F_LinearDamping)
mexPrintf ("hydroForces : F_Total:\n");
MEXSTREAM(F_Total)
mexPrintf ("leaving hydroForces\n");
#endif // DEBUG
    }


    void hydroBody::radiationForces(const double& t, const Matrix6Nd& vel, const Matrix6Nd& accel)
    {
        // matrix multiplication with acceleration
        double delay = 10e-8;
        int naccelvals = 0;

        if (_bodyToBody == true)
        {
            naccelvals = 6 * _numBodies;
        }
        else
        {
            naccelvals = 6;
        }
#if defined(DEBUG) && defined(MATLAB_MEX_FILE)
mexPrintf ("radiationForces 1\n");
#endif // DEBUG
        Array1Nd thisaccel(1, naccelvals);

        thisaccel = Array1Nd::Zero (1, naccelvals);
#if defined(DEBUG) && defined(MATLAB_MEX_FILE)
mexPrintf ("radiationForces 2\n");
#endif // DEBUG
        if (t > (_startTime + delay))
        {
            if (_stepCount > 2)
            {
#if defined(DEBUG) && defined(MATLAB_MEX_FILE)
mexPrintf ("radiationForces 3\n");
#endif // DEBUG
                linearInterp ( _timeStepHist(0,CONV_HIST_LEN-2),
                               _timeStepHist(0,CONV_HIST_LEN-1),
                               _accelHist.row(CONV_HIST_LEN-2),
                               _accelHist.row(CONV_HIST_LEN-1),
                               t - delay,
                               thisaccel );
            }
            else if (_stepCount == 2)
            {
#if defined(DEBUG) && defined(MATLAB_MEX_FILE)
mexPrintf ("radiationForces 4\n");
mexPrintf ("radiationForces : _timeStepHist:\n");
MEXSTREAM(_timeStepHist)
mexPrintf ("radiationForces : _accelHist:\n");
MEXSTREAM(_accelHist)
mexPrintf ("radiationForces : t - delay:\n");
MEXSTREAM(t - delay)
mexPrintf ("radiationForces : _accelHist.row(CONV_HIST_LEN - 2):\n");
MEXSTREAM(_accelHist.row(CONV_HIST_LEN - 2))
mexPrintf ("radiationForces : _accelHist.row(CONV_HIST_LEN - 1):\n");
MEXSTREAM(_accelHist.row(CONV_HIST_LEN - 1))
#endif // DEBUG
//                    thisaccel = interp1 ( obj.timeStepHist(end-_stepCount+1:end)',.
//                                          obj.accelHist(end-_stepCount+1:end,:),
//                                          t - delay,
//                                          'linear', 'extrap' );
                linearInterp ( _timeStepHist(0,CONV_HIST_LEN - 2),
                               _timeStepHist(0,CONV_HIST_LEN - 1),
                               _accelHist.row(CONV_HIST_LEN - 2),
                               _accelHist.row(CONV_HIST_LEN - 1),
                               t - delay,
                               thisaccel );
            }
            else if (_stepCount == 1)
            {
#if defined(DEBUG) && defined(MATLAB_MEX_FILE)
mexPrintf ("radiationForces 5\n");
mexPrintf ("radiationForces : t:\n");
MEXSTREAM(t)
mexPrintf ("radiationForces : _accelHist:\n");
MEXSTREAM(_accelHist)
mexPrintf ("radiationForces : _accelHist.row(CONV_HIST_LEN - 1):\n");
MEXSTREAM(_accelHist.row(CONV_HIST_LEN - 1))
#endif // DEBUG
                thisaccel = _accelHist.row (CONV_HIST_LEN-1);
            }

            //Map<const MatrixN1d> allaccels (thisaccel.data(), thisaccel.size());
            F_AddedMass = _hydroForce.fAddedMass * thisaccel.matrix ().transpose ();
#if defined(DEBUG) && defined(MATLAB_MEX_FILE)
mexPrintf ("radiationForces : _stepCount:\n");
MEXSTREAM(_stepCount)
mexPrintf ("radiationForces : thisaccel.transpose ():\n");
MEXSTREAM(thisaccel.transpose ())
mexPrintf ("radiationForces : _hydroForce.fAddedMass:\n");
MEXSTREAM(_hydroForce.fAddedMass)
mexPrintf ("radiationForces : F_AddedMass:\n");
MEXSTREAM(F_AddedMass)
#endif // DEBUG
        }
        else
        {
            F_AddedMass = Vector61d::Zero ();
        }
#if defined(DEBUG) && defined(MATLAB_MEX_FILE)
mexPrintf ("radiationForces 6\n");
#endif // DEBUG
        switch (_radiationMethod)
        {

            case (RadiationMethod::STATIC_COEFF):
            {
#if defined(DEBUG) && defined(MATLAB_MEX_FILE)
mexPrintf ("radiationForces 7\n");
#endif // DEBUG
                // simple static coefficients
                if (_bodyToBody == true)
                {
#if defined(DEBUG) && defined(MATLAB_MEX_FILE)
mexPrintf ("radiationForces 8\n");
#endif // DEBUG
                    Map<const MatrixN1d> allvels (vel.data(), vel.size());
                    F_RadiationDamping = _hydroForce.fDamping * allvels;
                }
                else
                {
#if defined(DEBUG) && defined(MATLAB_MEX_FILE)
mexPrintf ("radiationForces 9\n");
#endif // DEBUG
                    F_RadiationDamping = _hydroForce.fDamping * vel.col(_bodyNumber);
                }
                break;
            }
            case (RadiationMethod::CONVOLUTION_INTEGRAL):
            {
#if defined(DEBUG) && defined(MATLAB_MEX_FILE)
mexPrintf ("radiationForces 10\n");
#endif // DEBUG
                // convolution
                F_RadiationDamping = radiationConvolutionIntegral (t, vel);
                break;
            }
            case (RadiationMethod::STATESPACE):
            {
#if defined(DEBUG) && defined(MATLAB_MEX_FILE)
mexPrintf ("radiationForces 11\n");
#endif // DEBUG
                // state space
//                     if (_bodyToBody == true)
//                        F_RadiationDamping = _radForceSS.outputs (vel(:));
//                     else
//                         F_RadiationDamping = obj.radForceSS.outputs (vel(:,obj.bodyNumber));
//                     }
                break;
            }
            case (RadiationMethod::EXTERNAL):
            {
#if defined(DEBUG) && defined(MATLAB_MEX_FILE)
mexPrintf ("radiationForces 12\n");
#endif // DEBUG
                // state space, but calculated by some external solver,
                // e.g. by supplying MBDyn with the state-space matrices
                F_RadiationDamping = Vector61d::Zero ();
                break;
            }

        }

    }

    void hydroBody::linearExcitationForces(const double& t)
    {
        switch (_excitationMethod)
        {

            case (ExcitMethod::NO_WAVE):
            {
#if defined(DEBUG) && defined(MATLAB_MEX_FILE)
mexPrintf ("linearExcitationForces 1\n");
#endif // DEBUG
                // no wave
                F_ExcitLin = Vector61d::Zero ();
                break;
            }
            case (ExcitMethod::REGULAR_WAVE):
            {
#if defined(DEBUG) && defined(MATLAB_MEX_FILE)
mexPrintf ("linearExcitationForces 2\n");
#endif // DEBUG
                // regular wave

                // Calculates the wave force, F_wave, for the case of Regular Waves.

                // F_wave =   A * cos(w * t) * Re{Fext}
                //            -  A * sin(w * t) * Im{Fext}

                double wt;
                wt = _waves.w(0,0) * t;
#if defined(DEBUG) && defined(MATLAB_MEX_FILE)
mexPrintf ("linearExcitationForces : wt:\n");
MEXSTREAM(wt)
mexPrintf ("linearExcitationForces : _hydroForce.fExt.re.row(0).array ():\n");
MEXSTREAM(_hydroForce.fExt.re.row(0).array ())
mexPrintf ("linearExcitationForces : _hydroForce.fExt.im.row(0).array ():\n");
MEXSTREAM(_hydroForce.fExt.im.row(0).array ())
#endif // DEBUG
                F_ExcitLin = ( _waves.A(0,0) * (
                                    (cos (wt) * _hydroForce.fExt.re.row(0).array ())
                                    - (sin (wt) * _hydroForce.fExt.im.row(0).array ())
                                                 )
                             ).transpose ().matrix ();
#if defined(DEBUG) && defined(MATLAB_MEX_FILE)
mexPrintf ("linearExcitationForces : F_ExcitLin:\n");
MEXSTREAM(F_ExcitLin)
#endif // DEBUG
                break;
            }
            case (ExcitMethod::IRREGULAR_WAVE):
            {
#if defined(DEBUG) && defined(MATLAB_MEX_FILE)
mexPrintf ("linearExcitationForces 3\n");
#endif // DEBUG
                // irregular wave

                // Calculates the wave force, F_wave, for the case of Irregular Waves.
                //
                // F_wave = sum( F_wave(i))
                //
                // where i = each frequency bin.

                // TODO: check correct dimension/orientation of force output
                ArrayN1d B1 = sin (_waves.w * t + M_PI/2 + _waves.phase);

                ArrayN1d B11 = sin (_waves.w * t + _waves.phase);

                ArrayN1d C1 = sqrt (_waves.A * _waves.dw);

                // _hydroForce.fExt.re will be a 6 * Nfreqs matrix
                MatrixN6d D1 = C1.matrix ().asDiagonal () * _hydroForce.fExt.re;

                MatrixN6d D11 = C1.matrix ().asDiagonal () * _hydroForce.fExt.im;

                MatrixN6d E1 = B1.matrix ().asDiagonal () *  D1;

                MatrixN6d E11 = B11.matrix ().asDiagonal () * D11;

                F_ExcitLin = (E1 - E11).colwise ().sum ().transpose ();

                break;
            }
            case (ExcitMethod::USER_DEFINED_WAVE):
            {
#if defined(DEBUG) && defined(MATLAB_MEX_FILE)
mexPrintf ("linearExcitationForces 4\n");
#endif // DEBUG
                // user defined

                // Calculates the wave force, F_wave, for the case of User Defined Waves.
                //
                // F_wave = convolution calculation [1x6]

                //error ('not yet implemented')
                // TODO: make interpolation function for user defined waves, using ppval (C++ version)
                break;
            }
        }
    }

    void hydroBody::linearDampingForces(const Matrix6Nd& vel)
    {
#if defined(DEBUG) && defined(MATLAB_MEX_FILE)
mexPrintf ("linearDampingForces 1\n");
mexPrintf ("linearDampingForces : _linearDamping:\n");
MEXSTREAM(_linearDamping)
mexPrintf ("linearDampingForces : vel.col (_bodyNumber-1):\n");
MEXSTREAM(vel.col (_bodyNumber-1))
#endif // DEBUG
        F_LinearDamping = _linearDamping * vel.col (_bodyNumber-1);
    }

    void hydroBody::viscousDampingForces(const Matrix6Nd& vel)
    {
#if defined(DEBUG) && defined(MATLAB_MEX_FILE)
mexPrintf ("viscousDampingForces 1\n");
mexPrintf ("viscousDampingForces 1\n");
mexPrintf ("viscousDampingForces : _hydroForce.visDrag:\n");
MEXSTREAM(_hydroForce.visDrag)
mexPrintf ("viscousDampingForces : ( vel.col (_bodyNumber-1).array () * abs (vel.col (_bodyNumber-1).array ()) ).matrix ():\n");
MEXSTREAM(( vel.col (_bodyNumber-1).array () * abs (vel.col (_bodyNumber-1).array ()) ).matrix ())
#endif // DEBUG
        F_ViscousDamping = _hydroForce.visDrag * ( vel.col (_bodyNumber-1).array () * abs (vel.col (_bodyNumber-1).array ()) ).matrix ();
    }

    void hydroBody::morrisonElementForce(const double& t, const Vector61d& pos, const Matrix6Nd& vel, const Matrix6Nd& accel)
    {
        F_MorrisonElement = Vector61d::Zero ();
    }

    void hydroBody::nonlinearExcitationForces(const double& t, const Vector61d& pos, const Matrix6Nd& elv)
    {

    }

    void hydroBody::nonLinearBuoyancy(const double& t, const Vector61d& pos, const Matrix6Nd& elv)
    {
        F_Restoring = Vector61d::Zero ();
    }

    void hydroBody::waveElevationCalc(const double& t, const Vector61d& pos)
    {

//        if (_freeSurfaceMethod == FreeSurfaceMethod::NONLINEAR)
//        {
//            if isempty(_prev_waveElevation)
//            {
//
//                waveElevation = calc_elev (pos,t);
//
//                _prev_waveElevation = waveElevation;
//            }
//            else
//            {
//                if mod(t, _dtFeNonlin) < _dt/2
//                {
//
//                    waveElevation = calc_elev (pos,t);
//
//                    _prev_waveElevation = waveElevation;
//                }
//                else
//                {
//                    waveElevation = _prev_waveElevation;
//                }
//            }
//        }
//        else
//        {
            waveElevation = 0;
//        }

    }

    Vector61d hydroBody::radiationConvolutionIntegral(const double& t, const Matrix6Nd& vel)
    {
        Vector61d F_FM;

        // TODO: we could use advanceStep to ensure histories are
        // updated properly instead of the following test
        // TODO: allow nonuniform time spacing for convolutionIntegral
        if (std::abs (t - _radForceOldTime - _CIdt) < 1e-8)
        {
            // shift the old data one column to the left
            _radForceVelocity.leftCols (_radForceVelocity.cols ()-1) = _radForceVelocity.rightCols (_radForceVelocity.cols ()-1).eval ();

            // replace first column with new data
            Map<const MatrixN1d> allvels (vel.data(), vel.size());
            _radForceVelocity.col(_radForceVelocity.cols ()-1) = allvels;

            // integrate
            std::vector<Eigen::ArrayXXd> time_series;
            for (int i = 0; i < _radForce_IRKB_interp.size (); i++)
            {
                time_series.push_back(_radForce_IRKB_interp[i] * _radForceVelocity);
            }

            for (int i = 0; i < time_series.size (); i++)
            {
                F_FM(i,0) = trapz (_CTTime, time_series[i].rowwise ().sum ());
            }

            _radForceOldF_FM = F_FM;

            _radForceOldTime = t;
        }
        else
        {
            // use the old value which at the start of a simulation
            // is always zeros
            F_FM = _radForceOldF_FM;
        }

        return F_FM;
    }

    void hydroBody::applyRamp(const double& t, const Vector61d& nominal, Vector61d& ramped)
    {
        if (t < _rampT)
        {
#if defined(DEBUG) && defined(MATLAB_MEX_FILE)
mexPrintf ("applyRamp 1\n");
mexPrintf ("applyRamp : t:\n");
MEXSTREAM(t)
mexPrintf ("applyRamp : _rampT:\n");
MEXSTREAM(_rampT)
mexPrintf ("applyRamp : _hydroForce.fExt.im.row(0).array ():\n");
MEXSTREAM(nominal)
mexPrintf ("applyRamp : 0.5 * (1 + std::sin( M_PI * (t / _rampT) + 4.712388980384690)):\n");
MEXSTREAM(0.5 * (1 + std::sin( M_PI * (t / _rampT) + 4.712388980384690)))
mexPrintf ("applyRamp : ramped (before calc):\n");
MEXSTREAM(ramped)
#endif // DEBUG
            // (3 * pi/2) == 4.712388980384690
            ramped = nominal * 0.5 * (1 + std::sin( M_PI * (t / _rampT) + 4.712388980384690));
#if defined(DEBUG) && defined(MATLAB_MEX_FILE)
mexPrintf ("applyRamp : ramped (after calc):\n");
MEXSTREAM(ramped)
#endif // DEBUG
        }
        else
        {
#if defined(DEBUG) && defined(MATLAB_MEX_FILE)
mexPrintf ("applyRamp 2\n");
#endif // DEBUG
            ramped = nominal;
        }
    }

    void hydroBody::hydrostaticForces(const double& t, const Vector61d& pos)
    {

        Vector61d cgfull;
        cgfull << _cg, 0, 0, 0;

        Vector61d bodypos = pos - cgfull;

        double waveElevation = 0;

        switch (_freeSurfaceMethod)
        {

            case (FreeSurfaceMethod::LINEAR):
            {
                // linear hydrostatic restoring force

                //body_hspressure_out = [];

                F_Restoring = _hydroForce.linearHydroRestCoef * bodypos;
#if defined(DEBUG) && defined(MATLAB_MEX_FILE)
mexPrintf ("hydrostaticForces : bodypos:\n");
MEXSTREAM(bodypos)
mexPrintf ("hydrostaticForces : _hydroForce.linearHydroRestCoef * bodypos:\n");
MEXSTREAM(_hydroForce.linearHydroRestCoef * bodypos)
mexPrintf ("hydrostaticForces : F_Restoring (just after _hydroForce.linearHydroRestCoef * bodypos):\n");
MEXSTREAM(F_Restoring)
#endif // DEBUG
                double f_gravity = _g * _mass;

                double f_buoyancy = _rho * _g * _dispVol;

                // Add Net Buoyancy Force to Z-Direction
                F_Restoring(2,0) = F_Restoring(2,0) + (f_gravity - f_buoyancy);
#if defined(DEBUG) && defined(MATLAB_MEX_FILE)
mexPrintf ("hydrostaticForces : F_Restoring (just after F_Restoring(2,0) + (f_gravity - f_buoyancy)):\n");
MEXSTREAM(F_Restoring)
#endif // DEBUG
                Vector31d fbtemp;
                fbtemp << 0, 0, f_buoyancy;
                F_Restoring.block<3,1>(3,0) = F_Restoring.block<3,1>(3,0) + fbtemp.cross (_cb - _cg);
#if defined(DEBUG) && defined(MATLAB_MEX_FILE)
mexPrintf ("hydrostaticForces : _hydroForce.linearHydroRestCoef:\n");
MEXSTREAM(_hydroForce.linearHydroRestCoef)
mexPrintf ("hydrostaticForces : f_gravity:\n");
MEXSTREAM(f_gravity)
mexPrintf ("hydrostaticForces : f_buoyancy:\n");
MEXSTREAM(f_buoyancy)
mexPrintf ("hydrostaticForces : _mass:\n");
MEXSTREAM(_mass)
mexPrintf ("hydrostaticForces : _dispVol:\n");
MEXSTREAM(_dispVol)
mexPrintf ("hydrostaticForces : fbtemp.cross (_cb - _cg):\n");
MEXSTREAM(fbtemp.cross (_cb - _cg))
mexPrintf ("hydrostaticForces : F_Restoring:\n");
MEXSTREAM(F_Restoring)
#endif // DEBUG
                break;
            }
            case (FreeSurfaceMethod::NONLINEAR):
            {

//                waveElevationCalc (t, pos);
//
//                nonLinearBuoyancy (t, pos);
//
//                Vector61d fbtemp;
//                fbtemp << 0, 0, (_simu.g * _mass), 0, 0, 0;
//
//                // Add Net Buoyancy Force to Z-Direction
//                F_Restoring = -F_Restoring + fbtemp;

                break;

            }
        }

    }


    void hydroBody::advanceStep(const double& t, const Matrix6Nd& vel, const Matrix6Nd& accel)
    {
#if defined(DEBUG) && defined(MATLAB_MEX_FILE)
mexPrintf ("advanceStep 1\n");
mexPrintf ("advanceStep : t:\n");
MEXSTREAM(t)
mexPrintf ("advanceStep : _timeStepHist:\n");
MEXSTREAM(_timeStepHist)
mexPrintf ("advanceStep : _timeStepHist.leftCols (_timeStepHist.cols ()-1):\n");
MEXSTREAM(_timeStepHist.leftCols (_timeStepHist.cols ()-1))
mexPrintf ("advanceStep : _timeStepHist.rightCols (_timeStepHist.cols ()-1):\n");
MEXSTREAM(_timeStepHist.rightCols (_timeStepHist.cols ()-1))
#endif // DEBUG
        // update the time history
        _timeStepHist.leftCols (_timeStepHist.cols ()-1) = _timeStepHist.rightCols (_timeStepHist.cols ()-1).eval ();
#if defined(DEBUG) && defined(MATLAB_MEX_FILE)
mexPrintf ("advanceStep 2\n");
mexPrintf ("advanceStep : _timeStepHist:\n");
MEXSTREAM(_timeStepHist)
#endif // DEBUG
        _timeStepHist(0,_timeStepHist.cols()-1) = t;
#if defined(DEBUG) && defined(MATLAB_MEX_FILE)
mexPrintf ("advanceStep 3\n");
mexPrintf ("advanceStep : _timeStepHist:\n");
MEXSTREAM(_timeStepHist)
#endif // DEBUG
        // update the acceleration history
        _accelHist.topRows (_accelHist.rows ()-1) = _accelHist.bottomRows (_accelHist.rows ()-1).eval ();

        Map<const Matrix1Nd> allaccels (accel.data(), accel.size());
        _accelHist.bottomRows (1) = allaccels;

        if (_radiationMethod == RadiationMethod::STATESPACE)
        {
            // update the velocity history
            _velHist.topRows (_velHist.rows ()-1) = _velHist.bottomRows (_velHist.rows ()-1).eval ();

            Map<const Matrix1Nd> allvels (vel.data(), vel.size());
            _velHist.bottomRows (1) = allvels;

        }

        _stepCount = _stepCount + 1;
    }


    void linearInterp ( const Matrix11d x1,
                        const Matrix11d x2,
                        const Matrix1Nd y1,
                        const Matrix1Nd y2,
                        const double u,
                        Matrix1Nd &out )
    {
        Matrix1Nd m = (y2.array() - y1.array()) / (x2.array() - x1.array());
        Matrix1Nd c = y1.array() - m.array() * x1.array();

        out = m.array() * u + c.array();
    }

    void linearInterp ( const double x1,
                        const double x2,
                        const Array1Nd y1,
                        const Array1Nd y2,
                        const double u,
                        Array1Nd &out )
    {
        Array1Nd m = (y2 - y1) / (x2 - x1);
        Array1Nd c = y1 - m * x1;

        out = m * u + c;
    }

    double trapz (const Array1Nd x, const Array1Nd y)
    {
        Array11d two;
        two(0,0) = 2.0;

        Array1Nd z = diff (x) * ( y.segment(0,y.size ()-1) + y.segment(1,y.size ()-1) ) / two;

        return z.sum ();
    }

    double trapz (const ArrayN1d x, const ArrayN1d y)
    {
        Array11d two;
        two(0,0) = 2.0;
        Array1Nd z = diff (x) * ( y.segment(0,y.size ()-1) + y.segment(1,y.size ()-1) ) / two;

        return z.sum ();
    }

    Array1Nd diff (const Array1Nd x)
    {
        return x.segment(0,x.size ()-1) - x.segment(1,x.size ()-1);
    }

    ArrayN1d diff (const ArrayN1d x)
    {
        return x.segment(0,x.size ()-1) - x.segment(1,x.size ()-1);
    }

//    circShiftRight (ArrayXXd &arraytoshift)
//    {
//        Matrix<double, Dynamic, 1> tmp;
//
//        tmp = arraytoshift.col (arraytoshift.cols ()-1);
//
//        arraytoshift.rightCols(_radForceVelocity.cols-1) = arraytoshift.leftCols(_radForceVelocity.cols-1).eval ();
//
//        arraytoshift.col(0) = tmp;
//    }

} // namespace wsim
