#ifndef HYDROBODY_H
#define HYDROBODY_H

#include <vector>
#include <eigen3/Eigen/Dense>

#define EIGEN_NO_AUTOMATIC_RESIZING
#define CONV_HIST_LEN 3

typedef Eigen::Matrix<double, 6, Eigen::Dynamic> Matrix6Nd;
typedef Eigen::Array<double, CONV_HIST_LEN, Eigen::Dynamic> convArrayHlenNd;
typedef Eigen::Array<double, 1, CONV_HIST_LEN> convArray1Hlend;
typedef Eigen::Matrix<double, 1, Eigen::Dynamic> Matrix1Nd;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> MatrixN1d;
typedef Eigen::Matrix<double, Eigen::Dynamic, 6> MatrixN6d;
typedef Eigen::Matrix<double, 6, 6> Matrix66d;
typedef Eigen::Matrix<double, 1, 1> Matrix11d;
typedef Eigen::Matrix<double, 6, 1> Vector61d;
typedef Eigen::Matrix<double, 3, 1> Vector31d;
typedef Eigen::Matrix<double, 1, 3> Vector13d;
typedef Eigen::Array<double, 1, 1> Array11d;
typedef Eigen::Array<double, Eigen::Dynamic, 1> ArrayN1d;
typedef Eigen::Array<double, 1, Eigen::Dynamic> Array1Nd;
typedef Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> ArrayXXd;

namespace wsim
{
    enum class ExcitMethod { NO_WAVE,
                             REGULAR_WAVE,
                             IRREGULAR_WAVE,
                             USER_DEFINED_WAVE };

    enum class FreeSurfaceMethod { LINEAR,
                                   NONLINEAR };

    enum class RadiationMethod { STATIC_COEFF,
                                 CONVOLUTION_INTEGRAL,
                                 STATESPACE,
                                 EXTERNAL };

    struct excitFData {
        MatrixN6d re;
        MatrixN6d im;
    };

    struct hydroForceData {

        excitFData fExt;
        Matrix66d visDrag;
        Matrix66d linearHydroRestCoef;
        Eigen::MatrixXd fDamping;
        Eigen::MatrixXd fAddedMass;

    };

    struct waveData {

        ArrayN1d A;
        ArrayN1d w;
        ArrayN1d phase;
        ArrayN1d dw;
        ArrayN1d k;
        ArrayN1d Sf;

    };

    class hydroBody
    {
        public:

            hydroBody();
            virtual ~hydroBody();

            void init ( const unsigned int &bodyNumber,
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
                        const double &dispVol );

            void hydroForces ( const double &t,
                               const Vector61d pos,
                               const Matrix6Nd vel,
                               const Matrix6Nd accel );

            void advanceStep ( const double &t,
                               const Matrix6Nd &vel,
                               const Matrix6Nd &accel );

            Vector61d F_Total;
            Vector61d F_ExcitLin;
            Vector61d F_Excit;
            Vector61d F_ExcitRamp;
            Vector61d F_LinearDamping;
            Vector61d F_ViscousDamping;
            Vector61d F_AddedMass;
            Vector61d F_RadiationDamping;
            Vector61d F_Restoring;
            //Eigen::Matrix6Nd BodyHSPressure;
            Vector61d F_ExcitNonLin;
            //Matrix6Nd WaveNonLinearPressure;
            //Eigen::MatrixXd WaveLinearPressure;
            Vector61d F_MorrisonElement;
            double waveElevation;

        protected:

            void radiationForces ( const double &t,
                                   const Matrix6Nd &vel,
                                   const Matrix6Nd &accel );

            void linearExcitationForces ( const double &t );

            void linearDampingForces ( const Matrix6Nd &vel );

            void viscousDampingForces ( const Matrix6Nd &vel );

            void morrisonElementForce ( const double &t,
                                        const Vector61d &pos,
                                        const Matrix6Nd &vel,
                                        const Matrix6Nd &accel );

            void nonlinearExcitationForces ( const double &t,
                                             const Vector61d &pos,
                                             const Matrix6Nd &elv );

            void nonLinearBuoyancy ( const double &t,
                                     const Vector61d &pos,
                                     const Matrix6Nd &elv );

            void hydrostaticForces ( const double &t,
                                     const Vector61d &pos );

//            void triWaveElev ( const double &t,
//                            const Matrix6Nd &center );

            Vector61d radiationConvolutionIntegral( const double &t,
                                                    const Matrix6Nd &vel );

            void applyRamp ( const double &t,
                             const Vector61d &nominal,
                             Vector61d &ramped );

            void waveElevationCalc ( const double &t,
                                     const Vector61d &pos );



        private:

            bool _isInititialised;

            unsigned int _bodyNumber;
            bool _doViscousDamping;
            bool _doLinearDamping;
            bool _doNonLinearFKExcitation;
            bool _doMorrisonElementViscousDrag;
            ExcitMethod _excitationMethod;
            FreeSurfaceMethod _freeSurfaceMethod;
            RadiationMethod _radiationMethod;
            hydroForceData _hydroForce;
            waveData _waves;
            Matrix66d _linearDamping;
            Vector31d _cg;
            Vector31d _cb;
            double _mass;
            double _dispVol;

            // radiation related
            unsigned int _stepCount;
            convArray1Hlend _timeStepHist;
            convArrayHlenNd _accelHist;
            convArrayHlenNd _velHist;
            ArrayXXd _radForceVelocity;
            double _radForceOldTime;
            Vector61d _radForceOldF_FM;
            double _CIdt;
            std::vector<Eigen::ArrayXXd> _radForce_IRKB_interp;

            // simu
            double _g;
            double _rho;
            double _prev_waveElevation;
            double _dtFeNonlin;
            double _dt;
            double _startTime;
            bool _bodyToBody;
            unsigned int _numBodies;
            Array1Nd _CTTime;
            Matrix1Nd _lenJ;
            double _rampT;


    };



    void linearInterp ( const Matrix11d x1,
                        const Matrix11d x2,
                        const Matrix1Nd y1,
                        const Matrix1Nd y2,
                        const double u,
                        Matrix1Nd &out );

    void linearInterp ( const double x1,
                        const double x2,
                        const Array1Nd y1,
                        const Array1Nd y2,
                        const double u,
                        Array1Nd &out );

    double trapz (const Array1Nd x, const Array1Nd y);

    double trapz (const ArrayN1d x, const ArrayN1d y);

    Array1Nd diff (const Array1Nd x);

    ArrayN1d diff (const ArrayN1d x);

} // namespace wsim

#endif // HYDROBODY_H
