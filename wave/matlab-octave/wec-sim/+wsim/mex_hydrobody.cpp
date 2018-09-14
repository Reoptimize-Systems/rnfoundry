#include <vector>
#include <memory>
#include "hydrobody.h"
#include "mex.h"

// define unique signature for class BEFORE #including class_handle.hpp
// which  defines a default value if this is not supplied. The class
// handle signature is used to determine if two wrapped classes are of the
// same type
#define CLASS_HANDLE_SIGNATURE 0x3CE12CE8
#include "class_handle.hpp"

// A set of utility functions are provided in class_handle.hpp
// in the namespace mexutils. These can be to convert between
// some matlab and std data types, and ease various tasks for
// the mex interface
using namespace mexutils;
using namespace wsim;

Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>
wrapMxArrayDataInEigenMatrix (const mxArray* in_array)
{

    mxNumericArrayWrapper wrapped_in_array(in_array);

    if (wrapped_in_array.numDimensions () > 2)
    {
        mexErrMsgIdAndTxt(
            "wsim:wrapMxArrayDataInEigen",
            "Given array has > 2 dimensions. Can only create 2-dimensional matrices (and vectors).");
    }

    if (wrapped_in_array.numDimensions() == 1 || wrapped_in_array.numDimensions() == 0)
    {
        mexErrMsgIdAndTxt("wsim:wrapMxArrayDataInEigen", "Given array has 0 or 1 dimensions but we expected a 2-dimensional "
                                        "matrix (or row/column vector).");
        // Even when given a single value dimensionSize() is 2, with n=m=1. When does this happen?
    }

    // We can be sure now that the array is 2-dimensional (or 0, but then we're screwed anyway)
    int nrows = wrapped_in_array.getRows (); // or use array.rows()
    int ncols = wrapped_in_array.getColumns ();

    // get pointer to the underlying data array
    double* data = mxGetPr(in_array);

    Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> > eigen_map(data, nrows, ncols);

    return eigen_map;

}

Eigen::Map<Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>
wrapMxArrayDataInEigenArray (const mxArray* in_array)
{

    mxNumericArrayWrapper wrapped_in_array(in_array);

    if (wrapped_in_array.numDimensions () > 2)
    {
        mexErrMsgIdAndTxt(
            "wsim:wrapMxArrayDataInEigen",
            "Given array has > 2 dimensions. Can only create 2-dimensional matrices (and vectors).");
    }

    if (wrapped_in_array.numDimensions() == 1 || wrapped_in_array.numDimensions() == 0)
    {
        mexErrMsgIdAndTxt("wsim:wrapMxArrayDataInEigen", "Given array has 0 or 1 dimensions but we expected a 2-dimensional "
                                        "matrix (or row/column vector).");
        // Even when given a single value dimensionSize() is 2, with n=m=1. When does this happen?
    }

    // We can be sure now that the array is 2-dimensional (or 0, but then we're screwed anyway)
    int nrows = wrapped_in_array.getRows (); // or use array.rows()
    int ncols = wrapped_in_array.getColumns ();

    // get pointer to the underlying data array
    double* data = mxGetPr(in_array);

    Eigen::Map< Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> > eigen_map(data, nrows, ncols);

    return eigen_map;

}

void
copyEigenDataToMxArrayOuput (const mxArray* in_array, const Eigen::MatrixXd &eigen_mat)
{

    mxNumericArrayWrapper wrapped_in_array(in_array);

    if (wrapped_in_array.numDimensions () > 2)
    {
        mexErrMsgIdAndTxt(
            "wsim:wrapMxArrayDataInEigen",
            "Given array has > 2 dimensions. Can only create 2-dimensional matrices (and vectors).");
    }

    if (wrapped_in_array.numDimensions() == 1 || wrapped_in_array.numDimensions() == 0)
    {
        mexErrMsgIdAndTxt("wsim:wrapMxArrayDataInEigen", "Given array has 0 or 1 dimensions but we expected a 2-dimensional "
                                        "matrix (or row/column vector).");
        // Even when given a single value dimensionSize() is 2, with n=m=1. When does this happen?
    }

    // We can be sure now that the array is 2-dimensional (or 0, but then we're screwed anyway)
    int nrows = wrapped_in_array.getRows (); // or use array.rows()
    int ncols = wrapped_in_array.getColumns ();

    std::vector<mwSize> index;
    index.push_back(0);
    index.push_back(0);
    for (int c = 0; c < ncols; c++)
    {
        for (int r = 0; r < nrows; r++)
        {
            index[0] = r;
            index[1] = c;
            wrapped_in_array.setDoubleValue(index, eigen_mat(r,c));
        }
    }

}

// Wrapper class to convert matlab arguments to appropriate arguments
// for MBCNodal class
class hydroBody_wrapper
{
public:

    hydroBody_wrapper (void) : _hb (new hydroBody)
    {

    }

    ~hydroBody_wrapper (void)
    {
        // clean up
    }

    void Initialize (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {
        std::vector<int> nallowed;


        // eight or nine arguments must be supplied : refnode, refnoderot,
        //   nodes, getlabels, rot, accelerations, data_and_next, verboseflag, timeout, commethod, comstring, hostport
        nallowed.push_back (42);
        int nargin = mxnarginchk (nrhs, nallowed, 2);

        #ifdef DEBUG
        mexPrintf ("nargin: %d\n", nargin);
        #endif

        unsigned int bodyNumber;
        bool doViscousDamping;
        bool doLinearDamping;
        bool doNonLinearFKExcitation;
        bool doMorrisonElementViscousDrag;
        ExcitMethod excitationMethod;
        int excitationMethodNum;
        FreeSurfaceMethod freeSurfaceMethod;
        int freeSurfaceMethodNum;
        RadiationMethod radiationMethod;
        int radiationMethodNum;
        hydroForceData hydroForce;
        waveData waves;
        Matrix66d linearDamping;
        Vector31d cg;
        Vector31d cb;
        double CIdt;
        std::vector<Eigen::ArrayXXd> radForce_IRKB_interp;
        double g;
        double rho;
        double dtFeNonlin;
        double dt;
        double startTime;
        bool bodyToBody;
        unsigned int numBodies;
        Array1Nd CTTime;
        Matrix1Nd lenJ;
        double rampT;
        double mass;
        double dispVol;

        const mxArray* tmpmx;

        bodyNumber = (unsigned int)mxnthargscalar (nrhs, prhs, 1, 2);
        doViscousDamping = mxnthargscalarbool (nrhs, prhs, 2, 2);
        doLinearDamping = mxnthargscalarbool (nrhs, prhs, 3, 2);
        doNonLinearFKExcitation = mxnthargscalarbool (nrhs, prhs, 4, 2);
        doMorrisonElementViscousDrag = mxnthargscalarbool (nrhs, prhs, 5, 2);

        excitationMethodNum = (int)mxnthargscalar (nrhs, prhs, 6, 2);

        switch (excitationMethodNum)
        {
            case (0):
            {
                excitationMethod = ExcitMethod::NO_WAVE;
                break;
            }
            case (1):
            {
                excitationMethod = ExcitMethod::REGULAR_WAVE;
                break;
            }
            case (2):
            {
                excitationMethod = ExcitMethod::IRREGULAR_WAVE;
                break;
            }
            case (3):
            {
                excitationMethod = ExcitMethod::USER_DEFINED_WAVE;
                break;
            }
        }

        freeSurfaceMethodNum = (int)mxnthargscalar (nrhs, prhs, 7, 2);

        switch (freeSurfaceMethodNum)
        {
            case (0):
            {
                freeSurfaceMethod = FreeSurfaceMethod::LINEAR;
                break;
            }
            case (1):
            {
                freeSurfaceMethod = FreeSurfaceMethod::NONLINEAR;
                break;
            }
        }

        radiationMethodNum = (int)mxnthargscalar (nrhs, prhs, 8, 2);

        switch (radiationMethodNum)
        {
            case (0):
            {
                radiationMethod = RadiationMethod::STATIC_COEFF;
                break;
            }
            case (1):
            {
                radiationMethod = RadiationMethod::CONVOLUTION_INTEGRAL;
                break;
            }
            case (2):
            {
                radiationMethod = RadiationMethod::STATESPACE;
                break;
            }
            case (3):
            {
                radiationMethod = RadiationMethod::EXTERNAL;
                break;
            }
        }

        // hydroForce structure

        // fExt re
        tmpmx = mxnthargdoublemxArray (nrhs, prhs, 9, 2);
        hydroForce.fExt.re = wrapMxArrayDataInEigenMatrix (tmpmx);
        // fExt im
        tmpmx = mxnthargdoublemxArray (nrhs, prhs, 10, 2);
        hydroForce.fExt.im = wrapMxArrayDataInEigenMatrix (tmpmx);
        // visDrag
        tmpmx = mxnthargdoublemxArray (nrhs, prhs, 11, 2);
        hydroForce.visDrag = wrapMxArrayDataInEigenMatrix (tmpmx);
        // linearHydroRestCoef
        tmpmx = mxnthargdoublemxArray (nrhs, prhs, 12, 2);
        hydroForce.linearHydroRestCoef = wrapMxArrayDataInEigenMatrix (tmpmx);
        // fDamping
        tmpmx = mxnthargdoublemxArray (nrhs, prhs, 13, 2);
        hydroForce.fDamping = wrapMxArrayDataInEigenMatrix (tmpmx);
        // fAddedMass
        tmpmx = mxnthargdoublemxArray (nrhs, prhs, 14, 2);
        hydroForce.fAddedMass = wrapMxArrayDataInEigenMatrix (tmpmx);

        // waves structure

        // A
        tmpmx = mxnthargdoublemxArray (nrhs, prhs, 15, 2);
        waves.A = wrapMxArrayDataInEigenArray (tmpmx);
        // w
        tmpmx = mxnthargdoublemxArray (nrhs, prhs, 16, 2);
        waves.w = wrapMxArrayDataInEigenArray (tmpmx);
        // phase
        tmpmx = mxnthargdoublemxArray (nrhs, prhs, 17, 2);
        waves.phase = wrapMxArrayDataInEigenArray (tmpmx);
        // dw
        tmpmx = mxnthargdoublemxArray (nrhs, prhs, 18, 2);
        waves.dw = wrapMxArrayDataInEigenArray (tmpmx);
        // k
        tmpmx = mxnthargdoublemxArray (nrhs, prhs, 19, 2);
        waves.k = wrapMxArrayDataInEigenArray (tmpmx);
        // Sf
        tmpmx = mxnthargdoublemxArray (nrhs, prhs, 20, 2);
        waves.Sf = wrapMxArrayDataInEigenArray (tmpmx);

        // linearDamping
        tmpmx = mxnthargdoublemxArray (nrhs, prhs, 21, 2);
        linearDamping = wrapMxArrayDataInEigenMatrix (tmpmx);
        // cg
        tmpmx = mxnthargdoublemxArray (nrhs, prhs, 22, 2);
        cg = wrapMxArrayDataInEigenMatrix (tmpmx);
        // cb
        tmpmx = mxnthargdoublemxArray (nrhs, prhs, 23, 2);
        cb = wrapMxArrayDataInEigenMatrix (tmpmx);

        CIdt = mxnthargscalar(nrhs, prhs, 24, 2);

        // std::vector<Eigen::ArrayXXd> radForce_IRKB_interp; // 25--29
        int r = 0;
        int c = 0;
        for (int i = 0; i < 6; i++)
        {
            // get each slice of the matrix
            tmpmx = mxnthargdoublemxArray (nrhs, prhs, 25+i, 2);

            if (i == 0)
            {
                r = mxGetM (tmpmx);
                c = mxGetN (tmpmx);
            }

            radForce_IRKB_interp.push_back(Eigen::ArrayXXd::Zero (r,c));
            radForce_IRKB_interp[i] = wrapMxArrayDataInEigenArray (tmpmx);
        }

        g = mxnthargscalar(nrhs, prhs, 31, 2);
        rho = mxnthargscalar(nrhs, prhs, 32, 2);
        dtFeNonlin = mxnthargscalar(nrhs, prhs, 33, 2);
        dt = mxnthargscalar(nrhs, prhs, 34, 2);
        startTime = mxnthargscalar(nrhs, prhs, 35, 2);
        bodyToBody = mxnthargscalarbool (nrhs, prhs, 36, 2);
        numBodies = (unsigned int)mxnthargscalar(nrhs, prhs, 37, 2);

        // CTTime
        tmpmx = mxnthargdoublemxArray (nrhs, prhs, 38, 2);
        CTTime = wrapMxArrayDataInEigenArray (tmpmx);

        // lenJ
        tmpmx = mxnthargdoublemxArray (nrhs, prhs, 39, 2);
        lenJ = wrapMxArrayDataInEigenMatrix (tmpmx);

        // mass
        rampT = mxnthargscalar(nrhs, prhs, 40, 2);

        // mass
        mass = mxnthargscalar(nrhs, prhs, 41, 2);

        // dispVol
        dispVol = mxnthargscalar(nrhs, prhs, 42, 2);

        _hb->init ( bodyNumber,
                    doViscousDamping,
                    doLinearDamping,
                    doNonLinearFKExcitation,
                    doMorrisonElementViscousDrag,
                    excitationMethod,
                    freeSurfaceMethod,
                    radiationMethod,
                    hydroForce,
                    waves,
                    linearDamping,
                    cg,
                    cb,
                    CIdt,
                    radForce_IRKB_interp,
                    g,
                    rho,
                    dtFeNonlin,
                    dt,
                    startTime,
                    bodyToBody,
                    numBodies,
                    CTTime,
                    lenJ,
                    rampT,
                    mass,
                    dispVol );

    }


    void hydroForces (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {
        std::vector<int> nallowed;
        const mxArray* tmpmx;

        nallowed.push_back (4);
        int nargin = mxnarginchk (nrhs, nallowed, 2);

        double t = mxnthargscalar(nrhs, prhs, 1, 2);

        tmpmx = mxnthargdoublemxArray (nrhs, prhs, 2, 2);
        Vector61d pos = wrapMxArrayDataInEigenMatrix (tmpmx);

        tmpmx = mxnthargdoublemxArray (nrhs, prhs, 3, 2);
        Matrix6Nd vel = wrapMxArrayDataInEigenMatrix (tmpmx);

        tmpmx = mxnthargdoublemxArray (nrhs, prhs, 4, 2);
        Matrix6Nd accel = wrapMxArrayDataInEigenMatrix (tmpmx);

        // calculate the forces
        _hb->hydroForces ( t, pos, vel, accel );

#if defined(DEBUG) && defined(MATLAB_MEX)
mexPrintf ("just called _hb->hydroForces\n");
#endif // DEBUG


//         if (status)
//         {
//             mexErrMsgIdAndTxt ( "MBCNodal:getMotionFailure",
//                "GetMotion call failed.");
//         }

        if (nlhs > 0)
        {
            mxSetLHS (_hb->F_Total.data (), 1, nlhs, plhs, 6, 1);
        }

#if defined(DEBUG) && defined(MATLAB_MEX)
mexPrintf ("just called mxSetLHS(_hb->F_Total.data () ...\n");
#endif // DEBUG

        if (nlhs > 1)
        {
            const char *fieldnames[] = { "F_ExcitLin",
                                         "F_ViscousDamping",
                                         "F_LinearDamping",
                                         "F_AddedMass",
                                         "F_RadiationDamping",
                                         "F_Restoring",
                                         //"BodyHSPressure",
                                         "F_ExcitNonLin",
                                         //"WaveNonLinearPressure",
                                         //"WaveLinearPressure",
                                         "F_MorrisonElement",
                                         "F_Excit",
                                         "F_ExcitRamp" };

            mxArray* breakdown = mxCreateStructMatrix(1, 1, 10, fieldnames);

            mxArray* field_value_0 = mxCreateDoubleMatrix(6, 1, mxREAL);
            copyEigenDataToMxArrayOuput (field_value_0, _hb->F_ExcitLin);
            mxSetFieldByNumber (breakdown, 0, 0, field_value_0);
#if defined(DEBUG) && defined(MATLAB_MEX)
mexPrintf ("just set field_value_0\n");
#endif // DEBUG
            mxArray* field_value_1 = mxCreateDoubleMatrix(6, 1, mxREAL);
            copyEigenDataToMxArrayOuput (field_value_1, _hb->F_ViscousDamping);
            mxSetFieldByNumber (breakdown, 0, 1, field_value_1);
#if defined(DEBUG) && defined(MATLAB_MEX)
mexPrintf ("just set field_value_1\n");
#endif // DEBUG
            mxArray* field_value_2 = mxCreateDoubleMatrix(6, 1, mxREAL);
            copyEigenDataToMxArrayOuput (field_value_2, _hb->F_LinearDamping);
            mxSetFieldByNumber (breakdown, 0, 2, field_value_2);
#if defined(DEBUG) && defined(MATLAB_MEX)
mexPrintf ("just set field_value_2\n");
#endif // DEBUG
            mxArray* field_value_3 = mxCreateDoubleMatrix(6, 1, mxREAL);
            copyEigenDataToMxArrayOuput (field_value_3, _hb->F_AddedMass);
            mxSetFieldByNumber (breakdown, 0, 3, field_value_3);
#if defined(DEBUG) && defined(MATLAB_MEX)
mexPrintf ("just set field_value_3\n");
#endif // DEBUG
            mxArray* field_value_4 = mxCreateDoubleMatrix(6, 1, mxREAL);
            copyEigenDataToMxArrayOuput (field_value_4, _hb->F_RadiationDamping);
            mxSetFieldByNumber (breakdown, 0, 4, field_value_4);
#if defined(DEBUG) && defined(MATLAB_MEX)
mexPrintf ("just set field_value_4\n");
#endif // DEBUG
            mxArray* field_value_5 = mxCreateDoubleMatrix(6, 1, mxREAL);
            copyEigenDataToMxArrayOuput (field_value_5, _hb->F_Restoring);
            mxSetFieldByNumber (breakdown, 0, 5, field_value_5);
#if defined(DEBUG) && defined(MATLAB_MEX)
mexPrintf ("just set field_value_5\n");
#endif // DEBUG
            mxArray* field_value_6 = mxCreateDoubleMatrix(6, 1, mxREAL);
            copyEigenDataToMxArrayOuput (field_value_6, _hb->F_ExcitNonLin );
            mxSetFieldByNumber (breakdown, 0, 6, field_value_6);
#if defined(DEBUG) && defined(MATLAB_MEX)
mexPrintf ("just set field_value_6\n");
#endif // DEBUG
//            mxArray* field_value = mxCreateDoubleMatrix(6, 1, mxREAL);
//            copyEigenDataToMxArrayOuput (field_value, _hb->WaveNonLinearPressure);
//            mxSetFieldByNumber (breakdown, 0, 6, field_value);
            mxArray* field_value_7 = mxCreateDoubleMatrix(6, 1, mxREAL);
            copyEigenDataToMxArrayOuput (field_value_7, _hb->F_MorrisonElement);
            mxSetFieldByNumber (breakdown, 0, 7, field_value_6);
#if defined(DEBUG) && defined(MATLAB_MEX)
mexPrintf ("just set field_value_7\n");
#endif // DEBUG
            mxArray* field_value_8 = mxCreateDoubleMatrix(6, 1, mxREAL);
            copyEigenDataToMxArrayOuput (field_value_8, _hb->F_Excit);
            mxSetFieldByNumber (breakdown, 0, 8, field_value_7);
#if defined(DEBUG) && defined(MATLAB_MEX)
mexPrintf ("just set field_value_8\n");
#endif // DEBUG
            mxArray* field_value_9 = mxCreateDoubleMatrix(6, 1, mxREAL);
            copyEigenDataToMxArrayOuput (field_value_9, _hb->F_ExcitRamp);
            mxSetFieldByNumber (breakdown, 0, 9, field_value_8);
#if defined(DEBUG) && defined(MATLAB_MEX)
mexPrintf ("just set field_value_9\n");
#endif // DEBUG
            plhs[1] = breakdown;
        }

    }

    void advanceStep (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {
        std::vector<int> nallowed;

        nallowed.push_back (3);
        int nargin = mxnarginchk (nrhs, nallowed, 2);

        double t = mxnthargscalar(nrhs, prhs, 1, 2);

        const mxArray * tmpmx = mxnthargdoublemxArray (nrhs, prhs, 2, 2);
        Matrix6Nd vel = wrapMxArrayDataInEigenMatrix (tmpmx);

        tmpmx = mxnthargdoublemxArray (nrhs, prhs, 3, 2);
        Matrix6Nd accel = wrapMxArrayDataInEigenMatrix (tmpmx);

        _hb->advanceStep (t, vel, accel);

//         if (status)
//         {
//             mexErrMsgIdAndTxt ( "MBCNodal:getMotionFailure",
//                "GetMotion call failed.");
//         }
    }

private:

    // wrapped  hydroBody class. Object is created on the heap
    // using new in the constructor. Using unique_ptr ensures its
    // destruction when wrapper is done
    std::unique_ptr<hydroBody> _hb;

};


// mexfunction defintion, all the interaction with the class is done
// through this function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
     // use the macros provided by class_handle.hpp to create the interface
     // these macros assume the wrapper class can be constructed without
     // arguments. Note that the instance of the wrapper class is declared
     // with 'new' and created on the heap
     BEGIN_MEX_CLASS_WRAPPER(hydroBody_wrapper)
       REGISTER_CLASS_METHOD(hydroBody_wrapper,Initialize)
       REGISTER_CLASS_METHOD(hydroBody_wrapper,hydroForces)
       REGISTER_CLASS_METHOD(hydroBody_wrapper,advanceStep)
     END_MEX_CLASS_WRAPPER(hydroBody_wrapper)

}


