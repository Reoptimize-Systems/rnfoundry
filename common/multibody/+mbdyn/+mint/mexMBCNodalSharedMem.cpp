#include <vector>
#include <memory>
#include "mbcxxshared.h"
#include "mex.h"

// define unique signature for class BEFORE #including class_handle.hpp
// which  defines a default value if this is not supplied. The class
// handle signature is used to determine if two wrapped classes are of the
// same type
#define CLASS_HANDLE_SIGNATURE 0x3CE51CB2
#include "class_handle.hpp"

// A set of utility functions are provided in class_handle.hpp
// in the namespace mexutils. These can be to convert between
// some matlab and std data types, and ease various tasks for
// the mex interface
using namespace mexutils;

// Wrapper class to convert matlab arguments to appropriate arguments
// for MBCSharedMemNodal class
class MBCSharedMemNodal_wrapper
{
public:

    MBCSharedMemNodal_wrapper (void) : mbc (new MBCSharedMemNodal)
    {
        userefnode = false;
        refnoderot = MBC_ROT_NONE;
        nodes = 0;
        getlabels = true;
        accelerations = false;
        data_and_next = true;
        verboseflag = false;
        timeout = -1;
        rot = MBC_ROT_MAT;
    }

    ~MBCSharedMemNodal_wrapper (void)
    {
        // clean up
        #ifdef DEBUG
        mexPrintf ("Calling mbc->Close()\n");
        #endif
        mbc->Close ();
    }

    void Initialize (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {
        std::vector<int> nallowed;

        // Ten arguments must be supplied: 
        // refnode, 
        // refnoderot, 
        // nodes, 
        // getlabels, 
        // rot, 
        // accelerations, 
        // data_and_next, 
        // verboseflag, 
        // timeout, 
        // sharedmemname
        nallowed.push_back (10);
        int nargin = mxnarginchk (nrhs, nallowed, 2);

        #ifdef DEBUG
        mexPrintf ("nargin: %d\n", nargin);
        #endif

        userefnode = mxnthargscalarbool (nrhs, prhs, 1, 2);

        #ifdef DEBUG
        mexPrintf ("userefnode: %d\n", userefnode);
        #endif

        std::string refnoderotargstr = mxnthargstring (nrhs, prhs, 2, 2);

        #ifdef DEBUG
        mexPrintf ("refnoderotargstr: %s\n", refnoderotargstr.c_str());
        #endif

        if (userefnode)
        {
            // if there is supposed to be a reference node, look for the
            // string denoting the reference node orientation type to use
            if  (refnoderotargstr.compare ("none") == 0)
            {
                refnoderot = MBC_ROT_NONE;
            }
            else if (refnoderotargstr.compare ("orientation vector") == 0)
            {
                refnoderot = MBC_ROT_THETA;
            }
            else if (refnoderotargstr.compare ("orientation matrix") == 0)
            {
                refnoderot = MBC_ROT_MAT;
            }
            else if (refnoderotargstr.compare ("euler 123") == 0)
            {
                refnoderot = MBC_ROT_EULER_123;
            }
            else
            {
                mexErrMsgIdAndTxt ( "MBCSharedMemNodal:badrefrotmattype",
                   "Unrecgnised reference node rotation type.");
            }
        }
        else
        {
            // If userefnode is false, ignore the reference node
            // orientation type string and set to none to indicate the
            // reference node is not to be used
            refnoderot = MBC_ROT_NONE;
        }

        #ifdef DEBUG
        mexPrintf ("refnoderot: %d\n", refnoderot);
        #endif

        nodes = int (mxnthargscalar (nrhs, prhs, 3, 2));

        #ifdef DEBUG
        mexPrintf ("nodes: %d\n", nodes);
        #endif

        getlabels = mxnthargscalarbool (nrhs, prhs, 4, 2);

        #ifdef DEBUG
        mexPrintf ("getlabels: %d\n", getlabels);
        #endif

        std::string rotargstr = mxnthargstring (nrhs, prhs, 5, 2);

        #ifdef DEBUG
        mexPrintf ("rotargstr: %s\n", rotargstr.c_str());
        #endif

        if  (rotargstr.compare ("none") == 0)
        {
            rot = MBC_ROT_NONE;
        }
        else if (rotargstr.compare ("orientation vector") == 0)
        {
            rot = MBC_ROT_THETA;
        }
        else if (rotargstr.compare ("orientation matrix") == 0)
        {
            rot = MBC_ROT_MAT;
        }
        else if (rotargstr.compare ("euler 123") == 0)
        {
            rot = MBC_ROT_EULER_123;
        }
        else
        {
            mexErrMsgIdAndTxt ( "MBCSharedMemNodal:badrotmattype",
               "Unrecgnised node rotation type.");
        }
        
//         if (userefnode == false)
//         {
//             refnoderot = rot;
//             
//             #ifdef DEBUG
//             mexPrintf ("refnoderot changed to: %d\n", refnoderot);
//             #endif
//         }

        #ifdef DEBUG
        mexPrintf ("rot: %d\n", rot);
        #endif

        accelerations = mxnthargscalarbool (nrhs, prhs, 6, 2);

        #ifdef DEBUG
        mexPrintf ("accelerations: %d\n", accelerations);
        #endif

        data_and_next = mxnthargscalarbool (nrhs, prhs, 7, 2);

        #ifdef DEBUG
        mexPrintf ("data_and_next: %d\n", data_and_next);
        #endif

        verboseflag = mxnthargscalarbool (nrhs, prhs, 8, 2);

        #ifdef DEBUG
        mexPrintf ("verboseflag: %d\n", verboseflag);
        #endif

        timeout = int (mxnthargscalar (nrhs, prhs, 9, 2));

        #ifdef DEBUG
        mexPrintf ("timeout: %d\n", timeout);
        #endif

        std::string sharedmemname = mxnthargstring (nrhs, prhs, 10, 2);

        #ifdef DEBUG
        mexPrintf ("sharedmemname: %s\n", sharedmemname.c_str());
        #endif

        if (userefnode && refnoderot == MBC_ROT_NONE)
        {
            refnoderot = rot;
            
            #ifdef DEBUG
            mexPrintf ("refnoderot changed to: %d\n", refnoderot);
            #endif
        }

//         if (nomoments) {
//             rot = MBC_ROT_NONE;
//         }

        mbc->SetDataAndNext (data_and_next);

        mbc->SetVerbose (verboseflag);

        mbc->SetTimeout (timeout);
        
        #ifdef DEBUG
        mexPrintf ("refnoderot now: %d\n", refnoderot);
        #endif

        if (mbc->Initialize (refnoderot, nodes, getlabels, rot, accelerations, sharedmemname))
        {
            mexErrMsgIdAndTxt ("MBCSharedMemNodal:initfailed: ",
               "MBCSharedMemNodal::Initialize() failed");
        }

        mbc->Init ();

        /* "negotiate" configuration with MBDyn
         * errors out if configurations are inconsistent */
        if (mbc->Negotiate ())
        {
            mexErrMsgIdAndTxt ( "MBCSharedMemNodal:inconsistantConfig",
               "Negotiate call failed, indicating inconsistant configuration, check mbc file options etc. match options used here.");
        }

    }

    void GetStatus (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {
        std::vector<int> nallowed;

        nallowed.push_back (0);
        int nargin = mxnarginchk (nrhs, nallowed, 2);

        int status = mbc->GetStatus ();

        mxSetLHS (status, 1, nlhs, plhs);
    }

    void GetMotion (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {
        checkStatus ();

        std::vector<int> nallowed;

        nallowed.push_back (0);
        int nargin = mxnarginchk (nrhs, nallowed, 2);

        int status = mbc->GetMotion ();
//         if (status)
//         {
//             mexErrMsgIdAndTxt ( "MBCSharedMemNodal:getMotionFailure",
//                "GetMotion call failed.");
//         }

        mxSetLHS (status, 1, nlhs, plhs);
    }

    void GetNodes (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {
        checkStatus ();

        int n = mbc->GetNodes();

        mxSetLHS (n, 1, nlhs, plhs);

    }

    void KinematicsLabel (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {
        checkStatus ();

        if (!getlabels)
        {
            mexErrMsgIdAndTxt ( "MBCSharedMemNodal:nolabelsavailable",
               "You did not request labels at initialisation so they cannot be retreived.");
        }

        std::vector<int> nallowed;

        // one arg allowed, the node number
        nallowed.push_back (1);
        int nargin = mxnarginchk (nrhs, nallowed, 2);

        int nnodes = mbc->GetNodes();

        #ifdef DEBUG
        mexPrintf ("in KinematicsLabel, nnodes: %d\n", nnodes);
        #endif

        int n = int (mxnthargscalar (nrhs, prhs, 1, 2));

        #ifdef DEBUG
        mexPrintf ("in KinematicsLabel, n: %d\n", n);
        #endif

        if (n > nnodes)
        {
            mexErrMsgIdAndTxt ( "MBCSharedMemNodal:badNodeRequest",
               "Requested node number is greater than number of nodes.");
        }

        if (n < 1)
        {
            mexErrMsgIdAndTxt ( "MBCSharedMemNodal:badNodeRequest",
               "Requested node number was less than 1.");
        }

        #ifdef DEBUG
        mexPrintf ("About to call mbc->KinematicsLabel\n");
        #endif

        int label = mbc->GetKinematicsLabel (uint32_t(n));

        mxSetLHS (label, 1, nlhs, plhs);

    }

    void X (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {
        checkStatus ();

        std::vector<int> nallowed;

        // one arg, the node number
        nallowed.push_back (1);

        int nargin = mxnarginchk (nrhs, nallowed, 2);

        int n = int (mxnthargscalar (nrhs, prhs, 1, 2));

        std::vector<double> X;

        X.push_back (mbc->GetX(n, 1));
        X.push_back (mbc->GetX(n, 2));
        X.push_back (mbc->GetX(n, 3));

        mxSetLHS (X, 1, nlhs, plhs);
    }

    void XP (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {
        checkStatus ();

        std::vector<int> nallowed;

        // one arg, the node number
        nallowed.push_back (1);

        int nargin = mxnarginchk (nrhs, nallowed, 2);

        int n = int (mxnthargscalar (nrhs, prhs, 1, 2));

        std::vector<double> XP;

        XP.push_back (mbc->GetXP(n, 1));
        XP.push_back (mbc->GetXP(n, 2));
        XP.push_back (mbc->GetXP(n, 3));

        mxSetLHS (XP, 1, nlhs, plhs);
    }

    void XPP (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {
        checkStatus ();

        if (accelerations == false)
        {
           mexErrMsgIdAndTxt ("MBCSharedMemNodal:XPP:nouseaccelerations",
                    "You have set UseAccelerations to false, angular acceleration data is not available.");
        }

        std::vector<int> nallowed;

        // one arg, the node number
        nallowed.push_back (1);

        int nargin = mxnarginchk (nrhs, nallowed, 2);

        int n = int (mxnthargscalar (nrhs, prhs, 1, 2));

        std::vector<double> XPP;

        XPP.push_back (mbc->GetXPP(n, 1));
        XPP.push_back (mbc->GetXPP(n, 2));
        XPP.push_back (mbc->GetXPP(n, 3));

        mxSetLHS (XPP, 1, nlhs, plhs);
    }

    void Theta (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {
        checkStatus ();

        MBCType rottype = mbc->GetRot ();

        if (rottype != MBC_ROT_THETA)
        {
            mexErrMsgIdAndTxt ("MBCSharedMemNodal:Theta:wrongrottype",
                    "Rotation type is not set to orientation vector, so you cannot use Theta method.");
        }

        std::vector<int> nallowed;

        // one arg, the node number
        nallowed.push_back (1);

        int nargin = mxnarginchk (nrhs, nallowed, 2);

        int n = int (mxnthargscalar (nrhs, prhs, 1, 2));

        std::vector<double> Theta;

        Theta.push_back (mbc->GetTheta(n, 1));
        Theta.push_back (mbc->GetTheta(n, 2));
        Theta.push_back (mbc->GetTheta(n, 3));

        mxSetLHS (Theta, 1, nlhs, plhs);
    }

    void Omega (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {
        checkStatus ();

        std::vector<int> nallowed;

        // one arg, the node number
        nallowed.push_back (1);

        int nargin = mxnarginchk (nrhs, nallowed, 2);

        int n = int (mxnthargscalar (nrhs, prhs, 1, 2));

        std::vector<double> Omega;

        Omega.push_back (mbc->GetOmega(n, 1));
        Omega.push_back (mbc->GetOmega(n, 2));
        Omega.push_back (mbc->GetOmega(n, 3));

        mxSetLHS (Omega, 1, nlhs, plhs);
    }

    void OmegaP (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {
        checkStatus ();

        if (accelerations == false)
        {
           mexErrMsgIdAndTxt ("MBCSharedMemNodal:OmegaP:nouseaccelerations",
                    "You have set UseAccelerations to false, acceleration data is not available.");
        }

        std::vector<int> nallowed;

        // one arg, the node number
        nallowed.push_back (1);

        int nargin = mxnarginchk (nrhs, nallowed, 2);

        int n = int (mxnthargscalar (nrhs, prhs, 1, 2));

        std::vector<double> OmegaP;

        OmegaP.push_back (mbc->GetOmegaP(n, 1));
        OmegaP.push_back (mbc->GetOmegaP(n, 2));
        OmegaP.push_back (mbc->GetOmegaP(n, 3));

        mxSetLHS (OmegaP, 1, nlhs, plhs);
    }


    void F (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {
        checkStatus ();

        std::vector<int> nallowed;
        std::vector<mwSize> fmatind;

        // initialise the matrix index vectors
        fmatind.push_back(0);
        fmatind.push_back(0);

        nallowed.push_back (1);
        int nargin = mxnarginchk (nrhs, nallowed, 2);

        mxNumericArrayWrapper forces = mxnthargmatrix (nrhs, prhs, 1, 2);

        mwSize nrows = forces.getRows ();
        mwSize ncols = forces.getColumns ();

        unsigned nnodes = mbc->GetNodes();

        if ( (nrows != 3) | (ncols != nnodes) )
        {
            mexErrMsgIdAndTxt("MBCSharedMemNodal:F:badinputsize",
         "Input force matrix should have 3 rows and Nnodes columns, but is actually (%d x %d).",
                    nrows, ncols );
        }

        if (mbc->GetNodes() > 0)
        {
//            if (getlabels)
//            {
//                for (unsigned n = 1; n <= nnodes; n++)
//                {
//                    mbc->DynamicsLabel(n) = mbc->KinematicsLabel(n);
//                }
//            }

            #ifdef DEBUG
            mexPrintf ("In 'F', about to set forces\n");
            #endif

            for (unsigned n = 1; n <= nnodes; n++)
            {
                // note mxNumericArrayWrapper indices are zero based not 1-based

                fmatind[0] = (mwSize)0;
                fmatind[1] = (mwSize)(n-1);
                #ifdef DEBUG
                mexPrintf ("fmatind[0]: %d, fmatind[1]: %d, F(%d,1): %f\n", fmatind[0], fmatind[1], n, forces.getDoubleValue (fmatind));
                #endif
                mbc->SetF(n, 1, forces.getDoubleValue (fmatind));

                fmatind[0] = (mwSize)1;
                fmatind[1] = (mwSize)(n-1);
                #ifdef DEBUG
                mexPrintf ("fmatind[0]: %d, fmatind[1]: %d, F(%d,2): %f\n", fmatind[0], fmatind[1], n, forces.getDoubleValue (fmatind));
                #endif
                mbc->SetF(n, 2, forces.getDoubleValue (fmatind));

                fmatind[0] = (mwSize)2;
                fmatind[1] = (mwSize)(n-1);
                #ifdef DEBUG
                mexPrintf ("fmatind[0]: %d, fmatind[1]: %d, F(%d,3): %f\n", fmatind[0], fmatind[1], n, forces.getDoubleValue (fmatind));
                #endif
                mbc->SetF(n, 3, forces.getDoubleValue (fmatind));

            }

            #ifdef DEBUG
            mexPrintf ("In 'F', finished setting forces\n");
            #endif
        }
    }

    void M (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {
        checkStatus ();

        std::vector<int> nallowed;
        std::vector<mwSize> mmatind;

        // initialise the matrix index vectors
        mmatind.push_back(0);
        mmatind.push_back(0);

        nallowed.push_back (1);
        int nargin = mxnarginchk (nrhs, nallowed, 2);

        mxNumericArrayWrapper moments = mxnthargmatrix (nrhs, prhs, 1, 2);

        mwSize nrows = moments.getRows ();
        mwSize ncols = moments.getColumns ();

        unsigned nnodes = mbc->GetNodes();

        if ( (nrows != 3) | (ncols != nnodes))
        {
            mexErrMsgIdAndTxt("MBCSharedMemNodal:M:badinputsize",
         "Input moment matrix should have 3 rows and Nnodes columns, but is actually (%d x %d).",
                    nrows, ncols );
        }

        if (nnodes > 0)
        {
//            if (getlabels)
//            {
//                for (unsigned n = 1; n <= nnodes; n++)
//                {
//                    mbc->DynamicsLabel(n) = mbc->KinematicsLabel(n);
//                }
//            }

            #ifdef DEBUG
            mexPrintf ("In 'M', about to set moments\n");
            #endif

            for (unsigned n = 1; n <= nnodes; n++)
            {
                mmatind[0] = (mwSize)0;
                mmatind[1] = (mwSize)(n-1);
                #ifdef DEBUG
                mexPrintf ("mmatind[0]: %d, mmatind[1]: %d, F(%d,1): %f\n", mmatind[0], mmatind[1], n, moments.getDoubleValue (mmatind));
                #endif
                mbc->SetM(n, 1, moments.getDoubleValue (mmatind));
                mmatind[0] = (mwSize)1;
                mmatind[1] = (mwSize)(n-1);
                #ifdef DEBUG
                mexPrintf ("mmatind[0]: %d, mmatind[1]: %d, F(%d,1): %f\n", mmatind[0], mmatind[1], n, moments.getDoubleValue (mmatind));
                #endif
                mbc->SetM(n, 2, moments.getDoubleValue (mmatind));
                mmatind[0] = (mwSize)2;
                mmatind[1] = (mwSize)(n-1);
                #ifdef DEBUG
                mexPrintf ("mmatind[0]: %d, mmatind[1]: %d, F(%d,1): %f\n", mmatind[0], mmatind[1], n, moments.getDoubleValue (mmatind));
                #endif
                mbc->SetM(n, 3, moments.getDoubleValue (mmatind));
            }

            #ifdef DEBUG
            mexPrintf ("In 'M', finished setting moments\n");
            #endif
        }
    }


    void GetRot (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {
        checkStatus ();

        int nnodes = mbc->GetNodes ();
        std::vector<mwSize> index;
        std::vector<int> nallowed;

        nallowed.push_back (0);
        int nargin = mxnarginchk (nrhs, nallowed, 2);

        // This call gets the rotation matrix data into the mbc data
        // structure
        MBCType rottype = mbc->GetRot ();

        switch (rottype)
        {

            case MBC_ROT_THETA:
            {
                #ifdef DEBUG
                mexPrintf ("MBC_ROT_THETA\n");
                #endif

                // create the output matrix
                const mwSize dims[] = {3, nnodes};
                plhs[0] = mxCreateNumericArray (2, dims, mxDOUBLE_CLASS, mxREAL);

                // wrap it for easy indexing
                mxNumericArrayWrapper thetamat ( plhs[0] );
                // initialise the matrix index vector
                index.clear ();
                index.push_back (0);
                index.push_back (0);

                for (unsigned n = 1; n <= nnodes; n++)
                {

                    index[0] = (mwSize)0;
                    index[1] = (mwSize)(n-1);
                    thetamat.setDoubleValue (index, mbc->GetTheta (n, 1));

                    index[0] = (mwSize)1;
                    index[1] = (mwSize)(n-1);
                    thetamat.setDoubleValue (index, mbc->GetTheta (n, 2));

                    index[0] = (mwSize)2;
                    index[1] = (mwSize)(n-1);
                    thetamat.setDoubleValue (index, mbc->GetTheta (n, 3));

                }

                return;
            }
            case MBC_ROT_EULER_123:
            {
                #ifdef DEBUG
                mexPrintf ("MBC_ROT_EULER_123\n");
                #endif

                std::vector<double> euler123;
                // create the output matrix
                const mwSize dims[] = {3, nnodes};
                plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);

                // wrap it for easy indexing
                mxNumericArrayWrapper euler123mat ( plhs[0] );
                // initialise the matrix index vector
                index.clear ();
                index.push_back(0);
                index.push_back(0);

                for (unsigned n = 1; n <= nnodes; n++)
                {
                    index[0] = (mwSize)0;
                    index[1] = (mwSize)(n-1);
                    euler123mat.setDoubleValue (index, mbc->GetEuler123 (n, 1));

                    index[0] = (mwSize)1;
                    index[1] = (mwSize)(n-1);
                    euler123mat.setDoubleValue (index, mbc->GetEuler123 (n, 2));

                    index[0] = (mwSize)2;
                    index[1] = (mwSize)(n-1);
                    euler123mat.setDoubleValue (index, mbc->GetEuler123 (n, 3));
                }

                return;
            }
            default:
            {
                #ifdef DEBUG
                mexPrintf ("MBCSharedMemBase default %d (of %d, %d, %d, %d)\n", rottype, MBC_ROT_NONE, MBC_ROT_THETA, MBC_ROT_EULER_123, MBC_ROT_MAT);
                #endif

                // create the output matrix
                const mwSize dims[] = {3, 3, nnodes};
                plhs[0] = mxCreateNumericArray (3, dims, mxDOUBLE_CLASS, mxREAL);

                // wrap it for easy indexing
                mxNumericArrayWrapper rotmat ( plhs[0] );

                // initialise the matrix index vector
                index.clear ();
                index.push_back (0);
                index.push_back (0);
                if (nnodes > 1)
                {
                    // matlab collapses array dimensions of size 1
                    index.push_back (0);
                }

                // orientation matrix supplied as three column vectors:
                //
                //  r11,r21,r31,r12,r22,r32,r13,r23,r33;
                //
                // mbc->R(n, 1, 1)  mbc->R(n, 1, 2)  mbc->R(n, 1, 3)
                // mbc->R(n, 2, 1)  mbc->R(n, 2, 2)  mbc->R(n, 2, 3)
                // mbc->R(n, 3, 1)  mbc->R(n, 3, 2)  mbc->R(n, 3, 3)
                //

                for (unsigned n = 1; n <= nnodes; n++)
                {

                    if (nnodes > 1)
                    {
                        // matlab collapses array dimensions of size 1
                        index[2] = (mwSize)(n-1);
                    }


                    index[0] = (mwSize)0;
                    index[1] = (mwSize)0;
                    rotmat.setDoubleValue (index, mbc->GetR(n, 1, 1));

                    index[0] = (mwSize)1;
                    index[1] = (mwSize)0;
                    rotmat.setDoubleValue (index, mbc->GetR(n, 2, 1));

                    index[0] = (mwSize)2;
                    index[1] = (mwSize)0;
                    rotmat.setDoubleValue (index, mbc->GetR(n, 3, 1));


                    index[0] = (mwSize)0;
                    index[1] = (mwSize)1;
                    rotmat.setDoubleValue (index, mbc->GetR(n, 1, 2));

                    index[0] = (mwSize)1;
                    index[1] = (mwSize)1;
                    rotmat.setDoubleValue (index, mbc->GetR(n, 2, 2));

                    index[0] = (mwSize)2;
                    index[1] = (mwSize)1;
                    rotmat.setDoubleValue (index, mbc->GetR(n, 3, 2));


                    index[0] = (mwSize)0;
                    index[1] = (mwSize)2;
                    rotmat.setDoubleValue (index, mbc->GetR(n, 1, 3));

                    index[0] = (mwSize)1;
                    index[1] = (mwSize)2;
                    rotmat.setDoubleValue (index, mbc->GetR(n, 2, 3));

                    index[0] = (mwSize)2;
                    index[1] = (mwSize)2;
                    rotmat.setDoubleValue (index, mbc->GetR(n, 3, 3));

                }

                return;
            }

        }

    }

    void GetRefNodeRot (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {
        checkStatus ();

        std::vector<mwSize> index;
        std::vector<int> nallowed;

        if (!userefnode)
        {
            mexErrMsgIdAndTxt ( "MBCSharedMemNodal:getrefnoderot:norefnode",
               "GetRefNodeRot called, but there is no reference node.");
        }

        nallowed.push_back (0);
        int nargin = mxnarginchk (nrhs, nallowed, 2);

        // This call gets the rotation matrix data into the mbc data
        // structure
        MBCType rottype = mbc->GetRefNodeRot ();

        switch (rottype)
        {

            case MBC_ROT_THETA:
            {
                // create the output matrix
                const mwSize dims[] = {3, 1};
                plhs[0] = mxCreateNumericArray (2, dims, mxDOUBLE_CLASS, mxREAL);

                // wrap it for easy indexing
                mxNumericArrayWrapper thetamat ( plhs[0] );
                // initialise the matrix index vector
                index.push_back (0);
                index.push_back (0);

                index[0] = (mwSize)0;
                index[1] = (mwSize)0;
                thetamat.setDoubleValue (index, mbc->GetRefNodeTheta (1));

                index[0] = (mwSize)1;
                index[1] = (mwSize)0;
                thetamat.setDoubleValue (index, mbc->GetRefNodeTheta (2));

                index[0] = (mwSize)2;
                index[1] = (mwSize)0;
                thetamat.setDoubleValue (index, mbc->GetRefNodeTheta (3));

                return;
            }
            case MBC_ROT_EULER_123:
            {
                std::vector<double> euler123;
                // create the output matrix
                const mwSize dims[] = {3, 1};
                plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);

                // wrap it for easy indexing
                mxNumericArrayWrapper euler123mat ( plhs[0] );
                // initialise the matrix index vector
                index.push_back(0);
                index.push_back(0);

                index[0] = (mwSize)0;
                index[1] = (mwSize)0;
                euler123mat.setDoubleValue (index, mbc->GetRefNodeEuler123 (1));

                index[0] = (mwSize)1;
                index[1] = (mwSize)0;
                euler123mat.setDoubleValue (index, mbc->GetRefNodeEuler123 (2));

                index[0] = (mwSize)2;
                index[1] = (mwSize)0;
                euler123mat.setDoubleValue (index, mbc->GetRefNodeEuler123 (3));

                return;
            }

            default:
            {

                // create the output matrix
                const mwSize dims[] = {3, 3};
                plhs[0] = mxCreateNumericArray (2, dims, mxDOUBLE_CLASS, mxREAL);

                // wrap it for easy indexing
                mxNumericArrayWrapper rotmat ( plhs[0] );

                // initialise the matrix index vector
                index.push_back(0);
                index.push_back(0);

                // orientation matrix supplied as three column vectors:
                //
                //  r11,r21,r31,r12,r22,r32,r13,r23,r33;
                //
                // mbc->R(1, 1)  mbc->R(1, 2)  mbc->R(1, 3)
                // mbc->R(2, 1)  mbc->R(2, 2)  mbc->R(2, 3)
                // mbc->R(3, 1)  mbc->R(3, 2)  mbc->R(3, 3)
                //

                index[0] = (mwSize)0;
                index[1] = (mwSize)0;
                rotmat.setDoubleValue (index, mbc->GetRefNodeR (1, 1));

                index[0] = (mwSize)1;
                index[1] = (mwSize)0;
                rotmat.setDoubleValue (index, mbc->GetRefNodeR (2, 1));

                index[0] = (mwSize)2;
                index[1] = (mwSize)0;
                rotmat.setDoubleValue (index, mbc->GetRefNodeR (3, 1));


                index[0] = (mwSize)0;
                index[1] = (mwSize)1;
                rotmat.setDoubleValue (index, mbc->GetRefNodeR (1, 2));

                index[0] = (mwSize)1;
                index[1] = (mwSize)1;
                rotmat.setDoubleValue (index, mbc->GetRefNodeR (2, 2));

                index[0] = (mwSize)2;
                index[1] = (mwSize)1;
                rotmat.setDoubleValue (index, mbc->GetRefNodeR (3, 2));


                index[0] = (mwSize)0;
                index[1] = (mwSize)2;
                rotmat.setDoubleValue (index, mbc->GetRefNodeR (1, 3));

                index[0] = (mwSize)1;
                index[1] = (mwSize)2;
                rotmat.setDoubleValue (index, mbc->GetRefNodeR (2, 3));

                index[0] = (mwSize)2;
                index[1] = (mwSize)2;
                rotmat.setDoubleValue (index, mbc->GetRefNodeR (3, 3));

                return;
            }

        }

    }

    void PutForces (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {

        checkStatus ();

        std::vector<int> nallowed;

        // one input arg allowed, convergence flag
        nallowed.push_back (1);
        int nargin = mxnarginchk (nrhs, nallowed, 2);

        int iterflag = mxnthargscalarbool (nrhs, prhs, 1, 2);

        #ifdef DEBUG
        mexPrintf ("PutForces iterflag %d\n", iterflag);
        #endif

        int result = mbc->PutForces (iterflag);

        mxSetLHS (result, 1, nlhs, plhs);
    }

private:

    // wrapped  MBCSharedMemNodal class from mbcxxshared.h, class is created on the heap
    // using new in the constructor. Using unique_ptr ensures its
    // destruction when wrapper is done
    std::unique_ptr<MBCSharedMemNodal> mbc;

	bool userefnode;
	MBCType refnoderot;
    bool getlabels;
	int nodes;
	bool accelerations;
    bool data_and_next;
    bool verboseflag;
    int timeout;
	MBCType rot;

	void checkStatus (void)
	{

        switch (mbc->GetStatus ())
        {

            case MBCSharedMemNodal::FINISHED:
            {
                mexErrMsgIdAndTxt ( "MBCSharedMemNodal:mbdynfinished",
                   "MBDyn simulation is finished.");

                break;
            }
            case MBCSharedMemNodal::NOT_READY:
            {
                mexErrMsgIdAndTxt ( "MBCSharedMemNodal:mbdynfinished",
                   "Simulation is not ready (have you called Initialize?).");

                break;
            }
            default:
            {
                // do nothing
            }
        }

	}

};

// mexfunction defintion, all the interaction with the class is done
// through this function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
     // use the macros provided by class_handle.hpp to create the interface
     // these macros assume the wrapper class can be constructed without
     // arguments. Note that the instance of the wrapper class is declared
     // with 'new' and created on the heap
     BEGIN_MEX_CLASS_WRAPPER(MBCSharedMemNodal_wrapper)
       REGISTER_CLASS_METHOD(MBCSharedMemNodal_wrapper,Initialize)
       REGISTER_CLASS_METHOD(MBCSharedMemNodal_wrapper,GetMotion)
       REGISTER_CLASS_METHOD(MBCSharedMemNodal_wrapper,GetStatus)
       REGISTER_CLASS_METHOD(MBCSharedMemNodal_wrapper,GetNodes)
       REGISTER_CLASS_METHOD(MBCSharedMemNodal_wrapper,KinematicsLabel)
       REGISTER_CLASS_METHOD(MBCSharedMemNodal_wrapper,X)
       REGISTER_CLASS_METHOD(MBCSharedMemNodal_wrapper,XP)
       REGISTER_CLASS_METHOD(MBCSharedMemNodal_wrapper,XPP)
       REGISTER_CLASS_METHOD(MBCSharedMemNodal_wrapper,Theta)
       REGISTER_CLASS_METHOD(MBCSharedMemNodal_wrapper,Omega)
       REGISTER_CLASS_METHOD(MBCSharedMemNodal_wrapper,OmegaP)
       REGISTER_CLASS_METHOD(MBCSharedMemNodal_wrapper,PutForces)
       REGISTER_CLASS_METHOD(MBCSharedMemNodal_wrapper,F)
       REGISTER_CLASS_METHOD(MBCSharedMemNodal_wrapper,M)
       REGISTER_CLASS_METHOD(MBCSharedMemNodal_wrapper,GetRot)
       REGISTER_CLASS_METHOD(MBCSharedMemNodal_wrapper,GetRefNodeRot)
     END_MEX_CLASS_WRAPPER(MBCSharedMemNodal_wrapper)

}


