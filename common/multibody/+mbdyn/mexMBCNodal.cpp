#include <vector>
#include "mbcxx.h"
#include "mex.h"

// define unique signature for class BEFORE #including class_handle.hpp 
// which  defines a default value if this is not supplied. The class 
// handle signature is used to determine if two wrapped classes are of the 
// same type
#define CLASS_HANDLE_SIGNATURE 0x3CE51CB1
#include "class_handle.hpp"

// A set of utility functions are provided in class_handle.hpp
// in the namespace mexutils. These can be to convert between 
// some matlab and std data types, and ease various tasks for 
// the mex interface
using namespace mexutils;

// Wrapper class to convert matlab arguments to appropriate arguments
// for MBCNodal class
class MBCNodal_wrapper
{
public:
    
    MBCNodal_wrapper (void)
    {
        userefnode = false;
        refnoderot = MBCBase::NONE;
        nodes = 0;
        getlabels = true;
        accelerations = false;
        rot = MBCBase::MAT;
        
        // initialize data structure:
        mbc = new MBCNodal;
    }
    
    ~MBCNodal_wrapper (void)
    {
        // clean up
        mbc->Close ();
        delete mbc;
    }
    
    void Initialize (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
    {   
        std::vector<int> nallowed;
        
        // eight or nine arguments must be supplied : refnode, refnoderot, 
        //   nodes, getlabels, rot, accelerations, commethod, comstring, hostport
        nallowed.push_back (8);
        nallowed.push_back (9);
        int nargin = mxnarginchk (nrhs, nallowed, 2);
        
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
                refnoderot = MBCBase::NONE;
            }
            else if (refnoderotargstr.compare ("theta") == 0)
            {
                refnoderot = MBCBase::THETA;
            }
            else if (refnoderotargstr.compare ("mat") == 0)
            {
                refnoderot = MBCBase::MAT;
            }
            else if (refnoderotargstr.compare ("euler123") == 0)
            {
                refnoderot = MBCBase::EULER_123;
            }
            else
            {
                mexErrMsgIdAndTxt ( "MBCNodal:badrefrotmattype",
                   "Unrecgnised reference node rotation type.");
            }
        }
        else
        {
            // If userefnode is false, ignore the reference node 
            // orientation type string and set to none to indicate the 
            // reference node is not to be used
            refnoderot = MBCBase::NONE;
        }
        
        #ifdef DEBUG
        mexPrintf ("refnoderot: %d\n", userefnode);
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
            rot = MBCBase::NONE;
        }
        else if (rotargstr.compare ("theta") == 0)
        {
            rot = MBCBase::THETA;
        }
        else if (rotargstr.compare ("mat") == 0)
        {
            rot = MBCBase::MAT;
        }
        else if (rotargstr.compare ("euler123") == 0)
        {
            rot = MBCBase::EULER_123;
        }
        else
        {
            mexErrMsgIdAndTxt ( "MBCNodal:badrotmattype",
               "Unrecgnised node rotation type.");
        }
        
        #ifdef DEBUG
        mexPrintf ("rot: %d\n", rot);
        #endif
        
        accelerations = mxnthargscalarbool (nrhs, prhs, 6, 2);
        
        #ifdef DEBUG
        mexPrintf ("getlabels: %d\n", accelerations);
        #endif
        
        std::string commethod = mxnthargstring (nrhs, prhs, 7, 2);
        char* comstring = mxnthargchar (nrhs, prhs, 8, 2);
        
        #ifdef DEBUG
        mexPrintf ("commethod: %s\n", commethod.c_str());
        #endif
        
        #ifdef DEBUG
        mexPrintf ("comstring: %s\n", comstring);
        #endif
        
        if (userefnode && refnoderot == MBCBase::NONE) 
        {
            refnoderot = rot;
        }

//         if (nomoments) {
//             rot = MBCBase::NONE;
//         }
        
        if (mbc->Initialize (refnoderot, nodes, getlabels, rot, accelerations)) 
        {
            mexErrMsgIdAndTxt ("MBCNodal:initfailed: ",
               "MBCNodal::Initialize() failed");
        }

        if  (commethod.compare ("local") == 0)
        {
            /* initialize UNIX socket (path) */
            if (mbc->Init (comstring)) 
            {
                mexErrMsgIdAndTxt ( "MBCNodal:initlocalCommFailure",
                   "Starting local socket communication failed.");
            }
        }
        else if (commethod.compare ("inet") == 0)
        {
            if (nargin != 9)
            {
                mexErrMsgIdAndTxt ( "MBCNodal:noport",
                   "No inet port number was supplied.");
            }
            
            int port = int (mxnthargscalar (nrhs, prhs, 9, 2));
            
            /* initialize INET socket (host, port) */
            if (mbc->Init (comstring, port))
            {
                mexErrMsgIdAndTxt ( "MBCNodal:initInetCommFailure",
                "Starting inet socket communication failed.");
            }
        }
        else
        {
            mexErrMsgIdAndTxt ( "MBCNodal:badCommType",
               "Unrecognised communication method type (should be 'local' or 'inet').");
        }
        
        /* "negotiate" configuration with MBDyn
         * errors out if configurations are inconsistent */
        if (mbc->Negotiate ()) 
        {
            mexErrMsgIdAndTxt ( "MBCNodal:inconsistantConfig",
               "Negotiate call failed, indicating inconsistant configuration, check mbc file options etc. match options used here.");
        }
        
        // free stuff allocated by mx
        mxFree (comstring);
        
        
    }
    
    void GetMotion (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
    {
        std::vector<int> nallowed;
        
        nallowed.push_back (0);
        int nargin = mxnarginchk (nrhs, nallowed, 2);
        
        int status = mbc->GetMotion ();
//         if (status) 
//         {
//             mexErrMsgIdAndTxt ( "MBCNodal:getMotionFailure",
//                "GetMotion call failed.");
//         }
        
        mxSetLHS (status, 1, nlhs, plhs);
    }
    
    void GetNodes (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {
        
        int n = mbc->GetNodes();
        
        mxSetLHS (n, 1, nlhs, plhs);
        
    }
    
    void KinematicsLabel (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {
        if (!getlabels)
        {
            mexErrMsgIdAndTxt ( "MBCNodal:nolabelsavailable",
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
            mexErrMsgIdAndTxt ( "MBCNodal:badNodeRequest",
               "Requested node number is greater than number of nodes.");
        }
        
        if (n < 1)
        {
            mexErrMsgIdAndTxt ( "MBCNodal:badNodeRequest",
               "Requested node number was less than 1.");
        }
        
        #ifdef DEBUG
        mexPrintf ("About to call mbc->KinematicsLabel\n");
        #endif
        
        int label = mbc->KinematicsLabel (n);
        
        mxSetLHS (label, 1, nlhs, plhs);
        
    }
    
    void X (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {
        std::vector<int> nallowed;
        
        // one arg, the node number 
        nallowed.push_back (1);
        
        int nargin = mxnarginchk (nrhs, nallowed, 2);
        
        int n = int (mxnthargscalar (nrhs, prhs, 1, 2));
        
        std::vector<double> X;
        
        X.push_back (mbc->X(n, 1));
        X.push_back (mbc->X(n, 2));
        X.push_back (mbc->X(n, 3));
        
        mxSetLHS (X, 1, nlhs, plhs);
    }
    
    void XP (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {
        std::vector<int> nallowed;
        
        // one arg, the node number 
        nallowed.push_back (1);
        
        int nargin = mxnarginchk (nrhs, nallowed, 2);
        
        int n = int (mxnthargscalar (nrhs, prhs, 1, 2));
        
        std::vector<double> XP;
        
        XP.push_back (mbc->XP(n, 1));
        XP.push_back (mbc->XP(n, 2));
        XP.push_back (mbc->XP(n, 3));
        
        mxSetLHS (XP, 1, nlhs, plhs);
    }
    
    void XPP (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {
        if (accelerations == false)
        {
           mexErrMsgIdAndTxt ("MBCNodal:XPP:nouseaccelerations",
                    "You have set UseAccelerations to false, angular acceleration data is not available.");   
        }
        
        std::vector<int> nallowed;
        
        // one arg, the node number 
        nallowed.push_back (1);
        
        int nargin = mxnarginchk (nrhs, nallowed, 2);
        
        int n = int (mxnthargscalar (nrhs, prhs, 1, 2));
        
        std::vector<double> XPP;
        
        XPP.push_back (mbc->XPP(n, 1));
        XPP.push_back (mbc->XPP(n, 2));
        XPP.push_back (mbc->XPP(n, 3));
        
        mxSetLHS (XPP, 1, nlhs, plhs);
    }
    
    void Theta (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {
        std::vector<int> nallowed;
        
        // one arg, the node number 
        nallowed.push_back (1);
        
        int nargin = mxnarginchk (nrhs, nallowed, 2);
        
        int n = int (mxnthargscalar (nrhs, prhs, 1, 2));
        
        std::vector<double> Theta;
        
        Theta.push_back (mbc->Theta(n, 1));
        Theta.push_back (mbc->Theta(n, 2));
        Theta.push_back (mbc->Theta(n, 3));
        
        mxSetLHS (Theta, 1, nlhs, plhs);
    }
    
    void Omega (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {
        std::vector<int> nallowed;
        
        // one arg, the node number 
        nallowed.push_back (1);
        
        int nargin = mxnarginchk (nrhs, nallowed, 2);
        
        int n = int (mxnthargscalar (nrhs, prhs, 1, 2));
        
        std::vector<double> Omega;
        
        Omega.push_back (mbc->Omega(n, 1));
        Omega.push_back (mbc->Omega(n, 2));
        Omega.push_back (mbc->Omega(n, 3));
        
        mxSetLHS (Omega, 1, nlhs, plhs);
    }
    
    void OmegaP (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {
        
        if (accelerations == false)
        {
           mexErrMsgIdAndTxt ("MBCNodal:OmegaP:nouseaccelerations",
                    "You have set UseAccelerations to false, acceleration data is not available.");   
        }
        
        std::vector<int> nallowed;
        
        // one arg, the node number 
        nallowed.push_back (1);
        
        int nargin = mxnarginchk (nrhs, nallowed, 2);
        
        int n = int (mxnthargscalar (nrhs, prhs, 1, 2));
        
        std::vector<double> OmegaP;
        
        OmegaP.push_back (mbc->OmegaP(n, 1));
        OmegaP.push_back (mbc->OmegaP(n, 2));
        OmegaP.push_back (mbc->OmegaP(n, 3));
        
        mxSetLHS (OmegaP, 1, nlhs, plhs);
    }
    
    
    void F (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {
        std::vector<int> nallowed;
        std::vector<mwSize> fmatind;
        
        // initialise the matrix index vectors
        fmatind.push_back(0);
        fmatind.push_back(0);
                
        nallowed.push_back (1);
        int nargin = mxnarginchk (nrhs, nallowed, 2);
        
        /// \todo: add check for correct force matrix size (3 rows n cols)
        mxNumericArrayWrapper forces = mxnthargmatrix (nrhs, prhs, 1, 2);
        
        if (mbc->GetNodes() > 0)
        {
            if (getlabels)
            {
                for (unsigned n = 1; n <= mbc->GetNodes(); n++)
                {
                    mbc->DynamicsLabel(n) = mbc->KinematicsLabel(n);
                }
            }

            for (unsigned n = 1; n <= mbc->GetNodes(); n++)
            {
                // note mxNumericArrayWrapper indices are zero based not 1-based

                fmatind[0] = (mwSize)0;
                fmatind[1] = (mwSize)(n-1);
                #ifdef DEBUG
                mexPrintf ("fmatind[0]: %d, fmatind[1]: %d, F(%d,1): %f\n", fmatind[0], fmatind[1], n, forces.getDoubleValue (fmatind));
                #endif
                mbc->F(n, 1) = forces.getDoubleValue (fmatind);
                fmatind[0] = (mwSize)1;
                fmatind[1] = (mwSize)(n-1);
                #ifdef DEBUG
                mexPrintf ("fmatind[0]: %d, fmatind[1]: %d, F(%d,2): %f\n", fmatind[0], fmatind[1], n, forces.getDoubleValue (fmatind));
                #endif
                mbc->F(n, 2) = forces.getDoubleValue (fmatind);
                fmatind[0] = (mwSize)2;
                fmatind[1] = (mwSize)(n-1);
                #ifdef DEBUG
                mexPrintf ("fmatind[0]: %d, fmatind[1]: %d, F(%d,3): %f\n", fmatind[0], fmatind[1], n, forces.getDoubleValue (fmatind));
                #endif
                mbc->F(n, 3) = forces.getDoubleValue (fmatind);
            
            }   
        }
    }
    
    void M (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {
        std::vector<int> nallowed;
        std::vector<mwSize> mmatind;
        
        // initialise the matrix index vectors
        mmatind.push_back(0);
        mmatind.push_back(0);
                
        nallowed.push_back (1);
        int nargin = mxnarginchk (nrhs, nallowed, 2);
        
        /// \todo: add check for correct moment matrix size (3 rows n cols)
        mxNumericArrayWrapper moments = mxnthargmatrix (nrhs, prhs, 1, 2);
        
        if (mbc->GetNodes() > 0)
        {
            if (getlabels)
            {
                for (unsigned n = 1; n <= mbc->GetNodes(); n++)
                {
                    mbc->DynamicsLabel(n) = mbc->KinematicsLabel(n);
                }
            }
            
            for (unsigned n = 1; n <= mbc->GetNodes(); n++)
            {
                mmatind[0] = (mwSize)0;
                mmatind[1] = (mwSize)(n-1);
                mbc->M(n, 1) = moments.getDoubleValue (mmatind);
                mmatind[0] = (mwSize)1;
                mmatind[1] = (mwSize)(n-1);
                mbc->M(n, 2) = moments.getDoubleValue (mmatind);
                mmatind[0] = (mwSize)2;
                mmatind[1] = (mwSize)(n-1);
                mbc->M(n, 3) = moments.getDoubleValue (mmatind);
            }    
        }
    }
    
    
    void GetRot (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {
        
        int nnodes = mbc->GetNodes ();
        std::vector<mwSize> index;
        std::vector<int> nallowed;
        
        nallowed.push_back (0);
        int nargin = mxnarginchk (nrhs, nallowed, 2);
        
        // This call gets the rotation matrix data into the mbc data 
        // structure
        MBCBase::Rot rottype = mbc->GetRot ();
        
        switch (rottype) 
        {
                
            case MBCBase::THETA:
            {
                // create the output matrix
                const mwSize dims[] = {3, nnodes};
                plhs[0] = mxCreateNumericArray (2, dims, mxDOUBLE_CLASS, mxREAL);
                
                // wrap it for easy indexing
                mxNumericArrayWrapper thetamat ( plhs[0] );
                // initialise the matrix index vector
                index.push_back (0);
                index.push_back (0);
                
                for (unsigned n = 1; n <= nnodes; n++)
                {
                    
                    index[0] = (mwSize)0;
                    index[1] = (mwSize)(n-1);
                    thetamat.setDoubleValue (index, mbc->Theta (n, 1));

                    index[0] = (mwSize)1;
                    index[1] = (mwSize)(n-1);
                    thetamat.setDoubleValue (index, mbc->Theta (n, 2));

                    index[0] = (mwSize)2;
                    index[1] = (mwSize)(n-1);
                    thetamat.setDoubleValue (index, mbc->Theta (n, 3));
                
                }
                
                return;
            }
            case MBCBase::EULER_123:
            {
                std::vector<double> euler123;
                // create the output matrix
                const mwSize dims[] = {3, nnodes};
                plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
                
                // wrap it for easy indexing
                mxNumericArrayWrapper euler123mat ( plhs[0] );
                // initialise the matrix index vector
                index.push_back(0);
                index.push_back(0);
                
                for (unsigned n = 1; n <= nnodes; n++)
                {
                    index[0] = (mwSize)0;
                    index[1] = (mwSize)(n-1);
                    euler123mat.setDoubleValue (index, mbc->Euler123 (n, 1));

                    index[0] = (mwSize)1;
                    index[1] = (mwSize)(n-1);
                    euler123mat.setDoubleValue (index, mbc->Euler123 (n, 2));

                    index[0] = (mwSize)2;
                    index[1] = (mwSize)(n-1);
                    euler123mat.setDoubleValue (index, mbc->Euler123 (n, 3));
                }
                
                return;
            }  
            default:
            {
                // create the output matrix
                const mwSize dims[] = {3, 3, nnodes};
                plhs[0] = mxCreateNumericArray (3, dims, mxDOUBLE_CLASS, mxREAL);
                
                // wrap it for easy indexing
                mxNumericArrayWrapper rotmat ( plhs[0] );
                
                // initialise the matrix index vector
                index.push_back(0);
                index.push_back(0);
                index.push_back(0);
                
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
                    
                    index[2] = (mwSize)(n-1);
                    
                
                    index[0] = (mwSize)0;
                    index[1] = (mwSize)0;
                    rotmat.setDoubleValue (index, mbc->R(n, 1, 1));

                    index[0] = (mwSize)1;
                    index[1] = (mwSize)0;
                    rotmat.setDoubleValue (index, mbc->R(n, 2, 1));

                    index[0] = (mwSize)2;
                    index[1] = (mwSize)0;
                    rotmat.setDoubleValue (index, mbc->R(n, 3, 1));


                    index[0] = (mwSize)0;
                    index[1] = (mwSize)1;
                    rotmat.setDoubleValue (index, mbc->R(n, 1, 2));

                    index[0] = (mwSize)1;
                    index[1] = (mwSize)1;
                    rotmat.setDoubleValue (index, mbc->R(n, 2, 2));

                    index[0] = (mwSize)2;
                    index[1] = (mwSize)1;
                    rotmat.setDoubleValue (index, mbc->R(n, 3, 2));


                    index[0] = (mwSize)0;
                    index[1] = (mwSize)2;
                    rotmat.setDoubleValue (index, mbc->R(n, 1, 3));

                    index[0] = (mwSize)1;
                    index[1] = (mwSize)2;
                    rotmat.setDoubleValue (index, mbc->R(n, 2, 3));

                    index[0] = (mwSize)2;
                    index[1] = (mwSize)2;
                    rotmat.setDoubleValue (index, mbc->R(n, 3, 3));
                
                }
                
                return;
            }
                
        }
        
    }
    
    void GetRefNodeRot (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {
        std::vector<mwSize> index;
        std::vector<int> nallowed;
        
        if (!userefnode)
        {
            mexErrMsgIdAndTxt ( "MBCNodal:getrefnoderot:norefnode",
               "GetRefNodeRot called, but there is no reference node.");
        }
        
        nallowed.push_back (0);
        int nargin = mxnarginchk (nrhs, nallowed, 2);
        
        // This call gets the rotation matrix data into the mbc data 
        // structure
        MBCBase::Rot rottype = mbc->GetRefNodeRot ();
        
        switch (rottype) 
        {
                
            case MBCBase::THETA:
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
                thetamat.setDoubleValue (index, mbc->Theta (1));

                index[0] = (mwSize)1;
                index[1] = (mwSize)0;
                thetamat.setDoubleValue (index, mbc->Theta (2));

                index[0] = (mwSize)2;
                index[1] = (mwSize)0;
                thetamat.setDoubleValue (index, mbc->Theta (3));
                
                return;
            }
            case MBCBase::EULER_123:
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
                euler123mat.setDoubleValue (index, mbc->Euler123 (1));

                index[0] = (mwSize)1;
                index[1] = (mwSize)0;
                euler123mat.setDoubleValue (index, mbc->Euler123 (2));

                index[0] = (mwSize)2;
                index[1] = (mwSize)0;
                euler123mat.setDoubleValue (index, mbc->Euler123 (3));
                
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
                rotmat.setDoubleValue (index, mbc->R (1, 1));

                index[0] = (mwSize)1;
                index[1] = (mwSize)0;
                rotmat.setDoubleValue (index, mbc->R (2, 1));

                index[0] = (mwSize)2;
                index[1] = (mwSize)0;
                rotmat.setDoubleValue (index, mbc->R (3, 1));


                index[0] = (mwSize)0;
                index[1] = (mwSize)1;
                rotmat.setDoubleValue (index, mbc->R (1, 2));

                index[0] = (mwSize)1;
                index[1] = (mwSize)1;
                rotmat.setDoubleValue (index, mbc->R (2, 2));

                index[0] = (mwSize)2;
                index[1] = (mwSize)1;
                rotmat.setDoubleValue (index, mbc->R (3, 2));


                index[0] = (mwSize)0;
                index[1] = (mwSize)2;
                rotmat.setDoubleValue (index, mbc->R (1, 3));

                index[0] = (mwSize)1;
                index[1] = (mwSize)2;
                rotmat.setDoubleValue (index, mbc->R (2, 3));

                index[0] = (mwSize)2;
                index[1] = (mwSize)2;
                rotmat.setDoubleValue (index, mbc->R (3, 3));
                
                return;
            }
                
        }
        
    }            
            
    void PutForces (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {
        std::vector<int> nallowed;
        
        // one input arg allowed, convergence flag
        nallowed.push_back (1);
        int nargin = mxnarginchk (nrhs, nallowed, 2);
        
        int iterflag = mxnthargscalarbool (nrhs, prhs, 1, 2);
        
        int result = mbc->PutForces (iterflag);
        
        mxSetLHS (result, 1, nlhs, plhs);
    }

private:

    // wrapped  MBCNodal class from mbcxx.h, class is created on the heap
    // using new in the constructor
    MBCNodal *mbc;
    
	bool userefnode;
	MBCBase::Rot refnoderot;
    bool getlabels;
	int nodes;
	bool accelerations;
	MBCBase::Rot rot;

};

// mexfunction defintion, all the interaction with the class is done 
// through this function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{	
     // use the macros provided by class_handle.hpp to create the interface
     // these macros assume the wrapper class can be constructed without 
     // arguments. Note that the instance of the wrapper class is declared 
     // with 'new' and created on the heap
     BEGIN_MEX_CLASS_WRAPPER(MBCNodal_wrapper)
       REGISTER_CLASS_METHOD(MBCNodal_wrapper,Initialize)
       REGISTER_CLASS_METHOD(MBCNodal_wrapper,GetMotion)
       REGISTER_CLASS_METHOD(MBCNodal_wrapper,GetNodes)
       REGISTER_CLASS_METHOD(MBCNodal_wrapper,KinematicsLabel)
       REGISTER_CLASS_METHOD(MBCNodal_wrapper,X)
       REGISTER_CLASS_METHOD(MBCNodal_wrapper,XP)
       REGISTER_CLASS_METHOD(MBCNodal_wrapper,XPP)
       REGISTER_CLASS_METHOD(MBCNodal_wrapper,Theta)
       REGISTER_CLASS_METHOD(MBCNodal_wrapper,Omega)
       REGISTER_CLASS_METHOD(MBCNodal_wrapper,OmegaP)
       REGISTER_CLASS_METHOD(MBCNodal_wrapper,PutForces)
       REGISTER_CLASS_METHOD(MBCNodal_wrapper,F)
       REGISTER_CLASS_METHOD(MBCNodal_wrapper,M)
       REGISTER_CLASS_METHOD(MBCNodal_wrapper,GetRot)
       REGISTER_CLASS_METHOD(MBCNodal_wrapper,GetRefNodeRot)
     END_MEX_CLASS_WRAPPER(MBCNodal_wrapper)
    
}


