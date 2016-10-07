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
        labels = true;
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
        //   nodes, labels, rot, accelerations, commethod, comstring, hostport
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
        
        labels = mxnthargscalarbool (nrhs, prhs, 4, 2);
        
        #ifdef DEBUG
        mexPrintf ("labels: %d\n", labels);
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
        mexPrintf ("labels: %d\n", accelerations);
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
        
        if (mbc->Initialize (refnoderot, nodes, labels, rot, accelerations)) 
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
               "Unrecgnised communication method type (should be 'local' or 'inet').");
        }
        
        /* "negotiate" configuration with MBDyn
         * errors out if configurations are inconsistent */
        if (mbc->Negotiate ()) 
        {
            mexErrMsgIdAndTxt ( "MBCNodal:inconsistantConfig",
               "Negotiate call failed, indicating inconsistant mbc configuration.");
        }
        
        // free stuff allocated by mx
        mxFree (comstring);
        
        
    }
    
    void GetMotion (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
    {
        
        int status = mbc->GetMotion ();
        if (status) 
        {
            mexErrMsgIdAndTxt ( "MBCNodal:getMotionFailure",
               "GetMotion call failed.");
        }
        
        mxSetLHS (status, 1, nlhs, plhs);
    }
    
    void GetNodes (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {
        
        int n = mbc->GetNodes();
        
        mxSetLHS (n, 1, nlhs, plhs);
        
    }
    
    void KinematicsLabel (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {
        std::vector<int> nallowed;
        
        // one arg allowed, the node number
        nallowed.push_back (1);
        int nargin = mxnarginchk (nrhs, nallowed, 2);
        
        int nnodes = mbc->GetNodes();
        
        int n = int (mxnthargscalar (nrhs, prhs, 1, 2));
        
        if (n > (nnodes-1))
        {
            mexErrMsgIdAndTxt ( "MBCNodal:badNodeRequest",
               "Requested node number is greater than number of nodes.");
        }
        
        if (n < 0)
        {
            mexErrMsgIdAndTxt ( "MBCNodal:badNodeRequest",
               "Requested node number was negative.");
        }
        
        int label = mbc->KinematicsLabel(n);
        
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
            if (labels)
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
            if (labels)
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
	int nodes;
	bool labels;
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
       REGISTER_CLASS_METHOD(MBCNodal_wrapper,PutForces)
       REGISTER_CLASS_METHOD(MBCNodal_wrapper,F)
     END_MEX_CLASS_WRAPPER(MBCNodal_wrapper)
    
}


