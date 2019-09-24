#include <vector>
#include <memory>
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

    MBCNodal_wrapper (void) : mbc (new MBCNodal)
    {
        userefnode = false;
        refnoderot = MBCBase::NONE;
        nodes = 0;
        getlabels = true;
        accelerations = false;
        data_and_next = true;
        verboseflag = false;
        timeout = -1;
        rot = MBCBase::MAT;
        communicationInitialized = false;
    }

    ~MBCNodal_wrapper (void)
    {
        // clean up
        mbc->Close ();
    }

    void Initialize (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {
        std::vector<int> nallowed;

        // eight or nine arguments must be supplied : refnode, refnoderot,
        //   nodes, getlabels, rot, accelerations, data_and_next, verboseflag, timeout, commethod, comstring, hostport
        nallowed.push_back (11);
        nallowed.push_back (12);
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
                refnoderot = MBCBase::NONE;
            }
            else if (refnoderotargstr.compare ("orientation vector") == 0)
            {
                refnoderot = MBCBase::THETA;
            }
            else if (refnoderotargstr.compare ("orientation matrix") == 0)
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
        else if (rotargstr.compare ("orientation vector") == 0)
        {
            rot = MBCBase::THETA;
        }
        else if (rotargstr.compare ("orientation matrix") == 0)
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

        std::string commethod = mxnthargstring (nrhs, prhs, 10, 2);
        char* comstring = mxnthargchar (nrhs, prhs, 11, 2);

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

        mbc->SetDataAndNext (data_and_next);

        mbc->SetVerbose (verboseflag);

        mbc->SetTimeout (timeout);

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
            
            // if we reach here it means communication should be established
            communicationInitialized = true;
        }
        else if (commethod.compare ("inet") == 0)
        {
            if (nargin != 12)
            {
                mexErrMsgIdAndTxt ( "MBCNodal:noport",
                   "No inet port number was supplied.");
            }

            int port = int (mxnthargscalar (nrhs, prhs, 12, 2));

            #ifdef DEBUG
            mexPrintf ("port: %i\n", port);
            #endif

            /* initialize INET socket (host, port) */
            if (mbc->Init (comstring, port))
            {
                mexErrMsgIdAndTxt ( "MBCNodal:initInetCommFailure",
                "Starting inet socket communication failed.");
            }
            
            // if we reach here it means communication should be established
            communicationInitialized = true;
            
        }
        else
        {
            mexErrMsgIdAndTxt ( "MBCNodal:badCommType",
               "Unrecognised communication method type (should be 'local' or 'inet').");
        }

        // free stuff allocated by mx
        mxFree (comstring);


    }
    
    void Negotiate (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {
        std::vector<int> nallowed;
        
        nallowed.push_back (0);
        int nargin = mxnarginchk (nrhs, nallowed, 2);
        
        if (communicationInitialized == true)
        {
        
            /* "negotiate" configuration with MBDyn
             * errors out if configurations are inconsistent */
            if (mbc->Negotiate ())
            {
                mexErrMsgIdAndTxt ( "MBCNodal:inconsistantConfig",
                   "Negotiate call failed, which could indicate inconsistant configuration, check mbc file options etc. match options used here.");
            }
            
        }
        else
        {
            mexErrMsgIdAndTxt ( "MBCNodal:commsnotinitialized",
                   "Not negotiating as communication has not yet been established.");
        
        }
        
    }
    
//     void GetStatus (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
//     {
//         std::vector<int> nallowed;
//
//         nallowed.push_back (0);
//         int nargin = mxnarginchk (nrhs, nallowed, 2);
//
//         int status = mbc->GetStatus ();
//
//         mxSetLHS (status, 1, nlhs, plhs);
//     }

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
        std::vector<mwSize> out_index;
        std::vector<mwSize> in_index;

        // one arg, the node number(s)
        nallowed.push_back (1);

        int nargin = mxnarginchk (nrhs, nallowed, 2);

        //int n = int (mxnthargscalar (nrhs, prhs, 1, 2));
        mxNumericArrayWrapper nodenums = mxnthargmatrix (nrhs, prhs, 1, 2);

        mwSize nrows = nodenums.getRows ();
        mwSize ncols = nodenums.getColumns ();

        if ( (ncols != 1))
        {
            mexErrMsgIdAndTxt("MBCNodal:X:badinputsize",
         "Input node nums should be a column vector of node numbers, but is actually (%d x %d).",
                    nrows, ncols );
        }

        unsigned nnodes = mbc->GetNodes();

        // create the output matrix
        const mwSize dims[] = {3, nrows};
        plhs[0] = mxCreateNumericArray (2, dims, mxDOUBLE_CLASS, mxREAL);

        // wrap it for easy indexing
        mxNumericArrayWrapper outmat ( plhs[0] );

        // initialise the matrix index vectors
        in_index.clear ();
        in_index.push_back (0);
        in_index.push_back (0);
        in_index[1] = (mwSize)0; // always want col 1 of in index

        out_index.clear ();
        out_index.push_back (0);
        out_index.push_back (0);

        for (unsigned n = 1; n <= nrows; n++)
        {
            in_index[0] = (mwSize)n-1;

            int nodenum = int(nodenums.getInt32Value (in_index));

            if ( (nodenum > nnodes))
            {
                mexErrMsgIdAndTxt("MBCNodal:X:badnodenum",
             "Input node num %d is greater than the number of available nodes (%d).",
                        n, nnodes );
            }

            if ( (nodenum < 1))
            {
                mexErrMsgIdAndTxt("MBCNodal:X:badnodenum",
             "Input node num %d is less than 1 (it is %d).",
                        n, nodenum );
            }

            out_index[0] = (mwSize)0;
            out_index[1] = (mwSize)(n-1);
            outmat.setDoubleValue (out_index, mbc->X (nodenum, 1));

            out_index[0] = (mwSize)1;
            out_index[1] = (mwSize)(n-1);
            outmat.setDoubleValue (out_index, mbc->X (nodenum, 2));

            out_index[0] = (mwSize)2;
            out_index[1] = (mwSize)(n-1);
            outmat.setDoubleValue (out_index, mbc->X (nodenum, 3));

        }

    }

    void XP (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {
        std::vector<int> nallowed;
        std::vector<mwSize> out_index;
        std::vector<mwSize> in_index;

        // one arg, the node number(s)
        nallowed.push_back (1);

        int nargin = mxnarginchk (nrhs, nallowed, 2);

        //int n = int (mxnthargscalar (nrhs, prhs, 1, 2));
        mxNumericArrayWrapper nodenums = mxnthargmatrix (nrhs, prhs, 1, 2);

        mwSize nrows = nodenums.getRows ();
        mwSize ncols = nodenums.getColumns ();

        if ( (ncols != 1))
        {
            mexErrMsgIdAndTxt("MBCNodal:XP:badinputsize",
         "Input node nums should be a column vector of node numbers, but is actually (%d x %d).",
                    nrows, ncols );
        }

        unsigned nnodes = mbc->GetNodes();

        // create the output matrix
        const mwSize dims[] = {3, nrows};
        plhs[0] = mxCreateNumericArray (2, dims, mxDOUBLE_CLASS, mxREAL);

        // wrap it for easy indexing
        mxNumericArrayWrapper outmat ( plhs[0] );

        // initialise the matrix index vectors
        in_index.clear ();
        in_index.push_back (0);
        in_index.push_back (0);
        in_index[1] = (mwSize)0; // always want col 1 of in index

        out_index.clear ();
        out_index.push_back (0);
        out_index.push_back (0);

        for (unsigned n = 1; n <= nrows; n++)
        {
            in_index[0] = (mwSize)n-1;

            int nodenum = int(nodenums.getInt32Value (in_index));

            if ( (nodenum > nnodes))
            {
                mexErrMsgIdAndTxt("MBCNodal:XP:badnodenum",
             "Input node num %d is greater than the number of available nodes (%d).",
                        n, nnodes );
            }

            if ( (nodenum < 1))
            {
                mexErrMsgIdAndTxt("MBCNodal:XP:badnodenum",
             "Input node num %d is less than 1 (it is %d).",
                        n, nodenum );
            }

            out_index[0] = (mwSize)0;
            out_index[1] = (mwSize)(n-1);
            outmat.setDoubleValue (out_index, mbc->XP (nodenum, 1));

            out_index[0] = (mwSize)1;
            out_index[1] = (mwSize)(n-1);
            outmat.setDoubleValue (out_index, mbc->XP (nodenum, 2));

            out_index[0] = (mwSize)2;
            out_index[1] = (mwSize)(n-1);
            outmat.setDoubleValue (out_index, mbc->XP (nodenum, 3));

        }
    }

    void XPP (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {
        if (accelerations == false)
        {
           mexErrMsgIdAndTxt ("MBCNodal:XPP:nouseaccelerations",
                    "You have set UseAccelerations to false, angular acceleration data is not available.");
        }

        std::vector<int> nallowed;
        std::vector<mwSize> out_index;
        std::vector<mwSize> in_index;

        // one arg, the node number(s)
        nallowed.push_back (1);

        int nargin = mxnarginchk (nrhs, nallowed, 2);

        //int n = int (mxnthargscalar (nrhs, prhs, 1, 2));
        mxNumericArrayWrapper nodenums = mxnthargmatrix (nrhs, prhs, 1, 2);

        mwSize nrows = nodenums.getRows ();
        mwSize ncols = nodenums.getColumns ();

        if ( (ncols != 1))
        {
            mexErrMsgIdAndTxt("MBCNodal:XPP:badinputsize",
         "Input node nums should be a column vector of node numbers, but is actually (%d x %d).",
                    nrows, ncols );
        }

        unsigned nnodes = mbc->GetNodes();

        // create the output matrix
        const mwSize dims[] = {3, nrows};
        plhs[0] = mxCreateNumericArray (2, dims, mxDOUBLE_CLASS, mxREAL);

        // wrap it for easy indexing
        mxNumericArrayWrapper outmat ( plhs[0] );

        // initialise the matrix index vectors
        in_index.clear ();
        in_index.push_back (0);
        in_index.push_back (0);
        in_index[1] = (mwSize)0; // always want col 1 of in index

        out_index.clear ();
        out_index.push_back (0);
        out_index.push_back (0);

        for (unsigned n = 1; n <= nrows; n++)
        {
            in_index[0] = (mwSize)n-1;

            int nodenum = int(nodenums.getInt32Value (in_index));

            if ( (nodenum > nnodes))
            {
                mexErrMsgIdAndTxt("MBCNodal:XPP:badnodenum",
             "Input node num %d is greater than the number of available nodes (%d).",
                        n, nnodes );
            }

            if ( (nodenum < 1))
            {
                mexErrMsgIdAndTxt("MBCNodal:XPP:badnodenum",
             "Input node num %d is less than 1 (it is %d).",
                        n, nodenum );
            }

            out_index[0] = (mwSize)0;
            out_index[1] = (mwSize)(n-1);
            outmat.setDoubleValue (out_index, mbc->XPP (nodenum, 1));

            out_index[0] = (mwSize)1;
            out_index[1] = (mwSize)(n-1);
            outmat.setDoubleValue (out_index, mbc->XPP (nodenum, 2));

            out_index[0] = (mwSize)2;
            out_index[1] = (mwSize)(n-1);
            outmat.setDoubleValue (out_index, mbc->XPP (nodenum, 3));

        }
    }

    void Theta (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {

        MBCBase::Rot rottype = mbc->GetRot ();

        if (rottype != MBCBase::THETA)
        {
            mexErrMsgIdAndTxt ("MBCNodal:Theta:wrongrottype",
                    "Rotation type is not set to orientation vector, so you cannot use Theta method.");
        }

        std::vector<int> nallowed;
        std::vector<mwSize> out_index;
        std::vector<mwSize> in_index;

        // one arg, the node number(s)
        nallowed.push_back (1);

        int nargin = mxnarginchk (nrhs, nallowed, 2);

        //int n = int (mxnthargscalar (nrhs, prhs, 1, 2));
        mxNumericArrayWrapper nodenums = mxnthargmatrix (nrhs, prhs, 1, 2);

        mwSize nrows = nodenums.getRows ();
        mwSize ncols = nodenums.getColumns ();

        if ( (ncols != 1))
        {
            mexErrMsgIdAndTxt("MBCNodal:Theta:badinputsize",
         "Input node nums should be a column vector of node numbers, but is actually (%d x %d).",
                    nrows, ncols );
        }

        unsigned nnodes = mbc->GetNodes();

        // create the output matrix
        const mwSize dims[] = {3, nrows};
        plhs[0] = mxCreateNumericArray (2, dims, mxDOUBLE_CLASS, mxREAL);

        // wrap it for easy indexing
        mxNumericArrayWrapper outmat ( plhs[0] );

        // initialise the matrix index vectors
        in_index.clear ();
        in_index.push_back (0);
        in_index.push_back (0);
        in_index[1] = (mwSize)0; // always want col 1 of in index

        out_index.clear ();
        out_index.push_back (0);
        out_index.push_back (0);

        for (unsigned n = 1; n <= nrows; n++)
        {
            in_index[0] = (mwSize)n-1;

            int nodenum = int(nodenums.getInt32Value (in_index));

            if ( (nodenum > nnodes))
            {
                mexErrMsgIdAndTxt("MBCNodal:Theta:badnodenum",
             "Input node num %d is greater than the number of available nodes (%d).",
                        n, nnodes );
            }

            if ( (nodenum < 1))
            {
                mexErrMsgIdAndTxt("MBCNodal:Theta:badnodenum",
             "Input node num %d is less than 1 (it is %d).",
                        n, nodenum );
            }

            out_index[0] = (mwSize)0;
            out_index[1] = (mwSize)(n-1);
            outmat.setDoubleValue (out_index, mbc->Theta (nodenum, 1));

            out_index[0] = (mwSize)1;
            out_index[1] = (mwSize)(n-1);
            outmat.setDoubleValue (out_index, mbc->Theta (nodenum, 2));

            out_index[0] = (mwSize)2;
            out_index[1] = (mwSize)(n-1);
            outmat.setDoubleValue (out_index, mbc->Theta (nodenum, 3));

        }
    }

    void Eul (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {

        MBCBase::Rot rottype = mbc->GetRot ();

        if (rottype != MBCBase::EULER_123)
        {
            mexErrMsgIdAndTxt ("MBCNodal:Theta:wrongrottype",
                    "Rotation type is not set to euler123, so you cannot use Eul method.");
        }

        std::vector<int> nallowed;
        std::vector<mwSize> out_index;
        std::vector<mwSize> in_index;

        // one arg, the node number(s)
        nallowed.push_back (1);

        int nargin = mxnarginchk (nrhs, nallowed, 2);

        //int n = int (mxnthargscalar (nrhs, prhs, 1, 2));
        mxNumericArrayWrapper nodenums = mxnthargmatrix (nrhs, prhs, 1, 2);

        mwSize nrows = nodenums.getRows ();
        mwSize ncols = nodenums.getColumns ();

        if ( (ncols != 1))
        {
            mexErrMsgIdAndTxt("MBCNodal:Euler123:badinputsize",
         "Input node nums should be a column vector of node numbers, but is actually (%d x %d).",
                    nrows, ncols );
        }

        unsigned nnodes = mbc->GetNodes();

        // create the output matrix
        const mwSize dims[] = {3, nrows};
        plhs[0] = mxCreateNumericArray (2, dims, mxDOUBLE_CLASS, mxREAL);

        // wrap it for easy indexing
        mxNumericArrayWrapper outmat ( plhs[0] );

        // initialise the matrix index vectors
        in_index.clear ();
        in_index.push_back (0);
        in_index.push_back (0);
        in_index[1] = (mwSize)0; // always want col 1 of in index

        out_index.clear ();
        out_index.push_back (0);
        out_index.push_back (0);

        for (unsigned n = 1; n <= nrows; n++)
        {
            in_index[0] = (mwSize)n-1;

            int nodenum = int(nodenums.getInt32Value (in_index));

            if ( (nodenum > nnodes))
            {
                mexErrMsgIdAndTxt("MBCNodal:Euler123:badnodenum",
             "Input node num %d is greater than the number of available nodes (%d).",
                        n, nnodes );
            }

            if ( (nodenum < 1))
            {
                mexErrMsgIdAndTxt("MBCNodal:Euler123:badnodenum",
             "Input node num %d is less than 1 (it is %d).",
                        n, nodenum );
            }

            out_index[0] = (mwSize)0;
            out_index[1] = (mwSize)(n-1);
            outmat.setDoubleValue (out_index, mbc->Euler123 (nodenum, 1));

            out_index[0] = (mwSize)1;
            out_index[1] = (mwSize)(n-1);
            outmat.setDoubleValue (out_index, mbc->Euler123 (nodenum, 2));

            out_index[0] = (mwSize)2;
            out_index[1] = (mwSize)(n-1);
            outmat.setDoubleValue (out_index, mbc->Euler123 (nodenum, 3));

        }

    }

    void Omega (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {
        std::vector<int> nallowed;
        std::vector<mwSize> out_index;
        std::vector<mwSize> in_index;

        // one arg, the node number(s)
        nallowed.push_back (1);

        int nargin = mxnarginchk (nrhs, nallowed, 2);

        //int n = int (mxnthargscalar (nrhs, prhs, 1, 2));
        mxNumericArrayWrapper nodenums = mxnthargmatrix (nrhs, prhs, 1, 2);

        mwSize nrows = nodenums.getRows ();
        mwSize ncols = nodenums.getColumns ();

        if ( (ncols != 1))
        {
            mexErrMsgIdAndTxt("MBCNodal:Omega:badinputsize",
         "Input node nums should be a column vector of node numbers, but is actually (%d x %d).",
                    nrows, ncols );
        }

        unsigned nnodes = mbc->GetNodes();

        // create the output matrix
        const mwSize dims[] = {3, nrows};
        plhs[0] = mxCreateNumericArray (2, dims, mxDOUBLE_CLASS, mxREAL);

        // wrap it for easy indexing
        mxNumericArrayWrapper outmat ( plhs[0] );

        // initialise the matrix index vectors
        in_index.clear ();
        in_index.push_back (0);
        in_index.push_back (0);
        in_index[1] = (mwSize)0; // always want col 1 of in index

        out_index.clear ();
        out_index.push_back (0);
        out_index.push_back (0);

        for (unsigned n = 1; n <= nrows; n++)
        {
            in_index[0] = (mwSize)n-1;

            int nodenum = int(nodenums.getInt32Value (in_index));

            if ( (nodenum > nnodes))
            {
                mexErrMsgIdAndTxt("MBCNodal:Omega:badnodenum",
             "Input node num %d is greater than the number of available nodes (%d).",
                        n, nnodes );
            }

            if ( (nodenum < 1))
            {
                mexErrMsgIdAndTxt("MBCNodal:Omega:badnodenum",
             "Input node num %d is less than 1 (it is %d).",
                        n, nodenum );
            }

            out_index[0] = (mwSize)0;
            out_index[1] = (mwSize)(n-1);
            outmat.setDoubleValue (out_index, mbc->Omega (nodenum, 1));

            out_index[0] = (mwSize)1;
            out_index[1] = (mwSize)(n-1);
            outmat.setDoubleValue (out_index, mbc->Omega (nodenum, 2));

            out_index[0] = (mwSize)2;
            out_index[1] = (mwSize)(n-1);
            outmat.setDoubleValue (out_index, mbc->Omega (nodenum, 3));

        }
    }

    void OmegaP (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {

        if (accelerations == false)
        {
           mexErrMsgIdAndTxt ("MBCNodal:OmegaP:nouseaccelerations",
                    "You have set UseAccelerations to false, acceleration data is not available.");
        }

        std::vector<int> nallowed;
        std::vector<mwSize> out_index;
        std::vector<mwSize> in_index;

        // one arg, the node number(s)
        nallowed.push_back (1);

        int nargin = mxnarginchk (nrhs, nallowed, 2);

        //int n = int (mxnthargscalar (nrhs, prhs, 1, 2));
        mxNumericArrayWrapper nodenums = mxnthargmatrix (nrhs, prhs, 1, 2);

        mwSize nrows = nodenums.getRows ();
        mwSize ncols = nodenums.getColumns ();

        if ( (ncols != 1))
        {
            mexErrMsgIdAndTxt("MBCNodal:OmegaP:badinputsize",
         "Input node nums should be a column vector of node numbers, but is actually (%d x %d).",
                    nrows, ncols );
        }

        unsigned nnodes = mbc->GetNodes();

        // create the output matrix
        const mwSize dims[] = {3, nrows};
        plhs[0] = mxCreateNumericArray (2, dims, mxDOUBLE_CLASS, mxREAL);

        // wrap it for easy indexing
        mxNumericArrayWrapper outmat ( plhs[0] );

        // initialise the matrix index vectors
        in_index.clear ();
        in_index.push_back (0);
        in_index.push_back (0);
        in_index[1] = (mwSize)0; // always want col 1 of in index

        out_index.clear ();
        out_index.push_back (0);
        out_index.push_back (0);

        for (unsigned n = 1; n <= nrows; n++)
        {
            in_index[0] = (mwSize)n-1;

            int nodenum = int(nodenums.getInt32Value (in_index));

            if ( (nodenum > nnodes))
            {
                mexErrMsgIdAndTxt("MBCNodal:OmegaP:badnodenum",
             "Input node num %d is greater than the number of available nodes (%d).",
                        n, nnodes );
            }

            if ( (nodenum < 1))
            {
                mexErrMsgIdAndTxt("MBCNodal:OmegaP:badnodenum",
             "Input node num %d is less than 1 (it is %d).",
                        n, nodenum );
            }

            out_index[0] = (mwSize)0;
            out_index[1] = (mwSize)(n-1);
            outmat.setDoubleValue (out_index, mbc->OmegaP (nodenum, 1));

            out_index[0] = (mwSize)1;
            out_index[1] = (mwSize)(n-1);
            outmat.setDoubleValue (out_index, mbc->OmegaP (nodenum, 2));

            out_index[0] = (mwSize)2;
            out_index[1] = (mwSize)(n-1);
            outmat.setDoubleValue (out_index, mbc->OmegaP (nodenum, 3));

        }
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

        mxNumericArrayWrapper forces = mxnthargmatrix (nrhs, prhs, 1, 2);

        mwSize nrows = forces.getRows ();
        mwSize ncols = forces.getColumns ();

        unsigned nnodes = mbc->GetNodes();

        if ( (nrows != 3) | (ncols != nnodes) )
        {
            mexErrMsgIdAndTxt("MBCNodal:F:badinputsize",
         "Input force matrix should have 3 rows and Nnodes columns, but is actually (%d x %d).",
                    nrows, ncols );
        }

        if (mbc->GetNodes() > 0)
        {
            if (getlabels)
            {
                for (unsigned n = 1; n <= nnodes; n++)
                {
                    mbc->DynamicsLabel(n) = mbc->KinematicsLabel(n);
                }
            }

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

            #ifdef DEBUG
            mexPrintf ("In 'F', finished setting forces\n");
            #endif
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

        mxNumericArrayWrapper moments = mxnthargmatrix (nrhs, prhs, 1, 2);

        mwSize nrows = moments.getRows ();
        mwSize ncols = moments.getColumns ();

        unsigned nnodes = mbc->GetNodes();

        if ( (nrows != 3) | (ncols != nnodes))
        {
            mexErrMsgIdAndTxt("MBCNodal:M:badinputsize",
         "Input moment matrix should have 3 rows and Nnodes columns, but is actually (%d x %d).",
                    nrows, ncols );
        }

        if (nnodes > 0)
        {
            if (getlabels)
            {
                for (unsigned n = 1; n <= nnodes; n++)
                {
                    mbc->DynamicsLabel(n) = mbc->KinematicsLabel(n);
                }
            }

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
                mbc->M(n, 1) = moments.getDoubleValue (mmatind);
                mmatind[0] = (mwSize)1;
                mmatind[1] = (mwSize)(n-1);
                #ifdef DEBUG
                mexPrintf ("mmatind[0]: %d, mmatind[1]: %d, F(%d,1): %f\n", mmatind[0], mmatind[1], n, moments.getDoubleValue (mmatind));
                #endif
                mbc->M(n, 2) = moments.getDoubleValue (mmatind);
                mmatind[0] = (mwSize)2;
                mmatind[1] = (mwSize)(n-1);
                #ifdef DEBUG
                mexPrintf ("mmatind[0]: %d, mmatind[1]: %d, F(%d,1): %f\n", mmatind[0], mmatind[1], n, moments.getDoubleValue (mmatind));
                #endif
                mbc->M(n, 3) = moments.getDoubleValue (mmatind);
            }

            #ifdef DEBUG
            mexPrintf ("In 'M', finished setting moments\n");
            #endif
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
                #ifdef DEBUG
                mexPrintf ("MBCBase::THETA\n");
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
                #ifdef DEBUG
                mexPrintf ("MBCBase::EULER_123\n");
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
                #ifdef DEBUG
                mexPrintf ("MBCBase default %d (of %d, %d, %d, %d)\n", rottype, MBCBase::NONE, MBCBase::THETA, MBCBase::EULER_123, MBCBase::MAT);
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

        #ifdef DEBUG
        mexPrintf ("PutForces iterflag %d\n", iterflag);
        #endif

        int result = mbc->PutForces (iterflag);

        mxSetLHS (result, 1, nlhs, plhs);
    }

private:

    // wrapped  MBCNodal class from mbcxx.h, class is created on the heap
    // using new in the constructor. Using unique_ptr ensures its
    // destruction when wrapper is done
    std::unique_ptr<MBCNodal> mbc;

	bool userefnode;
	MBCBase::Rot refnoderot;
    bool getlabels;
	int nodes;
	bool accelerations;
    bool data_and_next;
    bool verboseflag;
    int timeout;
	MBCBase::Rot rot;
    
    bool communicationInitialized;
//
//    void Xreturn3ElVal (void (MBCNodal::*fpntr)(const int, const int), int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
//    {
//        std::vector<int> nallowed;
//        std::vector<mwSize> out_index;
//        std::vector<mwSize> in_index;
//
//        // one arg, the node number(s)
//        nallowed.push_back (1);
//
//        int nargin = mxnarginchk (nrhs, nallowed, 2);
//
//        //int n = int (mxnthargscalar (nrhs, prhs, 1, 2));
//        mxNumericArrayWrapper nodenums = mxnthargmatrix (nrhs, prhs, 1, 2);
//
//        mwSize nrows = nodenums.getRows ();
//        mwSize ncols = nodenums.getColumns ();
//
//        if ( (ncols != 1))
//        {
//            mexErrMsgIdAndTxt("MBCNodal:M:badinputsize",
//         "Input node nums should be a column vector of node numbers, but is actually (%d x %d).",
//                    nrows, ncols );
//        }
//
//        unsigned nnodes = mbc->GetNodes();
//
//        // create the output matrix
//        const mwSize dims[] = {3, nrows};
//        plhs[0] = mxCreateNumericArray (2, dims, mxDOUBLE_CLASS, mxREAL);
//
//        // wrap it for easy indexing
//        mxNumericArrayWrapper outmat ( plhs[0] );
//
//        // initialise the matrix index vectors
//        in_index.clear ();
//        in_index.push_back (0);
//        in_index.push_back (0);
//        in_index[1] = (mwSize)0; // always want col 1 of in index
//
//        out_index.clear ();
//        out_index.push_back (0);
//        out_index.push_back (0);
//
//        for (unsigned n = 1; n <= nrows; n++)
//        {
//            in_index[0] = (mwSize)n-1;
//
//            int nodenum = int(nodenums.getInt32Value (in_index));
//
//            if ( (nodenum > nnodes))
//            {
//                mexErrMsgIdAndTxt("MBCNodal:X:badnodenum",
//             "Input node num %d is greater than the number of available nodes (%d).",
//                        n, nnodes );
//            }
//
//            if ( (nodenum < 1))
//            {
//                mexErrMsgIdAndTxt("MBCNodal:X:badnodenum",
//             "Input node num %d is less than 1 (it is %d).",
//                        n, nodenum );
//            }
//
//            out_index[0] = (mwSize)0;
//            out_index[1] = (mwSize)(n-1);
//            outmat.setDoubleValue (out_index, mbc->X (nodenum, 1));
//
//            out_index[0] = (mwSize)1;
//            out_index[1] = (mwSize)(n-1);
//            outmat.setDoubleValue (out_index, mbc->X (nodenum, 2));
//
//            out_index[0] = (mwSize)2;
//            out_index[1] = (mwSize)(n-1);
//            outmat.setDoubleValue (out_index, mbc->X (nodenum, 3));
//
//        }
//
//        //mxSetLHS (outmat, 1, nlhs, plhs);
//    }

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
       REGISTER_CLASS_METHOD(MBCNodal_wrapper,Negotiate)
       REGISTER_CLASS_METHOD(MBCNodal_wrapper,GetMotion)
       //REGISTER_CLASS_METHOD(MBCNodal_wrapper,GetStatus)
       REGISTER_CLASS_METHOD(MBCNodal_wrapper,GetNodes)
       REGISTER_CLASS_METHOD(MBCNodal_wrapper,KinematicsLabel)
       REGISTER_CLASS_METHOD(MBCNodal_wrapper,X)
       REGISTER_CLASS_METHOD(MBCNodal_wrapper,XP)
       REGISTER_CLASS_METHOD(MBCNodal_wrapper,XPP)
       REGISTER_CLASS_METHOD(MBCNodal_wrapper,Theta)
       REGISTER_CLASS_METHOD(MBCNodal_wrapper,Eul)
       REGISTER_CLASS_METHOD(MBCNodal_wrapper,Omega)
       REGISTER_CLASS_METHOD(MBCNodal_wrapper,OmegaP)
       REGISTER_CLASS_METHOD(MBCNodal_wrapper,PutForces)
       REGISTER_CLASS_METHOD(MBCNodal_wrapper,F)
       REGISTER_CLASS_METHOD(MBCNodal_wrapper,M)
       REGISTER_CLASS_METHOD(MBCNodal_wrapper,GetRot)
       REGISTER_CLASS_METHOD(MBCNodal_wrapper,GetRefNodeRot)
     END_MEX_CLASS_WRAPPER(MBCNodal_wrapper)

}


