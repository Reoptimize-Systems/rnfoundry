#include <vector>
#include "m_phase_winding.h"
#include "coil.h"
#include "mex.h"
#include "class_handle.hpp"

// A set of utility functions are provided in class_handle.hpp
// in the namespace mexutils. These can be to convert between 
// some matlab and std data types, and ease various tasks for 
// the mex interface
using namespace mexutils;

using namespace Koil;

// 
// 
// // wrapper class to convert matlab arguments to appropriate arguments
// // for wrapped class
// class mPhaseWinding_wrapper
// {
// public:
    void ComputeWinding(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
    {   
        std::vector<int> nallowed;
        mPhaseWinding WL;
        bool SL = false;
        
        nallowed.push_back (4); // four arguments can be supplied
        int nargin = mxnarginchk (nrhs, nallowed);
        
        // get the input arguments
        int m = int (mxnthargscalar (nrhs, prhs, 1, 0));
        int Q = int (mxnthargscalar (nrhs, prhs, 2, 0));
        int p = int (mxnthargscalar (nrhs, prhs, 3, 0));
        int doSingleLayer = int (mxnthargscalar (nrhs, prhs, 4, 0));
        
        if (doSingleLayer == 0)
        {
            SL = false;
        }
        else
        {
            SL = true;
        }
        
        WL.ComputeWinding(m, Q, p, SL);
        
        //mexPrintf ("Computed winding, m: %d Q: %d p: %d SL: %d t: %d\n", m, Q, p, SL, WL.gett ());
        int nwindings = WL.windings.size ();
        std::vector<coil> coils = WL.windings[0].get_coils();
        int ncoilsperwinding = coils.size ();
        
        plhs[0] = mxCreateDoubleMatrix(ncoilsperwinding,2*nwindings,mxREAL);
        
        mxNumericArrayWrapper out = mxNumericArrayWrapper (plhs[0]);
        
        // populate the output matrix
        for (int windind=0; windind<nwindings; windind++)
        {
            coils = WL.windings[windind].get_coils ();
            
            for (int coilind=0; coilind<coils.size (); coilind++)
            {
                //mexPrintf ("%d, %d\n", coils[coilind].s (), coils[coilind].e ());
                
                std::vector<mwSize> index;
                index.push_back (coilind);
                index.push_back (2*windind);
                out.setDoubleValue (index, double (coils[coilind].s ()));
                
                index.clear ();
                index.push_back (coilind);
                index.push_back (2*windind+1);
                out.setDoubleValue (index, double (coils[coilind].e ()));
            }
        }
    }

//     void test(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
//     {
//         // return the value placed in _x
//         mxSetLHS (ex.test (), 1, nlhs, plhs);
//     }
// 
// private:
// 
//     // instance of the wrapped c++ class
//     mPhaseWinding WL;
// 
// };

// mexfunction defintion, all the interction with the class is done through
// this function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{	
//      // use the macros provided by class_handle.hpp to create the interface
//      // these macros assume the wrapper class can be constructed without arguments.
//      // Note that the instance of the wrapper class is declared with 'new' and
//      // created on the heap
//      BEGIN_MEX_CLASS_WRAPPER(mPhaseWinding_wrapper)
//        REGISTER_CLASS_METHOD(mPhaseWinding_wrapper,ComputeWinding)
//        REGISTER_CLASS_METHOD(mPhaseWinding_wrapper,test)
//      END_MEX_CLASS_WRAPPER(mPhaseWinding_wrapper)
    
    ComputeWinding(nlhs, plhs, nrhs, prhs);
    
}


