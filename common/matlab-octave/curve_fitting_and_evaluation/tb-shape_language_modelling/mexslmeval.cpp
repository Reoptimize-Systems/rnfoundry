#ifdef __cplusplus
extern "C"
{
#endif
    #include <gsl/gsl_histogram.h>
#ifdef __cplusplus
}
#endif  

#include <vector>
#include <cmath>
#include "mex.h"

/*
 * mexslmeval.cpp - mex version of slmeval.m
 *
 * This is a MEX-file for MATLAB.
 *
 */

double getcoef (const double * coef, const int nbins, const int bin, const int col)
{
    return *(coef + ((nbins+1)*(col-1) + bin));
}
        
void slmeval3dx0 (gsl_histogram * hist,
                  const int nbins, 
                  const double * bins, 
                  const double * coef, 
                  const int ndata, 
                  const double * x, 
                  const std::vector<double> dx,
                  double * y)
{
    size_t xbin = 0;
    size_t i = 0;

    for (i = 0; i < ndata; i++)
    {
        if (*(x+i) <= bins[0])
        {
            xbin = 0;
        }
        else if (*(x+i) >= bins[nbins-1])
        {
            xbin = nbins - 1;
        }
        else
        {
            gsl_histogram_find (hist, *(x+i), &xbin);
        }

        // f(x)
        double t = (*(x+i) - *(bins+xbin)) / dx[xbin];
        double t2 = t*t;
        double t3 = t*t*t;
        double s2 = std::pow((1.0-t),2);
        double s3 = std::pow((1.0-t),3);
                
        *(y+i) = (-getcoef(coef, nbins, xbin, 2) * (s3-s2) + getcoef(coef, nbins, xbin+1, 2) * (t3-t2)) * dx[xbin] 
                 + getcoef(coef, nbins, xbin, 1) * (3.0 * s2 - 2.0* s3)
                 + getcoef(coef, nbins, xbin+1, 1) * (3.0 * t2 - 2.0 * t3);

    }

}

void slmeval3dx1 (gsl_histogram * hist,
                  const int nbins, 
                  const double * bins, 
                  const double * coef, 
                  const int ndata, 
                  const double * x, 
                  const std::vector<double> dx,
                  double * y)
{
    size_t xbin = 0;
    size_t i = 0;

    for (i = 0; i < ndata; i++)
    {
        if (*(x+i) <= bins[0])
        {
            xbin = 0;
        }
        else if (*(x+i) >= bins[nbins-1])
        {
            xbin = nbins - 1;
        }
        else
        {
            gsl_histogram_find (hist, *(x+i), &xbin);
        }

        // first derivative for the cubic case
        double t = (*(x+i) - *(bins+xbin)) / dx[xbin];
        double t2 = t*t;
        double s = 1-t;
        double s2 = std::pow((1-t),2);
        *(y+i) = -getcoef(coef, nbins, xbin,2)*(-3*s2+2*s) +
                 getcoef(coef, nbins, xbin+1,2)*(3*t2-2*t) +
                 (getcoef(coef, nbins, xbin,1)*(-6*s+6*s2) +
                 getcoef(coef, nbins, xbin+1,1)*(6*t-6*t2)) / dx[xbin];

    }

}

void slmeval3dx2 (gsl_histogram * hist,
                  const int nbins, 
                  const double * bins, 
                  const double * coef, 
                  const int ndata, 
                  const double * x, 
                  const std::vector<double> dx,
                  double * y)
{
    size_t xbin = 0;
    size_t i = 0;

    for (i = 0; i < ndata; i++)
    {
        if (*(x+i) <= bins[0])
        {
            xbin = 0;
        }
        else if (*(x+i) >= bins[nbins-1])
        {
            xbin = nbins - 1;
        }
        else
        {
            gsl_histogram_find (hist, *(x+i), &xbin);
        }

        // second derivative of a cubic
        double t = (*(x+i) - *(bins+xbin)) / dx[xbin];

        double s = 1-t;
        
        *(y+i) = (-getcoef(coef, nbins, xbin,2)*(6*s - 2) +
                 getcoef(coef, nbins, xbin+1,2)*(6*t - 2))/dx[xbin] +
                 (getcoef(coef, nbins, xbin,1)*(6 - 12*s) +
                 getcoef(coef, nbins, xbin+1,1)*(6 - 12*t))/(std::pow(dx[xbin],2));

    }

}

void slmeval3dx3 (gsl_histogram * hist,
                  const int nbins, 
                  const double * bins, 
                  const double * coef, 
                  const int ndata, 
                  const double * x, 
                  const std::vector<double> dx,
                  double * y)
{
    size_t xbin = 0;
    size_t i = 0;

    for (i = 0; i < ndata; i++)
    {
        if (*(x+i) <= bins[0])
        {
            xbin = 0;
        }
        else if (*(x+i) >= bins[nbins-1])
        {
            xbin = nbins - 1;
        }
        else
        {
            gsl_histogram_find (hist, *(x+i), &xbin);
        }

        // second derivative of a cubic
        double t = (*(x+i) - *(bins+xbin)) / dx[xbin];

        double s = 1-t;
        
        *(y+i) = 6*(getcoef(coef, nbins, xbin,2) + getcoef(coef, nbins, xbin+1,2))/(std::pow(dx[xbin],2)) +
                        12*(getcoef(coef, nbins, xbin,1) - getcoef(coef, nbins, xbin+1,1))/(std::pow(dx[xbin],3));

    }

}


int slmeval3 (const int derivorder,
               const int nbins, 
               const double * bins, 
               const double* coef, 
               const int ndata, 
               const double * x, 
               double * y)
{
    size_t xbin = 0;
    size_t i = 0;
    std::vector<double> dx;
    int ret = 0;
    
    // resuerve enough space for the dx data in advance to avoid  
    // reallocating memory as we fill it
    dx.reserve (nbins-1);
    // now calculate dx
    for (i = 0; i < nbins; i++)
    {
        dx.push_back (bins[i+1] - bins[i]);
    }
    
    // create a histogram for binning the data
    gsl_histogram * hist = gsl_histogram_alloc (nbins);

    gsl_histogram_set_ranges (hist, bins, nbins+1);

    switch (derivorder)
    {
        case 0:
            
            slmeval3dx0 (hist, nbins, bins, coef, ndata, x, dx, y);
            
            break;
            
        case 1:
            
            slmeval3dx1 (hist, nbins, bins, coef, ndata, x, dx, y);
            
            break;
            
            
        case 2:
            
            slmeval3dx2 (hist, nbins, bins, coef, ndata, x, dx, y);
            
            break;
            
        case 3:
            
            slmeval3dx3 (hist, nbins, bins, coef, ndata, x, dx, y);
            
            break;
            
        case -1:
            
            ret = -2;
            
            break;
            
        default:
            
            ret = -1;
            
            break;
            
    }
                  
    gsl_histogram_free (hist);
    
    return ret;
}



// y = mexslmeval(x, knots, coef, evalmode)

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  /*  check for proper number of arguments */
  if(nrhs != 4)
    mexErrMsgTxt("Four inputs required.");
  if(nlhs != 1)
    mexErrMsgTxt("One output required.");

  /* check the input arguments */
  if( !mxIsNumeric(prhs[3]) || mxIsComplex(prhs[3]) || mxGetN(prhs[3])*mxGetM(prhs[3])!=1 ) {
    mexErrMsgTxt("Input derivorder must be a real valued scalar.");
  }
  
  if( !mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) ) {
    mexErrMsgTxt("Input x must be a real-valued matrix.");
  }
  
  if( !mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || mxGetN(prhs[1]) !=1 ) {
    mexErrMsgTxt("Input bins must be a real-valued column vector.");
  }
  
  if( !mxIsNumeric(prhs[2]) || mxIsComplex(prhs[2]) || mxGetN(prhs[2]) != 2 ) {
    mexErrMsgTxt("Input coefs must be a real-valued two column matrix.");
  }

  /*  get the vector inputs */
  double* x = mxGetPr(prhs[0]);
  int xrows = mxGetM(prhs[0]);
  int xcols = mxGetN(prhs[0]);
  int ndata =  xrows * xcols;
  
  double* bins = mxGetPr(prhs[1]);
  int nbins = mxGetM(prhs[1]) - 1;
  
  double* coef = mxGetPr(prhs[2]);
  int ncoef = mxGetM(prhs[1]);
  
  if (ncoef != (nbins+1))
  {
      mexErrMsgTxt("The number of coeffs does not match the number of knots.");
  }
  
  /*  get the scalar input evalmode */
  int evalmode = (int) mxGetScalar(prhs[3]);

  plhs[0] = mxCreateDoubleMatrix(xrows,xcols,mxREAL);
  double * y = mxGetPr(plhs[0]);
  
  /*  call the subroutine */
  int ret = slmeval3 (evalmode, nbins, bins, coef, ndata, x, y);

}
