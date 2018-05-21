#include "interpUtil.h"
#include <cmath>
#include <cstring>

#ifdef _MSC_VER
#include <ppl.h>
#define nParMin 1024
#endif

#define ERR_HEAD_PPUVAL "ppuval: "
#define ERR_HEAD_PPMVAL "ppmval: "

inline int findidx(double v, int l, int u, const double *a)
{
    if (a[l] >= v) return l;
    if (a[u] < v) return u;
    int m;
    while (l < u)
    {
        m = (l+u) / 2;
        (v >= a[m] ? l = m : u = m);
        if ((u - l) <= 1) return l;
    }
    // This return is spurious but dampens the compiler warning
    return l;
}

double nested_sum(unsigned int nest_levels,const double *x, const double **b,
	const int *ix, const int *p, const int *d,const double *c, const int *c_dim, int &ind)
{
    double sum = 0;

    if ( nest_levels == 1 )
    {
        // end of recursion
        double lx = x[0]-b[0][ix[0]];
        double pow_x = 1;

		ind += (ix[0] + p[0]*d[0]) * c_dim[0];
		int dec_x = p[0] * c_dim[0];
		int inc_x = p[0] * d[0] * c_dim[0];
        for ( int i = d[0]; i >= 1; i-- )
        {
            // Update the index from which retrieve coefficient
			ind -= dec_x;
			sum += pow_x*c[ind];
            pow_x *= lx;
        }
    }
    else
    {
        double lx = x[0]-b[0][ix[0]];
        double pow_x = 1;

		ind += (ix[0] + p[0]*d[0]) * c_dim[0];
		int dec_x = p[0] * c_dim[0];
        int inc_x = p[0] * d[0] * c_dim[0];
        for ( int i = d[0]; i > 0; i-- )
        {
			ind -= dec_x;
			sum += pow_x*nested_sum(nest_levels-1,x+1,b+1,ix+1,p+1,d+1,c,c_dim+1,ind);
            pow_x *= lx;
        }
    }
    return sum;
}

double eval_polynomial(const double *xPtr, int x_dims, const double **b, const int *ix,
        const int *p, const int *d, const double *coefs, const int* coefs_dim, int f_dim)
{
    switch (x_dims)
	{
	case 2 :
	{
		// optimize x_dims = 2 case to avoid recursion
		int dec_x = p[0] * coefs_dim[0];
        int	dec_y = p[1] * coefs_dim[1];
        int inc_y = p[1] * d[1] * coefs_dim[1];
		int off_x = f_dim + (ix[0] + p[0]*d[0]) * coefs_dim[0];
        int off_y = (ix[1] + p[1]*d[1]) * coefs_dim[1];

		double x = xPtr[0]-b[0][ix[0]];
        double y = xPtr[1]-b[1][ix[1]];
		double sum = 0;
        double pow_x = 1;
		double pow_y;
        for ( int i = d[0]; i > 0; i-- )
        {
            off_x -= dec_x;
			pow_y = 1;
            for ( int j = d[1]; j > 0; j-- )
            {
                off_y -= dec_y;
				sum += pow_x*pow_y*coefs[off_x + off_y];
                pow_y *= y;
            }
			off_y += inc_y;
            pow_x *= x;
        }
        return sum;
		break;
	}
	case 1 :
	{
		// optimize x_dims = 1 case to avoid recursion.
        // For one dimensional polynomial, it is easy to take
        // advantage of nested valuation. Saves few CPU cycles
        double sum = 0;
        double x = xPtr[0]-b[0][ix[0]];

        int off_x = f_dim+ix[0]*coefs_dim[0];
        int inc_x = p[0] * coefs_dim[0];

        sum = coefs[off_x];
        for (int l = 1; l < d[0]; l++)
        {
            off_x += inc_x;
            sum = sum * x + coefs[off_x];
        }

        return sum;
		break;
	}
	case 3 :
	{
        // optimize x_dims = 3 case to avoid recursion
        int dec_x = p[0] * coefs_dim[0];
        int	dec_y = p[1] * coefs_dim[1];
		int	dec_z = p[2] * coefs_dim[2];
        int inc_y = p[1] * d[1] * coefs_dim[1];
		int inc_z = p[2] * d[2] * coefs_dim[2];
		int off_x = f_dim + (ix[0] + p[0]*d[0]) * coefs_dim[0];
        int off_y = (ix[1] + p[1]*d[1]) * coefs_dim[1];
		int off_z = (ix[2] + p[2]*d[2]) * coefs_dim[2];

		double x = xPtr[0]-b[0][ix[0]];
        double y = xPtr[1]-b[1][ix[1]];
		double z = xPtr[2]-b[2][ix[2]];
		double sum = 0;
        double pow_x = 1;
		double pow_y;
		double pow_z;
        for (int i = d[0]; i > 0; i--  )
        {
			off_x -= dec_x;
			pow_y = 1;
            for (  int j = d[1]; j > 0; j-- )
            {
				off_y -= dec_y;
				pow_z = 1;
                for ( int k =  d[2]; k > 0; k--)
                {
					off_z -= dec_z;
					sum += pow_x*pow_y*pow_z*coefs[off_x + off_y+ off_z];
					pow_z *= z;
				}
				off_z += inc_z;
                pow_y *= y;
            }
			off_y += inc_y;
            pow_x *= x;
        }
        return sum;
		break;
	}
	default :
		// General case goes here
		// If you want, feel free to add more non recursive nested loops.
		// My personal need restricts to dim = 3
        return nested_sum(x_dims,xPtr,b,ix,p,d,coefs,coefs_dim,f_dim);
	}
}

void ppuval( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    // Check proper number of arguments:
    if ( nrhs < 2 )
    {
        mexErrMsgTxt(ERR_HEAD_PPUVAL "not enough input arguments");
        return;
    }
    else if (nrhs > 2 )
    {
        mexErrMsgTxt(ERR_HEAD_PPUVAL "too many input arguments");
    }
    const mxArray *X = prhs[0];
    const mxArray *pp_form = prhs[1];
    size_t n = mxGetNumberOfElements(X);
    if ( !mxIsDouble(X) )
    {
        mexErrMsgTxt(ERR_HEAD_PPUVAL "Evaluation sites must be numeric");
    }
    else if (!mxIsStruct(pp_form))
    {
        mexErrMsgTxt(ERR_HEAD_PPUVAL "Polynomial must be in 'pp'-form");
    }
    double *x = mxGetPr(X);
    mxArray *Form = mxGetField(pp_form,0,"form");
    mxArray *Breaks = mxGetField(pp_form,0,"breaks");
    mxArray *Coefs = mxGetField(pp_form,0,"coefs");
    mxArray *Order = mxGetField(pp_form,0,"order");
    mxArray *Dim = mxGetField(pp_form,0,"dim");
    mxArray *Pieces = mxGetField(pp_form,0,"pieces");
    if (Form == 0 || Breaks == 0 || Coefs == 0 || Order == 0 || Dim == 0 ||
        Pieces == 0 || strcmp(mxArrayToString(Form),"pp") != 0)
    {
        mexErrMsgTxt(ERR_HEAD_PPUVAL "Polynomial must be in 'pp'-form");
    }
    // assembly the pp-form struct
    double *breaks = mxGetPr(Breaks);
    double *pc = mxGetPr(Coefs);
    double *d_pieces = mxGetPr(Pieces);
    double *d_order = mxGetPr(Order);
    double *d_dim = mxGetPr(Dim);
    if (pc == 0 || d_pieces == 0 || d_order == 0 || d_dim == 0)
    {
        mexErrMsgTxt(ERR_HEAD_PPUVAL "Polynomial must be in 'pp'-form");
    }

    size_t pieces = (size_t)*d_pieces;
    size_t order = (size_t)*d_order;
    size_t dim = (size_t)*d_dim;

    if ( breaks == 0 || dim != 1 )
    {
        mexErrMsgTxt(ERR_HEAD_PPUVAL "Polynomial is not univariate."
            "Use ppmval-function for multivariate valuation");
    }

    plhs[0] = mxCreateNumericArray(mxGetNumberOfDimensions(X),
              mxGetDimensions(X) ,mxDOUBLE_CLASS, mxREAL);
    double *out = mxGetPr(plhs[0]);
#ifdef _MSC_VER
    if (n <= nParMin)
	{
#endif
		size_t j;
		double d;
		double lx;
		for (size_t i = 0; i < n; i++)
		{
			j = findidx(x[i], 0, pieces-1, breaks); /* breaks[j] <= X[i] < breaks[j+1] */
			lx = x[i] - breaks[j];
			d = pc[j];
			for (size_t l = 1; l < order; l++)
			{
				d = d * lx + pc[j+l*pieces];
			}
			out[i] = d;
		}
#ifdef _MSC_VER
	}
	else
	{
		Concurrency::parallel_for (size_t(0), n, [&](size_t i)
		{
			size_t j = findidx(x[i], 0, pieces-1, breaks); /* breaks[j] <= X[i] < breaks[j+1] */
			double lx = x[i] - breaks[j];
			double d = pc[j];
			for (size_t l = 1; l < order; l++)
			{
				d = d * lx + pc[j+l*pieces];
			}
			out[i] = d;
		});

	}
#endif
}

void ppmval( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    // Check proper number of arguments:
    if ( nrhs < 2 )
    {
        mexErrMsgTxt(ERR_HEAD_PPUVAL "not enough input arguments");
        return;
    }
    else if (nrhs > 2 )
    {
        mexErrMsgTxt(ERR_HEAD_PPUVAL "too many input arguments");
    }

    const mxArray *X = prhs[0];
    const mxArray *pp_form = prhs[1];

    if ( !mxIsDouble(X) )
    {
        mexErrMsgTxt(ERR_HEAD_PPUVAL "Evaluation sites must be numeric");
    }
    else if (mxGetNumberOfDimensions(X) != 2)
    {
        mexErrMsgTxt(ERR_HEAD_PPUVAL "Evaluation sites must be 2-d array");
    }
    else if (!mxIsStruct(pp_form))
    {
        mexErrMsgTxt(ERR_HEAD_PPUVAL "Polynomial must be in 'pp'-form");
    }
    size_t X_dim = mxGetM(X);
    size_t n_sites = mxGetN(X);

    double *x = mxGetPr(X);
    mxArray *Form = mxGetField(pp_form,0,"form");
    mxArray *Breaks = mxGetField(pp_form,0,"breaks");
    mxArray *Coefs = mxGetField(pp_form,0,"coefs");
    mxArray *Order = mxGetField(pp_form,0,"order");
    mxArray *Dim = mxGetField(pp_form,0,"dim");
    mxArray *Pieces = mxGetField(pp_form,0,"pieces");

    if (Form == 0 || Breaks == 0 || Coefs == 0 || Order == 0 || Dim == 0 ||
        Pieces == 0 || strcmp(mxArrayToString(Form),"pp") != 0)
    {
        mexErrMsgTxt(ERR_HEAD_PPUVAL "Polynomial must be in 'pp'-form");
    }
    // assembly the pp-form struct
    double *p = mxGetPr(Pieces);
    double *o = mxGetPr(Order);
    double *d_dim = mxGetPr(Dim);
    if (p == 0 || o == 0 || d_dim == 0)
    {
        mexErrMsgTxt(ERR_HEAD_PPUVAL "Polynomial must be in 'pp'-form");
    }

    if ( !mxIsCell(Breaks) )
    {
        mexErrMsgTxt(ERR_HEAD_PPMVAL "Polynomial is not multivariate."
            "Use ppval-function for univariate valuation");
    }
    // polynomial coefficients
    size_t coefs_dims = mxGetNumberOfDimensions(Coefs);
	const mwSize *coefs_dim = mxGetDimensions(Coefs);
    double *c = mxGetPr(Coefs);
    // polynomial is a mapping from R^(x_dims) -> R^(f_dims)
    size_t f_dims = (size_t) *d_dim;
    size_t x_dims = mxGetNumberOfElements(Breaks);

    const double **brk_pointers = new const double*[x_dims];
    const double *brk_pointer;
    int *pieces = new int[x_dims];  // number of intervals = number of breaks - 1
    int *order = new int[x_dims];   // polynomial order + 1
	int *prod_coefs_dim = new int[1+x_dims]; // preallocate arrays to store indices
    int *ix;
	prod_coefs_dim[0] = coefs_dim[0];

    for ( size_t i = 0; i < x_dims; i++ )
    {
        brk_pointer = mxGetPr(mxGetCell(Breaks,i));
        if (brk_pointer == 0)
        {
            mexErrMsgTxt(ERR_HEAD_PPMVAL "Polynomial breaks are corrupted in 'pp'-description");
        }
        else
        {
            brk_pointers[i] = brk_pointer;
			pieces[i] = (int) p[i];
			order[i] = (int) o[i];
			prod_coefs_dim[i+1] = prod_coefs_dim[i]*coefs_dim[i+1];
        }
    }

    // allocate memory to the result that is returned to Matlab
    plhs[0] = mxCreateDoubleMatrix(f_dims,n_sites, mxREAL);
    double *out = mxGetPr(plhs[0]);
	int off_x;
#ifdef _MSC_VER
	if (n_sites <= nParMin)
	{
#endif
		ix = new int[x_dims];
		/* serial execution route */
		for (size_t site = 0; site < n_sites; site++)
		{
			off_x = site*x_dims;
			for (size_t x_dim = 0; x_dim < x_dims; x_dim++)
			{
				//brk_pointer = mxGetPr(mxGetCell(breaks,x_dim));
				brk_pointer = brk_pointers[x_dim];
				// before polynomial can be evaluated on the given site, its all breaks needs to be found
				ix[x_dim] = findidx(x[x_dim+off_x],0,pieces[x_dim]-1,brk_pointer);
			}
			for (size_t f_dim = 0; f_dim < f_dims; f_dim++)
			{
				out[f_dim+site*f_dims] =  eval_polynomial(x+off_x, x_dims, brk_pointers,
                    ix, pieces, order,c, prod_coefs_dim, f_dim);
			}
		}
		delete[] ix;
#ifdef _MSC_VER
	}
	else
	{
		ix = new int[x_dims*n_sites]; // preallocate arrays to store indices
		/* parallel version */
		Concurrency::parallel_for (size_t(0), n_sites, [&](size_t site)
		{
			const double *brk_pointer;
			int off_x = site*x_dims;
			for (size_t x_dim = 0; x_dim < x_dims; x_dim++)
			{
				brk_pointer = brk_pointers[x_dim];
				// before polynomial can be evaluated on the given site, its all breaks needs to be found
				ix[x_dim+off_x] = findidx(x[x_dim+off_x],0,pieces[x_dim]-1,brk_pointer);
			}
		});
		Concurrency::parallel_for (size_t(0), n_sites, [&](size_t site)
		{
			int off_x = site*x_dims;
			for (size_t f_dim = 0; f_dim < f_dims; f_dim++)
			{
				out[f_dim+site*f_dims]  =  eval_polynomial(x+off_x, x_dims, brk_pointers,
                    ix, pieces, order,c, prod_coefs_dim, f_dim);
			}
		});
		delete[] ix;
	}
#endif
	delete[] prod_coefs_dim;
    delete[] brk_pointers;
    delete[] pieces;
    delete[] order;
}

