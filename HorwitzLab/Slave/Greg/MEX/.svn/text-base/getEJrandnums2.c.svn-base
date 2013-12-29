
/*****************************************************************/
/*  My attempt to implement EJ's random number generation        */
/*  routine as a MEX file.  Hopefully this will help in speeding */
/*  up getWhtnsStim.m.                                           */
/*                                                               */
/*  GDLH 7/23/01                                                 */
/*                                                               */
/*  Adding a third parameter to allow us to skip a bunch of      */
/*  the random numbers.  For instance, if the third parameter is */
/*  set to "192" then we return 192 numbers, skip the next 192,  */
/*  and so on until we've filled the array.                      */
/*                                                               */
/*  GDLH 4/7/04                                                  */
/*****************************************************************/

#include <stdio.h>
#include "mex.h"

/* Inputs */

#define N_NUMS        prhs[0]
#define SEED          prhs[1]
#define SKIP         prhs[2]

/* Outputs */

#define RAND_OUT      plhs[0]

/*******************/
/* EJ's code below */
/*******************/

/* This replicates the original Macintosh Toolbox Random() routine.
// The funciton is defined here explicitly for portability
// and independence from changes in the MacOS.  EJC 1999-12-22
// return value formerly 'short'*/
int StimRandShort(int *seed)
{
	int temp1, temp2, temp3;
	int result;

	temp1 = (*seed & 0xFFFF) * 0x41A7;
	temp2 = (*seed >> 16) * 0x41A7 + (temp1 >> 16);
	temp3 = temp2 * 2 >> 16;
	temp1 = (temp1 & 0xFFFF) - 0x7FFFFFFF;
	temp2 &= 0x7FFF;
	temp1 += (temp2 << 16) | (temp2 >> 16) + temp3;
	if (temp1 < 0)
		temp1 += 0x7FFFFFFF;
	*seed = temp1;
	result = temp1 & 0xFFFF;
	if (result == 0x8000)
		result = 0;
	return result;
}

void getnums(double op[], int iter, int seed, int skip)
{
   int i;
   int j;
   int skiptoggle;

   j = 0;
   skiptoggle = 0;
   i = 0;
   while(i<iter)
   {
     if (skiptoggle)
       {
          StimRandShort(&seed);
       }
     else
       {
         op[i] = StimRandShort(&seed);
         i++;
       }
     j++;
     if (j == skip)
       {
	 j = 0;
	 skiptoggle = !skiptoggle;
       }
   }
}

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
{
    double *op, *np, *sp, *skip;
    int intseed, intn, intskip;

    intskip = 0;

    if (nrhs < 2) 
    {
	mexErrMsgTxt("Two or three input arguments required.");      
    }
    np = mxGetPr(N_NUMS);
    sp = mxGetPr(SEED);
    intseed = (int)(*sp);
    intn = (int)(*np);
    if (nrhs == 3)
    {
       skip = mxGetPr(SKIP);
       intskip = (int)(*skip);
    }
    if (nlhs > 1) 
    {
       mexErrMsgTxt("Too many output arguments."); 
    } 

    /* Assign pointers to each input and output */
    
    /* Create matrix for the return argument */
    RAND_OUT = mxCreateDoubleMatrix (*np, 1,mxREAL);
    op = mxGetPr(RAND_OUT);
    /* Do the actual work */
    getnums(op,intn,intseed,intskip);
    return;
}
