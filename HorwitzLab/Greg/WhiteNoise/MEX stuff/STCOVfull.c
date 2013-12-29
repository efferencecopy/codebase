
/*****************************************************************/
/*  Function to calculate covariance matrices, efficiently       */
/*  and quickly.                                                 */
/*                                                               */
/*    This function takes two arguments.  Upon initialization,   */
/* the string 'init' as the first argument and the number of     */
/* stimulus dimensions as the second argument.  For subsequent   */
/* calls send the Nx3 RGB matrix as the first argument and a     */
/* Nx1 spikes matrix as the second argument.  Finally, send the  */
/* 'return' string as the last argument to get the output.       */
/* The output will be a structure containing STA, STCOV, and n.  */
/*                                                               */
/* Trying to get it to work with dynamic memory allocation.      */
/* GDLH 4/17/02                                                  */
/*                                                               */
/* Adding parameters for startframe and endframe (endframe       */
/* should always be the larger of the two).                      */
/* GDLH 10/17/02                                                 */
/*                                                               */
/* Now we can compute non-causal kernels.  So for instance       */
/* if you want to go from the time of the spike to 20 frames     */
/* preceding the spike you should make startframe = -18 and      */
/* endframe = 1.  Also, set frameskip = maxT = 0 to prevent      */
/* spikes at the start of each trial from getting omitted.       */
/* GDLH 8/9/04                                                   */
/*****************************************************************/

#include <stdio.h>
#include "mex.h"

/* Inputs */

#define RGBS            prhs[0]
#define SPIKES          prhs[1]

/* Outputs */

#define OUT            plhs[0]

/* Global variables */
double *STS;
double *STCross;
int ndims;

static void FreeMemory(void)
{

  int i;

  if (STS != NULL)
    {
      mxFree(STS);
    }
  if (STCross != NULL)
    {
      mxFree(STCross);
    }
}

void DoInit( int nrhs, const mxArray*prhs[], double *p_n, int *p_npixels, int *p_ncolors, int *p_startframe, int *p_endframe)
{

  int i;

  /* Initializing static variables to zero */
  *p_n = 0;

  /* Error checking  - GDLH include error check for number of elements in cell array */
  if ((nrhs != 2) | !mxIsCell(SPIKES))
  {
	mexErrMsgTxt("Argument after 'init' should be {npixels, ncolors, startframe, endframe}");
  } else {
      *p_npixels = (int) *mxGetPr(mxGetCell(SPIKES, 0));       		*p_ncolors = (int) *mxGetPr(mxGetCell(SPIKES, 1));
      *p_startframe = (int) *mxGetPr(mxGetCell(SPIKES, 2));
      *p_endframe = (int) *mxGetPr(mxGetCell(SPIKES, 3));
	if (*p_endframe < *p_startframe)
	  {
	    mexErrMsgTxt("Endframe must be greater than (or equal to) startframe\n");
	  }
        ndims = *p_npixels* *p_ncolors* (*p_endframe-*p_startframe+1);
        STS = mxCalloc(ndims, sizeof(double));
	STCross = mxCalloc(ndims*ndims, sizeof(double));
        if (STS && STCross)
	  {
	    mexMakeMemoryPersistent(STS);
	    mexMakeMemoryPersistent(STCross);
 	    mexAtExit(FreeMemory);
	  } else {
	    mexErrMsgTxt("Unable to allocate enough memory.  Aborting.");
	  }
    }

}

/* Come here once per spike */
void UpdateSTS(double *rgbarr, int whichframe, double nspikes, 
	       int npixels, int ncolors, int endframe)
{
  register int i;
  int ndimsperframe,count;

  ndimsperframe = npixels*ncolors;
  count = 0;

  for (i=0; i<ndims; i++)
    {
      STS[count] += nspikes*rgbarr[ndimsperframe*(whichframe-endframe+1)+i];
      count++;
    }
}

/* Come here once per spike */
void UpdateSTCross(double *rgbarr, int whichframe, double nspikes, 
	       int npixels, int ncolors, int endframe)
{
  int ndimsperframe;
  register int i,j;

  ndimsperframe = npixels*ncolors;
  for (i=0;i<ndims;i++)
    {
      for (j=0;j<=i;j++)
	{
	  STCross[i*ndims+j] += nspikes*
	    rgbarr[ndimsperframe*(whichframe-endframe+1)+i]*
	    rgbarr[ndimsperframe*(whichframe-endframe+1)+j]; 
	}
    }
}

/* Destructively modifies STS, STCross, and n */
void Update(const mxArray *rgbs, const mxArray *spikes, int npixels,
	    int ncolors, int startframe, int endframe, double *p_n) 
{
  int i, nstim;
  double *sparr, *rgbarr;
  
  sparr = mxGetPr(spikes);
  rgbarr = mxGetPr(rgbs);
  if (mxGetM(spikes) > mxGetN(spikes))
    {
      nstim = mxGetM(spikes);
    } else {
      nstim = mxGetN(spikes);
    }
  for(i=endframe-1;i<nstim;i++)  /* Ignoring startframe  */
    {
 	if (sparr[i] > 0)
	{
	  UpdateSTS(rgbarr,i,sparr[i],npixels,ncolors,endframe); 
	  UpdateSTCross(rgbarr,i,sparr[i],npixels,ncolors,endframe);
	  *p_n += (int) sparr[i];
	}
    }
}

void DupElements(double n)
{
  int i,j;

  for (i=0;i<ndims;i++)
    {
      for (j=0;j<=i;j++)
	{
	          /*    STCross[i*ndims+j] = (n*STCross[i*ndims+j]-(STS[i]*STS[j]))/denom;  */
		      STCross[j*ndims+i]  = STCross[i*ndims+j]; 
	}


    }
}

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
{
  static int npixels, ncolors, startframe, endframe;
  static double n;
  char *input_buf;
  int buflen, status, i;
  mxArray *ReturnVals[3];
 
  if (nlhs > 1) 
    {
      mexErrMsgTxt("Too many output arguments."); 
    } 

  if (mxIsChar(RGBS) == 1) 
    {
      buflen = (mxGetM(prhs[0]) * mxGetN(prhs[0])) + 1;
      input_buf = mxCalloc(buflen, sizeof(char));
      status = mxGetString(prhs[0], input_buf, buflen);
      if (status != 0) 
	mexWarnMsgTxt("Not enough space for string in argument 1");
      if (strcmp(input_buf,"init") == 0)
	{
	  DoInit(nrhs, prhs, &n, &npixels, &ncolors, &startframe, &endframe);
	}
      else if (strcmp(input_buf,"return") == 0)
	{
	  DupElements(n);
	  OUT = mxCreateCellMatrix(1,3);
	  ReturnVals[0] = mxCreateDoubleMatrix(ndims,1,mxREAL); 
	  memcpy(mxGetPr(ReturnVals[0]), STS, ndims*sizeof(double));
	  ReturnVals[1] = mxCreateDoubleMatrix(ndims,ndims,mxREAL); 
	  memcpy(mxGetPr(ReturnVals[1]), STCross, ndims*ndims*sizeof(double)); 
	  ReturnVals[2] = mxCreateDoubleMatrix(1,1,mxREAL); 
	  memcpy(mxGetPr(ReturnVals[2]), &n, sizeof(double)); 

        mxSetCell(OUT,0,ReturnVals[0]);
	  mxSetCell(OUT,1,ReturnVals[1]); 
        mxSetCell(OUT,2,ReturnVals[2]);
        return;
	}
      else
	{
	  mexErrMsgTxt("Unrecognized string argument");
	}
    } 
  else 
    {
      Update(RGBS,SPIKES,npixels,ncolors,startframe,endframe,&n);
    } 
}

