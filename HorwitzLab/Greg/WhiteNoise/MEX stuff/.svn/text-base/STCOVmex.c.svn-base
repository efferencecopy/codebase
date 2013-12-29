/*****************************************************************/
/*  Function to calculate covariance matrices, efficiently       */
/*  and quickly.                                                 */
/*                                                               */
/*    This function takes two arguments.  Upon initialization,   */
/* the string 'init' as the first argument and the number of     */
/* stimulus dimensions as the second argument in the format:     */
/* [npixs, ncols, nframes].                                      */
/*      In subsequent calls, send as the first argument a vector */
/* of RGBs (the total number of elements should be xpix*ncols*   */
/* nframes, and the order should be first all the elements for   */
/* frame 1, then all the elements for frame 2, and so on).       */
/* The second argument should be an nframesx1 matrix of spike    */
/* counts.                                                       */
/*       Send the 'return' string as a solo argument to get the  */
/* output - a structure containing {STS, STCross, and n.}        */
/*                                                               */
/* Trying to get it to work with dynamic memory allocation.      */
/* GDLH 4/17/02                                                  */
/*                                                               */
/* Based on STCOVmod3.c, this mex file will calculate STCross    */
/* matrices strictly on a frame-by-frame basis for use with      */
/* WhiteNoiseOnline.m                                            */
/* GDLH 2/24/08                                                  */
/*                                                               */
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

static void FreeMemory(void)
{

  if (STS != NULL)
    {
      mxFree(STS);
    }
  if (STCross != NULL)
    {
      mxFree(STCross);
    }
}

void DoInit( int nrhs, const mxArray*prhs[], double *p_n, int *p_npixels, int *p_ncolors, int *p_nframes)
{

  int ndims;
  int ndimsperframe;

  /* Initializing static variables to zero */
  *p_n = 0;

  /* Error checking */
  if ((nrhs != 2) | !mxIsCell(SPIKES))
  {
	mexErrMsgTxt("Argument after 'init' should be {npixels, ncolors, nframes}");
  } else {
      *p_npixels = (int) *mxGetPr(mxGetCell(SPIKES, 0));
      *p_ncolors = (int) *mxGetPr(mxGetCell(SPIKES, 1));
      *p_nframes = (int) *mxGetPr(mxGetCell(SPIKES, 2));
      ndimsperframe = *p_npixels* *p_ncolors;
      ndims = *p_npixels* *p_ncolors* *p_nframes;

      FreeMemory();
      STS = mxCalloc(ndims, sizeof(double));
	STCross = mxCalloc((ndimsperframe*ndimsperframe)* *p_nframes, sizeof(double)); /* no cross products between frames */
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
	       int npixels, int ncolors, int nframes)
{
  int i, j, count, ndims, ndimsperframe;

  ndimsperframe = npixels*ncolors;
  ndims = ndimsperframe*nframes;
  count = 0;

  
  /* We fill in the frame at the time of the spike first and */
  /* then take the previous frames in backward chronological order. */ 
  for (i=whichframe;i>whichframe-nframes;i--)
  {
    for (j=0; j<ndimsperframe; j++)
    {
       STS[count] += nspikes*rgbarr[i*ndimsperframe+j];
       count++;
    }
  }
}

/* Come here once per spike */
void UpdateSTCross(double *rgbarr, int whichframe, double nspikes, 
	       int npixels, int ncolors, int nframes)
{
  int ndimsperframe;
  register int i,j,k;

  ndimsperframe = npixels*ncolors;
  
  for (i=0;i<nframes;i++)
  {
     for (j=0;j<ndimsperframe;j++)
     {
        for (k=0;k<=j;k++)
	  {
	    STCross[i*(ndimsperframe*ndimsperframe)+j*ndimsperframe+k] += nspikes*
	         rgbarr[ndimsperframe*(whichframe-i)+j]*
	         rgbarr[ndimsperframe*(whichframe-i)+k]; 
	  }
     }
  }
}

/* Destructively modifies STS, STCross, and n */
void Update(const mxArray *rgbs, const mxArray *spikes, int npixels,
	    int ncolors, int nframes, double *p_n) 
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
  for(i=nframes-1;i<nstim;i++)
    {
      if (sparr[i] > 0)
	{
	  UpdateSTS(rgbarr,i,sparr[i],npixels,ncolors,nframes);
	  UpdateSTCross(rgbarr,i,sparr[i],npixels,ncolors,nframes);
	  *p_n += (int) sparr[i];
	}
    }
}

/* Not converting to STA and STCOV - just returning STS and STCross */ 
void DupElements(double n, int npixels, int ncolors, int nframes)
{
  int i,j,k,ndimsperframe,idx;
  double denom;
  
  denom = n*(n-1);
  ndimsperframe = npixels*ncolors;
  
  for (i=0;i<nframes;i++)
  {
     idx = i*(ndimsperframe*ndimsperframe);
     for (j=0;j<ndimsperframe;j++)
     {
        for (k=0;k<=j;k++)
	  {
		STCross[idx+k*ndimsperframe+j]  = STCross[idx+j*ndimsperframe+k]; 
	  }
     } 
  }
}

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
{
  static int npixels, ncolors, nframes;
  static double n;
  char *input_buf;
  int buflen, status, i, ndims, ndimsperframe;
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
	  DoInit(nrhs, prhs, &n, &npixels, &ncolors, &nframes);
	}
      else if (strcmp(input_buf,"return") == 0)
	{
	  ndimsperframe = npixels*ncolors;
        ndims = ndimsperframe*nframes;
        DupElements(n, npixels, ncolors, nframes);
	  OUT = mxCreateCellMatrix(1,3);

	  ReturnVals[0] = mxCreateDoubleMatrix(ndimsperframe,nframes,mxREAL);
	  memcpy(mxGetPr(ReturnVals[0]), STS, ndims*sizeof(double)); 


	  ReturnVals[1] = mxCreateDoubleMatrix(ndimsperframe*ndimsperframe,nframes,mxREAL); 
	  memcpy(mxGetPr(ReturnVals[1]), STCross, ndimsperframe*ndimsperframe*nframes*sizeof(double)); 
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
      Update(RGBS,SPIKES,npixels,ncolors,nframes,&n);
    } 
}
