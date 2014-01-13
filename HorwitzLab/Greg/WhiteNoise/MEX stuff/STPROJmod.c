
/*****************************************************************/
/*    Function to quickly find the projections of stimuli        */
/* onto a set of basis vectors.  These basis vectors may         */
/* have temporal, spatial, and color components.                 */    
/*                                                               */
/*    The first time this function is called, it should be with  */
/* the string, 'init', as the first argument and four subsequent */
/* arguments in a cell array.                                    */
/*     1) a matrix of basis vectors.                             */
/*     2) an offset (in frames) that we want to use when         */
/*        comparing the stimulus sequence to the basis vectors.  */
/*        This is basically the number of frames before the      */
/*        spike that we want to skip.                            */
/*     3) an estimate regarding the number of stimuli we'll be   */
/*        considering (This is important for memory allocation). */
/*     4) stimulus dimensions: [npixels ncolors nframes]         */
/*     At the moment, all basis vectors should be of length      */
/* n*64*3 (so that if a single element of the vector occurs on   */
/* frame 'i', then all of the elements on frame 'i' are likewise */
/* included).  The order of elements in the basis vectors should */
/* be such that color changes the fastest (from red to blue),    */
/* space changes second fastest (from 1:64) and time changes     */
/* last (from well before the spike to the time of the spike).   */
/*    On subsequent calls, send the Nx3 RGB matrix as the first  */
/* argument and a Nx1 spikes matrix as the second argument (see  */
/* getWhtnsStats3.m which is the calling function).              */
/*    Finally, send the 'return' string as the last argument to  */
/* get the output which will be a structure containing a matrix  */
/* of projections, and an Lspike vector which provides the       */
/* of spikes fired in response to each stimulus.                 */
/*                                                               */
/* GDLH 5/14/02                                                  */
/*                                                               */
/*****************************************************************/

#include <stdio.h>
#include "mex.h"

/* Inputs */

#define RGBS            prhs[0]
#define SPIKES          prhs[1]
#define BASISMAT		0
#define FRAMEOFFSET     1
#define NSTIMAPPROX     2
#define STIMDIMS        3

/* Outputs */

#define OUT            plhs[0]

/* Global variables */
double *BASIS; 
double *LSPIKE;
double *PROJ;

static void FreeMemory(void)
{
  if (BASIS != NULL)
    {
      mxFree(BASIS);
    }
  if (LSPIKE != NULL)
    {
      mxFree(LSPIKE);
    }
  if (PROJ != NULL)
    {
      mxFree(PROJ);
    }
}

void DoInit( int nrhs, const mxArray*prhs[], double *p_n, int *p_npixels, int *p_ncolors, int *p_nframes, int *p_frameoffset, int *p_nstimapprox, int *p_nBasis)
{

  int i;

  /* Initializing static variables to zero */
  *p_n = 0;

  /* Error checking */
  if (!mxIsCell(SPIKES))  /* Need better error checking here */
    {
	mexErrMsgTxt("Arguments after 'init' should be {[basis_matrix], frame_offset, nstim_approx, and [npixels, ncolors, nframes]}");
    } else {
      *p_nBasis = mxGetN(mxGetCell(SPIKES, 0));
      *p_frameoffset = (int) *mxGetPr(mxGetCell(SPIKES,1));
      *p_nstimapprox = (int) *mxGetPr(mxGetCell(SPIKES,2));
      *p_npixels = (int) *mxGetPr(mxGetCell(SPIKES,3));
      *p_ncolors = (int) *(mxGetPr(mxGetCell(SPIKES,3))+1);
      *p_nframes = (int) *(mxGetPr(mxGetCell(SPIKES,3))+2);
	/* printf("%d %d %d %d %d %d\n",*p_nBasis, *p_frameoffset, *p_nstimapprox, *p_npixels, *p_ncolors, *p_nframes); */
      if (mxGetM(mxGetCell(SPIKES, 0)) != *p_npixels* *p_ncolors* *p_nframes)
	{
	  mexErrMsgTxt("Dimensions implied by BasisMat does not match npixels*ncolors*nframes");
	}

      LSPIKE = mxCalloc(*p_nstimapprox, sizeof(double));
      PROJ = mxCalloc(*p_nstimapprox * *p_nBasis, sizeof(double));
      BASIS = mxCalloc(*p_npixels * *p_ncolors * *p_nframes * *p_nBasis, sizeof(double));
      for (i=0;i<*p_npixels * *p_ncolors * *p_nframes * *p_nBasis; i++)
      {
        BASIS[i] = *(mxGetPr(mxGetCell(SPIKES, 0))+i);
      }
        /* printf("i is %d, nBasis is %d, first element is %f\n",i, *p_nBasis, BASIS[0]); */
	  /* printf("last element is %f\n",BASIS[i-1]); */

      
      if (LSPIKE && PROJ && BASIS)
	{
	  mexMakeMemoryPersistent(LSPIKE);
	  mexMakeMemoryPersistent(PROJ);
	  mexMakeMemoryPersistent(BASIS);
	  mexAtExit(FreeMemory);
	} else {
	  mexErrMsgTxt("Unable to allocate enough memory.  Aborting.");
	}
    }
}

/* Update: Destructively modifies LSPIKE, PROJ, and n.                            */
/* The algorithm is to loop over stimuli (going through the passed in RGB matrix) */
/* to find individual stimulus chunks.  For each stimulus chunk, we loop over the */
/* number of basis vectors onto which we want to project the stimuli.  For each   */
/* combination of stimulus chunk and basis vector, we loop through, element by    */
/* element to take the dot-product.                                               */

void Update(const mxArray *rgbs, const mxArray *spikes, int npixels,
	    int ncolors, int nframes, int frameoffset, int nstimapprox, int nBasis, int *p_n) 
{
  int i, j, k;
  int nrgbs, ndimsperframe, ndimsperstimulus;
  double *sparr, *rgbarr;

  ndimsperframe = npixels*ncolors;
  ndimsperstimulus = ndimsperframe*nframes;

  sparr = mxGetPr(spikes);
  rgbarr = mxGetPr(rgbs);
  nrgbs = mxGetM(rgbs);

  for(i=0;i<nrgbs-ndimsperframe*(nframes+frameoffset-1);i=i+ndimsperframe)    /* Stimulus increment (in vector elements) */
  {
      for (j=0;j<nBasis;j++)  /* Looping over basis vectors */
      {
	for (k=0;k<ndimsperstimulus;k++)  /* Element increment (within stimulus) */
	  {
	    PROJ[(int)*p_n+j*nstimapprox] = PROJ[(int)*p_n+j*nstimapprox] +
                          BASIS[k+j*ndimsperstimulus]*rgbarr[k+i];
	  }
	 /*printf("rgb is %f, basis value is %f, first basis element is %f\n",rgbarr[k+i], BASIS[k+j*ndimsperstimulus], BASIS[0]);*/
	 /*printf("j is %d, idx is %d, last projection calculated was %f\n",j,(int)*p_n+j*nstimapprox,PROJ[(int)*p_n+j*nstimapprox]);*/  
      }
	LSPIKE[*p_n] = sparr[(i/ndimsperframe)+(nframes+frameoffset-1)];

      (*p_n)++;
	if (*p_n > nstimapprox)
	{
	  mexErrMsgTxt("Approximate number of stimuli too low.  Aborting.");
	}
    }
}

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
{
  static int npixels, ncolors, nframes;  /* Gives dimensions of the full stimulus space */
  static int frameoffset; /* How many frames we should skip before each spike */
  static int nstimapprox; /* The approximate length of the PROJ and LSPIKE arrays */
  static int nBasis; /* The number of basis vectors in the reduced space */
  /*  static int ndims; */  /* npixels*ncolors*nframes */
  static int n;   /* running count of how many stimuli we've looked at */
  char *input_buf;
  int buflen, status, i;
  mxArray *ReturnVals[2];
 
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
	  DoInit(nrhs, prhs, &n, &npixels, &ncolors, &nframes, &frameoffset, &nstimapprox, &nBasis);
	}
      else if (strcmp(input_buf,"return") == 0)
	{
	  OUT = mxCreateCellMatrix(1,2);
	  ReturnVals[0] = mxCreateDoubleMatrix(n,nBasis,mxREAL); 
          for (i=0; i<nBasis; i++)
	    {
	      memcpy(mxGetPr(ReturnVals[0])+(int)(i*n), &PROJ[i*nstimapprox], n*sizeof(double));
	    }
	  ReturnVals[1] = mxCreateDoubleMatrix(n,1,mxREAL); 
	  memcpy(mxGetPr(ReturnVals[1]), LSPIKE, n*sizeof(double)); 
        mxSetCell(OUT,0,ReturnVals[0]);
	  mxSetCell(OUT,1,ReturnVals[1]); 
	  return;
	}
      else
	{
	  mexErrMsgTxt("Unrecognized string argument");
	}
    } 
  else 
    {
      if (mxGetM(RGBS) != mxGetN(SPIKES)*npixels*ncolors)
	{
		mexErrMsgTxt("Number of RGB values not the correct multiple of number of elements in Lspike vector");
	}
      Update(RGBS,SPIKES,npixels,ncolors,nframes,frameoffset,nstimapprox,nBasis,&n);
    } 
}
