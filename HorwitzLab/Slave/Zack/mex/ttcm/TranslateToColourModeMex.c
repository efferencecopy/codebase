// Here's Zack's implementation of TranslateToColourMode. The comments do an okay job
// of guiding the reader through the logic of the various blocks of code. Please direct
// yourself to `help TranslateToColourMode` for more information about the color
// conversion modes and the legacy method on how we drive the Bits++ in colour mode.
//
// If any changes need to be made to this code, please confirm that the changes
// don't break edge cases. You can run the `test_ttcmmex` function which contains
// tests for _every_ valid permutation of the inputs. The tests in that function
// don't include cases for empty input arguments, but those are trivial to confirm
// (just search for `mxIsEmpty` in this file).

#include <math.h>
#include <string.h>
#include "mex.h"

#define IMAGE_IN prhs[0]
#define IMAGE_OUT plhs[0]

// human readable color conversion mode flags
enum { CCM_LEGACY = -1, CCM_HALFWIDTH, CCM_SUBSAMPLE, CCM_AVERAGE };

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray* prhs[])
{
    double *dp_imagein, *dp_imageout = NULL;
    uint8_T *bp_imagein, *bp_imageout = NULL;
    uint16_T *usp_imagein, *usp_imageout = NULL;
    
    mwSize insize[3], outsize[3], ndims, ninputelems, nrows, ncols, npages;
    mwIndex i, j, k, idxin, idxout;
    uint8_T hibyte, lobyte;
    // sane defaults
    int ccmode = CCM_LEGACY, useoldstyle16bit = 0;
    int scaling_factor = 1;
    bool found_int = false, need_to_free = false, need_scaling_factor = true;
    
    if (nrhs < 1 || nrhs > 3 || mxIsEmpty(IMAGE_IN))
        mexErrMsgTxt("You've specified too few or too many input arguments.\n");
    
    if (nlhs > 1)
        mexErrMsgTxt("There should be at most 1 output argument.\n");
    
    if (!(mxIsUint8(IMAGE_IN) || mxIsUint16(IMAGE_IN) || mxIsDouble(IMAGE_IN)))
        mexErrMsgTxt("Your image matrix must have class uint8, uint16, or double!\n");
    
    switch (nrhs) {
        case 3:
            if (!mxIsEmpty(prhs[2])) {
                ccmode = (int) *mxGetPr(prhs[2]);
                ccmode = (ccmode < CCM_LEGACY || ccmode > CCM_AVERAGE) ? CCM_LEGACY : ccmode;
            }
            // no break needed, I want to fall through
        case 2:
            if (mxIsUint8(IMAGE_IN)) {
                useoldstyle16bit = 0; // force this
            } else if (!mxIsEmpty(prhs[1])) {
                useoldstyle16bit = (int) *mxGetPr(prhs[1]);
                useoldstyle16bit = (useoldstyle16bit < 0 || useoldstyle16bit > 1) ? 0 : useoldstyle16bit;
            }
    }
    
    ndims = mxGetNumberOfDimensions(IMAGE_IN);
    memcpy(&insize, mxGetDimensions(IMAGE_IN), ndims * sizeof(mwSize));
    memcpy(&outsize, &insize, ndims * sizeof(mwSize));
    
    if (ndims < 3) { insize[2] = 1; outsize[2] = 1; }
    npages = insize[2]; ncols = insize[1]; nrows = insize[0];
    ninputelems = npages * ncols * nrows;
    if (ccmode != CCM_HALFWIDTH) outsize[1] *= 2; // only double the columns when the framebuffer is full width
    
    // Allocate memory for the output matrix with the same class as the input matrix when ccmode == CCM_LEGACY
    // Do all internal processing with the double precision `dp_imagein`
    if (mxIsUint8(IMAGE_IN) || mxIsUint16(IMAGE_IN)) {
        dp_imagein = mxCalloc(ninputelems, sizeof(double));
        need_to_free = true;
        need_scaling_factor = ccmode >= CCM_HALFWIDTH; // gotta scale when PsychImaging is running
        if (mxIsUint8(IMAGE_IN)) {
            bp_imagein = (uint8_T *) mxGetData(IMAGE_IN);
            for (i = 0; i < ninputelems; ++i) dp_imagein[i] = bp_imagein[i];
            if (ccmode == CCM_LEGACY) {
                IMAGE_OUT = mxCreateNumericArray(ndims, outsize, mxUINT8_CLASS, mxREAL);
                bp_imageout = (uint8_T *) mxGetData(IMAGE_OUT);
            }
        } else {
            usp_imagein = (uint16_T *) mxGetData(IMAGE_IN);
            for (i = 0; i < ninputelems; ++i) dp_imagein[i] = usp_imagein[i];
            if (ccmode == CCM_LEGACY) {
                IMAGE_OUT = mxCreateNumericArray(ndims, outsize, mxUINT16_CLASS, mxREAL);
                usp_imageout = (uint16_T *) mxGetData(IMAGE_OUT);
            }
        }
        if (ccmode >= CCM_HALFWIDTH) { // this accounts for the weird case of an input uint image matrix with PsychImaging
            IMAGE_OUT = mxCreateNumericArray(ndims, outsize, mxDOUBLE_CLASS, mxREAL);
            dp_imageout = mxGetPr(IMAGE_OUT); // gotta be doubles with PsychImaging
        }
    } else { // mxIsDouble
        dp_imagein = mxGetPr(IMAGE_IN);
        IMAGE_OUT = mxCreateNumericArray(ndims, outsize, mxDOUBLE_CLASS, mxREAL);
        dp_imageout = mxGetPr(IMAGE_OUT);
    }
    
    // Determine how to split up the high and low bytes in ccmode == CCM_LEGACY
    if (ccmode == CCM_LEGACY && !mxIsUint8(IMAGE_IN)) {
        for (i = 0; i < ninputelems; ++i) {
            if (dp_imagein[i] > 255.0) {
                useoldstyle16bit = 1;
                break;
            }
        }
        if (i == ninputelems) useoldstyle16bit = 0;
    }
    
    // Here we look for integers (different from 0,1) in the matrix of doubles. If one is found,
    // we know the matrix can't possibly have intensities ranging from 0 to 1.
    if (need_scaling_factor) {
        for (i = 0; i < ninputelems; ++i) {
            if (dp_imagein[i] != 1.0 && dp_imagein[i] != 0.0
                    && fabs(ceil(dp_imagein[i]) - floor(dp_imagein[i])) < 1E-5) {
                found_int = true;
                break;
            }
        }
        // figure out how to scale the matrix down to the 0.0 - 1.0 range
        if (ccmode >= CCM_HALFWIDTH && found_int) {
            for (j = i; j < ninputelems; ++j) {
                if (dp_imagein[j] > 255.0) {
                    scaling_factor = 65535;
                    break;
                }
            }
            if (j == ninputelems) scaling_factor = 255;
        }
    }
    // legacy Bits++ without PsychImaging (i.e., traditional TranslateToColourMode behavior)
    if ((found_int && ccmode == CCM_LEGACY) || (!mxIsDouble(IMAGE_IN) && ccmode == CCM_LEGACY)) {
        for (k = 0; k < npages; ++k) {
            for (j = 0; j < ncols; ++j) {
                for (i = 0; i < nrows; ++i) {
                    idxin = i + j*nrows + k*nrows*ncols;
                    idxout = i + j*2*nrows + k*nrows*outsize[1];
                    if (useoldstyle16bit) {
                        hibyte = (uint16_T) dp_imagein[idxin] >> 8;
                        lobyte = (uint16_T) dp_imagein[idxin] & 0xFF;
                    } else {
                        hibyte = (uint16_T) dp_imagein[idxin];
                        lobyte = 0;
                    }
                    // match the output's class to the input's
                    if (mxIsUint8(IMAGE_IN)) {
                        bp_imageout[idxout] = hibyte;
                        bp_imageout[idxout+outsize[0]] = lobyte;
                    } else if (mxIsUint16(IMAGE_IN)) {
                        usp_imageout[idxout] = hibyte;
                        usp_imageout[idxout+outsize[0]] = lobyte;
                    } else {
                        dp_imageout[idxout] = hibyte;
                        dp_imageout[idxout+outsize[0]] = lobyte;
                    }
                }
            }
        }
        // case 1: did the user give floating point intensities without PsychImaging? (issue warning)
        // case 2: half-width framebuffer without having to scale the image? pass back the input (no warning required)
    } else if ((!found_int && ccmode == CCM_LEGACY) || (ccmode == CCM_HALFWIDTH && scaling_factor == 1)) {
        mxDestroyArray(IMAGE_OUT);
        IMAGE_OUT = mxDuplicateArray(IMAGE_IN);
        if (ccmode == CCM_LEGACY) {
            mexWarnMsgTxt("You requested to translate to colour mode a matrix containing intensities\nranging "
                    "from 0 to 1 but indicated that PsychImaging isn't running. I'm simply\nreturning a copy "
                    "of the input matrix.");
        }
    } else {
        for (k = 0; k < npages; ++k) {
            for (j = 0; j < ncols; ++j) {
                for (i = 0; i < nrows; ++i) {
                    idxin = i + j*nrows + k*nrows*ncols;
                    idxout = i + j*2*nrows + k*nrows*outsize[1];
                    switch (ccmode) {
                        case CCM_HALFWIDTH:
                            dp_imageout[idxin] = dp_imagein[idxin] / scaling_factor;
                            break;
                        case CCM_AVERAGE: // assign columns 1:2:end (0-based indexing)
                            dp_imageout[idxout+nrows] = dp_imagein[idxin] / scaling_factor; // fall thru to assign columns 0:2:end
                        case CCM_SUBSAMPLE:
                            dp_imageout[idxout] = dp_imagein[idxin] / scaling_factor;
                            break;
                        default:
                            if (need_to_free) mxFree(dp_imagein);
                            mexErrMsgTxt("Will never get here\n");
                    }
                }
            }
        }
    }
    if (need_to_free) mxFree(dp_imagein);
}
