/*
 * File: _coder_removeHoles_api.h
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 01-Dec-2017 14:55:34
 */

#ifndef _CODER_REMOVEHOLES_API_H
#define _CODER_REMOVEHOLES_API_H

/* Include Files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_removeHoles_api.h"

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern void removeHoles(real_T x_data[], int32_T x_size[2], real_T y_data[],
  int32_T y_size[2]);
extern void removeHoles_api(const mxArray *prhs[1], const mxArray *plhs[1]);
extern void removeHoles_atexit(void);
extern void removeHoles_initialize(void);
extern void removeHoles_terminate(void);
extern void removeHoles_xil_terminate(void);

#endif

/*
 * File trailer for _coder_removeHoles_api.h
 *
 * [EOF]
 */
