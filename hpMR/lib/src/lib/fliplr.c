/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * fliplr.c
 *
 * Code generation for function 'fliplr'
 *
 */

/* Include files */
#include "fliplr.h"
#include "create2dEpi_codegen.h"

/* Function Definitions */
void fliplr(emxArray_real_T *x)
{
  int n;
  int nd2;
  int b_j1;
  int j2;
  double xtmp;
  n = x->size[1] - 1;
  nd2 = x->size[1] >> 1;
  for (b_j1 = 0; b_j1 < nd2; b_j1++) {
    j2 = n - b_j1;
    xtmp = x->data[b_j1];
    x->data[b_j1] = x->data[j2];
    x->data[j2] = xtmp;
  }
}

/* End of code generation (fliplr.c) */
