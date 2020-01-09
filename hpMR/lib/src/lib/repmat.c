/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * repmat.c
 *
 * Code generation for function 'repmat'
 *
 */

/* Include files */
#include "repmat.h"
#include "create2dEpi_codegen.h"
#include "create2dEpi_codegen_emxutil.h"

/* Function Definitions */
void repmat(const emxArray_real_T *a, double varargin_2, emxArray_real_T *b)
{
  int ncols;
  int i;
  int jtilecol;
  int ibtile;
  int jcol;
  ncols = b->size[0] * b->size[1];
  b->size[0] = 1;
  i = (int)varargin_2;
  b->size[1] = a->size[1] * i;
  emxEnsureCapacity_real_T(b, ncols);
  ncols = a->size[1];
  for (jtilecol = 0; jtilecol < i; jtilecol++) {
    ibtile = jtilecol * ncols;
    for (jcol = 0; jcol < ncols; jcol++) {
      b->data[ibtile + jcol] = a->data[jcol];
    }
  }
}

/* End of code generation (repmat.c) */
