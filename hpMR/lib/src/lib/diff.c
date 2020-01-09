/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * diff.c
 *
 * Code generation for function 'diff'
 *
 */

/* Include files */
#include "diff.h"
#include "create2dEpi_codegen.h"
#include "create2dEpi_codegen_emxutil.h"

/* Function Definitions */
void diff(const emxArray_real_T *x, emxArray_real_T *y)
{
  int dimSize;
  int ixLead;
  int iyLead;
  double work_data_idx_0;
  int m;
  double tmp1;
  double d;
  dimSize = x->size[1];
  if (x->size[1] == 0) {
    y->size[0] = 1;
    y->size[1] = 0;
  } else {
    ixLead = x->size[1] - 1;
    if (ixLead >= 1) {
      ixLead = 1;
    }

    if (ixLead < 1) {
      y->size[0] = 1;
      y->size[1] = 0;
    } else {
      ixLead = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = x->size[1] - 1;
      emxEnsureCapacity_real_T(y, ixLead);
      if (x->size[1] - 1 != 0) {
        ixLead = 1;
        iyLead = 0;
        work_data_idx_0 = x->data[0];
        for (m = 2; m <= dimSize; m++) {
          tmp1 = x->data[ixLead];
          d = tmp1;
          tmp1 -= work_data_idx_0;
          work_data_idx_0 = d;
          ixLead++;
          y->data[iyLead] = tmp1;
          iyLead++;
        }
      }
    }
  }
}

/* End of code generation (diff.c) */
