/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * interp1.c
 *
 * Code generation for function 'interp1'
 *
 */

/* Include files */
#include "interp1.h"
#include "create2dEpi_codegen.h"
#include "create2dEpi_codegen_emxutil.h"

/* Function Definitions */
void interp1(const emxArray_real_T *varargin_1, const emxArray_real_T
             *varargin_2, const emxArray_real_T *varargin_3, emxArray_real_T *Vq)
{
  emxArray_real_T *y;
  int j2;
  int nd2;
  emxArray_real_T *x;
  int nx;
  unsigned int outsize_idx_1;
  int low_i;
  double xtmp;
  int low_ip1;
  int mid_i;
  double r;
  emxInit_real_T(&y, 2);
  j2 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = varargin_2->size[1];
  emxEnsureCapacity_real_T(y, j2);
  nd2 = varargin_2->size[0] * varargin_2->size[1];
  for (j2 = 0; j2 < nd2; j2++) {
    y->data[j2] = varargin_2->data[j2];
  }

  emxInit_real_T(&x, 2);
  j2 = x->size[0] * x->size[1];
  x->size[0] = 1;
  x->size[1] = varargin_1->size[1];
  emxEnsureCapacity_real_T(x, j2);
  nd2 = varargin_1->size[0] * varargin_1->size[1];
  for (j2 = 0; j2 < nd2; j2++) {
    x->data[j2] = varargin_1->data[j2];
  }

  nx = varargin_1->size[1] - 1;
  outsize_idx_1 = (unsigned int)varargin_3->size[1];
  j2 = Vq->size[0] * Vq->size[1];
  Vq->size[0] = 1;
  Vq->size[1] = (int)outsize_idx_1;
  emxEnsureCapacity_real_T(Vq, j2);
  nd2 = (int)outsize_idx_1;
  for (j2 = 0; j2 < nd2; j2++) {
    Vq->data[j2] = 0.0;
  }

  if (varargin_3->size[1] != 0) {
    if (varargin_1->data[1] < varargin_1->data[0]) {
      j2 = varargin_1->size[1] >> 1;
      for (low_i = 0; low_i < j2; low_i++) {
        xtmp = x->data[low_i];
        nd2 = nx - low_i;
        x->data[low_i] = x->data[nd2];
        x->data[nd2] = xtmp;
      }

      nd2 = varargin_2->size[1] >> 1;
      for (low_i = 0; low_i < nd2; low_i++) {
        j2 = (varargin_2->size[1] - low_i) - 1;
        xtmp = y->data[low_i];
        y->data[low_i] = y->data[j2];
        y->data[j2] = xtmp;
      }
    }

    nd2 = varargin_3->size[1];
    for (j2 = 0; j2 < nd2; j2++) {
      xtmp = varargin_3->data[j2];
      if ((xtmp <= x->data[x->size[1] - 1]) && (xtmp >= x->data[0])) {
        nx = x->size[1];
        low_i = 1;
        low_ip1 = 2;
        while (nx > low_ip1) {
          mid_i = (low_i >> 1) + (nx >> 1);
          if (((low_i & 1) == 1) && ((nx & 1) == 1)) {
            mid_i++;
          }

          if (varargin_3->data[j2] >= x->data[mid_i - 1]) {
            low_i = mid_i;
            low_ip1 = mid_i + 1;
          } else {
            nx = mid_i;
          }
        }

        xtmp = x->data[low_i - 1];
        r = (varargin_3->data[j2] - xtmp) / (x->data[low_i] - xtmp);
        if (r == 0.0) {
          Vq->data[j2] = y->data[low_i - 1];
        } else if (r == 1.0) {
          Vq->data[j2] = y->data[low_i];
        } else {
          xtmp = y->data[low_i - 1];
          if (xtmp == y->data[low_i]) {
            Vq->data[j2] = xtmp;
          } else {
            Vq->data[j2] = (1.0 - r) * xtmp + r * y->data[low_i];
          }
        }
      }
    }
  }

  emxFree_real_T(&x);
  emxFree_real_T(&y);
}

/* End of code generation (interp1.c) */
