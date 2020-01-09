/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * nullAssignment.c
 *
 * Code generation for function 'nullAssignment'
 *
 */

/* Include files */
#include "nullAssignment.h"
#include "create2dEpi_codegen.h"
#include "create2dEpi_codegen_emxutil.h"

/* Function Definitions */
void b_nullAssignment(emxArray_boolean_T *x, const emxArray_int32_T *idx)
{
  emxArray_boolean_T *b;
  int nxin;
  int i;
  int k0;
  int k;
  int nxout;
  emxInit_boolean_T(&b, 2);
  nxin = x->size[1];
  i = b->size[0] * b->size[1];
  b->size[0] = 1;
  b->size[1] = x->size[1];
  emxEnsureCapacity_boolean_T(b, i);
  k0 = x->size[1];
  for (i = 0; i < k0; i++) {
    b->data[i] = false;
  }

  i = idx->size[1];
  for (k = 0; k < i; k++) {
    b->data[idx->data[k] - 1] = true;
  }

  k0 = 0;
  i = b->size[1];
  for (k = 0; k < i; k++) {
    k0 += b->data[k];
  }

  nxout = x->size[1] - k0;
  k0 = -1;
  for (k = 0; k < nxin; k++) {
    if ((k + 1 > b->size[1]) || (!b->data[k])) {
      k0++;
      x->data[k0] = x->data[k];
    }
  }

  emxFree_boolean_T(&b);
  i = x->size[0] * x->size[1];
  if (1 > nxout) {
    x->size[1] = 0;
  } else {
    x->size[1] = nxout;
  }

  emxEnsureCapacity_boolean_T(x, i);
}

void nullAssignment(emxArray_real_T *x, const emxArray_boolean_T *idx)
{
  int nxin;
  int k0;
  int i;
  int k;
  int nxout;
  nxin = x->size[1];
  k0 = 0;
  i = idx->size[1];
  for (k = 0; k < i; k++) {
    k0 += idx->data[k];
  }

  nxout = x->size[1] - k0;
  k0 = -1;
  for (k = 0; k < nxin; k++) {
    if ((k + 1 > idx->size[1]) || (!idx->data[k])) {
      k0++;
      x->data[k0] = x->data[k];
    }
  }

  i = x->size[0] * x->size[1];
  if (1 > nxout) {
    x->size[1] = 0;
  } else {
    x->size[1] = nxout;
  }

  emxEnsureCapacity_real_T(x, i);
}

/* End of code generation (nullAssignment.c) */
