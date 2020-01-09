/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * find.c
 *
 * Code generation for function 'find'
 *
 */

/* Include files */
#include "find.h"
#include "create2dEpi_codegen.h"
#include "create2dEpi_codegen_emxutil.h"

/* Function Definitions */
void eml_find(const emxArray_boolean_T *x, int i_data[], int i_size[2])
{
  emxArray_int32_T *i;
  int k;
  int idx;
  int b_i;
  int ii;
  boolean_T exitg1;
  emxInit_int32_T(&i, 2);
  k = (1 <= x->size[1]);
  idx = 0;
  b_i = i->size[0] * i->size[1];
  i->size[0] = 1;
  i->size[1] = k;
  emxEnsureCapacity_int32_T(i, b_i);
  ii = 0;
  exitg1 = false;
  while ((!exitg1) && (ii <= x->size[1] - 1)) {
    if (x->data[ii]) {
      idx++;
      i->data[idx - 1] = ii + 1;
      if (idx >= k) {
        exitg1 = true;
      } else {
        ii++;
      }
    } else {
      ii++;
    }
  }

  if (k == 1) {
    if (idx == 0) {
      i->size[0] = 1;
      i->size[1] = 0;
    }
  } else {
    b_i = i->size[0] * i->size[1];
    i->size[1] = (1 <= idx);
    emxEnsureCapacity_int32_T(i, b_i);
  }

  i_size[0] = 1;
  i_size[1] = i->size[1];
  ii = i->size[0] * i->size[1];
  for (b_i = 0; b_i < ii; b_i++) {
    i_data[b_i] = i->data[b_i];
  }

  emxFree_int32_T(&i);
}

/* End of code generation (find.c) */
