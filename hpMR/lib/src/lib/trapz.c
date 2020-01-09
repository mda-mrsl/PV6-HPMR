/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * trapz.c
 *
 * Code generation for function 'trapz'
 *
 */

/* Include files */
#include "trapz.h"
#include "create2dEpi_codegen.h"
#include "create2dEpi_codegen_emxutil.h"

/* Function Definitions */
double trapz(const emxArray_real_T *x)
{
  double z;
  emxArray_real_T *c;
  int i;
  int ix;
  int iac;
  int ia;
  z = 0.0;
  if (x->size[1] > 1) {
    emxInit_real_T(&c, 1);
    i = c->size[0];
    c->size[0] = x->size[1];
    emxEnsureCapacity_real_T(c, i);
    ix = x->size[1];
    for (i = 0; i < ix; i++) {
      c->data[i] = 1.0;
    }

    c->data[0] = 0.5;
    c->data[x->size[1] - 1] = 0.5;
    ix = 0;
    i = x->size[1];
    for (iac = 1; iac <= i; iac++) {
      for (ia = iac; ia <= iac; ia++) {
        z += x->data[ia - 1] * c->data[ix];
      }

      ix++;
    }

    emxFree_real_T(&c);
  }

  return z;
}

/* End of code generation (trapz.c) */
