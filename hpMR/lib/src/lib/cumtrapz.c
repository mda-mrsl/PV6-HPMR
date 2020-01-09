/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * cumtrapz.c
 *
 * Code generation for function 'cumtrapz'
 *
 */

/* Include files */
#include "cumtrapz.h"
#include "create2dEpi_codegen.h"
#include "create2dEpi_codegen_emxutil.h"

/* Function Definitions */
void b_cumtrapz(const emxArray_creal_T *x, emxArray_creal_T *z)
{
  int i;
  double s_re;
  double s_im;
  int iyz;
  double ylast_re;
  double ylast_im;
  int k;
  i = z->size[0] * z->size[1];
  z->size[0] = 1;
  z->size[1] = x->size[1];
  emxEnsureCapacity_creal_T(z, i);
  if (x->size[1] != 0) {
    s_re = 0.0;
    s_im = 0.0;
    iyz = 0;
    ylast_re = x->data[0].re;
    ylast_im = x->data[0].im;
    z->data[0].re = 0.0;
    z->data[0].im = 0.0;
    i = x->size[1];
    for (k = 0; k <= i - 2; k++) {
      iyz++;
      ylast_re += x->data[iyz].re;
      ylast_im += x->data[iyz].im;
      if (ylast_im == 0.0) {
        ylast_re /= 2.0;
        ylast_im = 0.0;
      } else if (ylast_re == 0.0) {
        ylast_re = 0.0;
        ylast_im /= 2.0;
      } else {
        ylast_re /= 2.0;
        ylast_im /= 2.0;
      }

      s_re += ylast_re;
      s_im += ylast_im;
      ylast_re = x->data[iyz].re;
      ylast_im = x->data[iyz].im;
      z->data[iyz].re = s_re;
      z->data[iyz].im = s_im;
    }
  }
}

void cumtrapz(const emxArray_real_T *x, emxArray_real_T *z)
{
  int i;
  double s;
  int iyz;
  double ylast;
  int k;
  i = z->size[0] * z->size[1];
  z->size[0] = 1;
  z->size[1] = x->size[1];
  emxEnsureCapacity_real_T(z, i);
  if (x->size[1] != 0) {
    s = 0.0;
    iyz = 0;
    ylast = x->data[0];
    z->data[0] = 0.0;
    i = x->size[1];
    for (k = 0; k <= i - 2; k++) {
      iyz++;
      s += (ylast + x->data[iyz]) / 2.0;
      ylast = x->data[iyz];
      z->data[iyz] = s;
    }
  }
}

/* End of code generation (cumtrapz.c) */
