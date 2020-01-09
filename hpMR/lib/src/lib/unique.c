/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * unique.c
 *
 * Code generation for function 'unique'
 *
 */

/* Include files */
#include "unique.h"
#include "create2dEpi_codegen.h"
#include "create2dEpi_codegen_emxutil.h"
#include <math.h>

/* Function Definitions */
void unique_vector(const emxArray_real_T *a, emxArray_real_T *b)
{
  emxArray_int32_T *idx;
  int na;
  int n;
  int i;
  int b_i;
  emxArray_int32_T *iwork;
  int k;
  double x;
  int i2;
  int j;
  int exitg1;
  int pEnd;
  double absx;
  int p;
  int q;
  int exponent;
  int qEnd;
  int kEnd;
  int i1;
  emxInit_int32_T(&idx, 2);
  na = a->size[1];
  n = a->size[1] + 1;
  i = idx->size[0] * idx->size[1];
  idx->size[0] = 1;
  idx->size[1] = a->size[1];
  emxEnsureCapacity_int32_T(idx, i);
  b_i = a->size[1];
  for (i = 0; i < b_i; i++) {
    idx->data[i] = 0;
  }

  if (a->size[1] != 0) {
    emxInit_int32_T(&iwork, 1);
    i = iwork->size[0];
    iwork->size[0] = a->size[1];
    emxEnsureCapacity_int32_T(iwork, i);
    i = a->size[1] - 1;
    for (k = 1; k <= i; k += 2) {
      if (a->data[k - 1] <= a->data[k]) {
        idx->data[k - 1] = k;
        idx->data[k] = k + 1;
      } else {
        idx->data[k - 1] = k + 1;
        idx->data[k] = k;
      }
    }

    if ((a->size[1] & 1) != 0) {
      idx->data[a->size[1] - 1] = a->size[1];
    }

    b_i = 2;
    while (b_i < n - 1) {
      i2 = b_i << 1;
      j = 1;
      for (pEnd = b_i + 1; pEnd < n; pEnd = qEnd + b_i) {
        p = j;
        q = pEnd;
        qEnd = j + i2;
        if (qEnd > n) {
          qEnd = n;
        }

        k = 0;
        kEnd = qEnd - j;
        while (k + 1 <= kEnd) {
          i = idx->data[q - 1];
          i1 = idx->data[p - 1];
          if (a->data[i1 - 1] <= a->data[i - 1]) {
            iwork->data[k] = i1;
            p++;
            if (p == pEnd) {
              while (q < qEnd) {
                k++;
                iwork->data[k] = idx->data[q - 1];
                q++;
              }
            }
          } else {
            iwork->data[k] = i;
            q++;
            if (q == qEnd) {
              while (p < pEnd) {
                k++;
                iwork->data[k] = idx->data[p - 1];
                p++;
              }
            }
          }

          k++;
        }

        for (k = 0; k < kEnd; k++) {
          idx->data[(j + k) - 1] = iwork->data[k];
        }

        j = qEnd;
      }

      b_i = i2;
    }

    emxFree_int32_T(&iwork);
  }

  i = b->size[0] * b->size[1];
  b->size[0] = 1;
  b->size[1] = a->size[1];
  emxEnsureCapacity_real_T(b, i);
  for (k = 0; k < na; k++) {
    b->data[k] = a->data[idx->data[k] - 1];
  }

  emxFree_int32_T(&idx);
  b_i = 0;
  k = 1;
  while (k <= na) {
    x = b->data[k - 1];
    do {
      exitg1 = 0;
      k++;
      if (k > na) {
        exitg1 = 1;
      } else {
        absx = fabs(x / 2.0);
        if (absx <= 2.2250738585072014E-308) {
          absx = 4.94065645841247E-324;
        } else {
          frexp(absx, &exponent);
          absx = ldexp(1.0, exponent - 53);
        }

        if (fabs(x - b->data[k - 1]) >= absx) {
          exitg1 = 1;
        }
      }
    } while (exitg1 == 0);

    b_i++;
    b->data[b_i - 1] = x;
  }

  i = b->size[0] * b->size[1];
  if (1 > b_i) {
    b->size[1] = 0;
  } else {
    b->size[1] = b_i;
  }

  emxEnsureCapacity_real_T(b, i);
}

/* End of code generation (unique.c) */
