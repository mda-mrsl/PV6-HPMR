/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * create2dEpi_codegen.h
 *
 * Code generation for function 'create2dEpi_codegen'
 *
 */

#ifndef CREATE2DEPI_CODEGEN_H
#define CREATE2DEPI_CODEGEN_H

/* Include files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "create2dEpi_codegen_types.h"

/* Function Declarations */
extern void create2dEpi_codegen(boolean_T epiType, double fov, int mtx, double
  dw, double Gmax, double Smax, double tGrad, double gmr, double pfy, const char
  saveFile[512], emxArray_real_T *G, int *NGRAD, int *NACQ, double *TE, double
  *PFY);

#endif

/* End of code generation (create2dEpi_codegen.h) */
