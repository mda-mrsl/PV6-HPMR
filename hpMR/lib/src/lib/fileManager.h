/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * fileManager.h
 *
 * Code generation for function 'fileManager'
 *
 */

#ifndef FILEMANAGER_H
#define FILEMANAGER_H

/* Include files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "create2dEpi_codegen_types.h"

/* Function Declarations */
extern int b_fileManager(void);
extern signed char cfopen(const char cfilename[512], const char * cpermission);
extern void fileManager(double varargin_1, FILE * *f, boolean_T *a);
extern void filedata_init(void);

#endif

/* End of code generation (fileManager.h) */
