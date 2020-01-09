/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * create2dEpi_codegen_initialize.c
 *
 * Code generation for function 'create2dEpi_codegen_initialize'
 *
 */

/* Include files */
#include "create2dEpi_codegen_initialize.h"
#include "create2dEpi_codegen.h"
#include "create2dEpi_codegen_data.h"
#include "fileManager.h"

/* Function Definitions */
void create2dEpi_codegen_initialize(void)
{
  filedata_init();
  isInitialized_create2dEpi_codegen = true;
}

/* End of code generation (create2dEpi_codegen_initialize.c) */
