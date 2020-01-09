/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * main.c
 *
 * Code generation for function 'main'
 *
 */

/*************************************************************************/
/* This automatically generated example C main file shows how to call    */
/* entry-point functions that MATLAB Coder generated. You must customize */
/* this file for your application. Do not modify this file directly.     */
/* Instead, make a copy of this file, modify it, and integrate it into   */
/* your development environment.                                         */
/*                                                                       */
/* This file initializes entry-point function arguments to a default     */
/* size and value before calling the entry-point functions. It does      */
/* not store or use any values returned from the entry-point functions.  */
/* If necessary, it does pre-allocate memory for returned values.        */
/* You can use this file as a starting point for a main function that    */
/* you can deploy in your application.                                   */
/*                                                                       */
/* After you copy the file, and before you deploy it, you must make the  */
/* following changes:                                                    */
/* * For variable-size function arguments, change the example sizes to   */
/* the sizes that your application requires.                             */
/* * Change the example values of function arguments to the values that  */
/* your application requires.                                            */
/* * If the entry-point functions return values, store these values or   */
/* otherwise use them as required by your application.                   */
/*                                                                       */
/*************************************************************************/
/* Include Files */
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <time.h>
#include "create2dEpi_codegen.h"
#include "create2dEpi_output.h"
#include "create2dEpi_codegen_terminate.h"
#include "create2dEpi_codegen_emxAPI.h"
#include "create2dEpi_codegen_initialize.h"


/* Function Definitions */

/*
 * Arguments    : lots (16)
 * Return Type  : int
 */
extern void create2dEpi_output(int a, double b, int c, double d, double e, double f, double g, double h, double i, char * j, double * A, double * B, int * C, int * D, double * E, double * F)
{
  emxArray_real_T *G;
  emxInitArray_real_T(&G, 2);
  
  bool useEpiType  = (a == 0) ? false : true;
  fflush(stdout);
  printf("Writing output files to %s\n", j);

  int outNGRAD, outNACQ;
  double outTE, outPFY;

  /* Initialize the application.
     You do not need to do this more than one time. */
  create2dEpi_codegen_initialize();

  /* Call the entry-point 'create2dEpi_codegen'. */
  create2dEpi_codegen(useEpiType, b, c, d, e, f, g, h, i, j, 
                    G, &outNGRAD, &outNACQ, &outTE, &outPFY);

  /* Terminate the application.
     You do not need to do this more than one time. */
  create2dEpi_codegen_terminate();

  *C = outNGRAD;
  *D = outNACQ;
  *E = outTE;
  *F = outPFY;
  printf("NGRAD = %d\nNACQ = %d\nTE = %g\nPFY = %g\n", outNGRAD, outNACQ, outTE, outPFY);

  if( outNGRAD < 65536) {
    for(int ii=0; ii<outNGRAD; ii++) {
      A[ii] = G->data[2*ii];
      B[ii] = G->data[2*ii+1];
    }
    emxDestroyArray_real_T(G);
    /*return 0;*/
  } else {
    printf("NGRAD = %d, create2dEpi_codegen_output will not return gradient waveform\n", outNGRAD);
    emxDestroyArray_real_T(G);
    /*return 1;*/
  }
  
}


/* End of code generation (main.c) */

