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
#include "main.h"


/* Function Definitions */

/*
 * Arguments    : int argc
 *                const char * const argv[]
 * Return Type  : int
 */
int main(int argc, const char * const argv[])
{
  /*(void)argc;
  (void)argv;*/

  int epiType, mtx, NGRAD, NACQ;
  double fov, dw, Gmax, Smax, tGrad, gmr, pfy, TE, PFY;
  char saveFile[512];
  double GradX[65535], GradY[65535] = {0.0};
  if( argc != 11){
    printf("Usage:\n create2dEpi_codegen epiType fov mtx dw Gmax Smax tGrad gmr pfy saveFile\n");
    return 1;
  } else {
    epiType  = atol(argv[1]);
    fov      = atof(argv[2]);
    mtx      = atol(argv[3]);
    dw       = atof(argv[4]);
    Gmax     = atof(argv[5]);
    Smax     = atof(argv[6]);
    tGrad    = atof(argv[7]);
    gmr      = atof(argv[8]);
    pfy      = atof(argv[9]);
    strcpy(saveFile, argv[10]);
    printf("saveFile = %s\n", saveFile);
  }

  /* Call create2dEpi_codegen_output */
  create2dEpi_output(epiType, fov, mtx, dw, Gmax, Smax, tGrad, 
          gmr, pfy, saveFile, GradX, GradY, &NGRAD, &NACQ, &TE, &PFY);

  /*if (retval != 0) {
    printf("Error encountered in create2dEpi_codegen_output\n");
    return 1;
  }*/
  
  printf("NGRAD = %d\nNACQ = %d\nTE = %g\nPFY = %g\n", NGRAD, NACQ, TE, PFY);
  if (NGRAD > 7) {
    int div_denom = 8;
    int cur_idx   = 0;
    for (int ii=0; ii<div_denom; ii++) {
      cur_idx = (ii+1)*NGRAD/div_denom;
      printf("Grad Point %d = [%g, %g] mT/m\n", cur_idx, GradX[cur_idx-1], GradY[cur_idx-1]);
    }
  }
  
  return 0;
}

/* End of code generation (main.c) */
