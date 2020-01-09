/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * fileManager.c
 *
 * Code generation for function 'fileManager'
 *
 */

/* Include files */
#include "fileManager.h"
#include "create2dEpi_codegen.h"
#include "create2dEpi_codegen_rtwutil.h"
#include <string.h>

/* Variable Definitions */
static FILE * eml_openfiles[20];
static boolean_T eml_autoflush[20];

/* Function Declarations */
static signed char filedata(void);

/* Function Definitions */
static signed char filedata(void)
{
  signed char f;
  int k;
  boolean_T exitg1;
  f = 0;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 20)) {
    if (eml_openfiles[k] == NULL) {
      f = (signed char)(k + 1);
      exitg1 = true;
    } else {
      k++;
    }
  }

  return f;
}

int b_fileManager(void)
{
  int f;
  int j;
  int cst;
  f = 0;
  for (j = 0; j < 20; j++) {
    if (eml_openfiles[j] != NULL) {
      cst = fclose(eml_openfiles[j]);
      if (cst == 0) {
        eml_openfiles[j] = NULL;
        eml_autoflush[j] = true;
      } else {
        f = -1;
      }
    }
  }

  return f;
}

signed char cfopen(const char cfilename[512], const char * cpermission)
{
  signed char fileid;
  signed char j;
  char ccfilename[513];
  FILE * filestar;
  int i;
  fileid = -1;
  j = filedata();
  if (j >= 1) {
    memcpy(&ccfilename[0], &cfilename[0], 512U * sizeof(char));
    ccfilename[512] = '\x00';
    filestar = fopen(&ccfilename[0], cpermission);
    if (filestar != NULL) {
      eml_openfiles[j - 1] = filestar;
      eml_autoflush[j - 1] = true;
      i = j + 2;
      if (i > 127) {
        i = 127;
      }

      fileid = (signed char)i;
    }
  }

  return fileid;
}

void fileManager(double varargin_1, FILE * *f, boolean_T *a)
{
  signed char fileid;
  fileid = (signed char)rt_roundd(varargin_1);
  if ((fileid < 0) || (varargin_1 != fileid)) {
    fileid = -1;
  }

  if (fileid >= 3) {
    *f = eml_openfiles[fileid - 3];
    *a = eml_autoflush[fileid - 3];
  } else if (fileid == 0) {
    *f = stdin;
    *a = true;
  } else if (fileid == 1) {
    *f = stdout;
    *a = true;
  } else if (fileid == 2) {
    *f = stderr;
    *a = true;
  } else {
    *f = NULL;
    *a = true;
  }
}

void filedata_init(void)
{
  FILE * a;
  int i;
  a = NULL;
  for (i = 0; i < 20; i++) {
    eml_autoflush[i] = false;
    eml_openfiles[i] = a;
  }
}

/* End of code generation (fileManager.c) */
