/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * create2dEpi_codegen.c
 *
 * Code generation for function 'create2dEpi_codegen'
 *
 */

/* Include files */
#include "create2dEpi_codegen.h"
#include "colon.h"
#include "create2dEpi_codegen_data.h"
#include "create2dEpi_codegen_emxutil.h"
#include "create2dEpi_codegen_initialize.h"
#include "create2dEpi_codegen_rtwutil.h"
#include "cumtrapz.h"
#include "diff.h"
#include "fileManager.h"
#include "find.h"
#include "fliplr.h"
#include "interp1.h"
#include "nullAssignment.h"
#include "print_processing.h"
#include "repmat.h"
#include "trapz.h"
#include "unique.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

/* Custom Source Code */
#include "create2dEpi_output.h"
#include "main.h"

/* Function Declarations */
static void findPattern2(const emxArray_real_T *array, const emxArray_real_T
  *pattern, emxArray_real_T *start);
static double rt_hypotd(double u0, double u1);
static double rt_remd(double u0, double u1);

/* Function Definitions */
static void findPattern2(const emxArray_real_T *array, const emxArray_real_T
  *pattern, emxArray_real_T *start)
{
  emxArray_boolean_T *locs;
  int i;
  double b_pattern;
  int loop_ub;
  emxArray_int32_T *ii;
  int nx;
  int idx;
  boolean_T exitg1;
  emxArray_boolean_T *b_locs;
  emxInit_boolean_T(&locs, 2);

  /* #################################################################### */
  /* findPattern2 Locate a pattern in an array. */
  /*  */
  /*    indices = findPattern2(array, pattern) finds the starting indices of */
  /*    pattern within array. */
  /*  */
  /*    Example: */
  /*    a = [0 1 4 9 16 4 9]; */
  /*    patt = [4 9]; */
  /*    indices = findPattern2(a,patt) */
  /*    indices = */
  /*         3     6 */
  /*  Let's assume for now that both the pattern and the array are non-empty */
  /*  VECTORS, but there's no checking for this.  */
  /*  For this algorithm, I loop over the pattern elements. */
  /*  First, find candidate locations; i.e., match the first element in the */
  /*  pattern. */
  i = locs->size[0] * locs->size[1];
  locs->size[0] = 1;
  locs->size[1] = array->size[1];
  emxEnsureCapacity_boolean_T(locs, i);
  b_pattern = pattern->data[0];
  loop_ub = array->size[0] * array->size[1];
  for (i = 0; i < loop_ub; i++) {
    locs->data[i] = (array->data[i] == b_pattern);
  }

  emxInit_int32_T(&ii, 2);
  nx = locs->size[1];
  idx = 0;
  i = ii->size[0] * ii->size[1];
  ii->size[0] = 1;
  ii->size[1] = locs->size[1];
  emxEnsureCapacity_int32_T(ii, i);
  loop_ub = 0;
  exitg1 = false;
  while ((!exitg1) && (loop_ub <= nx - 1)) {
    if (locs->data[loop_ub]) {
      idx++;
      ii->data[idx - 1] = loop_ub + 1;
      if (idx >= nx) {
        exitg1 = true;
      } else {
        loop_ub++;
      }
    } else {
      loop_ub++;
    }
  }

  if (locs->size[1] == 1) {
    if (idx == 0) {
      ii->size[0] = 1;
      ii->size[1] = 0;
    }
  } else {
    i = ii->size[0] * ii->size[1];
    if (1 > idx) {
      ii->size[1] = 0;
    } else {
      ii->size[1] = idx;
    }

    emxEnsureCapacity_int32_T(ii, i);
  }

  i = start->size[0] * start->size[1];
  start->size[0] = 1;
  start->size[1] = ii->size[1];
  emxEnsureCapacity_real_T(start, i);
  loop_ub = ii->size[0] * ii->size[1];
  for (i = 0; i < loop_ub; i++) {
    start->data[i] = ii->data[i];
  }

  emxFree_int32_T(&ii);

  /*  Next remove start values that are too close to the end to possibly match */
  /*  the pattern. */
  i = locs->size[0] * locs->size[1];
  locs->size[0] = 1;
  locs->size[1] = start->size[1];
  emxEnsureCapacity_boolean_T(locs, i);
  loop_ub = start->size[0] * start->size[1];
  for (i = 0; i < loop_ub; i++) {
    locs->data[i] = (((unsigned int)start->data[i] + pattern->size[1]) - 1U >
                     (unsigned int)array->size[1]);
  }

  nullAssignment(start, locs);

  /*  Next, loop over elements of pattern, usually much shorter than length of */
  /*  array, to check which possible locations are valid still. */
  i = pattern->size[1];
  emxInit_boolean_T(&b_locs, 2);
  for (nx = 0; nx <= i - 2; nx++) {
    /*  check viable locations in array */
    idx = locs->size[0] * locs->size[1];
    locs->size[0] = 1;
    locs->size[1] = start->size[1];
    emxEnsureCapacity_boolean_T(locs, idx);
    b_pattern = pattern->data[nx + 1];
    loop_ub = start->size[0] * start->size[1];
    for (idx = 0; idx < loop_ub; idx++) {
      locs->data[idx] = (b_pattern == array->data[(int)((start->data[idx] +
        ((double)nx + 2.0)) - 1.0) - 1]);
    }

    /*  delete false ones from indices */
    idx = b_locs->size[0] * b_locs->size[1];
    b_locs->size[0] = 1;
    b_locs->size[1] = locs->size[1];
    emxEnsureCapacity_boolean_T(b_locs, idx);
    loop_ub = locs->size[0] * locs->size[1];
    for (idx = 0; idx < loop_ub; idx++) {
      b_locs->data[idx] = !locs->data[idx];
    }

    nullAssignment(start, b_locs);
  }

  emxFree_boolean_T(&b_locs);
  emxFree_boolean_T(&locs);
}

static double rt_hypotd(double u0, double u1)
{
  double y;
  double a;
  double b;
  a = fabs(u0);
  b = fabs(u1);
  if (a < b) {
    a /= b;
    y = b * sqrt(a * a + 1.0);
  } else if (a > b) {
    b /= a;
    y = a * sqrt(b * b + 1.0);
  } else {
    y = a * 1.4142135623730951;
  }

  return y;
}

static double rt_remd(double u0, double u1)
{
  double y;
  double q;
  if ((u1 != 0.0) && (u1 != trunc(u1))) {
    q = fabs(u0 / u1);
    if (fabs(q - floor(q + 0.5)) <= DBL_EPSILON * q) {
      y = 0.0;
    } else {
      y = fmod(u0, u1);
    }
  } else {
    y = fmod(u0, u1);
  }

  return y;
}

void create2dEpi_codegen(boolean_T epiType, double fov, int mtx, double dw,
  double Gmax, double Smax, double tGrad, double gmr, double pfy, const char
  saveFile[512], emxArray_real_T *G, int *NGRAD, int *NACQ, double *TE, double
  *PFY)
{
  int i;
  int mty;
  char b_saveFile[513];
  char c_saveFile[513];
  char fileName[512];
  int i1;
  char cv[5];
  static const char cv1[5] = "tAcq";
  signed char fileid;
  FILE * b_NULL;
  FILE * filestar;
  boolean_T autoflush;
  double np;
  int i2;
  double dk;
  double gRead;
  emxArray_real_T *gBlip;
  int nIter;
  double fRamp;
  emxArray_creal_T *g;
  emxArray_creal_T *k;
  emxArray_real_T *gTmp;
  emxArray_real_T *gFlyb;
  emxArray_real_T *gPre;
  emxArray_real_T *gPreY;
  emxArray_real_T *gPreX;
  emxArray_real_T *t;
  emxArray_real_T *T;
  emxArray_real_T *KTmp;
  emxArray_real_T *idxADC;
  emxArray_boolean_T *idxAcq;
  emxArray_real_T *tAcq;
  emxArray_real_T *b_gRead;
  emxArray_real_T *kTmp;
  emxArray_int32_T *r;
  emxArray_real_T *r1;
  emxArray_int32_T *r2;
  emxArray_real_T *b_gTmp;
  emxArray_boolean_T *r3;
  emxArray_real_T *b_T;
  int exitg3;
  double nRamp;
  int loop_ub;
  double kRamp;
  double sRamp;
  double nRead;
  int b_loop_ub;
  double kBlip;
  double kRead;
  int n;
  double ex;
  int b_k;
  double d;
  double kMax;
  int exitg2;
  int b_i;
  int nEchoes;
  boolean_T guard1 = false;
  int tmp_data[1];
  int tmp_size[2];
  int exitg1;
  int iidx;
  int i3;
  unsigned int unnamed_idx_1;
  double validatedHoleFilling[5];
  if (isInitialized_create2dEpi_codegen == false) {
    create2dEpi_codegen_initialize();
  }

  /* CREATE2DEPI_CODEGEN designs 2D flyback or symmetric EPI trajectories for snapshot imaging */
  /*  */
  /*    Usage: [G, NGRAD, NACQ, TE, PFY] = create2dEpi_codegen(... */
  /*      epiType, fov, mtx, dw, Gmax, Smax, tGrad, gmr, pfy, saveFile) */
  /*  */
  /*    Inputs: */
  /*      epiType,  false for flyback, true for symmetric */
  /*      fov,      Field of view                [cm] */
  /*      mtx,      Matrix size  */
  /*      dw,       Acquisition dwell time       [us] */
  /*      Gmax,     Max gradient amplitude       [mT/m] */
  /*      Smax,     Max gradient slew rate       [T/m/s] */
  /*      tGrad,    Gradient timing resolution   [us] */
  /*      gmr,      Gyromagnetic ratio           [MHz/T] */
  /*      pfy,      Fraction of k-space sampled in blipped direction */
  /*      saveFile, String specifying file location to write text file 'tAcq'  */
  /*                containing acquisition time values for recon */
  /*                  */
  /*    Outputs: */
  /*      G,     2D snapshot EPI readout gradient in form Gx+iGy    [mT/m] */
  /*      NGRAD, Number of gradient points in readout waveform */
  /*      NACQ,  Number of acquisition points in readout waveform */
  /*      TE,    Echo time delay (from first gradient point) */
  /*      PFY,   Fraction of k-space sampled in blipped direction */
  /*    */
  /*  06/2019, Keith Michel */
  /*  Parse Inputs */
  printf("Entering %s (%s)\n\n", "create2dEpi", "v20191017");
  fflush(stdout);

  /*  roundoff tolerance for acquisition window k-vals */
  i = G->size[0] * G->size[1];
  G->size[0] = 2;
  G->size[1] = 1;
  emxEnsureCapacity_real_T(G, i);
  G->data[0] = 0.0;
  G->data[1] = 0.0;
  *NGRAD = 0;
  *NACQ = 0;
  *TE = 0.0;
  mty = (int)ceil((double)mtx * fmax(fmin(pfy, 1.0), 0.5));
  pfy = (double)mty / (double)mtx;
  memcpy(&b_saveFile[0], &saveFile[0], 512U * sizeof(char));
  b_saveFile[512] = '\x00';
  memcpy(&c_saveFile[0], &b_saveFile[0], 513U * sizeof(char));
  strcpy(fileName, c_saveFile);
  for (i1 = 0; i1 < 5; i1++) {
    cv[i1] = cv1[i1];
  }

  strcat(fileName, cv);

  /*  if ~saveFile */
  /*      fprintf(1, 'Not writing any output files\n'); */
  /*  else */
  memcpy(&b_saveFile[0], &fileName[0], 512U * sizeof(char));
  b_saveFile[512] = '\x00';
  memcpy(&c_saveFile[0], &b_saveFile[0], 513U * sizeof(char));
  printf("Saving acquisition times to %s\n", c_saveFile);
  fflush(stdout);
  fileid = cfopen(fileName, "wb");
  b_NULL = NULL;
  fileManager(fileid, &filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    fprintf(filestar,
            "# Acquisition timepoints in milliseconds from %s (%s). Positive values indicate recon samples\n",
            "create2dEpi", "v20191017");
    if (autoflush) {
      fflush(filestar);
    }
  }

  b_NULL = NULL;
  fileManager(fileid, &filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    i = mtx;
    if (mtx > 32767) {
      i = 32767;
    } else {
      if (mtx < -32768) {
        i = -32768;
      }
    }

    fprintf(filestar,
            "# Inputs: EPITYPE = %d, FOV = %g [cm], MTX = %d, DW = %g [us], GMAX = %g [mT/m], SMAX = %g [T/m/s], TGRAD = %g [us], GMR = %g [M"
            "Hz/T], PFY = %g, SAVEFILE = %s\n", (signed char)epiType, fov,
            (short)i, dw, Gmax, Smax, tGrad, gmr, pfy, b_saveFile);
    if (autoflush) {
      fflush(filestar);
    }
  }

  /*  end */
  if (rt_remd(dw, tGrad * 2.0) > 1.0E-6) {
    fprintf(stderr,
            "create2dEpi:dwellTime\nAcquisition dwell time must be even integer multiple of gradient timing resolution\n");
    fflush(stderr);
    b_fileManager();
  } else {
    np = rt_roundd(dw / tGrad);
    dw = tGrad * np;
    if (!epiType) {
      i = mtx;
      if (mtx > 32767) {
        i = 32767;
      } else {
        if (mtx < -32768) {
          i = -32768;
        }
      }

      i2 = mty;
      if (mty > 32767) {
        i2 = 32767;
      } else {
        if (mty < -32768) {
          i2 = -32768;
        }
      }

      printf("Flyback EPI, %g cm FOV, %d x %d mtx, %g Hz BW \n", fov, (short)i,
             (short)i2, 1.0E+6 / dw);
      fflush(stdout);
    } else {
      i = mtx;
      if (mtx > 32767) {
        i = 32767;
      } else {
        if (mtx < -32768) {
          i = -32768;
        }
      }

      i2 = mty;
      if (mty > 32767) {
        i2 = 32767;
      } else {
        if (mty < -32768) {
          i2 = -32768;
        }
      }

      printf("Symmetric EPI, %g cm FOV, %d x %d mtx, %g Hz BW \n", fov, (short)i,
             (short)i2, 1.0E+6 / dw);
      fflush(stdout);
    }

    /*  Design readout gradient */
    /*  force ramp duration to be integer multiple of dw/2 */
    /*  for symmetric waveforms, force phase blip duration to be same as ramps */
    dk = 1.0 / fov;

    /*  [cycles/cm] */
    gRead = 100000.0 * dk / gmr / dw;

    /*  [mT/m] */
    if (gRead > Gmax) {
      fprintf(stderr,
              "create2dEpi:maxGrad\nMaximum gradient amplitude exceeded on readout. Increase FOV or dwell time.\n");
      fflush(stderr);
      b_fileManager();
    } else {
      emxInit_real_T(&gBlip, 2);
      nIter = 0;
      fRamp = 0.5;
      i = gBlip->size[0] * gBlip->size[1];
      gBlip->size[0] = 1;
      gBlip->size[1] = 1;
      emxEnsureCapacity_real_T(gBlip, i);
      gBlip->data[0] = 0.0;
      if (!epiType) {
        printf("Designing readout ramp");
        fflush(stdout);
      } else {
        printf("Designing readout ramp and phase blip");
        fflush(stdout);
      }

      emxInit_creal_T(&g, 2);
      emxInit_creal_T(&k, 2);
      emxInit_real_T(&gTmp, 2);
      emxInit_real_T(&gFlyb, 2);
      emxInit_real_T(&gPre, 2);
      emxInit_real_T(&gPreY, 2);
      emxInit_real_T(&gPreX, 2);
      emxInit_real_T(&t, 2);
      emxInit_real_T(&T, 2);
      emxInit_real_T(&KTmp, 2);
      emxInit_real_T(&idxADC, 2);
      emxInit_boolean_T(&idxAcq, 2);
      emxInit_real_T(&tAcq, 2);
      emxInit_real_T(&b_gRead, 2);
      emxInit_real_T(&kTmp, 2);
      emxInit_int32_T(&r, 2);
      emxInit_real_T(&r1, 2);
      emxInit_int32_T(&r2, 2);
      emxInit_real_T(&b_gTmp, 2);
      emxInit_boolean_T(&r3, 2);
      emxInit_real_T(&b_T, 1);
      do {
        exitg3 = 0;
        printf(".");
        fflush(stdout);
        nIter++;
        if (nIter > 40) {
          fprintf(stderr,
                  "create2dEpi:rampIter\nMaximum iterations met in ramp design while loop\n");
          fflush(stderr);
          b_fileManager();
          exitg3 = 1;
        } else {
          nRamp = fRamp * np;
          if (nRamp < 0.0) {
            gTmp->size[0] = 1;
            gTmp->size[1] = 0;
          } else {
            i = gTmp->size[0] * gTmp->size[1];
            gTmp->size[0] = 1;
            loop_ub = (int)floor(nRamp);
            gTmp->size[1] = loop_ub + 1;
            emxEnsureCapacity_real_T(gTmp, i);
            for (i = 0; i <= loop_ub; i++) {
              gTmp->data[i] = i;
            }
          }

          i = gTmp->size[0] * gTmp->size[1];
          i2 = gTmp->size[0] * gTmp->size[1];
          gTmp->size[0] = 1;
          emxEnsureCapacity_real_T(gTmp, i2);
          loop_ub = i - 1;
          for (i = 0; i <= loop_ub; i++) {
            gTmp->data[i] = gTmp->data[i] * gRead / nRamp;
          }

          sRamp = (gTmp->data[1] - gTmp->data[0]) * 1000.0 / tGrad;
          if (!epiType) {
            if (sRamp > Smax) {
              fRamp += 0.5;
            } else {
              exitg3 = 2;
            }
          } else {
            i = r1->size[0] * r1->size[1];
            r1->size[0] = 1;
            r1->size[1] = gTmp->size[1];
            emxEnsureCapacity_real_T(r1, i);
            loop_ub = gTmp->size[0] * gTmp->size[1];
            for (i = 0; i < loop_ub; i++) {
              r1->data[i] = gTmp->data[i];
            }

            fliplr(r1);
            i = gBlip->size[0] * gBlip->size[1];
            gBlip->size[0] = 1;
            gBlip->size[1] = gTmp->size[1] + r1->size[1];
            emxEnsureCapacity_real_T(gBlip, i);
            loop_ub = gTmp->size[1];
            for (i = 0; i < loop_ub; i++) {
              gBlip->data[i] = gTmp->data[i];
            }

            loop_ub = r1->size[1];
            for (i = 0; i < loop_ub; i++) {
              gBlip->data[i + gTmp->size[1]] = r1->data[i];
            }

            kBlip = trapz(gBlip) * gmr * tGrad * 1.0E-5;
            i = gBlip->size[0] * gBlip->size[1];
            i2 = gBlip->size[0] * gBlip->size[1];
            gBlip->size[0] = 1;
            emxEnsureCapacity_real_T(gBlip, i2);
            loop_ub = i - 1;
            for (i = 0; i <= loop_ub; i++) {
              gBlip->data[i] = gBlip->data[i] * dk / kBlip;
            }

            if ((sRamp > Smax) || ((gBlip->data[1] - gBlip->data[0]) * 1000.0 /
                                   tGrad > Smax)) {
              fRamp += 0.5;
            } else {
              n = gBlip->size[1];
              if (gBlip->size[1] <= 2) {
                if (gBlip->size[1] == 1) {
                  ex = gBlip->data[0];
                } else if (gBlip->data[0] < gBlip->data[1]) {
                  ex = gBlip->data[1];
                } else {
                  ex = gBlip->data[0];
                }
              } else {
                ex = gBlip->data[0];
                for (b_k = 2; b_k <= n; b_k++) {
                  d = gBlip->data[b_k - 1];
                  if (ex < d) {
                    ex = d;
                  }
                }
              }

              if (ex > Gmax) {
                fRamp += 0.5;
              } else {
                exitg3 = 2;
              }
            }
          }
        }
      } while (exitg3 == 0);

      if (exitg3 != 1) {
        printf("\n\t%g us readout ramp, %g acquisition points\n", nRamp * tGrad,
               fRamp);
        fflush(stdout);
        kRamp = trapz(gTmp) * gmr * tGrad * 1.0E-5;

        /*  [cycles/cm] */
        sRamp = dw * (double)mtx;

        /*  [us] */
        nRead = rt_roundd(sRamp / tGrad);
        i = r1->size[0] * r1->size[1];
        r1->size[0] = 1;
        r1->size[1] = gTmp->size[1];
        emxEnsureCapacity_real_T(r1, i);
        loop_ub = gTmp->size[0] * gTmp->size[1];
        for (i = 0; i < loop_ub; i++) {
          r1->data[i] = gTmp->data[i];
        }

        fliplr(r1);
        i = b_gRead->size[0] * b_gRead->size[1];
        b_gRead->size[0] = 1;
        loop_ub = (int)nRead;
        b_gRead->size[1] = (gTmp->size[1] + loop_ub) + r1->size[1];
        emxEnsureCapacity_real_T(b_gRead, i);
        b_loop_ub = gTmp->size[1];
        for (i = 0; i < b_loop_ub; i++) {
          b_gRead->data[i] = gTmp->data[i];
        }

        for (i = 0; i < loop_ub; i++) {
          b_gRead->data[i + gTmp->size[1]] = gRead;
        }

        b_loop_ub = r1->size[1];
        for (i = 0; i < b_loop_ub; i++) {
          b_gRead->data[(i + gTmp->size[1]) + loop_ub] = r1->data[i];
        }

        kRead = trapz(b_gRead) * gmr * tGrad * 1.0E-5;

        /*  [cycles/cm] */
        n = b_gRead->size[1];
        if (b_gRead->size[1] <= 2) {
          if (b_gRead->size[1] == 1) {
            ex = b_gRead->data[0];
          } else if (b_gRead->data[0] < b_gRead->data[1]) {
            ex = b_gRead->data[1];
          } else {
            ex = b_gRead->data[0];
          }
        } else {
          ex = b_gRead->data[0];
          for (b_k = 2; b_k <= n; b_k++) {
            d = b_gRead->data[b_k - 1];
            if (ex < d) {
              ex = d;
            }
          }
        }

        printf("\t%g us readout, %g mT/m gradient\n", sRamp, ex);
        fflush(stdout);

        /*  Design flyback/prephasing gradient */
        /*  force duration to be integer multiple of dw */
        nIter = 0;
        kBlip = 0.5;
        fRamp = 1.0;
        kMax = (double)mtx / fov / 2.0;
        i = gFlyb->size[0] * gFlyb->size[1];
        gFlyb->size[0] = 1;
        gFlyb->size[1] = 1;
        emxEnsureCapacity_real_T(gFlyb, i);
        gFlyb->data[0] = 0.0;
        i = gPre->size[0] * gPre->size[1];
        gPre->size[0] = 1;
        gPre->size[1] = 1;
        emxEnsureCapacity_real_T(gPre, i);
        gPre->data[0] = 0.0;
        if (!epiType) {
          printf("Designing flyback gradient");
          fflush(stdout);
        } else {
          printf("Designing prephasing gradient");
          fflush(stdout);
        }

        do {
          exitg2 = 0;
          printf(".");
          fflush(stdout);
          nIter++;
          if (nIter > 40) {
            fprintf(stderr,
                    "create2dEpi:prepIter\nMaximum iterations met in flyback/prephasing lobe design while loop\n");
            fflush(stderr);
            b_fileManager();
            exitg2 = 2;
          } else {
            ex = kBlip * np;
            if (ex < 0.0) {
              gTmp->size[0] = 1;
              gTmp->size[1] = 0;
            } else {
              i = gTmp->size[0] * gTmp->size[1];
              gTmp->size[0] = 1;
              loop_ub = (int)floor(ex);
              gTmp->size[1] = loop_ub + 1;
              emxEnsureCapacity_real_T(gTmp, i);
              for (i = 0; i <= loop_ub; i++) {
                gTmp->data[i] = i;
              }
            }

            i = gTmp->size[0] * gTmp->size[1];
            i2 = gTmp->size[0] * gTmp->size[1];
            gTmp->size[0] = 1;
            emxEnsureCapacity_real_T(gTmp, i2);
            loop_ub = i - 1;
            for (i = 0; i <= loop_ub; i++) {
              gTmp->data[i] = gTmp->data[i] * Smax * tGrad / 1000.0;
            }

            b_k = gTmp->size[1];
            for (b_i = 0; b_i < b_k; b_i++) {
              if (gTmp->data[b_i] > Gmax) {
                gTmp->data[b_i] = Gmax;
              }
            }

            i = r1->size[0] * r1->size[1];
            r1->size[0] = 1;
            r1->size[1] = gTmp->size[1];
            emxEnsureCapacity_real_T(r1, i);
            loop_ub = gTmp->size[0] * gTmp->size[1];
            for (i = 0; i < loop_ub; i++) {
              r1->data[i] = gTmp->data[i];
            }

            fliplr(r1);
            d = gTmp->data[gTmp->size[1] - 1];
            i = b_gTmp->size[0] * b_gTmp->size[1];
            b_gTmp->size[0] = 1;
            loop_ub = (int)(fRamp * np);
            b_gTmp->size[1] = (gTmp->size[1] + loop_ub) + r1->size[1];
            emxEnsureCapacity_real_T(b_gTmp, i);
            b_loop_ub = gTmp->size[1];
            for (i = 0; i < b_loop_ub; i++) {
              b_gTmp->data[i] = gTmp->data[i];
            }

            for (i = 0; i < loop_ub; i++) {
              b_gTmp->data[i + gTmp->size[1]] = d;
            }

            b_loop_ub = r1->size[1];
            for (i = 0; i < b_loop_ub; i++) {
              b_gTmp->data[(i + gTmp->size[1]) + loop_ub] = r1->data[i];
            }

            i = gTmp->size[0] * gTmp->size[1];
            gTmp->size[0] = 1;
            gTmp->size[1] = b_gTmp->size[1];
            emxEnsureCapacity_real_T(gTmp, i);
            loop_ub = b_gTmp->size[0] * b_gTmp->size[1];
            for (i = 0; i < loop_ub; i++) {
              gTmp->data[i] = b_gTmp->data[i];
            }

            sRamp = trapz(gTmp) * gmr * tGrad * 1.0E-5;
            guard1 = false;
            if (!epiType) {
              if (sRamp >= kRead) {
                i = gFlyb->size[0] * gFlyb->size[1];
                gFlyb->size[0] = 1;
                gFlyb->size[1] = gTmp->size[1];
                emxEnsureCapacity_real_T(gFlyb, i);
                loop_ub = gTmp->size[0] * gTmp->size[1];
                for (i = 0; i < loop_ub; i++) {
                  gFlyb->data[i] = gTmp->data[i] * kRead / sRamp;
                }

                exitg2 = 1;
              } else {
                guard1 = true;
              }
            } else {
              ex = (kMax + kRamp) + dk;
              if (sRamp >= ex) {
                i = gPre->size[0] * gPre->size[1];
                gPre->size[0] = 1;
                gPre->size[1] = gTmp->size[1];
                emxEnsureCapacity_real_T(gPre, i);
                loop_ub = gTmp->size[0] * gTmp->size[1];
                for (i = 0; i < loop_ub; i++) {
                  gPre->data[i] = gTmp->data[i] * ex / sRamp;
                }

                exitg2 = 1;
              } else {
                guard1 = true;
              }
            }

            if (guard1) {
              n = gTmp->size[1];
              if (gTmp->size[1] <= 2) {
                if (gTmp->size[1] == 1) {
                  ex = gTmp->data[0];
                } else if (gTmp->data[0] < gTmp->data[1]) {
                  ex = gTmp->data[1];
                } else {
                  ex = gTmp->data[0];
                }
              } else {
                ex = gTmp->data[0];
                for (b_k = 2; b_k <= n; b_k++) {
                  d = gTmp->data[b_k - 1];
                  if (ex < d) {
                    ex = d;
                  }
                }
              }

              if (ex < Gmax) {
                kBlip += 0.5;
              }

              fRamp++;
            }
          }
        } while (exitg2 == 0);

        if (exitg2 == 1) {
          if (!epiType) {
            printf("\n\t%g us flyback gradient, %g acquisition points\n",
                   (double)gFlyb->size[1] * tGrad, 2.0 * kBlip + fRamp);
            fflush(stdout);
          } else {
            printf("\n\t%g us prephasing gradient, %g acquisition points\n",
                   (double)gPre->size[1] * tGrad, 2.0 * kBlip + fRamp);
            fflush(stdout);
          }

          /*  Design/validate prephasing gradients */
          /*  for flyback, scale flyback grad */
          /*  for symmetric, scale prephasing grad */
          if (!epiType) {
            ex = kMax - dk * ((double)mtx - (double)mty);
            i = gPreY->size[0] * gPreY->size[1];
            gPreY->size[0] = 1;
            gPreY->size[1] = gFlyb->size[1];
            emxEnsureCapacity_real_T(gPreY, i);
            loop_ub = gFlyb->size[0] * gFlyb->size[1];
            for (i = 0; i < loop_ub; i++) {
              gPreY->data[i] = gFlyb->data[i] * ex / kRead;
            }

            ex = (kMax + kRamp) + dk;
            i = gPreX->size[0] * gPreX->size[1];
            gPreX->size[0] = 1;
            gPreX->size[1] = gFlyb->size[1];
            emxEnsureCapacity_real_T(gPreX, i);
            loop_ub = gFlyb->size[0] * gFlyb->size[1];
            for (i = 0; i < loop_ub; i++) {
              gPreX->data[i] = gFlyb->data[i] * ex / kRead;
            }

            i = gBlip->size[0] * gBlip->size[1];
            gBlip->size[0] = 1;
            gBlip->size[1] = gFlyb->size[1];
            emxEnsureCapacity_real_T(gBlip, i);
            loop_ub = gFlyb->size[0] * gFlyb->size[1];
            for (i = 0; i < loop_ub; i++) {
              gBlip->data[i] = gFlyb->data[i] * dk / kRead;
            }

            nEchoes = mty;
          } else {
            ex = kMax - dk * ((double)mtx - (double)mty);
            i = gPreY->size[0] * gPreY->size[1];
            gPreY->size[0] = 1;
            gPreY->size[1] = gPre->size[1];
            emxEnsureCapacity_real_T(gPreY, i);
            sRamp = (kMax + kRamp) + dk;
            loop_ub = gPre->size[0] * gPre->size[1];
            for (i = 0; i < loop_ub; i++) {
              gPreY->data[i] = gPre->data[i] * ex / sRamp;
            }

            i = gPreX->size[0] * gPreX->size[1];
            gPreX->size[0] = 1;
            gPreX->size[1] = gPre->size[1];
            emxEnsureCapacity_real_T(gPreX, i);
            loop_ub = gPre->size[0] * gPre->size[1];
            for (i = 0; i < loop_ub; i++) {
              gPreX->data[i] = gPre->data[i];
            }

            nEchoes = (int)ceil((double)mty / 2.0);
          }

          /*  Adjust read prephasing gradient to correct acquisition k-values */
          /*  1 acq pt, is there a better value to use here? */
          if (2 > b_gRead->size[1]) {
            i = 0;
            i2 = 0;
          } else {
            i = 1;
            i2 = b_gRead->size[1];
          }

          b_i = gTmp->size[0] * gTmp->size[1];
          gTmp->size[0] = 1;
          loop_ub = (int)np;
          gTmp->size[1] = ((gPreX->size[1] + loop_ub) + i2) - i;
          emxEnsureCapacity_real_T(gTmp, b_i);
          b_loop_ub = gPreX->size[1];
          for (b_i = 0; b_i < b_loop_ub; b_i++) {
            gTmp->data[b_i] = -gPreX->data[b_i];
          }

          for (b_i = 0; b_i < loop_ub; b_i++) {
            gTmp->data[b_i + gPreX->size[1]] = 0.0;
          }

          b_loop_ub = i2 - i;
          for (i2 = 0; i2 < b_loop_ub; i2++) {
            gTmp->data[(i2 + gPreX->size[1]) + loop_ub] = b_gRead->data[i + i2];
          }

          if (gTmp->size[1] - 1 < 0) {
            t->size[0] = 1;
            t->size[1] = 0;
          } else {
            i = t->size[0] * t->size[1];
            t->size[0] = 1;
            t->size[1] = (int)((double)gTmp->size[1] - 1.0) + 1;
            emxEnsureCapacity_real_T(t, i);
            b_loop_ub = (int)((double)gTmp->size[1] - 1.0);
            for (i = 0; i <= b_loop_ub; i++) {
              t->data[i] = i;
            }
          }

          i = t->size[0] * t->size[1];
          i2 = t->size[0] * t->size[1];
          t->size[0] = 1;
          emxEnsureCapacity_real_T(t, i2);
          b_loop_ub = i - 1;
          for (i = 0; i <= b_loop_ub; i++) {
            t->data[i] *= tGrad;
          }

          ex = floor(t->data[t->size[1] - 1] / dw);
          if (ex < 0.0) {
            T->size[0] = 1;
            T->size[1] = 0;
          } else {
            i = T->size[0] * T->size[1];
            T->size[0] = 1;
            b_loop_ub = (int)ex;
            T->size[1] = b_loop_ub + 1;
            emxEnsureCapacity_real_T(T, i);
            for (i = 0; i <= b_loop_ub; i++) {
              T->data[i] = i;
            }
          }

          i = T->size[0] * T->size[1];
          i2 = T->size[0] * T->size[1];
          T->size[0] = 1;
          emxEnsureCapacity_real_T(T, i2);
          b_loop_ub = i - 1;
          for (i = 0; i <= b_loop_ub; i++) {
            T->data[i] *= dw;
          }

          cumtrapz(gTmp, kTmp);
          i = kTmp->size[0] * kTmp->size[1];
          i2 = kTmp->size[0] * kTmp->size[1];
          kTmp->size[0] = 1;
          emxEnsureCapacity_real_T(kTmp, i2);
          b_loop_ub = i - 1;
          for (i = 0; i <= b_loop_ub; i++) {
            kTmp->data[i] = kTmp->data[i] * gmr * tGrad * 1.0E-5;
          }

          if (!epiType) {
            /*  It works best to use 1/10 tol here for rounding KTmp. Not sure why... */
            interp1(t, kTmp, T, KTmp);
            i = KTmp->size[0] * KTmp->size[1];
            i2 = KTmp->size[0] * KTmp->size[1];
            KTmp->size[0] = 1;
            emxEnsureCapacity_real_T(KTmp, i2);
            b_loop_ub = i - 1;
            for (i = 0; i <= b_loop_ub; i++) {
              KTmp->data[i] = KTmp->data[i] / 0.001 * 10.0;
            }

            nIter = KTmp->size[1];
            for (b_k = 0; b_k < nIter; b_k++) {
              KTmp->data[b_k] = rt_roundd(KTmp->data[b_k]);
            }

            i = KTmp->size[0] * KTmp->size[1];
            i2 = KTmp->size[0] * KTmp->size[1];
            KTmp->size[0] = 1;
            emxEnsureCapacity_real_T(KTmp, i2);
            b_loop_ub = i - 1;
            for (i = 0; i <= b_loop_ub; i++) {
              KTmp->data[i] = KTmp->data[i] * 0.001 / 10.0;
            }
          } else {
            interp1(t, kTmp, T, KTmp);
          }

          i = r1->size[0] * r1->size[1];
          r1->size[0] = 1;
          r1->size[1] = KTmp->size[1];
          emxEnsureCapacity_real_T(r1, i);
          b_loop_ub = KTmp->size[0] * KTmp->size[1];
          for (i = 0; i < b_loop_ub; i++) {
            r1->data[i] = KTmp->data[i];
          }

          fliplr(r1);
          kBlip = kMax - dk;
          i = r3->size[0] * r3->size[1];
          r3->size[0] = 1;
          r3->size[1] = r1->size[1];
          emxEnsureCapacity_boolean_T(r3, i);
          b_loop_ub = r1->size[0] * r1->size[1];
          for (i = 0; i < b_loop_ub; i++) {
            r3->data[i] = (r1->data[i] <= kBlip);
          }

          eml_find(r3, tmp_data, tmp_size);
          i = gTmp->size[0] * gTmp->size[1];
          gTmp->size[0] = 1;
          gTmp->size[1] = tmp_size[1];
          emxEnsureCapacity_real_T(gTmp, i);
          b_loop_ub = tmp_size[0] * tmp_size[1];
          for (i = 0; i < b_loop_ub; i++) {
            gTmp->data[i] = (double)tmp_data[i] - 1.0;
          }

          ex = ((kMax + kRamp) + dk) + ((KTmp->data[(int)((double)KTmp->size[1]
            - gTmp->data[0]) - 1] - kMax) + dk);
          i = gPreX->size[0] * gPreX->size[1];
          i2 = gPreX->size[0] * gPreX->size[1];
          gPreX->size[0] = 1;
          emxEnsureCapacity_real_T(gPreX, i2);
          sRamp = (kMax + kRamp) + dk;
          b_loop_ub = i - 1;
          for (i = 0; i <= b_loop_ub; i++) {
            gPreX->data[i] = gPreX->data[i] * ex / sRamp;
          }

          /*  Assemble gradient waveform */
          /*  insert zeros between flybacks and readouts until k-values are correct */
          fRamp = 0.0;
          if ((dk == 0.0) || ((-kMax < kBlip) && (dk < 0.0)) || ((kBlip < -kMax)
               && (dk > 0.0))) {
            T->size[0] = 1;
            T->size[1] = 0;
          } else if ((floor(-kMax) == -kMax) && (floor(dk) == dk)) {
            i = T->size[0] * T->size[1];
            T->size[0] = 1;
            b_loop_ub = (int)floor((kBlip - (-kMax)) / dk);
            T->size[1] = b_loop_ub + 1;
            emxEnsureCapacity_real_T(T, i);
            for (i = 0; i <= b_loop_ub; i++) {
              T->data[i] = -kMax + dk * (double)i;
            }
          } else {
            eml_float_colon(-kMax, dk, kBlip, T);
          }

          printf("Assembling full gradient waveform");
          fflush(stdout);
          do {
            exitg1 = 0;
            printf(".");
            fflush(stdout);
            if (fRamp > dw / tGrad + 1.0) {
              fprintf(stderr,
                      "create2dEpi:idxAcq\nFailed to assemble gradient waveform. Can\'t find acquisition window.\n");
              fflush(stderr);
              b_fileManager();
              exitg1 = 1;
            } else {
              if (!epiType) {
                if (2 > b_gRead->size[1]) {
                  i = -2;
                  i2 = -1;
                } else {
                  i = -1;
                  i2 = b_gRead->size[1] - 1;
                }

                if (2 > gFlyb->size[1]) {
                  b_i = 0;
                  b_k = -1;
                } else {
                  b_i = 1;
                  b_k = gFlyb->size[1] - 1;
                }

                if (2 > gBlip->size[1]) {
                  iidx = -1;
                  i3 = -1;
                } else {
                  iidx = 0;
                  i3 = gBlip->size[1] - 1;
                }

                nIter = b_gTmp->size[0] * b_gTmp->size[1];
                b_gTmp->size[0] = 1;
                b_loop_ub = i2 - i;
                n = (int)fRamp;
                b_gTmp->size[1] = ((b_loop_ub + b_k) - b_i) + n;
                emxEnsureCapacity_real_T(b_gTmp, nIter);
                for (nIter = 0; nIter <= b_loop_ub - 2; nIter++) {
                  b_gTmp->data[nIter] = b_gRead->data[(i + nIter) + 2];
                }

                b_loop_ub = b_k - b_i;
                for (nIter = 0; nIter <= b_loop_ub; nIter++) {
                  b_gTmp->data[((nIter + i2) - i) - 1] = -gFlyb->data[b_i +
                    nIter];
                }

                for (nIter = 0; nIter < n; nIter++) {
                  b_gTmp->data[(((nIter + i2) - i) + b_k) - b_i] = 0.0;
                }

                repmat(b_gTmp, (double)mty - 1.0, r1);
                nIter = b_gRead->size[1] - 1;
                i = b_gTmp->size[0] * b_gTmp->size[1];
                b_gTmp->size[0] = 1;
                b_gTmp->size[1] = (((b_gRead->size[1] + i3) - iidx) + n) - 1;
                emxEnsureCapacity_real_T(b_gTmp, i);
                b_loop_ub = b_gRead->size[1] - 1;
                for (i = 0; i < b_loop_ub; i++) {
                  b_gTmp->data[i] = 0.0;
                }

                b_loop_ub = i3 - iidx;
                for (i = 0; i < b_loop_ub; i++) {
                  b_gTmp->data[i + nIter] = gBlip->data[(iidx + i) + 1];
                }

                for (i = 0; i < n; i++) {
                  b_gTmp->data[((i + nIter) + i3) - iidx] = 0.0;
                }

                repmat(b_gTmp, (double)mty - 1.0, KTmp);
                if (2 > b_gRead->size[1]) {
                  i = 0;
                  i2 = 0;
                } else {
                  i = 1;
                  i2 = b_gRead->size[1];
                }

                b_i = g->size[0] * g->size[1];
                g->size[0] = 1;
                g->size[1] = (((gPreX->size[1] + loop_ub) + r1->size[1]) + i2) -
                  i;
                emxEnsureCapacity_creal_T(g, b_i);
                b_loop_ub = gPreX->size[1];
                for (b_i = 0; b_i < b_loop_ub; b_i++) {
                  g->data[b_i].re = -gPreX->data[b_i];
                  g->data[b_i].im = 0.0 - gPreY->data[b_i];
                }

                for (b_i = 0; b_i < loop_ub; b_i++) {
                  g->data[b_i + gPreX->size[1]].re = 0.0;
                  g->data[b_i + gPreX->size[1]].im = 0.0;
                }

                b_loop_ub = r1->size[1];
                for (b_i = 0; b_i < b_loop_ub; b_i++) {
                  g->data[(b_i + gPreX->size[1]) + loop_ub].re = r1->data[b_i];
                  g->data[(b_i + gPreX->size[1]) + loop_ub].im = KTmp->data[b_i];
                }

                b_loop_ub = i2 - i;
                for (i2 = 0; i2 < b_loop_ub; i2++) {
                  g->data[((i2 + gPreX->size[1]) + loop_ub) + r1->size[1]].re =
                    b_gRead->data[i + i2];
                  g->data[((i2 + gPreX->size[1]) + loop_ub) + r1->size[1]].im =
                    0.0;
                }
              } else {
                nEchoes = (int)ceil((double)mty / 2.0);
                if (2 > b_gRead->size[1]) {
                  i = -1;
                  i2 = -1;
                } else {
                  i = 0;
                  i2 = b_gRead->size[1] - 1;
                }

                b_i = b_gTmp->size[0] * b_gTmp->size[1];
                b_gTmp->size[0] = 1;
                b_loop_ub = i2 - i;
                n = (int)fRamp;
                b_gTmp->size[1] = (b_loop_ub + b_gRead->size[1]) + n;
                emxEnsureCapacity_real_T(b_gTmp, b_i);
                for (b_i = 0; b_i < b_loop_ub; b_i++) {
                  b_gTmp->data[b_i] = b_gRead->data[(i + b_i) + 1];
                }

                b_loop_ub = b_gRead->size[1];
                for (b_i = 0; b_i < b_loop_ub; b_i++) {
                  b_gTmp->data[(b_i + i2) - i] = -b_gRead->data[b_i];
                }

                for (b_i = 0; b_i < n; b_i++) {
                  b_gTmp->data[((b_i + i2) - i) + b_gRead->size[1]] = 0.0;
                }

                repmat(b_gTmp, floor((double)mty / 2.0), gTmp);
                if (rt_remd(mty, 2.0) != 0.0) {
                  if (2 > b_gRead->size[1]) {
                    i = -1;
                    i2 = 0;
                  } else {
                    i = 0;
                    i2 = b_gRead->size[1];
                  }

                  b_i = gTmp->size[1];
                  b_loop_ub = (i2 - i) - 2;
                  i2 = gTmp->size[0] * gTmp->size[1];
                  gTmp->size[1] = (gTmp->size[1] + b_loop_ub) + 1;
                  emxEnsureCapacity_real_T(gTmp, i2);
                  for (i2 = 0; i2 <= b_loop_ub; i2++) {
                    gTmp->data[b_i + i2] = b_gRead->data[(i + i2) + 1];
                  }

                  nIter = (int)(nRead + fRamp);
                  i = b_gTmp->size[0] * b_gTmp->size[1];
                  b_gTmp->size[0] = 1;
                  b_loop_ub = (int)(nRead - 1.0);
                  b_gTmp->size[1] = ((nIter + gBlip->size[1]) + b_loop_ub) +
                    gBlip->size[1];
                  emxEnsureCapacity_real_T(b_gTmp, i);
                  for (i = 0; i < nIter; i++) {
                    b_gTmp->data[i] = 0.0;
                  }

                  n = gBlip->size[1];
                  for (i = 0; i < n; i++) {
                    b_gTmp->data[i + nIter] = gBlip->data[i];
                  }

                  for (i = 0; i < b_loop_ub; i++) {
                    b_gTmp->data[(i + nIter) + gBlip->size[1]] = 0.0;
                  }

                  n = gBlip->size[1];
                  for (i = 0; i < n; i++) {
                    b_gTmp->data[((i + nIter) + gBlip->size[1]) + b_loop_ub] =
                      gBlip->data[i];
                  }

                  repmat(b_gTmp, floor(((double)mty - 1.0) / 2.0) - 1.0, r1);
                  nIter = (int)(nRamp + nRead);
                  i = kTmp->size[0] * kTmp->size[1];
                  kTmp->size[0] = 1;
                  kTmp->size[1] = (nIter + gBlip->size[1]) + r1->size[1];
                  emxEnsureCapacity_real_T(kTmp, i);
                  for (i = 0; i < nIter; i++) {
                    kTmp->data[i] = 0.0;
                  }

                  n = gBlip->size[1];
                  for (i = 0; i < n; i++) {
                    kTmp->data[i + nIter] = gBlip->data[i];
                  }

                  n = r1->size[1];
                  for (i = 0; i < n; i++) {
                    kTmp->data[(i + nIter) + gBlip->size[1]] = r1->data[i];
                  }

                  i = kTmp->size[1];
                  n = gBlip->size[1];
                  i2 = kTmp->size[0] * kTmp->size[1];
                  kTmp->size[1] = (kTmp->size[1] + b_loop_ub) + gBlip->size[1];
                  emxEnsureCapacity_real_T(kTmp, i2);
                  for (i2 = 0; i2 < b_loop_ub; i2++) {
                    kTmp->data[i + i2] = 0.0;
                  }

                  for (i2 = 0; i2 < n; i2++) {
                    kTmp->data[(i + i2) + b_loop_ub] = gBlip->data[i2];
                  }
                } else {
                  nIter = (int)(nRead + fRamp);
                  i = b_gTmp->size[0] * b_gTmp->size[1];
                  b_gTmp->size[0] = 1;
                  b_loop_ub = (int)(nRead - 1.0);
                  b_gTmp->size[1] = ((nIter + gBlip->size[1]) + b_loop_ub) +
                    gBlip->size[1];
                  emxEnsureCapacity_real_T(b_gTmp, i);
                  for (i = 0; i < nIter; i++) {
                    b_gTmp->data[i] = 0.0;
                  }

                  n = gBlip->size[1];
                  for (i = 0; i < n; i++) {
                    b_gTmp->data[i + nIter] = gBlip->data[i];
                  }

                  for (i = 0; i < b_loop_ub; i++) {
                    b_gTmp->data[(i + nIter) + gBlip->size[1]] = 0.0;
                  }

                  n = gBlip->size[1];
                  for (i = 0; i < n; i++) {
                    b_gTmp->data[((i + nIter) + gBlip->size[1]) + b_loop_ub] =
                      gBlip->data[i];
                  }

                  repmat(b_gTmp, floor(((double)mty - 1.0) / 2.0), r1);
                  nIter = (int)(nRamp + nRead);
                  i = kTmp->size[0] * kTmp->size[1];
                  kTmp->size[0] = 1;
                  kTmp->size[1] = (nIter + gBlip->size[1]) + r1->size[1];
                  emxEnsureCapacity_real_T(kTmp, i);
                  for (i = 0; i < nIter; i++) {
                    kTmp->data[i] = 0.0;
                  }

                  b_loop_ub = gBlip->size[1];
                  for (i = 0; i < b_loop_ub; i++) {
                    kTmp->data[i + nIter] = gBlip->data[i];
                  }

                  b_loop_ub = r1->size[1];
                  for (i = 0; i < b_loop_ub; i++) {
                    kTmp->data[(i + nIter) + gBlip->size[1]] = r1->data[i];
                  }
                }

                nIter = gTmp->size[1] - kTmp->size[1];
                i = KTmp->size[0] * KTmp->size[1];
                KTmp->size[0] = 1;
                KTmp->size[1] = kTmp->size[1] + nIter;
                emxEnsureCapacity_real_T(KTmp, i);
                b_loop_ub = kTmp->size[1];
                for (i = 0; i < b_loop_ub; i++) {
                  KTmp->data[i] = kTmp->data[i];
                }

                for (i = 0; i < nIter; i++) {
                  KTmp->data[i + kTmp->size[1]] = 0.0;
                }

                i = g->size[0] * g->size[1];
                g->size[0] = 1;
                g->size[1] = (gPreX->size[1] + loop_ub) + gTmp->size[1];
                emxEnsureCapacity_creal_T(g, i);
                b_loop_ub = gPreX->size[1];
                for (i = 0; i < b_loop_ub; i++) {
                  g->data[i].re = -gPreX->data[i];
                  g->data[i].im = 0.0 - gPreY->data[i];
                }

                for (i = 0; i < loop_ub; i++) {
                  g->data[i + gPreX->size[1]].re = 0.0;
                  g->data[i + gPreX->size[1]].im = 0.0;
                }

                b_loop_ub = gTmp->size[1];
                for (i = 0; i < b_loop_ub; i++) {
                  g->data[(i + gPreX->size[1]) + loop_ub].re = gTmp->data[i];
                  g->data[(i + gPreX->size[1]) + loop_ub].im = KTmp->data[i];
                }
              }

              if (g->size[1] - 1 < 0) {
                t->size[0] = 1;
                t->size[1] = 0;
              } else {
                i = t->size[0] * t->size[1];
                t->size[0] = 1;
                t->size[1] = (int)((double)g->size[1] - 1.0) + 1;
                emxEnsureCapacity_real_T(t, i);
                b_loop_ub = (int)((double)g->size[1] - 1.0);
                for (i = 0; i <= b_loop_ub; i++) {
                  t->data[i] = i;
                }
              }

              sRamp = 0.001 * tGrad;
              i = t->size[0] * t->size[1];
              i2 = t->size[0] * t->size[1];
              t->size[0] = 1;
              emxEnsureCapacity_real_T(t, i2);
              b_loop_ub = i - 1;
              for (i = 0; i <= b_loop_ub; i++) {
                t->data[i] *= sRamp;
              }

              ex = floor(1000.0 * t->data[t->size[1] - 1] / dw);
              if (ex < 0.0) {
                gTmp->size[0] = 1;
                gTmp->size[1] = 0;
              } else {
                i = gTmp->size[0] * gTmp->size[1];
                gTmp->size[0] = 1;
                b_loop_ub = (int)ex;
                gTmp->size[1] = b_loop_ub + 1;
                emxEnsureCapacity_real_T(gTmp, i);
                for (i = 0; i <= b_loop_ub; i++) {
                  gTmp->data[i] = i;
                }
              }

              sRamp = 0.001 * dw;
              b_cumtrapz(g, k);
              i = k->size[0] * k->size[1];
              i2 = k->size[0] * k->size[1];
              k->size[0] = 1;
              emxEnsureCapacity_creal_T(k, i2);
              b_loop_ub = i - 1;
              for (i = 0; i <= b_loop_ub; i++) {
                k->data[i].re = 1.0E-5 * (tGrad * (gmr * k->data[i].re));
                k->data[i].im = 1.0E-5 * (tGrad * (gmr * k->data[i].im));
              }

              i = KTmp->size[0] * KTmp->size[1];
              KTmp->size[0] = 1;
              KTmp->size[1] = k->size[1];
              emxEnsureCapacity_real_T(KTmp, i);
              b_loop_ub = k->size[0] * k->size[1];
              for (i = 0; i < b_loop_ub; i++) {
                KTmp->data[i] = k->data[i].re;
              }

              i = b_gTmp->size[0] * b_gTmp->size[1];
              b_gTmp->size[0] = 1;
              b_gTmp->size[1] = gTmp->size[1];
              emxEnsureCapacity_real_T(b_gTmp, i);
              b_loop_ub = gTmp->size[0] * gTmp->size[1];
              for (i = 0; i < b_loop_ub; i++) {
                b_gTmp->data[i] = sRamp * gTmp->data[i];
              }

              interp1(t, KTmp, b_gTmp, kTmp);
              i = kTmp->size[0] * kTmp->size[1];
              i2 = kTmp->size[0] * kTmp->size[1];
              kTmp->size[0] = 1;
              emxEnsureCapacity_real_T(kTmp, i2);
              b_loop_ub = i - 1;
              for (i = 0; i <= b_loop_ub; i++) {
                kTmp->data[i] /= 0.001;
              }

              nIter = kTmp->size[1];
              for (b_k = 0; b_k < nIter; b_k++) {
                kTmp->data[b_k] = rt_roundd(kTmp->data[b_k]);
              }

              i = gTmp->size[0] * gTmp->size[1];
              gTmp->size[0] = 1;
              gTmp->size[1] = T->size[1];
              emxEnsureCapacity_real_T(gTmp, i);
              b_loop_ub = T->size[0] * T->size[1];
              for (i = 0; i < b_loop_ub; i++) {
                gTmp->data[i] = T->data[i] / 0.001;
              }

              nIter = gTmp->size[1];
              for (b_k = 0; b_k < nIter; b_k++) {
                gTmp->data[b_k] = rt_roundd(gTmp->data[b_k]);
              }

              i = KTmp->size[0] * KTmp->size[1];
              KTmp->size[0] = 1;
              KTmp->size[1] = kTmp->size[1];
              emxEnsureCapacity_real_T(KTmp, i);
              b_loop_ub = kTmp->size[0] * kTmp->size[1];
              for (i = 0; i < b_loop_ub; i++) {
                KTmp->data[i] = kTmp->data[i] * 0.001;
              }

              i = b_gTmp->size[0] * b_gTmp->size[1];
              b_gTmp->size[0] = 1;
              b_gTmp->size[1] = gTmp->size[1];
              emxEnsureCapacity_real_T(b_gTmp, i);
              b_loop_ub = gTmp->size[0] * gTmp->size[1];
              for (i = 0; i < b_loop_ub; i++) {
                b_gTmp->data[i] = gTmp->data[i] * 0.001;
              }

              findPattern2(KTmp, b_gTmp, gTmp);
              n = 0;
              i = gTmp->size[1];
              for (b_k = 0; b_k < i; b_k++) {
                if (gTmp->data[b_k] != 0.0) {
                  n++;
                }
              }

              if (n < nEchoes) {
                fRamp++;
              } else {
                b_k = gTmp->size[1] - 1;
                nIter = 0;
                for (b_i = 0; b_i <= b_k; b_i++) {
                  if (gTmp->data[b_i] > 0.0) {
                    nIter++;
                  }
                }

                i = idxADC->size[0] * idxADC->size[1];
                idxADC->size[0] = 1;
                idxADC->size[1] = nIter;
                emxEnsureCapacity_real_T(idxADC, i);
                n = 0;
                for (b_i = 0; b_i <= b_k; b_i++) {
                  d = gTmp->data[b_i];
                  if (d > 0.0) {
                    idxADC->data[n] = d;
                    n++;
                  }
                }

                if (!epiType) {
                  ex = k->data[k->size[1] - 1].re;
                  gRead = k->data[k->size[1] - 1].im;
                  i = g->size[1];
                  loop_ub = gFlyb->size[1];
                  i2 = g->size[0] * g->size[1];
                  g->size[1] += gFlyb->size[1];
                  emxEnsureCapacity_creal_T(g, i2);
                  for (i2 = 0; i2 < loop_ub; i2++) {
                    d = gFlyb->data[i2];
                    sRamp = d * gRead;
                    if (sRamp == 0.0) {
                      kBlip = 0.0 / kRead;
                      sRamp = 0.0;
                    } else {
                      kBlip = 0.0;
                      sRamp /= kRead;
                    }

                    b_i = i + i2;
                    g->data[b_i].re = -d * ex / kRead - kBlip;
                    g->data[b_i].im = 0.0 - sRamp;
                  }
                } else {
                  ex = k->data[k->size[1] - 1].re;
                  gRead = k->data[k->size[1] - 1].im;
                  fRamp = (kMax + kRamp) + dk;
                  i = g->size[1];
                  loop_ub = gPre->size[1];
                  i2 = g->size[0] * g->size[1];
                  g->size[1] += gPre->size[1];
                  emxEnsureCapacity_creal_T(g, i2);
                  kMax = (kMax + kRamp) + dk;
                  for (i2 = 0; i2 < loop_ub; i2++) {
                    d = gPre->data[i2];
                    sRamp = d * gRead;
                    if (sRamp == 0.0) {
                      kBlip = 0.0 / fRamp;
                      sRamp = 0.0;
                    } else {
                      kBlip = 0.0;
                      sRamp /= fRamp;
                    }

                    b_i = i + i2;
                    g->data[b_i].re = -d * ex / kMax - kBlip;
                    g->data[b_i].im = 0.0 - sRamp;
                  }
                }

                if (g->size[1] - 1 < 0) {
                  t->size[0] = 1;
                  t->size[1] = 0;
                } else {
                  i = t->size[0] * t->size[1];
                  t->size[0] = 1;
                  t->size[1] = (int)((double)g->size[1] - 1.0) + 1;
                  emxEnsureCapacity_real_T(t, i);
                  loop_ub = (int)((double)g->size[1] - 1.0);
                  for (i = 0; i <= loop_ub; i++) {
                    t->data[i] = i;
                  }
                }

                sRamp = 0.001 * tGrad;
                i = t->size[0] * t->size[1];
                i2 = t->size[0] * t->size[1];
                t->size[0] = 1;
                emxEnsureCapacity_real_T(t, i2);
                loop_ub = i - 1;
                for (i = 0; i <= loop_ub; i++) {
                  t->data[i] *= sRamp;
                }

                ex = floor(1000.0 * t->data[t->size[1] - 1] / dw);
                if (ex < 0.0) {
                  T->size[0] = 1;
                  T->size[1] = 0;
                } else {
                  i = T->size[0] * T->size[1];
                  T->size[0] = 1;
                  loop_ub = (int)ex;
                  T->size[1] = loop_ub + 1;
                  emxEnsureCapacity_real_T(T, i);
                  for (i = 0; i <= loop_ub; i++) {
                    T->data[i] = i;
                  }
                }

                sRamp = 0.001 * dw;
                i = T->size[0] * T->size[1];
                i2 = T->size[0] * T->size[1];
                T->size[0] = 1;
                emxEnsureCapacity_real_T(T, i2);
                loop_ub = i - 1;
                for (i = 0; i <= loop_ub; i++) {
                  T->data[i] *= sRamp;
                }

                i = b_gTmp->size[0] * b_gTmp->size[1];
                b_gTmp->size[0] = 1;
                b_gTmp->size[1] = g->size[1];
                emxEnsureCapacity_real_T(b_gTmp, i);
                loop_ub = g->size[0] * g->size[1];
                for (i = 0; i < loop_ub; i++) {
                  b_gTmp->data[i] = g->data[i].re;
                }

                cumtrapz(b_gTmp, r1);
                i = b_gTmp->size[0] * b_gTmp->size[1];
                b_gTmp->size[0] = 1;
                b_gTmp->size[1] = g->size[1];
                emxEnsureCapacity_real_T(b_gTmp, i);
                loop_ub = g->size[0] * g->size[1];
                for (i = 0; i < loop_ub; i++) {
                  b_gTmp->data[i] = g->data[i].im;
                }

                cumtrapz(b_gTmp, gTmp);
                i = k->size[0] * k->size[1];
                k->size[0] = 1;
                k->size[1] = r1->size[1];
                emxEnsureCapacity_creal_T(k, i);
                loop_ub = r1->size[0] * r1->size[1];
                for (i = 0; i < loop_ub; i++) {
                  k->data[i].re = 1.0E-5 * (tGrad * (gmr * r1->data[i]));
                  k->data[i].im = 1.0E-5 * (tGrad * (gmr * gTmp->data[i]));
                }

                unnamed_idx_1 = (unsigned int)T->size[1];
                i = idxAcq->size[0] * idxAcq->size[1];
                idxAcq->size[0] = 1;
                idxAcq->size[1] = (int)unnamed_idx_1;
                emxEnsureCapacity_boolean_T(idxAcq, i);
                loop_ub = (int)unnamed_idx_1;
                for (i = 0; i < loop_ub; i++) {
                  idxAcq->data[i] = false;
                }

                diff(idxADC, r1);
                unique_vector(r1, gTmp);
                if (gTmp->size[1] > 1) {
                  fprintf(stderr,
                          "create2dEpi:idxDiff\nEcho spacing is not uniform in time.");
                  fflush(stderr);
                  b_fileManager();
                } else {
                  for (nIter = 0; nIter < nEchoes; nIter++) {
                    d = idxADC->data[nIter];
                    sRamp = (d + (double)mtx) - 1.0;
                    if (d > sRamp) {
                      i = -1;
                      i2 = 0;
                    } else {
                      i = (int)d - 2;
                      i2 = (int)sRamp;
                    }

                    loop_ub = (i2 - i) - 1;
                    for (i2 = 0; i2 < loop_ub; i2++) {
                      idxAcq->data[(i + i2) + 1] = true;
                    }

                    if (epiType) {
                      d = idxADC->data[nIter] + ceil(gTmp->data[0] / 2.0);
                      sRamp = (d + (double)mtx) - 1.0;
                      if (d > sRamp) {
                        i = -1;
                        i2 = 0;
                      } else {
                        i = (int)d - 2;
                        i2 = (int)sRamp;
                      }

                      loop_ub = (i2 - i) - 1;
                      for (i2 = 0; i2 < loop_ub; i2++) {
                        idxAcq->data[(i + i2) + 1] = true;
                      }
                    }
                  }

                  if (epiType && (rt_remd(mty, 2.0) != 0.0)) {
                    d = idxADC->data[nEchoes - 1] + ceil(gTmp->data[0] / 2.0);
                    sRamp = (d + (double)mtx) - 1.0;
                    if (d > sRamp) {
                      i = -1;
                      i2 = 0;
                    } else {
                      i = (int)d - 2;
                      i2 = (int)sRamp;
                    }

                    loop_ub = (i2 - i) - 1;
                    for (i2 = 0; i2 < loop_ub; i2++) {
                      idxAcq->data[(i + i2) + 1] = false;
                    }

                    i = r2->size[0] * r2->size[1];
                    r2->size[0] = 1;
                    r2->size[1] = idxAcq->size[1] - T->size[1];
                    emxEnsureCapacity_int32_T(r2, i);
                    loop_ub = idxAcq->size[1] - T->size[1];
                    for (i = 0; i < loop_ub; i++) {
                      r2->data[i] = (T->size[1] + i) + 1;
                    }

                    b_nullAssignment(idxAcq, r2);
                  }

                  b_k = idxAcq->size[1] - 1;
                  nIter = 0;
                  for (b_i = 0; b_i <= b_k; b_i++) {
                    if (idxAcq->data[b_i]) {
                      nIter++;
                    }
                  }

                  i = tAcq->size[0] * tAcq->size[1];
                  tAcq->size[0] = 1;
                  tAcq->size[1] = nIter;
                  emxEnsureCapacity_real_T(tAcq, i);
                  n = 0;
                  for (b_i = 0; b_i <= b_k; b_i++) {
                    if (idxAcq->data[b_i]) {
                      tAcq->data[n] = T->data[b_i];
                      n++;
                    }
                  }

                  i = KTmp->size[0] * KTmp->size[1];
                  KTmp->size[0] = 1;
                  KTmp->size[1] = k->size[1];
                  emxEnsureCapacity_real_T(KTmp, i);
                  loop_ub = k->size[0] * k->size[1];
                  for (i = 0; i < loop_ub; i++) {
                    KTmp->data[i] = k->data[i].re;
                  }

                  interp1(t, KTmp, tAcq, r1);
                  i = KTmp->size[0] * KTmp->size[1];
                  KTmp->size[0] = 1;
                  KTmp->size[1] = k->size[1];
                  emxEnsureCapacity_real_T(KTmp, i);
                  loop_ub = k->size[0] * k->size[1];
                  for (i = 0; i < loop_ub; i++) {
                    KTmp->data[i] = k->data[i].im;
                  }

                  interp1(t, KTmp, tAcq, gTmp);
                  i = k->size[0] * k->size[1];
                  k->size[0] = 1;
                  k->size[1] = r1->size[1];
                  emxEnsureCapacity_creal_T(k, i);
                  loop_ub = r1->size[0] * r1->size[1];
                  for (i = 0; i < loop_ub; i++) {
                    k->data[i].re = r1->data[i];
                    k->data[i].im = gTmp->data[i];
                  }

                  nIter = k->size[1];
                  i = gTmp->size[0] * gTmp->size[1];
                  gTmp->size[0] = 1;
                  gTmp->size[1] = k->size[1];
                  emxEnsureCapacity_real_T(gTmp, i);
                  for (b_k = 0; b_k < nIter; b_k++) {
                    gTmp->data[b_k] = rt_hypotd(k->data[b_k].re, k->data[b_k].im);
                  }

                  n = gTmp->size[1];
                  if (gTmp->size[1] <= 2) {
                    if (gTmp->size[1] == 1) {
                      iidx = 0;
                    } else {
                      iidx = (gTmp->data[0] > gTmp->data[1]);
                    }
                  } else {
                    ex = gTmp->data[0];
                    iidx = 0;
                    for (b_k = 2; b_k <= n; b_k++) {
                      d = gTmp->data[b_k - 1];
                      if (ex > d) {
                        ex = d;
                        iidx = b_k - 1;
                      }
                    }
                  }

                  *TE = tAcq->data[iidx];
                  if (!epiType) {
                    sRamp = tAcq->data[mtx] - tAcq->data[0];
                  } else {
                    if (1 > mtx) {
                      i = 0;
                    } else {
                      i = mtx;
                    }

                    i2 = gTmp->size[0] * gTmp->size[1];
                    gTmp->size[0] = 1;
                    gTmp->size[1] = i;
                    emxEnsureCapacity_real_T(gTmp, i2);
                    for (b_k = 0; b_k < i; b_k++) {
                      gTmp->data[b_k] = rt_hypotd(k->data[b_k].re, k->data[b_k].
                        im);
                    }

                    n = gTmp->size[1];
                    if (gTmp->size[1] <= 2) {
                      if (gTmp->size[1] == 1) {
                        b_i = 1;
                      } else if (gTmp->data[0] > gTmp->data[1]) {
                        b_i = 2;
                      } else {
                        b_i = 1;
                      }
                    } else {
                      ex = gTmp->data[0];
                      b_i = 1;
                      for (b_k = 2; b_k <= n; b_k++) {
                        d = gTmp->data[b_k - 1];
                        if (ex > d) {
                          ex = d;
                          b_i = b_k;
                        }
                      }
                    }

                    d = (double)mtx * 2.0;
                    if ((double)mtx + 1.0 > d) {
                      i = 0;
                      i2 = 0;
                    } else {
                      i = mtx;
                      i2 = (int)d;
                    }

                    n = i2 - i;
                    i2 = gTmp->size[0] * gTmp->size[1];
                    gTmp->size[0] = 1;
                    gTmp->size[1] = n;
                    emxEnsureCapacity_real_T(gTmp, i2);
                    for (b_k = 0; b_k < n; b_k++) {
                      nIter = i + b_k;
                      gTmp->data[b_k] = rt_hypotd(k->data[nIter].re, k->
                        data[nIter].im);
                    }

                    n = gTmp->size[1];
                    if (gTmp->size[1] <= 2) {
                      if (gTmp->size[1] == 1) {
                        nIter = 1;
                      } else if (gTmp->data[0] > gTmp->data[1]) {
                        nIter = 2;
                      } else {
                        nIter = 1;
                      }
                    } else {
                      ex = gTmp->data[0];
                      nIter = 1;
                      for (b_k = 2; b_k <= n; b_k++) {
                        d = gTmp->data[b_k - 1];
                        if (ex > d) {
                          ex = d;
                          nIter = b_k;
                        }
                      }
                    }

                    sRamp = tAcq->data[(nIter + mtx) - 1] - tAcq->data[b_i - 1];
                  }

                  kBlip = (double)tAcq->size[1] / (double)T->size[1];
                  print_processing(t->data[t->size[1] - 1], tAcq->data[iidx],
                                   sRamp, 1000.0 / sRamp, 100.0 * kBlip,
                                   validatedHoleFilling);
                  printf("\nTotal duration %g ms, TE %g ms\n\tESP %g ms, Phase Bandwidth %g Hz\n\tEncoding Efficiency %g%%\n",
                         validatedHoleFilling[0], validatedHoleFilling[1],
                         validatedHoleFilling[2], validatedHoleFilling[3],
                         validatedHoleFilling[4]);
                  fflush(stdout);

                  /*  Acquisition time values > 0 are recon samples */
                  b_k = idxAcq->size[1] - 1;
                  nIter = 0;
                  for (b_i = 0; b_i <= b_k; b_i++) {
                    if (!idxAcq->data[b_i]) {
                      nIter++;
                    }
                  }

                  i = r->size[0] * r->size[1];
                  r->size[0] = 1;
                  r->size[1] = nIter;
                  emxEnsureCapacity_int32_T(r, i);
                  n = 0;
                  for (b_i = 0; b_i <= b_k; b_i++) {
                    if (!idxAcq->data[b_i]) {
                      r->data[n] = b_i + 1;
                      n++;
                    }
                  }

                  loop_ub = r->size[0] * r->size[1];
                  i = b_T->size[0];
                  b_T->size[0] = loop_ub;
                  emxEnsureCapacity_real_T(b_T, i);
                  for (i = 0; i < loop_ub; i++) {
                    b_T->data[i] = -T->data[r->data[i] - 1];
                  }

                  loop_ub = b_T->size[0];
                  for (i = 0; i < loop_ub; i++) {
                    T->data[r->data[i] - 1] = b_T->data[i];
                  }

                  /*  Save acquisition time values to text for recon */
                  /*  if saveFile */
                  print_processing(t->data[t->size[1] - 1], tAcq->data[iidx],
                                   sRamp, 1000.0 / sRamp, 100.0 * kBlip,
                                   validatedHoleFilling);
                  b_NULL = NULL;
                  fileManager(fileid, &filestar, &autoflush);
                  if (!(filestar == b_NULL)) {
                    fprintf(filestar,
                            "# Waveform Parameters: TOTALTIME = %g [ms], ECHOTIME = %g [ms], ECHOSPACING = %g [ms], PHASEBW = %g [Hz], ENCODINGEFFICIENCY = %"
                            "g %%\n", validatedHoleFilling[0],
                            validatedHoleFilling[1], validatedHoleFilling[2],
                            validatedHoleFilling[3], validatedHoleFilling[4]);
                    if (autoflush) {
                      fflush(filestar);
                    }
                  }

                  i = T->size[1];
                  for (nIter = 0; nIter < i; nIter++) {
                    b_NULL = NULL;
                    fileManager(fileid, &filestar, &autoflush);
                    if (!(filestar == b_NULL)) {
                      fprintf(filestar, "%.6f \n", T->data[nIter]);
                      if (autoflush) {
                        fflush(filestar);
                      }
                    }
                  }

                  b_fileManager();

                  /*  end */
                  /*  Convert gradient waveform to 2 x n real array for output */
                  i = G->size[0] * G->size[1];
                  G->size[0] = 2;
                  G->size[1] = g->size[1];
                  emxEnsureCapacity_real_T(G, i);
                  loop_ub = g->size[1];
                  for (i = 0; i < loop_ub; i++) {
                    G->data[2 * i] = g->data[i].re;
                  }

                  loop_ub = g->size[1];
                  for (i = 0; i < loop_ub; i++) {
                    G->data[2 * i + 1] = g->data[i].im;
                  }

                  *NGRAD = g->size[1];
                  *NACQ = T->size[1];
                  printf("\nLeaving %s (%s)\n", "create2dEpi", "v20191017");
                  fflush(stdout);
                }

                exitg1 = 1;
              }
            }
          } while (exitg1 == 0);
        }
      }

      emxFree_real_T(&b_T);
      emxFree_boolean_T(&r3);
      emxFree_real_T(&b_gTmp);
      emxFree_int32_T(&r2);
      emxFree_real_T(&r1);
      emxFree_int32_T(&r);
      emxFree_real_T(&kTmp);
      emxFree_real_T(&b_gRead);
      emxFree_real_T(&tAcq);
      emxFree_boolean_T(&idxAcq);
      emxFree_real_T(&idxADC);
      emxFree_real_T(&KTmp);
      emxFree_real_T(&T);
      emxFree_real_T(&t);
      emxFree_real_T(&gPreX);
      emxFree_real_T(&gPreY);
      emxFree_real_T(&gPre);
      emxFree_real_T(&gFlyb);
      emxFree_real_T(&gTmp);
      emxFree_real_T(&gBlip);
      emxFree_creal_T(&k);
      emxFree_creal_T(&g);
    }
  }

  *PFY = pfy;
}

/* End of code generation (create2dEpi_codegen.c) */
