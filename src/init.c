#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP _clampSeg_deconvolveJump(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _clampSeg_deconvolvePeak(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_clampSeg_deconvolveJump", (DL_FUNC) &_clampSeg_deconvolveJump,  8},
  {"_clampSeg_deconvolvePeak", (DL_FUNC) &_clampSeg_deconvolvePeak, 10},
  {NULL, NULL, 0}
};

void R_init_clampSeg(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
