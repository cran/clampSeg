#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP clampSeg_deconvolveJump(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP clampSeg_deconvolvePeak(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"clampSeg_deconvolveJump", (DL_FUNC) &clampSeg_deconvolveJump,  8},
  {"clampSeg_deconvolvePeak", (DL_FUNC) &clampSeg_deconvolvePeak, 10},
  {NULL, NULL, 0}
};

void R_init_clampSeg(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
