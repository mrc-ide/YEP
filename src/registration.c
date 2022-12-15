#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void SEIRV_Model_rhs_dde(void *);

/* .Call calls */
extern SEXP SEIRV_Model_contents(SEXP);
extern SEXP SEIRV_Model_create(SEXP);
extern SEXP SEIRV_Model_initial_conditions(SEXP, SEXP);
extern SEXP SEIRV_Model_metadata(SEXP);
extern SEXP SEIRV_Model_rhs_r(SEXP, SEXP, SEXP);
extern SEXP SEIRV_Model_set_initial(SEXP, SEXP, SEXP, SEXP);
extern SEXP SEIRV_Model_set_user(SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {"SEIRV_Model_rhs_dde", (DL_FUNC) &SEIRV_Model_rhs_dde, 1},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"SEIRV_Model_contents",           (DL_FUNC) &SEIRV_Model_contents,           1},
    {"SEIRV_Model_create",             (DL_FUNC) &SEIRV_Model_create,             1},
    {"SEIRV_Model_initial_conditions", (DL_FUNC) &SEIRV_Model_initial_conditions, 2},
    {"SEIRV_Model_metadata",           (DL_FUNC) &SEIRV_Model_metadata,           1},
    {"SEIRV_Model_rhs_r",              (DL_FUNC) &SEIRV_Model_rhs_r,              3},
    {"SEIRV_Model_set_initial",        (DL_FUNC) &SEIRV_Model_set_initial,        4},
    {"SEIRV_Model_set_user",           (DL_FUNC) &SEIRV_Model_set_user,           2},
    {NULL, NULL, 0}
};

void R_init_YEP(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
