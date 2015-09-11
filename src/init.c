#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "declarations.h"

static R_FortranMethodDef FortEntries[] = {
		{"fapproxinfmat", (DL_FUNC) &F77_SUB(approxinfmat), 4},
		{"farmainvcov", (DL_FUNC) &F77_SUB(armainvcov), 7},
		{"farcov", (DL_FUNC) &F77_SUB(arcov), 5},
		{NULL, NULL, 0}
};

void R_init_arimaPI(DllInfo *dll)
{
	R_registerRoutines(dll, NULL, NULL, FortEntries, NULL);
	R_useDynamicSymbols(dll, FALSE);
}
