#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern "C" {

    void R_init_spade(DllInfo *info) {
		static const R_CallMethodDef _CallFun[] = {
			// I don't think the functions are exported here because they're
			// all in extern "C" scopes.
			{NULL}
		};

    }

} // extern
