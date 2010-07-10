#ifndef UTIL_H
#define UTIL_H

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

// Replacement for assert
#define assert(X) { if (!(X)) { error("Assertion Failure %s [%s:%d]: %s\n",__FUNCTION__,__FILE__,__LINE__,#X); } }

#endif
