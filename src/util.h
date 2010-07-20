#ifndef UTIL_H
#define UTIL_H

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

// Replacement for assert
#define assert(X) { if (!(X)) { error("Assertion Failure %s [%s:%d]: %s\n",__FUNCTION__,__FILE__,__LINE__,#X); } }

namespace R {
    
    template<class T>
    class auto_ptr {
    public:
	typedef T element_type;

	explicit auto_ptr(T* p=0) : P(p) {}
	~auto_ptr() { if (P) Free(P); }

	T* get() const { return P; }
	T* release() { T* P_l = P; P = 0; return P_l; }

    private:
	T* P;
    };
}

#endif
