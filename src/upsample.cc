#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <cmath>
#include <stdint.h>
#include <limits>
#include <algorithm>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "util.h"

namespace {

    typedef double Data_t;
    typedef double Dist_t;
    typedef int    Idx_t;
    
    Dist_t 
    distance(const Data_t* a, const Data_t* b, size_t dim) {
	Dist_t d_l = 0.;
	for (size_t i=0; i<dim; i++)
	    d_l += fabs(a[i]-b[i]);
	return d_l;
    }

    void
    assign_cluster( const Data_t* data, const Data_t* cluster_data, Idx_t* cluster_assign,
		    size_t obs, size_t cls, size_t dim, Idx_t* assign) {
#ifdef _OPENMP
	#pragma omp parallel for shared(assign)
#endif
	for (size_t i=0; i<obs; i++) {
	    const Data_t* pt = &data[i*dim];
	   
	    // Compute nearest cluster with brute force  
	    Idx_t min_idx = 0; Dist_t min_dist = std::numeric_limits<Dist_t>::max();
	    for (size_t j=0; j<cls; j++) {
		Dist_t d = distance(pt,&cluster_data[j*dim],dim);
		if (d < min_dist) {
		    min_idx = j; 
		    min_dist = d;
		}
	    }

	    assign[i] = cluster_assign[min_idx];
	}
    }


} // anonymous namespace

extern "C" {

    SEXP
    SPADE_assign(SEXP tbl, SEXP cluster_data, SEXP cluster_assign) {
	// tbl,clusters is column major order, transposed to be row major
	size_t obs = static_cast<size_t>(INTEGER(GET_DIM(tbl))[1]), 
	       dim = static_cast<size_t>(INTEGER(GET_DIM(tbl))[0]),
	       cls = static_cast<size_t>(INTEGER(GET_DIM(cluster_data))[1]);
	int n_protected = 0;
	
	SEXP assign;
	
	PROTECT(assign = allocVector(INTSXP, obs)); ++n_protected;
	Idx_t *assign_l = static_cast<Idx_t*>(INTEGER(assign));

	assign_cluster(	static_cast<Data_t*>(REAL(tbl)), 
			static_cast<Data_t*>(REAL(cluster_data)), 
			static_cast<Idx_t*>(INTEGER(cluster_assign)), 
			obs, cls, dim, assign_l );	

	UNPROTECT(n_protected);
	return assign;
    }

} // extern
