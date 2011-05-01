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
    typedef int    Count_t;
    
    Dist_t 
    distance(const Data_t* a, const Data_t* b, size_t dim) {
	Dist_t d_l = 0.;
	for (size_t i=0; i<dim; i++)
	    d_l += fabs(a[i]-b[i]);
	return d_l;
    }

  Dist_t compute_median_min_dist(Data_t *data, size_t dim, size_t obs, size_t num_samples) {					
		Dist_t* min_dists = new Dist_t[num_samples];

		// Create an randomly shuffled array of indices for
		// randomly selected samples for computing minimum median distance
		size_t* idxs = new size_t[obs];
		for (size_t i=0; i<obs; i++)  			
			idxs[i] = i;
		std::random_shuffle(idxs, idxs+obs);

		#ifdef _OPENMP
		#pragma omp parallel for shared(min_dists)
		#endif
		for (size_t i=0; i<num_samples; i++) {  	    
	   	size_t idx = idxs[i];  // Uniformly sampled cell
	    
			// Compute distance to all other neighbors (skipping myself)
	    Data_t *point = &data[idx*dim];
	    Dist_t min_l = std::numeric_limits<Dist_t>::max();
	    for (size_t j=0; j<idx; j++)
				min_l = std::min(min_l,distance(point, &data[j*dim], dim));
	    for (size_t j=idx+1; j<obs; j++)
				min_l = std::min(min_l,distance(point, &data[j*dim], dim));
	    min_dists[i] = min_l;	
		}

		// Determine median
		std::nth_element(min_dists, min_dists+num_samples/2, min_dists+num_samples);
		Dist_t median_min_dist = min_dists[num_samples/2];
	
		delete idxs;
		delete min_dists;
	
		return median_min_dist;
  }

    void
    count_neighbors(Data_t* data, size_t dim, size_t obs, Dist_t kernel_width, Dist_t apprx_width, Count_t* densities) 
    {
		std::fill(densities, densities+obs, 0);
#ifdef _OPENMP
		#pragma omp parallel for shared(densities)	
#endif
		for (size_t i=0; i<obs; i++) {	    
			if (densities[i] > 0)
				continue;

			std::vector<size_t> apprxs;  // Keep track on observations we can approximate
			Data_t *point = &data[i*dim];
			Count_t c = 0;

			for (size_t j=0; j<obs; j++) {
				Dist_t d = distance(point, &data[j*dim], dim);
				if (d < kernel_width)
					c++;
				if (d < apprx_width)
					apprxs.push_back(j);
			}

			// Potential race condition on other density entries, use atomic
			// update to be safe
			for (size_t j=0; j<apprxs.size(); j++) {
#ifdef _OPENMP
				__sync_bool_compare_and_swap(densities+apprxs[j],0,c); //densities[apprxs[j]] = c;
#else
				densities[apprxs[j]] = c;
#endif
			
			}
			densities[i] = c;
		}

    }	
        
} // anonymous namespace

extern "C" {

    SEXP 
    SPADE_density(SEXP tbl, SEXP kernel_mult, SEXP apprx_mult, SEXP med_samples) {
	// tbl is column major order, transposed to be row major
	size_t obs = static_cast<size_t>(INTEGER(GET_DIM(tbl))[1]), 
	       dim = static_cast<size_t>(INTEGER(GET_DIM(tbl))[0]);
	Data_t kernel_mult_l = static_cast<Data_t>(asReal(kernel_mult)),
	       apprx_mult_l  = static_cast<Data_t>(asReal(apprx_mult));
	size_t med_samples_l = static_cast<size_t>(asInteger(med_samples));	

	SEXP densities;
	PROTECT(densities = allocVector(INTSXP, obs));
	Count_t* densities_l = static_cast<Count_t*>(INTEGER(densities));

	Dist_t median_min_dist = compute_median_min_dist(static_cast<Data_t*>(REAL(tbl)), dim, obs, std::min(obs,med_samples_l)),
	       kernel_width    = kernel_mult_l * median_min_dist,
	       apprx_width     = apprx_mult_l  * median_min_dist;
	
	count_neighbors(static_cast<Data_t*>(REAL(tbl)), dim, obs, kernel_width, apprx_width, densities_l);	
	
	UNPROTECT(1);
	return densities;
    }


} // extern
