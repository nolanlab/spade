#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <cmath>
#include <stdint.h>
#include <limits>
#include <algorithm>
#include <vector>
#include <utility>
#include <functional>
#include <set>

#ifdef HAVE_PARALLEL_STL
#include <parallel/algorithm>
#define PARALLEL_NAMESPACE std::__parallel
#else
#define PARALLEL_NAMESPACE std
#endif

#include "prng.h"

namespace {

    typedef double Data_t;
    typedef float  Dist_t;
    typedef int    Idx_t;
   
    typedef std::vector<bool> BV_t;

    Dist_t 
    distance(const Data_t* a, const Data_t* b, size_t dim) {
	Dist_t d_l = 0.;
	for (size_t i=0; i<dim; i++)
	    d_l += fabs(a[i]-b[i]);
	return d_l;
    }

    // Assumes no diagonal entries ...
    size_t tri_idx(size_t r, size_t c, size_t n) { 
	return (r < c) ? (n-1)*r - r*(r-1)/2 + c-r-1 : (n-1)*c - c*(c-1)/2 + r-c-1; 
    }

    struct TriIdxFixedRow {
	size_t r, n, p;
	TriIdxFixedRow(size_t r_a, size_t n_a) : r(r_a), n(n_a) { p = (n-1)*r - r*(r-1)/2 - r - 1; }
	inline size_t tri_idx(size_t c) const { return (r < c) ? p + c : (n-1)*c - c*(c-1)/2 + r-c-1; }
    };


    void
    pairwise_distances(const Data_t *data, size_t dim, size_t obs, Dist_t *dist) {
	#pragma omp parallel for shared(data)
	for (size_t i=0; i<obs; i++) {
	    const Data_t *point = &data[i*dim];
	    size_t t = tri_idx(i,i+1,obs);
	    for (size_t j=i+1; j<obs; j++) {
		dist[t++] = distance(point, &data[j*dim], dim);
	    }
	}
    }
 
    struct LT {
	const Dist_t  *pair;
	size_t         obs; 
	size_t         idx;
	const BV_t&    merged;
	
	const TriIdxFixedRow ti;
	
	LT(const Dist_t* pair_a, size_t obs_a, size_t idx_a, const BV_t& merged_a) : 
	    pair(pair_a), obs(obs_a), idx(idx_a), merged(merged_a), ti(idx, obs) {}
	
	bool operator()(size_t l, size_t r) const {
	    if      (l == idx)  return false;
	    else if (r == idx)  return true;
	    else if (merged[l]) return false;
	    else if (merged[r]) return true;
	    else {
		return pair[ti.tri_idx(l)] < pair[ti.tri_idx(r)];
	    }
	}
    }; 

    void
    possible_merges(const Dist_t *pair, size_t obs, size_t fold, const BV_t& merged, 
		    size_t idx, size_t* idxs_begin, size_t* idxs_end) {
	LT lt(pair, obs, idx, merged);
	PARALLEL_NAMESPACE::partial_sort(idxs_begin, idxs_begin+fold, idxs_end, lt);   
    }

    void
    single_linkage(Dist_t *pair, size_t obs, size_t into, size_t from) {	
	TriIdxFixedRow tii(into,obs), tif(from,obs);
	#pragma omp parallel for shared(pair, tii, tif)
	for (size_t i=0; i<obs; i++) {
	    if (i == into || i == from)
		continue;
	    size_t into_t = tii.tri_idx(i), from_t = tif.tri_idx(i);
	    pair[into_t] = std::min(pair[into_t], pair[from_t]);
	    pair[from_t] = std::numeric_limits<Dist_t>::max();
	}
    }

    struct TF {
	typedef const size_t argument_type;
	typedef bool         result_type;

	const std::vector<bool>& pred;
	
	TF(const std::vector<bool>& pred_a) : pred(pred_a) {}
	bool operator()(const size_t a) const { return pred[a]; }
    };

    void
    cluster(const Data_t* data, size_t obs, size_t dim, Idx_t* merge) {   
	Dist_t *pair = Calloc(obs*(obs-1)/2,Dist_t); // Declare a triangular matrix for pairwise distances
	pairwise_distances(data, dim, obs, pair); 

	std::vector<bool> valid(obs, true), merged(obs, false);
	size_t *idxs = Calloc(obs,size_t); for (size_t i=0; i<obs; i++) { idxs[i] = i; }
	Idx_t  *clst = Calloc(obs,Idx_t);  for (Idx_t i=0;  i<obs; i++) { clst[i] = -i-1; }

	size_t merge_idx = 0, merge_round = 0;

	TF tf_v(valid), tf_m(merged);

	size_t *idxs_end = idxs + obs; // End of valid indices
	while (merge_idx < (obs-1)) {

	    Rprintf("Merging round %zu ...\n", merge_round++);

	    // Prepare for a round of merging
	    idxs_end = std::partition(idxs, idxs_end, tf_v);	// Resest range of valid indices 	    
	    std::random_shuffle(idxs, idxs_end);		// Randomize order of cluster assembly	    
	    std::fill(merged.begin(), merged.end(), false);	// Reset merge tracking
    	    
	    size_t fold = std::max((size_t)1 /* fold 2x*/, (size_t)(idxs_end-idxs)/5000);

	    while (true) {

		// Find first un-merged cluster 
		size_t *into_p = PARALLEL_NAMESPACE::find_if(idxs, idxs_end, std::not1(tf_m));
		if (into_p == idxs_end)
		    goto END_ROUND;
		size_t into = *into_p;

		merged[into] = true; // Mark as merged in this round

		// Find other clusters that could be merged...
		possible_merges(pair, obs, fold, merged, into, idxs, idxs_end);	

		for (size_t i=0; i<fold; i++) {
		    size_t from = idxs[i];
		    if (!valid[from] || merged[from]) 
			goto END_ROUND;  // We have exhausted possible merges in this round

		    merged[from] = true; valid[from] = false; // Removed from further consideration
		    single_linkage(pair, obs, into, from);    // Update distances	

		    // Update history of cluster merges
		    merge[merge_idx*2] = clst[into]; merge[merge_idx*2+1] = clst[from];
		    clst[into] = merge_idx++ + 1;
		}
	    }
END_ROUND:;
	}	

	Free(idxs);
	Free(clst); 
	Free(pair);
    }

    void
    live_clusters(const Idx_t* merge, size_t obs, size_t k, std::set<Idx_t>& clusters) { 
	clusters.insert((Idx_t)obs-1); // Initialize with "final" cluster
	for (size_t i=obs-2; i>=0; i--) {
	    if (clusters.size() >= k)
		break;
	    clusters.erase((Idx_t)(i+1)); // Recall cluster Ids of 1-indexed
	    clusters.insert(merge[2*i]);
	    clusters.insert(merge[2*i+1]);
	}
    }

    void
    assign_observation(const Idx_t* merge, Idx_t cluster, Idx_t my_assgn, Idx_t* assgn) {
	if (cluster < 0) {
	    assgn[-(cluster+1)] = my_assgn;
	} else {
	    assign_observation(merge, merge[2*(cluster-1)], my_assgn, assgn);
	    assign_observation(merge, merge[2*(cluster-1)+1], my_assgn, assgn);
	}
    }

} // anonymous namespace


extern "C" {

    SEXP 
    FSPD_cluster(SEXP tbl, SEXP k) {
	// tbl is column major order, transposed to be row major
	size_t obs = static_cast<size_t>(INTEGER(GET_DIM(tbl))[1]), 
	       dim = static_cast<size_t>(INTEGER(GET_DIM(tbl))[0]);
	size_t k_l = static_cast<size_t>(asInteger(k));

	SEXP ans, names, merge, assgn; int n_protected = 0;
	
	// Cluster observations
	PROTECT(merge = allocMatrix(INTSXP, 2, obs-1)); n_protected++;
	Idx_t *merge_l = static_cast<Idx_t*>(INTEGER(merge));

	cluster(static_cast<Data_t*>(REAL(tbl)), obs, dim, merge_l);


	// Assigning observations to clusters	
	PROTECT(assgn = allocVector(INTSXP, obs)); n_protected++; 	
	
	std::set<Idx_t> clusters;
	live_clusters(merge_l, obs, k_l, clusters);  // Cut tree to find ~k clusters
	
	Idx_t cluster_index = 1; 
	for (std::set<Idx_t>::iterator i=clusters.begin(),e=clusters.end(); i!=e; ++i)
	    assign_observation(merge_l, *i, cluster_index++, static_cast<Idx_t*>(INTEGER(assgn)));
	   

	// Assemble return object
	PROTECT(ans = allocVector(VECSXP,2)); n_protected++;
	PROTECT(names = allocVector(STRSXP,2)); n_protected++;

	SET_VECTOR_ELT(ans, 0, merge);
	SET_STRING_ELT(names, 0, mkChar("merge"));
	SET_VECTOR_ELT(ans, 1, assgn);
	SET_STRING_ELT(names, 1, mkChar("assgn"));

	setAttrib(ans, R_NamesSymbol, names);
	UNPROTECT(n_protected);


	return ans;
    }

}

