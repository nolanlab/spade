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
    typedef double Dist_t;
    typedef int    Idx_t;
   
 
    Dist_t 
    distance(const Data_t* a, const Data_t* b, size_t dim) {
	Dist_t d_l = 0.;
	for (size_t i=0; i<dim; i++)
	    d_l += fabs(a[i]-b[i]);
	return d_l;
    }

    typedef std::vector<bool> BV_t;

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

    struct Merge_t {
	Idx_t from, into;
	void set(Idx_t f, Idx_t i) { from = f; into = i; } 
    };
    struct MergeCMP {
	const Dist_t*  height;
	MergeCMP(const Dist_t* height_a) : height(height_a) {}
	bool operator()(const Merge_t& a, const Merge_t& b) const {
	    return height[a.from] < height[b.from];
	}
    };
    

    void
    slink(const Data_t* data, size_t obs, size_t dim, Merge_t* merge) {
	R::auto_ptr<Idx_t>  P_ap(Calloc(obs, Idx_t));
	R::auto_ptr<Dist_t> L_ap(Calloc(obs, Dist_t));
	R::auto_ptr<Dist_t> M_ap(Calloc(obs, Dist_t));

	Idx_t*  P = P_ap.get();
	Dist_t* L = L_ap.get();
	Dist_t* M = M_ap.get();

	// Compute the pointer representation
	for (size_t i=0; i<obs; i++) {

	    // Steps corresponding to Sibson's 1972 paper "SLINK: An optimally
	    // efficient algorithmf or the single-link cluster method"

	    // Step 1: Initialize
	    P[i] = i;
	    L[i] = std::numeric_limits<Dist_t>::max();

	    // Step 2: Build out pairwise distances from objects in pointer
	    // represenation to the new object
	    const Data_t* pt = &data[i*dim];
	    #pragma omp parallel for shared(pt, M)
	    for (size_t j=0; j<i; j++) {
		M[j] = distance(pt, &data[j*dim], dim);
	    }

	    // Step 3: Update M, P, L
	    for (size_t j=0; j<i; j++) {
		Dist_t l = L[j], m = M[j];
		if (l >= m) {
		    M[P[j]] = std::min(M[P[j]], l);
		    L[j]    = m;
		    P[j]    = i;
		} else {
		    M[P[j]] = std::min(M[P[j]], m);
		}
	    }
    
	    // Step 4: Actualize the clusters
	    for (size_t j=0; j<i; j++) {
		if (L[j] >= L[P[j]])
		    P[j] = i;
	    }

	}

	// Convert the pointer representation to dendogram 
	for (size_t i=0; i<(obs-1); i++) {
	    merge[i].set(i, P[i]);
	}
	std::sort((Merge_t*)merge, (Merge_t*)merge + obs-1, MergeCMP(L)); 
	for (size_t i=0; i<obs; i++) {
	    P[i] = -(Idx_t)(i+1); // R is 1-indexed
	}
	for (size_t i=0; i<(obs-1); i++) {
	    Idx_t into = merge[i].into;
	    merge[i].set(P[merge[i].from],  P[merge[i].into]);
	    P[into] = i+1;
	}
    }

    void
    flatten_tree(Merge_t* slink, size_t obs) {
	Idx_t cur_idx = 0;
	
	std::vector<Idx_t> live_clusters; 
	live_clusters.reserve(2*obs); // Need to pre-allocate enough to prevent reallocation
	
	// Initialize live clusters with "single observation" clusters in
	// single-linkage sorted order
	{
	    Idx_t *part_beg = (Idx_t*)slink, 
		  *part_end = std::stable_partition((Idx_t*)slink, (Idx_t*)slink + 2*(obs-1), std::bind2nd(std::less<Idx_t>(),0));
	    live_clusters.insert(live_clusters.end(), part_beg, part_end);
	}

	// Merge clusters
	while (live_clusters.size() > 1) {	    
	    size_t num_live = live_clusters.size(), fold = std::max((size_t)1 /* fold 2x*/, num_live / 5000);	    
	    Idx_t *part_beg = &live_clusters[0], 
		  *part_end = part_beg + num_live; 
	    
	    // Merge clusters while at least two remaining in this round
	    while ((part_beg+1) < part_end) {
		slink[cur_idx++].set(*(part_beg++), *(part_beg++)); // Merge initial pair
		for (size_t i=1; i < fold && part_beg < part_end; i++) {
		    slink[cur_idx].set(*(part_beg++), cur_idx);  // Merge next nearest cluster until "fold" is reached
		    ++cur_idx;
		}
		live_clusters.push_back(cur_idx); // Only last formed cluster remains live
	    }
	    live_clusters.insert(live_clusters.end(), part_beg, part_end); // Re-insert any leftover clusters
	    live_clusters.erase(live_clusters.begin(), live_clusters.begin()+num_live);
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

	assert(sizeof(Merge_t) == 2*sizeof(Idx_t));

	SEXP ans, names, merge, assgn; int n_protected = 0;
	
	// Cluster observations
	PROTECT(merge = allocMatrix(INTSXP, 2, obs-1)); n_protected++;
	Idx_t *merge_l = static_cast<Idx_t*>(INTEGER(merge));
	
	slink(static_cast<Data_t*>(REAL(tbl)), obs, dim, (Merge_t*)merge_l);
	flatten_tree((Merge_t*)merge_l, obs);
	

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

