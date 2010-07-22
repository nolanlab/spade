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
    
    /* Maintain state for a single cluster
     */
    class ACluster {
    public:
	static size_t dim;
	static void init_global(size_t d) { dim = d; }

    public:
	Data_t* center;		// My cluster center
	bool valid;		// Cluster statuses
	bool merged;
	Idx_t cluster_id;	// Current ID in merging data structure
	const Data_t **members;	// Observations in my cluster
       	size_t num_members;


	// Satisfied with compiler generated constructors

	// Initialization and destruction
	void init_RM(size_t idx, const Data_t* data) {
	    center = Calloc(dim, Data_t);
	    memcpy(center, &data[idx*dim], dim*sizeof(Data_t));

	    valid = true; merged = false;
	    cluster_id = -(Idx_t)(idx+1);  // Recall R is 1-indexed

	    members = Calloc(1, const Data_t*);
	    members[0]  = &data[idx*dim];
	    num_members = 1;
	}
	void destroy() {
	    assert(!valid);
	    center  = (Free(center), (Data_t*)0);
	    members = (Free(members), (const Data_t**)0);
	    num_members = 0;
	}

	// Getters & Setters (many used in STL algorithms)
	bool get_valid() const { return valid; }
	void set_valid(bool v) { valid = v; }
	
	bool get_merged() const { return merged; }
	void set_merged(bool m) { merged = m; }

	Idx_t get_cluster_id() const { return cluster_id; }
	void set_cluster_id(Idx_t i) { cluster_id = i; }	

	// Helper structs
	
	// LessThan by data column (used in Median)
	struct ColLT_RM {
	    size_t col;
	    ColLT_RM(size_t c) : col(c) {}
	    bool operator()(const Data_t* l, const Data_t* r) { return l[col] < r[col]; }
	};

	Dist_t slink_distance(const ACluster& a) const {
	    Dist_t d = std::numeric_limits<Dist_t>::max();
	    for (size_t i=0; i<a.num_members; i++)
		d = std::min(d, distance(center, a.members[i], dim));
	    return d;
	}

	// LessThan by "single" link distance between my center and other cluster's members
	struct SLinkLT_RM {
	    const ACluster& base;
	    SLinkLT_RM(const ACluster& b) : base(b) {}
	    bool operator()(const ACluster& l, const ACluster& r) const {
		if      (l.merged) return false;
		else if (r.merged) return true;
		else return base.slink_distance(l) < base.slink_distance(r);
    	    }
	};

	// Reset cluster, including updating cluster center
	void reset_RM() {
	    merged = false;  // Reset merged for next round
	    for (size_t i=0; i<dim; i++) { // Update the cluster center as median of members
		ColLT_RM lt(i);
		std::nth_element(members, members + num_members/2, members + num_members, lt);
		center[i] = members[num_members/2][i];
	    }
	}

	// Merge in other cluster
	void merge_in(ACluster& rhs) {
	    rhs.merged = true; rhs.valid = false;

	    members = Realloc(members, num_members+rhs.num_members, const Data_t*);
	    memcpy(members+num_members, rhs.members, rhs.num_members*sizeof(Data_t*));
	    num_members += rhs.num_members;
	    
	    rhs.destroy();    
	}

    };

    size_t ACluster::dim;

    void
    cluster(const Data_t* data, size_t obs, size_t dim, Merge_t* merge) {
	ACluster::init_global(dim);
	
	Idx_t cur_merge = 0;  // Track current merge step

	R::auto_ptr<ACluster> c_ap(Calloc(obs, ACluster));
	ACluster *c_beg = c_ap.get(), *c_end = c_beg + obs;
	for (size_t i=0; i<obs; i++)  // Initialize clusters for row major data
	    c_beg[i].init_RM(i, data);
	
	while (cur_merge < (obs-1)) {
	    // Only looking at "valid" clusters
	    c_end = std::partition(c_beg, c_end, std::mem_fun_ref(&ACluster::get_valid));
	    std::random_shuffle(c_beg, c_end);
	    
	    for (ACluster *i=c_beg; i<c_end; i++) // Reset for current merging round
		i->reset_RM();	

	    size_t fold = std::max((size_t)1 /* 2x */, (size_t)(c_end - c_beg) / 5000);
	    while (true) {
		// Find first valid, un-merged cluster in this round
		ACluster* into = std::find_if(c_beg, c_end, std::not1(std::mem_fun_ref(&ACluster::get_merged)));
		if (into == c_end)
		    goto END_ROUND;		

		into->set_merged(true);

		// Find the "fold" nearest other clusters for row major data
		ACluster::SLinkLT_RM cmp(*into);
		std::partial_sort(into+1, std::min(into+1+fold,c_end), c_end, cmp);

		for (size_t i=1; i<=fold && (into+i)<c_end; i++) {
		    if (!into[i].get_valid() || into[i].get_merged())
			goto END_ROUND;

		    into->merge_in(into[i]);

		    // Record merging operation
		    merge[cur_merge].set(into->get_cluster_id(), into[i].get_cluster_id());
		    into->set_cluster_id(++cur_merge); // Recall merging is 1-indexed
		}
	    }
END_ROUND:;
	}

	// Clean up
	for (ACluster *i=c_beg; i<c_end; i++) {
	    if (i->get_valid()) {
		i->set_valid(false);
		i->destroy();
	    }
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

	cluster(static_cast<Data_t*>(REAL(tbl)), obs, dim, (Merge_t*)merge_l);
		

	// Assigning observations to clusters	
	PROTECT(assgn = allocVector(INTSXP, obs)); n_protected++; 	

	std::set<Idx_t> clusters;
	live_clusters((Idx_t*)merge_l, obs, k_l, clusters);  // Cut tree to find ~k clusters
	
	Idx_t cluster_index = 1; 
	for (std::set<Idx_t>::iterator i=clusters.begin(),e=clusters.end(); i!=e; ++i)
	    assign_observation((Idx_t*)merge_l, *i, cluster_index++, static_cast<Idx_t*>(INTEGER(assgn)));	   


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

