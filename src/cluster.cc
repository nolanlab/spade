#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Print.h>
#include <cmath>
#include <stdint.h>
#include <limits>
#include <algorithm>
#include <vector>
#include <utility>
#include <functional>
#include <set>
#include <queue>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef HAVE_PARALLEL_STL
#include <parallel/algorithm>
#define PARALLEL_NAMESPACE std::__parallel
#else
#define PARALLEL_NAMESPACE std
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

	class ACluster;

	class PQ {
	public:
		typedef std::pair<ACluster*,Dist_t> Value_t;
		struct CMP {
			bool operator()(const Value_t& l, const Value_t& r) { return l.second < r.second; }
		};
		std::priority_queue<Value_t,std::vector<Value_t>,CMP> Q;
		
		Dist_t MMD;  // Minimum merged distance (distance to nearest already merged neighbor)
		size_t Fold;

		PQ(size_t fold) : MMD(std::numeric_limits<Dist_t>::max()), Fold(fold) {}
	
		bool empty() const { return Q.empty(); }

		void normalize();	

		void push(ACluster* c, Dist_t d);    
		void push(PQ& pq);

		ACluster* furthestCluster() const { return Q.top().first; }
		void popFurthest() { Q.pop(); }
	};

    /* Maintain state for a single cluster
     */
	class ACluster {
	public:
		static size_t dim;
		static void init_global(size_t d) { dim = d; }

	public:
		Data_t* center;		// My cluster center
		bool valid, merged;  // Cluster statuses
		const Data_t **members;	// Observations in my cluster
		size_t num_members;


		// Satisfied with compiler generated constructors

		// Initialization and destruction
		void init_RM(const Data_t* data) {
			center = Calloc(dim, Data_t);
			memcpy(center, data, dim*sizeof(Data_t));

			valid = true; merged = false;

			members = Calloc(1, const Data_t*);
			members[0]  = data;
			num_members = 1;
		}
		void destroy() {
			center  = (Free(center), (Data_t*)0);
			members = (Free(members), (const Data_t**)0);
			num_members = 0;
		}

		// Getters & Setters (many used in STL algorithms)
		bool get_valid() const { return valid; }
		void set_valid(bool v) { valid = v; }

		bool get_merged() const { return merged; }
		void set_merged(bool m) { merged = m; }

		// Member observations
		size_t get_num_members() const { return num_members; }

		typedef const Data_t** member_iterator;
		member_iterator member_begin() const { return members; }
		member_iterator member_end() const { return members + num_members; }


		// Helper structs

		// LessThan by data column (used in Median)
		struct ColLT_RM {
			size_t col;
			ColLT_RM(size_t c) : col(c) {}
			bool operator()(const Data_t* l, const Data_t* r) const { return l[col] < r[col]; }
		};

		struct NMLT {
			bool operator()(const ACluster& l, const ACluster& r) const { return l.num_members < r.num_members; }
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

		void push_on_pq(ACluster* from, PQ& pq) const {
			for (member_iterator i=from->member_begin(), e=from->member_end(); i!=e; ++i)
				pq.push(from, distance(center, *i, dim)); 
		}

		void merge_in_pq(PQ& pq) {
			// Discard all potential merges further than nearest neighbor that
			// has already been merged in this round
			pq.normalize();

			while (!pq.empty()) {
				ACluster& rhs = *pq.furthestCluster();
				if (!rhs.get_merged())
					this->merge_in(rhs);
				pq.popFurthest();
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

	void PQ::normalize() {
		while (!Q.empty()) {
			if (Q.top().second >= MMD)
				Q.pop();
			else
				break;
		}
	}

	void PQ::push(ACluster* c, Dist_t d) {
		if (c->get_merged())
			MMD = std::min(MMD, d);
		else if (d < MMD) {
			Q.push(std::make_pair(c,d));
			if (Q.size() > Fold)
				Q.pop();
		}
	}

	void PQ::push(PQ& pq) {
		MMD = std::min(MMD, pq.MMD);	
		while (!pq.Q.empty()) {
			this->push(pq.Q.top().first, pq.Q.top().second);
			pq.Q.pop();
		}
	}
	
	void
    cluster(const Data_t* data, size_t obs, size_t dim, size_t k, Idx_t* assgn) {
		ACluster::init_global(dim);

		R::auto_ptr<ACluster> c_ap(Calloc(obs, ACluster));
		ACluster *c_beg = c_ap.get(), *c_end = c_beg + obs;
		for (size_t i=0; i<obs; i++)  // Initialize clusters for row major data
			c_beg[i].init_RM(&data[i*dim]);

		double max_rounds = log2((double)(obs)-1);
		for (size_t round = 0; ; round++) {
			Rprintf("  Estimated clustering progress: %2.0f%% ...\n",std::min(99.0,((double)round)/max_rounds * 100.0));

			// Only looking at "valid" clusters
			c_end = std::partition(c_beg, c_end, std::mem_fun_ref(&ACluster::get_valid));
						
			// Hueristic to strip out singleton clusters
			if (round == 5) {
				for (ACluster *i=c_beg; i<c_end; i++)
					if (i->get_num_members() == 1)
						i->set_valid(false);
				c_end = std::partition(c_beg, c_end, std::mem_fun_ref(&ACluster::get_valid));
			}

			// Stopping condition
			size_t num_valid = c_end - c_beg;
			if (num_valid < (size_t)(1.5 * k))
				break;

			if (round) {  // Clusters already initialized in "0" round
				ACluster::NMLT cmp; // Order by cluster size to preferentially merge "small" clusters
				std::sort(c_beg, c_end, cmp);
				for (ACluster *i=c_beg; i<c_end; i++)  // Reset clusters for current merging round
					i->reset_RM();	
			}

			size_t fold = std::max((size_t)1 /* 2x */, num_valid / 5000);
			for (ACluster* into=c_beg; into < c_end; ++into) {
				if (into->get_merged())
					continue;
				into->set_merged(true);

				PQ pq(fold);
#ifdef _OPENMP
				#pragma omp parallel shared(c_beg, c_end, into)
#endif
				{
					// Scan all "valid" clusters to fill in queue of nearest neighbors 
					PQ pq_p(fold);  // "private" priority queue

#ifdef _OPENMP
					#pragma omp for nowait
#endif
					for (size_t i=0; i < (c_end - c_beg); i++) {
						ACluster* from = c_beg + i;
						if (into != from)
							into->push_on_pq(from, pq_p);
					}

#ifdef _OPENMP
					#pragma omp critical
#endif
					{
						pq.push(pq_p);
					}
				}

				// Merge in clusters in queue
				into->merge_in_pq(pq);
			}
		}

		// Assignment and clean up
		Idx_t cur_Id = 1;  // Recall R is 1-indexed
		c_end = std::partition(c_beg, c_end, std::mem_fun_ref(&ACluster::get_valid));
		for (ACluster *i=c_beg; i<c_end; i++, cur_Id++) {
			for (ACluster::member_iterator b=i->member_begin(), e=i->member_end(); b!=e; ++b)
				assgn[(*b - data) / dim] = cur_Id;
			i->destroy();
		}	    
	
		Rprintf("  Estimated clustering progress: %3.0f%%\n",100.0);	
	}
		

   
} // anonymous namespace


extern "C" {

    SEXP 
    SPADE_cluster(SEXP tbl, SEXP k) {
		// tbl is column major order, transposed to be row major
		size_t obs = static_cast<size_t>(INTEGER(GET_DIM(tbl))[1]), 
			   dim = static_cast<size_t>(INTEGER(GET_DIM(tbl))[0]);
		size_t k_l = static_cast<size_t>(asInteger(k));

		SEXP ans, names, assgn; int n_protected = 0;

		// Cluster and assign observations to clusters	
		PROTECT(assgn = allocVector(INTSXP, obs)); n_protected++;
		Idx_t* assgn_l = static_cast<Idx_t*>(INTEGER(assgn));
		std::fill(assgn_l, assgn_l+obs, 0);
		cluster(static_cast<Data_t*>(REAL(tbl)), obs, dim, k_l, assgn_l);

		// Assemble return object
		PROTECT(ans = allocVector(VECSXP,1)); n_protected++;
		PROTECT(names = allocVector(STRSXP,1)); n_protected++;

		SET_VECTOR_ELT(ans, 0, assgn);
		SET_STRING_ELT(names, 0, mkChar("assgn"));

		setAttrib(ans, R_NamesSymbol, names);
		UNPROTECT(n_protected);


		return ans;
    }

}

