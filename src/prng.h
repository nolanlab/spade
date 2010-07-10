#ifndef PRNG_H
#define PRNG_H

#include <stdint.h>
#include "util.h"

// Produce a unique stream of pseudo random numbers using Galois LFSR 
class PRNG {
    static const uint32_t polynomial[];	
    uint32_t lfsr, mask;

public:
    PRNG(size_t min_period) : lfsr(1) {
	uint32_t bits = (uint32_t)ceil(log2(min_period+1)); 
	assert(bits >= 2 && bits < 20);
	mask = polynomial[bits];
    }

    void seed(uint32_t s) { assert(s); lfsr = s; }

    uint32_t sample() {
        lfsr = (lfsr >> 1) ^ (uint32_t)(0 - (lfsr & 1u) & mask);
        return lfsr;
    }
    uint32_t sample(uint32_t n) {
        uint32_t s;
        do { s = sample(); } while (s >= n); // Filter out entries above maximum value "n"
        return s;
    }
};

#endif
