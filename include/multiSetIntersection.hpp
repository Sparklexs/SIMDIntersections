/*
 * multiSetIntersectionSIMD.hpp
 *
 *  Created on: 2017Äê2ÔÂ5ÈÕ
 *      Author: John
 */

#ifndef INCLUDE_MULTISETINTERSECTION_HPP_
#define INCLUDE_MULTISETINTERSECTION_HPP_

#include "common.h"
#include "thomaswu.h"
#include "intersectionfactory.h"

namespace msis/*MultiSet InterSection*/{
typedef size_t (*setIntersectionFunction)(const uint32_t * set1,
		const size_t length1, const uint32_t * set2, const size_t length2,
		uint32_t *out);

// here adapt the range [start(0),end], different from __BSadvanceUntil
// whose range is [start+1,end-1]
long scalarBinarySearch(uint32_t min, const uint32_t * array, long end) {
	if (0 == end || *array >= min) {
		return 0;
	}

	size_t lower = 0;
	size_t upper = end;
	size_t mid;
	while (lower < upper) {
		mid = (lower + upper) / 2;
		if (array[mid] == min) {
			return mid;
		}

		if (array[mid] < min) {
			lower = mid + 1;
		} else {
			upper = mid;
		}
	}
	return upper;
}

// here adapt the range [start,end]
long scalargallop(uint32_t min, const uint32_t * array, long end) {
	// special handling for a possibly common sequential case
	if ((0 >= end) or (*array >= min)) {
		return 0;
	}

	size_t upper = 1; // could set larger
	// bootstrap an upper limit

	// sxs: here spansize is enlarged to the maximum where its corresponding
	// element is geq min
	while ((upper <= end) and (array[upper] < min))
		upper <<= 1;

	// we know that the next-smallest span was too small
	size_t lower = (upper / 2);

	if (upper > end)
		upper = end;

	if (array[upper] <= min) {    // means array has no item >= min
		return upper;
	}

	// else begin binary search
	size_t mid = 0;
	while (lower + 1 != upper) {
		mid = (lower + upper) / 2;
		if (array[mid] == min) {
			return mid;
		} else if (array[mid] < min)
			lower = mid;
		else
			upper = mid;
	}
	return upper;
}

/* return the exact position of intersection */
long simdgallop_v3_exact(UINT4 goal, const UINT4 *target, long ntargets) {
	const UINT4 *end_target, *stop_target, *init_target;
	const UINT4 V3_BLOCKSIZE = 4 * 32;
	const UINT4 SIMDWIDTH = 4;

	long low_offset, mid_offset, high_offset;

	init_target = target;
	stop_target = target + ntargets - V3_BLOCKSIZE;
	end_target = target + ntargets;

	if (target >= stop_target)
		return Intersection_find_scalar(goal, target, ntargets);

	if (target[V3_BLOCKSIZE - 1] < goal) {
		if (target + V3_BLOCKSIZE > stop_target)
			return V3_BLOCKSIZE
					+ Intersection_find_scalar(goal, target + V3_BLOCKSIZE,
							ntargets - V3_BLOCKSIZE);

		/* Galloping search */
		high_offset = 1;
		while (target[V3_BLOCKSIZE * high_offset + V3_BLOCKSIZE - 1] < goal) {
			if (target + (high_offset << 1) * V3_BLOCKSIZE <= stop_target) {
				low_offset = high_offset;
				high_offset <<= 1;
			} else if (target + V3_BLOCKSIZE * (high_offset + 1)
					<= stop_target) {
				high_offset = (stop_target - target) / V3_BLOCKSIZE;
				if (target[V3_BLOCKSIZE * high_offset + V3_BLOCKSIZE - 1]
						< goal) {
					target += V3_BLOCKSIZE * high_offset;
					return V3_BLOCKSIZE * high_offset
							+ Intersection_find_scalar(goal, target,
									(end_target - target));
				} else
					break;
			} else {
				target += V3_BLOCKSIZE * high_offset;
				return V3_BLOCKSIZE * high_offset
						+ Intersection_find_scalar(goal, target,
								(end_target - target));
			}
		}			// while-loop
		while (low_offset < high_offset) {
			mid_offset = (low_offset + high_offset) / 2;
			if (target[V3_BLOCKSIZE * mid_offset + V3_BLOCKSIZE - 1] < goal) {
				low_offset = mid_offset + 1;
			} else {
				high_offset = mid_offset;
			}
		}
		target += V3_BLOCKSIZE * high_offset;
	}
	// now search the block [target, target+block_size-1]
	__m128i Match;
	__m128i F0, Q0, Q1, Q2, Q3;
	UINT4 *hits;
	long pos, base;
	Match = _mm_set1_epi32(goal);

	if (target[SIMDWIDTH * 16 - 1] >= goal) {
		if (target[SIMDWIDTH * 8 - 1] >= goal) {
			base = 0;
			Q0 = _mm_or_si128(
					_mm_srli_epi32(
							_mm_cmpeq_epi32(
									_mm_lddqu_si128((__m128i *) target + 0),
									Match), 32 - 1),
					_mm_srli_epi32(
							_mm_cmpeq_epi32(
									_mm_lddqu_si128((__m128i *) target + 1),
									Match), 28 - 1));
			Q1 = _mm_or_si128(
					_mm_srli_epi32(
							_mm_cmpeq_epi32(
									_mm_lddqu_si128((__m128i *) target + 2),
									Match), 24 - 1),
					_mm_srli_epi32(
							_mm_cmpeq_epi32(
									_mm_lddqu_si128((__m128i *) target + 3),
									Match), 20 - 1));
			Q2 = _mm_or_si128(
					_mm_srli_epi32(
							_mm_cmpeq_epi32(
									_mm_lddqu_si128((__m128i *) target + 4),
									Match), 16 - 1),
					_mm_srli_epi32(
							_mm_cmpeq_epi32(
									_mm_lddqu_si128((__m128i *) target + 5),
									Match), 12 - 1));
			Q3 = _mm_or_si128(
					_mm_srli_epi32(
							_mm_cmpeq_epi32(
									_mm_lddqu_si128((__m128i *) target + 6),
									Match), 8 - 1),
					_mm_srli_epi32(
							_mm_cmpeq_epi32(
									_mm_lddqu_si128((__m128i *) target + 7),
									Match), 4 - 1));
		} else {
			base = 32;
			Q0 = _mm_or_si128(
					_mm_srli_epi32(
							_mm_cmpeq_epi32(
									_mm_lddqu_si128((__m128i *) target + 8),
									Match), 32 - 1),
					_mm_srli_epi32(
							_mm_cmpeq_epi32(
									_mm_lddqu_si128((__m128i *) target + 9),
									Match), 28 - 1));
			Q1 = _mm_or_si128(
					_mm_srli_epi32(
							_mm_cmpeq_epi32(
									_mm_lddqu_si128((__m128i *) target + 10),
									Match), 24 - 1),
					_mm_srli_epi32(
							_mm_cmpeq_epi32(
									_mm_lddqu_si128((__m128i *) target + 11),
									Match), 20 - 1));
			Q2 = _mm_or_si128(
					_mm_srli_epi32(
							_mm_cmpeq_epi32(
									_mm_lddqu_si128((__m128i *) target + 12),
									Match), 16 - 1),
					_mm_srli_epi32(
							_mm_cmpeq_epi32(
									_mm_lddqu_si128((__m128i *) target + 13),
									Match), 12 - 1));
			Q3 = _mm_or_si128(
					_mm_srli_epi32(
							_mm_cmpeq_epi32(
									_mm_lddqu_si128((__m128i *) target + 14),
									Match), 8 - 1),
					_mm_srli_epi32(
							_mm_cmpeq_epi32(
									_mm_lddqu_si128((__m128i *) target + 15),
									Match), 4 - 1));
		}
	} else {
		if (target[SIMDWIDTH * 24 - 1] >= goal) {
			base = 64;
			Q0 = _mm_or_si128(
					_mm_srli_epi32(
							_mm_cmpeq_epi32(
									_mm_lddqu_si128(
											(__m128i *) target + 0 + 16),
									Match), 32 - 1),
					_mm_srli_epi32(
							_mm_cmpeq_epi32(
									_mm_lddqu_si128(
											(__m128i *) target + 1 + 16),
									Match), 28 - 1));
			Q1 = _mm_or_si128(
					_mm_srli_epi32(
							_mm_cmpeq_epi32(
									_mm_lddqu_si128(
											(__m128i *) target + 2 + 16),
									Match), 24 - 1),
					_mm_srli_epi32(
							_mm_cmpeq_epi32(
									_mm_lddqu_si128(
											(__m128i *) target + 3 + 16),
									Match), 20 - 1));
			Q2 = _mm_or_si128(
					_mm_srli_epi32(
							_mm_cmpeq_epi32(
									_mm_lddqu_si128(
											(__m128i *) target + 4 + 16),
									Match), 16 - 1),
					_mm_srli_epi32(
							_mm_cmpeq_epi32(
									_mm_lddqu_si128(
											(__m128i *) target + 5 + 16),
									Match), 12 - 1));
			Q3 = _mm_or_si128(
					_mm_srli_epi32(
							_mm_cmpeq_epi32(
									_mm_lddqu_si128(
											(__m128i *) target + 6 + 16),
									Match), 8 - 1),
					_mm_srli_epi32(
							_mm_cmpeq_epi32(
									_mm_lddqu_si128(
											(__m128i *) target + 7 + 16),
									Match), 4 - 1));
		} else {
			base = 96;
			Q0 = _mm_or_si128(
					_mm_srli_epi32(
							_mm_cmpeq_epi32(
									_mm_lddqu_si128(
											(__m128i *) target + 8 + 16),
									Match), 32 - 1),
					_mm_srli_epi32(
							_mm_cmpeq_epi32(
									_mm_lddqu_si128(
											(__m128i *) target + 9 + 16),
									Match), 28 - 1));
			Q1 = _mm_or_si128(
					_mm_srli_epi32(
							_mm_cmpeq_epi32(
									_mm_lddqu_si128(
											(__m128i *) target + 10 + 16),
									Match), 24 - 1),
					_mm_srli_epi32(
							_mm_cmpeq_epi32(
									_mm_lddqu_si128(
											(__m128i *) target + 11 + 16),
									Match), 20 - 1));
			Q2 = _mm_or_si128(
					_mm_srli_epi32(
							_mm_cmpeq_epi32(
									_mm_lddqu_si128(
											(__m128i *) target + 12 + 16),
									Match), 16 - 1),
					_mm_srli_epi32(
							_mm_cmpeq_epi32(
									_mm_lddqu_si128(
											(__m128i *) target + 13 + 16),
									Match), 12 - 1));
			Q3 = _mm_or_si128(
					_mm_srli_epi32(
							_mm_cmpeq_epi32(
									_mm_lddqu_si128(
											(__m128i *) target + 14 + 16),
									Match), 8 - 1),
					_mm_srli_epi32(
							_mm_cmpeq_epi32(
									_mm_lddqu_si128(
											(__m128i *) target + 15 + 16),
									Match), 4 - 1));
		}
	}

	F0 = _mm_or_si128(_mm_or_si128(Q0, Q1), _mm_or_si128(Q2, Q3));
	if (_mm_testz_si128(F0, F0)) {
		return (target - init_target) + base;
	} else {
		hits = (UINT4*) &F0;
#ifdef __AVX2__
		F0 = _mm_sllv_epi32(F0,_mm_set_epi32(3,2,1,0));
		pos = __builtin_clz(hits[0] | hits[1] | hits[2] | hits[3]);
#else
		pos = __builtin_clz(
				hits[0] | (hits[1] << 1) | (hits[2] << 2) | (hits[3] << 3));
#endif
		return (target - init_target) + base + (31 - pos);
	}
}

/* return the head-position of mathced block */
long simdgallop_v3_rough(int *foundp, UINT4 goal, const UINT4 *target,
		long ntargets) {
	const UINT4 *end_target, *stop_target, *init_target;
	const UINT4 SIMDWIDTH = 4;
	const UINT4 V3_BLOCKSIZE = SIMDWIDTH * 32;

	long low_offset, mid_offset, high_offset;
	long pos, base;

	init_target = target;
	stop_target = target + ntargets - V3_BLOCKSIZE;
	end_target = target + ntargets;

	if (target >= stop_target) {
		if ((pos = Intersection_find_scalar(goal, target, ntargets)) <= ntargets
				&& target[pos] == goal)
			*foundp = 1;
		else
			*foundp = 0;
		return pos;
	}

	if (target[V3_BLOCKSIZE - 1] < goal) {
		if (target + V3_BLOCKSIZE > stop_target) {
			pos = V3_BLOCKSIZE
					+ Intersection_find_scalar(goal, target + V3_BLOCKSIZE,
							ntargets - V3_BLOCKSIZE);
			if (pos <= ntargets && target[pos] == goal)
				*foundp = 1;
			else
				*foundp = 0;
			return pos;
		}
		/* Galloping search */
		high_offset = 1;
		while (target[V3_BLOCKSIZE * high_offset + V3_BLOCKSIZE - 1] < goal) {
			if (target + (high_offset << 1) * V3_BLOCKSIZE <= stop_target) {
				low_offset = high_offset;
				high_offset <<= 1;
			} else if (target + V3_BLOCKSIZE * (high_offset + 1)
					<= stop_target) {
				high_offset = (stop_target - target) / V3_BLOCKSIZE;
				if (target[V3_BLOCKSIZE * high_offset + V3_BLOCKSIZE - 1]
						< goal) {
					target += V3_BLOCKSIZE * high_offset;
					pos = Intersection_find_scalar(goal, target,
							(end_target - target));
					if (pos + V3_BLOCKSIZE * high_offset <= ntargets
							&& target[pos] == goal)
						*foundp = 1;
					else
						*foundp = 0;
					return pos + V3_BLOCKSIZE * high_offset;
				} else
					break;
			} else {
				target += V3_BLOCKSIZE * high_offset;

				pos = Intersection_find_scalar(goal, target,
						(end_target - target));
				if (pos + V3_BLOCKSIZE * high_offset <= ntargets
						&& target[pos] == goal)
					*foundp = 1;
				else
					*foundp = 0;
				return pos + V3_BLOCKSIZE * high_offset;
			}
		}			// while-loop
		while (low_offset < high_offset) {
			mid_offset = (low_offset + high_offset) / 2;
			if (target[V3_BLOCKSIZE * mid_offset + V3_BLOCKSIZE - 1] < goal) {
				low_offset = mid_offset + 1;
			} else {
				high_offset = mid_offset;
			}
		}
		target += V3_BLOCKSIZE * high_offset;
	}
	// now search the block [target, target+block_size-1]
	__m128i Match;
	__m128i F0, Q0, Q1, Q2, Q3;
	Match = _mm_set1_epi32(goal);

	if (target[SIMDWIDTH * 16 - 1] >= goal) {
		if (target[SIMDWIDTH * 8 - 1] >= goal) {
			base = 0;
			Q0 = _mm_or_si128(
					_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 0),
							Match),
					_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 1),
							Match));
			Q1 = _mm_or_si128(
					_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 2),
							Match),
					_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 3),
							Match));
			Q2 = _mm_or_si128(
					_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 4),
							Match),
					_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 5),
							Match));
			Q3 = _mm_or_si128(
					_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 6),
							Match),
					_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 7),
							Match));
		} else {
			base = 32;
			Q0 = _mm_or_si128(
					_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 8),
							Match),
					_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 9),
							Match));
			Q1 = _mm_or_si128(
					_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 10),
							Match),
					_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 11),
							Match));
			Q2 = _mm_or_si128(
					_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 12),
							Match),
					_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 13),
							Match));
			Q3 = _mm_or_si128(
					_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 14),
							Match),
					_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 15),
							Match));
		}
	} else {
		if (target[SIMDWIDTH * 24 - 1] >= goal) {
			base = 64;
			Q0 = _mm_or_si128(
					_mm_cmpeq_epi32(
							_mm_lddqu_si128((__m128i *) target + 0 + 16),
							Match),
					_mm_cmpeq_epi32(
							_mm_lddqu_si128((__m128i *) target + 1 + 16),
							Match));
			Q1 = _mm_or_si128(
					_mm_cmpeq_epi32(
							_mm_lddqu_si128((__m128i *) target + 2 + 16),
							Match),
					_mm_cmpeq_epi32(
							_mm_lddqu_si128((__m128i *) target + 3 + 16),
							Match));
			Q2 = _mm_or_si128(
					_mm_cmpeq_epi32(
							_mm_lddqu_si128((__m128i *) target + 4 + 16),
							Match),
					_mm_cmpeq_epi32(
							_mm_lddqu_si128((__m128i *) target + 5 + 16),
							Match));
			Q3 = _mm_or_si128(
					_mm_cmpeq_epi32(
							_mm_lddqu_si128((__m128i *) target + 6 + 16),
							Match),
					_mm_cmpeq_epi32(
							_mm_lddqu_si128((__m128i *) target + 7 + 16),
							Match));
		} else {
			base = 96;
			Q0 = _mm_or_si128(
					_mm_cmpeq_epi32(
							_mm_lddqu_si128((__m128i *) target + 8 + 16),
							Match),
					_mm_cmpeq_epi32(
							_mm_lddqu_si128((__m128i *) target + 9 + 16),
							Match));
			Q1 = _mm_or_si128(
					_mm_cmpeq_epi32(
							_mm_lddqu_si128((__m128i *) target + 10 + 16),
							Match),
					_mm_cmpeq_epi32(
							_mm_lddqu_si128((__m128i *) target + 11 + 16),
							Match));
			Q2 = _mm_or_si128(
					_mm_cmpeq_epi32(
							_mm_lddqu_si128((__m128i *) target + 12 + 16),
							Match),
					_mm_cmpeq_epi32(
							_mm_lddqu_si128((__m128i *) target + 13 + 16),
							Match));
			Q3 = _mm_or_si128(
					_mm_cmpeq_epi32(
							_mm_lddqu_si128((__m128i *) target + 14 + 16),
							Match),
					_mm_cmpeq_epi32(
							_mm_lddqu_si128((__m128i *) target + 15 + 16),
							Match));
		}
	}

	F0 = _mm_or_si128(_mm_or_si128(Q0, Q1), _mm_or_si128(Q2, Q3));
	if (_mm_testz_si128(F0, F0))
		*foundp = 0;
	else
		*foundp = 1;
	return (target - init_target) + base;
}

long simdgallop_v3_greedy(UINT4 goal, const UINT4 *target, long ntargets) {
	const UINT4 *end_target, *stop_target, *init_target;
	const UINT4 SIMDWIDTH = 4;
	const UINT4 V3_BLOCKSIZE = SIMDWIDTH * 32;

	long low_offset, mid_offset, high_offset;
	long pos, base;

	init_target = target;
	stop_target = target + ntargets - V3_BLOCKSIZE;
	end_target = target + ntargets;

	if (target >= stop_target)
		return Intersection_find_scalar(goal, target, ntargets);

	if (target[V3_BLOCKSIZE - 1] < goal) {
		if (target + V3_BLOCKSIZE > stop_target)
			return V3_BLOCKSIZE
					+ Intersection_find_scalar(goal, target + V3_BLOCKSIZE,
							ntargets - V3_BLOCKSIZE);

		/* Galloping search */
		high_offset = 1;
		while (target[V3_BLOCKSIZE * high_offset + V3_BLOCKSIZE - 1] < goal) {
			if (target + (high_offset << 1) * V3_BLOCKSIZE <= stop_target) {
				low_offset = high_offset;
				high_offset <<= 1;
			} else if (target + V3_BLOCKSIZE * (high_offset + 1)
					<= stop_target) {
				high_offset = (stop_target - target) / V3_BLOCKSIZE;
				if (target[V3_BLOCKSIZE * high_offset + V3_BLOCKSIZE - 1]
						< goal) {
					target += V3_BLOCKSIZE * high_offset;
					return V3_BLOCKSIZE * high_offset
							+ Intersection_find_scalar(goal, target,
									(end_target - target));
				} else
					break;
			} else {
				target += V3_BLOCKSIZE * high_offset;
				return V3_BLOCKSIZE * high_offset
						+ Intersection_find_scalar(goal, target,
								(end_target - target));
			}
		}			// while-loop
		while (low_offset < high_offset) {
			mid_offset = (low_offset + high_offset) / 2;
			if (target[V3_BLOCKSIZE * mid_offset + V3_BLOCKSIZE - 1] < goal) {
				low_offset = mid_offset + 1;
			} else {
				high_offset = mid_offset;
			}
		}
		target += V3_BLOCKSIZE * high_offset;
	}
	// now search the block [target, target+block_size-1]
	__m128i Match;
	__m128i F0, Q0, Q1, Q2, Q3;
	Match = _mm_set1_epi32(goal);

	if (target[SIMDWIDTH * 16 - 1] >= goal) {
		if (target[SIMDWIDTH * 8 - 1] >= goal) {
			base = 0;
			Q0 = _mm_or_si128(
					_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 0),
							Match),
					_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 1),
							Match));
			Q1 = _mm_or_si128(
					_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 2),
							Match),
					_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 3),
							Match));
			Q2 = _mm_or_si128(
					_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 4),
							Match),
					_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 5),
							Match));
			Q3 = _mm_or_si128(
					_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 6),
							Match),
					_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 7),
							Match));
		} else {
			base = 32;
			Q0 = _mm_or_si128(
					_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 8),
							Match),
					_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 9),
							Match));
			Q1 = _mm_or_si128(
					_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 10),
							Match),
					_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 11),
							Match));
			Q2 = _mm_or_si128(
					_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 12),
							Match),
					_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 13),
							Match));
			Q3 = _mm_or_si128(
					_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 14),
							Match),
					_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 15),
							Match));
		}
	} else {
		if (target[SIMDWIDTH * 24 - 1] >= goal) {
			base = 64;
			Q0 = _mm_or_si128(
					_mm_cmpeq_epi32(
							_mm_lddqu_si128((__m128i *) target + 0 + 16),
							Match),
					_mm_cmpeq_epi32(
							_mm_lddqu_si128((__m128i *) target + 1 + 16),
							Match));
			Q1 = _mm_or_si128(
					_mm_cmpeq_epi32(
							_mm_lddqu_si128((__m128i *) target + 2 + 16),
							Match),
					_mm_cmpeq_epi32(
							_mm_lddqu_si128((__m128i *) target + 3 + 16),
							Match));
			Q2 = _mm_or_si128(
					_mm_cmpeq_epi32(
							_mm_lddqu_si128((__m128i *) target + 4 + 16),
							Match),
					_mm_cmpeq_epi32(
							_mm_lddqu_si128((__m128i *) target + 5 + 16),
							Match));
			Q3 = _mm_or_si128(
					_mm_cmpeq_epi32(
							_mm_lddqu_si128((__m128i *) target + 6 + 16),
							Match),
					_mm_cmpeq_epi32(
							_mm_lddqu_si128((__m128i *) target + 7 + 16),
							Match));
		} else {
			base = 96;
			Q0 = _mm_or_si128(
					_mm_cmpeq_epi32(
							_mm_lddqu_si128((__m128i *) target + 8 + 16),
							Match),
					_mm_cmpeq_epi32(
							_mm_lddqu_si128((__m128i *) target + 9 + 16),
							Match));
			Q1 = _mm_or_si128(
					_mm_cmpeq_epi32(
							_mm_lddqu_si128((__m128i *) target + 10 + 16),
							Match),
					_mm_cmpeq_epi32(
							_mm_lddqu_si128((__m128i *) target + 11 + 16),
							Match));
			Q2 = _mm_or_si128(
					_mm_cmpeq_epi32(
							_mm_lddqu_si128((__m128i *) target + 12 + 16),
							Match),
					_mm_cmpeq_epi32(
							_mm_lddqu_si128((__m128i *) target + 13 + 16),
							Match));
			Q3 = _mm_or_si128(
					_mm_cmpeq_epi32(
							_mm_lddqu_si128((__m128i *) target + 14 + 16),
							Match),
					_mm_cmpeq_epi32(
							_mm_lddqu_si128((__m128i *) target + 15 + 16),
							Match));
		}
	}

	F0 = _mm_or_si128(_mm_or_si128(Q0, Q1), _mm_or_si128(Q2, Q3));
	if (_mm_testz_si128(F0, F0))
		return (target - init_target) + base;
	else {
		target += base;
		Q0 = _mm_or_si128(
				_mm_srli_epi32(
						_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 0),
								Match), 32 - 1),
				_mm_srli_epi32(
						_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 1),
								Match), 28 - 1));
		Q1 = _mm_or_si128(
				_mm_srli_epi32(
						_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 2),
								Match), 24 - 1),
				_mm_srli_epi32(
						_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 3),
								Match), 20 - 1));
		Q2 = _mm_or_si128(
				_mm_srli_epi32(
						_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 4),
								Match), 16 - 1),
				_mm_srli_epi32(
						_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 5),
								Match), 12 - 1));
		Q3 = _mm_or_si128(
				_mm_srli_epi32(
						_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 6),
								Match), 8 - 1),
				_mm_srli_epi32(
						_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 7),
								Match), 4 - 1));
		F0 = _mm_or_si128(_mm_or_si128(Q0, Q1), _mm_or_si128(Q2, Q3));
		UINT4 *hits = (UINT4*) &F0;
#ifdef __AVX2__
		F0 = _mm_sllv_epi32(F0,_mm_set_epi32(3,2,1,0));
		pos = __builtin_clz(hits[0] | hits[1] | hits[2] | hits[3]);
#else
		pos = __builtin_clz(
				hits[0] | (hits[1] << 1) | (hits[2] << 2) | (hits[3] << 3));
#endif
		return (target - init_target) + (31 - pos);
	}
}

////////////////////////////////////////////////////////////////////
template<setIntersectionFunction FUNCTION>
void small_vs_small(const mySet &sets, std::vector<uint32_t> &out) {
	mySet::iterator it = sets.begin();
	// XXX: we'd like to use rvalue reference, however, it has conflicts
	// with "const", and it finally become an copy operation.
	std::vector<uint32_t> intersection(std::move(*it++));
	for (; it != sets.end(); it++) {
		// here we can change the intersection function to any regular scalar
		// or vector pair-set intersection algorithms.
		size_t inter_length = FUNCTION(intersection.data(), intersection.size(),
				it->data(), it->size(), intersection.data());
		intersection.resize(inter_length);
	}
	out.swap(intersection);
}

template<intersectionfindfunction FINDFUNCTION>
void set_vs_set(const mySet &sets, std::vector<uint32_t> &out) {
	size_t count = 0, currentset = 0;
	out.resize(sets.begin()->size());

	mySet::iterator it = sets.begin();
	std::vector<uint32_t> candidate(std::move(*it++));
	const mySet::iterator it_start = it;

	auto value = out.begin();
	auto eliminator = candidate.begin();
	std::vector<size_t> index(sets.size() - 1, 0);

	while (eliminator != candidate.end()) {
		index[currentset] += FINDFUNCTION(*eliminator,
				it->data() + index[currentset],
				it->size() - 1 - index[currentset]);

		if (it->at(index[currentset]) == *eliminator) {
			if (++currentset == sets.size() - 1) {
				*value++ = *eliminator;
				count++;
				// roll over
				currentset = 0;
				it = it_start;
				eliminator++;
			} else {
				++it;
			}
		} else {
			currentset = 0;
			it = it_start;
			eliminator++;
		}
	}
	out.resize(count);
}

// only used for find-functions that return rough position
template<flaggedintersectionfindfunction FINDFUNCTION>
void set_vs_set_flagged(const mySet &sets, std::vector<uint32_t> &out) {
	size_t count = 0, currentset = 0;
	out.resize(sets.begin()->size());

	mySet::iterator it = sets.begin();
	std::vector<uint32_t> candidate(std::move(*it++));
	const mySet::iterator it_start = it;

	int foundp;
	auto value = out.begin();
	auto eliminator = candidate.begin();
	std::vector<size_t> index(sets.size() - 1, 0);

	while (eliminator != candidate.end()) {
		index[currentset] += FINDFUNCTION(&foundp, *eliminator,
				it->data() + index[currentset],
				it->size() - 1 - index[currentset]);

		if (foundp == 1) {
			if (++currentset == sets.size() - 1) {
				*value++ = *eliminator;
				count++;
				// roll over
				currentset = 0;
				it = it_start;
				eliminator++;
			} else {
				++it;
			}
		} else {
			currentset = 0;
			it = it_start;
			eliminator++;
		}
	}
	out.resize(count);
}

template<intersectionfindfunction FINDFUNCTION>
void swapping_set_vs_set(const mySet& sets, std::vector<uint32_t>& out) {
	std::vector<std::pair<mySet::iterator, size_t /*current index*/>> vsets;
	out.resize(sets.begin()->size());

	auto sort_remaining = [&] {std::sort(vsets.begin(), vsets.end(),
				[](const std::pair<mySet::iterator, size_t>& lhs,
						const std::pair<mySet::iterator, size_t>& rhs) {
					return (lhs.first->size() - lhs.second) <
					(rhs.first->size() - rhs.second);
				});};

	for (auto it = sets.begin(); it != sets.end(); it++)
		vsets.emplace_back(it, 0);

	auto value = out.begin();
	size_t count = 0, currentset = 1;

	while (vsets[0].second != vsets[0].first->size()) {
		//simd gallop
		vsets[currentset].second += FINDFUNCTION(
				vsets[0].first->at(vsets[0].second),
				vsets[currentset].first->data() + vsets[currentset].second,
				vsets[currentset].first->size() - 1 - vsets[currentset].second);

		if (_LIKELY(
				vsets[currentset].first->at(vsets[currentset].second)
						> vsets[0].first->at(vsets[0].second))) {
			// greater
			vsets[0].second++;
			sort_remaining();
			currentset = 1;
		} else if (vsets[currentset].first->at(vsets[currentset].second)
				== vsets[0].first->at(vsets[0].second)) {
			// equal
			vsets[currentset].second++;
			if (++currentset == vsets.size()) {
				*value++ = vsets[0].first->at(vsets[0].second);
				count++;
				vsets[0].second++;
				// roll over
				sort_remaining();
				currentset = 1;
			}
		} else {
			// less
			break;
		}
	}
	out.resize(count);
}

/* for now we only implement its scalar version */
void adaptive_scalar(const mySet &sets, std::vector<uint32_t> &out) {
	out.resize(sets.begin()->size());
	auto value = out.begin();
	*value = sets.begin()->at(0);
	size_t count = 0, spansize = 1, elimset = 0, currentset = 1;
	mySet::iterator it = sets.begin();
	it++;
	std::vector<size_t> index(sets.size(), 0);

	// Here we do not adopt the following conditional statement
	// "while (index[currentset] <= it->size() - 1 || it->back() >= *value)"
	// because one condition is indispensable inside the loop, assign it
	// on the precise position will reduce trivial repetitive judgements.
	while (true) {
		if (it->at(index[currentset] + spansize) >= *value) {
			//binary search
			// FIXME: note spansize/2 may not be optimal since
			// spansize is finally clamped by the size of sequence
			// rather than simply doubled.
//			index[currentset] = binarySearch_wider(it->data(),
//					index[currentset] + spansize / 2,
//					index[currentset] + spansize, *value);
			index[currentset] += spansize / 2
					+ scalarBinarySearch(*value,
							it->data() + index[currentset] + spansize / 2,
							(spansize + 1) / 2);

			if (it->at(index[currentset]) == *value) {
				// found
				index[currentset]++;
				currentset = ++currentset % sets.size();
				it = ++it == sets.end() ? sets.begin() : it;
				if (currentset == elimset) {
					count++;
					// update eliminator and move to next set
					if (_UNLIKELY(index[currentset] == it->size() - 1))
						break;
					*++value = it->at(++index[currentset]);

					currentset = ++currentset % sets.size();
					it = ++it == sets.end() ? sets.begin() : it;
				}
			} else {
				// not found
				*value = it->at(index[currentset]);
				index[elimset]++;
				elimset = currentset;

				currentset = ++currentset % sets.size();
				it = ++it == sets.end() ? sets.begin() : it;
			}
			// here spansize is set to 0 to compare the last
			// element of the current set
			if (_LIKELY(index[currentset] + 1 < it->size()))
				spansize = 1;
			else if (index[currentset] + 1 == it->size())
				spansize = 0;
			else
				break;
			// here redirect to the beginning rather than
			// altering @spansize
			continue;
		} else if (_UNLIKELY(it->back() < *value))
			break;
		spansize =
				index[currentset] + (spansize << 1) < it->size() ?
						spansize << 1 : it->size() - index[currentset] - 1;
	}
	out.resize(count);
}

/* for now we only implement its scalar version */
void small_adaptive_scalar(const mySet &sets, std::vector<uint32_t> &out) {
	std::vector<std::pair<mySet::iterator, size_t /*current index*/>> vsets;
	out.resize(sets.begin()->size());

	auto sort_remaining = [&] {std::sort(vsets.begin(), vsets.end(),
				[](const std::pair<mySet::iterator, size_t>& lhs,
						const std::pair<mySet::iterator, size_t>& rhs) {
					return (lhs.first->size() - lhs.second) <
					(rhs.first->size() - rhs.second);
				});};

	for (auto it = sets.begin(); it != sets.end(); it++)
		vsets.emplace_back(it, 0);

	auto value = out.begin();
	size_t count = 0, spansize = 1, currentset = 1;

	while (vsets[0].second != vsets[0].first->size()) {
		if (vsets[currentset].first->at(vsets[currentset].second + spansize)
				>= vsets[0].first->at(vsets[0].second)) {
			// galloping overshot
			// FIXME: note spansize/2 may not be optimal since
			// spansize is finally clamped by the size of sequence
			// rather than simply doubled.
//			vsets[currentset].second = binarySearch_wider(
//					vsets[currentset].first->data(),
//					vsets[currentset].second + spansize / 2,
//					vsets[currentset].second + spansize,
//					vsets[0].first->at(vsets[0].second));
			vsets[currentset].second += spansize / 2
					+ scalarBinarySearch(vsets[0].first->at(vsets[0].second),
							vsets[currentset].first->data()
									+ vsets[currentset].second + spansize / 2,
							(spansize + 1) / 2);

			if (vsets[currentset].first->at(vsets[currentset].second)
					== vsets[0].first->at(vsets[0].second)) {
				// find the candidate
				vsets[currentset].second++;
				if (++currentset == vsets.size()) {
					// verify the full intersection
					*value++ = vsets[0].first->at(vsets[0].second);
					count++;
					vsets[0].second++;
					// roll over
					sort_remaining();
					currentset = 1;
				}
			} else {
				// doesn't find the candidate
				// reorder the sets and update candidate
//				vsets[currentset].second++;
				vsets[0].second++;
				sort_remaining();
				currentset = 1;
			}
			// here spansize is set to 0 to compare the last
			// element of the current set
			if (_LIKELY(
					vsets[currentset].second + 1
							< vsets[currentset].first->size()))
				spansize = 1;
			else if (vsets[currentset].second + 1
					== vsets[currentset].first->size())
				spansize = 0;
			else
				break;
			continue;
		} else if (_UNLIKELY(
				vsets[currentset].first->back()
						< vsets[0].first->at(vsets[0].second)))
			break;
		spansize =
				vsets[currentset].second + (spansize << 1)
						< vsets[currentset].first->size() ?
						spansize << 1 :
						vsets[currentset].first->size()
								- vsets[currentset].second - 1;
	}
	out.resize(count);
}

template<intersectionfindfunction FINDFUNCTION>
void sequential(const mySet &sets, std::vector<uint32_t> &out) {
	out.resize(sets.begin()->size());
	auto value = out.begin();
	*value = sets.begin()->at(0);
	size_t count = 0, elimset = 0, currentset = 1;
	mySet::iterator it = sets.begin();
	mySet::iterator it_elim = it++;
	std::vector<size_t> index(sets.size(), 0);

	while (_LIKELY(index[currentset] < it->size())) {
		//SIMD gallop
		index[currentset] += FINDFUNCTION(*value,
				it->data() + index[currentset],
				it->size() - 1 - index[currentset]);
		if (it->at(index[currentset]) == *value) {
			index[currentset]++;
			currentset = ++currentset % sets.size();
			it = ++it == sets.end() ? sets.begin() : it;
			if (currentset == elimset) {
				count++;
				// update eliminator and move to next set
				if (_UNLIKELY(++index[currentset] == it->size()))
					break;
				*++value = it->at(index[currentset]);

				currentset = ++currentset % sets.size();
				it = ++it == sets.end() ? sets.begin() : it;
			}
		} else if (_UNLIKELY(it->back() < *value)) {
			break;
		} else if (_UNLIKELY(++index[elimset] == it_elim->size()))
			break;
		else {
			// greater, however reaches the end.
			// XXX: following two sentences unreachable
//			if (_UNLIKELY(index[currentset] + 1 > it->size()))
//				break;
			// greater and doesn't reach the end
			// go on intersecting
			*value = it->at(index[currentset]);
			elimset = currentset;
			it_elim = it;

			currentset = ++currentset % sets.size();
			it = ++it == sets.end() ? sets.begin() : it;
		}
	}
	out.resize(count);
}

template<intersectionfindfunction FINDFUNCTION>
void max(const mySet &sets, std::vector<uint32_t> &out) {
	out.resize(sets.begin()->size());
	auto value = out.begin();
	*value = sets.begin()->at(0);
	size_t count = 0, intersect_count = 0, elimset = 0, currentset = 1;
	mySet::iterator it = sets.begin();
	mySet::iterator it_start = it++, it_elim = it_start;
	std::vector<size_t> index(sets.size(), 0);

	while (_LIKELY(index[currentset] < it->size())) {
		//SIMD gallop
		index[currentset] += FINDFUNCTION(*value,
				it->data() + index[currentset],
				it->size() - 1 - index[currentset]);

		if (it->at(index[currentset]) == *value) {
			intersect_count++;
			index[currentset]++;
			++currentset;
			++it;
//			currentset = ++currentset % sets.size();
//			it = ++it == sets.end() ? sets.begin() : it;
			if (intersect_count == sets.size() - 1) {
				// ensure positions of all the sets advance 1
				index[elimset]++;
				count++;
				// update eliminator from set 0 and move to set 1
				if (_UNLIKELY(index[0] == it_start->size()))
					break;
				intersect_count = 0;
				*++value = it_start->at(index[0]);
				currentset = 1;
				elimset = 0;
				it_elim = it = it_start;
				it++;
			} else if (currentset == elimset) {
				++currentset;
				++it;
//				currentset = ++currentset % sets.size();
//				it = ++it == sets.end() ? sets.begin() : it;
			}
		} else if (_UNLIKELY(it->back() < *value)) {
			break;
		} else if (_UNLIKELY(++index[elimset] == it_elim->size()))
			break;
		else {
			intersect_count = 0;
			if (currentset == 0
					|| it_start->at(index[0]) > it->at(index[currentset])) {
				*value = it_start->at(index[0]);
				elimset = 0;
				currentset = 1;

				it_elim = it = it_start;
				it++;
			} else {
				*value = it->at(index[currentset]);
				elimset = currentset;
				currentset = 0;

				it_elim = it;
				it = it_start;
			}
		}
	}
	out.resize(count);
}

/*
 * @description: basic methods intersecting two sorted sequences,
 * returned results are also sorted.
 */
// further optimization can be search @freq[0] and @freq[freq_end] in @rare,
// which helps target the real overlap of @rare and @freq
template<intersectionfindfunction FINDFUNCTION>
void BYintersect_sorted(const uint32_t *freq, const size_t &freq_end,
		const uint32_t *rare, const size_t &rare_end, uint32_t **out,
		uint32_t &count) {
	if (_LIKELY(
			freq_end == -1 || rare_end == -1 || freq[0] > rare[rare_end]
					|| rare[0] > freq[freq_end]))
		return;
	else {
		size_t rare_mid = rare_end / 2;
//		size_t freq_mid = msis::binarySearch_wider(freq, 0, freq_end,
//				rare[rare_mid]);
		size_t freq_mid = FINDFUNCTION(rare[rare_mid], freq, freq_end + 1);
		if (freq_mid > rare_mid)
			BYintersect_sorted<FINDFUNCTION>(freq, freq_mid - 1, rare,
					rare_mid - 1, out, count);
		else
			BYintersect_sorted<FINDFUNCTION>(rare, rare_mid - 1, freq,
					freq_mid - 1, out, count);
		if (freq[freq_mid] == rare[rare_mid++]) {
			*(*out)++ = freq[freq_mid++];
			count++;
		}
//		if (_UNLIKELY(freq_mid == freq_end || rare_mid == rare_end)) {
//			if (freq[freq_mid] == rare[rare_mid]) {
//				*(*out)++ = freq[freq_mid];
//				i++;
//			}
//		} else {
		if (freq_end - freq_mid > rare_end - rare_mid)
			BYintersect_sorted<FINDFUNCTION>(freq + freq_mid,
					freq_end - freq_mid, rare + rare_mid, rare_end - rare_mid,
					out, count);
		else
			BYintersect_sorted<FINDFUNCTION>(rare + rare_mid,
					rare_end - rare_mid, freq + freq_mid, freq_end - freq_mid,
					out, count);
//		}
	}
}

template<intersectionfindfunction FINDFUNCTION>
void BaezaYates(const mySet &sets, std::vector<uint32_t> &out) {
	mySet::iterator it = sets.begin();
	std::vector<uint32_t> intersection(std::move(*it++));

	uint32_t* out_init = intersection.data();
	uint32_t** pout = &out_init;

	for (; it != sets.end(); it++) {
		uint32_t count = 0;

		BYintersect_sorted<FINDFUNCTION>(it->data(), it->size() - 1,
				intersection.data(), intersection.size() - 1, pout, count);
		intersection.resize(count);

		out_init = intersection.data();
		pout = &out_init;
	}
	out.swap(intersection);
}

}
#endif /* INCLUDE_MULTISETINTERSECTION_HPP_ */
