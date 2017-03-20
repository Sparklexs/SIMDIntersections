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
const UINT4 SIMDWIDTH = 4;
const UINT4 V4_BLOCKSIZE = SIMDWIDTH * 1;
const UINT4 V8_BLOCKSIZE = SIMDWIDTH * 2;
const UINT4 V16_BLOCKSIZE = SIMDWIDTH * 4;
const UINT4 V32_BLOCKSIZE = SIMDWIDTH * 8;
const UINT4 V64_BLOCKSIZE = SIMDWIDTH * 16;
const UINT4 V128_BLOCKSIZE = SIMDWIDTH * 32;
const UINT4 V256_BLOCKSIZE = SIMDWIDTH * 64;
const UINT4 V512_BLOCKSIZE = SIMDWIDTH * 128;
const UINT4 V4_POS[] = { 0, 0, 1, 0, 2, 0, 0, 0, 3 };

typedef size_t (*setIntersectionFunction)(const uint32_t * set1,
		const size_t length1, const uint32_t * set2, const size_t length2,
		uint32_t *out); // used only for method `svs`

// here adapt the range [start(0),end], different from __BSadvanceUntil
// whose range is [start+1,end-1]
long scalarBinarySearch(UINT4 min, const UINT4 * array, long end) {
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
long scalarGallop(UINT4 min, const UINT4 * array, long end) {
	// special handling for a possibly common sql_exact case
	if ((0 >= end) or (*array >= min)) {
		return 0;
	}

	long upper = 1; // could set larger
	// bootstrap an upper limit

	// sxs: here spansize is enlarged to the maximum where its corresponding
	// element is geq min
	while ((upper <= end) and (array[upper] < min))
		upper <<= 1;

	// we know that the next-smallest span was too small
	long lower = (upper / 2);

	if (upper > end)
		upper = end;

	if (array[upper] <= min) {    // means array has no item >= min
		return upper;
	}

	// else begin binary search
	long mid = 0;
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

/*
 * some of the following functions can be furhter merged by using
 * templates, but not all of them. To keep the readability, templates
 * are discarded.
 */

long simdlinear_v4_exact(UINT4 goal, const UINT4 *target, long ntargets) {
	const UINT4 *end_target, *stop_target, *init_target;

	init_target = target;
	stop_target = &(target[ntargets - V4_BLOCKSIZE]);
	end_target = &(target[ntargets]);

	if (_UNLIKELY(target >= stop_target)) {
		return Intersection_find_scalar(goal, target, ntargets);
	}

	while (target[V4_BLOCKSIZE - 1] < goal) {
		target += V4_BLOCKSIZE;
		if (target >= stop_target) {
			return (target - init_target)
					+ Intersection_find_scalar(goal, target,/*ntargets*/
					(end_target - target));
		}
	}

	// now search the block [target, target+block_size-1]
	__m128i Match = _mm_set1_epi32(goal);
	UINT4 pos = V4_POS[_mm_movemask_ps(
			_mm_castsi128_ps(
					_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 0),
							Match)))];

	return (target - init_target) + pos;
}

long simdlinear_v4_rough(int*foundp, UINT4 goal, const UINT4 *target,
		long ntargets) {
	const UINT4 *stop_target, *init_target;
	long pos;

	init_target = target;
	stop_target = &(target[ntargets - V4_BLOCKSIZE]);

	if (_UNLIKELY(target >= stop_target)) {
		if ((pos = Intersection_find_scalar(goal, target, ntargets)) <= ntargets
				&& target[pos] == goal)
			*foundp = 1;
		else
			*foundp = 0;
		return pos;
	}

	while (target[V4_BLOCKSIZE - 1] < goal) {
		target += V4_BLOCKSIZE;
		if (target >= stop_target) {
			if ((pos = Intersection_find_scalar(goal, target, ntargets))
					<= ntargets && target[pos] == goal)
				*foundp = 1;
			else
				*foundp = 0;

			return (target - init_target) + pos;
		}
	}

	// now search the block [target, target+block_size-1]
	__m128i Match = _mm_set1_epi32(goal);

	__m128i F0 = _mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 0),
			Match);
	if (_mm_testz_si128(F0, F0))
		*foundp = 0;
	else
		*foundp = 1;
	return (target - init_target);
}

long simdlinear_v8_exact(UINT4 goal, const UINT4 *target, long ntargets) {
	const UINT4 *end_target, *stop_target, *init_target;

	init_target = target;
	stop_target = &(target[ntargets - V8_BLOCKSIZE]);
	end_target = &(target[ntargets]);

	if (_UNLIKELY(target >= stop_target)) {
		return Intersection_find_scalar(goal, target, ntargets);
	}

	while (target[V8_BLOCKSIZE - 1] < goal) {
		target += V8_BLOCKSIZE;
		if (target >= stop_target) {
			return (target - init_target)
					+ Intersection_find_scalar(goal, target,/*ntargets*/
					(end_target - target));
		}
	}

	// now search the block [target, target+block_size-1]
	__m128i Match = _mm_set1_epi32(goal);
	long pos;

	__m128i F0 = _mm_or_si128(
			_mm_srli_epi32(
					_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 0),
							Match), 32 - 1),
			_mm_srli_epi32(
					_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 1),
							Match), 28 - 1));

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

long simdlinear_v8_rough(int*foundp, UINT4 goal, const UINT4 *target,
		long ntargets) {
	const UINT4 *stop_target, *init_target;
	long pos;

	init_target = target;
	stop_target = &(target[ntargets - V8_BLOCKSIZE]);

	if (_UNLIKELY(target >= stop_target)) {
		if ((pos = Intersection_find_scalar(goal, target, ntargets)) <= ntargets
				&& target[pos] == goal)
			*foundp = 1;
		else
			*foundp = 0;
		return pos;
	}

	while (target[V8_BLOCKSIZE - 1] < goal) {
		target += V8_BLOCKSIZE;
		if (target >= stop_target) {
			if ((pos = Intersection_find_scalar(goal, target, ntargets))
					<= ntargets && target[pos] == goal)
				*foundp = 1;
			else
				*foundp = 0;

			return (target - init_target) + pos;
		}
	}

	// now search the block [target, target+block_size-1]
	__m128i Match = _mm_set1_epi32(goal);

	__m128i F0 = _mm_or_si128(
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 0), Match),
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 1), Match));
	if (_mm_testz_si128(F0, F0))
		*foundp = 0;
	else
		*foundp = 1;
	return (target - init_target);
}

long simdlinear_v8_rough_plow(int*foundp, UINT4 goal, const UINT4 *target,
		long ntargets) {
	const UINT4 *stop_target, *init_target;
	long pos;

	init_target = target;
	stop_target = &(target[ntargets - V8_BLOCKSIZE]);

	if (_UNLIKELY(target >= stop_target)) {
		if ((pos = Intersection_find_scalar(goal, target, ntargets)) <= ntargets
				&& target[pos] == goal)
			*foundp = 1;
		else
			*foundp = 0;
		return pos;
	}

	while (target[V8_BLOCKSIZE - 1] < goal) {
		target += V8_BLOCKSIZE;
		if (target >= stop_target) {
			if ((pos = Intersection_find_scalar(goal, target, ntargets))
					<= ntargets && target[pos] == goal)
				*foundp = 1;
			else
				*foundp = 0;

			return (target - init_target) + pos;
		}
	}

	// now search the block [target, target+block_size-1]
	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0 = _mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 0),
			Match);
	UINT4 i = 0;
	while (_mm_testz_si128(F0, F0) && ++i < V8_BLOCKSIZE / SIMDWIDTH) {
		F0 = _mm_or_si128(F0,
				_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + i),
						Match));
	}
	if (_mm_testz_si128(F0, F0)) {
		*foundp = 0;
		i--;
	} else
		*foundp = 1;
	return (target + i * SIMDWIDTH - init_target);
}

long simdlinear_v16_exact(UINT4 goal, const UINT4 *target, long ntargets) {
	const UINT4 *end_target, *stop_target, *init_target;

	init_target = target;
	stop_target = &(target[ntargets - V16_BLOCKSIZE]);
	end_target = &(target[ntargets]);

	if (_UNLIKELY(target >= stop_target)) {
		return Intersection_find_scalar(goal, target, ntargets);
	}

	while (target[V16_BLOCKSIZE - 1] < goal) {
		target += V16_BLOCKSIZE;
		if (target >= stop_target) {
			return (target - init_target)
					+ Intersection_find_scalar(goal, target,/*ntargets*/
					(end_target - target));
		}
	}

	// now search the block [target, target+block_size-1]
	__m128i Match = _mm_set1_epi32(goal);
	long pos;

	__m128i F0 = _mm_or_si128(
			_mm_or_si128(
					_mm_srli_epi32(
							_mm_cmpeq_epi32(
									_mm_lddqu_si128((__m128i *) target + 0),
									Match), 32 - 1),
					_mm_srli_epi32(
							_mm_cmpeq_epi32(
									_mm_lddqu_si128((__m128i *) target + 1),
									Match), 28 - 1)),
			_mm_or_si128(
					_mm_srli_epi32(
							_mm_cmpeq_epi32(
									_mm_lddqu_si128((__m128i *) target + 2),
									Match), 24 - 1),
					_mm_srli_epi32(
							_mm_cmpeq_epi32(
									_mm_lddqu_si128((__m128i *) target + 3),
									Match), 20 - 1)));

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

long simdlinear_v16_rough(int*foundp, UINT4 goal, const UINT4 *target,
		long ntargets) {
	const UINT4 *stop_target, *init_target;
	long pos;

	init_target = target;
	stop_target = &(target[ntargets - V16_BLOCKSIZE]);

	if (_UNLIKELY(target >= stop_target)) {
		if ((pos = Intersection_find_scalar(goal, target, ntargets)) <= ntargets
				&& target[pos] == goal)
			*foundp = 1;
		else
			*foundp = 0;
		return pos;
	}

	while (target[V16_BLOCKSIZE - 1] < goal) {
		target += V16_BLOCKSIZE;
		if (target >= stop_target) {
			if ((pos = Intersection_find_scalar(goal, target, ntargets))
					<= ntargets && target[pos] == goal)
				*foundp = 1;
			else
				*foundp = 0;

			return (target - init_target) + pos;
		}
	}

	// now search the block [target, target+block_size-1]
	__m128i Match = _mm_set1_epi32(goal);

	__m128i F0 = _mm_or_si128(
			_mm_or_si128(
					_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 0),
							Match),
					_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 1),
							Match)),
			_mm_or_si128(
					_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 2),
							Match),
					_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 3),
							Match)));
	if (_mm_testz_si128(F0, F0))
		*foundp = 0;
	else
		*foundp = 1;
	return (target - init_target);
}

long simdlinear_v16_rough_plow(int*foundp, UINT4 goal, const UINT4 *target,
		long ntargets) {
	const UINT4 *stop_target, *init_target;
	long pos;

	init_target = target;
	stop_target = &(target[ntargets - V16_BLOCKSIZE]);

	if (_UNLIKELY(target >= stop_target)) {
		if ((pos = Intersection_find_scalar(goal, target, ntargets)) <= ntargets
				&& target[pos] == goal)
			*foundp = 1;
		else
			*foundp = 0;
		return pos;
	}

	while (target[V16_BLOCKSIZE - 1] < goal) {
		target += V16_BLOCKSIZE;
		if (target >= stop_target) {
			if ((pos = Intersection_find_scalar(goal, target, ntargets))
					<= ntargets && target[pos] == goal)
				*foundp = 1;
			else
				*foundp = 0;

			return (target - init_target) + pos;
		}
	}

	// now search the block [target, target+block_size-1]
	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0 = _mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 0),
			Match);
	UINT4 i = 0;
	while (_mm_testz_si128(F0, F0) && ++i < V16_BLOCKSIZE / SIMDWIDTH) {
		F0 = _mm_or_si128(F0,
				_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + i),
						Match));
	}
	if (_mm_testz_si128(F0, F0)) {
		*foundp = 0;
		i--;
	} else
		*foundp = 1;
	return (target + i * SIMDWIDTH - init_target);
}

long simdlinear_v32_exact(UINT4 goal, const UINT4 *target, long ntargets) {
	const UINT4 *end_target, *stop_target, *init_target;

	init_target = target;
	stop_target = &(target[ntargets - V32_BLOCKSIZE]);
	end_target = &(target[ntargets]);

	if (_UNLIKELY(target >= stop_target)) {
		return Intersection_find_scalar(goal, target, ntargets);
	}

	while (target[V32_BLOCKSIZE - 1] < goal) {
		target += V32_BLOCKSIZE;
		if (target >= stop_target) {
			return (target - init_target)
					+ Intersection_find_scalar(goal, target,/*ntargets*/
					(end_target - target));
		}
	}

	// now search the block [target, target+block_size-1]
	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0, Q0, Q1, Q2, Q3;
	long pos;

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

long simdlinear_v32_rough(int*foundp, UINT4 goal, const UINT4 *target,
		long ntargets) {
	const UINT4 *stop_target, *init_target;
	long pos;

	init_target = target;
	stop_target = &(target[ntargets - V32_BLOCKSIZE]);

	if (_UNLIKELY(target >= stop_target)) {
		if ((pos = Intersection_find_scalar(goal, target, ntargets)) <= ntargets
				&& target[pos] == goal)
			*foundp = 1;
		else
			*foundp = 0;
		return pos;
	}

	while (target[V32_BLOCKSIZE - 1] < goal) {
		target += V32_BLOCKSIZE;
		if (target >= stop_target) {
			if ((pos = Intersection_find_scalar(goal, target, ntargets))
					<= ntargets && target[pos] == goal)
				*foundp = 1;
			else
				*foundp = 0;

			return (target - init_target) + pos;
		}
	}

	// now search the block [target, target+block_size-1]
	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0, Q0, Q1, Q2, Q3;

	Q0 = _mm_or_si128(
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 0), Match),
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 1), Match));
	Q1 = _mm_or_si128(
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 2), Match),
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 3), Match));
	Q2 = _mm_or_si128(
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 4), Match),
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 5), Match));
	Q3 = _mm_or_si128(
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 6), Match),
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 7), Match));

	F0 = _mm_or_si128(_mm_or_si128(Q0, Q1), _mm_or_si128(Q2, Q3));
	if (_mm_testz_si128(F0, F0))
		*foundp = 0;
	else
		*foundp = 1;
	return (target - init_target);
}

long simdlinear_v32_rough_plow(int*foundp, UINT4 goal, const UINT4 *target,
		long ntargets) {
	const UINT4 *stop_target, *init_target;
	long pos;

	init_target = target;
	stop_target = &(target[ntargets - V32_BLOCKSIZE]);

	if (_UNLIKELY(target >= stop_target)) {
		if ((pos = Intersection_find_scalar(goal, target, ntargets)) <= ntargets
				&& target[pos] == goal)
			*foundp = 1;
		else
			*foundp = 0;
		return pos;
	}

	while (target[V32_BLOCKSIZE - 1] < goal) {
		target += V32_BLOCKSIZE;
		if (target >= stop_target) {
			if ((pos = Intersection_find_scalar(goal, target, ntargets))
					<= ntargets && target[pos] == goal)
				*foundp = 1;
			else
				*foundp = 0;

			return (target - init_target) + pos;
		}
	}

	// now search the block [target, target+block_size-1]
	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0 = _mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 0),
			Match);
	UINT4 i = 0;
	while (_mm_testz_si128(F0, F0) && ++i < V32_BLOCKSIZE / SIMDWIDTH) {
		F0 = _mm_or_si128(F0,
				_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + i),
						Match));
	}
	if (_mm_testz_si128(F0, F0)) {
		*foundp = 0;
		i--;
	} else
		*foundp = 1;
	return (target + i * SIMDWIDTH - init_target);
}

long simdlinear_v64_exact(UINT4 goal, const UINT4 *target, long ntargets) {
	const UINT4 *end_target, *stop_target, *init_target;

	init_target = target;
	stop_target = &(target[ntargets - V64_BLOCKSIZE]);
	end_target = &(target[ntargets]);

	if (_UNLIKELY(target >= stop_target)) {
		return Intersection_find_scalar(goal, target, ntargets);
	}

	while (target[V64_BLOCKSIZE - 1] < goal) {
		target += V64_BLOCKSIZE;
		if (target >= stop_target) {
			return (target - init_target)
					+ Intersection_find_scalar(goal, target,/*ntargets*/
					(end_target - target));
		}
	}

	// now search the block [target, target+block_size-1]
	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0, Q0, Q1, Q2, Q3;
	long pos;

	if (target[V64_BLOCKSIZE / 2 - 1] < goal)
		target += 32;

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

long simdlinear_v64_rough(int*foundp, UINT4 goal, const UINT4 *target,
		long ntargets) {
	const UINT4 *stop_target, *init_target;
	long pos;

	init_target = target;
	stop_target = &(target[ntargets - V64_BLOCKSIZE]);

	if (_UNLIKELY(target >= stop_target)) {
		if ((pos = Intersection_find_scalar(goal, target, ntargets)) <= ntargets
				&& target[pos] == goal)
			*foundp = 1;
		else
			*foundp = 0;
		return pos;
	}

	while (target[V64_BLOCKSIZE - 1] < goal) {
		target += V64_BLOCKSIZE;
		if (target >= stop_target) {
			if ((pos = Intersection_find_scalar(goal, target, ntargets))
					<= ntargets && target[pos] == goal)
				*foundp = 1;
			else
				*foundp = 0;

			return (target - init_target) + pos;
		}
	}

	// now search the block [target, target+block_size-1]
	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0, Q0, Q1, Q2, Q3;

	if (target[V64_BLOCKSIZE / 2 - 1] < goal)
		target += 32;

	Q0 = _mm_or_si128(
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 0), Match),
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 1), Match));
	Q1 = _mm_or_si128(
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 2), Match),
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 3), Match));
	Q2 = _mm_or_si128(
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 4), Match),
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 5), Match));
	Q3 = _mm_or_si128(
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 6), Match),
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 7), Match));

	F0 = _mm_or_si128(_mm_or_si128(Q0, Q1), _mm_or_si128(Q2, Q3));
	if (_mm_testz_si128(F0, F0))
		*foundp = 0;
	else
		*foundp = 1;
	return (target - init_target);
}

long simdlinear_v64_rough_plow(int*foundp, UINT4 goal, const UINT4 *target,
		long ntargets) {
	const UINT4 *stop_target, *init_target;
	long pos;

	init_target = target;
	stop_target = &(target[ntargets - V64_BLOCKSIZE]);

	if (_UNLIKELY(target >= stop_target)) {
		if ((pos = Intersection_find_scalar(goal, target, ntargets)) <= ntargets
				&& target[pos] == goal)
			*foundp = 1;
		else
			*foundp = 0;
		return pos;
	}

	while (target[V64_BLOCKSIZE - 1] < goal) {
		target += V64_BLOCKSIZE;
		if (target >= stop_target) {
			if ((pos = Intersection_find_scalar(goal, target, ntargets))
					<= ntargets && target[pos] == goal)
				*foundp = 1;
			else
				*foundp = 0;

			return (target - init_target) + pos;
		}
	}

	// now search the block [target, target+block_size-1]
	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0 = _mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 0),
			Match);
	UINT4 i = 0;
	while (_mm_testz_si128(F0, F0) && ++i < V64_BLOCKSIZE / SIMDWIDTH) {
		F0 = _mm_or_si128(F0,
				_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + i),
						Match));
	}
	if (_mm_testz_si128(F0, F0)) {
		*foundp = 0;
		i--;
	} else
		*foundp = 1;
	return (target + i * SIMDWIDTH - init_target);
}

long simdlinear_v128_exact(UINT4 goal, const UINT4 *target, long ntargets) {
	const UINT4 *end_target, *stop_target, *init_target;

	init_target = target;
	stop_target = &(target[ntargets - V128_BLOCKSIZE]);
	end_target = &(target[ntargets]);

	if (_UNLIKELY(target >= stop_target)) {
		return Intersection_find_scalar(goal, target, ntargets);
	}

	while (target[V128_BLOCKSIZE - 1] < goal) {
		target += V128_BLOCKSIZE;
		if (target >= stop_target) {
			return (target - init_target)
					+ Intersection_find_scalar(goal, target,/*ntargets*/
					(end_target - target));
		}
	}

	// now search the block [target, target+block_size-1]
	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0, Q0, Q1, Q2, Q3;
	long pos;

	if (target[V128_BLOCKSIZE / 2 - 1] >= goal) {
		if (target[V128_BLOCKSIZE / 4 - 1] < goal)
			target += 32;
	} else {
		if (target[V128_BLOCKSIZE * 3 / 4 - 1] >= goal)
			target += 64;
		else
			target += 96;
	}

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

long simdlinear_v128_rough(int*foundp, UINT4 goal, const UINT4 *target,
		long ntargets) {
	const UINT4 *stop_target, *init_target;
	long pos;

	init_target = target;
	stop_target = &(target[ntargets - V128_BLOCKSIZE]);

	if (_UNLIKELY(target >= stop_target)) {
		if ((pos = Intersection_find_scalar(goal, target, ntargets)) <= ntargets
				&& target[pos] == goal)
			*foundp = 1;
		else
			*foundp = 0;
		return pos;
	}

	while (target[V128_BLOCKSIZE - 1] < goal) {
		target += V128_BLOCKSIZE;
		if (target >= stop_target) {
			if ((pos = Intersection_find_scalar(goal, target, ntargets))
					<= ntargets && target[pos] == goal)
				*foundp = 1;
			else
				*foundp = 0;

			return (target - init_target) + pos;
		}
	}

	// now search the block [target, target+block_size-1]
	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0, Q0, Q1, Q2, Q3;

	if (target[V128_BLOCKSIZE / 2 - 1] >= goal) {
		if (target[V128_BLOCKSIZE / 4 - 1] < goal)
			target += 32;
	} else {
		if (target[V128_BLOCKSIZE * 3 / 4 - 1] >= goal)
			target += 64;
		else
			target += 96;
	}

	Q0 = _mm_or_si128(
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 0), Match),
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 1), Match));
	Q1 = _mm_or_si128(
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 2), Match),
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 3), Match));
	Q2 = _mm_or_si128(
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 4), Match),
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 5), Match));
	Q3 = _mm_or_si128(
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 6), Match),
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 7), Match));

	F0 = _mm_or_si128(_mm_or_si128(Q0, Q1), _mm_or_si128(Q2, Q3));
	if (_mm_testz_si128(F0, F0))
		*foundp = 0;
	else
		*foundp = 1;
	return (target - init_target);
}

long simdlinear_v128_rough_plow(int*foundp, UINT4 goal, const UINT4 *target,
		long ntargets) {
	const UINT4 *stop_target, *init_target;
	long pos;

	init_target = target;
	stop_target = &(target[ntargets - V128_BLOCKSIZE]);

	if (_UNLIKELY(target >= stop_target)) {
		if ((pos = Intersection_find_scalar(goal, target, ntargets)) <= ntargets
				&& target[pos] == goal)
			*foundp = 1;
		else
			*foundp = 0;
		return pos;
	}

	while (target[V128_BLOCKSIZE - 1] < goal) {
		target += V128_BLOCKSIZE;
		if (target >= stop_target) {
			if ((pos = Intersection_find_scalar(goal, target, ntargets))
					<= ntargets && target[pos] == goal)
				*foundp = 1;
			else
				*foundp = 0;

			return (target - init_target) + pos;
		}
	}

	// now search the block [target, target+block_size-1]
	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0 = _mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 0),
			Match);
	UINT4 i = 0;
	while (_mm_testz_si128(F0, F0) && ++i < V128_BLOCKSIZE / SIMDWIDTH) {
		F0 = _mm_or_si128(F0,
				_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + i),
						Match));
	}
	if (_mm_testz_si128(F0, F0)) {
		*foundp = 0;
		i--;
	} else
		*foundp = 1;
	return (target + i * SIMDWIDTH - init_target);
}

long simdlinear_v256_exact(UINT4 goal, const UINT4 *target, long ntargets) {
	const UINT4 *end_target, *stop_target, *init_target;

	init_target = target;
	stop_target = &(target[ntargets - V256_BLOCKSIZE]);
	end_target = &(target[ntargets]);

	if (_UNLIKELY(target >= stop_target)) {
		return Intersection_find_scalar(goal, target, ntargets);
	}

	while (target[V256_BLOCKSIZE - 1] < goal) {
		target += V256_BLOCKSIZE;
		if (target >= stop_target) {
			return (target - init_target)
					+ Intersection_find_scalar(goal, target,/*ntargets*/
					(end_target - target));
		}
	}

	// now search the block [target, target+block_size-1]
	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0, Q0, Q1, Q2, Q3;
	long pos;

	if (target[V256_BLOCKSIZE / 2 - 1] >= goal) {
		if (target[V256_BLOCKSIZE / 4 - 1] >= goal) {
			if (target[V256_BLOCKSIZE / 8 - 1] < goal)
				target += 32;
		} else {
			if (target[V256_BLOCKSIZE / 8 * 3 - 1] >= goal) {
				target += 64;
			} else
				target += 96;
		}
	} else {
		if (target[V256_BLOCKSIZE / 4 + 128 - 1] >= goal) {
			if (target[V256_BLOCKSIZE / 8 + 128 - 1] >= goal)
				target += 128;
			else
				target += 160;

		} else {
			if (target[V256_BLOCKSIZE / 8 * 3 + 128 - 1] >= goal)
				target += 192;
			else
				target += 224;
		}
	}

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

long simdlinear_v256_rough(int*foundp, UINT4 goal, const UINT4 *target,
		long ntargets) {
	const UINT4 *stop_target, *init_target;
	long pos;

	init_target = target;
	stop_target = &(target[ntargets - V256_BLOCKSIZE]);

	if (_UNLIKELY(target >= stop_target)) {
		if ((pos = Intersection_find_scalar(goal, target, ntargets)) <= ntargets
				&& target[pos] == goal)
			*foundp = 1;
		else
			*foundp = 0;
		return pos;
	}

	while (target[V256_BLOCKSIZE - 1] < goal) {
		target += V256_BLOCKSIZE;
		if (target >= stop_target) {
			if ((pos = Intersection_find_scalar(goal, target, ntargets))
					<= ntargets && target[pos] == goal)
				*foundp = 1;
			else
				*foundp = 0;

			return (target - init_target) + pos;
		}
	}

	// now search the block [target, target+block_size-1]
	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0, Q0, Q1, Q2, Q3;

	if (target[V256_BLOCKSIZE / 2 - 1] >= goal) {
		if (target[V256_BLOCKSIZE / 4 - 1] >= goal) {
			if (target[V256_BLOCKSIZE / 8 - 1] < goal)
				target += 32;
		} else {
			if (target[V256_BLOCKSIZE / 8 * 3 - 1] >= goal)
				target += 64;
			else
				target += 96;
		}
	} else {
		if (target[V256_BLOCKSIZE / 4 + 128 - 1] >= goal) {
			if (target[V256_BLOCKSIZE / 8 + 128 - 1] >= goal)
				target += 128;
			else
				target += 160;
		} else {
			if (target[V256_BLOCKSIZE / 8 * 3 + 128 - 1] >= goal)
				target += 192;
			else
				target += 224;
		}
	}

	Q0 = _mm_or_si128(
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 0), Match),
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 1), Match));
	Q1 = _mm_or_si128(
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 2), Match),
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 3), Match));
	Q2 = _mm_or_si128(
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 4), Match),
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 5), Match));
	Q3 = _mm_or_si128(
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 6), Match),
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 7), Match));

	F0 = _mm_or_si128(_mm_or_si128(Q0, Q1), _mm_or_si128(Q2, Q3));
	if (_mm_testz_si128(F0, F0))
		*foundp = 0;
	else
		*foundp = 1;
	return (target - init_target);
}

long simdlinear_v256_rough_plow(int*foundp, UINT4 goal, const UINT4 *target,
		long ntargets) {
	const UINT4 *stop_target, *init_target;
	long pos;

	init_target = target;
	stop_target = &(target[ntargets - V256_BLOCKSIZE]);

	if (_UNLIKELY(target >= stop_target)) {
		if ((pos = Intersection_find_scalar(goal, target, ntargets)) <= ntargets
				&& target[pos] == goal)
			*foundp = 1;
		else
			*foundp = 0;
		return pos;
	}

	while (target[V256_BLOCKSIZE - 1] < goal) {
		target += V256_BLOCKSIZE;
		if (target >= stop_target) {
			if ((pos = Intersection_find_scalar(goal, target, ntargets))
					<= ntargets && target[pos] == goal)
				*foundp = 1;
			else
				*foundp = 0;

			return (target - init_target) + pos;
		}
	}

	// now search the block [target, target+block_size-1]
	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0 = _mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 0),
			Match);
	UINT4 i = 0;
	while (_mm_testz_si128(F0, F0) && ++i < V256_BLOCKSIZE / SIMDWIDTH) {
		F0 = _mm_or_si128(F0,
				_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + i),
						Match));
	}
	if (_mm_testz_si128(F0, F0)) {
		*foundp = 0;
		i--;
	} else
		*foundp = 1;
	return (target + i * SIMDWIDTH - init_target);
}

long simdlinear_v512_exact(UINT4 goal, const UINT4 *target, long ntargets) {
	const UINT4 *end_target, *stop_target, *init_target;

	init_target = target;
	stop_target = &(target[ntargets - V512_BLOCKSIZE]);
	end_target = &(target[ntargets]);

	if (_UNLIKELY(target >= stop_target)) {
		return Intersection_find_scalar(goal, target, ntargets);
	}

	while (target[V512_BLOCKSIZE - 1] < goal) {
		target += V512_BLOCKSIZE;
		if (target >= stop_target) {
			return (target - init_target)
					+ Intersection_find_scalar(goal, target,/*ntargets*/
					(end_target - target));
		}
	}

	// now search the block [target, target+block_size-1]
	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0, Q0, Q1, Q2, Q3;
	long pos;

	if (target[V512_BLOCKSIZE / 2 - 1] >= goal) {
		if (target[V512_BLOCKSIZE / 4 - 1] >= goal) {
			if (target[V512_BLOCKSIZE / 8 - 1] >= goal) {
				if (target[V512_BLOCKSIZE / 16 - 1] >= goal) {
				} else {
					target += 32;
				}
			} else {
				if (target[V512_BLOCKSIZE / 8 + 32 - 1] >= goal) {
					target += 64;
				} else {
					target += 96;
				}
			}
		} else {
			if (target[V512_BLOCKSIZE / 4 + 64 - 1] >= goal) {
				if (target[V512_BLOCKSIZE / 8 + 96 - 1] >= goal) {
					target += 128;
				} else {
					target += 160;
				}
			} else {
				if (target[V512_BLOCKSIZE / 4 + 96 - 1] >= goal) {
					target += 192;
				} else {
					target += 224;
				}
			}
		}
	} else {
		if (target[V512_BLOCKSIZE / 4 + 256 - 1] >= goal) {
			if (target[V512_BLOCKSIZE / 8 + 256 - 1] >= goal) {
				if (target[V512_BLOCKSIZE / 16 + 256 - 1] >= goal) {
					target += 256;
				} else {
					target += 288;
				}
			} else {
				if (target[V512_BLOCKSIZE / 8 + 288 - 1] >= goal) {
					target += 320;
				} else {
					target += 352;
				}
			}
		} else {
			if (target[V512_BLOCKSIZE / 4 + 320 - 1] >= goal) {
				if (target[V512_BLOCKSIZE / 8 + 352 - 1] >= goal) {
					target += 384;
				} else {
					target += 416;
				}
			} else {
				if (target[V512_BLOCKSIZE / 4 + 352 - 1] >= goal) {
					target += 448;
				} else {
					target += 480;
				}
			}
		}
	}

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

long simdlinear_v512_rough(int*foundp, UINT4 goal, const UINT4 *target,
		long ntargets) {
	const UINT4 *stop_target, *init_target;
	long pos;

	init_target = target;
	stop_target = &(target[ntargets - V512_BLOCKSIZE]);

	if (_UNLIKELY(target >= stop_target)) {
		if ((pos = Intersection_find_scalar(goal, target, ntargets)) <= ntargets
				&& target[pos] == goal)
			*foundp = 1;
		else
			*foundp = 0;
		return pos;
	}

	while (target[V512_BLOCKSIZE - 1] < goal) {
		target += V512_BLOCKSIZE;
		if (target >= stop_target) {
			if ((pos = Intersection_find_scalar(goal, target, ntargets))
					<= ntargets && target[pos] == goal)
				*foundp = 1;
			else
				*foundp = 0;

			return (target - init_target) + pos;
		}
	}

	// now search the block [target, target+block_size-1]
	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0, Q0, Q1, Q2, Q3;

	if (target[V512_BLOCKSIZE / 2 - 1] >= goal) {
		if (target[V512_BLOCKSIZE / 4 - 1] >= goal) {
			if (target[V512_BLOCKSIZE / 8 - 1] >= goal) {
				if (target[V512_BLOCKSIZE / 16 - 1] >= goal) {
				} else {
					target += 32;
				}
			} else {
				if (target[V512_BLOCKSIZE / 8 + 32 - 1] >= goal) {
					target += 64;
				} else {
					target += 96;
				}
			}
		} else {
			if (target[V512_BLOCKSIZE / 4 + 64 - 1] >= goal) {
				if (target[V512_BLOCKSIZE / 8 + 96 - 1] >= goal) {
					target += 128;
				} else {
					target += 160;
				}
			} else {
				if (target[V512_BLOCKSIZE / 4 + 96 - 1] >= goal) {
					target += 192;
				} else {
					target += 224;
				}
			}
		}
	} else {
		if (target[V512_BLOCKSIZE / 4 + 256 - 1] >= goal) {
			if (target[V512_BLOCKSIZE / 8 + 256 - 1] >= goal) {
				if (target[V512_BLOCKSIZE / 16 + 256 - 1] >= goal) {
					target += 256;
				} else {
					target += 288;
				}
			} else {
				if (target[V512_BLOCKSIZE / 8 + 288 - 1] >= goal) {
					target += 320;
				} else {
					target += 352;
				}
			}
		} else {
			if (target[V512_BLOCKSIZE / 4 + 320 - 1] >= goal) {
				if (target[V512_BLOCKSIZE / 8 + 352 - 1] >= goal) {
					target += 384;
				} else {
					target += 416;
				}
			} else {
				if (target[V512_BLOCKSIZE / 4 + 352 - 1] >= goal) {
					target += 448;
				} else {
					target += 480;
				}
			}
		}
	}

	Q0 = _mm_or_si128(
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 0), Match),
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 1), Match));
	Q1 = _mm_or_si128(
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 2), Match),
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 3), Match));
	Q2 = _mm_or_si128(
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 4), Match),
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 5), Match));
	Q3 = _mm_or_si128(
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 6), Match),
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 7), Match));

	F0 = _mm_or_si128(_mm_or_si128(Q0, Q1), _mm_or_si128(Q2, Q3));
	if (_mm_testz_si128(F0, F0))
		*foundp = 0;
	else
		*foundp = 1;
	return (target - init_target);
}

long simdlinear_v512_rough_plow(int*foundp, UINT4 goal, const UINT4 *target,
		long ntargets) {
	const UINT4 *stop_target, *init_target;
	long pos;

	init_target = target;
	stop_target = &(target[ntargets - V512_BLOCKSIZE]);

	if (_UNLIKELY(target >= stop_target)) {
		if ((pos = Intersection_find_scalar(goal, target, ntargets)) <= ntargets
				&& target[pos] == goal)
			*foundp = 1;
		else
			*foundp = 0;
		return pos;
	}

	while (target[V512_BLOCKSIZE - 1] < goal) {
		target += V512_BLOCKSIZE;
		if (target >= stop_target) {
			if ((pos = Intersection_find_scalar(goal, target, ntargets))
					<= ntargets && target[pos] == goal)
				*foundp = 1;
			else
				*foundp = 0;

			return (target - init_target) + pos;
		}
	}

	// now search the block [target, target+block_size-1]
	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0 = _mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 0),
			Match);
	UINT4 i = 0;
	while (_mm_testz_si128(F0, F0) && ++i < V512_BLOCKSIZE / SIMDWIDTH) {
		F0 = _mm_or_si128(F0,
				_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + i),
						Match));
	}
	if (_mm_testz_si128(F0, F0)) {
		*foundp = 0;
		i--;
	} else
		*foundp = 1;
	return (target + i * SIMDWIDTH - init_target);
}

long simdgallop_v4_exact(UINT4 goal, const UINT4 *target, long ntargets) {
	const UINT4 *end_target, *stop_target, *init_target;

	long low_offset = 0, mid_offset, high_offset = 1;

	init_target = target;
	stop_target = target + ntargets - V4_BLOCKSIZE;
	end_target = target + ntargets;

	if (target >= stop_target)
		return Intersection_find_scalar(goal, target, ntargets);

	if (target[V4_BLOCKSIZE - 1] < goal) {
		if (target + V4_BLOCKSIZE > stop_target)
			return V4_BLOCKSIZE
					+ Intersection_find_scalar(goal, target + V4_BLOCKSIZE,
							ntargets - V4_BLOCKSIZE);

		/* Galloping search */
		while (target[V4_BLOCKSIZE * high_offset + V4_BLOCKSIZE - 1] < goal) {
			if (target + (high_offset << 1) * V4_BLOCKSIZE <= stop_target) {
				low_offset = high_offset;
				high_offset <<= 1;
			} else if (target + V4_BLOCKSIZE * (high_offset + 1)
					<= stop_target) {
				high_offset = (stop_target - target) / V4_BLOCKSIZE;
				if (target[V4_BLOCKSIZE * high_offset + V4_BLOCKSIZE - 1]
						< goal) {
					target += V4_BLOCKSIZE * high_offset;
					return V4_BLOCKSIZE * high_offset
							+ Intersection_find_scalar(goal, target,
									(end_target - target));
				} else
					break;
			} else {
				target += V4_BLOCKSIZE * high_offset;
				return V4_BLOCKSIZE * high_offset
						+ Intersection_find_scalar(goal, target,
								(end_target - target));
			}
		}			// while-loop
		while (low_offset < high_offset) {
			mid_offset = (low_offset + high_offset) / 2;
			if (target[V4_BLOCKSIZE * mid_offset + V4_BLOCKSIZE - 1] < goal) {
				low_offset = mid_offset + 1;
			} else {
				high_offset = mid_offset;
			}
		}
		target += V4_BLOCKSIZE * high_offset;
	}
	// now search the block [target, target+block_size-1]
	__m128i Match = _mm_set1_epi32(goal);
	UINT4 pos = V4_POS[_mm_movemask_ps(
			_mm_castsi128_ps(
					_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 0),
							Match)))];

	return (target - init_target) + pos;
}

long simdgallop_v4_rough(int *foundp, UINT4 goal, const UINT4 *target,
		long ntargets) {
	const UINT4 *end_target, *stop_target, *init_target;

	long low_offset = 0, mid_offset, high_offset = 1;
	long pos;

	init_target = target;
	stop_target = target + ntargets - V4_BLOCKSIZE;
	end_target = target + ntargets;

	if (target >= stop_target) {
		if ((pos = Intersection_find_scalar(goal, target, ntargets)) <= ntargets
				&& target[pos] == goal)
			*foundp = 1;
		else
			*foundp = 0;
		return pos;
	}

	if (target[V4_BLOCKSIZE - 1] < goal) {
		if (target + V4_BLOCKSIZE > stop_target) {
			pos = V4_BLOCKSIZE
					+ Intersection_find_scalar(goal, target + V4_BLOCKSIZE,
							ntargets - V4_BLOCKSIZE);
			if (pos <= ntargets && target[pos] == goal)
				*foundp = 1;
			else
				*foundp = 0;
			return pos;
		}
		/* Galloping search */
		while (target[V4_BLOCKSIZE * high_offset + V4_BLOCKSIZE - 1] < goal) {
			if (target + (high_offset << 1) * V4_BLOCKSIZE <= stop_target) {
				low_offset = high_offset;
				high_offset <<= 1;
			} else if (target + V4_BLOCKSIZE * (high_offset + 1)
					<= stop_target) {
				high_offset = (stop_target - target) / V4_BLOCKSIZE;
				if (target[V4_BLOCKSIZE * high_offset + V4_BLOCKSIZE - 1]
						< goal) {
					target += V4_BLOCKSIZE * high_offset;
					pos = Intersection_find_scalar(goal, target,
							(end_target - target));
					if (pos + V4_BLOCKSIZE * high_offset <= ntargets
							&& target[pos] == goal)
						*foundp = 1;
					else
						*foundp = 0;
					return pos + V4_BLOCKSIZE * high_offset;
				} else
					break;
			} else {
				target += V4_BLOCKSIZE * high_offset;

				pos = Intersection_find_scalar(goal, target,
						(end_target - target));
				if (pos + V4_BLOCKSIZE * high_offset <= ntargets
						&& target[pos] == goal)
					*foundp = 1;
				else
					*foundp = 0;
				return pos + V4_BLOCKSIZE * high_offset;
			}
		}			// while-loop
		while (low_offset < high_offset) {
			mid_offset = (low_offset + high_offset) / 2;
			if (target[V4_BLOCKSIZE * mid_offset + V4_BLOCKSIZE - 1] < goal) {
				low_offset = mid_offset + 1;
			} else {
				high_offset = mid_offset;
			}
		}
		target += V4_BLOCKSIZE * high_offset;
	}
	// now search the block [target, target+block_size-1]
	__m128i Match = _mm_set1_epi32(goal);

	__m128i F0 = _mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 0),
			Match);
	if (_mm_testz_si128(F0, F0))
		*foundp = 0;
	else
		*foundp = 1;
	return (target - init_target);
}

long simdgallop_v8_exact(UINT4 goal, const UINT4 *target, long ntargets) {
	const UINT4 *end_target, *stop_target, *init_target;

	long low_offset = 0, mid_offset, high_offset = 1;

	init_target = target;
	stop_target = target + ntargets - V8_BLOCKSIZE;
	end_target = target + ntargets;

	if (target >= stop_target)
		return Intersection_find_scalar(goal, target, ntargets);

	if (target[V8_BLOCKSIZE - 1] < goal) {
		if (target + V8_BLOCKSIZE > stop_target)
			return V8_BLOCKSIZE
					+ Intersection_find_scalar(goal, target + V8_BLOCKSIZE,
							ntargets - V8_BLOCKSIZE);

		/* Galloping search */
		while (target[V8_BLOCKSIZE * high_offset + V8_BLOCKSIZE - 1] < goal) {
			if (target + (high_offset << 1) * V8_BLOCKSIZE <= stop_target) {
				low_offset = high_offset;
				high_offset <<= 1;
			} else if (target + V8_BLOCKSIZE * (high_offset + 1)
					<= stop_target) {
				high_offset = (stop_target - target) / V8_BLOCKSIZE;
				if (target[V8_BLOCKSIZE * high_offset + V8_BLOCKSIZE - 1]
						< goal) {
					target += V8_BLOCKSIZE * high_offset;
					return V8_BLOCKSIZE * high_offset
							+ Intersection_find_scalar(goal, target,
									(end_target - target));
				} else
					break;
			} else {
				target += V8_BLOCKSIZE * high_offset;
				return V8_BLOCKSIZE * high_offset
						+ Intersection_find_scalar(goal, target,
								(end_target - target));
			}
		}			// while-loop
		while (low_offset < high_offset) {
			mid_offset = (low_offset + high_offset) / 2;
			if (target[V8_BLOCKSIZE * mid_offset + V8_BLOCKSIZE - 1] < goal) {
				low_offset = mid_offset + 1;
			} else {
				high_offset = mid_offset;
			}
		}
		target += V8_BLOCKSIZE * high_offset;
	}
	// now search the block [target, target+block_size-1]
	__m128i Match = _mm_set1_epi32(goal);
	long pos;

	__m128i F0 = _mm_or_si128(
			_mm_srli_epi32(
					_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 0),
							Match), 32 - 1),
			_mm_srli_epi32(
					_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 1),
							Match), 28 - 1));

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

long simdgallop_v8_rough(int *foundp, UINT4 goal, const UINT4 *target,
		long ntargets) {
	const UINT4 *end_target, *stop_target, *init_target;

	long low_offset = 0, mid_offset, high_offset = 1;
	long pos;

	init_target = target;
	stop_target = target + ntargets - V8_BLOCKSIZE;
	end_target = target + ntargets;

	if (target >= stop_target) {
		if ((pos = Intersection_find_scalar(goal, target, ntargets)) <= ntargets
				&& target[pos] == goal)
			*foundp = 1;
		else
			*foundp = 0;
		return pos;
	}

	if (target[V8_BLOCKSIZE - 1] < goal) {
		if (target + V8_BLOCKSIZE > stop_target) {
			pos = V8_BLOCKSIZE
					+ Intersection_find_scalar(goal, target + V8_BLOCKSIZE,
							ntargets - V8_BLOCKSIZE);
			if (pos <= ntargets && target[pos] == goal)
				*foundp = 1;
			else
				*foundp = 0;
			return pos;
		}
		/* Galloping search */
		while (target[V8_BLOCKSIZE * high_offset + V8_BLOCKSIZE - 1] < goal) {
			if (target + (high_offset << 1) * V8_BLOCKSIZE <= stop_target) {
				low_offset = high_offset;
				high_offset <<= 1;
			} else if (target + V8_BLOCKSIZE * (high_offset + 1)
					<= stop_target) {
				high_offset = (stop_target - target) / V8_BLOCKSIZE;
				if (target[V8_BLOCKSIZE * high_offset + V8_BLOCKSIZE - 1]
						< goal) {
					target += V8_BLOCKSIZE * high_offset;
					pos = Intersection_find_scalar(goal, target,
							(end_target - target));
					if (pos + V8_BLOCKSIZE * high_offset <= ntargets
							&& target[pos] == goal)
						*foundp = 1;
					else
						*foundp = 0;
					return pos + V8_BLOCKSIZE * high_offset;
				} else
					break;
			} else {
				target += V8_BLOCKSIZE * high_offset;

				pos = Intersection_find_scalar(goal, target,
						(end_target - target));
				if (pos + V8_BLOCKSIZE * high_offset <= ntargets
						&& target[pos] == goal)
					*foundp = 1;
				else
					*foundp = 0;
				return pos + V8_BLOCKSIZE * high_offset;
			}
		}			// while-loop
		while (low_offset < high_offset) {
			mid_offset = (low_offset + high_offset) / 2;
			if (target[V8_BLOCKSIZE * mid_offset + V8_BLOCKSIZE - 1] < goal) {
				low_offset = mid_offset + 1;
			} else {
				high_offset = mid_offset;
			}
		}
		target += V8_BLOCKSIZE * high_offset;
	}
	// now search the block [target, target+block_size-1]
	__m128i Match = _mm_set1_epi32(goal);

	__m128i F0 = _mm_or_si128(
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 0), Match),
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 1), Match));
	if (_mm_testz_si128(F0, F0))
		*foundp = 0;
	else
		*foundp = 1;
	return (target - init_target);
}

long simdgallop_v8_rough_plow(int *foundp, UINT4 goal, const UINT4 *target,
		long ntargets) {
	const UINT4 *end_target, *stop_target, *init_target;

	long low_offset = 0, mid_offset, high_offset = 1;
	long pos;

	init_target = target;
	stop_target = target + ntargets - V8_BLOCKSIZE;
	end_target = target + ntargets;

	if (target >= stop_target) {
		if ((pos = Intersection_find_scalar(goal, target, ntargets)) <= ntargets
				&& target[pos] == goal)
			*foundp = 1;
		else
			*foundp = 0;
		return pos;
	}

	if (target[V8_BLOCKSIZE - 1] < goal) {
		if (target + V8_BLOCKSIZE > stop_target) {
			pos = V8_BLOCKSIZE
					+ Intersection_find_scalar(goal, target + V8_BLOCKSIZE,
							ntargets - V8_BLOCKSIZE);
			if (pos <= ntargets && target[pos] == goal)
				*foundp = 1;
			else
				*foundp = 0;
			return pos;
		}
		/* Galloping search */
		while (target[V8_BLOCKSIZE * high_offset + V8_BLOCKSIZE - 1] < goal) {
			if (target + (high_offset << 1) * V8_BLOCKSIZE <= stop_target) {
				low_offset = high_offset;
				high_offset <<= 1;
			} else if (target + V8_BLOCKSIZE * (high_offset + 1)
					<= stop_target) {
				high_offset = (stop_target - target) / V8_BLOCKSIZE;
				if (target[V8_BLOCKSIZE * high_offset + V8_BLOCKSIZE - 1]
						< goal) {
					target += V8_BLOCKSIZE * high_offset;
					pos = Intersection_find_scalar(goal, target,
							(end_target - target));
					if (pos + V8_BLOCKSIZE * high_offset <= ntargets
							&& target[pos] == goal)
						*foundp = 1;
					else
						*foundp = 0;
					return pos + V8_BLOCKSIZE * high_offset;
				} else
					break;
			} else {
				target += V8_BLOCKSIZE * high_offset;

				pos = Intersection_find_scalar(goal, target,
						(end_target - target));
				if (pos + V8_BLOCKSIZE * high_offset <= ntargets
						&& target[pos] == goal)
					*foundp = 1;
				else
					*foundp = 0;
				return pos + V8_BLOCKSIZE * high_offset;
			}
		}			// while-loop
		while (low_offset < high_offset) {
			mid_offset = (low_offset + high_offset) / 2;
			if (target[V8_BLOCKSIZE * mid_offset + V8_BLOCKSIZE - 1] < goal) {
				low_offset = mid_offset + 1;
			} else {
				high_offset = mid_offset;
			}
		}
		target += V8_BLOCKSIZE * high_offset;
	}
	// now search the block [target, target+block_size-1]
	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0 = _mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 0),
			Match);
	UINT4 i = 0;
	while (_mm_testz_si128(F0, F0) && ++i < V8_BLOCKSIZE / SIMDWIDTH) {
		F0 = _mm_or_si128(F0,
				_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + i),
						Match));
	}
	if (_mm_testz_si128(F0, F0)) {
		*foundp = 0;
		i--;
	} else
		*foundp = 1;
	return (target + i * SIMDWIDTH - init_target);
}

long simdgallop_v16_exact(UINT4 goal, const UINT4 *target, long ntargets) {
	const UINT4 *end_target, *stop_target, *init_target;

	long low_offset = 0, mid_offset, high_offset = 1;

	init_target = target;
	stop_target = target + ntargets - V16_BLOCKSIZE;
	end_target = target + ntargets;

	if (target >= stop_target)
		return Intersection_find_scalar(goal, target, ntargets);

	if (target[V16_BLOCKSIZE - 1] < goal) {
		if (target + V16_BLOCKSIZE > stop_target)
			return V16_BLOCKSIZE
					+ Intersection_find_scalar(goal, target + V16_BLOCKSIZE,
							ntargets - V16_BLOCKSIZE);

		/* Galloping search */
		while (target[V16_BLOCKSIZE * high_offset + V16_BLOCKSIZE - 1] < goal) {
			if (target + (high_offset << 1) * V16_BLOCKSIZE <= stop_target) {
				low_offset = high_offset;
				high_offset <<= 1;
			} else if (target + V16_BLOCKSIZE * (high_offset + 1)
					<= stop_target) {
				high_offset = (stop_target - target) / V16_BLOCKSIZE;
				if (target[V16_BLOCKSIZE * high_offset + V16_BLOCKSIZE - 1]
						< goal) {
					target += V16_BLOCKSIZE * high_offset;
					return V16_BLOCKSIZE * high_offset
							+ Intersection_find_scalar(goal, target,
									(end_target - target));
				} else
					break;
			} else {
				target += V16_BLOCKSIZE * high_offset;
				return V16_BLOCKSIZE * high_offset
						+ Intersection_find_scalar(goal, target,
								(end_target - target));
			}
		}			// while-loop
		while (low_offset < high_offset) {
			mid_offset = (low_offset + high_offset) / 2;
			if (target[V16_BLOCKSIZE * mid_offset + V16_BLOCKSIZE - 1] < goal) {
				low_offset = mid_offset + 1;
			} else {
				high_offset = mid_offset;
			}
		}
		target += V16_BLOCKSIZE * high_offset;
	}
	// now search the block [target, target+block_size-1]
	__m128i Match = _mm_set1_epi32(goal);
	long pos;

	__m128i F0 = _mm_or_si128(
			_mm_or_si128(
					_mm_srli_epi32(
							_mm_cmpeq_epi32(
									_mm_lddqu_si128((__m128i *) target + 0),
									Match), 32 - 1),
					_mm_srli_epi32(
							_mm_cmpeq_epi32(
									_mm_lddqu_si128((__m128i *) target + 1),
									Match), 28 - 1)),
			_mm_or_si128(
					_mm_srli_epi32(
							_mm_cmpeq_epi32(
									_mm_lddqu_si128((__m128i *) target + 2),
									Match), 24 - 1),
					_mm_srli_epi32(
							_mm_cmpeq_epi32(
									_mm_lddqu_si128((__m128i *) target + 3),
									Match), 20 - 1)));

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

long simdgallop_v16_rough(int *foundp, UINT4 goal, const UINT4 *target,
		long ntargets) {
	const UINT4 *end_target, *stop_target, *init_target;

	long low_offset = 0, mid_offset, high_offset = 1;
	long pos;

	init_target = target;
	stop_target = target + ntargets - V16_BLOCKSIZE;
	end_target = target + ntargets;

	if (target >= stop_target) {
		if ((pos = Intersection_find_scalar(goal, target, ntargets)) <= ntargets
				&& target[pos] == goal)
			*foundp = 1;
		else
			*foundp = 0;
		return pos;
	}

	if (target[V16_BLOCKSIZE - 1] < goal) {
		if (target + V16_BLOCKSIZE > stop_target) {
			pos = V16_BLOCKSIZE
					+ Intersection_find_scalar(goal, target + V16_BLOCKSIZE,
							ntargets - V16_BLOCKSIZE);
			if (pos <= ntargets && target[pos] == goal)
				*foundp = 1;
			else
				*foundp = 0;
			return pos;
		}
		/* Galloping search */
		while (target[V16_BLOCKSIZE * high_offset + V16_BLOCKSIZE - 1] < goal) {
			if (target + (high_offset << 1) * V16_BLOCKSIZE <= stop_target) {
				low_offset = high_offset;
				high_offset <<= 1;
			} else if (target + V16_BLOCKSIZE * (high_offset + 1)
					<= stop_target) {
				high_offset = (stop_target - target) / V16_BLOCKSIZE;
				if (target[V16_BLOCKSIZE * high_offset + V16_BLOCKSIZE - 1]
						< goal) {
					target += V16_BLOCKSIZE * high_offset;
					pos = Intersection_find_scalar(goal, target,
							(end_target - target));
					if (pos + V16_BLOCKSIZE * high_offset <= ntargets
							&& target[pos] == goal)
						*foundp = 1;
					else
						*foundp = 0;
					return pos + V16_BLOCKSIZE * high_offset;
				} else
					break;
			} else {
				target += V16_BLOCKSIZE * high_offset;

				pos = Intersection_find_scalar(goal, target,
						(end_target - target));
				if (pos + V16_BLOCKSIZE * high_offset <= ntargets
						&& target[pos] == goal)
					*foundp = 1;
				else
					*foundp = 0;
				return pos + V16_BLOCKSIZE * high_offset;
			}
		}			// while-loop
		while (low_offset < high_offset) {
			mid_offset = (low_offset + high_offset) / 2;
			if (target[V16_BLOCKSIZE * mid_offset + V16_BLOCKSIZE - 1] < goal) {
				low_offset = mid_offset + 1;
			} else {
				high_offset = mid_offset;
			}
		}
		target += V16_BLOCKSIZE * high_offset;
	}
	// now search the block [target, target+block_size-1]
	__m128i Match = _mm_set1_epi32(goal);

	__m128i F0 = _mm_or_si128(
			_mm_or_si128(
					_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 0),
							Match),
					_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 1),
							Match)),
			_mm_or_si128(
					_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 2),
							Match),
					_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 3),
							Match)));
	if (_mm_testz_si128(F0, F0))
		*foundp = 0;
	else
		*foundp = 1;
	return (target - init_target);
}

long simdgallop_v16_rough_plow(int *foundp, UINT4 goal, const UINT4 *target,
		long ntargets) {
	const UINT4 *end_target, *stop_target, *init_target;

	long low_offset = 0, mid_offset, high_offset = 1;
	long pos;

	init_target = target;
	stop_target = target + ntargets - V16_BLOCKSIZE;
	end_target = target + ntargets;

	if (target >= stop_target) {
		if ((pos = Intersection_find_scalar(goal, target, ntargets)) <= ntargets
				&& target[pos] == goal)
			*foundp = 1;
		else
			*foundp = 0;
		return pos;
	}

	if (target[V16_BLOCKSIZE - 1] < goal) {
		if (target + V16_BLOCKSIZE > stop_target) {
			pos = V16_BLOCKSIZE
					+ Intersection_find_scalar(goal, target + V16_BLOCKSIZE,
							ntargets - V16_BLOCKSIZE);
			if (pos <= ntargets && target[pos] == goal)
				*foundp = 1;
			else
				*foundp = 0;
			return pos;
		}
		/* Galloping search */
		while (target[V16_BLOCKSIZE * high_offset + V16_BLOCKSIZE - 1] < goal) {
			if (target + (high_offset << 1) * V16_BLOCKSIZE <= stop_target) {
				low_offset = high_offset;
				high_offset <<= 1;
			} else if (target + V16_BLOCKSIZE * (high_offset + 1)
					<= stop_target) {
				high_offset = (stop_target - target) / V16_BLOCKSIZE;
				if (target[V16_BLOCKSIZE * high_offset + V16_BLOCKSIZE - 1]
						< goal) {
					target += V16_BLOCKSIZE * high_offset;
					pos = Intersection_find_scalar(goal, target,
							(end_target - target));
					if (pos + V16_BLOCKSIZE * high_offset <= ntargets
							&& target[pos] == goal)
						*foundp = 1;
					else
						*foundp = 0;
					return pos + V16_BLOCKSIZE * high_offset;
				} else
					break;
			} else {
				target += V16_BLOCKSIZE * high_offset;

				pos = Intersection_find_scalar(goal, target,
						(end_target - target));
				if (pos + V16_BLOCKSIZE * high_offset <= ntargets
						&& target[pos] == goal)
					*foundp = 1;
				else
					*foundp = 0;
				return pos + V16_BLOCKSIZE * high_offset;
			}
		}			// while-loop
		while (low_offset < high_offset) {
			mid_offset = (low_offset + high_offset) / 2;
			if (target[V16_BLOCKSIZE * mid_offset + V16_BLOCKSIZE - 1] < goal) {
				low_offset = mid_offset + 1;
			} else {
				high_offset = mid_offset;
			}
		}
		target += V16_BLOCKSIZE * high_offset;
	}
	// now search the block [target, target+block_size-1]
	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0 = _mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 0),
			Match);
	UINT4 i = 0;
	while (_mm_testz_si128(F0, F0) && ++i < V16_BLOCKSIZE / SIMDWIDTH) {
		F0 = _mm_or_si128(F0,
				_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + i),
						Match));
	}
	if (_mm_testz_si128(F0, F0)) {
		*foundp = 0;
		i--;
	} else
		*foundp = 1;
	return (target + i * SIMDWIDTH - init_target);
}

long simdgallop_v32_exact(UINT4 goal, const UINT4 *target, long ntargets) {
	const UINT4 *end_target, *stop_target, *init_target;

	long low_offset = 0, mid_offset, high_offset = 1;

	init_target = target;
	stop_target = target + ntargets - V32_BLOCKSIZE;
	end_target = target + ntargets;

	if (target >= stop_target)
		return Intersection_find_scalar(goal, target, ntargets);

	if (target[V32_BLOCKSIZE - 1] < goal) {
		if (target + V32_BLOCKSIZE > stop_target)
			return V32_BLOCKSIZE
					+ Intersection_find_scalar(goal, target + V32_BLOCKSIZE,
							ntargets - V32_BLOCKSIZE);

		/* Galloping search */
		while (target[V32_BLOCKSIZE * high_offset + V32_BLOCKSIZE - 1] < goal) {
			if (target + (high_offset << 1) * V32_BLOCKSIZE <= stop_target) {
				low_offset = high_offset;
				high_offset <<= 1;
			} else if (target + V32_BLOCKSIZE * (high_offset + 1)
					<= stop_target) {
				high_offset = (stop_target - target) / V32_BLOCKSIZE;
				if (target[V32_BLOCKSIZE * high_offset + V32_BLOCKSIZE - 1]
						< goal) {
					target += V32_BLOCKSIZE * high_offset;
					return V32_BLOCKSIZE * high_offset
							+ Intersection_find_scalar(goal, target,
									(end_target - target));
				} else
					break;
			} else {
				target += V32_BLOCKSIZE * high_offset;
				return V32_BLOCKSIZE * high_offset
						+ Intersection_find_scalar(goal, target,
								(end_target - target));
			}
		}			// while-loop
		while (low_offset < high_offset) {
			mid_offset = (low_offset + high_offset) / 2;
			if (target[V32_BLOCKSIZE * mid_offset + V32_BLOCKSIZE - 1] < goal) {
				low_offset = mid_offset + 1;
			} else {
				high_offset = mid_offset;
			}
		}
		target += V32_BLOCKSIZE * high_offset;
	}
	// now search the block [target, target+block_size-1]
	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0, Q0, Q1, Q2, Q3;
	long pos;

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

long simdgallop_v32_rough(int *foundp, UINT4 goal, const UINT4 *target,
		long ntargets) {
	const UINT4 *end_target, *stop_target, *init_target;

	long low_offset = 0, mid_offset, high_offset = 1;
	long pos;

	init_target = target;
	stop_target = target + ntargets - V32_BLOCKSIZE;
	end_target = target + ntargets;

	if (target >= stop_target) {
		if ((pos = Intersection_find_scalar(goal, target, ntargets)) <= ntargets
				&& target[pos] == goal)
			*foundp = 1;
		else
			*foundp = 0;
		return pos;
	}

	if (target[V32_BLOCKSIZE - 1] < goal) {
		if (target + V32_BLOCKSIZE > stop_target) {
			pos = V32_BLOCKSIZE
					+ Intersection_find_scalar(goal, target + V32_BLOCKSIZE,
							ntargets - V32_BLOCKSIZE);
			if (pos <= ntargets && target[pos] == goal)
				*foundp = 1;
			else
				*foundp = 0;
			return pos;
		}
		/* Galloping search */
		while (target[V32_BLOCKSIZE * high_offset + V32_BLOCKSIZE - 1] < goal) {
			if (target + (high_offset << 1) * V32_BLOCKSIZE <= stop_target) {
				low_offset = high_offset;
				high_offset <<= 1;
			} else if (target + V32_BLOCKSIZE * (high_offset + 1)
					<= stop_target) {
				high_offset = (stop_target - target) / V32_BLOCKSIZE;
				if (target[V32_BLOCKSIZE * high_offset + V32_BLOCKSIZE - 1]
						< goal) {
					target += V32_BLOCKSIZE * high_offset;
					pos = Intersection_find_scalar(goal, target,
							(end_target - target));
					if (pos + V32_BLOCKSIZE * high_offset <= ntargets
							&& target[pos] == goal)
						*foundp = 1;
					else
						*foundp = 0;
					return pos + V32_BLOCKSIZE * high_offset;
				} else
					break;
			} else {
				target += V32_BLOCKSIZE * high_offset;

				pos = Intersection_find_scalar(goal, target,
						(end_target - target));
				if (pos + V32_BLOCKSIZE * high_offset <= ntargets
						&& target[pos] == goal)
					*foundp = 1;
				else
					*foundp = 0;
				return pos + V32_BLOCKSIZE * high_offset;
			}
		}			// while-loop
		while (low_offset < high_offset) {
			mid_offset = (low_offset + high_offset) / 2;
			if (target[V32_BLOCKSIZE * mid_offset + V32_BLOCKSIZE - 1] < goal) {
				low_offset = mid_offset + 1;
			} else {
				high_offset = mid_offset;
			}
		}
		target += V32_BLOCKSIZE * high_offset;
	}
	// now search the block [target, target+block_size-1]
	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0, Q0, Q1, Q2, Q3;

	Q0 = _mm_or_si128(
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 0), Match),
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 1), Match));
	Q1 = _mm_or_si128(
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 2), Match),
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 3), Match));
	Q2 = _mm_or_si128(
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 4), Match),
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 5), Match));
	Q3 = _mm_or_si128(
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 6), Match),
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 7), Match));

	F0 = _mm_or_si128(_mm_or_si128(Q0, Q1), _mm_or_si128(Q2, Q3));
	if (_mm_testz_si128(F0, F0))
		*foundp = 0;
	else
		*foundp = 1;
	return (target - init_target);
}

long simdgallop_v32_rough_plow(int *foundp, UINT4 goal, const UINT4 *target,
		long ntargets) {
	const UINT4 *end_target, *stop_target, *init_target;

	long low_offset = 0, mid_offset, high_offset = 1;
	long pos;

	init_target = target;
	stop_target = target + ntargets - V32_BLOCKSIZE;
	end_target = target + ntargets;

	if (target >= stop_target) {
		if ((pos = Intersection_find_scalar(goal, target, ntargets)) <= ntargets
				&& target[pos] == goal)
			*foundp = 1;
		else
			*foundp = 0;
		return pos;
	}

	if (target[V32_BLOCKSIZE - 1] < goal) {
		if (target + V32_BLOCKSIZE > stop_target) {
			pos = V32_BLOCKSIZE
					+ Intersection_find_scalar(goal, target + V32_BLOCKSIZE,
							ntargets - V32_BLOCKSIZE);
			if (pos <= ntargets && target[pos] == goal)
				*foundp = 1;
			else
				*foundp = 0;
			return pos;
		}
		/* Galloping search */
		while (target[V32_BLOCKSIZE * high_offset + V32_BLOCKSIZE - 1] < goal) {
			if (target + (high_offset << 1) * V32_BLOCKSIZE <= stop_target) {
				low_offset = high_offset;
				high_offset <<= 1;
			} else if (target + V32_BLOCKSIZE * (high_offset + 1)
					<= stop_target) {
				high_offset = (stop_target - target) / V32_BLOCKSIZE;
				if (target[V32_BLOCKSIZE * high_offset + V32_BLOCKSIZE - 1]
						< goal) {
					target += V32_BLOCKSIZE * high_offset;
					pos = Intersection_find_scalar(goal, target,
							(end_target - target));
					if (pos + V32_BLOCKSIZE * high_offset <= ntargets
							&& target[pos] == goal)
						*foundp = 1;
					else
						*foundp = 0;
					return pos + V32_BLOCKSIZE * high_offset;
				} else
					break;
			} else {
				target += V32_BLOCKSIZE * high_offset;

				pos = Intersection_find_scalar(goal, target,
						(end_target - target));
				if (pos + V32_BLOCKSIZE * high_offset <= ntargets
						&& target[pos] == goal)
					*foundp = 1;
				else
					*foundp = 0;
				return pos + V32_BLOCKSIZE * high_offset;
			}
		}			// while-loop
		while (low_offset < high_offset) {
			mid_offset = (low_offset + high_offset) / 2;
			if (target[V32_BLOCKSIZE * mid_offset + V32_BLOCKSIZE - 1] < goal) {
				low_offset = mid_offset + 1;
			} else {
				high_offset = mid_offset;
			}
		}
		target += V32_BLOCKSIZE * high_offset;
	}
	// now search the block [target, target+block_size-1]
	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0 = _mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 0),
			Match);
	UINT4 i = 0;
	while (_mm_testz_si128(F0, F0) && ++i < V32_BLOCKSIZE / SIMDWIDTH) {
		F0 = _mm_or_si128(F0,
				_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + i),
						Match));
	}
	if (_mm_testz_si128(F0, F0)) {
		*foundp = 0;
		i--;
	} else
		*foundp = 1;
	return (target + i * SIMDWIDTH - init_target);
}

long simdgallop_v64_exact(UINT4 goal, const UINT4 *target, long ntargets) {
	const UINT4 *end_target, *stop_target, *init_target;

	long low_offset = 0, mid_offset, high_offset = 1;

	init_target = target;
	stop_target = target + ntargets - V64_BLOCKSIZE;
	end_target = target + ntargets;

	if (target >= stop_target)
		return Intersection_find_scalar(goal, target, ntargets);

	if (target[V64_BLOCKSIZE - 1] < goal) {
		if (target + V64_BLOCKSIZE > stop_target)
			return V64_BLOCKSIZE
					+ Intersection_find_scalar(goal, target + V64_BLOCKSIZE,
							ntargets - V64_BLOCKSIZE);

		/* Galloping search */
		while (target[V64_BLOCKSIZE * high_offset + V64_BLOCKSIZE - 1] < goal) {
			if (target + (high_offset << 1) * V64_BLOCKSIZE <= stop_target) {
				low_offset = high_offset;
				high_offset <<= 1;
			} else if (target + V64_BLOCKSIZE * (high_offset + 1)
					<= stop_target) {
				high_offset = (stop_target - target) / V64_BLOCKSIZE;
				if (target[V64_BLOCKSIZE * high_offset + V64_BLOCKSIZE - 1]
						< goal) {
					target += V64_BLOCKSIZE * high_offset;
					return V64_BLOCKSIZE * high_offset
							+ Intersection_find_scalar(goal, target,
									(end_target - target));
				} else
					break;
			} else {
				target += V64_BLOCKSIZE * high_offset;
				return V64_BLOCKSIZE * high_offset
						+ Intersection_find_scalar(goal, target,
								(end_target - target));
			}
		}			// while-loop
		while (low_offset < high_offset) {
			mid_offset = (low_offset + high_offset) / 2;
			if (target[V64_BLOCKSIZE * mid_offset + V64_BLOCKSIZE - 1] < goal) {
				low_offset = mid_offset + 1;
			} else {
				high_offset = mid_offset;
			}
		}
		target += V64_BLOCKSIZE * high_offset;
	}
	// now search the block [target, target+block_size-1]
	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0, Q0, Q1, Q2, Q3;
	long pos;

	if (target[V64_BLOCKSIZE / 2 - 1] < goal)
		target += 32;

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

long simdgallop_v64_rough(int *foundp, UINT4 goal, const UINT4 *target,
		long ntargets) {
	const UINT4 *end_target, *stop_target, *init_target;

	long low_offset = 0, mid_offset, high_offset = 1;
	long pos;

	init_target = target;
	stop_target = target + ntargets - V64_BLOCKSIZE;
	end_target = target + ntargets;

	if (target >= stop_target) {
		if ((pos = Intersection_find_scalar(goal, target, ntargets)) <= ntargets
				&& target[pos] == goal)
			*foundp = 1;
		else
			*foundp = 0;
		return pos;
	}

	if (target[V64_BLOCKSIZE - 1] < goal) {
		if (target + V64_BLOCKSIZE > stop_target) {
			pos = V64_BLOCKSIZE
					+ Intersection_find_scalar(goal, target + V64_BLOCKSIZE,
							ntargets - V64_BLOCKSIZE);
			if (pos <= ntargets && target[pos] == goal)
				*foundp = 1;
			else
				*foundp = 0;
			return pos;
		}
		/* Galloping search */
		while (target[V64_BLOCKSIZE * high_offset + V64_BLOCKSIZE - 1] < goal) {
			if (target + (high_offset << 1) * V64_BLOCKSIZE <= stop_target) {
				low_offset = high_offset;
				high_offset <<= 1;
			} else if (target + V64_BLOCKSIZE * (high_offset + 1)
					<= stop_target) {
				high_offset = (stop_target - target) / V64_BLOCKSIZE;
				if (target[V64_BLOCKSIZE * high_offset + V64_BLOCKSIZE - 1]
						< goal) {
					target += V64_BLOCKSIZE * high_offset;
					pos = Intersection_find_scalar(goal, target,
							(end_target - target));
					if (pos + V64_BLOCKSIZE * high_offset <= ntargets
							&& target[pos] == goal)
						*foundp = 1;
					else
						*foundp = 0;
					return pos + V64_BLOCKSIZE * high_offset;
				} else
					break;
			} else {
				target += V64_BLOCKSIZE * high_offset;

				pos = Intersection_find_scalar(goal, target,
						(end_target - target));
				if (pos + V64_BLOCKSIZE * high_offset <= ntargets
						&& target[pos] == goal)
					*foundp = 1;
				else
					*foundp = 0;
				return pos + V64_BLOCKSIZE * high_offset;
			}
		}			// while-loop
		while (low_offset < high_offset) {
			mid_offset = (low_offset + high_offset) / 2;
			if (target[V64_BLOCKSIZE * mid_offset + V64_BLOCKSIZE - 1] < goal) {
				low_offset = mid_offset + 1;
			} else {
				high_offset = mid_offset;
			}
		}
		target += V64_BLOCKSIZE * high_offset;
	}
	// now search the block [target, target+block_size-1]
	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0, Q0, Q1, Q2, Q3;

	if (target[V64_BLOCKSIZE / 2 - 1] < goal)
		target += 32;

	Q0 = _mm_or_si128(
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 0), Match),
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 1), Match));
	Q1 = _mm_or_si128(
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 2), Match),
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 3), Match));
	Q2 = _mm_or_si128(
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 4), Match),
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 5), Match));
	Q3 = _mm_or_si128(
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 6), Match),
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 7), Match));

	F0 = _mm_or_si128(_mm_or_si128(Q0, Q1), _mm_or_si128(Q2, Q3));
	if (_mm_testz_si128(F0, F0))
		*foundp = 0;
	else
		*foundp = 1;
	return (target - init_target);
}

long simdgallop_v64_rough_plow(int *foundp, UINT4 goal, const UINT4 *target,
		long ntargets) {
	const UINT4 *end_target, *stop_target, *init_target;

	long low_offset = 0, mid_offset, high_offset = 1;
	long pos;

	init_target = target;
	stop_target = target + ntargets - V64_BLOCKSIZE;
	end_target = target + ntargets;

	if (target >= stop_target) {
		if ((pos = Intersection_find_scalar(goal, target, ntargets)) <= ntargets
				&& target[pos] == goal)
			*foundp = 1;
		else
			*foundp = 0;
		return pos;
	}

	if (target[V64_BLOCKSIZE - 1] < goal) {
		if (target + V64_BLOCKSIZE > stop_target) {
			pos = V64_BLOCKSIZE
					+ Intersection_find_scalar(goal, target + V64_BLOCKSIZE,
							ntargets - V64_BLOCKSIZE);
			if (pos <= ntargets && target[pos] == goal)
				*foundp = 1;
			else
				*foundp = 0;
			return pos;
		}
		/* Galloping search */
		while (target[V64_BLOCKSIZE * high_offset + V64_BLOCKSIZE - 1] < goal) {
			if (target + (high_offset << 1) * V64_BLOCKSIZE <= stop_target) {
				low_offset = high_offset;
				high_offset <<= 1;
			} else if (target + V64_BLOCKSIZE * (high_offset + 1)
					<= stop_target) {
				high_offset = (stop_target - target) / V64_BLOCKSIZE;
				if (target[V64_BLOCKSIZE * high_offset + V64_BLOCKSIZE - 1]
						< goal) {
					target += V64_BLOCKSIZE * high_offset;
					pos = Intersection_find_scalar(goal, target,
							(end_target - target));
					if (pos + V64_BLOCKSIZE * high_offset <= ntargets
							&& target[pos] == goal)
						*foundp = 1;
					else
						*foundp = 0;
					return pos + V64_BLOCKSIZE * high_offset;
				} else
					break;
			} else {
				target += V64_BLOCKSIZE * high_offset;

				pos = Intersection_find_scalar(goal, target,
						(end_target - target));
				if (pos + V64_BLOCKSIZE * high_offset <= ntargets
						&& target[pos] == goal)
					*foundp = 1;
				else
					*foundp = 0;
				return pos + V64_BLOCKSIZE * high_offset;
			}
		}			// while-loop
		while (low_offset < high_offset) {
			mid_offset = (low_offset + high_offset) / 2;
			if (target[V64_BLOCKSIZE * mid_offset + V64_BLOCKSIZE - 1] < goal) {
				low_offset = mid_offset + 1;
			} else {
				high_offset = mid_offset;
			}
		}
		target += V64_BLOCKSIZE * high_offset;
	}
	// now search the block [target, target+block_size-1]
	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0 = _mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 0),
			Match);
	UINT4 i = 0;
	while (_mm_testz_si128(F0, F0) && ++i < V64_BLOCKSIZE / SIMDWIDTH) {
		F0 = _mm_or_si128(F0,
				_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + i),
						Match));
	}
	if (_mm_testz_si128(F0, F0)) {
		*foundp = 0;
		i--;
	} else
		*foundp = 1;
	return (target + i * SIMDWIDTH - init_target);
}

/* return the exact position of intersection */
long simdgallop_v128_exact(UINT4 goal, const UINT4 *target, long ntargets) {
	const UINT4 *end_target, *stop_target, *init_target;

	long low_offset = 0, mid_offset, high_offset = 1;

	init_target = target;
	stop_target = target + ntargets - V128_BLOCKSIZE;
	end_target = target + ntargets;

	if (target >= stop_target)
		return Intersection_find_scalar(goal, target, ntargets);

	if (target[V128_BLOCKSIZE - 1] < goal) {
		if (target + V128_BLOCKSIZE > stop_target)
			return V128_BLOCKSIZE
					+ Intersection_find_scalar(goal, target + V128_BLOCKSIZE,
							ntargets - V128_BLOCKSIZE);

		/* Galloping search */
		while (target[V128_BLOCKSIZE * high_offset + V128_BLOCKSIZE - 1] < goal) {
			if (target + (high_offset << 1) * V128_BLOCKSIZE <= stop_target) {
				low_offset = high_offset;
				high_offset <<= 1;
			} else if (target + V128_BLOCKSIZE * (high_offset + 1)
					<= stop_target) {
				high_offset = (stop_target - target) / V128_BLOCKSIZE;
				if (target[V128_BLOCKSIZE * high_offset + V128_BLOCKSIZE - 1]
						< goal) {
					target += V128_BLOCKSIZE * high_offset;
					return V128_BLOCKSIZE * high_offset
							+ Intersection_find_scalar(goal, target,
									(end_target - target));
				} else
					break;
			} else {
				target += V128_BLOCKSIZE * high_offset;
				return V128_BLOCKSIZE * high_offset
						+ Intersection_find_scalar(goal, target,
								(end_target - target));
			}
		}			// while-loop
		while (low_offset < high_offset) {
			mid_offset = (low_offset + high_offset) / 2;
			if (target[V128_BLOCKSIZE * mid_offset + V128_BLOCKSIZE - 1]
					< goal) {
				low_offset = mid_offset + 1;
			} else {
				high_offset = mid_offset;
			}
		}
		target += V128_BLOCKSIZE * high_offset;
	}
	// now search the block [target, target+block_size-1]
	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0, Q0, Q1, Q2, Q3;
	long pos;

	if (target[V128_BLOCKSIZE / 2 - 1] >= goal) {
		if (target[V128_BLOCKSIZE / 4 - 1] < goal)
			target += 32;
	} else {
		if (target[V128_BLOCKSIZE * 3 / 4 - 1] >= goal)
			target += 64;
		else
			target += 96;
	}
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

/* return the head-position of mathced block */
long simdgallop_v128_rough(int *foundp, UINT4 goal, const UINT4 *target,
		long ntargets) {
	const UINT4 *end_target, *stop_target, *init_target;

	long low_offset = 0, mid_offset, high_offset = 1;
	long pos;

	init_target = target;
	stop_target = target + ntargets - V128_BLOCKSIZE;
	end_target = target + ntargets;

	if (target >= stop_target) {
		if ((pos = Intersection_find_scalar(goal, target, ntargets)) <= ntargets
				&& target[pos] == goal)
			*foundp = 1;
		else
			*foundp = 0;
		return pos;
	}

	if (target[V128_BLOCKSIZE - 1] < goal) {
		if (target + V128_BLOCKSIZE > stop_target) {
			pos = V128_BLOCKSIZE
					+ Intersection_find_scalar(goal, target + V128_BLOCKSIZE,
							ntargets - V128_BLOCKSIZE);
			if (pos <= ntargets && target[pos] == goal)
				*foundp = 1;
			else
				*foundp = 0;
			return pos;
		}
		/* Galloping search */
		while (target[V128_BLOCKSIZE * high_offset + V128_BLOCKSIZE - 1] < goal) {
			if (target + (high_offset << 1) * V128_BLOCKSIZE <= stop_target) {
				low_offset = high_offset;
				high_offset <<= 1;
			} else if (target + V128_BLOCKSIZE * (high_offset + 1)
					<= stop_target) {
				high_offset = (stop_target - target) / V128_BLOCKSIZE;
				if (target[V128_BLOCKSIZE * high_offset + V128_BLOCKSIZE - 1]
						< goal) {
					target += V128_BLOCKSIZE * high_offset;
					pos = Intersection_find_scalar(goal, target,
							(end_target - target));
					if (pos + V128_BLOCKSIZE * high_offset <= ntargets
							&& target[pos] == goal)
						*foundp = 1;
					else
						*foundp = 0;
					return pos + V128_BLOCKSIZE * high_offset;
				} else
					break;
			} else {
				target += V128_BLOCKSIZE * high_offset;

				pos = Intersection_find_scalar(goal, target,
						(end_target - target));
				if (pos + V128_BLOCKSIZE * high_offset <= ntargets
						&& target[pos] == goal)
					*foundp = 1;
				else
					*foundp = 0;
				return pos + V128_BLOCKSIZE * high_offset;
			}
		}			// while-loop
		while (low_offset < high_offset) {
			mid_offset = (low_offset + high_offset) / 2;
			if (target[V128_BLOCKSIZE * mid_offset + V128_BLOCKSIZE - 1]
					< goal) {
				low_offset = mid_offset + 1;
			} else {
				high_offset = mid_offset;
			}
		}
		target += V128_BLOCKSIZE * high_offset;
	}
	// now search the block [target, target+block_size-1]
	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0, Q0, Q1, Q2, Q3;

	if (target[V128_BLOCKSIZE / 2 - 1] >= goal) {
		if (target[V128_BLOCKSIZE / 4 - 1] < goal)
			target += 32;
	} else {
		if (target[V128_BLOCKSIZE * 3 / 4 - 1] >= goal)
			target += 64;
		else
			target += 96;
	}
	Q0 = _mm_or_si128(
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 0), Match),
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 1), Match));
	Q1 = _mm_or_si128(
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 2), Match),
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 3), Match));
	Q2 = _mm_or_si128(
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 4), Match),
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 5), Match));
	Q3 = _mm_or_si128(
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 6), Match),
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 7), Match));

	F0 = _mm_or_si128(_mm_or_si128(Q0, Q1), _mm_or_si128(Q2, Q3));
	if (_mm_testz_si128(F0, F0))
		*foundp = 0;
	else
		*foundp = 1;
	return (target - init_target);
}

long simdgallop_v128_rough_plow(int *foundp, UINT4 goal, const UINT4 *target,
		long ntargets) {
	const UINT4 *end_target, *stop_target, *init_target;

	long low_offset = 0, mid_offset, high_offset = 1;
	long pos;

	init_target = target;
	stop_target = target + ntargets - V128_BLOCKSIZE;
	end_target = target + ntargets;

	if (target >= stop_target) {
		if ((pos = Intersection_find_scalar(goal, target, ntargets)) <= ntargets
				&& target[pos] == goal)
			*foundp = 1;
		else
			*foundp = 0;
		return pos;
	}

	if (target[V128_BLOCKSIZE - 1] < goal) {
		if (target + V128_BLOCKSIZE > stop_target) {
			pos = V128_BLOCKSIZE
					+ Intersection_find_scalar(goal, target + V128_BLOCKSIZE,
							ntargets - V128_BLOCKSIZE);
			if (pos <= ntargets && target[pos] == goal)
				*foundp = 1;
			else
				*foundp = 0;
			return pos;
		}
		/* Galloping search */
		while (target[V128_BLOCKSIZE * high_offset + V128_BLOCKSIZE - 1] < goal) {
			if (target + (high_offset << 1) * V128_BLOCKSIZE <= stop_target) {
				low_offset = high_offset;
				high_offset <<= 1;
			} else if (target + V128_BLOCKSIZE * (high_offset + 1)
					<= stop_target) {
				high_offset = (stop_target - target) / V128_BLOCKSIZE;
				if (target[V128_BLOCKSIZE * high_offset + V128_BLOCKSIZE - 1]
						< goal) {
					target += V128_BLOCKSIZE * high_offset;
					pos = Intersection_find_scalar(goal, target,
							(end_target - target));
					if (pos + V128_BLOCKSIZE * high_offset <= ntargets
							&& target[pos] == goal)
						*foundp = 1;
					else
						*foundp = 0;
					return pos + V128_BLOCKSIZE * high_offset;
				} else
					break;
			} else {
				target += V128_BLOCKSIZE * high_offset;

				pos = Intersection_find_scalar(goal, target,
						(end_target - target));
				if (pos + V128_BLOCKSIZE * high_offset <= ntargets
						&& target[pos] == goal)
					*foundp = 1;
				else
					*foundp = 0;
				return pos + V128_BLOCKSIZE * high_offset;
			}
		}			// while-loop
		while (low_offset < high_offset) {
			mid_offset = (low_offset + high_offset) / 2;
			if (target[V128_BLOCKSIZE * mid_offset + V128_BLOCKSIZE - 1]
					< goal) {
				low_offset = mid_offset + 1;
			} else {
				high_offset = mid_offset;
			}
		}
		target += V128_BLOCKSIZE * high_offset;
	}
	// now search the block [target, target+block_size-1]
	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0 = _mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 0),
			Match);
	UINT4 i = 0;
	while (_mm_testz_si128(F0, F0) && ++i < V128_BLOCKSIZE / SIMDWIDTH) {
		F0 = _mm_or_si128(F0,
				_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + i),
						Match));
	}
	if (_mm_testz_si128(F0, F0)) {
		*foundp = 0;
		i--;
	} else
		*foundp = 1;
	return (target + i * SIMDWIDTH - init_target);
}

/* return the exact position of intersection */
/* return the head-position of mathced block */
long simdgallop_v128_greedy(UINT4 goal, const UINT4 *target, long ntargets) {
	const UINT4 *end_target, *stop_target, *init_target;

	long low_offset = 0, mid_offset, high_offset = 1;
	long pos;

	init_target = target;
	stop_target = target + ntargets - V128_BLOCKSIZE;
	end_target = target + ntargets;

	if (target >= stop_target)
		return Intersection_find_scalar(goal, target, ntargets);

	if (target[V128_BLOCKSIZE - 1] < goal) {
		if (target + V128_BLOCKSIZE > stop_target)
			return V128_BLOCKSIZE
					+ Intersection_find_scalar(goal, target + V128_BLOCKSIZE,
							ntargets - V128_BLOCKSIZE);

		/* Galloping search */
		while (target[V128_BLOCKSIZE * high_offset + V128_BLOCKSIZE - 1] < goal) {
			if (target + (high_offset << 1) * V128_BLOCKSIZE <= stop_target) {
				low_offset = high_offset;
				high_offset <<= 1;
			} else if (target + V128_BLOCKSIZE * (high_offset + 1)
					<= stop_target) {
				high_offset = (stop_target - target) / V128_BLOCKSIZE;
				if (target[V128_BLOCKSIZE * high_offset + V128_BLOCKSIZE - 1]
						< goal) {
					target += V128_BLOCKSIZE * high_offset;
					return V128_BLOCKSIZE * high_offset
							+ Intersection_find_scalar(goal, target,
									(end_target - target));
				} else
					break;
			} else {
				target += V128_BLOCKSIZE * high_offset;
				return V128_BLOCKSIZE * high_offset
						+ Intersection_find_scalar(goal, target,
								(end_target - target));
			}
		}			// while-loop
		while (low_offset < high_offset) {
			mid_offset = (low_offset + high_offset) / 2;
			if (target[V128_BLOCKSIZE * mid_offset + V128_BLOCKSIZE - 1]
					< goal) {
				low_offset = mid_offset + 1;
			} else {
				high_offset = mid_offset;
			}
		}
		target += V128_BLOCKSIZE * high_offset;
	}
	// now search the block [target, target+block_size-1]
	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0, Q0, Q1, Q2, Q3;

	if (target[V128_BLOCKSIZE / 2 - 1] >= goal) {
		if (target[V128_BLOCKSIZE / 4 - 1] < goal)
			target += 32;
	} else {
		if (target[V128_BLOCKSIZE * 3 / 4 - 1] >= goal)
			target += 64;

		else
			target += 96;
	}
	Q0 = _mm_or_si128(
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 0), Match),
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 1), Match));
	Q1 = _mm_or_si128(
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 2), Match),
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 3), Match));
	Q2 = _mm_or_si128(
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 4), Match),
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 5), Match));
	Q3 = _mm_or_si128(
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 6), Match),
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 7), Match));

	F0 = _mm_or_si128(_mm_or_si128(Q0, Q1), _mm_or_si128(Q2, Q3));
	if (_mm_testz_si128(F0, F0))
		return (target - init_target);
	else {
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

long simdgallop_v256_exact(UINT4 goal, const UINT4 *target, long ntargets) {
	const UINT4 *end_target, *stop_target, *init_target;

	long low_offset = 0, mid_offset, high_offset = 1;

	init_target = target;
	stop_target = target + ntargets - V256_BLOCKSIZE;
	end_target = target + ntargets;

	if (target >= stop_target)
		return Intersection_find_scalar(goal, target, ntargets);

	if (target[V256_BLOCKSIZE - 1] < goal) {
		if (target + V256_BLOCKSIZE > stop_target)
			return V256_BLOCKSIZE
					+ Intersection_find_scalar(goal, target + V256_BLOCKSIZE,
							ntargets - V256_BLOCKSIZE);

		/* Galloping search */
		while (target[V256_BLOCKSIZE * high_offset + V256_BLOCKSIZE - 1] < goal) {
			if (target + (high_offset << 1) * V256_BLOCKSIZE <= stop_target) {
				low_offset = high_offset;
				high_offset <<= 1;
			} else if (target + V256_BLOCKSIZE * (high_offset + 1)
					<= stop_target) {
				high_offset = (stop_target - target) / V256_BLOCKSIZE;
				if (target[V256_BLOCKSIZE * high_offset + V256_BLOCKSIZE - 1]
						< goal) {
					target += V256_BLOCKSIZE * high_offset;
					return V256_BLOCKSIZE * high_offset
							+ Intersection_find_scalar(goal, target,
									(end_target - target));
				} else
					break;
			} else {
				target += V256_BLOCKSIZE * high_offset;
				return V256_BLOCKSIZE * high_offset
						+ Intersection_find_scalar(goal, target,
								(end_target - target));
			}
		}			// while-loop
		while (low_offset < high_offset) {
			mid_offset = (low_offset + high_offset) / 2;
			if (target[V256_BLOCKSIZE * mid_offset + V256_BLOCKSIZE - 1]
					< goal) {
				low_offset = mid_offset + 1;
			} else {
				high_offset = mid_offset;
			}
		}
		target += V256_BLOCKSIZE * high_offset;
	}
	// now search the block [target, target+block_size-1]
	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0, Q0, Q1, Q2, Q3;
	long pos;

	if (target[V256_BLOCKSIZE / 2 - 1] >= goal) {
		if (target[V256_BLOCKSIZE / 4 - 1] >= goal) {
			if (target[V256_BLOCKSIZE / 8 - 1] < goal)
				target += 32;
		} else {
			if (target[V256_BLOCKSIZE / 8 * 3 - 1] >= goal) {
				target += 64;
			} else
				target += 96;
		}
	} else {
		if (target[V256_BLOCKSIZE / 4 + 128 - 1] >= goal) {
			if (target[V256_BLOCKSIZE / 8 + 128 - 1] >= goal)
				target += 128;
			else
				target += 160;

		} else {
			if (target[V256_BLOCKSIZE / 8 * 3 + 128 - 1] >= goal)
				target += 192;
			else
				target += 224;
		}
	}
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

long simdgallop_v256_rough(int *foundp, UINT4 goal, const UINT4 *target,
		long ntargets) {
	const UINT4 *end_target, *stop_target, *init_target;

	long low_offset = 0, mid_offset, high_offset = 1;
	long pos;

	init_target = target;
	stop_target = target + ntargets - V256_BLOCKSIZE;
	end_target = target + ntargets;

	if (target >= stop_target) {
		if ((pos = Intersection_find_scalar(goal, target, ntargets)) <= ntargets
				&& target[pos] == goal)
			*foundp = 1;
		else
			*foundp = 0;
		return pos;
	}

	if (target[V256_BLOCKSIZE - 1] < goal) {
		if (target + V256_BLOCKSIZE > stop_target) {
			pos = V256_BLOCKSIZE
					+ Intersection_find_scalar(goal, target + V256_BLOCKSIZE,
							ntargets - V256_BLOCKSIZE);
			if (pos <= ntargets && target[pos] == goal)
				*foundp = 1;
			else
				*foundp = 0;
			return pos;
		}
		/* Galloping search */
		while (target[V256_BLOCKSIZE * high_offset + V256_BLOCKSIZE - 1] < goal) {
			if (target + (high_offset << 1) * V256_BLOCKSIZE <= stop_target) {
				low_offset = high_offset;
				high_offset <<= 1;
			} else if (target + V256_BLOCKSIZE * (high_offset + 1)
					<= stop_target) {
				high_offset = (stop_target - target) / V256_BLOCKSIZE;
				if (target[V256_BLOCKSIZE * high_offset + V256_BLOCKSIZE - 1]
						< goal) {
					target += V256_BLOCKSIZE * high_offset;
					pos = Intersection_find_scalar(goal, target,
							(end_target - target));
					if (pos + V256_BLOCKSIZE * high_offset <= ntargets
							&& target[pos] == goal)
						*foundp = 1;
					else
						*foundp = 0;
					return pos + V256_BLOCKSIZE * high_offset;
				} else
					break;
			} else {
				target += V256_BLOCKSIZE * high_offset;

				pos = Intersection_find_scalar(goal, target,
						(end_target - target));
				if (pos + V256_BLOCKSIZE * high_offset <= ntargets
						&& target[pos] == goal)
					*foundp = 1;
				else
					*foundp = 0;
				return pos + V256_BLOCKSIZE * high_offset;
			}
		}			// while-loop
		while (low_offset < high_offset) {
			mid_offset = (low_offset + high_offset) / 2;
			if (target[V256_BLOCKSIZE * mid_offset + V256_BLOCKSIZE - 1]
					< goal) {
				low_offset = mid_offset + 1;
			} else {
				high_offset = mid_offset;
			}
		}
		target += V256_BLOCKSIZE * high_offset;
	}
// now search the block [target, target+block_size-1]

	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0, Q0, Q1, Q2, Q3;

	if (target[V256_BLOCKSIZE / 2 - 1] >= goal) {
		if (target[V256_BLOCKSIZE / 4 - 1] >= goal) {
			if (target[V256_BLOCKSIZE / 8 - 1] < goal)
				target += 32;
		} else {
			if (target[V256_BLOCKSIZE / 8 * 3 - 1] >= goal)
				target += 64;
			else
				target += 96;
		}
	} else {
		if (target[V256_BLOCKSIZE / 4 + 128 - 1] >= goal) {
			if (target[V256_BLOCKSIZE / 8 + 128 - 1] >= goal)
				target += 128;
			else
				target += 160;
		} else {
			if (target[V256_BLOCKSIZE / 8 * 3 + 128 - 1] >= goal)
				target += 192;
			else
				target += 224;
		}
	}
	Q0 = _mm_or_si128(
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 0), Match),
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 1), Match));
	Q1 = _mm_or_si128(
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 2), Match),
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 3), Match));
	Q2 = _mm_or_si128(
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 4), Match),
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 5), Match));
	Q3 = _mm_or_si128(
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 6), Match),
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 7), Match));

	F0 = _mm_or_si128(_mm_or_si128(Q0, Q1), _mm_or_si128(Q2, Q3));
	if (_mm_testz_si128(F0, F0))
		*foundp = 0;
	else
		*foundp = 1;
	return (target - init_target);
}

long simdgallop_v256_rough_plow(int *foundp, UINT4 goal, const UINT4 *target,
		long ntargets) {
	const UINT4 *end_target, *stop_target, *init_target;

	long low_offset = 0, mid_offset, high_offset = 1;
	long pos;

	init_target = target;
	stop_target = target + ntargets - V256_BLOCKSIZE;
	end_target = target + ntargets;

	if (target >= stop_target) {
		if ((pos = Intersection_find_scalar(goal, target, ntargets)) <= ntargets
				&& target[pos] == goal)
			*foundp = 1;
		else
			*foundp = 0;
		return pos;
	}

	if (target[V256_BLOCKSIZE - 1] < goal) {
		if (target + V256_BLOCKSIZE > stop_target) {
			pos = V256_BLOCKSIZE
					+ Intersection_find_scalar(goal, target + V256_BLOCKSIZE,
							ntargets - V256_BLOCKSIZE);
			if (pos <= ntargets && target[pos] == goal)
				*foundp = 1;
			else
				*foundp = 0;
			return pos;
		}
		/* Galloping search */
		while (target[V256_BLOCKSIZE * high_offset + V256_BLOCKSIZE - 1] < goal) {
			if (target + (high_offset << 1) * V256_BLOCKSIZE <= stop_target) {
				low_offset = high_offset;
				high_offset <<= 1;
			} else if (target + V256_BLOCKSIZE * (high_offset + 1)
					<= stop_target) {
				high_offset = (stop_target - target) / V256_BLOCKSIZE;
				if (target[V256_BLOCKSIZE * high_offset + V256_BLOCKSIZE - 1]
						< goal) {
					target += V256_BLOCKSIZE * high_offset;
					pos = Intersection_find_scalar(goal, target,
							(end_target - target));
					if (pos + V256_BLOCKSIZE * high_offset <= ntargets
							&& target[pos] == goal)
						*foundp = 1;
					else
						*foundp = 0;
					return pos + V256_BLOCKSIZE * high_offset;
				} else
					break;
			} else {
				target += V256_BLOCKSIZE * high_offset;

				pos = Intersection_find_scalar(goal, target,
						(end_target - target));
				if (pos + V256_BLOCKSIZE * high_offset <= ntargets
						&& target[pos] == goal)
					*foundp = 1;
				else
					*foundp = 0;
				return pos + V256_BLOCKSIZE * high_offset;
			}
		}			// while-loop
		while (low_offset < high_offset) {
			mid_offset = (low_offset + high_offset) / 2;
			if (target[V256_BLOCKSIZE * mid_offset + V256_BLOCKSIZE - 1]
					< goal) {
				low_offset = mid_offset + 1;
			} else {
				high_offset = mid_offset;
			}
		}
		target += V256_BLOCKSIZE * high_offset;
	}
// now search the block [target, target+block_size-1]

	__m128i Match = _mm_set1_epi32(goal);

	__m128i F0 = _mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 0),
			Match);
	UINT4 i = 0;
	while (_mm_testz_si128(F0, F0) && ++i < V256_BLOCKSIZE / SIMDWIDTH) {
		F0 = _mm_or_si128(F0,
				_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + i),
						Match));
	}
	if (_mm_testz_si128(F0, F0)) {
		*foundp = 0;
		i--;
	} else
		*foundp = 1;
	return (target + i * SIMDWIDTH - init_target);
}

long simdgallop_v512_exact(UINT4 goal, const UINT4 *target, long ntargets) {
	const UINT4 *end_target, *stop_target, *init_target;

	long low_offset = 0, mid_offset, high_offset = 1;

	init_target = target;
	stop_target = target + ntargets - V512_BLOCKSIZE;
	end_target = target + ntargets;

	if (target >= stop_target)
		return Intersection_find_scalar(goal, target, ntargets);

	if (target[V512_BLOCKSIZE - 1] < goal) {
		if (target + V512_BLOCKSIZE > stop_target)
			return V512_BLOCKSIZE
					+ Intersection_find_scalar(goal, target + V512_BLOCKSIZE,
							ntargets - V512_BLOCKSIZE);

		/* Galloping search */
		while (target[V512_BLOCKSIZE * high_offset + V512_BLOCKSIZE - 1] < goal) {
			if (target + (high_offset << 1) * V512_BLOCKSIZE <= stop_target) {
				low_offset = high_offset;
				high_offset <<= 1;
			} else if (target + V512_BLOCKSIZE * (high_offset + 1)
					<= stop_target) {
				high_offset = (stop_target - target) / V512_BLOCKSIZE;
				if (target[V512_BLOCKSIZE * high_offset + V512_BLOCKSIZE - 1]
						< goal) {
					target += V512_BLOCKSIZE * high_offset;
					return V512_BLOCKSIZE * high_offset
							+ Intersection_find_scalar(goal, target,
									(end_target - target));
				} else
					break;
			} else {
				target += V512_BLOCKSIZE * high_offset;
				return V512_BLOCKSIZE * high_offset
						+ Intersection_find_scalar(goal, target,
								(end_target - target));
			}
		}			// while-loop
		while (low_offset < high_offset) {
			mid_offset = (low_offset + high_offset) / 2;
			if (target[V512_BLOCKSIZE * mid_offset + V512_BLOCKSIZE - 1]
					< goal) {
				low_offset = mid_offset + 1;
			} else {
				high_offset = mid_offset;
			}
		}
		target += V512_BLOCKSIZE * high_offset;
	}
	// now search the block [target, target+block_size-1]
	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0, Q0, Q1, Q2, Q3;
	long pos;

	if (target[V512_BLOCKSIZE / 2 - 1] >= goal) {
		if (target[V512_BLOCKSIZE / 4 - 1] >= goal) {
			if (target[V512_BLOCKSIZE / 8 - 1] >= goal) {
				if (target[V512_BLOCKSIZE / 16 - 1] >= goal) {
				} else {
					target += 32;
				}
			} else {
				if (target[V512_BLOCKSIZE / 8 + 32 - 1] >= goal) {
					target += 64;
				} else {
					target += 96;
				}
			}
		} else {
			if (target[V512_BLOCKSIZE / 4 + 64 - 1] >= goal) {
				if (target[V512_BLOCKSIZE / 8 + 96 - 1] >= goal) {
					target += 128;
				} else {
					target += 160;
				}
			} else {
				if (target[V512_BLOCKSIZE / 4 + 96 - 1] >= goal) {
					target += 192;
				} else {
					target += 224;
				}
			}
		}
	} else {
		if (target[V512_BLOCKSIZE / 4 + 256 - 1] >= goal) {
			if (target[V512_BLOCKSIZE / 8 + 256 - 1] >= goal) {
				if (target[V512_BLOCKSIZE / 16 + 256 - 1] >= goal) {
					target += 256;
				} else {
					target += 288;
				}
			} else {
				if (target[V512_BLOCKSIZE / 8 + 288 - 1] >= goal) {
					target += 320;
				} else {
					target += 352;
				}
			}
		} else {
			if (target[V512_BLOCKSIZE / 4 + 320 - 1] >= goal) {
				if (target[V512_BLOCKSIZE / 8 + 352 - 1] >= goal) {
					target += 384;
				} else {
					target += 416;
				}
			} else {
				if (target[V512_BLOCKSIZE / 4 + 352 - 1] >= goal) {
					target += 448;
				} else {
					target += 480;
				}
			}
		}
	}

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

long simdgallop_v512_rough(int *foundp, UINT4 goal, const UINT4 *target,
		long ntargets) {
	const UINT4 *end_target, *stop_target, *init_target;

	long low_offset = 0, mid_offset, high_offset = 1;
	long pos;

	init_target = target;
	stop_target = target + ntargets - V512_BLOCKSIZE;
	end_target = target + ntargets;

	if (target >= stop_target) {
		if ((pos = Intersection_find_scalar(goal, target, ntargets)) <= ntargets
				&& target[pos] == goal)
			*foundp = 1;
		else
			*foundp = 0;
		return pos;
	}

	if (target[V512_BLOCKSIZE - 1] < goal) {
		if (target + V512_BLOCKSIZE > stop_target) {
			pos = V512_BLOCKSIZE
					+ Intersection_find_scalar(goal, target + V512_BLOCKSIZE,
							ntargets - V512_BLOCKSIZE);
			if (pos <= ntargets && target[pos] == goal)
				*foundp = 1;
			else
				*foundp = 0;
			return pos;
		}
		/* Galloping search */
		while (target[V512_BLOCKSIZE * high_offset + V512_BLOCKSIZE - 1] < goal) {
			if (target + (high_offset << 1) * V512_BLOCKSIZE <= stop_target) {
				low_offset = high_offset;
				high_offset <<= 1;
			} else if (target + V512_BLOCKSIZE * (high_offset + 1)
					<= stop_target) {
				high_offset = (stop_target - target) / V512_BLOCKSIZE;
				if (target[V512_BLOCKSIZE * high_offset + V512_BLOCKSIZE - 1]
						< goal) {
					target += V512_BLOCKSIZE * high_offset;
					pos = Intersection_find_scalar(goal, target,
							(end_target - target));
					if (pos + V512_BLOCKSIZE * high_offset <= ntargets
							&& target[pos] == goal)
						*foundp = 1;
					else
						*foundp = 0;
					return pos + V512_BLOCKSIZE * high_offset;
				} else
					break;
			} else {
				target += V512_BLOCKSIZE * high_offset;

				pos = Intersection_find_scalar(goal, target,
						(end_target - target));
				if (pos + V512_BLOCKSIZE * high_offset <= ntargets
						&& target[pos] == goal)
					*foundp = 1;
				else
					*foundp = 0;
				return pos + V512_BLOCKSIZE * high_offset;
			}
		}			// while-loop
		while (low_offset < high_offset) {
			mid_offset = (low_offset + high_offset) / 2;
			if (target[V512_BLOCKSIZE * mid_offset + V512_BLOCKSIZE - 1]
					< goal) {
				low_offset = mid_offset + 1;
			} else {
				high_offset = mid_offset;
			}
		}
		target += V512_BLOCKSIZE * high_offset;
	}
// now search the block [target, target+block_size-1]

	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0, Q0, Q1, Q2, Q3;

	if (target[V512_BLOCKSIZE / 2 - 1] >= goal) {
		if (target[V512_BLOCKSIZE / 4 - 1] >= goal) {
			if (target[V512_BLOCKSIZE / 8 - 1] >= goal) {
				if (target[V512_BLOCKSIZE / 16 - 1] >= goal) {
				} else {
					target += 32;
				}
			} else {
				if (target[V512_BLOCKSIZE / 8 + 32 - 1] >= goal) {
					target += 64;
				} else {
					target += 96;
				}
			}
		} else {
			if (target[V512_BLOCKSIZE / 4 + 64 - 1] >= goal) {
				if (target[V512_BLOCKSIZE / 8 + 96 - 1] >= goal) {
					target += 128;
				} else {
					target += 160;
				}
			} else {
				if (target[V512_BLOCKSIZE / 4 + 96 - 1] >= goal) {
					target += 192;
				} else {
					target += 224;
				}
			}
		}
	} else {
		if (target[V512_BLOCKSIZE / 4 + 256 - 1] >= goal) {
			if (target[V512_BLOCKSIZE / 8 + 256 - 1] >= goal) {
				if (target[V512_BLOCKSIZE / 16 + 256 - 1] >= goal) {
					target += 256;
				} else {
					target += 288;
				}
			} else {
				if (target[V512_BLOCKSIZE / 8 + 288 - 1] >= goal) {
					target += 320;
				} else {
					target += 352;
				}
			}
		} else {
			if (target[V512_BLOCKSIZE / 4 + 320 - 1] >= goal) {
				if (target[V512_BLOCKSIZE / 8 + 352 - 1] >= goal) {
					target += 384;
				} else {
					target += 416;
				}
			} else {
				if (target[V512_BLOCKSIZE / 4 + 352 - 1] >= goal) {
					target += 448;
				} else {
					target += 480;
				}
			}
		}
	}

	Q0 = _mm_or_si128(
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 0), Match),
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 1), Match));
	Q1 = _mm_or_si128(
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 2), Match),
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 3), Match));
	Q2 = _mm_or_si128(
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 4), Match),
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 5), Match));
	Q3 = _mm_or_si128(
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 6), Match),
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 7), Match));

	F0 = _mm_or_si128(_mm_or_si128(Q0, Q1), _mm_or_si128(Q2, Q3));
	if (_mm_testz_si128(F0, F0))
		*foundp = 0;
	else
		*foundp = 1;
	return (target - init_target);
}

long simdgallop_v512_rough_plow(int *foundp, UINT4 goal, const UINT4 *target,
		long ntargets) {
	const UINT4 *end_target, *stop_target, *init_target;

	long low_offset = 0, mid_offset, high_offset = 1;
	long pos;

	init_target = target;
	stop_target = target + ntargets - V512_BLOCKSIZE;
	end_target = target + ntargets;

	if (target >= stop_target) {
		if ((pos = Intersection_find_scalar(goal, target, ntargets)) <= ntargets
				&& target[pos] == goal)
			*foundp = 1;
		else
			*foundp = 0;
		return pos;
	}

	if (target[V512_BLOCKSIZE - 1] < goal) {
		if (target + V512_BLOCKSIZE > stop_target) {
			pos = V512_BLOCKSIZE
					+ Intersection_find_scalar(goal, target + V512_BLOCKSIZE,
							ntargets - V512_BLOCKSIZE);
			if (pos <= ntargets && target[pos] == goal)
				*foundp = 1;
			else
				*foundp = 0;
			return pos;
		}
		/* Galloping search */
		while (target[V512_BLOCKSIZE * high_offset + V512_BLOCKSIZE - 1] < goal) {
			if (target + (high_offset << 1) * V512_BLOCKSIZE <= stop_target) {
				low_offset = high_offset;
				high_offset <<= 1;
			} else if (target + V512_BLOCKSIZE * (high_offset + 1)
					<= stop_target) {
				high_offset = (stop_target - target) / V512_BLOCKSIZE;
				if (target[V512_BLOCKSIZE * high_offset + V512_BLOCKSIZE - 1]
						< goal) {
					target += V512_BLOCKSIZE * high_offset;
					pos = Intersection_find_scalar(goal, target,
							(end_target - target));
					if (pos + V512_BLOCKSIZE * high_offset <= ntargets
							&& target[pos] == goal)
						*foundp = 1;
					else
						*foundp = 0;
					return pos + V512_BLOCKSIZE * high_offset;
				} else
					break;
			} else {
				target += V512_BLOCKSIZE * high_offset;

				pos = Intersection_find_scalar(goal, target,
						(end_target - target));
				if (pos + V512_BLOCKSIZE * high_offset <= ntargets
						&& target[pos] == goal)
					*foundp = 1;
				else
					*foundp = 0;
				return pos + V512_BLOCKSIZE * high_offset;
			}
		}			// while-loop
		while (low_offset < high_offset) {
			mid_offset = (low_offset + high_offset) / 2;
			if (target[V512_BLOCKSIZE * mid_offset + V512_BLOCKSIZE - 1]
					< goal) {
				low_offset = mid_offset + 1;
			} else {
				high_offset = mid_offset;
			}
		}
		target += V512_BLOCKSIZE * high_offset;
	}
// now search the block [target, target+block_size-1]

	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0 = _mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 0),
			Match);
	UINT4 i = 0;
	while (_mm_testz_si128(F0, F0) && ++i < V512_BLOCKSIZE / SIMDWIDTH) {
		F0 = _mm_or_si128(F0,
				_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + i),
						Match));
	}
	if (_mm_testz_si128(F0, F0)) {
		*foundp = 0;
		i--;
	} else
		*foundp = 1;
	return (target + i * SIMDWIDTH - init_target);
}

////////////////////////////////////////////////////////////////////
/* small_vs_small */
template<setIntersectionFunction FUNCTION>
void svs(const mySet &sets, std::vector<uint32_t> &out) {
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

/* set_vs_set */
template<intersectionfindfunction FINDFUNCTION>
void SvS_exact(const mySet &sets, std::vector<uint32_t> &out) {
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
void SvS_rough(const mySet &sets, std::vector<uint32_t> &out) {
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

/* swapping_set_vs_set */
template<intersectionfindfunction FINDFUNCTION>
void s_SvS_exact(const mySet& sets, std::vector<uint32_t>& out) {
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
		vsets[currentset].second += FINDFUNCTION(
				vsets[0].first->at(vsets[0].second),
				vsets[currentset].first->data() + vsets[currentset].second,
				vsets[currentset].first->size() - 1 - vsets[currentset].second);

		if (vsets[currentset].first->at(vsets[currentset].second)
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
		} else if (_UNLIKELY(vsets[currentset].first->back() < *value)) {
			break;
		} else if (_UNLIKELY(++vsets[0].second == vsets[0].first->size()))
			break; // note here @vsets[0].second alreadly self-add 1
		else {
			// greater
			//vsets[0].second++;
			sort_remaining();
			currentset = 1;
		}
	}
	out.resize(count);
}

template<flaggedintersectionfindfunction FINDFUNCTION>
void s_SvS_rough(const mySet& sets, std::vector<uint32_t>& out) {
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
	int foundp;

	while (vsets[0].second != vsets[0].first->size()) {
		vsets[currentset].second += FINDFUNCTION(&foundp,
				vsets[0].first->at(vsets[0].second),
				vsets[currentset].first->data() + vsets[currentset].second,
				vsets[currentset].first->size() - 1 - vsets[currentset].second);

		if (foundp == 1) {
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
		} else if (_UNLIKELY(vsets[currentset].first->back() < *value)) {
			break;
		} else if (_UNLIKELY(++vsets[0].second == vsets[0].first->size()))
			break; // note here @vsets[0].second alreadly self-add 1
		else {
			// greater
			//vsets[0].second++;
			sort_remaining();
			currentset = 1;
		}
	}
	out.resize(count);
}

/* for now we only implement its scalar version */
/* adapitve */
void adp_scalar(const mySet &sets, std::vector<uint32_t> &out) {
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
/* small_adapitve */
void s_adp_scalar(const mySet &sets, std::vector<uint32_t> &out) {
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

/* sequential */
template<intersectionfindfunction FINDFUNCTION>
void sql_exact(const mySet &sets, std::vector<uint32_t> &out) {
	out.resize(sets.begin()->size());
	auto value = out.begin();
	*value = sets.begin()->at(0);
	size_t count = 0, elimset = 0, currentset = 1;
	mySet::iterator it = sets.begin();
	mySet::iterator it_elim = it++;
	std::vector<size_t> index(sets.size(), 0);

	while (_LIKELY(index[currentset] < it->size())) {
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
		} else if (_UNLIKELY(it->back() < *value))
			break;
		else if (_UNLIKELY(++index[elimset] == it_elim->size()))
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

template<flaggedintersectionfindfunction FINDFUNCTION>
void sql_rough(const mySet &sets, std::vector<uint32_t> &out) {
	out.resize(sets.begin()->size());
	auto value = out.begin();
	*value = sets.begin()->at(0);
	size_t count = 0, elimset = 0, currentset = 1;
	mySet::iterator it = sets.begin();
	mySet::iterator it_elim = it++;
	std::vector<size_t> index(sets.size(), 0);
	int foundp;

	while (_LIKELY(index[currentset] < it->size())) {
		index[currentset] += FINDFUNCTION(&foundp, *value,
				it->data() + index[currentset],
				it->size() - 1 - index[currentset]);
		if (foundp == 1) {
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
		} else if (_UNLIKELY(it->back() < *value))
			break;
		else if (_UNLIKELY(++index[elimset] == it_elim->size()))
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

/* small_sequential */
template<intersectionfindfunction FINDFUNCTION>
void s_sql_exact(const mySet &sets, std::vector<uint32_t> &out) {
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
	*value = sets.begin()->at(0);
	size_t count = 0, currentset = 1;

	while (vsets[0].second != vsets[0].first->size()) {
		vsets[currentset].second += FINDFUNCTION(
				vsets[0].first->at(vsets[0].second),
				vsets[currentset].first->data() + vsets[currentset].second,
				vsets[currentset].first->size() - 1 - vsets[currentset].second);
		if (vsets[currentset].first->at(vsets[currentset].second)
				== vsets[0].first->at(vsets[0].second)) {
			vsets[currentset].second++;
			if (++currentset == vsets.size()) {
				*value++ = vsets[0].first->at(vsets[0].second);
				count++;
				vsets[0].second++;
				sort_remaining();
				currentset = 1;
			}
		} else if (_UNLIKELY(
				vsets[currentset].first->back()
						< vsets[0].first->at(vsets[0].second)))
			break;
		else {
			vsets[0].second++;
			sort_remaining();
			currentset = 1;
		}
	}
	out.resize(count);
}

template<flaggedintersectionfindfunction FINDFUNCTION>
void s_sql_rough(const mySet &sets, std::vector<uint32_t> &out) {
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
	*value = sets.begin()->at(0);
	size_t count = 0, currentset = 1;
	int foundp;

	while (vsets[0].second != vsets[0].first->size()) {
		vsets[currentset].second += FINDFUNCTION(&foundp,
				vsets[0].first->at(vsets[0].second),
				vsets[currentset].first->data() + vsets[currentset].second,
				vsets[currentset].first->size() - 1 - vsets[currentset].second);
		if (foundp == 1) {
			vsets[currentset].second++;
			if (++currentset == vsets.size()) {
				*value++ = vsets[0].first->at(vsets[0].second);
				count++;
				vsets[0].second++;
				sort_remaining();
				currentset = 1;
			}
		} else if (_UNLIKELY(
				vsets[currentset].first->back()
						< vsets[0].first->at(vsets[0].second)))
			break;
		else {
			vsets[0].second++;
			sort_remaining();
			currentset = 1;
		}
	}
	out.resize(count);
}

template<intersectionfindfunction FINDFUNCTION>
void max_exact(const mySet &sets, std::vector<uint32_t> &out) {
	out.resize(sets.begin()->size());
	auto value = out.begin();
	*value = sets.begin()->at(0);
	size_t count = 0, intersect_count = 0, elimset = 0, currentset = 1;
	mySet::iterator it = sets.begin();
	mySet::iterator it_start = it++, it_elim = it_start;
	std::vector<size_t> index(sets.size(), 0);

	while (_LIKELY(index[currentset] < it->size())) {
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
		} else if (_UNLIKELY(it->back() < *value))
			break;
		else if (_UNLIKELY(++index[elimset] == it_elim->size()))
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

template<flaggedintersectionfindfunction FINDFUNCTION>
void max_rough(const mySet &sets, std::vector<uint32_t> &out) {
	out.resize(sets.begin()->size());
	auto value = out.begin();
	*value = sets.begin()->at(0);
	size_t count = 0, intersect_count = 0, elimset = 0, currentset = 1;
	mySet::iterator it = sets.begin();
	mySet::iterator it_start = it++, it_elim = it_start;
	std::vector<size_t> index(sets.size(), 0);
	int foundp;

	while (_LIKELY(index[currentset] < it->size())) {
		index[currentset] += FINDFUNCTION(&foundp, *value,
				it->data() + index[currentset],
				it->size() - 1 - index[currentset]);

		if (foundp == 1) {
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
		} else if (_UNLIKELY(it->back() < *value))
			break;
		else if (_UNLIKELY(++index[elimset] == it_elim->size()))
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
void BY(const mySet &sets, std::vector<uint32_t> &out) {
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
