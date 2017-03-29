/*
 * msis_searchingalg.hpp
 *
 *  Created on: 2017Äê3ÔÂ28ÈÕ
 *      Author: John
 */

#ifndef INCLUDE_MSIS_LINEAR_HPP_
#define INCLUDE_MSIS_LINEAR_HPP_

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

long simdlinear_4_exact(UINT4 goal, const UINT4 *target, long ntargets) {
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

long simdlinear_4_rough(int*foundp, UINT4 goal, const UINT4 *target,
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

long simdlinear_8_exact(UINT4 goal, const UINT4 *target, long ntargets) {
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

long simdlinear_8_rough(int*foundp, UINT4 goal, const UINT4 *target,
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

long simdlinear_8_rough_plow(int*foundp, UINT4 goal, const UINT4 *target,
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
		F0 = _mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + i), Match);
	}
	if (_mm_testz_si128(F0, F0)) {
		*foundp = 0;
		i = 0;
	} else
		*foundp = 1;
	return (target + i * SIMDWIDTH - init_target);
}

long simdlinear_16_exact(UINT4 goal, const UINT4 *target, long ntargets) {
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

long simdlinear_16_rough(int*foundp, UINT4 goal, const UINT4 *target,
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

long simdlinear_16_rough_plow(int*foundp, UINT4 goal, const UINT4 *target,
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
		F0 = _mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + i), Match);
	}
	if (_mm_testz_si128(F0, F0)) {
		*foundp = 0;
		i = 0;
	} else
		*foundp = 1;
	return (target + i * SIMDWIDTH - init_target);
}

long simdlinear_32_exact(UINT4 goal, const UINT4 *target, long ntargets) {
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

long simdlinear_32_rough(int*foundp, UINT4 goal, const UINT4 *target,
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

long simdlinear_32_rough_plow(int*foundp, UINT4 goal, const UINT4 *target,
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
		F0 = _mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + i), Match);
	}
	if (_mm_testz_si128(F0, F0)) {
		*foundp = 0;
		i = 0;
	} else
		*foundp = 1;
	return (target + i * SIMDWIDTH - init_target);
}

long simdlinear_64_exact(UINT4 goal, const UINT4 *target, long ntargets) {
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

long simdlinear_64_rough(int*foundp, UINT4 goal, const UINT4 *target,
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

long simdlinear_64_rough_plow(int*foundp, UINT4 goal, const UINT4 *target,
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
		F0 = _mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + i), Match);
	}
	if (_mm_testz_si128(F0, F0)) {
		*foundp = 0;
		i = 0;
	} else
		*foundp = 1;
	return (target + i * SIMDWIDTH - init_target);
}

long simdlinear_128_exact(UINT4 goal, const UINT4 *target, long ntargets) {
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

long simdlinear_128_rough(int*foundp, UINT4 goal, const UINT4 *target,
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

long simdlinear_128_rough_plow(int*foundp, UINT4 goal, const UINT4 *target,
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
		F0 = _mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + i), Match);
	}
	if (_mm_testz_si128(F0, F0)) {
		*foundp = 0;
		i = 0;
	} else
		*foundp = 1;
	return (target + i * SIMDWIDTH - init_target);
}

long simdlinear_256_exact(UINT4 goal, const UINT4 *target, long ntargets) {
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

long simdlinear_256_rough(int*foundp, UINT4 goal, const UINT4 *target,
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

long simdlinear_256_rough_plow(int*foundp, UINT4 goal, const UINT4 *target,
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
		F0 = _mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + i), Match);
	}
	if (_mm_testz_si128(F0, F0)) {
		*foundp = 0;
		i = 0;
	} else
		*foundp = 1;
	return (target + i * SIMDWIDTH - init_target);
}

long simdlinear_512_exact(UINT4 goal, const UINT4 *target, long ntargets) {
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

long simdlinear_512_rough(int*foundp, UINT4 goal, const UINT4 *target,
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

long simdlinear_512_rough_plow(int*foundp, UINT4 goal, const UINT4 *target,
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
		F0 = _mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + i), Match);
	}
	if (_mm_testz_si128(F0, F0)) {
		*foundp = 0;
		i = 0;
	} else
		*foundp = 1;
	return (target + i * SIMDWIDTH - init_target);
}
}

#endif /* INCLUDE_MSIS_LINEAR_HPP_ */
