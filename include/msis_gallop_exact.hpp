/*
 * msis_gallop.hpp
 *
 *  Created on: 2017Äê3ÔÂ28ÈÕ
 *      Author: John
 */

#ifndef INCLUDE_MSIS_GALLOP_EXACT_HPP_
#define INCLUDE_MSIS_GALLOP_EXACT_HPP_

#include "msis_linear.hpp"

namespace msis/*MultiSet InterSection*/{
long simdgallop_4_exact(UINT4 goal, const UINT4 *target, long ntargets) {
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

long simdgallop_8_exact(UINT4 goal, const UINT4 *target, long ntargets) {
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

long simdgallop_16_exact(UINT4 goal, const UINT4 *target, long ntargets) {
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

long simdgallop_32_exact(UINT4 goal, const UINT4 *target, long ntargets) {
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

long simdgallop_64_exact(UINT4 goal, const UINT4 *target, long ntargets) {
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

long simdgallop_128_exact(UINT4 goal, const UINT4 *target, long ntargets) {
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

long simdgallop_256_exact(UINT4 goal, const UINT4 *target, long ntargets) {
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

long simdgallop_512_exact(UINT4 goal, const UINT4 *target, long ntargets) {
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

long simdgallop_8_rough_plow(int *foundp, UINT4 goal, const UINT4 *target,
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
		F0 = _mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + i), Match);
	}
	if (_mm_testz_si128(F0, F0)) {
		*foundp = 0;
		i = 0;
	} else
		*foundp = 1;
	return (target + i * SIMDWIDTH - init_target);
}

long simdgallop_16_rough_plow(int *foundp, UINT4 goal, const UINT4 *target,
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
		F0 = _mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + i), Match);
	}
	if (_mm_testz_si128(F0, F0)) {
		*foundp = 0;
		i = 0;
	} else
		*foundp = 1;
	return (target + i * SIMDWIDTH - init_target);
}

long simdgallop_32_rough_plow(int *foundp, UINT4 goal, const UINT4 *target,
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
		F0 = _mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + i), Match);
	}
	if (_mm_testz_si128(F0, F0)) {
		*foundp = 0;
		i = 0;
	} else
		*foundp = 1;
	return (target + i * SIMDWIDTH - init_target);
}

long simdgallop_64_rough_plow(int *foundp, UINT4 goal, const UINT4 *target,
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
		F0 = _mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + i), Match);
	}
	if (_mm_testz_si128(F0, F0)) {
		*foundp = 0;
		i = 0;
	} else
		*foundp = 1;
	return (target + i * SIMDWIDTH - init_target);
}

/* return the exact position of intersection */
long simdgallop_128_rough_plow(int *foundp, UINT4 goal, const UINT4 *target,
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
		F0 = _mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + i), Match);
	}
	if (_mm_testz_si128(F0, F0)) {
		*foundp = 0;
		i = 0;
	} else
		*foundp = 1;
	return (target + i * SIMDWIDTH - init_target);
}

long simdgallop_128_rough_plow_2(int *foundp, UINT4 goal, const UINT4 *target,
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
	__m128i F0 = _mm_or_si128(
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 0), Match),
			_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 1), Match));
	UINT4 i = 0;
	while (_mm_testz_si128(F0, F0) && ++i < V128_BLOCKSIZE / (SIMDWIDTH * 2)) {
		F0 = _mm_or_si128(
				_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 2 * i),
						Match),
				_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 1 + 2 * i),
						Match));
	}
	if (_mm_testz_si128(F0, F0)) {
		*foundp = 0;
		i = 0;
	} else
		*foundp = 1;
	return (target + i * SIMDWIDTH * 2 - init_target);
}

long simdgallop_128_rough_plow_4(int *foundp, UINT4 goal, const UINT4 *target,
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
	UINT4 i = 0;
	while (_mm_testz_si128(F0, F0) && ++i < V128_BLOCKSIZE / (SIMDWIDTH * 4)) {
		F0 = _mm_or_si128(
				_mm_or_si128(
						_mm_cmpeq_epi32(
								_mm_lddqu_si128((__m128i *) target + 4 * i + 0),
								Match),
						_mm_cmpeq_epi32(
								_mm_lddqu_si128((__m128i *) target + 4 * i + 1),
								Match)),
				_mm_or_si128(
						_mm_cmpeq_epi32(
								_mm_lddqu_si128((__m128i *) target + 4 * i + 2),
								Match),
						_mm_cmpeq_epi32(
								_mm_lddqu_si128((__m128i *) target + 4 * i + 3),
								Match)));
	}
	if (_mm_testz_si128(F0, F0)) {
		*foundp = 0;
		i = 0;
	} else
		*foundp = 1;
	return (target + i * SIMDWIDTH * 4 - init_target);
}

long simdgallop_128_rough_plow_8(int *foundp, UINT4 goal, const UINT4 *target,
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
	__m128i F0 = _mm_or_si128(
			_mm_or_si128(
					_mm_or_si128(
							_mm_cmpeq_epi32(
									_mm_lddqu_si128((__m128i *) target + 0),
									Match),
							_mm_cmpeq_epi32(
									_mm_lddqu_si128((__m128i *) target + 1),
									Match)),
					_mm_or_si128(
							_mm_cmpeq_epi32(
									_mm_lddqu_si128((__m128i *) target + 2),
									Match),
							_mm_cmpeq_epi32(
									_mm_lddqu_si128((__m128i *) target + 3),
									Match))),
			_mm_or_si128(
					_mm_or_si128(
							_mm_cmpeq_epi32(
									_mm_lddqu_si128((__m128i *) target + 4),
									Match),
							_mm_cmpeq_epi32(
									_mm_lddqu_si128((__m128i *) target + 5),
									Match)),
					_mm_or_si128(
							_mm_cmpeq_epi32(
									_mm_lddqu_si128((__m128i *) target + 6),
									Match),
							_mm_cmpeq_epi32(
									_mm_lddqu_si128((__m128i *) target + 7),
									Match))));
	UINT4 i = 0;
	while (_mm_testz_si128(F0, F0) && ++i < V128_BLOCKSIZE / (SIMDWIDTH * 8)) {
		F0 = _mm_or_si128(
				_mm_or_si128(
						_mm_or_si128(
								_mm_cmpeq_epi32(
										_mm_lddqu_si128(
												(__m128i *) target + 8 * i + 0),
										Match),
								_mm_cmpeq_epi32(
										_mm_lddqu_si128(
												(__m128i *) target + 8 * i + 1),
										Match)),
						_mm_or_si128(
								_mm_cmpeq_epi32(
										_mm_lddqu_si128(
												(__m128i *) target + 8 * i + 2),
										Match),
								_mm_cmpeq_epi32(
										_mm_lddqu_si128(
												(__m128i *) target + 8 * i + 3),
										Match))),
				_mm_or_si128(
						_mm_or_si128(
								_mm_cmpeq_epi32(
										_mm_lddqu_si128(
												(__m128i *) target + 8 * i + 4),
										Match),
								_mm_cmpeq_epi32(
										_mm_lddqu_si128(
												(__m128i *) target + 8 * i + 5),
										Match)),
						_mm_or_si128(
								_mm_cmpeq_epi32(
										_mm_lddqu_si128(
												(__m128i *) target + 8 * i + 6),
										Match),
								_mm_cmpeq_epi32(
										_mm_lddqu_si128(
												(__m128i *) target + 8 * i + 7),
										Match))));
	}
	if (_mm_testz_si128(F0, F0)) {
		*foundp = 0;
		i = 0;
	} else
		*foundp = 1;
	return (target + i * SIMDWIDTH * 4 - init_target);
}

/* return the exact position of intersection */
/* return the head-position of mathced block */
long simdgallop_128_greedy(UINT4 goal, const UINT4 *target, long ntargets) {
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

long simdgallop_256_rough_plow(int *foundp, UINT4 goal, const UINT4 *target,
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
		F0 = _mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + i), Match);
	}
	if (_mm_testz_si128(F0, F0)) {
		*foundp = 0;
		i = 0;
	} else
		*foundp = 1;
	return (target + i * SIMDWIDTH - init_target);
}

long simdgallop_512_rough_plow(int *foundp, UINT4 goal, const UINT4 *target,
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

#endif /* INCLUDE_MSIS_GALLOP_EXACT_HPP_ */
