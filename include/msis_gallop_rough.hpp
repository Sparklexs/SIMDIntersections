/*
 * msis_gallop.hpp
 *
 *  Created on: 2017.3.28
 *      Author: John
 */

#ifndef INCLUDE_MSIS_GALLOP_ROUGH_HPP_
#define INCLUDE_MSIS_GALLOP_ROUGH_HPP_

#include "msis_linear.hpp"

namespace msis/*MultiSet InterSection*/{

long simdgallop_4_4_rough(int *foundp, UINT4 goal, const UINT4 *target,
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
		}
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
	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0 = _mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 0),
			Match);
	if (_mm_testz_si128(F0, F0))
		*foundp = 0;
	else
		*foundp = 1;
	return (target - init_target);
}
long simdgallop_8_4_rough(int *foundp, UINT4 goal, const UINT4 *target,
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
		}
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
	if (target[3] < goal) {
		target += 4;
	}
	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0 = _mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 0),
			Match);
	if (_mm_testz_si128(F0, F0))
		*foundp = 0;
	else
		*foundp = 1;
	return (target - init_target);
}
long simdgallop_8_8_rough(int *foundp, UINT4 goal, const UINT4 *target,
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
		}
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
long simdgallop_16_4_rough(int *foundp, UINT4 goal, const UINT4 *target,
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
		}
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
	if (target[7] >= goal) {
		if (target[3] < goal) {
			target += 4;
		}
	} else {
		if (target[11] >= goal) {
			target += 8;
		} else {
			target += 12;
		}
	}
	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0 = _mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 0),
			Match);
	if (_mm_testz_si128(F0, F0))
		*foundp = 0;
	else
		*foundp = 1;
	return (target - init_target);
}
long simdgallop_16_8_rough(int *foundp, UINT4 goal, const UINT4 *target,
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
		}
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
	if (target[7] < goal) {
		target += 8;
	}
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
long simdgallop_16_16_rough(int *foundp, UINT4 goal, const UINT4 *target,
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
		}
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
long simdgallop_32_4_rough(int *foundp, UINT4 goal, const UINT4 *target,
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
		}
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
	if (target[15] >= goal) {
		if (target[7] >= goal) {
			if (target[3] < goal) {
				target += 4;
			}
		} else {
			if (target[11] >= goal) {
				target += 8;
			} else {
				target += 12;
			}
		}
	} else {
		if (target[23] >= goal) {
			if (target[19] >= goal) {
				target += 16;
			} else {
				target += 20;
			}
		} else {
			if (target[27] >= goal) {
				target += 24;
			} else {
				target += 28;
			}
		}
	}
	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0 = _mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 0),
			Match);
	if (_mm_testz_si128(F0, F0))
		*foundp = 0;
	else
		*foundp = 1;
	return (target - init_target);
}
long simdgallop_32_8_rough(int *foundp, UINT4 goal, const UINT4 *target,
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
		}
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
	if (target[15] >= goal) {
		if (target[7] < goal) {
			target += 8;
		}
	} else {
		if (target[23] >= goal) {
			target += 16;
		} else {
			target += 24;
		}
	}
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
long simdgallop_32_16_rough(int *foundp, UINT4 goal, const UINT4 *target,
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
		}
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
	if (target[15] < goal) {
		target += 16;
	}
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
long simdgallop_32_32_rough(int *foundp, UINT4 goal, const UINT4 *target,
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
		}
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
	if (_mm_testz_si128(F0, F0))
		*foundp = 0;
	else
		*foundp = 1;
	return (target - init_target);
}
long simdgallop_64_4_rough(int *foundp, UINT4 goal, const UINT4 *target,
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
		}
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
	if (target[31] >= goal) {
		if (target[15] >= goal) {
			if (target[7] >= goal) {
				if (target[3] < goal) {
					target += 4;
				}
			} else {
				if (target[11] >= goal) {
					target += 8;
				} else {
					target += 12;
				}
			}
		} else {
			if (target[23] >= goal) {
				if (target[19] >= goal) {
					target += 16;
				} else {
					target += 20;
				}
			} else {
				if (target[27] >= goal) {
					target += 24;
				} else {
					target += 28;
				}
			}
		}
	} else {
		if (target[47] >= goal) {
			if (target[39] >= goal) {
				if (target[35] >= goal) {
					target += 32;
				} else {
					target += 36;
				}
			} else {
				if (target[43] >= goal) {
					target += 40;
				} else {
					target += 44;
				}
			}
		} else {
			if (target[55] >= goal) {
				if (target[51] >= goal) {
					target += 48;
				} else {
					target += 52;
				}
			} else {
				if (target[59] >= goal) {
					target += 56;
				} else {
					target += 60;
				}
			}
		}
	}
	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0 = _mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 0),
			Match);
	if (_mm_testz_si128(F0, F0))
		*foundp = 0;
	else
		*foundp = 1;
	return (target - init_target);
}
long simdgallop_64_8_rough(int *foundp, UINT4 goal, const UINT4 *target,
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
		}
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
	if (target[31] >= goal) {
		if (target[15] >= goal) {
			if (target[7] < goal) {
				target += 8;
			}
		} else {
			if (target[23] >= goal) {
				target += 16;
			} else {
				target += 24;
			}
		}
	} else {
		if (target[47] >= goal) {
			if (target[39] >= goal) {
				target += 32;
			} else {
				target += 40;
			}
		} else {
			if (target[55] >= goal) {
				target += 48;
			} else {
				target += 56;
			}
		}
	}
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
long simdgallop_64_16_rough(int *foundp, UINT4 goal, const UINT4 *target,
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
		}
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
	if (target[31] >= goal) {
		if (target[15] < goal) {
			target += 16;
		}
	} else {
		if (target[47] >= goal) {
			target += 32;
		} else {
			target += 48;
		}
	}
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
long simdgallop_64_32_rough(int *foundp, UINT4 goal, const UINT4 *target,
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
		}
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
	if (target[31] < goal) {
		target += 32;
	}
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
	if (_mm_testz_si128(F0, F0))
		*foundp = 0;
	else
		*foundp = 1;
	return (target - init_target);
}
long simdgallop_64_64_rough(int *foundp, UINT4 goal, const UINT4 *target,
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
		}
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
	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0 = _mm_or_si128(
			_mm_or_si128(
					_mm_or_si128(
							_mm_or_si128(
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 0),
											Match),
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 1),
											Match)),
							_mm_or_si128(
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 2),
											Match),
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 3),
											Match))),
					_mm_or_si128(
							_mm_or_si128(
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 4),
											Match),
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 5),
											Match)),
							_mm_or_si128(
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 6),
											Match),
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 7),
											Match)))),
			_mm_or_si128(
					_mm_or_si128(
							_mm_or_si128(
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 8),
											Match),
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 9),
											Match)),
							_mm_or_si128(
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 10),
											Match),
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 11),
											Match))),
					_mm_or_si128(
							_mm_or_si128(
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 12),
											Match),
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 13),
											Match)),
							_mm_or_si128(
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 14),
											Match),
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 15),
											Match)))));
	if (_mm_testz_si128(F0, F0))
		*foundp = 0;
	else
		*foundp = 1;
	return (target - init_target);
}
long simdgallop_128_4_rough(int *foundp, UINT4 goal, const UINT4 *target,
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
		}
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
	if (target[63] >= goal) {
		if (target[31] >= goal) {
			if (target[15] >= goal) {
				if (target[7] >= goal) {
					if (target[3] < goal) {
						target += 4;
					}
				} else {
					if (target[11] >= goal) {
						target += 8;
					} else {
						target += 12;
					}
				}
			} else {
				if (target[23] >= goal) {
					if (target[19] >= goal) {
						target += 16;
					} else {
						target += 20;
					}
				} else {
					if (target[27] >= goal) {
						target += 24;
					} else {
						target += 28;
					}
				}
			}
		} else {
			if (target[47] >= goal) {
				if (target[39] >= goal) {
					if (target[35] >= goal) {
						target += 32;
					} else {
						target += 36;
					}
				} else {
					if (target[43] >= goal) {
						target += 40;
					} else {
						target += 44;
					}
				}
			} else {
				if (target[55] >= goal) {
					if (target[51] >= goal) {
						target += 48;
					} else {
						target += 52;
					}
				} else {
					if (target[59] >= goal) {
						target += 56;
					} else {
						target += 60;
					}
				}
			}
		}
	} else {
		if (target[95] >= goal) {
			if (target[79] >= goal) {
				if (target[71] >= goal) {
					if (target[67] >= goal) {
						target += 64;
					} else {
						target += 68;
					}
				} else {
					if (target[75] >= goal) {
						target += 72;
					} else {
						target += 76;
					}
				}
			} else {
				if (target[87] >= goal) {
					if (target[83] >= goal) {
						target += 80;
					} else {
						target += 84;
					}
				} else {
					if (target[91] >= goal) {
						target += 88;
					} else {
						target += 92;
					}
				}
			}
		} else {
			if (target[111] >= goal) {
				if (target[103] >= goal) {
					if (target[99] >= goal) {
						target += 96;
					} else {
						target += 100;
					}
				} else {
					if (target[107] >= goal) {
						target += 104;
					} else {
						target += 108;
					}
				}
			} else {
				if (target[119] >= goal) {
					if (target[115] >= goal) {
						target += 112;
					} else {
						target += 116;
					}
				} else {
					if (target[123] >= goal) {
						target += 120;
					} else {
						target += 124;
					}
				}
			}
		}
	}
	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0 = _mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 0),
			Match);
	if (_mm_testz_si128(F0, F0))
		*foundp = 0;
	else
		*foundp = 1;
	return (target - init_target);
}
long simdgallop_128_8_rough(int *foundp, UINT4 goal, const UINT4 *target,
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
		}
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
	if (target[63] >= goal) {
		if (target[31] >= goal) {
			if (target[15] >= goal) {
				if (target[7] < goal) {
					target += 8;
				}
			} else {
				if (target[23] >= goal) {
					target += 16;
				} else {
					target += 24;
				}
			}
		} else {
			if (target[47] >= goal) {
				if (target[39] >= goal) {
					target += 32;
				} else {
					target += 40;
				}
			} else {
				if (target[55] >= goal) {
					target += 48;
				} else {
					target += 56;
				}
			}
		}
	} else {
		if (target[95] >= goal) {
			if (target[79] >= goal) {
				if (target[71] >= goal) {
					target += 64;
				} else {
					target += 72;
				}
			} else {
				if (target[87] >= goal) {
					target += 80;
				} else {
					target += 88;
				}
			}
		} else {
			if (target[111] >= goal) {
				if (target[103] >= goal) {
					target += 96;
				} else {
					target += 104;
				}
			} else {
				if (target[119] >= goal) {
					target += 112;
				} else {
					target += 120;
				}
			}
		}
	}
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
long simdgallop_128_16_rough(int *foundp, UINT4 goal, const UINT4 *target,
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
		}
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
	if (target[63] >= goal) {
		if (target[31] >= goal) {
			if (target[15] < goal) {
				target += 16;
			}
		} else {
			if (target[47] >= goal) {
				target += 32;
			} else {
				target += 48;
			}
		}
	} else {
		if (target[95] >= goal) {
			if (target[79] >= goal) {
				target += 64;
			} else {
				target += 80;
			}
		} else {
			if (target[111] >= goal) {
				target += 96;
			} else {
				target += 112;
			}
		}
	}
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
long simdgallop_128_32_rough(int *foundp, UINT4 goal, const UINT4 *target,
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
		}
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
	if (target[63] >= goal) {
		if (target[31] < goal) {
			target += 32;
		}
	} else {
		if (target[95] >= goal) {
			target += 64;
		} else {
			target += 96;
		}
	}
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
	if (_mm_testz_si128(F0, F0))
		*foundp = 0;
	else
		*foundp = 1;
	return (target - init_target);
}
long simdgallop_128_64_rough(int *foundp, UINT4 goal, const UINT4 *target,
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
		}
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
	if (target[63] < goal) {
		target += 64;
	}
	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0 = _mm_or_si128(
			_mm_or_si128(
					_mm_or_si128(
							_mm_or_si128(
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 0),
											Match),
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 1),
											Match)),
							_mm_or_si128(
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 2),
											Match),
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 3),
											Match))),
					_mm_or_si128(
							_mm_or_si128(
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 4),
											Match),
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 5),
											Match)),
							_mm_or_si128(
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 6),
											Match),
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 7),
											Match)))),
			_mm_or_si128(
					_mm_or_si128(
							_mm_or_si128(
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 8),
											Match),
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 9),
											Match)),
							_mm_or_si128(
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 10),
											Match),
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 11),
											Match))),
					_mm_or_si128(
							_mm_or_si128(
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 12),
											Match),
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 13),
											Match)),
							_mm_or_si128(
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 14),
											Match),
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 15),
											Match)))));
	if (_mm_testz_si128(F0, F0))
		*foundp = 0;
	else
		*foundp = 1;
	return (target - init_target);
}
long simdgallop_128_128_rough(int *foundp, UINT4 goal, const UINT4 *target,
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
		}
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
	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0 = _mm_or_si128(
			_mm_or_si128(
					_mm_or_si128(
							_mm_or_si128(
									_mm_or_si128(
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 0),
													Match),
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 1),
													Match)),
									_mm_or_si128(
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 2),
													Match),
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 3),
													Match))),
							_mm_or_si128(
									_mm_or_si128(
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 4),
													Match),
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 5),
													Match)),
									_mm_or_si128(
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 6),
													Match),
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 7),
													Match)))),
					_mm_or_si128(
							_mm_or_si128(
									_mm_or_si128(
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 8),
													Match),
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 9),
													Match)),
									_mm_or_si128(
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 10),
													Match),
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 11),
													Match))),
							_mm_or_si128(
									_mm_or_si128(
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 12),
													Match),
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 13),
													Match)),
									_mm_or_si128(
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 14),
													Match),
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 15),
													Match))))),
			_mm_or_si128(
					_mm_or_si128(
							_mm_or_si128(
									_mm_or_si128(
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 16),
													Match),
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 17),
													Match)),
									_mm_or_si128(
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 18),
													Match),
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 19),
													Match))),
							_mm_or_si128(
									_mm_or_si128(
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 20),
													Match),
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 21),
													Match)),
									_mm_or_si128(
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 22),
													Match),
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 23),
													Match)))),
					_mm_or_si128(
							_mm_or_si128(
									_mm_or_si128(
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 24),
													Match),
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 25),
													Match)),
									_mm_or_si128(
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 26),
													Match),
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 27),
													Match))),
							_mm_or_si128(
									_mm_or_si128(
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 28),
													Match),
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 29),
													Match)),
									_mm_or_si128(
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 30),
													Match),
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 31),
													Match))))));
	if (_mm_testz_si128(F0, F0))
		*foundp = 0;
	else
		*foundp = 1;
	return (target - init_target);
}
long simdgallop_256_4_rough(int *foundp, UINT4 goal, const UINT4 *target,
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
		}
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
	if (target[127] >= goal) {
		if (target[63] >= goal) {
			if (target[31] >= goal) {
				if (target[15] >= goal) {
					if (target[7] >= goal) {
						if (target[3] < goal) {
							target += 4;
						}
					} else {
						if (target[11] >= goal) {
							target += 8;
						} else {
							target += 12;
						}
					}
				} else {
					if (target[23] >= goal) {
						if (target[19] >= goal) {
							target += 16;
						} else {
							target += 20;
						}
					} else {
						if (target[27] >= goal) {
							target += 24;
						} else {
							target += 28;
						}
					}
				}
			} else {
				if (target[47] >= goal) {
					if (target[39] >= goal) {
						if (target[35] >= goal) {
							target += 32;
						} else {
							target += 36;
						}
					} else {
						if (target[43] >= goal) {
							target += 40;
						} else {
							target += 44;
						}
					}
				} else {
					if (target[55] >= goal) {
						if (target[51] >= goal) {
							target += 48;
						} else {
							target += 52;
						}
					} else {
						if (target[59] >= goal) {
							target += 56;
						} else {
							target += 60;
						}
					}
				}
			}
		} else {
			if (target[95] >= goal) {
				if (target[79] >= goal) {
					if (target[71] >= goal) {
						if (target[67] >= goal) {
							target += 64;
						} else {
							target += 68;
						}
					} else {
						if (target[75] >= goal) {
							target += 72;
						} else {
							target += 76;
						}
					}
				} else {
					if (target[87] >= goal) {
						if (target[83] >= goal) {
							target += 80;
						} else {
							target += 84;
						}
					} else {
						if (target[91] >= goal) {
							target += 88;
						} else {
							target += 92;
						}
					}
				}
			} else {
				if (target[111] >= goal) {
					if (target[103] >= goal) {
						if (target[99] >= goal) {
							target += 96;
						} else {
							target += 100;
						}
					} else {
						if (target[107] >= goal) {
							target += 104;
						} else {
							target += 108;
						}
					}
				} else {
					if (target[119] >= goal) {
						if (target[115] >= goal) {
							target += 112;
						} else {
							target += 116;
						}
					} else {
						if (target[123] >= goal) {
							target += 120;
						} else {
							target += 124;
						}
					}
				}
			}
		}
	} else {
		if (target[191] >= goal) {
			if (target[159] >= goal) {
				if (target[143] >= goal) {
					if (target[135] >= goal) {
						if (target[131] >= goal) {
							target += 128;
						} else {
							target += 132;
						}
					} else {
						if (target[139] >= goal) {
							target += 136;
						} else {
							target += 140;
						}
					}
				} else {
					if (target[151] >= goal) {
						if (target[147] >= goal) {
							target += 144;
						} else {
							target += 148;
						}
					} else {
						if (target[155] >= goal) {
							target += 152;
						} else {
							target += 156;
						}
					}
				}
			} else {
				if (target[175] >= goal) {
					if (target[167] >= goal) {
						if (target[163] >= goal) {
							target += 160;
						} else {
							target += 164;
						}
					} else {
						if (target[171] >= goal) {
							target += 168;
						} else {
							target += 172;
						}
					}
				} else {
					if (target[183] >= goal) {
						if (target[179] >= goal) {
							target += 176;
						} else {
							target += 180;
						}
					} else {
						if (target[187] >= goal) {
							target += 184;
						} else {
							target += 188;
						}
					}
				}
			}
		} else {
			if (target[223] >= goal) {
				if (target[207] >= goal) {
					if (target[199] >= goal) {
						if (target[195] >= goal) {
							target += 192;
						} else {
							target += 196;
						}
					} else {
						if (target[203] >= goal) {
							target += 200;
						} else {
							target += 204;
						}
					}
				} else {
					if (target[215] >= goal) {
						if (target[211] >= goal) {
							target += 208;
						} else {
							target += 212;
						}
					} else {
						if (target[219] >= goal) {
							target += 216;
						} else {
							target += 220;
						}
					}
				}
			} else {
				if (target[239] >= goal) {
					if (target[231] >= goal) {
						if (target[227] >= goal) {
							target += 224;
						} else {
							target += 228;
						}
					} else {
						if (target[235] >= goal) {
							target += 232;
						} else {
							target += 236;
						}
					}
				} else {
					if (target[247] >= goal) {
						if (target[243] >= goal) {
							target += 240;
						} else {
							target += 244;
						}
					} else {
						if (target[251] >= goal) {
							target += 248;
						} else {
							target += 252;
						}
					}
				}
			}
		}
	}
	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0 = _mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 0),
			Match);
	if (_mm_testz_si128(F0, F0))
		*foundp = 0;
	else
		*foundp = 1;
	return (target - init_target);
}
long simdgallop_256_8_rough(int *foundp, UINT4 goal, const UINT4 *target,
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
		}
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
	if (target[127] >= goal) {
		if (target[63] >= goal) {
			if (target[31] >= goal) {
				if (target[15] >= goal) {
					if (target[7] < goal) {
						target += 8;
					}
				} else {
					if (target[23] >= goal) {
						target += 16;
					} else {
						target += 24;
					}
				}
			} else {
				if (target[47] >= goal) {
					if (target[39] >= goal) {
						target += 32;
					} else {
						target += 40;
					}
				} else {
					if (target[55] >= goal) {
						target += 48;
					} else {
						target += 56;
					}
				}
			}
		} else {
			if (target[95] >= goal) {
				if (target[79] >= goal) {
					if (target[71] >= goal) {
						target += 64;
					} else {
						target += 72;
					}
				} else {
					if (target[87] >= goal) {
						target += 80;
					} else {
						target += 88;
					}
				}
			} else {
				if (target[111] >= goal) {
					if (target[103] >= goal) {
						target += 96;
					} else {
						target += 104;
					}
				} else {
					if (target[119] >= goal) {
						target += 112;
					} else {
						target += 120;
					}
				}
			}
		}
	} else {
		if (target[191] >= goal) {
			if (target[159] >= goal) {
				if (target[143] >= goal) {
					if (target[135] >= goal) {
						target += 128;
					} else {
						target += 136;
					}
				} else {
					if (target[151] >= goal) {
						target += 144;
					} else {
						target += 152;
					}
				}
			} else {
				if (target[175] >= goal) {
					if (target[167] >= goal) {
						target += 160;
					} else {
						target += 168;
					}
				} else {
					if (target[183] >= goal) {
						target += 176;
					} else {
						target += 184;
					}
				}
			}
		} else {
			if (target[223] >= goal) {
				if (target[207] >= goal) {
					if (target[199] >= goal) {
						target += 192;
					} else {
						target += 200;
					}
				} else {
					if (target[215] >= goal) {
						target += 208;
					} else {
						target += 216;
					}
				}
			} else {
				if (target[239] >= goal) {
					if (target[231] >= goal) {
						target += 224;
					} else {
						target += 232;
					}
				} else {
					if (target[247] >= goal) {
						target += 240;
					} else {
						target += 248;
					}
				}
			}
		}
	}
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
long simdgallop_256_16_rough(int *foundp, UINT4 goal, const UINT4 *target,
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
		}
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
	if (target[127] >= goal) {
		if (target[63] >= goal) {
			if (target[31] >= goal) {
				if (target[15] < goal) {
					target += 16;
				}
			} else {
				if (target[47] >= goal) {
					target += 32;
				} else {
					target += 48;
				}
			}
		} else {
			if (target[95] >= goal) {
				if (target[79] >= goal) {
					target += 64;
				} else {
					target += 80;
				}
			} else {
				if (target[111] >= goal) {
					target += 96;
				} else {
					target += 112;
				}
			}
		}
	} else {
		if (target[191] >= goal) {
			if (target[159] >= goal) {
				if (target[143] >= goal) {
					target += 128;
				} else {
					target += 144;
				}
			} else {
				if (target[175] >= goal) {
					target += 160;
				} else {
					target += 176;
				}
			}
		} else {
			if (target[223] >= goal) {
				if (target[207] >= goal) {
					target += 192;
				} else {
					target += 208;
				}
			} else {
				if (target[239] >= goal) {
					target += 224;
				} else {
					target += 240;
				}
			}
		}
	}
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
long simdgallop_256_32_rough(int *foundp, UINT4 goal, const UINT4 *target,
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
		}
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
	if (target[127] >= goal) {
		if (target[63] >= goal) {
			if (target[31] < goal) {
				target += 32;
			}
		} else {
			if (target[95] >= goal) {
				target += 64;
			} else {
				target += 96;
			}
		}
	} else {
		if (target[191] >= goal) {
			if (target[159] >= goal) {
				target += 128;
			} else {
				target += 160;
			}
		} else {
			if (target[223] >= goal) {
				target += 192;
			} else {
				target += 224;
			}
		}
	}
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
	if (_mm_testz_si128(F0, F0))
		*foundp = 0;
	else
		*foundp = 1;
	return (target - init_target);
}
long simdgallop_256_64_rough(int *foundp, UINT4 goal, const UINT4 *target,
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
		}
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
	if (target[127] >= goal) {
		if (target[63] < goal) {
			target += 64;
		}
	} else {
		if (target[191] >= goal) {
			target += 128;
		} else {
			target += 192;
		}
	}
	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0 = _mm_or_si128(
			_mm_or_si128(
					_mm_or_si128(
							_mm_or_si128(
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 0),
											Match),
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 1),
											Match)),
							_mm_or_si128(
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 2),
											Match),
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 3),
											Match))),
					_mm_or_si128(
							_mm_or_si128(
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 4),
											Match),
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 5),
											Match)),
							_mm_or_si128(
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 6),
											Match),
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 7),
											Match)))),
			_mm_or_si128(
					_mm_or_si128(
							_mm_or_si128(
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 8),
											Match),
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 9),
											Match)),
							_mm_or_si128(
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 10),
											Match),
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 11),
											Match))),
					_mm_or_si128(
							_mm_or_si128(
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 12),
											Match),
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 13),
											Match)),
							_mm_or_si128(
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 14),
											Match),
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 15),
											Match)))));
	if (_mm_testz_si128(F0, F0))
		*foundp = 0;
	else
		*foundp = 1;
	return (target - init_target);
}
long simdgallop_256_128_rough(int *foundp, UINT4 goal, const UINT4 *target,
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
		}
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
	if (target[127] < goal) {
		target += 128;
	}
	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0 = _mm_or_si128(
			_mm_or_si128(
					_mm_or_si128(
							_mm_or_si128(
									_mm_or_si128(
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 0),
													Match),
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 1),
													Match)),
									_mm_or_si128(
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 2),
													Match),
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 3),
													Match))),
							_mm_or_si128(
									_mm_or_si128(
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 4),
													Match),
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 5),
													Match)),
									_mm_or_si128(
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 6),
													Match),
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 7),
													Match)))),
					_mm_or_si128(
							_mm_or_si128(
									_mm_or_si128(
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 8),
													Match),
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 9),
													Match)),
									_mm_or_si128(
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 10),
													Match),
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 11),
													Match))),
							_mm_or_si128(
									_mm_or_si128(
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 12),
													Match),
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 13),
													Match)),
									_mm_or_si128(
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 14),
													Match),
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 15),
													Match))))),
			_mm_or_si128(
					_mm_or_si128(
							_mm_or_si128(
									_mm_or_si128(
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 16),
													Match),
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 17),
													Match)),
									_mm_or_si128(
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 18),
													Match),
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 19),
													Match))),
							_mm_or_si128(
									_mm_or_si128(
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 20),
													Match),
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 21),
													Match)),
									_mm_or_si128(
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 22),
													Match),
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 23),
													Match)))),
					_mm_or_si128(
							_mm_or_si128(
									_mm_or_si128(
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 24),
													Match),
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 25),
													Match)),
									_mm_or_si128(
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 26),
													Match),
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 27),
													Match))),
							_mm_or_si128(
									_mm_or_si128(
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 28),
													Match),
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 29),
													Match)),
									_mm_or_si128(
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 30),
													Match),
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 31),
													Match))))));
	if (_mm_testz_si128(F0, F0))
		*foundp = 0;
	else
		*foundp = 1;
	return (target - init_target);
}
long simdgallop_256_256_rough(int *foundp, UINT4 goal, const UINT4 *target,
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
		}
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
	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0 =
			_mm_or_si128(
					_mm_or_si128(
							_mm_or_si128(
									_mm_or_si128(
											_mm_or_si128(
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 0),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 1),
																	Match)),
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 2),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 3),
																	Match))),
											_mm_or_si128(
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 4),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 5),
																	Match)),
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 6),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 7),
																	Match)))),
									_mm_or_si128(
											_mm_or_si128(
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 8),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 9),
																	Match)),
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 10),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 11),
																	Match))),
											_mm_or_si128(
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 12),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 13),
																	Match)),
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 14),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 15),
																	Match))))),
							_mm_or_si128(
									_mm_or_si128(
											_mm_or_si128(
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 16),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 17),
																	Match)),
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 18),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 19),
																	Match))),
											_mm_or_si128(
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 20),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 21),
																	Match)),
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 22),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 23),
																	Match)))),
									_mm_or_si128(
											_mm_or_si128(
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 24),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 25),
																	Match)),
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 26),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 27),
																	Match))),
											_mm_or_si128(
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 28),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 29),
																	Match)),
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 30),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 31),
																	Match)))))),
					_mm_or_si128(
							_mm_or_si128(
									_mm_or_si128(
											_mm_or_si128(
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 32),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 33),
																	Match)),
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 34),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 35),
																	Match))),
											_mm_or_si128(
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 36),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 37),
																	Match)),
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 38),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 39),
																	Match)))),
									_mm_or_si128(
											_mm_or_si128(
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 40),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 41),
																	Match)),
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 42),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 43),
																	Match))),
											_mm_or_si128(
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 44),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 45),
																	Match)),
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 46),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 47),
																	Match))))),
							_mm_or_si128(
									_mm_or_si128(
											_mm_or_si128(
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 48),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 49),
																	Match)),
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 50),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 51),
																	Match))),
											_mm_or_si128(
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 52),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 53),
																	Match)),
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 54),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 55),
																	Match)))),
									_mm_or_si128(
											_mm_or_si128(
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 56),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 57),
																	Match)),
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 58),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 59),
																	Match))),
											_mm_or_si128(
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 60),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 61),
																	Match)),
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 62),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 63),
																	Match)))))));
	if (_mm_testz_si128(F0, F0))
		*foundp = 0;
	else
		*foundp = 1;
	return (target - init_target);
}
long simdgallop_512_4_rough(int *foundp, UINT4 goal, const UINT4 *target,
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
		}
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
	if (target[255] >= goal) {
		if (target[127] >= goal) {
			if (target[63] >= goal) {
				if (target[31] >= goal) {
					if (target[15] >= goal) {
						if (target[7] >= goal) {
							if (target[3] < goal) {
								target += 4;
							}
						} else {
							if (target[11] >= goal) {
								target += 8;
							} else {
								target += 12;
							}
						}
					} else {
						if (target[23] >= goal) {
							if (target[19] >= goal) {
								target += 16;
							} else {
								target += 20;
							}
						} else {
							if (target[27] >= goal) {
								target += 24;
							} else {
								target += 28;
							}
						}
					}
				} else {
					if (target[47] >= goal) {
						if (target[39] >= goal) {
							if (target[35] >= goal) {
								target += 32;
							} else {
								target += 36;
							}
						} else {
							if (target[43] >= goal) {
								target += 40;
							} else {
								target += 44;
							}
						}
					} else {
						if (target[55] >= goal) {
							if (target[51] >= goal) {
								target += 48;
							} else {
								target += 52;
							}
						} else {
							if (target[59] >= goal) {
								target += 56;
							} else {
								target += 60;
							}
						}
					}
				}
			} else {
				if (target[95] >= goal) {
					if (target[79] >= goal) {
						if (target[71] >= goal) {
							if (target[67] >= goal) {
								target += 64;
							} else {
								target += 68;
							}
						} else {
							if (target[75] >= goal) {
								target += 72;
							} else {
								target += 76;
							}
						}
					} else {
						if (target[87] >= goal) {
							if (target[83] >= goal) {
								target += 80;
							} else {
								target += 84;
							}
						} else {
							if (target[91] >= goal) {
								target += 88;
							} else {
								target += 92;
							}
						}
					}
				} else {
					if (target[111] >= goal) {
						if (target[103] >= goal) {
							if (target[99] >= goal) {
								target += 96;
							} else {
								target += 100;
							}
						} else {
							if (target[107] >= goal) {
								target += 104;
							} else {
								target += 108;
							}
						}
					} else {
						if (target[119] >= goal) {
							if (target[115] >= goal) {
								target += 112;
							} else {
								target += 116;
							}
						} else {
							if (target[123] >= goal) {
								target += 120;
							} else {
								target += 124;
							}
						}
					}
				}
			}
		} else {
			if (target[191] >= goal) {
				if (target[159] >= goal) {
					if (target[143] >= goal) {
						if (target[135] >= goal) {
							if (target[131] >= goal) {
								target += 128;
							} else {
								target += 132;
							}
						} else {
							if (target[139] >= goal) {
								target += 136;
							} else {
								target += 140;
							}
						}
					} else {
						if (target[151] >= goal) {
							if (target[147] >= goal) {
								target += 144;
							} else {
								target += 148;
							}
						} else {
							if (target[155] >= goal) {
								target += 152;
							} else {
								target += 156;
							}
						}
					}
				} else {
					if (target[175] >= goal) {
						if (target[167] >= goal) {
							if (target[163] >= goal) {
								target += 160;
							} else {
								target += 164;
							}
						} else {
							if (target[171] >= goal) {
								target += 168;
							} else {
								target += 172;
							}
						}
					} else {
						if (target[183] >= goal) {
							if (target[179] >= goal) {
								target += 176;
							} else {
								target += 180;
							}
						} else {
							if (target[187] >= goal) {
								target += 184;
							} else {
								target += 188;
							}
						}
					}
				}
			} else {
				if (target[223] >= goal) {
					if (target[207] >= goal) {
						if (target[199] >= goal) {
							if (target[195] >= goal) {
								target += 192;
							} else {
								target += 196;
							}
						} else {
							if (target[203] >= goal) {
								target += 200;
							} else {
								target += 204;
							}
						}
					} else {
						if (target[215] >= goal) {
							if (target[211] >= goal) {
								target += 208;
							} else {
								target += 212;
							}
						} else {
							if (target[219] >= goal) {
								target += 216;
							} else {
								target += 220;
							}
						}
					}
				} else {
					if (target[239] >= goal) {
						if (target[231] >= goal) {
							if (target[227] >= goal) {
								target += 224;
							} else {
								target += 228;
							}
						} else {
							if (target[235] >= goal) {
								target += 232;
							} else {
								target += 236;
							}
						}
					} else {
						if (target[247] >= goal) {
							if (target[243] >= goal) {
								target += 240;
							} else {
								target += 244;
							}
						} else {
							if (target[251] >= goal) {
								target += 248;
							} else {
								target += 252;
							}
						}
					}
				}
			}
		}
	} else {
		if (target[383] >= goal) {
			if (target[319] >= goal) {
				if (target[287] >= goal) {
					if (target[271] >= goal) {
						if (target[263] >= goal) {
							if (target[259] >= goal) {
								target += 256;
							} else {
								target += 260;
							}
						} else {
							if (target[267] >= goal) {
								target += 264;
							} else {
								target += 268;
							}
						}
					} else {
						if (target[279] >= goal) {
							if (target[275] >= goal) {
								target += 272;
							} else {
								target += 276;
							}
						} else {
							if (target[283] >= goal) {
								target += 280;
							} else {
								target += 284;
							}
						}
					}
				} else {
					if (target[303] >= goal) {
						if (target[295] >= goal) {
							if (target[291] >= goal) {
								target += 288;
							} else {
								target += 292;
							}
						} else {
							if (target[299] >= goal) {
								target += 296;
							} else {
								target += 300;
							}
						}
					} else {
						if (target[311] >= goal) {
							if (target[307] >= goal) {
								target += 304;
							} else {
								target += 308;
							}
						} else {
							if (target[315] >= goal) {
								target += 312;
							} else {
								target += 316;
							}
						}
					}
				}
			} else {
				if (target[351] >= goal) {
					if (target[335] >= goal) {
						if (target[327] >= goal) {
							if (target[323] >= goal) {
								target += 320;
							} else {
								target += 324;
							}
						} else {
							if (target[331] >= goal) {
								target += 328;
							} else {
								target += 332;
							}
						}
					} else {
						if (target[343] >= goal) {
							if (target[339] >= goal) {
								target += 336;
							} else {
								target += 340;
							}
						} else {
							if (target[347] >= goal) {
								target += 344;
							} else {
								target += 348;
							}
						}
					}
				} else {
					if (target[367] >= goal) {
						if (target[359] >= goal) {
							if (target[355] >= goal) {
								target += 352;
							} else {
								target += 356;
							}
						} else {
							if (target[363] >= goal) {
								target += 360;
							} else {
								target += 364;
							}
						}
					} else {
						if (target[375] >= goal) {
							if (target[371] >= goal) {
								target += 368;
							} else {
								target += 372;
							}
						} else {
							if (target[379] >= goal) {
								target += 376;
							} else {
								target += 380;
							}
						}
					}
				}
			}
		} else {
			if (target[447] >= goal) {
				if (target[415] >= goal) {
					if (target[399] >= goal) {
						if (target[391] >= goal) {
							if (target[387] >= goal) {
								target += 384;
							} else {
								target += 388;
							}
						} else {
							if (target[395] >= goal) {
								target += 392;
							} else {
								target += 396;
							}
						}
					} else {
						if (target[407] >= goal) {
							if (target[403] >= goal) {
								target += 400;
							} else {
								target += 404;
							}
						} else {
							if (target[411] >= goal) {
								target += 408;
							} else {
								target += 412;
							}
						}
					}
				} else {
					if (target[431] >= goal) {
						if (target[423] >= goal) {
							if (target[419] >= goal) {
								target += 416;
							} else {
								target += 420;
							}
						} else {
							if (target[427] >= goal) {
								target += 424;
							} else {
								target += 428;
							}
						}
					} else {
						if (target[439] >= goal) {
							if (target[435] >= goal) {
								target += 432;
							} else {
								target += 436;
							}
						} else {
							if (target[443] >= goal) {
								target += 440;
							} else {
								target += 444;
							}
						}
					}
				}
			} else {
				if (target[479] >= goal) {
					if (target[463] >= goal) {
						if (target[455] >= goal) {
							if (target[451] >= goal) {
								target += 448;
							} else {
								target += 452;
							}
						} else {
							if (target[459] >= goal) {
								target += 456;
							} else {
								target += 460;
							}
						}
					} else {
						if (target[471] >= goal) {
							if (target[467] >= goal) {
								target += 464;
							} else {
								target += 468;
							}
						} else {
							if (target[475] >= goal) {
								target += 472;
							} else {
								target += 476;
							}
						}
					}
				} else {
					if (target[495] >= goal) {
						if (target[487] >= goal) {
							if (target[483] >= goal) {
								target += 480;
							} else {
								target += 484;
							}
						} else {
							if (target[491] >= goal) {
								target += 488;
							} else {
								target += 492;
							}
						}
					} else {
						if (target[503] >= goal) {
							if (target[499] >= goal) {
								target += 496;
							} else {
								target += 500;
							}
						} else {
							if (target[507] >= goal) {
								target += 504;
							} else {
								target += 508;
							}
						}
					}
				}
			}
		}
	}
	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0 = _mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + 0),
			Match);
	if (_mm_testz_si128(F0, F0))
		*foundp = 0;
	else
		*foundp = 1;
	return (target - init_target);
}
long simdgallop_512_8_rough(int *foundp, UINT4 goal, const UINT4 *target,
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
		}
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
	if (target[255] >= goal) {
		if (target[127] >= goal) {
			if (target[63] >= goal) {
				if (target[31] >= goal) {
					if (target[15] >= goal) {
						if (target[7] < goal) {
							target += 8;
						}
					} else {
						if (target[23] >= goal) {
							target += 16;
						} else {
							target += 24;
						}
					}
				} else {
					if (target[47] >= goal) {
						if (target[39] >= goal) {
							target += 32;
						} else {
							target += 40;
						}
					} else {
						if (target[55] >= goal) {
							target += 48;
						} else {
							target += 56;
						}
					}
				}
			} else {
				if (target[95] >= goal) {
					if (target[79] >= goal) {
						if (target[71] >= goal) {
							target += 64;
						} else {
							target += 72;
						}
					} else {
						if (target[87] >= goal) {
							target += 80;
						} else {
							target += 88;
						}
					}
				} else {
					if (target[111] >= goal) {
						if (target[103] >= goal) {
							target += 96;
						} else {
							target += 104;
						}
					} else {
						if (target[119] >= goal) {
							target += 112;
						} else {
							target += 120;
						}
					}
				}
			}
		} else {
			if (target[191] >= goal) {
				if (target[159] >= goal) {
					if (target[143] >= goal) {
						if (target[135] >= goal) {
							target += 128;
						} else {
							target += 136;
						}
					} else {
						if (target[151] >= goal) {
							target += 144;
						} else {
							target += 152;
						}
					}
				} else {
					if (target[175] >= goal) {
						if (target[167] >= goal) {
							target += 160;
						} else {
							target += 168;
						}
					} else {
						if (target[183] >= goal) {
							target += 176;
						} else {
							target += 184;
						}
					}
				}
			} else {
				if (target[223] >= goal) {
					if (target[207] >= goal) {
						if (target[199] >= goal) {
							target += 192;
						} else {
							target += 200;
						}
					} else {
						if (target[215] >= goal) {
							target += 208;
						} else {
							target += 216;
						}
					}
				} else {
					if (target[239] >= goal) {
						if (target[231] >= goal) {
							target += 224;
						} else {
							target += 232;
						}
					} else {
						if (target[247] >= goal) {
							target += 240;
						} else {
							target += 248;
						}
					}
				}
			}
		}
	} else {
		if (target[383] >= goal) {
			if (target[319] >= goal) {
				if (target[287] >= goal) {
					if (target[271] >= goal) {
						if (target[263] >= goal) {
							target += 256;
						} else {
							target += 264;
						}
					} else {
						if (target[279] >= goal) {
							target += 272;
						} else {
							target += 280;
						}
					}
				} else {
					if (target[303] >= goal) {
						if (target[295] >= goal) {
							target += 288;
						} else {
							target += 296;
						}
					} else {
						if (target[311] >= goal) {
							target += 304;
						} else {
							target += 312;
						}
					}
				}
			} else {
				if (target[351] >= goal) {
					if (target[335] >= goal) {
						if (target[327] >= goal) {
							target += 320;
						} else {
							target += 328;
						}
					} else {
						if (target[343] >= goal) {
							target += 336;
						} else {
							target += 344;
						}
					}
				} else {
					if (target[367] >= goal) {
						if (target[359] >= goal) {
							target += 352;
						} else {
							target += 360;
						}
					} else {
						if (target[375] >= goal) {
							target += 368;
						} else {
							target += 376;
						}
					}
				}
			}
		} else {
			if (target[447] >= goal) {
				if (target[415] >= goal) {
					if (target[399] >= goal) {
						if (target[391] >= goal) {
							target += 384;
						} else {
							target += 392;
						}
					} else {
						if (target[407] >= goal) {
							target += 400;
						} else {
							target += 408;
						}
					}
				} else {
					if (target[431] >= goal) {
						if (target[423] >= goal) {
							target += 416;
						} else {
							target += 424;
						}
					} else {
						if (target[439] >= goal) {
							target += 432;
						} else {
							target += 440;
						}
					}
				}
			} else {
				if (target[479] >= goal) {
					if (target[463] >= goal) {
						if (target[455] >= goal) {
							target += 448;
						} else {
							target += 456;
						}
					} else {
						if (target[471] >= goal) {
							target += 464;
						} else {
							target += 472;
						}
					}
				} else {
					if (target[495] >= goal) {
						if (target[487] >= goal) {
							target += 480;
						} else {
							target += 488;
						}
					} else {
						if (target[503] >= goal) {
							target += 496;
						} else {
							target += 504;
						}
					}
				}
			}
		}
	}
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
long simdgallop_512_16_rough(int *foundp, UINT4 goal, const UINT4 *target,
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
		}
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
	if (target[255] >= goal) {
		if (target[127] >= goal) {
			if (target[63] >= goal) {
				if (target[31] >= goal) {
					if (target[15] < goal) {
						target += 16;
					}
				} else {
					if (target[47] >= goal) {
						target += 32;
					} else {
						target += 48;
					}
				}
			} else {
				if (target[95] >= goal) {
					if (target[79] >= goal) {
						target += 64;
					} else {
						target += 80;
					}
				} else {
					if (target[111] >= goal) {
						target += 96;
					} else {
						target += 112;
					}
				}
			}
		} else {
			if (target[191] >= goal) {
				if (target[159] >= goal) {
					if (target[143] >= goal) {
						target += 128;
					} else {
						target += 144;
					}
				} else {
					if (target[175] >= goal) {
						target += 160;
					} else {
						target += 176;
					}
				}
			} else {
				if (target[223] >= goal) {
					if (target[207] >= goal) {
						target += 192;
					} else {
						target += 208;
					}
				} else {
					if (target[239] >= goal) {
						target += 224;
					} else {
						target += 240;
					}
				}
			}
		}
	} else {
		if (target[383] >= goal) {
			if (target[319] >= goal) {
				if (target[287] >= goal) {
					if (target[271] >= goal) {
						target += 256;
					} else {
						target += 272;
					}
				} else {
					if (target[303] >= goal) {
						target += 288;
					} else {
						target += 304;
					}
				}
			} else {
				if (target[351] >= goal) {
					if (target[335] >= goal) {
						target += 320;
					} else {
						target += 336;
					}
				} else {
					if (target[367] >= goal) {
						target += 352;
					} else {
						target += 368;
					}
				}
			}
		} else {
			if (target[447] >= goal) {
				if (target[415] >= goal) {
					if (target[399] >= goal) {
						target += 384;
					} else {
						target += 400;
					}
				} else {
					if (target[431] >= goal) {
						target += 416;
					} else {
						target += 432;
					}
				}
			} else {
				if (target[479] >= goal) {
					if (target[463] >= goal) {
						target += 448;
					} else {
						target += 464;
					}
				} else {
					if (target[495] >= goal) {
						target += 480;
					} else {
						target += 496;
					}
				}
			}
		}
	}
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
long simdgallop_512_32_rough(int *foundp, UINT4 goal, const UINT4 *target,
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
		}
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
	if (target[255] >= goal) {
		if (target[127] >= goal) {
			if (target[63] >= goal) {
				if (target[31] < goal) {
					target += 32;
				}
			} else {
				if (target[95] >= goal) {
					target += 64;
				} else {
					target += 96;
				}
			}
		} else {
			if (target[191] >= goal) {
				if (target[159] >= goal) {
					target += 128;
				} else {
					target += 160;
				}
			} else {
				if (target[223] >= goal) {
					target += 192;
				} else {
					target += 224;
				}
			}
		}
	} else {
		if (target[383] >= goal) {
			if (target[319] >= goal) {
				if (target[287] >= goal) {
					target += 256;
				} else {
					target += 288;
				}
			} else {
				if (target[351] >= goal) {
					target += 320;
				} else {
					target += 352;
				}
			}
		} else {
			if (target[447] >= goal) {
				if (target[415] >= goal) {
					target += 384;
				} else {
					target += 416;
				}
			} else {
				if (target[479] >= goal) {
					target += 448;
				} else {
					target += 480;
				}
			}
		}
	}
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
	if (_mm_testz_si128(F0, F0))
		*foundp = 0;
	else
		*foundp = 1;
	return (target - init_target);
}
long simdgallop_512_64_rough(int *foundp, UINT4 goal, const UINT4 *target,
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
		}
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
	if (target[255] >= goal) {
		if (target[127] >= goal) {
			if (target[63] < goal) {
				target += 64;
			}
		} else {
			if (target[191] >= goal) {
				target += 128;
			} else {
				target += 192;
			}
		}
	} else {
		if (target[383] >= goal) {
			if (target[319] >= goal) {
				target += 256;
			} else {
				target += 320;
			}
		} else {
			if (target[447] >= goal) {
				target += 384;
			} else {
				target += 448;
			}
		}
	}
	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0 = _mm_or_si128(
			_mm_or_si128(
					_mm_or_si128(
							_mm_or_si128(
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 0),
											Match),
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 1),
											Match)),
							_mm_or_si128(
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 2),
											Match),
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 3),
											Match))),
					_mm_or_si128(
							_mm_or_si128(
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 4),
											Match),
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 5),
											Match)),
							_mm_or_si128(
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 6),
											Match),
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 7),
											Match)))),
			_mm_or_si128(
					_mm_or_si128(
							_mm_or_si128(
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 8),
											Match),
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 9),
											Match)),
							_mm_or_si128(
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 10),
											Match),
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 11),
											Match))),
					_mm_or_si128(
							_mm_or_si128(
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 12),
											Match),
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 13),
											Match)),
							_mm_or_si128(
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 14),
											Match),
									_mm_cmpeq_epi32(
											_mm_lddqu_si128(
													(__m128i *) target + 15),
											Match)))));
	if (_mm_testz_si128(F0, F0))
		*foundp = 0;
	else
		*foundp = 1;
	return (target - init_target);
}
long simdgallop_512_128_rough(int *foundp, UINT4 goal, const UINT4 *target,
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
		}
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
	if (target[255] >= goal) {
		if (target[127] < goal) {
			target += 128;
		}
	} else {
		if (target[383] >= goal) {
			target += 256;
		} else {
			target += 384;
		}
	}
	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0 = _mm_or_si128(
			_mm_or_si128(
					_mm_or_si128(
							_mm_or_si128(
									_mm_or_si128(
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 0),
													Match),
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 1),
													Match)),
									_mm_or_si128(
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 2),
													Match),
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 3),
													Match))),
							_mm_or_si128(
									_mm_or_si128(
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 4),
													Match),
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 5),
													Match)),
									_mm_or_si128(
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 6),
													Match),
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 7),
													Match)))),
					_mm_or_si128(
							_mm_or_si128(
									_mm_or_si128(
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 8),
													Match),
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 9),
													Match)),
									_mm_or_si128(
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 10),
													Match),
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 11),
													Match))),
							_mm_or_si128(
									_mm_or_si128(
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 12),
													Match),
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 13),
													Match)),
									_mm_or_si128(
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 14),
													Match),
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 15),
													Match))))),
			_mm_or_si128(
					_mm_or_si128(
							_mm_or_si128(
									_mm_or_si128(
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 16),
													Match),
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 17),
													Match)),
									_mm_or_si128(
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 18),
													Match),
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 19),
													Match))),
							_mm_or_si128(
									_mm_or_si128(
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 20),
													Match),
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 21),
													Match)),
									_mm_or_si128(
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 22),
													Match),
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 23),
													Match)))),
					_mm_or_si128(
							_mm_or_si128(
									_mm_or_si128(
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 24),
													Match),
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 25),
													Match)),
									_mm_or_si128(
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 26),
													Match),
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 27),
													Match))),
							_mm_or_si128(
									_mm_or_si128(
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 28),
													Match),
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 29),
													Match)),
									_mm_or_si128(
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 30),
													Match),
											_mm_cmpeq_epi32(
													_mm_lddqu_si128(
															(__m128i *) target
																	+ 31),
													Match))))));
	if (_mm_testz_si128(F0, F0))
		*foundp = 0;
	else
		*foundp = 1;
	return (target - init_target);
}
long simdgallop_512_256_rough(int *foundp, UINT4 goal, const UINT4 *target,
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
		}
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
	if (target[255] < goal) {
		target += 256;
	}
	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0 =
			_mm_or_si128(
					_mm_or_si128(
							_mm_or_si128(
									_mm_or_si128(
											_mm_or_si128(
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 0),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 1),
																	Match)),
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 2),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 3),
																	Match))),
											_mm_or_si128(
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 4),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 5),
																	Match)),
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 6),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 7),
																	Match)))),
									_mm_or_si128(
											_mm_or_si128(
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 8),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 9),
																	Match)),
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 10),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 11),
																	Match))),
											_mm_or_si128(
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 12),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 13),
																	Match)),
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 14),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 15),
																	Match))))),
							_mm_or_si128(
									_mm_or_si128(
											_mm_or_si128(
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 16),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 17),
																	Match)),
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 18),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 19),
																	Match))),
											_mm_or_si128(
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 20),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 21),
																	Match)),
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 22),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 23),
																	Match)))),
									_mm_or_si128(
											_mm_or_si128(
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 24),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 25),
																	Match)),
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 26),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 27),
																	Match))),
											_mm_or_si128(
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 28),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 29),
																	Match)),
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 30),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 31),
																	Match)))))),
					_mm_or_si128(
							_mm_or_si128(
									_mm_or_si128(
											_mm_or_si128(
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 32),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 33),
																	Match)),
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 34),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 35),
																	Match))),
											_mm_or_si128(
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 36),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 37),
																	Match)),
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 38),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 39),
																	Match)))),
									_mm_or_si128(
											_mm_or_si128(
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 40),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 41),
																	Match)),
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 42),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 43),
																	Match))),
											_mm_or_si128(
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 44),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 45),
																	Match)),
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 46),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 47),
																	Match))))),
							_mm_or_si128(
									_mm_or_si128(
											_mm_or_si128(
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 48),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 49),
																	Match)),
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 50),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 51),
																	Match))),
											_mm_or_si128(
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 52),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 53),
																	Match)),
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 54),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 55),
																	Match)))),
									_mm_or_si128(
											_mm_or_si128(
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 56),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 57),
																	Match)),
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 58),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 59),
																	Match))),
											_mm_or_si128(
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 60),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 61),
																	Match)),
													_mm_or_si128(
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 62),
																	Match),
															_mm_cmpeq_epi32(
																	_mm_lddqu_si128(
																			(__m128i *) target
																					+ 63),
																	Match)))))));
	if (_mm_testz_si128(F0, F0))
		*foundp = 0;
	else
		*foundp = 1;
	return (target - init_target);
}
long simdgallop_512_512_rough(int *foundp, UINT4 goal, const UINT4 *target,
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
		}
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
	__m128i Match = _mm_set1_epi32(goal);
	__m128i F0 =
			_mm_or_si128(
					_mm_or_si128(
							_mm_or_si128(
									_mm_or_si128(
											_mm_or_si128(
													_mm_or_si128(
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 0),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 1),
																			Match)),
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 2),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 3),
																			Match))),
													_mm_or_si128(
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 4),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 5),
																			Match)),
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 6),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 7),
																			Match)))),
											_mm_or_si128(
													_mm_or_si128(
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 8),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 9),
																			Match)),
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 10),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 11),
																			Match))),
													_mm_or_si128(
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 12),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 13),
																			Match)),
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 14),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 15),
																			Match))))),
									_mm_or_si128(
											_mm_or_si128(
													_mm_or_si128(
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 16),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 17),
																			Match)),
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 18),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 19),
																			Match))),
													_mm_or_si128(
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 20),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 21),
																			Match)),
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 22),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 23),
																			Match)))),
											_mm_or_si128(
													_mm_or_si128(
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 24),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 25),
																			Match)),
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 26),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 27),
																			Match))),
													_mm_or_si128(
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 28),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 29),
																			Match)),
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 30),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 31),
																			Match)))))),
							_mm_or_si128(
									_mm_or_si128(
											_mm_or_si128(
													_mm_or_si128(
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 32),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 33),
																			Match)),
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 34),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 35),
																			Match))),
													_mm_or_si128(
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 36),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 37),
																			Match)),
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 38),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 39),
																			Match)))),
											_mm_or_si128(
													_mm_or_si128(
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 40),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 41),
																			Match)),
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 42),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 43),
																			Match))),
													_mm_or_si128(
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 44),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 45),
																			Match)),
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 46),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 47),
																			Match))))),
									_mm_or_si128(
											_mm_or_si128(
													_mm_or_si128(
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 48),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 49),
																			Match)),
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 50),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 51),
																			Match))),
													_mm_or_si128(
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 52),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 53),
																			Match)),
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 54),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 55),
																			Match)))),
											_mm_or_si128(
													_mm_or_si128(
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 56),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 57),
																			Match)),
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 58),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 59),
																			Match))),
													_mm_or_si128(
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 60),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 61),
																			Match)),
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 62),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 63),
																			Match))))))),
					_mm_or_si128(
							_mm_or_si128(
									_mm_or_si128(
											_mm_or_si128(
													_mm_or_si128(
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 64),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 65),
																			Match)),
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 66),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 67),
																			Match))),
													_mm_or_si128(
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 68),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 69),
																			Match)),
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 70),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 71),
																			Match)))),
											_mm_or_si128(
													_mm_or_si128(
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 72),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 73),
																			Match)),
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 74),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 75),
																			Match))),
													_mm_or_si128(
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 76),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 77),
																			Match)),
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 78),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 79),
																			Match))))),
									_mm_or_si128(
											_mm_or_si128(
													_mm_or_si128(
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 80),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 81),
																			Match)),
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 82),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 83),
																			Match))),
													_mm_or_si128(
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 84),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 85),
																			Match)),
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 86),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 87),
																			Match)))),
											_mm_or_si128(
													_mm_or_si128(
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 88),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 89),
																			Match)),
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 90),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 91),
																			Match))),
													_mm_or_si128(
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 92),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 93),
																			Match)),
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 94),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 95),
																			Match)))))),
							_mm_or_si128(
									_mm_or_si128(
											_mm_or_si128(
													_mm_or_si128(
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 96),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 97),
																			Match)),
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 98),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 99),
																			Match))),
													_mm_or_si128(
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 100),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 101),
																			Match)),
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 102),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 103),
																			Match)))),
											_mm_or_si128(
													_mm_or_si128(
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 104),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 105),
																			Match)),
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 106),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 107),
																			Match))),
													_mm_or_si128(
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 108),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 109),
																			Match)),
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 110),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 111),
																			Match))))),
									_mm_or_si128(
											_mm_or_si128(
													_mm_or_si128(
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 112),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 113),
																			Match)),
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 114),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 115),
																			Match))),
													_mm_or_si128(
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 116),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 117),
																			Match)),
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 118),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 119),
																			Match)))),
											_mm_or_si128(
													_mm_or_si128(
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 120),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 121),
																			Match)),
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 122),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 123),
																			Match))),
													_mm_or_si128(
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 124),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 125),
																			Match)),
															_mm_or_si128(
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 126),
																			Match),
																	_mm_cmpeq_epi32(
																			_mm_lddqu_si128(
																					(__m128i *) target
																							+ 127),
																			Match))))))));
	if (_mm_testz_si128(F0, F0))
		*foundp = 0;
	else
		*foundp = 1;
	return (target - init_target);
}
}
#endif /* INCLUDE_MSIS_GALLOP_ROUGH_HPP_ */
