#ifndef MEDIUMINTERSECTION_H_
#define MEDIUMINTERSECTION_H_

#include "util.h"

#define VEC __m128i
#define VECLEN (sizeof(VEC)/sizeof(uint32_t))
#define VECMAX (VECLEN - 1)

size_t nate_count_scalar(const uint32_t *A, const size_t lenA,
		const uint32_t *B, const size_t lenB) {

	size_t count = 0;
	if (lenA == 0 || lenB == 0)
		return count;

	const uint32_t *endA = A + lenA;
	const uint32_t *endB = B + lenB;

	while (1) {
		while (*A < *B) {
			SKIP_FIRST_COMPARE: if (++A == endA)
				return count;
		}
		while (*A > *B) {
			if (++B == endB)
				return count;
		}
		if (*A == *B) {
			count++;
			if (++A == endA || ++B == endB)
				return count;
		} else {
			goto SKIP_FIRST_COMPARE;
		}
	}

	return count; // NOTREACHED
}

size_t nate_scalar(const uint32_t *A, const size_t lenA, const uint32_t *B,
		const size_t lenB, uint32_t * out) {
	const uint32_t * const initout(out);
	if (lenA == 0 || lenB == 0)
		return 0;

	const uint32_t *endA = A + lenA;
	const uint32_t *endB = B + lenB;

	while (1) {
		while (*A < *B) {
			SKIP_FIRST_COMPARE: if (++A == endA)
				return (out - initout);
		}
		while (*A > *B) {
			if (++B == endB)
				return (out - initout);
		}
		if (*A == *B) {
			*out++ = *A;
			if (++A == endA || ++B == endB)
				return (out - initout);
		} else {
			goto SKIP_FIRST_COMPARE;
		}
	}

	return (out - initout); // NOTREACHED
}

/**
 * Version of nate_scalar modified by D. Lemire
 * to have no goto.
 */
size_t nate_scalarwithoutgoto(const uint32_t *A, const size_t lenA,
		const uint32_t *B, const size_t lenB, uint32_t * out) {
	const uint32_t * const initout(out);
	if (lenA == 0 || lenB == 0)
		return 0;

	const uint32_t *endA = A + lenA;
	const uint32_t *endB = B + lenB;

	while (1) {
		while (*A < *B) {
			if (++A == endA)
				return (out - initout);
		}
		while (*A > *B) {
			if (++B == endB)
				return (out - initout);
		}
		if (*A == *B) {
			*out++ = *A;
			if (++A == endA || ++B == endB)
				return (out - initout);
		} else {
			if (++A == endA)
				return (out - initout);
		}
	}

	return (out - initout); // NOTREACHED
}

size_t nate_count_medium(const uint32_t *rare, const size_t lenRare,
		const uint32_t *freq, const size_t lenFreq) {

	// FUTURE: could swap freq and rare if inverted

	size_t count = 0;
	if (lenFreq == 0 || lenRare == 0)
		return count;

#define FREQSPACE (8 * VECLEN)
#define RARESPACE 1
	const uint32_t *stopFreq = freq + lenFreq - FREQSPACE;
	const uint32_t *stopRare = rare + lenRare - RARESPACE;

	if (freq > stopFreq) {
		return nate_count_scalar(freq, lenFreq, rare, lenRare);
	}

	uint32_t nextRare = *rare;
	uint32_t maxFreq = freq[7 * VECLEN + VECMAX];
	VEC M0, M1, M2, M3, M4, M5, M6, M7;
	VEC A0, A1, A2, A3;
	A0 = _mm_load_si128((VEC *) freq + 0);
	A1 = _mm_load_si128((VEC *) freq + 1);
	A2 = _mm_load_si128((VEC *) freq + 2);
	A3 = _mm_load_si128((VEC *) freq + 3);

	while (maxFreq < nextRare) { // advance freq to a possible match
		freq += VECLEN * 8; // NOTE: loop below requires this
		if (freq > stopFreq)
			goto FINISH_SCALAR;
		maxFreq = freq[VECLEN * 7 + VECMAX];
		A0 = _mm_load_si128((VEC *) freq + 0);
		A1 = _mm_load_si128((VEC *) freq + 1);
		A2 = _mm_load_si128((VEC *) freq + 2);
		A3 = _mm_load_si128((VEC *) freq + 3);
	}

	while (rare < stopRare) {
		uint32_t matchRare = nextRare;
		VEC Match = _mm_set1_epi32(matchRare);
		if (maxFreq < matchRare) { // if no match possible
			freq += VECLEN * 8; // advance 8 vectors
			if (freq > stopFreq)
				goto FINISH_SCALAR;
			maxFreq = freq[VECLEN * 7 + VECMAX];
			while (maxFreq < matchRare) { // if still no match possible
				freq += VECLEN * 8; // advance another 8 vectors
				if (freq > stopFreq)
					goto FINISH_SCALAR;
				maxFreq = freq[VECLEN * 7 + VECMAX];
				// printf("FREQ_ADVANCE2: %d %d\n", matchRare, maxFreq);
				M0 = _mm_load_si128((VEC *) freq + 0);
				M1 = _mm_load_si128((VEC *) freq + 1);
				M2 = _mm_load_si128((VEC *) freq + 2);
				M3 = _mm_load_si128((VEC *) freq + 3);
			}
			A0 = M0; // ratchet original vectors forward
			A1 = M1;
			A2 = M2;
			A3 = M3;
		} else { // if a match is still possible
			M0 = A0; // revert to original vectors
			M1 = A1;
			M2 = A2;
			M3 = A3;
		}

		rare += 1;
		nextRare = *rare;
		// NOTE: At this point M0, M1, M2, M3 must match freq + 0,1,2,3
		VEC Q0, Q1, Q2, Q3; // quarters
		VEC S0, S1; // semis
		VEC F0; // final

		M4 = _mm_load_si128((VEC *) freq + 4);
		M5 = _mm_load_si128((VEC *) freq + 5);
		M6 = _mm_load_si128((VEC *) freq + 6);
		M7 = _mm_load_si128((VEC *) freq + 7);

		M0 = _mm_cmpeq_epi32(M0, Match);
		M1 = _mm_cmpeq_epi32(M1, Match);
		Q0 = _mm_or_si128(M0, M1);
		M0 = _mm_load_si128((VEC *) freq + 8);
		M1 = _mm_load_si128((VEC *) freq + 9);

		M2 = _mm_cmpeq_epi32(M2, Match);
		M3 = _mm_cmpeq_epi32(M3, Match);
		Q1 = _mm_or_si128(M2, M3);
		M2 = _mm_load_si128((VEC *) freq + 10);
		M3 = _mm_load_si128((VEC *) freq + 11);

		M4 = _mm_cmpeq_epi32(M4, Match);
		M5 = _mm_cmpeq_epi32(M5, Match);
		Q2 = _mm_or_si128(M4, M5);

		M6 = _mm_cmpeq_epi32(M6, Match);
		M7 = _mm_cmpeq_epi32(M7, Match);
		Q3 = _mm_or_si128(M6, M7);

		S0 = _mm_or_si128(Q0, Q1);
		S1 = _mm_or_si128(Q2, Q3);
		F0 = _mm_or_si128(S0, S1);

		// PROFILE: check that this compiles to a CMOV
		if (_mm_testz_si128(F0, F0)) {
			// printf("NOT_FOUND: %d %d\n", matchRare, maxFreq);
		} else {
			// printf("FOUND: %d %d\n", matchRare, maxFreq);
			count += 1;
		}
	}

	FINISH_SCALAR: return count
			+ nate_count_scalar(freq, stopFreq + FREQSPACE - freq, rare,
					stopRare + RARESPACE - rare);
}

/**
 * Version hacked by D. Lemire, original by Nathan Kurz
 */
size_t natedan_count_medium(const uint32_t *rare, const size_t lenRare,
		const uint32_t *freq, const size_t lenFreq) {
	// FUTURE: could swap freq and rare if inverted
	size_t count = 0;
	if (lenFreq == 0 || lenRare == 0)
		return count;

	const size_t freqspace = 8 * VECLEN;
	const size_t rarespace = 1;

	const uint32_t *stopFreq = freq + lenFreq - freqspace;
	const uint32_t *stopRare = rare + lenRare - rarespace;
	if (freq > stopFreq) {
		return nate_count_scalar(freq, lenFreq, rare, lenRare);
//        goto FINISH_SCALAR;
	}
	uint32_t maxFreq = freq[7 * VECLEN + VECMAX];
	VEC M0, M1, M2, M3, M4, M5, M6, M7;

	while (maxFreq < *rare) { // advance freq to a possible match
		freq += VECLEN * 8; // NOTE: loop below requires this
		if (freq > stopFreq)
			goto FINISH_SCALAR;
		maxFreq = freq[VECLEN * 7 + VECMAX];
	}
	M0 = _mm_load_si128((VEC *) freq + 0);
	M1 = _mm_load_si128((VEC *) freq + 1);
	M2 = _mm_load_si128((VEC *) freq + 2);
	M3 = _mm_load_si128((VEC *) freq + 3);
	M4 = _mm_load_si128((VEC *) freq + 4);
	M5 = _mm_load_si128((VEC *) freq + 5);
	M6 = _mm_load_si128((VEC *) freq + 6);
	M7 = _mm_load_si128((VEC *) freq + 7);

	for (; rare < stopRare; ++rare) {
		const uint32_t matchRare = *rare; //nextRare;
		const VEC Match = _mm_set1_epi32(matchRare);
		// sxs: the following seems like can be further folded
		// but Lemire did this to avoid too much copy operation
		// of M0~M7, since if-branch can be skipped, but while
		// cannot.
		if (maxFreq < matchRare) { // if no match possible
			freq += VECLEN * 8; // advance 8 vectors
			if (freq > stopFreq)
				goto FINISH_SCALAR;
			maxFreq = freq[VECLEN * 7 + VECMAX];
			while (maxFreq < matchRare) { // if still no match possible
				freq += VECLEN * 8; // advance another 8 vectors
				if (freq > stopFreq)
					goto FINISH_SCALAR;
				maxFreq = freq[VECLEN * 7 + VECMAX];
			}
			M0 = _mm_load_si128((VEC *) freq + 0);
			M1 = _mm_load_si128((VEC *) freq + 1);
			M2 = _mm_load_si128((VEC *) freq + 2);
			M3 = _mm_load_si128((VEC *) freq + 3);
			M4 = _mm_load_si128((VEC *) freq + 4);
			M5 = _mm_load_si128((VEC *) freq + 5);
			M6 = _mm_load_si128((VEC *) freq + 6);
			M7 = _mm_load_si128((VEC *) freq + 7);

		}
		// sxs: found possible matchs in current range
		const VEC Q0 = _mm_or_si128(_mm_cmpeq_epi32(M0, Match),
				_mm_cmpeq_epi32(M1, Match));
		const VEC Q1 = _mm_or_si128(_mm_cmpeq_epi32(M2, Match),
				_mm_cmpeq_epi32(M3, Match));
		const VEC Q2 = _mm_or_si128(_mm_cmpeq_epi32(M4, Match),
				_mm_cmpeq_epi32(M5, Match));
		const VEC Q3 = _mm_or_si128(_mm_cmpeq_epi32(M6, Match),
				_mm_cmpeq_epi32(M7, Match));
		const VEC S0 = _mm_or_si128(Q0, Q1);
		const VEC S1 = _mm_or_si128(Q2, Q3);
		const VEC F0 = _mm_or_si128(S0, S1);
		if (_mm_testz_si128(F0, F0)) {
		} else {
			count += 1;
		}
	}

	FINISH_SCALAR: return count
			+ nate_count_scalar(freq, stopFreq + freqspace - freq, rare,
					stopRare + rarespace - rare);
}

/**
 * Version hacked by D. Lemire, original by Nathan Kurz
 * sxs: omit temporary variables in the former function `natedan_count_medium`
 */
size_t natedanalt_count_medium(const uint32_t *rare, const size_t lenRare,
		const uint32_t *freq, const size_t lenFreq) {
	// FUTURE: could swap freq and rare if inverted
	size_t count = 0;
	if (lenFreq == 0 || lenRare == 0)
		return count;

	const size_t freqspace = 8 * VECLEN;
	const size_t rarespace = 1;

	const uint32_t *stopFreq = freq + lenFreq - freqspace;
	const uint32_t *stopRare = rare + lenRare - rarespace;
	if (freq > stopFreq) {
		return nate_count_scalar(freq, lenFreq, rare, lenRare);
	}
	uint32_t maxFreq = freq[7 * VECLEN + VECMAX];

	while (maxFreq < *rare) { // advance freq to a possible match
		freq += VECLEN * 8; // NOTE: loop below requires this
		if (freq > stopFreq)
			goto FINISH_SCALAR;
		maxFreq = freq[VECLEN * 7 + VECMAX];
	}
	for (; rare < stopRare; ++rare) {
		const uint32_t matchRare = *rare; //nextRare;
		const VEC Match = _mm_set1_epi32(matchRare);
		while (maxFreq < matchRare) { // if no match possible
			freq += VECLEN * 8; // advance 8 vectors
			if (freq > stopFreq)
				goto FINISH_SCALAR;
			maxFreq = freq[VECLEN * 7 + VECMAX];
		}
		VEC F0;
		// sxs: that's how medium get its name
		if (freq[VECLEN * 3 + VECMAX] < matchRare) {
			const VEC Q2 = _mm_or_si128(
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 4), Match),
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 5), Match));
			const VEC Q3 = _mm_or_si128(
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 6), Match),
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 7), Match));
			F0 = _mm_or_si128(Q2, Q3);
		} else {
			const VEC Q0 = _mm_or_si128(
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 0), Match),
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 1), Match));
			const VEC Q1 = _mm_or_si128(
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 2), Match),
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 3), Match));
			F0 = _mm_or_si128(Q0, Q1);
		}
		if (_mm_testz_si128(F0, F0)) {
		} else {
			count += 1;
		}
	}

	FINISH_SCALAR: return count
			+ nate_count_scalar(freq, stopFreq + freqspace - freq, rare,
					stopRare + rarespace - rare);
}

/**
 * Version hacked by D. Lemire, original by Nathan Kurz
 */
size_t natedanalt_medium(const uint32_t *rare, const size_t lenRare,
		const uint32_t *freq, const size_t lenFreq, uint32_t * out) {
	// FUTURE: could swap freq and rare if inverted
	const uint32_t * const initout(out);
	size_t count = 0;
	if (lenFreq == 0 || lenRare == 0)
		return count;
	assert(lenRare <= lenFreq);

	const size_t freqspace = 8 * VECLEN;
	const size_t rarespace = 1;

	const uint32_t *stopFreq = freq + lenFreq - freqspace;
	const uint32_t *stopRare = rare + lenRare - rarespace;
	if (freq > stopFreq) {
		return nate_scalar(freq, lenFreq, rare, lenRare, out);
	}
	uint32_t maxFreq = freq[7 * VECLEN + VECMAX];

	while (maxFreq < *rare) { // advance freq to a possible match
		freq += VECLEN * 8; // NOTE: loop below requires this
		if (freq > stopFreq)
			goto FINISH_SCALAR;
		maxFreq = freq[VECLEN * 7 + VECMAX];
	}
	for (; rare < stopRare; ++rare) {
		const uint32_t matchRare = *rare; //nextRare;
		const VEC Match = _mm_set1_epi32(matchRare);
		while (maxFreq < matchRare) { // if no match possible
			freq += VECLEN * 8; // advance 8 vectors
			if (freq > stopFreq)
				goto FINISH_SCALAR;
			maxFreq = freq[VECLEN * 7 + VECMAX];
		}
		VEC F0;
		if (freq[VECLEN * 3 + VECMAX] < matchRare) {
			const VEC Q2 = _mm_or_si128(
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 4), Match),
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 5), Match));
			const VEC Q3 = _mm_or_si128(
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 6), Match),
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 7), Match));
			F0 = _mm_or_si128(Q2, Q3);
		} else {
			const VEC Q0 = _mm_or_si128(
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 0), Match),
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 1), Match));
			const VEC Q1 = _mm_or_si128(
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 2), Match),
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 3), Match));
			F0 = _mm_or_si128(Q0, Q1);
		}
		if (_mm_testz_si128(F0, F0)) {
		} else {
			*out++ = matchRare;
			// count += 1;
		}
	}

	FINISH_SCALAR: return (out - initout)
			+ nate_scalar(freq, stopFreq + freqspace - freq, rare,
					stopRare + rarespace - rare, out);
}

/**
 * Version hacked by D. Lemire, original by Nathan Kurz
 */
size_t danfar_count_medium(const uint32_t *rare, const size_t lenRare,
		const uint32_t *freq, const size_t lenFreq) {
	size_t count = 0;
	if (lenFreq == 0 || lenRare == 0)
		return count;

	const size_t freqspace = 16 * VECLEN;
	const size_t rarespace = 1;

	const uint32_t *stopFreq = freq + lenFreq - freqspace;
	const uint32_t *stopRare = rare + lenRare - rarespace;
	if (freq > stopFreq) {
		return nate_count_scalar(freq, lenFreq, rare, lenRare);
	}
	while (freq[VECLEN * 15 + VECMAX] < *rare) {
		freq += VECLEN * 16;
		if (freq > stopFreq)
			goto FINISH_SCALAR;
	}
	for (; rare < stopRare; ++rare) {
		const uint32_t matchRare = *rare; //nextRare;
		const VEC Match = _mm_set1_epi32(matchRare);
		while (freq[VECLEN * 15 + VECMAX] < matchRare) { // if no match possible
			freq += VECLEN * 16; // advance 8 vectors
			if (freq > stopFreq)
				goto FINISH_SCALAR;
		}
		VEC Q0, Q1, Q2, Q3;
		if (freq[VECLEN * 7 + VECMAX] < matchRare) {
			Q0 = _mm_or_si128(
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 12), Match),
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 13), Match));
			Q1 = _mm_or_si128(
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 14), Match),
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 15), Match));

			Q2 = _mm_or_si128(
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 8), Match),
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 9), Match));
			Q3 = _mm_or_si128(
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 10), Match),
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 11), Match));
		} else {
			Q0 = _mm_or_si128(
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 4), Match),
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 5), Match));
			Q1 = _mm_or_si128(
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 6), Match),
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 7), Match));
			Q2 = _mm_or_si128(
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 0), Match),
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 1), Match));
			Q3 = _mm_or_si128(
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 2), Match),
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 3), Match));
		}
		const VEC F0 = _mm_or_si128(_mm_or_si128(Q0, Q1), _mm_or_si128(Q2, Q3));
		if (_mm_testz_si128(F0, F0)) {
		} else {
			count += 1;
		}
	}

	FINISH_SCALAR: return count
			+ nate_count_scalar(freq, stopFreq + freqspace - freq, rare,
					stopRare + rarespace - rare);
}

/**
 * Version hacked by D. Lemire, original by Nathan Kurz
 */
size_t danfar_medium(const uint32_t *rare, const size_t lenRare,
		const uint32_t *freq, const size_t lenFreq, uint32_t * out) {
	if (lenFreq == 0 || lenRare == 0)
		return 0;
	assert(lenRare <= lenFreq);
	const uint32_t * const initout(out);

	const size_t freqspace = 16 * VECLEN;
	const size_t rarespace = 1;

	const uint32_t *stopFreq = freq + lenFreq - freqspace;
	const uint32_t *stopRare = rare + lenRare - rarespace;
	if (freq > stopFreq) {
		return nate_scalar(freq, lenFreq, rare, lenRare, out);
	}
	while (freq[VECLEN * 15 + VECMAX] < *rare) {
		freq += VECLEN * 16;
		if (freq > stopFreq)
			goto FINISH_SCALAR;
	}
	for (; rare < stopRare; ++rare) {
		const uint32_t matchRare = *rare; //nextRare;
		const VEC Match = _mm_set1_epi32(matchRare);
		while (freq[VECLEN * 15 + VECMAX] < matchRare) { // if no match possible
			freq += VECLEN * 16; // advance 8 vectors
			if (freq > stopFreq)
				goto FINISH_SCALAR;
		}
		VEC Q0, Q1, Q2, Q3;
		if (freq[VECLEN * 7 + VECMAX] < matchRare) {
			Q0 = _mm_or_si128(
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 8), Match),
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 9), Match));
			Q1 = _mm_or_si128(
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 10), Match),
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 11), Match));

			Q2 = _mm_or_si128(
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 12), Match),
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 13), Match));
			Q3 = _mm_or_si128(
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 14), Match),
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 15), Match));
		} else {
			Q0 = _mm_or_si128(
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 4), Match),
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 5), Match));
			Q1 = _mm_or_si128(
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 6), Match),
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 7), Match));
			Q2 = _mm_or_si128(
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 0), Match),
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 1), Match));
			Q3 = _mm_or_si128(
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 2), Match),
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 3), Match));
		}
		const VEC F0 = _mm_or_si128(_mm_or_si128(Q0, Q1), _mm_or_si128(Q2, Q3));
		if (_mm_testz_si128(F0, F0)) {
		} else {
			*out++ = matchRare;
		}
	}
	FINISH_SCALAR: return (out - initout)
			+ nate_scalar(freq, stopFreq + freqspace - freq, rare,
					stopRare + rarespace - rare, out);
}

/**
 * Version hacked by D. Lemire, original by Nathan Kurz
 * sxs: only difference is the last step which zero-tests the mask
 * _mm_testz_si128 or _mm_movemask_epi8
 */
size_t danfar_medium_mov(const uint32_t *rare, const size_t lenRare,
		const uint32_t *freq, const size_t lenFreq, uint32_t * out) {
	if (lenFreq == 0 || lenRare == 0)
		return 0;
	assert(lenRare <= lenFreq);
	const uint32_t * const initout(out);

	const size_t freqspace = 16 * VECLEN;
	const size_t rarespace = 1;

	const uint32_t *stopFreq = freq + lenFreq - freqspace;
	const uint32_t *stopRare = rare + lenRare - rarespace;
	if (freq > stopFreq) {
		return nate_scalar(freq, lenFreq, rare, lenRare, out);
	}
	while (freq[VECLEN * 15 + VECMAX] < *rare) {
		freq += VECLEN * 16;
		if (freq > stopFreq)
			goto FINISH_SCALAR;
	}
	for (; rare < stopRare; ++rare) {
		const uint32_t matchRare = *rare; //nextRare;
		const VEC Match = _mm_set1_epi32(matchRare);
		while (freq[VECLEN * 15 + VECMAX] < matchRare) { // if no match possible
			freq += VECLEN * 16; // advance 8 vectors
			if (freq > stopFreq)
				goto FINISH_SCALAR;
		}
		VEC Q0, Q1, Q2, Q3;
		if (freq[VECLEN * 7 + VECMAX] < matchRare) {
			Q0 = _mm_or_si128(
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 8), Match),
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 9), Match));
			Q1 = _mm_or_si128(
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 10), Match),
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 11), Match));

			Q2 = _mm_or_si128(
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 12), Match),
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 13), Match));
			Q3 = _mm_or_si128(
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 14), Match),
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 15), Match));
		} else {
			Q0 = _mm_or_si128(
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 4), Match),
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 5), Match));
			Q1 = _mm_or_si128(
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 6), Match),
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 7), Match));
			Q2 = _mm_or_si128(
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 0), Match),
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 1), Match));
			Q3 = _mm_or_si128(
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 2), Match),
					_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 3), Match));
		}
		const VEC F0 = _mm_or_si128(_mm_or_si128(Q0, Q1), _mm_or_si128(Q2, Q3));
		if (_mm_movemask_epi8(F0)) {
			*out++ = matchRare;
		}
	}

	FINISH_SCALAR: return (out - initout)
			+ nate_scalar(freq, stopFreq + freqspace - freq, rare,
					stopRare + rarespace - rare, out);
}

size_t danfarfar_medium(const uint32_t *rare, const size_t lenRare,
		const uint32_t *freq, const size_t lenFreq, uint32_t * out) {
	if (lenFreq == 0 || lenRare == 0)
		return 0;
	assert(lenRare <= lenFreq);
	const uint32_t * const initout(out);

	const size_t freqspace = 32 * VECLEN;
	const size_t rarespace = 1;

	const uint32_t *stopFreq = freq + lenFreq - freqspace;
	const uint32_t *stopRare = rare + lenRare - rarespace;
	if (freq > stopFreq) {
		return nate_scalar(freq, lenFreq, rare, lenRare, out);
	}
	while (freq[VECLEN * 31 + VECMAX] < *rare) {
		freq += VECLEN * 32;
		if (freq > stopFreq)
			goto FINISH_SCALAR;
	}
	for (; rare < stopRare; ++rare) {
		const uint32_t matchRare = *rare; //nextRare;
		const VEC Match = _mm_set1_epi32(matchRare);
		while (freq[VECLEN * 31 + VECMAX] < matchRare) { // if no match possible
			freq += VECLEN * 32; // advance 32 vectors
			if (freq > stopFreq)
				goto FINISH_SCALAR;
		}
		VEC Q0, Q1, Q2, Q3;
		if (freq[VECLEN * 15 + VECMAX] >= matchRare) {
			if (freq[VECLEN * 7 + VECMAX] < matchRare) {
				Q0 = _mm_or_si128(
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 8),
								Match),
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 9),
								Match));
				Q1 = _mm_or_si128(
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 10),
								Match),
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 11),
								Match));

				Q2 = _mm_or_si128(
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 12),
								Match),
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 13),
								Match));
				Q3 = _mm_or_si128(
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 14),
								Match),
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 15),
								Match));
			} else {
				Q0 = _mm_or_si128(
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 4),
								Match),
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 5),
								Match));
				Q1 = _mm_or_si128(
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 6),
								Match),
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 7),
								Match));
				Q2 = _mm_or_si128(
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 0),
								Match),
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 1),
								Match));
				Q3 = _mm_or_si128(
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 2),
								Match),
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 3),
								Match));
			}
		} else {
			if (freq[VECLEN * 23 + VECMAX] < matchRare) {
				Q0 = _mm_or_si128(
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 8 + 16),
								Match),
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 9 + 16),
								Match));
				Q1 = _mm_or_si128(
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 10 + 16),
								Match),
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 11 + 16),
								Match));

				Q2 = _mm_or_si128(
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 12 + 16),
								Match),
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 13 + 16),
								Match));
				Q3 = _mm_or_si128(
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 14 + 16),
								Match),
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 15 + 16),
								Match));
			} else {
				Q0 = _mm_or_si128(
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 4 + 16),
								Match),
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 5 + 16),
								Match));
				Q1 = _mm_or_si128(
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 6 + 16),
								Match),
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 7 + 16),
								Match));
				Q2 = _mm_or_si128(
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 0 + 16),
								Match),
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 1 + 16),
								Match));
				Q3 = _mm_or_si128(
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 2 + 16),
								Match),
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 3 + 16),
								Match));
			}

		}
		const VEC F0 = _mm_or_si128(_mm_or_si128(Q0, Q1), _mm_or_si128(Q2, Q3));
		if (_mm_testz_si128(F0, F0)) {
		} else {
			*out++ = matchRare;
		}
	}

	FINISH_SCALAR: return (out - initout)
			+ nate_scalar(freq, stopFreq + freqspace - freq, rare,
					stopRare + rarespace - rare, out);
}

//proxy for danfarfar_medium
size_t v3(const uint32_t *rare, const size_t lenRare, const uint32_t *freq,
		const size_t lenFreq, uint32_t * out) {
	return danfarfar_medium(rare, lenRare, freq, lenFreq, out);
}

#ifdef __AVX2__
#include <immintrin.h>

size_t danfarfar_medium_avx2(const uint32_t *rare, const size_t lenRare,
		const uint32_t *freq, const size_t lenFreq, uint32_t * out) {
	if (lenFreq == 0 || lenRare == 0)
	return 0;
	assert(lenRare <= lenFreq);
	const uint32_t * const initout (out);
	typedef __m256i vec;
	const uint32_t veclen = sizeof(vec) / sizeof(uint32_t);
	const size_t vecmax = veclen - 1;
	const size_t freqspace = 32 * veclen;
	const size_t rarespace = 1;

	const uint32_t *stopFreq = freq + lenFreq - freqspace;
	const uint32_t *stopRare = rare + lenRare - rarespace;
	if (freq > stopFreq) {
		return nate_scalar(freq, lenFreq, rare, lenRare, out);
	}
	while (freq[veclen * 31 + vecmax] < *rare) {
		freq += veclen * 32;
		if (freq > stopFreq)
		goto FINISH_SCALAR;
	}
	for (; rare < stopRare; ++rare) {
		const uint32_t matchRare = *rare; //nextRare;
		const vec Match = _mm256_set1_epi32(matchRare);
		while (freq[veclen * 31 + vecmax] < matchRare) { // if no match possible
			freq += veclen * 32;// advance 32 vectors
			if (freq > stopFreq)
			goto FINISH_SCALAR;
		}
		vec Q0,Q1,Q2,Q3;
		if(freq[veclen * 15 + vecmax] >= matchRare ) {
			if(freq[veclen * 7 + vecmax] < matchRare ) {
				Q0 = _mm256_or_si256(
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 8), Match),
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 9), Match));
				Q1 = _mm256_or_si256(
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 10), Match),
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 11), Match));

				Q2 = _mm256_or_si256(
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 12), Match),
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 13), Match));
				Q3 = _mm256_or_si256(
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 14), Match),
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 15), Match));
			} else {
				Q0 = _mm256_or_si256(
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 4), Match),
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 5), Match));
				Q1 = _mm256_or_si256(
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 6), Match),
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 7), Match));
				Q2 = _mm256_or_si256(
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 0), Match),
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 1), Match));
				Q3 = _mm256_or_si256(
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 2), Match),
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 3), Match));
			}
		}
		else
		{
			if(freq[veclen * 23 + vecmax] < matchRare ) {
				Q0 = _mm256_or_si256(
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 8 + 16), Match),
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 9 + 16), Match));
				Q1 = _mm256_or_si256(
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 10+ 16), Match),
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 11+ 16), Match));

				Q2 = _mm256_or_si256(
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 12+ 16), Match),
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 13+ 16), Match));
				Q3 = _mm256_or_si256(
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 14+ 16), Match),
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 15+ 16), Match));
			} else {
				Q0 = _mm256_or_si256(
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 4+ 16), Match),
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 5+ 16), Match));
				Q1 = _mm256_or_si256(
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 6+ 16), Match),
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 7+ 16), Match));
				Q2 = _mm256_or_si256(
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 0+ 16), Match),
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 1+ 16), Match));
				Q3 = _mm256_or_si256(
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 2+ 16), Match),
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 3+ 16), Match));
			}

		}
		const vec F0 = _mm256_or_si256(_mm256_or_si256(Q0, Q1),_mm256_or_si256(Q2, Q3));
		if (_mm256_testz_si256(F0, F0)) {
		} else {
			*out++ = matchRare;
		}
	}

	FINISH_SCALAR: return (out - initout) + nate_scalar(freq,
			stopFreq + freqspace - freq, rare, stopRare + rarespace - rare, out);
}

//proxy for danfarfar_medium
size_t v3avx2(const uint32_t *rare, const size_t lenRare,
		const uint32_t *freq, const size_t lenFreq, uint32_t * out) {
	return danfarfar_medium_avx2(rare,lenRare,freq,lenFreq,out);
}
#endif

size_t danfarfine_count_medium(const uint32_t *rare, const size_t lenRare,
		const uint32_t *freq, const size_t lenFreq) {
	size_t count = 0;
	if (lenFreq == 0 || lenRare == 0)
		return count;

	const size_t freqspace = 16 * VECLEN;
	const size_t rarespace = 1;

	const uint32_t *stopFreq = freq + lenFreq - freqspace;
	const uint32_t *stopRare = rare + lenRare - rarespace;
	if (freq > stopFreq) {
		return nate_count_scalar(freq, lenFreq, rare, lenRare);
	}
	while (freq[VECLEN * 15 + VECMAX] < *rare) {
		freq += VECLEN * 16;
		if (freq > stopFreq)
			goto FINISH_SCALAR;
	}
	for (; rare < stopRare; ++rare) {
		const uint32_t matchRare = *rare; //nextRare;
		const VEC Match = _mm_set1_epi32(matchRare);
		while (freq[VECLEN * 15 + VECMAX] < matchRare) { // if no match possible
			freq += VECLEN * 16; // advance 8 vectors
			if (freq > stopFreq)
				goto FINISH_SCALAR;
		}
		VEC Q0, Q1;
		if (freq[VECLEN * 7 + VECMAX] < matchRare) {
			if (freq[VECLEN * 11 + VECMAX] < matchRare) {
				Q0 = _mm_or_si128(
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 12),
								Match),
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 13),
								Match));
				Q1 = _mm_or_si128(
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 14),
								Match),
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 15),
								Match));

			} else {
				Q0 = _mm_or_si128(
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 8),
								Match),
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 9),
								Match));
				Q1 = _mm_or_si128(
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 10),
								Match),
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 11),
								Match));
			}
		} else {
			if (freq[VECLEN * 3 + VECMAX] < matchRare) {
				Q0 = _mm_or_si128(
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 4),
								Match),
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 5),
								Match));
				Q1 = _mm_or_si128(
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 6),
								Match),
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 7),
								Match));
			} else {
				Q0 = _mm_or_si128(
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 0),
								Match),
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 1),
								Match));
				Q1 = _mm_or_si128(
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 2),
								Match),
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 3),
								Match));
			}
		}
		const VEC F0 = _mm_or_si128(Q0, Q1);
		if (_mm_testz_si128(F0, F0)) {
		} else {
			count += 1;
		}
	}

	FINISH_SCALAR: return count
			+ nate_count_scalar(freq, stopFreq + freqspace - freq, rare,
					stopRare + rarespace - rare);
}

size_t SIMDgalloping(const uint32_t *rare, const size_t lenRare,
		const uint32_t *freq, const size_t lenFreq, uint32_t * out) {
	if (lenFreq == 0 || lenRare == 0)
		return 0;
	assert(lenRare <= lenFreq);
	const uint32_t * const initout(out);

	const size_t freqspace = 32 * VECLEN;
	const size_t rarespace = 1;

	const uint32_t *stopFreq = freq + lenFreq - freqspace;
	const uint32_t *stopRare = rare + lenRare - rarespace;
	if (freq > stopFreq) {
		return nate_scalar(freq, lenFreq, rare, lenRare, out);
	}
	for (; rare < stopRare; ++rare) {
		const uint32_t matchRare = *rare; //nextRare;
		const VEC Match = _mm_set1_epi32(matchRare);

		if (freq[VECLEN * 31 + VECMAX] < matchRare) { // if no match possible
			uint32_t offset = 1;
			if (freq + freqspace > stopFreq) {
				freq += freqspace;
				goto FINISH_SCALAR;
			}
			while (freq[freqspace * offset + VECLEN * 31 + VECMAX] < matchRare) {
				// if no match possible
				// sxs: keep on doubling the offset until we find the
				// range [offset/2,offset] where the match may fall into.
				// in case offset*2 exceeds the boundary, it will be fixed
				// to the last stopFreq
				if (freq + freqspace * (2 * offset) <= stopFreq) {
					offset *= 2;
				} else if (freq + freqspace * (offset + 1) <= stopFreq) {
					offset =
							static_cast<uint32_t>((stopFreq - freq) / freqspace);
					//offset += 1;
					if (freq[freqspace * offset + VECLEN * 31 + VECMAX]
							< matchRare) {
						freq += freqspace * offset;
						goto FINISH_SCALAR;
					} else {
						break;
					}
				} else {
					freq += freqspace * offset;
					goto FINISH_SCALAR;
				}
			}
			uint32_t lower = offset / 2;
			while (lower + 1 != offset) {
				const uint32_t mid = (lower + offset) / 2;
				if (freq[freqspace * mid + VECLEN * 31 + VECMAX] < matchRare)
					lower = mid;
				else
					offset = mid;
			}
			freq += freqspace * offset;
		}
		VEC Q0, Q1, Q2, Q3;
		if (freq[VECLEN * 15 + VECMAX] >= matchRare) {
			if (freq[VECLEN * 7 + VECMAX] < matchRare) {
				Q0 = _mm_or_si128(
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 8),
								Match),
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 9),
								Match));
				Q1 = _mm_or_si128(
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 10),
								Match),
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 11),
								Match));

				Q2 = _mm_or_si128(
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 12),
								Match),
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 13),
								Match));
				Q3 = _mm_or_si128(
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 14),
								Match),
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 15),
								Match));
			} else {
				Q0 = _mm_or_si128(
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 4),
								Match),
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 5),
								Match));
				Q1 = _mm_or_si128(
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 6),
								Match),
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 7),
								Match));
				Q2 = _mm_or_si128(
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 0),
								Match),
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 1),
								Match));
				Q3 = _mm_or_si128(
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 2),
								Match),
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 3),
								Match));
			}
		} else {
			if (freq[VECLEN * 23 + VECMAX] < matchRare) {
				Q0 = _mm_or_si128(
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 8 + 16),
								Match),
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 9 + 16),
								Match));
				Q1 = _mm_or_si128(
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 10 + 16),
								Match),
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 11 + 16),
								Match));

				Q2 = _mm_or_si128(
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 12 + 16),
								Match),
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 13 + 16),
								Match));
				Q3 = _mm_or_si128(
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 14 + 16),
								Match),
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 15 + 16),
								Match));
			} else {
				Q0 = _mm_or_si128(
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 4 + 16),
								Match),
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 5 + 16),
								Match));
				Q1 = _mm_or_si128(
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 6 + 16),
								Match),
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 7 + 16),
								Match));
				Q2 = _mm_or_si128(
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 0 + 16),
								Match),
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 1 + 16),
								Match));
				Q3 = _mm_or_si128(
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 2 + 16),
								Match),
						_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 3 + 16),
								Match));
			}

		}
		const VEC F0 = _mm_or_si128(_mm_or_si128(Q0, Q1), _mm_or_si128(Q2, Q3));
		if (_mm_testz_si128(F0, F0)) {
		} else {
			*out++ = matchRare;
		}
	}

	FINISH_SCALAR: return (out - initout)
			+ nate_scalar(freq, stopFreq + freqspace - freq, rare,
					stopRare + rarespace - rare, out);
}

#ifdef __AVX2__
#include <immintrin.h>
size_t SIMDgalloping_avx2(const uint32_t *rare, const size_t lenRare,
		const uint32_t *freq, const size_t lenFreq, uint32_t * out) {
	if (lenFreq == 0 || lenRare == 0)
	return 0;
	assert(lenRare <= lenFreq);
	const uint32_t * const initout (out);
	typedef __m256i vec;
	const uint32_t veclen = sizeof(vec) / sizeof(uint32_t);
	const size_t vecmax = veclen - 1;
	const size_t freqspace = 32 * veclen;
	const size_t rarespace = 1;

	const uint32_t *stopFreq = freq + lenFreq - freqspace;
	const uint32_t *stopRare = rare + lenRare - rarespace;
	if (freq > stopFreq) {
		return nate_scalar(freq, lenFreq, rare, lenRare, out);
	}
	for (; rare < stopRare; ++rare) {
		const uint32_t matchRare = *rare;            	//nextRare;
		const vec Match = _mm256_set1_epi32(matchRare);

		if (freq[veclen * 31 + vecmax] < matchRare) { // if no match possible
			uint32_t offset = 1;
			if (freq + veclen * 32 > stopFreq) {
				freq += veclen * 32;
				goto FINISH_SCALAR;
			}
			while (freq[veclen * offset * 32 + veclen * 31 + vecmax]
					< matchRare) { // if no match possible
				if (freq + veclen * (2 * offset ) * 32 <= stopFreq) {
					offset *= 2;
				} else if (freq + veclen * (offset + 1) * 32 <= stopFreq) {
					offset = static_cast<uint32_t>((stopFreq - freq ) / (veclen * 32));
					//offset += 1;
					if (freq[veclen * offset * 32 + veclen * 31 + vecmax]
							< matchRare) {
						freq += veclen * offset * 32;
						goto FINISH_SCALAR;
					} else {
						break;
					}
				} else {
					freq += veclen * offset * 32;
					goto FINISH_SCALAR;
				}
			}
			uint32_t lower = offset / 2;
			while (lower + 1 != offset) {
				const uint32_t mid = (lower + offset) / 2;
				if (freq[veclen * mid * 32 + veclen * 31 + vecmax]
						< matchRare)
				lower = mid;
				else
				offset = mid;
			}
			freq += veclen * offset * 32;
		}
		vec Q0,Q1,Q2,Q3;
		if (freq[veclen * 15 + vecmax] >= matchRare) {
			if (freq[veclen * 7 + vecmax] < matchRare) {
				Q0
				= _mm256_or_si256(
						_mm256_cmpeq_epi32(
								_mm256_loadu_si256((vec *) freq + 8), Match),
						_mm256_cmpeq_epi32(
								_mm256_loadu_si256((vec *) freq + 9), Match));
				Q1 = _mm256_or_si256(
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 10),
								Match),
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 11),
								Match));

				Q2 = _mm256_or_si256(
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 12),
								Match),
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 13),
								Match));
				Q3 = _mm256_or_si256(
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 14),
								Match),
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 15),
								Match));
			} else {
				Q0
				= _mm256_or_si256(
						_mm256_cmpeq_epi32(
								_mm256_loadu_si256((vec *) freq + 4), Match),
						_mm256_cmpeq_epi32(
								_mm256_loadu_si256((vec *) freq + 5), Match));
				Q1
				= _mm256_or_si256(
						_mm256_cmpeq_epi32(
								_mm256_loadu_si256((vec *) freq + 6), Match),
						_mm256_cmpeq_epi32(
								_mm256_loadu_si256((vec *) freq + 7), Match));
				Q2
				= _mm256_or_si256(
						_mm256_cmpeq_epi32(
								_mm256_loadu_si256((vec *) freq + 0), Match),
						_mm256_cmpeq_epi32(
								_mm256_loadu_si256((vec *) freq + 1), Match));
				Q3
				= _mm256_or_si256(
						_mm256_cmpeq_epi32(
								_mm256_loadu_si256((vec *) freq + 2), Match),
						_mm256_cmpeq_epi32(
								_mm256_loadu_si256((vec *) freq + 3), Match));
			}
		} else {
			if (freq[veclen * 23 + vecmax] < matchRare) {
				Q0 = _mm256_or_si256(
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 8 + 16),
								Match),
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 9 + 16),
								Match));
				Q1 = _mm256_or_si256(
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 10 + 16),
								Match),
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 11 + 16),
								Match));

				Q2 = _mm256_or_si256(
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 12 + 16),
								Match),
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 13 + 16),
								Match));
				Q3 = _mm256_or_si256(
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 14 + 16),
								Match),
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 15 + 16),
								Match));
			} else {
				Q0 = _mm256_or_si256(
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 4 + 16),
								Match),
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 5 + 16),
								Match));
				Q1 = _mm256_or_si256(
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 6 + 16),
								Match),
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 7 + 16),
								Match));
				Q2 = _mm256_or_si256(
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 0 + 16),
								Match),
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 1 + 16),
								Match));
				Q3 = _mm256_or_si256(
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 2 + 16),
								Match),
						_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 3 + 16),
								Match));
			}

		}
		const vec F0 = _mm256_or_si256(_mm256_or_si256(Q0, Q1),_mm256_or_si256(Q2, Q3));
		if (_mm256_testz_si256(F0, F0)) {
		} else {
			*out++ = matchRare;
		}
	}

	FINISH_SCALAR: return (out - initout) + nate_scalar(freq,
			stopFreq + freqspace - freq, rare, stopRare + rarespace - rare, out);
}

#endif
size_t SIMDgalloping2(const uint32_t *rare, const size_t lenRare,
		const uint32_t *freq, const size_t lenFreq, uint32_t * out) {
	if (lenFreq == 0 || lenRare == 0)
		return 0;
	assert(lenRare <= lenFreq);
	const uint32_t * const initout(out);

	const size_t freqspace = 8 * VECLEN;
	const size_t rarespace = 1;

	const uint32_t *stopFreq = freq + lenFreq - freqspace;
	const uint32_t *stopRare = rare + lenRare - rarespace;
	if (freq > stopFreq) {
		return nate_scalar(freq, lenFreq, rare, lenRare, out);
	}
	for (; rare < stopRare; ++rare) {
		const uint32_t matchRare = *rare; //nextRare;
		const VEC Match = _mm_set1_epi32(matchRare);

		if (freq[VECLEN * 7 + VECMAX] < matchRare) { // if no match possible
			uint32_t offset = 1;
			if (freq + VECLEN * 8 > stopFreq) {
				freq += VECLEN * 8;
				goto FINISH_SCALAR;
			}
			while (freq[VECLEN * offset * 8 + VECLEN * 7 + VECMAX] < matchRare) { // if no match possible
				if (freq + VECLEN * (2 * offset) * 8 <= stopFreq) {
					offset *= 2;
				} else if (freq + VECLEN * (offset + 1) * 8 <= stopFreq) {
					offset += 1;
				} else {
					freq += VECLEN * offset * 8;
					goto FINISH_SCALAR;
				}
			}
			uint32_t lower = offset / 2;
			while (lower + 1 != offset) {
				const uint32_t mid = (lower + offset) / 2;
				if (freq[VECLEN * mid * 8 + VECLEN * 7 + VECMAX] < matchRare)
					lower = mid;
				else
					offset = mid;
			}
			freq += VECLEN * offset * 8;
		}
		VEC Q0, Q1, Q2, Q3;
		Q0 = _mm_or_si128(
				_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 4), Match),
				_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 5), Match));
		Q1 = _mm_or_si128(
				_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 6), Match),
				_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 7), Match));
		Q2 = _mm_or_si128(
				_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 0), Match),
				_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 1), Match));
		Q3 = _mm_or_si128(
				_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 2), Match),
				_mm_cmpeq_epi32(_mm_load_si128((VEC *) freq + 3), Match));
		const VEC F0 = _mm_or_si128(_mm_or_si128(Q0, Q1), _mm_or_si128(Q2, Q3));
		if (_mm_testz_si128(F0, F0)) {
		} else {
			*out++ = matchRare;
		}
	}

	FINISH_SCALAR: return (out - initout)
			+ nate_scalar(freq, stopFreq + freqspace - freq, rare,
					stopRare + rarespace - rare, out);
}

#endif /* MEDIUMINTERSECTION_H_ */
