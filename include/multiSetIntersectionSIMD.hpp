/*
 * multiSetIntersectionSIMD.hpp
 *
 *  Created on: 2017Äê2ÔÂ5ÈÕ
 *      Author: John
 */

#ifndef INCLUDE_MULTISETINTERSECTIONSIMD_HPP_
#define INCLUDE_MULTISETINTERSECTIONSIMD_HPP_

#include "common.h"
#include "thomaswu.h"
#include "timer.h"

namespace msis/*MultiSet InterSection*/{

void testSIMD() {
	uint32_t array1[128];
	for (int i = 0; i < 128; i++)
		array1[i] = i;
	uint32_t match = 50;
	__m128i m, Match = _mm_set1_epi32(match), Result;
	for (int i = 0; i < 32; i++) {
		m = _mm_loadu_si128((__m128i *) array1 + i);
		Result = _mm_cmpeq_epi32(m, Match);
//		if (_mm_testz_si128(Result, Result)) {
////			std::cout << "i: " << std::dec << i * 4 << " not found: ";
////			int mask = _mm_movemask_ps((__m128 ) Result);
////			std::cout.setf(std::ios::showbase);
////			std::cout << std::hex << mask << std::endl;
////			std::cout.unsetf(std::ios::showbase);
//		} else {
//			std::cout << "i: " << std::dec << i * 4 << " found: ";
//			int mask = _mm_movemask_ps((__m128 ) Result);
//			std::cout.setf(std::ios::showbase);
//			std::cout << std::hex << mask << std::endl;
//			std::cout.unsetf(std::ios::showbase);
//		}
	}

	WallClockTimer timer;
	int bit = 8;
	while (true) {
		std::cout << bit << std::endl;

		Result = _mm_set1_epi32(bit);

		timer.reset();
		for (int i = 0; i < 1000000000; i++)
			_mm_testz_si128(Result, Result);
		size_t time = timer.split();
		std::cout << time << std::endl;

		timer.reset();
		for (int i = 0; i < 1000000000; i++)
			_mm_movemask_ps((__m128 ) Result);
		time = timer.split();
		std::cout << time << std::endl;

		bit >>= 1;
	}
}

template<intersectionfindfunction FINDFUNCTION>
void set_vs_set_SIMD(const mySet &sets, std::vector<uint32_t> &out) {
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

template<intersectionfindfunction FINDFUNCTION>
void swapping_set_vs_set_SIMD(const mySet& sets, std::vector<uint32_t>& out) {
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

template<intersectionfindfunction FINDFUNCTION>
void max_SIMD(const mySet &sets, std::vector<uint32_t> &out) {
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

template<intersectionfindfunction FINDFUNCTION>
void sequential_SIMD(const mySet &sets, std::vector<uint32_t> &out) {
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
}

#endif /* INCLUDE_MULTISETINTERSECTIONSIMD_HPP_ */
