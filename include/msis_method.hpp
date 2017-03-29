/*
 * multiSetIntersectionSIMD.hpp
 *
 *  Created on: 2017Äê2ÔÂ5ÈÕ
 *      Author: John
 */

#ifndef INCLUDE_MSIS_METHOD_HPP_
#define INCLUDE_MSIS_METHOD_HPP_

#include "common.h"
#include "thomaswu.h"
#include "intersectionfactory.h"
#include "msis_gallop_exact.hpp"
#include "msis_gallop_rough.hpp"

namespace msis/*MultiSet InterSection*/{
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
//	size_t invokes = 0;

	while (_LIKELY(index[currentset] < it->size())) {
		index[currentset] += FINDFUNCTION(&foundp, *value,
				it->data() + index[currentset],
				it->size() - 1 - index[currentset]);
//		invokes++;

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
//	std::cout << invokes << "  ";
	out.resize(count);
}

/*
 * @description: basic methods intersecting two sorted sequences,
 * returned results are also sorted.
 */
// further optimization can be search @freq[0] and @freq[freq_end] in @rare,
// which helps target the real overlap of @rare and @freq
template<intersectionfindfunction FINDFUNCTION>
void BYintersect_sorted_exact(const uint32_t *freq, const size_t &freq_end,
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
			BYintersect_sorted_exact<FINDFUNCTION>(freq, freq_mid - 1, rare,
					rare_mid - 1, out, count);
		else
			BYintersect_sorted_exact<FINDFUNCTION>(rare, rare_mid - 1, freq,
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
			BYintersect_sorted_exact<FINDFUNCTION>(freq + freq_mid,
					freq_end - freq_mid, rare + rare_mid, rare_end - rare_mid,
					out, count);
		else
			BYintersect_sorted_exact<FINDFUNCTION>(rare + rare_mid,
					rare_end - rare_mid, freq + freq_mid, freq_end - freq_mid,
					out, count);
//		}
	}
}

template<flaggedintersectionfindfunction FINDFUNCTION>
void BYintersect_sorted_rough(const uint32_t *freq, const size_t &freq_end,
		const uint32_t *rare, const size_t &rare_end, uint32_t **out,
		uint32_t &count) {
	if (_LIKELY(
			freq_end == -1 || rare_end == -1 || freq[0] > rare[rare_end]
					|| rare[0] > freq[freq_end]))
		return;
	else {
		size_t rare_mid = rare_end / 2;
		int foundp;
		size_t freq_mid = FINDFUNCTION(&foundp, rare[rare_mid], freq,
				freq_end + 1);
		if (freq_mid > rare_mid)
			BYintersect_sorted_rough<FINDFUNCTION>(freq, freq_mid - 1, rare,
					rare_mid - 1, out, count);
		else
			BYintersect_sorted_rough<FINDFUNCTION>(rare, rare_mid - 1, freq,
					freq_mid - 1, out, count);

		if (foundp == 1) {
			*(*out)++ = freq[freq_mid++];
			count++;
		}

		if (freq_end - freq_mid > rare_end - rare_mid)
			BYintersect_sorted_rough<FINDFUNCTION>(freq + freq_mid,
					freq_end - freq_mid, rare + rare_mid, rare_end - rare_mid,
					out, count);
		else
			BYintersect_sorted_rough<FINDFUNCTION>(rare + rare_mid,
					rare_end - rare_mid, freq + freq_mid, freq_end - freq_mid,
					out, count);
//		}
	}
}

template<intersectionfindfunction FINDFUNCTION>
void BY_exact(const mySet &sets, std::vector<uint32_t> &out) {
	mySet::iterator it = sets.begin();
	std::vector<uint32_t> intersection(std::move(*it++));

	uint32_t* out_init = intersection.data();
	uint32_t** pout = &out_init;

	for (; it != sets.end(); it++) {
		uint32_t count = 0;

		BYintersect_sorted_exact<FINDFUNCTION>(it->data(), it->size() - 1,
				intersection.data(), intersection.size() - 1, pout, count);
		intersection.resize(count);

		out_init = intersection.data();
		pout = &out_init;
	}
	out.swap(intersection);
}

template<flaggedintersectionfindfunction FINDFUNCTION>
void BY_rough(const mySet &sets, std::vector<uint32_t> &out) {
	mySet::iterator it = sets.begin();
	std::vector<uint32_t> intersection(std::move(*it++));

	uint32_t* out_init = intersection.data();
	uint32_t** pout = &out_init;

	for (; it != sets.end(); it++) {
		uint32_t count = 0;

		BYintersect_sorted_rough<FINDFUNCTION>(it->data(), it->size() - 1,
				intersection.data(), intersection.size() - 1, pout, count);
		intersection.resize(count);

		out_init = intersection.data();
		pout = &out_init;
	}
	out.swap(intersection);
}

}
#endif /* INCLUDE_MSIS_METHOD_HPP_ */
