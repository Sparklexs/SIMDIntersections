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

template<flaggedintersectionfindfunction FINDFUNCTION>
void svs_rough(const mySet &sets, std::vector<uint32_t> &out) {
	mySet::iterator it = sets.begin();
	// XXX: we'd like to use rvalue reference, however, it has conflicts
	// with "const", and it finally become an copy operation.
	std::vector<uint32_t> intersection(std::move(*it++));
	int foundp;

	for (; it != sets.end(); it++) {
		auto eliminator = intersection.begin();
		size_t index = 0;
		size_t inter = 0;
		while (eliminator != intersection.end()) {
			index += FINDFUNCTION(&foundp, *eliminator, it->data() + index,
					it->size() - 1 - index);
			if (foundp == 1) {
				intersection[inter++] = *eliminator;
				if (++index == it->size())
					break;
			}
			eliminator++;
		}
		intersection.resize(inter);
	}
	out.swap(intersection);
}

void svs_opt(const mySet &sets, std::vector<uint32_t> &out) {
	mySet::iterator it = sets.begin();
// XXX: we'd like to use rvalue reference, however, it has conflicts
// with "const", and it finally become an copy operation.
	std::vector<uint32_t> intersection(std::move(*it++));
//	size_t search = 0;
	int foundp;

	for (; it != sets.end(); it++) {
		msis::searchFUNC vec_search = optSearchFunc[32
				- __builtin_clz(it->size() / intersection.size())];

		auto eliminator = intersection.begin();
		size_t index = 0;
		size_t inter = 0;
		while (eliminator != intersection.end()) {
			index += vec_search(&foundp, *eliminator, it->data() + index,
					it->size() - 1 - index);
//			search++;
			if (foundp == 1) {
				intersection[inter++] = *eliminator;
				if (++index == it->size())
					break;
			}
			eliminator++;
		}
		intersection.resize(inter);
	}
	out.swap(intersection);
//	std::cout << "opt," << search << std::endl;
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
//	size_t search = 0, update = 0;

	while (eliminator != candidate.end()) {
		index[currentset] += FINDFUNCTION(*eliminator,
				it->data() + index[currentset],
				it->size() - 1 - index[currentset]);
//		search++;
		if (it->at(index[currentset]) == *eliminator) {
			index[currentset]++;
			it++;
			if (++currentset == sets.size() - 1) {
				*value++ = *eliminator;
				count++;
				// roll over
				currentset = 0;
				it = it_start;
				eliminator++;
			}
		} else {
			currentset = 0;
//			update++;
			it = it_start;
			eliminator++;
		}
	}
//	std::cout << "(" << search << ", " << update + count << ")  ";
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
//	size_t search = 0, update = 0;

	while (eliminator != candidate.end()) {
		index[currentset] += FINDFUNCTION(&foundp, *eliminator,
				it->data() + index[currentset],
				it->size() - 1 - index[currentset]);
//		search++;
		if (foundp == 1) {
			++index[currentset];
			it++;
			if (++currentset == sets.size() - 1) {
				*value++ = *eliminator;
				count++;
				// roll over
				currentset = 0;
				it = it_start;
				eliminator++;
			}
		} else {
			if (_UNLIKELY(index[currentset] == it->size()))
				break;
			currentset = 0;
//			update++;
			it = it_start;
			eliminator++;
		}
	}
//	std::cout << "(" << search << ", " << update + count << ")  ";
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
	*value = sets.begin()->at(0);
	size_t count = 0, currentset = 1;
//	size_t search = 0, update = 0;

	while (vsets[0].second != vsets[0].first->size()) {
		vsets[currentset].second += FINDFUNCTION(
				vsets[0].first->at(vsets[0].second),
				vsets[currentset].first->data() + vsets[currentset].second,
				vsets[currentset].first->size() - 1 - vsets[currentset].second);
//		search++;
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
//			update++;
			sort_remaining();
			currentset = 1;
		}
	}
//	std::cout << "(" << search << ", " << update + count << ")  ";
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
//	size_t search = 0, update = 0;

	while (vsets[0].second != vsets[0].first->size()) {
		vsets[currentset].second += FINDFUNCTION(&foundp,
				vsets[0].first->at(vsets[0].second),
				vsets[currentset].first->data() + vsets[currentset].second,
				vsets[currentset].first->size() - 1 - vsets[currentset].second);
//		search++;
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
//			update++;
			sort_remaining();
			currentset = 1;
		}
	}
//	std::cout << "(" << search << ", " << update + count << ")  ";
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
		} else if (_UNLIKELY(vsets[currentset].first->back() < *value))
			break;
		else if (_UNLIKELY(++vsets[0].second == vsets[0].first->size()))
			break;
		else {
//			vsets[0].second++;
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
//	size_t search = 0, update = 0;

	while (_LIKELY(index[currentset] < it->size())) {
		index[currentset] += FINDFUNCTION(*value,
				it->data() + index[currentset],
				it->size() - 1 - index[currentset]);
//		search++;
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
//			update++;
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
//	std::cout << "(" << search << ", " << update + count << ")  ";
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
//	size_t search = 0, update = 0;

	while (_LIKELY(index[currentset] < it->size())) {
		index[currentset] += FINDFUNCTION(&foundp, *value,
				it->data() + index[currentset],
				it->size() - 1 - index[currentset]);
//		search++;

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
//			update++;
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
//	std::cout << "max," << search
//			<< std::endl /*<< ", " << update + count << ")  "*/;
	out.resize(count);
}

void max_rough_opt(const mySet &sets, std::vector<uint32_t> &out) {
	out.resize(sets.begin()->size());
	auto value = out.begin();
	*value = sets.begin()->at(0);
	size_t count = 0, intersect_count = 0, elimset = 0, currentset = 1;
	mySet::iterator it = sets.begin();
	mySet::iterator it_start = it++, it_elim = it_start;
	std::vector<size_t> index(sets.size(), 0);
	int foundp;

	std::vector<msis::searchFUNC> vec_searches(sets.size());
	for (auto set : sets) {
		vec_searches[count++] = optSearchFunc[32
				- __builtin_clz(set.size() / sets.begin()->size())];
	}
	count = 0;

	while (_LIKELY(index[currentset] < it->size())) {
		index[currentset] += vec_searches[currentset](&foundp, *value,
				it->data() + index[currentset],
				it->size() - 1 - index[currentset]);

		if (foundp == 1) {
			intersect_count++;
			index[currentset]++;
			++currentset;
			++it;
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
void BYintersect_sorted(const uint32_t *freq, const size_t &freq_size,
		const uint32_t *rare, const size_t &rare_size, uint32_t **out,
		uint32_t &count) {
	if (_UNLIKELY(
			rare_size == 0 || freq[0] > rare[rare_size - 1]
					|| rare[0] > freq[freq_size - 1]))
		return;
	size_t rare_mid = rare_size / 2;
	if (_UNLIKELY(freq[freq_size - 1] < rare[rare_mid])) {
		// here freq_mid (actually freq_end) certainly less than rare_mid
		BYintersect_sorted<FINDFUNCTION>(rare, rare_mid, freq, freq_size, out,
				count);
		return;
	}
	size_t freq_mid = FINDFUNCTION(rare[rare_mid], freq, freq_size - 1);
	if (freq_mid > rare_mid)
		BYintersect_sorted<FINDFUNCTION>(freq, freq_mid, rare, rare_mid, out,
				count);
	else
		BYintersect_sorted<FINDFUNCTION>(rare, rare_mid, freq, freq_mid, out,
				count);
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
	if (freq_size - freq_mid > rare_size - rare_mid)
		BYintersect_sorted<FINDFUNCTION>(freq + freq_mid, freq_size - freq_mid,
				rare + rare_mid, rare_size - rare_mid, out, count);
	else
		BYintersect_sorted<FINDFUNCTION>(rare + rare_mid, rare_size - rare_mid,
				freq + freq_mid, freq_size - freq_mid, out, count);
//		}
}

template<intersectionfindfunction FINDFUNCTION>
void BY(const mySet &sets, std::vector<uint32_t> &out) {
	mySet::iterator it = sets.begin();
	std::vector<uint32_t> intersection(std::move(*it++));

	uint32_t* out_init = intersection.data();
	uint32_t** pout = &out_init;

	for (; it != sets.end(); it++) {
		uint32_t count = 0;
		BYintersect_sorted<FINDFUNCTION>(it->data(), it->size(),
				intersection.data(), intersection.size(), pout, count);
		intersection.resize(count);

		out_init = intersection.data();
		pout = &out_init;
	}
	out.swap(intersection);
}

/* deprecated */
template<intersectionfindfunction FINDFUNCTION>
void BYintersect_holistic(
		std::vector<std::pair<const uint32_t *, size_t /*size*/>> vsets,
		uint32_t **out, uint32_t &count) {
	if (_UNLIKELY(vsets[0].second == 0))
		return;

	size_t list_mids[vsets.size()];
	list_mids[0] = vsets[0].second / 2;

	bool isdivided = true;
	for (uint32_t i = 1; i < vsets.size(); i++) {
		if (_UNLIKELY(
				vsets[0].first[0] > vsets[i].first[vsets[i].second - 1]
						|| vsets[0].first[vsets[0].second - 1]
								< vsets[i].first[0]))
			return;
		// split from the end
		if (_UNLIKELY(
				vsets[i].first[vsets[i].second - 1]
						< vsets[0].first[list_mids[0]])) {
			list_mids[i] = vsets[i].second;
			isdivided = false;
			continue;
		}
		list_mids[i] = FINDFUNCTION(vsets[0].first[list_mids[0]],
				vsets[i].first, vsets[i].second - 1);
	}
	// first half
	std::vector<std::pair<const uint32_t *, size_t /*end*/>> firsthalf;
	for (uint32_t i = 0; i < vsets.size(); i++) {
		firsthalf.emplace_back(vsets[i].first, list_mids[i]);
	}
	std::sort(firsthalf.begin(), firsthalf.end(),
			[](const std::pair<const uint32_t *, size_t>& lhs,
					const std::pair<const uint32_t *, size_t>& rhs) {
				return lhs.second < rhs.second;});
	BYintersect_holistic<FINDFUNCTION>(firsthalf, out, count);

	// middle and second half
	if (isdivided) {
		bool found = true;
		for (uint32_t i = 1; i < vsets.size(); i++) {
			if (vsets[i].first[list_mids[i]] != vsets[0].first[list_mids[0]]) {
				found = false;
				continue;
			}
			list_mids[i]++;
		}
		if (found) {
			*(*out)++ = vsets[0].first[list_mids[0]];
			count++;
		}
		list_mids[0]++;

		for (uint32_t i = 0; i < vsets.size(); i++) {
			vsets[i].first += list_mids[i];
			vsets[i].second -= list_mids[i];
		}
		std::sort(vsets.begin(), vsets.end(),
				[](const std::pair<const uint32_t *, size_t>& lhs,
						const std::pair<const uint32_t *, size_t>& rhs) {
					return lhs.second < rhs.second;});
		BYintersect_holistic<FINDFUNCTION>(vsets, out, count);
	}
}

/* incomplete! */
template<intersectionfindfunction FINDFUNCTION>
void BYintersect_holistic_new(
		std::vector<std::pair<const uint32_t *, size_t /*size*/>> vsets,
		uint32_t **out, uint32_t &count) {
	if (_UNLIKELY(vsets[0].second == 0))
		return;

	size_t list_mids[vsets.size()];
	list_mids[0] = vsets[0].second / 2;
	for (size_t i = 1; i < vsets.size(); i++)
		list_mids[i] = 0;

	std::vector<std::pair<size_t, size_t /*size*/>> boundaries(vsets.size());
	boundaries[0].first = vsets[0].second / 2 - 1;
	boundaries[0].second = vsets[0].second / 2 + 1;
	for (size_t i = 1; i < vsets.size(); i++) {
		boundaries[i].first = vsets[i].second - 1;
		boundaries[i].second = 0;
	}

	uint32_t value = vsets[0].first[list_mids[0]];
	uint32_t currentset = 1, elimset = 0;
	bool left = true, right = true;

	while (currentset < vsets.size()) {
		if (_UNLIKELY(
				vsets[currentset].first[list_mids[currentset]]
						== vsets[elimset].first[list_mids[elimset]])) {
			if (boundaries[currentset].first > list_mids[currentset] - 1)
				boundaries[currentset].first = list_mids[currentset] - 1;
			if (boundaries[currentset].second < list_mids[currentset] + 1)
				boundaries[currentset].second = list_mids[currentset] + 1;
			currentset++;
			if (currentset == elimset)
				currentset++;
		} else if (vsets[currentset].first[list_mids[currentset]]
				< vsets[elimset].first[list_mids[elimset]]) {
			//TODO: boundary is negative
			if (boundaries[currentset].second >= vsets[currentset].second - 1) {
				right = false;
				break;
			}
			list_mids[currentset] = FINDFUNCTION(
					vsets[elimset].first[list_mids[elimset]],
					vsets[currentset].first + boundaries[currentset].second,
					vsets[currentset].second - 1);
			boundaries[currentset].second = list_mids[currentset] + 1;
			if (boundaries[currentset].first > list_mids[currentset] - 1)
				boundaries[currentset].first = list_mids[currentset] - 1;
			if (vsets[currentset].first[list_mids[currentset]]
					== vsets[elimset].first[list_mids[elimset]]) {
				currentset++;
				if (currentset == elimset)
					currentset++;
			} else {
				elimset = currentset;
				currentset = 0;
			}
		} else {
			//TODO: boundary is negative
			if (boundaries[currentset].first <= 0) {
				left = false;
				break;
			}
			list_mids[currentset] = FINDFUNCTION(
					vsets[elimset].first[list_mids[elimset]],
					vsets[currentset].first, boundaries[currentset].first);
			boundaries[currentset].first = list_mids[currentset] - 1;
			if (boundaries[currentset].second < list_mids[currentset] + 1)
				boundaries[currentset].second = list_mids[currentset] + 1;
			if (vsets[currentset].first[list_mids[currentset]]
					== vsets[elimset].first[list_mids[elimset]]) {
				currentset++;
				if (currentset == elimset)
					currentset++;
			} else {
				elimset = currentset;
				currentset = 0;
			}
		}
	}
	if (left) {
//
	}
	if (currentset = vsets.size()) {
		*(*out)++ = vsets[0].first[list_mids[0]];
		count++;
	}
	if (right) {
//
	}

}

template<intersectionfindfunction FINDFUNCTION>
void BYH(const mySet &sets, std::vector<uint32_t> &out) {
	out.resize(sets.begin()->size());
	std::vector<std::pair<const uint32_t *, size_t /*end*/>> vsets;
	uint32_t count = 0;
	uint32_t* out_init = out.data();
	uint32_t** pout = &out_init;

	for (auto it = sets.begin(); it != sets.end(); it++)
		vsets.emplace_back(it->data(), it->size());
	BYintersect_holistic<FINDFUNCTION>(vsets, pout, count);
	out.resize(count);
}
}
#endif /* INCLUDE_MSIS_METHOD_HPP_ */
