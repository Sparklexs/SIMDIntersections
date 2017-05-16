/*
 * multiSetIntersection.cpp
 * scalar algorithms for intersecting multiple sets.
 * Reference:
 * (1)  Culpepper, J. S., & Moffat, A. (2010).
 *      Efficient set intersection for inverted indexing.
 *      Acm Transactions on Information Systems, 29(1), 1.
 *
 * (2)  Barbay, J., L��pez-Ortiz, A., Lu, T., & Salinger, A.
 *              (2009). An experimental investigation of set
 *              intersection algorithms for text searching.
 *              Journal of Experimental Algorithmics, 14.
 *
 * (3)  Baezayates, R., & Salinger, A. (2010).
 *              Fast Intersection Algorithms for Sorted Sequences.
 *              Algorithms and Applications. Springer Berlin Heidelberg.
 *
 * (4)  Barbay, J., L��pezortiz, A., & Lu, T. (2006).
 *              Faster Adaptive Set Intersections for Text Searching.
 *              International Conference on Experimental Algorithms
 *              (Vol.4007, pp.146-157). Springer-Verlag.
 *
 *  Created on: 2016/12/20
 *      Author: SparkleXS
 */

#include "timer.h"
#include "synthetic.h"
#include <boost/preprocessor.hpp>
#include <msis_method.hpp>

class Set_collection {
	const static uint32_t HEADERSIZE = sizeof(float) + 3 * sizeof(uint32_t);
	const static uint32_t REPETITION = 20;
	const static uint32_t NUM = 10;
public:
	static mySet load_set(float _ir, uint32_t _sr, uint32_t _num,
			uint32_t _rank) {
		std::ifstream infile("sets.dat", std::ios::binary);
		if (!infile) {
			vector<float> intersectionsratios = { 0.10, 0.50, 1.00 };
			vector<uint32_t> sizeratios = { 1, 10, 100, 1000, 10000 };
			uint32_t logMinLength = 12; // log of minimal array size
			uint32_t ranks = REPETITION;
			uint32_t MaxBit = 31; // largest bit-length of element
			const uint32_t minlength = 1U << logMinLength;
			ClusteredDataGenerator cdg;

			std::ofstream outfile("sets.dat", std::ios::binary);
			for (float ir : intersectionsratios) {

				outfile.write((char*) &ir, sizeof(float));
				std::cout << "selectivity: " << ir << std::endl;

				for (uint32_t sr : sizeratios) {
					std::cout << "  size ratio: " << sr << std::endl;
					outfile.write((char*) &sr, sizeof(uint32_t));
					outfile.write((char*) &ranks, sizeof(uint32_t));

					for (uint32_t i = 2; i < 11; i++) {
						std::cout << "    num: " << i << std::endl;
						outfile.write((char*) &i, sizeof(uint32_t));
						long pos = outfile.tellp();

						outfile.write((char*) &i, sizeof(size_t)); // reserve
						size_t length = 0;

						for (uint32_t c = 0; c < ranks; c++) {

							mySet && sets = genMultipleSets(cdg, minlength, i,
									1U << MaxBit, static_cast<float>(sr), ir);
							for (auto set : sets) {

								uint32_t size = set.size();
								outfile.write((char*) &size, sizeof(uint32_t));
								outfile.write((char*) set.data(),
										set.size() * sizeof(uint32_t));
								length += 4 + set.size() * 4;
							}
						}
						outfile.seekp(pos);
						outfile.write((char*) &length, sizeof(size_t));
						outfile.seekp(0, std::ios_base::end);
					} //num
				} //sr
			} //ir
			outfile.flush();
			outfile.close();

			return load_set(_ir, _sr, _num, _rank);
		} else {
			size_t totallength;
			float selectivity;
			uint32_t ratio = 10000;
			uint32_t num = NUM;
			mySet multiset([](vector<uint32_t> lhs,vector<uint32_t> rhs) {
				return lhs.size() <= rhs.size();
			});

//			infile.seekg(pos);
			do {
				if (num == NUM) {
					if (ratio == 10000) {
						infile.read((char*) &selectivity, sizeof(float));
						infile.read((char*) &ratio, sizeof(uint32_t));
						infile.ignore(4); // total ranks must equal 20
						//infile.read((char*) &rank, sizeof(uint32_t));
						infile.read((char*) &num, sizeof(uint32_t));
						infile.read((char*) &totallength, sizeof(size_t));
					} else {
						infile.read((char*) &ratio, sizeof(uint32_t));
						infile.ignore(4); // total ranks must equal 20
						//infile.read((char*) &rank, sizeof(uint32_t));
						infile.read((char*) &num, sizeof(uint32_t));
						infile.read((char*) &totallength, sizeof(size_t));
					}
				} else {
					infile.read((char*) &num, sizeof(uint32_t));
					infile.read((char*) &totallength, sizeof(size_t));
				}

				infile.seekg(totallength, infile.cur);
			} while (!(selectivity == _ir and ratio == _sr and num == _num));

			infile.seekg(-totallength, infile.cur);
			uint32_t length;
			uint32_t rank = 0;

			while (rank != _rank) {
				for (uint32_t i = 0; i < num; i++) {
					infile.read((char*) &length, sizeof(uint32_t));
					infile.seekg(length * 4, infile.cur);
				}
				rank++;
			}
			for (uint32_t i = 0; i < num; i++) {
				infile.read((char*) &length, sizeof(uint32_t));
				vector<uint32_t> set(length);
				uint32_t* ptr = set.data();

				infile.read((char*) ptr, length * sizeof(uint32_t));
				multiset.insert(set);
			}
			return multiset;
		} // else
	} // method

};

void Zipper(const mySet &sets, std::vector<uint32_t> &out) {
	msis::svs<branchlessintersection>(sets, out);
}

void Schlegel(mySet &multiset, WallClockTimer timer,
		std::map<std::string, size_t>& times) {
	vector<vector<uint16_t> > pdata(multiset.size());
	mySet::iterator it = multiset.begin();

	for (uint32_t i = 0; i < multiset.size(); i++) {
		pdata[i].resize(it->size() * 4);
		const size_t c = partitioned::partition(it->data(), it->size(),
				&pdata[i][0], pdata[i].size());
		pdata[i].resize(c);
		vector<uint16_t>(pdata[i]).swap(pdata[i]);
		it++;
	}

	for (auto algo : partschemes) {
		timer.reset();
		for (uint32_t j = 0; j < 100; j++) {
			vector<uint16_t> intersection(pdata[0]);
			for (uint32_t k = 0; k < multiset.size(); k++) {
				size_t answer = algo.second(intersection.data(),
						pdata[k].data(), intersection.size(), pdata[k].size());
				intersection.resize(answer);
			}
		}
		times[algo.first] += timer.split();
	}
}

void Lemire_V1(const mySet &sets, std::vector<uint32_t> &out) {
	msis::svs<v1>(sets, out);
}

void Lemire_V3(const mySet &sets, std::vector<uint32_t> &out) {
	msis::svs<v3>(sets, out);
}

void Lemire_Gallop(const mySet &sets, std::vector<uint32_t> &out) {
	msis::svs<SIMDgalloping>(sets, out);
}

void Inoue(const mySet &sets, std::vector<uint32_t> &out) {
	msis::svs<SIMDIntersectWithPrefilter>(sets, out);
}

void set_vs_set_scalar(const mySet &sets, std::vector<uint32_t> &out) {
	msis::SvS_exact<msis::scalarBinarySearch>(sets, out);
}

void swapping_set_vs_set_scalar(const mySet& sets, vector<uint32_t>& out) {
	msis::s_SvS_exact<msis::scalarBinarySearch>(sets, out);
}

void sequential_scalar(const mySet &sets, std::vector<uint32_t> &out) {
	msis::sql_exact<msis::scalarGallop>(sets, out);
}

void small_sequential_scalar(const mySet &sets, std::vector<uint32_t> &out) {
	msis::s_sql_exact<msis::scalarGallop>(sets, out);
}

void max_scalar(const mySet &sets, std::vector<uint32_t> &out) {
	msis::max_exact<msis::scalarGallop>(sets, out);
}

void BaezaYates_scalar(const mySet &sets, std::vector<uint32_t> &out) {
	msis::BY_exact<msis::scalarBinarySearch>(sets, out);
}

typedef void (*intersectionFUNC)(const mySet &sets, std::vector<uint32_t> &out);
/* SIMD methods are template functions that can't be applied to typedef */
intersectionFUNC scalarFUNC[] = { Zipper, Lemire_V1, Lemire_V3, Lemire_Gallop,
		Inoue, set_vs_set_scalar, swapping_set_vs_set_scalar, /*msis::adp_scalar,
		 msis::s_adp_scalar,*/sequential_scalar, small_sequential_scalar,
		max_scalar, BaezaYates_scalar };

const std::string NAMESCALARFUNC[] = { "zip", "v1", "v3", "gallop", "inoue",
		"SvS", "s_SvS", /*"adp", "s_adp",*/"sql", "s_sql", "max", "bys" };

constexpr uint32_t NUMSCALARFUNC = sizeof(scalarFUNC) / sizeof(scalarFUNC[0]);

/////////////////////////////////////////
#define METHOD (msis::SvS)(msis::s_SvS)(msis::sql)(msis::s_sql)(msis::max)(msis::BY)(msis::svs)
#define HEAD (_exact<msis::simd)(_rough<msis::simd)
#define SEARCH (linear)(gallop)
#define SIZE (_4)(_8)(_16)(_32)(_64)(_128)(_256)(_512)
#define TAIL (_exact>)(_rough>)(_rough_plow>)
//////////////////////////////////////////

std::map<std::string, size_t> init_time_array() {
	std::map<std::string, size_t> times;
	for (uint32_t i = 0; i < NUMSCALARFUNC; i++) {
		times[NAMESCALARFUNC[i]] = 0;
	}

	for (auto algo : partschemes) {
		times[algo.first] = 0;
	}

	std::string name;
#define LOOP_BODY(R,PRODUCT)\
		name=(BOOST_PP_STRINGIZE(BOOST_PP_SEQ_ELEM(0,PRODUCT)));\
		name.append("_").append(BOOST_PP_STRINGIZE(BOOST_PP_SEQ_ELEM(2,PRODUCT)));\
		name.append(BOOST_PP_STRINGIZE(BOOST_PP_SEQ_ELEM(3,PRODUCT)));\
		name.append(BOOST_PP_STRINGIZE(BOOST_PP_SEQ_ELEM(4,PRODUCT)));\
		name.append(BOOST_PP_STRINGIZE(BOOST_PP_SEQ_ELEM(5,PRODUCT)));\
		name.resize(name.size()-1);\
		name=name.substr(6);\
	    times[name]=0;\

	BOOST_PP_SEQ_FOR_EACH_PRODUCT(LOOP_BODY, ((BOOST_PP_SEQ_ELEM(4,METHOD)))/*0:SvS,1:s_SvS,2:sql,3:s_sql,4:max*/
	((BOOST_PP_SEQ_ELEM(1,HEAD)))/*0:exact 1:rough*/
	((BOOST_PP_SEQ_ELEM(1,SEARCH)))/*0:linear 1:gallop*/
	((BOOST_PP_SEQ_ELEM(0,SIZE))) /*4,8,16,32,64,128,256,512*/
	(BOOST_PP_SEQ_FIRST_N(1,SIZE)) ((BOOST_PP_SEQ_ELEM(1,TAIL)))/*0:exact 1:rough 2:rough_plow*/
	);
	BOOST_PP_SEQ_FOR_EACH_PRODUCT(LOOP_BODY, ((BOOST_PP_SEQ_ELEM(4,METHOD)))/*0:SvS,1:s_SvS,2:sql,3:s_sql,4:max*/
	((BOOST_PP_SEQ_ELEM(1,HEAD)))/*0:exact 1:rough*/
	((BOOST_PP_SEQ_ELEM(1,SEARCH)))/*0:linear 1:gallop*/
	((BOOST_PP_SEQ_ELEM(1,SIZE))) /*4,8,16,32,64,128,256,512*/
	(BOOST_PP_SEQ_FIRST_N(2,SIZE)) ((BOOST_PP_SEQ_ELEM(1,TAIL)))/*0:exact 1:rough 2:rough_plow*/
	);
	BOOST_PP_SEQ_FOR_EACH_PRODUCT(LOOP_BODY, ((BOOST_PP_SEQ_ELEM(4,METHOD)))/*0:SvS,1:s_SvS,2:sql,3:s_sql,4:max*/
	((BOOST_PP_SEQ_ELEM(1,HEAD)))/*0:exact 1:rough*/
	((BOOST_PP_SEQ_ELEM(1,SEARCH)))/*0:linear 1:gallop*/
	((BOOST_PP_SEQ_ELEM(2,SIZE))) /*4,8,16,32,64,128,256,512*/
	(BOOST_PP_SEQ_FIRST_N(3,SIZE)) ((BOOST_PP_SEQ_ELEM(1,TAIL)))/*0:exact 1:rough 2:rough_plow*/
	);
	BOOST_PP_SEQ_FOR_EACH_PRODUCT(LOOP_BODY, ((BOOST_PP_SEQ_ELEM(4,METHOD)))/*0:SvS,1:s_SvS,2:sql,3:s_sql,4:max*/
	((BOOST_PP_SEQ_ELEM(1,HEAD)))/*0:exact 1:rough*/
	((BOOST_PP_SEQ_ELEM(1,SEARCH)))/*0:linear 1:gallop*/
	((BOOST_PP_SEQ_ELEM(3,SIZE))) /*4,8,16,32,64,128,256,512*/
	(BOOST_PP_SEQ_FIRST_N(4,SIZE)) ((BOOST_PP_SEQ_ELEM(1,TAIL)))/*0:exact 1:rough 2:rough_plow*/
	);
	BOOST_PP_SEQ_FOR_EACH_PRODUCT(LOOP_BODY, ((BOOST_PP_SEQ_ELEM(4,METHOD)))/*0:SvS,1:s_SvS,2:sql,3:s_sql,4:max*/
	((BOOST_PP_SEQ_ELEM(1,HEAD)))/*0:exact 1:rough*/
	((BOOST_PP_SEQ_ELEM(1,SEARCH)))/*0:linear 1:gallop*/
	((BOOST_PP_SEQ_ELEM(4,SIZE))) /*4,8,16,32,64,128,256,512*/
	(BOOST_PP_SEQ_FIRST_N(5,SIZE)) ((BOOST_PP_SEQ_ELEM(1,TAIL)))/*0:exact 1:rough 2:rough_plow*/
	);
	BOOST_PP_SEQ_FOR_EACH_PRODUCT(LOOP_BODY, ((BOOST_PP_SEQ_ELEM(4,METHOD)))/*0:SvS,1:s_SvS,2:sql,3:s_sql,4:max*/
	((BOOST_PP_SEQ_ELEM(1,HEAD)))/*0:exact 1:rough*/
	((BOOST_PP_SEQ_ELEM(1,SEARCH)))/*0:linear 1:gallop*/
	((BOOST_PP_SEQ_ELEM(5,SIZE))) /*4,8,16,32,64,128,256,512*/
	(BOOST_PP_SEQ_FIRST_N(6,SIZE)) ((BOOST_PP_SEQ_ELEM(1,TAIL)))/*0:exact 1:rough 2:rough_plow*/
	);
	BOOST_PP_SEQ_FOR_EACH_PRODUCT(LOOP_BODY, ((BOOST_PP_SEQ_ELEM(4,METHOD)))/*0:SvS,1:s_SvS,2:sql,3:s_sql,4:max*/
	((BOOST_PP_SEQ_ELEM(1,HEAD)))/*0:exact 1:rough*/
	((BOOST_PP_SEQ_ELEM(1,SEARCH)))/*0:linear 1:gallop*/
	((BOOST_PP_SEQ_ELEM(6,SIZE))) /*4,8,16,32,64,128,256,512*/
	(BOOST_PP_SEQ_FIRST_N(7,SIZE)) ((BOOST_PP_SEQ_ELEM(1,TAIL)))/*0:exact 1:rough 2:rough_plow*/
	);
	BOOST_PP_SEQ_FOR_EACH_PRODUCT(LOOP_BODY, ((BOOST_PP_SEQ_ELEM(4,METHOD)))/*0:SvS,1:s_SvS,2:sql,3:s_sql,4:max*/
	((BOOST_PP_SEQ_ELEM(1,HEAD)))/*0:exact 1:rough*/
	((BOOST_PP_SEQ_ELEM(1,SEARCH)))/*0:linear 1:gallop*/
	((BOOST_PP_SEQ_ELEM(7,SIZE))) /*4,8,16,32,64,128,256,512*/
	(BOOST_PP_SEQ_FIRST_N(8,SIZE)) ((BOOST_PP_SEQ_ELEM(1,TAIL)))/*0:exact 1:rough 2:rough_plow*/
	);
#undef LOOP_BODY
	return times;
}

void intersect_for_json(bool loadfromfile) {
	using namespace msis;

	vector<uint32_t> out;

	uint32_t MaxBit = 31; // largest bit-length of element
	uint32_t minlength;
	size_t REPETITION = 1000;
	size_t CASES = 20;

	WallClockTimer timer;
	vector<float> intersectionsratios = { 0.10, 0.30, 0.50, 0.70, 0.90 }; // odds
	vector<uint32_t> sizeratios = { 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024,
			2048, 4096, 8192 };

	std::map<std::string, size_t> times = init_time_array();
	for (uint32_t msb = 10; msb <= 13; ++msb) {
		minlength = 1U << msb;
		ostringstream name("time_", std::ios::ate);
		name << msb << ".json";
		ofstream of(name.str(), std::ios::binary);
		for (float ir : intersectionsratios) {
			printf("intersection ratio: \e[32m%3.0f%%\e[0m\n", ir * 100);
			for (uint32_t sr : sizeratios) {
				printf("  size ratio: \e[32m%4d\e[0m\n", sr);
				if (sr > 500)
					REPETITION = 200;
				else
					REPETITION = 1000;
				for (uint32_t num = 2; num < 11; num++) {
					times = init_time_array();
					time_t t = time(nullptr);
					tm* _tm = localtime(&t);
					printf("%02d:%02d:%02d> num: \e[32m%2d\e[0m  ",
							_tm->tm_hour, _tm->tm_min, _tm->tm_sec, num);

					for (uint32_t k = 0; k < CASES; k++) {
						mySet multiset;
						if (loadfromfile)
							multiset = Set_collection::load_set(ir, sr, num, k);
						else {
							ClusteredDataGenerator cdg;
							multiset = genMultipleSets(cdg, minlength, num,
									1U << MaxBit, static_cast<float>(sr), ir);
						}

//						 verification
//						auto it = multiset.begin();
//						vector<uint32_t> final_intersection = intersect(*it++,
//								*it++);
//						for (; it != multiset.end(); it++)
//							final_intersection = intersect(final_intersection,
//									*it);
//
//						msis::max_exact<msis::simdgallop_128_exact>(multiset,
//								out);
//						msis::max_rough<msis::simdgallop_128_4_rough>(multiset,
//								out);
////
////						msis::s_SvS_exact<msis::simdgallop_128_exact>(multiset,
////								out);
////						msis::s_SvS_rough<msis::simdgallop_128_32_rough>(
////								multiset, out);
////
////						msis::SvS_exact<msis::simdgallop_128_exact>(multiset,
////								out);
////						msis::SvS_rough<msis::simdgallop_128_32_rough>(multiset,
////								out);
//						std::cout << std::endl;
//						if (out != final_intersection) {
//							std::cerr << "bad result!  " << std::endl;
//							return;
//						} else
//							printf("good!  ");

						std::string name;
#define LOOP_BODY(R,PRODUCT)\
				timer.reset();\
				for (uint32_t howmany = 0; howmany < REPETITION;\
				++howmany) {\
					BOOST_PP_CAT(\
					BOOST_PP_CAT(\
					BOOST_PP_CAT(\
					BOOST_PP_CAT(\
					BOOST_PP_CAT(\
					BOOST_PP_SEQ_ELEM(0,PRODUCT),\
					BOOST_PP_SEQ_ELEM(1,PRODUCT)),\
					BOOST_PP_SEQ_ELEM(2,PRODUCT)),\
					BOOST_PP_SEQ_ELEM(3,PRODUCT)),\
					BOOST_PP_SEQ_ELEM(4,PRODUCT)),\
					BOOST_PP_SEQ_ELEM(5,PRODUCT))(multiset,out);\
				}\
				name=(BOOST_PP_STRINGIZE(BOOST_PP_SEQ_ELEM(0,PRODUCT)));\
				name.append("_").append(BOOST_PP_STRINGIZE(BOOST_PP_SEQ_ELEM(2,PRODUCT)));\
				name.append(BOOST_PP_STRINGIZE(BOOST_PP_SEQ_ELEM(3,PRODUCT)));\
				name.append(BOOST_PP_STRINGIZE(BOOST_PP_SEQ_ELEM(4,PRODUCT)));\
				name.append(BOOST_PP_STRINGIZE(BOOST_PP_SEQ_ELEM(5,PRODUCT)));\
				name.resize(name.size()-1);\
				name=name.substr(6);\
				times[name] += timer.split();\

						// test rough_gallop_rough
						BOOST_PP_SEQ_FOR_EACH_PRODUCT(LOOP_BODY,
								((BOOST_PP_SEQ_ELEM(6,METHOD)))/*0:SvS,1:s_SvS,2:sql,3:s_sql,4:max,5:BY,6:svs_me*/
								((BOOST_PP_SEQ_ELEM(1,HEAD)))/*0:exact 1:rough*/
								((BOOST_PP_SEQ_ELEM(1,SEARCH)))/*0:linear 1:gallop*/
								((BOOST_PP_SEQ_ELEM(0,SIZE))) /*4,8,16,32,64,128,256,512*/
								(BOOST_PP_SEQ_FIRST_N(1,SIZE)) ((BOOST_PP_SEQ_ELEM(1,TAIL)))/*0:exact 1:rough 2:rough_plow*/
								);
						BOOST_PP_SEQ_FOR_EACH_PRODUCT(LOOP_BODY,
								((BOOST_PP_SEQ_ELEM(6,METHOD)))/*0:SvS,1:s_SvS,2:sql,3:s_sql,4:max*/
								((BOOST_PP_SEQ_ELEM(1,HEAD)))/*0:exact 1:rough*/
								((BOOST_PP_SEQ_ELEM(1,SEARCH)))/*0:linear 1:gallop*/
								((BOOST_PP_SEQ_ELEM(1,SIZE))) /*4,8,16,32,64,128,256,512*/
								(BOOST_PP_SEQ_FIRST_N(2,SIZE)) ((BOOST_PP_SEQ_ELEM(1,TAIL)))/*0:exact 1:rough 2:rough_plow*/
								);
						BOOST_PP_SEQ_FOR_EACH_PRODUCT(LOOP_BODY,
								((BOOST_PP_SEQ_ELEM(6,METHOD)))/*0:SvS,1:s_SvS,2:sql,3:s_sql,4:max*/
								((BOOST_PP_SEQ_ELEM(1,HEAD)))/*0:exact 1:rough*/
								((BOOST_PP_SEQ_ELEM(1,SEARCH)))/*0:linear 1:gallop*/
								((BOOST_PP_SEQ_ELEM(2,SIZE))) /*4,8,16,32,64,128,256,512*/
								(BOOST_PP_SEQ_FIRST_N(3,SIZE)) ((BOOST_PP_SEQ_ELEM(1,TAIL)))/*0:exact 1:rough 2:rough_plow*/
								);
						BOOST_PP_SEQ_FOR_EACH_PRODUCT(LOOP_BODY,
								((BOOST_PP_SEQ_ELEM(6,METHOD)))/*0:SvS,1:s_SvS,2:sql,3:s_sql,4:max*/
								((BOOST_PP_SEQ_ELEM(1,HEAD)))/*0:exact 1:rough*/
								((BOOST_PP_SEQ_ELEM(1,SEARCH)))/*0:linear 1:gallop*/
								((BOOST_PP_SEQ_ELEM(3,SIZE))) /*4,8,16,32,64,128,256,512*/
								(BOOST_PP_SEQ_FIRST_N(4,SIZE)) ((BOOST_PP_SEQ_ELEM(1,TAIL)))/*0:exact 1:rough 2:rough_plow*/
								);
						BOOST_PP_SEQ_FOR_EACH_PRODUCT(LOOP_BODY,
								((BOOST_PP_SEQ_ELEM(6,METHOD)))/*0:SvS,1:s_SvS,2:sql,3:s_sql,4:max*/
								((BOOST_PP_SEQ_ELEM(1,HEAD)))/*0:exact 1:rough*/
								((BOOST_PP_SEQ_ELEM(1,SEARCH)))/*0:linear 1:gallop*/
								((BOOST_PP_SEQ_ELEM(4,SIZE))) /*4,8,16,32,64,128,256,512*/
								(BOOST_PP_SEQ_FIRST_N(5,SIZE)) ((BOOST_PP_SEQ_ELEM(1,TAIL)))/*0:exact 1:rough 2:rough_plow*/
								);
						BOOST_PP_SEQ_FOR_EACH_PRODUCT(LOOP_BODY,
								((BOOST_PP_SEQ_ELEM(6,METHOD)))/*0:SvS,1:s_SvS,2:sql,3:s_sql,4:max*/
								((BOOST_PP_SEQ_ELEM(1,HEAD)))/*0:exact 1:rough*/
								((BOOST_PP_SEQ_ELEM(1,SEARCH)))/*0:linear 1:gallop*/
								((BOOST_PP_SEQ_ELEM(5,SIZE))) /*4,8,16,32,64,128,256,512*/
								(BOOST_PP_SEQ_FIRST_N(6,SIZE)) ((BOOST_PP_SEQ_ELEM(1,TAIL)))/*0:exact 1:rough 2:rough_plow*/
								);
						BOOST_PP_SEQ_FOR_EACH_PRODUCT(LOOP_BODY,
								((BOOST_PP_SEQ_ELEM(6,METHOD)))/*0:SvS,1:s_SvS,2:sql,3:s_sql,4:max*/
								((BOOST_PP_SEQ_ELEM(1,HEAD)))/*0:exact 1:rough*/
								((BOOST_PP_SEQ_ELEM(1,SEARCH)))/*0:linear 1:gallop*/
								((BOOST_PP_SEQ_ELEM(6,SIZE))) /*4,8,16,32,64,128,256,512*/
								(BOOST_PP_SEQ_FIRST_N(7,SIZE)) ((BOOST_PP_SEQ_ELEM(1,TAIL)))/*0:exact 1:rough 2:rough_plow*/
								);
						BOOST_PP_SEQ_FOR_EACH_PRODUCT(LOOP_BODY,
								((BOOST_PP_SEQ_ELEM(6,METHOD)))/*0:SvS,1:s_SvS,2:sql,3:s_sql,4:max*/
								((BOOST_PP_SEQ_ELEM(1,HEAD)))/*0:exact 1:rough*/
								((BOOST_PP_SEQ_ELEM(1,SEARCH)))/*0:linear 1:gallop*/
								((BOOST_PP_SEQ_ELEM(7,SIZE))) /*4,8,16,32,64,128,256,512*/
								(BOOST_PP_SEQ_FIRST_N(8,SIZE)) ((BOOST_PP_SEQ_ELEM(1,TAIL)))/*0:exact 1:rough 2:rough_plow*/
								);
#undef LOOP_BODY
					} // for CASES
					for (auto time : times) {
						if (time.second != 0) {
							uint32_t type = 0;
							if (time.first.find("max") == std::string::npos)
								type = 1;
							of << "{\"type\": " << type << ", ";
							of << "\"time\": "
									<< (double) time.second / REPETITION / CASES
									<< ", ";
							size_t start = time.first.find('_');
							start = time.first.find('_', start + 1);
							size_t end = time.first.find('_', start + 1);
							uint32_t h;
							istringstream istr(
									time.first.substr(start + 1,
											end - start - 1));
							istr >> h;
							of << "\"h_exp\": " << h << ", ";
							h = std::log2(h);
							of << "\"h\": " << h << ", ";

							start = end;
							end = time.first.find('_', start + 1);

							uint32_t k;
							istr.clear();
							istr.str(
									time.first.substr(start + 1,
											end - start - 1));
							istr >> k;
							of << "\"k_exp\": " << k << ", ";
							k = std::log2(k);
							of << "\"k\": " << k << ", ";
							of << "\"diff\": " << h - k << ", ";
							of << "\"sr\": " << sr << ", ";
							of << "\"ir\": " << ir << ", ";
							of << "\"ratio\": "
									<< log2(minlength * sr / pow(2, h)) << "}"
									<< std::endl;
						} // if !=0
//							printf("%s: \e[31m%6.3f\e[0m  \n",
//									time.first.c_str(),
//									(double) time.second / REPETITION / CASES);
					} // for times
					printf("\n");
					fflush(stdout);
				} // for num
			} // for size-ratio
		} // for intersection-ratio
		of.flush();
		of.close();
	} // for minlength
}

void intersect_traditional_methods(bool loadfromfile) {
	using namespace msis;

	vector<uint32_t> out;

	uint32_t MaxBit = 31; // largest bit-length of element
	uint32_t minlength;
	size_t REPETITION = 100;
	size_t CASES = 20;

	WallClockTimer timer;
	vector<float> intersectionsratios = { 0.10, 0.50, 1.00 };
	vector<uint32_t> sizeratios = { 1, 10, 100, 1000, 10000 };

	FILE *pfile = fopen("output_traditional.csv", "w+");
	fprintf(pfile, "sr,ir,num,name,time\n");
	std::map<std::string, size_t> times = init_time_array();
	for (uint32_t msb = 10; msb <= 10; ++msb) {
		minlength = 1U << msb;

		for (float ir : intersectionsratios) {
			printf("intersection ratio: \e[32m%3.0f%%\e[0m\n", ir * 100);
			for (uint32_t sr : sizeratios) {
				printf("  size ratio: \e[32m%4d\e[0m\n", sr);
				if (sr > 100)
					REPETITION = 30;
				else
					REPETITION = 100;
				for (uint32_t num = 2; num < 11; num++) {
					times = init_time_array();
					time_t t = time(nullptr);
					tm* _tm = localtime(&t);
					printf("%02d:%02d:%02d> num: \e[32m%2d\e[0m  ",
							_tm->tm_hour, _tm->tm_min, _tm->tm_sec, num);

					for (uint32_t k = 0; k < CASES; k++) {
						mySet multiset;
						if (loadfromfile)
							multiset = Set_collection::load_set(ir, sr, num, k);
						else {
							ClusteredDataGenerator cdg;
							multiset = genMultipleSets(cdg, minlength, num,
									1U << MaxBit, static_cast<float>(sr), ir);
						}

						// start scalar intersection
						for (uint32_t j = 0; j < NUMSCALARFUNC; j++) {
							timer.reset();
							for (uint32_t howmany = 0; howmany < REPETITION;
									++howmany) {
								scalarFUNC[j](multiset, out);
							}

							times[NAMESCALARFUNC[j]] += timer.split();
						}
						Schlegel(multiset, timer, times);

					} // for CASES
					for (auto time : times) {
						if (time.second != 0) {
							printf("%s: \e[31m%6.0f\e[0m  ", time.first.c_str(),
									(double) time.second / REPETITION / CASES);
							fprintf(pfile, "%d,%.0f,%d,%s,%.0f\n", sr, ir * 100,
									num, time.first.c_str(),
									(double) time.second / REPETITION / CASES);
						} // if !=0
					} // for times
					printf("\n");
					fflush(pfile);
					fflush(stdout);
				} // for num
			} // for size-ratio
		} // for intersection-ratio
	} // for minlength
	fclose(pfile);
}

void intersect_my_methods(bool loadfromfile) {
	using namespace msis;

	vector<uint32_t> out;

	uint32_t MaxBit = 31; // largest bit-length of element
	uint32_t minlength;
	size_t REPETITION;
	size_t CASES = 20;

	WallClockTimer timer;
	vector<float> intersectionsratios = { /*0.10, 0.20, 0.30, 0.40, */0.50 /*, 0.60,
	 0.70, 0.80, 0.90, 1.00 */};
	vector<uint32_t> sizeratios = { /*1, 4, 8, 16, 32, 64, 128, 256,*/512, 1024,
			2048, 4096, 8192 };

	FILE *pfile = fopen("output_mine.csv", "w+");
	fprintf(pfile, "sr,ir,num,name,time\n");
	std::map<std::string, size_t> times;
	for (uint32_t msb = 10; msb <= 10; ++msb) {
		minlength = 1U << msb;

		for (float ir : intersectionsratios) {
			printf("intersection ratio: \e[32m%3.0f%%\e[0m\n", ir * 100);
			for (uint32_t sr : sizeratios) {
				printf("  size ratio: \e[32m%4d\e[0m\n", sr);
				if (sr > 1000)
					REPETITION = 100;
				else
					REPETITION = 200;
				for (uint32_t num = 2; num < 11; num++) {
//					times["SvS_gallop_exact"] = 0;
//					times["SvS_gallop_rough"] = 0;
//					times["s_SvS_gallop_exact"] = 0;
//					times["s_SvS_gallop_rough"] = 0;
//					times["max_rough_opt"] = 0;
//					times["max_gallop_exact"] = 0;
					times["max_gallop_rough"] = 0;
					times["lemire"] = 0;
					times["svs_opt"] = 0;
					times["max"] = 0;
					time_t t = time(nullptr);
					tm* _tm = localtime(&t);
					printf("%02d:%02d:%02d> num: \e[32m%2d\e[0m  \n",
							_tm->tm_hour, _tm->tm_min, _tm->tm_sec, num);

					for (uint32_t k = 0; k < CASES; k++) {
						mySet multiset;
						if (loadfromfile)
							multiset = Set_collection::load_set(ir, sr, num, k);
						else {
							ClusteredDataGenerator cdg;
							multiset = genMultipleSets(cdg, minlength, num,
									1U << MaxBit, static_cast<float>(sr), ir);
						}

						// verification
//						auto it = multiset.begin();
//						vector<uint32_t> final_intersection = intersect(*it++,
//								*it++);
//						for (; it != multiset.end(); it++)
//							final_intersection = intersect(final_intersection,
//									*it);
//
//						msis::max_exact<msis::simdgallop_128_exact>(multiset,
//								out);
//						msis::max_rough<msis::simdgallop_128_4_rough>(multiset,
//								out);
//						std::cout << std::endl;
//
//						msis::s_SvS_exact<msis::simdgallop_128_exact>(multiset,
//								out);
//						msis::s_SvS_rough<msis::simdgallop_128_32_rough>(
//								multiset, out);
//
//						msis::SvS_exact<msis::simdgallop_128_exact>(multiset,
//								out);
//						msis::SvS_rough<msis::simdgallop_128_32_rough>(multiset,
//								out);
//
//						svs_opt(multiset, out);
//						Lemire_Gallop(multiset, out);
//						if (out != final_intersection) {
//							std::cerr << "bad result!  " << std::endl;
//							return;
//						} else
//							printf("good!  ");

//						timer.reset();
//						for (uint32_t howmany = 0; howmany < REPETITION;
//								++howmany) {
//							SvS_exact<msis::simdgallop_128_exact>(multiset,
//									out);
//						}
//						times["SvS_gallop_exact"] += timer.split();
//
//
//						timer.reset();
//						for (uint32_t howmany = 0; howmany < REPETITION;
//								++howmany) {
//							SvS_rough<msis::simdgallop_128_32_rough>(multiset,
//									out);
//						}
//						times["SvS_gallop_rough"] += timer.split();
//
//						timer.reset();
//						for (uint32_t howmany = 0; howmany < REPETITION;
//								++howmany) {
//							s_SvS_exact<msis::simdgallop_128_exact>(multiset,
//									out);
//						}
//						times["s_SvS_gallop_exact"] += timer.split();
//
//						timer.reset();
//						for (uint32_t howmany = 0; howmany < REPETITION;
//								++howmany) {
//							s_SvS_rough<msis::simdgallop_128_32_rough>(multiset,
//									out);
//						}
//						times["s_SvS_gallop_rough"] += timer.split();
//
//						timer.reset();
//						for (uint32_t howmany = 0; howmany < REPETITION;
//								++howmany) {
//							max_exact<msis::simdgallop_128_exact>(multiset,
//									out);
//						}
//						times["max_gallop_exact"] += timer.split();
//
//						timer.reset();
//						for (uint32_t howmany = 0; howmany < REPETITION;
//								++howmany) {
//							max_rough_opt(multiset, out);
//						}
//						times["max_rough_opt"] += timer.split();

						timer.reset();
						for (uint32_t howmany = 0; howmany < REPETITION;
								++howmany) {
							max_scalar(multiset, out);
						}
						times["max"] += timer.split();

						timer.reset();
						for (uint32_t howmany = 0; howmany < REPETITION;
								++howmany) {
							max_rough<msis::simdgallop_128_32_rough>(multiset,
									out);
						}
						times["max_gallop_rough"] += timer.split();

						timer.reset();
						for (uint32_t howmany = 0; howmany < REPETITION;
								++howmany) {
							msis::svs<SIMDgalloping>(multiset, out);
						}
						times["lemire"] += timer.split();

						timer.reset();
						for (uint32_t howmany = 0; howmany < REPETITION;
								++howmany) {
							svs_opt(multiset, out);
						}
						times["svs_opt"] += timer.split();

					} // for CASES
					for (auto time : times) {
						/*if (time.second != 0)*/{
							printf("%s: \e[31m%6.0f\e[0m  ", time.first.c_str(),
									(double) time.second / REPETITION / CASES);
							fprintf(pfile, "%d,%.0f,%d,%s,%.0f\n", sr, ir * 100,
									num, time.first.c_str(),
									(double) time.second / REPETITION / CASES);
						} // if !=0
					} // for times
					printf("\n");
					fflush(pfile);
					fflush(stdout);
				} // for num
			} // for size-ratio
		} // for intersection-ratio
	} // for minlength
	fclose(pfile);
}

int main(int argc, char **argv) {
	int c;
	bool loadfromfile = false;
	while ((c = getopt(argc, argv, "f")) != -1) {
		switch (c) {
		case 'f':
			loadfromfile = true;
			break;
		}
	}
	if (loadfromfile)
		std::cout << "load sets from file!" << std::endl;
	else
		std::cout << "generate sets in runtime!" << std::endl;

//	intersect_for_json(loadfromfile);
	intersect_my_methods(loadfromfile);
//	intersect_traditional_methods(loadfromfile);
}
