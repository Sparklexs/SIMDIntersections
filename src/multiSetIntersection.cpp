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

#include "multiSetIntersection.hpp"
#include "timer.h"
#include "synthetic.h"
#include <boost/preprocessor.hpp>

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
	msis::BY<msis::scalarBinarySearch>(sets, out);
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
#define METHOD (msis::SvS)(msis::s_SvS)(msis::sql)(msis::s_sql)(msis::max)
#define HEAD (_exact<msis::simd)(_rough<msis::simd)
#define SEARCH (linear_v)(gallop_v)
#define SIZE (4)(8)(16)(32)(64)(128)(256)(512)
#define END (_exact>)(_rough>)(_rough_plow>)
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
	    name.resize(name.size()-2);\
	    name.append(BOOST_PP_STRINGIZE(BOOST_PP_SEQ_ELEM(4,PRODUCT)));\
	    name.resize(name.size()-1);\
	    name.append("_").append(BOOST_PP_STRINGIZE(BOOST_PP_SEQ_ELEM(3,PRODUCT)));\
	    name=name.substr(6);\
	    times[name]=0;\

// TEST ALL THE METHOD
	BOOST_PP_SEQ_FOR_EACH_PRODUCT(LOOP_BODY,
			(/*(BOOST_PP_SEQ_ELEM(3,METHOD))*/METHOD)/*SvS,s_SvS,sql,s_sql,max*/
			((BOOST_PP_SEQ_ELEM(0,HEAD)))/*0:exact 1:rough*/
			(/*(BOOST_PP_SEQ_ELEM(1,SEARCH))*/SEARCH)/*0:linear 1:gallop*/
			(SIZE/*BOOST_PP_SEQ_REST_N(5,SIZE)*/) /*4,8,16,32,64,128,256,512*/
			((BOOST_PP_SEQ_ELEM(0,END))));/*0:exact 1:rough 2:rough_plow*/

	BOOST_PP_SEQ_FOR_EACH_PRODUCT(LOOP_BODY,
			(/*(BOOST_PP_SEQ_ELEM(3,METHOD))*/METHOD)/*SvS,s_SvS,sql,s_sql,max*/
			((BOOST_PP_SEQ_ELEM(1,HEAD)))/*0:exact 1:rough*/
			(/*(BOOST_PP_SEQ_ELEM(1,SEARCH))*/SEARCH)/*0:linear 1:gallop*/
			(SIZE/*BOOST_PP_SEQ_REST_N(5,SIZE)*/) /*4,8,16,32,64,128,256,512*/
			((BOOST_PP_SEQ_ELEM(1,END))));/*0:exact 1:rough 2:rough_plow*/

	BOOST_PP_SEQ_FOR_EACH_PRODUCT(LOOP_BODY,
			(/*(BOOST_PP_SEQ_ELEM(3,METHOD))*/METHOD)/*SvS,s_SvS,sql,s_sql,max*/
			((BOOST_PP_SEQ_ELEM(1,HEAD)))/*0:exact 1:rough*/
			(/*(BOOST_PP_SEQ_ELEM(1,SEARCH))*/SEARCH)/*0:linear 1:gallop*/
			(BOOST_PP_SEQ_REST_N(1,SIZE)) /*4,8,16,32,64,128,256,512*/
			((BOOST_PP_SEQ_ELEM(2,END))));/*0:exact 1:rough 2:rough_plow*/
#undef LOOP_BODY
	return times;
}

void intersect_for_all(bool loadfromfile) {
	using namespace msis;

	vector<uint32_t> out;

	uint32_t logMinLength = 12; // log of minimal array size
	uint32_t MaxBit = 31; // largest bit-length of element
	const uint32_t minlength = 1U << logMinLength;
	size_t REPETITION = 100;
	size_t CASES = 20;

	WallClockTimer timer;
#ifdef __INTEL_COMPILER
// Intel's support for C++ sucks
	vector<float> intersectionsratios;
	intersectionsratios.push_back(1.00);
	intersectionsratios.push_back(0.80);
	intersectionsratios.push_back(0.60);
	intersectionsratios.push_back(0.20);
	intersectionsratios.push_back(0.10);
	intersectionsratios.push_back(0.05);
	intersectionsratios.push_back(0.01);
	vector < uint32_t > sizeratios;
	sizeratios.push_back(1);
	sizeratios.push_back(2);
	sizeratios.push_back(3);
	sizeratios.push_back(5);
	sizeratios.push_back(10);
	sizeratios.push_back(20);
	sizeratios.push_back(40);
	sizeratios.push_back(80);
	sizeratios.push_back(200);
	sizeratios.push_back(500);
	sizeratios.push_back(1000);
#else
// proper C++
//      vector<float> intersectionsratios = { 0.01, 0.05, 0.10, 0.20, 0.60, 0.80,
//                      1.00 };
//      vector<uint32_t> sizeratios = { 1, 10, 20, 50, 100, 500, 1000, 5000, 10000 };
	vector<float> intersectionsratios = { 0.10, 0.50, 1.00 };
	vector<uint32_t> sizeratios = { 1, 10, 100, 1000 /*, 10000*/};
#endif

	FILE *pfile = fopen("output_disk.csv", "w+");
	fprintf(pfile, "sr,ir,num,name,time\n");
	std::map<std::string, size_t> times = init_time_array();

	for (float ir : intersectionsratios) {
		printf("intersection ratio: \e[32m%3.0f%%\e[0m\n", ir * 100);
		for (uint32_t sr : sizeratios) {
			printf("  size ratio: \e[32m%4d\e[0m\n", sr);
			for (uint32_t num = 2; num < 11; num++) {
				time_t t = time(nullptr);
				tm* _tm = localtime(&t);
				printf("%02d:%02d:%02d> num: \e[32m%2d\e[0m  ", _tm->tm_hour,
						_tm->tm_min, _tm->tm_sec, num);
				if (num == 0 and sr > 100)
					REPETITION = 3;
				else
					REPETITION = 100;
				REPETITION = 1;

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
					auto it = multiset.begin();
					vector<uint32_t> final_intersection = intersect(*it++,
							*it++);
					for (; it != multiset.end(); it++)
						final_intersection = intersect(final_intersection, *it);

					// start scalar intersection
//					for (uint32_t j = 0; j < NUMSCALARFUNC; j++) {
//						timer.reset();
//						for (uint32_t howmany = 0; howmany < REPETITION;
//								++howmany) {
//							scalarFUNC[j](multiset, out);
//						}
//
////						if (out != final_intersection) {
////							std::cerr << "bad result!  " << std::endl;
////							return;
////						} else
////							printf("good!  ");
//
//						times[NAMESCALARFUNC[j]] += timer.split();
//					}
//					Schlegel(multiset, timer, times);

					std::string name;
#define LOOP_BODY(R,PRODUCT)\
	        timer.reset();\
	        for (uint32_t howmany = 0; howmany < REPETITION;\
	        	 ++howmany) {\
	               BOOST_PP_CAT(\
	               BOOST_PP_CAT(\
	               BOOST_PP_CAT(\
	               BOOST_PP_CAT(\
	               BOOST_PP_SEQ_ELEM(0,PRODUCT),\
	               BOOST_PP_SEQ_ELEM(1,PRODUCT)),\
	               BOOST_PP_SEQ_ELEM(2,PRODUCT)),\
	               BOOST_PP_SEQ_ELEM(3,PRODUCT)),\
	               BOOST_PP_SEQ_ELEM(4,PRODUCT))(multiset,out);\
	        }\
			if (out != final_intersection) {\
					std::cerr << "bad result!  " << std::endl;\
					return;\
			} else\
					printf("good!  ");\
	        name=(BOOST_PP_STRINGIZE(BOOST_PP_SEQ_ELEM(0,PRODUCT)));\
	        name.append("_").append(BOOST_PP_STRINGIZE(BOOST_PP_SEQ_ELEM(2,PRODUCT)));\
	        name.resize(name.size()-2);\
	        name.append(BOOST_PP_STRINGIZE(BOOST_PP_SEQ_ELEM(4,PRODUCT)));\
	        name.resize(name.size()-1);\
	        name.append("_").append(BOOST_PP_STRINGIZE(BOOST_PP_SEQ_ELEM(3,PRODUCT)));\
	        name=name.substr(6);\
	        times[name] += timer.split();\

//					BOOST_PP_SEQ_FOR_EACH_PRODUCT(LOOP_BODY,
//							((BOOST_PP_SEQ_ELEM(4,METHOD)))/*0:SvS,1:s_SvS,2:sql,3:s_sql,4:max*/
//							((BOOST_PP_SEQ_ELEM(1,HEAD)))/*0:exact 1:rough*/
//							(/*(BOOST_PP_SEQ_ELEM(1,SEARCH))*/SEARCH)/*0:linear 1:gallop*/
//							(BOOST_PP_SEQ_REST_N(1,SIZE)) /*4,8,16,32,64,128,256,512*/
//							(BOOST_PP_SEQ_POP_FRONT(END)));/*0:exact 1:rough 2:rough_plow*/
//
					// TEST ALL THE METHOD
//					BOOST_PP_SEQ_FOR_EACH_PRODUCT(LOOP_BODY,
//							(/*(BOOST_PP_SEQ_ELEM(3,METHOD))*/METHOD)/*SvS,s_SvS,sql,s_sql,max*/
//							((BOOST_PP_SEQ_ELEM(0,HEAD)))/*0:exact 1:rough*/
//							(/*(BOOST_PP_SEQ_ELEM(1,SEARCH))*/SEARCH)/*0:linear 1:gallop*/
//							(SIZE/*BOOST_PP_SEQ_REST_N(5,SIZE)*/) /*4,8,16,32,64,128,256,512*/
//							((BOOST_PP_SEQ_ELEM(0,END))));/*0:exact 1:rough 2:rough_plow*/
//
//					BOOST_PP_SEQ_FOR_EACH_PRODUCT(LOOP_BODY,
//							(/*(BOOST_PP_SEQ_ELEM(3,METHOD))*/METHOD)/*SvS,s_SvS,sql,s_sql,max*/
//							((BOOST_PP_SEQ_ELEM(1,HEAD)))/*0:exact 1:rough*/
//							(/*(BOOST_PP_SEQ_ELEM(1,SEARCH))*/SEARCH)/*0:linear 1:gallop*/
//							(SIZE/*BOOST_PP_SEQ_REST_N(5,SIZE)*/) /*4,8,16,32,64,128,256,512*/
//							((BOOST_PP_SEQ_ELEM(1,END))));/*0:exact 1:rough 2:rough_plow*/

					BOOST_PP_SEQ_FOR_EACH_PRODUCT(LOOP_BODY,
							(/*(BOOST_PP_SEQ_ELEM(3,METHOD))*/METHOD)/*SvS,s_SvS,sql,s_sql,max*/
							((BOOST_PP_SEQ_ELEM(1,HEAD)))/*0:exact 1:rough*/
							(/*(BOOST_PP_SEQ_ELEM(1,SEARCH))*/SEARCH)/*0:linear 1:gallop*/
							(BOOST_PP_SEQ_REST_N(1,SIZE)) /*4,8,16,32,64,128,256,512*/
							((BOOST_PP_SEQ_ELEM(2,END))));/*0:exact 1:rough 2:rough_plow*/
#undef LOOP_BODY
				}
				for (auto time : times) {
					if (time.second != 0) {
//						printf("%s: \e[31m%6.0f\e[0m  ", time.first.c_str(),
//								(double) time.second / REPETITION / CASES);
						fprintf(pfile, "%d,%.0f,%d,%s,%.0f\n", sr, ir * 100,
								num, time.first.c_str(),
								(double) time.second / REPETITION / CASES);
					}
				}
				printf("\n");
				fflush(pfile);
				fflush(stdout);
			} // for num
		} // for size-ratio
	} // for intersection-ratio
	fclose(pfile);
}

void intersect_for_10k(bool loadfromfile) {
	using namespace msis;

	vector<uint32_t> out;

	uint32_t logMinLength = 12; // log of minimal array size
	uint32_t MaxBit = 31; // largest bit-length of element
	const uint32_t minlength = 1U << logMinLength;
	size_t REPETITION = 100;
	size_t CASES = 20;

	WallClockTimer timer;
#ifdef __INTEL_COMPILER
// Intel's support for C++ sucks
	vector<float> intersectionsratios;
	intersectionsratios.push_back(1.00);
	intersectionsratios.push_back(0.80);
	intersectionsratios.push_back(0.60);
	intersectionsratios.push_back(0.20);
	intersectionsratios.push_back(0.10);
	intersectionsratios.push_back(0.05);
	intersectionsratios.push_back(0.01);
	vector < uint32_t > sizeratios;
	sizeratios.push_back(1);
	sizeratios.push_back(2);
	sizeratios.push_back(3);
	sizeratios.push_back(5);
	sizeratios.push_back(10);
	sizeratios.push_back(20);
	sizeratios.push_back(40);
	sizeratios.push_back(80);
	sizeratios.push_back(200);
	sizeratios.push_back(500);
	sizeratios.push_back(1000);
#else
// proper C++
//      vector<float> intersectionsratios = { 0.01, 0.05, 0.10, 0.20, 0.60, 0.80,
//                      1.00 };
//      vector<uint32_t> sizeratios = { 1, 10, 20, 50, 100, 500, 1000, 5000, 10000 };
	vector<float> intersectionsratios = { /*0.10,0.50,*/ 1.00 };
	vector<uint32_t> sizeratios = { 10000 };
#endif

	FILE *pfile = fopen("output_disk_10k.csv", "w+");
	fprintf(pfile, "sr,ir,num,name,time\n");
	std::map<std::string, size_t> times = init_time_array();

	for (float ir : intersectionsratios) {
		printf("intersection ratio: \e[32m%3.0f%%\e[0m\n", ir * 100);
		for (uint32_t sr : sizeratios) {
			printf("  size ratio: \e[32m%4d\e[0m\n", sr);
			for (uint32_t num = 2; num < 11; num++) {
				time_t t = time(nullptr);
				tm* _tm = localtime(&t);
				printf("%02d:%02d:%02d> num: \e[32m%2d\e[0m  ", _tm->tm_hour,
						_tm->tm_min, _tm->tm_sec, num);
				if (num == 0 and sr > 100)
					REPETITION = 3;
				else
					REPETITION = 100;

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
//					auto it = multiset.begin();
//					vector<uint32_t> final_intersection = intersect(*it++,
//							*it++);
//					for (; it != multiset.end(); it++)
//						final_intersection = intersect(final_intersection, *it);

					// start scalar intersection
					for (uint32_t j = 0; j < NUMSCALARFUNC; j++) {
						timer.reset();
						for (uint32_t howmany = 0; howmany < REPETITION;
								++howmany) {
							scalarFUNC[j](multiset, out);
						}

//						if (out != final_intersection) {
//							std::cerr << "bad result!  " << std::endl;
//							return;
//						} else
//							printf("good!  ");

						times[NAMESCALARFUNC[j]] += timer.split();
					}
					Schlegel(multiset, timer, times);

					std::string name;
#define LOOP_BODY(R,PRODUCT)\
	        timer.reset();\
	        for (uint32_t howmany = 0; howmany < REPETITION;\
	        	 ++howmany) {\
	               BOOST_PP_CAT(\
	               BOOST_PP_CAT(\
	               BOOST_PP_CAT(\
	               BOOST_PP_CAT(\
	               BOOST_PP_SEQ_ELEM(0,PRODUCT),\
	               BOOST_PP_SEQ_ELEM(1,PRODUCT)),\
	               BOOST_PP_SEQ_ELEM(2,PRODUCT)),\
	               BOOST_PP_SEQ_ELEM(3,PRODUCT)),\
	               BOOST_PP_SEQ_ELEM(4,PRODUCT))(multiset,out);\
	        }\
	        name=(BOOST_PP_STRINGIZE(BOOST_PP_SEQ_ELEM(0,PRODUCT)));\
	        name.append("_").append(BOOST_PP_STRINGIZE(BOOST_PP_SEQ_ELEM(2,PRODUCT)));\
	        name.resize(name.size()-2);\
	        name.append(BOOST_PP_STRINGIZE(BOOST_PP_SEQ_ELEM(4,PRODUCT)));\
	        name.resize(name.size()-1);\
	        name.append("_").append(BOOST_PP_STRINGIZE(BOOST_PP_SEQ_ELEM(3,PRODUCT)));\
	        name=name.substr(6);\
	        times[name] += timer.split();\

//					BOOST_PP_SEQ_FOR_EACH_PRODUCT(LOOP_BODY,
//							((BOOST_PP_SEQ_ELEM(4,METHOD)))/*0:SvS,1:s_SvS,2:sql,3:s_sql,4:max*/
//							((BOOST_PP_SEQ_ELEM(1,HEAD)))/*0:exact 1:rough*/
//							(/*(BOOST_PP_SEQ_ELEM(1,SEARCH))*/SEARCH)/*0:linear 1:gallop*/
//							(BOOST_PP_SEQ_REST_N(1,SIZE)) /*4,8,16,32,64,128,256,512*/
//							(BOOST_PP_SEQ_POP_FRONT(END)));/*0:exact 1:rough 2:rough_plow*/
//
					// TEST ALL THE METHOD
					BOOST_PP_SEQ_FOR_EACH_PRODUCT(LOOP_BODY,
							(/*(BOOST_PP_SEQ_ELEM(3,METHOD))*/METHOD)/*SvS,s_SvS,sql,s_sql,max*/
							((BOOST_PP_SEQ_ELEM(0,HEAD)))/*0:exact 1:rough*/
							(/*(BOOST_PP_SEQ_ELEM(1,SEARCH))*/SEARCH)/*0:linear 1:gallop*/
							(SIZE/*BOOST_PP_SEQ_REST_N(5,SIZE)*/) /*4,8,16,32,64,128,256,512*/
							((BOOST_PP_SEQ_ELEM(0,END))));/*0:exact 1:rough 2:rough_plow*/

					BOOST_PP_SEQ_FOR_EACH_PRODUCT(LOOP_BODY,
							(/*(BOOST_PP_SEQ_ELEM(3,METHOD))*/METHOD)/*SvS,s_SvS,sql,s_sql,max*/
							((BOOST_PP_SEQ_ELEM(1,HEAD)))/*0:exact 1:rough*/
							(/*(BOOST_PP_SEQ_ELEM(1,SEARCH))*/SEARCH)/*0:linear 1:gallop*/
							(SIZE/*BOOST_PP_SEQ_REST_N(5,SIZE)*/) /*4,8,16,32,64,128,256,512*/
							((BOOST_PP_SEQ_ELEM(1,END))));/*0:exact 1:rough 2:rough_plow*/

//					BOOST_PP_SEQ_FOR_EACH_PRODUCT(LOOP_BODY,
//							(/*(BOOST_PP_SEQ_ELEM(3,METHOD))*/METHOD)/*SvS,s_SvS,sql,s_sql,max*/
//							((BOOST_PP_SEQ_ELEM(1,HEAD)))/*0:exact 1:rough*/
//							(/*(BOOST_PP_SEQ_ELEM(1,SEARCH))*/SEARCH)/*0:linear 1:gallop*/
//							(BOOST_PP_SEQ_REST_N(1,SIZE)) /*4,8,16,32,64,128,256,512*/
//							((BOOST_PP_SEQ_ELEM(2,END))));/*0:exact 1:rough 2:rough_plow*/
#undef LOOP_BODY
				}
				for (auto time : times) {
					if (time.second != 0) {
//						printf("%s: \e[31m%6.0f\e[0m  ", time.first.c_str(),
//								(double) time.second / REPETITION / CASES);
						fprintf(pfile, "%d,%.0f,%d,%s,%.0f\n", sr, ir * 100,
								num, time.first.c_str(),
								(double) time.second / REPETITION / CASES);
					}
				}
				printf("\n");
				fflush(pfile);
				fflush(stdout);
			} // for num
		} // for size-ratio
	} // for intersection-ratio
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
	intersect_for_all(loadfromfile);
	intersect_for_10k(loadfromfile);
}
