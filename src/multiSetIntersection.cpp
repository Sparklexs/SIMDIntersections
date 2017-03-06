/*
 * multiSetIntersection.cpp
 * scalar algorithms for intersecting multiple sets.
 * Reference:
 * (1)	Culpepper, J. S., & Moffat, A. (2010).
 *  	Efficient set intersection for inverted indexing.
 *   	Acm Transactions on Information Systems, 29(1), 1.
 *
 * (2)	Barbay, J., L¨®pez-Ortiz, A., Lu, T., & Salinger, A.
 *		(2009). An experimental investigation of set
 *		intersection algorithms for text searching.
 *		Journal of Experimental Algorithmics, 14.
 *
 * (3)	Baezayates, R., & Salinger, A. (2010).
 *		Fast Intersection Algorithms for Sorted Sequences.
 *		Algorithms and Applications. Springer Berlin Heidelberg.
 *
 * (4)	Barbay, J., L¨®pezortiz, A., & Lu, T. (2006).
 *		Faster Adaptive Set Intersections for Text Searching.
 *		International Conference on Experimental Algorithms
 *		(Vol.4007, pp.146-157). Springer-Verlag.
 *
 *  Created on: 2016/12/20
 *      Author: SparkleXS
 */

#include "multiSetIntersection.hpp"
#include "timer.h"
#include "synthetic.h"
#include <boost/preprocessor.hpp>

void small_vs_small_scalar(const mySet &sets, std::vector<uint32_t> &out) {
	msis::svs<BSintersection>(sets, out);
}

void set_vs_set_scalar(const mySet &sets, std::vector<uint32_t> &out) {
	msis::SvS_exact<msis::scalarBinarySearch>(sets, out);
}

void swapping_set_vs_set_scalar(const mySet& sets, vector<uint32_t>& out) {
	msis::s_SvS_exact<msis::scalarBinarySearch>(sets, out);
}

void sequential_scalar(const mySet &sets, std::vector<uint32_t> &out) {
	msis::sql_exact<msis::scalargallop>(sets, out);
}

void small_sequential_scalar(const mySet &sets, std::vector<uint32_t> &out) {
	msis::s_sql_exact<msis::scalargallop>(sets, out);
}

void max__scalar(const mySet &sets, std::vector<uint32_t> &out) {
	msis::max_exact<msis::scalargallop>(sets, out);
}

void BaezaYates_scalar(const mySet &sets, std::vector<uint32_t> &out) {
	msis::BY<msis::scalarBinarySearch>(sets, out);
}

typedef void (*intersectionFUNC)(const mySet &sets, std::vector<uint32_t> &out);
/* SIMD methods are template functions that can't be applied to typedef */
intersectionFUNC scalarFUNC[] = { small_vs_small_scalar, set_vs_set_scalar,
		swapping_set_vs_set_scalar, msis::adp_scalar, msis::s_adp_scalar,
		sequential_scalar, small_sequential_scalar, max__scalar,
		BaezaYates_scalar };

constexpr int NUMSCALARFUNC = sizeof(scalarFUNC) / sizeof(scalarFUNC[0]);

const std::string NAMESCALARFUNC[] = { "small_vs_small", "set_vs_set",
		"swapping_set_vs_set", "adaptive", "small_adaptive", "sequential",
		"small_sequential", "max", "BaezaYates" };

/////////////////////////////////////////
#define METHOD (msis::SvS)(msis::s_SvS)(msis::sql)(msis::s_sql)(msis::max)
#define HEAD (_exact<msis::simd)(_rough<msis::simd)
#define SEARCH (linear_v)(gallop_v)
#define SIZE (4)(8)(16)(32)(64)(128)(256)(512)
#define END (_exact>)(_rough>)(_rough_plow>)
//////////////////////////////////////////

int main() {
	using namespace msis;

	uint32_t logMinLength = 12; // log of minimal array size
	uint32_t MaxBit = 31; // largest bit-length of element
	const uint32_t minlength = 1U << logMinLength;
	ClusteredDataGenerator cdg;

	mySet MultiSets;
	vector<uint32_t> out;

	WallClockTimer timer;
	size_t time = 0;
	const int REPETITION = 50;

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
	vector<float> intersectionsratios = { 0.01, 0.05, 0.10, 0.20, 0.60, 0.80,
			1.00 };

	vector<uint32_t> sizeratios = { 1, 10, 20, 50, 100, 500, 1000, 5000, 10000 };
#endif

	cout << "########### Intersection benchmark ###########" << endl;

	cout << "# columns are size ratios" << endl;
	cout << "# rows are intersection ratios" << endl;
// sxs: namely log(U/n)
	cout << "# average gaps in bits for smallest array: "

	<< std::setprecision(3) << log(1 + (1U << MaxBit) * 1.0 / minlength)
			<< endl;
	cout << "# average gaps in bits for last largest array: "

	<< std::setprecision(3)
			<< log(1 + (1U << MaxBit) * 1.0 / (sizeratios.back() * minlength))
			<< endl;
	cout << "#############################################" << endl << endl;

	for (float ir : intersectionsratios) {
		printf("intersection ratio: \e[32m%3.0f%%\e[0m\n", ir * 100);
		for (uint32_t sr : sizeratios) {
			printf("  size ratio: \e[32m%4d\e[0m\n", sr);
			for (int i = 3; i < 11; i++) {
				printf("    num: \e[32m%2d\e[0m  ", i);

				// generate sets
				MultiSets = genMultipleSets(cdg, minlength, i, 1U << MaxBit,
						static_cast<float>(sr), ir);
//				MultiSets = genMultipleSets(cdg, 8, i, 16,
//						static_cast<float>(sr), ir);

				// verification
				auto it = MultiSets.begin();
				vector<uint32_t> final_intersection = intersect(*it++, *it++);
				for (; it != MultiSets.end(); it++)
					final_intersection = intersect(final_intersection, *it);

				timer.reset();
				for (int howmany = 0; howmany < REPETITION; ++howmany) {
					msis::SvS_exact<msis::scalarBinarySearch>(MultiSets, out);
				}
				time = timer.split();
				printf("%s: \e[31m%6.0f\e[0m  ", "SvS_B",
						(double) time / REPETITION);

				timer.reset();
				for (int howmany = 0; howmany < REPETITION; ++howmany) {
					msis::SvS_exact<msis::scalargallop>(MultiSets, out);
				}
				time = timer.split();
				printf("%s: \e[31m%6.0f\e[0m  ", "SvS_G",
						(double) time / REPETITION);

				timer.reset();
				for (int howmany = 0; howmany < REPETITION; ++howmany) {
					msis::s_SvS_exact<msis::scalarBinarySearch>(MultiSets, out);
				}
				time = timer.split();
				printf("%s: \e[31m%6.0f\e[0m  ", "s_SvS_B",
						(double) time / REPETITION);

				timer.reset();
				for (int howmany = 0; howmany < REPETITION; ++howmany) {
					msis::s_SvS_exact<msis::scalargallop>(MultiSets, out);
				}
				time = timer.split();
				printf("%s: \e[31m%6.0f\e[0m  ", "s_SvS_G",
						(double) time / REPETITION);

//				std::string name;
//#define LOOP_BODY(R,PRODUCT)\
//		timer.reset();											\
//		for (int howmany = 0; howmany < REPETITION; ++howmany) {\
//			BOOST_PP_CAT(\
//					BOOST_PP_CAT(\
//					BOOST_PP_CAT(\
//					BOOST_PP_CAT(\
//					BOOST_PP_SEQ_ELEM(0,PRODUCT),\
//					BOOST_PP_SEQ_ELEM(1,PRODUCT)),\
//					BOOST_PP_SEQ_ELEM(2,PRODUCT)),\
//					BOOST_PP_SEQ_ELEM(3,PRODUCT)),\
//					BOOST_PP_SEQ_ELEM(4,PRODUCT))(MultiSets,out);\
//		}\
//		time = timer.split();\
//		name=(BOOST_PP_STRINGIZE(BOOST_PP_SEQ_ELEM(0,PRODUCT)));\
//		name.append(BOOST_PP_STRINGIZE(BOOST_PP_SEQ_ELEM(4,PRODUCT)));\
//		name.append("(").append(BOOST_PP_STRINGIZE(BOOST_PP_SEQ_ELEM(3,PRODUCT))).append(")");\
//		name=name.substr(6);\
//		printf("%s: \e[31m%6.0f\e[0m  ", name.c_str(),\
//					(double) time / REPETITION);\
//		if (out != final_intersection) {\
//			std::cerr << "bad result!  " << std::endl;\
//			return 1;\
//		} else\
//			printf("good!  ");\
//
//	BOOST_PP_SEQ_FOR_EACH_PRODUCT(LOOP_BODY,
//	((BOOST_PP_SEQ_ELEM(3,METHOD))/*METHOD*/)/*s_vs_s,s_s_vs_s,seq,max*/
//	((BOOST_PP_SEQ_ELEM(1,HEAD)))/*0:exact 1:rough*/
//	((BOOST_PP_SEQ_ELEM(1,SEARCH)))/*0:linear 1:gallop*/
//	(/*SIZE*/BOOST_PP_SEQ_REST_N(5,SIZE)) /*4,8,16,32,64,128,256,512*/
//	(/*(BOOST_PP_SEQ_ELEM(2,END))*/BOOST_PP_SEQ_POP_FRONT(END)));/*0:exact 1:rough 2:rough_plow*/
//#undef LOOP_BODY

				// start scalar intersection
//				for (int i = 0; i < NUMSCALARFUNC; i++) {
//					timer.reset();
//					for (int howmany = 0; howmany < REPETITION; ++howmany) {
//						scalarFUNC[i](MultiSets, out);
//					}
//					time = timer.split();
//					printf("%s: \e[31m%6.0f\e[0m  ", NAMESCALARFUNC[i].c_str(),
//							(double) time / REPETITION);
//
////					if (out != final_intersection) {
////						std::cerr << "bad result!  " << std::endl;
////						return 1;
////					} else
////						printf("good!  ");
//				} // for intersection
				printf("\n");
				fflush(stdout);
			} // for num
		} // for size-ratio
	} // for intersection-ratio
}
