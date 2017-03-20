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

void Zipper(const mySet &sets, std::vector<uint32_t> &out) {
	msis::svs<branchlessintersection>(sets, out);
}

void Schlegel(const std::vector<mySet> &sets, WallClockTimer timer, FILE* pfile,
		uint32_t sr, float ir, int i) {
	uint32_t num_set = sets[0].size();
	vector<vector<uint16_t> > pdata(sets.size() * num_set);
	for (auto multiset : sets) {
		mySet::iterator it = multiset.begin();
		for (uint32_t i = 0; i < multiset.size(); i++) {
			pdata[i].resize(it->size() * 4);
			const size_t c = partitioned::partition(it->data(), it->size(),
					&pdata[i][0], pdata[i].size());
			pdata[i].resize(c);
			vector<uint16_t>(pdata[i]).swap(pdata[i]);
			it++;
		}
	}
	for (auto algo : partschemes) {
		timer.reset();
		for (uint32_t i = 0; i < sets.size(); i++) {
			for (uint32_t j = 0; j < 100; j++) {
				vector<uint16_t> intersection(pdata[0 + i * num_set]);
				for (uint32_t k = 0; k < num_set; k++) {
					size_t answer = algo.second(intersection.data(),
							pdata[k + i * num_set].data(), intersection.size(),
							pdata[k + i * num_set].size());
					intersection.resize(answer);
				}
			}
		}
		size_t time = timer.split();
		fprintf(pfile, "%d,%.0f,%d,%s,%.0f\n", sr, ir, i, algo.first.c_str(),
				(double) time / 100 / sets.size());
//		printf("%s: \e[31m%6.0f\e[0m  ", algo.first.c_str(),
//				(double) time / 100 / sets.size());
	}
}

void Schlegel(const mySet &multiset, WallClockTimer timer, FILE* pfile,
		uint32_t sr, float ir, int i) {
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

		size_t time = timer.split();
		fprintf(pfile, "%d,%.0f,%d,%s,%.0f\n", sr, ir, i, algo.first.c_str(),
				(double) time / 100);
		printf("%s: \e[31m%6.0f\e[0m  ", algo.first.c_str(),
				(double) time / 100);
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

class Set_collection {
	const static uint32_t HEADERSIZE = sizeof(float) + 3 * sizeof(uint32_t);
	const static uint32_t REPETITION = 20;
	const static uint32_t NUM = 10;
public:
	static long pos;
	static float selectivity;
	static uint32_t ratio;
	static uint32_t rank;
	static uint32_t num;
	static void sequentially_load_file(mySet & multiset) {
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
//						outfile.flush();
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
//						outfile.flush();
					}
				}
			}
			outfile.flush();
			outfile.close();
			sequentially_load_file(multiset);
		} else {
			size_t totallength = 0;
			infile.seekg(pos);

			if (rank == REPETITION) {
				rank = 0;
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
			}
			if (ratio == 10000) {
				std::ostringstream ostr;
				ostr << selectivity << "_" << num;
				std::ofstream file(ostr.str(), std::ios::binary);

				uint32_t temp = 4096;
				file.write((char*) &temp, 4);

				size_t step = 1024 * 1024;
				size_t n = (totallength - 4) / step;
				char buffer[step];
				for (size_t i = 0; i < n; i++) {
					infile.read(buffer, step);
					file.write(buffer, step);
				}
				size_t last = (totallength - 4) % step;
				infile.read(buffer, last);
				file.write(buffer, last);
				file.flush();
				file.close();
			} else {
				infile.seekg(totallength - 4, std::ios_base::cur);
			}
			rank = 20;

//			uint32_t length;
//			for (uint32_t i = 0; i < num; i++) {
//				if (rank == 0 and i == 0)
//					length = 4096;
//				else
//					infile.read((char*) &length, sizeof(uint32_t));
//				vector<uint32_t> set(length);
//				uint32_t* ptr = set.data();
//
//				infile.read((char*) ptr, length * sizeof(uint32_t));
//				multiset.insert(set);
//			}
//			rank++;
			pos = infile.tellg();

		}
	}

};
long Set_collection::pos = 0;
float Set_collection::selectivity = 0;
uint32_t Set_collection::ratio = 10000;
uint32_t Set_collection::rank = REPETITION;
uint32_t Set_collection::num = NUM;

typedef void (*intersectionFUNC)(const mySet &sets, std::vector<uint32_t> &out);
/* SIMD methods are template functions that can't be applied to typedef */
intersectionFUNC scalarFUNC[] = { Zipper, Lemire_V1, Lemire_V3, Lemire_Gallop,
		Inoue, set_vs_set_scalar, swapping_set_vs_set_scalar, /*msis::adp_scalar,
		 msis::s_adp_scalar,*/sequential_scalar, small_sequential_scalar,
		max_scalar, BaezaYates_scalar };

constexpr int NUMSCALARFUNC = sizeof(scalarFUNC) / sizeof(scalarFUNC[0]);

const std::string NAMESCALARFUNC[] = { "zip", "v1", "v3", "gallop", "inoue",
		"SvS", "s_SvS", /*"adp", "s_adp",*/"sql", "s_sql", "max", "bys" };

/////////////////////////////////////////
#define METHOD (msis::SvS)(msis::s_SvS)(msis::sql)(msis::s_sql)(msis::max)
#define HEAD (_exact<msis::simd)(_rough<msis::simd)
#define SEARCH (linear_v)(gallop_v)
#define SIZE (4)(8)(16)(32)(64)(128)(256)(512)
#define END (_exact>)(_rough>)(_rough_plow>)
//////////////////////////////////////////

void intersect_realtime_data() {
	using namespace msis;
	uint32_t logMinLength = 12; // log of minimal array size
	uint32_t MaxBit = 31; // largest bit-length of element
	const uint32_t minlength = 1U << logMinLength;
	size_t REPETITION = 100;
	size_t CASES = 20;

	vector<uint32_t> out;

	WallClockTimer timer;
	size_t time = 0;
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
//	vector<float> intersectionsratios = { 0.01, 0.05, 0.10, 0.20, 0.60, 0.80,
//			1.00 };
//	vector<uint32_t> sizeratios = { 1, 10, 20, 50, 100, 500, 1000, 5000, 10000 };
	vector<float> intersectionsratios = { 0.10, 0.50, 1.00 };
	vector<uint32_t> sizeratios = { 1, 10, 100, 1000, 10000 };
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

	FILE *pfile = fopen("output_realtime.csv", "w+");
	fprintf(pfile, "sr,ir,num,name,time\n");

	for (float ir : intersectionsratios) {
		printf("intersection ratio: \e[32m%3.0f%%\e[0m\n", ir * 100);
		for (uint32_t sr : sizeratios) {
			printf("  size ratio: \e[32m%4d\e[0m\n", sr);
			for (int i = 2; i < 11; i++) {
				printf("    num: \e[32m%2d\e[0m  ", i);

				vector<mySet> multisets;
				ClusteredDataGenerator cdg;
				for (uint32_t c = 0; c < CASES; c++)
					multisets.push_back(
							genMultipleSets(cdg, minlength, i, 1U << MaxBit,
									static_cast<float>(sr), ir));

				// start scalar intersection
//				for (int j = 0; j < NUMSCALARFUNC; j++) {
//					if (i == 0 and sr > 100)
//						REPETITION = 3;
//					else
//						REPETITION = 100;
//
//					timer.reset();
//					for (uint32_t cases = 0; cases < CASES; cases++) {
//						// verification
////						auto it = sets[cases].begin();
////						vector<uint32_t> final_intersection = intersect(*it++,
////								*it++);
////						for (; it != sets[cases].end(); it++)
////							final_intersection = intersect(final_intersection,
////									*it);
////
//						for (uint32_t howmany = 0; howmany < REPETITION;
//								++howmany) {
//							scalarFUNC[j](multisets[cases], out);
//						}
////
////						if (out != final_intersection) {
////							std::cerr << "bad result!  " << std::endl;
////							return 1;
////						} else
////							printf("good!  ");
//					}
//					time = timer.split();
////					printf("%s: \e[31m%6.0f\e[0m  ", NAMESCALARFUNC[j].c_str(),
////							(double) time / REPETITION / CASES);
//					fprintf(pfile, "%d,%.0f,%d,%s,%.0f\n", sr, ir * 100, i,
//							NAMESCALARFUNC[j].c_str(),
//							(double) time / REPETITION / CASES);
//				} // for intersection
//				Schlegel(multisets, timer, pfile, sr, ir * 100, i);
//				printf("\n");

				std::string name;
#define LOOP_BODY(R,PRODUCT)\
		timer.reset();											\
		for (uint32_t cases = 0; cases < CASES; cases++) {\
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
				BOOST_PP_SEQ_ELEM(4,PRODUCT))(multisets[cases],out);\
			}\
		}\
		time = timer.split();\
		name=(BOOST_PP_STRINGIZE(BOOST_PP_SEQ_ELEM(0,PRODUCT)));\
		name.append("_").append(BOOST_PP_STRINGIZE(BOOST_PP_SEQ_ELEM(2,PRODUCT)));\
		name.resize(name.size()-2);\
		name.append(BOOST_PP_STRINGIZE(BOOST_PP_SEQ_ELEM(4,PRODUCT)));\
		name.resize(name.size()-1);\
		name.append("_").append(BOOST_PP_STRINGIZE(BOOST_PP_SEQ_ELEM(3,PRODUCT)));\
		name=name.substr(6);\
		fprintf(pfile, "%d,%.0f,%d,%s,%.0f\n", sr, ir * 100, i,\
				name.c_str(),(double) time / REPETITION/CASES);\
		printf("%s: \e[31m%6.0f\e[0m  \n", name.c_str(),\
				(double) time / REPETITION/CASES);

				BOOST_PP_SEQ_FOR_EACH_PRODUCT(LOOP_BODY,
						((BOOST_PP_SEQ_ELEM(4,METHOD)))/*0:SvS,1:s_SvS,2:sql,3:s_sql,4:max*/
						((BOOST_PP_SEQ_ELEM(1,HEAD)))/*0:exact 1:rough*/
						(/*(BOOST_PP_SEQ_ELEM(1,SEARCH))*/SEARCH)/*0:linear 1:gallop*/
						(BOOST_PP_SEQ_REST_N(1,SIZE)/*BOOST_PP_SEQ_REST_N(5,SIZE)*/) /*4,8,16,32,64,128,256,512*/
						(BOOST_PP_SEQ_POP_FRONT(END)));/*0:exact 1:rough 2:rough_plow*/

				// TEST ALL THE METHOD
//				BOOST_PP_SEQ_FOR_EACH_PRODUCT(LOOP_BODY,
//						(/*(BOOST_PP_SEQ_ELEM(3,METHOD))*/METHOD)/*SvS,s_SvS,sql,s_sql,max*/
//						((BOOST_PP_SEQ_ELEM(0,HEAD)))/*0:exact 1:rough*/
//						(/*(BOOST_PP_SEQ_ELEM(1,SEARCH))*/SEARCH)/*0:linear 1:gallop*/
//						(SIZE/*BOOST_PP_SEQ_REST_N(5,SIZE)*/) /*4,8,16,32,64,128,256,512*/
//						((BOOST_PP_SEQ_ELEM(0,END))/*BOOST_PP_SEQ_POP_FRONT(END)*/));/*0:exact 1:rough 2:rough_plow*/
//
//				BOOST_PP_SEQ_FOR_EACH_PRODUCT(LOOP_BODY,
//						(/*(BOOST_PP_SEQ_ELEM(3,METHOD))*/METHOD)/*SvS,s_SvS,sql,s_sql,max*/
//						((BOOST_PP_SEQ_ELEM(1,HEAD)))/*0:exact 1:rough*/
//						(/*(BOOST_PP_SEQ_ELEM(1,SEARCH))*/SEARCH)/*0:linear 1:gallop*/
//						(SIZE/*BOOST_PP_SEQ_REST_N(5,SIZE)*/) /*4,8,16,32,64,128,256,512*/
//						((BOOST_PP_SEQ_ELEM(1,END))));/*0:exact 1:rough 2:rough_plow*/
//
//				BOOST_PP_SEQ_FOR_EACH_PRODUCT(LOOP_BODY,
//						(/*(BOOST_PP_SEQ_ELEM(3,METHOD))*/METHOD)/*SvS,s_SvS,sql,s_sql,max*/
//						((BOOST_PP_SEQ_ELEM(1,HEAD)))s/*0:exact 1:rough*/
//						(/*(BOOST_PP_SEQ_ELEM(1,SEARCH))*/SEARCH)/*0:linear 1:gallop*/
//						(BOOST_PP_SEQ_REST_N(1,SIZE)/*BOOST_PP_SEQ_REST_N(5,SIZE)*/) /*4,8,16,32,64,128,256,512*/
//						((BOOST_PP_SEQ_ELEM(2,END))));/*0:exact 1:rough 2:rough_plow*/
#undef LOOP_BODY
				fflush(pfile);
				fflush(stdout);
			} // for num
		} // for size-ratio
	} // for intersection-ratio
	fclose(pfile);
}

void intersect_disk_data() {
	using namespace msis;

	mySet multiset([](vector<uint32_t> lhs,vector<uint32_t> rhs) {
		return lhs.size() <= rhs.size();
	});
	vector<uint32_t> out;

	size_t REPETITION;
	size_t CASES = 20;

	WallClockTimer timer;
	size_t time = 0;
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
//	vector<float> intersectionsratios = { 0.01, 0.05, 0.10, 0.20, 0.60, 0.80,
//			1.00 };
//	vector<uint32_t> sizeratios = { 1, 10, 20, 50, 100, 500, 1000, 5000, 10000 };
	vector<float> intersectionsratios = { 0.10, 0.50, 1.00 };
	vector<uint32_t> sizeratios = { 1, 100, 1000, 10000 };
#endif

	FILE *pfile = fopen("output_disk.csv", "w+");
	fprintf(pfile, "sr,ir,num,name,time\n");

	for (float ir : intersectionsratios) {
		printf("intersection ratio: \e[32m%3.0f%%\e[0m\n", ir * 100);
		for (uint32_t sr : sizeratios) {
			printf("  size ratio: \e[32m%4d\e[0m\n", sr);
			for (int i = 2; i < 11; i++) {
				printf("    num: \e[32m%2d\e[0m  ", i);

				for (uint32_t i = 0; i < CASES; i++) {
					multiset.clear();
					Set_collection::sequentially_load_file(multiset);
				}

//				// start scalar intersection
//				for (int j = 0; j < NUMSCALARFUNC; j++) {
//					if (i == 0 and sr > 100)
//						REPETITION = 3;
//					else
//						REPETITION = 100;
//					time = 0;
//					for (uint32_t i = 0; i < CASES; i++) {
//						multiset.clear();
//						Set_collection::sequentially_load_file(multiset);
//
//						// verification
////						auto it = multiset.begin();
////						vector<uint32_t> final_intersection = intersect(*it++,
////								*it++);
////						for (; it != multiset.end(); it++)
////							final_intersection = intersect(final_intersection,
////									*it);
//
//						timer.reset();
//						for (uint32_t howmany = 0; howmany < REPETITION;
//								++howmany) {
//							scalarFUNC[j](multiset, out);
//						}
////
////						if (out != final_intersection) {
////							std::cerr << "bad result!  " << std::endl;
////							return 1;
////						} else
////							printf("good!  ");
//
//						time += timer.split();
//						printf("%s: \e[31m%6.0f\e[0m  ",
//								NAMESCALARFUNC[j].c_str(),
//								(double) time / REPETITION / CASES);
//					}
//					fprintf(pfile, "%d,%.0f,%d,%s,%.0f\n", sr, ir * 100, i,
//							NAMESCALARFUNC[j].c_str(),
//							(double) time / REPETITION / CASES);
//				} // for intersection
//				Schlegel(multiset, timer, pfile, sr, ir * 100, i);
//				printf("\n");
//
//				std::string name;
//				Set_collection::pos = 0;
//#define LOOP_BODY(R,PRODUCT)\
//		time=0;\
//		for (uint32_t cases = 0; cases < CASES; cases++) {\
//			multiset.clear();\
//			Set_collection::sequentially_load_file(multiset);\
//			timer.reset();									\
//			for (uint32_t howmany = 0; howmany < REPETITION;\
//					++howmany) {\
//				BOOST_PP_CAT(\
//				BOOST_PP_CAT(\
//				BOOST_PP_CAT(\
//				BOOST_PP_CAT(\
//				BOOST_PP_SEQ_ELEM(0,PRODUCT),\
//				BOOST_PP_SEQ_ELEM(1,PRODUCT)),\
//				BOOST_PP_SEQ_ELEM(2,PRODUCT)),\
//				BOOST_PP_SEQ_ELEM(3,PRODUCT)),\
//				BOOST_PP_SEQ_ELEM(4,PRODUCT))(multiset,out);\
//			}\
//			time += timer.split();\
//		}\
//		name=(BOOST_PP_STRINGIZE(BOOST_PP_SEQ_ELEM(0,PRODUCT)));\
//		name.append("_").append(BOOST_PP_STRINGIZE(BOOST_PP_SEQ_ELEM(2,PRODUCT)));\
//		name.resize(name.size()-2);\
//		name.append(BOOST_PP_STRINGIZE(BOOST_PP_SEQ_ELEM(4,PRODUCT)));\
//		name.resize(name.size()-1);\
//		name.append("_").append(BOOST_PP_STRINGIZE(BOOST_PP_SEQ_ELEM(3,PRODUCT)));\
//		name=name.substr(6);\
//		fprintf(pfile, "%d,%.0f,%d,%s,%.0f\n", sr, ir * 100, i,\
//				name.c_str(),(double) time / REPETITION/ CASES);\
//
//				//printf("%s: \e[31m%6.0f\e[0m  \n", name.c_str(),\
//				//(double) time / REPETITION/CASES);
//
////				BOOST_PP_SEQ_FOR_EACH_PRODUCT(LOOP_BODY,
////				((BOOST_PP_SEQ_ELEM(0,METHOD)))/*SvS,s_SvS,sql,s_sql,max*/
////				((BOOST_PP_SEQ_ELEM(1,HEAD)))/*0:exact 1:rough*/
////				(/*(BOOST_PP_SEQ_ELEM(1,SEARCH))*/SEARCH)/*0:linear 1:gallop*/
////				(BOOST_PP_SEQ_REST_N(1,SIZE)/*BOOST_PP_SEQ_REST_N(5,SIZE)*/) /*4,8,16,32,64,128,256,512*/
////				(BOOST_PP_SEQ_POP_FRONT(END)));/*0:exact 1:rough 2:rough_plow*/
//
//				// TEST ALL THE METHOD
////				BOOST_PP_SEQ_FOR_EACH_PRODUCT(LOOP_BODY,
////						(/*(BOOST_PP_SEQ_ELEM(3,METHOD))*/METHOD)/*SvS,s_SvS,sql,s_sql,max*/
////						((BOOST_PP_SEQ_ELEM(0,HEAD)))/*0:exact 1:rough*/
////						(/*(BOOST_PP_SEQ_ELEM(1,SEARCH))*/SEARCH)/*0:linear 1:gallop*/
////						(SIZE/*BOOST_PP_SEQ_REST_N(5,SIZE)*/) /*4,8,16,32,64,128,256,512*/
////						((BOOST_PP_SEQ_ELEM(0,END))/*BOOST_PP_SEQ_POP_FRONT(END)*/));/*0:exact 1:rough 2:rough_plow*/
////
////				Set_collection::pos = 0;
////				BOOST_PP_SEQ_FOR_EACH_PRODUCT(LOOP_BODY,
////						(/*(BOOST_PP_SEQ_ELEM(3,METHOD))*/METHOD)/*SvS,s_SvS,sql,s_sql,max*/
////						((BOOST_PP_SEQ_ELEM(1,HEAD)))/*0:exact 1:rough*/
////						(/*(BOOST_PP_SEQ_ELEM(1,SEARCH))*/SEARCH)/*0:linear 1:gallop*/
////						(SIZE/*BOOST_PP_SEQ_REST_N(5,SIZE)*/) /*4,8,16,32,64,128,256,512*/
////						((BOOST_PP_SEQ_ELEM(1,END))));/*0:exact 1:rough 2:rough_plow*/
////
////				Set_collection::pos = 0;
////				BOOST_PP_SEQ_FOR_EACH_PRODUCT(LOOP_BODY,
////						(/*(BOOST_PP_SEQ_ELEM(3,METHOD))*/METHOD)/*SvS,s_SvS,sql,s_sql,max*/
////						((BOOST_PP_SEQ_ELEM(1,HEAD)))/*0:exact 1:rough*/
////						(/*(BOOST_PP_SEQ_ELEM(1,SEARCH))*/SEARCH)/*0:linear 1:gallop*/
////						(BOOST_PP_SEQ_REST_N(1,SIZE)/*BOOST_PP_SEQ_REST_N(5,SIZE)*/) /*4,8,16,32,64,128,256,512*/
////						((BOOST_PP_SEQ_ELEM(2,END))));/*0:exact 1:rough 2:rough_plow*/
//#undef LOOP_BODY
//				fflush(pfile);
//				fflush(stdout);
			} // for num
		} // for size-ratio
	} // for intersection-ratio
	fclose(pfile);
}

int main() {
	mySet multiset([](vector<uint32_t> lhs,vector<uint32_t> rhs) {
		return lhs.size() <= rhs.size();
	});
	while (Set_collection::pos != -1)
		Set_collection::sequentially_load_file(multiset);
//	intersect_disk_data();
}
