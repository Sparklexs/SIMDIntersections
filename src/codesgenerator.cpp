/*
 * codesgenerator.cpp
 *
 *  Created on: 2017Äê3ÔÂ28ÈÕ
 *      Author: John
 */
#include <stddef.h>
#include <stdlib.h>
#include <sys/_stdint.h>
#include <sys/time.h>
#include <cstring>
#include <iostream>

using namespace std;

void bisection(uint32_t, uint32_t, uint32_t);
void simdor(uint32_t, uint32_t);

void generate(uint32_t blocksize, uint32_t step) {
	cout << "long simdgallop_" << blocksize << "_" << step << "_rough"
			<< "(int *foundp, UINT4 goal, const UINT4 *target,\n\
		long ntargets) {\n\
	const UINT4 *end_target, *stop_target, *init_target;\n\
\n\
	long low_offset = 0, mid_offset, high_offset = 1;\n\
	long pos;\n\
\n\
	init_target = target;\n\
	stop_target = target + ntargets - V"
			<< blocksize
			<< "_BLOCKSIZE;\n\
	end_target = target + ntargets;\n\
\n\
	if (_UNLIKELY(target >= stop_target)) {\n\
		if ((pos = Intersection_find_scalar(goal, target, ntargets)) <= ntargets\n\
				&& target[pos] == goal)\n\
			*foundp = 1;\n\
		else\n\
			*foundp = 0;\n\
		return pos;\n\
	}\n\
\n\
	if (_UNLIKELY(target[V"
			<< blocksize << "_BLOCKSIZE - 1] < goal)) {\n\
		if (target + V"
			<< blocksize << "_BLOCKSIZE > stop_target) {\n\
			pos = V"
			<< blocksize
			<< "_BLOCKSIZE\n\
					+ Intersection_find_scalar(goal, target + V"
			<< blocksize << "_BLOCKSIZE,\n\
							ntargets - V" << blocksize
			<< "_BLOCKSIZE);\n\
			if (pos <= ntargets && target[pos] == goal)\n\
				*foundp = 1;\n\
			else\n\
				*foundp = 0;\n\
			return pos;\n\
		}\n\
		/* Galloping search */\n\
		while (target[V"
			<< blocksize << "_BLOCKSIZE * high_offset + V" << blocksize
			<< "_BLOCKSIZE - 1] < goal) {\n\
			if (target + (high_offset << 1) * V"
			<< blocksize
			<< "_BLOCKSIZE <= stop_target) {\n\
				//doubled\n\
				low_offset = high_offset;\n\
				high_offset <<= 1;\n\
			} else if (target + V"
			<< blocksize
			<< "_BLOCKSIZE * (high_offset + 1)\n\
					<= stop_target) {\n\
				//more than one block left\n\
				high_offset = (stop_target - target) / V"
			<< blocksize << "_BLOCKSIZE;\n\
				if (target[V" << blocksize
			<< "_BLOCKSIZE * high_offset + V" << blocksize
			<< "_BLOCKSIZE - 1]\n\
						< goal) {\n\
					target += V"
			<< blocksize
			<< "_BLOCKSIZE * high_offset;\n\
					pos = Intersection_find_scalar(goal, target,\n\
							(end_target - target));\n\
					if (pos + V"
			<< blocksize
			<< "_BLOCKSIZE * high_offset <= ntargets\n\
							&& target[pos] == goal)\n\
						*foundp = 1;\n\
					else\n\
						*foundp = 0;\n\
					return pos + V"
			<< blocksize
			<< "_BLOCKSIZE * high_offset;\n\
				} else\n\
					break;\n\
			} else {\n\
				//only one block left\n\
				target += V"
			<< blocksize
			<< "_BLOCKSIZE * high_offset;\n\
\n\
				pos = Intersection_find_scalar(goal, target,\n\
						(end_target - target));\n\
				if (pos + V"
			<< blocksize
			<< "_BLOCKSIZE * high_offset <= ntargets\n\
						&& target[pos] == goal)\n\
					*foundp = 1;\n\
				else\n\
					*foundp = 0;\n\
				return pos + V"
			<< blocksize
			<< "_BLOCKSIZE * high_offset;\n\
			}\n\
		}\n\
		while (low_offset < high_offset) {\n\
			mid_offset = (low_offset + high_offset) / 2;\n\
			if (target[V"
			<< blocksize << "_BLOCKSIZE * mid_offset + V" << blocksize
			<< "_BLOCKSIZE - 1]\n\
					< goal) {\n\
				low_offset = mid_offset + 1;\n\
			} else {\n\
				high_offset = mid_offset;\n\
			}\n\
		}\n\
		target += V"
			<< blocksize << "_BLOCKSIZE * high_offset;\n\
	}\n";
	bisection(0, blocksize, step);

	cout << "__m128i Match = _mm_set1_epi32(goal);\n";
	cout << "__m128i F0=";
	simdor(0, step / 4 - 1);
	cout << ";\n";
	cout
			<< "	if (_mm_testz_si128(F0, F0))\n\
		*foundp = 0;\n\
	else\n\
		*foundp = 1;\n\
	return (target - init_target);\n\
}"
			<< endl;
}

void bisection(uint32_t start, uint32_t end, uint32_t step) {
	if (end - start == step)
		return;
	if (end - start > 2 * step) {
		cout << "if(target[" << (end + start) / 2 - 1 << "] >= goal){" << endl;
		bisection(start, (end + start) / 2, step);
		cout << "}" << endl;

		cout << "else{" << endl;
		bisection((end + start) / 2, end, step);
		cout << "}" << endl;

	}
	// handle the last step inside if-loop
	else if (start == 0) {
		cout << "if(target[" << (end + start) / 2 - 1 << "] < goal){" << endl;
		cout << "target += " << (end + start) / 2 << ";}" << endl;
	} else {
		cout << "if(target[" << (end + start) / 2 - 1 << "] >= goal){" << endl;
		cout << "target += " << start << ";}" << endl;
		cout << "else{" << endl;
		cout << "target += " << (start + end) / 2 << ";}" << endl;
	}
}

void simdor(uint32_t start, uint32_t end) {
	if (start != end) {
		cout << "_mm_or_si128(";
		simdor(start, (start + end) / 2);
		cout << ",";
		simdor((start + end) / 2 + 1, end);
		cout << ")";
	} else {
		cout << "_mm_cmpeq_epi32(_mm_lddqu_si128((__m128i *) target + " << start
				<< "),Match)";
	}
}

int main() {
//	ifstream ifs("time.json");
//		ofstream ofs("time_4096.json");
//		if (!ifs.is_open())
//			return 1;
//		string line;
//		while (!ifs.eof()) {
//			std::getline(ifs, line);
//
//			if (line.find("\"minlength\": 4096") != std::string::npos) {
//				line.resize(line.find("minlength") - 3);
//				line.append("}");
//				ofs << line << std::endl;
//			}
//		}
//		ofs.flush();
//		ofs.close();

	cout
			<< "/*\n\
 * msis_gallop.hpp\n\
 *\n\
 *  Created on: 2017.3.28\n\
 *      Author: John\n\
 */\n\
\n\
#ifndef INCLUDE_MSIS_GALLOP_ROUGH_HPP_\n\
#define INCLUDE_MSIS_GALLOP_ROUGH_HPP_\n\
\n\
#include \"msis_linear.hpp\"\n\
\n\
namespace msis/*MultiSet InterSection*/{\n"
			<< endl;
	for (uint32_t i = 4; i <= 512; i *= 2) {
		for (uint32_t j = 4; j <= i; j *= 2) {
			generate(i, j);
		}
	}
	cout
			<< "typedef long (*searchFUNC)(int *foundp, UINT4 goal, const UINT4 *target,\n\
		long ntargets);\n\
searchFUNC optSearchFunc[32] = { nullptr/*0*/, simdgallop_4_4_rough/*1*/,\n\
		simdgallop_8_4_rough/*2*/, simdgallop_32_8_rough/*4*/,\n\
		simdgallop_64_8_rough/*8*/, simdgallop_256_16_rough/*16*/,\n\
		simdgallop_256_16_rough/*32*/, simdgallop_512_16_rough/*64*/,\n\
		simdgallop_512_32_rough/*128*/, simdgallop_512_32_rough/*256*/,\n\
		simdgallop_512_32_rough/*512*/, simdgallop_256_32_rough/*1024*/,\n\
		simdgallop_512_64_rough/*2048*/, simdgallop_512_64_rough/*4096*/,\n\
		simdgallop_512_64_rough/*8192*/, simdgallop_512_64_rough/*8192*/,\n\
		simdgallop_512_64_rough/*8192*/, simdgallop_512_64_rough/*8192*/,\n\
		simdgallop_512_64_rough/*8192*/, simdgallop_512_64_rough/*8192*/,\n\
		simdgallop_512_64_rough/*8192*/, simdgallop_512_64_rough/*8192*/,\n\
		simdgallop_512_64_rough/*8192*/, simdgallop_512_64_rough/*8192*/,\n\
		simdgallop_512_64_rough/*8192*/, simdgallop_512_64_rough/*8192*/,\n\
		simdgallop_512_64_rough/*8192*/, simdgallop_512_64_rough/*8192*/,\n\
		simdgallop_512_64_rough/*8192*/, simdgallop_512_64_rough/*8192*/,\n\
		simdgallop_512_64_rough/*8192*/, simdgallop_512_64_rough /*8192*/};"
			<< std::endl;
	cout << "}\n" << "#endif /* INCLUDE_MSIS_GALLOP_ROUGH_HPP_ */" << endl;
}

