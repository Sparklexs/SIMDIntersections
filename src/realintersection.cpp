/**
 * This code is released under the
 * Apache License Version 2.0 http://www.apache.org/licenses/.
 *
 * sxs: this method first generates synthetic data and makes some statistics,
 * then it tests skipping intersection and user-requested intersection
 * before every intersection, it pre-run all the data to ensure the correctness.
 * (c) Daniel Lemire, http://lemire.me/en/
 */
#include "common.h"
#include "intersectionfactory.h"
#include "timer.h"
#include "synthetic.h"
#include "util.h"
#include "skipping.h"

void printusage() {
	cout << " Try ./realintersection -r 40" << endl;
	cout << " Use the -s flag to specify just some scheme, choose from: "
			<< endl;
	for (string x : allRealNames())
		cout << " " << x << endl;
	cout << " Separate the schemes by a comma (e.g. -s schlegel,danscalar). "
			<< endl;
}

int main(int argc, char **argv) {
	size_t howmany = 0;
	bool natemode = false;
	bool safe = true;
	size_t loop = 1000;
	uint32_t S = 12;
	uint32_t ratio = 1;
	bool skipping = false;
	uint32_t skipgaplog = 0;
	map<string, intersectionfunction> myschemes(realschemes);
	int c;
	while ((c = getopt(argc, argv, "uns:S:m:l:r:hk:")) != -1)
		switch (c) {
		case 'k':
			skipping = true;
			skipgaplog = atoi(optarg);
			if ((S < 1) or (S > 31)) {
				cerr << "Skip param needs to be within [1,31]." << endl;
				printusage();
				return -1;
			}
			cout << "#using skip steps of " << (1 << skipgaplog) << endl;
			myschemes.clear();
			break;
		case 'u':
			safe = false;
			break;
		case 'n':
			natemode = true;
			break;
		case 'h':
			printusage();
			return 0;
		case 'S':
			S = atoi(optarg);
			if ((S < 1) or (S > 31)) {
				printusage();
				return -1;
			}
			break;
		case 'm':
			howmany = atoi(optarg);
			if (howmany < 1) {
				printusage();
				return -1;
			}
			break;
		case 's':
			myschemes.clear();
			{
				const string codecsstr(optarg);
				const vector<string> codecslst = split(codecsstr, ",:;");
				for (auto i = codecslst.begin(); i != codecslst.end(); ++i) {
					if (realschemes.find(*i) == realschemes.end()) {
						cerr << " Warning!!! Warning: unrecognized: " << *i
								<< endl;
						printusage();
						return -1;

					} else {
						const auto K = realschemes.find(*i);
						const std::string name = K->first;
						const intersectionfunction fn = K->second;
						myschemes[name] = fn;
					}
				}
			}
			break;
		case 'l':
			loop = atoi(optarg);
			if (loop < 1) {
				printusage();
				return -1;
			}
			break;
		case 'r':
			ratio = atoi(optarg);
			if (ratio < 1) {
				printusage();
				return -1;
			}
			break;
		default:
			abort();
		}
	if (howmany == 0) {
		if (natemode)
			howmany = 1;
		else
			howmany = 20;
	}
	uint32_t MaxBit = 31;
	if (natemode) {
		MaxBit = S + 7;
		if (MaxBit > 31)
			MaxBit = 31;
	}
//	cout << "# algo: ";
//	for (auto algo : myschemes) {
//		cout << algo.first << " ";
//	}
//	cout << endl;
	const uint32_t N = 1U << S;
	if (!natemode)
		if (ratio != 1) {
			cout << "# ratio = " << ratio << endl;
			cout << "# size of largest array = " << N << endl;
			cout << "# size of smallest array = " << N / ratio << endl;
		} else {
			cout << "# size of arrays = " << N << endl;
		}
	else
		cout << "# in natemode" << endl;

	ClusteredDataGenerator cdg;
	WallClockTimer z;
	size_t bogus = 0;
	size_t volume = 0;
	uint64_t time = 0;
	cout << "# first column is relative size of intersection" << endl;
	if (ratio > 1) {
		cout
				<< "# next two are estimated average bits per int for differential coding"
				<< endl;
	} else {
		cout
				<< "# next is estimated average bits per int for differential coding"
				<< endl;
	}

	cout
			<< "# other columns display speed in mis when computing the intersection"
			<< endl << endl;

	std::cout << "is\t" << "bpi\t\t";
//	for (auto ele : myschemes)
//		std::cout << ele.first << "\t";
	std::cout << std::endl;

	for (uint32_t gap = 0; gap + S <= MaxBit; gap += 1) {
		vector<vector<uint32_t> > data;
		for (size_t zz = 0; zz < howmany; ++zz) {
			if (natemode) {
				data.push_back(
						cdg.generateClustered((1U << (MaxBit - gap)) / ratio,
								1U << MaxBit));
				data.push_back(
						cdg.generateClustered((1U << (MaxBit - gap)),
								1U << MaxBit));
			} else {
				data.push_back(
						cdg.generateClustered(N / ratio, 1U << (gap + S)));
				data.push_back(cdg.generateClustered(N, 1U << (gap + S)));
			}
		}
		size_t intersize = 0;
		size_t smallestsize = 0;
		for (size_t k = 0; k < howmany; k++) {
			intersize += classicalintersectioncardinality(&data[2 * k][0],
					data[2 * k].size(), &data[2 * k + 1][0],
					data[2 * k + 1].size());
			smallestsize +=
					data[2 * k + 1].size() < data[2 * k].size() ?
							data[2 * k + 1].size() : data[2 * k].size();
		}

		cout << "\e[31m" << std::fixed << std::setprecision(3)
				<< static_cast<double>(intersize)
						/ static_cast<double>(smallestsize) << "\t";
		cout.flush();
		if (ratio > 1) {
			if (natemode) {
				cout
						<< log(
								1
										+ (1U << MaxBit) * 1.0
												/ ((1U << (MaxBit - gap))
														/ ratio)) / log(2)
						<< "\t";
				cout
						<< log(
								1
										+ (1U << MaxBit) * 1.0
												/ ((1U << (MaxBit - gap))))
								/ log(2) << "\t";
			} else {
				cout << log(1 + (1U << (gap + S)) * 1.0 / (N / ratio)) / log(2)
						<< "\t";
				cout << log(1 + (1U << (gap + S)) * 1.0 / (N)) / log(2) << "\t";
			}
		} else {
			if (natemode) {
				cout
						<< log(
								1
										+ (1U << MaxBit) * 1.0
												/ ((1U << (MaxBit - gap))))
								/ log(2) << "\t";
			} else {
				cout << log(1 + (1U << (gap + S)) * 1.0 / (N)) / log(2) << "\t";
			}

		}
		cout << "\n\e[0m";
		cout.flush();
		size_t maxsize = 0;
		for (auto x : data)
			if (x.size() > maxsize)
				maxsize = x.size();
		vector<uint32_t> buffer((maxsize + 15) / 16 * 16);
		/**
		 * Skipping is a standard technique in IR. We test it here.
		 */
		if (skipping) {
			vector<shared_ptr<Skipping> > sdata;
			for (vector<uint32_t> & x : data)
				sdata.emplace_back(
						shared_ptr<Skipping>(
								new Skipping(skipgaplog, x.data(),
										static_cast<uint32_t>(x.size()))));
			for (size_t k = 0; k < 2 * howmany; k += 2) {
				vector<uint32_t> out(buffer.size());
				size_t correctanswer = classicalintersection(&data[k][0],
						data[k].size(), &data[k + 1][0], data[k + 1].size(),
						&out[0]);
				out.resize(correctanswer);
				vector<uint32_t> out2(buffer.size());
				size_t thisschemesanswer = sdata[k]->intersect(*sdata[k + 1],
						&out2[0]);
				out2.resize(thisschemesanswer);
				if (out != out2) {
					if (thisschemesanswer != correctanswer) {
						cerr << "expecting cardinality of " << correctanswer;
						cerr << " got " << thisschemesanswer << "." << endl;
						int times = 0;
						for (size_t jj = 0;
								(jj < thisschemesanswer) && (jj < correctanswer)
										&& (times < 10); ++jj) {
							if (out[jj] != out2[jj]) {
								cout << " index = " << jj << " expected "
										<< out[jj] << " but got " << out2[jj]
										<< endl;
								++times;
							}
						}
						if (times == 0)
							cout
									<< "content is correct, but incomplete/too much"
									<< endl;
					} else {
						cerr << "Same cardinality " << correctanswer
								<< ". Good. " << endl;
						for (size_t jj = 0; jj < correctanswer; ++jj)
							if (out[jj] != out2[jj]) {
								cerr << "Differ at " << jj << " got "
										<< out2[jj] << " should find "
										<< out[jj] << endl;
								break;
							}
					}
					throw runtime_error("bug");
				}
			}
			volume = 0;
			z.reset();
			for (size_t L = 0; L < loop; ++L) {

				for (size_t k = 0; k < 2 * howmany; k += 2) {
					volume += data[k].size();
					volume += data[k + 1].size();
					bogus += sdata[k]->intersect(*sdata[k + 1], &buffer[0]);
				}

			}
			time = z.split();
			cout << std::setprecision(0)
					<< static_cast<double>(volume) / static_cast<double>(time)
					<< "\t";
			cout.flush();
		}
		/**
		 * End of Skipping
		 */
		for (auto algo : myschemes) {
			if (safe and buggyschemes.find(algo.first) == buggyschemes.end())
				for (size_t k = 0; k < 2 * howmany; k += 2) {
					vector<uint32_t> out(buffer.size());
					size_t correctanswer = classicalintersection(&data[k][0],
							data[k].size(), &data[k + 1][0], data[k + 1].size(),
							&out[0]);
					out.resize(correctanswer);
					vector<uint32_t> out2(buffer.size());
					size_t thisschemesanswer = algo.second(&data[k][0],
							data[k].size(), &data[k + 1][0], data[k + 1].size(),
							&out2[0]);
					out2.resize(thisschemesanswer);
					if (out != out2) {
						if (thisschemesanswer != correctanswer) {
							cerr << "expecting cardinality of "
									<< correctanswer;
							cerr << " got " << thisschemesanswer << "." << endl;
							if (correctanswer < 10)
								for (uint32_t x : out)
									cerr << x << endl;
						} else {
							cerr << "Same cardinality " << correctanswer
									<< ". Good. " << endl;
							for (size_t jj = 0; jj < correctanswer; ++jj)
								if (out[jj] != out2[jj]) {
									cerr << "Differ at " << jj << " got "
											<< out2[jj] << " should find "
											<< out[jj] << endl;
									break;
								}
						}
						throw runtime_error("bug");
					}
				}
			volume = 0;
			z.reset();
			for (size_t L = 0; L < loop; ++L) {

				for (size_t k = 0; k < 2 * howmany; k += 2) {
					volume += data[k].size();
					volume += data[k + 1].size();
					bogus += algo.second(&data[k][0], data[k].size(),
							&data[k + 1][0], data[k + 1].size(), &buffer[0]);
				}

			}
			time = z.split();
			cout << std::setw(10) << "(" + algo.first + ")" << std::setw(4)
					<< std::setprecision(0) << "\e[32m"
					<< static_cast<double>(volume) / static_cast<double>(time)
					<< "\e[0m  ";
			cout.flush();
		}

		cout << endl;

	}
	cout << "# bogus = " << bogus << endl;
}
