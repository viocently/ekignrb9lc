#include<iostream>
#include<fstream>
#include<vector>
#include<bitset>
#include<algorithm>
#include<thread>
#include<mutex>
#include<future>
#include<sstream>
#include<malloc.h>
#include<array>
#include"log.h"
#include"deg.h"
#include"thread_pool.h"
#include"trivium.h"
#include"triviumtrack.h"
#include"newtrivium.h"
#include"triviumcallback.h"
#include"triviumframework.h"

using namespace std;
using namespace thread_pool;
mutex stream_mutex;
mutex map_mutex;
mutex term_mutex;
int THREAD_NUM = 64;

array<vector<bitset<80>>, 288> conskCofs;

void preCompute(int rounds, bitset<80>& cube)
{
	for (int i = 0; i < 288; i++)
		conskCofs[i].clear();

	vector<bitset<288>> nonZeroTrackReps;
	triviumTrack(cube, nonZeroTrackReps);


	vector<bitset<288>> varsTrackReps;
	triviumVarsTrack(cube, varsTrackReps, nonZeroTrackReps);

	bitset<288> nonZeroFlag, varsFlag;
	setRoundFlags(rounds, cube, nonZeroFlag, varsFlag, nonZeroTrackReps, varsTrackReps);

	for (int i = 0; i < 288; i++)
		if (nonZeroFlag[i] == 1 && varsFlag[i] == 0)
		{
			map<bitset<80>, int, CMPS<80>> kcofMap;
			vector<bitset<80>> kcofs;
			bitset<288> lastState;
			lastState[i] = 1;
			bitset<80> cube0;
			MidSolutionCounter(rounds, cube0, lastState, kcofMap, 120, 2);
			for (auto& x : kcofMap)
				if (x.second % 2)
					kcofs.emplace_back(x.first);
			conskCofs[i] = kcofs;
		}

	logger("Precomputation for round " + to_string(rounds) + " Completed.");
}

void printSol(int rounds, const bitset<288>& vector, const map<bitset<80>,
	int, CMPS<80>> &solutions)
{
	logger(__func__ + string(": ") + to_string(rounds) + string("\t") +
		vector.to_string());

	string path = string("TERM/") + to_string(rounds) + string(".txt");
	ofstream os;
	os.open(path, ios::out | ios::app);
	os << rounds << "\t" << vector << endl;
	for (auto & it : solutions)
		os << it.first << "\t" << it.second << "\n";
	os << endl;
	os.close();
}

void resolveTerms(vector<bitset<80>>& kcofs, stringstream& ss)
{
	auto sep = "";
	ss << "(";
	for (auto& x : kcofs)
	{
		ss << sep;
		if (x.count() == 0)
			ss << "1";
		else
			for (int i = 0; i < 80; i++)
				if (x[i])
					ss << "k" << i;

		sep = "+";
	}
	ss << ")";
}


void filterTerms(const bitset<80>& cube, int rounds, vector<bitset<288>>& terms)
{
	int size0 = terms.size();
	vector<bitset<288>> tmp(terms.begin(), terms.end());
	terms.clear();
	for (auto& it : tmp)
	{
		auto d = computeDegree(cube, rounds, it);
		if (d >= cube.count())
			terms.push_back(it);
	}
	int size1 = terms.size();

	logger(__func__ + string(": ") + to_string(size0) + string("\t") +
		to_string(size1));
}

// this function is only called from SolutionFunc to get solutions from counterMap
void SecondSolutionFunc(int rounds,int midround, bitset<80>& cube, map<bitset<288>, int, CMPS<288> >& counterMap, map<bitset<80>, int, CMPS<80> >& solutions, vector<bitset<288>>& expandTerms, float time, int singlethread)
{
	solutions.clear();

	// first sort out midsolutions
	vector<bitset<288> > midsolutions;
	for (auto& x : counterMap)
	{
		if (x.second % 2 == 1)
			midsolutions.emplace_back(x.first);
		else
		{
			int midsols = 0;
			auto status = newMidSolutionCounter(midround, cube, x.first, midsols, time, singlethread);
			if (status == EXPAND || status == NOSOLUTION)
			{
				logger(__func__ + string(" : midsolution count error."));
				exit(-1);
			}

			int sols = x.second;
			if ((sols / midsols) % 2 == 1)
				midsolutions.emplace_back(x.first);
		}
	}

	// get final solutions from midsolutions
	// nonzero track
	vector<bitset<288>> nonZeroTrackReps;
	triviumTrack(cube, nonZeroTrackReps);

	// variables track
	vector<bitset<288>> varsTrackReps;
	triviumVarsTrack(cube, varsTrackReps, nonZeroTrackReps);

	bitset<288> nonZeroFlag;
	bitset<288> varsFlag;
	setRoundFlags(midround, cube, nonZeroFlag, varsFlag, nonZeroTrackReps, varsTrackReps);

	stringstream ss;

	for (auto& x : midsolutions)
	{
		/*
		bitset<288> midcube;
		vector<bitset<80>> midconskcofs = { bitset<80>(0) };

		for (int i = 0; i < 288; i++)
			if (varsFlag[i])
				midcube[i] = x[i];
			else if (nonZeroFlag[i] && x[i])
				midconskcofs = vecMul<80>(midconskcofs, conskCofs[i]);


		map<bitset<80>, int, CMPS<80>> midcubesols;
		auto status = MidSolutionCounter(midround, cube, midcube, midcubesols, 120, 2);
		if (status == NOSOLUTION)
		{
			logger(__func__ + string(" : get final solutions error."));
			exit(-1);
		}
		else if (status == EXPAND)
		{
			logger(__func__ + string(" expandTerm added : ") + x.to_string());

			lock_guard<mutex> guard(term_mutex);
			expandTerms.emplace_back(x);
		}

		vector<bitset<80>> midcubekcofs;
		for (auto& y : midcubesols)
			if (y.second % 2)
				midcubekcofs.emplace_back(y.first);

		vector<bitset<80>> midkcofs = vecMul<80>(midcubekcofs, midconskcofs);

		
		for (auto& y : midkcofs)
			solutions[y]++;
		
		*/
		bitset<288> midcube;
		for (int i = 0; i < 288; i++)
			if (varsFlag[i])
				midcube[i] = x[i];
		map<bitset<80>, int, CMPS<80>> midcubesols;
		auto status = MidSolutionCounter(midround, cube, midcube, midcubesols, 120, 2);
		if (status == NOSOLUTION)
		{
			logger(__func__ + string(" : get final solutions error."));
			exit(-1);
		}
		else if (status == EXPAND)
		{
			logger(__func__ + string(" expandTerm added : ") + x.to_string());

			lock_guard<mutex> guard(term_mutex);
			expandTerms.emplace_back(x);
		}

		vector<bitset<80>> midcubekcofs;
		for (auto& y : midcubesols)
			if (y.second % 2)
				midcubekcofs.emplace_back(y.first);

		if (midcubekcofs.size() > 0)
		{
			resolveTerms(midcubekcofs, ss);

			for (int i = 0; i < 288; i++)
				if (nonZeroFlag[i] && varsFlag[i] == 0 && x[i])
				{
					ss << "*";
					resolveTerms(conskCofs[i], ss);
				}

			ss << endl;
		}
				
	}

	//lock_guard<mutex> guard(stream_mutex);
	loggerStream(ss, rounds);

}


STATUS SolutionFunc(int rounds, int midround, bitset<80>& cube, bitset<288>& last, map<bitset<80>, int, CMPS<80>>& solutions, vector<bitset<288>>& expandTerms, float time, int singlethread)
{
	if (midround > rounds || rounds < 0 || midround < 0)
	{
		logger(__func__ + string(" : parameter error."));
		exit(-1);
	}

	map<bitset<288>, int, CMPS<288> > counterMap;
	auto status = newMidSolutionCounter(rounds, cube, last, counterMap, time, singlethread, midround);

	if (status == SOLUTION)
		SecondSolutionFunc(rounds,midround, cube, counterMap, solutions, expandTerms, time, singlethread);

	return status;
}

// below is threadworker
void SolutionSearcherWorker(bitset<288>& vec, int ROUND, bitset<80>& cube,
	vector<bitset<288>>& layerTerms, float time, int singlethread, int midround, map<bitset<80>, int, CMPS<80>>& totalSolutions, vector<bitset<288>>& expandTerms)
{
	map<bitset<80>, int, CMPS<80>> solutions;

	auto status = SolutionFunc(ROUND, midround, cube, vec, solutions, expandTerms, time, singlethread);

	/*stringstream ss;
	for (auto& x : solutions)
		ss << x.first << " " << x.second << endl;*/

	if (status == EXPAND)
	{
		lock_guard<mutex> guard(map_mutex);
		layerTerms.push_back(vec);
	}
	else if (status == SOLUTION)
	{
		/*lock_guard<mutex> guard(stream_mutex);
		for (auto& x : solutions)
			totalSolutions[x.first] = (totalSolutions[x.first] + x.second) % 2;*/
		//loggerStream(ss, ROUND);
	}
}



int main()
{

	vector<int> cubeIndex = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                  21, 22, 23, 24,25,26, 28,30, 32, 34, 36, 38, 40, 42, 45, 47, 49, 51, 53,
55, 57, 60, 62, 64, 66, 68, 70, 72, 77, 75, 79	};


	bitset<80> cube{ 0 };
	for (auto& it : cubeIndex)
		cube[it] = 1;

	cout << cube << endl;


	int midround = 90;
	ThreadPool thread_pool{};


	int ROUND = 848;
	int step = 300;
	int expandStep = 20;
	float expandTime = 3000;

	float time = 0;

	vector<bitset<288>> terms;
	vector<bitset<288>> layerTerms;

	BackExpandPolynomial(step, terms);

	ROUND = ROUND - step;
	filterTerms(cube, ROUND, terms);

	logger(__func__ + string("First expand"));

	int singlethread = 2;


	map<bitset<80>, int, CMPS<80>> totalSolutions;
	vector<bitset<288> > expandTerms;

	if (ROUND < CALLBACKBOUND)
		newExpandFunc(cube, terms, ROUND, expandStep, thread_pool, expandTime, singlethread);
	else
		expandFunc(cube, terms, ROUND, expandTime, singlethread);


	while (true)
	{

		if (ROUND > 600)
			time = 40;
		else if (ROUND > 500)
			time = 80;
		else if (ROUND > 400)
			time = 160;
		else if (ROUND > 300)
			time = 320;
		else if (ROUND > 250)
			time = 640;
		else if (ROUND > 100)
			time = 1200;
		else if (ROUND > 20)
			time = 3600;
		else
			time = 0;


		layerTerms.clear();
		preCompute(midround, cube);
		vector<std::future<void>> futures;
		for (auto& it : terms)
			futures.emplace_back(thread_pool.Submit(SolutionSearcherWorker, it,
				ROUND, ref(cube), ref(layerTerms), time, singlethread, midround, ref(totalSolutions), ref(expandTerms)));
		for (auto& it : futures)
			it.get();

		filterMap<80>(totalSolutions);
		filterExpandTerms<288>(expandTerms);

#ifdef _WIN32
#else
		showProcessMemUsage();
		malloc_trim(0);
		showProcessMemUsage();
#endif


		logger(string("layerTermsSize ") + to_string(layerTerms.size()));

		logger(string("expandTermsSize ") + to_string(expandTerms.size()));

		if (layerTerms.size() == 0 && expandTerms.size() == 0)
			break;
		else if (layerTerms.size() == 0)
		{
			logger("Start to handle expandTerms.");
			ROUND = midround;
			midround = ROUND / 2;
			layerTerms = expandTerms;
			filterTerms(cube, ROUND, layerTerms);
			expandTerms.clear();
		}

		if (ROUND < 10)
		{
			logger("ROUND is to small.");
			exit(-1);
		}

		terms = layerTerms;


		if (ROUND < FASTBOUND)
			expandStep = 20;

		if (ROUND - expandStep < midround)
			expandStep = ROUND - midround;

		// second Expand
		if (ROUND < CALLBACKBOUND)
			newExpandFunc(cube, terms, ROUND, expandStep, thread_pool, expandTime, singlethread);
		else
			expandFunc(cube, terms, ROUND, expandTime, singlethread);

		// if ROUND == midround,combine expandterms with terms since they are all in midround
		if (ROUND == midround)
		{
			logger("current ROUND : " + to_string(ROUND));
			logger("start combine expandTerms with terms.");
			map<bitset<288>, int, CMPS<288>> termCounter;
			for (auto& x : terms)
				termCounter[x] ++;
			for (auto& x : expandTerms)
				termCounter[x] ++;

			expandTerms.clear();
			terms.clear();
			for (auto& x : termCounter)
				if (x.second % 2 == 1)
					terms.emplace_back(x.first);

			midround = ROUND / 2;
		}


	}
	cout << "Success!" << endl;

	/*
	cout << "Output to file." << endl;

	// output result to file
	string path = string("TERM/") + string("res.txt");
	ofstream os;
	os.open(path, ios::out | ios::app);
	for (auto& it : totalSolutions)
	{
		bitset<80> term = it.first;
		if (term.count() == 0)
		{
			os << "1" << endl;
		}
		else
		{
			for (int i = 0; i < 80; i++)
				if (term[i])
					os << "k" << i;
			os << endl;
		}
	}
	os << endl;
	os.close();


	cout << "Output to file finished." << endl;*/
}
