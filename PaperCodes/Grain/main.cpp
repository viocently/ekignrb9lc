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
#include<sstream>
#include"log.h"
#include"thread_pool.h"
#include"grain.h"
#include"graintrack.h"
#include"newgrain.h"
#include"graincallback.h"
#include"grainframework.h"
#include"node.h"


using namespace std;
using namespace thread_pool;
mutex stream_mutex;
mutex map_mutex;
mutex term_mutex;
int THREAD_NUM = 64;

array<vector<bitset<128>>, 256> conskCofs;

void preCompute(int rounds, bitset<96> & cube)
{
	for (int i = 0; i < 256; i++)
		conskCofs[i].clear();

	bitset<256> flag1;
	flag1.set();



	vector<bitset<256> > nonZeroTrackReps;
	grainTrack(cube, nonZeroTrackReps);


	vector<bitset<256> > varsTrackReps;
	grainVarsTrack(cube, varsTrackReps, nonZeroTrackReps);

	bitset<256> nonZeroFlag = (rounds < nonZeroTrackReps.size()) ? nonZeroTrackReps[rounds] : flag1;
	bitset<256> varsFlag = (rounds < varsTrackReps.size()) ? varsTrackReps[rounds] : flag1;

	for (int i = 0; i < 256;i++)
		if (nonZeroFlag[i] == 1 && varsFlag[i] == 0)
		{
			map<bitset<128>, int, CMPS<128>> kcofMap;
			vector<bitset<128>> kcofs;
			bitset<256> lastState;
			lastState[i] = 1;
			bitset<96> cube0;
			MidSolutionCounter(rounds, cube0, lastState, kcofMap, 120, 2);
			for (auto& x : kcofMap)
				if (x.second % 2)
					kcofs.emplace_back(x.first);
			conskCofs[i] = kcofs;
		}

	logger("Precomputation for round " + to_string(rounds) + " Completed.");
}

void filterMap(map<bitset<256>, int, CMPS<256>>& mp)
{
	int size0 = mp.size();
	map< bitset<256>, int, CMPS<256>> tmp(mp);
	mp.clear();

	for (auto & it : tmp)
	{
		if (it.second % 2 == 1)
			mp[it.first] = it.second;
	}
	int size1 = mp.size();

	logger(__func__ + string(": ") + to_string(size0) + string("\t") +
		to_string(size1));
}

void printSol(int rounds, const bitset<256>& vector, const map<bitset<128>,
	int, CMPS<128>> &solutions)
{
	logger(__func__ + string(": ") + to_string(rounds) + string("\t") +
		vector.to_string());

	string path = string("TERM/") + to_string(rounds) + string(".txt");
	ofstream os;
	os.open(path, ios::out | ios::app);
	os << rounds << "\t" << vector << endl;
	for (auto &it : solutions)
		os << it.first << "\t" << it.second << "\n";
	os << endl;
	os.close();
}

void resolveTerms(vector<bitset<128>>& kcofs, stringstream& ss)
{
	auto sep = "";
	ss << "(";
	for (auto& x : kcofs)
	{
		ss << sep;
		if (x.count() == 0)
			ss << "1";
		else
			for (int i = 0; i < 128; i++)
				if (x[i])
					ss << "k" << i;

		sep = "+";
	}
	ss << ")";
}






// this function is only called from SolutionFunc to get solutions from counterMap
void SecondSolutionFunc(int rounds,int midround, bitset<96>& cube, map<bitset<256>, int, CMPS<256> >& counterMap, map<bitset<128>,int,CMPS<128> > & solutions,vector<bitset<256>> & expandTerms,float time, int singlethread)
{
	solutions.clear();

	// first sort out midsolutions
	vector<bitset<256> > midsolutions;
	for (auto& x : counterMap)
	{
		if (x.second % 2 == 1)
			midsolutions.emplace_back(x.first);
		else
		{
			bitset<256> midstate = x.first;
			int midsols = 0;
			auto status = newMidSolutionCounter(midround, cube, midstate, midsols, time, singlethread);
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

	bitset<256> flag1;
	flag1.set();



	vector<bitset<256> > nonZeroTrackReps;
	grainTrack(cube, nonZeroTrackReps);


	vector<bitset<256> > varsTrackReps;
	grainVarsTrack(cube, varsTrackReps, nonZeroTrackReps);

	bitset<256> nonZeroFlag = (midround < nonZeroTrackReps.size()) ? nonZeroTrackReps[midround] : flag1;
	bitset<256> varsFlag = (midround < varsTrackReps.size()) ? varsTrackReps[midround] : flag1;

	stringstream ss;

	for (auto& x : midsolutions)
	{
		/*
		bitset<256> midcube;
		for (int i = 0; i < 256; i++)
			if (varsFlag[i])
				midcube[i] = x[i];

		map<bitset<128>, int, CMPS<128>> midcubesols;
		auto status = MidSolutionCounter(midround, cube, x, midcubesols, time, singlethread);
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
			continue;
		}
		vector<bitset<128>> midcubekcofs;
		for (auto& y : midcubesols)
			if (y.second % 2)
				midcubekcofs.emplace_back(y.first);


		vector<bitset<128>> midconskcofs = { bitset<128>(0) };
		for(int i =0;i<256;i++)
			if(x[i] && nonZeroFlag[i] && varsFlag[i] == 0)
				midconskcofs = vecMul<128>(midconskcofs, conskCofs[i]);

		vector<bitset<128>> midkcofs = vecMul<128>(midcubekcofs, midconskcofs);


		for (auto& y : midkcofs)
			solutions[y]++;
			*/

		bitset<256> midcube;
		for (int i = 0; i < 256; i++)
			if (varsFlag[i])
				midcube[i] = x[i];

		map<bitset<128>, int, CMPS<128>> midcubesols;
		auto status = MidSolutionCounter(midround, cube, x, midcubesols, time, singlethread);
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
			continue;
		}
		vector<bitset<128>> midcubekcofs;
		for (auto& y : midcubesols)
			if (y.second % 2)
				midcubekcofs.emplace_back(y.first);


		if (midcubekcofs.size() > 0)
		{
			resolveTerms(midcubekcofs, ss);

			for (int i = 0; i < 256; i++)
				if (nonZeroFlag[i] && varsFlag[i] == 0 && x[i])
				{
					ss << "*";
					resolveTerms(conskCofs[i], ss);
				}

			ss << endl;
		}

	}

	lock_guard<mutex> guard(stream_mutex);
	loggerStream(ss, rounds);
}


STATUS SolutionFunc(int rounds, int midround, bitset<96>& cube, bitset<256>& last, map<bitset<128>, int, CMPS<128>>& solutions, vector<bitset<256>>& expandTerms,float time, int singlethread)
{
	if (midround > rounds || rounds < 0 || midround < 0)
	{
		logger(__func__ + string(" : parameter error."));
		exit(-1);
	}

	map<bitset<256>, int, CMPS<256> > counterMap;
	auto status = newMidSolutionCounter(rounds, cube, last, counterMap, time, singlethread, midround);

	if (status == SOLUTION)
		SecondSolutionFunc(rounds,midround, cube, counterMap, solutions, expandTerms,time, singlethread);
	
	return status;
}

// below is threadworker
void SolutionSearcherWorker(bitset<256>& vec, int ROUND, bitset<96>& cube,
	vector<bitset<256>>& layerTerms, float time, int singlethread,int midround,map<bitset<128>,int,CMPS<128>> & totalSolutions, vector<bitset<256>>& expandTerms)
{
	map<bitset<128>, int, CMPS<128>> solutions;

	auto status = SolutionFunc(ROUND, midround, cube, vec, solutions, expandTerms,time, singlethread);

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
	}
}




int main()
{

	bitset<96> cube;
	for (int i = 0; i < 96; i++)
		cube[i] = 1;
	cube[30] = 0;
		


	int midround = 40;
	ThreadPool thread_pool{};

	int ROUND = 191;
	int step = 60;
	int expandStep = 1;
	float expandTime = 3600;

	float time = 0;

	vector<bitset<256>> terms;
	vector<bitset<256>> layerTerms;

	BackExpandPolynomial(step, terms);

	ROUND = ROUND - step;

	logger(__func__ + string("First expand"));

	int singlethread = 2;

	map<bitset<128>, int, CMPS<128>> totalSolutions;
	vector<bitset<256> > expandTerms;

	// expand
	newExpandFunc(cube, terms, ROUND, expandStep, thread_pool, expandTime, singlethread,midround,expandTerms);
	logger("Second expand complete.");


	while (true)
	{
		if (ROUND > 120)
			time = 60;
		else if (ROUND > 110)
			time = 120;
		else if (ROUND > 100)
			time = 180;
		else if (ROUND > 33)
			time = 360;
		else if (ROUND > 10)
			time = 360;
		else
			time = 0;

		

		layerTerms.clear();
		preCompute(midround, cube);
		vector<std::future<void>> futures;
		for (auto& it : terms)
			futures.emplace_back(thread_pool.Submit(SolutionSearcherWorker, it,
				ROUND, ref(cube), ref(layerTerms), time, singlethread, midround, ref(totalSolutions),ref(expandTerms)));
		for (auto& it : futures)
			it.get();

#ifndef _WIN32
		showProcessMemUsage();
#endif

		filterMap<128>(totalSolutions);
		filterExpandTerms<256>(expandTerms);

#ifndef _WIN32
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
			expandTerms.clear();
			logger("Current ROUND : " + to_string(ROUND));
			logger("Current midround : " + to_string(midround));
		}

		if (ROUND < 10)
		{
			logger("ROUND is to small.");
			exit(-1);
		}

		terms = layerTerms;


		if (ROUND < FASTBOUND)
			expandStep = 1;

		if (ROUND - expandStep < midround)
			expandStep = ROUND - midround;

		newExpandFunc(cube, terms, ROUND, expandStep, thread_pool, expandTime, singlethread, midround,expandTerms);

		

		// back up
		/*for (auto& x : totalSolutions)
			loggerRes(x.first.to_string(), ROUND);*/
		for (auto& x : expandTerms)
			loggerExpandTerms(x.to_string(), ROUND);

	}
	cout << "Success!" << endl;


	cout << "Output to file." << endl;

	// output result to file
	string path = string("TERM/") + string("res.txt");
	ofstream os;
	os.open(path, ios::out | ios::app);
	for (auto& it : totalSolutions)
	{
		bitset<128> term = it.first;
		if (term.count() == 0)
		{
			os << "1" << endl;
		}
		else
		{
			for (int i = 0; i < 128; i++)
				if (term[i])
					os << "k" << i;
			os << endl;
		}
	}
	os << endl;
	os.close();


	cout << "Output to file finished." << endl; 
}
