#include<malloc.h>
#include<map>
#include<cmath>
#include<vector>
#include<bitset>
#include<string>
#include"trivium.h"
#include"gurobi_c++.h" 
#include"node.h"
#include"thread_pool.h"
#include"triviumcallback.h"
#include"log.h"
#include"triviumframework.h"


using namespace std;
using namespace thread_pool;
mutex threadMutex;
const int THREAD_NUM = 64;

// second Expand using divde-and-conquer
void ThreadWorkerExpand(const vector<bitset<288>>& layerTerm, const bitset<80>& cube,
	int start, int end, int step,
	map<bitset<288>, int, CMPS<288>>& layerMap, int singlethread)
{
	//map<bitset<80>, int, CMPS<80>> solutions;
	map<bitset<288>, int, CMPS<288>>  threadLayerMap;

	int END = end <= layerTerm.size() ? end : layerTerm.size();

	vector<bitset<288>> terms;


	for (int i = start; i < END; i++)
	{
		terms.clear();
		SecondBackExpandPolynomial(step, layerTerm[i], terms, singlethread);


		for (auto& it : terms)
			threadLayerMap[it] += 1;
	}

	lock_guard<mutex> guard(threadMutex);
	for (auto& it : threadLayerMap)
		layerMap[it.first] += it.second;
}

void expandFunc(bitset<80>& cube, vector<bitset<288> >& layerTerms, int& start, float time, int singlethread)
{
	vector<thread> threads;

	int step = 0;
	map<bitset<288>, int, CMPS<288> > layerMap;

	while (layerMap.size() < SIZEBOUND)
	{
		step++;
		layerMap.clear();

		int size = layerTerms.size();

		int base = (size / THREAD_NUM) >= 1 ? (size / THREAD_NUM) :
			1;
		int size_multiple = base * THREAD_NUM;

		threads.clear();

		// expand the terms 
		for (int i = 0; i < THREAD_NUM; i++)
			threads.push_back(thread(ThreadWorkerExpand, ref(
				layerTerms), ref(cube), base * i, base * (i + 1), step,
				ref(layerMap), singlethread));
		for (auto& th : threads) th.join();

		if (size > size_multiple)
		{
			threads.clear();
			for (int i = 0; i < THREAD_NUM; i++)
				threads.push_back(thread(ThreadWorkerExpand, ref(
					layerTerms), ref(cube), size_multiple + i, size_multiple + i + 1, step,
					ref(layerMap), singlethread));
			for (auto& th : threads) th.join();
		}
	}

	start = start - step;

	logger("Step New " + to_string(step));
	logger("Round New " + to_string(start));
	logger(string("layerMapSize ") + to_string(layerMap.size()));

	layerTerms.clear();
	filterMap<288>(layerMap);

	for (auto& x : layerMap)
		layerTerms.emplace_back(x.first);

	filterTerms(cube, start, layerTerms);

	for (auto& x : layerTerms)
		loggerState(x.to_string(), start);
}





// second Expand using callback function
void newThreadWorkerCallBack(int rounds, int midround, bitset<80>& cube, bitset<288>& last, map<bitset<288>, int, CMPS<288> >& resMap, float time, int threads)
{
	vector<bitset<288>> res;
	auto status = uMidSolutionCallBack(rounds, cube, last, res, time, threads, midround);

	if (status == EXPAND)
	{
		logger("timeout" + string(" : need expand."));
		res.clear();
		SecondBackExpandPolynomial(rounds - midround, last, res, threads);
	}

	{
		lock_guard<mutex> guard(threadMutex);
		for (auto& x : res)
			resMap[x] ++;
	}
}

void newExpandFunc(bitset<80>& cube, vector<bitset<288> >& states, int& start, int step, ThreadPool& threadpool, float time, int threads)
{
	if (step > start)
		return;

	int midround = start - step;
	int execCnt = 0;
	while ((states.size() < SIZEBOUND && start > FASTBOUND) || (start <= FASTBOUND && execCnt == 0))
	{
		map<bitset<288>, int, CMPS<288> > resMap;
		vector<std::future<void>> futures;

		for (auto& x : states)
			futures.emplace_back(threadpool.Submit(newThreadWorkerCallBack, start, midround, ref(cube), ref(x), ref(resMap), time, threads));

		for (auto& x : futures)
			x.get();

		logger("from " + to_string(start) + " to " + to_string(midround) + " complete.");

		// filterMap
		filterMap<288>(resMap);

#ifdef _WIN32
#else
		showProcessMemUsage();
		malloc_trim(0);
		showProcessMemUsage();
#endif

		// save resMap
		states.clear();
		for (auto& x : resMap)
		{
			loggerState(x.first.to_string(), midround);
			states.emplace_back(x.first);
		}

		// update start and midround and states
		start = midround;
		midround = start - step;

		filterTerms(cube, start, states);

		logger("next states Size : " + to_string(states.size()));
		execCnt++;
	}
}

