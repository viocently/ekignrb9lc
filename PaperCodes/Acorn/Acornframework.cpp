#include<map>
#include<malloc.h>
#include<cmath>
#include<vector>
#include<bitset>
#include<string>
#include"Acorn.h"
#include"gurobi_c++.h" 
#include"node.h"
#include"thread_pool.h"
#include"Acorncallback.h"
#include"log.h"
#include"Acornframework.h"


using namespace std;
using namespace thread_pool;
mutex threadMutex;
const int THREAD_NUM = 64;

// divide and conquer back expand thread worker
void ThreadWorkerBackExpandDivide(const vector<bitset<549>>& layerTerm, 
	int start, int end, int current, int step,
	map<bitset<549>, int, CMPS<549>>& layerMap, int singlethread)
{
	map<bitset<549>, int, CMPS<549>>  threadLayerMap;

	int END = end <= layerTerm.size() ? end : layerTerm.size();

	vector<bitset<549>> terms;


	for (int i = start; i < END; i++)
	{
		terms.clear();
		SecondBackExpandPolynomial(current,step,layerTerm[i], terms, singlethread);

		for (auto& it : terms)
			threadLayerMap[it] += 1;
	}

	lock_guard<mutex> guard(threadMutex);
	for (auto& it : threadLayerMap)
		layerMap[it.first] += it.second;
}

void BackExpandDivideFunc(vector<bitset<549> >& layerTerms, int& start, float time, int singlethread)
{
	vector<thread> threads;

	int step = 0;
	map<bitset<549>, int, CMPS<549> > layerMap;

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
			threads.push_back(thread(ThreadWorkerBackExpandDivide, ref(
				layerTerms),  base * i, base * (i + 1), start,step,
				ref(layerMap), singlethread));
		for (auto& th : threads) th.join();

		if (size > size_multiple)
		{
			threads.clear();
			for (int i = 0; i < THREAD_NUM; i++)
				threads.push_back(thread(ThreadWorkerBackExpandDivide, ref(
					layerTerms),  size_multiple + i, size_multiple + i + 1, start,step,
					ref(layerMap), singlethread));
			for (auto& th : threads) th.join();
		}
	}

	start = start - step;

	logger("Step New " + to_string(step));
	logger("Round New " + to_string(start));
	logger(string("layerMapSize ") + to_string(layerMap.size()));

	layerTerms.clear();
	filterMap<549>(layerMap);

	for (auto& x : layerMap)
		layerTerms.emplace_back(x.first);

	for (auto& x : layerTerms)
		loggerState(x.to_string(), start);
}


//  back expand thread worker using callback functions
void ThreadWorkerBackExpandCallback(int rounds, int midround, bitset<293>& startcube, bitset<549>& last, map<bitset<549>, int, CMPS<549> >& resMap, float time, int threads)
{
	vector<bitset<549>> res;
	auto status = uMidSolutionCallBack(rounds, startcube, last, res, time, threads, midround, rounds - 30 > 350 ? rounds - 30 : 350);

	if (status == EXPAND)
	{
		logger("timeout" + string(" : need expand."));
		logger("timeout : " + last.to_string());
		res.clear();
		SecondBackExpandPolynomial(rounds, rounds - midround, last, res, threads);
	}

	{
		lock_guard<mutex> guard(threadMutex);
		for (auto& x : res)
			resMap[x] ++;
	}
}

void BackExpandCallbackFunc(bitset<293>& startcube, vector<bitset<549> >& states, int& start, int step, ThreadPool& threadpool, float time, int threads)
{
	if (step > start)
		return;

	int midround = start - step;
	int execCnt = 0;
	while ((states.size() < SIZEBOUND && start > FASTBOUND) || (start <= FASTBOUND && execCnt == 0))
	{
		map<bitset<549>, int, CMPS<549> > resMap;
		vector<std::future<void>> futures;

		for (auto& x : states)
			futures.emplace_back(threadpool.Submit(ThreadWorkerBackExpandCallback, start, midround, ref(startcube), ref(x), ref(resMap), time, threads));

		for (auto& x : futures)
			x.get();

		logger("from " + to_string(start) + " to " + to_string(midround) + " complete.");

		// filterMap
		filterMap<549>(resMap);

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

		logger("next states Size : " + to_string(states.size()));
		execCnt++;
	}
}