#include<map>
#include<cmath>
#include<vector>
#include<bitset>
#include<string>
#include<malloc.h>
#include"grain.h"
#include"gurobi_c++.h" 
#include"node.h"
#include"thread_pool.h"
#include"graincallback.h"
#include"log.h"
#include"grainframework.h"


using namespace std;
using namespace thread_pool;
mutex threadMutex;

// second Expand
void newThreadWorker(int rounds, int midround, bitset<96>& cube, bitset<256>& last, map<bitset<256>, int, CMPS<256> >& resMap, float time, int threads)
{
	vector<bitset<256>> res;
	//auto status = newMidSolutionCallBack(rounds, cube, last, midround, res, time, threads);
	int coRound = midround;
	auto status = uMidSolutionCallBack(rounds, cube, last, midround, res, time, threads,coRound);

	if (status == EXPAND)
	{
		logger("timeout" + string(" : need expand."));
		logger("timeout : " + last.to_string());
		res.clear();
		SecondBackExpandPolynomial(rounds - midround, last, res, threads);
	}

	
	{
		lock_guard<mutex> guard(threadMutex);
		for (auto& x : res)
			resMap[x] ++;
	}
}

void newExpandFunc(bitset<96>& cube, vector<bitset<256> >& states, int& start, int step, ThreadPool& threadpool, float time, int threads,int &midround, vector<bitset<256> > & expandTerms)
{
	if (step > start)
		return;
	int nextround = start - step;
	int execCnt = 0;
	while ((states.size() < SIZEBOUND && start > ROUNDSBOUND) || (start <= ROUNDSBOUND && execCnt == 0))
	{

		map<bitset<256>, int, CMPS<256> > resMap;
		vector<std::future<void>> futures;
		for (auto& x : states)
			futures.emplace_back(threadpool.Submit(newThreadWorker, start, nextround, ref(cube), ref(x), ref(resMap), time, threads));
		for (auto& x : futures)
			x.get();

		logger("from " + to_string(start) + " to " + to_string(nextround) + " complete.");

		// if ROUND == midround,combine expandterms with terms since they are all in midround
		if (nextround == midround)
		{
			logger("current ROUND : " + to_string(nextround));
			logger("start combine expandTerms with terms.");
			for (auto& x : expandTerms)
				resMap[x]++;
			expandTerms.clear();
			midround = midround / 2;
		}

		// filterMap
		filterMap(resMap);

#ifndef _WIN32
		showProcessMemUsage();
		malloc_trim(0);
		showProcessMemUsage();
#endif

		// save resMap
		states.clear();
		for (auto& x : resMap)
		{
			loggerState(x.first.to_string(), nextround);
			states.emplace_back(x.first);
		}

		// update start and midround and states
		start = nextround;
		if(start - step >= midround)
			nextround = start - step;
		else
			nextround = midround;

		logger("next states Size : " + to_string(states.size()));
		execCnt++;
	}
}
