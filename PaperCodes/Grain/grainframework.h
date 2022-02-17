#ifndef __GRAIN_FRAMEWORK_H__
#define  __GRAIN_FRAMEWORK_H__

#define ROUNDSBOUND 25
#define SIZEBOUND 15000
#define FASTBOUND 45

#include<vector>
#include<bitset>
#include<map>
#include"node.h"
#include"thread_pool.h"

using namespace thread_pool;

void newExpandFunc(bitset<96>& cube, vector<bitset<256> >& states, int& start, int step, ThreadPool& threadpool, float time, int threads, int& midround, vector<bitset<256> >& expandTerms);
#endif

