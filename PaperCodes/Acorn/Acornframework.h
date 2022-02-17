#ifndef __ACORN_FRAMEWORK_H__
#define  __ACORN_FRAMEWORK_H__


#define FASTBOUND 420
#define SIZEBOUND 15000
#define CALLBACKBOUND 540

#include<vector>
#include<bitset>
#include<map>
#include"node.h"
#include"thread_pool.h"

using namespace thread_pool;
void BackExpandDivideFunc(vector<bitset<549> >& layerTerms, int& start, float time, int singlethread);
void BackExpandCallbackFunc(bitset<293>& startcube, vector<bitset<549> >& states, int& start, int step, ThreadPool& threadpool, float time, int threads);
#endif