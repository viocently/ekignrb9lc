#ifndef __TRIVIUM_FRAMEWORK_H__
#define  __TRIVIUM_FRAMEWORK_H__

// make sure CALLBACKBOUND is bigger than ROUDNSBOUND
#define FASTBOUND 220
#define SIZEBOUND 15000
#define CALLBACKBOUND 600

#include<vector>
#include<bitset>
#include<map>
#include"node.h"
#include"thread_pool.h"

using namespace thread_pool;
void expandFunc(bitset<80>& cube, vector<bitset<288> >& layerTerms, int& start, float time, int singlethread);
void newExpandFunc(bitset<80>& cube, vector<bitset<288> >& states, int& start, int step, ThreadPool& threadpool, float time, int threads);

#endif