#ifndef __KREYVIUM_FRAMEWORK_H__
#define  __KREYVIUM_FRAMEWORK_H__

#define FASTBOUND 240
#define SIZEBOUND 10000
#define CALLBACKBOUND 420

#include<vector>
#include<bitset>
#include<map>
#include"node.h"
#include"thread_pool.h"

using namespace thread_pool;
void expandFunc(bitset<128>& cube, vector<bitset<256 + 288> >& layerTerms, int& start, float time, int singlethread);
void newExpandFunc(bitset<128>& cube, vector<bitset<256 + 288> >& states, int& start, int step, ThreadPool& threadpool, float time, int threads);

#endif
